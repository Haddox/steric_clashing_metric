"""
This script contains various modules for protein design

Hugh Haddox, January 29, 2018
"""

import os
import re
import random
import math
import time
import numpy
import pandas
import subprocess
import glob
import Bio.PDB
from Bio.Seq import Seq
import gzip
import bz2


def WriteSbatchFile(
    sbatch_file_name, command_file_name=None, command=None, queue_type='short',
    memory='2g', array=None, group_size=None, n_commands=None
    ):
    """
    Write a file to submit a job via `sbatch`

    Args:
        `sbatch_file_name`: The name of the sbatch file. This will also
            serve as a prefix for the name of the output and error files.
        `command_file_name`: The name of a file with a list of commands
            to execute (string). The default is `None`, but if one is given,
            then this function will write an `sbatch` file that is specific
            for carying out an array of commands. If this is not given,
            the `command` argument must be specified. Only one can be
            specified at once.
        `command`: A command-line argument to be carried out (string).
            The default is `None`. If this is not given, the
            `command_file_name` argument must be specified. Only one can
            be specified at once.
        `queue_type`: The queue type ('short', 'medium', 'long')
        `array`: The size of the array (default: None)
        `memory`: The amount of memory in megabytes, unless other unit is specified. Default is '2g', i.e., 2 GB.
        `group_size`: The number of commands per group in a group run
        `n_commands`: The total number of commands. When `group_size` is
            specified, this number will be used to compute the size of the
            array that will be called.

    Retruns:
        A file for submitting a job via `sbatch`
    """

    # If a file with an array of commands is provided, write an `sbatch`
    # file that is suited for this task
    if command_file_name:

        # Make sure the `command` argument hasn't been called
        assert command==None

        # Write the file
        with open(sbatch_file_name, 'w') as f:

            # ... to execute the commands in groups
            if group_size:
                array = math.ceil(n_commands/group_size)
                f.write('#!/bin/bash\n')
                f.write('#SBATCH -p {0}\n'.format(queue_type))
                f.write('#SBATCH -a 1-{0}\n'.format(array)) # total tasks/$GROUP_SIZE
                f.write('#SBATCH --mem={0}\n'.format(memory))
                f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
                f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
                f.write('GROUP_SIZE={0}\n'.format(group_size))
                f.write('for I in $(seq 1 $GROUP_SIZE)\n')
                f.write('do\n')
                f.write('\tJ=$(($SLURM_ARRAY_TASK_ID * $GROUP_SIZE + $I - $GROUP_SIZE))\n')
                f.write('\tCMD=$(sed -n "${J}p" %s)\n'%command_file_name)
                f.write('\t${CMD}\n')
                f.write('done\n')

            # ... to execute the commands individually
            else:
                f.write('#!/bin/bash\n')
                f.write('#SBATCH -p {0}\n'.format(queue_type))
                if array:
                    f.write('#SBATCH -a {0}\n'.format(array))
                f.write('#SBATCH --mem={0}\n'.format(memory))
                f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
                f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
                f.write('CMD=$(head -n $SLURM_ARRAY_TASK_ID {0} | tail -1)\n'.format(command_file_name))
                f.write('exec ${CMD}')

    # Write an `sbatch` file to carry out the specified `command`
    else:
        with open(sbatch_file_name, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -p {0}\n'.format(queue_type))
            if array:
                f.write('#SBATCH -a {0}\n'.format(array))
            f.write('#SBATCH --mem={0}\n'.format(memory))
            f.write('#SBATCH -o {0}.out\n'.format(sbatch_file_name))
            f.write('#SBATCH -e {0}.err\n'.format(sbatch_file_name))
            f.write(command)


def parse_hbnet_metadata_from_pdb(pdb_path, verbose=False):
    """
    Parse residues in an HBNet given an input PDB file

    Args:
        `pdb_path`: the path to a PDB file with a hydrogen-bond network made
            using `HBNet`, as well as information about that network in the PDB
            file, requires adding the flag `-out::file::pdb_comments` to the
            design run

    Returns:
        `hbnet_dict`: a dictionary with the following keys:
            `hbnet_network_residues`: a list of tupples, where each tupple is
                made up of two objects: i) a string giving the amino-acid
                identity using single-letter abbreviations and ii) an integer
                giving the position of the amino acid, with positions indexed
                starting at one relative to the input sequence. Alternatively, a
                tupple may be replaced by the string `backbone` if the network
                includes a backbone residue
            `hbnet_size`: the size of the HBNet
            `hbnet_score`: the score of the network according to HBNet
            `hbnet_num_hbonds`: number of hydrogen bonds in the network
            `hbnet_percent_hbond_capacity`: same as in HBNet
            `hbnet_num_unsat_Hpol`: the number of unsatisfied polar hydrogen atoms
                in the network
    """

    # Compile a regular-expressions pattern to parse residues in the network
    net_res_pattern = re.compile(r'(?P<chain>\w)_(?P<aa>\w)_(?P<res_n>\d+)')

    # Find the PDB file of interest
    if pdb_path[-4:] != '.pdb':
        pdb_path += '.pdb'
    starting_backbone_prefix = pdb_path[:]
    assert os.path.isfile(pdb_path), "Could not find the PDB file: {0}".format(pdb_path)

    # Parse the PDB file for which residues are in the HBNet
    network_residues = []
    with open(pdb_path) as f:
        found_network = False
        for line in f:
            if '#network' not in line:
                continue
            line_elements = line.strip().split('\t')
            assert len(line_elements) == 7
            net_res = line_elements[1].split(',')
            size = float(line_elements[2])
            score = float(line_elements[3])
            num_hbonds = float(line_elements[4])
            percent_hbond_capacity = float(line_elements[5])
            num_unsat_Hpol = float(line_elements[6])
            for res_i in net_res:
                match = re.search(net_res_pattern, res_i)
                if match:
                    (chain, aa, res_n) = match.group('chain', 'aa', 'res_n')
                    network_residues.append((aa, int(res_n)))
                elif res_i in ('backbone', 'backbonecycle_'):
                    network_residues.append(res_i)
                else:
                    raise ValueError("Could not parse the entry: {0}".format(res_i))
            found_network = True
            break
    if not found_network:
        raise ValueError("Could not find the expected network information in the PDB file: {0}".format(pdb_path))

    # Store the information in a dictionary and then return the dictionary
    hbnet_dict = {
        'hbnet_network_residues': network_residues,
        'hbnet_size': size,
        'hbnet_score': score,
        'hbnet_num_hbonds': num_hbonds,
        'hbnet_percent_hbond_capacity': percent_hbond_capacity,
        'hbnet_num_unsat_Hpol': num_unsat_Hpol
    }
    return (hbnet_dict)


def check_residues_are_intact(residues, sequence, verbose=False):
    """
    Make sure an input list of residues are present in an input sequence

    Args:
        `residues`: a list of tupples, where each tupple is made up of
            two objects: i) a string giving the amino-acid identity using
            single-letter abbreviations and ii) an integer giving the
            position of the amino acid, with positions indexed starting
            at one relative to the input sequence. Alternatively, a tupple
            may be replaced by the string `backbone`, in which case this
            residue is considered to be intact by default.
        `sequence`: a string giving the amino-acid sequence

    Returns:
        `True` if each of the input residues is present in the input sequence,
            both in terms of amino-acid identity and position. `False` if
            at least one residue is NOT intact. `None` if the `residues`
            argument is also `None`.

    Code for `doctest`:
    >>> residues = [('A', 1), ('C', 3), 'backbone']
    >>> sequence = 'ARCQHYP'
    >>> check_residues_are_intact(residues, sequence)
    True
    >>> residues = [('A', 1), ('C', 4)]
    >>> check_residues_are_intact(residues, sequence)
    False
    """

    # In a trivial case where `None` is passed as the `residues` variable,
    # return `None`
    if residues == None:
        return None

    # Otherwise, examine the residues
    residues_intact = []
    for residue in residues:
        if residue in ['backbone', 'backbonecycle_']:
            residues_intact.append(True)
        else:
            (aa, res_n) = residue
            if not isinstance(res_n, int):
                raise ValueError("The following residue number is not an integer: {0}".format(res_n))
            if sequence[res_n - 1].upper() == aa.upper():
                residues_intact.append(True)
            else:
                residues_intact.append(False)
                if verbose:
                    raise ValueError("Residue number {0} is not the amino acid {1} in the sequence:\n{2}".format(
                        res_n, aa, sequence
                    ))
    if sum(residues_intact) == len(residues):
        return True
    else:
        return False


def extract_total_score_from_pdb_file(pdb_file):
    """
    Extract the total score from the Rosetta score table in an input PDB file

    Args:
        `pdb_file`: a path to a PDB file with a Rosetta score table (string)

    Returns:
        The total score as computed by Rosetta(float)
    """
    with open(pdb_file) as f:
        in_energy_table = False
        for line in f:
            if '#BEGIN_POSE_ENERGIES_TABLE' in line:
                in_energy_table = True
            if not in_energy_table:
                continue
            elements = line.strip().split(' ')
            if elements[0] == 'label':
                n_labels = len(elements)
            elif elements[0] == 'pose':
                assert len(elements) == n_labels
                total_score = float(elements[-1])
                break
    return total_score


def compute_score_per_res(
    rosettapath, xml_file, file_listing_pdbs, weights_file, output_dir,
    scores_file_prefix, extra_script_vars=None, extra_args=None,
    flags_file=None, submit_sbatch_job=True, queue_type='medium', memory='2g'
    ):
    """
    Compute the `score_per_res` for input proteins using a given energy function

    Args:
        `rosettapath`: a path to Rosetta (string)
        `xml_file`: a path to an input XML file for scoring (string)
        `file_listing_pdbs`: a path to a file giving paths of PDB files to analyze (string)
        `weights_file`: a path to an input weights file or a string that specifies the
            weights in a way that is recognized by Rosetta
        `output_dir`: the path to the directory in which to store output designs
        `scores_file_prefix`: the prefix of the name of the out scores file in the `output_dir`
        `extra_script_vars`: a list of strings that will be added to the command-line
            argument following `-parser:script_vars` flag, allowing additional variables
            to be specified here and passed to the XML. Default value is `None`.
        `extra_args`: a list of strings that will be added as command-line arguments. Default
            value is `None`.
        `flags_file`: a path to an input file with Rosetta flags. Default value is `None`.
            But, if a file is provided, this file will be referenced in the command-line
            argument for making designs.
        `submit_sbatch_job`: a boolean specifying whether to submit the sbatch job
            (default=`True`)
        `queue_type`: the sbatch queue to use (string; default: medium)
        `memory`: the amount of memory to use in the sbatch job submission

    Returns:
        A scores file in `output_dir`
    """

    # Start putting together the command-line argument
    cmd = ' '.join([
        rosettapath,
        '-l {0}'.format(file_listing_pdbs),
        '-parser:protocol {0}'.format(xml_file),
        '-out:prefix {0}'.format(output_dir),
        '-out:file:scorefile {0}.sc'.format(scores_file_prefix),
        '-parser:script_vars',
        'weights_file={0}'.format(weights_file)
    ])

    # Add extra `script_vars` and aguments if appropriate, including ones that
    # are provided or a file prefix to use for dumping PDBs
    if extra_script_vars:
        for arg in extra_script_vars:
            cmd += ' {0}'.format(arg)
    if extra_args:
        for arg in extra_args:
            cmd += ' {0}'.format(arg)

    # If `flags_file` is provided, add it to the command
    if flags_file:
        cmd += ' @{0}'.format(flags_file)

    # Write an sbatch file to carry out the command
    sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)

    WriteSbatchFile(
        sbatch_file_name,
        command=cmd,
        queue_type=queue_type,
        memory=memory
    )
    if submit_sbatch_job:
        subprocess.check_call(['sbatch', sbatch_file_name])


def compute_score_using_rosettascripts(
    score_app_path, weights_file, output_dir,
    scores_file_prefix, file_listing_pdbs=False, silent_file=False,
    extra_args=None, flags_file=None,
    submit_sbatch_job=True, queue_type='medium', memory='2g'
    ):
    """
    Compute the `score_per_res` for input proteins using a given energy function

    Args:
        `score_app_path`: a path to the "score" application in Rosetta (string)
        `file_listing_pdbs`: a path to a file giving paths of PDB files to analyze (string)
        `weights_file`: a path to an input weights file or a string that specifies the
            weights in a way that is recognized by Rosetta
        `output_dir`: the path to the directory in which to store output designs
        `scores_file_prefix`: the prefix of the name of the out scores file in the `output_dir`
        `extra_args`: a list of strings that will be added as command-line arguments. Default
            value is `None`.
        `flags_file`: a path to an input file with Rosetta flags. Default value is `None`.
            But, if a file is provided, this file will be referenced in the command-line
            argument for making designs.
        `submit_sbatch_job`: a boolean specifying whether to submit the sbatch job
            (default=`True`)
        `queue_type`: the sbatch queue to use (string; default: medium)
        `memory`: the amount of memory to use in the sbatch job submission

    Returns:
        A scores file in `output_dir`
    """
    
    # Start putting together the command-line argument
    if file_listing_pdbs:
        infile_arg = '-l {0}'.format(file_listing_pdbs)
    elif silent_file:
        infile_arg = f'-in:file:silent {silent_file}'
    else:
        raise ValueError("Must specify `file_listing_pdbs` or `silent_file`")
    
    cmd = ' '.join([
        score_app_path,
        infile_arg,
        '-out:prefix {0}'.format(output_dir),
        '-out:file:scorefile {0}.sc'.format(scores_file_prefix),
        '-score:weights {0}'.format(weights_file)
    ])

    # Add extra aguments if appropriate
    if extra_args:
        for arg in extra_args:
            cmd += ' {0}'.format(arg)

    # If `flags_file` is provided, add it to the command
    if flags_file:
        cmd += ' @{0}'.format(flags_file)

    # Write an sbatch file to carry out the command
    sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)

    WriteSbatchFile(
        sbatch_file_name,
        command=cmd,
        queue_type=queue_type,
        memory=memory
    )
    if submit_sbatch_job:
        subprocess.check_call(['sbatch', sbatch_file_name])


def compute_score_per_res_after_relax(
    rosettapath, databasepath, xml_file, relaxscript, pdbs, weights_file,
    output_dir, scores_file_prefix, extra_script_vars=None, extra_args=None,
    flags_file=None, submit_sbatch_job=True, queue_type='medium', memory='2g',
    sbatch_group_size=None, subsample_pdbs=False, n_repeats=1
    ):
    """
    Compute the `score_per_res` for input proteins using a given energy function

    Args:
        `rosettapath`: a path to Rosetta (string)
        `databasepath`: a path to a Rosetta database (string)
        `xml_file`: a path to an input XML file for scoring (string)
        `relaxscript`: a path to a relax script to use (string)
        `pdbs`: a list of paths of PDB files to analyze
        `weights_file`: a path to an input weights file or a string that specifies the
            weights in a way that is recognized by Rosetta
        `output_dir`: the path to the directory in which to store output designs
        `scores_file_prefix`: the prefix of the name of the out scores file in the `output_dir`
        `extra_script_vars`: a list of strings that will be added to the command-line
            argument following `-parser:script_vars` flag, allowing additional variables
            to be specified here and passed to the XML. Default value is `None`.
        `extra_args`: a list of strings that will be added as command-line arguments. Default
            value is `None`.
        `flags_file`: a path to an input file with Rosetta flags. Default value is `None`.
            But, if a file is provided, this file will be referenced in the command-line
            argument for making designs.
        `submit_sbatch_job`: a boolean specifying whether to submit the sbatch job
            (default=`True`)
        `queue_type`: the sbatch queue to use (string; default: medium)
        `memory`: the amount of memory to use in the sbatch job submission
        `sbatch_group_size`: the number of commands to group together when computing
            scores.
        `subsample_pdbs`: default is "False". If an integer is provided, it is the number of
            PDBs that will be subsampled in a pseudo-random fashion from the list of input
            PDBs.
        `n_repeats`: number of `FastRelax` repeats (integer; default: 1)

    Returns:
        A scores file in `output_dir`
    """

    # Subsample PDBs if specified
    if subsample_pdbs:
        random.seed(a=3)
        pdbs = random.sample(pdbs, subsample_pdbs)

    # Start putting together the command-line argument
    cmds = []
    for pdb in pdbs:
        cmd = ' '.join([
            rosettapath,
            '-database {0}'.format(databasepath),
            '-in:file:s {0}'.format(pdb),
            '-parser:protocol {0}'.format(xml_file),
            '-out:prefix {0}'.format(output_dir),
            '-out:file:scorefile {0}.sc'.format(scores_file_prefix),
            '-parser:script_vars',
            'weights_file={0}'.format(weights_file),
            'relaxscript={0}'.format(relaxscript),
            'n_repeats={0}'.format(n_repeats)
        ])

        # Add extra `script_vars` and aguments if appropriate, including ones that
        # are provided or a file prefix to use for dumping PDBs
        if extra_script_vars:
            for arg in extra_script_vars:
                cmd += ' {0}'.format(arg)
        if extra_args:
            for arg in extra_args:
                cmd += ' {0}'.format(arg)

        # If `flags_file` is provided, add it to the command
        if flags_file:
            cmd += ' @{0}'.format(flags_file)

        cmds.append(cmd)

    # Write the commands to a file
    command_file = '{0}_commands.txt'.format(scores_file_prefix)
    with open(command_file, 'w') as f:
        for cmd in cmds:
            f.write('{0}\n'.format(cmd))

    # Write a file that uses sbatch to execute the command in the above file
    sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)
    if not sbatch_group_size:
        WriteSbatchFile(
            sbatch_file_name,
            command_file_name=command_file,
            queue_type=queue_type,
            memory=memory
        )
        if submit_sbatch_job:
            subprocess.check_call(['sbatch', '-a', '1-{0}'.format(len(cmds)), sbatch_file_name])

    else:
        WriteSbatchFile(
            sbatch_file_name,
            command_file_name=command_file,
            queue_type=queue_type,
            memory=memory,
            group_size=sbatch_group_size,
            n_commands=len(cmds)
        )
        if submit_sbatch_job:
            subprocess.check_call(['sbatch', sbatch_file_name])


def forward_fold_design(
    pdb_path, forward_folding_dir, minirosetta_path,
    e_function_name, weights_file, flags_file, extra_args=[],
    run_jobs=False, n_abinitio_trajectories=9, n_abinitio_parallel_pools=3,
    abinitio_starting_pool_index=0, relax_starting_pool_index=0,
    n_relax_trajectories=10, n_relax_parallel_pools=1, verbose=False,
    queue_type='medium'
    ):
    """
    Generate command files to forward fold a sequence from an input PDB file

    Args:
        `pdb_path`: the path to a PDB file (string)
        `forward_folding_dir`: the directory in which to store the results
        `n_trajectories`


        `n_abinitio_trajectories`: the number of ab initio jobs to run (int)
        `n_relax_trajectories`: the number of relax jobs to run (int)


    Returns:
        - A directory with fragments for use in forward folding
    """

    # Make a directory in which to store designs of forward folding and
    # subsequent results if the directory doesn't already exist
    if not os.path.isdir(forward_folding_dir):
        os.makedirs(forward_folding_dir)
        if verbose:
            print("Making the following directory for storing all results of forward foldign: {0}".format(forward_folding_dir))

    # Convert paths to absolute paths to enable the below scripts to be carried
    # out in directories that differ from where the PDBs are stored
    cwd = os.getcwd()
    pdb_path = os.path.join(cwd, pdb_path)

    #----------------------------------------------------------
    # Make fragments for forwardward folding, unless a fragment folder already
    # exists
    #----------------------------------------------------------
    pdb_basename = os.path.basename(pdb_path)
    fragment_dir = os.path.join(
        forward_folding_dir, '{0}_fragments'.format(pdb_basename.replace('.pdb', ''))
    )
    if not os.path.isdir(fragment_dir):
        if verbose:
            print("Making fragments")
        cmd = [
            'python',
            '/home/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/make_fragments.py',
            '-pdbs', pdb_path
        ]
        make_fragments_outfile = os.path.join(forward_folding_dir, 'make_fragments.out')
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=forward_folding_dir)
        (out, err) = process.communicate()
        with open(make_fragments_outfile, 'wb') as f:
            f.write(out)
        if err:
            print(err)
            raise ValueError('Error in making fragments. The error message is:\n{0}'.format(err))

    #----------------------------------------------------------
    # Run `AbinitioRelax` using `minirosetta` using the exact same command-line
    # arguments as the standard protocol for doing this on boinc
    #----------------------------------------------------------
    # Define input variables
    frag3_file = os.path.join(fragment_dir, '00001.200.3mers')
    frag9_file = os.path.join(fragment_dir, '00001.200.9mers')
    assert os.path.isfile(frag3_file)
    assert os.path.isfile(frag9_file)

    # Determine how many parallel jobs to submit and how many structures to
    # generate per job as a function of the relevant input variables
    if (n_abinitio_trajectories % n_abinitio_parallel_pools) != 0:
        raise ValueError("The number of abinitio trajectories ({0}) is not divisible by the number of jobs ({1})".format(
            n_abinitio_trajectories, n_abinitio_parallel_pools
        ))
    abinitio_nstruct = int(n_abinitio_trajectories / n_abinitio_parallel_pools)

    # Create a directory in which to store the abinitio results and then write
    # all commands to a file
    abinitio_resultsdir = os.path.join(fragment_dir, e_function_name, 'abinitio/')
    if not os.path.isdir(abinitio_resultsdir):
        if verbose:
            print("Commencing ab initio runs")
        os.makedirs(abinitio_resultsdir)
        abinitio_commands_file = os.path.join(abinitio_resultsdir, 'commands.txt')
        with open(abinitio_commands_file, 'w') as f:
            for pool_i in list(range(abinitio_starting_pool_index, abinitio_starting_pool_index + n_abinitio_parallel_pools)):
                abinitio_silent_file = os.path.join(
                    abinitio_resultsdir, 'pool_{0}.out'.format(pool_i)
                )
                cmd = ' '.join([
                    minirosetta_path,
                    '-in:file:native', pdb_path,
                    '-frag3', frag3_file,
                    '-frag9', frag9_file,
                    '-nstruct', str(abinitio_nstruct),
                    '-out:file:silent', abinitio_silent_file,
                    '-silent_gz 1',
                    '-score:weights', weights_file,
                    '-abinitio::increase_cycles 10',
                    '-abinitio::fastrelax 1',
                    '-abinitio::rsd_wt_helix 0.5',
                    '-abinitio::rsd_wt_loop 0.5',
                    '-abinitio::use_filters false',
                    '-abinitio::rg_reweight 0.5',
                    '-relax::default_repeats 5',
                    '-ex1 1',
                    '-ex2aro 1',
                    '-mute all'
                ])
                if extra_args:
                    for arg in extra_args:
                        cmd += ' {0}'.format(arg)
                if flags_file:
                    cmd += ' @{0}'.format(flags_file)
                f.write('{0}\n'.format(cmd))

        # Write an `sbatch` input file that can be used to carry out the commands
        # written to the above file
        abinitio_sbatch_file = os.path.join(abinitio_resultsdir, 'run.sbatch')
        queue_type = queue_type
        WriteSbatchFile(
            abinitio_sbatch_file,
            command_file_name=abinitio_commands_file,
            queue_type=queue_type
        )
    else:
        abinitio_commands_file = None
        abinitio_sbatch_file = None

    #----------------------------------------------------------
    # Run `relax` using `minirosetta` using the exact same command-line
    # arguments as the standard protocol for doing this on boinc
    #----------------------------------------------------------

    # Determine how many parallel jobs to submit and how many structures to
    # generate per job as a function of the relevant input variables
    if (n_relax_trajectories % n_relax_parallel_pools) != 0:
        raise ValueError("The number of relax trajectories ({0}) is not divisible by the number of jobs ({1})".format(
            n_relax_trajectories, n_relax_parallel_pools
        ))
    relax_nstruct = int(n_relax_trajectories / n_relax_parallel_pools)

    # Create a directory in which to store the relax results and then write
    # all commands to a file
    relax_resultsdir = os.path.join(fragment_dir, e_function_name, 'relax/')
    if not os.path.isdir(relax_resultsdir):
        if verbose:
            print("Commencing relax runs")
        os.makedirs(relax_resultsdir)
        relax_commands_file = os.path.join(relax_resultsdir, 'commands.txt')
        with open(relax_commands_file, 'w') as f:
            for pool_i in list(range(relax_starting_pool_index, relax_starting_pool_index + n_relax_parallel_pools)):
                relax_silent_file = os.path.join(
                    relax_resultsdir, 'pool_{0}.out'.format(pool_i)
                )
                cmd = ' '.join([
                    minirosetta_path,
                    '-run:protocol relax',
                    '-in:file:native', pdb_path, #
                    '-in:file:s', pdb_path, #
                    '-frag3', frag3_file, #
                    '-frag9', frag9_file, #
                    '-nstruct', str(relax_nstruct),
                    '-out:file:silent', relax_silent_file,
                    '-silent_gz 1', #
                    '-score:weights', weights_file,
                    '-abinitio::rg_reweight 0.5', #
                    '-relax::default_repeats 5', #
                    '-ex1 1', #
                    '-ex2aro 1', #
                    '-mute all', #
                    '-in:file:fullatom 1'
                ])
                if extra_args:
                    for arg in extra_args:
                        cmd += ' {0}'.format(arg)
                if flags_file:
                    cmd += ' @{0}'.format(flags_file)
                f.write('{0}\n'.format(cmd))

        # Write an `sbatch` input file that can be used to carry out the commands
        # written to the above file
        relax_sbatch_file = os.path.join(relax_resultsdir, 'run.sbatch')
        queue_type = queue_type
        WriteSbatchFile(
            relax_sbatch_file,
            command_file_name=relax_commands_file,
            queue_type=queue_type
        )
    else:
        relax_commands_file = None
        relax_sbatch_file = None

    return (abinitio_commands_file, abinitio_sbatch_file, relax_commands_file, relax_sbatch_file)


def parse_silent_file(file):
    """
    Parse a silent file for the name, score, and rms of all structures
    """

    # Initiate a dictionary to keep track of relevant info
    # to turn into a dataframe
    scores_dict = {
        key : []
        for key in ['description', 'score', 'rms']
    }

    # Go through the file and record the above info for each
    # entry, identifying column names from a header line and
    # column values from subsequent rows, recognizing relevant
    # lines through unique strings.
    if file[-3:] == '.gz':
        gzip_file = True
        f = gzip.open(file, 'rb')
    else:
        gzip_file = False
        f = open(file)
    header_line_found = False
    for line in f:
        if gzip_file:
            line = line.decode('ascii')
        line = line.strip().split()
        if line[:2] == ['SCORE:', 'score']:
            header_line = line
            header_line_found = True
        elif 'SCORE:' in line:
            assert header_line_found == True, "Found a line of scores before the header line"
            assert len(line) == len(header_line)
            score_dict_i = {
                score_term : value
                for (score_term, value) in zip(header_line, line)
            }
            scores_dict['description'].append(score_dict_i['description'])
            scores_dict['score'].append(score_dict_i['score'])
            scores_dict['rms'].append(score_dict_i['rms'])
    f.close()

    # Turn the scores dictionary into a dataframe and return
    # the dataframe
    df = pandas.DataFrame.from_dict(scores_dict)
    return df

def cart_or_dualspace_relax(
    relax_app_path, input_pdbs, weights_file, results_dir, scores_file_prefix,
    relax_space='cartesian', nstruct=5, relax_script=False, default_repeats=5,
    extra_args=None, flags_file=None, submit_sbatch_job=True, queue_type='short',
    memory='2g', sbatch_group_size=None
    ):
    """Relax a design in cartesian or dual space"""

    # Assemble commands for all input pdbs
    cmds = []
    for pdb in input_pdbs:

        # Start assembling the command
        cmd = ' '.join([
            relax_app_path,
            '-s {0}'.format(pdb),
            '-use_input_sc',
            '-constrain_relax_to_start_coords',
            '-nstruct {0}'.format(nstruct),
            '-relax:coord_constrain_sidechains',
            '-relax:{0}'.format(relax_space),
            '-relax:default_repeats {0}'.format(default_repeats),
            '-relax:min_type lbfgs_armijo_nonmonotone',
            '-score:weights {0}'.format(weights_file),
            '-out:prefix {0}'.format(results_dir),
            '-out:file:scorefile {0}.sc'.format(scores_file_prefix),
        ])

        # Add additional flags relating to the relax protocol
        if relax_script:
            cmd += ' -relax:script {0}'.format(relax_script)
        else:
            cmd += ' -optimization:default_max_cycles 200'

        # Add any additional arguments or flags
        if extra_args:
            for arg in extra_args:
                cmd += ' {0}'.format(arg)
        if flags_file:
            cmd += ' @{0}'.format(flags_file)
        cmds.append(cmd)

    # Write the commands to a file
    command_file = '{0}_commands.txt'.format(scores_file_prefix)
    with open(command_file, 'w') as f:
        for cmd in cmds:
            f.write('{0}\n'.format(cmd))

    # Write a file that uses sbatch to execute the command in the above file
    sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)
    if not sbatch_group_size:
        WriteSbatchFile(
            sbatch_file_name,
            command_file_name=command_file,
            queue_type=queue_type,
            memory=memory
        )
        if submit_sbatch_job:
            subprocess.check_call(['sbatch', '-a', '1-{0}'.format(len(cmds)), sbatch_file_name])
    else:
        WriteSbatchFile(
            sbatch_file_name,
            command_file_name=command_file,
            queue_type=queue_type,
            memory=memory,
            group_size=sbatch_group_size,
            n_commands=len(cmds)
        )
        if submit_sbatch_job:
            subprocess.check_call(['sbatch', sbatch_file_name])


def relax_design(
    relax_app_path, input_pdbs, weights_file, results_dir, scores_file_prefix,
    relax_space=False, relax_script=False, extra_args=None,
    flags_file=None, submit_sbatch_job=True, queue_type='short', memory='2g',
    sbatch_group_size=None
    ):
    """Relax a design using the Rosetta relax application"""

    # Assemble commands for all input pdbs
    cmds = []
    for pdb in input_pdbs:

        # Start assembling the command
        cmd = ' '.join([
            relax_app_path,
            '-s {0}'.format(pdb),
            '-use_input_sc',
            '-relax:min_type lbfgs_armijo_nonmonotone',
            '-score:weights {0}'.format(weights_file),
            '-out:prefix {0}'.format(results_dir),
            '-out:file:scorefile {0}.sc'.format(scores_file_prefix),
        ])

        # Add additional flags relating to the relax protocol
        if relax_space:
            cmd += ' -relax:{0}'.format(relax_space)
        if relax_script:
            cmd += ' -relax:script {0}'.format(relax_script)

        # Add any additional arguments or flags
        if extra_args:
            for arg in extra_args:
                cmd += ' {0}'.format(arg)
        if flags_file:
            cmd += ' @{0}'.format(flags_file)
        cmds.append(cmd)

    # If there's just one command, then carry it out. Otherwise, if there's
    # more than one command, then write a file with the commands and then
    # carry them out in a single or multiple jobs
    sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)
    if len(cmds) == 1:
        WriteSbatchFile(
            sbatch_file_name,
            command=cmds[0],
            queue_type=queue_type,
            memory=memory
        )
        if submit_sbatch_job:
            subprocess.check_call(['sbatch', sbatch_file_name])
    else:
        # Write the commands to a file
        command_file = '{0}_commands.txt'.format(scores_file_prefix)
        with open(command_file, 'w') as f:
            for cmd in cmds:
                f.write('{0}\n'.format(cmd))

        # Write a file that uses sbatch to execute the command in the above file
        sbatch_file_name = '{0}.sbatch'.format(scores_file_prefix)
        if not sbatch_group_size:
            WriteSbatchFile(
                sbatch_file_name,
                command_file_name=command_file,
                queue_type=queue_type,
                memory=memory
            )
            if submit_sbatch_job:
                subprocess.check_call([
                    'sbatch', '-a', '1-{0}'.format(len(cmds)), sbatch_file_name
                ])
        else:
            WriteSbatchFile(
                sbatch_file_name,
                command_file_name=command_file,
                queue_type=queue_type,
                memory=memory,
                group_size=sbatch_group_size,
                n_commands=len(cmds)
            )
            if submit_sbatch_job:
                subprocess.check_call(['sbatch', sbatch_file_name])


def pre_ddG_cart_relax(
    relax_app_path, cart_relax_script, input_pdb, weights_file, relax_results_dir,
    nstruct=20, extra_args=None, flags_file=None, use_sbatch=True, run_sbatch=True,
    queue_type='short'
):
    """Relax a design in cartesian space for downstream use in the `cartesian_ddg` application"""

    # Assemble the command
    cmd = ' '.join([
        relax_app_path,
        '-s {0}'.format(input_pdb),
        '-use_input_sc',
        '-constrain_relax_to_start_coords',
        '-ignore_unrecognized_res',
        '-nstruct {0}'.format(nstruct),
        '-relax:coord_constrain_sidechains',
        '-relax:cartesian',
        '-score:weights {0}'.format(weights_file),
        '-relax:min_type lbfgs_armijo_nonmonotone',
        '-relax:script {0}'.format(cart_relax_script),
        '-out:prefix {0}'.format(relax_results_dir)
    ])
    if extra_args:
        for arg in extra_args:
            cmd += ' {0}'.format(arg)
    if flags_file:
        cmd += ' @{0}'.format(flags_file)

    # Carry out the command using `sbatch`
    if use_sbatch:
        sbatch_file_name = os.path.join(relax_results_dir, 'relax.sbatch')
        WriteSbatchFile(
            sbatch_file_name,
            command=cmd,
            queue_type=queue_type,
            memory='1g'
        )
        if run_sbatch:
            sbatch_cmd = ['sbatch', sbatch_file_name]
            subprocess.check_call(sbatch_cmd)
        return None
    else:
        return cmd


def perform_computational_ssm(
    input_pdb, weights_file, mutfilesdir, cartesian_ddg_path, databasepath, cartesian_ddg_resultsdir,
    extra_args=None, flags_file=None, mutation_subset=None, use_sbatch=True, run_sbatch=True, queue_type="short",
    verbose=False
):
    """
    Perform a computational SSM of an input protein using the `cartesian_ddg` application

    Must use *absolute* paths when specifying the below files, unless otherwise specified

    Args:
        `input_pdb`: the absolute path to the input PDB to use in the SSM
        `mutfilesdir`: the absolute path to an output directory in which to store mutation
            files that are used to run the `cartesian_ddg` applicaiton. This directory must
            already exist.
        `cartesian_ddg_path`: the absolute path to the `cartesian_ddg` application
        `databasepath`: the absolute path to the Rosetta `database`
        `cartesian_ddg_resultsdir`: the absolute path to a directory where the results of the
            `cartesian_ddg` application will be stored, including mutant PDBs and files
            with the estimated ddG values

    Returns:
        Mutation files in `mutfilesdir` and results of the `cartesian_ddg` application in
            `cartesian_ddg_resultsdir`
    """

    # First, read in the wildtype amino-acid sequence of the protein from the input PDB
    parser = Bio.PDB.PDBParser()
    ppb = Bio.PDB.PPBuilder()
    structure = parser.get_structure('test', input_pdb)
    sequence = ppb.build_peptides(structure)[0].get_sequence()

    # Then, generate input files with each single amino-acid mutation to test in the SSM,
    # with one file per mutation
    amino_acids = IUPAC.IUPACProtein.letters
    amino_acids = list(amino_acids)
    assert len(amino_acids) == 20
    mutations = []
    for (site, wt_aa) in enumerate(sequence, 1):
        for mut_aa in amino_acids:
            if wt_aa == mut_aa:
                continue
            mutation = '{0}{1}{2}'.format(wt_aa, site, mut_aa)
            mut_file = os.path.join(mutfilesdir, '{0}.txt'.format(mutation))
            mutations.append((mutation, mut_file))
            if not os.path.isfile(mut_file):
                with open(mut_file, 'w') as f:
                    f.write('total 1\n')
                    f.write('1\n')
                    f.write('{0} {1} {2}\n'.format(wt_aa, site, mut_aa))

    # Then, run the `cartesian_ddg` application one mutation at a time with each mutation file
    # as input
    for (mutation, mutfile) in mutations:

        # Don't rerun the `cartesian_ddg` application if there is already a results file
        cartesian_ddg_results_file = os.path.join(cartesian_ddg_resultsdir, '{0}.ddg'.format(mutation))
        if os.path.isfile(cartesian_ddg_results_file):
            continue

        # If `mutation_subset` is specified, only cary out the protocol for the specified subset
        # of mutations
        if mutation_subset:
            if mutation not in mutation_subset:
                continue

        # Assemble the command
        cmd = ' '.join([
            cartesian_ddg_path,
            '-database {0}'.format(databasepath),
            '-s {0}'.format(input_pdb),
            '-score:weights {0}'.format(weights_file),
            '-ddg:mut_file {0}'.format(mutfile),
            '-ddg:iterations 3', # can be flexible; 3 is fast and reasonable
            '-ddg::cartesian',
            '-ddg::dump_pdbs true',
            '-bbnbr 1', # bb dof, suggestion: i-1, i, i+1
            '-fa_max_dis 9.0', # modify fa_atr and fa_sol behavior, really important for protein stability (default: 6)
        ])
        if extra_args:
            for arg in extra_args:
                cmd += ' {0}'.format(arg)
        if flags_file:
            cmd += ' @{0}'.format(flags_file)

        # Carry out the command using `sbatch`
        if use_sbatch:
            sbatch_file_name = os.path.join(cartesian_ddg_resultsdir, '{0}.sbatch'.format(mutation))
            WriteSbatchFile(
                sbatch_file_name,
                command=cmd,
                queue_type=queue_type,
                memory='1g'
            )
            sbatch_cmd = ['sbatch', '{0}.sbatch'.format(mutation)]
            if run_sbatch:
                process = subprocess.Popen(
                    sbatch_cmd,
                    cwd=cartesian_ddg_resultsdir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE
                )
                out, err = process.communicate()
                if verbose:
                    print("out:\n{0}".format(out))
                    print("err:\n{0}".format(err))

    pass # Don't return anything


def read_cartesian_ddg_outfile(file):
    """
    Read in scores from an output file generated by the `cartesian_ddg` application
    """

    # Identify the positions of columns with strings that describes what the numerical
    # value is in the column that immediately follows
    header_column_indices = range(4, 56, 2)

    # Extract data line by line using the header columns as guides and storing the results
    # in a dictionary
    assert os.path.isfile(file), "Could not find the file: {0}".format(file)
    with open(file) as f:
        for (i, line) in enumerate(f):
            line = line.strip().split()
            assert len(line) == 56

            # If it is the first line, then initiate a dictionary to track the results
            if i == 0:
                header_columns = [line[header_index] for header_index in header_column_indices]
                scores_dict = {
                    key : []
                    for key in ['round', 'mut_or_wt', 'total_score'] + header_columns
                }
            scores_dict['round'].append(line[1][-2])
            mut_or_wt = line[2][:-1]
            if 'WT' in mut_or_wt:
                assert 'MUT' not in mut_or_wt
                scores_dict['mut_or_wt'].append('wt')
            elif 'MUT' in mut_or_wt:
                scores_dict['mut_or_wt'].append('mut')
            else:
                raise ValueError('Could not parse `WT` or `MUT` from the string: {0}'.format(mut_or_wt))
            scores_dict['total_score'].append(float(line[3]))
            for header_index in header_column_indices:
                scores_dict[line[header_index]].append(line[header_index + 1])

    # Return a dataframe of the scores
    scores_df = pandas.DataFrame.from_dict(scores_dict)
    return scores_df


def forward_fold_designs_on_boinc(
    pdb_paths, forward_folding_dir, convert_rel_to_abs_path=True,
    boinc_id_prefix=None, run_jobs=False, nj_relax=10, nj_abinitio=200
    ):
    """
    Forward fold a set of designs
    """

    # If indicated by input, as is the default, convert relative paths to
    # absolute paths
    if convert_rel_to_abs_path:
        cwd = os.getcwd()
        pdb_paths = [os.path.join(cwd, pdb_path) for pdb_path in pdb_paths]

    # Make a list of PDB files that do NOT already have a directory with the
    # results of fragment picking
    pdbs_for_fragment_picking = [
        pdb_path for pdb_path in pdb_paths
        if not os.path.isdir(os.path.join(
            forward_folding_dir,
            os.path.basename(pdb_path).replace('.pdb', '') + '_fragments'
        ))
    ]

    # Make fragments for PDBs that passed the above step
    if pdbs_for_fragment_picking:
        print("Making fragments")
        cmd = [
            'python',
            '/home/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/make_fragments.py',
            '-pdbs'
        ] + pdbs_for_fragment_picking
        make_fragments_outfile = os.path.join(
            forward_folding_dir, 'make_fragments.out'
        )
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, cwd=forward_folding_dir
        )
        (out, err) = process.communicate()
        with open(make_fragments_outfile, 'wb') as f:
            f.write(out)
        if err:
            print(err)
            raise ValueError('Error in making fragments. The error message is:\n{0}'.format(err))

    # Prepare `boinc` job for each set of fragments
    fragment_dirs = [
        os.path.join(
            forward_folding_dir,
            os.path.basename(pdb_path).replace('.pdb', '') + '_fragments'
        ) for pdb_path in pdb_paths
    ]
    for fragment_dir in fragment_dirs:

        if not os.path.isdir(fragment_dir):
            print("Missing the expected fragments directory: {0}".format(fragment_dir))
            continue

        # Specify a unique ID for the job
        boinc_run_id = os.path.split(fragment_dir)[1].replace('_fragments', '')
        if boinc_id_prefix:
            boinc_run_id = boinc_id_prefix + boinc_run_id
        assert len(boinc_run_id) >= 8

        # Run the command to prepare the job if the target outfile doesn't exist
        cmd = [
            'python2',
            '/home/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/prepare_fold.py',
            '-nj_relax', '{0}'.format(nj_relax),
            '-nj_abinitio', '{0}'.format(nj_abinitio),
            boinc_run_id
        ]
        cmd = ' '.join(cmd)
        run_job_file = os.path.join(fragment_dir, 'run.fold.boinc.job')
        if not os.path.isfile(run_job_file):
            print("Preparing boinc jobs for the fragment directory: {0}".format(fragment_dir))
            prepare_job_outfile = os.path.join(fragment_dir, 'prepare_job.out')
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=fragment_dir, shell=True)
            (out, err) = process.communicate()
            with open(prepare_job_outfile, 'wb') as f:
                f.write(out)
            if err:
                print(err)
                raise ValueError('Error in preparing boinc job. The error message is:\n{0}'.format(err))

        # Submit the job to boinc
        if run_jobs:
            submit_job_outfile = os.path.join(fragment_dir, 'submit_job.out')
            if not os.path.isfile(submit_job_outfile):
                print("Submitting boinc job for the fragment directory: {0}".format(fragment_dir))
                cmd = '/projects/boinc/bin/boinc_submit run.fold.boinc.job'
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=fragment_dir, shell=True)
                (out, err) = process.communicate()
                with open(submit_job_outfile, 'wb') as f:
                    f.write(out)
                if err:
                    print(err)
                    raise ValueError('Error in submitting boinc job. The error message is:\n{0}'.format(err))
            #else:
                #print("Job already submitted for the fragment directory: {0}".format(fragment_dir))


def copy_boinc_results(
    run_name, boinc_results_dir, new_results_dir, verbose=False, overwrite=False, file_suffix='.sc.bz2',
    file_suffixes_to_ignore=None
):
    "Copy results from `boinc` to a target directory"

    # Define the directory where the output files are stored
    run_results_dir = os.path.join(
        boinc_results_dir,
        run_name[:8]
    )

    # Copy output scores files to the target destination, one
    # file for the "fold" job, one for the "relax" job
    for run_type in ['fold', 'relax', 'abinitio']:

        # Find the output file
        output_file = glob.glob(os.path.join(
            run_results_dir,
            '{0}_{1}_*{2}'.format(run_name, run_type, file_suffix)
        ))
        if file_suffixes_to_ignore:
            for suffix in file_suffixes_to_ignore:
                output_file = [
                    f for f in output_file
                    if suffix not in f
                ]
        if len(output_file) == 0:
            if verbose:
                print("Could not find an output file for {0} and {1}".format(
                    run_name, run_type
                ))
            continue
        elif len(output_file) > 1:
            raise ValueError("Found {0} {1} output files: {2}".format(
                len(output_file),
                run_type,
                ' '.join(output_file)
            ))
        output_file = output_file[0]

        # Copy the file to the target destination if the file
        # doesn't already exist
        new_output_file = os.path.join(
            new_results_dir,
            os.path.basename(output_file)
        )
        if not os.path.isfile(new_output_file) or overwrite:
            cmd = ['cp', output_file, new_output_file]
            subprocess.check_call(cmd)

    pass


def copy_boinc_results2(
    batch_id, boinc_results_dir, new_results_dir, verbose=False, overwrite=False,
    file_suffix='.sc.bz2', file_suffixes_to_ignore=None
):
    "Copy results from `boinc` to a target directory"

    # Define the directory where the output files are stored
    run_results_dir = os.path.join(
        boinc_results_dir,
        batch_id
    )

    # Copy output scores files to the target destination, one
    # file for the "fold" job, one for the "relax" job
    for run_type in ['fold', 'relax', 'abinitio']:

        # Find the output file
        output_file = glob.glob(os.path.join(
            run_results_dir,
            '*_{0}_*{1}'.format(run_type, file_suffix)
        ))
        if file_suffixes_to_ignore:
            for suffix in file_suffixes_to_ignore:
                output_file = [
                    f for f in output_file
                    if suffix not in f
                ]
        if len(output_file) == 0:
            if verbose:
                print("Could not find an output file for {0}".format(run_type))
            continue
        elif len(output_file) > 1:
            raise ValueError("Found {0} {1} output files: {2}".format(
                len(output_file),
                run_type,
                ' '.join(output_file)
            ))
        output_file = output_file[0]

        # Copy the file to the target destination if the file
        # doesn't already exist
        new_output_file = os.path.join(
            new_results_dir,
            os.path.basename(output_file)
        )
        if not os.path.isfile(new_output_file) or overwrite:
            cmd = ['cp', output_file, new_output_file]
            subprocess.check_call(cmd)

    pass


def parse_boinc_silent_scores_file(file, max_datapoints=None):
    """
    Parse a silent scores file from a forward-folding run on boinc
    for the name, score, and rms of all structures
    """

    # Initiate a dictionary to keep track of relevant info
    # to turn into a dataframe
    scores_dict = {
        key : []
        for key in ['description', 'score', 'rms']
    }

    # Go through the file and record the above info for each
    # entry, identifying column names from a header line and
    # column values from subsequent rows, recognizing relevant
    # lines through unique strings.
    if file[-3:] == '.gz':
        binary_data = True
        f = gzip.open(file, 'rb')
    elif file[-4:] == '.bz2':
        binary_data = True
        f = bz2.BZ2File(file, 'r')
    else:
        binary_data = False
        f = open(file, 'r')
    header_line_found = False
    for line in f:
        if binary_data:
            line = line.decode('ascii')
        line = line.strip().split()
        if line[:2] == ['SCORE:', 'score']:
            header_line = line
            header_line_found = True
        elif 'SCORE:' in line:
            assert header_line_found == True, "Found a line of scores before the header line"
            # Continue if the line with data has a fewer number of elements than
            # the header line
            if len(line) < len(header_line):
                continue
            # ... otherwise, record data for that line
            score_dict_i = {
                score_term : value
                for (score_term, value) in zip(header_line, line)
            }
            if 'description' not in score_dict_i.keys():
                print("Failed to find expected columns in the following line from the file: {0}".format(file))
                print(line)
            scores_dict['description'].append(score_dict_i['description'])
            scores_dict['score'].append(float(score_dict_i['score']))
            if 'rms' in score_dict_i:
                scores_dict['rms'].append(float(score_dict_i['rms']))
            else:
                scores_dict['rms'].append(numpy.nan)
        if max_datapoints:
            if len(scores_dict['description']) > max_datapoints:
                break
    f.close()

    # Turn the scores dictionary into a dataframe and return
    # the dataframe
    df = pandas.DataFrame.from_dict(scores_dict)
    return df


def copy_pdbs(pdbs, new_pdb_dir):
    """Copy a list of PDBs to a target directory"""
    if not os.path.isdir(new_pdb_dir):
        os.makedirs(new_pdb_dir)
    for pdb in pdbs:
        assert os.path.isfile(pdb)
        new_pdb_path = os.path.join(
            new_pdb_dir,
            os.path.basename(pdb)
        )
        if not os.path.isfile(new_pdb_path):
            subprocess.check_call(['cp', pdb, new_pdb_path])


def compute_hydropathy_of_sequence(seq, only_consider_AFILMVWY=False):
    """
    Compute the hydrophathy for a given sequence

    Specifically, sum the hydropathy of each amino acid across a
    sequence, where values of hydropathy are taken from Table 2 of
    Kyte, J. and Doolittle, R. 1982. A simple method for displaying
    the hydropathic character of a protein. J. Mol. Biol. 157: 105-132.

    Args:
        `seq`: a string giving the sequence of the protein
        `only_consider_AFILMVWY`: a bool. If True, this function
            will only sum amino-acid hydropathies over hydrophobic
            amino acids (AFILMVWY)
    Returns:
        `hydropathy`: a float giving the hydropathy of the input
            sequence
    """

    # Define the hydropathy of each amino acid according to Table 2
    # from Kyte and Doolittle, 1982, J. Mol. Biol.
    aa_hydropathy = {
        "I": 4.5,
        "V": 4.2,
        "L": 3.8,
        "F": 2.8,
        "C": 2.5,
        "M": 1.9,
        "A": 1.8,
        "G": -0.4,
        "T": -0.7,
        "W": -0.9,
        "S": -0.8,
        "Y": -1.3,
        "P": -1.6,
        "H": -3.2,
        "E": -3.5,
        "Q": -3.5,
        "D": -3.5,
        "N": -3.5,
        "K": -3.9,
        "R": -4.5
    }

    # Sum the hydropathy of each amino acid in the sequence, only
    # considering hydrophobic amino acids if specified above
    if only_consider_AFILMVWY:
        return sum([
            aa_hydropathy[aa] for aa in seq
            if aa in list('AFILMVWY')
        ])
    else:
        return sum([aa_hydropathy[aa] for aa in seq])

def divJensenShannon(p1, p2):
    """Jensen-Shannon divergence between two distributions.

    From: https://jbloomlab.github.io/dms_tools2/_modules/dms_tools2/compareprefs.html#divJensenShannon

    The logarithms are taken to base 2, so the result will be
    between 0 and 1.

    Args:
        `p1` and `p2` (array-like)
            The two distributions for which we compute divergence.

    Returns:
        The Jensen-Shannon divergence as a float.

    >>> p1 = [0.5, 0.2, 0.2, 0.1]
    >>> p2 = [0.4, 0.1, 0.3, 0.2]
    >>> p3 = [0.0, 0.2, 0.2, 0.6]
    >>> numpy.allclose(divJensenShannon(p1, p1), 0, atol=1e-5)
    True
    >>> numpy.allclose(divJensenShannon(p1, p2), 0.035789, atol=1e-5)
    True
    >>> numpy.allclose(divJensenShannon(p1, p3), 0.392914, atol=1e-5)
    True
    """
    p1 = numpy.array(p1)
    p2 = numpy.array(p2)

    def _kldiv(a, b):
        with numpy.errstate(all='ignore'):
            kl = numpy.sum([
                v for v in a * numpy.log2(a / b)
                if not numpy.isnan(v)
            ])
        return kl

    m = 0.5 * (p1 + p2)

    return 0.5 * (_kldiv(p1, m) +_kldiv(p2, m))


def count_number_of_codons_in_sequence(seq, codons):
    """
    Count the number of codons in an input sequence

    Code for doctest:
    >>> seq = 'GCCATGCAAGCTTTTGCT'
    >>> codons = ['GCT', 'GCC', 'GCA', 'GCG']
    >>> counts = count_number_of_codons_in_sequence(seq, codons)
    >>> (counts['GCT'], counts['GCC'], counts['GCA'], counts['GCG'])
    (2, 1, 0, 0)
    >>> seq = 'AGCCATGCAAGCTTTTGC'
    >>> counts = count_number_of_codons_in_sequence(seq, codons)
    >>> (counts['GCT'], counts['GCC'], counts['GCA'], counts['GCG'])
    (0, 0, 1, 0)
    """

    codon_counts_dict = {
        codon : 0
        for codon in codons
    }
    seq = iter(list(seq))
    while True:
        try:
            codon = ''.join([
                next(seq), next(seq), next(seq)
            ])
            if codon in codons:
                codon_counts_dict[codon] += 1
        except StopIteration:
            break
    return codon_counts_dict

def mutate_protein_sequence(seq, site_n, wt_aa, mut_aa):
    """
    Introduce an amino-acid mutation into a protein sequence

    Code for doctest:
    >>> seq = 'MTREIPLLG'
    >>> site_n = 2
    >>> wt_aa = 'T'
    >>> mut_aa = 'V'
    >>> mutate_protein_sequence(seq, site_n, wt_aa, mut_aa)
    'MVREIPLLG'
    """

    # Go through each site in a sequence and mutate the target site
    mut_seq = []
    for (site_i, aa_i) in enumerate(seq, 1):
        if site_i == int(site_n):
            if aa_i != wt_aa:
                raise ValueError("Expected the wildtype amino-acid {0}, but found {1}".format(wt_aa, mut_aa))
            mut_seq.append(mut_aa)
        else:
            mut_seq.append(aa_i)
    mut_seq = ''.join(mut_seq)
    assert len(seq) == len(mut_seq)
    return mut_seq

def mutate_codon_in_sequence(seq, codon_n, wt_aa, mut_codon, mut_aa):
    """
    Introduce an codon mutation into a DNA sequence

    Code for doctest:
    >>> seq = 'ATGACCACC'
    >>> codon_n = 2
    >>> wt_aa = 'T'
    >>> mut_codon = 'GTC'
    >>> mut_aa = 'V'
    >>> mutate_codon_in_sequence(seq, codon_n, wt_aa, mut_codon, mut_aa)
    'ATGGTCACC'
    """

    # Make sure that the mutant codon matches the expected mutant
    # amino acid
    translated_mut_aa = str(
        Seq(mut_codon, IUPAC.unambiguous_dna).translate()
    )
    if translated_mut_aa != mut_aa:
        raise ValueError("The tranlated mutated codon gives {0}, but found {1}".format(translated_mut_aa, mut_aa))

    # Iterate through codons in a sequence and mutate the target one
    mut_dna_seq = []
    dna_iter = iter(list(seq))
    codon_i = 1
    found_codon_to_mutate = False
    while True:

        # Try getting the next codon by iterating on the
        # next three nucleotides in the sequence. If there
        # are no codons left, except the error and break the loop
        try:
            # Get the codon
            codon = ''.join([
                next(dna_iter), next(dna_iter), next(dna_iter)
            ])

            # See if this codon is the one to mutate. If so,
            # then check that it translates to the expected
            # amino acid, and make the codon mutation
            if codon_i == int(codon_n):
                assert not found_codon_to_mutate
                aa_i = str(
                    Seq(codon, IUPAC.unambiguous_dna).translate()
                )
                if aa_i != wt_aa:
                    raise ValueError("Expected the wildtype amino-acid {0}, but found {1}".format(wt_aa, aa_i))
                mut_dna_seq.append(mut_codon)
                found_codon_to_mutate = True
            else:
                mut_dna_seq.append(codon)
            codon_i += 1

        except StopIteration:
            break

    # Return the full-length mutated sequence
    assert found_codon_to_mutate, "Did not find codon to mutate"
    mut_dna_seq = ''.join(mut_dna_seq)
    assert len(seq) == len(mut_dna_seq)
    return mut_dna_seq

def compute_distance_between_codons(codon_x, codon_y):
    """
    Compute the edit distance between two codons

    Code for doctest:
    >>> codon_x = 'ACG'
    >>> codon_y = 'AGG'
    >>> codon_z = 'GCC'
    >>> compute_distance_between_codons(codon_x, codon_x)
    0
    >>> compute_distance_between_codons(codon_x, codon_y)
    1
    >>> compute_distance_between_codons(codon_x, codon_z)
    2
    >>> compute_distance_between_codons(codon_y, codon_z)
    3
    """
    assert len(codon_x) == len(codon_y) == 3
    dist = sum([i!=j for (i, j) in zip(list(codon_x), list(codon_y))])
    return dist

def get_codon_from_dna_seq(dna_seq, codon_n):
    """
    Get a codon from a DNA sequence based on the codon number, with
    numbering indexed starting at 1.

    Code for doctest:
    >>> dna_seq = 'ACGTTTCAC'
    >>> get_codon_from_dna_seq(dna_seq, 1)
    'ACG'
    >>> get_codon_from_dna_seq(dna_seq, 2)
    'TTT'
    >>> get_codon_from_dna_seq(dna_seq, 3)
    'CAC'
    """
    starting_nt_index = (codon_n - 1) * 3
    codon = dna_seq[starting_nt_index : starting_nt_index + 3]
    return codon


def forward_fold_designs_from_sequence_on_boinc(
    fastas_folder, forward_folding_dir, setup_boinc_alone_script,
    convert_rel_to_abs_path=True, boinc_id_prefix=None, run_jobs=False,
    nj_relax=10, nj_abinitio=200
    ):
    """
    Forward fold a set of designs starting from the sequence (no PDB considered)
    """

    # If indicated by input, as is the default, convert relative paths to
    # absolute paths
    if convert_rel_to_abs_path:
        cwd = os.getcwd()
        fastas_folder = os.path.join(cwd, fastas_folder)

    # Make fragments
    make_fragments_outfile = os.path.join(
        forward_folding_dir, 'make_fragments.out'
    )
    if not os.path.isfile(make_fragments_outfile):
        print("Making fragments in the directory: {0}".format(forward_folding_dir))
        cmd = ' '.join([
            'python',
            '/home/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/make_fragments_fasta.py',
            '-fastas_folder {0}'.format(fastas_folder)
        ])
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=forward_folding_dir, shell=True)
        (out, err) = process.communicate()
        with open(make_fragments_outfile, 'wb') as f:
            f.write(out)
        if err:
            print(err)
            raise ValueError('Error in making fragments. The error message is:\n{0}'.format(err))

    # Prepare `boinc` job for each set of fragments
    prepare_job_outfile = os.path.join(forward_folding_dir, 'prepare_job.out')
    if not os.path.isfile(prepare_job_outfile):
        print("Preparing boinc jobs in the directory: {0}".format(forward_folding_dir))
        cmd = ['ruby', setup_boinc_alone_script]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=forward_folding_dir)
        (out, err) = process.communicate()
        with open(prepare_job_outfile, 'wb') as f:
            f.write(out)
        if err:
            print(err)
            raise ValueError('Error in preparing boinc job. The error message is:\n{0}'.format(err))

    # Submit the job to boinc
    if run_jobs:
        fragment_dirs = [
            d for d in glob.glob(os.path.join(forward_folding_dir, '*_fragments'))
            if os.path.isdir(d)
        ]
        for fragment_dir in fragment_dirs:
            submit_job_outfile = os.path.join(fragment_dir, 'submit_job.out')
            if not os.path.isfile(submit_job_outfile):
                print("Submitting boinc job for the fragment directory: {0}".format(fragment_dir))
                cmd = '/projects/boinc/bin/boinc_submit run.fold.boinc.job'
                process = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=fragment_dir, shell=True)
                (out, err) = process.communicate()
                with open(submit_job_outfile, 'wb') as f:
                    f.write(out)
                if err:
                    print(err)
                    raise ValueError('Error in submitting boinc job. The error message is:\n{0}'.format(err))
            #else:
                #print("Job already submitted for the fragment directory: {0}".format(fragment_dir))

                
def extract_pdb_from_boinc_output(
    working_dir, silent_file_basename, tags, output_dir=False
):

    # Unzip silent file
    silent_file = os.path.join(working_dir, silent_file_basename)
    if (os.path.isfile(silent_file)) and ('.bz2' in silent_file):
        cmd = ' '.join(['bzip2', '-d', silent_file_basename])
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, cwd=working_dir, shell=True
        )
        (out, err) = process.communicate()
        if err:
            print(err)
            raise ValueError('Error in unzipping file. The error message is:\n{0}'.format(err))
        time.sleep(10)

    # Extract structure
    unzipped_silent_file = silent_file_basename.replace('.bz2', '')
    if not os.path.isfile(
        os.path.join(working_dir, unzipped_silent_file)
    ):
        return "No file: {0}".format(unzipped_silent_file)
    cmd = ' '.join([
        '/software/rosetta/versions/v2018.31-dev60339/main/source/bin/extract_pdbs',
        '-database /software/rosetta/versions/v2018.31-dev60339/database/',
        '-in:file:silent {0}'.format(unzipped_silent_file),
        '-in:file:tags {0}'.format(' '.join(tags)),
        '-in:file:silent_struct_type binary',
        '-silent_read_through_errors'
    ])
    if output_dir:
        cmd += ' -out:prefix {0}'.format(output_dir)
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, cwd=working_dir, shell=True
    )
    (out, err) = process.communicate()
    if err:
        print(err)
        raise ValueError('Error in extracting PDBs. The error message is:\n{0}'.format(err))

                
                
if __name__ == "__main__":
    import doctest
    doctest.testmod()
