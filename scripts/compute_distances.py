"""Compute energies between all pairs of atoms in a structure"""

# Import `Python` modules
import os
import sys
import argparse
import glob
import subprocess
import pandas
import numpy as np
import math
import pyrosetta
import pyrosetta.rosetta
flags = [
    '-beta_nov16', '-corrections:beta_nov16',
    '-mute all', '-mute core', '-mute protocols',
]
pyrosetta.init(extra_options=' '.join(flags))
sf = pyrosetta.get_fa_scorefxn()
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)
fix_scorefxn(sf)

# Define functions
def compute_distances(
    pdb,
    use_tenA_neighbor_residues=True,
    analyze_all_sites=True,
    sites_to_analyze=[],
):
    """
    Compute distances between all pairs of atoms in a protein
    
    Note, this function does not consider:
        * atom pairs separated by distances greater than
            0.5 A above the sum of the VDW radii of the
            two atoms
        * atom pairs from the same residue, or atom from
            adjacent residues where at least one of the
            atoms is a backbone atom
    
    Args:
        *pdb*: the path to an input PDB file
        *use_tenA_neighbor_residues*: boolean (default: True). If
            True, then this function only considers atom pairs
            between residues that are near each other in 3D space,
            defined as the two Cb atoms being within 10A of each
            other
        *analyze_all_sites*: boolean (default: True). If False,
            only compute distances across the subset of sites
            in the list passed to the argument *sites_to_analyze*
        *sites_to_analyze*: If *analyze_all_sites* is False, then
            only analyze the site numbers given in this list.
    Returns:
        A dataframe where each row corresponds to an atom pair
            and columns provide metadata, such as inter-atomic
            distance (d) and the names of the residues and atoms
            involved
    """
    
    # Read in PDB as pose and score the pose
    pose = pyrosetta.pose_from_pdb(pdb)
    sf(pose)

    # Initiate a dictionary for storing inter-atomic
    # distances
    energies_dict = {
        key : []
        for key in [
            'd',

            'res_i_n', 'res_i_pdb_n', 'res_i_name', 'res_i_chain',
            'atom_i_n', 'atom_i_name', 'atom_i_type_name',
            'atom_i_bb', 'atom_i_lj_radius',

            'res_j_n', 'res_j_pdb_n', 'res_j_name', 'res_j_chain',
            'atom_j_n', 'atom_j_name', 'atom_j_type_name',
            'atom_j_bb', 'atom_j_lj_radius',
        ]
    }

    # Loop over all residues in the protein
    pose_size = pose.size()
    res_i_ns = list(range(1, pose_size+1))
    for res_i_n in res_i_ns:

        # If you specified a subset of sites to analyze,
        # then skip over sites that aren't in this subset
        if not analyze_all_sites:
            if res_i_n not in sites_to_analyze:
                continue

        # Make a list of neighboring residues to loop over
        # when computing inter-atomic distances
        res_i = pose.residue(res_i_n)
        if use_tenA_neighbor_residues:
            neighbors = pyrosetta.rosetta.core.select.get_tenA_neighbor_residues(
                pose,
                pyrosetta.Vector1([
                    i == res_i_n for i in range(1, pose.size()+1)
                ])
            )
            res_j_ns = [
                res_n for (bool_n, res_n)
                in zip(neighbors, range(1, pose_size+1))
                if bool_n
            ]
        else:
            res_j_ns = list(range(1, pose.size()+1))

        # Loop over all neighbors to residue i, computing
        # all inter-atomic distances between residue pairs.
        for res_j_n in res_j_ns:
            
            # Skip over a residue pair of res_i_n => res_j_n so
            # as to avoid double counting and to skip distances
            # between atoms in the same residues. If specified,
            # also skip over sites that shouldn't be analyzed.
            if res_i_n >= res_j_n:
                continue
            if not analyze_all_sites:
                if res_j_n not in sites_to_analyze:
                    continue
            res_j = pose.residue(res_j_n)
            for atom_i_n in list(range(1, res_i.natoms()+1)):
                
                # Get atom name and coordinates
                atom_type_i = res_i.atom_type(atom_i_n)
                atom_i_xyz = res_i.xyz(atom_i_n)

                for atom_j_n in list(range(1, res_j.natoms()+1)):

                    # Get atom name and corrdinates
                    atom_type_j = res_j.atom_type(atom_j_n)
                    atom_j_xyz = res_j.xyz(atom_j_n)

                    # Record distance and other metadata
                    energies_dict['d'].append(
                        (atom_i_xyz - atom_j_xyz).norm()
                    )
                    energies_dict['res_i_n'].append(res_i_n)
                    energies_dict['res_i_pdb_n'].append(
                        pose.pdb_info().number(res_i_n)
                    )
                    energies_dict['res_i_name'].append(
                        res_i.name3()
                    )
                    energies_dict['res_i_chain'].append(
                        pose.pdb_info().chain(res_i_n)
                    )
                    energies_dict['atom_i_n'].append(atom_i_n)
                    energies_dict['atom_i_name'].append(
                        res_i.atom_name(atom_i_n).strip()
                    )
                    energies_dict['atom_i_type_name'].append(
                        atom_type_i.name().strip()
                    )
                    energies_dict['atom_i_bb'].append(
                        res_i.atom_is_backbone(atom_i_n)
                    )
                    energies_dict['atom_i_lj_radius'].append(
                        atom_type_i.lj_radius()
                    )

                    energies_dict['res_j_n'].append(res_j_n)
                    energies_dict['res_j_pdb_n'].append(
                        pose.pdb_info().number(res_j_n)
                    )
                    energies_dict['res_j_name'].append(res_j.name3())
                    energies_dict['res_j_chain'].append(
                        pose.pdb_info().chain(res_j_n)
                    )
                    energies_dict['atom_j_n'].append(atom_j_n)
                    energies_dict['atom_j_name'].append(
                        res_j.atom_name(atom_j_n).strip()
                    )
                    energies_dict['atom_j_type_name'].append(
                        atom_type_j.name().strip()
                    )
                    energies_dict['atom_j_bb'].append(
                        res_j.atom_is_backbone(atom_j_n)
                    )
                    energies_dict['atom_j_lj_radius'].append(
                        atom_type_j.lj_radius()
                    )

    energies_df = pandas.DataFrame(energies_dict)
    
    # Subset to inter-atomic distances of d-o <= 0.5 A,
    # where o is the sum of the VDW radii
    energies_df['o'] = \
        energies_df['atom_i_lj_radius'] + \
        energies_df['atom_j_lj_radius']
    energies_df['d-o'] = energies_df['d'] - energies_df['o']
    energies_df = energies_df[
        energies_df['d-o'] <= 0.5
    ]
    
    # Compute the distance between residues in primary
    # sequence. Then, drop rows where the two atoms
    # are from two adjacent residues and at least one
    # atom is a backbone atom
    energies_df['seq_dist'] = \
        energies_df['res_i_n'] - energies_df['res_j_n']
    energies_df['seq_dist'] = energies_df['seq_dist'].abs()
    energies_df['drop'] = \
        (energies_df['seq_dist'] < 2) & \
        (energies_df['atom_i_bb'] | energies_df['atom_j_bb'])
    energies_df = energies_df[~energies_df['drop']]
    
    # Add a column that gives the atom pair, with the
    # two atom atom names sorted in alphabetical order
    energies_df['atom_pair'] = energies_df.apply(
        lambda row: ':'.join(sorted([
            row['atom_i_type_name'], row['atom_j_type_name']
        ])),
        axis=1
    )
    
    # Add columns indicating rows with certain combinations
    # of hydrophobic atoms
    hydrophobic_residues = [
        'ALA', 'VAL', 'LEU', 'ILE', 'MET',
        'PHE', 'TRP', 'TYR',
    ]
    hydrophobic_carbons = ['CH1', 'CH2', 'CH3', 'CH0', 'aroC']
    hydrophobic_hydrogens = ['Hapo', 'Haro']
    for x in ['i', 'j']:
        energies_df[f'hydrophobic_carbon_{x}'] = \
            (energies_df[f'res_{x}_name'].isin(hydrophobic_residues) &
            energies_df[f'atom_{x}_type_name'].isin(hydrophobic_carbons))
        energies_df[f'hydrophobic_hydrogen_{x}'] = \
            (energies_df[f'res_{x}_name'].isin(hydrophobic_residues) &
            energies_df[f'atom_{x}_type_name'].isin(hydrophobic_hydrogens))

    energies_df['hydrophobic_C_C'] = \
        (energies_df['hydrophobic_carbon_i'] &
        energies_df['hydrophobic_carbon_j'])
    energies_df['hydrophobic_C_H'] = (
        (
            energies_df['hydrophobic_carbon_i'] &
            energies_df['hydrophobic_hydrogen_j']
        ) |
        (
            energies_df['hydrophobic_carbon_j'] &
            energies_df['hydrophobic_hydrogen_i']
        ))
    energies_df['hydrophobic_H_H'] = \
        ((energies_df['hydrophobic_hydrogen_i']) &
        (energies_df['hydrophobic_hydrogen_j']))
    energies_df['C_Obb'] = energies_df['atom_pair'].isin([
        'CH1:OCbb', 'CH2:OCbb', 'CH3:OCbb',
        'OCbb:aroC', 'CH0:OCbb',
    ])
    
    return energies_df

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Run the main code
def main():
    """Read in command-line arguments and execute the main code"""

    # Read in command-line arguments using `argparse`
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", help="the path to a directory with input PDBs")
    parser.add_argument("--file_with_thresholds", help="the path to a file with thresholds used to determine if a pair of atoms is clashing. The file `results/natives/thresholds.csv` has thresholds computed from a set of ~80 high-resolution crystal structures of native proteins", default='')
    parser.add_argument("--use_tenA_neighbor_residues", type=str2bool, default=True, help="a boolean (default: True). Only compute distances between atoms from residues with C-beta atoms within 10 Angstroms. This helps with speed and appears to capture nearly all clashes.")
    parser.add_argument("--output_file_prefix", help="a prefix that that will be used as the start of the path of output files. If `--report_clashes` is set to True (default), the script will output a file with the suffix `n_clashes.csv` that quantifies the number of clashes that surpass thresholds from the above input file. If `--report_distances` is set to True (default is False), then the script will also output a file with all inter-atomic distances used to compute the number of clashes.")
    parser.add_argument("--report_clashes", type=str2bool, default=True, help="a boolean (default: True). Write an output file with the number of clashes per structure surpassing thresholds")
    parser.add_argument("--report_distances", type=str2bool, default=False, help="a boolean (default: False). Write an output file with inter-atomic distances")
    args = parser.parse_args()
    
    # Read in thresholds
    if args.report_clashes:
        thresholds_df = pandas.read_csv(args.file_with_thresholds)
        atom_pairs = sorted(list(set(thresholds_df['atom_pair'])))
        threshold_quantiles = [
            col for col in thresholds_df.columns.values
            if 'q' in col
        ]

        # Initialize a dictionary for recording the number of
        # clashes per structure surpassing the thresholds
        n_clashes_dict = {
            'pdb' : []
        }
        for q in threshold_quantiles:
            n_clashes_dict.update({
                f'{atom_pair}_{q}' : []
                for atom_pair in atom_pairs
        })
    
    # Cycle over PDBs, compute inter-atomic distances for each one,
    # and compute the number of clashes that surpass thresholds
    dfs = []
    pdbs = glob.glob(os.path.join(args.pdb_dir, '*.pdb'))
    for pdb in pdbs:
        
        print(f"Analyzing {pdb}")
        
        # Compute inter-atomic distances
        df = compute_distances(pdb, args.use_tenA_neighbor_residues)
        df['pdb'] = pdb
        if args.report_distances:
            dfs.append(df)
        
        # Iterate over atom pairs and record the number of
        # clashes that surpass thresholds specific to the
        # atom pair
        if args.report_clashes:
            n_clashes_dict['pdb'].append(pdb)
            for atom_pair in atom_pairs:
                for q in threshold_quantiles:

                    # Get the threshold for a specific atom pair
                    # and quantile
                    threshold = float(thresholds_df[
                        thresholds_df['atom_pair'] == atom_pair
                    ][q])

                    # Get data for all instances of this atom
                    # pair in the input structure
                    if atom_pair in df.columns.values:
                        df_i = df[df[atom_pair] == True]
                    else:
                        df_i = df[df['atom_pair'] == atom_pair]

                    # Compute the number of clashes surpassing
                    # the threshold
                    n_clashes_dict[f'{atom_pair}_{q}'].append(
                        sum(df_i['d'] < threshold)
                    )

    # Write an output file with the number of clashes
    # per structure
    if args.report_clashes:
        output_file = f'{args.output_file_prefix}n_clashes.csv'
        n_clashes_df = pandas.DataFrame(n_clashes_dict)
        n_clashes_df.to_csv(output_file, index=False)
    
    # Write an output file with inter-atomic distances
    if args.report_distances:
        output_file = f'{args.output_file_prefix}distances.csv'
        df = pandas.concat(dfs)
        df.to_csv(output_file, index=False)

if __name__ == '__main__':
    main()
