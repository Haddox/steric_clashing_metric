# Quantify steric clashing between atoms

This repository uses `PyRosetta` to compute the level of clashing between atoms in a set of input molecules.

## How the metric works

The code in this repo computes the number of "bad" clashes per input PDB, quantified separately for different atom types. "Bad" clashes are defined as those that exceed thresholds of clashing seen in native proteins. Specifically, for a given atom pair (e.g., methyl carbons and backbone oxygens), I pre-computed the distribution of inter-atomic distances seen in a set of high-resolution crystal structures of native proteins, exluding extremely long distances. I next computed the shortest inter-atomic distances (i.e., biggest clashes) in the distribution, drawing thresholds that corresponded to the 1st and 5th percentile. I computed thresholds for all combinations of atom types with at least 100 counts in the input natives.

For each input structure, the code computes the number of clashes that surpass each of the thresholds defined above. For instance, an input protein may have three pairs of methyl carbons and backbone oxygens that have shorter distances (bigger clashes) than the 1st percentile of distances in natives. The script reports this kind of number for each combination of atoms and each threshold.

The clashing metrics that I would suggest looking at first are:

* `hydrophobic_C_C_q1`: the number of carbon-carbon interactions between hydrophobic sidechains where the inter-atomic distance is less than the 1st percentile of distances seen in natives.
* `hydrophobic_C_H_q1`: same as above, but for carbon-hydrogen interactions between hydrophobic sidechains.
* `hydrophobic_H_H_q1`: same as above, but for hydrogen-hydrogen interactions between hydrophobic sidechains.
* `C_Obb_q1`: same as above, but for interactions between a carbon from a hydrophobic sidechain and a backbone oxygen.

I typically see the largest numbers of clashes in these categories of atoms. But, the code reports clashing between all pairs of atom types observed in the set of native proteins. Below, I describe in more detail how to run the code and what the output files report.

## How to run the code

**Example command:** Below is a simple example of how to run the code on a set of input proteins. Each flag is described below, along with a set of optional flags. The output of the code is also described below.

```
python scripts/compute_distances.py --pdb_dir PDB_DIR --file_with_thresholds results/natives/thresholds.csv --output_file_prefix OUTPUT_FILE_PREFIX
```

**Required flags:**
* `--pdb_dir`: the path to a directory with input PDBs
* `--file_with_thresholds`: the path to a file with thresholds used to determine if a pair of atoms is clashing. The file `results/natives/thresholds.csv` has thresholds computed from a set of ~80 high-resolution crystal structures of native proteins
* `--output_file_prefix`: a prefix that that will be used as the start of the path of output files. If `--report_clashes` is set to True (default), the script will output a file with the suffix `n_clashes.csv` that quantifies the number of clashes that surpass thresholds from the above input file. If `--report_distances` is set to True (default is False), then the script will also output a file with all inter-atomic distances used to compute the number of clashes.

**Optional flags:**
* `--report_clashes`: a boolean (default: True). Write an output file with the number of clashes per structure surpassing thresholds. The file will have the suffix `n_clashes.csv`.
* `--report_distances`: a boolean (default: False). Write an output file with inter-atomic distances used to compute the number of clashes. The file will have the suffix `distances.csv`.
* `--use_tenA_neighbor_residues`: a boolean (default: True). Only compute distances between atoms from residues with C-beta atoms within 10 Angstroms. This helps with speed and appears to capture nearly all clashes.

**How I pre-computed clashing thresholds from native proteins**
The Jupyter notebook `analysis_code.ipynb` has the code I used for this step.

## Output files

* `n_clashes.csv`: A dataframe with one row per protein and columns with the following data:
    * `pdb`: the path to an input PDB
    * `X_q1` and `X_q5`: there are several columns that follow this naming convention, where `X` specifies a specific combination of atoms (e.g., hydrophobic carbons), and `q1` and `q5` specify which threshold was used to compute the number of clashes.
        * `X` = `hydrophobic_C_C`: the atoms are carbons from hydrophobic sidechains (ALIVMFYW).
        * `X` = `hydrophobic_C_H`: the atoms are one carbon and one hydrogen, both from hydrophobic sidechains (ALIVMFYW).
        * `X` = `hydrophobic_H_H`: the atoms are carbons from hydrophobic sidechains (ALIVMFYW).
        * `X` = `C_Obb`: the atoms are one carbon from a hydrophobic sidechain (ALIVMFYW) and one backbone oxygen.
        * `q1`: threshold corresponds to the 1st percentile of the distance distribution from native proteins.
        * `q5`: threshold corresponds to the 1st percentile of the distance distribution from native proteins.
    * `X:Y_q1` and `X:Y_q5`: there are several other columns that follow this naming convention, where `X` and `Y` specify the two Rosetta atom types involved in the clash (e.g., `CH3:OCbb`), and where `q1` and `q5`specify which threshold was used to compute the number of clashes.

* `distances.csv`: This file is only generated if the `--report_distances` flag is set to True. It is a dataframe with one row per atom pair, consisting of all inter-atomic distances used to compute the number of clashes in `n_clashes.csv`. This dataframe contains many columns. Here are a few important ones:
    * `pdb`: the path to the input PDB with the atom pair
    * `d`: the distance between the atom pair
    * `o`: the sum of the atomic radii of the two atoms, according to Rosetta.
    * there are several other columns that provide the name and number of the pair of atoms and residues involved, where `i` and `j` are used to distinguish between atoms/residues on either side of the clash.