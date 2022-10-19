[![Build Status](https://travis-ci.org/abradle/duck.svg?branch=master)](https://travis-ci.org/abradle/duck)
[![stable](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.1.0-blue.svg?style=flat)](https://github.com/abradle/duck)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/abradle/duck/blob/master/LICENSE.txt)

# Installation

## Conda

Make a fresh Conda environment, notice conda might take long time to resolve it I don't know why.
```
git clone git@github.com:AlvaroSmorras/openduck.git
cd duck
conda env create -f environment.yaml 
conda activate openduck_latest
python setup.py install
```

#### Running

Activate conda and run like this:
```
conda activate openduck_latest
python ./scripts/run_full_duck_pond.py -i inputs_file.txt
```

where inputs_file.txt contains the paths to directories with openduck input (1 directory = 1 interaction).
The input directories need to contain:
- prepared apo protein file
- protonated ligand file
- run.yaml file with DUck parameters


where run.yaml is a file like the following:

```
prot_code: '1n2v'
prot_int: 'A_ASP_156_OD2'
lig_id: 'BDI'
cutoff: 9
md_len: 0.5
distance: 2.5
init_velocity: 0.00001
num_smd_cycles: 1
gpu_id: 0
apo_pdb_file: '1n2v_apo.pdb'
mol_file: ligand.mol
```

##### Preparing system for amber

Script to prepare the prot-chunk system using open duck, but to launch the simulations in amber
If you want to keep structural waters, have a file with them called "waters_to_retain.pdb"

```{bash}
conda activate openduck_latest
python ./scripts/duck_prepare_sys_for_amber.py -p receptor.pdb -l ligand.mol -i interaction_string(e.g. A_ASN_140_ND2) -q Slurm|SGE --HMR True -r nreplicas -w wqb_threshold -f SMIRNOFF|GAFF2
```

if you have multiple ligands and want to prepare them in batch. There is also this other script. Notice that the script will generate a folder for each ligand, but only one queue file, that will simulate all the ligands through an array.

```{bash}
conda activate openduck_latest
python ./scripts/batch_duck_prepare_for_amber.py -p receptor.pdb -l ligands.sdf -i interaction_string(e.g. A_ASN_140_ND2) -q Slurm|SGE --HMR -r nreplicas -w wqb_threshold -f SMIRNOFF|GAFF2 -n number_of_threads
```
