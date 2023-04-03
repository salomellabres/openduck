[![stable](http://badges.github.io/stability-badges/dist/experimental.svg)](http://github.com/badges/stability-badges)
[![Version](http://img.shields.io/badge/version-0.1.0-blue.svg?style=flat)](https://github.com/CBDD/openduck)
[![License](http://img.shields.io/badge/license-Apache%202.0-blue.svg?style=flat)](https://github.com/CBDD/openduck/LICENSE)

# Welcome to OpenDUck

OpenDUck is an opensource implementation of Dynamic Undocking. Dynamic Undocking (DUck) is a particular implementation of steered molecular dynamics developed as a fast algorithm to calculate the work necessary to reach a quasi-bound state at which a certain netive contact is broken between the ligand and receptor. Due to its cheap and fast nature, DUck is mainly employed as a post-docking filter during virtual screening campaings.

Ruiz-Carmona, S., Schmidtke, P., Luque, F. et al. Dynamic undocking and the quasi-bound state as tools for drug discovery. Nature Chem 9, 201–206 (2017). https://doi.org/10.1038/nchem.2660 

## Installation

### Conda

Make a fresh Conda environment, notice conda might take long time to resolve it I don't know why.
```{bash}
$ git clone git@github.com:CBDD/openduck.git
$ cd openduck
$ conda env create -f environment.yaml 
$ conda activate openduck
$ python setup.py install
```

## OpenDUck usage

Openduck can be used as a python library or as an executable. The executable provided is openduck, which is set up an entry point in the path. Alternatively you can use the openduck.py in scripts.

Openduck comes with several submodules to launch the different steps of Dynamic Undocking

```{bash}
$ conda activate openduck

$ openduck -h

usage: openduck [-h]  ...

Open Dynamic Undocking

optional arguments:
  -h, --help            show this help message and exit

Open Dynamic Undocking toolkit. Choose one of the following actions:
  
    OpenMM_prepare      
                        Preparation of systems for OpenMM simulations
    OpenMM_full-protocol
                        OpenDuck OpenMM full protocol either with or without chunking the protein.
    OpenMM_from-equilibrated
                        OpenDuck openMM protocol starting from a pre-equilibrated system (e.g. from duck_prepare_sys.py)
    OpenMM_from-amber   
                        OpenDuck openMM protocol starting from an amber topology and coordinates (prmtop and inpcrd).
    AMBER_prepare       
                        Preparation of systems, inputs and queue files for AMBER simulations
    report              
                        Generate report for openduck results.
    chunk               
                        Chunk a protein for Dynamic Undocking.
```

Each of the OpenDUck modules accepts the input either as command line arguments or in a yaml file.
The arguments in the input.yaml are expected to follow the command-line nomenclature.

### Running OpenDUck in OpenMM

The openMM implementation of OpenDUck is divided in two steps, the preparation using openmm_forcefields and the production (MD/SMD) following the Dynamic Undocking protocol.
During the preparation step, the protein is reduced to the pocket during the _chunking_ step, then the ligand and receptor are parametrized with the specified forcefields and the solvation box is generated.
After equilibration, the sequential MD/SMD simulations, sampling and steering the main interaction from the specified distances, are produced in OpenMM.

The openMM submodules are setup to allow independent usage of the previously mentioned steps:

    * OpenMM_full-protocol: Executes the full OpenDuck protocol, from the protein-ligand files to the Steered simulations
    * OpenMM_prepare: Executes exclusively the preparation and equilibration of the protein-ligand complex
    * OpenMM_from-equilibrated: Executes exclusively the production of MD/SMD runs from an equilibrated input.
    * OpenMM_from-amber: Starts the openDuck simulations from an amber topology.

#### Example full-protocol

The complete OpenDUck protocol through OpenMM can be executed using the OpenMM_full-protocol submodule. The usage of the module can be obtained through the help.

```{bash}
$ openduck OpenMM_full-protocol -h

usage: openduck OpenMM_full-protocol [-h] [-y YAML_INPUT] [-l LIGAND] [-i INTERACTION] [-r RECEPTOR] [-g GPU_ID] [--do-chunk] [-c CUTOFF] [-b] [-f {SMIRNOFF,GAFF2}]
                                     [-w {tip3p,spce}] [-ff {amber99sb,amber14-all}] [-ion IONIC_STRENGTH] [-s SOLVENT_BUFFER_DISTANCE] [-water WATERS_TO_RETAIN] [-fl]
                                     [-F FORCE_CONSTANT_EQ] [-n SMD_CYCLES] [-m MD_LENGTH] [-W WQB_THRESHOLD] [-v INIT_VELOCITIES] [-d INIT_DISTANCE]

optional arguments:
  -h, --help            show this help message and exit

Main arguments:
  -y YAML_INPUT, --yaml-input YAML_INPUT
                        Input yaml file with the all the arguments for the full openDUck protocol.
  -l LIGAND, --ligand LIGAND
                        ligand mol file to use as reference for interaction.
  -i INTERACTION, --interaction INTERACTION
                        Protein atom to use for ligand interaction.
  -r RECEPTOR, --receptor RECEPTOR
                        Protein pdb file to chunk, or chunked protein if mode is "for_chunk".
  -g GPU_ID, --gpu-id GPU_ID
                        GPU ID, if not specified, runs on CPU only.

Chunking arguments:
  --do-chunk            Chunk initial receptor based on the interaction with ligand and add cappings.
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff distance to define chunk.
  -b, --ignore-buffers  Do not remove buffers (solvent, ions etc.)

Parametrization arguments:
  -f {SMIRNOFF,GAFF2}, --small_molecule_forcefield {SMIRNOFF,GAFF2}
                        Small Molecules forcefield.
  -w {tip3p,spce}, --water-model {tip3p,spce}
                        Water model to parametrize the solvent with.
  -ff {amber99sb,amber14-all}, --protein-forcefield {amber99sb,amber14-all}
                        Protein forcefield to parametrize the chunked protein.
  -ion IONIC_STRENGTH, --ionic-strength IONIC_STRENGTH
                        Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M
  -s SOLVENT_BUFFER_DISTANCE, --solvent-buffer-distance SOLVENT_BUFFER_DISTANCE
                        Buffer distance between the periodic box and the protein. Default = 10 A
  -water WATERS_TO_RETAIN, --waters-to-retain WATERS_TO_RETAIN
                        PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.
  -fl, --fix-ligand     Some simple fixes for the ligand: ensure tetravalent nitrogens have the right charge assigned and add hydrogens to carbon atoms.

MD/SMD Production arguments:
  -F FORCE_CONSTANT_EQ, --force-constant_eq FORCE_CONSTANT_EQ
                        Force Constant for equilibration
  -n SMD_CYCLES, --smd-cycles SMD_CYCLES
                        Number of MD/SMD cycles to perfrom
  -m MD_LENGTH, --md-length MD_LENGTH
                        Lenght of md sampling between smd runs in ns.
  -W WQB_THRESHOLD, --wqb-threshold WQB_THRESHOLD
                        Minimum WQB threshold to stop simulations.
  -v INIT_VELOCITIES, --init-velocities INIT_VELOCITIES
                        Set initial velocities when heating
  -d INIT_DISTANCE, --init-distance INIT_DISTANCE
                        Set initial HB distance for SMD
```

A valid example is provided in the test subfolder, and can be executed using the command-line arguments like the following:

```{bash}
$ openduck OpenMM_full-protocol -l brd4_lig.mol -r 4LR6_aligned_chunk_nowat.pdb -i A_ASN_140_ND2 -g 0 -f GAFF2 -ff amber14-all -w tip3p -w 6 -n 20
```

Similarly, an execution with the input parameters in a yaml is also possible (the example has slightly different parameters as the command line, just for ilustration purposes ).
```{bash}
openduck OpenMM_full-protocol -y input.yaml
```
where the _input.yaml_ file has the following
```{yaml}
# Main Arguments
interaction : A_ASN_140_ND2
receptor_pdb : 4LR6_aligned_chunk_nowat.pdb
ligand_mol : brd4_lig.mol
gpu_id : 0

# Chunking Arguments
do_chunk : True
cutoff : 10
ignore_buffers : False

# Preparation Arguments
small_molecule_forcefield : smirnoff
protein_forcefield : amber14-all
water_model : tip3p
ionic_strength : 0.05
waters_to_retain : waters_to_retain.pdb
solvent_buffer_distance : 15
force_constant_eq : 1

# Production Arguments
smd_cycles : 20
md_length : 0.5
wqb_threshold : 6
init_velocities : 0.00001
init_distance : 2.5
fix_ligand : False
```

### Running OpenDUck in Amber

Originally Dynamic Undocking was designed to run in Amber, in OpenDUck we have maintained the same protocol for AMBER, replacing the parametrization and file-generation previously done by MOE with an open implementation. Similarly as with the old protocol, freedom on the quequeing headers and preparation are integrated both in the commandline arguments and in the input yaml.
Preparation flags are very similar to those for OpenMM, as most of the implementation is shared.

```{bash}
$ openduck AMBER_prepare -l ../1a28_lig.mol -r ../1a28_prot.pdb -i A_ARG_766_NH2 --do-chunk -c 10 -b False -f GAFF2 -ff amber99sb -w SPCE -ion 0.05 -s 11 --HMR -n 4 -W 4 -q Slurm
```

Tipically, we Dynamic Undocking is employed as a post-docking filter in a virtual screening campaign. As such, multiple ligands might be generated at the time. A batch execution of the DUck preparation can be executed specifying the _--batch_ argument and giving an sdf file with multiple ligands as input. Each ligand and the associated file structure will be generated in the LIG_target_{n} subfolder where n is the ligand's position in the sdf.

```{bash}
$ openduck AMBER_prepare -l ../brd4_ligands.sdf -r ../4LR6_aligned_chunk_nowat.pdb -i A_ASN_140_ND2 --waters-to-retain ../waters_to_retain.pdb -f gaff2 -ff amber14-all -w spce -ion 1 -s 30 --HMR --smd-cycles 10 -wqb-threshold 8 -q SGE --seed 1235467890 --batch --threads 8
```

The queueing template is stored in duck/templates/queueing_templates. To customize a template, a new file can be added in the directory with the following format: {queue_name}_array.q or {queue_name}.q either if the expected execution is in an array or not. For a local execution, the _local_ argument can be given for a plain list of commands.

The templates have the following format:
```
#!/bin/bash
#HEADER

{functions}

{commands}
```

### Analysis

After the successfully running openDUck in either openMM or Amber, compiling the results can be performed using the report submodule.
By default, the lowest $W_{QB}$ is used as a reporter for the robustness of the studied interaction, however, a better (and more formally appropiate) descriptor has been shown to be the $\Delta_{QB}$ obtained using the Jarzynski equality. When specified, a bootstrapping is performed to assess the convergence of the reported value.

**important**

Due the different formating of Amber and OpenMM steered molecular dynamics output, the _--format_ argument is needed when analyzing either output.
When analyzing a single result, the default pattern is the current directory (.). If multiple directories are ment to be analysed, specify the wildcard pattern.

```{bash}
$ openduck report -h

usage: openduck report [-h] [-p PATTERN] [-d {min,single,avg,jarzynski,all}] [-o OUTPUT] [-of {csv,sdf,tbl}] [--plot] [-s SUBSAMPLE_SIZE] [-i ITERATIONS] [-t STEP_THRESHOLD] [-f {amber,openmm}]

optional arguments:
  -h, --help            show this help message and exit
  -p PATTERN, --pattern PATTERN
                        Wildcard pattern to find folders with DUck data
  -d {min,single,avg,jarzynski,all}, --data {min,single,avg,jarzynski,all}
                        Mode to compile the report [min | single | avg | jarzynski | all]
  -o OUTPUT, --output OUTPUT
                        Output file, default is printing report to stdout.
  -of {csv,sdf,tbl}, --output-format {csv,sdf,tbl}
                        Output format, [csv | sdf | tbl].
  --plot                Plot work or energy values to file.
  -s SUBSAMPLE_SIZE, --subsample-size SUBSAMPLE_SIZE
                        Subsample size for jarzynski bootstrapping.
  -i ITERATIONS, --iterations ITERATIONS
                        Number of bootstrapping iterations for jarzynski analysis.
  -t STEP_THRESHOLD, --step-threshold STEP_THRESHOLD
                        Steps_treshold to find the minima.
  -f {amber,openmm}, --format {amber,openmm}
                        Data where the results come from. Amber by default.
```

#### Report Example

Parsing $W_{QB}$ in a tabulted output format from an example from the benchmark dataset.
```{bash}
$ cd Iridium/1a28/ARG766_HN2
$ openduck report -p '*/'

System	WQB
SMIRNOFF/	5.6297999999999995
GAFF2/	5.74429
SMIRNOFF_nohmr/	6.922040000000001
GAFF2_nohmr/	6.264200000000001
GAFF2_bigbox/	5.180540000000001
GAFF2_tip4pew/	7.60925
GAFF2_spce/	6.17347
```

Calculating and reporting the $\Delta_{QB}$ of the same simulations as before
```{bash}
$ openduck report -p '*/' -d jarzynski

System	Jarzynski	Jarzynski_SD	Jarzynski_SEM
SMIRNOFF/	6.391165243930965	0.12701922010484387	0.02033935321315065
GAFF2/  6.990007044361647	0.23871054814271678	0.038224279367893604
SMIRNOFF_nohmr/	7.690209327251983	0.10126872662893496	0.01621597423328342
GAFF2_nohmr/    7.075924262571702   0.10628038471267455 0.01701848178973498
GAFF2_bigbox/   6.670905906961162   0.3167545344338863  0.05072131880828811
GAFF2_tip4pew/  9.007291835327141   0.33968295997031067 0.05439280525909316
GAFF2_spce/ 7.31123937991228    0.24521160219077848 0.03926528115039682

``` 