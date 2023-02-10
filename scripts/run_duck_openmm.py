import argparse
import pickle
import shutil
import sys
from pathlib import Path
from os import chdir
from os.path import isfile
import yaml


# load functions for duck_steps
from duck.steps.chunk import duck_chunk
from duck.steps.parametrize import prepare_system
from duck.utils.cal_ints import find_interaction, find_atom, is_lig
from duck.steps.equlibrate import do_equlibrate
from duck.utils.check_system import check_if_equlibrated
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md

def args_sanitation(parser):
    args = parser.parse_args()
    # check if everything is ok
    if args.mode == 'full-protocol':
        if (args.yaml_input is None) and (args.ligand is None or args.interaction is None or args.receptor is None):
            parser.error('The input needs to be either the input yaml or specified in the command line (ligand, receptor interaction).')
        elif args.yaml_input:
            input_arguments = yaml.load(open(args.yaml_input), Loader=yaml.FullLoader)
            if all(item in list(input_arguments.keys()) for item in ['receptor_pdb', 'interaction', 'ligand_mol']):
                #transfer all required items
                args.ligand = str(input_arguments['ligand_mol'])
                args.receptor = str(input_arguments['receptor_pdb'])
                args.interaction = str(input_arguments['interaction'])
                
                # overwrite the defaults from command line args
                if 'gpu_id' in input_arguments: args.gpu_id =  input_arguments['gpu_id']
                if 'do_chunk' in input_arguments: args.do_chunk =  bool(input_arguments['do_chunk'])
                if 'cutoff' in input_arguments: args.cutoff =  float(input_arguments['cutoff'])
                if 'ignore_buffers' in input_arguments: args.ignore_buffers =  bool(input_arguments['ignore_buffers'])
                if 'small_molecule_forcefield' in input_arguments: args.small_molecule_forcefield =  str(input_arguments['small_molecule_forcefield'])
                if 'water_model' in input_arguments: args.water_model =  str(input_arguments['water_model'])
                if 'protein_forcefield' in input_arguments: args.protein_forcefield =  str(input_arguments['protein_forcefield'])
                if 'ionic_strength' in input_arguments: args.ionic_strength =  float(input_arguments['ionic_strength'])
                if 'solvent_buffer_distance' in input_arguments: args.solvent_buffer_distance =  float(input_arguments['solvent_buffer_distance'])
                if 'waters_to_retain' in input_arguments: args.waters_to_retain =  str(input_arguments['waters_to_retain'])
                if 'num_smd_cycles' in input_arguments: args.smd_cycles =  int(input_arguments['num_smd_cycles'])
                if 'md_length' in input_arguments: args.md_length =  float(input_arguments['md_length'])
                if 'init_velocities' in input_arguments: args.init_velocities =  float(input_arguments['init_velocities'])
                if 'init_distance' in input_arguments: args.init_distance =  float(input_arguments['init_distance'])
            else:
                parser.error('You need to specify at least "ligand_mol", "receptor_pdb" and "interaction" in the yaml file.')
        elif (args.ligand is None or args.interaction is None or args.receptor is None):
            parser.error('The parameters --ligand, --interaction and --receptor are required.')
        else:
            # all good
            pass
    elif args.mode == 'from-equilibration':
        if (args.yaml_input is None) and (args.input is None or args.pickle is None):
            parser.error('The input needs to be either the input yaml or specified in the command line equilibrated input and pickle from parametrization.')
        elif args.yaml_input:
            input_arguments = yaml.load(open(args.yaml_input), Loader=yaml.FullLoader)
            if all(item in list(input_arguments.keys()) for item in ['pickle', 'equilibrated_system']):
                #transfer required items
                args.pickle = input_arguments['pickle']
                args.equilibrated_system = input_arguments['equilibrated_system']

                # overwrite the defaults from command line args
                if 'gpu_id' in input_arguments: args.gpu_id =  input_arguments['gpu_id']
                if 'num_smd_cycles' in input_arguments: args.smd_cycles =  int(input_arguments['num_smd_cycles'])
                if 'md_length' in input_arguments: args.md_length =  float(input_arguments['md_length'])
                if 'init_velocities' in input_arguments: args.init_velocities =  float(input_arguments['init_velocities'])
                if 'init_distance' in input_arguments: args.init_distance =  float(input_arguments['init_distance'])

    return args

                
                
                



def parse_input():
    parser = argparse.ArgumentParser(description='Run Dynamic Undocking using openMM')
    # run modes full | from equil
    modes = parser.add_subparsers(title='OpenDuck starting mode', help='Modes to run OpenMM, "full-protocol" with or without chunking and "from-equilibrated from a prepared and equilibrated system.')
    full = modes.add_parser('full-protocol', help='OpenDuck OpenMM full protocol either with or without chunking the protein.')
    full.set_defaults(mode='full-protocol')
    full.add_argument('-y', '--yaml-input', type=str, default = None, help='Input yaml file with the all the parameters for the full openDUck protocol.')
    full.add_argument('-l', '--ligand', type=str, default = None, help='ligand mol file to use as reference for interaction.')
    full.add_argument('-i', '--interaction', type=str, default = None, help='Protein atom to use for ligand interaction.')
    full.add_argument('-r', '--receptor', type=str, default = None, help='Protein pdb file to chunk, or chunked protein if mode is "for_chunk".')
    full.add_argument('-g', '--gpu-id', type=int, default=None, help='GPU ID, if not specified, runs on CPU only.')
    # chunking args
    chunk = full.add_argument_group('Chunking parameters')
    chunk.add_argument('--do-chunk', action='store_true', help='Chunk initial receptor based on the interaction with ligand and add cappings.')
    chunk.add_argument('-c', '--cutoff', type=float, default = 9, help='Cutoff distance to define chunk.')
    chunk.add_argument('-b', '--ignore-buffers', action='store_true', help='Do not remove buffers (solvent, ions etc.)')
    # preparation args
    prep = full.add_argument_group('Parametrization parameters')
    prep.add_argument('-f', '--small_molecule_forcefield', type=str,  default = 'SMIRNOFF', choices=('SMIRNOFF', 'GAFF2'), help='Small Molecules forcefield.')
    prep.add_argument('-w', '--water-model', default='tip3p', type=str.lower, choices = ('TIP3P', 'SPCE'), help='Water model to parametrize the solvent with.')
    prep.add_argument('-pf','--protein-forcefield', default='amber99sb', type=str.lower, choices=('amber99sb', 'amber14-all'), help='Protein forcefield to parametrize the chunked protein.')
    prep.add_argument('-ion','--ionic-strength', default=0.1, type=float, help='Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M')
    prep.add_argument('-s','--solvent-buffer-distance', default=10, type=float, help='Buffer distance between the periodic box and the protein. Default = 10 A')
    prep.add_argument('-water','--waters-to-retain', default='waters_to_retain.pdb', type=str, help='PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.')
    # MD/SMD args
    prod = full.add_argument_group('MD/SMD Production parameters')
    prod.add_argument('-n', '--smd_cycles', type=int, default = 20, help='Number of MD/SMD cycles to perfrom')
    prod.add_argument('-m', '--md-length', type=float, default=0.5, help='Lenght of md sampling between smd runs in ns.')
    prod.add_argument('-v', '--init-velocities', type=float, default=0.00001, help='Set initial velocities when heating')
    prod.add_argument('-d', '--init-distance', type=float, default=2.5, help='Set initial HB distance for SMD')
    

    #run from equil
    equil = modes.add_parser('from-equilibrated', help='OpenDuck openMM protocol starting from a pre-equilibrated system (e.g. from duck_prepare_sys.py)')
    equil.add_argument('-y', '--yaml-input', type=str, default = None, help='Input yaml file with the all the parameters for the full openDUck protocol.')
    equil.add_argument('-s', '--equilibrated_system', default=None, help='Equilibrated system as input (*.chk).')
    equil.add_argument('-p', '--pickle', default=None, help='Pickle output from preparation.')
    equil.add_argument('-n', '--smd_cycles', type=int, default = 20, help='Number of MD/SMD cycles to perfrom')
    equil.add_argument('-m', '--md-length', type=float, default=0.5, help='Lenght of md sampling between smd runs in ns.')
    equil.add_argument('-v', '--init-velocities', type=float, default=0.00001, help='Set initial velocities when heating')
    equil.add_argument('-d' '--init-distance', type=float, default=2.5, help='Set initial HB distance for SMD')
    equil.add_argument('-g', '--gpu-id', type=int, default=None, help='GPU ID, if not specified, runs on CPU only.')
    equil.set_defaults(mode='from-equilibrated')

    args = args_sanitation(parser)

    return args
if __name__ == '__main__':
    
    args = parse_input()
    
    