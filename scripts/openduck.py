import argparse
import pickle
import shutil
import sys
from pathlib import Path
import os
import yaml

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
                if 'force_constant_eq' in input_arguments: args.force_constant_eq =  float(input_arguments['force_constant_eq'])

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
            else:
                parser.error('You need to specify at least "pickle" and "equilibrated_system" in the yaml file.')
        elif (args.pickle is None or args.equilibrated_system is None):
            parser.error('The parameters --pickle and --equilibrated_system are required.')
        else:
            #all good
            pass
    elif args.mode == 'Amber-prep':
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
                if 'do_chunk' in input_arguments: args.do_chunk =  bool(input_arguments['do_chunk'])
                if 'cutoff' in input_arguments: args.cutoff =  float(input_arguments['cutoff'])
                if 'ignore_buffers' in input_arguments: args.ignore_buffers =  bool(input_arguments['ignore_buffers'])
                if 'small_molecule_forcefield' in input_arguments: args.small_molecule_forcefield =  str(input_arguments['small_molecule_forcefield'])
                if 'water_model' in input_arguments: args.water_model =  str(input_arguments['water_model'])
                if 'protein_forcefield' in input_arguments: args.protein_forcefield =  str(input_arguments['protein_forcefield'])
                if 'ionic_strength' in input_arguments: args.ionic_strength =  float(input_arguments['ionic_strength'])
                if 'solvent_buffer_distance' in input_arguments: args.solvent_buffer_distance =  float(input_arguments['solvent_buffer_distance'])
                if 'waters_to_retain' in input_arguments: args.waters_to_retain =  str(input_arguments['waters_to_retain'])
                if 'seed' in input_arguments: args.seed =  str(input_arguments['seed'])
                if 'queue_template' in input_arguments: args.queue_template = str(input_arguments['queue_template'])
                if 'HMR' in input_arguments: args.HMR = bool(input_arguments['HMR'])
                if 'wqb_threshold' in input_arguments: args.wqb_threshold = bool(input_arguments['wqb_threshold'])
                if 'smd_cycles' in input_arguments: args.smd_cycles = bool(input_arguments['smd_cycles'])
                if 'batch' in input_arguments: args.batch = int(input_arguments['batch'])
                if 'number_of_cpus' in input_arguments: args.number_of_cpus = int(input_arguments['number_of_cpus'])

            else:
                parser.error('You need to specify at least "ligand_mol", "receptor_pdb" and "interaction" in the yaml file.')
        elif (args.ligand is None or args.interaction is None or args.receptor is None):
            parser.error('The parameters --ligand, --interaction and --receptor are required.')
        else:
            # all good
            pass

    return args

def parse_input():
    parser = argparse.ArgumentParser(description='Open Dynamic Undocking')
    # run modes full | from equil
    modes = parser.add_subparsers(title='OpenDuck starting mode', help='Modes to run OpenDuck, "full-protocol" with or without chunking and "from-equilibrated from a prepared and equilibrated system.')
    full = modes.add_parser('OpenMM_full-protocol', help='OpenDuck OpenMM full protocol either with or without chunking the protein.')
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
    prep.add_argument('-ff','--protein-forcefield', default='amber99sb', type=str.lower, choices=('amber99sb', 'amber14-all'), help='Protein forcefield to parametrize the chunked protein.')
    prep.add_argument('-ion','--ionic-strength', default=0.1, type=float, help='Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M')
    prep.add_argument('-s','--solvent-buffer-distance', default=10, type=float, help='Buffer distance between the periodic box and the protein. Default = 10 A')
    prep.add_argument('-water','--waters-to-retain', default='waters_to_retain.pdb', type=str, help='PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.')
    # MD/SMD args
    prod = full.add_argument_group('MD/SMD Production parameters')
    prod.add_argument('-F', '--force-constant_eq', type=float, default = 1, help='Force Constant for equilibration')
    prod.add_argument('-n', '--smd-cycles', type=int, default = 20, help='Number of MD/SMD cycles to perfrom')
    prod.add_argument('-m', '--md-length', type=float, default=0.5, help='Lenght of md sampling between smd runs in ns.')
    prod.add_argument('-v', '--init-velocities', type=float, default=0.00001, help='Set initial velocities when heating')
    prod.add_argument('-d', '--init-distance', type=float, default=2.5, help='Set initial HB distance for SMD')
    
    #run from equil
    equil = modes.add_parser('OpenMM_from-equilibrated', help='OpenDuck openMM protocol starting from a pre-equilibrated system (e.g. from duck_prepare_sys.py)')
    equil.add_argument('-y', '--yaml-input', type=str, default = None, help='Input yaml file with the all the parameters for the openMM simulations from the equilibrated system.')
    equil.add_argument('-s', '--equilibrated-system', default=None, help='Equilibrated system as input (*.chk).')
    equil.add_argument('-p', '--pickle', default=None, help='Pickle output from preparation.')
    equil.add_argument('-n', '--smd-cycles', type=int, default = 20, help='Number of MD/SMD cycles to perfrom')
    equil.add_argument('-m', '--md-length', type=float, default=0.5, help='Lenght of md sampling between smd runs in ns.')
    equil.add_argument('-v', '--init-velocities', type=float, default=0.00001, help='Set initial velocities when heating')
    equil.add_argument('-d' '--init-distance', type=float, default=2.5, help='Set initial HB distance for SMD')
    equil.add_argument('-g', '--gpu-id', type=int, default=None, help='GPU ID, if not specified, runs on CPU only.')
    equil.set_defaults(mode='from-equilibrated')

    #Preparation for amber
    amber = modes.add_parser('preparation_for_AMBER', help='Preparation of systems, inputs and queue files for AMBER simulations')
    amber.add_argument('-y', '--yaml-input', type=str, default = None, help='Input yaml file with the all the parameters for the system preparation and inputs/queueing for AMBER.')
    amber.add_argument('-l', '--ligand', type=str, default = None, help='ligand mol file to use as reference for interaction.')
    amber.add_argument('-i', '--interaction', type=str, default = None, help='Protein atom to use for ligand interaction.')
    amber.add_argument('-r', '--receptor', type=str, default = None, help='Protein pdb file to chunk, or chunked protein if mode is "for_chunk".')
    equil.set_defaults(mode='Amber-preparation')

    # Chunking args
    amber_chunk = amber.add_argument_group('Chunking parameters')
    amber_chunk.add_argument('--do-chunk', action='store_true', help='Chunk initial receptor based on the interaction with ligand and add cappings.')
    amber_chunk.add_argument('-c', '--cutoff', type=float, default = 9, help='Cutoff distance to define chunk.')
    amber_chunk.add_argument('-b', '--ignore-buffers', action='store_true', help='Do not remove buffers (solvent, ions etc.)')

    # preparation args
    amber_prep = amber.add_argument_group('Parametrization parameters')
    amber_prep.add_argument('-f', '--small_molecule_forcefield', type=str,  default = 'SMIRNOFF', choices=('SMIRNOFF', 'GAFF2'), help='Small Molecules forcefield.')
    amber_prep.add_argument('-w', '--water-model', default='tip3p', type=str.lower, choices = ('TIP3P', 'SPCE'), help='Water model to parametrize the solvent with.')
    amber_prep.add_argument('-q', '--queue-template', type=str, default = 'local', help='Write out a queue file from templates.')
    amber_prep.add_argument('-H','--HMR', action='store_true', help ='Perform Hydrogen Mass Repartition on the topology and use it for the input files')
    amber_prep.add_argument('-n', '--smd-cycles', type=int, default=5, help='Ammount of SMD replicas to perform')
    amber_prep.add_argument('-w', '--wqb-threshold', type=float, default=7.0, help='WQB threshold to stop the simulations')
    amber_prep.add_argument('-ff','--protein-forcefield', default='amber99sb', type=str.lower, choices=('amber99sb', 'amber14-all'), help='Protein forcefield to parametrize the chunked protein.')
    amber_prep.add_argument('-ion','--ionic-strength', default=0.1, type=float, help='Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M')
    amber_prep.add_argument('-s','--solvent-buffer-distance', default=10, type=float, help='Buffer distance between the periodic box and the protein. Default = 10 A')
    amber_prep.add_argument('-water','--waters-to-retain', default='waters_to_retain.pdb', type=str, help='PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.')
    amber_prep.add_argument('--seed', default='-1', type=str, help='Specify seed for amber inputs.')
    amber_prep.add_argument('-B', '--batch', default=False, action='store_true', help='Batch processing for multi-ligand sdf')
    amber_prep.add_argument('-ncpu', '--number-of-cpus', default=1, type=int, help='Define number of cpus for batch processing.')

    args = args_sanitation(parser)

    return args

def duck_smd_runs(input_checkpoint, pickle, num_runs, md_len, gpu_id, start_dist, init_velocity, save_dir):
    from duck.steps.normal_md import perform_md
    from duck.steps.steered_md import run_steered_md
    from duck.utils.check_system import check_if_equlibrated

    
    shutil.copyfile(input_checkpoint, "equil.chk")
    shutil.copyfile(pickle, "complex_system.pickle")

    # Now do the MD
    # remember start_dist
    if not Path(save_dir).exists(): save_dir.mkdir()
    for i in range(num_runs):
        if i == 0:
            md_start = "equil.chk"
        else:
            md_start = str(Path(save_dir, "md_" + str(i - 1) + ".chk"))
        log_file = str(Path(save_dir, "md_" + str(i) + ".csv"))
        perform_md(
            checkpoint_in_file=md_start,
            checkpoint_out_file=str(Path(save_dir, "md_" + str(i) + ".chk")),
            csv_out_file=log_file,
            pdb_out_file=str(Path(save_dir, "md_" + str(i) + ".pdb")),
            dcd_out_file=str(Path(save_dir, "md_" + str(i) + ".dcd")),
            md_len=md_len,
            gpu_id=gpu_id,
        )
        # Open the file and check that the potential is stable and negative
        if not check_if_equlibrated(log_file, 3):
            print("SYSTEM NOT EQUILIBRATED")
            sys.exit()
        # Now find the interaction and save to a file


        run_steered_md(
            300,
            str(Path(save_dir, "md_" + str(i) + ".chk")),
            str(Path(save_dir, "smd_" + str(i) + "_300.csv")),
            str(Path(save_dir, "smd_" + str(i) + "_300.dat")),
            str(Path(save_dir, "smd_" + str(i) + "_300.pdb")),
            str(Path(save_dir, "smd_" + str(i) + "_300.dcd")),
            start_dist,
            init_velocity=init_velocity,
            gpu_id=gpu_id,
        )
        run_steered_md(
            325,
            str(Path(save_dir, "md_" + str(i) + ".chk")),
            str(Path(save_dir, "smd_" + str(i) + "_325.csv")),
            str(Path(save_dir, "smd_" + str(i) + "_325.dat")),
            str(Path(save_dir, "smd_" + str(i) + "_325.pdb")),
            str(Path(save_dir, "smd_" + str(i) + "_325.dcd")),
            start_dist,
            init_velocity=init_velocity,
            gpu_id=gpu_id,
        )

def prepare_sys_for_amber(ligand_file, protein_file, chunk_file, interaction, HMR,  small_molecule_forcefield='SMIRNOFF', water_ff_str = 'tip3p.xml', forcefield_str='amber99sb.xml', ionic_strength = 0.1, box_buffer_distance = 10, waters_to_retain="waters_to_retain.pdb", seed='-1'):
    from duck.steps.parametrize import prepare_system
    from duck.utils.cal_ints import find_interaction
    from duck.utils.amber_inputs import Amber_templates
    # Parameterize the ligand
    prepare_system(ligand_file, chunk_file, forcefield_str=forcefield_str,
                   hmr=HMR, small_molecule_ff=small_molecule_forcefield, water_ff_str = water_ff_str,
                   box_buffer_distance = box_buffer_distance, ionicStrength = ionic_strength, waters_to_retain="waters_to_retain.pdb")
    
    # Now find the interaction and save to a file
    results = find_interaction(interaction, protein_file)
    print(results) # what happens to these?
    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + results
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)
    #p = (parmed_structure, prot_index, ligand_index, pairmeandistance)
    p[0].save('system_complex.inpcrd', overwrite=True)

    
    #do_equlibrate(force_constant_equilibrate=args.force_constant_eq, gpu_id=args.gpu_id, keyInteraction=p[1:])
    amber = Amber_templates(structure=p[0], interaction=p[1:],hmr=HMR, seed=seed)
    amber.write_all_inputs()

def AMBER_prepare_ligand_in_folder(ligand_string, lig_indx, protein, chunk, interaction, HMR, base_dir, small_molecule_forcefield = 'SMIRNOFF', water_model = 'tip3p', forcefield = 'amber99sb', ion_strength = 0.1, box_buffer_distance = 10, waters_to_retain='waters_to_retain.pdb', seed='-1'):
    from duck.utils.amber_inputs import write_string_to_file
    from contextlib import redirect_stdout,redirect_stderr

    os.chdir(base_dir)
    
    #Create the ligand folder
    if os.path.isdir(f'LIG_target_{lig_indx}'):
        print(f'WARNING: LIG_target_{lig_indx} already exist and it will be overwritten.')
        shutil.rmtree(f'./LIG_target_{lig_indx}', ignore_errors=True)
    os.mkdir(f'LIG_target_{lig_indx}')
    os.chdir(f'LIG_target_{lig_indx}')
    print(f'Working on LIG_target_{lig_indx}')
    with open('preparation.out', 'w') as o:
        with redirect_stdout(o):
            # Copying files to ligand foldef; ligand and prot
            write_string_to_file(string=ligand_string, file=f'lig_{lig_indx}.mol')
            shutil.copyfile(f'../{protein}', f'./{protein}', follow_symlinks=True)
            if os.path.isfile(f'../{waters_to_retain}'):
                shutil.copyfile(f'../{waters_to_retain}', f'./{waters_to_retain}', follow_symlinks=True)

            prepare_sys_for_amber(f'lig_{lig_indx}.mol', protein, chunk, interaction, HMR,
                                  small_molecule_forcefield=small_molecule_forcefield, water_ff_str=f'{water_model}',
                                  forcefield_str=f'{forcefield}.xml', ion_strenght = ion_strength,
                                  box_buffer_distance = box_buffer_distance, waters_to_retain=f"{waters_to_retain}", seed=seed)

    #os.chdir(f'..')
    return(f'Lig_target_{lig_indx} prepared correctly')

#### main functions
def do_full_openMM_protocol(args):
    # adapted from run_full_duck_pipeline.py
    from duck.utils.check_system import check_if_equlibrated
    from duck.steps.equlibrate import do_equlibrate
    from duck.steps.parametrize import prepare_system
    from duck.utils.cal_ints import find_interaction
    # create chunk
    if args.do_chunk:
        print('Chunking protein')
        from duck.steps.chunk import duck_chunk
        chunked_file = duck_chunk(args.receptor,args.ligand,args.interaction,args.cutoff, ignore_buffers=args.ignore_buffers)
    else: chunked_file = args.protein
    # prepare system
    prepare_system(args.ligand, chunked_file, forcefield_str=f'{args.protein_forcefield}.xml', water_ff_str = f'{args.water_model}',
            small_molecule_ff=args.small_molecule_forcefield, waters_to_retain=args.waters_to_retain,
            box_buffer_distance = args.solvent_buffer_distance, ionicStrength = args.ionic_strength)
    results = find_interaction(args.interaction, args.protein)
    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + results
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)
    p[0].save('system_complex.inpcrd', overwrite=True)

    # Equlibration
    do_equlibrate(force_constant_equilibrate=args.force_constant_eq, gpu_id=args.gpu_id, keyInteraction=results)
    if not check_if_equlibrated("density.csv", 1):
        raise EquilibrationError("System is not equilibrated.") # Does this exist?
    
    # set up phase, I don't know why are the names changed here. Might be better to ommit it
    pickle_path = Path('complex_system.pickle')
    #new_pickle_path = Path('cs.pickle')
    #pickle_path.rename(new_pickle_path)
    equil_path = Path('equil.chk')
    #new_equil_path = Path('eql.chk')
    #equil_path.rename(new_equil_path)
    #print('checkpoint_path', equil_path)
    save_dir = Path('duck_runs')
    if not save_dir.exists(): save_dir.mkdir()

    # Now production
    duck_smd_runs(input_checkpoint=equil_path,
                pickle=pickle_path,
                num_runs=args.smd_cycles,
                md_len=args.md_length,
                gpu_id=args.gpu_id,
                start_dist=args.init_distance,
                init_velocity=args.init_velocity,
                save_dir=save_dir)

def do_openMM_from_equil(args):
    save_dir = Path('duck_runs')
    if not save_dir.exists(): save_dir.mkdir()
    #only need to do production
    duck_smd_runs(input_checkpoint=args.equilibrated_system,
            pickle=args.pickle,
            num_runs=args.smd_cycles,
            md_len=args.md_length,
            gpu_id=args.gpu_id,
            start_dist=args.init_distance,
            init_velocity=args.init_velocity,
            save_dir=save_dir)

def do_AMBER_preparation(args):
    # adapted from duck_prepare_sys_for_amber.py and batch_duck_prepare_for_amber.py
    from duck.utils.amber_inputs import Queue_templates
     
    # create chunk if needed
    if args.do_chunk:
        from duck.steps.chunk import duck_chunk
        print('Chunking protein')
        chunk_file = duck_chunk(args.receptor,args.ligand,args.interaction,args.cutoff, ignore_buffers=args.ignore_buffers)
    else: chunk_file = args.protein

    if args.batch:
        import multiprocessing as mp
        from duck.utils.amber_inputs import log_result, handle_error, ligand_string_generator
        pool = mp.Pool(args.number_of_cpus)
        base_dir = os.getcwd()

        r = [pool.apply_async(AMBER_prepare_ligand_in_folder,
                          args=(ligand_string, j+1, args.protein, chunk_file,
                                args.interaction, args.HMR, base_dir,
                                args.small_molecule_forcefield, args.water_model, args.protein_forcefield,
                                args.ionic_strength, args.solvent_buffer_distance, args.waters_to_retain, args.seed),
                          callback=log_result,
                          error_callback=handle_error) for j, ligand_string in enumerate(ligand_string_generator(args.ligands))]
        pool.close()
        pool.join()

        queue = Queue_templates(wqb_threshold=args.wqb_threshold, replicas=args.replicas, array_limit=len(r), hmr=args.HMR)
    else:
        if args.number_of_cpus != 1:
            print('WARNING: The number of threads does not have an impact if the batch mode is not enabled.')
        prepare_sys_for_amber(args.ligand, args.protein, chunk_file, args.interaction, args.HMR, 
        small_molecule_forcefield=args.small_molecule_forcefield, water_ff_str = args.water_model,
        forcefield_str=f'{args.protein_forcefield}.xml', ionic_strength = args.ionic_strength,
        box_buffer_distance = args.solvent_buffer_distance, waters_to_retain=args.waters_to_retain, seed=args.seed)
        queue = Queue_templates(wqb_threshold=args.wqb_threshold, replicas=args.replicas, hmr=args.HMR)
    queue.write_queue_file(kind=args.queue_template)

    
if __name__ == '__main__':
    
    args = parse_input()
    if args.mode == 'full-protocol':
        do_full_openMM_protocol(args)
    elif args.mode == 'from-equilibration':
        do_openMM_from_equil(args)
    elif args.mode == 'Amber-preparation':
        do_AMBER_preparation(args)
