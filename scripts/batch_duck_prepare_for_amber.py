import argparse
from cgitb import small
import pickle
import os
import shutil
import multiprocessing as mp
from contextlib import redirect_stdout,redirect_stderr

from duck.steps.parametrize import prepare_system
from duck.utils.cal_ints import find_interaction
from duck.utils.amber_inputs import Queue_templates, Amber_templates

def ligand_string_generator(file):
    with open(file) as fh:
        mol = []
        for line in fh:
            line=line.rstrip()
            mol.append(line)
            if line == '$$$$':
                new_mol = mol
                mol = []
                yield '\n'.join(new_mol)

def write_string_to_file(self, file,string):
    with open(file, 'w') as fh:
        fh.write(string)

def prepare_sys_for_amber(ligand_file, protein_file, chunk_file, interaction, HMR,  small_molecule_forcefield='SMIRNOFF', water_ff_str = 'tip3p.xml', forcefield_str='amber99sb.xml', ionic_strength = 0.1, box_buffer_distance = 10, waters_to_retain="waters_to_retain.pdb"):
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
    amber = Amber_templates(structure=p[0], interaction=p[1:],hmr=HMR)
    amber.write_all_inputs()

def prepare_ligand_in_folder(ligand_string, lig_indx, protein, chunk, interaction, HMR, base_dir, small_molecule_forcefield = 'SMIRNOFF', water_model = 'tip3p', forcefield = 'amber99sb', ion_strength = 0.1, box_buffer_distance = 10, waters_to_retain='waters_to_retain.pdb'):

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
                                  box_buffer_distance = box_buffer_distance, waters_to_retain="waters_to_retain.pdb")

    #os.chdir(f'..')
    return(f'Lig_target_{lig_indx} prepared correctly')

#global?
result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)
# handle raised errors
def handle_error(error):
	print(error, flush=True)

def main():
    parser = argparse.ArgumentParser(description='Prepare system for dynamic undocking')
    parser.add_argument('-p', '--protein', help='chunk protein in PDB format')
    parser.add_argument('-l', '--ligands', help='Ligands in sdf format')
    parser.add_argument('-i', '--interaction', help='Protein atom to use for ligand interaction.')
    parser.add_argument('-q', '--queue-template', type=str, default = None, help='Write out a queue template from the following: [Slurm | SGE | local]')
    parser.add_argument('-H', '--HMR', action='store_true', help ='Perform Hydrogen Mass Repartition on the topology and use it for the input files')
    parser.add_argument('-r', '--replicas', type=int, default=5, help='Ammount of SMD replicas to perform')
    parser.add_argument('-w', '--wqb_threshold', type=float, default=7.0, help='WQB threshold to stop the simulations')
    parser.add_argument('-n', '--n-threads', type=int, default=None, help='Ammount of CPU to use, default will be all available CPU')
    parser.add_argument('-f', '--small_molecule_forcefield', type=str, default='SMIRNOFF', help='Small Molecules forcefield to employ from the following: [SMIRNOFF | GAFF2 | ESPALOMA]')
    parser.add_argument('-c', '--chunk', default = None, help='Chunked protein')
    parser.add_argument('-s', '--water-model', default='tip3p', type=str.lower, help='Water model to parametrize the solvent with. Chose from the following: [TIP3P | TIP4PEW | SPCE] ')
    parser.add_argument('-pf','--protein-forcefield', default='amber99sb', type=str.lower, help='Protein forcefield to parametrize the chunked protein. Chose form the following: [amber99sb | amber14-all]')
    parser.add_argument('-ion','--ionic-strength', default=0.1, type=float, help='Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M')
    parser.add_argument('-b','--solvent-buffer-distance', default=10, type=float, help='Buffer distance between the periodic box and the protein. Default = 10 A')
    parser.add_argument('-water','--waters-to-retain', default='waters_to_retain.pdb', type=str, help='PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.')

    
    args = parser.parse_args()

    # Initializing pool of cpus
    if not args.n_threads:
        pool = mp.Pool(mp.cpu_count())
        print(f'Number of Threads to use not specified, using {mp.cpu_count()}')
    else:
        pool = mp.Pool(args.n_threads)

    #If not chunk given, use protein as chunk (Only relevant for purposes of identifying the interaction)
    if not args.chunk: args.chunk = args.protein

    #Where am I?
    base_dir = os.getcwd()
    
    # Iterate_ligands
    r = [pool.apply_async(prepare_ligand_in_folder,
                          args=(ligand_string, j+1, args.protein, args.chunk,
                                args.interaction, args.HMR, base_dir,
                                args.small_molecule_forcefield, args.water_model, args.protein_forcefield,
                                args.ionic_strength, args.solvent_buffer_distance, args.waters_to_retain),
                          callback=log_result,
                          error_callback=handle_error) for j, ligand_string in enumerate(ligand_string_generator(args.ligands))]
    pool.close()
    pool.join()

    # write queue array
    if args.queue_template:
        queue = Queue_templates(wqb_threshold=args.wqb_threshold, replicas=args.replicas, array_limit=len(r), hmr=args.HMR)
        queue.write_queue_file(kind=args.queue_template)
    #handle exceptions and results to see if everything went well
    for result in r:
        value = result.get()
        print(value)
    

if __name__ == "__main__":
    main()
