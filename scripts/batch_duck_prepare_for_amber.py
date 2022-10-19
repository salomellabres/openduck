import argparse
from cgitb import small
import pickle
import os
import shutil
import multiprocessing as mp
from contextlib import redirect_stdout,redirect_stderr
try:
    from duck.steps.parametrize import prepare_system
    from duck.utils.cal_ints import find_interaction
    from duck.utils.amber_inputs import write_all_inputs, write_queue_template, write_string_to_file, write_getWqbValues
except ModuleNotFoundError:
    print('Dependencies missing; check openmm, pdbfixer, and yank are installed from Omnia.')

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

def prepare_sys_for_amber(ligand_file, protein_file, interaction, HMR,  small_molecule_forcefield='SMIRNOFF'):
    # Parameterize the ligand
    prepare_system(ligand_file, protein_file, forcefield_str="amber99sb.xml", hmr=HMR, small_molecule_ff=small_molecule_forcefield)
    
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
    
    write_all_inputs(p[0], p[1:], hmr = HMR)
    write_getWqbValues()

def prepare_ligand_in_folder(ligand_string, lig_indx, protein, interaction, HMR, base_dir, small_molecule_forcefield = 'SMIRNOFF'):

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
            if os.path.isfile('../waters_to_retain.pdb'):
                shutil.copyfile(f'../waters_to_retain.pdb', f'./waters_to_retain.pdb', follow_symlinks=True)

            prepare_sys_for_amber(f'lig_{lig_indx}.mol', protein, interaction, HMR, small_molecule_forcefield=small_molecule_forcefield)

    #os.chdir(f'..')
    return(f'Lig_target_{lig_indx} prepared correctly')

#global?
result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

def main():
    parser = argparse.ArgumentParser(description='Prepare system for dynamic undocking')
    parser.add_argument('-p', '--protein', help='chunk protein in PDB format')
    parser.add_argument('-l', '--ligands', help='Ligands in sdf format')
    parser.add_argument('-i', '--interaction', help='Protein atom to use for ligand interaction.')
    parser.add_argument('-q', '--queue-template', type=str, default = None, help='Write out a queue template from the following: [Slurm | SGE]')
    parser.add_argument('-H','--HMR', action='store_true', help ='Perform Hydrogen Mass Repartition on the topology and use it for the input files')
    parser.add_argument('-r', '--replicas', type=int, default=5, help='Ammount of SMD replicas to perform')
    parser.add_argument('-w', '--wqb_threshold', type=float, default=7.0, help='WQB threshold to stop the simulations')
    parser.add_argument('-n', '--n-threads', type=int, default=None, help='Ammount of CPU to use, default will be all available CPU')
    parser.add_argument('-f', '--small_molecule_forcefield', type=str, default='SMIRNOFF', help='Ŝmall Molecules forcefield to employ from the following: [SMIRNOFF | GAFF2 | ESPALOMA]')
    args = parser.parse_args()

    if not args.n_threads:
        pool = mp.Pool(mp.cpu_count())
        print(f'Number of Threads to use not specified, using {mp.cpu_count()}')
    else:
        pool = mp.Pool(args.n_threads)
    base_dir = os.getcwd()
    # Iterate_ligands
    r = [pool.apply_async(prepare_ligand_in_folder, args=(ligand_string, j+1, args.protein, args.interaction, args.HMR, base_dir, args.small_molecule_forcefield), callback=log_result) for j, ligand_string in enumerate(ligand_string_generator(args.ligands))]
    pool.close()
    pool.join()

    #handle exceptions and results to see if everything went well
    for result in r:
        value = result.get()
        print(value)
    if args.queue_template:
        write_queue_template(args.queue_template, hmr = args.HMR, replicas=args.replicas, wqb_threshold=args.wqb_threshold, array_limit=len(r))
    

if __name__ == "__main__":
    main()
