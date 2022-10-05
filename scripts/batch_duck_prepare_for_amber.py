import argparse
import pickle
import os
import shutil
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


def main():
    parser = argparse.ArgumentParser(description='Prepare system for dynamic undocking')
    parser.add_argument('-p', '--protein', help='chunk protein in PDB format')
    parser.add_argument('-l', '--ligands', help='Ligands in sdf format')
    parser.add_argument('-i', '--interaction', help='Protein atom to use for ligand interaction.')
    parser.add_argument('-q', '--queue-template', type=str, default = None, help='Write out a queue template from the following: [Slurm | SGE]')
    parser.add_argument('--HMR', type=bool, default=True, help ='Perform Hydrogen Mass Repartition on the topology and use it for the input files')
    parser.add_argument('-r', '--replicas', type=int, default=5, help='Ammount of SMD replicas to perform')
    parser.add_argument('-w', '--wqb_threshold', type=float, default=7.0, help='WQB threshold to stop the simulations')
    args = parser.parse_args()

    # Iterate_ligands
    for j, ligand_string in enumerate(ligand_string_generator(args.ligands)):
        i=j+1
        #Create the ligand folder
        if os.path.isdir(f'LIG_target_{i}'):
            print(f'WARNING: LIG_target_{i} already exist and it will be overwritten.')
            shutil.rmtree(f'./LIG_target_{i}', ignore_errors=True)
        os.mkdir(f'LIG_target_{i}')
        os.chdir(f'LIG_target_{i}')

        # Copying files to ligand foldef; ligand and prot
        write_string_to_file(string=ligand_string, file=f'lig_{i}.mol')
        shutil.copyfile(f'../{args.protein}', f'./{args.protein}', follow_symlinks=True)
        if os.path.isfile('../waters_to_retain.pdb'):
            shutil.copyfile(f'../waters_to_retain.pdb', f'./waters_to_retain.pdb', follow_symlinks=True)

        # Parameterize the ligand    
        prepare_system(f'lig_{i}.mol', args.protein, forcefield_str="amber99sb.xml", hmr=args.HMR)
        
        # Now find the interaction and save to a file
        results = find_interaction(args.interaction, args.protein)
        print(results) # what happens to these?
        with open('complex_system.pickle', 'rb') as f:
            p = pickle.load(f) + results
        with open('complex_system.pickle', 'wb') as f:
            pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)
        #p = (parmed_structure, prot_index, ligand_index, pairmeandistance)
        p[0].save('system_complex.inpcrd', overwrite=True)

        
        #do_equlibrate(force_constant_equilibrate=args.force_constant_eq, gpu_id=args.gpu_id, keyInteraction=p[1:])
        
        write_all_inputs(p[0], p[1:], hmr = args.HMR)
        write_getWqbValues()
        os.chdir(f'..')

    if args.queue_template:
        write_queue_template(args.queue_template, hmr = args.HMR, replicas=args.replicas, wqb_threshold=args.wqb_threshold, array_limit=i)
    

if __name__ == "__main__":
    main()
