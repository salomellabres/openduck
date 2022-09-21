import argparse
import pickle

try:
    from duck.steps.parametrize import prepare_system
    from duck.utils.cal_ints import find_interaction
    from duck.utils.amber_inputs import write_all_inputs, write_queue_template
    from duck.steps.equlibrate import do_equlibrate
except ModuleNotFoundError:
    print('Dependencies missing; check openmm, pdbfixer, and yank are installed from Omnia.')


def main():
    parser = argparse.ArgumentParser(description='Prepare system for dynamic undocking')
    parser.add_argument('-p', '--protein', help='Apoprotein in PDB format')
    parser.add_argument('-l', '--ligand', help='Ligand in mol format')
    # parser.add_argument('-o', '--output', help="PDB output")
    parser.add_argument('-c', '--chunk', help='Chunked protein')
    parser.add_argument('-i', '--interaction', help='Protein atom to use for ligand interaction.')
    parser.add_argument('-s', '--seed', type=int, help='Random seed.')
    parser.add_argument('--gpu-id', type=int, help='GPU ID (optional); if not specified, runs on CPU only.')
    parser.add_argument('--force-constant-eq', type=float, default=1.0, help='Force constant for equilibration.')
    parser.add_argument('--queue-template', type=str, default = None, help='Write out a queue template from the following: Slurm, CTEP')

    args = parser.parse_args()
    # Parameterize the ligand
    prepare_system(args.ligand, args.chunk, forcefield_str="amber99sb.xml")
    # Now find the interaction and save to a file
    results = find_interaction(args.interaction, args.protein)
    print(results) # what happens to these?
    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + results
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)
    #p = (parmed_structure, prot_index, ligand_index, pairmeandistance)
    p[0].save('system_complex.inpcrd', overwrite=True)
    
    do_equlibrate(force_constant_equilibrate=args.force_constant_eq, gpu_id=args.gpu_id, keyInteraction=p[1:])
    
    write_all_inputs(p[0], p[1:])

    if args.queue_template:
        write_queue_template(args.queue_template)
    


if __name__ == "__main__":
    main()
