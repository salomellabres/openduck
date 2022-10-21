import argparse
import pickle

try:
    from duck.steps.parametrize import prepare_system
    from duck.utils.cal_ints import find_interaction
    from duck.utils.amber_inputs import write_all_inputs, write_queue_template
except ModuleNotFoundError:
    print('Dependencies missing; check openmm, pdbfixer, and yank are installed from Omnia.')


def main():
    parser = argparse.ArgumentParser(description='Prepare system for dynamic undocking')
    parser.add_argument('-p', '--protein', help='chunk protein in PDB format')
    parser.add_argument('-l', '--ligand', help='Ligands in sdf format')
    parser.add_argument('-i', '--interaction', help='Protein atom to use for ligand interaction.')
    parser.add_argument('-q', '--queue-template', type=str, default = None, help='Write out a queue template from the following: [Slurm | SGE]')
    parser.add_argument('-H','--HMR', action='store_true', help ='Perform Hydrogen Mass Repartition on the topology and use it for the input files')
    parser.add_argument('-r', '--replicas', type=int, default=5, help='Ammount of SMD replicas to perform')
    parser.add_argument('-w', '--wqb_threshold', type=float, default=7.0, help='WQB threshold to stop the simulations')
    parser.add_argument('-f', '--small_molecule_forcefield', type=str, default='SMIRNOFF', help='Small Molecules forcefield to employ from the following: [SMIRNOFF | GAFF2 | ESPALOMA]')
    parser.add_argument('-c', '--chunk', default = None, help='Chunked protein')

    args = parser.parse_args()
    
    #If not chunk given, use protein as chunk (Only relevant for purposes of identifying the interaction)
    if not args.chunk: args.chunk = args.protein

    # Parameterize the ligand
    prepare_system(args.ligand, args.chunk, forcefield_str="amber99sb.xml", hmr=args.HMR, small_molecule_ff=args.small_molecule_forcefield)
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

    if args.queue_template:
        write_queue_template(args.queue_template, hmr = args.HMR, replicas=args.replicas, wqb_threshold=args.wqb_threshold)
    


if __name__ == "__main__":
    main()
