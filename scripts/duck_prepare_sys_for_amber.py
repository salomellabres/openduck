#!/usr/bin/env python

import argparse
import pickle
from duck.steps.parametrize import prepare_system
from duck.utils.cal_ints import find_interaction
from duck.utils.amber_inputs import Queue_templates, Amber_templates


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
    parser.add_argument('-s', '--water-model', default='tip3p', type=str.lower, help='Water model to parametrize the solvent with. Chose from the following: [TIP3P | TIP4PEW | SPCE] ')
    parser.add_argument('-pf','--protein-forcefield', default='amber99sb', type=str.lower, help='Protein forcefield to parametrize the chunked protein. Chose form the following: [amber99sb | amber14-all]')
    parser.add_argument('-ion','--ionic-strength', default=0.1, type=float, help='Ionic strength (concentration) of the counter ion salts (Na+/Cl+). Default = 0.1 M')
    parser.add_argument('-b','--solvent-buffer-distance', default=10, type=float, help='Buffer distance between the periodic box and the protein. Default = 10 A')
    parser.add_argument('-water','--waters-to-retain', default='waters_to_retain.pdb', type=str, help='PDB File with structural waters to retain water moleules. Default is waters_to_retain.pdb.')

    args = parser.parse_args()
    
    #If not chunk given, use protein as chunk (Only relevant for purposes of identifying the interaction)
    if not args.chunk: args.chunk = args.protein

    # Parameterize the ligand
    prepare_system(args.ligand, args.chunk, forcefield_str=f'{args.protein_forcefield}.xml', water_ff_str = f'{args.water_model}',
                   hmr=args.HMR, small_molecule_ff=args.small_molecule_forcefield, waters_to_retain=args.waters_to_retain,
                   box_buffer_distance = args.solvent_buffer_distance, ionicStrength = args.ionic_strength)
    # Now find the interaction and save to a file
    results = find_interaction(args.interaction, args.protein)
    print(results) # what happens to these?
    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + results
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)
    p[0].save('system_complex.inpcrd', overwrite=True)
        
    amber = Amber_templates(structure=p[0], interaction=p[1:],hmr=args.HMR)
    amber.write_all_inputs()

    if args.queue_template:
        queue = Queue_templates(wqb_threshold=args.wqb_threshold, replicas=args.replicas, hmr=args.HMR)
        queue.write_queue_file(kind=args.queue_template)
    


if __name__ == "__main__":
    main()
