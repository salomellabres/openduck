from duck.steps.chunk import chunk_with_amber,do_tleap,remove_prot_buffers_alt_locs,find_disulfides
import argparse

def main(prot_file, mol_file, interaction, cutoff, output_name = 'protein_out.pdb', ignore_buffers=False):
    # A couple of file name
    orig_file = prot_file

    
    chunk_protein_prot = f'protonated_{output_name}'
    # Do the removal of buffer mols and alt locs
    if not ignore_buffers:
        prot_file = remove_prot_buffers_alt_locs(prot_file)
    # Do the chunking and the protonation
    # Chunk
    chunk_with_amber(mol_file,prot_file,interaction,output_name,cutoff,orig_file)
    # Protontate
    disulfides = find_disulfides(output_name)
    do_tleap(output_name, chunk_protein_prot, disulfides)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare chunk from protein for Dynamic Undocking')
    parser.add_argument('-y', '--input', type=str, default = None, help='Input yaml with the protein interaction (prot_int), distance selection cutoff (cutoff),the receptor pdb (apo_pdf_file) and the ligand mol file (mol_file).')
    parser.add_argument('-l', '--ligand', type=str, default = None, help='ligand mol file to use as reference for interaction.')
    parser.add_argument('-i', '--interaction', type=str, default = None, help='Protein atom to use for ligand interaction.')
    parser.add_argument('-p', '--protein', type=str, default = None, help='Protein pdb file to chunk.')
    parser.add_argument('-c', '--cutoff', type=float, default = None, help='Cutoff distance to define chunk.')
    parser.add_argument('-o', '--output', type=float, default = 'protein_out.pdb', help='Output file for chunked receptor, default: protein_out.pdb. Protonated output will be stored in protonated_{output}')
    parser.add_argument('-b', '--ignore-buffers', action='store_true', help='Do not remove buffers (solvent, ions etc.)')

    args = parser.parse_args()
    if args.input is None and (args.lig is None or args.interaction is None or args.protein is None or args.cutoff is None):
        parser.error('The input needs to be either the input yaml with the specified items, or the protein, ligand, interaction and cutoff arguments.')
    elif args.input:
        import yaml
        out_data = yaml.load(open(args.input))
        prot_int = out_data["prot_int"]
        cutoff = out_data["cutoff"]
        prot_file = out_data["apo_pdb_file"]
        mol_file = out_data["mol_file"]
    else:
        prot_file = args.protein
        cutoff = args.cutoff
        mol_file = args.lig
        prot_int = args.prot_int
    main(prot_file,mol_file,prot_int,cutoff, output_name=args.output, ignore_buffers=args.ignore_buffers)