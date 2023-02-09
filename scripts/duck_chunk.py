#!/usr/bin/env python
from duck.steps.chunk import duck_chunk
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Prepare chunk from protein for Dynamic Undocking')
    parser.add_argument('-y', '--input', type=str, default = None, help='Input yaml with the protein interaction (prot_int), distance selection cutoff (cutoff),the receptor pdb (apo_pdf_file) and the ligand mol file (mol_file).')
    parser.add_argument('-l', '--ligand', type=str, default = None, help='ligand mol file to use as reference for interaction.')
    parser.add_argument('-i', '--interaction', type=str, default = None, help='Protein atom to use for ligand interaction.')
    parser.add_argument('-p', '--protein', type=str, default = None, help='Protein pdb file to chunk.')
    parser.add_argument('-c', '--cutoff', type=float, default = None, help='Cutoff distance to define chunk.')
    parser.add_argument('-o', '--output', type=str, default = 'protein_out.pdb', help='Output file for chunked receptor, default: protein_out.pdb. Protonated output will be stored in protonated_{output}')
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
        ignore_buffers=False
        if 'ignore_buffers' in out_data:
            ignore_buffers = out_data['ignore_buffers']
        output = 'protein_out.pdb'
        if 'output' in out_data:
            output = out_data['output']
    else:
        prot_file = args.protein
        cutoff = args.cutoff
        mol_file = args.lig
        prot_int = args.prot_int
        ignore_buffers = args.ignore_buffers
        output = args.output
    duck_chunk(prot_file,mol_file,prot_int,cutoff, output_name=output, ignore_buffers=ignore_buffers)