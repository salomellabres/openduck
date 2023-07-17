import json, pickle, sys, os
from parmed.geometry import distance2
from parmed.topologyobjects import Atom
import operator
import parmed
import math

def check_same(atom, chain, res_name, res_number, atom_name):
    '''
    Check if the atom object corresponds to the chain, residue_name, residue number and atom name
    '''
    if atom.residue.name == res_name:
        if atom.residue.number == res_number:
            if atom.name == atom_name:
                if atom.residue.chain == chain:
                    return True
    return False

def is_lig(atom):
    # Non-hydrogen
    if atom.residue.name == "UNL" and atom.atomic_number > 1:
        return True
def is_good_atomtype(atom, elements):
    if type(elements) != list:
        elements = list(elements)
    for e in elements:
        if type(e) != int:
            raise ValueError(f'Element {e} specified in the accepted elements is not an integrer.')
    # Only oxygen and nitrogen
    if atom.atomic_number in elements:
        return True
def find_atom(res_atom=None, prot_file=None, combined_pmd=None):
    '''
    Parse the combined_system pickle and protein file to find the atom id based on the interaction string in this format 'chain_resname_resid_atomname' (i.e. "A_LYS_311_N")
    '''
    # Parse the input data like this -> "A_LYS_311_N"
    chain = res_atom.split("_")[0]
    res_name = res_atom.split("_")[1]
    res_number = int(res_atom.split("_")[2])
    atom_name = res_atom.split("_")[3]
    # Read the original PDB File and find the atom coords
    protein = parmed.load_file(prot_file, structure=True)
    for atom in protein.atoms:
        if check_same(atom, chain, res_name, res_number, atom_name):
            prot_atom = atom
            break
    else:
        raise ValueError("The specified interaction atom was not found in the protein file.")

    distance_atom_1 = [(x.idx, distance2(x, prot_atom)) for x in combined_pmd.atoms]
    distance_atom_1.sort(key=operator.itemgetter(1))
    return distance_atom_1, prot_atom

def find_result(res_atom=None, prot_file=None, combined_pmd=None, accepted_lig_elements=[7,8]):
    # Find the
    distance_atom_1, prot_atom = find_atom(res_atom, prot_file, combined_pmd)
    # Now find the one nearest
    distance_atom_2 = [
        (x.idx, distance2(x, prot_atom)) for x in combined_pmd.atoms if (is_lig(x) and is_good_atomtype(x,accepted_lig_elements ))
    ]
    distance_atom_2.sort(key=operator.itemgetter(1))
    # These are the interactions to find
    index_one = distance_atom_1[0][0]
    # The ligand one
    index_two = distance_atom_2[0][0]
    out_res = [index_one, index_two, math.sqrt(distance_atom_2[0][1])]
    return index_one, index_two, out_res, distance_atom_2[0][1]

def find_interaction(res_atom=None, prot_file=None, lig_HB_elements=[7,8]):
    '''
    Find interaction atom based on the protein_file and the interaction string in this format 'chain_resname_resid_atomname' (i.e. "A_LYS_311_N")
    '''
    output_file = "indice_prueba.text"
    # Read files
    print("loading pickle")
    pickle_in = open("complex_system.pickle", "rb")
    combined_pmd = pickle.load(pickle_in)[0]
    pickle_in.close()
    index_one, index_two, out_res, dist = find_result(res_atom, prot_file, combined_pmd, accepted_lig_elements=lig_HB_elements)
    out_f = open(output_file, "w")
    out_f.write(json.dumps(out_res))
    out_f.close()
    return [index_one, index_two, math.sqrt(dist)]

def find_interaction_amber_input(combined_pmd, chunk_file, res_atom):
    '''
    Find interaction atom based on the protein_file and amber topology and the interaction string in this format 'chain_resname_resid_atomname' (i.e. "A_LYS_311_N")
    '''
    chunk = parmed.load_file(chunk_file, structure=True)
    if chunk_file.split('.')[-1]== 'mol2':
        rename_mol2_residues(chunk)

    chain = '' # chain information gets lost in the chunk
    res_name = res_atom.split("_")[1]
    res_number = int(res_atom.split("_")[2])
    atom_name = res_atom.split("_")[-1]

    for atom in chunk.atoms:
        # chain information gets lost in the chunking process, which is expected. Maybe need a check to make sure
        if check_same(atom, chain, res_name, res_number, atom_name):
            #print(atom.name, atom.residue.name)
            chunk_atom = atom
            break
    if 'chunk_atom' not in locals():
        raise ValueError('Cannot find the interaction atom in the chunk. Check residue labelling')
    else:
        # Find the positions of the protein atoms in the combined_system atoms and the chunk
        # chunk_dict = chunk_residue_number: combined_pmd_residue_number
        chunk_dict = {}
        for cd, cr in zip(combined_pmd.residues, chunk.residues):
            chunk_dict[cr.number] = cd.number
        for atom in combined_pmd.atoms:
            if check_same(atom, chain, res_name, chunk_dict[chunk_atom.residue.number], atom_name):
                #print(atom.name, atom.residue.name)
                target_protein_atom = atom
                break
        if 'target_protein_atom' not in locals():
            raise ValueError('Cannot find the interaction atom in the prepared system. Check residue labelling')
        
        index_one = target_protein_atom.idx
        distance_atom_2 = [(x.idx, distance2(x, target_protein_atom)) for x in combined_pmd.atoms if is_lig(x)]
        distance_atom_2.sort(key=operator.itemgetter(1))
        index_two = distance_atom_2[0][0]
        dist = distance_atom_2[0][1]

        output_file='indice.txt'
        out_f = open(output_file, "w")
        out_f.write(json.dumps([index_one, index_two, math.sqrt(dist)]))
        out_f.close()

        return [index_one, index_two, math.sqrt(dist)]

def clean_up_files(path='./', files_to_delete=[]):
    for file in os.listdir(path):
        # how to check what we want to keep?
        if file in files_to_delete:
            os.remove(path+file)

def rename_mol2_residues(parmed_protein):
    """
    Crashes with altlocs currently...
    :param parmed_protein:
    :return:
    """
    for r in parmed_protein.residues:
        res_name = r.name[:3]
        true_res_number = int(r.name[3:])
        r.name = res_name
        r.number += true_res_number - r.number

if __name__ == "__main__":
    # Define the input
    res_atom = sys.argv[1]
    prot_file = sys.argv[2]
    protid, ligid, dist = find_interaction(res_atom, prot_file)
    print('receptor ID %s\nligand ID %s\ndistance %s'%(protid+1, ligid+1, dist))
