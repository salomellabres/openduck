"""
A collection of functions for chunking that allow for more control over the starting chunk.
The idea is to input any waters and ions we want to keep/remove aftre visual inspection of the initial chunk.
"""
import os, pkg_resources
import argparse
from rdkit import Chem
import parmed
from parmed.geometry import distance2
import operator
import yaml

def check_same(atom, chain, res_name, res_number, atom_name):
    if atom.residue.name == res_name:
        if atom.residue.number == res_number:
            if atom.name == atom_name:
                if atom.residue.chain == chain:
                    return True
    return False

def get_corresponding_residue(residue, reference_protein):
    prot_atom = residue.atoms[0]
    distance_atom_1 = [(x.idx, distance2(x, prot_atom)) for x in reference_protein.atoms]
    distance_atom_1.sort(key=operator.itemgetter(1))
    return reference_protein.atoms[distance_atom_1[0][0]].residue

def find_residue(prot_file, chain, resname, resnum):
    parmed_prot = parmed.load_file(prot_file)
    for residue in parmed_prot.residues:
        if residue.chain == chain:
            if residue.name == resname:
                if residue.number == int(resnum):
                    return residue

def find_atom(res_atom=None, prot_file=None, combined_pmd=None):
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
    distance_atom_1 = [(x.idx, distance2(x, prot_atom)) for x in combined_pmd.atoms]
    distance_atom_1.sort(key=operator.itemgetter(1))
    return distance_atom_1, prot_atom

def return_tleap(prot_protein_chunk, out_save, disulfides=[]):
    param_f_path = pkg_resources.resource_filename(
        "duck", "parameters/tleap/leaprc.ff14SB.redq"
    )
    bond_list = [
        "bond mol." + str(num[0]) + ".SG mol." + str(num[1]) + ".SG "
        for num in disulfides
    ]
    return (
        """source """
        + param_f_path
        + """
mol = loadpdb """
        + prot_protein_chunk
        + """
"""
        + "\n".join(bond_list)
        + """
savepdb mol """
        + out_save
        + """
quit"""
    )


def do_tleap(prot_protein_chunk, out_save, disulfides=[]):
    # Now do tleap
    out_f = open("run.tleap", "w")
    out_f.write(return_tleap(prot_protein_chunk, out_save, disulfides))
    out_f.close()
    os.system("tleap -f run.tleap")


def add_cap(x, atom_set):
    # Add the cap
    atom_set.add(x.idx)
    # Add the capping atoms
    for y in x.bond_partners:
        atom_set.add(y.idx)
    return atom_set


def find_neighbour_residues(residues):
    single_joins = {}
    double_joins = {}
    atom_set = set()
    for resid in residues:
        residue_set = set()
        double_res_set = set()
        for atom in resid.atoms:
            for bond_1 in atom.bonds:
                if bond_1.measure() > 3.0:
                    print(
                        "Skipping too long bond (" + str(bond_1.measure()) + ") :  ",
                        bond_1.atom1,
                        bond_1.atom2,
                    )
                    continue
                if bond_1.atom1 == atom:
                    x = bond_1.atom2
                else:
                    x = bond_1.atom1
                residue_set.add(x.residue)
                for atom_two in x.residue.atoms:
                    for bond_2 in atom_two.bonds:
                        if bond_2.measure() > 3.0:
                            print(
                                "Skipping too long bond ("
                                + str(bond_2.measure())
                                + ") :  ",
                                bond_2.atom1,
                                bond_2.atom2,
                            )
                            continue
                        if bond_2.atom1 == atom:
                            conn_two = bond_2.atom2
                        elif bond_2.atom1 == atom_two:
                            conn_two = bond_2.atom2
                        else:
                            conn_two = bond_2.atom1
                        double_res_set.add(conn_two.residue)
                atom_set = add_cap(x, atom_set)
        double_res_set.difference_update(residue_set)
        single_joins[resid] = residue_set
        double_joins[resid] = double_res_set
    return single_joins, double_joins, atom_set


def find_neighbours(residues):
    single_joins, double_joins, atom_set = find_neighbour_residues(residues)
    new_residues = set()
    for resid_one in single_joins:
        for resid_two in single_joins:
            if resid_one == resid_two:
                continue
            single_join = single_joins[resid_two].intersection(single_joins[resid_one])
            # If it's just a terminal residue
            if resid_one not in single_joins[resid_two]:
                double_join = double_joins[resid_two].intersection(
                    single_joins[resid_one]
                )
            else:
                double_join = set()
            new_set = single_join.union(double_join)
            for new_res in new_set:
                residues.add(new_res)
                new_residues.add(new_res)
    return new_residues, atom_set


def convert_to_ace_nme(subset):
    remove_res_ids = []
    remove_atom_ids = []
    for residue in subset.residues:
        if len(residue) == 3:
            if set([x.name for x in residue.atoms]) == set(["CA", "C", "O"]):
                residue.name = "ACE"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
            # If it's a proline
            if set([x.name for x in residue.atoms]) == set(["CA", "CD", "N"]):
                residue.name = "NME"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
                    if atom.name == "CD":
                        remove_atom_ids.append(str(atom.idx + 1))
        elif len(residue) == 2:
            if set([x.name for x in residue.atoms]) == set(["CA", "N"]):
                residue.name = "NME"
                for atom in residue.atoms:
                    if atom.name == "CA":
                        atom.name = "CH3"
            elif set([x.name for x in residue.atoms]) == set(["CB", "SG"]):
                remove_res_ids.append(str(residue.idx + 1))
    if remove_atom_ids != []:
        subset = subset["!(@" + ",".join(remove_atom_ids) + ")"]
    if remove_res_ids != []:
        subset = subset["!(:" + ",".join(remove_res_ids) + ")"]
    return subset


def remove_prot_buffers_alt_locs(prot_file):
    output_file = "no_buffer_altlocs.pdb"
    solvents = ["NA", "CL", "SO4", "EDO", "FMT", "P04", "DMS", "EPE"]
    # Remove hydrogens and solvents and buffers
    protein = parmed.load_file(prot_file)["!(:HOH," + ",".join(solvents) + ")"][
        "!(:=@H=)"
    ]
    protein.write_pdb(output_file, altlocs="first")
    return output_file

def remove_prot_altlocs(prot_file):
    output_file = "no_buffer_altlocs.pdb"
    # Remove hydrogens and solvents and buffers
    protein = parmed.load_file(prot_file)[
        "!((:=@H=) & !:HOH)"
    ]
    waters = protein[":HOH"]
    waters.write_pdb("waters_to_retain.pdb")
    protein = protein["!(:=@H=)"]
    protein.write_pdb(output_file, altlocs="first")
    return output_file


def find_disulfides(input_file, threshold=6.2):
    structure = parmed.load_file(input_file)
    sulfurs = [x for x in structure.atoms if x.residue.name == "CYS" and x.name == "SG"]
    #sulfurs = [x for x in structure.atoms if (x.residue.name == "CYS" or x.residue.name == "CYX") and x.name == "SG"]
    disulfides = []
    for atom_one in sulfurs:
        for atom_two in sulfurs:
            if atom_one.idx >= atom_two.idx:
                continue
            dist = distance2(atom_one, atom_two)
            if dist < threshold:
                atom_one.residue.name = "CYX"
                atom_two.residue.name = "CYX"
                disulfides.append((atom_one.residue.number, atom_two.residue.number))
    structure.write_pdb(input_file)
    return disulfides


def find_res_idx(protein, chain, res_name, res_num):
    for residue in protein.residues:
        if residue.chain == chain:
            if residue.name == res_name:
                if residue.number == res_num:
                    return residue.idx + 1



def prot_with_pdb_fixer(chunk_protein, chunk_prot_protein):
    os.system(
        "pdbfixer "
        + chunk_protein
        + " --replace-nonstandard --output="
        + chunk_prot_protein
    )
    return [chunk_prot_protein]


def add_ter_records(input_file, output_file):
    lines = open(input_file).readlines()
    output_f = open(output_file, "w")
    for line in lines:
        output_f.write(line)
        if "NME" in line:
            output_f.write("TER\n")
    return [output_file]

def chunk_with_amber(
    mol_file="MURD-x0349.mol",
    prot_file="MURD-x0349_apo.pdb",
    interaction="A_LYS_311_N",
    out_save="protein_out.pdb",
    cutoff=9.0,
    orig_prot="MURD-x0349_apo.pdb",
    residues_to_add=None,
    residues_to_remove=None
):
    # Load up the topology
    mol = Chem.MolFromMolFile(mol_file)
    pdb_mol_file = mol_file.replace(".mol", ".pdb") # not sure why?
    Chem.MolToPDBFile(mol, pdb_mol_file)
    protein = parmed.load_file(prot_file)
    # get these details
    atom_idx, prot_atom = find_atom(interaction, orig_prot, protein)
    mask = parmed.amber.AmberMask(
        protein, "@" + str(atom_idx[0][0] + 1) + "<:" + str(cutoff)
    )
    residues = set(
        [protein.atoms[i].residue for i, x in enumerate(mask.Selection()) if x == 1]
    )

        
#     residues = set([x for x in residues if x.name != "UNL"])
    # # Find all the residues that are connected to two residues in this list of residues.
    new_residues, atom_set = find_neighbours(residues)

    if residues_to_remove is not None:
        parmed_residues_to_remove = []
        # the residues are input as a list of strings, in the same format as the interaction residues:
        for rcode in residues_to_remove:
            chain, rname, rnum = rcode.split('_')
            print(chain, rname, rnum)
            tar_res = find_residue(prot_file, chain, rname, rnum)
            print(tar_res)
            if tar_res is not None:
                parmed_residues_to_remove.append(get_corresponding_residue(tar_res, protein))
        all_remove, to_remove_set = find_terminal(set(parmed_residues_to_remove))
        all_remove = list((set([protein.atoms[a].residue for a in all_remove]))) + parmed_residues_to_remove
        print(all_remove)
        print(parmed_residues_to_remove)
        residues = set([x for x in list(residues) if x not in all_remove])
        atom_set = set([a for a in atom_set if protein.atoms[a].residue not in all_remove])
#         residues = set([x for x in list(residues) if x not in parmed_residues_to_remove])
#         atom_set = set([a for a in atom_set if protein.atoms[a].residue not in parmed_residues_to_remove])
        
#     if residues_to_add is not None:
#         res_list = list(residues)
#         parmed_residues_to_add = []
#         for rcode in residues_to_add:
#             chain, rname, rnum = rcode.split('_')
#             tar_res = find_residue(protein, chain, rname, rnum)
            
#             if tar_res is not None:
#                 parmed_residues_to_add.append(tar_res)
#         for res in parmed_residues_to_add:
#             if res not in res_list:
#                 res_list.append(res)
#         residues = set(res_list)
            
    # Collect the atoms
    atom_idx = []
    for res in residues:
        atom_idx.extend([x.idx for x in res.atoms if x.altloc in ["", "A"]])
    
    atom_idx.extend(atom_set)
    subset = protein[atom_idx]
    subset.write_pdb(out_save.replace(".pdb", "_no_ace_nme.pdb"))
    subset = convert_to_ace_nme(
        parmed.load_file(out_save.replace(".pdb", "_no_ace_nme.pdb"))
    )
    subset.write_pdb(out_save)
    add_ter_records(out_save, out_save)
    return [out_save]

def duck_chunk(protein, ligand, interaction, chunk_name, cutoff, ignore_buffers=False, to_keep=None, to_remove=None):
    """
    Same as Simon's duck_chunk.py script, but reads in yaml rather than commandline
    :return:
    """
    chunk_protein = chunk_name
    chunk_protein_prot = chunk_name.replace(".pdb", "_prot.pdb")
    # Do the removal of buffer mols and alt locs
    if not ignore_buffers:
        prot_file = remove_prot_buffers_alt_locs(protein)
    else:
        prot_file = remove_prot_altlocs(protein)
    # Do the chunking and the protonation
    chunk_with_amber(
        ligand, prot_file, interaction, chunk_protein, cutoff, protein, residues_to_add=to_keep, residues_to_remove=to_remove
    )
    # Protonate
    disulfides = find_disulfides(chunk_protein)
    do_tleap(chunk_protein, chunk_protein_prot, disulfides)

    return chunk_protein_prot

def main(direc):
    os.chdir(direc)
    yaml_file = Path('run.yaml')
    out_data = yaml.load(yaml_file.read_text(), Loader=yaml.FullLoader)
        # Set all the params
    protein_file = out_data["apo_pdb_file"]
    ligand_file = out_data["mol_file"]
    chunk_file = "protein_out.pdb"
    protein_interaction = out_data["prot_int"]
    cutoff = float(out_data["cutoff"]) # doesn't matter here, as we're supplying the chunk
    
    try:
        buffers = int(out_data["ignore_buffers"])
        if buffers == 0:
            buffers = False
        else:
            buffers = True
    except KeyError:
        buffers = False
    print(buffers)
    try:
        to_keep = out_data['residues_to_keep'].split(' ')
        if len(to_keep) == 0:
            to_keep = None
            
        to_remove = out_data['residues_to_remove'].split(' ')
        if len(to_remove) == 0:
            to_remove = None

    except KeyError:
        to_keep = None
        to_remove = None
    duck_chunk(protein_file, ligand_file, protein_interaction, chunk_file, cutoff, ignore_buffers=buffers, to_keep=to_keep, to_remove=to_remove)
    
#     md_len = float(out_data["md_len"])
#     start_distance = float(out_data["distance"])
#     init_velocity = float(out_data["init_velocity"])
#     num_smd_cycles = int(out_data["num_smd_cycles"])
#     gpu_id = int(out_data["gpu_id"])

if __name__ == "__main__":
    from pathlib import Path
    
    parser = argparse.ArgumentParser(description='Perform dynamic undocking in the duck_pond')
    parser.add_argument('-i', '--input_file', help='input file containing the paths to directories with DUck input')

    args = parser.parse_args()
    dir_list = Path(args.input_file).read_text().strip().split('\n')

#     dir_list = [x.resolve() for x in Path("fragment_waters").glob('*') if x.is_dir()]
    
    for i, d in enumerate(dir_list):
        main(d)
