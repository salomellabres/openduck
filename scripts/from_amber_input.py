import argparse
import pickle
import shutil
import sys
from pathlib import Path
from os import chdir
import yaml
import requests
import json
import socket

# from duck.steps.parametrize import prepare_system
from duck.steps.equlibrate import do_equlibrate
from duck.utils.check_system import check_if_equlibrated
from duck.utils import duck_stuff
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md

import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import parmed
from parmed.geometry import distance2
import math
import operator


def check_same(atom, chain, res_name, res_number, atom_name):
    if atom.residue.name == res_name:
        if atom.residue.number == res_number:
            if atom.name == atom_name:
                if atom.residue.chain == chain:
                    return True
    return False


def is_lig(atom):
    # Non-hydrogen
    if (atom.residue.name == "UNL" or atom.residue.name=="LIG") and atom.atomic_number > 1:
        return True


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


def find_interaction_amber_input(combined_pmd, chunk_file, res_atom):

    chunk = parmed.load_file(chunk_file, structure=True)
    if chunk_file.split('.')[-1]== 'mol2':
        rename_mol2_residues(chunk)

    chain = '' # chain information gets lost in the chunk
    res_name = res_atom.split("_")[0][:3]
    res_number = int(res_atom.split("_")[0][3:6])
    atom_name = res_atom.split("_")[-1]

    for atom in chunk.atoms:
        # chain information gets lost in the chunking process, which is expected. Maybe need a check to make sure
        if check_same(atom, chain, res_name, res_number, atom_name):
            print(atom.name, atom.residue.name)
            chunk_atom = atom
            break

    if 'chunk_atom' not in locals():
        print('Cannot find the interaction atom in the chunk. Check residue labelling')
        return

    else:
        # Find the positions of the protein atoms in the combined_system atoms and the chunk

        # chunk_dict = chunk_residue_number: combined_pmd_residue_number
        chunk_dict = {}
        for cd, cr in zip(combined_pmd.residues, chunk.residues):
            chunk_dict[cr.number] = cd.number

        for atom in combined_pmd.atoms:
            if check_same(atom, chain, res_name, chunk_dict[chunk_atom.residue.number], atom_name):
                print(atom.name, atom.residue.name)
                target_protein_atom = atom
                break

        if 'target_protein_atom' not in locals():
            print('Cannot find the interaction atom in the prepared system. Check residue labelling')
            return

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



def equilibrate_from_amber_prep(interaction, prmtop_file, inpcrd_file, chunk,  gpu_id, force_constant_eq=1.0):
    # Load the prepared system:
    c_pm = parmed.load_file(prmtop_file, inpcrd_file)
    # Find the interations
    keyInteraction = find_interaction_amber_input(combined_pmd=c_pm, chunk_file=chunk, res_atom=interaction)
    # Platform definition
    platformProperties = {}
    if gpu_id != None:
        platform = mm.Platform_getPlatformByName("CUDA")
        platformProperties["CudaPrecision"] = "double"
    else:
        platform = mm.Platform_getPlatformByName("CPU")
    platformProperties["DeterministicForces"] = 'true'

    complex = "./complex_system.pickle"
    pickle_out = open(complex, "wb")
    pickle.dump([c_pm], pickle_out)
    pickle_out.close()

    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + keyInteraction
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)

    do_equlibrate(keyInteraction=keyInteraction, force_constant_equilibrate=force_constant_eq, gpu_id=gpu_id)
    if not check_if_equlibrated("density.csv", 1):
        raise EquilibrationError("System is not equilibrated.")


def duck_smd_runs(input_checkpoint, pickle, num_runs, md_len, gpu_id, start_dist, init_velocity, save_dir):
    shutil.copyfile(input_checkpoint, "equil.chk")
    shutil.copyfile(pickle, "complex_system.pickle")

    # Now do the MD
    # remember start_dist
    if not Path(save_dir).exists(): save_dir.mkdir()
    for i in range(num_runs):
        if i == 0:
            md_start = "equil.chk"
        else:
            md_start = str(Path(save_dir, "md_" + str(i - 1) + ".chk"))
        log_file = str(Path(save_dir, "md_" + str(i) + ".csv"))
        perform_md(
            md_start,
            str(Path(save_dir, "md_" + str(i) + ".chk")),
            log_file,
            str(Path(save_dir, "md_" + str(i) + ".pdb")),
            md_len=md_len,
            gpu_id=gpu_id,
        )
        # Open the file and check that the potential is stable and negative
        if not check_if_equlibrated(log_file, 3):
            print("SYSTEM NOT EQUILIBRATED")
            sys.exit()
        # Now find the interaction and save to a file


        run_steered_md(
            300,
            str(Path(save_dir, "md_" + str(i) + ".chk")),
            str(Path(save_dir, "smd_" + str(i) + "_300.csv")),
            str(Path(save_dir, "smd_" + str(i) + "_300.dat")),
            str(Path(save_dir, "smd_" + str(i) + "_300.pdb")),
            str(Path(save_dir, "smd_" + str(i) + "_300.dcd")),
            start_dist,
            init_velocity=init_velocity,
            gpu_id=gpu_id,
        )
        run_steered_md(
            325,
            str(Path(save_dir, "md_" + str(i) + ".chk")),
            str(Path(save_dir, "smd_" + str(i) + "_325.csv")),
            str(Path(save_dir, "smd_" + str(i) + "_325.dat")),
            str(Path(save_dir, "smd_" + str(i) + "_325.pdb")),
            str(Path(save_dir, "smd_" + str(i) + "_325.dcd")),
            start_dist,
            init_velocity=init_velocity,
            gpu_id=gpu_id,
        )

def run_single_direc(direc):
    """

    :param direc: Path to directory holding the data and yaml file
    :return:
    """
    chdir(direc)
    yaml_file = Path('run.yaml')
    out_data = yaml.load(yaml_file.read_text(), Loader=yaml.FullLoader)

    # Set all the params
    prmtop = out_data["prmtop"]
    inpcrd = out_data["inpcrd"]
    chunk = out_data["chunk"]
    protein_code = out_data["prot_code"]
    protein_interaction = out_data["prot_int"]
    cutoff = float(out_data["cutoff"])
    md_len = float(out_data["md_len"])
    start_distance = float(out_data["distance"])
    init_velocity = float(out_data["init_velocity"])
    num_smd_cycles = int(out_data["num_smd_cycles"])
    gpu_id = int(out_data["gpu_id"])

    # Fix the gpu_id
    if socket.gethostname() != 'gandalf':
        gpu_id = 0

    save_dir = Path(Path.cwd(), 'duck_runs')
    if not save_dir.exists():
        save_dir.mkdir()


    equilibrate_from_amber_prep(interaction=protein_interaction,
                                prmtop_file=prmtop,
                                inpcrd_file=inpcrd,
                                chunk=chunk,
                                gpu_id=gpu_id,
                                force_constant_eq=1.0)

    pickle_path = Path('complex_system.pickle')
    new_pickle_path = Path('cs.pickle')
    pickle_path.rename(new_pickle_path)

    equil_path = Path('equil.chk')
    new_equil_path = Path('eql.chk')
    equil_path.rename(new_equil_path)
    print('checkpoint_path', equil_path)

    duck_smd_runs(input_checkpoint=new_equil_path,
                  pickle=new_pickle_path,
                  num_runs=num_smd_cycles,
                  md_len=md_len,
                  gpu_id=gpu_id,
                  start_dist=start_distance,
                  init_velocity=init_velocity,
                  save_dir=save_dir)

def main():
    """
    Runs from an input file containing the paths to directories with DUck input
    :return:
    """

    parser = argparse.ArgumentParser(description='Perform dynamic undocking in the duck_pond')
    parser.add_argument('-i', '--input_file', help='input file containing the paths to directories with DUck input')

    args = parser.parse_args()
    dir_list = Path(args.input_file).read_text().strip().split('\n')

    requests.post('https://notify.run/jJveF2DXBuQgHNJr', data=f"Starting run")
    for i, d in enumerate(dir_list):
        try:
            run_single_direc(Path(d))
            requests.post('https://notify.run/jJveF2DXBuQgHNJr',
                          data=f"Finished run {Path(d).name}")
        except Exception as e:
            with open(Path(d, 'error.log'), 'w') as f:
                f.write(str(e))


if __name__=='__main__':
    from multiprocessing import Pool
    import sys

    input_file = sys.argv[1]
    dir_list = Path(Path.cwd(),input_file).read_text().strip().split('\n')

    with Pool(2) as p:
        p.map(run_single_direc, dir_list, chunksize=1)

