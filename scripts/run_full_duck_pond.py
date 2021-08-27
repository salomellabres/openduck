import argparse
import pickle
import shutil
import sys
from pathlib import Path
from os import chdir
import yaml
import requests
import operator
import math

from parmed.geometry import distance2
from duck.steps.chunk import (
    chunk_with_amber,
    do_tleap,
    remove_prot_buffers_alt_locs,
    find_disulfides,
)
from duck.steps.parametrize import prepare_system
from duck.utils.cal_ints import find_interaction, find_atom, is_lig
from duck.steps.equlibrate import do_equlibrate
from duck.utils.check_system import check_if_equlibrated
from duck.steps.normal_md import perform_md
from duck.steps.steered_md import run_steered_md

def find_ligand_interaction(res_atom=None, protein_file=None, ligand_coords=None, ligand_atomnum=None):
    output_file = "indice.text"
    if not res_atom or prot_file:
        if os.path.isfile(output_file):
            return json.load(open(output_file))
    # Read files
    print("loading pickle")
    pickle_in = open("complex_system.pickle", "rb")
    combined_pmd = pickle.load(pickle_in)[0]
    pickle_in.close()
    distance_atom_1, prot_atom = find_atom(res_atom, prot_file, combined_pmd)
    distance_atom_2 = [
        (x.idx, distance2(x, ligand_coords)) for x in combined_pmd.atoms if (is_lig(x) and x.atomic_number == ligand_atomnum)
    ]
    distance_atom_2.sort(key=operator.itemgetter(1))
    # These are the interactions to find
    index_one = distance_atom_1[0][0]
    # The ligand one
    index_two = distance_atom_2[0][0]
    out_res = [index_one, index_two, math.sqrt(distance_atom_2[0][1])]
    return index_one, index_two, out_res, distance_atom_2[0][1]

def duck_chunk(protein, ligand, interaction, cutoff, ignore_buffers=False):
    """
    Same as Simon's duck_chunk.py script, but reads in yaml rather than commandline
    :return:
    """
    # A couple of file name
    orig_file = prot_file = protein
    chunk_protein = "protein_out.pdb"
    chunk_protein_prot = "protein_out_prot.pdb"
    # Do the removal of buffer mols and alt locs
    if not ignore_buffers:
        prot_file = remove_prot_buffers_alt_locs(prot_file)
    # Do the chunking and the protonation
    chunk_with_amber(
        ligand, prot_file, interaction, chunk_protein, cutoff, orig_file
    )
    # Protonate
    disulfides = find_disulfides(chunk_protein)
    do_tleap(chunk_protein, chunk_protein_prot, disulfides)

    return chunk_protein_prot


def prepare_sys(protein, ligand, interaction, chunk, gpu_id, force_constant_eq=1.0, lig_coords=None, lig_atomnum=None):
    """
    Same as Simon's duck_prepare_sys.py script
    :return:
    """
    prepare_system(ligand, chunk)
    # Now find the interaction and save to a file
    
    if (lig_coords is not None) and (lig_atomnum is not None):
        results = find_ligand_interaction()
    else:
        results = find_interaction(interaction, protein)
    print(results)  # what happens to these?
    
    with open('complex_system.pickle', 'rb') as f:
        p = pickle.load(f) + results
    with open('complex_system.pickle', 'wb') as f:
        pickle.dump(p, f, protocol=pickle.HIGHEST_PROTOCOL)

    #     pickle.dump(l, 'complex_system.pickle')
    # Now do the equlibration
    do_equlibrate(force_constant_equilibrate=force_constant_eq, gpu_id=gpu_id)
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
            checkpoint_in_file=md_start,
            checkpoint_out_file=str(Path(save_dir, "md_" + str(i) + ".chk")),
            csv_out_file=log_file,
            pdb_out_file=str(Path(save_dir, "md_" + str(i) + ".pdb")),
            dcd_out_file=str(Path(save_dir, "md_" + str(i) + ".dcd")),
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
    protein_file = out_data["apo_pdb_file"]
    ligand_file = out_data["mol_file"]
    protein_code = out_data["prot_code"]
    protein_interaction = out_data["prot_int"]
    cutoff = float(out_data["cutoff"])
    md_len = float(out_data["md_len"])
    start_distance = float(out_data["distance"])
    init_velocity = float(out_data["init_velocity"])
    num_smd_cycles = int(out_data["num_smd_cycles"])
    gpu_id = int(out_data["gpu_id"])
    
    try:
        lig_coords = [float(x) for x in out_data['lig_coords'].split(' ')]
        lig_atomnum = int(out_data['lig_atomnum'])
        
    except IndexError:
        lig_coords = None
        lig_atomnum = None

    save_dir = Path(direc, 'duck_runs')
    if not save_dir.exists(): save_dir.mkdir()

    chunk_prot_fname = duck_chunk(protein=protein_file,
                                 ligand=ligand_file,
                                 interaction=protein_interaction,
                                 cutoff=cutoff)


    prepare_sys(protein=protein_file,
                ligand=ligand_file,
                interaction=protein_interaction,
                chunk=str(Path(direc,chunk_prot_fname)),
                gpu_id=gpu_id)

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
        #try:
        run_single_direc(Path(d))
        #requests.post('https://notify.run/jJveF2DXBuQgHNJr', data=f"Ran {i} out of {len(d)} runs, name: {str(d.name)}")
        #except Exception as e:
            #with open(Path(d, 'error.log'), 'w') as f:
                #f.write(str(e))

if __name__=='__main__':
    main()
