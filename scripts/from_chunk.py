import argparse
from pathlib import Path
from os import chdir
import yaml
import requests

from run_full_duck_pond import prepare_sys, duck_smd_runs

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
    chunk_file = out_data["chunk"]
    protein_code = out_data["prot_code"] # This doesn't seem to get used anywhere generally
    protein_interaction = out_data["prot_int"]
    cutoff = float(out_data["cutoff"]) # doesn't matter here, as we're supplying the chunk
    md_len = float(out_data["md_len"])
    start_distance = float(out_data["distance"])
    init_velocity = float(out_data["init_velocity"])
    num_smd_cycles = int(out_data["num_smd_cycles"])
    gpu_id = int(out_data["gpu_id"])

    try:
        buffers = int(out_data["ignore_buffers"])
        if buffers == 0:
            buffers = False
        else:
            buffers = True
    except KeyError:
        buffers = False

    try:
        lig_coords = [float(x) for x in out_data['lig_coords'].split(' ')]
        lig_atomnum = int(out_data['lig_atomnum'])
    except Exception as e:
        print(e)
        lig_coords = None
        lig_atomnum = None

    save_dir = Path(direc, 'duck_runs')
    if not save_dir.exists(): save_dir.mkdir()

    prepare_sys(protein=protein_file,
                ligand=ligand_file,
                interaction=protein_interaction,
                chunk=chunk_file,
                gpu_id=gpu_id,
		lig_coords=lig_coords,
		lig_atomnum=lig_atomnum)
    save_dir = Path(direc, 'duck_runs')

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
            requests.post('https://notify.run/jJveF2DXBuQgHNJr', data=f"Ran {i} out of {len(dir_list)} runs, name: {Path(d).name}")
        except Exception as e:
            with open(Path(d, 'error.log'), 'w') as f:
                f.write(str(e))

if __name__=='__main__':
    main()
