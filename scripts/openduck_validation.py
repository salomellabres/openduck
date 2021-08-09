"""
Script to compare the results of Openduck runs to what has been previously published, as well as runs from Maciej
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from duck_analysis import get_Wqb_value

def make_duck_dict(to_sim):

    #dir_paths = [list(Path(p).glob('duck_runs/smd_*_300.dat')) for p in to_sim]
    dir_paths = [list(Path(p).glob('duck_runs/smd_*_300.dat')) for p in to_sim]

    duck_dict = {}
    for i, dp in enumerate(dir_paths):
        all_data = [get_Wqb_value(p) for p in dp]

        for it, d in enumerate(all_data):
            duck_dict[Path(to_sim[i]).name] = round(np.min([ad[0] for ad in all_data]), 2)
    return duck_dict

def make_wqb_df(data_direcs, sample_size=20, num_samples=10):
    """

    :param data_direcs: list of paths to folders with duck results
    :param sample_size:
    :param num_samples:
    :return:
    """
    duck_dict = make_duck_dict(data_direcs)
    print(duck_dict)
    dir_paths = [list(Path(p).glob('duck_runs/smd_*_300.dat')) for p in data_direcs]

    np.random.seed(3)
    min_mean_dict = {}
    stdev_dict = {}

    for i, dp in enumerate(dir_paths):
        # check that we only include folders with actual duck results calculated
        if len(dp) > 0:
            all_data = np.array([get_Wqb_value(p)[0] for p in dp])
            p_name = dp[0].parents[1].stem
            print(p_name, all_data.shape)
            wqb_mins = []

            for q in range(num_samples):
                sub = np.random.choice(all_data, sample_size, replace=False)
                wqb_mins.append(np.min(sub))
            mean_min = np.mean(wqb_mins)
            min_mean_dict[p_name] = mean_min
            stdev_dict[p_name] = np.array(wqb_mins).std()

    pdb_list = []
    interactions = []
    openduck_2021_list = []
    mean_min_list = []
    stdev_list = []

    for key, wqb in duck_dict.items():
        try:
            mean_min = min_mean_dict[key]
            mean_min_list.append(mean_min)
            stdev = stdev_dict[key]
            stdev_list.append(stdev)
            openduck_2021_list.append(wqb)
            pdb_id = key.split('_')[0]
            pdb_list.append(pdb_id)
            interaction = '_'.join(key.split('_')[1:])
            interactions.append(interaction)
        except KeyError as e:
            print(e)

    wqb_df = pd.DataFrame()
    wqb_df['PDBID'] = pdb_list
    wqb_df['Interaction'] = interactions
    wqb_df['openduck_2021'] = openduck_2021_list
    wqb_df['openduck_mean_min'] = np.round(mean_min_list, 2)
    wqb_df['mean_min_stdev'] = np.round(stdev_list, 2)

    return wqb_df

def compare_to_paper(validation_csv_path, wqb_df):
    paper_df = pd.read_csv(validation_csv_path)
    paper_wqb_list = []
    paper_mean_wqb_list = []
    paper_wqb_stdev = []
    paper_target = []

    for idx, row in wqb_df.iterrows():
        pdb_id = row['PDBID'].split('-')[0]
        if pdb_id == '1fhd-2':
            pdb_id = '1fhd'
        inter_res = row['Interaction'].split('_')[0]
        inter_atom = row['Interaction'].split('_')[1]
        paper_sub_df = paper_df[paper_df['PDB'] == pdb_id]
        paper_row = paper_sub_df[(paper_sub_df['Residue'] == inter_res) & (paper_sub_df['Prot_atom'] == inter_atom)]
        if len(paper_row)>1:
            print(paper_row)
            paper_row=pd.DataFrame(paper_row.iloc[0], columns=paper_sub_df.columns)
        print(pdb_id, inter_res)
        try:
            paper_wqb_list.append(paper_row["WQB [kcal/mol]"].values[0])
            paper_mean_wqb_list.append(paper_row["mean_wqb [kcal/mol]"].values[0])
            paper_wqb_stdev.append(paper_row["std [kcal/mol]"].values[0])
            paper_target.append(paper_row["target_name"].values[0])
        except IndexError:
            paper_wqb_list.append(paper_row["WQB [kcal/mol]"])
            paper_mean_wqb_list.append(paper_row["mean_wqb [kcal/mol]"])
            paper_wqb_stdev.append(paper_row["std [kcal/mol]"])
            paper_target.append(paper_row["target_name"])

    wqb_df['reference_Wqb'] = paper_wqb_list
    wqb_df['reference_mean_wqb'] = paper_mean_wqb_list
    wqb_df['reference_std'] = paper_wqb_stdev
    wqb_df['target_name'] = paper_target
    return wqb_df

def find_iridium_target(validation_csv_path):
    import requests
    pdb_url = "https://data.rcsb.org/rest/v1/core/uniprot/"

    target_names = []
    target_uniprot = []

    paper_df = pd.read_csv(str(validation_csv_path))

    for i, r, in paper_df.iterrows():
        pdb_code = r['PDB'].split('_')[0]
        if '1hgh' in pdb_code:
            pdb_code = '1hgh'
        # Get the data for the first entity
        resp = requests.get(pdb_url + pdb_code+ "/1/")
        if len(resp.json()) == 1:
            target_names.append(str(resp.json()[0]['rcsb_uniprot_entry_name']))
            target_uniprot.append(resp.json()[0]['rcsb_id'])
        else:
            print(f'something went wrong with {pdb_code}')
            print(resp.json())
            target_names.append(" ")
            target_uniprot.append(" ")

    paper_df['target_name'] = target_names
    paper_df['target_uniprot'] = target_uniprot
    paper_df.to_csv(validation_csv_path)
    return paper_df

if __name__ == "__main__":
    # res_dir = Path('/home/jin76872/Desktop/Mih/Data/duck_stuff/openduck_validation_march2021/openduck_validation_march2021')
    # d_direcs = Path(res_dir, 'duck_pond_inputs_all.txt').read_text().strip().split('\n')
    # d_direcs = [Path(res_dir, Path(dd).name) for dd in d_direcs]

   
    res_dir = Path(
        '/home/jin76872/Desktop/Mih/Data/duck_stuff/openduck_iridium_april_2021/iridium_ducks')
    all_runs = list(Path('iridium_ducks').glob('iridium_dat_files/*.dat'))
    to_sim = list(set([(x).name.split("_smd")[0] for x in all_runs]))

    d_direcs = [list(res_dir.glob(f'iridium_dat_files/{p}_smd_*.dat')) for p in to_sim]
    #d_direcs = list(res_dir.glob('*/*'))
    #d_direcs = [x for x in d_direcs if (x.is_dir() and Path(x, 'duck_runs').exists())]

    #openduck_df = make_wqb_df(data_direcs=d_direcs, sample_size=10, num_samples=4)
    openduck_df = pd.read_csv(str(Path(res_dir.parent, "iridium_openduck_values.csv")), index_col=0)
    # print(openduck_df)
    val_csv_path = Path('/dls/science/groups/i04-1/software/mihaela/Data/duck_stuff/duck_paper_validation_values.csv')
    #vali_df = find_iridium_target(val_csv_path)
    final_df = compare_to_paper(val_csv_path, openduck_df)
    interaction_names = [f"{k}_{i}" for k, i in zip(final_df['PDBID'].values, final_df['Interaction'].values)]
    final_df['interaction_id'] = interaction_names
    final_df = final_df.sort_values(by='PDBID')
    final_df.to_csv(Path(res_dir, 'openduck_to_paper_comparison_full_iridum.csv'))
    """

    val_csv_path = Path('/dls/science/groups/i04-1/software/mihaela/Data/duck_stuff/seraphic_paper_values.csv')
    # vali_df = find_iridium_target(val_csv_path)
    res_dir = Path(
        '/home/jin76872/Desktop/Mih/Data/duck_stuff/openduck_seraphic_may_2021/seraphic_ducks')
    openduck_df = pd.read_csv(str(Path(res_dir.parent, "seraphic_openduck_values.csv")), index_col=0)
    final_df = compare_to_paper(val_csv_path, openduck_df)
    interaction_names = [f"{k}_{i}" for k, i in zip(final_df['PDBID'].values, final_df['Interaction'].values)]
    final_df['interaction_id'] = interaction_names
    final_df = final_df.sort_values(by='PDBID')
    final_df.to_csv(Path(res_dir, 'openduck_to_paper_comparison_full_seraphic.csv'))
    """
