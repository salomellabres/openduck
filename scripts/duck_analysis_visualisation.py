import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

def get_Wqb_value(file_duck_dat):
    """
    Taken from Maciej's MOD repo: https://github.com/MaciejMajew/OpenDUck/blob/master/duck_final/getWqbValue.py
    :param file_duck_dat:
    :return:
    """
    f = open(file_duck_dat, 'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:, 3]
    # split it into segments of 200 points
    num_segments = int(len(data) / 200)
    # alayze each segment to see if minimum in the segment is the local minimum
    # local minimum is the point with the lowest value of 200 neighbouring points
    # first local minumum is miminum used later to duck analysis
    for segment in range(num_segments):
        # detecting minium inthe segment
        sub_data = data[segment * 200: (segment + 1) * 200]
        sub_Work = sub_data[:, 3]
        index_local = np.argmin(sub_Work)
        # segment of 200 points arround detected minimum
        index_global = index_local + segment * 200
        if index_global > 100:
            sub2_data = data[index_global - 100: index_global + 101]
        else:
            sub2_data = data[0: index_global + 101]
        sub2_Work = sub2_data[:, 3]
        index_local2 = np.argmin(sub2_Work)
        if index_global < 100:
            if index_local2 == index_global:
                Wqb_min_index = index_global
            break
        else:
            if index_local2 == 100:
                Wqb_min_index = index_global
                break
    # For some of my traces, the local minimum doesn't work. IN that case, try defaulting to the global minimum.
    if not 'Wqb_min_index' in locals():
        print('Cannot find minimum Wqb index for file {} \n Setting it to the minimium work observed in the trace'.format(file_duck_dat))
        Wqb_min_index = np.argmin(Work)

    Wqb_min = Work[Wqb_min_index]
    sub_max_data = data[Wqb_min_index:]
    sub_max_Work = sub_max_data[:, 3]

    Wqb_max = max(sub_max_Work)

    Wqb_value = Wqb_max - Wqb_min
    return (Wqb_value, data, Wqb_min)

def make_all_data_dict(dir_paths, num_runs=40):
    dat_dict = {}
    to_sim = list(set([(x).name.split("_smd")[0] for x in all_runs]))
    for i, dp in enumerate(dir_paths):
        print(to_sim[i])
        try:
            all_data = [get_Wqb_value(p) for p in dp]
            dat_dict[to_sim[i]] = [(p, i) for p, i in zip(dp, all_data)]
        except IndexError:
            continue
        if len(all_data) < num_runs:
            continue
        
    return dat_dict

def make_wqb_images(dat_dict, save_dir):
    if not Path(save_dir).exists():
        Path(save_dir).mkdir()
        
    for interaction, run in dat_dict.items():
        all_data=[run[i][1] for i in range(len(run))]
        plt.figure()
        for it, d in enumerate(all_data):
            dists = 1*d[1][:, 0]
            Work = d[1][:,3]-d[2]
            plt.plot(dists, Work)
            plt.text(x=dists[-1], y=Work[-1], s=str(round(d[0], 2)) + ' ' + str(it))
            plt.title(f'{interaction} Wqb: {round(np.min([ad[0] for ad in all_data]), 2)}')
            plt.xlabel('Distance (Angstrom)')
            plt.ylabel('Wqb (kcal/mol)')
        plt.savefig(str(Path(save_dir, f"{interaction}.png")))
        plt.close()
            
def make_wqb_df(dat_dict, sample_size=20, num_samples=10):
    np.random.seed(3)
    min_mean_dict = {}
    stdev_dict = {}
    
    for interaction, run in dat_dict.items():
        all_data=[run[i][1][0] for i in range(len(run))]
        print(all_data)
        wqb_mins = []
        try:
            for q in range(num_samples):
                sub = np.random.choice(all_data, sample_size, replace=False)
                wqb_mins.append(np.min(sub))
            mean_min = np.mean(wqb_mins)
            min_mean_dict[interaction] = mean_min
            stdev_dict[interaction] = np.array(wqb_mins).std()
        except Exception as e:
            print(e)
                
    pdb_list = []
    interactions = []
    openduck_2021_list = []
    mean_min_list = []
    stdev_list = []

    for key, run1 in dat_dict.items():
        all_dat = [run1[i][1][0] for i in range(len(run1))]
        try:
            mean_min = min_mean_dict[key]
            mean_min_list.append(mean_min)
            stdev = stdev_dict[key]
            stdev_list.append(stdev)
            openduck_2021_list.append(round(np.min([all_dat]), 2))
            pdb_id = key.split('_')[0]
            pdb_list.append(pdb_id)
            interaction = '_'.join(key.split('_')[1:])
            interactions.append(interaction)
        except KeyError as e:
            print('EXCEPTION IN SECOND LOOP')
            print(e)

    wqb_df = pd.DataFrame()
    wqb_df['PDBID'] = pdb_list
    wqb_df['Interaction'] = interactions
    wqb_df['openduck_2021'] = openduck_2021_list
    wqb_df['openduck_mean_min'] = np.round(mean_min_list, 2)
    wqb_df['mean_min_stdev'] = np.round(stdev_list, 2)

    return wqb_df

def make_all_wqb_runs_df(dat_dict):
    interact_list = []
    dat_file_list = []
    all_wqb_list = []

    for inter in dat_dict.keys():
        for di in dat_dict[inter]:
            interact_list.append(inter)
            dat_file_list.append(di[0].name)
            all_wqb_list.append(di[1][0])
    all_df = pd.DataFrame(data=np.array([interact_list, dat_file_list, all_wqb_list]).T, columns= ['Interaction', 'dat_file','Wqb'])
    return all_df

def get_chunk_residue(original_pdb, chunk_pdb, target_residue_name, target_residue_number):
    pdb_univ = mda.Universe(original_pdb, format='PDB')
    chunk_univ = mda.Universe(chunk_pdb, format='PDB')

    # Get the calpha atom of the pdb target residue:
    target_res_calpha = pdb_univ.select_atoms(f"resname {target_residue_name} and resnum {target_residue_number} and name CA")

    # Get all the chunk calphas of the same residue type
    chunk_resis = chunk_univ.select_atoms(f"resname {target_residue_name} and name CA")

    # Get the distances between the trget calpha and those of the chunk
    distances = np.array([mda.analysis.rms.rmsd(target_res_calpha[0].position, chunk_resis[i].position) for i in range(len(chunk_resis))])

    smallest_distance = distances[np.argmin(distances)]
    if smallest_distance > 0.5:
        print(f'Alignment for {original_pdb} displaced')

    # Target chunk atom should be the closest
    target_chunk_atom = chunk_resis[np.argmin(distances)]

    return (target_chunk_atom.resname, target_chunk_atom.resid)


def vmd_template(pdb_file, dcd_file, ligand_selection, residue_selection):
    lig_and_res_sel = "{" + ligand_selection + " and (" + residue_selection + ")}"
    ligand_selection = "{"+ ligand_selection + "}"
    residue_selection = "{"+ residue_selection + "}"

    print(lig_and_res_sel)

    vmd_str = """
    #!/usr/local/bin/vmd
    # VMD script written by save_state $Revision: 1.47 $
    # VMD version: 1.9.3
    set viewplist {}
    set fixedlist {}
    proc vmdrestoremymaterials {} {
      set mlist { Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy BlownGlass GlassBubble RTChrome }
      set mymlist [material list]
      foreach mat $mlist {
        if { [lsearch $mymlist $mat] == -1 } { 
          material add $mat
        }
      }
    }
    vmdrestoremymaterials
    """

    vmd_str += """
    # Display settings

    mol new {0} type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol addfile {1} type dcd first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

    mol delrep 0 top
    mol representation Licorice 0.300000 12.000000 12.000000
    mol color Name
    mol selection {2}
    mol material Opaque
    mol addrep top
    mol selupdate 0 top 0
    mol colupdate 0 top 0
    mol scaleminmax top 0 0.000000 0.000000
    mol smoothrep top 0 0
    """.format(pdb_file, dcd_file, ligand_selection)

    vmd_str += """
    mol drawframes top 0 {now}

    mol representation Lines 1.000000
    mol color Name
    mol selection {protein}
    mol material Opaque
    mol addrep top
    mol selupdate 1 top 0
    mol colupdate 1 top 0
    mol scaleminmax top 1 0.000000 0.000000
    mol smoothrep top 1 0
    mol drawframes top 1 {now}
    """

    vmd_str += """
    mol representation Licorice 0.300000 12.000000 12.000000
    mol color Name
    mol selection {0}
    mol material Opaque
    mol addrep top
    mol selupdate 2 top 0
    mol colupdate 2 top 0
    mol scaleminmax top 2 0.000000 0.000000
    mol smoothrep top 2 0
    """.format(residue_selection)

    vmd_str += """
    mol drawframes top 2 {now}
    mol representation HBonds 3.000000 20.000000 10.000000
    mol color Name
    """

    vmd_str += """
    mol selection { not water}
    mol material Opaque
    mol addrep top
    mol selupdate 3 top 0
    mol colupdate 3 top 0
    mol scaleminmax top 3 0.000000 0.000000
    mol smoothrep top 3 0
    """

    vmd_str += """
    mol drawframes top 3 {now}
    """

    vmd_str += "mol rename top {}".format(Path(pdb_file).name)
    vmd_str += """
    lappend viewplist [molinfo top]
    set topmol [molinfo top]
    # done with molecule 0
    foreach v $viewplist {
      molinfo $v set {center_matrix rotate_matrix scale_matrix global_matrix} $viewpoints($v)
    }
    foreach v $fixedlist {
      molinfo $v set fixed 1
    }
    unset viewplist
    unset fixedlist
    mol top $topmol
    unset topmol

    label textsize 1.0
    pbc wrap -center com -centersel "protein" -compound residue -all 
    pbc wrap -center com -centersel "protein" -compound residue -all 
    """
    #return vmd_str + vmd_str1 + vmd_str2 + vmd_str3 + vmd_str4 + vmd_str5 + vmd_str6
    return vmd_str
    

if __name__ == "__main__":
    #First, get the openduck values 
    
    all_runs = list(Path('.').glob('dat_files/*.dat'))
    to_sim = list(set([(x).name.split("_smd")[0] for x in all_runs]))
    dir_paths = [list(Path('.').glob(f'dat_files/{p}_smd_*.dat')) for p in to_sim]
    duck_dict = make_all_data_dict(dir_paths, num_runs=10)
    
    # Save the calculated values
    all_runs_df = make_all_wqb_runs_df(duck_dict)
    all_runs_df.to_csv('openduck_all_runs.csv')
    iridium_df = make_wqb_df(duck_dict, sample_size=10, num_samples=4)
    iridium_df.to_csv('openduck_values.csv')
    
    image_dir = Path('duck_images')
    make_wqb_images(duck_dict, image_dir)

    tar_dir = Path('.').resolve()
    for tp in [x for x in tar_dir.glob('*') if (x.is_dir() and '.ipy' not in x.name)]:
        try:
            with open(Path(tp, 'run.yaml')) as f:
                yaml_dict = yaml.load(f, Loader=yaml.FullLoader)
            original_pdb = Path(tp, yaml_dict["apo_pdb_file"])
            chunk_pdb = Path(tp, "complex.pdb")
            resname = yaml_dict["prot_int"].split('_')[1]
            resnum = yaml_dict["prot_int"].split('_')[2]
            print(resname, resnum)
            chunk_resname, chunk_resnum = get_chunk_residue(original_pdb, chunk_pdb,resname, resnum)
            print(chunk_resname, chunk_resnum)

            for traj in tp.glob('duck_runs/*.dcd'):
                pdb_file = str(traj).replace('.dcd', '.pdb')
                vmd_resname = "resname " + chunk_resname + ' and  resid ' + str(chunk_resnum)
                lig_name = 'resname UNL'
                vmd_script = vmd_template(pdb_file, str(traj), lig_name, vmd_resname)
                vmd_script_path = Path(str(traj).replace('.dcd', '.vmd'))
                vmd_script_path.write_text(vmd_script)

        except Exception as e:
            print(e)
            continue
    

    
    





