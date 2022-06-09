import matplotlib
#matplotlib.use('Agg')
from pathlib import Path
import yaml

from scipy.signal import find_peaks, peak_widths
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import traceback

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

def make_all_data_dict(dir_paths, to_sim, num_runs=40):
    dat_dict = {}
    
    for i, dp in enumerate(dir_paths):
        print(to_sim[i])
        try:
            all_data = [get_Wqb_value(str(p)) for p in dp]
            dat_dict[to_sim[i]] = [(p, i) for p, i in zip(dp, all_data)]
        except IndexError:
            continue
        if len(all_data) < num_runs:
            continue
        
    return dat_dict

def find_outliers_median(otype_vals):
    med = np.median(otype_vals)
    diffs = otype_vals - med
    outs = np.where((np.abs(diffs) > 2*otype_vals.std()) & (np.abs(diffs) > 2.0))
    return outs[0], diffs[outs[0]], med

def find_outliers(otype_vals):
    # First, map the values onto a histogram that spans possible work trace space
    hdata = np.histogram(otype_vals, bins=100, range=(-5.0, 40))
    peaks = find_peaks(hdata[0], threshold=2, distance=6)
    if len(peaks[0]) > 1:
        print(peaks)
        widths = peak_widths(hdata[0], peaks[0])
        print(widths)
        num_samples = []
        for i in range(len(peaks[0])):
            #get the peak bins
            bottom_idx = int(np.floor(widths[2][i]))
            top_idx = int(np.ceil(widths[3][i]))
            vals = sum(hdata[0][bottom_idx:top_idx+1])
            num_samples.append(vals)
        #find the most populated peak
        pp_idx = peaks[0][np.argmax(np.array(num_samples))]
        pp_val = hdata[1][pp_idx]
        # Anything over 2std away is an ouitlier
        diffs = otype_vals - pp_val
        print(diffs.shape)
        outs = np.where((np.abs(diffs) > 2*otype_vals.std()) & (np.abs(diffs) > 2.0))
        #outs = np.where(((otype_vals < hdata[1][bottom_idx]) | (otype_vals > hdata[1][top_idx])) & (np.abs(diffs) > 2.0))
        return outs[0], diffs[outs[0]], pp_val
    else:
        return(find_outliers_median(otype_vals))
        
    

def get_outlying_trajectories(tar_interact):
   
    tar_paths = [x for x in Path("../iridium_ducks/iridium_dat_files").glob(f"{tar_interact}*")]
    direcs_to_sim = [tar_interact]
    int_dict = make_all_data_dict([tar_paths], to_sim=direcs_to_sim)
    m = int_dict[tar_interact]
    
    df = pd.DataFrame()
    df['trajectory'] = [m[i][0].name for i in range(len(m))]
    df['starts'] = np.array([m[i][1][1][:, 3][0]- m[i][1][2] for i in range(len(m))])
    df['minima'] = np.array([m[i][1][2] for i in range(len(m))])
    df['final'] = np.array([m[i][1][1][:, 3][-1]- m[i][1][2] for i in range(len(m))])
    
    outliers = []
    outlier_types = []
    outlier_distance = []
    peak_vals = []
    notes = []
    
    for otype in ['starts', 'minima', 'final']:
        #med = df[otype].median()
        #outlying = df[(diffs > 2*df[otype].std()) & (diffs > 1.0)]
        try:
            out, diffs, peak_val = find_outliers(df[otype].values)
            notes.extend([""]*len(out))
        except Exception as e:
            print(traceback.format_exc())
            out, diffs, peak_val = find_outliers_median(df[otype].values)
            notes.extend(["Peak detection failed. Peak value set to median"]*len(out))
            
        outlying=df.iloc[out]
        outliers.extend(outlying['trajectory'])
        outlier_types.extend([otype] * len(outlying))
        outlier_distance.extend(diffs)
        peak_vals.extend([peak_val]*len(outlying))
        
    out_df = pd.DataFrame()
    out_df['interaction'] = [tar_interact] * len(outliers)
    out_df['trajectories'] = outliers
    out_df['type'] = outlier_types
    out_df['distance_to_peak'] = [round(i,2) for i in outlier_distance]
    out_df['peak_val'] = [round(i,2) for i in peak_vals]
    out_df['notes'] = notes
    return out_df

if __name__ == "__main__":
    all_df = pd.read_csv("../openduck_to_paper_comparison_full_iridum.csv", index_col=0)
    #interactions = all_df['interaction_id'].values
    interactions = ['1yot_GLY216_N', '1l7f_GLU227_OE1', '1l7f_ARG118_NH1', '1n2j_GLN72_OE1', '1qhi_GLN125_NE2', '1tt1_THR143_OG1']
    out_trajs = []
    for i in interactions:
        try:
            out_trajs.append(get_outlying_trajectories(i))
        except Exception as e:
            print(traceback.format_exc())
            continue
        
    all_df = pd.concat(out_trajs)
    all_df.to_csv("outlying_trajectories_iridium.csv")