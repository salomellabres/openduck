from pathlib import Path
import numpy as np
from scipy.stats import ttest_ind, ks_2samp
import pandas as pd
import matplotlib.pyplot as plt


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

def get_Wqb_value_all(dat_files_list, protein_name, plot_data=True):
    """
    Again taken with minor modifications from Maciej's repo
    :param dat_files_list:
    :param protein_name:
    :return:
    """

    Wqb_values = []

    for fil in dat_files_list:
        Wqb_data = get_Wqb_value(fil)
        Wqb_values.append(Wqb_data[0])

    Wqb = min(Wqb_values)
    if plot_data:
        fig = plt.figure(figsize=(7, 7))
        plt.plot(Wqb_data[1][:,0], Wqb_data[1][:,3]-Wqb_data[2])
        plt.xlabel('HB Distance (A)')
        plt.ylabel('Work (kcal/mol)')
        fig.suptitle(f"{protein_name}, Wqb: {round(Wqb,2)}")
        plt.savefig(Path(f"{protein_name}_wqb_plot.png"))
        plt.close()

    return(Wqb)

def wqb_comparison(dat_files_list_300, dat_files_list_325, protein_name, plot_data=True):
    """
    Plots a comparison between the Wqb value distributions at 300 and 325 and calculates a
    probablity that they come from the same distribution.
    Try t-test and Kolmogorov-Smirnov?
    :param dat_files_list_300:
    :param dat_files_list_325:
    :param protein_name:
    :return:
    """
    wqbs_300 = [get_Wqb_value(f)[0] for f in dat_files_list_300]
    wqbs_325 = [get_Wqb_value(i)[0] for i in dat_files_list_325]

    # Try the t-test with equal and unequal variances, as well as Kolmogorov-Smirnov
    ev_a, ev_b = ttest_ind(wqbs_300, wqbs_325, equal_var=True)
    un_a, un_b = ttest_ind(wqbs_300, wqbs_325, equal_var=False)
    ks_s, ks_p = ks_2samp(wqbs_300,  wqbs_325)

    if plot_data:
        fig, ax = plt.figure((7,7))
        ax.hist(wqbs_300, c='b', alpha=0.5)
        ax.hist(wqbs_325, c='r', alpha=0.5)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
        textstr = f""" t-test (equal variances): {ev_a}, p-value: {ev_b} \n t-test(unequal variances: {un_a}, p-value: {un_b} \n K-S 2samples: {ks_s}, p-value:{ks_p}"""
        ax.text(0.7, 0.7, textstr, bbox=props)
        fig.suptitle(protein_name)
        plt.show()
    return [[ev_a, ev_b], [un_a, un_b], [ks_s, ks_p]]

