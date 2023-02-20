#!/usr/bin/env python

import glob
import pandas as pd
import numpy as np
from scipy.stats import shapiro, sem
import matplotlib.pyplot as plt
import seaborn as sns
from math import log, exp
import os
import pickle

def read_DUckdat(temp=None):
	# read SMD reports from amber (4th line is the work)
	# If temp not None, it only takes the specified temperature works
	# its important in order to accurately calculate the jarzynski dG

    work_values = {}
    pattern = 'DUCK*/duck.dat' 
    for dat_file in glob.glob(pattern):
        if temp == '325K' and temp not in dat_file:
            continue
        elif temp == '300K' and '325K' in dat_file:
            continue
        with open(dat_file) as fh:
            work_values[dat_file] = []
            RC = []
            for line in fh:
                line = line.split()
                work_values[dat_file].append(float(line[3]))
                RC.append(line[0])
    return pd.DataFrame(work_values, index = RC)

def normalize_by_series(df, index_threshold = 2500):
	# Normalize the work values to have the minima at 0 kcal/mol
	# Find the minima in the first half of steps (i.e 2500).
	# Should be tweaked but its the normal behaviour.
    import copy
    norm_df = copy.deepcopy(df)
    for col in norm_df.columns:
        norm_df[col] = norm_df[col]-min(norm_df[col][:index_threshold])
    return norm_df

def get_expavg_FD_df(work_df, T=300, calculate_FD=True):
    #Calcular Jar // adapted from Luca's script
    kB=0.0019872  #Units: kcal/mol/K
    B=1/(kB*T)
    final_df = pd.DataFrame(columns=['expavg', 'sqrtMSE', 'MSE', 'Bj', 'v', 'av', 'ai', 'Wdis', 'avg', 'variance', 'FD', 'sqrtMSE_FD'],
                            index = work_df.index)
    t_work = work_df.T
    expBwork  = t_work.applymap(lambda x : exp(x*(-B)))
    N = len(work_df.columns) # Should be 5000 for DUck

    for rc, work_step in work_df.iterrows():
        #initialize variables to prevent errors in saving empty values
        expavg, sqrtMSE, MSE, Bj, v, av, ai, Wdis, average, variance, FD, sqrtMSE_FD = (np.NAN for x in range(12))
        
        #Normal statistics
        average  = np.mean(work_step.array)
        variance = np.var(work_step)

        #Boltzman average for the step
        expavg   = (-log(np.mean(expBwork.loc[:,rc])))/B
        if calculate_FD: # I have added the option to skip it, as it was giving some overflow float problems with very big dissipations during sampling
            # Other parameter for the fluctuation dissipation formula (applied to the Boltzmann avg)
            # Dissipated work due to heat loss, Wdis = T*dS-dQ if Wdis >=0
            Wdis = 0.5*B*variance
            if np.isnan(Wdis) or Wdis < 0.02: # for very small Wdis we can neglect it
                ai, av = 1, 1
            else:
                ai = float(log(30*Wdis*B))/float(log(15*(exp(2*B*Wdis)-1)))
                av=float(log(100*Wdis*B))/float(log(50*(exp(2*B*Wdis)-1)))
            try:
                Bj = Wdis/N**ai
                v = variance/N**av
                MSE = Bj**2+v
                sqrtMSE = MSE**0.5
                # Proper fluctuation dissipation (FD) using Kullback-leibler divergence
                FD = average-Wdis
                MSE_FD = 2*Wdis/(B*N)+(2*Wdis**2)/(N-1)
                sqrtMSE_FD = MSE_FD**0.5
            except:
                #print(ai, av)
                pass



        # write out on the final_df
        final_df.loc[rc] = [expavg, sqrtMSE, MSE, Bj,v, av, ai, Wdis, average, variance, FD, sqrtMSE_FD]
    return final_df

def plot_expavg_FD(raw_data_df, FD_df):
    RC = raw_data_df.index.values
    fig, ax = plt.subplots()
    # Plot raw Wqb
    for i, duck_file in raw_data_df.T.iterrows():
        ax.plot(RC,duck_file.values, 'b', linewidth=0.5)

    # Plot expavg (Boltzmann avg) & FD (Fluctuation dissipation)	
    ax.plot(RC, FD_df['expavg'].values, 'g', linewidth=4)
    ax.plot(RC, FD_df['FD'].values,'r', linewidth=4)
    ax.set_xticks(np.append(RC[::1000], RC[-1]))


    ax.set_ylabel('Free Energy (kcal/mol)')
    ax.set_xlabel('Distance (\u212B)')
    fig.savefig('Wqb_plot_jarzynski.png') 

def shapiro_test(work_df):
	fig, ax = plt.subplots()
	#values = work_df.tail(1).values[0]
	values = np.array([work_df[coll].max() for coll in work_df.columns]) # get the real maximum for the shapiro, as it will be the one used
	statistic, p = shapiro(values)
	y, x, _ = ax.hist(values, bins = len(values)//2)
	ax.text(y = max(y)-1, x = max(x)-1.5, s ='Shapiro p = %s'%round(p,4))
	sns.despine()
	ax.set_xlabel('Wqb (kcal/mol)')
	fig.savefig('wqb_hist.pdf')
	if p >= 0.05:
		return True
	else:
		return False

def bootstrap_df(norm_data_list, sample_size = 20, samples=20, plot=True, temps=[300, 325]):
	# Sample by bootstrapping the works in order to see the convergence of jarzynski.
	# Each temperature is sampled separatedly and then merged when the jarzynski is already calculated.
	sample_dfs = []
	for temp, norm_df in zip(temps, norm_data_list):
		for i in range(samples):
			new_df = norm_df.sample(n=sample_size, replace=True,axis='columns')
			new_df.columns = [f'Sample_{x}'for x in range(sample_size)]
			sample_dfs.append(get_expavg_FD_df(new_df, T=temp, calculate_FD=False))
	flat_dict = {'CV':[],'Jarzynski':[]}
	for sample in sample_dfs:
		for i,row in sample.iterrows():
			flat_dict['CV'].append(float(i))
			flat_dict['Jarzynski'].append(row['expavg'])
	flat_df = pd.DataFrame(flat_dict)
	if plot:
		fig, ax = plt.subplots()
		sns.lineplot(data=flat_df, x='CV',y='Jarzynski', errorbar='sd', ax=ax)
		ax.set_ylabel('Free Energy (kcal/mol)')
		ax.set_xlabel('Distance (\u212B)')
		fig.savefig('bootstraped_WQB_plot.png')
	return flat_df, sample_dfs

def get_stats_from_bootstrapping(flat_df):
	CVs = list(norm_df.index)
	stats_dict = {'CV':[],'AvgJarzynski':[], 'sd':[], 'sem':[]}
	for step in CVs:
		sub_df = flat_df[flat_df['CV'] == float(step)]
		values = np.array(sub_df['Jarzynski'])
		stats_dict['CV'].append(float(step))
		stats_dict['AvgJarzynski'].append(np.mean(values))
		stats_dict['sd'].append(np.std(values))
		stats_dict['sem'].append(sem(values))
	stats_df = pd.DataFrame(stats_dict)
	return stats_df

def average_dataframes(dataf1, dataf2):
	new_df = pd.DataFrame(index=dataf1.index, columns=dataf1.columns)
	for i,r in dataf1.iterrows():
		for col,v300 in r.items():
			v325 = dataf2.loc[i,col]
			new_df.at[i,col] = (v300+v325)/2
	return new_df

def get_real_jarzynski_from_bootstrapping(sample_dfs, save=None, splitting_point=2500):
	# real calculated wqb should come from the maximums and not from the end points
	# as such, the reported values come from the following
	max_values = [max(list(sample_df['expavg'])[list(sample_df['expavg'])[:splitting_point].index(min(list(list(sample_df['expavg'])[:splitting_point]))):])
              for sample_df in sample_dfs]
	avg = np.mean(max_values)
	sd = np.std(max_values)
	sem_v = sem(max_values)
	if save:
		with open(save, 'w') as f:
			f.write('\t'.join([str(x) for x in [avg,sd,sem_v]]))
	return avg,sd,sem_v, max_values
def save_pickle(file, object):
    with open(file, 'wb') as fh:
        pickle.dump(object, fh)
def read_pickle(file):
    with open(file, 'r') as fh:
        obj = pickle.load(fh)
    return obj
def main():
    temperatures = [300, 325]
    norm_datas, FD_datas = [],[]
    for T in temperatures:
        WQB_df = read_DUckdat(temp=f'{T}K')
        norm_df = normalize_by_series(WQB_df)
        FD_df = get_expavg_FD_df(norm_df,T=T)
        norm_datas.append(norm_df)
        FD_datas.append(FD_df)
    norm_df = pd.concat(norm_datas, axis=1)
    CVs = list(norm_df.index())
    normality = shapiro_test(norm_df)
    FD_df = average_dataframes(*FD_datas)
    plot_expavg_FD(norm_df, FD_df)
    flat_bootstrapped_df, sampled_dfs = bootstrap_df(norm_datas)
    save_pickle('resampling.pickle',(flat_bootstrapped_df, sampled_dfs))
    stats_df = get_stats_from_bootstrapping(flat_bootstrapped_df, CVs)
    avg,sd,sem_v, max_stats = get_real_jarzynski_from_bootstrapping(sampled_dfs, save='jarz_sd_sem.tbl')
    print(max_stats)
	#save data in csv
    FD_df.to_csv('fluctuation_dissipation.csv')
    stats_df.to_csv('bootstrapped_stats.csv')
if __name__=='__main__':  
    # usage: launch script in the LIG_target folder to obtain jarzynski report
    #test_folder = '/home/aserrano/Documents/openduck_validation/Iridium/1ml1/GLY173_N/SMIRNOFF'
    #os.chdir(test_folder)
    main()