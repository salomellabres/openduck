#!/usr/bin/env python
import glob
import os
import numpy as np
import pandas as pd
import sys
import argparse
from rdkit import Chem

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#### Script from Maciej
def get_Wqb_value_AMBER(file_duck_dat):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[0]), float(a[1]), float(a[2]), float(a[3])])
    f.close()
    data = np.array(data)
    Work = data[:,3]
    #split it into segments of 200 points 
    num_segments = int(len(data)/200) 
    #alayze each segment to see if minimum in the segment is the local minimum
    #local minimum is the point with the lowest value of 200 neighbouring points
    #first local minumum is miminum used later to duck analysis
    Wqb_min_index = False
    for segment in range(num_segments):
        #detecting minium inthe segment
        sub_data = data[segment * 200 : (segment + 1) * 200]
        sub_Work = sub_data[:,3] 
        index_local = np.argmin(sub_Work)
        #segment of 200 points arround detected minimum
        index_global = index_local + segment * 200
        if index_global > 100:
            sub2_data = data[index_global - 100 : index_global + 101]
        else:
            sub2_data = data[0 : index_global + 101]
        sub2_Work = sub2_data[:,3]
        index_local2 = np.argmin(sub2_Work)
        if index_global < 100:
            if index_local2 == index_global:
                
                Wqb_min_index = index_global
            break
        else:
            if index_local2 == 100:
                Wqb_min_index = index_global
                break
    if not Wqb_min_index:
        eprint('Warning, could not find minima in {file_duck_dat}')
        return (0, data, 0)
    Wqb_min = Work[Wqb_min_index]
    sub_max_data = data[Wqb_min_index:]
    sub_max_Work = sub_max_data[:,3]
    Wqb_max_index = np.argmax(sub_max_Work)
    Wqb_max_index_global = Wqb_max_index + Wqb_min_index
    
    Wqb_max = max(sub_max_Work)
    
    Wqb_value = Wqb_max - Wqb_min
    
    return Wqb_value, data, Wqb_min

def get_Wqb_value_AMBER_all(prefix = 'DUCK', file = 'duck.dat'):
    folder = []
    for fol in os.listdir(os.getcwd()):
        if fol.startswith(prefix):
            folder.append(fol)
            
    Wqb_values = []
    for fol in folder:
        if os.path.isfile(fol+'/'+file):
            Wqb_data = get_Wqb_value_AMBER(fol+'/'+file)
            Wqb_values.append((fol, Wqb_data[0]))
    if len(Wqb_values):
        wqb = min([x[1] for x in Wqb_values])
    else:
        wqb = 'Nan'
    return wqb, Wqb_values

#######

def read_file(file):
    with open(file) as f:
        wqb = f.read()
        try: 
            wqb = float(wqb)
        except:
            wqb = calculate_wqb(os.path.join(*file.split('/')[:-1]))
    return wqb

def calculate_wqb(folder):
    #know where we are
    currdir = os.getcwd()

    #go to duck folder
    os.chdir(folder)

    wqb = get_Wqb_value_AMBER_all()
    os.chdir(currdir)
    return wqb

def iterate_systems(pattern):
    folders = glob.glob(pattern)
    r = {}
    for folder in folders:
        wqb = calculate_wqb(folder)
        r[folder] = wqb
    return r

def flatten_wqb_dict(info_dict):
    values = list(info_dict.values())
    k = ['System', *[x[0] for x in values[0][1]]]
    v = [list(info_dict.keys())]
    v.extend([[] for i in range(len(values[0][1]))])
    for i in range(len(values)):
        for j in range(len(values[0][1])):
            v[j+1].append(values[i][1][j][1])
    return {ks: vs for ks, vs in zip(k, v)}

def build_report_df(info_dict, mode='min'):
    if mode == 'min':
        df = pd.DataFrame({'System': list(info_dict.keys()), 'Wqb':[wqb[0] for wqb in info_dict.values()]})
    elif mode == 'all':
        flat_dict = flatten_wqb_dict(info_dict)
        df = pd.DataFrame(flat_dict)
    elif mode == 'avg':
        df = pd.DataFrame({'System': list(info_dict.keys()), 'Wqb':[wqb[0] for wqb in info_dict.values()],
                            'Average': [np.mean([x[1] for x in wqb[1]]) for wqb in info_dict.values()],
                            'SD': [np.std([x[1] for x in wqb[1]]) for wqb in info_dict.values()]})
    else:
        raise ValueError(f'{mode} is not a valid report mode. Try with min, all or avg.')
    return df

def get_mols_and_format(data_df, mode='min'):
    mols = []
    for i, row in data_df.iterrows():
        mol = Chem.MolFromPDBFile(os.path.join(row['System'], 'ligand.pdb'))
        mol.SetProp('_Name', row['System'])
        if mode == 'min':
            mol.SetProp('Wqb', str(row['Wqb']) )
        elif mode == 'all':
            data = list(row.values)
            data.pop(0) # pop system
            mol.SetProp('All Wqb' ,','.join([str(x) for x in data]))
        elif mode == 'avg':
            mol.SetProp('Wqb', str(row['Wqb']) )
            mol.SetProp('Wqb Avg', str(row['Average']) )
            mol.SetProp('Wqb sd', str(row['SD']) )
        mols.append(mol)
    return mols
def main():
    parser = argparse.ArgumentParser(description='Collect data from duck to report')
    parser.add_argument('-p', '--pattern', default='.', type=str, help='Bash wildcard pattern to find folders with duck data')
    parser.add_argument('-m', '--mode', type=str, default='min', help='Mode to compile the report [min | all | avg]')
    parser.add_argument('-o', '--output', default='stdout', help = 'Output file, default is printing report to stdout.')
    parser.add_argument('-of', '--output_format', default='tbl', type=str, help='Output format, [csv | sdf | tbl].')
    args = parser.parse_args()

    if args.output == 'stdout':
        args.output = sys.stdout
        
    wqb_info = iterate_systems(args.pattern)
    df = build_report_df(wqb_info, mode=args.mode)
    if args.output_format == 'csv':
        df.to_csv(args.output, index=False)
    elif args.output_format == 'tbl':
        df.to_csv(args.output, index=False, sep='\t')
    elif args.output_format == 'sdf':
        mols = get_mols_and_format(df, mode=args.mode)
        with Chem.SDWriter(args.output) as w:
            [w.write(mol) for mol in mols]
if __name__ =='__main__':
    main()
    