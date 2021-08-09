import numpy as np
import os
import argparse

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

def get_wqb_simple(file_duck_dat):
    f = open(file_duck_dat,'r')
    data = []
    for line in f:
        a = line.split()
        data.append([float(a[1]), float(a[3]), float(a[5]), float(a[8])])
    f.close()
    data = np.array(data[1:])
    Work = data[:,3]
    Wqb_max = max(Work[400:])
    Wqb_min = min(Work[:400])
    Wqb_value = Wqb_max - Wqb_min
    return(Wqb_value, data, Wqb_min)


def get_Wqb_value_all(input_dir):
    file_list = []
    for fil in os.listdir(input_dir):
        if fil[-3:] == 'dat':
            file_list.append(fil)

    Wqb_values = []
    for fil in file_list:
        Wqb_data = get_wqb_simple(fil)
        Wqb_values.append(Wqb_data[0])

    Wqb = min(Wqb_values)
    return(Wqb)

def main():
    parser = argparse.ArgumentParser(description='Get WQB score from OpenDUck data')
    parser.add_argument('-d', '--dir', help='Directory with location of OpenDUck data')
    parser.add_argument('-l', '--ligand', help='Ligand in mol format')
    parser.add_argument('-o', '--output', help='Ligand output in mol forma, with wqb value')

    args = parser.parse_args()

    if args.dir:
        input_dir = args.dir
    else:
        input_dir = os.getcwd()
    
    wqb_val = get_Wqb_value_all(input_dir)

    if args.ligand:
        with open(args.ligand) as f:
            records = f.read().split('$$$$')
            print(records)
        if (len(records) > 2) or (len(records) == 2 and not records[1].isspace()):
            # if there is more than 1 record; 2 is ok if the second is whitespace
            raise IOError('The mol file contains multiple records.')
        else:
            wqb_str = "> <SCORE.DUCK_WQB>\n{}\n".format(wqb_val)
            records[0] += wqb_str
            print(records)
            with open(args.output, 'w') as f:
                f.write('$$$$'.join(records))
    else:
        print(wqb_val)


if __name__ == '__main__':
    main()

