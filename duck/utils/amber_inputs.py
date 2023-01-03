import os

def write_string_to_file(file,string):
    with open(file, 'w') as fh:
        fh.write(string)

def write_min_and_equil_inputs(chunk_residues, interaction, hmr=False):
    # defining strings to write
    min_str = f"""&cntrl
imin=1, maxcyc=2000,
ntpr=100,
ntr=1,
restraintmask=':{chunk_residues} & !@H=',
restraint_wt=1.0,
/
    """
    heating_init=f"""&cntrl
imin=0,
ntx=1,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntp=0, ntc=2,
ntb=1, ntf=2, cut=9.0,
ntt=3, temp0=150, tempi=100, ig=-1, gamma_ln=4.0,
nstlim= 50000, dt=0.002,
ntr=1,
restraintmask=':{chunk_residues} & !@H=',
restraint_wt=1.0,
nmropt=1,
&end
&wt type='END' /
DISANG=dist_md.rst
"""
    heating = """&cntrl
imin=0,
ntx=5, irest=1,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntp=0, ntc=2,
ntb=1, ntf=2, cut=9.0, 
ntt=3, temp0={temp}, ig=-1,  gamma_ln=4.0,
nstlim= 50000, dt=0.002,
ntr=1,
restraintmask=':{chunk_residues} & !@H=', 
restraint_wt=1.0,
nmropt=1,
&end
&wt type='END' /
DISANG=dist_md.rst
"""
    time_step = '0.002'
    iterations='500000'
    if hmr:
        time_step = '0.004'
        iterations= '250000'
    eq_str = f"""&cntrl
imin=0,
ntx=5, irest=1,
iwrap=0,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntp=1, ntc=2, taup=2.0, 
ntb=2, ntf=2, cut=9.0,
ntt=3, temp0=300.0, ig=-1,  gamma_ln=4.0,
nstlim={iterations}, dt={time_step},
ntr=1,
restraintmask=':{chunk_residues} & !@H=', 
restraint_wt=1.0,
nmropt=1,
&end
&wt type='END' /
DISANG=dist_md.rst
    """
    prot_idx, lig_idx, pairmeandistance_i = interaction
    dist_str = f"""#Prevent dissociation by penalizing interaction break (>3.0A)
&rst iat={prot_idx+1},{lig_idx+1}, r2=2.00, r3=3.00, r4=4.00, rk2=1.0, rk3=10.0, /
    """

    # writing inputs for minimization, equilibration and heating
    print('Writing min and eq inputs')
    write_string_to_file('1_min.in', min_str)
    write_string_to_file('2_heat150.in', heating_init)
    write_string_to_file('2_heat200.in', heating.format(temp='200', chunk_residues=chunk_residues))
    write_string_to_file('2_heat250.in', heating.format(temp='250', chunk_residues=chunk_residues))
    write_string_to_file('2_heat300.in', heating.format(temp='300', chunk_residues=chunk_residues))
    write_string_to_file('3_eq.in', eq_str)
    write_string_to_file('dist_md.rst', dist_str)

def write_md_inputs(chunk_residues, interaction, hmr=False):
    time_step = '0.002'
    iterations = '500000'
    top = '{top}'
    if hmr:
        time_step = '0.004'
        iterations= '250000'
        top = 'HMR_'+top
    md_str =f"""&cntrl
ntx=5, irest=1,
iwrap=0,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntc=2, ntf=2,
ntb=1, cut=9.0,
ntt=3, temp0=300.0, gamma_ln=4.0, ig=-1,
nstlim={iterations}, dt={time_step},
ntr=1,
restraintmask=':{chunk_residues} & !@H=', 
restraint_wt=1.0,
nmropt=1,
&end
&wt type='END' /
DISANG=dist_md.rst
    """
    if not os.path.isfile('dist_md.rst'):
        prot_idx, lig_idx, pairmeandistance_i = interaction
        dist_str = f"""#Prevent dissociation by penalizing interaction break (>3.0A)\n&rst iat={prot_idx+1},{lig_idx+1}, r2=2.00, r3=3.00, r4=4.00, rk2=1.0, rk3=10.0, /"""
        write_string_to_file('dist_md.rst', dist_str)
    write_string_to_file('md.in', md_str)

def write_smd_inputs(chunk_residues, interaction, hmr=False):
    time_step = '0.002'
    iterations = '250000'
    savefreq = '50'
    if hmr:
        time_step = '0.004'
        iterations = '125000'
        savefreq = '25'
    smd_str=f"""&cntrl
ntx = 5, irest=1,
iwrap=0,
ntb=1,
ntt=3, temp0=300.0, gamma_ln=4.0,
nstlim={iterations}, dt={time_step},
ntc=2, ntf=2, cut=9.0,
ntxo=1, ntpr=2000, ntwx=0, ntwe=1000, ntwr=0, ioutfm=1,
jar=1,
ntr=1, restraintmask=':{chunk_residues} & !@H=', restraint_wt=1.0,
/
&wt type='DUMPFREQ', istep1={savefreq} /
&wt type='END'   /
DISANG=../dist_duck.rst
DUMPAVE=duck.dat
LISTIN=POUT
LISTOUT=POUT
    """
    smd_325_str=f"""&cntrl
ntx = 5, irest=1,
iwrap=0,
ntb=1,
ntt=3, temp0=325.0, gamma_ln=4.0,
nstlim={iterations}, dt={time_step},
ntc=2, ntf=2, cut=9.0,
ntxo=1, ntpr=2000, ntwx=0, ntwe=1000, ntwr=0, ioutfm=1,
jar=1,
ntr=1, restraintmask=':{chunk_residues} & !@H=', restraint_wt=1.0,
/
&wt type='DUMPFREQ', istep1={savefreq} /
&wt type='END'   /
DISANG=../dist_duck.rst
DUMPAVE=duck.dat
LISTIN=POUT
LISTOUT=POUT
    """
    prot_idx, lig_idx, pairmeandistance_i = interaction
    dist_duck_str=f"""#change distance from (2.50) to unbound (5.00)
&rst iat={prot_idx+1},{lig_idx+1}, r2=2.50, rk2=50.00, r2a=5.00, /
    """
    write_string_to_file('duck.in', smd_str)
    write_string_to_file('duck_325K.in', smd_325_str)
    write_string_to_file('dist_duck.rst', dist_duck_str)

def extract_residuenumbers(structure):
    r_id = []
    for r in structure.residues:
        if r.name != 'WAT' and r.name != 'UNL' and r.name != 'NA' and r.name != 'CL':
            r_id.append(int(r.number)+1)
    chunk_residues = re_range(r_id)
    return chunk_residues
def re_range(lst):
    # from here https://stackoverflow.com/questions/9847601/convert-list-of-numbers-to-string-ranges
    # not sure if it works super well!!!!
    n = len(lst)
    result = []
    scan = 0
    while n - scan > 2:
        step = 1
        if lst[scan + 2] - lst[scan + 1] != step:
            result.append(str(lst[scan]))
            scan += 1
            continue

        for j in range(scan+2, n-1):
            if lst[j+1] - lst[j] != step:
                result.append('{}-{}'.format(lst[scan], lst[j]+1))
                scan = j+1
                break
        else:
            result.append('{}-{}'.format(lst[scan], lst[n-1]+1))
            return ','.join(result)

    if n - scan == 1:
        result.append(str(lst[scan]))
    elif n - scan == 2:
        result.append(','.join(map(str, lst[scan:])))
    return ','.join(result)

def write_getWqbValues():
    script_string = """#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import os
import sys

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
        print('The minima could not be found, setting the Wqb to 0')
        exit(0)
    Wqb_min = Work[Wqb_min_index]
    sub_max_data = data[Wqb_min_index:]
    sub_max_Work = sub_max_data[:,3]
    Wqb_max_index = np.argmax(sub_max_Work)
    Wqb_max_index_global = Wqb_max_index + Wqb_min_index
    
    Wqb_max = max(sub_max_Work)
    
    Wqb_value = Wqb_max - Wqb_min
    
    return(Wqb_value, data, Wqb_min)

def get_Wqb_value_AMBER_all(prefix = 'DUCK', file = 'duck.dat'):
    folder = []
    for fol in os.listdir(os.getcwd()):
        if fol.startswith(prefix):
            folder.append(fol)
            
    Wqb_values = []
    for fol in folder:
        if os.path.isfile(fol+'/'+file):
            Wqb_data = get_Wqb_value_AMBER(fol+'/'+file)
            Wqb_values.append(Wqb_data[0])

    Wqb = min(Wqb_values)
    return(Wqb)
    

if __name__ == '__main__':
    if len(sys.argv) > 2:
        print(sys.argv)
        print(get_Wqb_value_AMBER_all(sys.argv[1],sys.argv[2]))
    elif len(sys.argv) > 1:
        print(get_Wqb_value_AMBER_all(sys.argv[1]))
    else:
        print(get_Wqb_value_AMBER_all())"""
    write_string_to_file('getWqbValues.py', script_string)

def write_all_inputs(structure,interaction, hmr=False):
    chunk_residues = extract_residuenumbers(structure)
    write_min_and_equil_inputs(chunk_residues, interaction, hmr=hmr)
    write_md_inputs(chunk_residues, interaction, hmr=hmr)
    write_smd_inputs(chunk_residues, interaction,hmr=hmr)

def write_queue_template(template, hmr = False, replicas=5, wqb_threshold=7, array_limit=False):
    top = 'system_complex.prmtop'
    if hmr: top = 'HMR_'+top
    functions_str="""##### FUNCTIONS ####
#Function adapted from 'submit_duck_smd_gpu.csh' of the DUck std pipeline
prepare_duck_and_launch(){{
   nustart=$1
   nuend=$2
   temp=$3
   nu=$nustart
   while (($nu <= $nuend)); do
      if [ "$temp" == '300K' ]; then
         dir=DUCK_${nu}
         mkdir $dir
         cd $dir
         if [ "$nu" == "0" ]; then
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../{top} -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../{top} -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../{top} -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../{top} -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi
         cd ..  
      elif [ "$temp" == '325K' ]; then
         dir=DUCK_325K_${nu}
         mkdir $dir
         cd $dir
	 if [ "$nu" == "0" ]; then
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../{top} -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../{top} -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../{top} -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../{top} -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi 
         cd ..
      fi
      nu=$((nu+1))
   done

}}

# Function to check if WQB is lower than 0.1 using getWqbValues.py
# getWqbValues.py is a script from Maciej modified.
# We use this one instead of the R version, as R is not available in the IQTC

check_WQB(){{
   wqb_limit=$1
   lowest_wqb=$(python getWqbValues.py)
   echo $lowest_wqb > wqb.log
   are_we_above_the_limit=$(echo "$lowest_wqb < $wqb_limit" | bc )
   if [ "$are_we_above_the_limit" == "1" ]; then
      echo "Wqb lower than ${wqb_limit}, stoping DUck run"
      cp -r ./* $LIG_TARGET/
      exit 
   fi
}}
    """.format(top=top, nu='{nu}',wqb_limit='{wqb_limit}')

    commands_str="""#### Runing Duck ####
# Minimization&Equilibration
pmemd.cuda -O -i 1_min.in -o min.out -p {top} -c system_complex.inpcrd -r min.rst -ref system_complex.inpcrd
pmemd.cuda -O -i 2_heat150.in -o 2_heat150.out -p {top} -c min.rst -r  2_heat150.rst -x 2_heat150.nc -ref system_complex.inpcrd
pmemd.cuda -O -i 2_heat200.in -o 2_heat200.out -p {top} -c 2_heat150.rst -r 2_heat200.rst -x 2_heat200.nc -ref 2_heat150.rst
pmemd.cuda -O -i 2_heat250.in -o 2_heat250.out -p {top} -c 2_heat200.rst -r 2_heat250.rst -x 2_heat250.nc -ref 2_heat200.rst
pmemd.cuda -O -i 2_heat300.in -o 2_heat300.out -p {top} -c 2_heat250.rst -r 2_heat300.rst -x 2_heat300.nc -ref 2_heat250.rst
pmemd.cuda -O -i 3_eq.in -o 3_eq.out -p {top} -c 2_heat300.rst -r 3_eq.rst -x 3_eq.nc -ref 2_heat300.rst -e 3_eq.ene

#Launch DUck 0 and check wqb
prepare_duck_and_launch 0 0 300K
check_WQB $min_wqb

#Launch DUck_325K 0 and check wqb
prepare_duck_and_launch 0 0 325K
check_WQB $min_wqb

#For each replica wanted do: MD, prepare SMD & launch SMD
for ((i=1;i<=$replicas;++i)); do
   if [ "$i" == "1" ]; then
      pmemd.cuda -O -i md.in -o md1.out -p {top} -c 3_eq.rst -r md1.rst -x md1.nc -ref 3_eq.rst
   else
      pmemd.cuda -O -i md.in -o md${i}.out -p {top} -c md$((i-1)).rst -r md${i}.rst -x md${i}.nc -ref 3_eq.rst
   fi

   prepare_duck_and_launch $i $i 300K
   prepare_duck_and_launch $i $i 325K
   check_WQB $min_wqb

done
    """.format(top=top, i='{i}')
    params_str = f"""
#### PARAMS ####
replicas={replicas}
min_wqb={wqb_threshold}
    """
    if not array_limit:
        queue_strings = {'Slurm': f"""#!/bin/bash
#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1      

{functions_str}

#### Modules ####
#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB
module load amber/20

{params_str}

{commands_str}
exit
        """,
        'SGE':f"""#!/bin/bash
#$ -N DUck_queue          # The name of the job, can be whatever makes sense to you
#$ -S /bin/bash          # Force sh if not Sun Grid Engine default shell
#$ -cwd                 # The batchsystem should use the current directory as working directory.
#$ -q fartorgpu.q            # Queue name where the job should be placed into.
#$ -o DUck.q.o             # Redirect output stream to this file.
#$ -e DUck.q.e             # Redirect error stream to this file.
#$ -l h_rt=15:00:00 # Time limit
#$ -pe gpu 1
#$ -m e

{functions_str}

#### Modules ####
#How to check modules in the queue
#module_fartor av

#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB)
. /etc/profile
module load amber/20_cuda9.0_ompi

#Necessary to use a free GPU
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`

#### Coping the files to node ####
#Things will need to run in $TMPDIR
LIG_TARGET=$PWD
cp -r $LIG_TARGET/* $TMPDIR
cd $TMPDIR

#Remove local output files, as it removes the queue ones when copying back
rm DUck.q.o DUck.q.e

#Where is the calculation being done?
echo "TMPDIR is $TMPDIR"

{params_str}

{commands_str}

#### Coping the files back to local ####
cp -r ./* $LIG_TARGET/
cd $LIG_TARGET

exit
        """}
    else:
        queue_strings = {'Slurm': f"""#!/bin/bash
#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1
#SBATCH --array=1-{array_limit}      

{functions_str}

#### Modules ####
#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB
module load amber/20
cd LIG_target_$SLURM_ARRAY_TASK_ID
{params_str}

{commands_str}
cd ..
exit
        """,
        'SGE':f"""#!/bin/bash
#$ -N DUck_queue          # The name of the job, can be whatever makes sense to you
#$ -S /bin/bash          # Force sh if not Sun Grid Engine default shell
#$ -cwd                 # The batchsystem should use the current directory as working directory.
#$ -q fartorgpu.q            # Queue name where the job should be placed into.
#$ -o DUck.q.o             # Redirect output stream to this file.
#$ -e DUck.q.e             # Redirect error stream to this file.
#$ -l h_rt=15:00:00 # Time limit
#$ -pe gpu 1
#$ -m e
#$ -t 1-{array_limit}

{functions_str}

#### Modules ####
#How to check modules in the queue
#module_fartor av

#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB)
. /etc/profile
module load amber/20_cuda9.0_ompi

#Necessary to use a free GPU
export CUDA_VISIBLE_DEVICES=`cat $TMPDIR/.gpus`

#### Coping the files to node ####
#Things will need to run in $TMPDIR
LIG_TARGET=$PWD/LIG_target_$SGE_TASK_ID
cp -r $LIG_TARGET/* $TMPDIR
cd $TMPDIR

#Remove local output files, as it removes the queue ones when copying back
rm DUck.q.o DUck.q.e

#Where is the calculation being done?
echo "TMPDIR is $TMPDIR"

{params_str}

{commands_str}

#### Coping the files back to local ####
cp -r ./* $LIG_TARGET/
cd $LIG_TARGET

exit
        """}

    print('Writing queue files')
    if template not in queue_strings:
        print('Warning wrong queue template. Only {} accepted'.format(list(queue_strings.keys())))
    if array_limit: write_string_to_file('array_duck_queue.q', queue_strings[template])
    else: write_string_to_file('duck_queue.q', queue_strings[template])
    write_getWqbValues()

if __name__=='__main__':
    import sys
    if len(sys.argv) != 3:
        print(f'usage: python {sys.argv[0]} template hmr_boolean')
    template = sys.argv[1]
    hmr = sys.argv[2]
    write_queue_template(template, hmr=hmr)

