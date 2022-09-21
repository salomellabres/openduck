import os

def write_string_to_file(file,string):
    with open(file, 'w') as fh:
        fh.write(string)

def write_min_and_equil_inputs(chunk_residues, interaction):
    # defining strings to write
    min_str = f"""&cntrl
imin=1, maxcyc=10000,
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
iwrap=1,
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
iwrap=1,
ntr=1,
restraintmask=':{chunk_residues} & !@H=', 
restraint_wt=1.0,
nmropt=1,
&end
&wt type='END' /
DISANG=dist_md.rst
"""
    eq_str = f"""&cntrl
imin=0,
ntx=5, irest=1,
iwrap=1,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntp=1, ntc=2, taup=2.0, 
ntb=2, ntf=2, cut=9.0,
ntt=3, temp0=300.0, ig=-1,  gamma_ln=4.0,
nstlim=500000, dt=0.002,
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

def write_md_inputs(chunk_residues, interaction):
    md_str =f"""&cntrl
ntx=5, irest=1,
iwrap=1,
ntxo=1, ntpr=2000, ntwx=0, ntwv=0, ntwe=0, ntwr=0, ioutfm=1,
ntc=2, ntf=2,
ntb=1, cut=9.0,
ntt=3, temp0=300.0, gamma_ln=4.0, ig=-1,
nstlim=250000, dt=0.002,
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

def write_smd_inputs(chunk_residues, interaction):
    smd_str=f"""ntx = 5, irest=1,
iwrap=1,
ntb=1,
ntt=3, temp0=300.0, gamma_ln=4.0,
nstlim=250000, dt=0.002,
ntc=2, ntf=2, cut=9.0,
ntxo=1, ntpr=2000, ntwx=0, ntwe=1000, ntwr=0, ioutfm=1,
jar=1,
ntr=1, restraintmask=':{chunk_residues} & !@H=', restraint_wt=1.0,
/
&wt type='DUMPFREQ', istep1=50 /
&wt type='END'   /
DISANG=../dist_duck.rst
DUMPAVE=duck.dat
LISTIN=POUT
LISTOUT=POUT
    """
    smd_325_str=f"""ntx = 5, irest=1,
iwrap=1,
ntb=1,
ntt=3, temp0=325.0, gamma_ln=4.0,
nstlim=250000, dt=0.002,
ntc=2, ntf=2, cut=9.0,
ntxo=1, ntpr=2000, ntwx=0, ntwe=1000, ntwr=0, ioutfm=1,
jar=1,
ntr=1, restraintmask=':{chunk_residues} & !@H=', restraint_wt=1.0,
/
&wt type='DUMPFREQ', istep1=50 /
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
            result.append('{}-{}'.format(lst[scan], lst[j]+1))
            return ','.join(result)

    if n - scan == 1:
        result.append(str(lst[scan]))
    elif n - scan == 2:
        result.append(','.join(map(str, lst[scan:])))
    return ','.join(result)

def write_all_inputs(structure,interaction):
    chunk_residues = extract_residuenumbers(structure)
    write_min_and_equil_inputs(chunk_residues, interaction)
    write_md_inputs(chunk_residues, interaction)
    write_smd_inputs(chunk_residues, interaction)

def write_queue_template(template):
    queue_strings = {'Slurm': """#!/bin/bash
#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1      
           
##### FUNCTIONS ####
#Function adapted from 'submit_duck_smd_gpu.csh' of the DUck std pipeline
prepare_duck_and_launch(){
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
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi
         cd ..  
      elif [ "$temp" == '325K' ]; then
         dir=DUCK_325K_${nu}
         mkdir $dir
         cd $dir
	 if [ "$nu" == "0" ]; then
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../3_eq.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst" > cmd${nu}
         else
            pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.e -ref ../3_eq.rst
            echo "pmemd.cuda -O -i ../duck_325K.in -o duck_${nu}.o -p ../system_complex.prmtop -c ../md${nu}.rst -r duck_${nu}.rst -x duck_${nu}.nc -e duck_${nu}.nc -e duck_${nu}.e -e ../3_eq.rst" > cmd${nu}
         fi 
         cd ..
      fi
      nu=$((nu+1))
   done

}

# Function to check if WQB is lower than 0.1 using getWqbValues.py
# getWqbValues.py is a script from Maciej modified.
# We use this one instead of the R version, as R is not available in the IQTC

check_WQB(){
   wqb_limit=$1
   lowest_wqb=$(python getWqbValues.py)
   echo $lowest_wqb > wqb.log
   are_we_above_the_limit=$(echo "$lowest_wqb < $wqb_limit" | bc )
   if [ "$are_we_above_the_limit" == "1" ]; then
      echo "Wqb lower than ${wqb_limit}, stoping DUck run"
      cp -r ./* $LIG_TARGET/
      exit 
   fi
}



#### Modules ####
#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB
module load amber/20


#### PARAMS ####
replicas=5
min_wqb=7

#### Runing Duck ####
# Minimization&Equilibration
pmemd.cuda -O -i 1_min.in -o min.out -p system_complex.prmtop -c system_complex.inpcrd -r min.rst -ref system_complex.inpcrd
pmemd.cuda -O -i 2_heat150.in -o 2_heat150.out -p system_complex.prmtop -c min.rst -r  2_heat150.rst -x 2_heat150.nc -ref system_complex.inpcrd
pmemd.cuda -O -i 2_heat200.in -o 2_heat200.out -p system_complex.prmtop -c 2_heat150.rst -r 2_heat200.rst -x 2_heat200.nc -ref 2_heat150.rst
pmemd.cuda -O -i 2_heat250.in -o 2_heat250.out -p system_complex.prmtop -c 2_heat200.rst -r 2_heat250.rst -x 2_heat250.nc -ref 2_heat200.rst
pmemd.cuda -O -i 2_heat300.in -o 2_heat300.out -p system_complex.prmtop -c 2_heat250.rst -r 2_heat300.rst -x 2_heat300.nc -ref 2_heat_250.rst
pmemd.cuda -O -i 3_eq.in -o 3_eq.out -p system_complex.prmtop -c 2_heat300.rst -r 3_eq.rst -x 3_eq.nc -ref 2_heat300.rst -e 3_eq.ene

#Launch DUck 0
prepare_duck_and_launch 0 0 300K
prepare_duck_and_launch 0 0 325K

#Check if WQB is not lower than the min
check_WQB $min_wqb

#For each replica wanted do: MD, prepare SMD & launch SMD
for ((i=1;i<=$replicas;++i)); do
   if [ "$i" == "1" ]; then
      pmemd.cuda -O -i md.in -o md1.out -p system_complex.prmtop -c 3_eq.rst -r md1.rst -x md1.nc -ref 3_eq.rst
   else
      pmemd.cuda -O -i md.in -o md${i}.out -p system_complex.prmtop -c md$((i-1)).rst -r md${i}.rst -x md${i}.nc -ref 3_eq.rst
   fi

   prepare_duck_and_launch $i $i 300K
   prepare_duck_and_launch $i $i 325K
   check_WQB $min_wqb

done

exit
    """}
    print('Writing queue files')
    if template not in queue_strings:
        print('Warning wrong queue template. Only {} accepted'.format(list(queue_strings.keys())))
    write_string_to_file('merged_duck_template.q', queue_strings[template])
