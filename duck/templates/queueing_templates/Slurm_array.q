#!/bin/bash

#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1  
#SBATCH --array=1-{array_limit}   

#### Modules ####
module load amber

{functions}

{mmpbsa}

cd LIG_target_$SLURM_ARRAY_TASK_ID

{commands}

exit
