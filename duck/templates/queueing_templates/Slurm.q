#!/bin/bash

#SBATCH --job-name=DUck   
#SBATCH -D .                       
#SBATCH --time=72:00:00            
#SBATCH --output=DUck.q.o         
#SBATCH --error=DUck.q.e          
#SBATCH --ntasks=1                 
#SBATCH --gres=gpu:1               
#SBATCH --cpus-per-task=1     

#### Modules ####
#Load modules (we are missing R in here, but python is installed, so we could use Maciej's scripts to check the WQB
module load amber/20

{functions}

{commands}

exit
