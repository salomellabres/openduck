#!/bin/bash
#$ -N DUck_queue          # The name of the job, can be whatever makes sense to you
#$ -S /bin/bash          # Force sh if not Sun Grid Engine default shell
#$ -cwd                 # The batchsystem should use the current directory as working directory.
#$ -q fartorgpu.q            # Queue name where the job should be placed into.
#$ -o DUck.q.o             # Redirect output stream to this file.
#$ -e DUck.q.e             # Redirect error stream to this file.
#$ -l h_rt=15:00:00 # Time limit
#$ -pe gpu 1
#$ -m e

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

{functions}

{commands}

#### Coping the files back to local ####
cp -r ./* $LIG_TARGET/
cd $LIG_TARGET

exit