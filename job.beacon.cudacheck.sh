#!/bin/bash
#PBS -S /bin/bash
#PBS -A ACF-UTK0011
#PBS -m bae
#PBS -M wphan@vols.utk.edu
#PBS -j oe
#PBS -v PATH,LD_LIBRARY_PATH
#PBS -N checkcuda
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:01:00
#PBS -l feature=beacon_gpu

export EXE=$HOME/exciting-plus-gpu/cudacheck-beacon.sh

# Load modules
module load cuda/10.0

# Start the job
echo "`date` Job $PBS_JOBID launched from `hostname`"
cd $SCRATCHDIR
exec $EXE
echo "`date` Done"
