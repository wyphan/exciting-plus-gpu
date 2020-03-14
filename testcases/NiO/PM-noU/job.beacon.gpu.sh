#!/bin/bash
#PBS -S /bin/bash
#PBS -A ACF-UTK0011
#PBS -m bae
#PBS -M wphan@vols.utk.edu
#PBS -j oe
#PBS -v PATH,LD_LIBRARY_PATH
#PBS -N crpa
#PBS -l nodes=1:ppn=16
#PBS -l walltime=00:20:00
#PBS -l feature=beacon_gpu

export EXE="elk-gpu"
export EXEDIR="${HOME}/exciting-plus-gpu/bin"
export JOB="${PBS_JOBNAME}-nk4-pd-ngsh10"
export JOBNUM=`echo "${PBS_JOBID}" | awk -F. '{print $1}'`

export LOGFILE="${PBS_JOBNAME}.log"
export CPUSPECFILE="lscpu.log"
export GPUSPECFILE="CUDAdeviceQuery.log"

export RESDIR="${HOME}/exciting+/NiO-GGA-noU/${JOB}"

# Load modules
module load mkl
module load hdf5
module load cuda/10.0
export CHKCUDA=${CUDA_DIR}/extras/demo_suite/deviceQuery

# Prepare the job
echo "`date` Job $PBS_JOBID launched from `hostname`"
cd $PBS_O_WORKDIR/gpu
echo "Workdir is `pwd`"
cp -L ${EXEDIR}/${EXE} ./
lscpu > $CPUSPECFILE
exec $CHKCUDA > $GPUSPECFILE
tail -n 2 $GPUSPECFILE 

echo "`date` Launching $EXE with $PBS_NNODES ranks, $OMP_NUM_THREADS threads"
mpirun -n $PBS_NNODES -f $PBS_NODEFILE ./$EXE

# Compress output
cp ${PBS_JOBNAME}.o${JOBNUM} $LOGFILE
tar cJf ${JOB}-gpu.tar.xz q *.OUT *.hdf5 $LOGFILE $CPUSPECFILE $GPUSPECFILE
mv ${JOB}-gpu.tar.xz $RESDIR/
rm -f ./$EXE
cd $RESDIR

echo "`date` Done"
