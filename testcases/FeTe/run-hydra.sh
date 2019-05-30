#!/bin/bash

# Environment variables
OMPI_ROOT=/usr/lib64/openmpi
CUDA_PATH=/usr/local/cuda
MKL_PATH=/home/wphan/intel/mkl
export PATH="${OMPI_ROOT}/bin:$PATH"
export LD_LIBRARY_PATH="${OMPI_ROOT}/lib:${MKL_PATH}/lib/intel64:${CUDA_PATH}/lib64:$LD_LIBRARY_PATH"
WDIR=$(pwd)

# Command to execute
EXE="/home/wphan/exciting-plus-gpu-dbg/bin/elk-cpu"
#EXE="/home/wphan/exciting-plus-gpu-dbg/bin/elk-gpu"

# Processes per node
PPN=8

# Set up hostfile
cat <<__EOF__ > my_hosts
hydra2 slots=8
__EOF__

# Run
mpirun -x LD_LIBRARY_PATH=$LD_LIBRARY_PATH --prefix $OMPI_ROOT -wdir $WDIR \
       -hostfile my_hosts -npernode $PPN $EXE
