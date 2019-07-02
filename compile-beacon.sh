#!/bin/bash

# Load modules
module load mkl
module load hdf5

# Compilers
export MAKE=make
export F90=mpiifort

# Make the binary
cp make.inc.beacon.intel make.inc
make clean
make

# CPU-only version
mv src/elk src/elk-cpu

# Enable CUDA
module load cuda
export CUDA_PATH=${CUDA_DIR}
export NVCC=${CUDA_PATH}/bin/nvcc

# Compile the necessary codes.
$NVCC -c -g -G cublas_fortran.cu
$F90 -cpp -c -g cublas_fortran_iso.f90
$F90 -cpp -g -D_MPI_ -c -I./src/ genmegqblh_cublas.f90

# Move the appropriate files over
cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
cp cublas_fortran_iso.o  cublas_fortran.o *.mod src/addons/expigqr/

# re-Make the binary.
cp make.inc.beacon.intel.gpu make.inc
make

# Copy the hybrid CPU+GPU version
mv src/elk src/elk-gpu

# Keep the two different versions
make install
rm bin/elk
cd bin
ln -s -T ../src/elk-cpu elk-cpu
ln -s -T ../src/elk-gpu elk-gpu
cd ..
