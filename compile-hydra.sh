#!/bin/bash

# Compilers
export MAKE=make
export F90=mpif90
export CUDA_PATH=/usr/local/cuda
export NVCC=$CUDA_PATH/bin/nvcc

# Make the binary
cp make.inc.hydra.gcc.mkl.gpu make.inc
make

# Copy the CPU-only version
cp src/elk src/elk-cpu

# Compile the necessary codes.
$NVCC -c -g -G cublas_fortran.cu
$F90 -cpp -c -g cublas_fortran_iso.f90
$F90 -cpp -g -D_MPI_ -c -I./src/ genmegqblh_cublas.f90

# Move the appropriate files over
cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
cp cublas_fortran_iso.o  cublas_fortran.o *.mod src/addons/expigqr/

# re-Make the binary.
cp make.inc.hydra.gcc.mkl.gpu make.inc
make

# Copy the hybrid CPU+GPU version
cp src/elk src/elk-gpu

# Keep the two different versions
make install
rm src/elk
rm bin/elk
cd bin
ln -s -T ../src/elk-cpu elk-cpu
ln -s -T ../src/elk-gpu elk-gpu
cd ..

