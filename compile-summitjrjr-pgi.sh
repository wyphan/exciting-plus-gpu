#!/bin/bash

# Compilers
export MAKE=make
export F90=${PGI}/linux86-64/2019/mpi/openmpi-3.1.3/bin/mpif90
export CUDA_PATH=/usr/local/cuda
export NVCC=$CUDA_PATH/bin/nvcc

# Clean up
make clean
rm *.o *.mod
rm bin/elk-cpu bin/elk-gpu

# Make the binary
cp make.inc.local.pgi.cpu make.inc
make

# CPU-only version
mv src/elk src/elk-cpu

# Compile the necessary codes.
$NVCC -c -g -G cublas_fortran.cu
$F90 -Mpreprocess -c -g cublas_fortran_iso.f90
$F90 -Mpreprocess -g -D_MPI_ -c -I./src/ genmegqblh_cublas.f90
#$F90 -Mpreprocess -g -D_MPI_ -c -I./src/ genvscrn_cublas.f90

# Move the appropriate files over
cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
#cp genvscrn_cublas.o   src/addons/genvscrn.o
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/expigqr/
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/addons/
cp cublas_fortran_iso.o cublas_fortran.o cublas_f.mod src/

# Remove main and mod_mpi_grid so they will be recompiled
rm src/main.o
#rm src/mod_mpi_grid.mod # Keep old module for pstop()
rm src/addons/mod_mpi_grid.o

# re-Make the binary.
cp make.inc.local.pgi.gpu make.inc
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

