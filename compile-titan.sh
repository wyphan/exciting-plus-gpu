#Notes from Junqi's initial implementation readme:

# Compile the necessary codes.
nvcc -c -g -G cublas_fortran.cu
ftn -cpp -c -g cublas_fortran_iso.f90
ftn -cpp -g -D_MPI_ -c -I/ccs/home/ssawyer1/exciting-plus-rgvw-mod/src/ genmegqblh_cublas.f90

# Move the appropriate files over
cp genmegqblh_cublas.o src/addons/expigqr/genmegqblh.o
cp cublas_fortran_iso.o  cublas_fortran.o *.mod src/addons/expigqr/

# re-Make the binary.
cp make.inc.titan.intel.gpu make.inc
make

