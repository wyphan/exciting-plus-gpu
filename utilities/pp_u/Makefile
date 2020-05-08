include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c -o $(*F).o $<

all: pp_u4

mod_hdf5.o:
ifdef HDF5_LIB
	$(F90SERIAL) -c $(F90_OPTS_SERIAL) $(HDF5_OPTS) ../../src/addons/mod_hdf5.f90 -o mod_hdf5.o
else
	echo "Please set HDF5_LIB, HDF5_INC, and HDF5_OPTS in make.inc"
endif

pp_u4: mod_hdf5.o pp_u4.o 
ifdef HDF5_LIB
	$(F90SERIAL) -o $@ mod_hdf5.o pp_u4.o $(HDF5_LIB)
else
	echo "Please set HDF5_LIB, HDF5_INC, and HDF5_OPTS in make.inc"
endif

install: pp_u4
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./pp_u4 ../../bin/pp_u4-${COMPILER}

clean:
	-rm pp_u4 *.o *.mod

lsutil:
	echo "pp_u:  pp_u4"

utilhelp: | lsutil
	echo "Post processing utility for CRPA when HDF5 is enabled"
	echo "  For use with output from Exciting-Plus task 801 (crpa)."
	echo "  Obviously, requires HDF5 to build and run."
	echo "TODO: write input documentation"
