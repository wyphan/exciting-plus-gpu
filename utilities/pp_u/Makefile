include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c -o $(*F).o $<

all: pp_u4

mod_hdf5.o:
	$(F90) -c $(F90_OPTS) $(HDF5_OPTS) ../../src/addons/mod_hdf5.f90 -o mod_hdf5.o

pp_u4: mod_hdf5.o pp_u4.o 
	$(F90SERIAL) -o $@ mod_hdf5.o pp_u4.o $(HDF5_LIB)

install: pp_u4
	cd ../../bin
	ln -s -T ../utilities/pp_u/pp_u4 ../../bin/pp_u4

clean:
	rm pp_u4 *.o *.mod
