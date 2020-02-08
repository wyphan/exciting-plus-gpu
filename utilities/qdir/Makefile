include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -I../src -c -o $(*F).o $<

all: qdir

qdir: qdir.o
	$(F90SERIAL) -o $@ qdir.o ../../src/libelk.a

install: qdir
	ln -s -T ../utilities/qdir/qdir ../../bin/qdir

clean:
	rm *.o *.mod
