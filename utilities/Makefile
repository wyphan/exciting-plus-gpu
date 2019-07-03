include ../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -I../src -c -o $(*F).o $<

all: qdir 

qdir: qdir.o
	$(F90SERIAL) -o $@ qdir.o  $(LIBS) ../src/libelk.a
		
clean:
	rm *.o *.mod
