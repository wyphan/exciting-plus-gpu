include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c -o $(*F).o $<

all: pdos bndchr

pdos: addorb.o pdos.o
	$(F90) -o $@ addorb.o pdos.o
bndchr: addorb.o bndchr.o
	$(F90) -o $@ addorb.o bndchr.o

install: pdos bndchr
	cd ../../bin
	ln -s -T ../utilities/pp/pdos ../../bin/pdos
	ln -s -T ../utilities/pp/bndchr ../../bin/bndchr

clean:
	rm pdos bndchr *.o
