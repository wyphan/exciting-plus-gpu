include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c -o $(*F).o $<

all: pdos bndchr

pdos: addorb.o pdos.o
	$(F90SERIAL) -o $@ addorb.o pdos.o

bndchr: addorb.o bndchr.o
	$(F90SERIAL) -o $@ addorb.o bndchr.o

install: pdos bndchr
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./pdos ../../bin/pdos-${COMPILER}
	cp ./bndchr ../../bin/bndchr-${COMPILER}

clean:
	-rm pdos bndchr *.o

lsutil:
	echo "pp:  pdos bndchr"

utilhelp: | lsutil
	echo "Post processing utilities for ground state"
	echo "- pdos"
	echo "  For use with output from Exciting-Plus task 811 (dosrlm)."
	echo "  Extracts the projected density of states from PDOS.OUT"
	echo "  for Bloch and/or Wannier basis."
	echo "- bndchr"
	echo "  For use with output from Exciting-Plus task 822 (bandrlm)."
	echo "  Extracts the l-resolved band structure from BNDCHR.OUT"
	echo "  for Bloch and/or Wannier basis."
	echo "TODO: write input documentation"
