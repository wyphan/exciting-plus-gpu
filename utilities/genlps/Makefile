include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) -c $<

all: genlps

genlps: main.o factnm.o ylmroty.o rlmrot.o euler.o
	$(F90SERIAL) -o $@ $^ 

install:
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./genlps ../../bin/genlps-${COMPILER}

clean:
	-rm genlps *.o

lsutil:
	echo "genlps:  genlps"

utilhelp: | lsutil
	echo "Pre processing utility to generate local rotations"
	echo "  Generates the 'lps' (local point symmetry) block"
	echo "  See README.txt for more details"
