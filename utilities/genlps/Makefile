include ../../make.inc

.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) -c $<

# prevent invocation of m2c
%.o: %.mod

all: genlps

genlps: genlps.o factnm.o ylmroty.o rlmrot.o euler.o
	$(F90SERIAL) -o $@ $^ 

axes3d: mod_axes3d.o axes3d.o
	$(F90SERIAL) -o $@ $^

gramschmidt: mod_axes3d.o gramschmidt.o
	$(F90SERIAL) -o $@ $^ $(LAPACK_LIB)

install:
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./genlps ../../bin/genlps-${COMPILER}
	cp ./axes3d ../../bin/axes3d-${COMPILER}
	cp ./gramschmidt ../../bin/gramschmidt-${COMPILER}

clean:
	-rm genlps axes3d gramschmidt *.o *.mod

lsutil:
	echo "genlps:  genlps axes3d gramschmidt"

utilhelp: | lsutil
	echo "Pre processing utility to generate local rotations"
	echo "- genlps: Generates the 'lps' (local point symmetry) block"
	echo "          See README.txt for more details"
	echo "- axes3d: Orthogonalizes the three Cartesian axes that are input"
	echo "          into genlps. Input and output are both genlps.in"
	echo "- gramschmidt: Performs Gram-Schmidt orthogonalization using the"
	echo "               modified Gram-Schmidt algorithm."
	echo "               Takes in two vectors from stdin and outputs"
	echo "               three vectors to stdout."
