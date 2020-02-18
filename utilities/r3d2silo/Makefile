include ../../make.inc

# Most modern Fortran compilers can also process Fortran 77 source files
F77SERIAL=${F90SERIAL}
F77_OPTS_SERIAL=${F90_OPTS_SERIAL}

# Include directories 
IDIR=/usr/local/include

# Directory where libsilo is installed, if different from the usual place
SILODIR=/opt/gcclibs

LIBS=-L${SILODIR}/lib -lsilo -lm

.SUFFIXES: .o .f90 .f
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -I$(IDIR) -I$(SILODIR)/include -c $<
.f.o:
	$(F77SERIAL) $(F77_OPTS_SERIAL) -I$(IDIR) -I$(SILODIR)/include -c $<

all: r3d2silo

r3d2silo: main.o readr3d.o
	$(F90) -o $@ $^ $(LIBS)

install: r3d2silo
	ln -s -T ../utilities/r3d2silo/r3d2silo ../../bin/r3d2silo

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out core r3d2silo
