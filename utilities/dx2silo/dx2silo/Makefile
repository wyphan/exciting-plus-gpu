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

all: dx2silo

dx2silo: main.o readbin.o readdxh.o
	$(F90SERIAL) -o $@ $^ $(LIBS)

install: dx2silo
	cp  ../utilities/dx2silo/dx2silo ../../bin/dx2silo

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out core dx2silo
