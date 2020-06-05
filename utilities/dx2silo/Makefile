include ../../make.inc

# Most modern Fortran compilers can also process Fortran 77 source files
F77SERIAL=$(F90SERIAL)
F77_OPTS_SERIAL=$(F90_OPTS_SERIAL)

.SUFFIXES: .o .f90 .f
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c $<

# Only the main file depends on silo.inc,
# and it's also the only one with this suffix
.f.o:
ifdef SILO_INC
	$(F77SERIAL) $(F77_OPTS_SERIAL) $(SILO_INC) -c $<
else
	echo "Please set SILO_LIB and SILO_INC in make.inc"
endif

all: dx2silo

# Make sure to supply SILO_INC and SILO_LIB through the make.inc file
dx2silo: main.o readbin.o readdxh.o
ifdef SILO_LIB
	$(F90SERIAL) $(F90_OPTS_SERIAL) -o $@ $^ $(SILO_LIB)
else
	echo "Please set SILO_LIB and SILO_INC in make.inc"
endif

install: dx2silo
	cp ./dx2silo ../../bin/dx2silo-${COMPILER}

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out core dx2silo

lsutil:
	echo "dx2silo:  dx2silo"

utilhelp: | lsutil
	echo "Post processing utility for converting OpenDX to Silo format"
	echo "  For use with output from Exciting-Plus task 863 (wann_plot_3d)."
	echo "  Obviously, requires Silo library to build and run."
