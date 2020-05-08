include ../../make.inc

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90SERIAL) $(F90_OPTS_SERIAL) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------

SRC_modules = modmain.f90

SRC_main = main.f90

SRC_routines = readinput.f90 gencrystal.f90 sgsymb.f90 seitzgen.f90 \
 gengroup.f90 seitzeq.f90 seitzmul.f90 writegeom.f90 geomplot.f90

SRC_lib = findprim.f90 r3frac.f90 r3taxi.f90 r3mv.f90 r3cross.f90 r3minv.f90 \
 r3mm.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_lib)

OBJ = $(SRC:.f90=.o)
EXE = spacegroup

spacegroup: $(OBJ)
	$(F90SERIAL) $(F90_OPTS_SERIAL) -o $(EXE) $(OBJ)

clean:
	-rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT *.xsf *.ascii gmon.out \
         *.aux *.dvi *.log *.pdf *.tex *.toc $(EXE)

# Maintain compatibility with Elk
copy: install

install: spacegroup
	[[ -d ../../bin ]] || mkdir -p ../../bin
	cp ./spacegroup ../../bin/spacegroup-${COMPILER}

# TODO: Update documentation to reflect "ANGSTROM" field
doc:
	../../src/protex -s \
           $(SRC_main) $(SRC_modules) $(SRC_routines) $(SRC_lib) \
           > spacegroup.tex
	pdflatex spacegroup; pdflatex spacegroup; pdflatex spacegroup

clean-doc:
	-rm -f $(addprefix spacegroup.,aux log dvi toc pdf tex)

lsutil:
	echo "spacegroup:  spacegroup"

utilhelp: | lsutil
	echo "  Produces GEOMETRY.OUT for inclusion into elk.in"
	echo "  Documentation available at docs/spacegroup.pdf"
	echo "  (requires pdflatex to build)"
