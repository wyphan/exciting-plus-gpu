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
SRC = modmain.f90 eos.f90 eveos.f90 pveos.f90 readinput.f90 getedata.f90 \
      fitdata.f90 output.f90 minf_nm.f90 fmin_nm.f90

OBJ = $(SRC:.f90=.o)
EXE = eos

eos:	$(OBJ)
	$(F90SERIAL) $(F90_OPTS_SERIAL) -o $(EXE) $(OBJ)

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno *.OUT gmon.out $(EXE)
