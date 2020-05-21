###############################################################################
# Main Makefile for Exciting-Plus
# Last edited: May 20, 2020 (WYP)
###############################################################################

# Note (WYP): this Makefile is heavily skewed towards GNU Make, 
#             especially the '-C' option
ifndef $MAKE
MAKE = make
endif

# Don't forget to supply a 'make.inc' file
# COMPILER and EXE_SFX should be set inside make.inc
include make.inc

# List of default utilities (ignored if UTILS already set somewhere else)
DEFUTILS := spacegroup pp

#------------------------------------------------------------------------------

.PHONY: all
all: clean-elk elk install-elk

clean: clean-elk clean-utils

distclean: clean-elk clean-allutils clean-docs rmdirs

install: install-elk install-utils

# Note: Only runs elk's standard tests/test-{001..016}
# TODO: add Summit-specific tests and streamline job script generation
test:
	cd tests; ./tests.sh

docs: docs-elk docs-spacegroup # docs-crpa

install-docs: docs | mkdir-docs
	cp src/elk.pdf docs/
	cp utilities/spacegroup/spacegroup.pdf docs/

#------------------------------------------------------------------------------

ifndef NVCC

elk:
	cd src; $(MAKE) gensrc; $(MAKE) elk; cd ..

clean-elk:
	cd src; $(MAKE) clean

install-elk: elk | mkdir-bin
	cp src/elk bin/elk-$(COMPILER)-$(EXE_SFX)

else

elk: elk-cpu elk-gpu

elk-cpu:
	cd src; $(MAKE) gensrc; $(MAKE) elk; cd ..
	mv src/elk elk-cpu

elk-gpu: cudareplace cudawrap | elk-cpu
	cd src; $(MAKE) elk; cd ..
	mv src/elk src/elk-gpu

clean-elk: | clean-cuda
	cd src; $(MAKE) clean
	-rm src/elk-cpu
	-rm src/elk-gpu

install-elk: install-elk-cpu install-elk-gpu

install-elk-cpu: elk-cpu | mkdir-bin
	cp src/elk-cpu bin/elk-cpu-$(COMPILER)-$(EXE_SFX)

install-elk-gpu: elk-gpu | mkdir-bin
	cp src/elk-gpu bin/elk-gpu-$(COMPILER)-$(EXE_SFX)

endif

#------------------------------------------------------------------------------

docs-elk: | elk
	$(MAKE) -C src doc

docs-spacegroup: | spacegroup
	$(MAKE) -C utilities/spacegroup doc

#docs-crpa: | mkdir-docs
#	cp src/addons/CRPA-Calculation.pdf docs/

clean-docs: clean-docs-elk clean-docs-spacegroup # clean-docs-crpa

clean-docs-elk:
	$(MAKE) -C src clean-doc

clean-docs-spacegroup:
	$(MAKE) -C utilities/spacegroup clean-doc

#------------------------------------------------------------------------------

ifndef $(UTILS)
UTILS := $(DEFUTILS)
endif

utils: $(UTILS)

lsutils:
	echo "$(UTILS)"

clean-utils:
	$(foreach dir,$(UTILS),$(MAKE) -C utilities/$(dir) clean;)

install-utils: utils | mkdir-bin
	$(foreach dir,$(UTILS),$(MAKE) -C utilities/$(dir) install;)

UTILDIRS := $(subst /,,$(subst utilities,,$(shell ls -d utilities/*/)))

lsutildirs:
	echo "$(UTILDIRS)"

$(%UTILDIRS):
	$(MAKE) -C utilities/$@;

$(foreach util,$(UTILDIRS),lsutil-$(util)):
	$(MAKE) -C utilities/$(subst lsutil-,,$@) utilhelp;

$(foreach util,$(UTILDIRS),clean-$(util)):
	$(MAKE) -C utilities/$(subst clean-,,$@) clean;

$(foreach util,$(UTILDIRS),install-$(util)): $@
	$(MAKE) -C utilities/$(subst install-,,$@) install;

allutils: $(%UTILDIRS)

lsallutils:
	$(foreach dir,$(UTILDIRS),$(MAKE) -C utilities/$(dir) lsutil;)

clean-allutils:
	$(foreach dir,$(UTILDIRS),$(MAKE) -C utilities/$(dir) clean;)

install-allutils: allutils | mkdir-bin
	$(foreach dir,$(UTILDIRS),$(MAKE) -C utilities/$(dir) install;)

#------------------------------------------------------------------------------

mkdirs: mkdir-bin mkdir-docs

mkdir-bin:
	[[ -d bin ]] || mkdir bin

mkdir-docs:
	[[ -d docs ]] || mkdir docs

rmdirs: rmdir-bin rmdir-docs

rmdir-bin:
	rm -rf bin/

rmdir-docs:
	rm -rf docs/

# Function to dump unexpanded variable
showvar = $(info $(1)='$(value $(1))')

showvars:
	@touch .showvars
	$(info ) $(call showvar,MAKE)
	$(foreach var,COMPILER F90 CC CXX,$(call showvar,$(var))) $(info )
	$(foreach var,CPP_OPTS F90_OPTS F90_LINK_OPTS,$(call showvar,$(var)))
	$(info ) $(call showvar,EXE_SFX) $(info )
	$(call showvar,OMP_OPTS) $(info )
ifdef F90SERIAL
	$(call showvar,F90SERIAL) $(call showvar,F90_OPTS_SERIAL) $(info )
else
	$(info F90SERIAL not set) $(info )
endif
	$(call showvar,LIBS) $(info )
	$(call showvar,LAPACK_LIB) $(info )
ifdef HDF5_LIB
	$(foreach var,HDF5_CPP_OPTS HDF5_INC HDF5_LIB,$(call showvar,$(var)))
	$(info )
else
	$(info HDF5 vars not set) $(info )
endif
	@rm .showvars

# Workaround for TAU -OptLinking
lsvars:
	echo "export LIBS=\"$(LIBS)\"" > libs.sh; chmod +x libs.sh
#	echo "export F90_OPTS=\"$(F90_OPTS)\"" > f90opts.sh; chmod +x f90opts.sh

#------------------------------------------------------------------------------

ifdef NVCC

# TODO: There should be a better way than this
cudawrap: cublas_fortran.o cublas_fortran_iso.o
cublas_fortran.o:
	$(NVCC) $(NVCC_OPTS) -c -o $(*F).o cublas_fortran.cu
cublas_fortran_iso.o:
	$(F90) $(F90_OPTS) -c -o $(*F).o cublas_fortran_iso.f90

# TODO: Adapt for more than one replaced subroutine
genmegqblh_cublas.o: cudawrap
	$(F90) $(F90_OPTS) -I./src/ -c -o $(*F).o genmegqblh_cublas.f90

cudareplace: genmegqblh_cublas.o
	cp ./genmegqblh_cublas.o ./src/addons/expigqr/genmegqblh.o
	cp ./cublas_f.mod ./src/addons/expigqr/

endif

clean-cuda:
	-rm ./cublas_f.mod ./cublas_fortran_iso.o ./cublas_fortran.o
	-rm ./genmegqblh_cublas.o 
