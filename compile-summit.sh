#!/bin/bash

about() {
  echo "Exciting-Plus compile script for Summit (ORNL)"
  echo "Last edited: May 6, 2020 (WYP)"
}

# Check whether script is executed from Summit login node
curnode=`hostname --fqdn | awk 'BEGIN { FS ="." } ; { print $2 }'`
if [ "x$curnode" != "xsummit" ]; then
  echo "ERROR: script not executed on Summit"
  exit 42 # Don't panic
fi

usage() { echo "Usage: $0 [compiler] [task]"; }

tasklist() {
  echo "Available tasks:"
  echo "  help,"
  echo "  elk, tau,"
  echo "  pp, pp_u, pp_u4, spacegroup, utils"
  return 0
} 

# TODO: accomodate multiple compiler versions and extract them automatically
IBMVER="IBM XL 16.1.1-5 (default compiler)"
PGIVER="PGI 19.9" # ", 20.1"
compilers() {
  echo "On Summit, Exciting-Plus has been tested with the following compilers:"
  echo "  ibm   ${IBMVER} (default compiler)"
  echo "  pgi   ${PGIVER}"
#  echo "  gcc   GCC 6.4.0"
#  echo "  llvm  Clang/Flang 8.0.0+git"
  return 0
}

helptext() {
  echo "Available tasks:"
  echo
  echo "  help       Show this help text"
  echo
  echo "  elk        Compile Exciting-Plus"
  echo "  tau        Compile Exciting-Plus with TAU 2.29.1 + chosen compiler"
  echo
  echo "  pp         Compile 'bndchr' and 'pdos' utilities"
  echo "  pp_u       Compile 'pp_u4' utility"
  echo "  spacegroup Compile 'spacegroup' utility"
#  echo "  eos        Compile 'eos' utility"
#  echo "  plot3d     Compile 'sicvlm' and 'plot_wan_dens' utilities"
  echo "  utils      Compile all of the above utilities"
  echo
  echo "If no compiler choice is given, then the default compiler will be used."
  echo "By default, these are turned on: MPI, OpenMP, ESSL, HDF5"
  echo "Modify the appropriate 'make.inc' files for finer-grained control"
  echo "For now, please don't supply two compilers or two tasks"
  echo "TODO: improve compile script"
}

# Default choices (can be overriden through environment variables)
if [ "x$MAKE"     == "x"  ]; then MAKE=make; fi
if [ "x$COMPILER" == "x"  ]; then COMPILER=ibm; fi
if [ "x$USEESSL"  != "x0" ]; then export USEESSL=1; fi
if [ "x$USEHDF5"  != "x0" ]; then export USEHDF5=1; fi

# Default choices
export BUILDELK=1
export BUILDUTILS=0
export USETAU=0

# Debugging shortcuts
export EXCDIR=`pwd`
export WANN="${EXCDIR}/src/addons/wann/"
export EXPI="${EXCDIR}/src/addons/expigqr/"

# Function to print '=' 80 times, adapted from this link
# https://stackoverflow.com/questions/5349718/how-can-i-repeat-a-character-in-bash
hline() { printf '=%.0s' {1..80}; printf '\n'; }

# Array to store list of utilities to compile
declare -a UTILS
declare -a DEFUTILS

# Function to parse a single task
# TODO: account for multiple compilers and/or utils
parsetask() {
  case "$1" in

  # Show full help text
    help | -h | --help )
      about; echo; usage;
      echo; hline; echo;
      compilers;
      echo; hline; echo;
      helptext; echo;
      return 0
      ;;

  # Build Exciting-Plus, CPU-only version
    elk )
      export BUILDELK=1
      return 0
      ;;

  # Build instrumented Exciting-Plus for profiling with TAU
    tau )
      export USETAU=1
      export COMPILER="tau-${COMPILER}"
      return 0
      ;;

  # Compiler choice
    ibm | pgi | gcc | llvm )
      export BUILDELK=1
      export COMPILER="$1"
      return 0
      ;;

  # Utilities choice
    pp | pp_u | spacegroup | eos | plot3d )
      export BUILDELK=0
      export BUILDUTILS=1
      UTILS+=("$1")
      return 0
      ;;

  # Alias for pp_u4 -> pp_u
    pp_u4 )
      export BUILDELK=0
      export BUILDUTILS=1
      UTILS+=("pp_u")
      return 0
      ;;

  # Default set of utilities
    utils )
      export BUILDELK=0
      export BUILDUTILS=1
      UTILS=("pp" "pp_u" "spacegroup")
      return 0
      ;;

  # Invalid input
    *)
      echo "Unknown task $1"; return 1 ;;

  esac
}

# Parse arguments
# TODO: rewrite for any number of arguments
if [ "x$2" != "x" ]; then
  # argc = 2
  parsetask "$1"; if [ "x$?" != "x0" ]; then tasklist; exit 1; fi
  parsetask "$2"; if [ "x$?" != "x0" ]; then tasklist; exit 2; fi
elif [ "x$1" != "x" ]; then
  # argc = 1
  parsetask "$1"; if [ "x$?" != "x0" ]; then tasklist; exit 1; fi
fi

# ESSL depends on libxlf90_r
# This function extracts IBM XL compiler paths and saves them to xlvars.sh
getxlvars() {
  module load xl
  cat > summit-xlvars.sh << __EOF__
export OLCF_XL_ROOT=${OLCF_XL_ROOT}
export OLCF_XLF_ROOT=${OLCF_XLF_ROOT}
export OLCF_XLC_ROOT=${OLCF_XLC_ROOT}
export OLCF_XLMASS_ROOT=${OLCF_XLMASS_ROOT}
export OLCF_XLSMP_ROOT=${OLCF_XLSMP_ROOT}
__EOF__
  chmod +x summit-xlvars.sh
}

# Note: no longer needed after Mar 10 2020 stack upgrade
# PGI's OpenMP implementation relies on GCC's libatomic
# This function extracts the GCC compiler path and saves it to xlvars.sh
#getgccvars() {
#  module load gcc
#  echo "export OLCF_GCC_ROOT=${OLCF_GCC_ROOT}" > summit-gccvars.sh
#  chmod +x summit-gccvars.sh
#}

# TODO: decouple tau options from compiler
case ${COMPILER} in

  ibm)
    module load xl
    export COMPILERVER="${IBMVER}"
    ;;

  pgi)
    getxlvars
    #getgccvars
    module load pgi
    export COMPILERVER="${PGIVER}"
    #source ./summit-gccvars.sh
    ;;

  gcc)
    echo "Compiler not tested yet (TODO: rewrite make.inc.summit.gcc.cpu)"
    exit 1
    #getxlvars
    #module load gcc
    #export COMPILERVER="${GCCVER}"
    ;;

  llvm)
    echo "Compiler not tested yet (TODO: write make.inc.summit.llvm.cpu)"
    exit 1
    #getxlvars
    #module load llvm
    #export COMPILERVER="${LLVMVER}"
    ;;

  tau-ibm)
    # TODO: Resolve ticket #419691
    module load xl/16.1.1-5
    export COMPILERVER="${IBMVER}"
    export TAUVER="2.29.1"
    module load tau/${TAUVER}
    export TAU_MAKEFILE=${TAU_DIR}/lib/Makefile.tau-xl_16.1.1-5-papi-mpi
    module load papi
    ;;

  tau-pgi)
    # TODO: Resolve ticket #419691 and test PGI 20.1
    getxlvars # for ESSL
    #getgccvars
    module load pgi/19.9
    export COMPILERVER="${PGIVER}"
    export TAUVER="2.29.1"
    module load tau/${TAUVER}
    export TAU_MAKEFILE=${TAU_DIR}/lib/Makefile.tau-pgi_19.9-papi-mpi-pgi
    module load papi
    #source ./summit-gccvars.sh
    ;;

  tau-gcc)
    echo "Compiler not yet tested (TODO: write make.inc.summit.tau-gcc.cpu)"
    exit 1
    ;;

  tau-llvm)
    echo "Compiler not yet tested (TODO: write make.inc.summit.tau-llvm.cpu)"
    exit 1
    ;;

  *)
    echo "Unsupported compiler"
    exit 1
esac

# Build Exciting-Plus CPU-only version
if [ "x${BUILDELK}" == "x1" ]; then

  clear; hline; echo;

  if [ "x${USETAU}" == "x1" ]; then
    echo "`date` Building elk-cpu with TAU ${TAUVER} and ${COMPILERVER}"
    echo "Using TAU_MAKEFILE:"
    echo "  ${TAU_MAKEFILE##*/}"
  else
    echo "`date` Building elk-cpu with ${COMPILERVER}"
  fi

  echo; hline; echo

  # Load IBM ESSL
  if [ "x${USEESSL}" == "x1" ]; then
    module load essl
    echo "Using IBM ESSL"
    if [ "${COMPILER:(-3)}" != "ibm" ]; then source ./summit-xlvars.sh; fi
  fi

  # Load HDF5
  if [ "x${USEHDF5}" == "x1" ]; then
    module load hdf5
    echo "Using HDF5"
  fi

  # Copy the appropriate make.inc
  # TODO: Write the unavailable make.inc files
  cp make.inc.summit.${COMPILER}.cpu make.inc

  # Extract link line from make.inc
  if [ "x${USETAU}" == "x1" ]; then
    make lsvars
    source ./libs.sh
    # Apply options
    export TAU_OPTIONS="-optCompInst -optRevert -optTrackIO -optLinking=\"${LIBS}\""
  fi

  # Clean build directory
  ${MAKE} clean
  #rm *.o *.mod

  # Build elk-cpu and check error code
  ${MAKE}
  RETVAL=$?
  if [ $RETVAL != 0 ]; then
    # Build failed
    echo; hline; echo;
    echo "`date` Build failed for elk-cpu with error code ${RETVAL}"
    echo; hline; echo;
    exit $RETVAL
  else
    # Build completed, install elk-cpu
    ${MAKE} install-elk
    echo; hline; echo;
    echo "`date` Success! Built and installed elk-cpu to ./bin/"
    echo; hline; echo;
  fi
  
fi # BUILDELK

# Build and install the utilities
if [ "x${BUILDUTILS}" == "x1" ]; then
  for util in ${UTILS[@]}; do

    echo; hline; echo;
    echo "`date` Building ${util}"
    echo; hline; echo;

    dir="utilities/${util}"

    # pp_u needs HDF5
    if [ "${util}" == "pp_u" ] && [ "x${USEHDF5}" != "x1" ]; then
      echo "pp_u requires HDF5"
      exit 1
    else
      module load hdf5 || echo "Using HDF5"
    fi
    ${MAKE} -C "${dir}" clean;

    # Build the utility and catch error code
    ${MAKE} -C ${dir}
    RETVAL=$?
    if [ $RETVAL != 0 ]; then
      # Build failed
      echo; hline; echo;
      echo "`date` Build failed for ${util} with error code ${RETVAL}"
      echo; hline; echo;
      exit $RETVAL
    else
      # Build completed, install this utility
      ${MAKE} -C ${dir} install
      files="$(${MAKE} -s -C ${dir} lsutil)"
      echo; hline; echo;
      echo "`date` Installed ${files#${util}:  } to ./bin/"
      echo; hline; echo;
    fi

  done
fi

# Clean up variables
unset COMPILER
unset COMPILERVER
unset BUILDELK
unset BUILDUTILS
unset UTILS
unset USEESSL
unset USEHDF5
unset USETAU
unset TAUVER

echo; hline; echo;
echo " Done! "
echo; hline; echo;

exit 0
