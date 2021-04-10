#!/bin/bash

about() {
  echo "Exciting-Plus compile script for Andes (ORNL)"
  echo "Last edited: Apr 10, 2021 (WYP)"
}

# Check whether script is executed from Rhea login node
curnode=`hostname --fqdn | awk 'BEGIN { FS ="-" } ; { print $1 }'`
if [ "x$curnode" != "xandes" ]; then
  echo "ERROR: script not executed on Andes"
  exit 42 # Don't panic
fi

usage() { echo "Usage: $0 [compiler] [task]"; }

tasklist() {
  echo "Available tasks:"
  echo "  help,"
  echo "  elk,"
  echo "  pp, pp_u, pp_u4, spacegroup, dx2silo, utils"
  return 0
} 

# TODO: accomodate multiple compiler versions and extract them automatically
INTELVER="Intel 19.0"
PGIVER="PGI 19.10"
GCCVER="GCC 9.3.0"
LLVMVER="AOCC 2.2.0 (based on LLVM 10.0)"
compilers() {
  echo "On Andes, Exciting-Plus has been tested with the following compilers:"
  echo "  llvm ${LLVMVER} (default compiler)"
#  echo "  intel ${INTELVER}"
#  echo "  pgi   ${PGIVER}"
#  echo "  gcc   ${GCCVER}"
  return 0
}

helptext() {
  echo "Available tasks:"
  echo
  echo "  help       Show this help text"
  echo
  echo "  elk        Compile Exciting-Plus"
  echo
  echo "  pp         Compile 'bndchr' and 'pdos' utilities"
  echo "  pp_u       Compile 'pp_u4' utility"
  echo "  spacegroup Compile 'spacegroup' utility"
#  echo "  eos        Compile 'eos' utility"
#  echo "  plot3d     Compile 'sicvlm' and 'plot_wan_dens' utilities"
  echo "  dx2silo    Compile 'dx2silo' utility"
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
if [ "x$COMPILER" == "x"  ]; then COMPILER=llvm; fi
if [ "x$USEOBLAS" != "x0" ]; then export USEOBLAS=1; fi
#if [ "x$USEHDF5"  != "x0" ]; then export USEHDF5=1; fi
#if [ "x$USESILO"  != "x0" ]; then export USESILO=1; fi
if [ "x$USEFFTW"  != "x0" ]; then export USEFFTW=1; fi

# Default choices
export BUILDELK=1
export BUILDUTILS=1

# Debugging shortcuts
export EXCDIR=`pwd`

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

  # Compiler choice
    llvm | pgi | gcc | intel )
      export BUILDELK=1
      export COMPILER="$1"
      return 0
      ;;

  # Utilities choice
    pp | spacegroup | eos | plot3d )
      export BUILDELK=0
      export BUILDUTILS=1
      UTILS+=("$1")
      return 0
      ;;

  # Alias for pp_u4 -> pp_u
    pp_u | pp_u4 )
      export BUILDELK=0
      export BUILDUTILS=1
      export USEHDF5=1
      UTILS+=("pp_u")
      return 0
      ;;

    dx2silo )
      export BUILDELK=0
      export BUILDUTILS=1
      export USESILO=1
      UTILS+=("dx2silo")
      return 0
      ;;

  # Default set of utilities
    utils )
      export BUILDELK=0
      export BUILDUTILS=1
      UTILS=("pp" "pp_u" "spacegroup" "dx2silo")
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

case ${COMPILER} in

  llvm)
    module load aocc
    export COMPILERVER="${LLVMVER}"
    #export USEAOCL=1
    ;;

  intel)
    echo "Compiler not tested yet (TODO: write make.inc.andes.intel.cpu)"
    exit 1
    #module load intel
    #export COMPILERVER="${INTELVER}"
    #export USEMKL=1
    ;;

  pgi)
    echo "Compiler not tested yet (TODO: write make.inc.andes.pgi.cpu)"
    exit 1
    #module load intel
    #getintelvars
    #module load pgi
    #export COMPILERVER="${PGIVER}"
    ;;

  gcc)
    module load gcc
    export COMPILERVER="${GCCVER}"
    ;;

  *)
    echo "Unsupported compiler"
    exit 1
esac

# Copy the appropriate make.inc
# TODO: Write the unavailable make.inc files
cp make.inc.andes.${COMPILER}.cpu make.inc

# Build Exciting-Plus CPU-only version
if [ "x${BUILDELK}" == "x1" ]; then

  clear; hline; echo;
  echo "`date` Building elk-cpu with ${COMPILERVER}"
  echo; hline; echo

  # Load OpenBLAS
  if [ "x${USEOBLAS}" == "x1" ]; then
    module load openblas
    echo "Using OpenBLAS"
  fi

  # Load reference BLAS and LAPACK
  if [ "x${USEREFBLAS}" == "x1" ]; then
    module load reflapack
    echo "Using reference BLAS and LAPACK"
  fi

  # Load AMD AOCL
  if [ "x${USEAOCL}" == "x1" ]; then
    module load aocl
    echo "Using AMD AOCL"
  fi

  # Load Intel MKL
  if [ "x${USEMKL}" == "x1" ]; then
    echo "Using Intel MKL"
    if [ "${COMPILER}" != "intel" ]; then source ./rhea-intelvars.sh; fi
  fi

  # Load FFTW 3
  if [ "x${USEFFTW}" == "x1" ]; then
    module load fftw
    echo "Using FFTW 3"
  fi

  # Load HDF5
  if [ "x${USEHDF5}" == "x1" ]; then
    module load hdf5
    echo "Using HDF5"
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

    # dx2silo needs Silo (duh)
    if [ "${util}" == "dx2silo" ] && [ "x${USESILO}" != "x1" ]; then
      echo "dx2silo requires Silo"
      exit 1
    else
      module load silo || echo "Using Silo"
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
unset USEMKL
unset USEHDF5
unset USESILO

echo; hline; echo;
echo " Done! "
echo; hline; echo;

exit 0
