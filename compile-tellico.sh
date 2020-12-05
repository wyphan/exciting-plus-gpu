#!/bin/bash

about() {
  echo "Exciting-Plus compile script for Tellico (ICL)"
  echo "Last edited: Oct 18, 2020 (WYP)"
}

# Check whether script is executed from Summit login node
curnode=`hostname --fqdn | awk 'BEGIN { FS ="-" } ; { print $1 }'`
if [ "x$curnode" != "xtellico" ]; then
  echo "ERROR: script not executed on Tellico"
  exit 42 # Don't panic
fi

usage() { echo "Usage: $0 [compiler] [task]"; }

tasklist() {
  echo "Available tasks:"
  echo "  help,"
  echo "  elk, tau, acc, tau-acc, scorep"
  echo "  pp, pp_u, pp_u4, spacegroup, utils"
  return 0
} 

# TODO: accomodate multiple compiler versions and extract them automatically
#IBMVER="IBM XL 16.1.1-5"
NVVER="NVIDIA HPC SDK 20.9"
#GCCVER="GCC 10.2.0"

compilers() {
  echo "On Tellico, Exciting-Plus has been tested with the following compilers:"
  echo "  nv   ${NVVER}"
#  echo "  ibm   ${IBMVER}"
#  echo "  gcc   ${GCCVER}"
  return 0
}

helptext() {
  echo "Available tasks:"
  echo
  echo "  help       Show this help text"
  echo
  echo "  elk        Compile Exciting-Plus"
  echo "  acc        Compile Exciting-Plus with OpenACC (requires PGI compiler)"
  echo
  echo "  pp         Compile 'bndchr' and 'pdos' utilities"
  echo "  pp_u       Compile 'pp_u4' utility"
  echo "  spacegroup Compile 'spacegroup' utility"
#  echo "  eos        Compile 'eos' utility"
#  echo "  plot3d     Compile 'sicvlm' and 'plot_wan_dens' utilities"
#  echo "  dx2silo    Compile 'dx2silo' utility"
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
if [ "x$COMPILER" == "x"  ]; then COMPILER=nv; fi
if [ "x$USEESSL"  != "x0" ]; then export USEESSL=1; fi
if [ "x$USEHDF5"  != "x0" ]; then export USEHDF5=1; fi
if [ "x$USEFFTW"  != "x0" ]; then export USEFFTW=1; fi
if [ "x$USEACC"   == "x"  ]; then export USEACC=none; fi

# Default choices
export BUILDELK=1
export BUILDUTILS=0
export USECUDA=1
export USEREFBLAS=0
export MAKEJOBS=8 # Reasonable enough5A

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
      export USEACC=none
      return 0
      ;;

  # Build Exciting-Plus, OpenACC version
    acc )
      export BUILDELK=1
      export USEACC=volta
      if [ "x$COMPILER" == "x" ]; then export COMPILER=nv; fi
      return 0
      ;;

  # Compiler choice
    ibm | nv | gcc )
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

    #dx2silo )
      #export BUILDELK=0
      #export BUILDUTILS=1
      #export USESILO=1
      #UTILS+=("dx2silo")
      #return 0
      #;;

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

case ${COMPILER} in

  nv)
    module load nvhpc/20.9
    export COMPILERVER="${PGIVER}"
    ;;

  ibm)
    echo "Compiler not tested yet (TODO: write make.inc.ascent.ibm.cpu)"
    exit 1
    #module load xl
    #export COMPILERVER="${IBMVER}"
    ;;

  gcc)
    echo "Compiler not tested yet (TODO: write make.inc.ascent.gcc.cpu)"
    exit 1
    #module load gcc
    #export COMPILERVER="${GCCVER}"
    ;;

  *)
    echo "Unsupported compiler"
    exit 1
    ;;

esac

# Copy the appropriate make.inc
# TODO: Write the unavailable make.inc files
case ${USEACC} in
  none)
    cp -L make.inc.tellico.${COMPILER}.cpu make.inc
    ;;
  volta)
    cp -L make.inc.tellico.nv.acc make.inc
    module load cuda
    if [ "x$USEESSL" == "x1" ]; then 
      # IBM ESSL isn't complete, add reference LAPACK
      module load lapack
    fi
    ;;
  *)
    echo "Error USEACC=$USEACC"
    exit 1
    ;;
esac

# Build Exciting-Plus
if [ "x${BUILDELK}" == "x1" ]; then

  clear; hline; echo;

  if [ "x${USETAU}" == "x1" ]; then
    echo "`date` Building elk-cpu with TAU ${TAUVER} and ${COMPILERVER}"
    echo "Using TAU_MAKEFILE:"
    echo "  ${TAU_MAKEFILE##*/}"
  elif [ "x${USESCOREP}" == "x1" ]; then
    echo "Note: Score-P is available only for select compilers on Summit; make sure ${COMPILERVER} is included."
    echo "`date` Building elk-cpu with Score-P 6.0 and ${COMPILERVER}"
    module load scorep/6.0
  else
    echo "`date` Building elk-cpu with ${COMPILERVER}"
  fi

  echo; hline; echo

  # Load IBM ESSL
  if [ "x${USEESSL}" == "x1" ]; then
    module load essl
    echo "Using IBM ESSL"
    #if [ "${COMPILER:(-3)}" != "ibm" ]; then source ./summit-xlvars.sh; fi
  fi

  # Load reference BLAS and LAPACK
  if [ "x${USEREFBLAS}" == "x1" ]; then
    module load lapack
    echo "Using reference BLAS and LAPACK"
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

  # Load CUDA
  if [ "x${USECUDA}" == "x1" ]; then
    module load cuda
    echo "Using CUDA (for nvTX)"
  fi

  # Clean build directory
  ${MAKE} clean
  #rm *.o *.mod

  ${MAKE} -j ${MAKEJOBS}
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
unset USEFFTW
unset USECUDA
unset USEACC

echo; hline; echo;
echo " Done! "
echo; hline; echo;

exit 0
