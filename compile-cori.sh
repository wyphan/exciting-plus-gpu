#!/bin/bash

about() {
  echo "Exciting-Plus compile script for Cori and Cori-GPU (NERSC)"
  echo "Last edited: Jun 22, 2020 (WYP)"
}

# Check whether script is executed from Cori login node or Cori-GPU node
curnode=`hostname`
if [ [ "x${curnode::4}" != "xcori" ] || [ "x${curnode::4}" != "xcgpu" ] ]; then
  echo "ERROR: script not executed on Cori"
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
GCCVER="GCC 8.3.0"
compilers() {
  echo "On Cori, Exciting-Plus has been tested with the following compilers:"
  echo "  intel ${INTELVER} through Cray compiler wrappers"
  echo "  pgi   ${PGIVER} through Cray compiler wrappers"
#  echo "  gcc   ${GCCVER}"
  return 0
}

helptext() {
  echo "Available tasks:"
  echo
  echo "  help       Show this help text"
  echo
  echo "  elk        Compile Exciting-Plus for Cori-Haswell nodes"
  echo "  knl        Compile Exciting-Plus for Cori-KNL nodes"
  echo "  acc        Compile Exciting-Plus for Cori-GPU nodes"
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
if [ "x$COMPILER" == "x"  ]; then COMPILER=intel; fi
if [ "x$USEMKL"   != "x0" ]; then export USEMKL=1; fi
if [ "x$USEHDF5"  != "x0" ]; then export USEHDF5=1; fi
if [ "x$USEACC"   == "x"  ]; then export USEACC=none; fi

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
      export USEACC=none
      return 0
      ;;

  # Build Exciting-Plus, KNL version
    knl )
      export BUILDELK=1
      export USEACC=knl
      export COMPILER=intel
      ;;

  # Build Exciting-Plus, OpenACC version
    acc )
      export BUILDELK=1
      export USEACC=tesla
      export COMPILER=pgi
      ;;
    
  # Compiler choice
    intel | pgi | gcc )
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

# MKL is loaded through the Intel compiler module
# This function extracts environment variables set through the module
# and saves them to rhea-intelvars.sh
getintelvars() {
  module load intel
  cat > cori-intelvars.sh << __EOF__
export INTEL_PATH=${INTEL_PATH}
export MKLROOT=${MKLROOT}
__EOF__
  chmod +x cori-intelvars.sh
}

case ${COMPILER} in

  intel)
    export COMPILERVER="${INTELVER}"
    ;;

  pgi)
    module swap pgi intel
    getintelvars
    module load intel pgi
    export COMPILERVER="${PGIVER}"
    ;;

  gcc)
    echo "Compiler not tested yet (TODO: write make.inc.cori.gcc.cpu)"
    exit 1
    #module load gcc
    #export COMPILERVER="${GCCVER}"
    ;;

  *)
    echo "Unsupported compiler"
    exit 1
esac

  # Copy the appropriate make.inc and perform necessary module swaps
  # TODO: Write the unavailable make.inc files
  case ${USEACC} in
    none )
      cp make.inc.cori.${COMPILER}.cpu make.inc
      ;;
    knl )
      cp make.inc.cori.intel.knl make.inc
      module swap cray-haswell cray-mic-knl
      ;;
    tesla )
      if [ "x${curnode::4}" != "xcgpu" ]; then
	echo "Warning: per Cori-GPU documentation, cannot cross-compile from Cori login node"
	exit 42
      fi
      cp make.inc.cori.pgi.acc make.inc
      module load cuda
      module load mvapich2
      ;;
    *)
      echo "Error USEACC=$USEACC"
      exit 1
      ;;
  esac

# Build Exciting-Plus CPU-only version
if [ "x${BUILDELK}" == "x1" ]; then

  clear; hline; echo;
  echo "`date` Building elk-cpu with ${COMPILERVER}"
  echo; hline; echo

  # Load Intel MKL
  if [ "x${USEMKL}" == "x1" ]; then
    echo "Using Intel MKL"
    if [ "${COMPILER}" != "intel" ]; then source ./cori-intelvars.sh; fi
  fi

  # Load HDF5
  if [ "x${USEHDF5}" == "x1" ]; then
    module load cray-hdf5
    echo "Using HDF5 (serial)"
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
      module load cray-hdf5 || echo "Using HDF5 (serial)"
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
