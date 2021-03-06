#!/bin/bash

about() {
  echo "Exciting-Plus compile script for wyp-ThinkPad13"
  echo "Last edited: Dec 23, 2020 (WYP)"
}

usage() { echo "Usage: $0 [compiler] [task]"; }

tasklist() {
  echo "Available tasks:"
  echo "  help,"
  echo "  elk, acc,"
  echo "  pp, pp_u, pp_u4, spacegroup, utils"
  return 0
} 

GCCVER="GCC 9.3.0"
PGIVER="PGI 19.10"
NVVER="NVIDIA HPC SDK 20.11"
INTELVER="Intel 2020 Update 2"
compilers() {
  echo "On ThinkPad13, Exciting-Plus has been tested with the following compilers:"
  echo "  gcc   ${GCCVER} (default compiler)"
  echo "  pgi   ${PGIVER}"
  echo "  nv    ${NVVER}"
  echo "  intel ${INTELVER}"
  return 0
}

helptext() {
  echo "Available tasks:"
  echo
  echo "  help       Show this help text"
  echo
  echo "  elk        Compile Exciting-Plus"
#  echo "  tau        Compile Exciting-Plus with TAU 2.29.1 + chosen compiler"
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
  echo "By default, these are turned on: MPI, OpenMP, OpenBLAS"
  echo "Modify the appropriate 'make.inc' files for finer-grained control"
  echo "For now, please don't supply two compilers or two tasks"
  echo "TODO: improve compile script"
}

# Default choices (can be overriden through environment variables)
if [ "x$MAKE"     == "x"  ]; then MAKE=make; fi
if [ "x$COMPILER" == "x"  ]; then COMPILER=gcc; fi
if [ "x$USEOBLAS" != "x1" ]; then export USEOBLAS=0; export USEMKL=1; fi
if [ "x$USEREFBLAS" != "x1" ]; then export USEREFBLAS=0; export USEMKL=1; fi
if [ "x$USEMKL"   == "x0" ]; then export USEREFBLAS=1; fi
if [ "x$USEHDF5"  != "x0" ]; then export USEHDF5=1; fi
if [ "x$USEFFTW"  != "x0" ]; then export USEFFTW=1; fi
if [ "x$USEACC"   == "x"  ]; then export USEACC=none; fi

# Default choices
export BUILDELK=1
export BUILDUTILS=0
export USETAU=0
export MAKEJOBS=4

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
    help | "-h" | "--help" )
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
    #tau )
      #export USETAU=1
      #export COMPILER="tau-${COMPILER}"
      #return 0
      #;;

  # Build Exciting-Plus, OpenACC version
    acc )
      export BUILDELK=1
      export USEACC=cpu
      ;;

  # Compiler choice
    gcc | pgi | nv | intel )
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

# TODO: decouple tau options from compiler
case ${COMPILER} in

  gcc)
    module load gcc/9.3.0
    module load openmp/4.0.3-gcc9.3.0
    export COMPILERVER="${GCCVER}"

    ;;
    
  pgi)
    module load pgi/19.10-nollvm
    module load pgilibs/19.10
    module load openmpi/3.1.3-pgi19.10+cuda10.1
    export COMPILERVER="${PGIVER}"
    ;;

  nv)
    module load nvhpc/20.11
    module load nvlibs/20.11
    module load openmpi/4.0.5-nv20.11+cuda11.0
    export COMPILERVER="${NVVER}"
    ;;

  intel)
    module load intel/19.1
    module load intellibs/19.1
    module load intelmpi/19.1
    export COMPILERVER="${INTELVER}"
    ;;

  tau-gcc)
    echo "Compiler not yet tested (TODO: write make.inc.basecamp.tau-gcc.cpu)"
    exit 1
    ;;

  tau-pgi)
    echo "Compiler not yet tested (TODO: write make.inc.basecamp.tau-pgi.cpu)"
    exit 1
    #module load pgi
    #export COMPILERVER="${PGIVER}"
    #export TAUVER="2.29.1"
    #module load tau
    #export TAU_MAKEFILE="${TAU_DIR}/lib/Makefile.tau-pgi-papi-mpi-pgi"
    #module load papi
    ;;

  tau-nv)
    echo "Compiler not yet tested (TODO: write make.inc.basecamp.tau-nv.cpu)"
    exit 1
    ;;

  tau-intel)
    echo "Compiler not yet tested (TODO: write make.inc.basecamp.tau-intel.cpu)"
    exit 1
    ;;

  *)
    echo "Unsupported compiler"
    exit 1
esac

# Copy the appropriate make.inc
# TODO: Write the unavailable make.inc files
case ${USEACC} in
  none)
    cp make.inc.tp13.${COMPILER}.cpu make.inc
    ;;
  cpu)
    echo "OpenACC (on cpu)"
    cp make.inc.tp13.${COMPILER}.acc make.inc
    module load cuda
    module load magma
    ;;
  *)
    echo "Error USEACC=$USEACC"
    exit 1
    ;;
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

  # Load OpenBLAS
  if [ "x${USEOBLAS}" == "x1" ]; then
    echo "Using OpenBLAS"
  fi

  # Load reference BLAS and LAPACK
  if [ "x${USEREFBLAS}" == "x1" ]; then
    echo "Using reference BLAS and LAPACK"
  fi

  # Load Intel MKL
  if [ "x${USEMKL}" == "x1" ]; then
    echo "Using Intel MKL"
    if [ "${COMPILER}" != "intel" ]; then module load mkl; fi
  fi

  # Load FFTW 3
  if [ "x${USEFFTW}" == "x1" ]; then
    echo "Using FFTW 3"
  fi

  # Load HDF5
  if [ "x${USEHDF5}" == "x1" ]; then
    echo "Using HDF5"
  fi

  # Extract link line from make.inc
  if [ "x${USETAU}" == "x1" ]; then
    ${MAKE} lsvars
    source ./libs.sh
    # Apply options
    export TAU_OPTIONS="-optCompInst -optRevert -optTrackIO -optLinking=\"${LIBS}\""
  fi

  # Clean build directory
  ${MAKE} clean
  #rm *.o *.mod

  # Build elk-cpu and check error code
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
      echo "Using HDF5"
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
unset USEHDF5
unset USETAU
unset TAUVER
unset USEOBLAS
unset USEREFBLAS
unset USEMKL

echo; hline; echo;
echo " Done! "
echo; hline; echo;

exit 0
