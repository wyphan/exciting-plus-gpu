#!/bin/bash

# Environment variables
OMPI_ROOT=/usr/lib64/openmpi
WDIR=`pwd`

# Command to execute
EXE="elk"

# Processes per node
PPN=8

# Set up hostfile
cat <<__EOF__ > my_hosts
hydra2
hydra3
__EOF__

# Run
mpirun --prefix $OMPI_ROOT -wdir $WDIR -hostfile my_hosts -npernode $PPN $EXE
