#!/bin/bash

#export exe="elk-pgi-prof-bare"
#export exe="elk-pgi-dbg-acc"
export exe="elk-pgi-prof-acc"

export pfx="crpa-NiO-PM"

# Load NSight Systems 2020.4
#export NSYSDIR=/ccsopen/proj/gen148/eecm/nsight-systems-cli/2020.4.1/target-linux-ppc64le
export NSYSDIR=/ccs/proj/mat201/nsight-systems-cli/2020.4.1/target-linux-ppc64le
export nsys=${NSYSDIR}/nsys

# Only profile rank 0
if [ "x${OMPI_COMM_WORLD_RANK}" == "x0" ]; then
  ${nsys} profile --stats=true -t openacc -o "${pfx}_%q{OMPI_COMM_WORLD_RANK}" ./${exe}
else
  ./${exe}
fi
