#!/bin/bash
. /opt/intel/mkl/bin/mklvars.sh intel64
export LD_RUN_PATH="${PWD}/lib:$LD_LIBRARY_PATH"
make