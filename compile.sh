#!/bin/bash
. /opt/intel/mkl/bin/mklvars.sh intel64
export LD_RUN_PATH="${PWD}/lib:$LD_LIBRARY_PATH"
g++ -w -std=c++11 -DEIGEN_USE_MKL_ALL -fPIC -shared -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp -m64 -I${MKLROOT}/include -I. -Ofast -DNDEBUG source/model.cpp -o lib/libmodel.so
echo model compiled
g++ -w -std=c++11 -Wl,--no-as-needed -Llib -lmodel -I. -Ofast -DNDEBUG source/mLearn.cpp -o mLearn
echo mLearn compiled