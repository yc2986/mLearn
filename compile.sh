#!/bin/bash
. /opt/intel/mkl/bin/mklvars.sh intel64
<<<<<<< HEAD
export LD_RUN_PATH="${PWD}/lib:$LD_LIBRARY_PATH"
make
=======
export LD_RUN_PATH="${PWD}/lib:/opt/intel/compilers_and_libraries_2016.2.181/linux/tbb/lib/intel64/gcc4.4:/opt/intel/compilers_and_libraries_2016.2.181/linux/compiler/lib/intel64:/opt/intel/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64"
g++ -w -std=c++11 -DEIGEN_USE_MKL_ALL -fPIC -shared -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -ldl -fopenmp -m64 -I${MKLROOT}/include -I. -Ofast -DNDEBUG source/model.cpp -o lib/libmodel.so
echo model compiled
g++ -w -std=c++11 -Wl,--no-as-needed -Llib -lmodel -I. -Ofast -DNDEBUG source/mLearn.cpp -o mLearn
echo mLearn compiled
>>>>>>> 94c80b15fd61ddbdfe80b47587ec7ee858df16cd
