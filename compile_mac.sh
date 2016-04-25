#!/bin/bash
. /opt/intel/mkl/bin/mklvars.sh intel64
export LD_DYLIB_INSTALL_NAME="${PWD}/lib:/opt/intel//compilers_and_libraries_2016.2.146/mac/tbb/lib:/opt/intel//compilers_and_libraries_2016.2.146/mac/compiler/lib:/opt/intel//compilers_and_libraries_2016.2.146/mac/mkl/lib"
clang++ -std=c++11 -DEIGEN_USE_MKL_ALL -dynamiclib -std=c++11  -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -m64 -I${MKLROOT}/include -I. -Ofast -DNDEBUG source/model.cpp -o lib/libmodel.dylib

install_name_tool -change @rpath/libmkl_intel_lp64.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_intel_lp64.dylib lib/libmodel.dylib
install_name_tool -change @rpath/libmkl_core.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_core.dylib lib/libmodel.dylib
install_name_tool -change @rpath/libmkl_intel_thread.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_intel_thread.dylib lib/libmodel.dylib
install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/compiler/lib/libiomp5.dylib lib/libmodel.dylib

echo Done

clang++ -std=c++11 -DEIGEN_USE_MKL_ALL -Llib/ -lmodel -Iinclude -L${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -m64 -I${MKLROOT}/include -I. -Ofast -DNDEBUG source/mLearn.cpp -o mLearn

install_name_tool -change @rpath/libmkl_intel_lp64.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_intel_lp64.dylib ./mLearn
install_name_tool -change @rpath/libmkl_core.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_core.dylib ./mLearn
install_name_tool -change @rpath/libmkl_intel_thread.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/mkl/lib/libmkl_intel_thread.dylib ./mLearn
install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2016.2.146/mac/compiler/lib/libiomp5.dylib ./mLearn

echo Done