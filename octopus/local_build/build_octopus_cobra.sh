#! /bin/bash -l
# compilation script for octopus
# needs to be executed in a subfolder (e.g. _build) of the root directory
#
# if you want to skip the configure step, call as ./build_octopus.sh noconfig
# if you want to skip the configure step and the dependency checking,
#   call as ./build_octopus.sh noconfig nodep

if [[ ! -f ../configure.ac ]]; then
  echo "Error! Please execute this script in a subfolder of the root directory."
  exit 1
fi

compiler=intel
#compiler=gnu
#cuda=yes

mpi=impi
#mpi=openmpi

echo Using $compiler compiler.
[[ $cuda == yes ]] && echo Building with CUDA support.

module purge
if [[ $compiler == intel ]]; then
  # newest intel version
  module load intel/19.1.3 impi/2019.9 mkl/2020.4 

  #export I_MPI_LINK=dbg
  export I_MPI_LINK=opt

  if [[ $mpi == impi ]]; then
    export CC=mpiicc
    export FC=mpiifort
    export CXX=mpiicpc
    export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
  else
    export CC=icc
    export FC=ifort
    export CXX=icpc
    export MKL="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"
  fi
  export CFLAGS="-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fma -ip -g -traceback"
  export FCFLAGS="$CFLAGS"
  export CXXFLAGS="$CFLAGS"
elif [[ $compiler == gnu ]]; then
  if [[ $mpi == impi ]]; then
    module load gcc/9 intel/19.1.3 impi/2019.9 mkl/2020.4 
    export CC=mpicc
    export FC=mpif90
    export CXX=mpicxx
    export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"
  elif [[ $mpi == openmpi ]]; then
    module load gcc/9 openmpi/4 mkl/2020.4
    export CC=mpicc
    export FC=mpif90
    export CXX=mpicxx
    export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl"
  else
    export CC=gcc
    export FC=gfortran
    export CXX=g++
    export MKL="-L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"
  fi
  #export CFLAGS="-O3 -march=skylake-avx512 -g"
  export CFLAGS="-O3 -march=native -g"
  #export FCFLAGS="$CFLAGS -Wall -Wno-maybe-uninitialized  -Wno-unused-dummy-argument -Wno-c-binding-type -fbacktrace -fcheck=all -fbounds-check -finit-real=snan -ffpe-trap=zero,invalid"
  export FCFLAGS="$CFLAGS -Wall -Wno-maybe-uninitialized  -Wno-unused-dummy-argument -Wno-c-binding-type -fallow-argument-mismatch -fallow-invalid-boz"
  export CXXFLAGS="$CFLAGS"
else
  echo "Compiler $compiler unknown."
  exit 1
fi

module load gsl/2.4 hdf5-serial/1.10.9 netcdf-serial/4.4.1 libxc/5.1.2 metis/5.1 parmetis/4.0 cgal/5.2 boost/1.74
if [[ $compiler == intel ]]; then
  # gcc 10 needed for newer c++ libraries
  module load gcc/9
fi
export CXXFLAGS="$CXXFLAGS -std=c++14"
export LDFLAGS="-Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$ELPA_HOME/lib:$METIS_HOME/lib:$PARMETIS_HOME/lib"

if [[ $cuda == yes ]]; then
  module load cuda/11.2
  # for CUDA-aware MPI:
  if [[ $mpi == openmpi ]]; then
    module load openmpi_gpu/4
  fi
  CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
  LDFLAGS="$LDFLAGS:$CUDA_HOME/lib64"
else
  CUDA_FLAGS=""
fi
module list

export INSTALLDIR=$PWD/installed

if [[ $mpi == impi ]]; then
  ENABLE_MPI=--enable-mpi
  SCALAPACK=$MKL
elif [[ $mpi == openmpi ]]; then
  ENABLE_MPI=--enable-mpi
  SCALAPACK=$MKL
  #export SCALAPACK="/u/sohlmann/sw/src/scalapack-2.1.0/libscalapack.a"
else
  ENABLE_MPI=
  SCALAPACK=no
  PARMETIS_HOME=no
fi

if [[ "$1" != "noconfig" ]]; then
  pushd .. && autoreconf -i && popd
  ../configure $CUDA_FLAGS \
    FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
    FCCPP="cpp -ffreestanding" \
    --prefix=$INSTALLDIR \
    $ENABLE_MPI \
    --enable-openmp \
    --disable-gdlib \
    --enable-shared --disable-static \
    --with-gsl-prefix="$GSL_HOME" \
    --with-libxc-prefix="$LIBXC_HOME" \
    --with-blas="$MKL" \
    --with-lapack="$MKL" \
    --with-blacs="$SCALAPACK" \
    --with-scalapack="$SCALAPACK" \
    --with-netcdf-prefix="$NETCDF_HOME" \
    --with-metis-prefix="$METIS_HOME" \
    --with-parmetis-prefix="$PARMETIS_HOME" \
    --with-boost="$BOOST_HOME" \
    || exit 1
fi

echo "\n\nBuilding octopus...\n"

if [[ "$2" == "nodep" ]]; then
  NODEP=1
else
  NODEP=0
fi
make NODEP=$NODEP -j20 && make -j20 NODEP=$NODEP install || exit 1

mkdir -p $INSTALLDIR/.build.doc/
cp -f config.log $INSTALLDIR/.build.doc/
cp -f $0 $INSTALLDIR/.build.doc/

echo "... done"
