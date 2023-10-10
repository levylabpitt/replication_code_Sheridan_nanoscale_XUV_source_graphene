#!/bin/bash
#
# Compilation script for octopus. If this script is not executed in a subdirectory (e.g. _build) of the root octopus directory,
# you have to pass the path to octopus directory to the compilation script using the --root-path flag
# See ./build_octopus.sh --help for complete list of available options.
#
# Usage examples:
#
# build octopus for CPU, using gnu compilers and openmpi, with octopus at "$PWD/.." and install it to "$PWD/installed"
# > ./build_octopus_mpcdf.sh
#
# build octopus for CPU, using intel compilers and intel-mpi, install it to "/u/user/octopus/"
# > ./build_octopus_mpcdf.sh --compiler=intel --mpi=impi --prefix=/u/user/octopus/
#
# build octopus for CPU without MPI (with OpenMP support only)
# > ./build_octopus_mpcdf.sh --no-mpi
#
# build octopus for CPU + GPU, using gnu compilers and cuda-aware openmpi (if possible)
# > ./build_octopus_mpcdf.sh --cuda
#
# build octopus for CPU + GPU, using gnu compilers and cuda-unaware(!) openmpi
# > ./build_octopus_mpcdf.sh --cuda --no-cuda-aware-mpi
#
# rebuild without running the configuration script
# > ./build_octopus_mpcdf.sh --no-config
#
# rebuild without compiling the dependencies (WARNING: DO NOT USE THIS AFTER CHANGING THE FUNCTION/SUBROUTINE INTERFACES)
# > ./build_octopus_mpcdf.sh --no-dep
#
#
# NOTE: You must load the required runtime modules in your job submission script before
#       executing the octopus for your calculations. The list of these modules is printed
#       at the end of compilation and also is saved in your installation path at
#       .build.doc/required_runtime_modules file. You can also generate this list by passing
#       '--show-modules-only' flag to this script without compiling the code.

set -e

# version of this script
VERSION="1.4"


# default value of cli arguments
_arg_config="on"
_arg_dep="on"
_arg_cuda="off"
_arg_cuda_aware_mpi="on"
_arg_debug="off"
_arg_compiler="gnu"
_arg_mpi="openmpi"
_arg_prefix="$PWD/installed"
_arg_root_path="${PWD%/*}"
_arg_show_modules_only="off"


# default modules
# Latest tests:
## by MFT on COBRA, RAVEN & ADA clusters @01.2023 using octopus: fd2983011b37c086e7fd979e99afef197beab7a7
# other clusters may have different version of these modules
# update the versions if necessary!
#
## for OpenMP only
modules_gcc="gcc/11"
modules_intel="intel/21.5.0"
modules_libs_openmp="mkl/2022.0 gsl hdf5-serial netcdf-serial etsf_io libxc/5.1.5 libvdwxc cgal boost/1.74"
## for MPI
modules_intel_impi="intel/21.5.0 impi/2021.5"
modules_gcc_impi="gcc/11 impi/2021.5"
modules_gcc_openmpi="gcc/11 openmpi/4"
modules_libs_mpi="mkl/2022.0 gsl hdf5-serial netcdf-serial etsf_io libxc/5.1.5 libvdwxc-mpi cgal boost/1.74 elpa/mpi/standard/2021.05.001 metis-64 parmetis-64 dftbplus/22.1"
## For GPU
modules_cuda="cuda/11.4"
modules_cuda_aware_mpi="openmpi_gpu/4"


# configuration features to verify
## common features we expect to detect in all cases
common_features_list=('CGAL' 'NETCDF' 'ETSF_IO' 'FFTW3' 'LIBXC5' 'LIBXC_FXC' 'LIBXC_KXC' 'LIBVDWXC' 'OPENMP')
## features of mpi version
mpi_features_list=('MPI' 'ELPA' 'PARMETIS' 'LIBVDWXC_MPI') # TODO: add "DFTBPLUS" to this list after https://gitlab.com/octopus-code/octopus/-/issues/583 is solved
## features of cuda version
cuda_features_list=('CUDA')

expected_features_list=( "${common_features_list[@]}" )


die()
{
	_die_ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_die_ret}"
}

print_help()
{
	printf '%s\n' "Compilation script for octopus. This script needs to be executed in a subdirectory (e.g. _build) of the root directory or the path to octopus root directory must be set using --root-path flag."
	print_usage
}

print_usage()
{
	printf '\nUsage: %s [--(no-)config] [--(no-)dep] [--(no-)cuda] [--(no-)cuda-aware-mpi] [--(no-)debug] [--compiler <arg>] [--no-mpi] [--mpi <arg>] [--prefix <arg>] [--show-modules-only] [-v|--version] [-h|--help]\n' "$0"
	printf '\t%s\n' "--config, --no-config: run configure script (on by default)"
	printf '\t%s\n' "--dep, --no-dep: compile dependencies (on by default)"
	printf '\t%s\n' "--cuda, --no-cuda: compile with CUDA support (off by default)"
	printf '\t%s\n' "--cuda-aware-mpi, --no-cuda-aware-mpi: compile with CUDA aware MPI support (on by default)"
	printf '\t%s\n' "--debug, --no-debug: compile with debug flags [-O0, I_MPI_LINK=dbg] (off by default)"
	printf '\t%s\n' "--compiler: compiler type. available options: 'intel' or 'gnu' (default: 'gnu')"
	printf '\t%s\n' "--no-mpi: compile without the MPI support (equivalent to '--mpi=off')"
	printf '\t%s\n' "--mpi: MPI implementation. available options: 'impi', 'openmpi', or 'off' (default: 'openmpi')"
	printf '\t%s\n' "--prefix: installation path (default: './installed/')"
	printf '\t%s\n' "--root-path: path to the root of octopus repository (default: '../')"
	printf '\t%s\n' "--show-modules-only: print the required runtime modules and exit"
	printf '\t%s\n' "-v, --version: print compilation script version"
	printf '\t%s\n' "-h, --help: print help (this list!)"
}


parse_commandline()
{
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			--no-config|--config)
				test "${1:0:5}" = "--no-" && _arg_config="off"
				;;
			--no-dep|--dep)
				test "${1:0:5}" = "--no-" && _arg_dep="off"
				;;
			--no-cuda|--cuda)
				test "${1:0:5}" = "--no-" || _arg_cuda="on"
				;;
			--no-cuda-aware-mpi|--cuda-aware-mpi)
				test "${1:0:5}" = "--no-" && _arg_cuda_aware_mpi="off"
				;;
			--no-debug|--debug)
				test "${1:0:5}" = "--no-" || _arg_debug="on"
				;;
			--compiler)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_compiler="$2"
				shift
				;;
			--compiler=*)
				_arg_compiler="${_key##--compiler=}"
				;;
			--no-mpi)
				_arg_mpi="off"
				;;
			--mpi)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_mpi="$2"
				shift
				;;
			--mpi=*)
				_arg_mpi="${_key##--mpi=}"
				;;
			--prefix)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_prefix="$2"
				shift
				;;
			--prefix=*)
				_arg_prefix="${_key##--prefix=}"
				;;
			--root-path)
				test $# -lt 2 && die "Missing value for the optional argument '$_key'." 1
				_arg_root_path="$2"
				shift
				;;
			--root-path=*)
				_arg_root_path="${_key##--root-path=}"
				;;
			--show-modules-only)
				_arg_show_modules_only="on"
				;;
			-v|--version)
				echo Octopus compilation script v$VERSION
				exit 0
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_PRINT_HELP=yes die "FATAL ERROR: Unexpected argument '$1'" 1
				;;
		esac
		shift
	done
}

# verify features are detected in config.h as "#define HAVE_*** 1"
verify_configuration()
{
	read -ra _expected_features <<< "$@"

	for feature in "${_expected_features[@]}"; do
		printf "checking feature: %s ..." "$feature"
		_definition="#define HAVE_$feature 1"
		grep -q "$_definition" config.h || ( echo "NOT FOUND!" && echo "ERROR: Configuration script failed to detect $feature" && exit 1 )
		echo "OK!"
	done

}

parse_commandline "$@"

# expand ~ in path arguments
_arg_prefix="${_arg_prefix/#\~/$HOME}"
_arg_root_path="${_arg_root_path/#\~/$HOME}"

echo "====================================="
[ "$_arg_show_modules_only" = "off" ] && echo "run configure: $_arg_config"
[ "$_arg_show_modules_only" = "off" ] && echo "build dependencies: $_arg_dep"
echo "compiler: $_arg_compiler"
echo "mpi: $_arg_mpi"
echo "cuda support: $_arg_cuda"
[ "$_arg_mpi" = "openmpi" ] && [ "$_arg_cuda" = "on" ] && echo "cuda-aware-mpi: $_arg_cuda_aware_mpi"
[ "$_arg_show_modules_only" = "off" ] && echo "debug mode: $_arg_debug"
[ "$_arg_show_modules_only" = "off" ] && echo "installation path: $_arg_prefix"
[ "$_arg_show_modules_only" = "off" ] && echo "octopus repository path: $_arg_root_path"
echo "====================================="

if [ ! -f "$_arg_root_path/configure.ac" ]; then
	echo "Error! Could not find configure.ac in octopus repository path ($_arg_root_path). Please pass the right path to the compilation script using the \"--root-path\" flag"
	exit 1
fi


# validating script arguments

if [ "$_arg_compiler" != "intel" ] && [ "$_arg_compiler" != "gnu" ]; then
	echo "invalid --compiler option: $_arg_compiler"
	print_usage
	exit 1
fi

if [ "$_arg_mpi" != "impi" ] && [ "$_arg_mpi" != "openmpi" ] && [ "$_arg_mpi" != "off" ]; then
	echo "invalid --mpi option: $_arg_mpi"
	print_usage
	exit 1
fi

if [ "$_arg_compiler" = "intel" ] && [ "$_arg_mpi" = "openmpi" ]; then
	echo "WARNING: Compiling Octopus with intel compilers and linking to OpenMPI library is not supported!"
	echo "Will use IntelMPI (impi) instead!"
	_arg_mpi="impi"
fi



module purge

if [ "$_arg_compiler" = "intel" ]; then
	if [ "$_arg_mpi" = "off" ]; then

		module load $modules_intel
		required_modules+=" ${modules_intel}"

		export CC=icc
		export FC=ifort
		export CXX=icpc
		MKL_FLAGS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl"

	elif [ "$_arg_mpi" = "impi" ]; then

		module load $modules_intel_impi
		required_modules+=" ${modules_intel_impi}"

		if [ "$_arg_debug" = "on" ]; then
			export I_MPI_LINK=dbg
		else
			export I_MPI_LINK=opt
		fi

		export CC=mpiicc
		export FC=mpiifort
		export CXX=mpiicpc
		MKL_FLAGS="-lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"

	elif [ "$_arg_mpi" = "openmpi" ]; then
		# we should never reach here! Adapt the script if we support Intel compilers + OpenMPI in future.
		echo "ERROR: Compiling Octopus with intel compilers and linking to OpenMPI library is not currently supported!"
		exit 1
	fi

	if [ "$_arg_debug" = "on" ]; then
		export CFLAGS="-O0 -g -traceback"
	else
		export CFLAGS="-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fma -ip -g -traceback"
	fi
	export FCFLAGS="$CFLAGS"
	export CXXFLAGS="$CFLAGS"

elif [ "$_arg_compiler" = "gnu" ]; then
	if [ "$_arg_mpi" = "off" ]; then

		module load $modules_gcc
		required_modules+=" ${modules_gcc}"

		export CC=gcc
		export FC=gfortran
		export CXX=g++
		MKL_FLAGS="-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl"

	elif [ "$_arg_mpi" = "impi" ]; then

		module load $modules_gcc_impi
		required_modules+=" ${modules_gcc_impi}"

		export CC=mpigcc
		export FC=mpigfortran
		export CXX=mpig++
		MKL_FLAGS="-lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"

	elif [ "$_arg_mpi" = "openmpi" ]; then

		module load $modules_gcc_openmpi
		required_modules+=" ${modules_gcc_openmpi}"

		export CC=mpicc
		export FC=mpif90
		export CXX=mpicxx
		MKL_FLAGS="-lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -lgomp -lpthread -lm -ldl"
	fi

	if [ $_arg_debug = on ]; then
		export CFLAGS="-O0 -g"
	else
		export CFLAGS="-O3 -march=native -g"
	fi

	export FCFLAGS="$CFLAGS -Wall -Wno-maybe-uninitialized -Wno-unused-dummy-argument -Wno-c-binding-type -fallow-argument-mismatch -fallow-invalid-boz"
	export CXXFLAGS="$CFLAGS"
fi

if [ "$_arg_mpi" = "off" ]; then
	module load $modules_libs_openmp
	export MKL="-L${MKLROOT}/lib/intel64 ${MKL_FLAGS}"
	ENABLE_MPI=
	SCALAPACK=no
	PARMETIS_HOME=no
	METIS_HOME=no
	ELPA_HOME=no
	DFTBPLUS_HOME=no
else
	module load $modules_libs_mpi
	export MKL="-L${MKLROOT}/lib/intel64 ${MKL_FLAGS}"
	expected_features_list+=( "${mpi_features_list[@]}" )
	ENABLE_MPI=--enable-mpi
	SCALAPACK=$MKL

fi


if [ "$_arg_compiler" = "intel" ]; then
	# gcc 10 or later is needed for newer C++ libraries. This needs to be loaded after $modules_libs_...
	module load gcc/11
fi

if [ "$_arg_cuda" = "on" ]; then
	module load $modules_cuda
	required_modules+=" ${modules_cuda}"
	expected_features_list+=( "${cuda_features_list[@]}" )

	if [ "$_arg_mpi" = "openmpi" ] && [ "$_arg_cuda_aware_mpi" = "on" ] ; then
		# on COBRA the underlying network drivers do not support CUDA-aware MPI
		module -s load $modules_cuda_aware_mpi && required_modules+=" ${modules_cuda_aware_mpi}"\
		|| echo "WARNING: CUDA-aware MPI module does not seem to be available in the current environment. Proceeding without the CUDA-aware MPI support!"
	fi
	CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
	LDFLAGS="$LDFLAGS:$CUDA_HOME/lib64"
else
	CUDA_FLAGS=""
fi

if [ "$_arg_show_modules_only" = "on" ]; then
	echo -e "\nRequired runtime modules for the selected configuration would be: ${required_modules}\n"
	exit 0
fi

module list

export CXXFLAGS="$CXXFLAGS -std=c++14"
export LDFLAGS="-Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$ELPA_HOME/lib:$METIS_HOME/lib:$PARMETIS_HOME/lib:$LIBXC_HOME/lib"

if [ "$_arg_config" != "off" ]; then
	pushd "$_arg_root_path" && autoreconf -i && popd
	"$_arg_root_path"/configure $CUDA_FLAGS \
	FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
	--prefix="$_arg_prefix" \
	$ENABLE_MPI \
	--enable-openmp \
	--disable-gdlib \
	--enable-shared --disable-static \
	--with-gsl-prefix="$GSL_HOME" \
	--with-cgal="$CGAL_HOME" \
	--with-libxc-prefix="$LIBXC_HOME" \
	--with-libvdwxc-prefix="$LIBVDWXC_HOME" \
	--with-blas="$MKL" \
	--with-lapack="$MKL" \
	--with-blacs="$SCALAPACK" \
	--with-scalapack="$SCALAPACK" \
	--with-netcdf-prefix="$NETCDF_HOME" \
	--with-etsf-io-prefix="$ETSF_IO_HOME" \
	--with-metis-prefix="$METIS_HOME" \
	--with-parmetis-prefix="$PARMETIS_HOME" \
	--with-boost="$BOOST_HOME" \
	--with-elpa-prefix="$ELPA_HOME" \
	--with-dftbplus-prefix="$DFTBPLUS_HOME"

	verify_configuration "${expected_features_list[@]}"
fi


if [ "$_arg_dep" = "off" ]; then
	NODEP=1
else
	NODEP=0
fi

make NODEP=$NODEP -j20 && make -j20 NODEP=$NODEP install

# keep the compilation script, its cli arguments, config.log, and the list of required runtime modules
# at .build.doc/ in installation path for future reference
mkdir -p "$_arg_prefix"/.build.doc/
cp -f config.log "$_arg_prefix"/.build.doc/
cp -f "$0" "$_arg_prefix"/.build.doc/
echo "$0 $*" > "$_arg_prefix"/.build.doc/compilation_command
echo "$required_modules" > "$_arg_prefix"/.build.doc/required_runtime_modules
echo "====================================="
echo "Compilation done!"
echo "A copy of this compilation script, its CLI arguments, config.log, and the list of required runtime modules are all stored in $_arg_prefix/.build.doc/ for future reference"
echo -e "\nRequired runtime modules: ${required_modules}\n"
