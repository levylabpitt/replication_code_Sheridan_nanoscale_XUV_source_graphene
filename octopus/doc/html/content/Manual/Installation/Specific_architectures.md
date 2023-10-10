---
title: "Specific architectures"
#series: "Manual"
---

{{% notice note %}}
This page needs updating!
{{% /notice %}}

Here the {{< manual "Install" "general instructions" >}} for building Octopus are supplemented by instructions for some specific large supercomputers. Included also is how to configure the optional ETSF_IO library, and how you can run the testsuite in parallel.

### Generic Ubuntu 10.10

Install the following packages: 
```bash
 * Required: gfortran, gcc, make, automake, m4, libtool, libgsl0-dev, libblas-dev, liblapack-dev, libfftw3-dev
 * Optional: mpi-default-dev, libgd2-xpm-dev, libsparskit-dev, libnetcdf-dev, libtrilinos-dev
```

[Unfortunately you should not use <tt>libblacs-mpi-dev, libscalapack-mpi-dev</tt> packages because they can give incorrect results. To use BLACS and ScaLAPACK (optional) you must build them yourself.]

First [download libxc](https://www.tddft.org/programs/libxc/download), extract with {{< code "tar xzf" >}}, enter the resulting directory, and execute:

```bash
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3
 make install
 make check
```

For octopus, download and extract, and execute the following (inserting the path where you put libxc):

```bash
 ./configure --prefix=`pwd` CC=gcc CXX=g++ FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3 --with-libxc-prefix=''LIBXC_PATH''
 make install
```

If you are using MPI, replace the compilers by <tt>CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx FC=/usr/bin/mpif90</tt> and add <tt>--enable-mpi</tt>.

To build ETSF_IO (optional), <tt>libnetcdf-dev</tt> is required.

```bash
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCFLAGS="-O3 -ffree-line-length-none" CFLAGS="-O3" --with-netcdf-ldflags="-lnetcdff"
 make install
```

Add <tt>--with-etsf-io-prefix="$DIR/etsf_io-1.0.4"</tt> to the Octopus configure line, where <tt>$DIR/etsf_io-1.0.4</tt> is where you installed ETSF_IO.

(Aug 2018)

<!--
### Generic MacOS

Octopus is now available in MacPorts.
For more info, and how to build by hand with MacPorts libraries, see [[Compiling_octopus_in_OS_X]].
-->

### United States


{{% expand "Sequoia/Vulcan (IBM Blue Gene/Q)" %}}
Supercomputers at Lawrence Livermore National Laboratory based on the IBM Blue Gene/Q architecture. This was tested in Vulcan, but it should work on Sequoia as well.

https://computing.llnl.gov/tutorials/bgq/

```bash
#include_input scripts/build/build_octopus_vulcan.sh
```

(This build does not use Scalapack)
{{% /expand %}}


{{% expand "Stampede" %}}
10 PFLOPS (PF) Dell Linux Cluster based on 6,400+ Dell PowerEdge server nodes, each outfitted with 2 Intel Xeon E5 (Sandy Bridge) processors and an Intel Xeon Phi Coprocessor (MIC Architecture). Texas Advanced Computing Center, University of Texas, Austin, USA.

https://portal.tacc.utexas.edu/user-guides/stampede

```bash
#include_input scripts/build/build_octopus_stampede.sh
```

testsuite script:

```bash
#include_input scripts/test/test_octopus_stampede.sh
```

(13 September 2016)
{{% /expand %}}

{{% expand "Edison" %}}
Cray XC30 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

https://www.nersc.gov/nusers/systems/

```bash
#include_input scripts/build/build_octopus_edison.sh
```

To run the testsuite in parallel:

```bash
#include_input scripts/test/test_octopus_edison.sh
```

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```bash
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo" \
 --with-netcdf-module-path="$NETCDF_DIR/include" --with-netcdf-prefix="$NETCDF_DIR"
 make -j 6 install
 make check
```

(version 8.1, July 2018)
{{% /expand %}}

{{% expand "Cori" %}}
Cray XC40 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

https://www.nersc.gov/users/computational-systems/cori/

```bash
 module load cray-fftw gsl cray-netcdf-hdf5parallel parpack
```

For libxc (you can use the installed version instead)
```bash
 ./configure --prefix=`pwd` CC=cc FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

For octopus (note: using parpack here causes scalapack problems):
```bash
 ./configure --prefix=`pwd` --with-libxc-prefix=/usr/common/software/libxc/4.2.3/intel/haswell --enable-mpi \
 CC=cc CXX=CC FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo" \
 FCCPP="/lib/cpp -ffreestanding" --with-fftw-prefix=$FFTW_DIR/.. \
 --with-arpack="$ARPACK" --with-parpack=no --with-berkeleygw-prefix=/usr/common/software/berkeleygw/2.0
 make -j 6 install
```

To run the testsuite in parallel:

```bash
#include_input scripts/test/test_octopus_cori.sh
```

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```bash
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

(30 July 2018)
{{% /expand %}}


{{% expand "Lonestar" %}}
Dell Linux cluster at Texas Advanced Computing Center (TACC), University of Texas, Austin, USA

Part of National Science Foundation's [https://teragrid.org/ TeraGrid]

https://services.tacc.utexas.edu/index.php/lonestar-user-guide

```bash
#include_input scripts/build/build_octopus_lonestar.sh
```

To build ''ETSF_IO'':

```bash
 module load netcdf
 ./configure --prefix=`pwd` FC=mpif90 --with-netcdf-module-path=$TACC_NETCDF_INC --with-netcdf-ldflags=-L$TACC_NETCDF_LIB FCFLAGS=-O3
 make install
 make check
```

To run the testsuite in parallel:

```bash
#include_input scripts/test/test_octopus_lonestar.sh
```

(2 May 2011)
{{% /expand %}}

{{% expand "Lawrencium" %}}
Dell Linux cluster at Lawrence Berkeley National Laboratory, USA

Cluster available only to Lawrence Berkeley National Laboratory researchers

https://lrc.lbl.gov/html/Lawrencium.html

```bash
#include_input scripts/build/build_octopus_lawrencium.sh
```

To build ''ETSF_IO'':

```bash
 module load netcdf/4.2-intel-p
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3" --with-netcdf-module-path=$NETCDF_DIR/include/ --with-netcdf-ldflags="-L$NETCDF_DIR/lib/ -lnetcdf"
```

(21 Aug 2012)
{{% /expand %}}

### Europe

{{% expand "Curie" %}}

##### PFFT 1.0.5 and BigDFT 1.6.0

The Curie supercomputer, owned by GENCI and operated into the TGCC by CEA, is the first French Tier0 system open to scientists through the French participation into the PRACE research infrastructure.

Curie is offering 3 different fractions of x86-64 computing resources for addressing a wide range of scientific challenges and offering an aggregate peak performance of 2 PetaFlops.

General information: https://www-hpc.cea.fr/en/complexe/tgcc-curie.htm

Compilation with PFFT version 1.0.5 and LibISF taken from BigDFT 1.6.0 (previously compiled):

```bash
#include_input scripts/build/build_octopus_curie_pfft1.0.5.sh
```

BigDFT was configured like this:

```bash
 #!/bin/bash 
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
   --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
   --prefix=$WORKDIR/software/bigdft
```


##### PFFT 1.0.7 and BigDFT 1.7.0

```bash
#include_input scripts/build/build_octopus_curie_pfft1.0.7.sh
```

LibISF compilation:
```bash
 #!/bin/bash
 
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
    --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
    --prefix=$WORKDIR/poisson/local/libisf \
    --disable-libbigdft --disable-binaries --enable-dynamic-libraries
```

FFTW 3.3.3 compilation script is in https://www-user.tu-chemnitz.de/~mpip/software/install_fftw-3.3.3_gcc.sh and it is changed as follow:

```bash
 < myprefix=$HOME/local
 --- 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 23c29,30
 < wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 ---
 > - wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 > cp ../fftw-$FFTW_VERSION.tar.gz .
 876c883
 < ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared \
 ---
 > ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared --disable-shared \ 
```


PFFTW 1.0.7 compilation script is in https://www-user.tu-chemnitz.de/~mpip/software/install_pfft-1.0.7-alpha_gcc.sh and it is changed as follow:

```bash
 < myprefix=$HOME/local
 ---
 > 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 24c31,32
 < wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 ---
 > -wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 > cp ../pfft-$PFFT_VERSION.tar.gz .
```

{{% /expand %}}

{{% expand "Fermi/Juqueen" %}}
IBM Blue Gene/Q, Cineca, Italy

General information: https://www.hpc.cineca.it/content/ibm-fermi-user-guide

The exactly the same script has been tested in the Juqueen, Forschungszentrum Jülich, Jülich Supercomputing Centre (JSC), Jülich, Germany. General information: https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUQUEEN/JUQUEEN_node.html

```bash
#include_input scripts/build/build_octopus_juqueen.sh
```

To build FFTW3. Small changes are made to the script taken from M. Pippig web site https://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
It has a patch and is build with MPI and multithreaded (multithreaded may not work in Juqueen):

```bash
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw-threads-${FFTW_VERSION}"
 TMP="tmp-mine-fftw-${FFTW_VERSION}" 
 
 # bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)" 
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
  fi 
 
 mkdir $TMP && cd $TMP 
 
 wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 # fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
 	  sbo[pe] = (int) (pe * (b * p->tblock) * vn);
 	  rbs[pe] = (int) (db * bt * vn);
 	  rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -	  if (sbs[pe] != (b * p->tblock) * vn
 -	      || rbs[pe] != (p->block * bt) * vn)
 +	  if (dbt != p->tblock
 +	      || db != p->block)
 	       equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC=mpixlc_r F77=mpixlf90_r MPICC=mpixlc_r MPICXX=mpixlcxx_r MPIFC=mpixlf90_r MPILIBS=" " \
 CFLAGS='-O3 -g' \
 FCFLAGS='-O3 -g' \ 
 
 make -j 4
 make install
```

To build PFFT, also taken from https://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
```bash
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft-threads-$PFFT_VERSION
 FFTWDIR=$myprefix/fftw-threads-$FFTW_VERSION
 TMP="tmp-mine-pfft-$PFFT_VERSION" 
 
 # bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)"
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTDIR --with-fftw3=$FFTWDIR \
 MPICC=mpif90 MPICXX=mpicxx MPIFC=mpif90 CC=mpicc FC=mpif90 \
 CFLAGS='-O3 -g -qmaxmem=-1' \
 FCFLAGS='-O3 -g -qmaxmem=-1' \
 --disable-shared
 
 make -j 4
 make install
```
{{% /expand %}}

{{% expand "Jugene" %}}

FZJ-JSC IBM Blue Gene/P, JÜLICH SUPERCOMPUTING CENTRE (JSC), Germany

May also work on other Blue Gene/P computers

General information: https://www.fz-juelich.de/jsc/jugene

Usage information: https://www.fz-juelich.de/jsc/jugene/usage/
```bash
#include_input scripts/build/build_octopus_jugene.sh
```

To run you have to change to $WORK directory and run LoadLeveler. Below is an example of the tutorial (https://www.tddft.org/programs/octopus/wiki/index.php/Tutorial:Benzene_molecule) executing in 32 nodes in the SMP mode (4 threads per chip). The execution mode is hybrid; OpenMP/MPI. 

Here is a example job input script:
```bash
 # @ job_name = Octopus_Sample_1
 # @ comment = "BGP Job by Size"
 # @ error = $(job_name).$(jobid).out
 # @ output = $(job_name).$(jobid).out
 # @ environment = COPY_ALL;
 # @ wall_clock_limit = 00:30:00
 # @ notification = error
 # @ notify_user = foo@gmail.com
 # @ job_type = bluegene
 # @ bg_size = 32
 # @ queue
 mpirun -exe $OCTOPUS_HOME/bin/octopus_mpi -mode SMP -verbose 1 
```
To execute the above script:
```bash
 llsubmit executing_script
```
{{% /expand %}}

{{% expand "MareNostrum II" %}}

MareNostrum is a supercomputer in Europe, the most powerful in Spain. This is one of the seven supercomputers of the Spanish Supercomputing Network. The supercomputer consists of 2,560 JS21 blade computing nodes, each with 2 dual-core IBM 64-bit PowerPC 970MP processors running at 2.3 GHz for 10,240 CPUs in total. The computing nodes of MareNostrum communicate primarily through Myrinet.

You may need to compile GSL and FFTW libraries before starting then installation.

Here is a script to build Octopus:
```bash
#include_input scripts/build/build_octopus_marenostrum2.sh
```

As you can see in the script the --with-gsl-prefix has to point to your GSL and --with-fft-lib to your FFTW. You should compile those with the same CFLAGS and FCFLAGS. Here is an example of configuring GSL (we use GSL-1.11):
```bash
 ./configure --prefix=$outdirGSL CC=gcc -m64 --enable-static --disable-shared
```
We use FFTW-3.1.3. To configure FFTW3 just write:
```bash
 ./configure --prefix=$outdirFFT
```

Recently Marenostrum has changed to the SLURM manager. To execute a job you have to submit it with jobsubmit (more details in the [https://www.bsc.es/media/859.pdf User's Guide]).

{{% /expand %}}

{{% expand "MareNostrum III" %}}

MareNostrum III is a supercomputer based on Intel SandyBridge processors, iDataPlex Compute Racks, a Linux Operating System and an Infiniband interconnection. More information: https://www.bsc.es/marenostrum-support-services/mn3

```bash
#include_input scripts/build/build_octopus_marenostrum3.sh
```
{{% /expand %}}

{{% expand "Corvo" %}}

Corvo is the computational cluster of the Nano-bio Spectroscopy Group. It consists of 1504 processor cores and 5284 GiB of RAM connected by an Infiniband network.

Script to build LIBFM library included in the Scafacos library (modified Scafacos svn version 1920, this is valid since version 9887 of Octopus):

```bash
#include_input scripts/build/build_octopus_corvo.sh
```

The scripts to compile FFTW 3.3.2 and PFFT 1.0.5

```sh
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw/${FFTW_VERSION}"
 TMP="tmp-fftw-${FFTW_VERSION}"
 
 # bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 # fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
           sbo[pe] = (int) (pe * (b * p->tblock) * vn);
           rbs[pe] = (int) (db * bt * vn);
           rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -         if (sbs[pe] != (b * p->tblock) * vn
 -             || rbs[pe] != (p->block * bt) * vn)
 +         if (dbt != p->tblock
 +             || db != p->block)
               equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC="mpicc" F77="mpif90" MPICC="mpicc" \
 CFLAGS="-O3" FFLAGS="-O3" \
 MPILIBS=" "
 
 make -j 4
 make install
```


```bash
 #!/bin/sh -e
 
 myprefix=/home/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft/$PFFT_VERSION
 FFTWDIR=$myprefix/fftw/$FFTW_VERSION
 TMP="tmp-pfft-$PFFT_VERSION"
 module unload fftw/3.2.2
 
 # bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 export CC=mpicc
 export FC=mpif90 
 
 wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR --disable-shared
 
 make install
```

This is the script to compile Octopus with PFFT and FMM support (through Scafacos library):

```bash
#include_input scripts/build/build_octopus_corvo_pfft_fmm.sh
```
{{% /expand %}}

{{% expand "LaPalma" %}}

LaPalma is a supercomputer of the IAC. It is a part of the old MareNostrum II computer. It was not possible to compile with XL compilers, so GNU 4.6.1 version is used.

Previously to Octopus, Libxc (2.2), FFTW (3.3.4), GSL (1.16) and OpenBLAS (0.2.13) were required to compile.

GSL
```bash
 export LD_LIBRARY_PATH=/usr/lib:/gpfs/apps/MPICH2/slurm-2.5.6/64/lib/:/gpfs/apps/GCC/4.9.2/lib64/
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.9.2/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.9.2/bin/gfortran"
 
 ./configure --prefix=$HOME/local/gsl/gnu --enable-static --disable-shared
```

OpenBLAS
```bash
 make PREFIX=$HOME/local/openblas/4.6.1 BINARY=64 CC=/gpfs/apps/GCC/4.6.1/bin/gcc FC=/gpfs/apps/GCC/4.6.1/bin/gfortran
```
Libxc
```bash
 export PATH=/gpfs/apps/MPICH2/mx/default/64/bin:$PATH
 export LD_LIBRARY_PATH=/gpfs/apps/MPICH2/mx/default/64/lib:/gpfs/apps/MPICH2/slurm/64/lib
 
 export LD_LIBRARY_PATH=/gpfs/apps/GCC/4.6.1/lib64/:$HOME/local/openblas/lib:$LD_LIBRARY_PATH
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.6.1/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.6.1/bin/gfortran"
 
 export CFLAGS="-m64 -O3"
 export FCFLAGS=$CFLAGS
 
 ./configure --prefix=$HOME/local/libxc/gnu/4.6.1/2.2
```

Octopus 

```bash
#include_input scripts/build/build_octopus_lapalma.sh
```
{{% /expand %}}

{{% expand "XBES" %}}

Located at the University of Hamburg in Hamburg. It has a login node and a 24 nodes with Intel Xeon. Here a basic configuration script:

```bash
#include_input scripts/build/build_octopus_xbes.sh
```
{{% /expand %}}

{{% expand "Magerit (CeSViMa)" %}}

Magerit is a cluster located in "Centro de Supercomputación y Visualización de Madrid". It consists on 245 nodes eServer BladeCenter PS702 with 2 Power7 processors each of 8 cores at 3'3 GHz (422,4 GFlops) and with 32 GB de RAM.
Here the "optimal" configuration script which uses GNU compilers with ESSL IBM library.

```bash
#include_input scripts/build/build_octopus_magerit.sh
```
{{% /expand %}}


{{% expand "MPCDF systems (cobra, raven, ada)" %}}

Cobra and Raven are the HPC machines of the Max-Planck Compute and Data Facility in Garching. On these clusters, as well as ADA, Octopus can be built with the following script. A copy of this compilation script, its CLI arguments, config.log, and the list of required runtime modules will be all stored at .build.doc/ in installation path for future reference.

```bash
#include_input scripts/build/build_octopus_mpcdf.sh
```

{{% /expand %}}

{{< manual-foot prev="Tutorial:Running Octopus on Graphical Processing Units (GPUs)" next="Manual:Appendix:Porting Octopus and Platform Specific Instructions" >}}
---------------------------------------------
