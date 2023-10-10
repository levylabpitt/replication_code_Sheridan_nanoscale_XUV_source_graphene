---
title: "MacOS installation"
#series: "Manual"
hidden: True
---

### MacOS Monterey 12.6

Here the {{< manual "Install" "general instructions" >}} for building Octopus on MacOS. 
This guide assumes you are using the latest MacOS version, i.e. Monterey 12.6.

This guide was tested on an Apple M1 Pro chip, but it should be compatible with an Inter processor. However, you should check the installation folder for Homebrew.

#### Before Installing

Install [HomeBrew](https://brew.sh/). To do that, just follow the procedure on the website. 

You may also need to install Xcode Command Line Tools. In case you need it, it is enough to open the Terminal app and execute {{<code "xcode-select --install">}} 

Find out the architecture of your CPU: To do that it’s enough to click on the Apple logo in the top left corner and select the first item in the drop-down (About this Mac…). Then look at the information for the chip, if it says Apple M1/M1 Pro/M2 then the architecture is arm64, otherwise, it’s x86_64 (Intel CPU)

Finally, **remember to check on Homebrew website for the compatibility of the packages with your architecture** To do it, got to https://formulae.brew.sh/ and type the formula you want to install (e.g. libxc).

Now open the Terminal app and follow the procedure below:

#### Install required and optinal libraries

##### Install compilers

 * GCC

   Execute:
   ```bash
   brew install gcc
   ```

 * OpenMPI

   Follow this step only if you wish to run Octopus in parallel (which is recommended). Execute:
   ```bash
   brew install open-mpi
   ```

 * Autotools

   Execute:
   ```bash
   brew install autoconf
   brew install automake
   brew install libtool
   ```

##### Install libraries

 * LibXC
  
   This library contains the exchange and correlations functionals. This library is required. Execute:
   ```bash
   brew install libxc
   ```

 * GSL

   This is a numerical library for C and C++. This library is required. Execute:
   ```bash
   brew install gsl
   ```

 * FFTW
 
   This is a library for computing discrete Fourier transform. This library is required. Execute:
   ```bash
   brew install fftw
   ```

 * Lapack and BLAS

   These are numerical libraries for linear algebra (Linear Algebra PACKage and Basic Linear Algebra Subprograms). These library are required. Execute:
   ```bash
   brew install lapack
   ```

 * Scalapack and BLACS

   These are numerical libraries for parallel linear algebra operations (SCAlable Linear Algebra PACKage and Basic Linear Algebra Communication Subprograms). 
   These library are required to run Octopus in parallel. Install them only if you installed OpenMPI. Execute:
   ```bash
   brew install lapack
   ```

#### Clone and compile Octopus

Now we shall start compiling Octopus. You may want to restart your Mac. Then, create a folder in your favorite location and open a Terminal window in that folder. To do this, type `cd ` then drag and drop the folder from the Finder to the Terminal. <br>

##### Get the source code

**GIT CLONE:** To use git clone you will need to create an SSH key and a GitLab account. Then execute: 
```bash
git clone https://gitlab.com/octopus-code/octopus.git
```

**DOWNLOAD SOURCE:** Go to the {{<versioned-link "releases" "releases">}} and download the lastest version. 
Place the downloaded tsr.gz file and place it in the folder you created. Then extract all files by double clicking or using the terminal command `tar xzf octopus-{{<octopus-version>}}.tar.gz`

##### Compile the code

Execute:
```bash
mkdir local_build
cd local_build
cp ../octopus/scripts/build/build_octopus_macos.sh ./
source build_octopus_macos.sh
```

Now Octopus was installed in {{<file "../installed/bin">}}

#### Execute Octopus
Create a folder and a file called “inp” in it. Put all options (https://www.octopus-code.org/documentation/main/) required for your calculation.

{{%notice warning "IMPORTANT"%}}
Due to some error with rapidxml library, please always set: `PseudopotentialSet = 0` <br>
You can still use a pseudopotential by: 

 1. Using the species block:
 ```text
 %Species
  'C' | species_pseudo | file | “path_to_file”
 %
 ```
 2. Setting
 ```text
 PseudopotentialSet = hscv_lda
 ```
{{%/notice%}}

Execute Octopus in serial by running:
```bash
<path_to_octopus_installation>/octopus
```

Execute Octopus in parallel by running 
```bash
/opt/homebrew/bin/mpirun -n 4 <path_to_octopus_installation>/octopus
```
where `-n 4` tells the number of processors you wish to use (in this case 4).

