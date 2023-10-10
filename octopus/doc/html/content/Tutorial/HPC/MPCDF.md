---
title: "Using the MPCDF Systems"
weight: 1
section: "Tutorial MPCDF_Systems"
series: "Tutorials"
tutorials: ["HPC"]
author: "Sebastian Ohlmann"
description: "Start to use the MPCDF systems"
---


For this course, resources are provided at the MPCDF (Max Planck Computing and
Data Facility), the central computing center of the Max Planck Society.

##  Connecting to MPCDF systems via a gateway

The compute resources at the MPCDF run Linux, and login is only possible via
SSH in combination with two-factor authentication (see below). Any modern
operating system (Linux, MacOS, Windows 10) provides an SSH client which is
typically invoked via the ssh command from a terminal. To access a system at
the MPCDF from the public Internet it is necessary to log into one of the
gateway machines first, and then log from the gateway system into the target
system.

You can connect to the gateway machine `gatezero.mpcdf.mpg.de` and then further
to HPC systems. If you connect for the very first time, you will get a warning
that the host’s authenticity cannot be verified. Provided that you can find
the presented “fingerprint” in the list below, it is safe to answer with “yes”
to continue connecting:
* ED25519 key fingerprint is {{<code "SHA256:qjBJoqcJcCM0LyTqtj09BAxS74u81SizY9zob+XwEOA">}}.
* RSA key fingerprint is {{<code "SHA256:zF/sNLAYqwwRlY3/lhb1A805pGiQiF3GhGP1bBCpvik">}}.

In order to connect to the gateway machine via ssh, you need to setup
two-factor authentication (2FA, see https://docs.mpcdf.mpg.de/faq/2fa.html).
For this, the following steps are needed:
* Visit https://selfservice.mpcdf.mpg.de and log in
* In the menu bar at the top of the page, click “My account > Security”
* Select “Configure 2FA” and provide your password
* Choose a primary token type, recommended is an OTP app (OTP one time password)
* Scan the QR code with the app on your phone
* Validate the token by providing a valid OTP
* Choose a secondary token type

##  ssh config made easy

In order to avoid typing in your password repeatedly, you can configure a
ControlMaster setup for ssh (see
https://docs.mpcdf.mpg.de/faq/2fa.html-do-i-have-to-type-in-an-otp-every-time-i-access-the-secured-systems).
The following snippet can be added to {{<file "~/.ssh/config">}} and should work for Linux
and MacOS (replace YOUR_USER_NAME with your user name):

```bash
# Correctly resolve short names of gateway machines and HPC nodes
Match originalhost gate*,cobra,raven
    CanonicalDomains mpcdf.mpg.de
    CanonicalizeFallbackLocal no
    CanonicalizeHostname yes

# Keep a tunnel open for the day when accessing the gate machines
Match canonical host gate*
    User YOUR_USER_NAME
    Compression yes
    ServerAliveInterval 120
    ControlMaster auto
    ControlPersist 10h
    ControlPath ~/.ssh/master-%C

# Keep a tunnel open for the day when accessing the HPC nodes
Match canonical host cobra*,raven*
    User YOUR_USER_NAME
    Compression yes
    ControlMaster auto
    ControlPersist 10h
    ControlPath ~/.ssh/master-%C
    - OpenSSH >=7.3
    ProxyJump gatezero
    - OpenSSH <=7.2
    -ProxyCommand ssh -W %h:%p gatezero

```

##  Connect to the supercomputer

We will use the supercomputer cobra in this tutorial, so you can use `ssh cobra` to connect via gatezero to the login node of cobra.

More information on cobra can be found at
https://docs.mpcdf.mpg.de/doc/computing/cobra-user-guide.html

So please login to cobra. You will find your home directory at
{{<file "/u/YOUR_USER_NAME">}}, which is on one of two fast parallel file systems (GPFS).
For the purposes of this course, you can store data and run simulations under
this folder. For production runs with serious I/O, please use the ptmp file
system under {{<file "/ptmp/YOUR_USER_NAME">}} because it is larger and more powerful. Be
aware that files that have not been accessed for more than 12 weeks are deleted
on ptmp.

##  Software environment

The software stack on cobra is available via environment modules
(https://docs.mpcdf.mpg.de/doc/computing/software/environment-modules.html).
* `module avail` shows available software packages
* `module load package_name` loads a module, such that binaries are in the path; for many packages, a variable named <PKG>_HOME (with the package name) is exported to be used for compiling and linking codes
* `module unload package_name` unloads a module and cleans the environment variables
* `module purge` unloads all modules and creates a clean environment

To manage the plethora of software packages resulting from all the relevant
combinations of compilers and MPI libraries, we have decided to organize the
environment module system for accessing these packages in a natural hierarchical
manner. Compilers (gcc, intel) are located on the uppermost level, depending
libraries (e.g., MPI) on the second level, more depending libraries on a third
level. This means that not all the modules are visible initially: only after
loading a compiler module, will the modules depending on this become available.
Similarly, loading an MPI module in addition will make the modules depending on
the MPI library available.


In case you know the name of the module you wish to load, but you are not sure
about the available versions or what dependencies need to be loaded first, you
can try to use the ‘find-module’ command. This tool searches for the MODULENAME
string through a list of all installed modules:
```bash
find-module MODULENAME
```

Many software packages, such as octopus, are available as modules. Please run
`find-module octopus` to find all octopus modules and the compiler and MPI
modules you might need to load in order to make a specific version available.

{{<notice note>}}
From version 12 on, {{<octopus>}} is self-contained. For previous versions, you need to load the compiler, etc. beforehand.
{{</notice>}}

Three modules of octopus are available for combinations of compiler and MPI
packages (here for version {{<octopus-version>}}):
* octopus/{{<octopus-version>}}: standard package, compiled against MKL and using FFTs from MKL
* octopus-pfft/{{<octopus-version>}}: compiled against PFFT and FFTW to make the better-scaling PFFT backend for the FFT Poisson solver available (only on raven)
* octopus-gpu/{{<octopus-version>}}: version compiled to run on GPUs

To load the default version of octopus/{{<octopus-version>}}, run:
```bash
module purge
module load octopus/{{<octopus-version>}}
```

Now, octopus is available in the path, so you can confirm the octopus version with
```bash
oct-help -v
```

This should print out the version of octopus:
```
octopus {{<octopus-version>}} (git commit )
```

{{< tutorial-footer >}}
