#!/usr/bin/env bash
#
# Copyright (C) 2020 S. Ohlmann
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#

# Paths.
prefix=/u/davismw/TDDFT-real-time/octopus/local_build/installed
datarootdir=${prefix}/share
pkgdatadir=${datarootdir}/octopus
testsuite="$pkgdatadir/performance_testsuite"

# Usage.
function usage() {
    cat <<EOF

 Copyright (C) 2020 by Sebastian Ohlmann
 (based on script by H. Appel and others)

Usage: oct-run_testsuite.sh [options] [arguments]
    
     -h            this message
     -s            use slurm with salloc
     -l            local run
     -t TESTS      tests to run [default: all]
     -d DIRECTORY  directory where to look for testsuite files
     -p PREFIX     installation prefix [default: $prefix]
     -c            do comparison; two arguments are needed:
                     the reference file and the current file with performance data

Report bugs to <octopus-devel@tddft.org>.
EOF
 exit 0;
}


# Parse command line.

# Some default settings.
local_run="no"
SLURM=${SLURM:-false}
SALLOC=${SALLOC:-true}
TESTS=${TESTS:-all}
NODES=${NODES:-1}
TASKS=${TASKS:-1}
NPROCS=${NPROCS:-1}
comparison="no"

while getopts "hlscp:d:t:" opt ; do
    case "$opt" in
        h) usage;;
        l) local_run="yes";;
        s) SLURM=true;;
        c) comparison="yes";;
        p) prefix="$OPTARG";;
        d) directory="$OPTARG";;
        t) TESTS="$OPTARG";;
        ?) echo "Error parsing arguments" 1>&2; exit 1;;
    esac
done
shift $[ OPTIND - 1 ]

if [ "${local_run}" == "yes" ]; then
    bin_directory=$(pwd)/../../src
    if [ -n "$directory" ]; then
        testsuite="$directory"
    else
        testsuite=$(pwd)
    fi
    bin_testsuite=$testsuite
else
    bin_directory=$prefix/bin
    if [ -n "$directory" ]; then
        testsuite="$directory"
    else
        testsuite=$prefix/share/octopus/testsuite/performance
    fi
    bin_testsuite=$bin_directory
fi

run_directory=${run_directory:-.}

# run comparison script?
if [ "${comparison}" == "yes" ]; then
    exec python $testsuite/compare_results.py $@
fi

# check for binary dir
if [ ! -d "${bin_directory}" ]; then
    echo "Specified binary directory '$bin_directory' does not exist." 1>&2
    exit 1
fi

echo "*****************************************"
echo "  Running octopus performance testsuite  "
echo "*****************************************"

if [[ "$SLURM" == "true" && "$SALLOC" == "true" ]]; then
  cd $run_directory
  salloc --nodes $NODES --ntasks-per-node $(($TASKS/$NODES)) -K15 make -f $testsuite/Makefile.performance -j $(($TASKS/$NPROCS)) bin_directory=$bin_directory bin_testsuite=$bin_testsuite testsuite_dir=$testsuite
else
  cd $run_directory
  make -f $testsuite/Makefile.performance -j $(($TASKS/$NPROCS)) bin_directory=$bin_directory bin_testsuite=$bin_testsuite testsuite_dir=$testsuite
fi
