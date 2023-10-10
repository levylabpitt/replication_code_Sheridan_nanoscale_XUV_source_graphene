#!/usr/bin/env bash

# these variables need to be defined:
# HELPER_DIR=~/HUGO/octopus-documentation/scripts/
# OCTOPUS_TOP=~/Octopus/octopus/

cp $OCTOPUS_TOP/testsuite/tutorials/06-octopus_basics-periodic_systems.02-silicon_converged.inp inp

. ./kpts.sh

cp tutorial.sh kpts.sh ktps.log $OCTOPUS_TOP/doc/tutorials/periodic_systems/getting_started_with_periodic_systems/2.k-points
