#!/bin/bash

SEEDNUM=$1  # first argument given to the script
numEvents=$2

if [ -z $2 ]
then
echo 'Need a $RANDOM_SEED and $NUM_EVENTS as argument'
exit
fi

## end of variables
source /afs/cern.ch/user/j/jung/work/public/rapgap/setup-rapgap.sh
wait
export RASEED=$SEEDNUM

# run RAPGAP based on a unique steer_file with different seed numbers to get different test.hepmc files
/afs/cern.ch/user/j/jung/work/public/rapgap/local/bin/rapgap_hepmc < steer_test
wait

alias alienv="/cvmfs/alice.cern.ch/bin/alienv"
alienv enter VO_ALICE@AliPhysics::v5-09-05-01-1,VO_ALICE@AliDPG::prod-201811-01-1
$ALIDPG_ROOT/bin/aliroot_dpgsim.sh  --run 266025 --mode Muon,full --uid 1 --generator Custom --nevents $numEvents
