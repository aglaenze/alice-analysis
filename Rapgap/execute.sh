#!/bin/bash

numEvents=$1

if [ -z $1 ]
then
echo 'Needs $NUM_EVENTS as argument'
exit
fi

if test -f test.hepmc
then
echo "test.hepmc found, continuing"
else
echo "test.hepmc does not exist, exit"
exit
fi

$ALIDPG_ROOT/bin/aliroot_dpgsim.sh  --run 266025 --mode Muon,full --uid 1 --generator Custom --nevents $numEvents
#./aliroot_dpgsim.sh  --run 266025 --mode Muon,full --uid 1 --generator Custom --nevents $numEvents
