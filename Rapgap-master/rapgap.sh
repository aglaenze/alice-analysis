#!/bin/bash

if test -f test.hepmc
then
rm test.hepmc
fi

#source /afs/cern.ch/user/j/jung/work/public/rapgap/setup-rapgap.sh

source /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_96 x86_64-centos7-gcc8-opt

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/afs/cern.ch/user/j/jung/work/public/hepmc/lib:/afs/cern.ch/user/j/jung/work/public/pythia6/pythia6428/lib:/cvmfs/sft.cern.ch/lcg/releases/LCG_96/GSL/2.5/x86_64-centos7-gcc8-opt/lib

export LHAPDF_DATA_PATH=${LHAPDF_DATA_PATH}:/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:/cvmfs/sft.cern.ch/lcg/releases/LCG_96/MCGenerators/lhapdf/6.2.3/x86_64-centos7-gcc8-opt/share/LHAPDF/

#export RASEED=12314
export HEPMCOUT=test.hepmc


if [ -z $1 ]
then
export RASEED=$RANDOM
else
export RASEED=$1
fi

echo "RASEED = $RASEED"
echo
/afs/cern.ch/user/j/jung/work/public/rapgap/local/bin/rapgap_hepmc < steer_test
