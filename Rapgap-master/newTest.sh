#!/bin/bash

if [ -z $1 ]
then
echo 'Please type ./newTest.sh $k'
exit
fi

testFolder=test${1}
if [ -d $testFolder ]
then
echo "folder $testFolder already exists"
exit
fi

mkdir $testFolder
cd $testFolder

# copy all useful scripts
cp ../steer_test .
cp ../launchJobs.sh .
cp ../clean.sh .
cp ../execute.sh .
cp ../MergeAODs.C .
cp ../GeneratorCustom.C .

ln -s ../*root .

mkdir output log error
mkdir hepmcFolder
cp ../rapgap.sh .

#source /afs/cern.ch/user/j/jung/work/public/rapgap/setup-rapgap.sh


