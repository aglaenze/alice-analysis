#!/bin/bash

numberOfJobs=2
numEvents=1	# this number must be < 700

#SEEDNUM='$(( 100 * $ProcId + $RANDOM ))'  # varies with the job
SEEDNUM='$(ProcId)'  # varies with the job

phraseIn=AliAOD.root
phraseOut=AliAOD'-$(ProcId).root'
dataFolder='.'

# We will need:
# - GeneratorCustom.C
# - OCDBsim.root
# - OCDBrec.root
inputFile1=$dataFolder/GeneratorCustom.C
inputFile2=$dataFolder/OCDBsim.root
inputFile3=$dataFolder/OCDBrec.root

#For ALICE
#alienv enter VO_ALICE@AliPhysics::v5-09-05-01-1,VO_ALICE@AliDPG::prod-201811-01-1
#ALIDPG_ROOT=/cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/AliDPG/prod-201811-01-1
#echo $ALIDPG_ROOT

# First: create all test.hepmc files needed
#source /afs/cern.ch/user/j/jung/work/public/rapgap/setup-rapgap.sh

#for (( i=1; i<=$numberOfJobs; i++ ))
#do
#if [ ! -d hepmcFolder/hepmc$i ];
#then
#mkdir hepmcFolder/hepmc$i
#fi

# run RAPGAP based on a unique steer_file with different seed numbers to get different test.hepmc files
#export RASEED=$(( $i +3 ))
#/afs/cern.ch/user/j/jung/work/public/rapgap/local/bin/rapgap_hepmc < steer_test
#wait
#mv test.hepmc hepmcFolder/hepmc$i
#done
#export RASEED=$SEEDNUM

#wait

if test -f job.sub
then
rm job.sub
fi

# Create the job.sub file
touch job.sub
echo 'executable        = execute.sh' >> job.sub
#echo 'executable        = rapgap.sh' >> job.sub
echo 'arguments         = '$numEvents  >> job.sub
echo 'output            = output/ex.$(ClusterId).$(ProcId).out' >> job.sub
echo 'should_transfer_files = YES' >> job.sub
echo 'input             = '$dataFolder >> job.sub
#echo 'transfer_input_files = test.hepmc' >> job.sub
#echo 'transfer_input_files = rapgap.sh',$ALIDPG_ROOT/,$ALICE_ROOT,$ALIDPG_ROOT/bin/aliroot_dpgsim.sh,$inputFile1,$inputFile2,$inputFile3 >> job.sub
echo 'transfer_input_files = rapgap.sh' >> job.sub
echo '+PreCmd		 = "rapgap.sh"' >> job.sub
echo 'error             = error/ex.$(ClusterId).$(ProcId).err' >> job.sub
echo 'log               = log/ex.$(ClusterId).$(ProcId).log' >> job.sub
echo 'getenv            = true' >> job.sub
echo '+MaxRuntime       = 460000' >> job.sub
echo 'request_cpus      = 4' >> job.sub
echo 'transfer_output_remaps  = "'$phraseIn'='$phraseOut'" ' >> job.sub
#echo 'transfer_output_files = sim.log' >> job.sub
echo 'queue' $numberOfJobs >> job.sub

wait
#exit
condor_submit job.sub

