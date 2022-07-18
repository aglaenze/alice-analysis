#!/bin/bash

numberOfJobs=10
SEEDNUM=$((100*$ProcId+$RANDOM))  # varies with the job
numEvents=1


alias alienv="/cvmfs/alice.cern.ch/bin/alienv"
alienv enter VO_ALICE@AliPhysics::v5-09-05-01-1,VO_ALICE@AliDPG::prod-201811-01-1
$ALIDPG_ROOT/bin/aliroot_dpgsim.sh  --run 266025 --mode Muon,full --uid 1 --generator Custom --nevents $numEvents

phraseIn=AliAOD.root
phraseOut=AliAOD'-$(ProcId).root'
dataFolder='.'

# We well need:
# - GeneratorCustom.C
# - OCDBsim.root
# - OCDBrec.root
#inputFile1=$dataFolder/GeneratorCustom.C
#inputFile2=$dataFolder/OCDBsim.root
#inputFile3=$dataFolder/OCDBrec.root

if test -f job.sub
then
rm job.sub
fi

# Create the job.sub file
touch job.sub
echo 'executable        = execute.sh' >> job.sub
echo 'arguments         = '$SEEDNUM $numEvents  >> job.sub
echo 'output            = output/ex.$(ClusterId).$(ProcId).out' >> job.sub
echo 'input             = '$dataFolder >> job.sub
#transfer_input_files = input/input-74.sh
#+PreCmd                 = "input-74.sh"
echo 'error             = error/ex.$(ClusterId).$(ProcId).err' >> job.sub
echo 'log               = log/ex.$(ClusterId).$(ProcId).log' >> job.sub
echo 'getenv            = true' >> job.sub
echo '+MaxRuntime       = 460000' >> job.sub
echo 'request_cpus      = 4' >> job.sub
echo 'transfer_output_remaps  = "'$phraseIn'='$phraseOut'" ' >> job.sub
echo 'queue' $numberOfJobs >> job.sub


condor_submit job.sub
