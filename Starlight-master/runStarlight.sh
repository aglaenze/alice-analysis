#!/bin/bash

## variables
N_EVENTS=1000000
#N_EVENTS=0
runStarlight=1	## 0 for testing, 1 to run Starlight


start=`date +%s`

clean()
{
filesToDelete="*.so *.d *.pcm *ACLiC* *.in *.out *dict*"
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}

errorMessage()
{
echo
echo 'Error: please type ./execute.sh $year $prodMode $config'
echo 'with'
echo '--> year = 2013 or 2016'
echo -e '--> production mode =\t 1 (kIncohJpsiToMu)'
echo -e '\t\t\t 2 (kCohJpsiToMu)'
echo -e '\t\t\t 3 (kIncohPsi2sToMu)'
echo -e '\t\t\t 4 (kTwoGammaToMu)'
echo -e '\t\t\t 5 (kCohJpsiToMubis) --> production mode = 3 instead of 2'
echo '--> config = 1 (p-Pb) or 2 (Pb-p) or 3 (Pb-Pb)'
echo
exit
}

write()
{
echo 'baseFileName =' $baseFileName >> $slightIn
#echo 'OUTPUT_HEADER = 1' >> $slightIn
echo 'BEAM_1_Z =' $BEAM_1_Z'    #Z of projectile' >> $slightIn
echo 'BEAM_1_A =' $BEAM_1_A'   #A of projectile' >> $slightIn
echo 'BEAM_2_Z =' $BEAM_2_Z'   #Z of target' >> $slightIn
echo 'BEAM_2_A =' $BEAM_2_A'   #A of target' >> $slightIn
echo 'BEAM_1_GAMMA = '$BEAM_1_GAMMA'	 #Gamma of the colliding ion 1' >> $slightIn
echo 'BEAM_2_GAMMA = '$BEAM_2_GAMMA'	#Gamma of the colliding ion 2' >> $slightIn
echo 'W_MAX =' $W_MAX'   #Max value of w' >> $slightIn
echo 'W_MIN =' $W_MIN'    #Min value of w' >> $slightIn
echo 'W_N_BINS =' $W_N_BINS'    #Bins i w' >> $slightIn
echo 'RAP_MAX = 5    #max y' >> $slightIn
echo 'RAP_N_BINS = 1001    #Bins in y' >> $slightIn
#echo 'RAP_N_BINS = 900    #Bins in y' >> $slightIn
echo 'CUT_PT = 0 #Cut in pT? 0 = (no, 1 = yes)' >> $slightIn
echo 'PT_MIN = 0.0 #Minimum pT in GeV' >> $slightIn
echo 'PT_MAX = 20.0 #Maximum pT in GeV' >> $slightIn
echo 'CUT_ETA = 0 #Cut in pseudorapidity? (0 = no, 1 = yes)' >> $slightIn
echo 'ETA_MIN = -4 #Minimum pseudorapidity' >> $slightIn
echo 'ETA_MAX = -2.5 #Maximum pseudorapidity' >> $slightIn
echo 'PROD_MODE =' $PROD_MODE'     	#gg or gP switch (1 = 2-photon, 2 = coherent vector meson (narrow), 3 = coherent vector meson (wide), 4 = incoherent vector meson)' >> $slightIn
echo 'N_EVENTS =' $N_EVENTS'   	#Number of events' >> $slightIn
echo 'PROD_PID =' $PROD_PID'   	#Channel of interest' >> $slightIn
echo 'RND_SEED =' $RND_SEED' 	#Random number seed' >> $slightIn
#echo 'OUTPUT_FORMAT = 2     #Form of the output' >> $slightIn
echo 'BREAKUP_MODE = 5     #Controls the nuclear breakup' >> $slightIn
echo 'INTERFERENCE = 0     #Interference (0 = off, 1 = on)' >> $slightIn
echo 'IF_STRENGTH = 1.    #% of interference (0.0 - 0.1)' >> $slightIn
#echo 'COHERENT = 1     #Coherent=1,Incoherent=0' >> $slightIn
#echo 'INCO_FACTOR = 1.    #percentage of incoherence' >> $slightIn
#echo 'BFORD = 9.5' >> $slightIn
echo 'INT_PT_MAX = 0.24  #Maximum pt considered, when interference is turned on' >> $slightIn
echo 'INT_PT_N_BINS = 120   #Number of pt bins when interference is turned on' >> $slightIn
echo 'XSEC_METHOD = 1' >> $slightIn
echo 'N_THREADS = 8' >> $slightIn
#echo 'PYTHIA_FULL_EVENTRECORD = 1' >> $slightIn
echo 'PRINT_VM = 2' >> $slightIn
}

## Tests parameters

if [ -z $3 ]
then
errorMessage
fi

# tests on the first parameter: the year
if [ $1 = 2013 ]
then
BEAM_1=4263.2
BEAM_2=1680.7
elif [ $1 = 2016 ]
then
#BEAM_1=6927.0
#BEAM_1=6929.6
BEAM_1=6928		# p
#BEAM_2=2731.1
#BEAM_2=2751.9
BEAM_2=2745		# Pb
else
errorMessage
fi

# tests on the second parameter: the process
if [ $2 = 1 ]
then
process=kIncohJpsiToMu
PROD_MODE=4
PROD_PID=443013
W_MAX=3.09738
W_MIN=3.09645
W_N_BINS=20
RND_SEED=578537
elif [ $2 = 2 ]
then
process=kCohJpsiToMu
PROD_MODE=2
PROD_PID=443013
W_MAX=3.09738
W_MIN=3.09645
W_N_BINS=20
RND_SEED=578537
elif [ $2 = 3 ]
then
process=kIncohPsi2sToMu
PROD_MODE=4
PROD_PID=444013
W_MAX=3.6876
W_MIN=3.68461
W_N_BINS=20
RND_SEED=578537
elif [ $2 = 4 ]
then
process=kTwoGammaToMu
PROD_MODE=1
PROD_PID=13
W_MAX=5.0
W_MIN=0.8
W_N_BINS=169
#W_MIN=1.2
#W_MAX=15.0
#W_N_BINS=$(bc <<< "($W_MAX - $W_MIN)*10+1")
#W_N_BINS=500
#RND_SEED=912665125
RND_SEED=5574533
#echo $W_N_BINS
elif [ $2 = 5 ]
then
process=kCohJpsiToMuBis
PROD_MODE=3
PROD_PID=443013
W_MAX=3.09738
W_MIN=3.09645
W_N_BINS=20
RND_SEED=578537
else
errorMessage
fi

# tests on the third parameter: the configuration (p-Pb or Pb-p)
if [ $3 = 1 ]
then
config="p-Pb"
BEAM_1_GAMMA=$BEAM_1
BEAM_2_GAMMA=$BEAM_2
BEAM_1_Z=1
BEAM_1_A=1
BEAM_2_Z=82
BEAM_2_A=208
elif [ $3 = 2 ]
then
config="Pb-p"
BEAM_1_GAMMA=$BEAM_2
BEAM_2_GAMMA=$BEAM_1
BEAM_1_Z=82
BEAM_1_A=208
BEAM_2_Z=1
BEAM_2_A=1
elif [ $3 = 3 ]
then
config="Pb-Pb"
BEAM_1_GAMMA=$BEAM_2
BEAM_2_GAMMA=$BEAM_2
BEAM_1_Z=82
BEAM_1_A=208
BEAM_2_Z=82
BEAM_2_A=208
else
errorMessage
fi

# Action starts here
echo
echo Process: $process
echo 'Configuration' $config
echo

year=$1
baseFileName=slight-$year-$process-$config
slightIn2="/Users/aglaenzer/Softwares/Starlight/build/$baseFileName.in"
slightIn="/Users/aglaenzer/Softwares/Starlight/build/slight.in"
slightOut="/Users/aglaenzer/Softwares/Starlight/build/$baseFileName.out"
outputFile="/Users/aglaenzer/Softwares/Starlight/build/output-$year-$config.txt"

if test -f $slightIn;
then
    rm $slightIn
fi

dir=$PWD
write
outputFolder=$dir/files/$process
if ! test -d $outputFolder ; then
mkdir $outputFolder
fi


if [ $runStarlight = 1 ]; then
echo 'Writing in' $outputFile
cd /Users/aglaenzer/Softwares/Starlight/build
./starlight >& $outputFile&
echo 'Running Starlight...'
wait
echo 'Done'
cd $dir
mv $slightIn2 "$outputFolder/slight-$year-$config.in"
mv $slightOut "$outputFolder/slight-$year-$config.out"
mv $outputFile "${outputFolder}/output-$year-$config.txt"
#mv $slightIn2 $outputFolder/
#mv $slightOut $outputFolder/
fi

clean

end=`date +%s`

runtime=$((end-start))

echo "Run time = $runtime"
