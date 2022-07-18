#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"
periods="{\"LHC16r\"}"

#mMin=2.0
#mMax=2.8

mMin=2.0
mMax=2.5

ptMin=0
ptMax=3


## end of variables

## beginning of function definitions

initiateParameters() {
if test -f $inputfile
then
rm $inputfile
fi

touch $inputfile
echo "
bExc = $bExc
gammaPbYield = $gammaPbYield

mMin = $mMin
mMax = $mMax

ptMin = $ptMin
ptMax = $ptMax

rapMin = $rapMin
rapMax = $rapMax

muLandau = $mu
sigmaLandau = $sigma

pt0 = $pt0
nInc = $nInc

useCuts = $useCuts
logScale = $logScale
drawPulls = $drawPulls		# draws graphs data-fit, in any case the graphs (data-fit)/sigma are plotted
exp = $exp	# if true, use the exponential for the dissociative contribution. if false, use the power law
exclusiveOnly = $exclusiveOnly
" >> $inputfile
}

init() {
# first initiate LHC16r
inputfile="input-LHC16r.txt"
bExc=$bExc16r
gammaPbYield=$gammaPbYield16r
mu=$mu16r
sigma=$sigma16r
pt0=$binc16r
nInc=$ninc16r
initiateParameters

# then LHC16s
inputfile="input-LHC16s.txt"
bExc=$bExc16s
gammaPbYield=$gammaPbYield16s
mu=$mu16s
sigma=$sigma16s
pt0=$binc16s
nInc=$ninc16s
initiateParameters
}


clean() {
filesToDelete="*.so *.d *.pcm *ACLiC* *.in *.out Include/*.so Include/*.d Include/*.pcm"
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}

run() {
clean
init
root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
}

errorMessage()
{
echo
echo 'Error: please type ./Fit.sh $rapRange'
echo 'with'
echo -e '--> rapRange =\t 1 (-4.0 < y < -2.5)'
echo -e '\t\t 2 (-4.0 < y < -3.25)'
echo -e '\t\t 3 (-3.25 < y < -2.5)'
echo
exit
}

if [ -z $1 ]
then
errorMessage
fi


# tests on the first parameter: the rapidity
if [ $1 = 1 ]
then
rapMin=-4.0
rapMax=-2.5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245


elif [ $1 = 2 ]
then
rapMin=-4.0
rapMax=-3.25
mu16r=0.0832
sigma16r=0.0366
mu16s=0.0708
sigma16s=0.0059


elif [ $1 = 3 ]
then
rapMin=-3.25
rapMax=-2.5
mu16r=0.0752
sigma16r=0.0315
mu16s=0.0612
sigma16s=0.0297

else
errorMessage
fi


## end of function definitions


muonfilter=std
#muonfilter=nopXdca
#muonfilter=noLpt

#I have a diffrent folder inside "rootFiles/" for each of the muon filter processed in the analysis
path="~/alice/analysis-alice/p-Pb-2016/rootFiles"
path_to_rootfiles_data=\"${path}/$muonfilter\"
path_to_rootfiles_MC=\"${path}/MC-$muonfilter\"

# It is assumed that the format of the root file names is
# - AnalysisResults_LHC16r_MC_kIncohJpsiToMu.root for MC data
# - AnalysisResults_LHC16r.root for real data

echo $muonfilter
echo $periods


## beginning of action
init
run
exit

