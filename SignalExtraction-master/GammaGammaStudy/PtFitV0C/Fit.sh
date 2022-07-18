#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"
periods="{\"LHC16r\"}"

#mMin=2.0
#mMax=2.8

mMin=1.0
mMax=1.5

ptMin=0.000001
#ptMin=1
ptMax=3
exclusiveOnly=1



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
clean
root -l -q "PtPlotV0C.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
#exit
}

errorMessage()
{
echo
echo 'Error: please type ./Fit.sh $rapRange $massRange'
echo 'with'
echo -e '--> rapRange =\t 1 (-4.0 < y < -2.5)'
echo -e '\t\t 2 (-4.0 < y < -3.25)'
echo -e '\t\t 3 (-3.25 < y < -2.5)'
echo 'and'
echo -e '--> massRange =\t 1 (1.0 < M < 1.5)'
echo -e '\t\t 2 (1.5 < M < 2.0)'
echo -e '\t\t 3 (2.0 < M < 2.5)'
echo
exit
}

# tests on the first parameter: the rapidity

defineRanges() {
if [ $rapRange = 1 ]
then
rapMin=-4.0
rapMax=-2.5
elif [ $rapRange = 2 ]
then
rapMin=-4.0
rapMax=-3.25
elif [ $rapRange = 3 ]
then
rapMin=-3.25
rapMax=-2.5
else
errorMessage
fi

# tests on the second parameter: the mass
if [ $massRange = 1 ]
then
mMin=1.0
mMax=1.5
elif [ $massRange = 2 ]
then
mMin=1.5
mMax=2.0
elif [ $massRange = 3 ]
then
mMin=2.0
mMax=2.5
else
errorMessage
fi
}

getParam() {
# initialisation of parameters
if [ $rapRange = 1 ]
then
if [ $massRange = 1 ]
then
mu16r=0.067
mu16rErr=0.003
sigma16r=0.030
sigma16rErr=0.002
elif [ $massRange = 2 ]
then
mu16r=0.077
mu16rErr=0.004
sigma16r=0.035
sigma16rErr=0.003
elif [ $massRange = 3 ]
then
mu16r=0.089
mu16rErr=0.006
sigma16r=0.036
sigma16rErr=0.004
fi

elif [ $rapRange = 2 ]
then
if [ $massRange = 1 ]
then
mu16r=0.070
mu16rErr=0.003
sigma16r=0.032
sigma16rErr=0.002
elif [ $massRange = 2 ]
then
mu16r=0.076
mu16rErr=0.005
sigma16r=0.034
sigma16rErr=0.003
elif [ $massRange = 3 ]
then
mu16r=0.093
mu16rErr=0.008
sigma16r=0.032
sigma16rErr=0.004
fi

elif [ $rapRange = 3 ]
then
if [ $massRange = 1 ]
then
mu16r=0.054
mu16rErr=0.005
sigma16r=0.023
sigma16rErr=0.003
elif [ $massRange = 2 ]
then
mu16r=0.078
mu16rErr=0.008
sigma16r=0.039
sigma16rErr=0.006
elif [ $massRange = 3 ]
then
mu16r=0.083
mu16rErr=0.009
sigma16r=0.038
sigma16rErr=0.006
fi

fi
}

## end of function definitions


muonfilter=std
#muonfilter=nopXdca
#muonfilter=noLpt

#I have a diffrent folder inside "rootFiles/" for each of the muon filter processed in the analysis
#path="~/alice/analysis-alice/p-Pb-2016/rootFiles"
path="/Volumes/Transcend/rootFiles-pPb"
path_to_rootfiles_data=\"${path}/$muonfilter\"
path_to_rootfiles_MC=\"${path}/MC-$muonfilter\"

# It is assumed that the format of the root file names is
# - AnalysisResults_LHC16r_MC_kIncohJpsiToMu.root for MC data
# - AnalysisResults_LHC16r.root for real data

## beginning of action
if [ -z $2 ]
then
errorMessage
fi

echo $muonfilter
echo $periods

rapRange=$1
massRange=$2

defineRanges
getParam

clean
init
run
exit


