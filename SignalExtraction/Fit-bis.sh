#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"
periods="{\"LHC16r\"}"

#mMin=2.0
#mMax=2.8

mMin=2.5
mMax=3.5

ptMin=0
ptMax=3
useCuts=1
logScale=0
drawPulls=0			# draws graphs data-fit, in any case the graphs (data-fit)/sigma are plotted
exp=0				# if true, use the exponential for the dissociative contribution. if false, use the power law
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
#root -l -q "SplotZ.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
#root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
#exit
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
gammaPbYield16r=56
gammaPbYield16s=84
bExc16r=3.62
bExc16s=5.49
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi

if [ $exclusiveOnly = 1 ]
then
binc16s=2.43
ninc16s=6.19
else
binc16s=1.71
ninc16s=7.00
fi

elif [ $1 = 2 ]
then
rapMin=-4.0
rapMax=-3.25
gammaPbYield16r=25
gammaPbYield16s=43
bExc16r=3.38
bExc16s=5.15
mu16r=0.0832
sigma16r=0.0366
mu16s=0.0708
sigma16s=0.0059

if [ $exclusiveOnly = 1 ]
then
binc16r=3.88
ninc16r=7.0
else
binc16r=8.00
ninc16r=1.00
fi

if [ $exclusiveOnly = 1 ]
then
binc16s=1.78
ninc16s=7.0
else
binc16s=2.3452
ninc16s=8
fi

elif [ $1 = 3 ]
then
rapMin=-3.25
rapMax=-2.5
gammaPbYield16r=33
gammaPbYield16s=43
bExc16r=3.86
bExc16s=5.78
mu16r=0.0752
sigma16r=0.0315
mu16s=0.0612
sigma16s=0.0297

if [ $exclusiveOnly = 1 ]
then
binc16r=3.03
ninc16r=6.90
else
binc16r=2.16
ninc16r=7.00
fi

if [ $exclusiveOnly = 1 ]
then
binc16s=4.59
ninc16s=2.96
else
binc16s=1.4270
ninc16s=4.7255
fi

elif [ $1 = 11 ]
then
rapMin=-4.0
rapMax=-3.75
gammaPbYield16r=5
gammaPbYield16s=3.3
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


elif [ $1 = 12 ]
then
rapMin=-3.75
rapMax=-3.5
gammaPbYield16r=5
gammaPbYield16s=12.3
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


elif [ $1 = 13 ]
then
rapMin=-3.5
rapMax=-3.25
gammaPbYield16r=5
gammaPbYield16s=23.9
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


if [ $exclusiveOnly = 1 ]
then
binc16s=4.59
ninc16s=2.96
else
binc16s=1.4270
ninc16s=4.7255
fi

elif [ $1 = 14 ]
then
rapMin=-3.25
rapMax=-3.0
gammaPbYield16r=5
gammaPbYield16s=26.9
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


elif [ $1 = 15 ]
then
rapMin=-3.0
rapMax=-2.75
gammaPbYield16r=5
gammaPbYield16s=17.0
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


elif [ $1 = 16 ]
then
rapMin=-2.75
rapMax=-2.5
gammaPbYield16r=5
gammaPbYield16s=5.1
bExc16r=3
bExc16s=5
mu16r=0.0799
sigma16r=0.0351
mu16s=0.0660
sigma16s=0.0245

if [ $exclusiveOnly = 1 ]
then
binc16r=4.02
ninc16r=4.20

else
binc16r=5.92
ninc16r=1.77
fi


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



rm *txt


rapMin=-4.0
rapMax=-3.25
gammaPbYield16r=25
gammaPbYield16s=43.4
bExc16r=3.38
bExc16s=4.51
mu16r=0.0832
sigma16r=0.0366
init
clean
#root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
bExc16r=$(echo "scale=2; 2.4+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.0+${k}*0.04" | bc)
init
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
done

rapMin=-3.25
rapMax=-2.5
gammaPbYield16r=33
gammaPbYield16s=39.6
bExc16r=3.86
bExc16s=5.36
mu16r=0.0752
sigma16r=0.0315
init
clean
#root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
bExc16r=$(echo "scale=2; 2.4+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.0+${k}*0.04" | bc)
init
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
done

rapMin=-4.0
rapMax=-2.5
gammaPbYield16r=56
gammaPbYield16s=78.5
bExc16r=3.62
bExc16s=5.72
mu16r=0.0799
sigma16r=0.0351
init
clean
#root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
clean

bExc16r=$(echo "scale=2; 2.4+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.0+${k}*0.04" | bc)
init
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
done


exit




clean
#root -l -q "GetTemplates.C+($path_to_rootfiles_data, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
#exit


root -l -q "NakedPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax, $useCuts)"
clean


clean
mMin=1.3
mMax=2.3
root -l -q "Background.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
exit



