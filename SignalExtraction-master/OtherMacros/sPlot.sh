#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"

#mMin=2.0
#mMax=2.8

mMin=2.5
mMax=3.5

ptMin=0
ptMax=5
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
initiateParameters

# then LHC16s
inputfile="input-LHC16s.txt"
bExc=$bExc16s
gammaPbYield=$gammaPbYield16s
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


# tests on the first parameter: the year
if [ $1 = 1 ]
then
rapMin=-4.0
rapMax=-2.5
gammaPbYield16r=49
gammaPbYield16s=78.5
bExc16r=3.80
bExc16s=5.52
elif [ $1 = 2 ]
then
rapMin=-4.0
rapMax=-3.25
gammaPbYield16r=23
gammaPbYield16s=43.4
bExc16r=3.49
bExc16s=5.52
elif [ $1 = 3 ]
then
rapMin=-3.25
rapMax=-2.5
gammaPbYield16r=27
gammaPbYield16s=39.6
bExc16r=3.66
bExc16s=5.52
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
clean
root -l -q "SplotRap.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
exit



rm *txt
exp=0


rapMin=-4
rapMax=-3.25
gammaPbYield16r=19.1
gammaPbYield16s=28.9
init
clean
root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
bExc16r=$(echo "scale=2; 2.8+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.2+${k}*0.04" | bc)
init
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
done

rapMin=-3.25
rapMax=-2.5
gammaPbYield16r=48.4
gammaPbYield16s=32.5
init
clean
root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
bExc16r=$(echo "scale=2; 2.8+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.2+${k}*0.04" | bc)
init
clean
root -l -q "TwoDPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean
done

rapMin=-4
rapMax=-2.5
gammaPbYield16r=48.4
gammaPbYield16s=75.1
init
clean
root -l -q "Splot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods)"
clean

for ((k=0;k<70;k++)); do
clean

bExc16r=$(echo "scale=2; 2.8+${k}*0.04" | bc)
bExc16s=$(echo "scale=2; 4.2+${k}*0.04" | bc)
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



