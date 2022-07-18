#!/bin/bash

## variables
periods="{\"LHC16r\", \"LHC16s\"}"
#periods="{\"LHC16r\"}"

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
exclusiveOnly=0





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

muonfilter=std
#muonfilter=nopXdca
#muonfilter=noLpt

#I have a diffrent folder inside "rootFiles/" for each of the muon filter processed in the analysis
path="~/alice/analysis-alice/p-Pb-2016/rootFiles"
path_to_rootfiles_data=\"${path}/$muonfilter\"
path_to_rootfiles_MC=\"${path}/MC-$muonfilter\"



echo $muonfilter
echo $periods

clean
root -l -q "TailParameters.C+($path_to_rootfiles_MC, $periods, $logScale)"
clean
exit


clean
root -l -q "GetTemplates.C+($path_to_rootfiles_data, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
exit


root -l -q "NakedPlot.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax, $useCuts)"
clean


clean
mMin=1.3
mMax=2.3
root -l -q "Background.C+($path_to_rootfiles_data, $path_to_rootfiles_MC, $periods, $mMin, $mMax, $ptMin, $ptMax)"
clean
exit



