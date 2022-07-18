#!/bin/bash

## Variables

mMin=2.5
mMax=3.5
drawRunByRun=false
periods="{\"LHC16r\", \"LHC16s\"}"
periods="{\"LHC16r\"}"

## end of variable

clean() {
filesToDelete="*.d *.so *.pcm *ACLiC* *.in *.out"
for entry in $filesToDelete
do
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}

errorMessage() {
echo "Please type ./ExtractEff i j"
echo "where i can be"
echo -e "\t1 = kIncohJpsiToMu"
echo -e "\t2 = kCohJpsiToMu"
echo -e "\t3 = kTwoGammaToMuLow"
echo -e "\t4 = kTwoGammaToMuMedium"
echo -e "\t5 = kIncohPsi2sToMuPi"
echo -e "\t6 = kIncohPsi2sToMu"
echo "and j can be"
echo -e "\t1 = -4.0 < y < -2.5"
echo -e "\t2 = -4.0 < y < -3.25"
echo -e "\t3 = -3.25 < y < -2.5"
exit
}

if [ -z $2 ]
then
errorMessage
fi

if [ $1 = 1 ]
then
process=kIncohJpsiToMu
mMin=2.5
mMax=3.5
elif [ $1 = 2 ]
then
process=kCohJpsiToMu
mMin=2.5
mMax=3.5
elif [ $1 = 3 ]
then
process=kTwoGammaToMuLow
elif [ $1 = 4 ]
then
process=kTwoGammaToMuMedium
elif [ $1 = 5 ]
then
process=kIncohPsi2sToMuPi
mMin=2.5
mMax=3.5
elif [ $1 = 6 ]
then
process=kIncohPsi2sToMu
mMin=2.5
mMax=4.0
else
errorMessage
fi

if [ $2 = 1 ]
then
rapMin=2.5
rapMax=4.0
elif [ $2 = 2 ]
then
rapMin=3.25
rapMax=4.0
elif [ $2 = 3 ]
then
rapMin=2.5
rapMax=3.25

elif [ $2 = 11 ]
then
rapMin=3.75
rapMax=4.0
elif [ $2 = 12 ]
then
rapMin=3.5
rapMax=3.75

elif [ $2 = 13 ]
then
rapMin=3.25
rapMax=3.5
elif [ $2 = 14 ]
then
rapMin=3.0
rapMax=3.25

elif [ $2 = 15 ]
then
rapMin=2.75
rapMax=3.0
elif [ $2 = 16 ]
then
rapMin=2.5
rapMax=2.75


else
errorMessage
fi

clean
root -l -q "Efficiency.C(\"$process\", $periods, $mMin, $mMax, $rapMin, $rapMax)"
clean
#root -l -q "PlotEfficiency.C+(\"$process\", $periods, $mMin, $mMax, $rapMin, $rapMax, true)"
#root -l -q "CompareGenToLumi.C+(\"$process\", $periods)"
clean
