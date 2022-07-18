#!/bin/bash


clean() {
filesToDelete="*.so *.d *.pcm *ACLiC* *.in *.out"
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
echo 'Error: please type ./Fit.sh $rapMode $massMode'
echo 'with'
echo -e '--> rapMode =\t 1 (-4.0 < y < -2.5)'
echo -e '\t\t 2 (-4.0 < y < -3.25)'
echo -e '\t\t 3 (-3.25 < y < -2.5)'
echo -e 'and massMode =\t 1 (1.0 < M < 1.5)'
echo -e '\t\t 2 (1.5 < M < 2.0)'
echo -e '\t\t 3 (2.0 < M < 2.5)'
echo
exit
}


## end of function definitions

if [ -z $2 ]
then
errorMessage
fi

rapMode=$1
massMode=$2
clean
root -l -q "FitResultAnalysis.C($rapMode, $massMode)"
clean
exit



