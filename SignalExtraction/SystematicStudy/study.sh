#!/bin/bash

## variables

rapMin=-4
#rapMax=-3.2
rapMax=-2.5

exp=0				# if true, use the exponential for the dissociative contribution. if false, use the power law

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
echo 'Error: please type ./Fit.sh $rapRange'
echo 'with'
echo -e '--> rapRange =\t 1 (-4.0 < y < -2.5)'
echo -e '\t\t 2 (-4.0 < y < -3.25)'
echo -e '\t\t 3 (-3.25 < y < -2.5)'
echo
exit
}


## end of function definitions

if [ -z $1 ]
then
errorMessage
fi

rapMode=$1
clean
root -l -q "FitResultAnalysis.C+($exp, $rapMode)"
clean
exit



