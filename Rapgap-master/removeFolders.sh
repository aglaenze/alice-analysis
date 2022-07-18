#!/bin/bash

for (( i=30; i<=80; i++ ))
do
testFolder=test${i}
if [ -d $testFolder ]
then
rm -rf $testFolder
echo "$testFolder removed"
fi
done



