#!/bin/bash

delete() {
for entry in $filesToDelete
do
#echo $entry
if test -f $entry;
then
#echo Deleting $entry
rm $entry
fi
done
}

filesToDelete="*.xml *.pcm *.d *.so *LHC16*Task* *LHC17*Task* outputs_valid *merge* stderr *validation* myAnalysis.C temp* *plugin*"
delete

#alien.py rm -f *
alien.py rm -r -f LHC16r_ana/
alien.py rm -r -f LHC16s_ana/
alien.py rm -r -f LHC16r_MC*/
alien.py rm -r -f LHC16s_MC*/
alien.py rm -r -f myWorkingDir/
alien.py rm -r -f recycle/
alien.py rm -f Stage_1
alien.py rm -f .plugin_test_copy~
echo " "
alien.py ls -a
y=$(eval "alien.py ls -a | wc -l")
if [ $y = 0 ]
#if [ -z 'alien.py ls -a' ]
then
echo "all good - everything was deleted in alien"
echo " "
fi
