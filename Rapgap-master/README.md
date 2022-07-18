

The idea is to use the rapgap software, installed on lxplus by Hannes Jung, in order to simulate dissociative J/Psi photoproduction off protons. The full code is meant to run on lxplus.

# How to use

## Preliminary

1) Connect to lxplus
```ssh -Y aglaenze@lxplus.cern.ch ```

2) Go to your rapgap working directory

3) Get OCDB snapshots matching the run you need, and put them in the working directory. I am using run 266025, so these files: 

- /alice/sim/2017/LHC17e4/OCDB/266025/OCDBsim.root

- /alice/sim/2017/LHC17e4/OCDB/266025/OCDBrec.root

## Then

4) Load ALICE environment

```alias loadenv="alienv enter VO_ALICE@AliPhysics::v5-09-05-01-1,VO_ALICE@AliDPG::prod-201811-01-1" ; ```
```loadenv ```

5) Initiate a new test folder

``` ./newTest $k ```

6) Go in the new test directory

``` cd test$k ```

7) If need be, modify the file steer_test

## To run "locally" on lxplus

``` ./rapgap.sh ```


Optional: ``` screen -S $name ```


``` ./execute.sh $numberOfEvents ```

## With Condor and parallel jobs

For info on Condor: https://batchdocs.web.cern.ch/tutorial/introduction.html

First, if need be, modify in launchJobs.sh the number of events and / or the number of parallel jobs


``` ./launchJobs.sh ```

When jobs are finished, the AliAOD.root files are merged by doing the following:

``` /cvmfs/alice.cern.ch/bin/alienv enter AliPhysics/vAN-20181121_ROOT6-1 ```


``` root -l -q MergeAODs.C  ``` 
