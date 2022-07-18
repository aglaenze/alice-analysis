## This is my Starlight directory for J/Psi photoproduction in p-Pb collisions

## How to use:

- First, generate a Starlight simulation that will run over one specific process
```./runStarlight.sh $year $processNumber $config```
where 

$year can be 2013 or 2016 (it will change the gamma of the colliding ions, they are written to match ALICE data)

$processNumber can be 1 to 4

$config is 1 (p-Pb) or 2 (Pb-p)

In case of doubt, juste type ```./runStarlight.sh``` and it will display the possible options.

- Once Stalight has run, it wrote in a folder files/$process the input and output files. There is 1 input file (something like slight.in) and 2 output files (slight.out and output.txt)

- By executing ```./getXsections.sh $year $processNumber $config```, this prints in the terminal the cross section for the given process AND IN A GIVEN RAPIDITY RANGE AND MASS RANGE (to be specified in getXsections.sh).
