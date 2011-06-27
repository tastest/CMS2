To run the code in this directory, you need to do:

ln -s ../../../CORE CORE
ln -s /data/tmp/cms2 data

doQCDFREstimate will run the QCDFRestimator code. This will extract the FR
from QCD and apply it to W+Jets. Do:

root [2] .L doQCDFREstimate.C 
root [3] doAll()

To apply the FR to all the samples and estimate the number of fakes, one needs
to run EstimateFakes.C:

root [0] .L EstimateFakes.C
root [1] EstimateFakes(1893376) //whatever bitmask you want

The second step is still not ready.....working on porting that part of the
code.
