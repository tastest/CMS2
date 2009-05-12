To run the code in this directory, one will have to have already done eval `scramv1 runtime -sh`
in a CMSSW area (preferably in the same release used to make your ntuples)

Then, export CMS2_LOCATION to wherever your CMS2 directory is.

All the macros you need are in the above directory (printHists, etc)

doQCDFREstimate will run the QCDFRestimator code. This will extract the FR
from QCD and apply it to W+Jets. Do:

root [0] .L ../setup.C
roto [1] .L setup()
root [2] .L doQCDFREstimate.C 
root [3] doAll()

To apply the FR to all the samples and estimate the number of fakes, one needs
to run EstimateFakes.C:

root [0] .L EstimateFakes.C
root [1] EstimateFakes(1893376) //whatever bitmask you want 
