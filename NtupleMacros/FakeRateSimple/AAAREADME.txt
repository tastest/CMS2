export CMS2_LOCATION to wherever your CMS2 directory is.


doQCDFREstimate will run the QCDFRestimator code. This will extract the FR
from QCD and apply it to W+Jets. Do:

root [0] .L doQCDFREstimate.C
root [1] doAll()

To apply the FR to all the samples and estimate the number of fakes, one needs
to run EstimateFakes.C:

