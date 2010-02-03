************Old instructions*****************
To run the code in this directory, one will have to have already done eval `scramv1 runtime -sh`
in a CMSSW area (preferably in the same release used to make your ntuples)

Then, export CMS2_LOCATION to wherever your CMS2 directory is.

doQCDFREstimate will run the QCDFRestimator code. This will extract the FR
from QCD and apply it to W+Jets. Do:

root [0] .L doQCDFREstimate.C 
root [1] doAll()

To apply the FR to all the samples and estimate the number of fakes, one needs
to run EstimateFakes.C:

root [0] .L EstimateFakes.C
root [1] EstimateFakes(1893376) //whatever bitmask you want 


********************New instructions*********************
I tried to run this code 1 year later (2/1/10) on uaf-6.t2.ucsd.edu (since 
this is where the ntuples are) but ran into several problems (possibly due 
to different root versions, slc5 vs slc4, different libFWLite because I'm 
trying to use CMSSW_3_X whatever. Anyway, to get it to work after you've 
checked out the code from CVS, do:

export SCRAM_ARCH=slc5_ia32_gcc434
export CMS2_LOCATION=<wherever the directory that contains CMS2 is>/CMS2
ln -s /data/tmp/cms2-V01-02-06/
use: /code/osgcode/cmssoft/cms/slc5_ia32_gcc434/lcg/root/5.22.00d-cms6/bin/root

The rest of the steps are the same as the old instructions. The major differences are that I've added the utility files to this directory. I got a lot of errors while compiling from the directory above. Probably could have figured out how to get it to work, but this was easiest :)
