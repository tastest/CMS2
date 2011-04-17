#!/bin/bash

#
# This script makes the plots
# it should be run as a cron job
#

cd /tas03/home/dlevans/express_gather/CMSSW_4_1_2/src/CMS2/NtupleMacros/HuntGather2011/plotter
export SCRAM_ARCH=slc5_amd64_gcc434
eval `scramv1 runtime -sh`

echo root -q -b makeGatherPlots.C\(\"/tas/cms2/gather/\",\EXPRESS\)
scp /tas03/home/dlevans/express_gather/CMSSW_4_1_2/src/CMS2/NtupleMacros/HuntGather2011/output/*.png uaf-6.t2.ucsd.edu:/home/users/dlevans/public_html/gather2011/express_gather/
ssh uaf-6.t2.ucsd.edu -q "touch /home/users/dlevans/public_html/gather2011/express_gather/index.php"

