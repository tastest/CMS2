#!/bin/bash

#
# This script makes the plots
# it should be run as a cron job
#

cd /home/users/dlevans/gather/production/CMSSW_4_1_2_patch1/src/CMS2/NtupleMacros/HuntGather2011/plotter
export SCRAM_ARCH=slc5_amd64_gcc434
source /code/osgcode/cmssoft/cms/cmsset_default.sh
eval `scramv1 runtime -sh`

root -q -b makeGatherPlots.C\(\"/nfs-3/userdata/cms2/gather/\",\EXPRESS\)
cp /home/users/dlevans/gather/production/CMSSW_4_1_2_patch1/src/CMS2/NtupleMacros/HuntGather2011/output/post_techstop*.png /home/users/dlevans/public_html/gather2011/express_gather/
touch /home/users/dlevans/public_html/gather2011/express_gather/index.php

