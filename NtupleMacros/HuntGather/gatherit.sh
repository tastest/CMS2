#!/bin/bash
export ROOTSYS=/code/osgcode/UCSD_root/root_v5.24.00
export PATH=$PATH:$ROOTSYS/bin
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
source /data/vdt/setup.sh

#while [ 1 ] 
#do
# gather it
  cd /nfs-3/userdata/yanjuntu/hunt
  root -q -b gather_doAll.C

# copy plots to web dir
#  rm -f ~jribnik/public_html/gather/tmp/*.png
 # mv -f ~jribnik/public_html/gather/gather_yesterday/*.png  ~jribnik/public_html/gather/tmp
 # mv -f ~jribnik/public_html/gather/gather_yesterday/*total.out ~jribnik/public_html/gather/tmp
 # mv -f ~jribnik/public_html/gather/*.png ~jribnik/public_html/gather/gather_yesterday/.
 # mv -f ~jribnik/public_html/gather/*total.out ~jribnik/public_html/gather/gather_yesterday/.
 # mv *.png ~jribnik/public_html/gather/
 # source eventListMaker.sh
 # mv *total.out ~jribnik/public_html/gather/ 
 # 
 # sleep 86400
#done
