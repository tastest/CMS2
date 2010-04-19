#!/bin/bash

# Run basic plotting and scanning script
# and transfer plots and dumps to uaf
# Note that because cms-tas0? is not slc5
# we are executing processCands.sh on lxplus5
# and transferring everything back over
cp wcands.txt wcands.txt_previous
root -b -l -q makePlots.C
# check if there is a change in the list of run, lumisection, event (we will also send a mail
# about a new W candidate if a candidate disappers right now).
cat wcands.txt          | awk '{print $4" "$6" "$8}' | grep -v run | grep ^[^" "] | sort > wcands.txt_short
cat wcands.txt_previous | awk '{print $4" "$6" "$8}' | grep -v run | grep ^[^" "] | sort > wcands.txt_previous_short
#
foundmorecandidates=`diff wcands.txt_short wcands.txt_previous_short`
 if [ `echo -n $foundmorecandidates | wc -c` -gt 0 ];
   then
   echo "**********ALERT***ALERT***ALERT************"
   echo "There are new W candidates! sending mail..."
   echo "*******************************************"
# prepare and send ALERT email:
   echo "Check out: http://uaf-2.t2.ucsd.edu/~jribnik/whunt/" > wcand_new.txt 
   echo "" >> wcand_new.txt 
   head -n 3 wcands.txt >> wcand_new.txt 
   diff wcands.txt wcands.txt_previous >> wcand_new.txt
   cat wcand_new.txt | mail -s "[WHunt ALERT!] New Candidate found" ingo.bloch@cern.ch
   cat wcand_new.txt | mail -s "[WHunt ALERT!] New Candidate found" jribnik@cern.ch
   cat wcand_new.txt | mail -s "[WHunt ALERT!] New Candidate found" ayagil@physics.ucsd.edu
   cat wcand_new.txt | mail -s "[WHunt ALERT!] New Candidate found" fgolf@physics.ucsd.edu
   rm wcand_new.txt 
 else
   echo "No new W candidates"
 fi
rm wcands.txt_short
rm wcands.txt_previous_short
#
# Run Dave / Yanyan data MC comparison:
echo "Starting data/MC comparison looper"
cd ../HybridLooper/oldstyle/
make build
# was makeWPlots.py before
python makeWSummaryPlot.py wfinder
#
scp wcands.txt lxplus303:~/scratch0/whunt
ssh lxplus303 /afs/cern.ch/user/j/jribnik/scratch0/whunt/processCands.sh wcands.txt /tas03/disk01/whunt
scp plots/* uaf-4.t2.ucsd.edu:~/public_html/whunt/plots
scp dumps/* uaf-4.t2.ucsd.edu:~/public_html/whunt/dumps
scp picks/* uaf-4.t2.ucsd.edu:~/devel/fireShot/picks
ssh uaf-4.t2.ucsd.edu "touch ~/devel/fireShot/.fireShot.bang"
