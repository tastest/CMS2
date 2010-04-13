#!/bin/bash

# Run basic plotting and scanning script
# and transfer plots and dumps to uaf
# Note that because cms-tas0? is not slc5
# we are executing wcands.sh on lxplus5
# and transferring everything back over
cp wcands.txt wcands.txt_previous
root -b -l -q makePlots.C
foundmorecandidates=`diff wcands.txt wcands.txt_previous`
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
scp wcands.txt lxplus303:~/scratch0/whunt
ssh lxplus303 /afs/cern.ch/user/j/jribnik/scratch0/whunt/wcands.sh
scp plots/* uaf-2.t2.ucsd.edu:~/public_html/whunt/plots
scp dumps/* uaf-2.t2.ucsd.edu:~/public_html/whunt/dumps
scp picks/* uaf-2.t2.ucsd.edu:~/public_html/whunt/picks
