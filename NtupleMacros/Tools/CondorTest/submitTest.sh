#! /bin/bash

echo "Submitting 5 test jobs."
echo "Use   condor_q $USER   to check your jobs" 
echo "When your jobs finish you can find their std out/err files in /data/tmp/${USER}/TestJob/std_logs"

./lib_sh/submit.sh -e ./test/test.sh -a 1 -i ./test/test.txt -u "TestJob"
./lib_sh/submit.sh -e ./test/test.sh -a 2 -i ./test/test.txt -u "TestJob"
./lib_sh/submit.sh -e ./test/test.sh -a 3 -i ./test/test.txt -u "TestJob"
./lib_sh/submit.sh -e ./test/test.sh -a 4 -i ./test/test.txt -u "TestJob"
./lib_sh/submit.sh -e ./test/test.sh -a 5 -i ./test/test.txt -u "TestJob"

echo 
echo "Done submitting test jobs"