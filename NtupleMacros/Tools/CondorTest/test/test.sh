#! /bin/bash
echo "Arg1     = $1"
echo "Hostname = $HOSTNAME"
echo "Date     = `date`"
echo "`ls -l .`"
cat test.txt
echo "Who am I = `whoami`"
echo `env`

sleep 5m
