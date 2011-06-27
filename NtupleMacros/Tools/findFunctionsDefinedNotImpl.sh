#!/bin/bash

if [ $# -ne 2 ]; then
    echo "USAGE: ./findFunctionsDefinedNotImpl.sh headerfile.h implementation.cc"
    exit 1;
fi

HEADER=$1
IMPL=$2

echo "The following functions are defined in " $HEADER
echo "... but not implemented in " $IMPL

ARGS=`cat $HEADER | grep '(' | awk -F '(' '{print $1}' | awk '{print $2}'`
for ARG in $ARGS; do
    RESULTS=`grep $ARG $IMPL`
    if [ "$RESULTS" == "" ]; then
        echo $ARG
    fi
done

exit 0

