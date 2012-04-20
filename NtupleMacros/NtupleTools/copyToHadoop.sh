#! /bin/bash

: ${1?"No hadoop directory specified in which to copy the ntuples. Exiting."}
if [ ! -d "$1" ]; then
	echo "Directory, $1, does not exist. Attempting to make."
	mkdir -p "$1"
	if [ ! -d "$1" ]; then
		echo "Failed to make directory, $1. Exiting"
		exit 1
	fi
fi

if [ ! -r "$1" ]; then
	echo "User does not have read permisions on directory, $1. Exiting."
	ls -ld $1
	exit 1
fi

if [ ! -w "$1" ]; then
	echo "User does not have write permisions on directory, $1. Exiting."
	ls -ld $1
	exit 1
fi

if [ ! -x "$1" ]; then
	echo "User does not have execute permisions on directory, $1. Exiting."
	ls -ld $1
	exit 1
fi

counter=0
errorcounter=0
for rootfile in `ls merge*.root`; do
	echo "Copying file $rootfile to dir $1."
	hadoop fs -copyFromLocal $rootfile ${1#/hadoop}
	if [ $? != 0 ]; then
		echo " => Error copying file $rootfile to dir ${1}!"
		errorcounter=$((errorcounter+1))
	fi
	counter=$((counter+1))
done

num=`ls -1 merge*.root | wc -l`
copied=`ls -1 ${1}/merge*.root | wc -l`

echo
echo
echo "Done copying files."
echo "Summary:  Number to copy    $num files."
echo "          Attempted to copy $counter files."
echo "          Copied            $copied files."
echo "          Error copying     $errorcounter files."