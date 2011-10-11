#! /bin/bash
## check the arguments
: ${1?"Error: Supply the configuration file as the first argument. Exiting"}
if [ ! -e $1 ]; then
	echo "Error: config file $1 does not exist. Exiting"
	exit 1
fi
time=${2-0}

## load the config file
. libsh/loadConfig.sh $1
: ${input_dir:?"Error: input_dir not set, check config file $1. Exiting"}
: ${output_dir:?"Error: output_dir not set, check config file $1. Exiting"}
: ${reskim_list:?"Error: reskim_list not set, check config file $1. Exiting"}

## remove the reskim_list file
[ -e ${output_dir}${reskim_list} ] && rm ${output_dir}${reskim_list}

## decide on how to look for files
command=
if [ $time -gt 0 ]; then
	command="find $input_dir -maxdepth 1 -mtime -$time -name '*.root'"
else
	command="find $input_dir -maxdepth 1 -name '*.root'"
fi
files=`eval $command`

## find the files that need to be resubmitted
resubmit_files=
for file in $files; do
	short_file=${file##*/}
	skim_file=${short_file%%.root}_skim.root
	if [ ! `find $output_dir -maxdepth 1 -name $skim_file` ]; then
		resubmit_files="${resubmit_files}${short_file} "
		#echo $resubmit_files
	fi
done

## make the list of resubmit files
printf "%s\n" $resubmit_files > ${input_dir}${reskim_list}

## tell the user how many files to resubmit and where to find this list
sleep 2 #hadoop is slow, wait a 2 seconds before doing anything with the file above
echo "`cat ${input_dir}${reskim_list} | wc -w` files need to be submitted. Submit list can be found at ${input_dir}${reskim_list}"