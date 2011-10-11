#! /bin/bash
## load the configuration file
: ${1?"Error: No configuration file. To skim you must supply a configuration file as argument 1. Exiting."}
if [ ! -e $1 ]; then
	echo "Error: Could not find configuration file $1. Exiting"
fi
. libsh/loadConfig.sh $1
##


## Check that the proper variables are loaded
: ${input_dir:?"Error: var input_dir not set. Check the config file $1. Exiting."}
: ${output_dir:?"Error: var output_dir not set. Check the config file $1. Exiting."}
: ${skim_C:?"Error: skim_C not set. Check the config file $1. Exiting."}
: ${makeSkim_C:?"Error: makeSkim_C not set. Check the config file $1. Exiting."}
: ${cms2_C:?"Error: cms2_C not set. Check the config file $1. Exiting."}
: ${cms2_h:?"Error: cms2_h not set. Check the config file $1. Exiting."}
: ${libminifwlite:?"Error: libminifwlite not set. Check the config file $1. Exiting."}
: ${core_tgz:?"Error: core_tgz not set. Check the config file $1. Exiting."}
: ${dataset_name:?"Error: dataset_name not set. Check the config file $1. Exiting"}
: ${isData:?"Error: isData not set. Check the config file $1. Exiting"}
##


## Check that all the variables point somewhere
if [ ! -d "$input_dir" ]; then
	echo "Error: input_dir $input_dir does not exist. Exiting."
	exit 1
fi
if [[ "$input_dir" != /hadoop/* ]]; then
	echo "Error: input_dir path $input_dir does not start with /hadoop/. Files to be skimmed must be on hadoop. Exiting."
	exit 1
fi
if [[ "$output_dir" != /hadoop/* ]]; then
	echo "Error: output_dir path $output_dir does not start with /hadoop/. Files to be skimmed must end on hadoop. Exiting."
	exit 1
fi
if [ ! -d "$output_dir" ]; then
	echo "Warning: output_dir $output_dir does not exist. Making now."
	mkdir -p $output_dir || exit 1
fi
if [ ! -e "$skim_C" ]; then
	echo "Error:  skim_C code $skim_C does not exist. Exiting."
	exit 1
fi
if [ ! -e "$cms2_C" ]; then
	echo "Error:  cms2_C code $cms2_C does not exist. Exiting."
	exit 1
fi
if [ ! -e "$cms2_h" ]; then
	echo "Error:  cms2_h code $cms2_h does not exist. Exiting."
	exit 1
fi
if [ ! -e "$libminifwlite" ]; then
	echo "Error:  libminifwlite library $libminifwlite does not exist. Exiting."
	exit 1
fi
if [ ! -e "$core_tgz" ]; then
	echo "Error:  core_tgz file $core_tgz does not exist. Exiting."
	exit 1
fi
if [[ $isData != false &&  $isData != true ]]; then
	echo "Error: isData=$isData. This is not a bool. Should be true or false. Exiting."
	exit 1
fi
##

#remember where the submit scripts are and then move into the dir with the ntuples
submit_dir=$PWD
cd $input_dir
##


## figure out where to get the list of files and then submit
source=
if [ "$reskim_list" = "" ]; then
	source='ls *.root'
else
    source='cat $reskim_list'
fi
for ntuple in `eval $source`; do
	$submit_dir/libsh/submit.sh -e $submit_dir/skim/wrapper.sh -a "$input_dir/$ntuple $output_dir ${skim_C##*/} ${libminifwlite##*/} $isData"  -i "$cms2_h,$cms2_C,$core_tgz,$makeSkim_C,$libminifwlite,$skim_C" -u "$dataset_name"
done



