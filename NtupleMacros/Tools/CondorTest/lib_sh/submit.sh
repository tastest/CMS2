#! /bin/bash
#? ##########################################################
#? SCRIPT: submit.sh - a script to create condor submit files
#? AUTHOR: Ian MacNeill
#? DESCRIPTION: For use on uaf machines. Creates a condor
#?              submit file to submit jobs to condor batch
#?              system. Supports glideinwms submissions and 
#?              condor g direct submissions. Does all the 
#?              work in the current directory.
#? ##########################################################


#TO DO
#1) separate and log std out and err



## FUNCTIONS
die(){
#@ DESCRIPTION: print error message and exit with supplied return code
#@ USAGE: die STATUS [MESSAGE]
	error=$1
	shift
	[ -n "$*" ] && printf "\n%s\n" "$*" >&2
	exit "$error"
}

usage(){
#@ DESCRIPTION: print usage information
#@ USAGE: usage	
	printf "\n*-------------------------------------- USAGE for %s ----------------------------------------*\n" "$scriptname" 
	printf "%s\n" "$usage" 
	printf "\n\n\n" 
}

## INITIALIZE VARIABLES
scriptname=${0##*/}
usage="    -h Displays help
    -v Sets verbose mode (not implemented yet)
    -d Use direct condor submission to UCSD from UAF instead of glideinWMS.
    -t Test mode. Does everything but submit.
    -e [executable]
          The executable to run on the worker node.
    -a [\"args\"]
          List of arguments for the executable. Must be separated by spaces and the whole list
          enclosed in double quotes.
    -i [\"input files\"]
          Comma separated list of input files with full path to be transfered for job.
    -s [\"desired sites\"]
          Defaults to UCSD.
    -c [\"config file\"]
          Unused option for now.
    -u [\"unique_identifier\"] 
          This will be the name of the dir where your output is stored. Defaults to \"UserJob\".

    REMEMBER: 
     * Please use full paths to all files that you pass as arguments to this script. You may have
       weird behavior if you don't.
     * If you use -d for direct submission, you must be at UCSD and have given your DN to Terrence
       to run on the tier 3. Otherwise you'll just get t2 nodes."

help=0
verbose=0
test=0
use_glidein=1
unique_identifier="UserJob"
executable=
desired_sites="\"UCSD\""
arguments=
input_files=
universe="vanilla"


# if the user set up their cmssw environment, X509_USER_PROXY should be defined, if not, look for it
[ -n "${X509_USER_PROXY+x}" ] && user_proxy=$X509_USER_PROXY || user_proxy=`find /tmp/x509up_u* -user $USER`
: ${user_proxy:?"Error: Could not find a user proxy for the user. Exiting."} #check to see if the variable is set but empty, if it is, a proxy couldn't be found, quit

## GET INPUT OPTIONS
#options executable, args, input files, desired sites, verbose?
while getopts hvdte:a:i:s:c:u: opt; do
	case $opt in
		h) help=1;;
		v) verbose=1;;
		d) use_glidein=0;;
		t) test=1;;
		e) executable=$OPTARG;;
		a) arguments="$OPTARG";;
		i) input_files="$OPTARG";;
		s) desired_sites="$OPTARG";;
		c) config_file="$OPTARG";;
		u) unique_identifier="$OPTARG";;
		*) printf "Invalid option(s)\n" && exit 1;;
	esac
done

shift "$(( OPTIND - 1 ))"


## DISPLAY HELP IF ASKED FOR
if [ "$help" = 1 ]; then # if they ask for the help, display the usage and quit
	usage
	exit 0
fi

## CHECK FOR EXISTENCE OF Executable
if [ ! -f "$executable" ]; then
	die 1 "ERROR: executable \"$executable\" doesn't exist."

fi

## CHECK FOR EXISTENCE OF Files
file_exists=1
OLDIFS=$IFS
IFS=$','
for file in $input_files; do
	if [ -n $file ] && [ ! -f $file ]; then		
		printf "\nERROR: file \"%s\" doesn't exist.\n" $file  >&2
		file_exists=0
	fi
done
IFS=$OLDIFS

if [ "$file_exists" = 0 ]; then
	die 1 "ERROR: Exiting do to previous non-existent file(s)."
fi

##if the user wants direct submission, make sure they're at ucsd
if [ $use_glidein = 0 ]; then
	if [[ $HOSTNAME != *uaf*ucsd* ]]; then
		die 1 "ERROR: Attempting to make a direct submission from a machine($HOSTNAME) other than uaf. Don't use -d option."
	else
		universe=grid
	fi
fi

if [ -n "$config_file" ] && [ ! -f "$config_file" ]; then
	die 1 "ERROR: Config file $config_file doesn't exist."
fi


## MORE VARIABLES BASED ON INPUT
unique_identifier="${unique_identifier:-condor_job}"
log_dir=/data/tmp/${USER}/${unique_identifier}/submit_logs
submit_log=$log_dir/submit_`date "+%m_%d_%Y"`.log
condor_log=$log_dir/condor_`date "+%m_%d_%Y"`.log
std_dir=/data/tmp/${USER}/${unique_identifier}/std_logs
submit_file=/data/tmp/${USER}/${unique_identifier}/${unique_identifier}.cmd


## MAKE SOME DIRECTORIES IF THEY DON't EXIST
[ -d /data/tmp/${USER} ] || mkdir -p /data/tmp/${USER}
[ -d $log_dir ] || mkdir -p $log_dir
[ -d $std_dir ] || mkdir -p $std_dir


## WRITE SUBMIT FILE

printf "universe=$universe\n" > $submit_file
if [ $use_glidein = 0 ]; then
	printf "Grid_Resource=gt2 osg-gw-4.t2.ucsd.edu:2119/jobmanager-condor\n" >> $submit_file
fi
cat >> $submit_file <<@EOF
executable=$executable
arguments=$arguments
transfer_executable=True
when_to_transfer_output = ON_EXIT
#the actual executable to run is not transfered by its name.
#In fact, some sites may do weird things like renaming it and such.
transfer_input_files=$input_files
+DESIRED_Sites=$desired_sites 
+Owner = undefined 
log=$condor_log
output=${std_dir}/1e.\$(Cluster).\$(Process).out
error =${std_dir}/1e.\$(Cluster).\$(Process).err
notification=Never
#x509userproxy=$ENV(X509_USER_PROXY)	
x509userproxy=$user_proxy
queue
	
@EOF


## PRINT INFO FOR LOGGING
line="-------------------------------------"
printf "\n\n" 
printf "%s\n" $line 
printf "Dir:  %s\n" $unique_identifier 
printf "Exe:  %s\n" $executable 
printf "Arg:  %s\n" $arguments 
printf "In :  %s\n" $input_files
printf "Site: %s\n" $desired_sites 



## SUBMIT THE JOB
error_code=0
if [ $test = 0 ]; then
	#cp $tmp_submit_file ./$submit_file
	condor_submit $submit_file
	error_code=$?
else
	printf "In test mode, will not submit %s\n" $submit_file
fi

printf "%s\n" $line 


if [ "$error_code" != 0 ]; then
	printf "Error submitting condor job. Exit Code %s.\n" $error_code >&2
fi


exit $error_code
#### figure out how to make this submit log work
#>& $submitting_log

