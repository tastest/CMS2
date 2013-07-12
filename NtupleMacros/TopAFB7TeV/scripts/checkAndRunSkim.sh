dataset_name=$1
input_dir=$2 
output_dir=$3 
pattern=$4
skim_log_dir=$3/skim_log
echo $1
echo $3
[ ! -d "$skim_log_dir" ] && echo Create $skim_log_dir && mkdir -p  $skim_log_dir 
echo "Start skims at "`date`
dateS=`date '+%Y.%m.%d-%H.%M.%S'`

echo "                             "
echo "Start  TPrime skimming "
echo "                             "
fIn=$input_dir/${pattern}
fOut=`echo ${output_dir}/skimmed_ntuples.root` 
root -l -b -q "makeSkim.C(\"${fIn}\",\"${fOut}_ready\",\"TPrime\")"  >&  $skim_log_dir/Skim.log.${dateS}
[ -f "${fOut}_ready" ] && echo ${fOut}_ready && mv -f ${fOut}_ready ${fOut}

echo "Done skimming "`date`
