dataset_name=$1
input_dir=$2 
output_dir=$3 
skim_log_dir=$3/skim_log
echo The input directory is $1
echo The output directory is $3
[ ! -d "$skim_log_dir" ] && echo Create $skim_log_dir && mkdir -p  $skim_log_dir 
dateS=`date '+%Y.%m.%d-%H.%M.%S'`
echo "Start skims at "`date`
echo "Done skimming "`date`

for di in ${input_dir}; do
    find ${di} -name merged\*.root | while read -r f; do
	fO=`echo ${f} | sed -e "s?${input_dir}/merged?${output_dir}/skimmed?g" `
	[ ! -f "${fO}" -o "${f}" -nt ${fO} ]   && echo ${fO} && root -l -b -q "makeSkim.C(\"${f}\",\"${fO}_ready\")"
	[ -f "${fO}_ready" ] && echo ${fO}_ready && mv -f ${fO}_ready ${fO}
    done
done   >&  skim_log_dir/Skim.log.${dateS}

echo "Done skimming "`date`