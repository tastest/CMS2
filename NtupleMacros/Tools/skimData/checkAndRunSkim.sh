dataset_name=$1
input_dir=$2 
output_dir=$3 
skim_log_dir=$3/skim_log
pattern=$4
miniFWlib=$5
skim_C=$6

echo The input directory is $2
echo The output directory is $3
[ ! -d "$skim_log_dir" ] && echo Create $skim_log_dir && mkdir -p  $skim_log_dir 
dateS=`date '+%Y.%m.%d-%H.%M.%S'`
echo "Start skims at "`date`

for di in ${input_dir}; do
    ls ${input_dir}/${pattern} 
    find ${di} -name ${pattern} | while read -r f; do
	fO=`echo ${f} | sed -e "s?${input_dir}/?${output_dir}/?g" `
	[ ! -f "${fO}" -o "${f}" -nt ${fO} ]   && echo ${fO} && root -l -b -q "makeSkim.C(\"${f}\",\"${fO}_ready\",\"${skim_C}\",\"${miniFWlib}\")"
	[ -f "${fO}_ready" ] && echo ${fO}_ready && mv -f ${fO}_ready ${fO}
    done
done   >&  ${skim_log_dir}/Skim.log.${dateS}

echo "Done skimming "`date`