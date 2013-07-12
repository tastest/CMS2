dataset_name=$1
input_dir_sub=$2/${dataset_name} 
output_dir_sub=$3/${dataset_name} 
skim_log_dir=${output_dir_sub}/skim_log
pattern=$4
miniFWlib=$5
skim_C=$6

echo The input directory is ${input_dir_sub}
echo The output directory is ${output_dir_sub}
[ ! -d "$skim_log_dir" ] && echo Create $skim_log_dir && mkdir -p  $skim_log_dir 
dateS=`date '+%Y.%m.%d-%H.%M.%S'`
echo "Start skims at "`date`

for di in ${input_dir_sub}; do
    ls ${input_dir_sub}/${pattern} 
    find ${di} -name ${pattern} | while read -r f; do
	fO=`echo ${f} | sed -e "s?${input_dir_sub}?${output_dir_sub}?g" `
	[ ! -f "${fO}" -o "${f}" -nt ${fO} ]   && echo ${fO} && root -l -b -q "makeSkim.C(\"${f}\",\"${fO}_ready\",\"${skim_C}\",\"${miniFWlib}\")"
	[ -f "${fO}_ready" ] && echo ${fO}_ready && mv -f ${fO}_ready ${fO}
    done
done   >&  ${skim_log_dir}/Skim.log.${dateS}

echo "Done skimming "`date`