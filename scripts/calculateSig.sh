bk_file=$1
sig_file_1=$2
sig_file_2=$3
count=1

echo "\begin{tabular}{l| c  c  c  c c c}"
echo "\hline"
echo "Var & & \$W'_{400}$ & Significance & AxigluonR & Significance & SM\\\\"
cat $bk_file |cut -d"&" -f7 | awk '{print $1,$3}'  | while read -r fr err; do
    a=` echo $fr `
    b=` cat ${sig_file_1} |cut -d"&" -f7 | awk '{print $1 }' | sed -n "${count}p" `
    c=` cat ${sig_file_2} |cut -d"&" -f7 | awk '{print $1 }' | sed -n "${count}p" `
    berr=` cat ${sig_file_1} |cut -d"&" -f7 | awk '{print $3 }' | sed -n "${count}p" `
    cerr=` cat ${sig_file_2} |cut -d"&" -f7 | awk '{print $3 }' | sed -n "${count}p" `
    bsig=`echo $a $b $berr | awk '{print (($2-$1)/$3 > 0) ? ($2-$1)/$3 : ($1-$2)/$3 }' `
    csig=`echo $a $c $cerr | awk '{print (($2-$1)/$3 > 0) ? ($2-$1)/$3 : ($1-$2)/$3 }' `
    prefix=` cat ${sig_file_1} |cut -d"&" -f1-2  | sed -n "${count}p" `
    printf "%s %s %6.2f %s %6.2f %s %6.1f%s %s %6.2f %s %6.2f %s %6.1f%s %s %6.2f %s %6.2f %s \n" "$prefix" "&" $b "$\pm$" $berr "&" $bsig "$\sigma$" "&" $c "$\pm$" $cerr "&" $csig "$\sigma$" "&" $a "$\pm$" $err "\\\\"  
    (( count++ ))  
done 

echo "\hline"
echo "\end{tabular}"