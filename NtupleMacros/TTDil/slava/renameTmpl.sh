finN=$1
[ "$finN" = "" ] && finN=GIVEMENAME
tmpl=CMS2_SW2_looper 
cp $tmpl.h $finN.h
cp $tmpl.C $finN.C
echo Making $finN.h and $finN.C from $tmpl.h and $tmpl.C
for f in `ls $finN.{h,C}`; do cat $f | sed -e "s/${tmpl}/${finN}/g" > a.tmp && mv a.tmp $f; done
echo done
