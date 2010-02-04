#!/bin/bash
skimInCard=$1
[ "x$1" == "x" ] && echo "Specify the skim card file" && exit 21
[ ! -f "${skimInCard}" ] && echo "Skim input card does not exist " && exit 22
grep [A-Z] ${skimInCard} | while read -r ib ds ob x y rest
 do cL="$ib $ds $ob $x $y $rest"
 [ "x${y}" =="x" ] && "Input line $cL with less than 5 words " && exit 220
 [ "x${rest}" != "x" ] && "Card line $cL longer than 5 words: trailing input ignored "
 [ ! -d "${ib}" ] && echo "Base dir ${ib} does not exist or input in wrong format" && exit 23
 [ ! -d "${ib}/${ds}" ] && echo "Dataset dir ${ib}/${ds} does not exist or input in wrong format" && exit 24
 [ ! -d "${ob}" ] && echo "Output base dir ${ob} does not exist or input in wrong format " && exit 25
 echo ${x} | grep [A-Za-z] >& /dev/null && echo "Bad number format ${x}" && exit 26
 echo ${y} | grep [A-Za-z] >& /dev/null && echo "Bad number format ${x}" && exit 26
done
# now do the job
cat ${skimInCard} | while read -e ib ds ob x y
 do [ ! -d "$ob/$ds" ] && mkdirhier $ob/$ds
 [ -f "$ob/$ds/.lock" ] && echo "Got lock in $ob/$ds : remove to continue" && continue
 touch $ob/$ds/.lock
 find $ib/$ds -name merged\*.root | grep merged |while read -r ln
   do sn=`echo $ln | sed -e "s?$ib/$ds/??g"`
   rmd=" skims/skimDilPtXY.C(\"$ib/$ds\",\"$ob/$ds\",\"$sn\",20,10)"
   ofl=$ob/$ds/${sn}_2010.log
   echo "At "`date`" : $rmd" > $ofl
   root -l -b -q $rmd >> $ofl 2>&1 
   echo done at `date`>> $ofl
 done
done 

