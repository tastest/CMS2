#!/bin/bash

#
# Write html
#

TAG=$1
HEIGHT=320
WIDTH=320

echo "
<html>
	<head>
	<title>
	Validation results for $TAG
	</title>
	</head>

	<body>	
	<h1>Semi-Automated validation results for $TAG</h1>
	<br>
	<br>
"

for FILE_S in `ls results/*.png | grep $TAG | grep -v _b_`; do

        FILE_B=`echo $FILE_S | sed 's/\(.*\)_s_\(.*\)/\1_b_\2/g'`
	TYPE=`echo $FILE_S | sed 's/.*_\(.*\)_s_.*/\1/g'`
	VAR=`echo $FILE_S | sed 's/.*_\(.*\)_E.*\.png/\1/g'`
	DET=`echo $FILE_S | sed 's/.*_\(.*\)\.png/\1/g'`

	if [ $TYPE == "overlay" ]; then
		echo "
		<h2>$VAR ($DET): Signal (left), Background (right) <br>
			- before all selections (solid), after all selections (points)</h2>
		<img src=$FILE_S HEIGHT=$HEIGHT WIDTH=$WIDTH>
	        <img src=$FILE_B HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
		"
		if [  $VAR == "pdgid" ]; then
                echo "
                	<h2>$VAR ($DET): Signal (left), Background (right) <br>
                        	- after all selections (points)</h2>
                	<img src="results/val01_s_h1_hyp_debug_after_cand01_pdgid_$DET.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
               	 	<img src="results/val01_b_h1_hyp_debug_after_cand01_pdgid_$DET.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                	"
		fi

	fi
	if [ $TYPE == "eff" ]; then
	        if [ $VAR == "pt" ] || [ $VAR == "eta" ]; then
                	echo "
                	<h2>Eff($VAR, $DET): Signal (left), Background (right) <br>
                        	- Efficiency as a function of $VAR to pass all selections</h2>
                	<img src=$FILE_S HEIGHT=$HEIGHT WIDTH=$WIDTH>
                	<img src=$FILE_B HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                	"
		fi
	fi

done

echo "
	</body>
</html>
"

