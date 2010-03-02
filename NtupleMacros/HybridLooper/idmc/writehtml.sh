#!/bin/bash

#
# Write html
#

VER=$1
TAG=$2
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
	echo "
        <h2>All Selections (EE): Signal (left), Background (right) <br>
                - All selections w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_cand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_cand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_cand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_cand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

	<h2>Isolation (EE): Signal (left), Background (right) <br>
        	- Isolation w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_isocand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_isocand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_isocand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_isocand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>ID (EE): Signal (left), Background (right) <br>
                - ID w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV (EE): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_convcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_convcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_convcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_convcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV W.R.T Iso+ID (EE): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. Iso+ID</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconvcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconvcand01_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconvcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconvcand01_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>


        <h2>All Selections (EB): Signal (left), Background (right) <br>
                - All selections w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_cand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_cand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_cand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_cand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>Isolation (EB): Signal (left), Background (right) <br>
                - Isolation w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_isocand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_isocand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_isocand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_isocand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>ID (EB): Signal (left), Background (right) <br>
                - ID w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV (EB): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_convcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_convcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_convcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_convcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV W.R.T Iso+ID (EB): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. Iso+ID</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconvcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconvcand01_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconvcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconvcand01_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
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
                	<img src="results/$VER\_s_h1_hyp_debug_after_cand01_pdgid_$DET.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
               	 	<img src="results/$VER\_b_h1_hyp_debug_after_cand01_pdgid_$DET.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                	"
		fi

	fi
#	if [ $TYPE == "eff" ]; then
#            EFF=`echo $FILE_S | sed 's/.*after_\(.*\).png/\1/g'`
#
#	        if [ $VAR == "pt" ] || [ $VAR == "eta" ]; then
#                	echo "
#                	<h2>Eff($VAR, $DET): Signal (left), Background (right) <br>
#                        	- Efficiency as a function of $VAR to pass all selections ($EFF)</h2>
#                	<img src=$FILE_S HEIGHT=$HEIGHT WIDTH=$WIDTH>
#                	<img src=$FILE_B HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
#                	"
#		fi
#	fi

done

echo "
	</body>
</html>
"

