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
	Validation results for $TAG $VER
	</title>
	</head>

	<body>	
	<h1>Semi-Automated validation results for $TAG $VER</h1>
	<br>
	<br>
"
	echo "
    <h1> Definitions </h1>
        In this validation:<br>
        <b>Signal:</b><br>
    chain_ttbar->Add(ntuple_location + "/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");<br>
    select EMU events, where muon and electron are truth matched and have W as mother<br>      

        <b>Background:</b><br>
    chain_wmunu->Add(ntuple_location + "/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");<br>
    select EMU events, where the muon is truth match and from a W<br>
    <br>
    <br>
    <b>MC pdgid catagories using elFakeMCCategory in mcSelections.cc</b><br>
    // 0: real electron from W or Z<br>
    // 1: electrons from gamma (conversion) <br>
    // 2: electron candidate or its mother is a light hadron<br>
    // 3: real electron from heavy hadron decay<br>
    // 4: the rest<br>

    <br>
    <br>
    <A HREF=#EFF>Efficiencies of various parts of selection</A><br>
    <A HREF=#EFFALL>Efficiency of full selection</A><br>
    <A HREF=#DIST>Distributions of all variables before and after application of full selection</A><br>

    <br>
    <br>

    <A NAME=EFF>
    <h1>Efficiencies of various parts of selection</h1>

	<h2>Isolation (EE): Signal (left), Background (right) <br>
        	- Isolation w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_iso$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_iso$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_iso$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_iso$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>ID (EE): Signal (left), Background (right) <br>
                - ID w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_id$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_id$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_id$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_id$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV (EE): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_conv$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_conv$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_conv$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_conv$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV W.R.T Iso+ID (EE): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. Iso+ID</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconv$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconv$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconv$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconv$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>


        <h2>Isolation (EB): Signal (left), Background (right) <br>
                - Isolation w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_iso$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_iso$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_iso$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_iso$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>ID (EB): Signal (left), Background (right) <br>
                - ID w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_id$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_id$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_id$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_id$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV (EB): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_conv$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_conv$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_conv$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_conv$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>CONV W.R.T Iso+ID (EB): Signal (left), Background (right) <br>
                - Conversion rejection w.r.t. Iso+ID</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconv$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconv$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_idisoconv$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_idisoconv$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

    <br>
    <br>
    <A NAME=EFFALL>
    <h1>Efficiency of full selection</h1>

        <h2>All Selections (EE): Signal (left), Background (right) <br>
                - All selections w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_$TAG\_pt_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_$TAG\_eta_EE.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>

        <h2>All Selections (EB): Signal (left), Background (right) <br>
                - All selections w.r.t. basic denominator</h2>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_$TAG\_pt_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
                <img src="results/$VER\_eff_s_h1_hyp_debug_after_$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH>
                <img src="results/$VER\_eff_b_h1_hyp_debug_after_$TAG\_eta_EB.png" HEIGHT=$HEIGHT WIDTH=$WIDTH><br>


    <br>
    <br>
    <A NAME=DIST>
    <h1>Distributions of all variables before and after application of full selection</h1>
	"


for FILE_S in `ls results/*.png | grep $TAG | grep $VER | grep _s_`; do

    FILE_B=`echo $FILE_S | sed 's/\(.*\)_s_\(.*\)/\1_b_\2/g'`
    FILE_SB=`echo $FILE_S | sed 's/\(.*\)_s_\(.*\)/\1_sb_\2/g'`
	TYPE=`echo $FILE_S | sed 's/.*_\(.*\)_s_.*/\1/g'`
	VAR=`echo $FILE_S | sed 's/.*_\(.*\)_E.*\.png/\1/g'`
	DET=`echo $FILE_S | sed 's/.*_\(.*\)\.png/\1/g'`

	if [ $TYPE == "overlay" ]; then
		echo "
		<h2>$VAR ($DET): Signal (left), Background (right) <br>
            - comparing S before and after all selections and likewise for B</h2>
		<img src=$FILE_S HEIGHT=$HEIGHT WIDTH=$WIDTH>
        <img src=$FILE_B HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
        "
        if [ $VAR != "pdgid" ]; then
                echo "
                <h2>- comparing S and B after all selections</h2><br>
            <img src=$FILE_SB HEIGHT=$HEIGHT WIDTH=$WIDTH><br>
            "
		fi

        if [ $VAR == "pdgid" ]; then
            echo "
                <h2>$VAR ($DET): Signal (left), Background (right) <br>
                    	- comparing S before and after all selections and likewise for B</h2>
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

