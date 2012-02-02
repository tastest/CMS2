#!/bin/bash

REV=1.2

# make include dir
if [ ! -d MuonAnalysisTools/data ]; then
    mkdir -p MuonAnalysisTools/data
fi

# IDIsoCombined uses ID+IP+Iso
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml > MuonAnalysisTools/data/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml
