#!/bin/bash

REV=1.2

# make include dir
if [ ! -d EgammaAnalysisTools/data ]; then
    mkdir -p EgammaAnalysisTools/data
fi

# likelihood id
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet0HighPtWithLHV3_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet0LowPtWithLHV3_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet1HighPtWithLHV3_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet1LowPtWithLHV3_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet2HighPtWithLHV3_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPtWithLHV3_BDTG.weights.xml > EgammaAnalysisTools/data/Subdet2LowPtWithLHV3_BDTG.weights.xml

