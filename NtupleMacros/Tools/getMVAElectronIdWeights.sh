#!/bin/bash

REV=1.1

# make include dir
if [ ! -d EgammaAnalysisTools/data ]; then
    mkdir -p EgammaAnalysisTools/data
fi

# V1 uses no IP information
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0LowPt_V1_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1LowPt_V1_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2LowPt_V1_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0HighPt_V1_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1HighPt_V1_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_V1_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2HighPt_V1_BDTG.weights.xml

# V2 uses d0, IP3D, IP3D sig 
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0LowPt_V2_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1LowPt_V2_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2LowPt_V2_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0HighPt_V2_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1HighPt_V2_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_V2_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2HighPt_V2_BDTG.weights.xml

