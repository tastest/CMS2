#!/bin/bash

REV=1.1

# make include dir
if [ ! -d EgammaAnalysisTools/data ]; then
    mkdir -p EgammaAnalysisTools/data
fi

# NoIPInfo uses no IP information
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0LowPt_NoIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1LowPt_NoIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2LowPt_NoIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0HighPt_NoIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1HighPt_NoIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2HighPt_NoIPInfo_BDTG.weights.xml

# WithIPInfo uses d0, IP3D, IP3D sig 
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0LowPt_WithIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1LowPt_WithIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2LowPt_WithIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0HighPt_WithIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1HighPt_WithIPInfo_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2HighPt_WithIPInfo_BDTG.weights.xml

# IDIsoCombined uses ID+IP+Iso
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml
cvs co -r $REV -p UserCode/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml  > EgammaAnalysisTools/data/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml

