#!/bin/bash

REV=edm-20Oct2010
PDFS=/afs/cern.ch/user/e/emanuele/public/4Likelihood/PDFsSQLite/CMSSW_3_5_X/pdfs_MC.root

# likelihood id
cvs co -d EgammaAnalysisTools/src/ -r $REV UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodPdf.cc
cvs co -d EgammaAnalysisTools/include/ -r $REV UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodPdf.h
cvs co -d EgammaAnalysisTools/src/ -r $REV UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodSpecies.cc
cvs co -d EgammaAnalysisTools/include/ -r $REV UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodSpecies.h
cvs co -d EgammaAnalysisTools/src/ -r $REV UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodPdfProduct.cc
cvs co -d EgammaAnalysisTools/include/ -r $REV UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodPdfProduct.h
cvs co -d EgammaAnalysisTools/src/ -r $REV UserCode/emanuele/EgammaAnalysisTools/src/ElectronLikelihood.cc
cvs co -d EgammaAnalysisTools/include/ -r $REV UserCode/emanuele/EgammaAnalysisTools/include/ElectronLikelihood.h

# other includes needed
cvs co -d EgammaAnalysisTools/include -r $REV UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodMeasurements.h
cvs co -d EgammaAnalysisTools/include -r $REV UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodSwitches.h

# my example and makefile
cvs co -p UserCode/DLEvans/LikelihoodElectronId/Makefile > EgammaAnalysisTools/Makefile
cvs co -p UserCode/DLEvans/LikelihoodElectronId/Makefile.arch > EgammaAnalysisTools/Makefile.arch
cvs co -p UserCode/DLEvans/LikelihoodElectronId/include/LikelihoodUtil.h > EgammaAnalysisTools/include/LikelihoodUtil.h
cvs co -p UserCode/DLEvans/LikelihoodElectronId/src/LikelihoodUtil.cc > EgammaAnalysisTools/src/LikelihoodUtil.cc

# make directory to hold libs
if [ ! -d EgammaAnalysisTools/lib ]; then
    mkdir EgammaAnalysisTools/lib
fi

# make directory to hold pdfs
if [ ! -d EgammaAnalysisTools/PDFs ]; then
    mkdir EgammaAnalysisTools/PDFs
fi

# get pdfs
scp lxplus.cern.ch:$PDFS EgammaAnalysisTools/PDFs

