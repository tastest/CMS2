#!/bin/bash

#REV=edm-Feb2011
#PDFS=/afs/cern.ch/user/e/emanuele/public/4Likelihood/PDFsSQLite/CMSSW_3_5_X/pdfs_MC.root
REV=edm06May11
PDFS=/afs/cern.ch/user/e/emanuele/public/4Likelihood/PDFsSQLite/CMSSW_4_1_X/pdfs_MC.root

# make include dir
if [ ! -d EgammaAnalysisTools/include ]; then
    mkdir -p EgammaAnalysisTools/include
fi
# make source dir
if [ ! -d EgammaAnalysisTools/src ]; then
    mkdir -p EgammaAnalysisTools/src
fi

# likelihood id
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodPdf.cc             > EgammaAnalysisTools/src/LikelihoodPdf.cc
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodPdf.h          > EgammaAnalysisTools/include/LikelihoodPdf.h
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodSpecies.cc         > EgammaAnalysisTools/src/LikelihoodSpecies.cc
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodSpecies.h      > EgammaAnalysisTools/include/LikelihoodSpecies.h
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/src/LikelihoodPdfProduct.cc      > EgammaAnalysisTools/src/LikelihoodPdfProduct.cc
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodPdfProduct.h   > EgammaAnalysisTools/include/LikelihoodPdfProduct.h
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/src/ElectronLikelihood.cc        > EgammaAnalysisTools/src/ElectronLikelihood.cc
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/ElectronLikelihood.h     > EgammaAnalysisTools/include/ElectronLikelihood.h

# other includes needed
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodMeasurements.h > EgammaAnalysisTools/include/LikelihoodMeasurements.h
cvs co -r $REV -p UserCode/emanuele/EgammaAnalysisTools/include/LikelihoodSwitches.h     > EgammaAnalysisTools/include/LikelihoodSwitches.h

# my example and makefile
cvs co -p UserCode/DLEvans/LikelihoodElectronId/Makefile >                  EgammaAnalysisTools/Makefile
cvs co -p UserCode/DLEvans/LikelihoodElectronId/Makefile.arch >             EgammaAnalysisTools/Makefile.arch
cvs co -p UserCode/DLEvans/LikelihoodElectronId/include/LikelihoodUtil.h >  EgammaAnalysisTools/include/LikelihoodUtil.h
cvs co -p UserCode/DLEvans/LikelihoodElectronId/src/LikelihoodUtil.cc >     EgammaAnalysisTools/src/LikelihoodUtil.cc

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

