#!/bin/bash

#TAG="V00-02-06"
#root -b -q skimMCLeptonTree.C+\(0,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_ee.root\"\)
#root -b -q skimMCLeptonTree.C+\(1,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_mm.root\"\)

TAG="V00-02-07"
DATASET="DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1_AODSIM"
root -b -q skimMCLeptonTree.C+\(0,\"/smurf/dlevans/LeptonTree/${TAG}/${DATASET}\",\"merged_ee.root\"\)
root -b -q skimMCLeptonTree.C+\(1,\"/smurf/dlevans/LeptonTree/${TAG}/${DATASET}\",\"merged_mm.root\"\)

