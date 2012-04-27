#!/bin/bash

###
CERT="../runlists/HWW_2011.jmu"
TAG="V00-01-05"
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2011BPromptV1/\",\"merged_HWW_2011.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2011BPromptV1/\",\"merged_HWW_2011.root\"\)

###
CERT="../runlists/Cert_190456-191276_8TeV_PromptReco_Collisions12_JSON.jmu"
TAG="V00-01-06"
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)
root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/PhotonsRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)

##
CERT="../runlists/DCSOnly_2012A_230412.jmu"
TAG="V00-01-05"
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_DCSOnly_2012A_230412.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_DCSOnly_2012A_230412.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/\",\"merged_DCSOnly_2012A_230412.root\"\)

