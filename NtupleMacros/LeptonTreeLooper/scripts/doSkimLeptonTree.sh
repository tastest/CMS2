#!/bin/bash

###
CERT="../runlists/Cert_190456-191276_8TeV_PromptReco_Collisions12_JSON.jmu"
TAG="V00-01-05"
root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)
root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_Cert_190456-191276.root\"\)

##
CERT="../runlists/DCSOnly_2012A_230412.jmu"
TAG="V00-01-05"
root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_DCSOnly_2012A_230412.root\"\)
root -b -q skimLeptonTree.C+\(\"${CERT}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_DCSOnly_2012A_230412.root\"\)

