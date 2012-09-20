#!/bin/bash

### 28th april good run list
CERT="Cert_190456-193336_8TeV_PromptReco_Collisions12_JSON"
CERTNAME="../runlists/${CERT}.jmu"
TAG="V00-02-00"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/PhotonsRun2012APromptV1/\",\"merged_${CERT}.root\"\)

TAG="V00-02-00_frtests2"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/\",\"merged_${CERT}.root\"\)

TAG="forBen"
CERT="Cert_190456-193557_8TeV_PromptReco_Collisions12_JSON_goodruns"
CERTNAME="/tas/benhoob/home/jsons/${CERT}.txt"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/PhotonsRun2012APromptV1_OldTrig/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/PhotonsRun2012APromptV1_NewTrig/\",\"merged_${CERT}.root\"\)

# May 18th
TAG="V00-02-02"
CERT="Cert_190456-194076_8TeV_PromptReco_Collisions12_JSON"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012APromptV1/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012APromptV1/\",\"merged_${CERT}.root\"\)

# May 25th

TAG="V00-02-03"
CERT="Cert_190456-194479_8TeV_PromptReco_Collisions12_JSON"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012/\",\"merged_${CERT}.root\"\)

#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\"\)

# 920 pb for SS
TAG="V00-02-03"
CERT="runlist_920pb_for_ss"
CERTNAME="../runlists/${CERT}.txt"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\"\)

############

# 2.96 /fb
TAG="V00-02-04"
CERT="Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\"\)


TAG="V00-02-05"
CERT="Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\",1\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\",1\)

#########

# V00-02-04 skimmed with 3.5 /fb runlist...
TAG="V00-02-04"
CERT="Cert_190456-195658_8TeV_PromptReco_Collisions12_JSON"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectronRun2012/\",\"merged_${CERT}.root\",0\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMuRun2012/\",\"merged_${CERT}.root\",0\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\",1\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\",1\)

# with bit bug fix
TAG="V00-02-05"
CERT="Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\"\)

############

# V00-02-07
# first attempt with 53X for HCP

TAG="V00-02-07"
CERT="Cert_190456-202305_8TeV_PromptReco_Collisions12_JSON"
CERTNAME="../runlists/${CERT}.jmu"
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectron_Run2012A-13Jul2012-v1_AOD_190456_193621/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectron_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectron_Run2012B-13Jul2012-v1_AOD_193834_196531/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectron_Run2012C-24Aug2012-v1_AOD_198022_198523/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleElectron_Run2012C-PromptReco-v2_AOD_198934_202950/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMu_Run2012A-13Jul2012-v1_AOD_190456_193621/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMu_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMu_Run2012B-13Jul2012-v4_AOD_193834_196531/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/DoubleMu_Run2012C-PromptReco-v2_AOD_198934_202950/\",\"merged_${CERT}.root\",0\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectron_Run2012A-13Jul2012-v1_AOD_190456_193621/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectron_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectron_Run2012B-13Jul2012-v1_AOD_193834_196531/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectron_Run2012C-24Aug2012-v1_AOD_198022_198523/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectron_Run2012C-PromptReco-v2_AOD_198934_202950/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMu_Run2012A-13Jul2012-v1_AOD_190456_193621/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMu_Run2012A-recover-06Aug2012-v1_AOD_190782_190949/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMu_Run2012B-13Jul2012-v1_AOD_193834_196531/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMu_Run2012C-24Aug2012-v1_AOD_198022_198523/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMu_Run2012C-PromptReco-v2_AOD_198934_202950/\",\"merged_${CERT}.root\",1\)

