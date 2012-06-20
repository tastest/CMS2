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
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\",1\)
root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\",1\)

# with bit bug fix
TAG="V00-02-05"
CERT="Cert_190456-195396_8TeV_PromptReco_Collisions12_JSON_v2"
CERTNAME="../runlists/${CERT}.jmu"
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleElectronRun2012/\",\"merged_${CERT}.root\"\)
#root -b -q skimLeptonTree.C+\(\"${CERTNAME}\",\"/smurf/dlevans/LeptonTree/${TAG}/SingleMuRun2012/\",\"merged_${CERT}.root\"\)



