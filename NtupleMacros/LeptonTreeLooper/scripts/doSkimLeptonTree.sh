#!/bin/bash
root -b -q skimLeptonTree.C+\(\"../runlists/DCSOnly_2012A_160412.jmu\",\"/smurf/dlevans/LeptonTree/V00-01-02/DoubleElectronRun2012APromptV1/merged.root\",\"/smurf/dlevans/LeptonTree/V00-01-02/DoubleElectronRun2012APromptV1/merged_dcs.root\"\)
root -b -q skimLeptonTree.C+\(\"../runlists/DCSOnly_2012A_160412.jmu\",\"/smurf/dlevans/LeptonTree/V00-01-02/DoubleElectronRun2012APromptV1_egammaMVAFix/merged.root\",\"/smurf/dlevans/LeptonTree/V00-01-02/DoubleElectronRun2012APromptV1_egammaMVAFix/merged_dcs.root\"\)

root -b -q skimLeptonTree.C+\(\"../runlists/DCSOnly_2012A_160412.jmu\",\"/smurf/dlevans/LeptonTree/V00-01-02/SingleElectronRun2012APromptV1/merged.root\",\"/smurf/dlevans/LeptonTree/V00-01-02/SingleElectronRun2012APromptV1/merged_dcs.root\"\)

root -b -q skimLeptonTree.C+\(\"../runlists/DCSOnly_2012A_160412.jmu\",\"/smurf/dlevans/LeptonTree/V00-01-02/SingleMuRun2012APromptV1/merged.root\",\"/smurf/dlevans/LeptonTree/V00-01-02/SingleMuRun2012APromptV1/merged_dcs.root\"\)

