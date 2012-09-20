#!/bin/bash

TAG="V00-02-06"
root -b -q skimMCLeptonTree.C+\(0,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_ee.root\"\)
root -b -q skimMCLeptonTree.C+\(1,\"/smurf/dlevans/LeptonTree/${TAG}/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/\",\"merged_mm.root\"\)


