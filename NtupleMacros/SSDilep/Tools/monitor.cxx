{
gROOT->LoadMacro("monitor.C+");

TChain* chain1 = new TChain("Events","TTbar-madgraph");
chain1->Add("/data/tmp/dmytro/cms2-V01-02-01/TTJets-madgraph/*.root");

TChain* chain2 = new TChain("Events","WW-2l");
chain2->Add("/data/tmp/dmytro/cms2-V01-02-01/WW_2l-Pythia/merged_ntuple.root");

TChain* chain3 = new TChain("Events","DY");
chain3->Add("/data/tmp/dmytro/cms2-V01-02-01/ZJets-madgraph/*.root");

TChain* chain4 = new TChain("Events","WJets");
chain4->Add("/data/tmp/dmytro/cms2-V01-02-01/WJets-madgraph/*.root");

drawAll(chain1,chain2,chain3,chain4);
}
