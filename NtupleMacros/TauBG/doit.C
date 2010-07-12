void doit(){

gROOT->LoadMacro("textFileMaker.C+");

TChain* chain1 = new TChain("Events");
chain1->Add("/tas/cms2/Ztautau_Spring10-START3X_V26_S09-v1/V03-04-08-01/*.root");
// chain1->Add("/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/*.root");

textFileMaker * b = new textFileMaker();
b->ScanChain(chain1);
}
