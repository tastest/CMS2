{

using namespace std;
gROOT->ProcessLine(".x setup.C(true)");
gROOT->LoadMacro("ScanChain.C++");
//gROOT->ProcessLine(".x runLooper.C");
TChain *chain = new TChain("Events");
//chain->Add("/store/disk02/gutsche/cms2-V1-02-06/SingleElectronFlatPt5To100/merged_ntuple*.root");
//chain->Add("/store/disk01/yanjuntu/SingleElectron83055667ceda31da68d5aa532f72919d/preprocessing/ntuple*.root");
chain->Add("/store/disk01/yanjuntu/SingleGamma4016bf049eded18e8e9dac4c09773936/preprocessing/ntuple*.root");
ScanChain(chain);

const char* outFile = "myHist.root";
hist::saveHist(outFile);
hist::deleteHistos();
}
