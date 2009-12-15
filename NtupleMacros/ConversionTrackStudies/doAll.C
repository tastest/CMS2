
{

gROOT->ProcessLine(".L MyScanChain.C++");

TFile f_data_paolo("histos_data_paolo_124120.root", "RECREATE");
f_data_paolo->cd();
TChain *chain_paolo = new TChain("Events");
chain_paolo->Add("/store/disk00/EgammaIso/V00-00-10/meridian/124120/ntuple_data.root");
ScanChain(true, chain_paolo);
f_data_paolo->Write();
f_data_paolo->Close();

TFile f_mc("histos_mc_V8D_2360GeV.root", "RECREATE");
f_mc->cd();
TChain *chain_mc = new TChain("Events");
chain_mc->Add("/store/disk00/EgammaIso/V00-00-10/mc/MinBias_Summer09-STARTUP3X_V8D_2360GeV-v2/ntuple_mc_*.root");
ScanChain(false, chain_mc);
f_mc->Write();
f_mc->Close();

}

