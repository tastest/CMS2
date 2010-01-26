
{

	gROOT->ProcessLine(".L MyScanChain.C++");

	TFile f_mc("histos_mc.root", "RECREATE");
	f_mc->cd();
	TChain *chain_mc = new TChain("Events");
	chain_mc->Add("/store/disk02/cms2/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/ntuple_mc_*.root");
	ScanChain(false, chain_mc);
	f_mc->Write();
	f_mc->Close();

}

