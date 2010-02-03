
void doAll3X() {

	// from slavas code
	gSystem->Load("libGui.so");
	gSystem->Load("libPhysics.so");

	// Load and compile something to allow proper treatment of vectors
	// Not clear that it is needed
	gSystem->CompileMacro("loader.C", "++k", "libloader");
	// end from slava

	//
	// utilities
	//
	gROOT->ProcessLine(".L CMS2.cc+");

	gROOT->ProcessLine(".L ../../CORE/electronSelections.cc+");
	gROOT->ProcessLine(".L ../../CORE/utilities.cc+");
	gROOT->ProcessLine(".L ../../CORE/selections.cc+");
	gROOT->ProcessLine(".L ../histtools.C+");

	//
	// the looper
	//
	gROOT->ProcessLine(".L MyScanChain.cc+");

	//
	// output file for histograms
	//

	MyScanChain *looper = new MyScanChain();

	//
	// chains for input files
	// 
	// SM
	// ttbar
	TChain *chain_ttbar = new TChain("Events");
	chain_ttbar->Add("/store/disk02/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");
	// ww
	TChain *chain_ww = new TChain("Events");
	chain_ww->Add("/store/disk02/cms2/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    // wz
    TChain *chain_wz = new TChain("Events");
    chain_wz->Add("/store/disk02/cms2/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    // zz
    TChain *chain_zz = new TChain("Events");
    chain_zz->Add("/store/disk02/cms2/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// dyee
    TChain *chain_dyee = new TChain("Events");
    chain_dyee->Add("/store/disk02/cms2/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merged_ntuple*.root");
	// dymm
    TChain *chain_dymm = new TChain("Events");
    chain_dymm->Add("/store/disk02/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");


	// BSM
	// LM0
	TChain *chain_lm0 = new TChain("Events");
	chain_lm0->Add("/store/disk02/cms2/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

	// 
	// do looping
	//

	looper->ScanChain(false, "ttbar", chain_ttbar);
	looper->ScanChain(false, "ww", chain_ww);
	looper->ScanChain(false, "wz", chain_wz);
	looper->ScanChain(false, "zz", chain_zz);
	looper->ScanChain(false, "dyee", chain_dyee);
	looper->ScanChain(false, "dymm", chain_dymm);

	//
	// write histograms
	// 

	const char* outFile = "histos_mc_3x.root";
	hist::saveHist(outFile); 
	hist::deleteHistos();

	//
	// tidy up
	//
	delete looper;

	delete chain_ttbar;
	delete chain_ww;
	delete chain_wz;
	delete chain_zz;
	delete chain_dyee;
	delete chain_dymm;
	delete chain_lm0;
}

