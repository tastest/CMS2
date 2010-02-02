
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

	// BSM
	// LM0
	TChain *chain_lm0 = new TChain("Events");
	chain_lm0->Add("/store/disk02/cms2/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

	// 
	// do looping
	//

	looper->ScanChain(false, "ttbar", chain_ttbar);

	looper->ScanChain(false, "lm0", chain_lm0);

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
	delete chain_lm0;
}

