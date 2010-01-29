

void doAll2X() {

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

	//
	// the looper
	//
	gROOT->ProcessLine(".L MyScanChain.C+");

	//
	// output file for histograms
	//
	TFile f_mc("histos_mc_2x.root", "RECREATE");
	f_mc.cd();

	//
	// chains for input files
	//
	// ttbar
	TChain *chain_ttbar = new TChain("Events");
	//chain_ttbar->Add("/store/disk02/cms2-V02-00-12/TTBar_Summer09-MC_31X_V3_7TeV-v1/merged_ntuple*.root");
        chain_ttbar->Add("/store/disk01/cms2-V01-03-01/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged_ntuple*.root");

	// 
	// do looping
	//
	// ttbar
	ScanChain(false, "ttbar", chain_ttbar);

	//
	// write histograms
	// 

	f_mc.Write();
	f_mc.Close();

}

