
{

	//
	// utilities
	//
	gROOT->ProcessLine(".L CMS2.cc+");
	gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");

	//
	// the looper
	//
        gROOT->ProcessLine(".L MyScanChain.C+");

	//
	// output file for histograms
	//
	TFile f_mc("histos_mc.root", "RECREATE");
	f_mc.cd();

	//
	// chains for input files
	//
	// ttbar
        TChain *chain_ttbar = new TChain("Events");
        chain_ttbar->Add("/store/disk02/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");

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

