
void doAll() {

	// from slavas code
	gSystem->Load("libGui.so");
	gSystem->Load("libPhysics.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

	// Load and compile something to allow proper treatment of vectors
	// Not clear that it is needed
	gSystem->CompileMacro("loader.C", "++k", "libloader");
	// end from slava

	//
	// utilities
	//
	gROOT->ProcessLine(".L CMS2.cc+");

    gROOT->ProcessLine(".L ../../CORE/mcSelections.cc+");
    gROOT->ProcessLine(".L ../../CORE/trackSelections.cc+");
    gROOT->ProcessLine(".L ../../CORE/jetSelections.cc+");
    gROOT->ProcessLine(".L ../../CORE/metSelections.cc+");
    gROOT->ProcessLine(".L ../../CORE/eventSelections.cc+");
	gROOT->ProcessLine(".L ../../CORE/electronSelections.cc+");
	gROOT->ProcessLine(".L ../../CORE/utilities.cc+");
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
	//TString ntuple_location = "/data/tmp/";
    TString ntuple_location = "/store/disk02/";


	// SM
	// ttbar
	TChain *chain_ttbar = new TChain("Events");
	chain_ttbar->Add(ntuple_location + "/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");
	// ww
	TChain *chain_ww = new TChain("Events");
	chain_ww->Add(ntuple_location + "/cms2/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// wz
	TChain *chain_wz = new TChain("Events");
	chain_wz->Add(ntuple_location + "/cms2/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// zz
	TChain *chain_zz = new TChain("Events");
	chain_zz->Add(ntuple_location + "/cms2/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// dyee
	TChain *chain_dyee = new TChain("Events");
	chain_dyee->Add(ntuple_location + "/cms2/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merged_ntuple*.root");
	// dymm
	TChain *chain_dymm = new TChain("Events");
	chain_dymm->Add(ntuple_location + "/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// dytt
    TChain *chain_dytt = new TChain("Events");
    chain_dytt->Add(ntuple_location + "/cms2/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// wjets
    TChain *chain_wjets = new TChain("Events");
    chain_wjets->Add(ntuple_location + "/cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	//chain_wjets->Add(ntuple_location + "/cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/dilep-skim/wjets_skim.root");
	// qcd pt30
	TChain *chain_qcd30 = new TChain("Events");
	chain_qcd30->Add(ntuple_location + "/cms2/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// wmunu
    TChain *chain_wmunu = new TChain("Events");
    chain_wmunu->Add(ntuple_location + "/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// photonjets
	TChain *chain_photonjets = new TChain("Events");
    chain_photonjets->Add(ntuple_location + "/cms2/PhotonJet_Pt170to300_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    chain_photonjets->Add(ntuple_location + "/cms2/PhotonJet_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    chain_photonjets->Add(ntuple_location + "/cms2/PhotonJet_Pt30to50_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    chain_photonjets->Add(ntuple_location + "/cms2/PhotonJet_Pt50to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    chain_photonjets->Add(ntuple_location + "/cms2/PhotonJet_Pt80to120_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");


	// BSM
	// LM0
	TChain *chain_lm0 = new TChain("Events");
	chain_lm0->Add(ntuple_location + "/cms2/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

    TChain *chain_lm4 = new TChain("Events");
    chain_lm4->Add(ntuple_location + "/cms2/LM4_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");

	// Technical
	// single particle gun electrons
	TChain *chain_elegunstartup = new TChain("Events");
	chain_elegunstartup->Add(ntuple_location + "/cms2/SingleElectronPt5to100_336patch4/V03-00-35/merged_ntuple*.root");
    TChain *chain_elegunideal = new TChain("Events");
    chain_elegunideal->Add(ntuple_location + "/cms2/SingleElectronPt5to100_336patch4MC31XV9/V03-00-35/merged_ntuple*.root");

    // data
    //
    TChain *chain_v0 = new TChain("Events");
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/*.root.v0");



	// 
	// do looping
	//

	//looper->ScanChain(false, "ttbar", chain_ttbar);
//	looper->ScanChain(false, "ww", chain_ww);
//	looper->ScanChain(false, "wz", chain_wz);
//	looper->ScanChain(false, "zz", chain_zz);
//	looper->ScanChain(false, "dyee", chain_dyee);
//	looper->ScanChain(false, "dymm", chain_dymm);
//    looper->ScanChain(false, "dytt", chain_dytt);
//    looper->ScanChain(false, "wjets", chain_wjets);


//	looper->ScanChain(false, "elegunstartup", chain_elegunstartup);
//  looper->ScanChain(false, "elegunideal", chain_elegunideal);

//	looper->ScanChain(false, "QCDpt30", chain_qcd30);
//	looper->ScanChain(false, "wm", chain_wmunu);
//	looper->ScanChain(false, "photonjets", chain_photonjets);

//    looper->ScanChain(false, "lm0", chain_lm0);
//    looper->ScanChain(false, "lm4", chain_lm4);

    looper->ScanChain(true, "v0", chain_v0);

	//
	// write histograms
	// 

	const char* outFile = "histos_mc.root";
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
	delete chain_dytt;
	delete chain_wjets;
	delete chain_wmunu;

	delete chain_qcd30;
	delete chain_photonjets;

	delete chain_lm0;
    delete chain_lm4;

	delete chain_elegunstartup;
	delete chain_elegunideal;

    delete chain_v0;

}

