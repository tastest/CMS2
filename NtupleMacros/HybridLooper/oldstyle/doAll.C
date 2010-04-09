
void doAll() {

	//
	// the looper
	//
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");

    gROOT->ProcessLine(".L histtools.C+");

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

	// Technical
	// single particle gun electrons
	TChain *chain_elegunstartup = new TChain("Events");
	chain_elegunstartup->Add(ntuple_location + "/cms2/SingleElectronPt5to100_336patch4/V03-00-35/merged_ntuple*.root");
    TChain *chain_elegunideal = new TChain("Events");
    chain_elegunideal->Add(ntuple_location + "/cms2/SingleElectronPt5to100_336patch4MC31XV9/V03-00-35/merged_ntuple*.root");

    // 7 TeV minbias
    TChain *chain_minbias = new TChain("Events");
    chain_minbias->Add("/tas01/disk01/cms2/MinBias_Spring10-START3X_V25B_356ReReco-v1/V03-03-07/merged_ntuple*.root");

    // data
    //
    // the slava-tuple
    TChain *chain_v0 = new TChain("Events");
    // run 440?
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/_cms__store_express_Commissioning10_ExpressPhysics_FEVT_v7_*_*_440_*.root.v0");
    // run 442?
    chain_v0->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/_cms__store_express_Commissioning10_ExpressPhysics_FEVT_v7_*_*_442_*.root.v0");

    // run 132569
    TChain *chain_132569 = new TChain("Events");
    chain_132569->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/_cms__store_express_Commissioning10_ExpressPhysics_FEVT_v*_*_*_569_*.root.v0");

    // run 132579
    TChain *chain_132572 = new TChain("Events");
    chain_132572->Add("/tas03/disk02/slava77/reltestdata/CMSSW_3_5_6-cms2-data/_cms__store_express_Commissioning10_ExpressPhysics_FEVT_v*_*_*_572_*.
root.v0");


    // the goodcoll-tuple
    TChain *chain_goodcoll = new TChain("Events");
    chain_goodcoll->Add("/tas01/disk01/cms2/MinimumBias_Commissioning10-GOODCOLL-v7_r132440_r132442/V03-03-07/merged*.root");

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
	//looper->ScanChain(false, "wm", chain_wmunu);
//	looper->ScanChain(false, "photonjets", chain_photonjets);

    //looper->ScanChain(true, "v0", chain_v0);

    //looper->ScanChain(true, "run132569", chain_132569);
    //looper->ScanChain(true, "run132572", chain_132572);

    looper->ScanChain(true, "goodcoll", chain_goodcoll);

    //looper->ScanChain(false, "minbias", chain_minbias);

	//
	// write histograms
	// 

	const char* outFile = "histos_data.root";
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

	delete chain_elegunstartup;
	delete chain_elegunideal;

    delete chain_v0;
    delete chain_goodcoll;

    delete chain_minbias;

    delete chain_132569;
    delete chain_132572;

}

