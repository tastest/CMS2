
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

    gROOT->ProcessLine(".L ../histtools.C+");

	//
	// output file for histograms
	//

//
// danger!  keep this up to date
//
enum ElectronSelection {
    ELEPASS_PT10,
    ELEPASS_PT20,
    ELEPASS_PT10NOT20,
    ELEPASS_DPHI,
    ELEPASS_DETA,
    ELEPASS_HOE,
    ELEPASS_LSHAPE,
    ELEPASS_ISO,
    ELEPASS_EXTRA,
    ELEPASS_D0,
    ELEPASS_ID,
    ELEPASS_NOTCONV,
    ELEPASS_NOMUON,
    ELEPASS_TYPE,
    ELEPASS_FIDUCIAL,
    ELEPASS_FULLSELECTION
};


    TString fileNameString = "pt20up";
    //TString fileNameString = "pt10to20";
    //TString fileNameString = "pt10up";
    elecuts_t configured_cuts = (1<<ELEPASS_PT20);
    //elecuts_t configured_cuts = (1<<ELEPASS_PT10NOT20);
    //elecuts_t configured_cuts = (1<<ELEPASS_PT10);
	MyScanChain *looper = new MyScanChain(configured_cuts);

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
    //chain_wjets->Add(ntuple_location + "/cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	chain_wjets->Add(ntuple_location + "/cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/dilep-skim/wjets_skim.root");
	// qcd pt30
	TChain *chain_qcd30 = new TChain("Events");
	chain_qcd30->Add(ntuple_location + "/cms2/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
	// wmunu
    TChain *chain_wmunu = new TChain("Events");
    chain_wmunu->Add(ntuple_location + "/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
    // inclusive mu15 dilep filt
    TChain *chain_inclmu15 = new TChain("Events");
    chain_inclmu15->Add(ntuple_location + "/cms2/InclusiveMu15_Summer09-MC_31X_V3_7TeV-v1_dilepfilt/V03-00-35/merged_ntuple*.root");

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

    // test ntuple with new ele id computed in cmssw
    // from
    // '/store/relval/CMSSW_3_5_2/RelValZEE/GEN-SIM-RECO/MC_3XY_V21-v1/0016/*.root',
    //
    TChain *chain_eleidval = new TChain("Events");
    chain_eleidval->Add(ntuple_location + "/dlevans/relval_zee_ntuple.root");

	// 
	// do looping
	//
//	looper->ScanChain(false, "ww", chain_ww);
//	looper->ScanChain(false, "wz", chain_wz);
//	looper->ScanChain(false, "zz", chain_zz);
//	looper->ScanChain(false, "dyee", chain_dyee);
//	looper->ScanChain(false, "dymm", chain_dymm);
//    looper->ScanChain(false, "dytt", chain_dytt);
//    looper->ScanChain(false, "wjets", chain_wjets);
//	looper->ScanChain(false, "elegunstartup", chain_elegunstartup);
//  looper->ScanChain(false, "elegunideal", chain_elegunideal);
	//looper->ScanChain(false, "QCDpt30", chain_qcd30);
    //looper->ScanChain(false, "InclusiveMuPt15", chain_inclmu15);


    // ele id studies
    //

    //looper->ScanChain(false, "ttbar", chain_ttbar);
	//looper->ScanChain(false, "wm", chain_wmunu);

    // ele id sanity check with sani id
    looper->ScanChain(false, "eleidval", chain_eleidval);


	//
	// write histograms
	// 

	const char* outFile = "histos_eleid_saniv02_" + fileNameString + ".root";
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
    delete chain_inclmu15;
	delete chain_qcd30;
	delete chain_lm0;
	delete chain_elegunstartup;
	delete chain_elegunideal;
    delete chain_eleidval;

}

