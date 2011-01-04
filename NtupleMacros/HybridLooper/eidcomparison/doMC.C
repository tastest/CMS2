
void doMC() {

    //
    // the looper
    //
    gSystem->SetIncludePath(Form("%s -I../../Tools", gSystem->GetIncludePath()));

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

    // for likelihood id
    gSystem->Load("../../Tools/EgammaAnalysisTools/lib/libElectronLikelihoodId.so");

    gROOT->ProcessLine(".L histtools.C+");

    //
    // output file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    //
    // chains for input files
    TString ntuple_location = "/nfs-3/userdata/cms2/";
    //TString ntuple_location = "/tas/cms2/";

    // zee
    TChain *chain_zee = new TChain("Events");
    //chain_zee->Add(ntuple_location +"Zee_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");
    chain_zee->Add("/tas/dlevans/data_skim/zee_skim_2020.root");
    chain_zee->Add(ntuple_location +"DYee_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_zee->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    
    // zmumu
    TChain *chain_zmumu = new TChain("Events");
    //chain_zmumu->Add(ntuple_location +"Zmumu_Spring10-START3X_V26_S09-v1/V03-04-13-07/*.root");
    chain_zmumu->Add("/tas/dlevans/data_skim/zmumu_skim_2020.root");
    chain_zmumu->Add(ntuple_location +"DYmumu_M10to20_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");
    chain_zmumu->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");

    // ztautau
    TChain *chain_ztautau = new TChain("Events");
    //chain_ztautau->Add(ntuple_location +"Ztautau_Spring10-START3X_V26_S09-v1/V03-04-13-07/*.root");
    chain_ztautau->Add("/tas/dlevans/data_skim/ztautau_skim_2020.root");
    chain_ztautau->Add(ntuple_location +"ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/dilepPt2010Skim/skimmed_ntuple.root");

    // ttbar
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add(ntuple_location + "TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/diLepPt2010Skim/skimmed_ntuple.root");

    // ttbar (38X)
    TChain *chain_ttbar38X = new TChain("Events");
    chain_ttbar38X->Add(ntuple_location + "TT_TuneZ2_7TeV-pythia6-tauola_Fall10-START38_V12-v1_GEN-SIM-RECO/V03-06-14/merged_ntuple*.root");

   // wjets 
    TChain *chain_wjets = new TChain("Events");
    chain_wjets->Add(ntuple_location + "WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-13-07/diLepPt2010Skim/skimmed_ntuple.root");

    // ZZ
    TChain *chain_zz = new TChain("Events");
    chain_zz->Add(ntuple_location + "ZZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    // WZ

    TChain *chain_wz = new TChain("Events");
    chain_wz->Add(ntuple_location + "WZ_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple.root");

    // WW
    TChain *chain_ww = new TChain("Events");
    chain_ww->Add(ntuple_location + "WW_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple.root");

    // LM0
    TChain *chain_lm0 = new TChain("Events");
    chain_lm0->Add(ntuple_location + "LM0_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged_ntuple*.root");

    // ttbar Mtt1TeV
    TChain *chain_ttbarmtt1tev = new TChain("Events");
    chain_ttbarmtt1tev->Add(ntuple_location + "TTbarJets_Mtt1TeV-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    // LM10
    TChain *chain_lm10 = new TChain("Events");
    chain_lm10->Add(ntuple_location + "LM10_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged_ntuple*.root");

    // LM6
    TChain *chain_lm6 = new TChain("Events");
    chain_lm6->Add(ntuple_location + "LM6_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged_ntuple*.root");

    // zprime
    TChain *chain_zprime = new TChain("Events");
    chain_zprime->Add(ntuple_location + "Zprime_M1TeV_W10GeV-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    // z5jet
    TChain *chain_z5jet = new TChain("Events");
    chain_z5jet->Add(ntuple_location + "Z5Jets_Pt800to1600-alpgen_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged_ntuple*.root");

    //
    // run it
    //
    // k factor for pythia samples
    //looper->ScanChain(false, "dyee", chain_zee, 1.0);
    //looper->ScanChain(false, "dymm", chain_zmumu, 1.0);
    //looper->ScanChain(false, "dytt", chain_ztautau, 1.0);
    //looper->ScanChain(false, "zz", chain_zz, 1.0);
    //looper->ScanChain(false, "wz", chain_zz, 1.0);
    //looper->ScanChain(false, "ww", chain_zz, 1.0);
    // k factor for mg sample
    looper->ScanChain(false, "wjets", chain_wjets, (31314./28049.));
    //looper->ScanChain(false, "ttbar", chain_ttbar, (157.5/165.0));
//    looper->ScanChain(false, "ttbar38X", chain_ttbar38X, (157.5/165.0));


    // BSM or signal like
    //looper->ScanChain(false, "ttbarmtt1tev", chain_ttbarmtt1tev);
    looper->ScanChain(false, "lm0", chain_lm0, (1.0));
    //looper->ScanChain(false, "lm10", chain_lm10, (1.0));
    //looper->ScanChain(false, "zprime", chain_zprime, (1.0));
    //looper->ScanChain(false, "z5jet", chain_z5jet, (1.0));
    //looper->ScanChain(false, "lm6", chain_lm6, (1.0));

    //
    // write histograms
    // 

    //const char* outFile = "histos_mc.root";
    //hist::saveHist(outFile); 
    //hist::deleteHistos();

    //
    // tidy up
    //
    delete looper;
    delete chain_zee;
    delete chain_zmumu;
    delete chain_ztautau;
    delete chain_wjets;
    delete chain_ttbar;
    delete chain_zz;
    delete chain_wz;
    delete chain_ww;
    delete chain_ttbarmtt1tev;
    delete chain_lm0;
    delete chain_lm10;
    delete chain_lm6;
    delete chain_zprime;
    delete chain_z5jet;

}

