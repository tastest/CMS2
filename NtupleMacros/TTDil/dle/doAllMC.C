
void doAllMC(unsigned int electronId) {

    //
    // the looper
    //
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libGenVector.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");
    gROOT->ProcessLine(".L ../histtools.C+");

    MyScanChain *looper = new MyScanChain(electronId);

    //
    // chains for input files
    TString ntuple_location = "/tas07/disk00/cms2/";

    // SM
    // ttbar
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add(ntuple_location + "TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/merged_ntuple*.root");
    //chain_ttbar->Add("/tas01/disk02/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root");
    // ww
/*
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
*/
    // wjets
    TChain *chain_wjets = new TChain("Events");
    chain_wjets->Add(ntuple_location + "WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/dilep-skim.root");
    
    // zjets
    TChain *chain_zjets = new TChain("Events");
    chain_zjets->Add(ntuple_location + "ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root");

    // qcd pt30
    TChain *chain_qcd30 = new TChain("Events");
    chain_qcd30->Add(ntuple_location + "QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root");
/*
    // wmunu
    TChain *chain_wmunu = new TChain("Events");
    chain_wmunu->Add(ntuple_location + "/cms2/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root");
*/

    // 
    // do looping
    //

    //looper->ScanChain(false, "ttbar", chain_ttbar);
    //looper->ScanChain(false, "ww", chain_ww);
    //looper->ScanChain(false, "wz", chain_wz);
    //looper->ScanChain(false, "zz", chain_zz);
    //looper->ScanChain(false, "dyee", chain_dyee);
    //looper->ScanChain(false, "dymm", chain_dymm);
    //looper->ScanChain(false, "dytt", chain_dytt);
    //looper->ScanChain(false, "wjets", chain_wjets);
    //looper->ScanChain(false, "zjets", chain_zjets);
    looper->ScanChain(false, "QCDpt30", chain_qcd30);
    //looper->ScanChain(false, "wm", chain_wmunu);

    //
    // write histograms
    // 

    const char* outFile = "histos_mc_spike.root";
    hist::saveHist(outFile); 
    hist::deleteHistos();

    //
    // tidy up
    //
    delete looper;

    delete chain_ttbar;
/*    delete chain_ww;
    delete chain_wz;
    delete chain_zz;
    delete chain_dyee;
    delete chain_dymm;
    delete chain_dytt;
*/    delete chain_wjets;
    delete chain_zjets;
/*    delete chain_wmunu;
    delete chain_qcd30;
*/

}

