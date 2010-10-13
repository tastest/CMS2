
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
    gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

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
    };


    TString fileNameString = "pt20up";
    //TString fileNameString = "pt10to20";
    //TString fileNameString = "pt10up";
    cuts_t configured_cuts = (1<<ELEPASS_PT20);
    //cuts_t configured_cuts = (1<<ELEPASS_PT10NOT20);
    //cuts_t configured_cuts = (1<<ELEPASS_PT10);
    MyScanChain *looper = new MyScanChain(configured_cuts);

    //
    // chains for input files
    TString ntuple_location = "/store/disk00/cms2/";

    // SM
    // ttbar
    TChain *chain_ttbar = new TChain("Events");
    chain_ttbar->Add(ntuple_location + "TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/merged_ntuple*.root");
    // wjets
    TChain *chain_wjets = new TChain("Events");
    chain_wjets->Add(ntuple_location + "doesnotexist/wjets_skim.root");
    // qcd pt30
    TChain *chain_qcd30 = new TChain("Events");
    chain_qcd30->Add(ntuple_location + "QCD_Pt30_Spring10-START3X_V26_S09-v1/V03-04-08/merged_ntuple*.root");

    //
    // Validation
    TChain *chain_eleidval = new TChain("Events");
    chain_eleidval->Add("/tmp/dlevans/ntuple.root");
    

    // 
    // do looping
    //

    //looper->ScanChain(false, "ttbar", chain_ttbar);
    //looper->ScanChain(false, "wjets", chain_wjets);
    //looper->ScanChain(false, "QCDpt30", chain_qcd30);

    looper->ScanChain(false, "eleidval", chain_eleidval);

    //
    // write histograms
    // 
    const char* outFile = "histos_eleid_" + fileNameString + ".root";
    hist::saveHist(outFile); 
    hist::deleteHistos();

    //
    // tidy up
    //
    //delete looper;
    delete chain_ttbar;
    delete chain_wjets;
    delete chain_qcd30;
    delete chain_eleidval;

}

