
void doData() {

    //
    // the looper
    //
    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gSystem->Load("libCMS2NtupleMacrosCORE.so");
    gSystem->Load("libCMS2NtupleMacrosLooper.so");
    gSystem->Load("../Tools/MiniFWLite/libMiniFWLite.so");
    gROOT->ProcessLine(".L histtools.C+");

    //
    // output file for histograms
    //

    MyScanChain *looper = new MyScanChain();

    //
    // Where is the data?
    //
    //TString ntuple_location = "/tas/cms2/";
    TString ntuple_location = "/nfs-3/userdata/";

    //
    // latest data
    //
    TChain *chain_data = new TChain("Events");

    // Re reco
    chain_data->Add("/nfs-3/userdata/cms2/Mu_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/skimmed_ntuple_*.root");
    chain_data->Add("/nfs-3/userdata/cms2/EG_Run2010A-Sep17ReReco_v2_RECO/V03-06-09/diLepPt1020Skim/skimmed_ntuple_*.root");

    // Prompt reco CMS2 version V03-06-09
    chain_data->Add("/nfs-3/userdata/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-09/diLepPt1020Skim/skimmed_ntuple_*.root");
    chain_data->Add("/nfs-3/userdata/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-09/diLepPt1020Skim/skimmed_ntuple_*.root");

    // Prompt reco CMS2 version V03-06-14
    chain_data->Add("/nfs-3/userdata/cms2/Mu_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple_*.root");
    chain_data->Add("/nfs-3/userdata/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/diLepPt1020Skim/skimmed_ntuple_*.root");


    TChain *chain_minbias = new TChain("Events");
    chain_minbias->Add("/home/users/jmuelmen/CMSSW_3_6_1_patch4/src/CMS2/NtupleMacros/NtupleTools/dilep_skim_2.root");
    chain_minbias->Add("/nfs-3/userdata/fgolf/SSskims/data/skimmed_ntuple*.root");


    //  
    // run it
    //

    looper->setGoodRunList("runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt");

    looper->ScanChain(true, "data", chain_data);
    //looper->ScanChain(true, "data", chain_skim);
    //looper->ScanChain(true, "minbias", chain_minbias);

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
    delete chain_data;

}

