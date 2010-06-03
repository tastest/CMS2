
void doAllData() {

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

    MyScanChain *looper = new MyScanChain(0);

    std::string ntuple_location = "/tas07/disk00/cms2";

    // DATA

    // may 6th egamma
    TChain *chain_may6eg = new TChain("Events");
    chain_may6eg->Add("/tas07/disk00/cms2/MinimumBias_Commissioning10-May6thPDSkim2_SD_EG-v1/V03-04-09/merged*.root");


    TChain *chain_whunt_skim = new TChain("Events");
    //chain_whunt_skim->Add("/tas03/disk01/whunt/skim/emuskim_*.root");
    // on dle laptop
    //chain_whunt_skim->Add("/Users/dlevans/tas03/disk01/whunt/skim/emuskim_merged.root");
    // skim merged
    //chain_whunt_skim->Add("/tmp/emuskim_merged_1400_160510.root");
    //chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim*.root");

    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132440_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132596_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132598_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132599_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132601_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132602_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132605_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132606_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132656_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132658_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132659_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132661_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132662_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132716_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132959_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132960_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132961_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132965_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_132968_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133029_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133031_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133034_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133035_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133036_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133046_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133158_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133874_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133875_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133876_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133877_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133881_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133885_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133927_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_133928_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135059_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135149_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135175_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135445_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135521_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135523_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135525_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135534_*.root");
    chain_whunt_skim->Add("/tas03/disk01/dilephunt/skim/dilepskim_135535_*.root");

    //
    // Frank Golf
    //
    TChain *ch_golf = new TChain("Events");
    ch_golf->Add(ntuple_location + "/MinimumBias_Commissioning10-May6thPDSkim2_SD_Mu-v1/V03-04-09/merged_ntuple*.root");
    ch_golf->Add(ntuple_location + "/MinimumBias_Commissioning10-SD_Mu-v9_RECO/singleLepPt5Skim/*.root"  );
    ch_golf->Add(ntuple_location + "/Mu_Run2010A-PromptReco-v1_RECO/V03-04-09-01/merged_ntuple*.root"      );
    ch_golf->Add(ntuple_location + "/Mu_Run2010A-PromptReco-v1_RECO/V03-04-10-01/merged_ntuple*.root"      );
    ch_golf->Add(ntuple_location + "/Mu_Run2010A-PromptReco-v2_RECO/V03-04-10-01/merged_ntuple*.root"      );
    ch_golf->Add(ntuple_location + "/MinimumBias_Commissioning10-May6thPDSkim2_SD_EG-v1/V03-04-09/merged_ntuple*.root");
    ch_golf->Add(ntuple_location + "/MinimumBias_Commissioning10-SD_EG-v9_RECO/singleLepPt5Skim/*.root"  );
    ch_golf->Add(ntuple_location + "/EG_Run2010A-PromptReco-v1_RECO/V03-04-09-01/merged_ntuple*.root"      );
    ch_golf->Add(ntuple_location + "/EG_Run2010A-PromptReco-v1_RECO/V03-04-10-01/merged_ntuple*.root"      );
    ch_golf->Add(ntuple_location + "/EG_Run2010A-PromptReco-v2_RECO/V03-04-10-01/merged_ntuple*.root"      );


    // do the looping
    //looper->ScanChain(true, "whunt", chain_whunt_skim);
    looper->ScanChain(true, "may6eg", chain_may6eg);

    //
    // write histograms
    // 
    const char* outFile = "histos_data_spike.root";
    hist::saveHist(outFile); 
    hist::deleteHistos();

    //
    // tidy up
    //
    delete looper;
    delete chain_whunt_skim;
    delete chain_may6eg;
}

