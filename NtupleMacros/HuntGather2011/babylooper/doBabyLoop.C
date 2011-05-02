{

    gROOT->ProcessLine(".L BabyLooper.cc+");
    gROOT->ProcessLine(".L histtools.C+");

    //
    // the baby looper
    //

    BabyLooper *babyLooper = new BabyLooper();

    //
    // the baby mc
    //

    TChain *ch_baby_dyeemm = new TChain("tree");
    ch_baby_dyeemm->Add("/tas/cms2/gather/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root");

    //
    // the baby data
    //


    // 2011 double electron re-reco
    TChain *ch_baby_El2011AApr22ReReco = new TChain("tree");
    ch_baby_El2011AApr22ReReco->Add("/tas/cms2/gather/data/CMSSW_4_1_2_patch1_V04-01-05/DoubleElectron_Run2011A-Apr22ReReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-05_merged/V04-01-05/baby_gather.root");

    // 2011 double electron prompt
    TChain *ch_baby_El2011APrompt = new TChain("tree");
    ch_baby_El2011APrompt->Add("/tas/cms2/gather/data/CMSSW_4_1_2_patch1_V04-00-13/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather.root");
    ch_baby_El2011APrompt->Add("/tas/cms2/gather/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleElectron_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby*.root");

    // 2011 double muon prompt
    TChain *ch_baby_Mu2011APrompt = new TChain("tree");
    ch_baby_Mu2011APrompt->Add("/tas/cms2/gather/data/CMSSW_4_1_2_patch1_V04-00-13/DoubleMu_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather.root");
    ch_baby_Mu2011APrompt->Add("/tas/cms2/gather/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleMu_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby*.root");

    // 2011 EXPRESS
    TChain *ch_baby_Express2011 = new TChain("tree");
    ch_baby_Express2011->Add("/tas/cms2/gather/data/ExpressPhysics_Run2011A-Express-v2_FEVT/V04-01-02/baby*.root");
    ch_baby_Express2011->Add("/tas/cms2/gather/data/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/baby_gather.root");

    // 2010B Muons
    TChain *ch_baby_mu2010B = new TChain("tree");
    ch_baby_mu2010B->Add("/tas/cms2/gather/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root");

    // 2010B Electrons
    TChain *ch_baby_el2010B = new TChain("tree");
    ch_baby_el2010B->Add("/tas/cms2/gather/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");

    //
    // do the looping
    //

    // DATA

    babyLooper->setGoodRunList("../runlists/dcs_jmu.txt");
    babyLooper->Loop("El2011AApr22ReReco", ch_baby_El2011AApr22ReReco);
    babyLooper->Loop("El2011APromptR1", ch_baby_El2011APrompt, 160329, 163333);
    babyLooper->Loop("El2011APromptR2", ch_baby_El2011APrompt, 163334, 999999);
    babyLooper->Loop("Mu2011APrompt", ch_baby_Mu2011APrompt); 
    babyLooper->Loop("Express2011", ch_baby_Express2011);

    babyLooper->setGoodRunList("../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt");
    babyLooper->Loop("Mu2010BNov4", ch_baby_mu2010B);
    babyLooper->Loop("El2010BNov4", ch_baby_el2010B);


    // MC
    babyLooper->Loop("dyeemm", ch_baby_dyeemm);     

    //
    // write histograms
    // 

    const char* outFile = "histos_baby.root";
    hist::saveHist(outFile);

    delete babyLooper;

}




