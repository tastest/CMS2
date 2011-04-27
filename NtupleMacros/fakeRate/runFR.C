#include "ChainFromText.cc"
void runFR(){

gROOT->LoadMacro("myBabyMaker.C++");

//////////
// 2011 //
//////////

  // Data

    // Double Electron
    TChain *chain1 = new TChain("Events");
    //chain1->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-01-02/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-01-02_merged/V04-01-02/*.root");
    chain1->Add("/home/users/dbarge/ntuple_production/devel_CMSSW_4_1_2_patch1_V04-01-03/crab/DoubleElectron/DoubleElectron_vOriginal.root");
    myBabyMaker* baby1 = new myBabyMaker();
    baby1->ScanChain(chain1, "triggerStudy/DoubleElectron_triggerStudy_27Apr11_v1.root", true, -1);
    return;

    // Single Muon
    TChain *chain3 = new TChain("Events");
    //chain3->Add("/nfs-4/userdata/cms2/SingleMu_Run2011A-PromptReco-v1_AOD/V04-00-13/*.root");
    chain3->Add("/home/users/dbarge/ntuple_production/devel_CMSSW_4_1_2_patch1_V04-01-03/crab/SingleMu/F41C615E-5C50-E011-99CA-0030487CD716.root");
    myBabyMaker* baby3 = new myBabyMaker();
    baby3->ScanChain(chain3, "triggerStudy/SingleMu_triggerStudy_27Apr11_v1.root", true, -1);

    // Double Muon
    TChain *chain2 = new TChain("Events");
    //chain2->Add("/hadoop/cms/store/user/yanjuntu/CMSSW_4_1_2_patch1_V04-00-13/DoubleMu_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/*.root");
    chain2->Add("/home/users/dbarge/ntuple_production/devel_CMSSW_4_1_2_patch1_V04-01-03/crab/DoubleMu/E4DB0E07-E355-E011-BFE6-003048F11C58.root");
    myBabyMaker* baby2 = new myBabyMaker();
    baby2->ScanChain(chain2, "triggerStudy/DoubleMu_triggerStudy_27Apr11_v1.root", true, -1);

    return;

  // MC

//////////
// 2010 //
//////////

// Data

  // EG
  TChain *chain1 = ChainFromText( "input_data/eg_uaf.txt" );
  //TChain *chain1 = new TChain("Events");
  //chain1->Add("/nfs-3/userdata/cms2/Electron_Run2010B-PromptReco-v2_RECO/V03-06-14/merged_ntuple_149294_1.root");
  myBabyMaker * baby1 = new myBabyMaker();
  baby1->ScanChain(chain1, "EG_2010.root", true, -1);

  // 2011 testing
  //TChain *chain1 = new TChain("Events");
  //chain1->Add("/nfs-3/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/*.root");
  //myBabyMaker * baby1 = new myBabyMaker();
  //baby1->ScanChain(chain1, "ttbar.root", false, -1);
  //return;

  // Mu
  TChain *chain2 = ChainFromText( "input_data/mu_uaf.txt" );
  myBabyMaker * baby2 = new myBabyMaker();
  baby2->ScanChain(chain2, "Mu.root", true, -1);

  // EGMon
  TChain *chain3 = ChainFromText( "input_data/egmon_uaf.txt" );
  myBabyMaker * baby3 = new myBabyMaker();
  baby3->ScanChain(chain3, "EGMon.root", true, -1);

// Monte Carlo

  // QCD 30to50
  TChain *chain4 = ChainFromText( "input_data/qcd_pt_30to50_fall10_uaf.txt" );
  myBabyMaker * baby4 = new myBabyMaker();
  baby4->ScanChain(chain4, "qcd_pt_30to50_fall10.root", false, -1);

  // QCD 50to80
  TChain *chain5 = ChainFromText( "input_data/qcd_pt_50to80_fall10_uaf.txt" );
  myBabyMaker * baby5 = new myBabyMaker();
  baby5->ScanChain(chain5, "qcd_pt_50to80_fall10.root", false, -1);

  // QCD 80to120
  TChain *chain6 = ChainFromText( "input_data/qcd_pt_80to120_fall10_uaf.txt" );
  myBabyMaker * baby6 = new myBabyMaker();
  baby6->ScanChain(chain6, "qcd_pt_80to120_fall10.root", false, -1);

  /*
  // QCD 30
  TChain *chain4 = ChainFromText( "input_data/qcd30_uaf.txt" );
  myBabyMaker * baby4 = new myBabyMaker();
  baby4->ScanChain(chain4, "qcd30.root", false, -1);

  // QCD 50
  TChain *chain5 = ChainFromText( "input_data/qcd50_uaf.txt" );
  myBabyMaker * baby5 = new myBabyMaker();
  baby5->ScanChain(chain5, "qcd50.root", false, -1);

  // QCD 80
  TChain *chain6 = ChainFromText( "input_data/qcd80_uaf.txt" );
  myBabyMaker * baby6 = new myBabyMaker();
  baby6->ScanChain(chain6, "qcd80.root", false, -1);
  */

  // MuEnriched 10
  TChain *chain7 = ChainFromText( "input_data/mu10_uaf.txt" );
  myBabyMaker * baby7 = new myBabyMaker();
  baby7->ScanChain(chain7, "mu10.root", false, -1);

  // MuEnriched 15
  TChain *chain8 = ChainFromText( "input_data/mu15_uaf.txt" );
  myBabyMaker * baby8 = new myBabyMaker();
  baby8->ScanChain(chain8, "mu15.root", false, -1);

  // Wmunu
//  TChain *chain6 = ChainFromText( "input_data/wmunu_uaf.txt" );
//  myBabyMaker * baby6 = new myBabyMaker();
//  baby6->ScanChain(chain6, "Wmunu.root", false, -1);

  // Wenu
//  TChain *chain7 = ChainFromText( "input_data/wenu_uaf.txt" );
//  myBabyMaker * baby7 = new myBabyMaker();
//  baby7->ScanChain(chain7, "Wenu.root", false, -1);

  // InclusiveMu 15
//  TChain *chain8 = ChainFromText( "input_data/inclMu_uaf.txt" );
//  myBabyMaker * baby8 = new myBabyMaker();
//  baby8->ScanChain(chain8, "inclMu.root", false, -1);

}
