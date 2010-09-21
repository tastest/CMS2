#ifndef __CINT__
#include "TChain.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include "histtools.h"
#include "ossusy_looper.h"

#include <iostream>
#endif

void pickSkimIfExists( TChain *ch, const std::string& base, const std::string& skimPrefix )
{
  TChain *dummy = new TChain("Events");

  /*
    if (skimPrefix != "") {
    std::string skimName = Form("data/skim_v2/%s_skimmednTuple.root", skimPrefix.c_str());
    if (dummy->Add(skimName.c_str())) {
    int nFiles = ch->Add(skimName.c_str());
    std::cout << "Skim " << skimName.c_str() << " exists: use it. Loaded "<< nFiles << " files" << std::endl;
    return;
    } else
    std::cout << "Skim " << skimName.c_str() << " does not exist: use " << base.c_str() << std::endl;
    }
  */
  if (dummy->Add(base.c_str())) {
    int nFiles = ch->Add(base.c_str());
    std::cout << "Main " <<base.c_str() << " exists: use it. Loaded " 
              << nFiles << " files" << std::endl;
  } else
    std::cout << "FUCK SHIT DAMN!" << std::endl;

  // be paranoid
  if (nFiles == 0) {
    std::cout << "ERROR: expected to read files " 
              << base.c_str() << "  but found none" << std::endl;
    assert(0);
  }

  return;
}

void doAll_ossusy_looper(bool skipFWLite = true)
{

  //Load CORE stuff
  gROOT->ProcessLine(".L CORE/CMS2.cc+");
  gROOT->ProcessLine(".L CORE/trackSelections.cc+");
  gROOT->ProcessLine(".L CORE/metSelections.cc+");
  gROOT->ProcessLine(".L CORE/eventSelections.cc+");
  gROOT->ProcessLine(".L CORE/electronSelectionsParameters.cc+");
  gROOT->ProcessLine(".L CORE/electronSelections.cc+");
  gROOT->ProcessLine(".L CORE/muonSelections.cc+");
  gROOT->ProcessLine(".L CORE/SimpleFakeRate.cc+");
  gROOT->ProcessLine(".L CORE/mcSelections.cc+");
  gROOT->ProcessLine(".L CORE/MT2/MT2.cc+");


  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)", skipFWLite));

  // Load FWLite
  gSystem->Load("Tools/MiniFWLite/libMiniFWLite.so");

  // Load and compile the looping code
  gSystem->CompileMacro("ossusy_looper.C","++k", "libossusy_looper");

  ossusy_looper* looper = new ossusy_looper();
  //use OS/SS baseline selection as documented in:
  //http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/SusyStudies3x
  looper->set_susybaseline(0);
  //make baby ntuple
  looper->set_createTree(1);
  //use bitmask selection
  looper->set_useBitMask(0);


  // K-factors
  // these have been k-factors NLO/LO before
  // now using them as sample normalizations to NLO

  // these two are taken from Ceballos's pdf. 
  // It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.;  // 375pb, 127000 events processed
  float kttotr    = 1.;  // 375pb, 127000 events processed
  float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjets    = 1.;  // 11850 pb, 980000 events processed
  float kWcharm   = 1.1;
  float kZjets    = 1.;  
  float kDYee     = 1.;  // 1230 pb,  970360 events processed
  float kDYmm     = 1.;  // 1230 pb,  970360 events processed
  float kDYtautau = 1.;  // 1230 pb,  970360 events processed
  float kppMuX    = 1.;  // xsec/nevents
  float kEM       = 1.;
  float ktW       = 1.;  // the evtScale is all negative for some reason
  float kVQQ      = 1.;
  float kLM0      = 1.;
  float kLM1      = 1.;
  float kLM2      = 1.;
  float kLM3      = 1.;
  float kLM4      = 1.;
  float kLM5      = 1.;
  float kLM6      = 1.;
  float kLM7      = 1.;
  float kLM8      = 1.;
  float kLM9      = 1.;
  float kLM10     = 1.;
  float kLM11     = 1.;
  float kLM12     = 1.;
  float kLM13     = 1.;
  float kLMscan   = 1.;
  float kML1      = 1.;
  float kML2      = 1.;
  float kML3      = 1.;
  float kML4      = 1.;
  float kML5      = 1.;
  float kML6      = 1.;
  float kML7      = 1.;
  float kML8      = 1.;

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;
  int preWW       = 1;
  int preWZ       = 1;
  int preZZ       = 1;
  int preWjets    = 1;
  int preWcharm   = 1;
  int preZjets    = 1;
  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;
  int preppMuX    = 1;
  int preEM       = 1;
  int pretW       = 1;
  int preVQQ      = 1;
  int preLM0      = 1;
  int preLM1      = 1;
  int preLM2      = 1;
  int preLM3      = 1;
  int preLM4      = 1;
  int preLM5      = 1;
  int preLM6      = 1;
  int preLM7      = 1;
  int preLM8      = 1;
  int preLM9      = 1;
  int preLM10     = 1;
  int preLM11     = 1;
  int preLM12     = 1;
  int preLM13     = 1;
  int preML1      = 1;
  int preML2      = 1;
  int preML3      = 1;
  int preML4      = 1;
  int preML5      = 1;
  int preML6      = 1;
  int preML7      = 1;
  int preML8      = 1;
  int preLMscan = 1;
  
  //Flags for files to run over
  bool rundata     = 1;
  bool runttdil    = 1;
  bool runttotr    = 1;
  bool runWW       = 1;
  bool runWZ       = 1;
  bool runZZ       = 1;
  bool runWjets    = 1;
  bool runWcharm   = 0;
  bool runZjets    = 1;
  bool runDYee     = 0;
  bool runDYmm     = 0;
  bool runDYtautau = 0;
  bool runppMuX    = 0;
  bool runEM       = 0;
  bool runtW       = 1;
  bool runVQQ      = 0;
  bool runLM0      = 1;
  bool runLM1      = 1;
  bool runLM2      = 1;
  bool runLM3      = 1;
  bool runLM4      = 1;
  bool runLM5      = 1;
  bool runLM6      = 1;
  bool runLM7      = 1;
  bool runLM8      = 1;
  bool runLM9      = 1;
  bool runLM10     = 1;
  bool runLM11     = 1;
  bool runLM12     = 1;
  bool runLM13     = 1;
  bool runML1      = 0;
  bool runML2      = 0;
  bool runML3      = 0;
  bool runML4      = 0;
  bool runML5      = 0;
  bool runML6      = 0;
  bool runML7      = 0;
  bool runML8      = 0;
  bool runLMscan   = 0; 

/*  
  //Flags for files to run over
  bool rundata     = 1;
  bool runttdil    = 0;
  bool runttotr    = 0;
  bool runWW       = 0;
  bool runWZ       = 0;
  bool runZZ       = 0;
  bool runWjets    = 0;
  bool runWcharm   = 0;
  bool runZjets    = 0;
  bool runDYee     = 0;
  bool runDYmm     = 0;
  bool runDYtautau = 0;
  bool runppMuX    = 0;
  bool runEM       = 0;
  bool runtW       = 0;
  bool runVQQ      = 0;
  bool runLM0      = 0;
  bool runLM1      = 0;
  bool runLM2      = 0;
  bool runLM3      = 0;
  bool runLM4      = 0;
  bool runLM5      = 0;
  bool runLM6      = 0;
  bool runLM7      = 0;
  bool runLM8      = 0;
  bool runLM9      = 0;
  bool runLM10     = 0;
  bool runLM11     = 0;
  bool runLM12     = 0;
  bool runLM13     = 0;
  bool runML1      = 0;
  bool runML2      = 0;
  bool runML3      = 0;
  bool runML4      = 0;
  bool runML5      = 0;
  bool runML6      = 0;
  bool runML7      = 0;
  bool runML8      = 0;
  bool runLMscan   = 0; 
*/

  TChain* chdata = new  TChain("Events");
  if(rundata){
    pickSkimIfExists(chdata, 
                     "/tas/fgolf/SSskims/data/dilep_skim_2.root",
                     "data");
    pickSkimIfExists(chdata, 
                     "/tas/fgolf/SSskims/data/skimmed_ntuple*.root",
                     "data");
    pickSkimIfExists(chdata, 
                     "/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed_ntuple*.root",
                     "data");
    pickSkimIfExists(chdata, 
                     "/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-12/diLepPt1020Skim/skimmed_ntuple*.root",
                     "data");
    pickSkimIfExists(chdata, 
                     "/tas/cms2/EG_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139783_0.root",
                     "data");
    pickSkimIfExists(chdata, 
                     "/tas/cms2/Mu_Run2010A-PromptReco-v4_RECO/V03-04-26-02/singleLepPt10Skim/skimmed_ntuple_139783_0.root",
                     "data");
  }
  
  TChain* chZjets = new  TChain("Events");
  if(runZjets){
    pickSkimIfExists(chZjets, 
                     "/tas/cms2/ZJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-08/merged*root",
                     "Zjets");
  }

  TChain* chtopdil = new TChain("Events");
  if (runttdil) {
    pickSkimIfExists(chtopdil, 
                     "/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/merged*.root",
                     "TTJets");
  }

  TChain* chtopotr = new TChain("Events");
  if(runttotr){
    pickSkimIfExists(chtopotr, 
                     "/tas/cms2/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/merged*.root",
                     "TTJets");
  }

  TChain* chww = new TChain("Events");
  if(runWW){
    pickSkimIfExists(chww, 
                     "/tas/cms2/WW_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/merged*root",
                     "WW");
  }

  TChain* chWZ = new TChain("Events");
  if(runWZ){
    pickSkimIfExists(chWZ, 
                     "/tas/cms2/WZ_Spring10-START3X_V26_S09-v1/V03-04-08/merged*root",
                     "WZ_incl"); // can try WZ_3l-Pythia
  }

  TChain* chZZ = new TChain("Events");
  if(runZZ){
    pickSkimIfExists(chZZ, 
                     "/tas/cms2/ZZ_Spring10-START3X_V26_S09-v1_DiLep/V03-04-08/merged*.root", 
                     "ZZ");
  }

  TChain* chWjets = new  TChain("Events");
  if(runWjets){
    pickSkimIfExists(chWjets, 
                     "/tas/cms2/WJets-madgraph_Spring10-START3X_V26_S09-v1_SingleLep/V03-04-08/merged*.root",
                     "WJets");
  }

  TChain* chWcharm = new TChain("Events");
  if(runWcharm){
    pickSkimIfExists(chWcharm, 
                     //is this the correct sample????
   		     "/tas/cms2/WCJets_7TeV-madgraph_Spring10-START3X_V26-v1/V03-04-13-01/merged*root",
                     "Wc");
  }
  
  TChain* chDYtautau = new  TChain("Events");
  if(runDYtautau){
    pickSkimIfExists(chDYtautau, 
                     "/tas/cms2/Ztautau_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged*.root", 
                     "DYtautau");
  }
  
  // the low-mass splice has no choice other than the skim
  //    pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Ztautau_M20");

  TChain* chDYee = new  TChain("Events");
  if(runDYee){
    pickSkimIfExists(chDYee, 
                     "/tas/cms2/Zee_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged*.root", 
                     "DYee");
  }
  // the low-mass splice has no choice other than the skim
  //    pickSkimIfExists(chDYee, "data/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Zee_M20");
  
  if(runDYmm){
    TChain* chDYmm = new  TChain("Events");
    pickSkimIfExists(chDYmm, 
                     "/tas/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged*.root", 
                     "DYmm");
  }
  // the low-mass splice has no choice other than the skim
  //    pickSkimIfExists(chDYmm, "data/Zmumu_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root", "Zmumu_M20");
  
  
  // ppMuX
  TChain* chppMuX = new  TChain("Events");
  if (runppMuX) {
    pickSkimIfExists(chppMuX, 
                     "/tas/cms2/InclusiveMu15_Spring10-START3X_V26_S09-v1/V03-04-13-07/",
                     "InclusiveMuPt15"); 
    // can try InclusiveMu5Pt50 .. figure out how to merge later
  }
  
  // ppEM
  //only pt80to170 sample is currently available!!!
  TChain* chEM =  new  TChain("Events");
  if (runEM) {
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
    pickSkimIfExists(chEM, 
                     "/tas/cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
                     "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
//     pickSkimIfExists(chEM, 
//                      "data3x/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", 
//                      "_skimSimple2020");
  }

  // tW
  TChain* chtW = new  TChain("Events");
  if (runtW) {
    pickSkimIfExists(chtW, 
                     "/tas/cms2/SingleTop_sChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",
                     "SingleTop_sChannel"); 
    pickSkimIfExists(chtW, 
                     "/tas/cms2/SingleTop_tChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-13-07/merged*.root",
                     "SingleTop_tChannel"); 
    pickSkimIfExists(chtW, 
                     "/tas/cms2/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/V03-04-07/merged*.root",
                     "SingleTop_tWChannel"); 
  }

  // VQQ
  TChain* chVQQ = new TChain("Events");
  if (runVQQ) {
    pickSkimIfExists(chVQQ, 
                     //is this the correct sample???
                     "/tas/cms2/VqqJets-madgraph_Spring10-START3X_V26_S09-v1/merged*.root",
  		     //"data/VQQ-madgraph_Summer08_IDEAL_V11_redigi_v2/merged*.root", 
  		     "VQQ");
  }
  
  // LM points currently being processed!!!!
  // LM0
  TChain *chLM0 = new TChain("Events");
  if (runLM0) {
    pickSkimIfExists(chLM0, 
                     "/tas/cms2/LM0_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM0");
  }
  
  // LM1
  TChain *chLM1 = new TChain("Events");
  if (runLM1) {
    pickSkimIfExists(chLM1, 
                     "/tas/cms2/LM1_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM1");
  }

  // LM2
  TChain *chLM2 = new TChain("Events");
  if (runLM2) {
    pickSkimIfExists(chLM2, 
                     "/tas/cms2/LM2_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM2");
  }

  // LM3
  TChain *chLM3 = new TChain("Events");
  if (runLM3) {
    pickSkimIfExists(chLM3, 
                     "/tas/cms2/LM3_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM3");
  }

  // LM4
  TChain *chLM4 = new TChain("Events");
  if (runLM4) {
    pickSkimIfExists(chLM4, 
                     "/tas/cms2/LM4_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM4");
  }

  // LM5
  TChain *chLM5 = new TChain("Events");
  if (runLM5) {
    pickSkimIfExists(chLM5, 
                     "/tas/cms2/LM5_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM5");
  }

  // LM6
  TChain *chLM6 = new TChain("Events");
  if (runLM6) {
    pickSkimIfExists(chLM6, 
                     "/tas/cms2/LM6_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM6");
  }

  // LM7
  TChain *chLM7 = new TChain("Events");
  if (runLM7) {
    pickSkimIfExists(chLM7, 
                     "/tas/cms2/LM7_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM7");
  }

  // LM8
  TChain *chLM8 = new TChain("Events");
  if (runLM8) {
    pickSkimIfExists(chLM8, 
                     "/tas/cms2/LM8_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM8");
  }

  // LM9
  TChain *chLM9 = new TChain("Events");
  if (runLM9) {
    pickSkimIfExists(chLM9, 
                     "/tas/cms2/LM9_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM9");
  }

  // LM10
  TChain *chLM10 = new TChain("Events");
  if (runLM10) {
    pickSkimIfExists(chLM10, 
                     "/tas/cms2/LM10_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM10");
  }

  // LM11 
  TChain *chLM11 = new TChain("Events");
  if (runLM11) {
    pickSkimIfExists(chLM11, 
                     "/tas/cms2/LM11_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM11");
  }

  // LM12
  TChain *chLM12 = new TChain("Events");
  if (runLM12) {
    pickSkimIfExists(chLM12, 
                     "/tas/cms2/LM12_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM12");
  }

  // LM13
  TChain *chLM13 = new TChain("Events");
  if (runLM13) {
    pickSkimIfExists(chLM13, 
                     "/tas/cms2/LM13_Spring10-START3X_V26_S09-v1/V03-04-13-01/merged*.root", 
                     "SUSY_LM13");
  }
  
  // ML1
  TChain *chML1 = new TChain("Events");
  if (runML1) {
    pickSkimIfExists(chML1, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML01_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML1");
  }

  // ML2
  TChain *chML2 = new TChain("Events");
  if (runML2) {
    pickSkimIfExists(chML2, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML02_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML2");
  }

  // ML3
  TChain *chML3 = new TChain("Events");
  if (runML3) {
    pickSkimIfExists(chML3, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML03_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML3");
  }

  // ML4
  TChain *chML4 = new TChain("Events");
  if (runML4) {
    pickSkimIfExists(chML4, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML04_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML4");
  }

  // ML5
  TChain *chML5 = new TChain("Events");
  if (runML5) {
    pickSkimIfExists(chML5, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML05_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML5");
  }

  // ML6
  TChain *chML6 = new TChain("Events");
  if (runML6) {
    pickSkimIfExists(chML6, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML06_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML6");
  }

  // ML7
  TChain *chML7 = new TChain("Events");
  if (runML7) {
    pickSkimIfExists(chML7, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML07_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML7");
  }
  
  // ML8
  TChain *chML8 = new TChain("Events");
  if (runML8) {
    pickSkimIfExists(chML8, 
                     "/tas/cms2/PhysicsProcess_PYTHIA6_SUSY_GMSM_SC_ML08_7TeV_v0/V03-04-13-01-gmsb/merged*root",
                     "SUSY_ML8");
  }

  // LMscan
  TChain *chLMscan = new TChain("Events");
  if (runLMscan) {
    pickSkimIfExists(chLMscan, 
                     "data3x/TANB3_CMSW336FASTv3/V03-00-37/merged*.root",
                     "LMscan");
  }


  //--------------------------------
  //set luminosity to scale to
  //--------------------------------
  float lumi              = 2.7945e-3;     //2.8 pb-1
  bool  calculateTCMET    = true; //redo tcmet calculation on the fly

  char* jetTypeStrings[3] = {"JPT", "calo","pfjet"};
  char* metTypeStrings[4] = {"tcmet", "muon", "muonjes","pfmet"};
  char* zvetoStrings[4]   = {"", "_allzveto", "_nozveto","_selectz"};
  bool doFakeApp          = false;

  // Process files one at a time, and color them as needed
  for (int jetTypeIdx = 2; jetTypeIdx < 3; ++jetTypeIdx)
    {
      for (int metTypeIdx = 0; metTypeIdx < 1; ++metTypeIdx)
        {
          for (int zvetoIdx = 0; zvetoIdx < 1; ++zvetoIdx)
            {

              ossusy_looper::JetTypeEnum  jetType(jetTypeIdx);
              ossusy_looper::MetTypeEnum  metType(metTypeIdx);
              ossusy_looper::ZVetoEnum    zveto(zvetoIdx);

              if (rundata) {
                cout << "Processing data" << endl;
                looper->ScanChain(chdata,"data", 1, 1, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing data" << endl;
                hist::color("data", kBlack);
              }
              if (runZjets) {
                cout << "Processing Zjets" << endl;
                looper->ScanChain(chZjets,"Zjets", kZjets, preZjets, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing Zjets" << endl;
                hist::color("Zjets", kBlack);
              }
              if (runttdil) {
                cout << "Processing ttbar dileptonic.. " << endl;
                looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ttbar dileptonic.. " << endl;
                hist::color("ttdil", kYellow);
              }
              if (runttotr) {
                cout << "Processing ttbar no-dileptons.. " << endl;
                looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ttbar no-dileptons.. " << endl;
                hist::color("ttotr", 30);
              }
              if (runWW) {
                cout << "Processing WW.." << endl;
                looper->ScanChain(chww,"ww", kWW, preWW, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing WW.." << endl;
                hist::color("ww", kRed);
              }
              if (runWZ) {
                cout << "Processing WZ.." << endl;
                looper->ScanChain(chWZ,"wz", kWZ, preWZ, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing WZ.." << endl;
                hist::color("wz", kBlue);
              }
              if (runZZ) {
                cout << "Processing ZZ.." << endl;
                looper->ScanChain(chZZ,"zz", kZZ, preZZ, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ZZ.." << endl;
                hist::color("zz", kGreen);
              }
              if (runWjets) {
                cout << "Processing Wjets.." << endl;
                looper->ScanChain(chWjets,"wjets", kWjets, preWjets, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing Wjets.." << endl;
                hist::color("wjets", 40);
              }
              if (runWcharm) {
                cout << "Processing Wcharm.." << endl;
                looper->ScanChain(chWcharm, "wcharm", kWcharm, preWcharm, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing Wcharm.." << endl;
                hist::color("wcharm", 50);
              }
              if (runDYtautau) {
                cout << "Processing DY->tautau" << endl;
                looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing DY->tautau" << endl;
                hist::color("DYtautau", kBlack);
              }
              if (runDYee) {
                cout << "Processing DY->ee" << endl;
                looper->ScanChain(chDYee,"DYee", kDYee, preDYee, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done rocessing DY->ee" << endl;
                hist::color("DYee", kMagenta);
              }
              if (runDYmm) {
                cout << "Processing DY->mm" << endl;
                looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing DY->mm" << endl;
                hist::color("DYmm", kCyan);
              }
              if (runppMuX) {
                cout << "Processing ppMuX" << endl;
                looper->ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ppMuX" << endl;
                hist::color("ppMuX", 51);
              }
              if (runEM) {
                cout << "Processing EM" << endl;
                looper->ScanChain(chEM,"EM", kEM, preEM, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing EM" << endl;
                hist::color("EM", 49);
              }
              if (runtW) {
                cout << "Processing tW" << endl;
                looper->ScanChain(chtW,"tW", ktW, pretW, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing tW" << endl;
                hist::color("tW", 63);
              }
              if (runVQQ) { 
                cout << "Processing VQQ" << endl;
                looper->ScanChain(chVQQ,"VQQ", kVQQ, preVQQ, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing VQQ" << endl;
                hist::color("VQQ", 45);
              }
              if (runLM0) {
                cout << "Processing LM0" << endl;
                looper->ScanChain(chLM0, "LM0", kLM0, preLM0, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM0" << endl;
                hist::color("LM0", kOrange);
              }
              if (runLM1) {
                cout << "Processing LM1" << endl;
                looper->ScanChain(chLM1, "LM1", kLM1, preLM1, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM1" << endl;
                hist::color("LM1", kOrange+1);
              }
              if (runLM2) {
                cout << "Processing LM2" << endl;
                looper->ScanChain(chLM2, "LM2", kLM2, preLM2, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM2" << endl;
                hist::color("LM2", kOrange+2);
              }
              if (runLM3) {
                cout << "Processing LM3" << endl;
                looper->ScanChain(chLM3, "LM3", kLM3, preLM3, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM3" << endl;
                hist::color("LM3", kOrange+3);
              }
              if (runLM4) {
                cout << "Processing LM4" << endl;
                looper->ScanChain(chLM4, "LM4", kLM4, preLM4, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM4" << endl;
                hist::color("LM4", kOrange+4);
              }
              if (runLM5) {
                cout << "Processing LM5" << endl;
                looper->ScanChain(chLM5, "LM5", kLM5, preLM5, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM5" << endl;
                hist::color("LM5", kOrange+5);
              }
              if (runLM6) {
                cout << "Processing LM6" << endl;
                looper->ScanChain(chLM6, "LM6", kLM6, preLM6, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM6" << endl;
                hist::color("LM6", kOrange+6);
              }
              if (runLM7) {
                cout << "Processing LM7" << endl;
                looper->ScanChain(chLM7, "LM7", kLM7, preLM7, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM7" << endl;
                hist::color("LM7", kOrange+7);
              }
              if (runLM8) {
                cout << "Processing LM8" << endl;
                looper->ScanChain(chLM8, "LM8", kLM8, preLM8, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM8" << endl;
                hist::color("LM8", kOrange+8);
              }
              if (runLM9) {
                cout << "Processing LM9" << endl;
                looper->ScanChain(chLM9, "LM9", kLM9, preLM9, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM9" << endl;
                hist::color("LM9", kOrange+9);
              }
              if (runLM10) {
                cout << "Processing LM10" << endl;
                looper->ScanChain(chLM10, "LM10", kLM10, preLM10, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM10" << endl;
                hist::color("LM10", kOrange+10);
              }
              if (runLM11) { 
                cout << "Processing LM11" << endl;
                looper->ScanChain(chLM11, "LM11", kLM11, preLM11, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM11" << endl;
                hist::color("LM11", kOrange-7);
              }
              if (runLM12) {
                cout << "Processing LM12" << endl;
                looper->ScanChain(chLM12, "LM12", kLM12, preLM12, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM12" << endl;
                hist::color("LM12", kOrange-7);
              }
              if (runLM13) {
                cout << "Processing LM13" << endl;
                looper->ScanChain(chLM13, "LM13", kLM13, preLM13, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LM13" << endl;
                hist::color("LM13", kOrange-7);
              }
              if (runML1) {
                cout << "Processing ML1" << endl;
                looper->ScanChain(chML1, "ML1", kML1, preML1, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML1" << endl;
              }
              if (runML2) {
                cout << "Processing ML2" << endl;
                looper->ScanChain(chML2, "ML2", kML2, preML2, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML2" << endl;
              }
              if (runML3) {
                cout << "Processing ML3" << endl;
                looper->ScanChain(chML3, "ML3", kML3, preML3, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML3" << endl;
              }
              if (runML4) {
                cout << "Processing ML4" << endl;
                looper->ScanChain(chML4, "ML4", kML4, preML4, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML4" << endl;
              }
              if (runML5) {
                cout << "Processing ML5" << endl;
                looper->ScanChain(chML5, "ML5", kML5, preML5, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML5" << endl;
              }
              if (runML6) {
                cout << "Processing ML6" << endl;
                looper->ScanChain(chML6, "ML6", kML6, preML6, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML6" << endl;
              }
              if (runML7) {
                cout << "Processing ML7" << endl;
                looper->ScanChain(chML7, "ML7", kML7, preML7, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML7" << endl;
              }
              if (runML8) {
                cout << "Processing ML8" << endl;
                looper->ScanChain(chML8, "ML8", kML8, preML8, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing ML8" << endl;
              }
              if (runLMscan) {
                cout << "Processing LMscan" << endl;
                looper->ScanChain(chLMscan, "LMscan", kLMscan, preLMscan, lumi, jetType, metType, zveto, doFakeApp, calculateTCMET);
                cout << "Done processing LMscan" << endl;
                hist::color("LMscan", kOrange-7);
              }


              // save all the histograms
              if(doFakeApp) {
                const char* outFile = Form("output/ossusy_%s_%s%s_FakeApp.root", 
                                           jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx],zvetoStrings[zvetoIdx]);
              }
              else {
                const char* outFile = Form("output/ossusy_%s_%s%s.root", 
                                           jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx],zvetoStrings[zvetoIdx]);
              }

              //const char* outFile = Form("victory_baseline_genmetgt50_nosumjetptcut_%s_%s_pleasework_varbins.root", 
              //jetTypeStrings[jetTypeIdx], metTypeStrings[metTypeIdx]);
              TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
              rootdir->cd();
              saveHist(outFile);
              deleteHistos();

            }//zvetoIdx
        } // metTypeIdx
    } // jetTypeIdx

  gSystem->Exit(0);
}
