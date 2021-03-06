//==============================================================
//
// This runs over the skimmed files (see AAREADME.txt)
//
// To run on unskimmed files, change the file names.
//
//
//
//==============================================================
#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include <iostream>
#include "histtools.h"
#include "ttDilCounts_looper.h"
#endif //__CINT__

void pickSkimIfExists(TChain* ch, const std::string& base, const std::string& skimExt){
  TChain* dummy = new TChain("Events");
  if (skimExt != ""){
    std::string skimName = base+skimExt;
    if (dummy->Add(skimName.c_str())){
      int nFiles = ch->Add(skimName.c_str());
      std::cout<<"Skim "<<skimName.c_str()<<" exists: use it. Loaded "<<nFiles<<" files"<<std::endl;
      return;
    } else {
      std::cout<<"Skim "<<skimName.c_str()<<" does not exist ==> will use "<<base.c_str()<<std::endl;
    }
  }
  int nFiles = ch->Add(base.c_str());
  std::cout<<"Main "<<base.c_str()<<" exists: use it. Loaded "<<nFiles<<" files"<<std::endl;
  //be a bit paranoid here
  if (nFiles == 0) {
    std::cout<<"ERROR: expected to read files "<<base.c_str()<<" \n\t but found none"<<std::endl;
    assert(0);
  }
  return;
  
}

void doAll(unsigned int bitmask, bool skipFWLite = false){
  //here is a list to the combinations of cuts useful for the analysis:
  // 1957888 -- baseline
  // 35512320 -- baseline using tcmet
  // 1695744 -- baseline without MET
  // 1433600 -- baseline without zveto
  // 1926144 -- baseline without tight iso (only loose iso)
  // 18735104 -- baseline JES-up
  // 10346496 -- baseline JES-down
  // 1171456 -- baseline without MET, without zveto
  // 2220032 -- baseline without MET, without zveto, using AN09/047 (v<=4) trigger selection

  // 1941504 -- baseline without duplicate removal

  
  //cut <-> bit mask
  //ID cuts               -> 2**0 (1)
  //Isolation cuts        -> 2**1 (2) (default is both legs are isolated. Using relative isolation, TRK+CALO)
  //                       2**8+2**1 (Require one-only hyp lepton to be isolated. In emu, that will mean that the muon
  //                                  will have rel iso > 0.92 and the el will have 0.6 < relIso < 0.92)
  //                       2**9+2**1 (Require one-only hyp lepton to be isolated. In emu, that will mean that the el
  //                                  will have rel iso > 0.92 and the mu will have 0.6 < relIso < 0.92)
  //                       2**9+2**8+2**1 (require both leptons to have 0.6 < relIso < 0.92)
  //DileptonMassVeto      -> 2**2 (4)
  //METcut                -> 2**3 (8)
  //nJets                 -> 2**4 (16)
  //extra MuTag           -> 2**5 (32)
  //METveto               -> 2**6 (64)
  //Extra MuTag (pt>5)    -> 2**7 (128)
  //looseDileptonSelection, TTdil note 2008
  //                      -> 2**10 (1024)
  //fullMultipleHypsOnly  -> 2**11 (2048) !!!!! Not implemented, so does nothing right now !!!!!!!
  //applyZWindow cut      -> 2**12 (4096)
  //Opp. Sign Selection   -> 2**13 (8192)
  //fillMaxWeightDilOnly  -> 2**14 (16284
  //leptonIsolationDilSelectionTTDil08 -> 2**15 (uses trk and calo isolation seperately, reltrkIso > 0.9, 
  //                                             relCaloIso > 0.9. NO OTHER CUTS BUT ISOLATION APPLIED)
  //looseDilSelectionNoIsoTTDil08      -> 2**16 (basic muon preselection cuts, no isolation)
  //lepton20Eta2p4DilSelection         -> 2**17 (only pt and eta cuts applied to leptons)
  //metBaselineSelectionTTDil108       -> 2**18 (calls passPatMet_OF20_SF30 -> corrMET  > 20 (emu),
  //                                             corrMET > 30 (ee, mumu))
  //dilepMassVetoCutTTDil08            -> 2**19 pretty much what it sounds like
  //applyTriggersMu9orLisoE15          -> 2**20 mm -> HLT_Mu9 
  //                                            ee -> HLT_Ele15_SW_L1R
  //                                            em -> HLT_Mu9 or HLT_Ele15_SW_L1R
  //applyTriggersTTDil08JanTrial       -> 2**21 mm -> HLT_Mu15_L1Mu7 || HLT_DoubleMu3
  //                                            ee -> HLT_IsoEle18_L1R || HLT_DoubleIsoEle12_L1R
  //                                            em -> HLT_IsoEle18_L1R || HLT_Mu15_L1Mu7 || HLT_IsoEle10_Mu10_L1R
  // dilepAdditionalMassVetoCutTTDil08 -> 2**22
  // corJES10ptDn                      -> 2**23 rescale JES by 10% Down
  // corJES10ptUp                      -> 2**24 rescale JES by 10% up
  // useTcMet                          -> 2**25 use tcmet
  // useJPT                            -> 2**26 use JPT
  // muJetClean                        -> 2**27 do not count jets within 0.4 of muons 
  // dilTruthMatch                     -> 2**28 require the lepton to be MCtruth matched (coming off a hard-scattering lepton)
 
  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)",skipFWLite));

  // Load and compile the looping code
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  

  // Flag for jet selection
  // true  = hyp_jet selection (15 GeV uncorrectedm eta<3)
  // false = 30 GeV corrected, eta<2.4
  bool oldjet=false;

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO
  
  //these two are taken from Ceballos's pdf. 
  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; //375pb, 127000 events processed
  float kttotr    = 1.; //375pb, 127000 events processed

  float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjets    = 1.; //11850 pb, 980000 events processed
  float kWcharm   = 1.1;
  float kDYee     = 1.;  //1230 pb,  970360 events processed
  float kDYmm     = 1.;  //1230 pb,  970360 events processed
  float kDYtautau = 1.;  //1230 pb,  970360 events processed
  float kppMuX    = 1.; //xsec/nevents
  float kEM       = 1.;
  float ktW       = 1.; //the evtScale is all negative for some reason
  float kVQQ      = 1;

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;
  int preWW       = 1;
  int preWZ       = 1;
  int preZZ       = 1;
  int preWjets    = 1;
  int preWcharm   = 1;
  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;
  int preppMuX    = 1;
  int preEM       = 1;
  int pretW       = 1;
  int preVQQ      = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = true;
  bool runWW       = true;
  bool runWZ       = true;
  bool runZZ       = true;
  bool runWjets    = true;
  bool runWcharm   = false;
  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;
  bool runppMuX    = true;
  bool runEM       = true;
  bool runtW       = true;
  bool runVQQ      = true;

  TChain* chtopdil = new TChain("Events");
  pickSkimIfExists(chtopdil, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020anydil");

  TChain* chtopotr = new TChain("Events");
  pickSkimIfExists(chtopotr, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020nodil");

  TChain* chww = new TChain("Events");
  pickSkimIfExists(chww, "data/WW_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");

  TChain* chWZ = new TChain("Events");
  pickSkimIfExists(chWZ, "data/WZ_incl_Summer08_IDEAL_V11_redigi_v1/merged*.root", ""); // can try WZ_3l-Pythia

  TChain* chZZ = new TChain("Events");
  pickSkimIfExists(chZZ, "data/ZZ_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");
  
  TChain* chWjets = new  TChain("Events");
  pickSkimIfExists(chWjets, "data/WJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");

  TChain* chWcharm = new TChain("Events");
  pickSkimIfExists(chWcharm, "data/Wc-madgraph_Fall08_IDEAL_V11_redigi_v1/merged*.root", "");

  TChain* chDYtautau = new  TChain("Events");
  pickSkimIfExists(chDYtautau, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020tautau");
  //the low-mass splice has no choice other than the skim
  pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");
  
  TChain* chDYee = new  TChain("Events");
  pickSkimIfExists(chDYee, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020ee");
  //the low-mass splice has no choice other than the skim
  pickSkimIfExists(chDYee, "data/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");

  TChain* chDYmm = new  TChain("Events");
  pickSkimIfExists(chDYmm, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020mm");
  //the low-mass splice has no choice other than the skim
  pickSkimIfExists(chDYmm, "data/Zmumu_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");
  
  //ppMuX
  TChain* chppMuX = new  TChain("Events");
  if (runppMuX) {
    pickSkimIfExists(chppMuX, "data/InclusiveMuPt15_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020"); 
    //can try InclusiveMu5Pt50 .. figure out how to merge later
  }
  
  //ppEM
  TChain* chEM =  new  TChain("Events");
  if (runEM) {
    pickSkimIfExists(chEM, "data/QCD_EMenriched_Pt20to30_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
    pickSkimIfExists(chEM, "data/QCD_EMenriched_Pt30to80_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
    pickSkimIfExists(chEM, "data/QCD_EMenriched_Pt80to170_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
    pickSkimIfExists(chEM, "data/QCD_BCtoE_Pt20to30_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");
    pickSkimIfExists(chEM, "data/QCD_BCtoE_Pt30to80_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");
    pickSkimIfExists(chEM, "data/QCD_BCtoE_Pt80to170_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");
  }

  //tW
  TChain* chtW = new  TChain("Events");
  if (runtW) {
    pickSkimIfExists(chtW, "data/SingleTop_sChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 
    pickSkimIfExists(chtW, "data/SingleTop_tChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 
    pickSkimIfExists(chtW, "data/SingleTop_tWChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 
  }

  //VQQ
  TChain* chVQQ = new TChain("Events");
  if (runVQQ) {
    pickSkimIfExists(chVQQ, "data/VQQ-madgraph_Fall08_IDEAL_V9_v1/merged*.root", "");
  }


  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    cout << "Done Processing ttbar dileptonic.. "<<endl;
    hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }
  if (runWW) {
    cout << "Processing WW.."<<endl;
    looper->ScanChain(chww,"ww", kWW, preWW, oldjet, bitmask);
    hist::color("ww", kRed);
  }
  if (runWZ) {
    cout << "Processing WZ.."<<endl;
    looper->ScanChain(chWZ,"wz", kWZ, preWZ, oldjet, bitmask);
    hist::color("wz", kBlue);
  }
  if (runZZ) {
    cout << "Processing ZZ.."<<endl;
    looper->ScanChain(chZZ,"zz", kZZ, preZZ, oldjet, bitmask);
    hist::color("zz", kGreen);
  }
  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
    looper->ScanChain(chWjets,"wjets", kWjets, preWjets, oldjet, bitmask);
    hist::color("wjets", 40);
  }

  if (runWcharm) {
    cout << "Processing Wcharm.." << endl;
    looper->ScanChain(chWcharm, "wcharm", kWcharm, preWcharm, oldjet, bitmask);
    hist::color("wcharm", 50);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    looper->ScanChain(chDYee,"DYee", kDYee, preDYee, oldjet, bitmask);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, oldjet, bitmask);
    hist::color("DYmm", kCyan);
  }
  if (runppMuX) {
    cout << "Processing ppMuX"<<endl;
    looper->ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, oldjet, bitmask);
    hist::color("ppMuX", 51);
  }
  if (runEM) {
    cout << "Processing EM"<<endl;
    looper->ScanChain(chEM,"EM", kEM, preEM, oldjet, bitmask);
    hist::color("EM", 49);
  }
  if (runtW) {
    cout << "Processing tW"<<endl;
    looper->ScanChain(chtW,"tW", ktW, pretW, oldjet, bitmask);
    hist::color("tW", 63);
  }
    
  if (runVQQ) {
    cout << "Processing VQQ"<<endl;
    looper->ScanChain(chVQQ,"VQQ", kVQQ, preVQQ, oldjet, bitmask);
    hist::color("VQQ", 45);
  }

  //save all the histograms
    
  const char* outFile = Form("myHist_testfix_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
void doAllCombined(unsigned int bitmask, bool skipFWLite = false){
  //here is a list to the combinations of cuts useful for the analysis:
  // 1957888 -- baseline
  // 35512320 -- baseline using tcmet
  // 1695744 -- baseline without MET
  // 1433600 -- baseline without zveto
  // 1926144 -- baseline without tight iso (only loose iso)
  // 18735104 -- baseline JES-up
  // 10346496 -- baseline JES-down
  // 1171456 -- baseline without MET, without zveto
  // 2220032 -- baseline without MET, without zveto, using AN09/047 (v<=4) trigger selection

  // 1941504 -- baseline without duplicate removal
  // 1139712 -- baseline without MET, without zveto, no tight-iso (loose only) == "loose leptons"
  // 1172480 -- baseline without MET, without zveto, with tight iso
  // 538828800 -- baseline with dil dispatch by the highest mass
  // 538042368 -- baseline without MET, without zveto, with dil dispatch by the highest mass
  // 538010624 -- baseline without MET, without zveto, no tight-iso (loose only), with dil dispatch by the highest mass
  // 538043392 -- baseline without MET, without zveto, with tight iso, with dil dispatch by the highest mass
  // 1075699712 -- baseline with dil dispatch by the highest mass
  // 1074913280 -- baseline without MET, without zveto, with dil dispatch by the highest pt
  // 1074881536 -- baseline without MET, without zveto, no tight-iso (loose only), with dil dispatch by the highest pt
  // 1074914304 -- baseline without MET, without zveto, with tight iso, with dil dispatch by the highest pt


  //cut <-> bit mask
  //ID cuts               -> 2**0 (1)
  //Isolation cuts        -> 2**1 (2) (default is both legs are isolated. Using relative isolation, TRK+CALO)
  //                       2**8+2**1 (Require one-only hyp lepton to be isolated. In emu, that will mean that the muon
  //                                  will have rel iso > 0.92 and the el will have 0.6 < relIso < 0.92)
  //                       2**9+2**1 (Require one-only hyp lepton to be isolated. In emu, that will mean that the el
  //                                  will have rel iso > 0.92 and the mu will have 0.6 < relIso < 0.92)
  //                       2**9+2**8+2**1 (require both leptons to have 0.6 < relIso < 0.92)
  //DileptonMassVeto      -> 2**2 (4)
  //METcut                -> 2**3 (8)
  //nJets                 -> 2**4 (16)
  //extra MuTag           -> 2**5 (32)
  //METveto               -> 2**6 (64)
  //Extra MuTag (pt>5)    -> 2**7 (128)
  //looseDileptonSelection, TTdil note 2008
  //                      -> 2**10 (1024)
  //fullMultipleHypsOnly  -> 2**11 (2048) !!!!! Not implemented, so does nothing right now !!!!!!!
  //applyZWindow cut      -> 2**12 (4096)
  //Opp. Sign Selection   -> 2**13 (8192)
  //fillMaxWeightDilOnly  -> 2**14 (16284
  //leptonIsolationDilSelectionTTDil08 -> 2**15 (uses trk and calo isolation seperately, reltrkIso > 0.9, 
  //                                             relCaloIso > 0.9. NO OTHER CUTS BUT ISOLATION APPLIED)
  //looseDilSelectionNoIsoTTDil08      -> 2**16 (basic muon preselection cuts, no isolation)
  //lepton20Eta2p4DilSelection         -> 2**17 (only pt and eta cuts applied to leptons)
  //metBaselineSelectionTTDil108       -> 2**18 (calls passPatMet_OF20_SF30 -> corrMET  > 20 (emu),
  //                                             corrMET > 30 (ee, mumu))
  //dilepMassVetoCutTTDil08            -> 2**19 pretty much what it sounds like
  //applyTriggersMu9orLisoE15          -> 2**20 mm -> HLT_Mu9 
  //                                            ee -> HLT_Ele15_SW_L1R
  //                                            em -> HLT_Mu9 or HLT_Ele15_SW_L1R
  //applyTriggersTTDil08JanTrial       -> 2**21 mm -> HLT_Mu15_L1Mu7 || HLT_DoubleMu3
  //                                            ee -> HLT_IsoEle18_L1R || HLT_DoubleIsoEle12_L1R
  //                                            em -> HLT_IsoEle18_L1R || HLT_Mu15_L1Mu7 || HLT_IsoEle10_Mu10_L1R
  // dilepAdditionalMassVetoCutTTDil08 -> 2**22
  // corJES10ptDn                      -> 2**23 rescale JES by 10% Down
  // corJES10ptUp                      -> 2**24 rescale JES by 10% up
  // useTcMet                          -> 2**25 use tcmet
  // useJPT                            -> 2**26 use JPT
  // muJetClean                        -> 2**27 do not count jets within 0.4 of muons 
  // dilTruthMatch                     -> 2**28 require the lepton to be MCtruth matched (coming off a hard-scattering lepton)
 
  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)",skipFWLite));

  // Load and compile the looping code
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  

  // Flag for jet selection
  // true  = hyp_jet selection (15 GeV uncorrectedm eta<3)
  // false = 30 GeV corrected, eta<2.4
  bool oldjet=false;

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO
  
  //these two are taken from Ceballos's pdf. 
  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; //375pb, 127000 events processed
  float kttotr    = 1.; //375pb, 127000 events processed

  float kVV       = 1.;
  float kWjets    = 1.;
  float kDYeemm   = 1.;
  float kDYtautau = 1.;
  float kQCD      = 1.;
  float kt        = 1.;
  float kVgamma   = 1.;
  float kLM0      = 1.;

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;
  int preVV       = 1;
  int preWjets    = 1;
  int preDYeemm   = 1;
  int preDYtautau = 1;
  int preQCD      = 1;
  int pret        = 1;
  int preVgamma   = 1;
  int preLM0      = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = true;
  bool runVV       = true;
  bool runWjets    = true;
  bool runDYeemm   = true;
  bool runDYtautau = true;
  bool runQCD      = true;
  bool runt        = true;
  bool runVgamma   = true;
  bool runLM0      = true;

  TChain* chtopdil = new TChain("Events");
  pickSkimIfExists(chtopdil, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020anydil");

  TChain* chtopotr = new TChain("Events");
  pickSkimIfExists(chtopotr, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020nodil");

  TChain* chVV = new TChain("Events");
  pickSkimIfExists(chVV, "data/WW_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");
  pickSkimIfExists(chVV, "data/WZ_incl_Summer08_IDEAL_V11_redigi_v1/merged*.root", ""); // can try WZ_3l-Pythia
  pickSkimIfExists(chVV, "data/ZZ_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");
  
  TChain* chWjets = new  TChain("Events");
  pickSkimIfExists(chWjets, "data/WJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");

  TChain* chDYtautau = new  TChain("Events");
  pickSkimIfExists(chDYtautau, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020tautau");
  //the low-mass splice has no choice other than the skim
  //  pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");
  pickSkimIfExists(chDYtautau, "data/VQQ-madgraph_Fall08_IDEAL_V9_v1/merged*.root", "");
  
  TChain* chDYeemm = new  TChain("Events");
  pickSkimIfExists(chDYeemm, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "");
  //the low-mass splice has no choice other than the skim
  //  pickSkimIfExists(chDYeemm, "data/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");
  //  pickSkimIfExists(chDYeemm, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020mm");
  //the low-mass splice has no choice other than the skim
  //  pickSkimIfExists(chDYeemm, "data/Zmumu_M20_Summer08_IDEAL_V11_redigi_v1/merged*.root_skimSimple2020_20m50", "");
  pickSkimIfExists(chDYeemm, "data/VQQ-madgraph_Fall08_IDEAL_V9_v1/merged*.root", "");
  
  //ppMuX
  TChain* chQCD = new  TChain("Events");
  pickSkimIfExists(chQCD, "data/InclusiveMuPt15_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020"); 
  //  pickSkimIfExists(chQCD, "data/QCD_EMenriched_Pt20to30_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
  //  pickSkimIfExists(chQCD, "data/QCD_EMenriched_Pt30to80_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
  //  pickSkimIfExists(chQCD, "data/QCD_EMenriched_Pt80to170_Summer08_IDEAL_V11_redigi_v2-SingleLepton/merged*.root", "_skimSimple2020");
  //  pickSkimIfExists(chQCD, "data/QCD_BCtoE_Pt20to30_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");
  //  pickSkimIfExists(chQCD, "data/QCD_BCtoE_Pt30to80_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");
  //  pickSkimIfExists(chQCD, "data/QCD_BCtoE_Pt80to170_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged*.root", "_skimSimple2020");

  //tW
  TChain* cht = new  TChain("Events");
  pickSkimIfExists(cht, "data/SingleTop_sChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 
  pickSkimIfExists(cht, "data/SingleTop_tChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 
  pickSkimIfExists(cht, "data/SingleTop_tWChannel_Summer08_IDEAL_V11_redigi_v3/merged*.root", ""); 

  //Vgamma
  TChain* chVgamma = new TChain("Events");
  //  pickSkimIfExists(chVgamma, "data/AVJets-madgraph_Fall08_IDEAL_V9_v3/merged*.root", "_skimSimple2020");

  //LMs are run and loaded at the same time

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    cout << "Done Processing ttbar dileptonic.. "<<endl;
    hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }
  if (runVV) {
    cout << "Processing VV.."<<endl;
    looper->ScanChain(chVV,"VV", kVV, preVV, oldjet, bitmask);
    hist::color("VV", kRed);
  }
  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
    looper->ScanChain(chWjets,"wjets", kWjets, preWjets, oldjet, bitmask);
    hist::color("wjets", 40);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYeemm) {
    cout << "Processing DY->ee/mm" << endl;
    looper->ScanChain(chDYeemm,"DYeemm", kDYeemm, preDYeemm, oldjet, bitmask);
    hist::color("DYeemm", kMagenta);
  }
  if (runQCD) {
    cout << "Processing ppMuX and EM"<<endl;
    looper->ScanChain(chQCD,"QCD", kQCD, preQCD, oldjet, bitmask);
    hist::color("QCD", 51);
  }
  if (runt) {
    cout << "Processing t"<<endl;
    looper->ScanChain(cht,"t", kt, pret, oldjet, bitmask);
    hist::color("t_", 63);
  }

  if (runVgamma){
    cout << "Processing Vgamma ... "<<endl;
    looper->ScanChain(chVgamma, "Vgamma", kVgamma, preVgamma, oldjet, bitmask);
  }
    
  if (runLM0){
    std::vector<TString> lmEs;
    std::vector<TString> lmEds;
    lmEs.push_back("LM0"); lmEds.push_back("SUSY_LM0-sftsht_Summer08_IDEAL_V11_v1");
    lmEs.push_back("LM1"); lmEds.push_back("SUSY_LM1-sftsht_Summer08_IDEAL_V11_redigi_v1");
    lmEs.push_back("LM2"); lmEds.push_back("SUSY_LM2-sftsht_Summer08_IDEAL_V11_redigi_v1");
    lmEs.push_back("LM3"); lmEds.push_back("SUSY_LM3-sftsht_Summer08_IDEAL_V11_redigi_v1");
    lmEs.push_back("LM4"); lmEds.push_back("SUSY_LM4-sftsht_Summer08_IDEAL_V11_redigi_v1");
    lmEs.push_back("LM5"); lmEds.push_back("SUSY_LM5-sftsht_Summer08_IDEAL_V11_redigi_v1");
    for(unsigned int iLm=0; iLm < lmEs.size(); ++iLm){
      cout << "Processing  ... "<<lmEs[iLm].Data()<<endl;
      TChain* chLM = new TChain("Events");
      pickSkimIfExists(chLM, Form("data/%s/merged*.root",lmEds[iLm].Data()), "");
      looper->ScanChain(chLM, lmEs[iLm].Data(), kLM0, preLM0, oldjet, bitmask);
    }
  }
    
  //save all the histograms
    
  const char* outFile = 0;
  if (!runVgamma) outFile = Form("myHistComb_testfix_%d_%s.root", bitmask, looper->compactConfig.c_str());
  else outFile = Form("myHistComb_wExtras_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  

void doDYandTT_PY(unsigned int bitmask, bool skipFWLite = false){
  //see cuts written up above in the doAll()
 
  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)",skipFWLite));

  // Load and compile the looping code
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  
  // Flag for jet selection
  // true  = hyp_jet selection (15 GeV uncorrectedm eta<3)
  // false = 30 GeV corrected, eta<2.4
  bool oldjet=false;

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO
  
  //these two are taken from Ceballos's pdf. 
  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; 
  float kttotr    = 1.; 

  float kDYee     = 1.;
  float kDYmm     = 1.; 
  float kDYtautau = 1.; 

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;

  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = true;

  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;

  TChain* chtopdil = new TChain("Events");
  pickSkimIfExists(chtopdil, "data/TauolaTTbar-Pythia/merged*.root", "_skimSimple2020anydil");
  TChain* chtopotr = new TChain("Events");
  pickSkimIfExists(chtopotr, "data/TauolaTTbar-Pythia/merged*.root", "_skimSimple2020nodil");


  //Need to include the same mass range
  TChain* chDYtautau = new  TChain("Events");
  pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V9_v1/merged*.root_skimSimple2020_m50", "");
  TChain* chDYee = new  TChain("Events");
  pickSkimIfExists(chDYee, "data/Zee_M20_Summer08_IDEAL_V9_reco-v3/merged*.root_skimSimple2020_m50", "");
  TChain* chDYmm = new  TChain("Events");
  pickSkimIfExists(chDYmm, "data/Zmumu_M20_Summer08_IDEAL_V9_reco-v2/merged*.root_skimSimple2020_m50", "");

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    cout << "Done Processing ttbar dileptonic.. "<<endl;
    hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }

  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    looper->ScanChain(chDYee,"DYee", kDYee, preDYee, oldjet, bitmask);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, oldjet, bitmask);
    hist::color("DYmm", kCyan);
  }
  //save all the histograms
    
  const char* outFile = Form("myHist_DYandTT_PY_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
void doDYandTT_MG(unsigned int bitmask, bool skipFWLite = false){
  //see cuts written up above in the doAll()
 
  // Load various tools  
  gROOT->ProcessLine(Form(".x setup.C(%d)",skipFWLite));

  // Load and compile the looping code
  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  
  // Flag for jet selection
  // true  = hyp_jet selection (15 GeV uncorrectedm eta<3)
  // false = 30 GeV corrected, eta<2.4
  bool oldjet=false;

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO
  
  //these two are taken from Ceballos's pdf. 
  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; 
  float kttotr    = 1.; 

  float kDYee     = 1.;
  float kDYmm     = 1.; 
  float kDYtautau = 1.; 

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;

  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = true;

  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;

  TChain* chtopdil = new TChain("Events");
  pickSkimIfExists(chtopdil, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020anydil");
  TChain* chtopotr = new TChain("Events");
  pickSkimIfExists(chtopotr, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020nodil");

  TChain* chDYtautau = new  TChain("Events");
  pickSkimIfExists(chDYtautau, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020tautau");
  TChain* chDYee = new  TChain("Events");
  pickSkimIfExists(chDYee, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020ee");
  TChain* chDYmm = new  TChain("Events");
  pickSkimIfExists(chDYmm, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020mm");
  

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    looper->ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    cout << "Done Processing ttbar dileptonic.. "<<endl;
    hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    looper->ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }

  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper->ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    looper->ScanChain(chDYee,"DYee", kDYee, preDYee, oldjet, bitmask);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    looper->ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, oldjet, bitmask);
    hist::color("DYmm", kCyan);
  }
  //save all the histograms
    
  const char* outFile = Form("myHist_DYandTT_MG_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
