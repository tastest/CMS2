//==============================================================
//
// This runs over the skimmed files (see AAREADME.txt)
//
// To run on unskimmed files, change the file names.
//
//
//
//==============================================================
{
  
  //No cuts - Cut0
  //Bitmask = 00000

  //ID cuts only - Cut1:
  //Bitmask = 00001
  
  //ID+isolation - Cut2
  //Bitmask = 00011
  
  //ID+isolation+dileptonMass veto - Cut3
  //Bitmask = 00111
  
  //ID+isolation+dileptonMass veto+MET (all of it) - Cut4
  //Bitmask 01111
  
  //ID+isolation+dileptonMass veto+MET (all of it)+NJETS>=2 - Cut5
  //Bitmask 11111
  
  //all except ID and nJet Cut - Cut6 
  //Bitmask 01110
  
  //all except iso and NJetCut - Cut7
  //Bitmask 01101
  
  //all except massVeto and NJet Cut - Cut 8
  //Bitmask 01011
  
  //all except MET and NJet  - Cut 9
  //Bitmask 00111
 
  //ID+isolation+dileptonMass veto + dima MET (only dima met!) + NJets - cuts 10
  //BitMask - 11111

  // Flag for jet selection
  // true  = hyp_jet selection (15 GeV uncorrectedm eta<3)
  // false = 30 GeV corrected, eta<2.4
  bool oldjet=false;

  // K-factors
  float kttdil    = 1.85;
  float kttotr    = 1.85;
  float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjets    = 1.12;
  float kDYee     = 1.12;
  float kDYmm     = 1.12;
  float kDYtautau = 1.12;
  float kppMuX    = 1.;
  float kEM       = 1.;
  float ktW       = -1.; //the evtScale is all negative for some reason
  float kWQQ      = 1;

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;
  int preWW       = 1;
  int preWZ       = 1;
  int preZZ       = 1;
  int preWjets    = 1;
  int preDYee     = 1;
  int preDYmm     = 1;
  int preDYtautau = 1;
  int preppMuX    = 1;
  int preEM       = 1;
  int pretW       = 1;
  int preWQQ      = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = false;
  bool runWW       = false;
  bool runWZ       = false;
  bool runZZ       = false;
  bool runWjets    = false;
  bool runDYee     = false;
  bool runDYmm     = false;
  bool runDYtautau = false;
  bool runppMuX    = false;
  bool runEM       = false;
  bool runtW       = false;
  bool runWQQ      = false;

  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gROOT->LoadMacro("loader.C+");

  // Load various tools  
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  //gROOT->ProcessLine(".L myLoopingFunctionFlags.C+");
  gROOT->ProcessLine(".L CMS2.C++");

  TChain* chtopdil = new TChain("event");
  //chtopdil->Add("data/Chowder/10pb/skim/Chowder-BahMeh_ttbar_dil_2020.root");
  //chtopdil->Add("data/cms2_electron_soup_postprocessed_split_ttbar/ntuple_ttbar_1.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_1.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_2.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_3.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_4.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_6.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_7.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_8.root");
  chtopdil->Add("/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08//ntuple_diet_9.root");

  TChain* chtopotr = new TChain("event");
  chtopotr->Add("data/Chowder/10pb/skim/Chowder-BahMeh_ttbar_nodil_2020.root");

  TChain* chww = new TChain("event");
  chww->Add("data/signal/skim/ntuplemaker_WW_incl_2020.root");

  TChain* chWZ = new TChain("event");
  chWZ->Add("data/signal/skim/ntuplemaker_WZ_incl_2020.root");

  TChain* chZZ = new TChain("event");
  chZZ->Add("data/signal/skim/ntuplemaker_ZZ_incl_2020.root");
  
  TChain* chWjets = new  TChain("event");
  chWjets->Add("data/Chowder/10pb/skim/Chowder-BahMeh_wjets_2020.root");

  TChain* chDYtautau = new  TChain("event");
  chDYtautau->Add("data/Chowder/10pb/skim/Chowder-BahMeh_ztautau_2020.root");
  chDYtautau->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DrellYan_ll_40/skim/DrellYan_tautau_200_inf.root");
  chDYtautau->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DY_ll_10-40/skim/DrellYan_tautau_10_40.root");

  TChain* chDYee = new  TChain("event");
  chDYee->Add("data/Chowder/10pb/skim/Chowder-BahMeh_zee_2020.root");
  chDYee->Add("data/Chowder/10pb/skim/Chowder-BahMeh_zee_2020_1.root");
  chDYee->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DrellYan_ll_40/skim/DrellYan_ee_200_inf.root");
  chDYee->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DY_ll_10-40/skim/DrellYan_ee_10_40.root");

  TChain* chDYmm = new  TChain("event");
  chDYmm->Add("data/Chowder/10pb/skim/Chowder-BahMeh_zmumu_2020.root");
  chDYmm->Add("data/Chowder/10pb/skim/Chowder-BahMeh_zmumu_2020_1.root");
  chDYmm->Add("data/Chowder/10pb/skim/Chowder-BahMeh_zmumu_2020_2.root");
  chDYmm->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DrellYan_ll_40/skim/DrellYan_mumu_200_inf.root");
  chDYmm->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/DY_ll_10-40/skim/DrellYan_mumu_10_40.root");

  //ppMuX
  TChain* chppMuX = new  TChain("event");
  if (runppMuX) {
    chppMuX->Add("data/signal/skim/ntuplemaker_ppMuPt20_2020.root");
  }
  
  //ppEM
  TChain* chEM =  new  TChain("event");
  if (runEM) {
    chEM->Add("data/signal/skim/ntuplemaker_QCDJetsEnriched_eehyp_2020.root");
  }

  //tW
  TChain* chtW = new  TChain("event");
  if (runtW) {
    chtW->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/TopRex_tWinclusive/ntuple_merged_skim2020.root");
  }

  //WQQ
  TChain* chWQQ = new TChain("event");
  if (runWQQ) {
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wbb_0jets/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wbb_1jet/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wbb_2jets/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wbb_3jets/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wcc_0jets/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wcc_1jet/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wcc_2jets/ntuple_merged.root");
    chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wcc_3jets/ntuple_merged.root");
  }


  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  //do no cuts first:
  
  int bitmaskR = atoi(getenv("FFFbitmask"));
  unsigned int bitmask = (unsigned int)bitmaskR;
  
  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    CMS2 tt;
    tt.ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    //hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }
  if (runWW) {
    cout << "Processing WW.."<<endl;
    ScanChain(chww,"ww", kWW, preWW, oldjet, bitmask);
    hist::color("ww", kRed);
  }
  if (runWZ) {
    cout << "Processing WZ.."<<endl;
    ScanChain(chWZ,"wz", kWZ, preWZ, oldjet, bitmask);
    hist::color("wz", kBlue);
  }
  if (runZZ) {
    cout << "Processing ZZ.."<<endl;
    ScanChain(chZZ,"zz", kZZ, preZZ, oldjet, bitmask);
    hist::color("zz", kGreen);
  }
  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
    ScanChain(chWjets,"wjets", kWjets, preWjets, oldjet, bitmask);
    hist::color("wjets", 40);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    ScanChain(chDYee,"DYee", kDYee, preDYee, oldjet, bitmask);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, oldjet, bitmask);
    hist::color("DYmm", kCyan);
  }
  if (runppMuX) {
    cout << "Processing ppMuX"<<endl;
    ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, oldjet, bitmask);
    hist::color("ppMuX", 51);
  }
  if (runEM) {
    cout << "Processing EM"<<endl;
    ScanChain(chEM,"EM", kEM, preEM, oldjet, bitmask);
    hist::color("EM", 49);
  }
  if (runtW) {
    cout << "Processing tW"<<endl;
    ScanChain(chtW,"tW", ktW, pretW, oldjet, bitmask);
    hist::color("tW", 63);
  }
    
  if (runWQQ) {
    cout << "Processing WQQ"<<endl;
    ScanChain(chWQQ,"WQQ", kWQQ, preWQQ, oldjet, bitmask);
    hist::color("WQQ", 45);
  }

  //save all the histograms
    
  const char* outFile = Form("myHist_%d.root", bitmask);
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
