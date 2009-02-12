//==============================================================
//
// This runs over the skimmed files (see AAREADME.txt)
//
// To run on unskimmed files, change the file names.
//
//
//
//==============================================================
void doAll(unsigned int bitmask){
  
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
  float kDYee     = 1.;  //1230 pb,  970360 events processed
  float kDYmm     = 1.;  //1230 pb,  970360 events processed
  float kDYtautau = 1.;  //1230 pb,  970360 events processed
  float kppMuX    = 1.; //xsec/nevents
  float kEM       = 1.;
  float ktW       = 1.; //the evtScale is all negative for some reason
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
  bool runttotr    = true;
  bool runWW       = true;
  bool runWZ       = true;
  bool runZZ       = true;
  bool runWjets    = true;
  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;
  bool runppMuX    = false;
  bool runEM       = false;
  bool runtW       = true;
  bool runWQQ      = true;

  // Load and compile something to allow proper treatment of vectors
  // Not clear that it is needed
  gROOT->LoadMacro("loader.C+");

  // Load various tools  
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L ttDilCounts_looper.C++");

  TChain* chtopdil = new TChain("Events");
//   chtopdil->Add("/net/stau/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08/ntuple_diet_*.root");
//  chtopdil->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/TTJets-madgraph/merged*.root");
  chtopdil->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v2/merged*.root");

  TChain* chtopotr = new TChain("Events");
//   chtopotr->Add("/net/stau/cdf26/dietcms2/V01-00-04/PYTHIA_TauolaTTbar_Summer08/ntuple_diet_*.root");
//   chtopotr->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/TTJets-madgraph/merged*.root");
  chtopotr->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v2/merged*.root");

  TChain* chww = new TChain("Events");
//   chww->Add("data/signal/skim/ntuplemaker_WW_incl_2020.root");
//   chww->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/WW_2l-Pythia/merged*.root");
  chww->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/merged*.root");

  TChain* chWZ = new TChain("Events");
//   chWZ->Add("data/signal/skim/ntuplemaker_WZ_incl_2020.root");
//   chWZ->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/WZ_incl-Pythia/merged*.root"); // can try WZ_3l-Pythia
  chWZ->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/WZ_3l_Summer08_IDEAL_V9_v2/merged*.root"); // can try WZ_3l-Pythia

  TChain* chZZ = new TChain("Events");
//   chZZ->Add("data/signal/skim/ntuplemaker_ZZ_incl_2020.root");
//   chZZ->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/ZZ_2l2n-Pythia/merged*.root");
  chZZ->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/merged*.root");
  
  TChain* chWjets = new  TChain("Events");
  //  chWjets->Add("/net/stau/cdf26/dietcms2/V01-00-04/MadGraph_Wjets_Summer08/ntuple_diet_*.root");
//   chWjets->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/WJets-madgraph/merged*.root");
  chWjets->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged*.root");

  TChain* chDYtautau = new  TChain("Events");
  //  chDYtautau->Add("/net/stau/cdf26/dietcms2/V01-00-04/MadGraph_Zjets_Summer08/ntuple_diet_*.root");
//   chDYtautau->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/ZJets-madgraph/merged*.root");
  chDYtautau->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_reco-v2/merged*.root");
  
  TChain* chDYee = new  TChain("Events");
  //  chDYee->Add("/net/stau/cdf26/dietcms2/V01-00-04/MadGraph_Zjets_Summer08/ntuple_diet_*.root");
//   chDYee->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/ZJets-madgraph/merged*.root");
  chDYee->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_reco-v2/merged*.root");

  TChain* chDYmm = new  TChain("Events");
//  chDYmm->Add("/net/stau/cdf26/dietcms2/V01-00-04/MadGraph_Zjets_Summer08/ntuple_diet_*.root");
//   chDYmm->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/ZJets-madgraph/merged*.root");
  chDYmm->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_reco-v2/merged*.root");
  
  //ppMuX
  TChain* chppMuX = new  TChain("Events");
  if (runppMuX) {
    //    chppMuX->Add("/data3/slava77/cms/mc/Mupt15Inclusive_Summer08/V01-00-04/diet/ntuple_*.root");
//     chppMuX->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/InclusiveMuPt15/merged*.root"); 
    chppMuX->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/InclusiveMuPt15/merged*.root"); 
    //can try InclusiveMu5Pt50 .. figure out how to merge later
  }
  
  //ppEM
  TChain* chEM =  new  TChain("Events");
  if (runEM) {
//     chEM->Add("data/signal/skim/ntuplemaker_QCDJetsEnriched_eehyp_2020.root");
//     chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/QCD_EMenriched_Pt20to30/merged*.root");
//     chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/QCD_EMenriched_Pt30to80/merged*.root");
//     chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/QCD_EMenriched_Pt80to170/merged*.root");
    chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/QCD_EMenriched_Pt20to30/merged*.root");
    chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/QCD_EMenriched_Pt30to80/merged*.root");
    chEM->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/QCD_EMenriched_Pt80to170/merged*.root");
  }

  //tW
  TChain* chtW = new  TChain("Events");
  if (runtW) {
//     chtW->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/TopRex_tWinclusive/ntuple_merged_skim2020.root");
//     chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/SingleTop_sChannel-madgraph-LHE/merged*.root"); 
//     chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/SingleTop_tChannel-madgraph-LHE/merged*.root"); 
//     chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-01/SingleTop_tWChannel-madgraph-LHE/merged*.root"); 
    chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/SingleTop_sChannel-madgraph-LHE/merged*.root"); 
    chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/SingleTop_tChannel-madgraph-LHE/merged*.root"); 
    chtW->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/SingleTop_tWChannel-madgraph-LHE/merged*.root"); 
  }

  //WQQ
  TChain* chWQQ = new TChain("Events");
  if (runWQQ) {
//     chWQQ->Add("/net/bquark/data2/slava77/cms/cms1/V04-03-05/Alpgen_Wbb_0jets/ntuple_merged.root");
//    chWQQ->Add("");
    chWQQ->Add("/net/bquark/data3/slava77/cms/mc/cms2/V01-02-06/VQQ-madgraph_Fall08_IDEAL_V9_v1/merged*.root");
  }


  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  //do no cuts first:
//   int bitmaskR = atoi(getenv("FFFbitmask"));
//   unsigned int bitmask = (unsigned int)bitmaskR;

//  looperTTDisambiguation looper;
  ttDilCounts_looper looper;

  // Process files one at a time, and color them as needed
  if (runttdil) {
    cout << "Processing ttbar dileptonic.. "<<endl;
    looper.ScanChain(chtopdil,"ttdil", kttdil, prettdil, oldjet, bitmask);
    hist::color("ttdil", kYellow);
  }
  if (runttotr) {
    cout << "Processing ttbar no-dileptons.. "<<endl;
    looper.ScanChain(chtopotr,"ttotr", kttotr, prettotr, oldjet, bitmask);
    hist::color("ttotr", 30);
  }
  if (runWW) {
    cout << "Processing WW.."<<endl;
    looper.ScanChain(chww,"ww", kWW, preWW, oldjet, bitmask);
    hist::color("ww", kRed);
  }
  if (runWZ) {
    cout << "Processing WZ.."<<endl;
    looper.ScanChain(chWZ,"wz", kWZ, preWZ, oldjet, bitmask);
    hist::color("wz", kBlue);
  }
  if (runZZ) {
    cout << "Processing ZZ.."<<endl;
    looper.ScanChain(chZZ,"zz", kZZ, preZZ, oldjet, bitmask);
    hist::color("zz", kGreen);
  }
  if (runWjets) {
    cout << "Processing Wjets.."<<endl;
    looper.ScanChain(chWjets,"wjets", kWjets, preWjets, oldjet, bitmask);
    hist::color("wjets", 40);
  }
  if (runDYtautau) {
    cout << "Processing DY->tautau" << endl;
    looper.ScanChain(chDYtautau,"DYtautau", kDYtautau, preDYtautau, oldjet, bitmask);
    hist::color("DYtautau", kBlack);
  }
  if (runDYee) {
    cout << "Processing DY->ee" << endl;
    looper.ScanChain(chDYee,"DYee", kDYee, preDYee, oldjet, bitmask);
    hist::color("DYee", kMagenta);
  }
  if (runDYmm) {
    cout << "Processing DY->mm" << endl;
    looper.ScanChain(chDYmm,"DYmm", kDYmm, preDYmm, oldjet, bitmask);
    hist::color("DYmm", kCyan);
  }
  if (runppMuX) {
    cout << "Processing ppMuX"<<endl;
    looper.ScanChain(chppMuX,"ppMuX", kppMuX, preppMuX, oldjet, bitmask);
    hist::color("ppMuX", 51);
  }
  if (runEM) {
    cout << "Processing EM"<<endl;
    looper.ScanChain(chEM,"EM", kEM, preEM, oldjet, bitmask);
    hist::color("EM", 49);
  }
  if (runtW) {
    cout << "Processing tW"<<endl;
    looper.ScanChain(chtW,"tW", ktW, pretW, oldjet, bitmask);
    hist::color("tW", 63);
  }
    
  if (runWQQ) {
    cout << "Processing WQQ"<<endl;
    looper.ScanChain(chWQQ,"WQQ", kWQQ, preWQQ, oldjet, bitmask);
    hist::color("WQQ", 45);
  }

  //save all the histograms
    
  const char* outFile = Form("myHist_%d.root", bitmask);
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
