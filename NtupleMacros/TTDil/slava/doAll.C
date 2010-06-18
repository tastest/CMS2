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
#include "histscripts/histtools.h"
#include "ttDilCounts_looper.h"
#include "ttDilRefS10_looper.h"
#endif //__CINT__

#include "ProcDSS.h"

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
    gSystem->Exit(33);
    //    assert(0);
  }
  return;
  
}


void pickSkimIfExists(ProcDSS& pds, const std::string& name, 
		      const std::string& base, const std::string& skimExt,
                      float scale = 1., bool useWeigtFromBranch = true, bool checkDups = false){
  TChain* dummy = new TChain("Events");
  TChain* ch = new TChain("Events");
  ProcDSChain tmpDs(ch, name, scale, useWeigtFromBranch, checkDups );
  pds.add(tmpDs);
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
void pickSkimIfExists(ProcDSS& pds, const std::string& base, const std::string& skimExt,
                      float scale = 1., bool useWeigtFromBranch = true, bool checkDups = false){
  pickSkimIfExists(pds, pds.name, base, skimExt, scale, useWeigtFromBranch, checkDups);
}

void doAllCombined(unsigned long int bitmask, bool runTTbarOnly = false){
  using namespace std;
  //See AAA_ListOfBits.txt for the list of bits and their meaning

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO  

  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; //375pb, 127000 events processed
  float kttotr    = 1.; //375pb, 127000 events processed

  float kVV       = 1.;
  float kWjets    = 1.;
  float kDYeemm   = 1.;
  float kDYtautau = 1.;
  float kQCD      = 1.;
  float kt        = 1.;
  float ktotr     = 1.;
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
  int pretotr     = 1;
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
  bool runtotr     = true;
  bool runVgamma   = false;
  bool runLMs      = true;

  if (runTTbarOnly){
    runttdil=true;
    runttotr=true;
    runVV=false;
    runWjets=false;
    runDYeemm=false;
    runDYtautau=false;
    runQCD=false;
    runt=false;
    runtotr=false;
    runVgamma=false;
    runLMs=false;
  }

  ProcDSS chtopdil("ttdil",kttdil, prettdil);
  if (runttdil){
    pickSkimIfExists(chtopdil, "/data/tmp/slava77/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged*.root", "", 165.E+3*1./626610., false, false);
  }    
  
  ProcDSS chtopotr("ttotr",kttotr, prettotr);
  if (runttotr){
    pickSkimIfExists(chtopotr, "/data/tmp/slava77/cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged*.root", "", 165.E+3*1./626610., false, false);
  }    
  
  ProcDSS chVV("VV",kVV, preVV, bitmask);
  if (runVV){
    pickSkimIfExists(chVV, "/data/tmp/slava77/cms2/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "", 43.E+3*1./120280., false, false);
    pickSkimIfExists(chVV, "/data/tmp/slava77/cms2/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "", 18.E+3*1./114070., false, false); // can try WZ_3l-Pythia
    pickSkimIfExists(chVV, "/data/tmp/slava77/cms2/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "", 5.9E+3*1./145368., false, false);
  }    

  ProcDSS chWjets("wjets",kWjets, preWjets);
  if (runWjets){
    pickSkimIfExists(chWjets, "/data/tmp/slava77/cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 28.E+6*1./11217412., false, false);
  }
  
  ProcDSS chDYtautau("DYtautau", kDYtautau, preDYtautau);
  if(runDYtautau){
    pickSkimIfExists(chDYtautau, "/data/tmp/cms2/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/skimDilepton-20-10/skimZtautau.root", "", 1482.E+3*1./2193025., false, false);
  }

  ProcDSS chDYeemm("DYeemm", kDYeemm, preDYeemm);
  if (runDYeemm){
    pickSkimIfExists(chDYeemm, "/data/tmp/cms2/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/skimDilepton-20-10/skimZmumu.root", "", 1482.E+3*1./2078932., false, false);
    pickSkimIfExists(chDYeemm, "/data/tmp/cms2/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/skimDilepton-20-10/skimZee.root", "", 1482.E+3*1./2513855., false, false);
  }
  
  //ppMuX
  ProcDSS chQCD("QCD", kQCD, preQCD);
  if(runQCD){
    //ppMuX here  
    pickSkimIfExists(chQCD, "/data/tmp/cms2/InclusiveMu15_Summer09-MC_31X_V3_7TeV-v1_dilepfilt/V03-00-35/merged*.root", "", 1., true, false); 
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_BCtoE_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 1., true, false);
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_BCtoE_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 56.700382, false, false);
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_BCtoE_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 8.0148010, false, false);
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_EMEnriched_Pt20to30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 235.5E9/33545633.*0.0073, false, false);
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_EMEnriched_Pt30to80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 1., true, false);
    pickSkimIfExists(chQCD, "/data/tmp/slava77/cms2/QCD_EMEnriched_Pt80to170_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "_skimDil20.10", 24.373495, false, false);
  }

  //tW
  ProcDSS cht("t",kt, pret);
  if (runt){
    pickSkimIfExists(cht, "/data/tmp/slava77/cms2/SingleTop_tWChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root", "", 10.6E+3*1./473237., false, false); 
  }    

  ProcDSS chtotr("totr", ktotr, pretotr);
  if (runttotr){
    pickSkimIfExists(chtotr, "/data/tmp/slava77/cms2/SingleTop_sChannel-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged*.root", "", 4.6E+3*0.32442*1./409787., false, false); 
    pickSkimIfExists(chtotr, "/data/tmp/slava77/cms2/SingleTop_tChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged*.root", "", 63.E+3*0.32442*1./515167., false, false); 
  }

  //Vgamma
  ProcDSS chVgamma("Vgamma", kVgamma, preVgamma);
  if(runVgamma){
  }

  //LMs are run and loaded at the same time
  std::vector<TString> lmEs;
  std::vector<TString> lmEds;
  if (runLMs){
    lmEs.push_back("LM0"); lmEds.push_back("LM0_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM1"); lmEds.push_back("LM1_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM2"); lmEds.push_back("LM2_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM2m"); lmEds.push_back("LM2mhfeq360_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM3"); lmEds.push_back("LM3_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM4"); lmEds.push_back("LM4_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM5"); lmEds.push_back("LM5_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM6"); lmEds.push_back("LM6_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM7"); lmEds.push_back("LM7_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM8"); lmEds.push_back("LM8_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM9"); lmEds.push_back("LM9_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM9p"); lmEds.push_back("LM9p_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM9t"); lmEds.push_back("LM9t175_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM10"); lmEds.push_back("LM10_Summer09-MC_31X_V3_7TeV-v1");
    lmEs.push_back("LM12"); lmEds.push_back("LM12_Summer09-MC_31X_V3_7TeV-v1");
    for(unsigned int iLm=0; iLm < lmEs.size(); ++iLm){
      cout << "Loading to test presence of files  ... "<<lmEs[iLm].Data()<<endl;
      ProcDSS chLM(lmEs[iLm].Data(),kLM0, preLM0);
      pickSkimIfExists(chLM, Form("/data/tmp/slava77/cms2/%s/V03-00-35/merged*.root",lmEds[iLm].Data()), "", 1., true, false);
    }
  }

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (chtopdil.size()>0) {
    cout << "Processing ... "<<chtopdil.name<<endl;
    looper->ScanChain(chtopdil, bitmask);
    hist::color(chtopdil.name.c_str(), kYellow);
  }
  if (chtopotr.size()>0) {
    cout << "Processing  "<<chtopotr.name<<endl;
    looper->ScanChain(chtopotr, bitmask);
    hist::color(chtopotr.name.c_str(), 30);
  }
  if (chVV.size()>0) {
    cout << "Processing .."<<chVV.name<<endl;
    looper->ScanChain(chVV, bitmask);
    hist::color(chVV.name.c_str(), kRed);
  }
  if (chWjets.size()>0) {
    cout << "Processing "<<chWjets.name<<endl;
    looper->ScanChain(chWjets, bitmask);
    hist::color(chWjets.name.c_str(), 40);
  }
  if (chDYtautau.size()>0) {
    cout << "Processing "<<chDYtautau.name << endl;
    looper->ScanChain(chDYtautau, bitmask);
    hist::color(chDYtautau.name.c_str(), kBlack);
  }
  if (chDYeemm.size()>0) {
    cout << "Processing "<<chDYeemm.name << endl;
    looper->ScanChain(chDYeemm, bitmask);
    hist::color(chDYeemm.name.c_str(), kMagenta);
  }
  if (chQCD.size()>0) {
    cout << "Processing "<<chQCD.name<<endl;
    looper->ScanChain(chQCD, bitmask);
    hist::color(chQCD.name.c_str(), 51);
  }
  if (cht.size()>0) {
    cout << "Processing "<<cht.name<<endl;
    looper->ScanChain(cht, bitmask);
    hist::color(std::string(cht.name+"_").c_str(), 63);
  }
  if (chtotr.size()>0) {
    cout << "Processing "<<chtotr.name<<endl;
    looper->ScanChain(chtotr, bitmask);
    hist::color(std::string(chtotr.name+"_").c_str(), 67);
  }

  if (chVgamma.size()>0){
    cout << "Processing ... "<<chVgamma.name<<endl;
    looper->ScanChain(chVgamma, bitmask);
  }
    
  if (runLMs){
    for(unsigned int iLm=0; iLm < lmEs.size(); ++iLm){
      ProcDSS chLM(lmEs[iLm].Data(), kLM0, preLM0);
      cout << "Processing  ... "<<chLM.name<<endl;
      pickSkimIfExists(chLM, Form("/data/tmp/slava77/cms2/%s/V03-00-35/merged*.root",lmEds[iLm].Data()), "", 1., true, false);
      looper->ScanChain(chLM, bitmask);
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
  

void doAllRefS10(unsigned long int bitmask, bool runTTbarOnly = false){
  using namespace std;
  //See AAA_ListOfBits.txt for the list of bits and their meaning

  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to NLO  

  //It looks like the top x-section is for mtop = 175 GeV
  float kttdil    = 1.; //375pb, 127000 events processed
  float kttotr    = 1.; //375pb, 127000 events processed

  // Prescales
  int prettdil    = 1;
  int prettotr    = 1;

  // Flags for files to run over
  bool runttdil    = true;
  bool runttotr    = true;

  if (runTTbarOnly){
    runttdil=true;
    runttotr=true;
  }

  ProcDSS chtopdil("ttdil",kttdil, prettdil);
  if (runttdil){
    //    pickSkimIfExists(chtopdil, "/tas05/disk00/slava77/Reference*ntupleRefS10.root", "", 100., false, false);
    pickSkimIfExists(chtopdil, "/tas05/disk00/slava77/Reference*-V03-03-13-01.root", "", 100., false, false);
  }    
  
  ProcDSS chtopotr("ttotr",kttotr, prettotr);
  if (runttotr){
    //    pickSkimIfExists(chtopotr, "/tas05/disk00/slava77/Reference*ntupleRefS10.root", "", 100., false, false);
    pickSkimIfExists(chtopotr, "/tas05/disk00/slava77/Reference*-V03-03-13-01.root", "", 100., false, false);
  }    
  
  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilRefS10_looper* looper = new ttDilRefS10_looper();
  
  // Process files one at a time, and color them as needed
  if (chtopdil.size()>0) {
    cout << "Processing ... "<<chtopdil.name<<endl;
    looper->ScanChain(chtopdil, bitmask);
    hist::color(chtopdil.name.c_str(), kYellow);
  }
  if (chtopotr.size()>0) {
    cout << "Processing  "<<chtopotr.name<<endl;
    looper->ScanChain(chtopotr, bitmask);
    hist::color(chtopotr.name.c_str(), 30);
  }

  //save all the histograms    
  const char* outFile = 0;
  outFile = Form("histComb_refS10_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;
  gSystem->Exit(0);

}
  

void doDYandTT_PY(unsigned int bitmask){
  using namespace std;
  //see cuts written up above in the doAll()
 
  // Load various tools  
  //  gROOT->ProcessLine(Form(".x setup.C(%d)",skipFWLite));

  // Load and compile the looping code
  //  gSystem->CompileMacro("ttDilCounts_looper.C", "++k", "libttDilCounts_looper");
  
  // K-factors
  //these have been k-factors NLO/LO before
  //now using them as sample normalizations to something other than in the evt_scale1fb
  
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

  ProcDSS chtopdil("ttdil", kttdil, prettdil);
  pickSkimIfExists(chtopdil, "data/TauolaTTbar-Pythia/merged*.root", "_skimSimple2020anydil");
  ProcDSS chtopotr("ttotr", kttotr, prettotr);
  pickSkimIfExists(chtopotr, "data/TauolaTTbar-Pythia/merged*.root", "_skimSimple2020nodil");


  //Need to include the same mass range
  ProcDSS chDYtautau("DYtautau", kDYtautau, preDYtautau);
  pickSkimIfExists(chDYtautau, "data/Ztautau_M20_Summer08_IDEAL_V9_v1/merged*.root_skimSimple2020_m50", "");
  ProcDSS chDYee("DYee", kDYee, preDYee);
  pickSkimIfExists(chDYee, "data/Zee_M20_Summer08_IDEAL_V9_reco-v3/merged*.root_skimSimple2020_m50", "");
  ProcDSS chDYmm("DYmm", kDYmm, preDYmm);
  pickSkimIfExists(chDYmm, "data/Zmumu_M20_Summer08_IDEAL_V9_reco-v2/merged*.root_skimSimple2020_m50", "");

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (chtopotr.size()>0) {
    cout << "Processing "<<chtopdil.name<<endl;
    looper->ScanChain(chtopdil, bitmask);
    cout << "Done Processing ttbar dileptonic.. "<<endl;
    hist::color(chtopdil.name.c_str(), kYellow);
  }
  if (chtopdil.size()>0) {
    cout << "Processing .. "<< chtopotr.name<<endl;
    looper->ScanChain(chtopotr, bitmask);
    hist::color(chtopotr.name.c_str(), 30);
  }

  if (chDYtautau.size()>0) {
    cout << "Processing "<<chDYtautau.name << endl;
    looper->ScanChain(chDYtautau,bitmask);
    hist::color(chDYtautau.name.c_str(), kBlack);
  }
  if (chDYee.size()) {
    cout << "Processing "<<chDYee.name << endl;
    looper->ScanChain(chDYee,bitmask);
    hist::color(chDYee.name.c_str(), kMagenta);
  }
  if (chDYmm.size()) {
    cout << "Processing "<<chDYmm.name << endl;
    looper->ScanChain(chDYmm,bitmask);
    hist::color(chDYmm.name.c_str(), kCyan);
  }
  //save all the histograms
    
  const char* outFile = Form("myHist_DYandTT_PY_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
void doDYandTT_MG(unsigned int bitmask){
  using namespace std;
  //see cuts written up above in the doAll()

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

  ProcDSS chtopdil("ttdil", kttdil, prettdil);
  pickSkimIfExists(chtopdil, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020anydil");
  ProcDSS chtopotr("ttotr", kttotr, prettotr);
  pickSkimIfExists(chtopotr, "data/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10/merged*.root", "_skimSimple2020nodil");

  ProcDSS chDYtautau("DYtautau", kDYtautau, preDYtautau);
  pickSkimIfExists(chDYtautau, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020tautau");
  ProcDSS chDYee("DYee", kDYee, preDYee);
  pickSkimIfExists(chDYee, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020ee");
  ProcDSS chDYmm("DYmm", kDYmm, preDYmm);
  pickSkimIfExists(chDYmm, "data/ZJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root", "_skimSimple2020mm");
  

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };
  


  ttDilCounts_looper* looper = new ttDilCounts_looper();
  
  // Process files one at a time, and color them as needed
  if (chtopdil.size()) {
    cout << "Processing  "<<chtopdil.name<<endl;
    looper->ScanChain(chtopdil, bitmask);
    hist::color(chtopdil.name.c_str(), kYellow);
  }
  if (chtopotr.size()) {
    cout << "Processing .. "<<chtopotr.name<<endl;
    looper->ScanChain(chtopotr, bitmask);
    hist::color(chtopotr.name.c_str(), 30);
  }

  if (chDYtautau.size()) {
    cout << "Processing "<<chDYtautau.name << endl;
    looper->ScanChain(chDYtautau,bitmask);
    hist::color(chDYtautau.name.c_str(), kBlack);
  }
  if (chDYee.size()) {
    cout << "Processing "<<chDYee.name << endl;
    looper->ScanChain(chDYee,bitmask);
    hist::color(chDYee.name.c_str(), kMagenta);
  }
  if (chDYmm.size()) {
    cout << "Processing "<<chDYmm.name << endl;
    looper->ScanChain(chDYmm,bitmask);
    hist::color(chDYmm.name.c_str(), kCyan);
  }
  //save all the histograms
    
  const char* outFile = Form("myHist_DYandTT_MG_%d_%s.root", bitmask, looper->compactConfig.c_str());
  hist::saveHist(outFile);
  hist::deleteHistos();

  cout << "Finished with bitmask "<<bitmask << endl;


  gSystem->Exit(0);

}
  
