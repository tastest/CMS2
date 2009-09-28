//void doAll(unsigned int bitmask, bool skipFWLite = false){
doAll(){
 float kWW       = 1.;
  float kWZ       = 1.;
  float kZZ       = 1.;
  float kWjets    = 1.; 
  float kWcharm   = 1;
  float kDYee     = 1.;  
  float kDYmm     = 1.;  
  float kDYtautau = 1.; 
  float kppMuX    = 1.; 
  float kEM       = 1.;
  float ktW       = 1.; 
  float kVQQ      = 1;
  float kZJets     = 1.; 
  float kWJets     = 1.;
  float kTTBar     = 1.; 
  float kQCD       = 1.; 
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
  int preZJets     = 1;
  int preWJets     = 1;
  int preTTBar     = 1;
  int preQCD       = 1;

  bool runttdil    = true;
  bool runttotr    = true;
  bool runWW       = true;
  bool runWZ       = true;
  bool runZZ       = true;
  bool runWjets    = true;
  bool runWcharm   = true;
  bool runDYee     = true;
  bool runDYmm     = true;
  bool runDYtautau = true;
  bool runppMuX    = true;
  bool runEM       = true;
  bool runtW       = true;
  bool runVQQ      = true;
  bool runZJets    = true;
  bool runWJets    = true;
  bool runTTBar    = true;
  bool runQCD      = true;
  
using namespace std;
 gROOT->ProcessLine(".x setup.C(true)");
 gSystem->CompileMacro("Ana_looper.C", "++k", "libAna_looper");
 
 string dataset;
 if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
   cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
   return;
 }
 dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");

 
 TChain *chDYtautau = new TChain("Events");
 chDYtautau->Add((dataset+"/cms2-V01-03-01/Ztautau_M20_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root").c_str());
 TChain *chDYee = new TChain("Events");
 chDYee->Add((dataset+"/cms2-V01-03-01/Zee_M20_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root").c_str());
 TChain *chZJets = new TChain("Events");
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chZJets->Add((dataset+"/cms2-V01-03-01/Z_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chZJets->Add((dataset+"/cms2-V01-03-01/Z_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chZJets->Add((dataset+"/cms2-V01-03-01/Z_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chZJets->Add((dataset+"/cms2-V01-03-01/Z_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());


 TChain *chWJets = new TChain("Events");
 chWJets->Add((dataset+"/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 

 chWJets->Add((dataset+"/cms2-V01-03-01/W_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chWJets->Add((dataset+"/cms2-V01-03-01/W_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chWJets->Add((dataset+"/cms2-V01-03-01/W_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());

 chWJets->Add((dataset+"/cms2-V01-03-01/W_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 
 chWJets->Add((dataset+"/cms2-V01-03-01/W_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());
 chWJets->Add((dataset+"/cms2-V01-03-01/W_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple*.root").c_str());




 TChain *chTTBar = new TChain("Events");
 chTTBar->Add((dataset+"/cms2-V01-03-01/TTJets-madgraph_Fall08_IDEAL_V11_redigi_v10-SingleLepton/merged_ntuple*.root").c_str());
 TChain *chQCD = new TChain("Events");
 chQCD->Add((dataset+"/cms2-V01-03-01/QCDpt30_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root").c_str());
 
 Ana_looper* looper = new Ana_looper();
//  if (runDYtautau) {
//     cout << "Processing DY->tautau" << endl;
//     looper->ScanChain(chDYtautau, -1, "DYtautau",kDYtautau,  preDYtautau);
//     hist::color("DYtautau", kBlack);
//   }
//   if (runDYee) {
//     cout << "Processing DY->ee" << endl;
//     looper->ScanChain(chDYee, -1, "DYee",  kDYee,  preDYee);
//     hist::color("DYee", kMagenta);
//   }

  if (runZJets) {
    cout << "Processing Z+0jets" << endl;
    looper->ScanChain(chZJets, -1, "ZJets",  kZJets,  preZJets);
    hist::color("ZJets", kBlue);
  }
  if (runWJets) {
    cout << "Processing W+0jets" << endl;
    looper->ScanChain(chWJets, -1, "WJets",  kWJets,  preWJets);
    hist::color("WJets", kRed);
  }

  if (runTTBar) {
    cout << "Processing TTBar" << endl;
    looper->ScanChain(chTTBar, -1, "TTBar",  kTTBar,  preTTBar);
    hist::color("TTBar", kRed);
  }

  if (runQCD) {
    cout << "Processing QCD" << endl;
    looper->ScanChain(chQCD, -1, "QCD",  kQCD,  preQCD);
    hist::color("QCD", kRed);
  }
 const char* outFile = "myHist.root";
 hist::saveHist(outFile);
 hist::deleteHistos();
}
