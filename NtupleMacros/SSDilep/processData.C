//==============================================================
//
// This runs over the skimmed files (see AAREADME.txt)
//
// To run on unskimmed files, change the file names, but be
// careful about Drell Yan.
//
// The DY skims have separated out the three generated final states,
// while the unskimmed file has them all together, so it is a bit 
// more complicated.  You should:
// (1) Uncomment the block following the "Full Drell Yan file" comment
// (2) Optionally comment out blocks where the skimmed DY files are opened
// (3) Replace these three statements
//        ScanTree(tDYtautau,"DYtautau", -1, 1.2);
//        ScanTree(tDYmm,"DYmm", -1, 1.2);
//        ScanTree(tDYee,"DYee", -1, 1.2);
//     by
//        ScanTree(tDY,"DYtautau", 2, 1.2);
//        ScanTree(tDY,"DYmm",     1, 1.2);
//        ScanTree(tDY,"DYee",     0, 1.2);
//     Note the change in the 3rd calling parameter!
//
//==============================================================
{
// Output file
  const char* outFile = "processed_data_tag.root";
  const bool identifyVVEvents = false; // careful with this option. You don't need it for not mixed samples
  const bool identifyDYEvents = false;

// Flags for files to run over
bool runWW    = false;
bool runWZ    = false;
bool runZZ    = false;
bool runWjets = false;
bool runDYee  = false;
bool runDYmm  = false;
bool runDYtt  = false;
bool runDY  = false;
bool runttbar = true;
bool runtW    = false;
bool runLM0x    = false;
bool runLM1x    = false;
bool runLM2x    = false;
bool runLM3x    = false;
bool runLM4x    = false;
bool runLM5x    = false;
bool runLM6x    = false;
bool runLM7x    = false;
bool runLM8x    = false;
bool runLM9x    = false;

bool runWjetBackground1 = false;
bool runWjetBackground2 = false;

// Load various tools
// // string old_path = gROOT->GetMacroPath();
// // gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

// // gROOT->ProcessLine(".x setup.C");
// // gROOT->SetMacroPath(old_path.c_str());
// // gROOT->ProcessLine(".L fitWjets.C+");

 gROOT->ProcessLine(Form(".x setup.C(%d)", 1));
 gSystem->CompileMacro("doAnalysis.C", "++k", "libsusyosdiltp");
// gSystem->CompileMacro("doFlipAnalysis.C", "++k", "libsusyosdiltp");
// gSystem->CompileMacro("doAnalysisttbar.C", "++k", "libsusyttdiltp");
// gSystem->CompileMacro("doFakeAnalysis.C", "++k", "libsusyosdiltp");
// gSystem->CompileMacro("doAnalysisCompare.C", "++k", "libsusyosdiltp");

// read dataset prefix
 string dataset;
 if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
   cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
   return;
 }
 dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");
 
//WW file
TChain *fWW = new TChain("Events");
if (runWW) {
  fWW->Add((dataset+"/WW_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//WZ file
TChain *fWZ = new TChain("Events");
if (runWZ) {
  fWZ->Add((dataset+"/WZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//ZZ file
TChain *fZZ = new TChain("Events");
if (runZZ) {
  fZZ->Add((dataset+"/ZZ_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//Wjets file
TChain *fWjets = new TChain("Events");
if (runWjets) {
  fWjets->Add((dataset+"/Wenu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
  fWjets->Add((dataset+"/Wmunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
  fWjets->Add((dataset+"/Wtaunu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

// DY file

TChain *fDY = new TChain("Events");
if (runDY) {
  fDY->Add((dataset+"/Zee_Summer09-MC_31X_V3_7TeV_TrackingParticles-v1/V03-00-35/merged_ntuple*.root").c_str());
  fDY->Add((dataset+"/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
  fDY->Add((dataset+"/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}


//DYee file
TChain *fDYee = new TChain("Events");
if (runDYee) {
  fDYee->Add((dataset+"/cms2-V01-03-01/Zee_M20_Summer08_IDEAL_V11_redigi_v1/merged_ntuple*.root").c_str());
}

//DYmm file
TChain *fDYmm = new TChain("Events");
if (runDYmm) {
  fDYmm->Add((dataset+"/Zmumu_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//DYtt file
TChain *fDYtt = new TChain("Events");
if (runDYtt) {
  fDYtt->Add((dataset+"/Ztautau_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//ttbar file
TChain *fttbar = new TChain("Events");
if (runttbar) {
  fttbar->Add((dataset+"/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root").c_str());
}

//tW file
TChain *ftW = new TChain("Events");
if (runtW) {
  ftW->Add((dataset+"/SingleTop_sChannel-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple.root").c_str());
  ftW->Add((dataset+"/SingleTop_tChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple.root").c_str());
  ftW->Add((dataset+"/SingleTop_tWChannel-madgraph_Summer09-MC_31X_V3_7TeV-v2/V03-00-35/merged_ntuple.root").c_str());
}

//LM0
TChain *flm0x = new TChain("Events");
if (runLM0x) {
   flm0x->Add((dataset+"/LM0_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM1
TChain *flm1x = new TChain("Events");
if (runLM1x) {
  flm1x->Add((dataset+"/LM1_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM2
TChain *flm2x = new TChain("Events");
if (runLM2x) {
  flm2x->Add((dataset+"/LM2_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM3
TChain *flm3x = new TChain("Events");
if (runLM3x) {
  flm3x->Add((dataset+"/LM3_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM4
TChain *flm4x = new TChain("Events");
if (runLM4x) {
  flm4x->Add((dataset+"/LM4_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM5
TChain *flm5x = new TChain("Events");
if (runLM5x) {
  flm5x->Add((dataset+"/LM5_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM6
TChain *flm6x = new TChain("Events");
if (runLM6x) {
  flm6x->Add((dataset+"/LM6_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM7
TChain *flm7x = new TChain("Events");
if (runLM7x) {
  flm7x->Add((dataset+"/LM7_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM8
TChain *flm8x = new TChain("Events");
if (runLM8x) {
  flm8x->Add((dataset+"/LM8_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}

//LM9
TChain *flm9x = new TChain("Events");
if (runLM9x) {
  flm9x->Add((dataset+"/LM9_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root").c_str());
}



// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

 RooDataSet *fullDataSet(0);

// Process files one at a time, and color them as needed
if (runWW) {
  cout << "Processing WW.."<< endl;
  RooDataSet* data = ScanChain(fWW, WW, identifyVVEvents);
  if( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("ww", kRed);
}

if (runWZ) {
  cout << "Processing WZ.."<< endl;
  RooDataSet* data = ScanChain(fWZ, WZ, identifyVVEvents);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("wz", kBlue);
}

if (runZZ) {
  cout << "Processing ZZ.."<< endl;
  RooDataSet* data = ScanChain(fZZ, ZZ, identifyVVEvents);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("zz", kGreen);
}

if (runWjets) {
  cout << "Processing Wjets.."<<endl;
  RooDataSet* data = ScanChain(fWjets, Wjets, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("wjets", 40);
}

if (runDY) {
  cout << "Processing DY.."<<endl;
  RooDataSet* data = ScanChain(fDY, DY, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("dy", kMagenta);
}

if (runDYee) {
  cout << "Processing DYee.."<<endl;
  RooDataSet* data = ScanChain(fDYee, DYee, identifyDYEvents);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("dyee", kMagenta);
}

if (runDYmm) {
  cout << "Processing DYmm.."<<endl;
  RooDataSet* data = ScanChain(fDYmm, DYmm, identifyDYEvents);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("dymm", kCyan);
}

if (runDYtt) {
  cout << "Processing DYtt.."<<endl;
  RooDataSet* data = ScanChain(fDYtt, DYtt, identifyDYEvents);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("dytt", kBlack);
}

if (runttbar) {
  cout << "Processing ttbar.."<<endl;
  RooDataSet* data = ScanChain(fttbar, ttbar, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("ttbar", kYellow);
}

if (runtW) {
  cout << "Processing tW.."<<endl;
  RooDataSet* data = ScanChain(ftW, tW, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("tw", 63);
}

if (runLM0x) {
  cout << "Processing LM0.."<<endl;
  RooDataSet* data = ScanChain(flm0x, LM0x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("lm0x", kBlack);
}

if (runLM1x) {
   cout << "Processing LM1.."<<endl;
  RooDataSet* data = ScanChain(flm1x, LM1x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm1x", kBlack);
}

if (runLM2x) {
   cout << "Processing LM2.."<<endl;
  RooDataSet* data = ScanChain(flm2x, LM2x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm2x", kBlack);
}

if (runLM3x) {
   cout << "Processing LM3.."<<endl;
  RooDataSet* data = ScanChain(flm3x, LM3x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm3x", kBlack);
}

if (runLM4x) {
   cout << "Processing LM4.."<<endl;
  RooDataSet* data = ScanChain(flm4x, LM4x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm4x", kBlack);
}

if (runLM5x) {
   cout << "Processing LM5.."<<endl;
  RooDataSet* data = ScanChain(flm5x, LM5x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm5x", kBlack);
}

if (runLM6x) {
   cout << "Processing LM6.."<<endl;
  RooDataSet* data = ScanChain(flm6x, LM6x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm6x", kBlack);
}

if (runLM7x) {
   cout << "Processing LM7.."<<endl;
  RooDataSet* data = ScanChain(flm7x, LM7x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm7x", kBlack);
}

if (runLM8x) {
   cout << "Processing LM8.."<<endl;
  RooDataSet* data = ScanChain(flm8x, LM8x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm8x", kBlack);
}

if (runLM9x) {
   cout << "Processing LM9.."<<endl;
  RooDataSet* data = ScanChain(flm9x, LM9x, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
   hist::color("lm9x", kBlack);
}

//save all the histograms
//saveHist(outFile);
 cout << "got up to here" << endl;
 TList* list = gDirectory->GetList() ;
 TIterator* iter = list->MakeIterator();
 
 TRegexp re("*",kTRUE) ;
 
 TFile outf(outFile,"RECREATE") ;
 while(obj=iter->Next()) {
   if (TString(obj->GetName()).Index(re)>=0) {
     obj->Write() ;
     cout << "." ;
     cout.flush() ;
   }
 }
 cout << endl ;
 outf.Close() ;
 
 delete iter ;

 if ( runWjetBackground1 ){ 
   TCanvas* c1 = new TCanvas("wjetsBackgroundEstimates_sidebandfit","",800,800);
   c1->Divide(2,2);
   c1->cd(1);
   fit_isolation(fullDataSet,0,2,"Wjets e-fake background (pdf2)");
   c1->cd(2);
   fit_isolation(fullDataSet,0,1,"Wjets e-fake background (pdf1)");
   c1->cd(3);
   fit_isolation(fullDataSet,1,2,"Wjets mu-fake background (pdf2)");
   c1->cd(4);
   fit_isolation(fullDataSet,1,1,"Wjets mu-fake background (pdf1)");
 }
 if ( runWjetBackground2 ){ 
   TFile* fcs = TFile::Open("fakeIsoControlSamples.root");
   if ( fcs ){
     TCanvas* c2 = new TCanvas("wjetsBackgroundEstimates_qcd_sideband","",600,900);
     c2->Divide(2,3);
     c2->cd(1);
     fit_isolation(fullDataSet,0,3,"Wjets e-fake background (QCD30)",h_electron_qcd30);
     c2->cd(2);
     fit_isolation(fullDataSet,1,3,"Wjets mu-fake background (QCD30)",h_muon_qcd30);
     c2->cd(3);
     fit_isolation(fullDataSet,0,3,"Wjets e-fake background (QCD80)",h_electron_qcd80);
     c2->cd(4);
     fit_isolation(fullDataSet,1,3,"Wjets mu-fake background (QCD80)",h_muon_qcd80);
     c2->cd(5);
     fit_isolation(fullDataSet,0,3,"Wjets e-fake background (QCD170)",h_electron_qcd170);
     c2->cd(6);
     fit_isolation(fullDataSet,1,3,"Wjets mu-fake background (QCD170)",h_muon_qcd170);
   }
 }
}
