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

// Flags for files to run over
bool runWW    = true;
bool runWZ    = true;
bool runZZ    = true;
bool runWjets = true;
bool runDYee  = true;
bool runDYmm  = true;
bool runDYtt  = true;
bool runttbar = true;
bool runtW    = true;

bool runWjetBackground1 = false;
bool runWjetBackground2 = false;

// Load various tools
 string old_path = gROOT->GetMacroPath();
 gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

 gROOT->ProcessLine(".x setup.C");
 gROOT->SetMacroPath(old_path.c_str());
 gROOT->ProcessLine(".L fitWjets.C+");

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
  fWW->Add((dataset+"/cms2-V01-02-06/WW_2l_Summer08_IDEAL_V9_v2/merged_ntuple.root").c_str());
}

//WZ file
TChain *fWZ = new TChain("Events");
if (runWZ) {
  fWZ->Add((dataset+"/cms2-V01-02-06/WZ_3l_Summer08_IDEAL_V9_v2/merged_ntuple.root").c_str());
}

//ZZ file
TChain *fZZ = new TChain("Events");
if (runZZ) {
  fZZ->Add((dataset+"/cms2-V01-02-06/ZZ_2l2n_Summer08_IDEAL_V9_v2/merged_ntuple.root").c_str());
}

//Wjets file
TChain *fWjets = new TChain("Events");
if (runWjets) {
  fWjets->Add((dataset+"/cms2-V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//DYee file
TChain *fDYee = new TChain("Events");
if (runDYee) {
  fDYee->Add((dataset+"/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//DYmm file
TChain *fDYmm = new TChain("Events");
if (runDYmm) {
  fDYmm->Add((dataset+"/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//DYtt file
TChain *fDYtt = new TChain("Events");
if (runDYtt) {
  fDYtt->Add((dataset+"/cms2-V01-02-06/ZJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//ttbar file
TChain *fttbar = new TChain("Events");
if (runttbar) {
  fttbar->Add((dataset+"/cms2-V01-02-06/TTJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//tW file
TChain *ftW = new TChain("Events");
if (runtW) {
  ftW->Add((dataset+"/cms2-V01-02-06/SingleTop_tWChannel-madgraph-LHE/merged_ntuple.root").c_str());
}

// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

 RooDataSet *fullDataSet(0);

// Process files one at a time, and color them as needed
if (runWW) {
  cout << "Processing WW.."<< endl;
  RooDataSet* data = ScanChain(fWW, WW);
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
  RooDataSet* data = ScanChain(fWZ, WZ);
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
  RooDataSet* data = ScanChain(fZZ, ZZ);
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
  RooDataSet* data = ScanChain(fWjets, Wjets);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("wjets", 40);
}

if (runDYee) {
  cout << "Processing DYee.."<<endl;
  RooDataSet* data = ScanChain(fDYee, DYee);
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
  RooDataSet* data = ScanChain(fDYmm, DYmm);
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
  RooDataSet* data = ScanChain(fDYtt, DYtt);
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
  RooDataSet* data = ScanChain(fttbar, ttbar);
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
  RooDataSet* data = ScanChain(ftW, tW);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("tw", 63);
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
