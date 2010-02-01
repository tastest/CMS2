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
  using namespace std;
// Output file
  const char* outFile = "processed_data.root";
  const bool identifyVVEvents = false; // careful with this option. You don't need it for not mixed samples
  const bool identifyDYEvents = false;

// Flags for files to run over
  bool runWW    = true;
  bool runWZ    = false;
  bool runZZ    = false;
  bool runWjets = false;
  bool runDYee  = false;
  bool runDYmm  = false;
  bool runDYtt  = false;
  bool runttbar = false;
  bool runtW    = true;
  bool runQCD   = false; 

// Load various tools
 gROOT->ProcessLine(".x init.C");
 gSystem->Load("libCMS2NtupleMacrosCORE");
 gSystem->Load("libRooFit.so");
 gSystem->Load("libCMS2NtupleMacrosLooper");
 gROOT->LoadMacro("../Tools/getMyHistosNames.C");
 gROOT->LoadMacro("../Tools/histtools.C+");
 gROOT->LoadMacro("../Tools/browseStacks.C");

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
  fWW->Add((dataset+"/cms2-V01-02-06/WW_Summer08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
  // fWW->Add("/data/tmp/cms2-V01-02-06/VVJets/merged_vvjets.root");
}

//WZ file
TChain *fWZ = new TChain("Events");
if (runWZ) {
  fWZ->Add((dataset+"/cms2-V01-02-06/WZ_incl_Summer08_IDEAL_V9_v2/merged_ntuple*.root").c_str());
  // fWZ->Add("/data/tmp/cms2-V01-02-06/VVJets/merged_vvjets.root");
}

//ZZ file
TChain *fZZ = new TChain("Events");
if (runZZ) {
  fZZ->Add((dataset+"/cms2-V01-02-06/ZZ_Summer08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
  // fZZ->Add("/data/tmp/cms2-V01-02-06/VVJets/merged_vvjets.root");
}

//Wjets file
TChain *fWjets = new TChain("Events");
if (runWjets) {
  fWjets->Add((dataset+"/cms2-V01-02-06/WJets-madgraph_Fall08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
}

//DYee file
TChain *fDYee = new TChain("Events");
if (runDYee) {
  fDYee->Add((dataset+"/cms2-V01-02-06/Zee_M20_Summer08_IDEAL_V9_reco-v3/merged_ntuple*.root").c_str());
}

//DYmm file
TChain *fDYmm = new TChain("Events");
if (runDYmm) {
  fDYmm->Add((dataset+"/cms2-V01-02-06/Zmumu_M20_Summer08_IDEAL_V9_reco-v2/merged_ntuple*.root").c_str());
}

//DYtt file
TChain *fDYtt = new TChain("Events");
if (runDYtt) {
  fDYtt->Add((dataset+"/cms2-V01-02-06/Ztautau_M20_Summer08_IDEAL_V9_v1/merged_ntuple*.root").c_str());
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

//QCD file
TChain *fqcd = new TChain("Events");
if (runQCD) {
  fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMuPt15/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/InclusiveMu5Pt50/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt20to30/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt30to80/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_EMenriched_Pt80to170/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt20to30/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt30to80/merged_ntuple*.root").c_str());
  fqcd->Add((dataset+"/cms2-V01-02-06/QCD_BCtoE_Pt80to170/merged_ntuple*.root").c_str());
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

if (runQCD) {
  cout << "Processing QCD.."<<endl;
  RooDataSet* data = ScanChain(fWjets, qcd, false);
  if ( data ){
    if ( fullDataSet )
      fullDataSet->append(*data);
    else
      fullDataSet=data;
  }
  hist::color("qcd", 40);
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
 if (fullDataSet) {
   std::string description("Full N-1 dataset:");
   if (runWW)     description+=" WW";
   if (runWZ)     description+=" WZ";
   if (runZZ)     description+=" ZZ";
   if (runWjets)  description+=" Wjets";
   if (runDYee)   description+=" DYee";
   if (runDYmm)   description+=" DYmm";
   if (runDYtt)   description+=" DYtt";
   if (runttbar)  description+=" ttbar";
   if (runtW)     description+=" tW";
   if (runQCD)    description+=" QCD";
   fullDataSet->SetName("fulldataset");
   fullDataSet->SetTitle(description.c_str());
   fullDataSet->Write();
 }
 outf.Close() ;
 
 delete iter ;

}
