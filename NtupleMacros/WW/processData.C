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
bool runWZ    = false;
bool runZZ    = false;
bool runWjets = true;
bool runDYee  = true;
bool runDYmm  = true;
bool runDYtt  = true;
bool runttbar = true;
bool runtW    = true;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
//gROOT->ProcessLine(".L fkwLoopingFunctionFast.C+");
//gROOT->ProcessLine(".L ClaudioLoopingFunctionFast.C+");

// read dataset prefix
 ifstream f;
 f.open("dataset.txt");
 if ( ! f.is_open() ) {
   cout << "Dataset location is not set. Please create dataset.txt file that has dataset prefix" <<endl;
   return;
 }
 string dataset;
 f >> dataset;
 f.close();

//WW file
TChain *fWW = new TChain("Events");
if (runWW) {
  fWW->Add((dataset+"/ww.root").c_str());
}

//WZ file
TChain *fWZ = new TChain("Events");
if (runWZ) {
  fWZ->Add("/data/tmp/dietcms2/cms2_WZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
}

//ZZ file
TChain *fZZ = new TChain("Events");
if (runZZ) {
  fZZ->Add("/data/tmp/dietcms2/cms2_ZZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
}

//Wjets file
TChain *fWjets = new TChain("Events");
if (runWjets) {
  fWjets->Add((dataset+"/wjets.root").c_str());
}

//DYee file
TChain *fDYee = new TChain("Events");
if (runDYee) {
  fDYee->Add((dataset+"/dy.root").c_str());
}

//DYmm file
TChain *fDYmm = new TChain("Events");
if (runDYmm) {
  fDYmm->Add((dataset+"/dy.root").c_str());
}

//DYtt file
TChain *fDYtt = new TChain("Events");
if (runDYtt) {
  fDYtt->Add((dataset+"/dy.root").c_str());
}

//ttbar file
TChain *fttbar = new TChain("Events");
if (runttbar) {
  fttbar->Add((dataset+"/ttbar.root").c_str());
}

//tW file
TChain *ftW = new TChain("Events");
if (runtW) {
  ftW->Add((dataset+"/tw.root").c_str());
}

// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

// Process files one at a time, and color them as needed
if (runWW) {
  cout << "Processing WW.."<< endl;
  ScanChain(fWW, WW);
  hist::color("ww", kRed);
}

if (runWZ) {
  cout << "Processing WZ.."<< endl;
  ScanChain(fWZ, WZ);
  hist::color("wz", kBlue);
}

if (runZZ) {
  cout << "Processing ZZ.."<< endl;
  ScanChain(fZZ, ZZ);
  hist::color("zz", kGreen);
}

if (runWjets) {
  cout << "Processing Wjets.."<<endl;
  ScanChain(fWjets, Wjets);
  hist::color("wjets", 40);
}

if (runDYee) {
  cout << "Processing DYee.."<<endl;
  ScanChain(fDYee, DYee);
  hist::color("dyee", kMagenta);
}

if (runDYmm) {
  cout << "Processing DYmm.."<<endl;
  ScanChain(fDYmm, DYmm);
  hist::color("dymm", kCyan);
}

if (runDYtt) {
  cout << "Processing DYtt.."<<endl;
  ScanChain(fDYtt, DYtt);
  hist::color("dytt", kBlack);
}

if (runttbar) {
  cout << "Processing ttbar.."<<endl;
  ScanChain(fttbar, ttbar);
  hist::color("ttbar", kYellow);
}

if (runtW) {
  cout << "Processing tW.."<<endl;
  ScanChain(ftW, tW);
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

}
