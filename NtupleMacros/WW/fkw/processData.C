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
const char* outFile = "processed_data_tag_boostedZstudy.root";

// Flags for files to run over
bool runWW    = false;
bool runWZ    = false;
bool runZZ    = false;
bool runWjets = false;
bool runDYee  = true;
bool runDYmm  = true;
bool runDYtt  = true;
bool runttbar = false;
bool runtW    = false;
bool runLM1   = false;
bool runLM4   = false;
bool runLM8   = false;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// read dataset prefix
 string dataset;
 if ( ! gSystem->Getenv("CMS2_NTUPLE_LOCATION") ){
   cout << "ERROR: Dataset location is not set. Please set CMS2_NTUPLE_LOCATION." <<endl;
   return;
 }
 dataset = gSystem->Getenv("CMS2_NTUPLE_LOCATION");
 
//LM1 file
TChain *fLM1 = new TChain("Events");
if (runLM1) {
  fLM1->Add((dataset+"/cms2-V01-02-06/SUSY_LM1-sftsht/merged_ntuple*.root").c_str());
}
//LM4 file
TChain *fLM4 = new TChain("Events");
if (runLM4) {
  fLM4->Add((dataset+"/cms2-V01-02-06/SUSY_LM4-sftsht/merged_ntuple*.root").c_str());
}
//LM8 file
TChain *fLM8 = new TChain("Events");
if (runLM8) {
  fLM8->Add((dataset+"/cms2-V01-02-06/SUSY_LM8-sftsht/merged_ntuple*.root").c_str());
}

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

// Process files one at a time, and color them as needed
if (runLM1) {
  cout << "Processing LM1.."<< endl;
  ScanChain(fLM1, LM1);
  hist::color("LM1", kGreen);
}
if (runLM4) {
  cout << "Processing LM4.."<< endl;
  ScanChain(fLM4, LM4);
  hist::color("LM4", 50);
}
if (runLM8) {
  cout << "Processing LM8.."<< endl;
  ScanChain(fLM8, LM8);
  hist::color("LM8", 60);
}

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
