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
char* outFile = "myHist_WW.root";

// Flags for files to run over
bool runWW    = true;
bool runWZ    = false;
bool runZZ    = false;
bool runWjets = false;
bool runDYee  = false;
bool runDYmm  = false;
bool runDYtt  = false;
bool runttbar = false;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
gROOT->ProcessLine(".L ClaudioLoopingFunction.C+");

//WW file
if (runWW) {
  TChain *fWW = new TChain("Events");
  fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1.root");
  fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_2.root");
  fWW->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/jmuelmen/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_3.root");
}

//WZ file
if (runWZ) {
  TChain *fWZ = new TChain("Events");
  fWZ->Add("data/wz/ntuple*.root");
}

//ZZ file
if (runZZ) {
  TChain *fZZ = new TChain("Events");
  fZZ->Add("data/zz/ntuple*.root");
}

//Wjets file
if (runWjets) {
  TChain *fWjets = new TChain("Events");
//   fWjets->Add("data/electron_soup/ntuple*.root");
  fWjets->Add("data/muon_soup/muon_soup.root");
}

//DYee file
if (runDYee) {
  TChain *fDYee = new TChain("Events");
//   fDYee->Add("data/electron_soup/ntuple*.root");
  fDYee->Add("data/muon_soup/muon_soup.root");
}

//DYmm file
if (runDYmm) {
  TChain *fDYmm = new TChain("Events");
//   fDYmm->Add("data/electron_soup/ntuple*.root");
  fDYmm->Add("data/muon_soup/muon_soup.root");
}

//DYtt file
if (runDYtt) {
  TChain *fDYtt = new TChain("Events");
//   fDYtt->Add("data/electron_soup/ntuple*.root");
  fDYtt->Add("data/muon_soup/muon_soup.root");
}

//ttbar file
if (runttbar) {
  TChain *fttbar = new TChain("Events");
//   fttbar->Add("data/electron_soup/ntuple*.root");
  fttbar->Add("data/muon_soup/muon_soup.root");
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
  TChain *t = TreePipe_glob("data/wz/ntuple*.root", "Events"); 
  ScanChain(t, WZ);
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

//save all the histograms
saveHist(outFile);
}
