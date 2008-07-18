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
char* outFile = "myHist.root";

// Flags for files to run over
bool runWW = true;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
gROOT->ProcessLine(".L myLoopingFunction.C+");

//WW file
if (runWW) {
  TChain *fWW = new TChain("Events");
  fWW->Add("ntuple*.root");
}

// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

// Process files one at a time, and color them as needed
if (runWW) {
  cout << "Processing WW.."<<endl;
  ScanChain(fWW);
  hist::color("ww", kBlue);
}

//save all the histograms
saveHist(outFile);

}
