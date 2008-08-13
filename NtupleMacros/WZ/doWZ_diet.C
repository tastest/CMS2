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
bool runWW    = true;
bool runWZ    = true;
bool runZZ    = true;
bool runWjets = true;
bool runDYee  = true;
bool runDYmm  = true;
bool runDYtt  = true;
bool runttbar = true;
bool runtW    = true;

// Load various tools
gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

// Load various tools
gROOT->ProcessLine(".x setup.C");

// Load and compile the looping code
gROOT->ProcessLine(".L CMS2.C+");

//WW file
if (runWW) {
     TChain *fWW = new TChain("Events");
     fWW->Add("/data/tmp/dietcms2/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
     fWW->Add("/data/tmp/dietcms2/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_2*.root");
     fWW->Add("/data/tmp/dietcms2/cms2_WW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_3*.root");
}

//WZ file
if (runWZ) {
     TChain *fWZ = new TChain("Events");
     fWZ->Add("/data/tmp/dietcms2/cms2_WZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
}

//ZZ file
if (runZZ) {
     TChain *fZZ = new TChain("Events");
     fZZ->Add("/data/tmp/dietcms2/cms2_ZZ_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
}

//Wjets file
if (runWjets) {
     TChain *fWjets = new TChain("Events");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_1*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_2*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_3*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_4*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_5*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_6*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_7*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_1*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_2*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_3*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_4*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_5*.root");
     fWjets->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_Wjet_1ac3e5ee47ec18c52404771fc74ace5d/ntuple_Wjet_6*.root");
}

//DYee file
if (runDYee) {
     TChain *fDYee = new TChain("Events");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYee->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7*.root");
}

//DYmm file
if (runDYmm) {
  TChain *fDYmm = new TChain("Events");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYmm->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7*.root");
}

//DYtt file
if (runDYtt) {
  TChain *fDYtt = new TChain("Events");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_1*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_2*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_3*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_4*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_5*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_6*.root");
     fDYtt->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_DY_095a1b92de7be271c0725c75b81c98ea/ntuple_DY_7*.root");
}

//ttbar file
if (runttbar) {
     TChain *fttbar = new TChain("Events");
     fttbar->Add("/data/tmp/dietcms2/cms2_muon_soup_postprocessed_split_ttbar_a8faccabfc12f0755896310bdda19928/ntuple_ttbar_1*.root");
     fttbar->Add("/data/tmp/dietcms2/cms2_electron_soup_postprocessed_split_ttbar_a8faccabfc12f0755896310bdda19928/ntuple_ttbar_1*.root");
}

if (runtW) {
     TChain *ftW = new TChain("Events");
     ftW->Add("/data/tmp/dietcms2/cms2_tW_signal_postprocessed_a38953977b4f365d80e08a78eaaff932/ntuple_signal_1*.root");
}

// Define colors numbers:
gStyle->SetPalette(1);
enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

// Process files one at a time, and color them as needed
if (runttbar) {
  cout << "Processing ttbar.."<<endl;
  ScanChain(fttbar,"tt",-1,1.85);
  hist::color("tt", kYellow);
}

if (runWW) {
  cout << "Processing WW.."<< endl;
  ScanChain(fWW, "ww",-1,1.0);
  hist::color("ww", kRed);
}

if (runWZ) {
  cout << "Processing WZ.."<< endl;
  ScanChain(fWZ, "wz",-1,1.0);
  hist::color("wz", kBlue);
}

if (runZZ) {
  cout << "Processing ZZ.."<< endl;
  ScanChain(fZZ, "zz",-1,1.0);
  hist::color("zz", kGreen);
}

if (runWjets) {
  cout << "Processing Wjets.."<<endl;
  ScanChain(fWjets, "wjets",-1,1.0);
  hist::color("wjets", 40);
}

if (runDYtt) {
  cout << "Processing DYtt.."<<endl;
  ScanChain(fDYtt, "DYtautau",2,1.12);
  hist::color("DYtautau", kBlack);
}

if (runDYee) {
  cout << "Processing DYee.."<<endl;
  ScanChain(fDYee, "DYee",0,1.12);
  hist::color("DYee", kMagenta);
}

if (runDYmm) {
  cout << "Processing DYmm.."<<endl;
  ScanChain(fDYmm, "DYmm",1,1.12);
  hist::color("DYmm", kCyan);
}

if (runtW) {
  cout << "Processing tW.."<<endl;
  ScanChain(ftW, "tW",-1,1.0);
  hist::color("tW", 63);
}

//save all the histograms
saveHist(outFile);

}
