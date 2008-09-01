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

  // ntuple base directory
  const char* base = "/home/gutsche/data/dietcms2";

  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

  // Load various tools
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L CMS2.C+");

  //WW file
  if (runWW) {
    TChain *fWW = new TChain("Events");
    fWW->Add(Form("%s/WW_signal/*.root",base));
  }

  //WZ file
  if (runWZ) {
    TChain *fWZ = new TChain("Events");
    fWZ->Add(Form("%s/WZ_signal/*.root",base));
  }

  //ZZ file
  if (runZZ) {
    TChain *fZZ = new TChain("Events");
    fZZ->Add(Form("%s/ZZ_signal/*.root",base));
  }

  //Wjets file
  if (runWjets) {
    TChain *fWjets = new TChain("Events");
    fWjets->Add(Form("%s/cms2_muon_soup_postprocessed_split_Wjet/*.root",base));
    fWjets->Add(Form("%s/cms2_electron_soup_postprocessed_split_Wjet/*.root",base));
  }

  //DYee file
  if (runDYee) {
    TChain *fDYee = new TChain("Events");
    fDYee->Add(Form("%s/cms2_muon_soup_postprocessed_split_DY/*.root",base));
    fDYee->Add(Form("%s/cms2_electron_soup_postprocessed_split_DY/*.root",base));
  }

  //DYmm file
  if (runDYmm) {
    TChain *fDYmm = new TChain("Events");
    fDYmm->Add(Form("%s/cms2_muon_soup_postprocessed_split_DY/*.root",base));
    fDYmm->Add(Form("%s/cms2_electron_soup_postprocessed_split_DY/*.root",base));
  }

  //DYtt file
  if (runDYtt) {
    TChain *fDYtt = new TChain("Events");
    fDYtt->Add(Form("%s/cms2_muon_soup_postprocessed_split_DY/*.root",base));
    fDYtt->Add(Form("%s/cms2_electron_soup_postprocessed_split_DY/*.root",base));
  }

  //ttbar file
  if (runttbar) {
    TChain *fttbar = new TChain("Events");
    fttbar->Add(Form("%s/cms2_muon_soup_postprocessed_split_ttbar/*.root",base));
    fttbar->Add(Form("%s/cms2_electron_soup_postprocessed_split_ttbar/*.root",base));
  }

  if (runtW) {
    TChain *ftW = new TChain("Events");
    ftW->Add(Form("%s/tW_signal/*.root",base));
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
