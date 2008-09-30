{

  // Output file
  char* outFile = "myHist_V00_04_00.root";

  //
  // ATTENTION
  //
  // V00-04-0x samples not complete
  // take WZ,ZZ,LM1 from V00-05-00
  //

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
  bool runLM1   = true;

  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

  // Load various tools
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L CMS2_V00_04_00.C+");

  //WW file
  if (runWW) {
    TChain *fWW = new TChain("Events");
    fWW->Add("/data/tmp/cms2-V00-04-0x/ww.root");
  }

  //WZ file
  if (runWZ) {
    TChain *fWZ = new TChain("Events");
    fWZ->Add("/data/tmp/cms2-V00-05-00/merge_WZ.root");
  }

  //ZZ file
  if (runZZ) {
    TChain *fZZ = new TChain("Events");
    fZZ->Add("/data/tmp/cms2-V00-05-00/merge_ZZ.root");
  }

  //Wjets file
  if (runWjets) {
    TChain *fWjets = new TChain("Events");
    fWjets->Add("/data/tmp/cms2-V00-04-0x/wjets.root");
  }

  //DYee file
  if (runDYee) {
    TChain *fDYee = new TChain("Events");
    fDYee->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //DYmm file
  if (runDYmm) {
    TChain *fDYmm = new TChain("Events");
    fDYmm->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //DYtt file
  if (runDYtt) {
    TChain *fDYtt = new TChain("Events");
    fDYtt->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //ttbar file
  if (runttbar) {
    TChain *fttbar = new TChain("Events");
    fttbar->Add("/data/tmp/cms2-V00-04-0x/ttbar.root");
  }

  //tW file
  if (runtW) {
    TChain *ftW = new TChain("Events");
    ftW->Add("/data/tmp/cms2-V00-04-0x/tw.root");
  }

  //LM1 file
  if (runLM1) {
    TChain *fLM1 = new TChain("Events");
    fLM1->Add("/data/tmp/cms2-V00-05-00/mergeLM1.root");
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

  if (runLM1) {
    cout << "Processing LM1.."<<endl;
    ScanChain(fLM1, "LM1",-1,1.0);
    hist::color("LM1", 98);
  }

  //save all the histograms
  saveHist(outFile);

}
