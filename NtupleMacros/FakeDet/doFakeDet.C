{

  // Output file
  char* outFile = "myHist.root";

  // Flags for files to run over
  bool run_qcd_0_15 	= true;
  bool run_qcd_15_20 	= true;
  bool run_qcd_20_30 	= true;
  bool run_qcd_30_50 	= true;
  bool run_qcd_50_80 	= true;
  bool run_qcd_80_120 	= true;
  bool run_qcd_120_170 	= true;
  bool run_qcd_170_230 	= true;
  bool run_qcd_230_300 	= true;
  bool run_qcd_300_380 	= true;
  bool run_qcd_380_470 	= true;
  bool run_qcd_470_600 	= true;

  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

  // Load various tools
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L CMS2.C+");

  // 0_15
  if ( run_qcd_0_15 ) {
    TChain *fQCD_0_15 = new TChain("Events");
    fQCD_0_15->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_0_15.root");
  }

  // 15_20
  if ( run_qcd_15_20 ) {
    TChain *fQCD_15_20 = new TChain("Events");
    fQCD_15_20->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_15_20.root");
  }

  // 20_30
  if ( run_qcd_20_30 ) {
    TChain *fQCD_20_30 = new TChain("Events");
    fQCD_20_30->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_20_30.root");
  }

  // 30_50
  if ( run_qcd_30_50 ) {
    TChain *fQCD_30_50 = new TChain("Events");
    fQCD_30_50->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_30_50.root");
  }

  // 50_80
  if ( run_qcd_50_80 ) {
    TChain *fQCD_50_80 = new TChain("Events");
    fQCD_50_80->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_50_80.root");
  }

  // 80_120
  if ( run_qcd_80_120 ) {
    TChain *fQCD_80_120 = new TChain("Events");
    fQCD_80_120->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_80_120.root");
  }

  // 120_170
  if ( run_qcd_120_170 ) {
    TChain *fQCD_120_170 = new TChain("Events");
    fQCD_120_170->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_120_170.root");
  }

  // 170_230
  if ( run_qcd_170_230 ) {
    TChain *fQCD_170_230 = new TChain("Events");
    fQCD_170_230->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_170_230.root");
//     fQCD_170_230->Add("/tmp/gutsche/qcd_pt_170_230.root");
  }

  // 230_300
  if ( run_qcd_230_300 ) {
    TChain *fQCD_230_300 = new TChain("Events");
    fQCD_230_300->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_230_300.root");
  }

  // 300_380
  if ( run_qcd_300_380 ) {
    TChain *fQCD_300_380 = new TChain("Events");
    fQCD_300_380->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_300_380.root");
  }

  // 380_470
  if ( run_qcd_380_470 ) {
    TChain *fQCD_380_470 = new TChain("Events");
    fQCD_380_470->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_380_470.root");
  }

  // 470_600
  if ( run_qcd_470_600 ) {
    TChain *fQCD_470_600 = new TChain("Events");
    fQCD_470_600->Add("/uscms_data/d1/cms1/cms2-V00-05-01/qcd_pt_470_600.root");
  }

  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

  // 0_15
  if ( run_qcd_0_15 ) {
    cout << "Processing QCD 0_15.." << endl;
    ScanChain(fQCD_0_15,"qcd_0_15");
    hist::color("qcd_0_15",1);
  }

  // 15_20
  if ( run_qcd_15_20 ) {
    cout << "Processing QCD 15_20.." << endl;
    ScanChain(fQCD_15_20,"qcd_15_20");
    hist::color("qcd_15_20",2);
  }

  // 20_30
  if ( run_qcd_20_30 ) {
    cout << "Processing QCD 20_30.." << endl;
    ScanChain(fQCD_20_30,"qcd_20_30");
    hist::color("qcd_20_30",3);
  }

  // 30_50
  if ( run_qcd_30_50 ) {
    cout << "Processing QCD 30_50.." << endl;
    ScanChain(fQCD_30_50,"qcd_30_50");
    hist::color("qcd_30_50",4);
  }

  // 50_80
  if ( run_qcd_50_80 ) {
    cout << "Processing QCD 50_80.." << endl;
    ScanChain(fQCD_50_80,"qcd_50_80");
    hist::color("qcd_50_80",6);
  }

  // 80_120
  if ( run_qcd_80_120 ) {
    cout << "Processing QCD 80_120.." << endl;
    ScanChain(fQCD_80_120,"qcd_80_120");
    hist::color("qcd_80_120",7);
  }

  // 120_170
  if ( run_qcd_120_170 ) {
    cout << "Processing QCD 120_170.." << endl;
    ScanChain(fQCD_120_170,"qcd_120_170");
    hist::color("qcd_120_170",8);
  }

  // 170_230
  if ( run_qcd_170_230 ) {
    cout << "Processing QCD 170_230.." << endl;
    ScanChain(fQCD_170_230,"qcd_170_230");
    hist::color("qcd_170_230",9);
  }

  // 230_300
  if ( run_qcd_230_300 ) {
    cout << "Processing QCD 230_300.." << endl;
    ScanChain(fQCD_230_300,"qcd_230_300");
    hist::color("qcd_230_300",28);
  }

  // 300_380
  if ( run_qcd_300_380 ) {
    cout << "Processing QCD 300_380.." << endl;
    ScanChain(fQCD_300_380,"qcd_300_380");
    hist::color("qcd_300_380",46);
  }

  // 380_470
  if ( run_qcd_380_470 ) {
    cout << "Processing QCD 380_470.." << endl;
    ScanChain(fQCD_380_470,"qcd_380_470");
    hist::color("qcd_380_470",42);
  }

  // 470_600
  if ( run_qcd_470_600 ) {
    cout << "Processing QCD 470_600.." << endl;
    ScanChain(fQCD_470_600,"qcd_470_600");
    hist::color("qcd_470_600",31);
  }

  //save all the histograms
  saveHist(outFile);

}
