  {

  // Output file
  //char* outFile = "myHist_V00_05_00.root";
  char* outFile = "lm1_stacks.root";
  //char* holder = "";

  //char tag[100];
  char* tag;

  //Flag for naming outputfile based on cuts
  bool outFile_name_cuts = false;

  // Flags for files to run over
  bool runWW    = false;
  bool runWZ    = false;
  bool runZZ    = false;
  bool runWjets = false;
  bool runDYee  = false;
  bool runDYmm  = false;
  bool runDYtt  = false;
  bool runttbar = false;
  bool runtW    = false;
  bool runLM1   = true;

  gStyle->SetOptStat(110011);
  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../../Tools/").c_str());

  // Load various tools
  gROOT->ProcessLine(".x setup.C");
  
  // Load and compile the looping code
  // The old (original, or original with my edits) is in CMS2_V00_05_00.C
  //gROOT->ProcessLine(".L CMS2_V00_05_00.C+");
  gROOT->ProcessLine(".L warren_looper.C+");
  //gROOT->ProcessLine(".L warren_hists.C+");
  gROOT->ProcessLine(".L warren_functions.C+");

  
  //WW file
  if (runWW) {
    TChain *fWW = new TChain("Events");
    //fWW->Add("/data/tmp/cms2-V00-05-00/mergeWW.root");
	fWW->Add("/data/tmp/cms2-V00-04-0x/ww.root");
  }

  //WZ file
  if (runWZ) {
    TChain *fWZ = new TChain("Events");
//     fWZ->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/WZ_incl_P3-V00-05-00/ntuple_diet_merge.root");
    fWZ->Add("/data/tmp/cms2-V00-05-00/merge_WZ.root");
  }

  //ZZ file
  if (runZZ) {
    TChain *fZZ = new TChain("Events");
//     fZZ->Add("dcap://dcap-2.t2.ucsd.edu:22139//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/ZZ_incl_P3-V00-05-00/ntuple_diet_merge.root");
    fZZ->Add("/data/tmp/cms2-V00-05-00/merge_ZZ.root");
  }

  //Wjets file
  if (runWjets) {
    TChain *fWjets = new TChain("Events");
    //fWjets->Add("/data/tmp/cms2-V00-05-00/electronsWjet.root");
    //fWjets->Add("/data/tmp/cms2-V00-05-00/muonsWjet.root");
	fWjets->Add("/data/tmp/cms2-V00-04-0x/wjets.root");
  }

  //DYee file
  if (runDYee) {
    TChain *fDYee = new TChain("Events");
    //fDYee->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
	//Why add muonsDY.root to fDYee (?)
    //fDYee->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
	fDYee->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //DYmm file
  if (runDYmm) {
    TChain *fDYmm = new TChain("Events");
    //fDYmm->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
    //fDYmm->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
	fDYmm->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //DYtt file
  if (runDYtt) {
    TChain *fDYtt = new TChain("Events");
    //fDYtt->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
    //fDYtt->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
	fDYtt->Add("/data/tmp/cms2-V00-04-0x/dy.root");
  }

  //ttbar file
  if (runttbar) {
    TChain *fttbar = new TChain("Events");
    //fttbar->Add("/data/tmp/cms2-V00-05-00/electronsttbar.root");
    //fttbar->Add("/data/tmp/cms2-V00-05-00/muonsttbar.root");
	fttbar->Add("/data/tmp/cms2-V00-04-0x/ttbar.root");
  }

  //tW file
  if (runtW) {
    TChain *ftW = new TChain("Events");
    //ftW->Add("/data/tmp/cms2-V00-05-00/mergetW.root");
    //ftW->Add("/data/tmp/cms2-V00-04-00/mergetW.root");
	ftW->Add("/data/tmp/cms2-V00-04-0x/tw.root");
  }

  //LM1 file
  //Adding multiple files because merge isn't around atm.
  if (runLM1) {
    TChain *fLM1 = new TChain("Events");
    //fLM1->Add("/data/tmp/cms2-V00-05-00/mergeLM1.root");
	//fLM1->Add("/data/tmp/spadhi/V00-05-00/LM1/ntuple_post_signal_*.root");
	//fLM1->Add("/data/tmp/spadhi/V00-05-00/LM1/ntuple_post_signal_1.root");
	//These are on uaf-6 (?)
	fLM1->Add("/data/tmp/wandrews/V00-05-002/LM1/ntuple_post_signal_*.root");

	tag = "05-002";
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
    //holder = ScanChain(fLM1, "LM1",-1,1.0);
	ScanChain(fLM1, "LM1",tag,-1,1.0);
    hist::color("LM1", 98);
  }

  if( outFile_name_cuts ) {
	strcpy(outFile, holder);
	strcat(outFile, ".root");
  }
  //outFile = "lm1_stacks.root";
  //save all the histograms
  //saveHist(outFile);

}
