{

  // Flags for files to run over
  bool runWW    = false;
  bool runWZ    = false;
  bool runZZ    = false;
  bool runWjets = true;
  bool runDYee  = false;
  bool runDYmm  = false;
  bool runDYtt  = false;
  bool runttbar = false;
  bool runtW    = false;
  bool runLM1   = false;

  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

  // Load various tools
  gROOT->ProcessLine(".x setup.C");

  // Load and compile the looping code
  gROOT->ProcessLine(".L FakeTest.C+");

  // load fake rate template histograms
  TFile *fakeRateFile = new TFile("fakeRates.root"); 
  TH2F *theFakeRate = (TH2F*)fakeRateFile->Get("fakeRateTemplate_wo_leading_elt_fakeRatesFull");
  theFakeRate->SetDirectory(0);

  gDirectory->cd("Rint:") ;

  //WW file
  if (runWW) {
    TChain *fWW = new TChain("Events");
    fWW->Add("/data/tmp/cms2-V00-05-00/mergeWW.root");
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
    fWjets->Add("/uscms_data/d1/cms1/cms2-V00-04-0x/wjets.root");
  }

  //DYee file
  if (runDYee) {
    TChain *fDYee = new TChain("Events");
    fDYee->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
    fDYee->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
  }

  //DYmm file
  if (runDYmm) {
    TChain *fDYmm = new TChain("Events");
    fDYmm->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
    fDYmm->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
  }

  //DYtt file
  if (runDYtt) {
    TChain *fDYtt = new TChain("Events");
    fDYtt->Add("/data/tmp/cms2-V00-05-00/electronsDY.root");
    fDYtt->Add("/data/tmp/cms2-V00-05-00/muonsDY.root");
  }

  //ttbar file
  if (runttbar) {
    TChain *fttbar = new TChain("Events");
    fttbar->Add("/data/tmp/cms2-V00-05-00/electronsttbar.root");
    fttbar->Add("/data/tmp/cms2-V00-05-00/muonsttbar.root");
  }

  //tW file
  if (runtW) {
    TChain *ftW = new TChain("Events");
    ftW->Add("/data/tmp/cms2-V00-05-00/mergetW.root");
  }

  //LM1 file
  if (runLM1) {
    TChain *fLM1 = new TChain("Events");
    fLM1->Add("/data/tmp/cms2-V00-05-00/mergeLM1.root");
  }


  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

  //Process one file at a time
  //flag = 0:  normal processing
  //flag = 1:  lepton + fakeable object
  //flag = 2:  lepton + fakeable object weighted by fake probability
  for (int flag=0; flag<3; flag++) {
    //for (int flag=2; flag<3; flag++) {

    // do not prepare the histos 
    // where the FOratio is not applied
    // to the FOs:
    if(flag == 1) continue;

    // Process files one at a time, and color them as needed
    if (runttbar) {
      cout << "Processing ttbar.."<<endl;
      ScanChain(fttbar,"tt",-1,1.85, flag, theFakeRate);
      //hist::color("tt", kYellow);
    }

    if (runWW) {
      cout << "Processing WW.."<< endl;
      ScanChain(fWW, "ww",-1,1.0, flag, theFakeRate);
      //hist::color("ww", kRed);
    }

    if (runWZ) {
      cout << "Processing WZ.."<< endl;
      ScanChain(fWZ, "wz",-1,1.0, flag, theFakeRate);
      //hist::color("wz", kBlue);
    }

    if (runZZ) {
      cout << "Processing ZZ.."<< endl;
      ScanChain(fZZ, "zz",-1,1.0, flag, theFakeRate);
      //hist::color("zz", kGreen);
    }

    if (runWjets) {
      cout << "Processing Wjets.."<<endl;
      ScanChain(fWjets, "wjets",-1,1.0, flag, theFakeRate);
      //hist::color("wjets", 40);
    }

    if (runDYtt) {
      cout << "Processing DYtt.."<<endl;
      ScanChain(fDYtt, "DYtautau",2,1.12, flag, theFakeRate);
      //hist::color("DYtautau", kBlack);
    }

    if (runDYee) {
      cout << "Processing DYee.."<<endl;
      ScanChain(fDYee, "DYee",0,1.12, flag, theFakeRate);
      //hist::color("DYee", kMagenta);
    }

    if (runDYmm) {
      cout << "Processing DYmm.."<<endl;
      ScanChain(fDYmm, "DYmm",1,1.12, flag, theFakeRate);
      //hist::color("DYmm", kCyan);
    }

    if (runtW) {
      cout << "Processing tW.."<<endl;
      ScanChain(ftW, "tW",-1,1.0, flag, theFakeRate);
      //hist::color("tW", 63);
    }

    if (runLM1) {
      cout << "Processing LM1.."<<endl;
      ScanChain(fLM1, "LM1",-1,1.0, flag, theFakeRate);
      //hist::color("LM1", 98);
    }

    //    gDirectory->ls();

    // save all the histograms
    saveHist(Form("myWHist_elt_%d.root",flag));

    // to be safe delete stuff 
    deleteHistos();

  }

}
