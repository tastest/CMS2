  {
  char* tag;

  // Flags for files to run over
  bool runLM1 = false;
  bool runLM3 = true;
  bool runLM4 = true;
  bool runLM5 = true;
  bool runLM6 = true;
  bool runLM7 = true;
  bool runLM8 = true;
  bool runLM9 = true;
  bool runLM10= true;
  bool runLM2 = false;
  bool runLM11= false;

  gStyle->SetOptStat(1110111);
  // Load various tools
  gROOT->SetMacroPath((string(gROOT->GetMacroPath()) + ":" + "../Tools/").c_str());

  // Load various tools
  //gROOT->ProcessLine(".x setup.C");
  
  // Load and compile the looping code
  gROOT->ProcessLine(".L warren_looper.C+");
  gROOT->ProcessLine(".L warren_functions.C+");

  //LM1 file
  if (runLM1) {
    TChain *fLM1 = new TChain("Events");
    //fLM1->Add("/data/tmp/cms2-V00-05-00/mergeLM1.root");
	//fLM1->Add("/data/tmp/spadhi/V00-05-00/LM1/ntuple_post_signal_*.root");
	//fLM1->Add("/data/tmp/spadhi/V00-05-00/LM1/ntuple_post_signal_1.root");
	//These are on uaf-6 
	fLM1->Add("/data/tmp/wandrews/V00-05-002/LM1/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM2 IS NOT DONE NTUPLIZING YET
  //LM2 file
  if (runLM2) {
    TChain *fLM2 = new TChain("Events");
	//fLM2->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM2_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	fLM3->Add("/data/tmp/wandrews/V00-05-002/LM3/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM3 file
  if (runLM3) {
    TChain *fLM3 = new TChain("Events");
	//fLM3->Add("dcache:dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM3_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	//fLM3->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM3_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_1.root");
	fLM3->Add("/data/tmp/wandrews/V00-05-002/LM3/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM4 file
  if (runLM4) {
    TChain *fLM4 = new TChain("Events");
	fLM4->Add("/data/tmp/wandrews/V00-05-002/LM4/ntuple_post_signal_*.root");
	//fLM4->AddFile("dcap://dcap-2.t2.ucsd.edu:22142//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM4_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_2.root");
	tag = "V00-05-002";
  }

  //LM5 file
  if (runLM5) {
    TChain *fLM5 = new TChain("Events");
	fLM5->Add("/data/tmp/wandrews/V00-05-002/LM5/ntuple_post_signal_*.root");
	//fLM5->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM5_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM6 file
  if (runLM6) {
    TChain *fLM6 = new TChain("Events");
	fLM6->Add("/data/tmp/wandrews/V00-05-002/LM6/ntuple_post_signal_*.root");
	//fLM6->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data6/cms/phedex/store/user/spadhi/CMS2/LM6_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM7 file
  if (runLM7) {
    TChain *fLM7 = new TChain("Events");
	fLM7->Add("/data/tmp/wandrews/V00-05-002/LM7/ntuple_post_signal_*.root");
	//fLM7->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM7_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM8 file
  if (runLM8) {
    TChain *fLM8 = new TChain("Events");
	fLM8->Add("/data/tmp/wandrews/V00-05-002/LM8/ntuple_post_signal_*.root");
	//fLM8->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM8_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM9 file
  if (runLM9) {
    TChain *fLM9 = new TChain("Events");
	fLM9->Add("/data/tmp/wandrews/V00-05-002/LM9/ntuple_post_signal_*.root");
	//fLM9->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM9_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM10 file
  if (runLM10) {
    TChain *fLM10 = new TChain("Events");
	fLM10->Add("/data/tmp/wandrews/V00-05-002/LM10/ntuple_post_signal_*.root");
	//fLM10->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM10_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  //LM11 IS NOT DONE NTUPLIZING YET
  //LM11 file
  if (runLM11) {
    TChain *fLM11 = new TChain("Events");
	fLM11->Add("/data/tmp/wandrews/V00-05-002/LM11/ntuple_post_signal_*.root");
	//fLM11->Add("dcap://dcap-2.t2.ucsd.edu:22138//pnfs/t2.ucsd.edu/data4/cms/phedex/store/user/spadhi/CMS2/LM11_sftsdkpyt_P3-V00-05-002/ntuple_post_signal_*.root");
	tag = "V00-05-002";
  }

  
  // Define colors numbers:
  gStyle->SetPalette(1);
  enum EColor { kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan };

  // Process files one at a time
  if (runLM1) {
    cout << "\nProcessing LM1.." << endl;
	ScanChain(fLM1, "LM1",tag,-1,1.0);
  }

  if (runLM2) {
    cout << "\nProcessing LM2.." << endl;
	ScanChain(fLM2, "LM2",tag,-1,1.0);
  }

  if (runLM3) {
    cout << "\nProcessing LM3.." << endl;
	ScanChain(fLM3, "LM3",tag,-1,1.0);
  }

  if (runLM4) {
    cout << "\nProcessing LM4.." << endl;
	ScanChain(fLM4, "LM4",tag,-1,1.0);
  }

  if (runLM5) {
    cout << "\nProcessing LM5.." << endl;
	ScanChain(fLM5, "LM5",tag,-1,1.0);
  }

  if (runLM6) {
    cout << "\nProcessing LM6.." << endl;
	ScanChain(fLM6, "LM6",tag,-1,1.0);
  }

  if (runLM7) {
    cout << "\nProcessing LM7.." << endl;
	ScanChain(fLM7, "LM7",tag,-1,1.0);
  }

  if (runLM8) {
    cout << "\nProcessing LM8.." << endl;
	ScanChain(fLM8, "LM8",tag,-1,1.0);
  }

  if (runLM9) {
    cout << "\nProcessing LM9.." << endl;
	ScanChain(fLM9, "LM9",tag,-1,1.0);
  }

  if (runLM10) {
    cout << "\nProcessing LM10.." << endl;
	ScanChain(fLM10, "LM10",tag,-1,1.0);
  }

  if (runLM11) {
    cout << "\nProcessing LM11.." << endl;
	ScanChain(fLM11, "LM11",tag,-1,1.0);
  }

}
