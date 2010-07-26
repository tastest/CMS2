{

  gROOT->ProcessLine(".L ScanChain.C++");

  TChain *ch = new TChain("Events"); 
  ch->Add("/tas/cms2/Zmumu_Spring10-START3X_V26_S09-v1/V03-04-08-01/merged_ntuple*.root");
  
  myBabyMaker *looper = new myBabyMaker();
  looper->ScanChain(ch);
}
