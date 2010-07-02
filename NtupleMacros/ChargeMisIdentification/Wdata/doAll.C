{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("Events"); 
  ch->Add("/home/users/jaehyeok/CMSSW_3_6_1_patch3_V03-04-25/src/ntuple.root");
  ScanChain(ch); 
}