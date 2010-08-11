{

  gROOT->ProcessLine(".L doit.C++");

  doit();

  //TChain *ch = new TChain("Events"); 
  //ch->Add("/tas/cms2/Mu_Run2010A-Jun9thReReco_v1_RECO/V03-05-01/singleLepPt5Skim/skimmed_ntuple_1.root");
  //ScanChain(ch); 

}
