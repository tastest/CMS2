{

  gROOT->LoadMacro("CMS2.C+");

  TChain *chain = new TChain("Events");
  chain->Add("cms2_notqaf.root");

  ScanChain(chain);

}
