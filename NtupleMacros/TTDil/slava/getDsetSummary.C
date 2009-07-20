void getDsetSummary(const char* name){
  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();
  e = new TChain("Events");
  e->Add(Form("%s/merged_ntuple.root", name));
  e->Scan("evtscale1fb:evtxsecexcl:evtkfactor:evtnEvts:evtfilteff:evtnEvts*evtscale1fb*1e-3:evtnEvts:evtnEvts/(evtnEvts*evtscale1fb*1e-3)", "", "", 3);
  

}
