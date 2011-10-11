void makeSkim(const char* ntupleI, const char* ntupleO, const char* isData, const char *skim_C, const char* libMiniFWLite_so){
  gSystem->Load(libMiniFWLite_so);
  TTree::SetMaxTreeSize(39000000000ULL);
  //e = new TChain("Events");    ////// may not need the next 3 lines
  //e->Add(nmI);
  //e->SetBranchStatus("EventAuxiliary", 0);
  //  if (e->GetEntries()!=0){   ////want empty files to come back for better accounting
  gROOT->Macro(Form("%s++(\"%s\",\"%s\",false,%s)",skim_C,ntupleI,ntupleO,isData));
  
}
