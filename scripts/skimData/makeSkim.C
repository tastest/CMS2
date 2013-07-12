void makeSkim(const char* nmI, const char* nmO,  const char* skim_C, const char* miniFWlib){
  gSystem->Load(miniFWlib);
  TTree::SetMaxTreeSize(39000000000ULL);
  e = new TChain("Events");
  e->Add(nmI);
  if(e->GetBranch("EventAuxiliary") != 0)e->SetBranchStatus("EventAuxiliary", 0);
  if (e->GetEntries()!=0){
    gROOT->Macro(Form("%s++(\"%s\",\"%s\")",skim_C,nmI,nmO));
  }
}
