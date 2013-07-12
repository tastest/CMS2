void makeSkim(const char* nmI, const char* nmO, const char* expr=""){
  gSystem->Load("~/MiniFWlib/libMiniFWLite_CMSSW_5_3_2_patch4_V05-03-13.so");
  TTree::SetMaxTreeSize(39000000000ULL);
  e = new TChain("Events");
  e->Add(nmI);
  if(e->GetBranch("EventAuxiliary") != 0)e->SetBranchStatus("EventAuxiliary", 0);
  if (e->GetEntries()!=0){
    if (expr=="TPrime"){
      gROOT->LoadMacro("ntupleFilterTPrime.cc++");
      ntupleFilterTPrime(nmI, nmO);
    }
  }
  //  
  //e = new TChain("Events");
  //e->Add(nmI);
  //e->SetBranchStatus("EventAuxiliary", 0);
  //if (e->GetEntries()!=0){
  //  fCp = new TFile(nmO,"recreate");
 
  //  chCp = e->CopyTree(expr);
  //  chCp->Write();
  //  fCp = gFile; fCp->Write();
  //  fCp->Close();
  //}
}
