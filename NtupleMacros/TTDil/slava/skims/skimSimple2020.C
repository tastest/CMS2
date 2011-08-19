void skimSimple2020(char* fname){
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  TFile* tf = new TFile(fname);
  TTree* tr = (TTree*)tf->Get("Events");

  TFile* tfCp = new TFile(Form("%s_skimSimple2020",fname), "recreate");

  std::string selection = "Sum$(hyp_ll_p4.pt()>20&&hyp_lt_p4.pt()>20)>0";
  int nEvents = tr->GetEntries(); 
  int skimEvents = tr->GetEntries(selection.c_str()); 
  std::cout<<selection.c_str()<<" will select "<<skimEvents<<" out of "<<nEvents<<std::endl;   

  trCp = tr->CopyTree(selection.c_str(), "");

  trCp->Write();

  tfCp->Write();

  std::cout<<"Closing file "<<tfCp->GetName()<<" with "<<trCp->GetEntries()<<" entries in "<<trCp->GetName()<<std::endl;
  tfCp->Close();
}
