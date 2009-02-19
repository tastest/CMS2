void skimSimple2020genpNoDil(char* fname){
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  TFile* tf = new TFile(fname);
  TTree* tr = (TTree*)tf->Get("Events");

  TFile* tfCp;
  tfCp = new TFile(Form("%s_skimSimple2020nodil",fname), "recreate");

  TTree* trCp;
  std::string selection = "Sum$(hyp_ll_p4.pt()>20&&hyp_lt_p4.pt()>20)>0&&Sum$(abs(genps_id)==11||abs(genps_id)==13||abs(genps_id)==15)<2";
  int nEvents = tr->GetEntries();
  int skimEvents = tr->GetEntries(selection.c_str());
  std::cout<<selection.c_str()<<" will select "<<skimEvents<<" out of "<<nEvents<<std::endl;  
  trCp = tr->CopyTree(selection.c_str(), "");

  trCp->Write();

  tfCp->Write();

  std::cout<<"Closing file "<<tfCp->GetName()<<" with "<<trCp->GetEntries()<<" entries in "<<trCp->GetName()<<std::endl;
  tfCp->Close();
}
