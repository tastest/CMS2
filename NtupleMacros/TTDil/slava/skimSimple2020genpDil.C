void skimSimple2020genpDil(char* fname, int type){
  if (type != 11 && type != 13 && type != 15) return;
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();
  
  TFile* tf = new TFile(fname);
  TTree* tr = (TTree*)tf->Get("Events");

  TFile* tfCp;
  if (type==11) tfCp = new TFile(Form("%s_skimSimple2020ee",fname), "recreate");
  if (type==13) tfCp = new TFile(Form("%s_skimSimple2020mm",fname), "recreate");
  if (type==15) tfCp = new TFile(Form("%s_skimSimple2020tautau",fname), "recreate");

  TTree* trCp;
  std::string selection = Form("Sum$(hyp_ll_p4.pt()>20&&hyp_lt_p4.pt()>20)>0&&Sum$(abs(genps_id)==%d)==2",type);
  int nEvents = tr->GetEntries();
  int skimEvents = tr->GetEntries(selection.c_str());
  std::cout<<selection.c_str()<<" will select "<<skimEvents<<" out of "<<nEvents<<std::endl;  
  trCp = tr->CopyTree(selection.c_str(), "");

  trCp->Write();

  tfCp->Write();

  std::cout<<"Closing file "<<tfCp->GetName()<<" with "<<trCp->GetEntries()<<" entries in "<<trCp->GetName()<<std::endl;
  tfCp->Close();
}
