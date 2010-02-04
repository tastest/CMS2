void skimDilPtXY(char* inD, char* outD, char* fname, float x=20, float y=10){
  gSystem->Load("libMiniFWLite");
  //  AutoLibraryLoader::enable();
  
  std::string fnameS = Form("%s/%s",inD,fname);
  TFile* tf = new TFile(fnameS.c_str());
  TTree* tr = (TTree*)tf->Get("Events");
  if (tr ==0){
    
  }

  std::string fonameS = Form("%s/%s_skimDil%d.%d",outD,fname,(int)x,(int)y);
  TFile* tfCp = new TFile(fonameS.c_str());
  if (tfCp->GetListOfKeys()->GetSize()==0){
    tfCp  = new TFile(fonameS.c_str(), "recreate");
  } else {
    std::cout<<fonameS<<" file exists: remove it first"<<std::endl;
    gSystem->Exit(55);
  }


  std::string selection = Form("Sum$((hyp_ll_p4.pt()>%f&&hyp_lt_p4.pt()>%f)||(hyp_lt_p4.pt()>%f&&hyp_ll_p4.pt()>%f))>0",x,y,x,y);
  int nEvents = tr->GetEntries(); 
  int skimEvents = tr->GetEntries(selection.c_str()); 
  std::cout<<selection.c_str()<<" will select "<<skimEvents<<" out of "<<nEvents<<std::endl;   

  trCp = tr->CopyTree(selection.c_str(), "");

  trCp->Write();

  tfCp->Write();

  std::cout<<"Closing file "<<tfCp->GetName()<<" with "<<trCp->GetEntries()<<" entries in "<<trCp->GetName()<<std::endl;
  tfCp->Close();
}
