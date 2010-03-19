#include <iomanip>

void getYields(){

  //TFile *f=TFile::Open("victory_baseline_metgt50_sumjetptgt200_cand01_calo_tcmet_3x.root");
  TFile *f=TFile::Open("victory_baseline_metgt50_sumjetptgt200_calo_tcmet_3x.root");

  vector<char*> samples;
  samples.push_back("ttdil"); 
  samples.push_back("ttotr");
  samples.push_back("Zjets"); 
  samples.push_back("wjets"); 
  samples.push_back("tW");   
  samples.push_back("ww"); 
  samples.push_back("wz");
  samples.push_back("zz"); 
  samples.push_back("DYee");
  samples.push_back("LM0");
  samples.push_back("LM1");
  samples.push_back("LM2");
  samples.push_back("LM3");
  samples.push_back("LM4");
  samples.push_back("LM5");
  samples.push_back("LM6");
  samples.push_back("LM7");
  samples.push_back("LM8");
  samples.push_back("LM9");
  samples.push_back("LM10");
  samples.push_back("LM11");
  samples.push_back("LM12");
     
  makeTable(f,samples,100);
  makeTable(f,samples,175);

}

void makeTable(TFile* f, vector<char*> samples, float metcut){
  
  const unsigned int nsamples=samples.size();
  char* leptype[4]={"all","ee","mm","em"};
  TH1F *h;
  
  cout<<endl<<"Making table for tcmet > "<<metcut<<" GeV"<<endl;
  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  cout<<"|    Sample    |            all    |             ee    |             mm    |             em   |"<<endl;
  
  for(int isample = 0 ; isample < nsamples ; isample++){

    h = (TH1F*) f->Get(Form("%s_htcmet_allj_all",samples.at(isample),"all"));
    if(h==0) continue;
    delete h;
    
    cout<<"|"<<setw(10)<<samples[isample];

    for(int ilep = 0 ; ilep < 4 ; ilep++){
      h = (TH1F*) f->Get(Form("%s_htcmet_allj_%s",samples.at(isample),leptype[ilep]));
      
      int minbin = h->FindBin(metcut);
      int maxbin = h->GetNbinsX()+1;
      cout<<setprecision(3)<<"    |"<<setw(15)<<h->Integral(minbin,maxbin);
      
      delete h;
    }

    cout<<"   |"<<endl;
  }

  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  
}
