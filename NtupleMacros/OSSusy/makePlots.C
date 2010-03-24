#include "TROOT.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <string>
#include "TLegendEntry.h"
//#include "histtools.h"

namespace hist {
  void add(const char* outHistName, const char* patORpfx0);
  TLegend* legend(TCanvas* canvas, Option_t* option = "lpf", Bool_t addColor = kFALSE, Int_t token = -1, 
		  Float_t xmin = 0.75, Float_t ymin = 0.75, Float_t xmax = 0.99, Float_t ymax = 0.99);
}

////void loadHist(const char* filename, const char* pfx = 0, const char* pat = "*", Bool_t doAdd = kFALSE);
void loadHist(const char* filename, const char* directory = 0, const char* pfx = 0, 
	      const char* pat = "*", Bool_t doAdd = kFALSE);

void makeTable(char* filename, vector<char*> samples, float metcut,
               TH1F* hall, TH1F* hee, TH1F* hmm, TH1F* hem);

void makePlots(char* filename) {

  gROOT->SetStyle("Plain");
  gROOT->ProcessLine("gStyle->SetOptStat(0)");

  //Input samples--------------------------------------------------
  vector<char*> samples;
  samples.push_back("ttdil"); 
  samples.push_back("ttotr"); 
  samples.push_back("Zjets"); 
  samples.push_back("wjets"); 
  samples.push_back("tW"); 
  samples.push_back("ww"); 
  samples.push_back("wz"); 
  samples.push_back("zz"); 
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
   
  //Load and add histos--------------------------------------------
  loadHist(filename, 0, 0, "*dilPt_allj*", kTRUE);
  loadHist(filename, 0, 0, "*dilPtSmeared_allj*", kTRUE);
  loadHist(filename, 0, 0, "*tcmet_allj*", kTRUE);
  //loadHist(filename.c_str(), 0, 0, "*metmuonjes_allj*", kTRUE);
  std::cout << "adding histograms..." << std::endl;

  hist::add("sm_dilPt_allj_all", "^[^L].*_hdilPt_allj_all$");
  hist::add("sm_dilPtSmeared_allj_all", "^[^L].*_hdilPtSmeared_allj_all$");
  //hist::add("sm_dilPt_2j_all", "^LM11_hdilPt_2j_all$");
  hist::add("sm_tcmet_allj_all", "^[^L].*_htcmet_allj_all$");
  hist::add("sm_tcmet_allj_ee",  "^[^L].*_htcmet_allj_ee$");
  hist::add("sm_tcmet_allj_mm",  "^[^L].*_htcmet_allj_mm$");
  hist::add("sm_tcmet_allj_em",  "^[^L].*_htcmet_allj_em$");
  //hist::add("sm_metmuonjes_2j_all", "^[^L].*_hmetmuonjes_2j_all$");

  TH1F* sm_dilPt =     (TH1F*)gDirectory->Get("sm_dilPt_allj_all");
  TH1F* sm_dilPtSmeared = (TH1F*)gDirectory->Get("sm_dilPtSmeared_allj_all");
  TH1F* sm_tcmet_all = (TH1F*)gDirectory->Get("sm_tcmet_allj_all");
  TH1F* sm_tcmet_ee =  (TH1F*)gDirectory->Get("sm_tcmet_allj_ee");
  TH1F* sm_tcmet_mm =  (TH1F*)gDirectory->Get("sm_tcmet_allj_mm");
  TH1F* sm_tcmet_em =  (TH1F*)gDirectory->Get("sm_tcmet_allj_em");

  if(sm_dilPt == 0){
    cout<<"Error unable to find histos sm_dilPt_allj_all"<<endl;
    return;
  }
  if(sm_dilPtSmeared == 0){
    cout<<"Error unable to find histos sm_dilPtSmeared_allj_all"<<endl;
    return;
  }
  if(sm_tcmet_all == 0){
    cout<<"Error unable to find histos sm_tcmet_allj_all"<<endl;
    return;
  }
  if(sm_tcmet_ee == 0){
    cout<<"Error unable to find histos sm_tcmet_allj_ee"<<endl;
    return;
  }
  if(sm_tcmet_mm == 0){
    cout<<"Error unable to find histos sm_tcmet_allj_mm"<<endl;
    return;
  }
  if(sm_tcmet_em == 0){
    cout<<"Error unable to find histos sm_tcmet_allj_em"<<endl;
    return;
  }
  
  //SM+LM0 dilPt hist
  TH1F* sm_LM0_dilPt = (TH1F*) sm_dilPt->Clone();
  TH1F* hLM0_dilPt   = (TH1F*) gDirectory->Get("LM0_hdilPt_allj_all");
  sm_LM0_dilPt->Add(hLM0_dilPt);

  //SM+LM1 dilPt hist
  TH1F* sm_LM1_dilPt = (TH1F*) sm_dilPt->Clone();
  TH1F* hLM1_dilPt   = (TH1F*) gDirectory->Get("LM1_hdilPt_allj_all");
  sm_LM1_dilPt->Add(hLM1_dilPt);

  //SM+LM0 dilPtSmeared hist
  TH1F* sm_LM0_dilPtSmeared = (TH1F*) sm_dilPtSmeared->Clone();
  TH1F* hLM0_dilPtSmeared   = (TH1F*) gDirectory->Get("LM0_hdilPtSmeared_allj_all");
  sm_LM0_dilPtSmeared->Add(hLM0_dilPtSmeared);

  //SM+LM1 dilPtSmeared hist
  TH1F* sm_LM1_dilPtSmeared = (TH1F*) sm_dilPtSmeared->Clone();
  TH1F* hLM1_dilPtSmeared   = (TH1F*) gDirectory->Get("LM1_hdilPtSmeared_allj_all");
  sm_LM1_dilPtSmeared->Add(hLM1_dilPtSmeared);

  //SM+LM0 tcmet hist
  TH1F* sm_LM0_tcmet = (TH1F*) sm_tcmet_all->Clone();
  TH1F* hLM0_tcmet   = (TH1F*) gDirectory->Get("LM0_htcmet_allj_all");
  sm_LM0_tcmet->Add(hLM0_tcmet);

  //SM+LM1 tcmet hist
  TH1F* sm_LM1_tcmet = (TH1F*) sm_tcmet_all->Clone();
  TH1F* hLM1_tcmet   = (TH1F*) gDirectory->Get("LM1_htcmet_allj_all");
  sm_LM1_tcmet->Add(hLM1_tcmet);


  //SM only prediction
  std::cout << "making predictions..." << std::endl;

  int bin50  = sm_dilPt->FindBin(50);
  int bin100 = sm_dilPt->FindBin(100);
  int bin175 = sm_dilPt->FindBin(175);
  double scale = sm_dilPt->Integral(bin50, 101) / sm_dilPt->Integral(0, 101);
  sm_dilPt->Scale( 1. / scale );
  double dilpt_prediction100 = sm_dilPt->Integral(bin100, 101);
  double dilpt_prediction175 = sm_dilPt->Integral(bin175, 101);

  int bin50_smeared  = sm_dilPtSmeared->FindBin(50);
  int bin100_smeared = sm_dilPtSmeared->FindBin(100);
  int bin175_smeared = sm_dilPtSmeared->FindBin(175);
  double scale_smeared = sm_dilPtSmeared->Integral(bin50_smeared, 101) / sm_dilPtSmeared->Integral(0, 101);
  sm_dilPtSmeared->Scale( 1. / scale_smeared );
  double dilptSmeared_prediction100 = sm_dilPtSmeared->Integral(bin100_smeared, 101);
  double dilptSmeared_prediction175 = sm_dilPtSmeared->Integral(bin175_smeared, 101);

  std::cout << "the scale factor is: " << scale << std::endl;
  std::cout << "the smeared scale factor is: " << scale_smeared << std::endl;

  //SM+LM0 prediction
  double scaleLM0 = sm_LM0_dilPt->Integral(bin50, 101) / sm_LM0_dilPt->Integral(0, 101);
  sm_LM0_dilPt->Scale(1. / scaleLM0);
  double dilpt_LM0_prediction100 = sm_LM0_dilPt->Integral(bin100, 101);
  double dilpt_LM0_prediction175 = sm_LM0_dilPt->Integral(bin175, 101);

  //SM+LM1 prediction
  double scaleLM1 = sm_LM1_dilPt->Integral(bin50, 101) / sm_LM1_dilPt->Integral(0, 101);
  sm_LM1_dilPt->Scale(1. / scaleLM1);
  double dilpt_LM1_prediction100 = sm_LM1_dilPt->Integral(bin100, 101);
  double dilpt_LM1_prediction175 = sm_LM1_dilPt->Integral(bin175, 101);

  //SM+LM0 smeared prediction
  double scaleLM0_smeared = sm_LM0_dilPtSmeared->Integral(bin50_smeared, 101) / sm_LM0_dilPtSmeared->Integral(0, 101);
  sm_LM0_dilPtSmeared->Scale(1. / scaleLM0_smeared);
  double dilptSmeared_LM0_prediction100 = sm_LM0_dilPtSmeared->Integral(bin100_smeared, 101);
  double dilptSmeared_LM0_prediction175 = sm_LM0_dilPtSmeared->Integral(bin175_smeared, 101);

  //SM+LM1 smeared prediction
  double scaleLM1_smeared = sm_LM1_dilPtSmeared->Integral(bin50_smeared, 101) / sm_LM1_dilPtSmeared->Integral(0, 101);
  sm_LM1_dilPtSmeared->Scale(1. / scaleLM1_smeared);
  double dilptSmeared_LM1_prediction100 = sm_LM1_dilPtSmeared->Integral(bin100_smeared, 101);
  double dilptSmeared_LM1_prediction175 = sm_LM1_dilPtSmeared->Integral(bin175_smeared, 101);

  //draw plot
  TCanvas* c1 = new TCanvas;
  c1->cd();
  sm_dilPt->SetMarkerColor(kBlue);
  sm_dilPt->SetLineColor(kBlue);
  sm_dilPt->SetMarkerStyle(22);
  sm_dilPt->Draw();

  TH1F* craphisto = new TH1F();
  craphisto->SetFillColor(0);
  craphisto->SetMarkerColor(0);
  craphisto->SetLineColor(0);
  craphisto->Draw("sames");
  
  sm_tcmet_all->SetMarkerColor(kRed);
  sm_tcmet_all->SetLineColor(kRed);
  sm_tcmet_all->SetMarkerStyle(24);
  sm_tcmet_all->Draw("sames");
  
  sm_dilPt->SetXTitle("p_{T} (GeV)");
  sm_dilPt->SetTitle("");
  
  TLegend* legend = hist::legend(c1, "p");
  TList* list = legend->GetListOfPrimitives();
  TLegendEntry* leg1 = (TLegendEntry*)list->At(0);

  leg1->SetLabel("Predicted SM MET distribution");
  TLegendEntry* leg3 = (TLegendEntry*)list->At(1);
  leg3->SetLabel("from dilepton p_{T} distribution");
  TLegendEntry* leg2 = (TLegendEntry*)list->At(2);
  leg2->SetLabel("Observed SM MET distribution");
  
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  legend->Draw();

  std::cout << dilpt_prediction100 << " events predicted from dilepton pt distribution for " 
            << " tcmet > 100 GeV" << std::endl;
  std::cout << dilptSmeared_prediction100 << " events predicted from smeared dilepton pt distribution for " 
            << " tcmet > 100 GeV" << std::endl;
  std::cout << dilpt_prediction175 << " events predicted from dilepton pt distribution for " 
            << " tcmet > 175 GeV" << std::endl;
  std::cout << dilptSmeared_prediction175 << " events predicted from smeared dilepton pt distribution for " 
            << " tcmet > 175 GeV" << std::endl;
  
  //get observed number of events above 100, 175 GeV
  std::cout << "working on tcmet/jpts..." << std::endl;
  int metbin100 = sm_tcmet_all->FindBin(100);
  int metbin175 = sm_tcmet_all->FindBin(175);

  double met_observed100 = sm_tcmet_all->Integral(metbin100, 101);
  double met_observed175 = sm_tcmet_all->Integral(metbin175, 101);

  double met_LM0_observed100 = sm_LM0_tcmet->Integral(metbin100, 101);
  double met_LM0_observed175 = sm_LM0_tcmet->Integral(metbin175, 101);

  double met_LM1_observed100 = sm_LM1_tcmet->Integral(metbin100, 101);
  double met_LM1_observed175 = sm_LM1_tcmet->Integral(metbin175, 101);

  std::cout << met_observed100 << " events observed in " 
            << " tcmet distribution above 100 GeV for SM." << std::endl;  
  std::cout << met_observed175 << " events observed in " 
            << " tcmet distribution above 175 GeV for SM." << std::endl;  

  //print out predictions and observations
  cout<<endl<<endl<<endl;
  cout<<"---------------------------------------------------------------"<<endl;
  cout<<"|                         |  met > 100 GeV  |  met > 175 GeV  |"<<endl;
  cout<<"|"<<setw(10)<<"SM Observed              |"<<setprecision(3)<<setw(15)<<met_observed100<<"  |"<<setw(15)<<met_observed175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM Prediction            |"<<setprecision(3)<<setw(15)<<dilpt_prediction100<<"  |"<<setw(15)<<dilpt_prediction175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM Smeared Prediction    |"<<setprecision(3)<<setw(15)<<dilptSmeared_prediction100<<"  |"<<setw(15)<<dilptSmeared_prediction175<<"  |"<<endl;
  cout<<"|                         |"<<setw(18)<<"|"<<setw(18)<<"|"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM0 Observed          |"<<setprecision(3)<<setw(15)<<met_LM0_observed100<<"  |"<<setw(15)<<met_LM0_observed175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM0 Prediction        |"<<setprecision(3)<<setw(15)<<dilpt_LM0_prediction100<<"  |"<<setw(15)<<dilpt_LM0_prediction175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM0 Smeared Prediction|"<<setprecision(3)<<setw(15)<<dilptSmeared_LM0_prediction100<<"  |"<<setw(15)<<dilptSmeared_LM0_prediction175<<"  |"<<endl;
  cout<<"|                         |"<<setw(18)<<"|"<<setw(18)<<"|"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM1 Observed          |"<<setprecision(3)<<setw(15)<<met_LM1_observed100<<"  |"<<setw(15)<<met_LM1_observed175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM1 Prediction        |"<<setprecision(3)<<setw(15)<<dilpt_LM1_prediction100<<"  |"<<setw(15)<<dilpt_LM1_prediction175<<"  |"<<endl;
  cout<<"|"<<setw(10)<<"SM+LM1 Smeared Prediction|"<<setprecision(3)<<setw(15)<<dilptSmeared_LM1_prediction100<<"  |"<<setw(15)<<dilptSmeared_LM1_prediction175<<"  |"<<endl;
  cout<<"---------------------------------------------------------------"<<endl;

  //make yield tables
  //makeTable(filename, samples, 100, sm_tcmet_all, sm_tcmet_ee, sm_tcmet_mm, sm_tcmet_em);
  //makeTable(filename, samples, 175, sm_tcmet_all, sm_tcmet_ee, sm_tcmet_mm, sm_tcmet_em);
  
  //make yield tables (side-by-side)
  TFile *f=TFile::Open(filename);

  const unsigned int nsamples=samples.size();
  char* leptype[4]={"all","ee","mm","em"};
  TH1F *h;
  TH1F* hall = (TH1F*)sm_tcmet_all->Clone();
  TH1F* hee  = (TH1F*)sm_tcmet_ee->Clone();
  TH1F* hmm  = (TH1F*)sm_tcmet_mm->Clone();
  TH1F* hem  = (TH1F*)sm_tcmet_em->Clone();

  cout<<endl<<endl<<endl;
  cout<<"------------------------------------------------------------------------------------------------"
      <<"-------------------------------------------------------------------------------"<<endl;

  cout<<"|              |  met > 100 GeV    |  met > 100 GeV    |  met > 100 GeV    |  met > 100 GeV    |"
      <<"  met > 175 GeV    |  met > 175 GeV    |  met > 175 GeV    |  met > 175 GeV   |"<<endl;

  cout<<"|    Sample    |            all    |             ee    |             mm    |             em    |"
      <<"            all    |             ee    |             mm    |             em   |"<<endl;
  
  bin100 = hall->FindBin(100);
  bin175 = hall->FindBin(175);
  int maxbin = hall->GetNbinsX()+1;
  cout<<"|"<<setw(10)<<"SM TOT";
  cout<<setprecision(3)<<"    |"<<setw(15)<<hall->Integral(bin100,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hee-> Integral(bin100,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hmm-> Integral(bin100,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hem-> Integral(bin100,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hall->Integral(bin175,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hee-> Integral(bin175,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hmm-> Integral(bin175,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hem-> Integral(bin175,maxbin)<<"   |"<<endl;

  for(unsigned int isample = 0 ; isample < nsamples ; isample++){

    h = (TH1F*) f->Get(Form("%s_htcmet_allj_all",samples.at(isample),"all"));
    if(h==0) continue;
    delete h;
    
    cout<<"|"<<setw(10)<<samples[isample];

    for(int ilep = 0 ; ilep < 4 ; ilep++){
      h = (TH1F*) f->Get(Form("%s_htcmet_allj_%s",samples.at(isample),leptype[ilep]));
      
      cout<<setprecision(3)<<"    |"<<setw(15)<<h->Integral(bin100,maxbin);
      
      delete h;
    }

    for(int ilep = 0 ; ilep < 4 ; ilep++){
      h = (TH1F*) f->Get(Form("%s_htcmet_allj_%s",samples.at(isample),leptype[ilep]));
      
      cout<<setprecision(3)<<"    |"<<setw(15)<<h->Integral(bin175,maxbin);
      
      delete h;
    }

    cout<<"   |"<<endl;
  }

  cout<<"------------------------------------------------------------------------------------------------"
      <<"-------------------------------------------------------------------------------"<<endl;

 
  delete f;

}


void makeTable(char* filename, vector<char*> samples, float metcut,
               TH1F* hall, TH1F* hee, TH1F* hmm, TH1F* hem){
  
  TFile *f=TFile::Open(filename);
  
  const unsigned int nsamples=samples.size();
  char* leptype[4]={"all","ee","mm","em"};
  TH1F *h;
  
  cout<<endl<<"Making yield table for tcmet > "<<metcut<<" GeV"<<endl;
  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
  cout<<"|    Sample    |            all    |             ee    |             mm    |             em   |"<<endl;
  
  int minbin = hall->FindBin(metcut);
  int maxbin = hall->GetNbinsX()+1;
  cout<<"|"<<setw(10)<<"SM TOT";
  cout<<setprecision(3)<<"    |"<<setw(15)<<hall->Integral(minbin,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hee->Integral(minbin,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hmm->Integral(minbin,maxbin);
  cout<<setprecision(3)<<"    |"<<setw(15)<<hem->Integral(minbin,maxbin)<<"   |"<<endl;

  for(unsigned int isample = 0 ; isample < nsamples ; isample++){

    h = (TH1F*) f->Get(Form("%s_htcmet_allj_all",samples.at(isample),"all"));
    if(h==0) continue;
    delete h;
    
    cout<<"|"<<setw(10)<<samples[isample];

    for(int ilep = 0 ; ilep < 4 ; ilep++){
      h = (TH1F*) f->Get(Form("%s_htcmet_allj_%s",samples.at(isample),leptype[ilep]));
      
      cout<<setprecision(3)<<"    |"<<setw(15)<<h->Integral(minbin,maxbin);
      
      delete h;
    }

    cout<<"   |"<<endl;
  }

  cout<<"-----------------------------------------------------------------------------------------------"<<endl;
 
  delete f;
}


































//   for( unsigned int i = nSMSamples; i < 19; ++i ) {

//     TH1F* hmet   = (TH1F*)gDirectory->Get(Form("%s%s%s%s", samples[i], "_h", met.c_str(), "_2j_all"));

//     if (hmet   == 0) return;

//     int lmbin100 = hmet->FindBin(100);
//     int lmbin175 = hmet->FindBin(175);
//     double lm_observed100 = hmet->Integral(lmbin100, 101);
//     double lm_observed175 = hmet->Integral(lmbin175, 101);

//     std::cout << setw(15) << samples[i] << setw(25) << setprecision(3) << lm_observed100+met_observed100 
//               << setw(25) << setprecision(3) << dilpt_prediction100 << setw(25) << setprecision(3) 
//               << lm_observed175+met_observed175 << setw(25) << setprecision(3) << dilpt_prediction175 << std::endl;

//   }

//   std::cout << "\n\n" << std::endl;  
//   std::cout << setw(15) << "Sample" << setw(25) << "Observed (100)" << setw(25) << "Observed (175)" << std::endl;

//   for( unsigned int j = 0; j < nSMSamples; ++j ) {

//     TH1F* hmet   = (TH1F*)gDirectory->Get(Form("%s%s%s%s", samples[j], "_h", met.c_str(), "_2j_all"));

//     if (hmet   == 0) return;

//     int smbin100 = hmet->FindBin(100);
//     int smbin175 = hmet->FindBin(175);
//     double sm_observed100 = hmet->Integral(smbin100, 101);
//     double sm_observed175 = hmet->Integral(smbin175, 101);

//     std::cout << setw(15) << samples[j] << setw(25) << setprecision(3) << sm_observed100 << setw(25) << setprecision(3) << sm_observed175 << std::endl;


//   }
