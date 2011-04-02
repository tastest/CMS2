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
#include "THStack.h"
#include "TStyle.h"

using namespace std;

//variables to plot
vector<char*> vars;
vector<char*> xtitles; 
vector<char*> smprefix;
vector<char*> susyprefix;

int colors[]={5,4,2,3,6,7,8,9,10,11,12};
int susycolors[]={1,2,4,6};

void initialize(){

  //variables to plot
  //vars.push_back("htcmet");    xtitles.push_back("tcmet (GeV)");
  vars.push_back("hdilMass");  xtitles.push_back("dilepton mass (GeV)");
  vars.push_back("hdilPt");    xtitles.push_back("dilepton p_{T} (GeV)");
  //vars.push_back("hnJet");     xtitles.push_back("Jet Multiplicity");
  //vars.push_back("hdildR");    xtitles.push_back("#DeltaR(ll)");
  //vars.push_back("hmlj3");     xtitles.push_back("M(lj)3");
  //vars.push_back("hdphiLep");  xtitles.push_back("#Delta#phi(ll)");
  //vars.push_back("hmt2j");     xtitles.push_back("MT2Jets");
  
  //SM samples to include
  smprefix.push_back("ttdil");
  smprefix.push_back("ww"); 
  smprefix.push_back("wz");
  smprefix.push_back("zz");
  smprefix.push_back("wjets");
  smprefix.push_back("Zjets");   
  smprefix.push_back("tW"); 
  
  //SUSY samples to include
  susyprefix.push_back("LM0");
  susyprefix.push_back("LM1");
  susyprefix.push_back("LM4");
  susyprefix.push_back("LM9");
  susyprefix.push_back("LM10");
}

TH1F* getCloneHist(TH1F* hin, int color);
TH1F* getFSHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color);
TH1F* getSFHist(TH1F* hee, TH1F* hmm, int color);
TH1F* getSUSYCloneHist(TH1F* hin, int color);
TH1F* getFSSUSYHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color);
void formatHist(THStack *stack,char* var, char* histtype);
TLegend *getLegend();
THStack* getSMStack  (TFile *file, char* varname, char* xtitle, char* histtype);
void drawSUSYHists(TFile *file, char* varname, char* histtype);


void drawPlots(char* filename , bool printgif = false){
  
  initialize();

  const unsigned int nVars = vars.size();
  assert(vars.size() == xtitles.size());
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open(filename);
  file->cd();
  
  TCanvas *canArray[nVars];
  TCanvas *canvas;

  TLegend *leg = getLegend();
  
  //loop over variables-----------------------------------------------------
  for(unsigned int iVar = 0 ; iVar < nVars ; iVar++ ){

    char* allj="allj_";
    if(strcmp(vars[iVar],"hnJet")==0)allj="";

    canArray[iVar] = new TCanvas(Form("%s_can",vars.at(iVar)),Form("%s_can",vars.at(iVar)),1600,600);
    canArray[iVar]->Divide(2,1);
    
    //draw hists with all dilepton flavors combined
    canArray[iVar]->cd(1);
    THStack *stack=getSMStack(file,vars.at(iVar),xtitles.at(iVar),"all");
    stack->Draw();
    if(iVar==0)leg->Draw();
    drawSUSYHists(file, vars.at(iVar),"all");
   
    //draw flavor-subtracted hists
    canArray[iVar]->cd(2);
    THStack *stackFS = getSMStack(file,vars.at(iVar),xtitles.at(iVar),"fs");
    stackFS->Draw();
    drawSUSYHists(file, vars.at(iVar),"fs");
    
    if(printgif) canArray[iVar]->Print(Form("plots/%s.gif",vars.at(iVar)));
    //if(histmin < stackFS->GetMinimum() ) stackFS->SetMinimum(histmin);
    //if(histmax > stackFS->GetMaximum() ) stackFS->SetMaximum(histmax);
        
  }

}

void formatHist(THStack *stack,char* var, char* histtype){
  if(strcmp(var,"hnJet")==0 && strcmp(histtype,"fs")==0){
    stack->SetMinimum(-10);
    stack->SetMaximum(15);
  }
  if(strcmp(var,"hdilPt")==0 && strcmp(histtype,"fs")==0){
    stack->SetMinimum(-1.5);
    stack->SetMaximum(4);
  }
  if(strcmp(var,"hmt2j")==0 && strcmp(histtype,"fs")==0){
    stack->SetMinimum(-1.2);
    stack->SetMaximum(2.7);
  }
  //if(strcmp(var,"hdilMass")==0 && strcmp(histtype,"fs")==0){
  //  stack->SetMinimum(-20);
  //  stack->SetMaximum(20);
  //}
  if(strcmp(var,"hdilMass")==0 && (strcmp(histtype,"sf")==0 || strcmp(histtype,"em")==0)){
    stack->SetMaximum(20);
  }
  if(strcmp(var,"hdildR")==0 && (strcmp(histtype,"all")==0)){
    stack->SetMaximum(10);
  }
  if(strcmp(var,"hdildR")==0 && (strcmp(histtype,"fs")==0)){
    stack->SetMaximum(6);
    stack->SetMinimum(-1.5);
  }
  if(strcmp(var,"hdphiLep")==0 && (strcmp(histtype,"all")==0)){
    stack->SetMaximum(15);
  }
  if(strcmp(var,"hdphiLep")==0 && (strcmp(histtype,"fs")==0)){
    stack->SetMinimum(-2);
    stack->SetMaximum(7);
  }
}

TH1F* getFSHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color){
  TH1F* h1 = getCloneHist(hee,color);
  TH1F* h2 = getCloneHist(hmm,color);
  TH1F* h3 = getCloneHist(hem,color);

  h1->Add(h2,1);
  h1->Add(h3,-1);

  return h1;

}

TH1F* getSFHist(TH1F* hee, TH1F* hmm, int color){
  TH1F* h1 = getCloneHist(hee,color);
  TH1F* h2 = getCloneHist(hmm,color);
  
  h1->Add(h2,1);
  
  return h1;

}


TH1F* getFSSUSYHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color){
  TH1F* h1 = getSUSYCloneHist(hee,color);
  TH1F* h2 = getSUSYCloneHist(hmm,color);
  TH1F* h3 = getSUSYCloneHist(hem,color);
  
  h1->Add(h2,1);
  h1->Add(h3,-1);
  
  return h1;
  
}

TH1F* getSFSUSYHist(TH1F* hee, TH1F* hmm, int color){
  TH1F* h1 = getSUSYCloneHist(hee,color);
  TH1F* h2 = getSUSYCloneHist(hmm,color);
  
  h1->Add(h2,1);
  
  return h1;
  
}

TH1F* getCloneHist(TH1F* hin, int color){
  TH1F* hout=new TH1F(hin->GetName(),hin->GetTitle(),
		      hin->GetNbinsX(),hin->GetXaxis()->GetXmin(),hin->GetXaxis()->GetXmax());

  for(int ibin=1;ibin<=hin->GetNbinsX();ibin++)
    hout->SetBinContent(ibin,hin->GetBinContent(ibin));

    hout->SetLineColor(1);
    hout->SetMarkerColor(color);
    hout->SetFillColor(color);

  return hout;
}

TH1F* getSUSYCloneHist(TH1F* hin, int color){
  TH1F* hout=new TH1F(hin->GetName(),hin->GetTitle(),
		      hin->GetNbinsX(),hin->GetXaxis()->GetXmin(),hin->GetXaxis()->GetXmax());
  
  for(int ibin=1;ibin<=hin->GetNbinsX();ibin++){
    hout->SetBinContent(ibin,hin->GetBinContent(ibin));
    hout->SetBinError  (ibin,hin->GetBinError(ibin));
  }
  hout->SetMarkerColor(color);
  hout->SetLineColor(color);
  return hout;
}

TLegend *getLegend(){

  const unsigned int nSM = smprefix.size();
  const unsigned int nSUSY = susyprefix.size();

  TLegend *leg =new TLegend(0.7,0.5,0.85,0.85);
  
  for(unsigned int i = 0 ; i < nSM ; i++){
    TH1* hdummy=new TH1("hSMdummy","",1,0,1);
    hdummy->SetLineColor(1);
    hdummy->SetMarkerColor(colors[i]);
    hdummy->SetFillColor(colors[i]);
    leg->AddEntry(hdummy,smprefix.at(i),"f");
    
  }
  for(unsigned int i=0;i<nSUSY;i++){
    TH1* hdummy=new TH1("hdummy","",1,0,1);
    hdummy->SetLineColor(i+1);
    hdummy->SetMarkerColor(i+1);
    leg->AddEntry(hdummy,susyprefix.at(i));
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  return leg;
}

int rebin(char* prefix){
  int r = 1;
  //if(strcmp(prefix,"hmt2j")==0) r = 10;

  return r;

}
THStack* getSMStack(TFile *file, char* varname, char* xtitle,char* histtype){
  
  const unsigned int nSM = smprefix.size();
  
  THStack *stack=new THStack("stack","");
  TH1F* h[nSM];

  char* allj="allj_";
  if(strcmp(varname,"hnJet")==0)allj="";



  for(unsigned int iSM=0;iSM<nSM;iSM++){
    if(strcmp(histtype,"all")==0 || strcmp(histtype,"ee")==0 || 
       strcmp(histtype,"mm")==0  || strcmp(histtype,"em")==0){
      
      cout<<"Getting "<<Form("%s_%s_%s%s",smprefix.at(iSM),varname,allj,histtype)<<endl;
      h[iSM]  = getCloneHist((TH1F*)(file->Get(Form("%s_%s_%s%s",smprefix.at(iSM),varname,allj,histtype))), 
			     colors[iSM]);
      
    }
    if(strcmp(histtype,"fs")==0){

      cout<<"Getting FS "<<Form("%s_%s_allj_all",smprefix.at(iSM),varname)<<endl;
      h[iSM]  = getFSHist((TH1F*)(file->Get(Form("%s_%s_%see",smprefix.at(iSM),varname,allj)))->Clone(),
			  (TH1F*)(file->Get(Form("%s_%s_%smm",smprefix.at(iSM),varname,allj)))->Clone(),
			  (TH1F*)(file->Get(Form("%s_%s_%sem",smprefix.at(iSM),varname,allj)))->Clone(),
			  colors[iSM]);

     
    }
    if(strcmp(histtype,"sf")==0){
      
      cout<<"Getting SF "<<Form("%s_%s_allj_all",smprefix.at(iSM),varname)<<endl;
      h[iSM]  = getSFHist((TH1F*)(file->Get(Form("%s_%s_%see",smprefix.at(iSM),varname,allj)))->Clone(),
			  (TH1F*)(file->Get(Form("%s_%s_%smm",smprefix.at(iSM),varname,allj)))->Clone(),
			  colors[iSM]);


      
    }
    
    
    if(rebin(varname)>1) h[iSM]->Rebin(rebin(varname));
    stack->Add(h[iSM]);
  }
  
  if(strcmp(histtype,"all")==0) stack->SetTitle(Form("%s (ee + #mu#mu + e#mu)",xtitle));
  if(strcmp(histtype,"fs")==0)  stack->SetTitle(Form("%s (ee + #mu#mu - e#mu)",xtitle));
  if(strcmp(histtype,"sf")==0)  stack->SetTitle(Form("%s (ee + #mu#mu)",xtitle));
  if(strcmp(histtype,"ee")==0)  stack->SetTitle(Form("%s (ee)",xtitle));
  if(strcmp(histtype,"mm")==0)  stack->SetTitle(Form("%s (#mu#mu)",xtitle));
  if(strcmp(histtype,"em")==0)  stack->SetTitle(Form("%s (e#mu)",xtitle));

  formatHist(stack,varname,histtype);

  return stack;
}



void drawSUSYHists(TFile *file, char* varname, char* histtype){

  //histmin = 0.;
  //histmax = 0.;
  const unsigned int nSUSY = susyprefix.size();


  TH1F* hsusy[nSUSY];
  char* allj="allj_";
  if(strcmp(varname,"hnJet")==0)allj="";
  
  for(unsigned int iSUSY=0;iSUSY<nSUSY;iSUSY++){
    if(strcmp(histtype,"all")==0 || strcmp(histtype,"ee")==0 || 
       strcmp(histtype,"mm")==0  || strcmp(histtype,"em")==0){
      hsusy[iSUSY]  = getSUSYCloneHist((TH1F*)file->Get(Form("%s_%s_%s%s",susyprefix.at(iSUSY),varname,allj,histtype))->Clone(),susycolors[iSUSY]);
    }
    if(strcmp(histtype,"fs")==0){ 
      hsusy[iSUSY]  = getFSSUSYHist((TH1F*)(file->Get(Form("%s_%s_%see",susyprefix.at(iSUSY),varname,allj)))->Clone(),
				    (TH1F*)(file->Get(Form("%s_%s_%smm",susyprefix.at(iSUSY),varname,allj)))->Clone(),
				    (TH1F*)(file->Get(Form("%s_%s_%sem",susyprefix.at(iSUSY),varname,allj)))->Clone(),
				    susycolors[iSUSY]);
    }
    if(strcmp(histtype,"sf")==0){ 
      hsusy[iSUSY]  = getSFSUSYHist((TH1F*)(file->Get(Form("%s_%s_%see",susyprefix.at(iSUSY),varname,allj)))->Clone(),
				    (TH1F*)(file->Get(Form("%s_%s_%smm",susyprefix.at(iSUSY),varname,allj)))->Clone(),
				    susycolors[iSUSY]);
    }

    //cout<<"integral "<<varname<<" "<<histtype<<" "<<iSUSY<<" "<<hsusy[iSUSY]->Integral()<<endl;
    if(rebin(varname)>1) hsusy[iSUSY]->Rebin(rebin(varname));
    hsusy[iSUSY]->Draw("sameE1");
    
    //if( hsusy[iSUSY]->GetMinimum() < histmin ) histmin = hsusy[iSUSY]->GetMinimum();
    //if( hsusy[iSUSY]->GetMaximum() > histmax ) histmax = hsusy[iSUSY]->GetMaximum();
  }
}

