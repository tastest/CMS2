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
#include "TLine.h"

using namespace std;

//variables to plot
vector<char*> vars;
vector<char*> xtitles; 
vector<char*> smprefix;
vector<char*> susyprefix;

int colors[]={9,4,2,3,6,7,8,5,10,11,12};
int susycolors[]={1,2,4,6};

void initialize(){

  //------------------------
  //variables to plot
  //------------------------
  
  vars.push_back("htcmet");      xtitles.push_back("tcmet (GeV)");
  vars.push_back("hsumJetPt");   xtitles.push_back("scalar sum jet p_{T} (GeV)");
  vars.push_back("hdilMass");   xtitles.push_back("dilepton mass (GeV)");
  vars.push_back("hdilPt");     xtitles.push_back("dilepton p_{T} (GeV)");
  //vars.push_back("hminLepPt");  xtitles.push_back("min lepton p_{T} (GeV)");
  //vars.push_back("hmaxLepPt");  xtitles.push_back("max lepton p_{T} (GeV)");
  //vars.push_back("hptJet1");    xtitles.push_back("leading jet p_{T} (GeV)");
  //vars.push_back("hptJet2");    xtitles.push_back("2nd leading jet p_{T} (GeV)");
  //vars.push_back("hptJet3");    xtitles.push_back("3rd leading jet p_{T} (GeV)");
  //vars.push_back("hptJet4");    xtitles.push_back("4th leading jet p_{T} (GeV)");
  //vars.push_back("hmeffJet");   xtitles.push_back("effective mass (GeV)");
  //vars.push_back("hmt");        xtitles.push_back("transverse mass (GeV)");
  vars.push_back("hnJet");     xtitles.push_back("Jet Multiplicity");
  //vars.push_back("hdildR");    xtitles.push_back("#DeltaR(ll)");
  //vars.push_back("hmlj3");     xtitles.push_back("M(lj)3");
  //vars.push_back("hdphiLep");  xtitles.push_back("#Delta#phi(ll)");
  //vars.push_back("hmt2j");     xtitles.push_back("MT2Jets");
  
  //------------------------
  //SM samples to include
  //------------------------
  
  smprefix.push_back("ttotr");
  smprefix.push_back("Zjets");     
  smprefix.push_back("wjets");     
  smprefix.push_back("ww"); 
  smprefix.push_back("wz");
  smprefix.push_back("zz");
  smprefix.push_back("tW");
  smprefix.push_back("ttdil");
  
  //------------------------
  //SUSY samples to include
  //------------------------
  
  //susyprefix.push_back("LM0");
  //susyprefix.push_back("LM1");
  //susyprefix.push_back("LM4");
  //susyprefix.push_back("LM9");
  //susyprefix.push_back("LM10");
}

TH1F* getCloneHist(TH1F* hin, int color);
TH1F* getFSHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color);
TH1F* getSFHist(TH1F* hee, TH1F* hmm, int color);
TH1F* getSUSYCloneHist(TH1F* hin, int color);
TH1F* getFSSUSYHist(TH1F* hee, TH1F* hmm, TH1F* hem, int color);
void formatHist(TH1F *hist,     char* var, char* histtype);
TLegend *getLegend();
THStack* getSMStack  (TFile *file, char* varname, char* xtitle, char* histtype);
void drawSUSYHists(TFile *file, char* varname, char* histtype);
TH1F* getDataHist(TFile *file, char* varname, char* xtitle,char* histtype);
void drawOverlayPlot( TH1F* h , THStack* stack );

void compareDataMC(char* filename , bool printgif = false){
  
  initialize();

  const unsigned int nVars = vars.size();
  assert(vars.size() == xtitles.size());
  
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open(filename);
  file->cd();
  
  TCanvas *canArray[nVars];
  //TCanvas *canvas;

  TLegend *leg = getLegend();
  
  //loop over variables-----------------------------------------------------
  for(unsigned int iVar = 0 ; iVar < nVars ; iVar++ ){

    char* allj="allj_";
    if(strcmp(vars[iVar],"hnJet")==0)allj="";

    canArray[iVar] = new TCanvas(Form("%s_can",vars.at(iVar)),Form("%s_can",vars.at(iVar)),1100,650);
    canArray[iVar]->Divide(2,2);
      
    //tot
    canArray[iVar]->cd(4);
    TH1F*    hist  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "all");
    THStack* stack = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "all");
    drawOverlayPlot( hist , stack );
    formatHist( hist , vars.at(iVar) , "all" );
    leg->Draw();

    //ee
    canArray[iVar]->cd(1);
    TH1F*    histee  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "ee");
    THStack* stackee = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "ee");
    drawOverlayPlot( histee , stackee );
    formatHist( histee , vars.at(iVar) , "ee" );
    leg->Draw();

    //em
    canArray[iVar]->cd(3);
    TH1F*    histem  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "em");
    THStack* stackem = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "em");
    drawOverlayPlot( histem , stackem );
    formatHist( histem , vars.at(iVar) , "em" );
    leg->Draw();
    
    //mm
    canArray[iVar]->cd(2);
    TH1F*    histmm  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "mm");
    THStack* stackmm = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "mm");
    drawOverlayPlot( histmm , stackmm );
    formatHist( histmm , vars.at(iVar) , "mm" );
    leg->Draw();
    
    if(printgif) canArray[iVar]->Print(Form("plots/%s.gif",vars.at(iVar)));
    //if(histmin < stackFS->GetMinimum() ) stackFS->SetMinimum(histmin);
    //if(histmax > stackFS->GetMaximum() ) stackFS->SetMaximum(histmax);
        
  }

}

void drawOverlayPlot( TH1F* h , THStack* stack ){
  h->Draw("E1");
  stack->Draw("same");
  h->Draw("sameE1");
  h->Draw("sameaxis");
}

void formatHist(TH1F *hist,char* var, char* histtype){
  TLine line;
  line.SetLineColor(2);
  line.SetLineWidth(2);
  line.SetLineStyle(2);

  if(strcmp(var,"hnJet")==0 && strcmp(histtype,"fs")==0){
    hist->SetMinimum(-10);
    hist->SetMaximum(15);
  }
  if(strcmp(var,"hdilPt")==0 && strcmp(histtype,"fs")==0){
    hist->SetMinimum(-1.5);
    hist->SetMaximum(4);
  }
  if(strcmp(var,"hmt2j")==0 && strcmp(histtype,"fs")==0){
    hist->SetMinimum(-1.2);
    hist->SetMaximum(2.7);
  }
  if(strcmp(var,"hdildR")==0 && (strcmp(histtype,"all")==0)){
    hist->SetMaximum(10);
  }
  if(strcmp(var,"hdildR")==0 && (strcmp(histtype,"fs")==0)){
    hist->SetMaximum(6);
    hist->SetMinimum(-1.5);
  }
  if(strcmp(var,"hdphiLep")==0 && (strcmp(histtype,"all")==0)){
    hist->SetMaximum(15);
  }
  if(strcmp(var,"hdphiLep")==0 && (strcmp(histtype,"fs")==0)){
    hist->SetMinimum(-2);
    hist->SetMaximum(7);
  }
  if( strcmp(var,"htcmet") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 50 , 0 , 50 , 1.05 * max );
  }
  if( strcmp(var,"hsumJetPt") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 100 , 0 , 100 , 1.05 * max );
    hist->GetXaxis()->SetRangeUser(0,500);
  }
  if( strcmp(var,"hdilMass") == 0 && ( histtype == "ee" || histtype == "mm" ) ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 76  , 0 , 76  , 1.05 * max );
    line.DrawLine( 106 , 0 , 106 , 1.05 * max );
  }
  if( strcmp(var,"hnJet") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 1.5 , 0 , 1.5 , 1.05 * max );
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

  TLegend *leg =new TLegend(0.6,0.25,0.8,0.85);
  
  TH1F* hdatadummy = new TH1F("hdatadummy","",1,0,1);
  leg->AddEntry(hdatadummy,"DATA!!!");

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
  if(strcmp(prefix,"hmt")==0)       r = 10;
  if(strcmp(prefix,"htcmet")==0)    r = 5;
  if(strcmp(prefix,"hsumJetPt")==0) r = 5;
  if(strcmp(prefix,"hdilMass")==0)  r = 5;
  if(strcmp(prefix,"hdilPt")==0)    r = 5;

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
  
  //if(strcmp(histtype,"all")==0) stack->SetTitle(Form("%s (ee + #mu#mu + e#mu)",xtitle));
  //stack->GetXaxis()->SetTitle(xtitle);
  if(strcmp(histtype,"all")==0) stack->SetTitle(xtitle);
  if(strcmp(histtype,"fs")==0)  stack->SetTitle(Form("%s (ee + #mu#mu - e#mu)",xtitle));
  if(strcmp(histtype,"sf")==0)  stack->SetTitle(Form("%s (ee + #mu#mu)",xtitle));
  if(strcmp(histtype,"ee")==0)  stack->SetTitle(Form("%s (ee)",xtitle));
  if(strcmp(histtype,"mm")==0)  stack->SetTitle(Form("%s (#mu#mu)",xtitle));
  if(strcmp(histtype,"em")==0)  stack->SetTitle(Form("%s (e#mu)",xtitle));
  
  //formatHist(stack,varname,histtype);
  


  return stack;
}

TH1F* getDataHist(TFile *file, char* varname, char* xtitle,char* histtype){
  
  TH1F* hist = new TH1F();

  char* allj="allj_";
  if(strcmp(varname,"hnJet")==0)allj="";

  
  if(strcmp(histtype,"all")==0 || strcmp(histtype,"ee")==0 || 
     strcmp(histtype,"mm")==0  || strcmp(histtype,"em")==0){
    
    cout<<"Getting "<<Form("%s_%s_%s%s","data",varname,allj,histtype)<<endl;
    hist  = getCloneHist((TH1F*)(file->Get(Form("%s_%s_%s%s","data",varname,allj,histtype))), 
                         1);
    
  }
  if(strcmp(histtype,"fs")==0){
    
    cout<<"Getting FS "<<Form("%s_%s_allj_all","data",varname)<<endl;
    hist  = getFSHist((TH1F*)(file->Get(Form("%s_%s_%see","data",varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%smm","data",varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%sem","data",varname,allj)))->Clone(),
                      1);
    
    
  }
  if(strcmp(histtype,"sf")==0){
    
    cout<<"Getting SF "<<Form("%s_%s_allj_all","data",varname)<<endl;
    hist  = getSFHist((TH1F*)(file->Get(Form("%s_%s_%see","data",varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%smm","data",varname,allj)))->Clone(),
                      1);
    
    
    
  }
  
  if(rebin(varname)>1) hist->Rebin(rebin(varname));
  
  //if(strcmp(histtype,"all")==0) hist->SetTitle(Form("%s (ee + #mu#mu + e#mu)",xtitle));
  hist->GetXaxis()->SetTitle(xtitle);
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(1);
  if(strcmp(histtype,"all")==0) hist->SetTitle("tot");
  if(strcmp(histtype,"fs")==0)  hist->SetTitle("ee + #mu#mu - e#mu");
  if(strcmp(histtype,"sf")==0)  hist->SetTitle("ee + #mu#mu)");
  if(strcmp(histtype,"ee")==0)  hist->SetTitle("ee");
  if(strcmp(histtype,"mm")==0)  hist->SetTitle("#mu#mu");
  if(strcmp(histtype,"em")==0)  hist->SetTitle("e#mu");

  //formatHist(hist,varname,histtype);

  return hist;
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

