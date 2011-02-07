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
#include "TLatex.h"

using namespace std;

bool drawSUSY    = true;
const char* data = "dataskim";
//variables to plot
vector<char*> vars;
vector<char*> xtitles; 
vector<char*> smprefix;
vector<char*> smleg;
vector<char*> susyprefix;

int colors[]={0,3,7,4,2,5,-1,11,12};
int susycolors[]={2,4,4,6};

void initialize(){

  //------------------------
  //variables to plot
  //------------------------

  vars.push_back("hsumJptPt");      xtitles.push_back("H_{T} (GeV)");
  vars.push_back("htcmet_sqrtht");  xtitles.push_back("y (GeV^{1/2})"); 
  vars.push_back("hdilMass");       xtitles.push_back("M(ll) (GeV)");
  vars.push_back("hdilPt");         xtitles.push_back("p_{T}(ll) (GeV)");

  /*
  vars.push_back("htcmet");         xtitles.push_back("#slash{E}_{T} (GeV)");
  
  vars.push_back("hsumJptPt");      xtitles.push_back("H_{T} (GeV)");
  vars.push_back("htcmet_sqrtht");  xtitles.push_back("#slash{E}_{T}  /  #sqrt{H_{T}} (GeV^{1/2})"); 
  //vars.push_back("htcmet_sqrtht");  xtitles.push_back("y (GeV^{1/2})"); 
  vars.push_back("hnJpt");          xtitles.push_back("Jet Multiplicity");  
  vars.push_back("hdilMass");       xtitles.push_back("M(ll) (GeV)");
  vars.push_back("hdilPt");         xtitles.push_back("p_{T}(ll) (GeV)");

  
  vars.push_back("hmt2core");       xtitles.push_back("MT2 (GeV)");
  vars.push_back("hmt2jcore");      xtitles.push_back("MT2J (GeV)");
  vars.push_back("htopMass");       xtitles.push_back("top mass estimate (GeV)");
  vars.push_back("hmeffJPT");       xtitles.push_back("effective mass (GeV)");
                                     
  vars.push_back("hmaxLepPt");      xtitles.push_back("leading lepton p_{T} (GeV)"); 
  vars.push_back("hminLepPt");      xtitles.push_back("2nd leading lepton p_{T} (GeV)");
  vars.push_back("hmaxLepEta");     xtitles.push_back("leading lepton #eta");
  vars.push_back("hminLepEta");     xtitles.push_back("2nd leading lepton #eta"); 

  vars.push_back("hptJpt1");        xtitles.push_back("leading jet p_{T} (GeV)");
  vars.push_back("hptJpt2");        xtitles.push_back("2nd leading jet p_{T} (GeV)");
  vars.push_back("hptJpt3");        xtitles.push_back("3rd leading jet p_{T} (GeV)");
  vars.push_back("hptJpt4");        xtitles.push_back("4th leading jet p_{T} (GeV)");

  vars.push_back("hetaJpt1");       xtitles.push_back("leading jet #eta");
  vars.push_back("hetaJpt2");       xtitles.push_back("2nd leading jet #eta");
  vars.push_back("hetaJpt3");       xtitles.push_back("3rd leading jet #eta");
  vars.push_back("hetaJpt4");       xtitles.push_back("4th leading jet #eta");
  
  vars.push_back("hnBtagJpt");      xtitles.push_back("nBtags"); 
  vars.push_back("hptBtagJpt1");    xtitles.push_back("leading b-jet p_{T} (GeV)");
  vars.push_back("hptBtagJpt2");    xtitles.push_back("2nd leading b-jet p_{T} (GeV)");
  vars.push_back("hptBtagJpt3");    xtitles.push_back("3rd leading b-jet p_{T} (GeV)");
  vars.push_back("hptBtagJpt4");    xtitles.push_back("4th leading b-jet p_{T} (GeV)");
  
  vars.push_back("hdrLep");         xtitles.push_back("#DeltaR(lep_{1}lep_{2})"); 
  vars.push_back("hdrJ1J2");        xtitles.push_back("#DeltaR(jet_{1}jet_{2})"); 
  vars.push_back("hdphiLep");       xtitles.push_back("#Delta#phi(ll)");
  */  

  




  //vars.push_back("hmt");          xtitles.push_back("transverse mass (GeV)");
 
  //vars.push_back("hptJet1");      xtitles.push_back("leading jet p_{T} (GeV)");
  //vars.push_back("hptJet2");      xtitles.push_back("2nd leading jet p_{T} (GeV)");
  //vars.push_back("hptJet3");      xtitles.push_back("3rd leading jet p_{T} (GeV)");
  //vars.push_back("hptJet4");      xtitles.push_back("4th leading jet p_{T} (GeV)");
    
  //vars.push_back("hetaJet1");     xtitles.push_back("leading jet #eta");
  //vars.push_back("hetaJet2");     xtitles.push_back("2nd leading jet #eta");
  //vars.push_back("hetaJet3");     xtitles.push_back("3rd leading jet #eta");
  //vars.push_back("hetaJet4");     xtitles.push_back("4th leading jet #eta");


  //vars.push_back("hdildR");    xtitles.push_back("#DeltaR(ll)");
  //vars.push_back("hmlj3");     xtitles.push_back("M(lj)3");
  
  //------------------------
  //SM samples to include
  //------------------------
  

  //smprefix.push_back("Zjets");   smleg.push_back("Z+jets");
  smprefix.push_back("wjets");    smleg.push_back("#font[12]{W} + jets");
  smprefix.push_back("VV");       smleg.push_back("#font[12]{VV}");
  //smprefix.push_back("ww");       smleg.push_back("WW");
  //smprefix.push_back("wz");       smleg.push_back("WZ");
  //smprefix.push_back("zz");       smleg.push_back("ZZ");
  smprefix.push_back("tW");       smleg.push_back("single top");
  smprefix.push_back("DYall");    smleg.push_back("#font[12]{Z^{0} #rightarrow l^{+}l^{-}}");
  smprefix.push_back("ttotr");    smleg.push_back("#font[12]{t#bar{t}#rightarrow}other");
  //smprefix.push_back("ttdil");    smleg.push_back("t#bar{t}#rightarrow l^{+}l^{-}");
  smprefix.push_back("ttdil");    smleg.push_back("#font[12]{t}#bar{#font[12]{t}}#rightarrow #font[12]{l}^{+}#font[12]{l}^{-}");


  //smprefix.push_back("DYee");     smleg.push_back("DYee");
  //smprefix.push_back("DYmm");     smleg.push_back("DYmm");
  //smprefix.push_back("DYtautau"); smleg.push_back("DYtautau");

  //smprefix.push_back("LM0");     smleg.push_back("LM0");
  smprefix.push_back("LM1");     smleg.push_back("LM1");
  
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
char* caption(char* var);

void compareDataMC(char* filename , int print = 0){
  
  initialize();

  const unsigned int nVars = vars.size();
  assert(vars.size() == xtitles.size());
  
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  TFile *file = TFile::Open(filename);
  file->cd();
  
  TCanvas* canArray[nVars];
  TPad*    plotPadArray[nVars];
  TPad*    legPadArray[nVars];
  TPad*    captionPadArray[nVars];

  TCanvas* canvas;

  TLatex t;
  t.SetNDC();

  if( print == 2 ){
    canvas = new TCanvas("canvas","canvas",1100,650);
    gStyle->SetPaperSize(22,28);
    canvas->Print("plots/datamc.ps[");
  }


  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.07);

  TLegend *leg = getLegend();
  
  //loop over variables-----------------------------------------------------
  //for(unsigned int iVar = 0 ; iVar < nVars ; iVar++ ){
  for(unsigned int iVar = 0 ; iVar < 1 ; iVar++ ){

    char* allj="allj_";
    if(strcmp(vars[iVar],"hnJet")==0)     allj="";
    if(strcmp(vars[iVar],"hnJpt")==0)     allj="";
    if(strcmp(vars[iVar],"hnBtagJpt")==0) allj="";

    if( print == 2 ){
      canvas->cd();
    }else{
      canArray[iVar] = new TCanvas(Form("%s_can",vars.at(iVar)),Form("%s_can",vars.at(iVar)),1100,650);
      canArray[iVar]->cd();
    }


    //plotPadArray[iVar] = new TPad(Form("pad_%i",iVar),Form("pad_%i",iVar),0,0.1,0.8,1);
    plotPadArray[iVar] = new TPad(Form("pad_%i",iVar),Form("pad_%i",iVar),0,0,0.8,1);
    plotPadArray[iVar]->Draw();
    plotPadArray[iVar]->cd();
    
    //plotPadArray[iVar]->Divide(2,2);
    
    //plotPadArray[iVar]->cd(1);
    TH1F*    hist0  = getDataHist( file , vars.at(0) , xtitles.at(0) , "all");
    THStack* stack0 = getSMStack(  file , vars.at(0) , xtitles.at(0) , "all");
    drawOverlayPlot( hist0 , stack0 );
    formatHist( hist0 , vars.at(0) , "all" );
    text->DrawLatex(0.8,0.5,"(a)");

    /*
    plotPadArray[iVar]->cd(2);
    TH1F*    hist1  = getDataHist( file , vars.at(1) , xtitles.at(1) , "all");
    THStack* stack1 = getSMStack(  file , vars.at(1) , xtitles.at(1) , "all");
    drawOverlayPlot( hist1 , stack1 );
    formatHist( hist1 , vars.at(1) , "all" );
    text->DrawLatex(0.8,0.5,"(b)");

    plotPadArray[iVar]->cd(3);
    TH1F*    hist2  = getDataHist( file , vars.at(2) , xtitles.at(2) , "all");
    THStack* stack2 = getSMStack(  file , vars.at(2) , xtitles.at(2) , "all");
    drawOverlayPlot( hist2 , stack2 );
    formatHist( hist2 , vars.at(2) , "all" );
    text->DrawLatex(0.8,0.5,"(c)");

    plotPadArray[iVar]->cd(4);
    TH1F*    hist3  = getDataHist( file , vars.at(3) , xtitles.at(3) , "all");
    THStack* stack3 = getSMStack(  file , vars.at(3) , xtitles.at(3) , "all");
    drawOverlayPlot( hist3 , stack3 );
    formatHist( hist3 , vars.at(3) , "all" );
    text->DrawLatex(0.8,0.5,"(d)");    
    */

    /*
    plotPadArray[iVar] = new TPad(Form("pad_%i",iVar),Form("pad_%i",iVar),0,0.1,0.8,1);
    plotPadArray[iVar]->Draw();
    plotPadArray[iVar]->cd();
    plotPadArray[iVar]->Divide(2,2);

    //tot
    //canArray[iVar]->cd(4);
    plotPadArray[iVar]->cd(4);
    TH1F*    hist  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "all");
    THStack* stack = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "all");
    drawOverlayPlot( hist , stack );
    if( drawSUSY ) drawSUSYHists( file, vars.at(iVar) , "all" );
    formatHist( hist , vars.at(iVar) , "all" );
    //leg->Draw();


    //ee
    //canArray[iVar]->cd(1);
    plotPadArray[iVar]->cd(1);
    TH1F*    histee  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "ee");
    THStack* stackee = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "ee");
    drawOverlayPlot( histee , stackee );
    if( drawSUSY ) drawSUSYHists( file, vars.at(iVar) , "ee" );
    formatHist( histee , vars.at(iVar) , "ee" );

    //leg->Draw();

    //em
    //canArray[iVar]->cd(3);
    plotPadArray[iVar]->cd(3);
    TH1F*    histem  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "em");
    THStack* stackem = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "em");
    drawOverlayPlot( histem , stackem );
    if( drawSUSY ) drawSUSYHists( file, vars.at(iVar) , "em" );
    formatHist( histem , vars.at(iVar) , "em" );
    //leg->Draw();
 
    //mm
    //canArray[iVar]->cd(2);
    plotPadArray[iVar]->cd(2);
    TH1F*    histmm  = getDataHist( file , vars.at(iVar) , xtitles.at(iVar) , "mm");
    THStack* stackmm = getSMStack(  file , vars.at(iVar) , xtitles.at(iVar) , "mm");
    drawOverlayPlot( histmm , stackmm );
    if( drawSUSY ) drawSUSYHists( file, vars.at(iVar) , "mm" );
    formatHist( histmm , vars.at(iVar) , "mm" );
    //leg->Draw();
    */

    if( print == 2 ){
      canvas->cd();
    }else{
      canArray[iVar]->cd();
    }
   
    legPadArray[iVar] = new TPad(Form("legpad_%i",iVar),Form("legpad_%i",iVar),0.8,0.1,1,1);
    legPadArray[iVar]->Draw();
    legPadArray[iVar]->cd();
    leg->Draw();

    if( print == 2 ){
      canvas->cd();
    }else{
      canArray[iVar]->cd();
    }
    
    /*
    captionPadArray[iVar] = new TPad(Form("captionpad_%i",iVar),Form("captionpad_%i",iVar),0,0,0.8,0.1);
    captionPadArray[iVar]->Draw();
    captionPadArray[iVar]->cd();
    t.SetTextSize(0.35);
    t.SetTextAlign(22);
    t.DrawLatex(0.5,0.5,Form("Fig. C%i: %s",iVar+1,caption(vars.at(iVar))));
    */
    
    if     ( print == 1 ) canArray[iVar]->Print(Form("plots/%s.png",vars.at(iVar)));
    else if( print == 2 ) {
      canvas->Print("plots/datamc.ps");
      canvas->Clear();
    }
 
    //if(histmin < stackFS->GetMinimum() ) stackFS->SetMinimum(histmin);
    //if(histmax > stackFS->GetMaximum() ) stackFS->SetMaximum(histmax);
        
  }

  if( print == 2 ) {
    canvas->Print("plots/datamc.ps]");
    canvas->Clear();
  }

}

char* caption(char* var){
  if( var == "htcmet" ){
    return "Track-corrected missing transverse energy.";
  }
  if( var == "htcmet_sqrtht" ){
    return "Track-corrected MET divided by square root of scalar sum of jet p_{T}'s.";
  }
  if( var == "hsumJptPt" ){   
    return "Scalar sum of jet transverse momenta.";
  }
  if( var == "hnJpt" ){       
    return "Number of jets.";
  }
  if( var == "hdilMass" ){ 
    return "Dilepton mass.";
  }   
  if( var == "hdilPt" ){  
    return "Transverse momentum of dilepton system.";

  }    
  if( var == "htopMass" ){ 
    return "Reconstructed top mass.";
  }   
  if( var == "hmaxLepPt" ){   
    return "Transverse momentum of leading lepton.";
  }
  if( var == "hminLepPt" ){   
    return "Transverse momentum of non-leading lepton.";
  }
  if( var == "hmaxLepEta" ){  
    return "Pseudo-rapidity of leading lepton.";
  }
  if( var == "hminLepEta" ){  
    return "Pseudo-rapidity of non-leading lepton.";
  }
  if( var == "hptJpt1" ){     
    return "Transverse momentum of leading jet.";
  }
  if( var == "hptJpt2" ){     
    return "Transverse momentum of 2^{nd} leading jet.";
  }
  if( var == "hptJpt3" ){     
    return "Transverse momentum of 3^{rd} leading jet.";
  }
  if( var == "hptJpt4" ){     
    return "Transverse momentum of 4^{th} leading jet.";
  }
  if( var == "hetaJpt1" ){    
    return "Pseudo-rapidity of leading jet.";
  }
  if( var == "hetaJpt2" ){    
    return "Pseudo-rapidity of 2^{nd} leading jet.";
  }
  if( var == "hetaJpt3" ){    
    return "Pseudo-rapidity of 3^{rd} leading jet.";
  }
  if( var == "hetaJpt4" ){    
    return "Pseudo-rapidity of 4^{th} leading jet.";
  }
  if( var == "hnBtagJpt" ){
    return "Number of b-tagged jets (simpleSecondarVertexTagger HE medium WP).";
  }
  if( var == "hptBtagJpt1" ){ 
    return "Transverse momentum of leading b-tagged jet.";
  }
  if( var == "hptBtagJpt2" ){ 
    return "Transverse momentum of 2^{nd} leading b-tagged jet.";
  }
  if( var == "hptBtagJpt3" ){ 
    return "Transverse momentum of 3^{rd} leading b-tagged jet.";
  }
  if( var == "hptBtagJpt4" ){ 
    return "Transverse momentum of 4^{th} leading b-tagged jet.";
  }
  if( var == "hdrLep" ){  
    return "#DeltaR between leptons.";
  }
  if( var == "hdrJ1J2" ){     
    return "#DeltaR between 2 leading jets.";
  }
  if( var == "hdphiLep" ){    
    return "#Delta#phi between leptons.";
  }
  if( var == "hmeffJPT" ){    
    return "Scalar sum of jet and lepton transverse momenta and missing transverse energy.";
  }
  if( var == "hmt" ){         
    return "Transverse mass of leading lepton and missing tranverse energy.";
  }
  if( var == "hmt2core" ){    
    return "MT2";
  }
  if( var == "hmt2jcore" ){   
    return "MT2J";
  }
  return "Error, could not find caption";


}

void drawOverlayPlot( TH1F* h , THStack* stack ){
  h->Draw("E1");
  stack->Draw("same");
  h->Draw("sameE1");
  h->Draw("sameaxis");



  if( stack->GetMaximum() > h->GetMaximum() + h->GetBinError(h->GetMaximumBin() ) ) 
    h->SetMaximum( 1.1 * stack->GetMaximum() );
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
    //hist->GetXaxis()->SetRangeUser(0,500);
  }
  if( strcmp(var,"hsumJptPt") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 100 , 0 , 100 , 1.05 * max );
    //hist->GetXaxis()->SetRangeUser(0,500);
  }
  if( strcmp(var,"hdilMass") == 0 && ( histtype == "ee" || histtype == "mm" ) ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 76  , 0 , 76  , 1.05 * max );
    line.DrawLine( 106 , 0 , 106 , 1.05 * max );
  }

  if( strcmp(var,"htopMass") == 0 && ( histtype == "ee" || histtype == "mm" ) ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    //    line.DrawLine( 76  , 0 , 76  , 1.05 * max );
    //    line.DrawLine( 106 , 0 , 106 , 1.05 * max );
  }

  if( strcmp(var,"hnJet") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 2 , 0 , 2 , 1.05 * max );
  }
  if( strcmp(var,"hnJpt") == 0 ){
    float max = hist->GetBinContent( hist->GetMaximumBin() );
    max += hist->GetBinError( hist->GetMaximumBin() );
    line.DrawLine( 2 , 0 , 2 , 1.05 * max );
  }
  TLatex *t = new TLatex();
  t->SetNDC();

  t->SetTextSize(0.05);
  t->DrawLatex(0.6,0.86,"CMS");
  t->DrawLatex(0.6,0.78,"34.0 pb^{-1} at #sqrt{s} = 7 TeV");
  if( histtype == "all" ) t->DrawLatex(0.6,0.7,"Events with ee/#mu#mu/e#mu");
  if( histtype == "ee"  ) t->DrawLatex(0.6,0.7,"Events with ee");
  if( histtype == "mm"  ) t->DrawLatex(0.6,0.7,"Events with #mu#mu");
  if( histtype == "em"  ) t->DrawLatex(0.6,0.7,"Events with e#mu");
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

  if( color < 0 ){
    hout->SetLineColor(1);
    hout->SetMarkerColor(1);
    hout->SetFillColor(0);
    hout->SetLineStyle(2);
  }else{
    hout->SetLineColor(1);
    hout->SetMarkerColor(color);
    hout->SetFillColor(color);
  }

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

  const unsigned int nSM   = smprefix.size();
  const unsigned int nSUSY = susyprefix.size();

  //TLegend *leg =new TLegend(0.6,0.35,0.8,0.85);
  TLegend *leg =new TLegend(0.1,0.3,0.7,0.8);
  
  TH1F* hdatadummy = new TH1F("hdatadummy","",1,0,1);
  leg->AddEntry(hdatadummy,"data");

  //SM MC
  for(int i = nSM - 1 ; i >= 0 ; i--){
    if( TString( smprefix.at(i) ).Contains("LM") ) continue;
    TH1* hdummy=new TH1("hSMdummy","",1,0,1);
    hdummy->SetLineColor(1);
    hdummy->SetMarkerColor(colors[i]);
    hdummy->SetFillColor(colors[i]);
    leg->AddEntry(hdummy,smleg.at(i),"f");
    
  }

  //LM MC
  for(int i = nSM - 1 ; i >= 0 ; i--){
    if( !TString( smprefix.at(i) ).Contains("LM") ) continue;
    TH1* hdummy=new TH1("hSMdummy","",1,0,1);
    hdummy->SetLineColor(1);
    hdummy->SetLineStyle(2);;
    hdummy->SetMarkerColor(0);
    hdummy->SetFillColor(0);
    leg->AddEntry(hdummy,smleg.at(i),"f");
    
  }

  for(unsigned int i=0;i<nSUSY;i++){
    TH1* hdummy=new TH1("hdummy","",1,0,1);
    hdummy->SetLineColor(i+1);
    hdummy->SetMarkerColor(i+1);
    hdummy->SetLineStyle(2);
    leg->AddEntry(hdummy,susyprefix.at(i),"l");
  }
  leg->SetFillColor(0);
  leg->SetBorderSize(1);

  leg->SetTextSize(0.1);
  return leg;
}

int rebin(char* prefix){
  int r = 1;

  if(strcmp(prefix,"htcmet")==0)    r = 5;
  if(strcmp(prefix,"htcmet_sqrtht")==0)    r = 20;
  if(strcmp(prefix,"hsumJetPt")==0) r = 10;
  if(strcmp(prefix,"hsumJptPt")==0) r = 10;
  if(strcmp(prefix,"hdilMass")==0)  r = 5;

  if(strcmp(prefix,"htopMass")==0)  r = 20;

  if(strcmp(prefix,"hdilPt")==0)    r = 5;
  if(strcmp(prefix,"hminLepPt")==0) r = 4;
  if(strcmp(prefix,"hmaxLepPt")==0) r = 4;
  if(strcmp(prefix,"hminLepEta")==0) r = 6;
  if(strcmp(prefix,"hmaxLepEta")==0) r = 6;
  if(strcmp(prefix,"hptJpt1")==0)   r = 4;
  if(strcmp(prefix,"hptJpt2")==0)   r = 4;
  if(strcmp(prefix,"hptJpt3")==0)   r = 4;
  if(strcmp(prefix,"hptJpt4")==0)   r = 4;
  if(strcmp(prefix,"hptJet1")==0)   r = 4;
  if(strcmp(prefix,"hptJet2")==0)   r = 4;
  if(strcmp(prefix,"hptJet3")==0)   r = 4;
  if(strcmp(prefix,"hptJet4")==0)   r = 4;
  if(strcmp(prefix,"hmeffJPT")==0)  r = 5;
  if(strcmp(prefix,"hmt")==0)       r = 50;
  if(strcmp(prefix,"hmt2core")==0)  r = 5;
  if(strcmp(prefix,"hmt2jcore")==0) r = 20;
  if(strcmp(prefix,"hdphiLep")==0)  r = 5;
  if(strcmp(prefix,"hetaJpt1")==0)  r = 5;
  if(strcmp(prefix,"hetaJpt2")==0)  r = 5;
  if(strcmp(prefix,"hetaJpt3")==0)  r = 5;
  if(strcmp(prefix,"hetaJpt4")==0)  r = 5;
  if(strcmp(prefix,"hetaJet1")==0)  r = 5;
  if(strcmp(prefix,"hetaJet2")==0)  r = 5;
  if(strcmp(prefix,"hetaJet3")==0)  r = 5;
  if(strcmp(prefix,"hetaJet4")==0)  r = 5;
  if(strcmp(prefix,"hptBtagJpt1")==0)   r = 4;
  if(strcmp(prefix,"hptBtagJpt2")==0)   r = 4;
  if(strcmp(prefix,"hptBtagJpt3")==0)   r = 4;
  if(strcmp(prefix,"hptBtagJpt4")==0)   r = 4;
  if(strcmp(prefix,"hetaBtagJpt1")==0)  r = 5;
  if(strcmp(prefix,"hetaBtagJpt2")==0)  r = 5;
  if(strcmp(prefix,"hetaBtagJpt3")==0)  r = 5;
  if(strcmp(prefix,"hetaBtagJpt4")==0)  r = 5;
  if(strcmp(prefix,"hdrLep")==0)        r = 5;
  if(strcmp(prefix,"hdrJ1J2")==0)       r = 5;
  
  return r;

}
THStack* getSMStack(TFile *file, char* varname, char* xtitle,char* histtype){
  
  const unsigned int nSM = smprefix.size();
  
  THStack *stack=new THStack("stack","");
  TH1F* h[nSM];

  char* allj="allj_";
  if(strcmp(varname,"hnJet")==0)    allj="";
  if(strcmp(varname,"hnJpt")==0)    allj="";
  if(strcmp(varname,"hnBtagJpt")==0)allj="";


  
  for(unsigned int iSM=0;iSM<nSM;iSM++){
    if(strcmp(histtype,"all")==0 || strcmp(histtype,"ee")==0 || 
       strcmp(histtype,"mm")==0  || strcmp(histtype,"em")==0){
      
      cout<<"Getting "<<Form("%s_%s_%s%s",smprefix.at(iSM),varname,allj,histtype)<<endl;

      if( strcmp( smprefix.at(iSM) , "DYall" ) == 0 ){
        TH1F *temp1 = getCloneHist((TH1F*)(file->Get(Form("DYee_%s_%s%s",    varname,allj,histtype))), colors[iSM]);
        TH1F *temp2 = getCloneHist((TH1F*)(file->Get(Form("DYmm_%s_%s%s",    varname,allj,histtype))), colors[iSM]);
        TH1F *temp3 = getCloneHist((TH1F*)(file->Get(Form("DYtautau_%s_%s%s",varname,allj,histtype))), colors[iSM]);
        h[iSM] = (TH1F*) temp1->Clone();
        h[iSM]->Add(temp2);
        h[iSM]->Add(temp3);
      }

      else if( strcmp( smprefix.at(iSM) , "VV" ) == 0 ){
        TH1F *temp1 = getCloneHist((TH1F*)(file->Get(Form("ww_%s_%s%s",    varname,allj,histtype))), colors[iSM]);
        TH1F *temp2 = getCloneHist((TH1F*)(file->Get(Form("wz_%s_%s%s",    varname,allj,histtype))), colors[iSM]);
        TH1F *temp3 = getCloneHist((TH1F*)(file->Get(Form("zz_%s_%s%s",    varname,allj,histtype))), colors[iSM]);
        h[iSM] = (TH1F*) temp1->Clone();
        h[iSM]->Add(temp2);
        h[iSM]->Add(temp3);
      }

      else{
        h[iSM]  = getCloneHist((TH1F*)(file->Get(Form("%s_%s_%s%s",smprefix.at(iSM),varname,allj,histtype))), colors[iSM]);
      }
    }
    if(strcmp(histtype,"fs")==0){
      cout << "FS NOT SUPPORTED" << endl;
      exit(0);

      cout<<"Getting FS "<<Form("%s_%s_allj_all",smprefix.at(iSM),varname)<<endl;
      h[iSM]  = getFSHist((TH1F*)(file->Get(Form("%s_%s_%see",smprefix.at(iSM),varname,allj)))->Clone(),
			  (TH1F*)(file->Get(Form("%s_%s_%smm",smprefix.at(iSM),varname,allj)))->Clone(),
			  (TH1F*)(file->Get(Form("%s_%s_%sem",smprefix.at(iSM),varname,allj)))->Clone(),
			  colors[iSM]);

     
    }
    if(strcmp(histtype,"sf")==0){
      cout << "SF NOT SUPPORTED" << endl;
      exit(0);
      
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
  if(strcmp(varname,"hnJet")==0)    allj="";
  if(strcmp(varname,"hnJpt")==0)    allj="";
  if(strcmp(varname,"hnBtagJpt")==0)allj="";

  
  if(strcmp(histtype,"all")==0 || strcmp(histtype,"ee")==0 || 
     strcmp(histtype,"mm")==0  || strcmp(histtype,"em")==0){
    
    cout<<"Getting "<<Form("%s_%s_%s%s",data,varname,allj,histtype)<<endl;
    hist  = getCloneHist((TH1F*)(file->Get(Form("%s_%s_%s%s",data,varname,allj,histtype))), 
                         1);
    
  }
  if(strcmp(histtype,"fs")==0){
    
    cout<<"Getting FS "<<Form("%s_%s_allj_all",data,varname)<<endl;
    hist  = getFSHist((TH1F*)(file->Get(Form("%s_%s_%see",data,varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%smm",data,varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%sem",data,varname,allj)))->Clone(),
                      1);
    
    
  }
  if(strcmp(histtype,"sf")==0){
    
    cout<<"Getting SF "<<Form("%s_%s_allj_all",data,varname)<<endl;
    hist  = getSFHist((TH1F*)(file->Get(Form("%s_%s_%see",data,varname,allj)))->Clone(),
                      (TH1F*)(file->Get(Form("%s_%s_%smm",data,varname,allj)))->Clone(),
                      1);
    
    
    
  }
  
  if(rebin(varname)>1) hist->Rebin(rebin(varname));
  
  //if(strcmp(histtype,"all")==0) hist->SetTitle(Form("%s (ee + #mu#mu + e#mu)",xtitle));
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetYaxis()->SetTitle("Events   ");
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetTitleOffset(0.8);
  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetTitleOffset(0.8);
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
  if(strcmp(varname,"hnJpt")==0)allj="";
  
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
    hsusy[iSUSY]->SetLineStyle(2);
    hsusy[iSUSY]->Draw("samehist");
    
    //if( hsusy[iSUSY]->GetMinimum() < histmin ) histmin = hsusy[iSUSY]->GetMinimum();
    //if( hsusy[iSUSY]->GetMaximum() > histmax ) histmax = hsusy[iSUSY]->GetMaximum();
  }
}

