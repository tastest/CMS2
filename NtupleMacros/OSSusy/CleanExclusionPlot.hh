#ifndef MAKE_EXCLUSIONPLOT_HH
#define MAKE_EXCLUSIONPLOT_HH

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraphErrors.h"


#include "vector.h"

bool removePoints = false;

void CleanExclusionPlot( char* filename );
  

void CommandMSUGRA(TString plotName , char* filename);

TGraph* set_sneutrino_d0_1(){
  double sn_m0[9]={0,0,55,90,100,110,100,55,0};
  double sn_m12[9]={0,140,220,240,240,230,210,150,0};

  TGraph* sn_d0_gr = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(1001);


  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(){
  double sn_m0[5]={0,50,105,190,0};
  double sn_m12[5]={0,140,205,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(5,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(1001);


  return sn_d0_gr_2;
}



TGraph* set_lep_ch(){

  double ch_m0[5];
  double ch_m12[5];

 
    ch_m0[0] = 0;
    ch_m0[1] = 0;
    ch_m0[2] = 600;
    ch_m0[3] = 1000;
    ch_m0[4] = 1000;

    ch_m12[0] = 0;
    ch_m12[1] = 130;
    ch_m12[2] = 120;
    ch_m12[3] = 113;
    ch_m12[4] = 0;

  
  TGraph* ch_gr = new TGraph(5,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}


TGraph* set_lep_sl(){



  double sl_m0[] = {0,0,30,50,60,75,80,90,100};
  double sl_m12[] = {0,245,240,220,200,150,100,50,0}; 
  
  TGraph* lep_sl = new TGraph(9,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetFillStyle(1001);
  
  return lep_sl;
}

TGraph* set_tev_sg_cdf(){

  double sg_m0[] = {0,50,100,150,200,250,300,350,400,450,500,550,600,600};
  double sg_m12[] = {0,170,160,155,150,122,116,112,110,106,105,100,98,0};
  TGraph* sg_gr = new TGraph(14,sg_m0,sg_m12);

  sg_gr->SetFillColor(2);
  sg_gr->SetFillStyle(1001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(){
  double sgd_m0[] = {0,50,100,150,200,250,300,350,400,450,500,550,600,600};
  double sgd_m12[] = {0,173,170,168,160,150,140,130,125,120,120,120,120,0};
  TGraph* sgd_gr = new TGraph(14,sgd_m0,sgd_m12);

  sgd_gr->SetFillColor(41);
  sgd_gr->SetFillStyle(1001);

  return sgd_gr;

}

TGraph* set_tev_tlp_cdf(){
  double tlp1_m0[] = {0,20,40,60,70,80,90,80,70,60};
  double tlp1_m12[] = {170,185,200,215,220,215,210,190,175,160};
  TGraph* tlp1_gr = new TGraph(10,tlp1_m0,tlp1_m12);

  tlp1_gr->SetFillColor(4);
  tlp1_gr->SetFillStyle(1001);

  return tlp1_gr;
}

TGraph* set_tev_tlp_d0(){
  double tlp2_m0[] = {70,80,90,100,105,110,120,130,140};
  double tlp2_m12[] = {160,172,184,196,205,195,185,173,160};
  TGraph* tlp2_gr = new TGraph(9,tlp2_m0,tlp2_m12);

  tlp2_gr->SetFillColor(4);
  tlp2_gr->SetFillStyle(1001); 

  return tlp2_gr;

}

TGraph* set_tev_stau(){
  double st_m0[] = {0,30,200,0,0};
  double st_m12[] = {230,240,1000,1000,230};
  TGraph* st_gr = new TGraph(5,st_m0,st_m12);

  st_gr->SetFillColor(40);
  st_gr->SetFillStyle(1001);


  return st_gr;

}

TF1* constant_squark(int i){
//---lines of constant gluino/squark
  double coef1 = 0.35;
  double coef2[] = {5,5,4.6,4.1};

  char hname[200];

  sprintf(hname,"lnsq_%i",i); 

  
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]-[2])",0,1000);
  lnsq->SetParameter(0,(500+150*(i-1))*(500+150*(i-1))/coef2[i]);
  lnsq->SetParameter(1,1./coef2[i]);
  lnsq->SetParameter(2,-coef1*91*91*99./101/coef2[i]);//--tanbeta=10 --> cos2beta = -99/101
  lnsq->SetLineWidth(1);


  lnsq->SetLineColor(kBlack);

  return lnsq;
}


TF1* constant_gluino(int i){
//---lines of constant gluino/squark
  double coef1 = 0.35;
  double coef2[] = {5,5,4.6,4.1};

  char hname[200];

  sprintf(hname,"lngl_%i",i); 
    
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,1000);
  lngl->SetParameter(0,(500+150.*(i-1))/2.4);
  lngl->SetParameter(1,-40./1400);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kBlack);

  return lngl;
}


TLatex* constant_squark_text(Int_t it,TF1& lnsq){
  char legnm[200];

  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",500+150*(it-1));
  TLatex *t3;
 
  if( it == 1 ){
    t3 = new TLatex(220+10*(it-1),lnsq.Eval(220+10*(it-1))+10,legnm);
    t3->SetTextAngle(-12);
  }
  else{
    t3 = new TLatex(180+10*(it-1),lnsq.Eval(180+10*(it-1))-12,legnm);
    t3->SetTextAngle(-10);
  }
  t3->SetTextSize(0.04);

  
  return t3;
}

TLatex* constant_gluino_text(Int_t it,TF1& lngl){
  char legnm[200];

  sprintf(legnm,"#font[92]{#tilde{g}(%i)GeV}",500+150*(it-1));
  TLatex* t4 = new TLatex(423,lngl.Eval(400),legnm);
  t4->SetTextSize(0.04);
  t4->SetTextAlign(13);

  return t4;
}





TLegend* makeStauLegend(Double_t txtsz){
  TLegend* legst = new TLegend(0.22,0.83,0.2,0.85);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);

  return legst;


}

TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz){
  //TLegend* legexp = new TLegend(0.6796657,0.6718346,0.9791086,0.9509044,NULL,"brNDC");
  //TLegend* legexp = new TLegend(0.61,0.65,0.93,0.92,NULL,"brNDC");
  //TLegend* myleg = new TLegend(0.35,0.81,0.6,0.91,NULL,"brNDC");
  TLegend* legexp = new TLegend(0.61,0.67,0.9,0.92,NULL,"brNDC");

  legexp->SetFillColor(0);
  legexp->AddEntry(&sg_gr,"CDF   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, tan#beta = 5,  2 fb^{-1}","f");  
  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, tan#beta = 3,  2.1 fb^{-1}","f");  
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  //legexp->AddEntry(&tev_sn,"D0  #tilde{#nu}","f");  
  legexp->AddEntry(&tev_sn,"D0  #tilde{#chi}_{1}^{#pm}, #tilde{#chi}_{2}^{0}","f");  
  legexp->SetShadowColor(0);
  txtsz-=0.015;
  legexp->SetTextSize(txtsz);

  return legexp;

}

/*
TGraphErrors* getLO_signalCont(){



  Int_t nl = 10;
  Double_t xl[nl];
  Double_t yl[nl];
  Double_t exl[nl];
  Double_t eyl[nl];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 0;
  yl[0] = 270;
  xl[1] = 100;
  yl[1] = 263;
  xl[2] = 200;
  yl[2] = 254;
  xl[3] = 250;
  yl[3] = 246;
  xl[4] = 300;
  yl[4] = 212;
  xl[5] = 340;
  yl[5] = 176;
  xl[6] = 400;
  yl[6] = 143;
  xl[7] = 430;
  yl[7] = 134;
  xl[8] = 450;
  yl[8] = 130;
  xl[9] = 490;
  yl[9] = 50;
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;
}
*/

/*
TGraphErrors* getExpected_NLOunc(){

 Int_t nl = 11;
  Double_t xl[nl];
  Double_t yl[nl];
  Double_t exl[nl];
  Double_t eyl[nl];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 35;
  yl[0] = 285;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 276;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 230;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 155;
  xl[10] = 550;
  yl[10] = 50;
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kWhite);
  //gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  //  TSpline3 *s = new TSpline3("grs",gr1);
  // s->SetLineColor(kBlue);
  // s->SetLineStyle(2);
  // s->SetLineWidth(3);
  

  return gr1;





}
*/

/*
TGraphErrors* getObserved_NLOunc(){

  Int_t nl = 11;
  Double_t xl[nl];
  Double_t yl[nl];
  Double_t exl[nl];
  Double_t eyl[nl];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 35;
  yl[0] = 275;
  xl[1] = 100;
  yl[1] = 270;
  xl[2] = 150;
  yl[2] = 268;
  xl[3] = 200;
  yl[3] = 265;
  xl[4] = 250;
  yl[4] = 255;
  xl[5] = 300;
  yl[5] = 236;
  xl[6] = 350;
  yl[6] = 200;
  xl[7] = 400;
  yl[7] = 172;
  xl[8] = 450;
  yl[8] = 158;
  xl[9] = 500;
  yl[9] = 132;
  xl[10] = 510;
  yl[10] = 120;
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kWhite);
  // gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;

}
*/

TGraphErrors* getObserved_NLOunc( char* filename ){

  cout << "Retrieving observed NLO exclusion curve from " << filename << endl;
  TFile *f = TFile::Open(filename);
  TGraphErrors *gre = (TGraphErrors*) f->Get("limitgraph_NLO_obs");
  gre->SetMarkerColor(kWhite);

  if( removePoints ){
    Double_t xp;
    Double_t yp;
    
    for( int i = 12 ; i < 50 ; i++ ){
      gre->GetPoint(i,xp,yp);
      gre->SetPoint(i,xp,100);
    }
  }

  return gre;
}

TGraphErrors* getExpected_NLOunc( char* filename ){

  //cout get observed contour, reduce reach by 10 GeV
  cout << "Retrieving expected NLO exclusion curve from " << filename << endl;
  TFile *f = TFile::Open(filename);
  TGraphErrors *gre = (TGraphErrors*) f->Get("limitgraph_NLO_obs");
  gre->SetMarkerColor(kWhite);

  Double_t xp;
  Double_t yp;
    
  for( int i = 0 ; i < 50 ; i++ ){
    gre->GetPoint(i,xp,yp);
    gre->SetPoint(i,xp,yp-10);
  }

  return gre;

  /*
  cout << "Retrieving expected NLO exclusion curve from " << filename << endl;
  TFile *f = TFile::Open(filename);
  TGraphErrors *gre = (TGraphErrors*) f->Get("limitgraph_NLO_exp");
  gre->SetMarkerColor(kWhite);

  if( removePoints ){
    
    Double_t xp;
    Double_t yp;
    
    for( int i = 12 ; i < 50 ; i++ ){
      gre->GetPoint(i,xp,yp);
      gre->SetPoint(i,xp,100);
    }
  }

  return gre;
  */

}

TGraphErrors* getLO_signalCont( char* filename ){

  cout << "Retrieving observed LO exclusion curve from " << filename << endl;
  TFile *f = TFile::Open(filename);
  TGraphErrors *gre = (TGraphErrors*) f->Get("limitgraph_LO_obs");
  gre->SetMarkerColor(kGreen+2);
  gre->SetMarkerStyle(21);

  if( removePoints ){

    Double_t xp;
    Double_t yp;
    
    for( int i = 12 ; i < 50 ; i++ ){
      gre->GetPoint(i,xp,yp);
      gre->SetPoint(i,xp,100);
    }
  }

  return gre;

}


#endif
