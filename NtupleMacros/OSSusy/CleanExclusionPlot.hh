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
void CleanExclusionPlot();
  

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

  
//   cout << "Retrieving expected NLO exclusion curve from " << filename << endl;
//   TFile *f = TFile::Open(filename);
//   TGraphErrors *gre = (TGraphErrors*) f->Get("limitgraph_NLO_exp");
//   gre->SetMarkerColor(kWhite);

//   if( removePoints ){
    
//     Double_t xp;
//     Double_t yp;
    
//     for( int i = 12 ; i < 50 ; i++ ){
//       gre->GetPoint(i,xp,yp);
//       gre->SetPoint(i,xp,100);
//     }
//   }

//   return gre;
  

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


/*
TGraphErrors* getObserved_NLOunc(char* filename){

  Int_t n=67;
  Double_t x[n],y[n];
  Double_t xerr[n],yerr[n];
  for( unsigned int ierr = 0 ; ierr < n ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
  }

  Int_t i = -1;

  x[++i]=10;
  y[i]=255.;

  x[++i]=20;
  y[i]=255.;

  x[++i]=30.;
  y[i]=265.;

  x[++i]=40.;
  y[i]=265.;

  x[++i]=70.;
  y[i]=265.;
  
  x[++i]=80.;
  y[i]=275.;

  x[++i]=100.;
  y[i]=275.;

  x[++i]=120.;
  y[i]=275.;

  x[++i]=130.;
  y[i]=275.;

  x[++i]=140.;
  y[i]=275.;

  x[++i]=145.;
  y[i]=265.;

  x[++i]=150.;
  y[i]=255.;

  x[++i]=160.;
  y[i]=235.;

  x[++i]=170.;
  y[i]=225.;

  x[++i]=180.;
  y[i]=215.;

  x[++i]=190.;
  y[i]=208.;

  x[++i]=200.;
  y[i]=205.;

  x[++i]=210.;
  y[i]=195.;

  x[++i]=220.;
  y[i]=170.;

  x[++i]=230.;
  y[i]=185.;

  x[++i]=240.;
  y[i]=175.;

  x[++i]=250.;
  y[i]=175.;

  x[++i]=260.;
  y[i]=155.;

  x[++i]=270.;
  y[i]=165.;

  x[++i]=280.;
  y[i]=162.;

  x[++i]=290.;
  y[i]=165.;

  x[++i]=300.;
  y[i]=165.;

  x[++i]=310.;
  y[i]=160.;

  x[++i]=300.;
  y[i]=155.;

  x[++i]=290.;
  y[i]=155.;

  x[++i]=280.;
  y[i]=160.;

  x[++i]=275.;
  y[i]=150.;

  x[++i]=280.;
  y[i]=135.;

  x[++i]=290.;
  y[i]=125.;;

  x[++i]=300.;
  y[i]=125.;;

  x[++i]=310.;
  y[i]=120.;

  x[++i]=320.;
  y[i]=125.;

  x[++i]=330.;
  y[i]=125.;

  x[++i]=340.;
  y[i]=135.;

  x[++i]=350.;
  y[i]=135.;

  x[++i]=360.;
  y[i]=135.;

  x[++i]=370.;
  y[i]=135.;

  x[++i]=380.;
  y[i]=135.;

  x[++i]=390.;
  y[i]=135.;

  x[++i]=400.;
  y[i]=145.;

  x[++i]=420.;
  y[i]=145.;

  x[++i]=430.;
  y[i]=145.;

  x[++i]=440.;
  y[i]=145.;

  x[++i]=450.;
  y[i]=145.;

  x[++i]=460.;
  y[i]=145.;

  x[++i]=465.;
  y[i]=140.;

  x[++i]=460.;
  y[i]=135.;

  x[++i]=450.;
  y[i]=135.;

  x[++i]=440.;
  y[i]=135.;

  x[++i]=430.;
  y[i]=135.;

  x[++i]=420.;
  y[i]=135.;

  x[++i]=410.;
  y[i]=125.;

  x[++i]=400.;
  y[i]=125.;

  x[++i]=390.;
  y[i]=125.;

  x[++i]=380.;
  y[i]=125.;

  x[++i]=370.;
  y[i]=125.;

  x[++i]=360.;
  y[i]=115.;

  x[++i]=350.;
  y[i]=115.;

  x[++i]=340.;
  y[i]=115.;

  x[++i]=330.;
  y[i]=115.;

  x[++i]=320.;
  y[i]=115.;

  x[++i]=310.;
  y[i]=115.;

  x[++i]=300.;
  y[i]=100.;

  TGraphErrors* gr = new TGraphErrors(n,x,y,xerr,yerr);
  gr->SetMarkerColor(kWhite);
  return gr;

}



TGraphErrors* getExpected_NLOunc(char* filename){


  Int_t nexp=71;
  Double_t xexp[nexp],yexp[nexp];
  Double_t xerr[nexp],yerr[nexp];
  for( unsigned int ierr = 0 ; ierr < nexp ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
  }

  Int_t i = -1;

  xexp[++i]=10;
  yexp[i]=255.;

  xexp[++i]=20;
  yexp[i]=245.;

  xexp[++i]=30.;
  yexp[i]=255.;

  xexp[++i]=40.;
  yexp[i]=255.;

  xexp[++i]=50.;
  yexp[i]=265.;

  xexp[++i]=60.;
  yexp[i]=265.;

  xexp[++i]=70.;
  yexp[i]=265.;
  
  xexp[++i]=80.;
  yexp[i]=273.;

  xexp[++i]=90.;
  yexp[i]=265.;

  xexp[++i]=100.;
  yexp[i]=265.;

  xexp[++i]=110.;
  yexp[i]=265.;

  xexp[++i]=120.;
  yexp[i]=275.;

  xexp[++i]=130.;
  yexp[i]=275.;

  xexp[++i]=140.;
  yexp[i]=235.;

  xexp[++i]=150.;
  yexp[i]=235.;

  xexp[++i]=160.;
  yexp[i]=235.;

  xexp[++i]=170.;
  yexp[i]=225.;

  xexp[++i]=180.;
  yexp[i]=205.;

  xexp[++i]=190.;
  yexp[i]=205.;

  xexp[++i]=200.;
  yexp[i]=205.;

  xexp[++i]=210.;
  yexp[i]=195.;

  xexp[++i]=220.;
  yexp[i]=155.;

  xexp[++i]=230.;
  yexp[i]=155.;

  xexp[++i]=240.;
  yexp[i]=165.;

  xexp[++i]=250.;
  yexp[i]=155.;

  xexp[++i]=260.;
  yexp[i]=155.;

  xexp[++i]=270.;
  yexp[i]=155.;

  //
  xexp[++i]=280.;
  yexp[i]=160.;

  xexp[++i]=290.;
  yexp[i]=165.;

  xexp[++i]=300.;
  yexp[i]=165.;

  xexp[++i]=310.;
  yexp[i]=160.;

  xexp[++i]=300.;
  yexp[i]=155.;

  xexp[++i]=290.;
  yexp[i]=159.;

  xexp[++i]=280.;
  yexp[i]=160.;

  xexp[++i]=275.;
  yexp[i]=150.;

  xexp[++i]=255.;
  yexp[i]=135.;

  xexp[++i]=290.;
  yexp[i]=125.;;

  xexp[++i]=300.;
  yexp[i]=125.;;

  xexp[++i]=310.;
  yexp[i]=120.;

  xexp[++i]=320.;
  yexp[i]=125.;

  xexp[++i]=330.;
  yexp[i]=125.;

  xexp[++i]=340.;
  yexp[i]=125.;

  xexp[++i]=350.;
  yexp[i]=125.;

  xexp[++i]=360.;
  yexp[i]=135.;

  xexp[++i]=370.;
  yexp[i]=135.;

  xexp[++i]=380.;
  yexp[i]=135.;

  xexp[++i]=390.;
  yexp[i]=135.;

  xexp[++i]=400.;
  yexp[i]=145.;

  xexp[++i]=420.;
  yexp[i]=145.;

  xexp[++i]=430.;
  yexp[i]=145.;

  xexp[++i]=440.;
  yexp[i]=145.;

  xexp[++i]=450.;
  yexp[i]=145.;

  xexp[++i]=460.;
  yexp[i]=145.;

  xexp[++i]=460.;
  yexp[i]=140.;

  xexp[++i]=460.;
  yexp[i]=135.;

  xexp[++i]=450.;
  yexp[i]=135.;

  xexp[++i]=440.;
  yexp[i]=135.;

  xexp[++i]=430.;
  yexp[i]=135.;

  xexp[++i]=420.;
  yexp[i]=135.;

  xexp[++i]=410.;
  yexp[i]=130.;

  xexp[++i]=400.;
  yexp[i]=130.;

  xexp[++i]=390.;
  yexp[i]=125.;

  xexp[++i]=380.;
  yexp[i]=125.;

  xexp[++i]=370.;
  yexp[i]=125.;

  xexp[++i]=360.;
  yexp[i]=115.;

  xexp[++i]=350.;
  yexp[i]=115.;

  xexp[++i]=340.;
  yexp[i]=115.;

  xexp[++i]=330.;
  yexp[i]=115.;

  xexp[++i]=320.;
  yexp[i]=115.;

  xexp[++i]=310.;
  yexp[i]=115.;

  xexp[++i]=300.;
  yexp[i]=100.;

  //cout << i<< endl;

  TGraphErrors* grexp = new TGraphErrors(nexp,xexp,yexp,xerr,yerr);
  grexp->SetMarkerColor(kWhite);
  return grexp;

}

TGraphErrors* getLO_signalCont(char* filename){
  //TGraphErrors* getLO(char* filename){
  
  TFile *f = TFile::Open(filename);

  Int_t nLO=39;
  Double_t xLO[nLO],yLO[nLO];
  Double_t xerr[nLO],yerr[nLO];
  for( unsigned int ierr = 0 ; ierr < 39 ; ++ierr ){
    xerr[ierr] = 0.;
    yerr[ierr] = 0.;
    xLO[ierr] = 0.;
    yLO[ierr] = 0.;
  }

  Int_t i = -1;

  xLO[++i]=10;
  yLO[i]=245.;

  xLO[++i]=20;
  yLO[i]=245.;

  xLO[++i]=30.;
  yLO[i]=245.;

  xLO[++i]=40.;
  yLO[i]=255.;

  xLO[++i]=50.;
  yLO[i]=255.;

  xLO[++i]=60.;
  yLO[i]=255.;
  
  xLO[++i]=70.;
  yLO[i]=263.;
  
  xLO[++i]=80.;
  yLO[i]=265.;

  xLO[++i]=100.;
  yLO[i]=265.;

  xLO[++i]=120.;
  yLO[i]=265.;

  xLO[++i]=130.;
  yLO[i]=245.;

  xLO[++i]=140.;
  yLO[i]=235.;

  xLO[++i]=150.;
  yLO[i]=225.;

  xLO[++i]=160.;
  yLO[i]=215.;

  xLO[++i]=170.;
  yLO[i]=205.;

  xLO[++i]=180.;
  yLO[i]=195.;

  xLO[++i]=190.;
  yLO[i]=195.;

  xLO[++i]=200.;
  yLO[i]=180.;

  xLO[++i]=210.;
  yLO[i]=175.;

  xLO[++i]=220.;
  yLO[i]=155.;

  xLO[++i]=230.;
  yLO[i]=150.;

  xLO[++i]=240.;
  yLO[i]=145.;

  xLO[++i]=250.;
  yLO[i]=155.;

  xLO[++i]=255.;
  yLO[i]=140.;

  xLO[++i]=255.;
  yLO[i]=135.;

  xLO[++i]=240.;
  yLO[i]=125.;

  xLO[++i]=230.;
  yLO[i]=115.;

  xLO[++i]=240.;
  yLO[i]=115.;

  xLO[++i]=250.;
  yLO[i]=115.;

  xLO[++i]=260.;
  yLO[i]=115.;

  xLO[++i]=270.;
  yLO[i]=125.;

  xLO[++i]=280.;
  yLO[i]=125.;

  xLO[++i]=290.;
  yLO[i]=125.;;

  xLO[++i]=300.;
  yLO[i]=125.;;

  xLO[++i]=310.;
  yLO[i]=120.;

  xLO[++i]=320.;
  yLO[i]=125.;

  xLO[++i]=320.;
  yLO[i]=120.;

  xLO[++i]=310.;
  yLO[i]=110.;

  xLO[++i]=305.;
  yLO[i]=100.;
  

  //cout << i<< endl;
  TGraphErrors* grLO  = new TGraphErrors(nLO,xLO, yLO,xerr,yerr);
  grLO->SetMarkerColor(kGreen+2);
  grLO->SetMarkerStyle(21);
  return grLO;
  
//     Double_t xxLO[9],yyLO[9];
//     xxLO[0]=345.;
//     yyLO[0]=120.;

//     xxLO[1]=350.;
//     yyLO[1]=125.;

//     xxLO[2]=360.;
//     yyLO[2]=135.;

//     xxLO[3]=370.;
//     yyLO[3]=135.;

//     xxLO[4]=375.;
//     yyLO[4]=130.;

//     xxLO[5]=370.;
//     yyLO[5]=125.;

//     xxLO[6]=360.;
//     yyLO[6]=115.;

//     xxLO[7]=350.;
//     yyLO[7]=115.;

//     xxLO[8]=345.;
//     yyLO[8]=120.;
 
//     grLO2 = new TGraph(9,  xxLO,yyLO);

//     // Draw the curves
//     grLO.Draw("C");
//     grLO2.Draw("C");

}
*/


#endif
