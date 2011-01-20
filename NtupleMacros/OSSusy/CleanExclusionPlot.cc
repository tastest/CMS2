#include <sstream>
#include <string>
#include <iomanip>
#include "CleanExclusionPlot.hh"
 
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "vector.h"
#include "TMath.h"

void CleanExclusionPlot( char* filename ){
  gStyle->SetPalette(1);

  CommandMSUGRA("35pb_expected_11.root" , filename );
  

}

void CommandMSUGRA(TString plotName_ , char* filename ){

  // Output file
  TFile* output = new TFile( plotName_, "RECREATE" );
  if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }

  //set old exclusion Limits
  TGraph* LEP_ch = set_lep_ch();
  TGraph* LEP_sl = set_lep_sl();//slepton curve
  TGraph* TEV_sg_cdf = set_tev_sg_cdf();//squark gluino cdf
  TGraph* TEV_sg_d0 = set_tev_sg_d0();//squark gluino d0
  TGraph* TEV_tlp_cdf = set_tev_tlp_cdf();//trilepton cdf
  TGraph* TEV_tlp_d0 = set_tev_tlp_d0();//trilepton d0
  TGraph* stau = set_tev_stau();//stau 

  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1();
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2();

  //constant ssqquark and gluino lines
  TF1* lnsq[4];
  TF1* lngl[4];

  TLatex* sq_text[4];
  TLatex* gl_text[4];

  for(int i = 0; i < 4; i++){
    lnsq[i] = constant_squark(i);
    sq_text[i] = constant_squark_text(i,*lnsq[i]);
    lngl[i] = constant_gluino(i);
    gl_text[i] = constant_gluino_text(i,*lngl[i]);
  }

  //Legends
  TLegend* legst = makeStauLegend(0.05);
  TLegend* legexp = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.04);

  //set Addition text in plot
  TLatex* mytext = new TLatex(135.,458.,"CMS preliminary 2010, L_{int} = 34 pb^{-1}, #sqrt{s} = 7 TeV");
  mytext->SetTextSize(0.03);

  TLatex* mytext_two = new TLatex(80,355,"tan#beta = 3, A_{0} = 0, sign(#mu)>0");
  //TLatex* mytext_two = new TLatex(80,355,"tan#beta = 10, A_{0} = 0, sign(#mu)>0");
  mytext_two->SetTextSize(0.03);

  //make Canvas
  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,1900,1000);
  output->cd();

 
  //the exclusion plots  
  TGraphErrors* First  = getObserved_NLOunc(  filename );
  TGraphErrors* Second = getExpected_NLOunc( filename );//getLO_jetMultis();
  TGraphErrors* Third  = getLO_signalCont(    filename );

  TH2F* hdummy = new TH2F("hdummy","dummy hist",100,0,500,100,80,450);

  First->GetXaxis()->SetRangeUser(0,500);
  First->GetYaxis()->SetRangeUser(80,450);
  First->GetXaxis()->SetTitle("m_{0} (GeV)");
  First->GetYaxis()->SetTitle("m_{1/2} (GeV)");

  TSpline3 *sFirst = new TSpline3("sFirst",First);
  sFirst->SetLineColor(kRed);
  sFirst->SetLineWidth(3);
  First->SetLineColor(kRed);
  First->SetLineWidth(3);

  TSpline3 *sSecond = new TSpline3("sSecond",Second);
  sSecond->SetLineColor(kBlue);
  sSecond->SetLineStyle(2);
  sSecond->SetLineWidth(3);
  Second->SetLineColor(kBlue);
  Second->SetLineStyle(2);
  Second->SetLineWidth(3);
  
  TSpline3 *sThird = new TSpline3("sThird",Third);
  sThird->SetLineColor(kGreen+2);
  sThird->SetLineStyle(4);
  sThird->SetLineWidth(3);
  Third->SetLineColor(kGreen+2);
  Third->SetLineStyle(4);
  Third->SetLineWidth(3);
  

  //Set Legend for Exclusion lines
  //TLegend* myleg = new TLegend(0.3966245,0.8271605,0.6666667,0.9485597,NULL,"brNDC");
  TLegend* myleg = new TLegend(0.35,0.81,0.6,0.92,NULL,"brNDC");
  myleg->SetFillColor(0);
  myleg->SetShadowColor(0);
  myleg->SetTextSize(0.03);
  myleg->AddEntry(sSecond,"NLO Expected Limit","l");
  myleg->AddEntry(sFirst,"NLO Observed Limit","l"); 
  myleg->AddEntry(sThird,"LO Observed Limit","l");

  //Now start drawing-------------------------------------------------
  First->Draw("AP");//graph with white dots just for setting the axis right
  
  bool smooth = false;

  if( smooth ){
    //the exclusion lines (smoothed)
    sFirst->Draw("same");
    sSecond->Draw("same");
    sThird->Draw("same");
  }else{
    //the exclusion lines (un-smoothed)
    First->Draw("samec");
    Second->Draw("samec");
    Third->Draw("samec");
  }

  //the squark and gluino lines plus text
  for (int it=1;it<4;it++) {   
    lngl[it]->Draw("same");   
    lnsq[it]->Draw("same");
    sq_text[it]->Draw();
    gl_text[it]->Draw();
  }

  //exclusion plots of other experiments
  TEV_sn_d0_1->Draw("fsame");
  TEV_sn_d0_2->Draw("fsame");
  TEV_sg_d0->Draw("fsame");
  TEV_sg_cdf->Draw("fsame");
  LEP_ch->Draw("fsame");
  LEP_sl->Draw("fsame");

  //stau is lightest stable particle region
  stau->Draw("fsame");

  //exclusion plot legend
  legexp->Draw();
  legst->Draw();
  myleg->Draw();

  //draw the additional text in histogram
  mytext->Draw("same");
  mytext_two->Draw("same");

  //LM0, LM1
//   float m0[2]  = {200,60};
//   float m12[2] = {160,250};
//   TGraph *gr = new TGraph(2,m0,m12);
//   gr->Draw("sameP");

  hdummy->Draw("axissame");
  //First->Draw("sameP");
  
  //write out canvas
  cvsSys->Write("exclusion.pdf");
  

  output->Write();
  output->Close();
  delete output; 
}
