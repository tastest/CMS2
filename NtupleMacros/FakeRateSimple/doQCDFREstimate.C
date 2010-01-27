#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include <iostream>
#include "histools.h"
#include "QCDFRestimator.h"
#endif __CINT__


void TestPrediction() {
  
//  gROOT->ProcessLine(".L setup.C");
//  gROOT->ProcessLine("setup()");
//  gROOT->ProcessLine(Form(".x setup.C(%d)", 1));
  hist::deleteHistos();
  hist::loadHist("QCDFRplots.root");
   
  TH1F *QCD_FRpt[2];
  TH1F *QCD_FReta[2];
  
  TH1F *WJets_FOpt[2];
  TH1F *WJets_FOeta[2];

  TH1F *WJets_numpt[2];
  TH1F *WJets_numeta[2];

  TH3F *WJets_nJets3D[2];

  TH1F *WJets_actualnJets[2];
  TH1F *WJets_predictednJets[2];

  TH1F *WJets_predictedpt[2];
  TH1F *WJets_predictedeta[2];
  
  
  char *flavor[2] = {"el", "mu"};
  for(unsigned int i = 0; i < 2; i++) {
    QCD_FRpt[i] = (TH1F*)gDirectory->Get(Form("QCD_FRpt_%s", flavor[i]));
    QCD_FReta[i] = (TH1F*)gDirectory->Get(Form("QCD_FReta_%s", flavor[i]));
    
    WJets_FOpt[i] = (TH1F*)gDirectory->Get(Form("WJets_FOpt_%s", flavor[i]));
    WJets_FOeta[i] = (TH1F*)gDirectory->Get(Form("WJets_FOeta_%s", flavor[i]));

    WJets_numpt[i] = (TH1F*)gDirectory->Get(Form("WJets_numpt_%s", flavor[i]));
    WJets_numeta[i] = (TH1F*)gDirectory->Get(Form("WJets_numeta_%s", flavor[i]));

    WJets_nJets3D[i] = (TH3F*)gDirectory->Get(Form("WJets_nJets3D_%s", flavor[i]));

    WJets_actualnJets[i] = (TH1F*)gDirectory->Get(Form("WJets_actualnJets_%s", flavor[i]));
    WJets_predictednJets[i] = (TH1F*)gDirectory->Get(Form("WJets_predictednJets_%s", flavor[i]));
  }      
  
  Float_t pt[4] = {10,20,60,150};
  Float_t eta[3] = {0, 1.479, 2.4};
  WJets_predictedpt[0] = new TH1F(Form("WJets_predictedpt_%s", flavor[0]), 
				  "Predicted Pt Distribution in WJets", 3, pt);
  WJets_predictedeta[0] = new TH1F(Form("WJets_predictedeta_%s", flavor[0]), 
				   "Predicted Eta Distribution in WJets", 2, eta);
  
  WJets_predictedpt[1] = new TH1F(Form("WJets_predictedpt_%s", flavor[1]), 
				  "Predicted Pt Distribution in WJets", 3, pt);
  WJets_predictedeta[1] = new TH1F(Form("WJets_predictedeta_%s", flavor[1]), 
				   "Predicted Eta Distribution in WJets", 2, eta);
  
       
    
  for(unsigned int i = 0; i < 2 ; i++) {

    for(int j = 1; j < QCD_FRpt[i]->GetNbinsX() + 1; j++) {
      Float_t qcd       = QCD_FRpt[i]->GetBinContent(j);
      Float_t qcd_err   = QCD_FRpt[i]->GetBinError(j);
      Float_t wjets     = WJets_FOpt[i]->GetBinContent(j);
      Float_t wjets_err = WJets_FOpt[i]->GetBinError(j);
      Float_t temp = qcd*wjets;
      float err2 = pow(qcd*wjets_err,2) + pow(qcd_err*wjets,2);
      WJets_predictedpt[i]->SetBinContent(j, temp);
      WJets_predictedpt[i]->SetBinError(j, sqrt(err2));
    }
    for(int j = 1; j < QCD_FReta[i]->GetNbinsX() + 1; j++) {
      Float_t qcd       = QCD_FReta[i]->GetBinContent(j);
      Float_t qcd_err   = QCD_FReta[i]->GetBinError(j);
      Float_t wjets     = WJets_FOeta[i]->GetBinContent(j);
      Float_t wjets_err = WJets_FOeta[i]->GetBinError(j);
      Float_t temp = qcd*wjets;
      float err2 = pow(qcd*wjets_err,2) + pow(qcd_err*wjets,2);
      WJets_predictedeta[i]->SetBinContent(j, temp);
      WJets_predictedeta[i]->SetBinError(j, sqrt(err2));
    }
    
  }//flavor loop
  hist::color("WJets_actual", kRed);
  hist::color("WJets_num", kRed);

    
  TCanvas *c1 = new TCanvas();
  WJets_numeta[0]->Draw("e");
  WJets_predictedeta[0]->Draw("samese");
  hist::setrangey(c1);
  c1->SaveAs("elPredictedeta.png");
  
  TCanvas *c2 = new TCanvas();
  WJets_numpt[0]->Draw("e");
  WJets_predictedpt[0]->Draw("sames");
  hist::setrangey(c2);
  c2->SaveAs("elPredictedpt.png");


  TCanvas *c21 = new TCanvas();
  WJets_numeta[1]->Draw("e");
  WJets_predictedeta[1]->Draw("samese");
  hist::setrangey(c21);
  c21->SaveAs("muPredictedeta.png");

  TCanvas *c22 = new TCanvas();
  WJets_numpt[1]->Draw("e");
  WJets_predictedpt[1]->Draw("sames");
  hist::setrangey(c22);
  c22->SaveAs("muPredictedpt.png");

  TCanvas *c3 = new TCanvas();
  ((TH1F*)gDirectory->Get("WJets_actualnJets_el"))->SetMarkerStyle(3);
  (TH1F*)gDirectory->Get("WJets_actualnJets_el")->Draw("e");
  ((TH1F*)gDirectory->Get("WJets_predictednJets_el"))->SetMarkerStyle(25);
  (TH1F*)gDirectory->Get("WJets_predictednJets_el")->Draw("samese");
  hist::setrangey(c3);
  c3->SaveAs("elnJetsPredicted.png");

  TCanvas *c4 = new TCanvas();
  TH2F *QCDFRptvseta_el = (TH2F*)gDirectory->Get("QCD_FRptvseta_el");
  QCDFRptvseta_el->SetMinimum(0.);
  QCDFRptvseta_el->SetMaximum(0.2);
  QCDFRptvseta_el->Draw("LEGO2");
  c4->SaveAs("elFR2D.png");

  TCanvas *c5 = new TCanvas();
  TH2F *QCDFRErrptvseta_el = (TH2F*)gDirectory->Get("QCD_FRErrptvseta_el");
  QCDFRErrptvseta_el->SetMinimum(0.);
  QCDFRErrptvseta_el->SetMaximum(0.2);
  QCDFRErrptvseta_el->Draw("LEGO2");
  c5->SaveAs("elFRErr2D.png");
  
  TCanvas *c6 = new TCanvas();
  TH1F *QCD_FReta_el = (TH1F*)gDirectory->Get("QCD_FReta_el");
  QCD_FReta_el->SetMinimum(0.);
  QCD_FReta_el->SetMaximum(0.2);
  QCD_FReta_el->Draw();
  c6->SaveAs("elFReta.png");


  TCanvas *c7 = new TCanvas();
  TH1F *QCD_FRpt_el = gDirectory->Get("QCD_FRpt_el");
  QCD_FRpt_el->SetMinimum(0.);
  QCD_FRpt_el->SetMaximum(0.2);
  QCD_FRpt_el->Draw();
  c7->SaveAs("elFRpt.png");

  TCanvas *c8 = new TCanvas();
  WJets_numeta[0]->Draw("e");
  WJets_predictedeta[0]->Draw("samese");
  c8->SaveAs("elPredictedeta.png");
  
  TCanvas *c9 = new TCanvas();
  WJets_numpt[0]->Draw("e");
  WJets_predictedpt[0]->Draw("sames");
  c9->SaveAs("elPredictedpt.png");

  TCanvas *c10 = new TCanvas();
  ((TH1F*)gDirectory->Get("WJets_actualnJets_mu"))->SetMarkerStyle(3);
  (TH1F*)gDirectory->Get("WJets_actualnJets_mu")->Draw("e");
  ((TH1F*)gDirectory->Get("WJets_predictednJets_mu"))->SetMarkerStyle(25);
  (TH1F*)gDirectory->Get("WJets_predictednJets_mu")->Draw("samese");
  hist::setrangey(c3);
  c10->SaveAs("munJetsPredicted.png");

  TCanvas *c11 = new TCanvas();
  (TH1F*)gDirectory->Get("QCD_FRptvseta_mu")->Draw("LEGO2");
  c11->SaveAs("muFR2D.png");

  TCanvas *c12 = new TCanvas();
  TH2F *QCDFRErrptvseta_mu = (TH2F*)gDirectory->Get("QCD_FRErrptvseta_mu");
  QCDFRErrptvseta_mu->SetMaximum(0.2);
  QCDFRErrptvseta_mu->Draw("LEGO2");
  c12->SaveAs("muFRErr2D.png");
  
  TCanvas *c13 = new TCanvas();
  TH1F *QCD_FReta_mu = (TH1F*)gDirectory->Get("QCD_FReta_mu");
  QCD_FReta_mu->Draw();
  c13->SaveAs("muFReta.png");


  TCanvas *c14 = new TCanvas();
  TH1F *QCD_FRpt_mu = gDirectory->Get("QCD_FRpt_mu");
  QCD_FRpt_mu->Draw();
  c14->SaveAs("muFRpt.png");


  //Print the actual and predicted numbers with errors
  //do the errors for the WJets
  float predictionError[2] = {0.0, 0.0};
  for (unsigned int i=0; i < 2; i++) {
    Float_t totalErr = 0.;
    for(unsigned int ieta = 1; ieta < WJets_nJets3D[i]->GetNbinsY() + 1; ieta++) {
      for(unsigned int ipt = 1; ipt < WJets_nJets3D[i]->GetNbinsZ() + 1; ipt++) {
	Float_t temp23 = 0.;  
	for(unsigned int iJet = 1; iJet < WJets_nJets3D[i]->GetNbinsX() + 1; iJet++) {
	  temp23 = temp23 + WJets_nJets3D[i]->GetBinContent(iJet, ieta, ipt);
	}
	totalErr = pow(temp23,2) + totalErr;
      }
    }
    predictionError[i] = sqrt(totalErr);
  }
  //output summary
  cout << "Actual, predicted electron fakes, WJets: " << WJets_actualnJets[0]->GetEntries() 
       << ", " << WJets_predictednJets[0]->Integral() << " +/- " << predictionError[0] << endl;
  cout << "Actual, predicted muon fakes, WJets: " << WJets_actualnJets[1]->GetEntries() 
       << ", " << WJets_predictednJets[1]->Integral() << " +/- " << predictionError[1] << endl;

}
    
void doAll(){
  
  gROOT->ProcessLine(Form(".x setup.C(%d)", 1));
  gSystem->CompileMacro("QCDFRestimator.C","++k", "libQCDFRestimator");

  const char* location = gSystem->Getenv("CMS2_NTUPLE_LOCATION");

  TChain *ch_pthat30to80 = new TChain("Events");
  ch_pthat30to80->Add(Form("%s/%s",location,"cms2-V01-03-01/QCDpt30_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root"));
  //  ch_pthat30to80->Add(Form("%s/%s",location,"cms2-V01-03-01/QCDpt80_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root"));

//  TChain *ch_pthat80 = new TChain("Events");
//  ch_pthat80->Add("/data/tmp/cms2-V01-03-01/QCDpt80_Summer08_IDEAL_V11_redigi_v1-SingleLepton/merged_ntuple*.root");
  
  //WJets chain
  TChain *ch_WJets = new TChain("Events");
  ch_WJets->Add(Form("%s/%s",location,"cms2-V01-03-01/WJets-madgraph_Summer08_IDEAL_V11_redigi_v1/merged*.root"));
/*
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_10.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_11.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_12.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_13.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_14.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_1.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_2.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_3.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_4.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_5.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_6.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_7.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_8.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_9.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_0jet-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_1.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_2.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_1jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple_1.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_2jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_3jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_3jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_3jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_3jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_3jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_4jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_4jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_4jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_4jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_4jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_5jet_Pt0to100-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_5jet_Pt100to300-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_5jet_Pt1600toInf-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_5jet_Pt300to800-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
ch_WJets->Add("/data/tmp/cms2-V01-03-01/W_5jet_Pt800to1600-alpgen_Summer08_IDEAL_V12_RECOSIM_v1-SingleLepton/merged_ntuple.root ");
*/
  QCDFRestimator *looper = new QCDFRestimator();
  
  looper->ScanChainQCD(ch_pthat30to80, "QCD", 1.0, 1.0, 0.0, 99999999999999999999999999999.);
  cout << "Done QCD" << endl;

  //do WJets
  looper->ScanChainWJets(ch_WJets, "WJets", 1.0, 1.0);
  cout << "Done WJets" << endl;
  hist::saveHist("QCDFRplots.root");
  cout << "Done saving hists" << endl;
  delete looper;

  TestPrediction();


  TCanvas* c = new TCanvas();
  c->Divide(2,2);
  
  c->cd(1);
  QCD_FRpt_el->SetMinimum(0.0);
  QCD_FRpt_el->SetMaximum(0.2);
  QCD_FRpt_el->Draw();
    
  c->cd(2);
  QCD_FReta_el->SetMinimum(0.0);
  QCD_FReta_el->SetMaximum(0.2);
  QCD_FReta_el->Draw();
  
  c->cd(3);
  QCD_FRpt_mu->SetMinimum(0.0);
  QCD_FRpt_mu->SetMaximum(0.2);
  QCD_FRpt_mu->Draw();
  
  c->cd(4);
  QCD_FReta_mu->SetMinimum(0.0);
  QCD_FReta_mu->SetMaximum(0.2);
  QCD_FReta_mu->Draw();
    
}
    
