#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include <iostream>
#include "histools.h"
#include "QCDFRestimator.h"
#endif __CINT__
#include "overlay.C"

void TestPrediction() {
  
  //  gROOT->ProcessLine(".L setup.C");
  //  gROOT->ProcessLine("setup()");
  //  gROOT->ProcessLine(Form(".x setup.C(%d)", 1));
  hist::deleteHistos();
  hist::loadHist("QCDFRplots.root");
   
  TH3F *WJets_nJets3D[2];

  TH1F *WJets_actualnJets[2];
  TH1F *WJets_predictednJets[2];
  
  char *flavor[2] = {"el", "mu"};
  for(unsigned int i = 0; i < 2; i++) {
    WJets_nJets3D[i] = (TH3F*)gDirectory->Get(Form("WJets_nJets3D_%s", flavor[i]));

    WJets_actualnJets[i] = (TH1F*)gDirectory->Get(Form("WJets_actualnJets_%s", flavor[i]));
    WJets_predictednJets[i] = (TH1F*)gDirectory->Get(Form("WJets_predictednJets_%s", flavor[i]));
  }      
  
  hist::color("WJets_actual", kRed);
  hist::color("WJets_num", kRed);


  TH3F *WJets_TrueCat3D[2];

  TH1F *WJets_actualTrueCat[2];
  TH1F *WJets_predictedTrueCat[2];
  
  char *flavor[2] = {"el", "mu"};
  for(unsigned int i = 0; i < 2; i++) {
    WJets_TrueCat3D[i] = (TH3F*)gDirectory->Get(Form("WJets_TrueCat3D_%s", flavor[i]));

    WJets_actualTrueCat[i] = (TH1F*)gDirectory->Get(Form("WJets_actualTrueCat_%s", flavor[i]));
    WJets_predictedTrueCat[i] = (TH1F*)gDirectory->Get(Form("WJets_predictedTrueCat_%s", flavor[i]));
  }      
  
  hist::color("WJets_actual", kRed);
  hist::color("WJets_num", kRed);
  hist::color("TTbar_actual", kRed);
  hist::color("TTbar_num", kRed);

    
  //TCanvas *c1 = new TCanvas();
  //WJets_numeta[0]->Draw("e");
  //WJets_predictedeta[0]->Draw("samese");
  //hist::setrangey(c1);
  //c1->SaveAs("elPredictedeta.png");
  
  //TCanvas *c2 = new TCanvas();
  //WJets_numpt[0]->Draw("e");
  //WJets_predictedpt[0]->Draw("sames");
  //hist::setrangey(c2);
  //c2->SaveAs("elPredictedpt.png");


  //TCanvas *c21 = new TCanvas();
  //WJets_numeta[1]->Draw("e");
  //WJets_predictedeta[1]->Draw("samese");
  //hist::setrangey(c21);
  //c21->SaveAs("muPredictedeta.png");

  //TCanvas *c22 = new TCanvas();
  //WJets_numpt[1]->Draw("e");
  //WJets_predictedpt[1]->Draw("sames");
  //hist::setrangey(c22);
  //c22->SaveAs("muPredictedpt.png");

  // dbarge
  overlay2( (TH1F*)gDirectory->Get("WJets_actualnJets_el"), 
            (TH1F*)gDirectory->Get("WJets_predictednJets_el"), "elWJetsnJetsPredicted.png");
  overlay2( (TH1F*)gDirectory->Get("WJets_actualTrueCat_el"), 
            (TH1F*)gDirectory->Get("WJets_predictedTrueCat_el"), "elWJetTrueCatPredicted.png");
  overlay2( (TH1F*)gDirectory->Get("TTbar_actualnJets_el"), 
            (TH1F*)gDirectory->Get("TTbar_predictednJets_el"), "elTTbarnJetsPredicted.png");
  overlay2( (TH1F*)gDirectory->Get("TTbar_actualTrueCat_el"), 
            (TH1F*)gDirectory->Get("TTbar_predictedTrueCat_el"), "elTTbarTrueCatPredicted.png");

  TCanvas *c4 = new TCanvas("c4","c4",1280,960);
  TH2F *QCDFRptvseta_el = (TH2F*)gDirectory->Get("QCD_FRptvseta_el");
  QCDFRptvseta_el->SetMinimum(0.);
  QCDFRptvseta_el->SetMaximum(0.2);
  QCDFRptvseta_el->Draw("LEGO2");
  c4->SaveAs("elFR2D.png");

  TCanvas *c5 = new TCanvas("c5","c5",1280,960);
  TH2F *QCDFRErrptvseta_el = (TH2F*)gDirectory->Get("QCD_FRErrptvseta_el");
  QCDFRErrptvseta_el->SetMinimum(0.);
  QCDFRErrptvseta_el->SetMaximum(0.2);
  QCDFRErrptvseta_el->Draw("LEGO2");
  c5->SaveAs("elFRErr2D.png");
  
  TCanvas *c6 = new TCanvas("c6","c6",1280,960);
  TH1F *QCD_FReta_el = (TH1F*)gDirectory->Get("QCD_FReta_el");
  QCD_FReta_el->SetMinimum(0.);
  QCD_FReta_el->SetMaximum(0.2);
  QCD_FReta_el->Draw();
  c6->SaveAs("elFReta.png");


  TCanvas *c7 = new TCanvas("c7","c7",1280,960);
  TH1F *QCD_FRpt_el = gDirectory->Get("QCD_FRpt_el");
  QCD_FRpt_el->SetMinimum(0.);
  QCD_FRpt_el->SetMaximum(0.2);
  QCD_FRpt_el->Draw();
  c7->SaveAs("elFRpt.png");

  //TCanvas *c8 = new TCanvas();
  //WJets_numeta[0]->Draw("e");
  //WJets_predictedeta[0]->Draw("samese");
  //c8->SaveAs("elPredictedeta.png");
  
  //TCanvas *c9 = new TCanvas();
  //WJets_numpt[0]->Draw("e");
  //WJets_predictedpt[0]->Draw("sames");
  //c9->SaveAs("elPredictedpt.png");

  TCanvas *c10 = new TCanvas("c10","c10",1280,960);
  ((TH1F*)gDirectory->Get("WJets_actualnJets_mu"))->SetMarkerStyle(3);
  (TH1F*)gDirectory->Get("WJets_actualnJets_mu")->Draw("e");
  ((TH1F*)gDirectory->Get("WJets_predictednJets_mu"))->SetMarkerStyle(25);
  (TH1F*)gDirectory->Get("WJets_predictednJets_mu")->Draw("samese");
  //hist::setrangey(c3);
  c10->SaveAs("munJetsPredicted.png");

  TCanvas *c11 = new TCanvas("c11","c11",1280,960);
  TH2F *QCDFRptvseta_mu = (TH2F*)gDirectory->Get("QCD_FRptvseta_mu");
  QCDFRptvseta_mu->SetMinimum(0.);
  QCDFRptvseta_mu->SetMaximum(0.2);
  QCDFRptvseta_mu->Draw("LEGO2");
  c11->SaveAs("muFR2D.png");

  TCanvas *c12 = new TCanvas("c12","c12",1280,960);
  TH2F *QCDFRErrptvseta_mu = (TH2F*)gDirectory->Get("QCD_FRErrptvseta_mu");
  QCDFRErrptvseta_mu->SetMaximum(0.2);
  QCDFRErrptvseta_mu->Draw("LEGO2");
  c12->SaveAs("muFRErr2D.png");
  
  TCanvas *c13 = new TCanvas("c13","c13",1280,960);
  TH1F *QCD_FReta_mu = (TH1F*)gDirectory->Get("QCD_FReta_mu");
  QCD_FReta_mu->Draw();
  c13->SaveAs("muFReta.png");

  TCanvas *c14 = new TCanvas("c14","c14",1280,960);
  TH1F *QCD_FRpt_mu = gDirectory->Get("QCD_FRpt_mu");
  QCD_FRpt_mu->Draw();
  c14->SaveAs("muFRpt.png");

  TCanvas *c15 = new TCanvas("c15","c15",1280,960);
  c15->Divide(3,3);
  c15->cd(1);
  TH1F *h_el_WJnum = gDirectory->Get("WJets_truecomposition_num_el");
  h_el_WJnum->SetMaximum(1.1);
  h_el_WJnum->DrawNormalized("9");
  c15->cd(2);
  TH1F *h_el_WJdenom = gDirectory->Get("WJets_truecomposition_denom_el");
  h_el_WJdenom->SetMaximum(1.1);
  h_el_WJdenom->DrawNormalized("9");
  c15->cd(3);
  TH1F *h_el_WJratio = gDirectory->Get("WJets_truecomposition_ratio_el");
  h_el_WJratio->Draw("9");
  c15->cd(4);
  TH1F *h_el_TTnum = gDirectory->Get("TTbar_truecomposition_num_el");
  h_el_TTnum->SetMaximum(1.1);
  h_el_TTnum->DrawNormalized("9");
  c15->cd(5);
  TH1F *h_el_TTdenom = gDirectory->Get("TTbar_truecomposition_denom_el");
  h_el_TTdenom->SetMaximum(1.1);
  h_el_TTdenom->DrawNormalized("9");
  c15->cd(6);
  TH1F *h_el_TTratio = gDirectory->Get("TTbar_truecomposition_ratio_el");
  h_el_TTratio->Draw("9");
  c15->cd(7);
  TH1F *h_el_QCDnum = gDirectory->Get("QCD_truecomposition_num_el");
  h_el_QCDnum->SetMaximum(1.1);
  h_el_QCDnum->DrawNormalized("9");
  c15->cd(8);
  TH1F *h_el_QCDdenom = gDirectory->Get("QCD_truecomposition_denom_el");
  h_el_QCDdenom->SetMaximum(1.1);
  h_el_QCDdenom->DrawNormalized("9");
  c15->cd(9);
  TH1F *h_el_QCDratio = gDirectory->Get("QCD_truecomposition_ratio_el");
  h_el_QCDratio->Draw("9");
  c15->Draw("9");
  c15->SaveAs("eFakeComposition.png");




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
  ch_pthat30to80->Add(Form("%s/%s",location,"cms2/QCD_Pt30_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root"));
  ch_pthat30to80->Add(Form("%s/%s",location,"cms2/QCD_Pt80_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root"));

  //WJets chain
  TChain *ch_WJets = new TChain("Events");
  ch_WJets->Add(Form("%s/%s",location,"cms2/WJets-madgraph_Summer09-MC_31X_V3_7TeV-v1/V03-00-35/merged_ntuple*.root"));

  //TTbar chain
  TChain *ch_TTbar = new TChain("Events");
  ch_TTbar->Add(Form("%s/%s",location,"cms2/TTbar_Summer09-MC_31X_V3_7TeV-v1/V03-00-34/merged_ntuple*.root"));

  QCDFRestimator *looper = new QCDFRestimator();
  
  looper->ScanChainQCD(ch_pthat30to80, "QCD", 1.0, 1.0, 0.0, 99999999999999999999999999999.);
  cout << "Done QCD" << endl;

  //do WJets
  looper->ScanChainAppTest(ch_WJets, "WJets", 1.0, 1.0);
  cout << "Done WJets" << endl;

  //do TTbar
  looper->ScanChainAppTest(ch_TTbar, "TTbar", 1.0, 1.0);
  cout << "Done TTbar" << endl;

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
    
