{

  gStyle->SetOptStat(1111111);
  gStyle->SetPalette(1,0);
  gROOT->UseCurrentStyle();

  bool subtract_realE = false;

  int rebinvalue_pt = 1;
  int rebinvalue_eta = 1;

  TFile *_file0 = TFile::Open("Results.root");
  TFile *_file2 = TFile::Open("Wjets_Fakerate.root");

  TCanvas * elept_em = new TCanvas("elept_em","elept_em");
  TCanvas * eleeta_em = new TCanvas("eleeta_em","eleeta_em");
  TCanvas * nJet_em = new TCanvas("nJet_em","nJet_em");
  TCanvas * MET_em = new TCanvas("MET_em","MET_em");

  TCanvas * elept_em_true = new TCanvas("elept_em_true","elept_em_true");

  // plot Pt comparison
  
  elept_em->cd();
  TH1F* wjets_elPt_em_observed = (TH1F*) ((_file0->Get("wjets_elPt_em"))->Clone("wjets_elPt_em_observed"));
  wjets_elPt_em_observed->Rebin(rebinvalue_pt);
  wjets_elPt_em_observed->SetLineColor(kRed);
  wjets_elPt_em_observed->SetLineWidth(2.);
  wjets_elPt_em_observed->Draw();

  if(subtract_realE){
    // also subtract true ele from observed
    TH1F* wjets_elPt_em_true_trueele_observed = (TH1F*) ((_file0->Get("wjets_elPtTrue_em"))->Clone("wjets_elPt_em_true_trueele_observed"));
    wjets_elPt_em_true_trueele_observed->Rebin(rebinvalue_pt);
    wjets_elPt_em_true_trueele_observed->SetLineColor(8);
    wjets_elPt_em_true_trueele_observed->SetLineWidth(2.);
    //    wjets_elPt_em_true_trueele_observed->Draw("sames");

    TH1F* observed_minusTrueele = wjets_elPt_em_observed->Clone("observed_minusTrueele");
    observed_minusTrueele->Add(wjets_elPt_em_true_trueele_observed,-1);
    observed_minusTrueele->SetLineColor(kYellow);
    observed_minusTrueele->SetLineWidth(2.);
    observed_minusTrueele->Draw("sames");
  }

  TH1F* wjets_elPt_em_predicted = (TH1F*) ((_file2->Get("wjets_elPt_em"))->Clone("wjets_elPt_em_predicted"));
  wjets_elPt_em_predicted->Rebin(rebinvalue_pt);
  wjets_elPt_em_predicted->SetLineColor(kBlue);
  wjets_elPt_em_predicted->SetLineWidth(2.);
  wjets_elPt_em_predicted->Draw("sames");

  if(subtract_realE){
    // subtract from prediction
    TH1F* wjets_elPt_em_true_trueele_predicted = (TH1F*) ((_file2->Get("wjets_elPtTrue_em"))->Clone("wjets_elPt_em_true_trueele_predicted"));
    wjets_elPt_em_true_trueele_predicted->Rebin(rebinvalue_pt);
    wjets_elPt_em_true_trueele_predicted->SetLineColor(8);
    wjets_elPt_em_true_trueele_predicted->SetLineWidth(2.);
    //    wjets_elPt_em_true_trueele_predicted->Draw("sames");

    TH1F* predicted_minusTrueele = wjets_elPt_em_predicted->Clone("predicted_minusTrueele");
    predicted_minusTrueele->Add(wjets_elPt_em_true_trueele_predicted,-1);
    predicted_minusTrueele->SetLineColor(kMagenta);
    predicted_minusTrueele->SetLineWidth(2.);
    predicted_minusTrueele->Draw("sames");
  }

   TLatex *   lable_ref = new TLatex(0.55,0.75,"Red: observed");
   lable_ref->SetNDC();
   lable_ref->SetTextSize(0.04);
   lable_ref->SetTextColor(kRed);
   TLatex *   lable_new = new TLatex(0.55,0.7,"Blue: predicted");
   lable_new->SetNDC();
   lable_new->SetTextSize(0.04);
   lable_new->SetTextColor(kBlue);
   lable_new->Draw();
   lable_ref->Draw();

  // Plot eta comparison
  eleeta_em->cd();
  TH1F* wjets_elEta_em_observed = (TH1F*) ((_file0->Get("wjets_elEta_em"))->Clone("wjets_elEta_em_observed"));
  wjets_elEta_em_observed->Rebin(rebinvalue_eta);
  wjets_elEta_em_observed->SetLineColor(kRed);
  wjets_elEta_em_observed->SetLineWidth(2.);
  wjets_elEta_em_observed->Draw();

  if(subtract_realE && 42 != 42){
    // also need to subtract realE from prediction!
    // eta plot not implemented for true ele!!!
    TH1F* wjets_elEta_em_true_trueele = (TH1F*) ((_file0->Get("wjets_elEtaTrue_em"))->Clone("wjets_elEta_em_true_trueele"));
    wjets_elEta_em_true_trueele->Rebin(rebinvalue_eta);
    wjets_elEta_em_true_trueele->SetLineColor(8);
    wjets_elEta_em_true_trueele->SetLineWidth(2.);
    wjets_elEta_em_true_trueele->Draw("sames");

    TH1F* subtract = wjets_elEta_em_observed->Clone("subtract");
    subtract->Add(wjets_elEta_em_true_trueele,-1);
    subtract->SetLineColor(kYellow);
    subtract->SetLineWidth(2.);
    subtract->Draw("sames");
  }

  TH1F* wjets_elEta_em_predicted = (TH1F*) ((_file2->Get("wjets_elEta_em"))->Clone("wjets_elEta_em_predicted"));
  wjets_elEta_em_predicted->Rebin(rebinvalue_eta);
  wjets_elEta_em_predicted->SetLineColor(kBlue);
  wjets_elEta_em_predicted->SetLineWidth(2.);
  wjets_elEta_em_predicted->Draw("sames");

   lable_new->Draw();
   lable_ref->Draw();

  // plot nJet comparison
  nJet_em->cd();
  nJet_em->SetLogy();
  TH1F* wjets_nJet_em_observed = (TH1F*) ((_file0->Get("wjets_nJet_em"))->Clone("wjets_nJet_em_observed"));
  //  wjets_nJet_em_observed->Rebin(5);
  wjets_nJet_em_observed->SetLineColor(kRed);
  wjets_nJet_em_observed->SetLineWidth(2.);
  wjets_nJet_em_observed->Draw();

  TH1F* wjets_nJet_em_predicted = (TH1F*) ((_file2->Get("wjets_nJet_em"))->Clone("wjets_nJet_em_predicted"));
  //  wjets_nJet_em_predicted->Rebin(5);
  wjets_nJet_em_predicted->SetLineColor(kBlue);
  wjets_nJet_em_predicted->SetLineWidth(2.);
  wjets_nJet_em_predicted->Draw("sames");
   lable_new->Draw();
   lable_ref->Draw();

  // plot MET comparison
  MET_em->cd();
  TH1F* wjets_met_em_observed = (TH1F*) ((_file0->Get("wjets_met_em"))->Clone("wjets_met_em_observed"));
  //  wjets_met_em_observed->Rebin(5);
  wjets_met_em_observed->SetLineColor(kRed);
  wjets_met_em_observed->SetLineWidth(2.);
  wjets_met_em_observed->Draw();

  TH1F* wjets_met_em_predicted = (TH1F*) ((_file2->Get("wjets_met_em"))->Clone("wjets_met_em_predicted"));
  //  wjets_met_em_predicted->Rebin(5);
  wjets_met_em_predicted->SetLineColor(kBlue);
  wjets_met_em_predicted->SetLineWidth(2.);
  wjets_met_em_predicted->Draw("sames");
   lable_new->Draw();
   lable_ref->Draw();

  // Draw the electron / fake contributions of the plotted fake rate:
  elept_em_true->cd();
  TH1F* wjets_elPt_em_true_allfakeprod = (TH1F*) ((_file2->Get("wjets_elPt_em"))->Clone("wjets_elPt_em_true_allfakeprod"));
  //  wjets_elPt_em_true_allfakeprod->Rebin(5);
  wjets_elPt_em_true_allfakeprod->SetLineColor(kBlue);
  wjets_elPt_em_true_allfakeprod->SetLineWidth(2.);
  wjets_elPt_em_true_allfakeprod->Draw();

}
