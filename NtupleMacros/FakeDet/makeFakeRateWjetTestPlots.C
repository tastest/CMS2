{

  gStyle->SetOptStat(1111111);
  gStyle->SetPalette(1,0);
  gROOT->UseCurrentStyle();

  bool subtract_realE = false;

  int rebinvalue_pt = 4;
  int rebinvalue_eta = 2;

  TFile *_file0 = TFile::Open("myWHist_elt_0.root");
  TFile *_file2 = TFile::Open("myWHist_elt_2.root");

  TCanvas * elept_em = new TCanvas("elept_em","elept_em");
  TCanvas * eleeta_em = new TCanvas("eleeta_em","eleeta_em");
  TCanvas * nJet_em = new TCanvas("nJet_em","nJet_em");
  TCanvas * MET_em = new TCanvas("MET_em","MET_em");

  TCanvas * elept_em_true = new TCanvas("elept_em_true","elept_em_true");

  // plot Pt comparison
  
  elept_em->cd();
  TH1F* wjets_helePt_em_observed = (TH1F*) ((_file0->Get("wjets_helePt_em"))->Clone("wjets_helePt_em_observed"));
  wjets_helePt_em_observed->Rebin(rebinvalue_pt);
  wjets_helePt_em_observed->SetLineColor(kRed);
  wjets_helePt_em_observed->SetLineWidth(2.);
  wjets_helePt_em_observed->Draw();

  if(subtract_realE){
    // also subtract true ele from observed
    TH1F* wjets_helePt_em_true_trueele_observed = (TH1F*) ((_file0->Get("wjets_helePtTrue_em"))->Clone("wjets_helePt_em_true_trueele_observed"));
    wjets_helePt_em_true_trueele_observed->Rebin(rebinvalue_pt);
    wjets_helePt_em_true_trueele_observed->SetLineColor(8);
    wjets_helePt_em_true_trueele_observed->SetLineWidth(2.);
    //    wjets_helePt_em_true_trueele_observed->Draw("sames");

    TH1F* observed_minusTrueele = wjets_helePt_em_observed->Clone("observed_minusTrueele");
    observed_minusTrueele->Add(wjets_helePt_em_true_trueele_observed,-1);
    observed_minusTrueele->SetLineColor(kYellow);
    observed_minusTrueele->SetLineWidth(2.);
    observed_minusTrueele->Draw("sames");
  }

  TH1F* wjets_helePt_em_predicted = (TH1F*) ((_file2->Get("wjets_helePt_em"))->Clone("wjets_helePt_em_predicted"));
  wjets_helePt_em_predicted->Rebin(rebinvalue_pt);
  wjets_helePt_em_predicted->SetLineColor(kBlue);
  wjets_helePt_em_predicted->SetLineWidth(2.);
  wjets_helePt_em_predicted->Draw("sames");

  if(subtract_realE){
    // subtract from prediction
    TH1F* wjets_helePt_em_true_trueele_predicted = (TH1F*) ((_file2->Get("wjets_helePtTrue_em"))->Clone("wjets_helePt_em_true_trueele_predicted"));
    wjets_helePt_em_true_trueele_predicted->Rebin(rebinvalue_pt);
    wjets_helePt_em_true_trueele_predicted->SetLineColor(8);
    wjets_helePt_em_true_trueele_predicted->SetLineWidth(2.);
    //    wjets_helePt_em_true_trueele_predicted->Draw("sames");

    TH1F* predicted_minusTrueele = wjets_helePt_em_predicted->Clone("predicted_minusTrueele");
    predicted_minusTrueele->Add(wjets_helePt_em_true_trueele_predicted,-1);
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
  TH1F* wjets_heleEta_em_observed = (TH1F*) ((_file0->Get("wjets_heleEta_em"))->Clone("wjets_heleEta_em_observed"));
  wjets_heleEta_em_observed->Rebin(rebinvalue_eta);
  wjets_heleEta_em_observed->SetLineColor(kRed);
  wjets_heleEta_em_observed->SetLineWidth(2.);
  wjets_heleEta_em_observed->Draw();

  if(subtract_realE && 42 != 42){
    // also need to subtract realE from prediction!
    // eta plot not implemented for true ele!!!
    TH1F* wjets_heleEta_em_true_trueele = (TH1F*) ((_file0->Get("wjets_heleEtaTrue_em"))->Clone("wjets_heleEta_em_true_trueele"));
    wjets_heleEta_em_true_trueele->Rebin(rebinvalue_eta);
    wjets_heleEta_em_true_trueele->SetLineColor(8);
    wjets_heleEta_em_true_trueele->SetLineWidth(2.);
    wjets_heleEta_em_true_trueele->Draw("sames");

    TH1F* subtract = wjets_heleEta_em_observed->Clone("subtract");
    subtract->Add(wjets_heleEta_em_true_trueele,-1);
    subtract->SetLineColor(kYellow);
    subtract->SetLineWidth(2.);
    subtract->Draw("sames");
  }

  TH1F* wjets_heleEta_em_predicted = (TH1F*) ((_file2->Get("wjets_heleEta_em"))->Clone("wjets_heleEta_em_predicted"));
  wjets_heleEta_em_predicted->Rebin(rebinvalue_eta);
  wjets_heleEta_em_predicted->SetLineColor(kBlue);
  wjets_heleEta_em_predicted->SetLineWidth(2.);
  wjets_heleEta_em_predicted->Draw("sames");

   lable_new->Draw();
   lable_ref->Draw();

  // plot nJet comparison
  nJet_em->cd();
  nJet_em->SetLogy();
  TH1F* wjets_hnJet_em_observed = (TH1F*) ((_file0->Get("wjets_hnJet_em"))->Clone("wjets_hnJet_em_observed"));
  //  wjets_hnJet_em_observed->Rebin(5);
  wjets_hnJet_em_observed->SetLineColor(kRed);
  wjets_hnJet_em_observed->SetLineWidth(2.);
  wjets_hnJet_em_observed->Draw();

  TH1F* wjets_hnJet_em_predicted = (TH1F*) ((_file2->Get("wjets_hnJet_em"))->Clone("wjets_hnJet_em_predicted"));
  //  wjets_hnJet_em_predicted->Rebin(5);
  wjets_hnJet_em_predicted->SetLineColor(kBlue);
  wjets_hnJet_em_predicted->SetLineWidth(2.);
  wjets_hnJet_em_predicted->Draw("sames");
   lable_new->Draw();
   lable_ref->Draw();

  // plot MET comparison
  MET_em->cd();
  TH1F* wjets_hMET_em_observed = (TH1F*) ((_file0->Get("wjets_hMET_em"))->Clone("wjets_hMET_em_observed"));
  //  wjets_hMET_em_observed->Rebin(5);
  wjets_hMET_em_observed->SetLineColor(kRed);
  wjets_hMET_em_observed->SetLineWidth(2.);
  wjets_hMET_em_observed->Draw();

  TH1F* wjets_hMET_em_predicted = (TH1F*) ((_file2->Get("wjets_hMET_em"))->Clone("wjets_hMET_em_predicted"));
  //  wjets_hMET_em_predicted->Rebin(5);
  wjets_hMET_em_predicted->SetLineColor(kBlue);
  wjets_hMET_em_predicted->SetLineWidth(2.);
  wjets_hMET_em_predicted->Draw("sames");
   lable_new->Draw();
   lable_ref->Draw();

  // Draw the electron / fake contributions of the plotted fake rate:
  elept_em_true->cd();
  TH1F* wjets_helePt_em_true_allfakeprod = (TH1F*) ((_file2->Get("wjets_helePt_em"))->Clone("wjets_helePt_em_true_allfakeprod"));
  //  wjets_helePt_em_true_allfakeprod->Rebin(5);
  wjets_helePt_em_true_allfakeprod->SetLineColor(kBlue);
  wjets_helePt_em_true_allfakeprod->SetLineWidth(2.);
  wjets_helePt_em_true_allfakeprod->Draw();

  TH1F* wjets_helePt_em_true_trueele = (TH1F*) ((_file2->Get("wjets_helePtTrue_em"))->Clone("wjets_helePt_em_true_trueele"));
  //  wjets_helePt_em_true_trueele->Rebin(5);
  wjets_helePt_em_true_trueele->SetLineColor(8);
  wjets_helePt_em_true_trueele->SetLineWidth(2.);
  wjets_helePt_em_true_trueele->Draw("sames");

  TH1F* wjets_helePt_em_trueW_trueWele = (TH1F*) ((_file2->Get("wjets_helePtTrueEW_em"))->Clone("wjets_helePt_em_trueW_trueWele"));
  //  wjets_helePt_em_trueW_trueWele->Rebin(5);
  wjets_helePt_em_trueW_trueWele->SetLineColor(kYellow);
  wjets_helePt_em_trueW_trueWele->SetLineWidth(2.);
  wjets_helePt_em_trueW_trueWele->Draw("sames");

  TH1F* wjets_helePt_em_fake_fakeele = (TH1F*) ((_file2->Get("wjets_helePtFake_em"))->Clone("wjets_helePt_em_fake_fakeele"));
  //  wjets_helePt_em_fake_fakeele->Rebin(5);
  wjets_helePt_em_fake_fakeele->SetLineColor(kRed);
  wjets_helePt_em_fake_fakeele->SetLineWidth(2.);
  wjets_helePt_em_fake_fakeele->Draw("sames");

}
