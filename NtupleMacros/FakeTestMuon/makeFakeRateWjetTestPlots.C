{

  gStyle->SetOptStat(1111111);
  gStyle->SetPalette(1,0);
  gROOT->UseCurrentStyle();

  int rebinvalue_pt = 1;
  int rebinvalue_eta = 1;

  TFile *_file0 = TFile::Open("Results.root");
  TFile *_file2 = TFile::Open("Wjets_Fakerate.root");

  TCanvas * muopt_em = new TCanvas("muopt_em","muopt_em");
  TCanvas * muoeta_em = new TCanvas("muoeta_em","muoeta_em");
  TCanvas * nJet_em = new TCanvas("nJet_em","nJet_em");
  TCanvas * MET_em = new TCanvas("MET_em","MET_em");

  TCanvas * muopt_em_true = new TCanvas("muopt_em_true","muopt_em_true");

  // plot Pt comparison
  
  muopt_em->cd();
  TH1F* wjets_muPt_em_observed = (TH1F*) ((_file0->Get("wjets_muPt_em"))->Clone("wjets_muPt_em_observed"));
  wjets_muPt_em_observed->Rebin(rebinvalue_pt);
  wjets_muPt_em_observed->SetLineColor(kRed);
  wjets_muPt_em_observed->SetLineWidth(2.);
  wjets_muPt_em_observed->Draw();
  TH1F* wjets_muPt_em_predicted = (TH1F*) ((_file2->Get("wjets_muPt_em"))->Clone("wjets_muPt_em_predicted"));
  wjets_muPt_em_predicted->Rebin(rebinvalue_pt);
  wjets_muPt_em_predicted->SetLineColor(kBlue);
  wjets_muPt_em_predicted->SetLineWidth(2.);
  wjets_muPt_em_predicted->Draw("sames");
  muopt_em->Update();
  TPaveStats *stats_mupt = (TPaveStats*)(wjets_muPt_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_mupt->SetX1NDC(0.78);
  stats_mupt->SetY1NDC(0.43);
  stats_mupt->SetX2NDC(0.98);
  stats_mupt->SetY2NDC(0.71);
  muopt_em->Update();
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
  muopt_em->SaveAs("muo_pt.png");

  // Plot eta comparison
  muoeta_em->cd();
  TH1F* wjets_muEta_em_observed = (TH1F*) ((_file0->Get("wjets_muEta_em"))->Clone("wjets_muEta_em_observed"));
  wjets_muEta_em_observed->Rebin(rebinvalue_eta);
  wjets_muEta_em_observed->SetLineColor(kRed);
  wjets_muEta_em_observed->SetLineWidth(2.);
  wjets_muEta_em_observed->Draw();
  TH1F* wjets_muEta_em_predicted = (TH1F*) ((_file2->Get("wjets_muEta_em"))->Clone("wjets_muEta_em_predicted"));
  wjets_muEta_em_predicted->Rebin(rebinvalue_eta);
  wjets_muEta_em_predicted->SetLineColor(kBlue);
  wjets_muEta_em_predicted->SetLineWidth(2.);
  wjets_muEta_em_predicted->Draw("sames");
  muoeta_em->Update();
  TPaveStats *stats_muoeta = (TPaveStats*)(wjets_muEta_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_muoeta->SetX1NDC(0.78);
  stats_muoeta->SetY1NDC(0.43);
  stats_muoeta->SetX2NDC(0.98);
  stats_muoeta->SetY2NDC(0.71);
  muoeta_em->Update();
  lable_new->Draw();
  lable_ref->Draw();
  muoeta_em->SaveAs("muo_eta.png");

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
  nJet_em->Update();
  TPaveStats *stats_nJet = (TPaveStats*)(wjets_nJet_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_nJet->SetX1NDC(0.78);
  stats_nJet->SetY1NDC(0.43);
  stats_nJet->SetX2NDC(0.98);
  stats_nJet->SetY2NDC(0.71);
  nJet_em->Update();
  lable_new->Draw();
  lable_ref->Draw();
  nJet_em->SaveAs("muo_nJet.png");

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
  MET_em->Update();
  TPaveStats *stats_MET = (TPaveStats*)(wjets_met_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_MET->SetX1NDC(0.78);
  stats_MET->SetY1NDC(0.43);
  stats_MET->SetX2NDC(0.98);
  stats_MET->SetY2NDC(0.71);
  MET_em->Update();
  lable_new->Draw();
  lable_ref->Draw();
  MET_em->SaveAs("muo_MET.png");

}
