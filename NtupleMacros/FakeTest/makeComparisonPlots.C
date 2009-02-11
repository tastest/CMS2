{

  gStyle->SetOptStat(1111111);
  gStyle->SetPalette(1,0);
  gROOT->UseCurrentStyle();

  int rebinvalue_pt = 1;
  int rebinvalue_eta = 1;

  TFile *_file0 = TFile::Open("090211_03/Results.root");
  TFile *_file2 = TFile::Open("090211_02/Results.root");

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
  TH1F* wjets_elPt_em_predicted = (TH1F*) ((_file2->Get("wjets_elPt_em"))->Clone("wjets_elPt_em_predicted"));
  wjets_elPt_em_predicted->Rebin(rebinvalue_pt);
  wjets_elPt_em_predicted->SetLineColor(kBlue);
  wjets_elPt_em_predicted->SetLineWidth(2.);
  wjets_elPt_em_predicted->Draw("sames");
  elept_em->Update();
  TPaveStats *stats_elpt = (TPaveStats*)(wjets_elPt_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_elpt->SetX1NDC(0.78);
  stats_elpt->SetY1NDC(0.43);
  stats_elpt->SetX2NDC(0.98);
  stats_elpt->SetY2NDC(0.71);
  elept_em->Update();
  TLatex *   lable_ref = new TLatex(0.55,0.75,"Red: with std::abs");
  lable_ref->SetNDC();
  lable_ref->SetTextSize(0.04);
  lable_ref->SetTextColor(kRed);
  TLatex *   lable_new = new TLatex(0.55,0.7,"Blue: with TMath::Abs");
  lable_new->SetNDC();
  lable_new->SetTextSize(0.04);
  lable_new->SetTextColor(kBlue);
  lable_new->Draw();
  lable_ref->Draw();
  elept_em->SaveAs("ele_pt.png");

  // Plot eta comparison
  eleeta_em->cd();
  TH1F* wjets_elEta_em_observed = (TH1F*) ((_file0->Get("wjets_elEta_em"))->Clone("wjets_elEta_em_observed"));
  wjets_elEta_em_observed->Rebin(rebinvalue_eta);
  wjets_elEta_em_observed->SetLineColor(kRed);
  wjets_elEta_em_observed->SetLineWidth(2.);
  wjets_elEta_em_observed->Draw();
  TH1F* wjets_elEta_em_predicted = (TH1F*) ((_file2->Get("wjets_elEta_em"))->Clone("wjets_elEta_em_predicted"));
  wjets_elEta_em_predicted->Rebin(rebinvalue_eta);
  wjets_elEta_em_predicted->SetLineColor(kBlue);
  wjets_elEta_em_predicted->SetLineWidth(2.);
  wjets_elEta_em_predicted->Draw("sames");
  eleeta_em->Update();
  TPaveStats *stats_eleeta = (TPaveStats*)(wjets_elEta_em_predicted->GetListOfFunctions()->FindObject("stats"));
  stats_eleeta->SetX1NDC(0.78);
  stats_eleeta->SetY1NDC(0.43);
  stats_eleeta->SetX2NDC(0.98);
  stats_eleeta->SetY2NDC(0.71);
  eleeta_em->Update();
  lable_new->Draw();
  lable_ref->Draw();
  eleeta_em->SaveAs("ele_eta.png");

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
  nJet_em->SaveAs("ele_nJet.png");

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
  MET_em->SaveAs("ele_MET.png");

}
