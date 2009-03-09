{
  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111111);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1,0);

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

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

  //  TCanvas * elept_em_true = new TCanvas("elept_em_true","elept_em_true");

  // plot Pt comparison
  
  elept_em->cd();
  TH1F* wjets_helePt_em_observed = (TH1F*) ((_file0->Get("wjets_elPt_em"))->Clone("wjets_helePt_em_observed"));
  wjets_helePt_em_observed->Rebin(rebinvalue_pt);
  wjets_helePt_em_observed->SetLineColor(kRed);
  wjets_helePt_em_observed->SetFillColor(kWhite);
  wjets_helePt_em_observed->SetLineWidth(2.);
  wjets_helePt_em_observed->SetName("observed");
  wjets_helePt_em_observed->GetYaxis()->SetTitle("Events");
  wjets_helePt_em_observed->GetYaxis()->SetTitleOffset(1.2);
  wjets_helePt_em_observed->GetYaxis()->SetTitleSize(0.04);
  wjets_helePt_em_observed->GetXaxis()->SetTitle("p_{T}^{electron} (GeV)");
  wjets_helePt_em_observed->GetXaxis()->SetTitleOffset(1.2);
  wjets_helePt_em_observed->GetXaxis()->SetTitleSize(0.04);
  wjets_helePt_em_observed->SetMarkerStyle(20);
  wjets_helePt_em_observed->SetMarkerColor(kRed);
  wjets_helePt_em_observed->SetMarkerSize(1.1);
  wjets_helePt_em_observed->Draw();


  TH1F* wjets_helePt_em_predicted = (TH1F*) ((_file2->Get("wjets_elPt_em"))->Clone("wjets_helePt_em_predicted"));
  wjets_helePt_em_predicted->Rebin(rebinvalue_pt);
  wjets_helePt_em_predicted->SetLineColor(kBlue);
  wjets_helePt_em_predicted->SetFillColor(kWhite);
  wjets_helePt_em_predicted->SetLineWidth(2.);
  wjets_helePt_em_predicted->SetName("predicted");
  wjets_helePt_em_predicted->SetMarkerStyle(28);
  wjets_helePt_em_predicted->SetMarkerColor(kBlue);
  wjets_helePt_em_predicted->SetMarkerSize(1.2);
  wjets_helePt_em_predicted->Draw("sames");

  leg = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  //  leg->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leg->AddEntry(wjets_helePt_em_observed,       "Observed","lpf");
  leg->AddEntry(wjets_helePt_em_predicted,      "Predicted","lpf");
  leg->Draw();


  elept_em->Update();

  TPaveStats *fake_rate_stats = (TPaveStats*)(wjets_helePt_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats != 0 ) {
    fake_rate_stats->SetX1NDC(0.8);
    fake_rate_stats->SetY1NDC(0.75);
    fake_rate_stats->SetX2NDC(0.99);
    fake_rate_stats->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats1 = (TPaveStats*)(wjets_helePt_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats1 != 0 ) {
    fake_rate_stats1->SetX1NDC(0.8);
    fake_rate_stats1->SetY1NDC(0.5);
    fake_rate_stats1->SetX2NDC(0.99);
    fake_rate_stats1->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);
  elept_em->Update();


   TLatex *   lable_ref = new TLatex(0.55,0.75,"");
   lable_ref->SetNDC();
   lable_ref->SetTextSize(0.04);
   lable_ref->SetTextColor(kRed);
   TLatex *   lable_new = new TLatex(0.55,0.7,"");
   lable_new->SetNDC();
   lable_new->SetTextSize(0.04);
   lable_new->SetTextColor(kBlue);
   lable_new->Draw();
   lable_ref->Draw();

  // Plot eta comparison
  eleeta_em->cd();
  TH1F* wjets_heleEta_em_observed = (TH1F*) ((_file0->Get("wjets_elEta_em"))->Clone("wjets_heleEta_em_observed"));
  wjets_heleEta_em_observed->Rebin(rebinvalue_eta);
  wjets_heleEta_em_observed->SetLineColor(kRed);
  wjets_heleEta_em_observed->SetFillColor(kWhite);
  wjets_heleEta_em_observed->SetLineWidth(2.);
  wjets_heleEta_em_observed->SetName("observed");
  wjets_heleEta_em_observed->GetYaxis()->SetTitle("Events");
  wjets_heleEta_em_observed->GetYaxis()->SetTitleOffset(1.2);
  wjets_heleEta_em_observed->GetYaxis()->SetTitleSize(0.04);
  wjets_heleEta_em_observed->GetXaxis()->SetTitle("#eta^{electron}");
  wjets_heleEta_em_observed->GetXaxis()->SetTitleOffset(1.2);
  wjets_heleEta_em_observed->GetXaxis()->SetTitleSize(0.04);
  wjets_heleEta_em_observed->SetMarkerStyle(20);
  wjets_heleEta_em_observed->SetMarkerColor(kRed);
  wjets_heleEta_em_observed->SetMarkerSize(1.1);
  wjets_heleEta_em_observed->Draw();


  TH1F* wjets_heleEta_em_predicted = (TH1F*) ((_file2->Get("wjets_elEta_em"))->Clone("wjets_heleEta_em_predicted"));
  wjets_heleEta_em_predicted->Rebin(rebinvalue_eta);
  wjets_heleEta_em_predicted->SetLineColor(kBlue);
  wjets_heleEta_em_predicted->SetFillColor(kWhite);
  wjets_heleEta_em_predicted->SetLineWidth(2.);
  wjets_heleEta_em_predicted->SetName("predicted");
  wjets_heleEta_em_predicted->SetMarkerStyle(28);
  wjets_heleEta_em_predicted->SetMarkerColor(kBlue);
  wjets_heleEta_em_predicted->SetMarkerSize(1.2);
  wjets_heleEta_em_predicted->Draw("sames");

  eleeta_em->Update();
//    lable_new->Draw();
//    lable_ref->Draw();
  leteleEta = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leteleEta->SetLineColor(1);
  leteleEta->SetLineStyle(1);
  leteleEta->SetLineWidth(1);
  leteleEta->SetFillColor(10);
  leteleEta->SetBorderSize(1);
  //  leteleEta->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leteleEta->AddEntry(wjets_heleEta_em_observed,       "Observed","lpf");
  leteleEta->AddEntry(wjets_heleEta_em_predicted,      "Predicted","lpf");
  leteleEta->Draw();


  TPaveStats *fake_rate_stats5 = (TPaveStats*)(wjets_heleEta_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats5 != 0 ) {
    fake_rate_stats5->SetX1NDC(0.8);
    fake_rate_stats5->SetY1NDC(0.75);
    fake_rate_stats5->SetX2NDC(0.99);
    fake_rate_stats5->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats6 = (TPaveStats*)(wjets_heleEta_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats6 != 0 ) {
    fake_rate_stats6->SetX1NDC(0.8);
    fake_rate_stats6->SetY1NDC(0.5);
    fake_rate_stats6->SetX2NDC(0.99);
    fake_rate_stats6->SetY2NDC(0.74);
  }
  eleeta_em->Update();
  leg->Draw();

  // plot nJet comparison
  nJet_em->cd();
  nJet_em->SetLogy();
  TH1F* wjets_hnJet_em_observed = (TH1F*) ((_file0->Get("wjets_nJet_em"))->Clone("wjets_hnJet_em_observed"));
  //  wjets_hnJet_em_observed->Rebin(5);
  wjets_hnJet_em_observed->SetLineColor(kRed);
  wjets_hnJet_em_observed->SetFillColor(kWhite);
  wjets_hnJet_em_observed->SetLineWidth(2.);
  wjets_hnJet_em_observed->SetName("observed");
  wjets_hnJet_em_observed->GetYaxis()->SetTitle("Events");
  wjets_hnJet_em_observed->GetYaxis()->SetTitleOffset(1.2);
  wjets_hnJet_em_observed->GetYaxis()->SetTitleSize(0.04);
  wjets_hnJet_em_observed->GetXaxis()->SetTitle("Number of Jets");
  wjets_hnJet_em_observed->GetXaxis()->SetTitleOffset(1.2);
  wjets_hnJet_em_observed->GetXaxis()->SetTitleSize(0.04);
  wjets_hnJet_em_observed->SetMarkerStyle(20);
  wjets_hnJet_em_observed->SetMarkerColor(kRed);
  wjets_hnJet_em_observed->SetMarkerSize(1.1);
  wjets_hnJet_em_observed->Draw();

  TH1F* wjets_hnJet_em_predicted = (TH1F*) ((_file2->Get("wjets_nJet_em"))->Clone("wjets_hnJet_em_predicted"));
  //  wjets_hnJet_em_predicted->Rebin(5);
  wjets_hnJet_em_predicted->SetLineColor(kBlue);
  wjets_hnJet_em_predicted->SetFillColor(kWhite);
  wjets_hnJet_em_predicted->SetLineWidth(2.);
  wjets_hnJet_em_predicted->SetName("predicted");
  wjets_hnJet_em_predicted->SetMarkerStyle(28);
  wjets_hnJet_em_predicted->SetMarkerColor(kBlue);
  wjets_hnJet_em_predicted->SetMarkerSize(1.2);
  wjets_hnJet_em_predicted->Draw("sames");
//    lable_new->Draw();
//    lable_ref->Draw();
  nJet_em->Update();
  letnJet = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letnJet->SetLineColor(1);
  letnJet->SetLineStyle(1);
  letnJet->SetLineWidth(1);
  letnJet->SetFillColor(10);
  letnJet->SetBorderSize(1);
  //  letnJet->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letnJet->AddEntry(wjets_hnJet_em_observed,       "Observed","lpf");
  letnJet->AddEntry(wjets_hnJet_em_predicted,      "Predicted","lpf");
  letnJet->Draw();


  nJet_em->Update();
  //  nJet_em->SetNDC();
  TPaveStats *fake_rate_stats2 = (TPaveStats*)(wjets_hnJet_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2 != 0 ) {
    fake_rate_stats2->SetX1NDC(0.8);
    fake_rate_stats2->SetY1NDC(0.75);
    fake_rate_stats2->SetX2NDC(0.99);
    fake_rate_stats2->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3 = (TPaveStats*)(wjets_hnJet_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3 != 0 ) {
    fake_rate_stats3->SetX1NDC(0.8);
    fake_rate_stats3->SetY1NDC(0.5);
    fake_rate_stats3->SetX2NDC(0.99);
    fake_rate_stats3->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  nJet_em->Update();
   TLatex *   lable_ref3 = new TLatex(0.55,0.75,"");
   lable_ref3->SetNDC();
   lable_ref3->SetTextSize(0.04);
   lable_ref3->SetTextColor(kRed);
   TLatex *   lable_new3 = new TLatex(0.55,0.7,"");
   lable_new3->SetNDC();
   lable_new3->SetTextSize(0.04);
   lable_new3->SetTextColor(kBlue);
   lable_new3->Draw();
   lable_ref3->Draw();

  // plot MET comparison
  MET_em->cd();
  TH1F* wjets_hMET_em_observed = (TH1F*) ((_file0->Get("wjets_met_em"))->Clone("wjets_hMET_em_observed"));
  //  wjets_hMET_em_observed->Rebin(5);
  wjets_hMET_em_observed->SetLineColor(kRed);
  wjets_hMET_em_observed->SetFillColor(kWhite);
  wjets_hMET_em_observed->SetLineWidth(2.);
  wjets_hMET_em_observed->SetName("observed");
  wjets_hMET_em_observed->GetYaxis()->SetTitle("Events");
  wjets_hMET_em_observed->GetYaxis()->SetTitleOffset(1.2);
  wjets_hMET_em_observed->GetYaxis()->SetTitleSize(0.04);
  wjets_hMET_em_observed->GetXaxis()->SetTitle("MET (GeV)");
  wjets_hMET_em_observed->GetXaxis()->SetTitleOffset(1.2);
  wjets_hMET_em_observed->GetXaxis()->SetTitleSize(0.04);
  wjets_hMET_em_observed->SetMarkerStyle(20);
  wjets_hMET_em_observed->SetMarkerColor(kRed);
  wjets_hMET_em_observed->SetMarkerSize(1.1);
  wjets_hMET_em_observed->Draw();

  TH1F* wjets_hMET_em_predicted = (TH1F*) ((_file2->Get("wjets_met_em"))->Clone("wjets_hMET_em_predicted"));
  //  wjets_hMET_em_predicted->Rebin(5);
  wjets_hMET_em_predicted->SetLineColor(kBlue);
  wjets_hMET_em_predicted->SetFillColor(kWhite);
  wjets_hMET_em_predicted->SetLineWidth(2.);
  wjets_hMET_em_predicted->SetName("predicted");
  wjets_hMET_em_predicted->SetMarkerStyle(28);
  wjets_hMET_em_predicted->SetMarkerColor(kBlue);
  wjets_hMET_em_predicted->SetMarkerSize(1.2);
  wjets_hMET_em_predicted->Draw("sames");
//    lable_new->Draw();
//    lable_ref->Draw();

  letMET = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letMET->SetLineColor(1);
  letMET->SetLineStyle(1);
  letMET->SetLineWidth(1);
  letMET->SetFillColor(10);
  letMET->SetBorderSize(1);
  //  letMET->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letMET->AddEntry(wjets_hnJet_em_observed,       "Observed","lpf");
  letMET->AddEntry(wjets_hnJet_em_predicted,      "Predicted","lpf");
  letMET->Draw();


  MET_em->Update();
  //  MET_em->SetNDC();
  TPaveStats *fake_rate_stats2 = (TPaveStats*)(wjets_hMET_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2 != 0 ) {
    fake_rate_stats2->SetX1NDC(0.8);
    fake_rate_stats2->SetY1NDC(0.75);
    fake_rate_stats2->SetX2NDC(0.99);
    fake_rate_stats2->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3 = (TPaveStats*)(wjets_hMET_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3 != 0 ) {
    fake_rate_stats3->SetX1NDC(0.8);
    fake_rate_stats3->SetY1NDC(0.5);
    fake_rate_stats3->SetX2NDC(0.99);
    fake_rate_stats3->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  MET_em->Update();
   TLatex *   lable_ref2 = new TLatex(0.55,0.75,"");
   lable_ref2->SetNDC();
   lable_ref2->SetTextSize(0.04);
   lable_ref2->SetTextColor(kRed);
   TLatex *   lable_new2 = new TLatex(0.55,0.7,"");
   lable_new2->SetNDC();
   lable_new2->SetTextSize(0.04);
   lable_new2->SetTextColor(kBlue);
   lable_new2->Draw();
   lable_ref2->Draw();


   elept_em->Print("elept_em.pdf");
   eleeta_em->Print("eleeta_em.pdf");
   nJet_em->Print("nJet_em.pdf");
   MET_em->Print("MET_em_em.pdf");
   

}
