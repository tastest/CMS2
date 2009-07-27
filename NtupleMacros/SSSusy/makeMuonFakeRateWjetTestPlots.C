{
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(1111111);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1,0);

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);

  gROOT->UseCurrentStyle();

  bool subtract_realE = false;

  // use SameSign?
  bool doSS = true;

  //  TString sample="wjetsAlpgen_";
  //  TString sample="wjets_";
  TString sample="ttbar_";
  
  TString observed = "Observed (Numerator)";
  //  TString predicted = "FO not Numerator";
  TString predicted = "Predicted";
  
  int rebinvalue_pt = 1;
  int rebinvalue_eta = 1;

  if( sample.Contains("wjets") ) {
    if(doSS) {
      TFile *_file0 = TFile::Open("Wjets_SS_MuonNumerator.root");
      TFile *_file2 = TFile::Open("Wjets_SS_MuonFakerate.root");
    }
    else {
      TFile *_file0 = TFile::Open("Wjets_MuonNumerator.root");
      TFile *_file2 = TFile::Open("Wjets_MuonFakerate.root");
    }
  }
  else {
    if(doSS) {
      TFile *_file0 = TFile::Open("Ttbar_SS_MuonNumerator.root");
      TFile *_file2 = TFile::Open("Ttbar_SS_MuonFakerate.root");
    }
    else {
      TFile *_file0 = TFile::Open("Ttbar_MuonNumerator.root");
      TFile *_file2 = TFile::Open("Ttbar_MuonFakerate.root");
    }

  }
  //   TFile *_file0 = TFile::Open("Wjets_SS_FOs_Not_Numerator.root");
  //  TFile *_file2 = TFile::Open("Wjets_FOs_Not_Numerator.root");

//   TFile *_file0 = TFile::Open("Wjets_SS_Numerator.root");
//   TFile *_file2 = TFile::Open("Wjets_Numerator.root");

  TCanvas * muopt_em = new TCanvas("muopt_em","muopt_em");
  TCanvas * muoeta_em = new TCanvas("muoeta_em","muoeta_em");
  TCanvas * muoEtaPtBelow20_em = new TCanvas("muoEtaPtBelow20_em","muoEtaPtBelow20_em");
  TCanvas * muoEtaPtAbove20_em = new TCanvas("muoEtaPtAbove20_em","muoEtaPtAbove20_em");
  TCanvas * nJet_em = new TCanvas("nJet_em","nJet_em");
  TCanvas * cmuRelIso_em = new TCanvas("cmuRelIsoN1_em","cmuRelIsoN1_em");
  TCanvas * MET_em = new TCanvas("MET_em","MET_em");
  TCanvas * muPdgId_em = new TCanvas("muPdgId_em","muPdgId_em");
  TCanvas * muMoPdgId_em = new TCanvas("muMoPdgId_em","muMoPdgId_em");
  TCanvas * muPdgIdCat_em = new TCanvas("muPdgIdCat_em","muPdgIdCat_em");

  //  TCanvas * muopt_em_true = new TCanvas("muopt_em_true","muopt_em_true");

  // plot Pt comparison
  
  muopt_em->cd();
  muopt_em->SetLogy();
  TH1F* hmuoPt_em_observed = (TH1F*) ((_file0->Get(sample+"muPt_em"))->Clone("hmuoPt_em_observed"));
  hmuoPt_em_observed->Rebin(rebinvalue_pt);
  hmuoPt_em_observed->SetLineColor(kRed);
  hmuoPt_em_observed->SetFillColor(kWhite);
  hmuoPt_em_observed->SetLineWidth(2.);
  hmuoPt_em_observed->SetName(observed);
  hmuoPt_em_observed->GetYaxis()->SetTitle("Events");
  hmuoPt_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuoPt_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuoPt_em_observed->GetXaxis()->SetTitle("p_{T}^{muon} (GeV)");
  hmuoPt_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuoPt_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuoPt_em_observed->SetMarkerStyle(20);
  hmuoPt_em_observed->SetMarkerColor(kRed);
  hmuoPt_em_observed->SetMarkerSize(1.1);


  TH1F* hmuoPt_em_predicted = (TH1F*) ((_file2->Get(sample+"muPt_em"))->Clone("hmuoPt_em_predicted"));
  hmuoPt_em_predicted->Rebin(rebinvalue_pt);
  hmuoPt_em_predicted->SetLineColor(kBlue);
  hmuoPt_em_predicted->SetFillColor(kWhite);
  hmuoPt_em_predicted->SetLineWidth(2.);
  hmuoPt_em_predicted->GetYaxis()->SetTitle("Events");
  hmuoPt_em_predicted->GetXaxis()->SetTitle("p_{T}^{muon} (GeV)");
  hmuoPt_em_predicted->SetName(predicted);
  hmuoPt_em_predicted->SetMarkerStyle(28);
  hmuoPt_em_predicted->SetMarkerColor(kBlue);
  hmuoPt_em_predicted->SetMarkerSize(1.2);
  if(hmuoPt_em_observed->GetMaximum() >= hmuoPt_em_predicted->GetMaximum() ) {
    hmuoPt_em_observed->Draw();
    hmuoPt_em_predicted->Draw("hist P sames");
  }
  else {
    hmuoPt_em_predicted->Draw("hist P");
    hmuoPt_em_observed->Draw("sames");
  }

  leg = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  //  leg->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leg->AddEntry(hmuoPt_em_observed,       observed,"lpf");
  leg->AddEntry(hmuoPt_em_predicted,      predicted,"lpf");
  leg->Draw();


  muopt_em->Update();

  TPaveStats *fake_rate_stats = (TPaveStats*)(hmuoPt_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats != 0 ) {
    fake_rate_stats->SetX1NDC(0.8);
    fake_rate_stats->SetY1NDC(0.75);
    fake_rate_stats->SetX2NDC(0.99);
    fake_rate_stats->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats1 = (TPaveStats*)(hmuoPt_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats1 != 0 ) {
    fake_rate_stats1->SetX1NDC(0.8);
    fake_rate_stats1->SetY1NDC(0.5);
    fake_rate_stats1->SetX2NDC(0.99);
    fake_rate_stats1->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);
  muopt_em->Update();


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
  muoeta_em->cd();
  TH1F* hmuoEta_em_observed = (TH1F*) ((_file0->Get(sample+"muEta_em"))->Clone("hmuoEta_em_observed"));
  hmuoEta_em_observed->Rebin(rebinvalue_eta);
  hmuoEta_em_observed->SetLineColor(kRed);
  hmuoEta_em_observed->SetFillColor(kWhite);
  hmuoEta_em_observed->SetLineWidth(2.);
  hmuoEta_em_observed->SetName(observed);
  hmuoEta_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuoEta_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuoEta_em_observed->GetYaxis()->SetTitle("Events");
  hmuoEta_em_observed->GetXaxis()->SetTitle("#eta^{muon}");
  hmuoEta_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuoEta_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuoEta_em_observed->SetMarkerStyle(20);
  hmuoEta_em_observed->SetMarkerColor(kRed);
  hmuoEta_em_observed->SetMarkerSize(1.1);
  //  hmuoEta_em_observed->Draw();

  TH1F* hmuoEta_em_predicted = (TH1F*) ((_file2->Get(sample+"muEta_em"))->Clone("hmuoEta_em_predicted"));
  hmuoEta_em_predicted->Rebin(rebinvalue_eta);
  hmuoEta_em_predicted->SetLineColor(kBlue);
  hmuoEta_em_predicted->SetFillColor(kWhite);
  hmuoEta_em_predicted->SetLineWidth(2.);
  hmuoEta_em_predicted->GetYaxis()->SetTitle("Events");
  hmuoEta_em_predicted->GetXaxis()->SetTitle("#eta^{muon}");
  hmuoEta_em_predicted->SetName(predicted);
  hmuoEta_em_predicted->SetMarkerStyle(28);
  hmuoEta_em_predicted->SetMarkerColor(kBlue);
  hmuoEta_em_predicted->SetMarkerSize(1.2);
  //  hmuoEta_em_predicted->Draw("sames");
  if(hmuoEta_em_observed->GetMaximum() >= hmuoEta_em_predicted->GetMaximum() ) {
    hmuoEta_em_observed->Draw();
    hmuoEta_em_predicted->Draw("hist P sames");
  }
  else {
    hmuoEta_em_predicted->Draw("hist P");
    hmuoEta_em_observed->Draw("sames");
  }

  muoeta_em->Update();
//    lable_new->Draw();
//    lable_ref->Draw();
  letmuoEta = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmuoEta->SetLineColor(1);
  letmuoEta->SetLineStyle(1);
  letmuoEta->SetLineWidth(1);
  letmuoEta->SetFillColor(10);
  letmuoEta->SetBorderSize(1);
  //  letmuoEta->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmuoEta->AddEntry(hmuoEta_em_observed,       observed,"lpf");
  letmuoEta->AddEntry(hmuoEta_em_predicted,      predicted,"lpf");
  letmuoEta->Draw();


  TPaveStats *fake_rate_stats5 = (TPaveStats*)(hmuoEta_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats5 != 0 ) {
    fake_rate_stats5->SetX1NDC(0.8);
    fake_rate_stats5->SetY1NDC(0.75);
    fake_rate_stats5->SetX2NDC(0.99);
    fake_rate_stats5->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats6 = (TPaveStats*)(hmuoEta_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats6 != 0 ) {
    fake_rate_stats6->SetX1NDC(0.8);
    fake_rate_stats6->SetY1NDC(0.5);
    fake_rate_stats6->SetX2NDC(0.99);
    fake_rate_stats6->SetY2NDC(0.74);
  }
  muoeta_em->Update();
  leg->Draw();


  // plot nJet comparison
  nJet_em->cd();
  nJet_em->SetLogy();
  TH1F* hnJet_em_observed = (TH1F*) ((_file0->Get(sample+"nJet_em"))->Clone("hnJet_em_observed"));
  //  hnJet_em_observed->Rebin(5);
  hnJet_em_observed->SetLineColor(kRed);
  hnJet_em_observed->SetFillColor(kWhite);
  hnJet_em_observed->SetLineWidth(2.);
  hnJet_em_observed->SetName(observed);
  hnJet_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hnJet_em_observed->GetYaxis()->SetTitleSize(0.04);
  hnJet_em_observed->GetYaxis()->SetTitle("Events");
  hnJet_em_observed->GetXaxis()->SetTitle("Number of Jets");
  hnJet_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hnJet_em_observed->GetXaxis()->SetTitleSize(0.04);
  hnJet_em_observed->SetMarkerStyle(20);
  hnJet_em_observed->SetMarkerColor(kRed);
  hnJet_em_observed->SetMarkerSize(1.1);
  //  hnJet_em_observed->Draw();

  TH1F* hnJet_em_predicted = (TH1F*) ((_file2->Get(sample+"nJet_em"))->Clone("hnJet_em_predicted"));
  //  hnJet_em_predicted->Rebin(5);
  hnJet_em_predicted->SetLineColor(kBlue);
  hnJet_em_predicted->SetFillColor(kWhite);
  hnJet_em_predicted->SetLineWidth(2.);
  hnJet_em_predicted->GetYaxis()->SetTitle("Events");
  hnJet_em_predicted->GetXaxis()->SetTitle("Number of Jets");
  hnJet_em_predicted->SetName(predicted);
  hnJet_em_predicted->SetMarkerStyle(28);
  hnJet_em_predicted->SetMarkerColor(kBlue);
  hnJet_em_predicted->SetMarkerSize(1.2);
  //  hnJet_em_predicted->Draw("sames");
  if(hnJet_em_observed->GetMaximum() >= hnJet_em_predicted->GetMaximum() ) {
    hnJet_em_observed->Draw();
    hnJet_em_predicted->Draw("hist P sames");
  }
  else {
    hnJet_em_predicted->Draw("hist P");
    hnJet_em_observed->Draw("sames");
  }

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
  letnJet->AddEntry(hnJet_em_observed,       observed,"lpf");
  letnJet->AddEntry(hnJet_em_predicted,      predicted,"lpf");
  letnJet->Draw();


  nJet_em->Update();
  //  nJet_em->SetNDC();
  TPaveStats *fake_rate_stats2 = (TPaveStats*)(hnJet_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2 != 0 ) {
    fake_rate_stats2->SetX1NDC(0.8);
    fake_rate_stats2->SetY1NDC(0.75);
    fake_rate_stats2->SetX2NDC(0.99);
    fake_rate_stats2->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3 = (TPaveStats*)(hnJet_em_predicted->GetListOfFunctions()->FindObject("stats"));
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

   //BBBB
  // plot muPdgId comparison
  muPdgId_em->cd();
  //  muPdgId_em->SetLogy();
  TH1F* hmuPdgId_em_observed = (TH1F*) ((_file0->Get(sample+"muPdgId_em"))->Clone("hmuPdgId_em_observed"));
  //  hmuPdgId_em_observed->Rebin(5);
  hmuPdgId_em_observed->SetLineColor(kRed);
  hmuPdgId_em_observed->SetFillColor(kWhite);
  hmuPdgId_em_observed->SetLineWidth(2.);
  hmuPdgId_em_observed->SetName(observed);
  hmuPdgId_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuPdgId_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuPdgId_em_observed->GetYaxis()->SetTitle("Events");
  hmuPdgId_em_observed->GetXaxis()->SetTitle("muon Pdg ID");
  hmuPdgId_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuPdgId_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuPdgId_em_observed->SetMarkerStyle(20);
  hmuPdgId_em_observed->SetMarkerColor(kRed);
  hmuPdgId_em_observed->SetMarkerSize(1.1);
  //  hmuPdgId_em_observed->Draw();

  TH1F* hmuPdgId_em_predicted = (TH1F*) ((_file2->Get(sample+"muPdgId_em"))->Clone("hmuPdgId_em_predicted"));
  //  hmuPdgId_em_predicted->Rebin(5);
  hmuPdgId_em_predicted->SetLineColor(kBlue);
  hmuPdgId_em_predicted->SetFillColor(kWhite);
  hmuPdgId_em_predicted->SetLineWidth(2.);
  hmuPdgId_em_predicted->SetName(predicted);
  hmuPdgId_em_predicted->GetYaxis()->SetTitle("Events");
  hmuPdgId_em_predicted->GetXaxis()->SetTitle("muon Pdg ID");
  hmuPdgId_em_predicted->SetMarkerStyle(28);
  hmuPdgId_em_predicted->SetMarkerColor(kBlue);
  hmuPdgId_em_predicted->SetMarkerSize(1.2);
  //  hmuPdgId_em_predicted->Draw("sames");
  if(hmuPdgId_em_observed->GetMaximum() >= hmuPdgId_em_predicted->GetMaximum() ) {
    hmuPdgId_em_observed->Draw();
    hmuPdgId_em_predicted->Draw("P sames");
  }
  else {
    hmuPdgId_em_predicted->Draw("P");
    hmuPdgId_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  muPdgId_em->Update();
  letmuPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmuPdgId->SetLineColor(1);
  letmuPdgId->SetLineStyle(1);
  letmuPdgId->SetLineWidth(1);
  letmuPdgId->SetFillColor(10);
  letmuPdgId->SetBorderSize(1);
  //  letmuPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmuPdgId->AddEntry(hmuPdgId_em_observed,       observed,"lpf");
  letmuPdgId->AddEntry(hmuPdgId_em_predicted,      predicted,"lpf");
  letmuPdgId->Draw();


  muPdgId_em->Update();
  //  muPdgId_em->SetNDC();
  TPaveStats *fake_rate_muIdstats1 = (TPaveStats*)(hmuPdgId_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muIdstats1 != 0 ) {
    fake_rate_muIdstats1->SetX1NDC(0.8);
    fake_rate_muIdstats1->SetY1NDC(0.75);
    fake_rate_muIdstats1->SetX2NDC(0.99);
    fake_rate_muIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_muPdgIdstats2 = (TPaveStats*)(hmuPdgId_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muPdgIdstats2 != 0 ) {
    fake_rate_muPdgIdstats2->SetX1NDC(0.8);
    fake_rate_muPdgIdstats2->SetY1NDC(0.5);
    fake_rate_muPdgIdstats2->SetX2NDC(0.99);
    fake_rate_muPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  muPdgId_em->Update();
   TLatex *   lable_refmuPdgId = new TLatex(0.55,0.75,"");
   lable_refmuPdgId->SetNDC();
   lable_refmuPdgId->SetTextSize(0.04);
   lable_refmuPdgId->SetTextColor(kRed);
   TLatex *   lable_newmuPdgId = new TLatex(0.55,0.7,"");
   lable_newmuPdgId->SetNDC();
   lable_newmuPdgId->SetTextSize(0.04);
   lable_newmuPdgId->SetTextColor(kBlue);
   lable_newmuPdgId->Draw();
   lable_refmuPdgId->Draw();

   //EEEE
  // plot muMoPdgId comparison
  muMoPdgId_em->cd();
  //  muMoPdgId_em->SetLogy();
  TH1F* hmuMoPdgId_em_observed = (TH1F*) ((_file0->Get(sample+"muMoPdgId_em"))->Clone("hmuMoPdgId_em_observed"));
  //  hmuMoPdgId_em_observed->Rebin(5);
  hmuMoPdgId_em_observed->SetLineColor(kRed);
  hmuMoPdgId_em_observed->SetFillColor(kWhite);
  hmuMoPdgId_em_observed->SetLineWidth(2.);
  hmuMoPdgId_em_observed->SetName(observed);
  hmuMoPdgId_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuMoPdgId_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuMoPdgId_em_observed->GetYaxis()->SetTitle("Events");
  hmuMoPdgId_em_observed->GetXaxis()->SetTitle("muon Mother Pdg ID");
  hmuMoPdgId_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuMoPdgId_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuMoPdgId_em_observed->SetMarkerStyle(20);
  hmuMoPdgId_em_observed->SetMarkerColor(kRed);
  hmuMoPdgId_em_observed->SetMarkerSize(1.1);
  //  hmuMoPdgId_em_observed->Draw();

  TH1F* hmuMoPdgId_em_predicted = (TH1F*) ((_file2->Get(sample+"muMoPdgId_em"))->Clone("hmuMoPdgId_em_predicted"));
  //  hmuMoPdgId_em_predicted->Rebin(5);
  hmuMoPdgId_em_predicted->SetLineColor(kBlue);
  hmuMoPdgId_em_predicted->SetFillColor(kWhite);
  hmuMoPdgId_em_predicted->SetLineWidth(2.);
  hmuMoPdgId_em_predicted->SetName(predicted);
  hmuMoPdgId_em_predicted->GetYaxis()->SetTitle("Events");
  hmuMoPdgId_em_predicted->GetXaxis()->SetTitle("muon Mother Pdg ID");
  hmuMoPdgId_em_predicted->SetMarkerStyle(28);
  hmuMoPdgId_em_predicted->SetMarkerColor(kBlue);
  hmuMoPdgId_em_predicted->SetMarkerSize(1.2);
  //  hmuMoPdgId_em_predicted->Draw("sames");
  if(hmuMoPdgId_em_observed->GetMaximum() >= hmuMoPdgId_em_predicted->GetMaximum() ) {
    hmuMoPdgId_em_observed->Draw();
    hmuMoPdgId_em_predicted->Draw("P sames");
  }
  else {
    hmuMoPdgId_em_predicted->Draw("P");
    hmuMoPdgId_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  muMoPdgId_em->Update();
  letmuMoPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmuMoPdgId->SetLineColor(1);
  letmuMoPdgId->SetLineStyle(1);
  letmuMoPdgId->SetLineWidth(1);
  letmuMoPdgId->SetFillColor(10);
  letmuMoPdgId->SetBorderSize(1);
  //  letmuMoPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmuMoPdgId->AddEntry(hmuMoPdgId_em_observed,       observed,"lpf");
  letmuMoPdgId->AddEntry(hmuMoPdgId_em_predicted,      predicted,"lpf");
  letmuMoPdgId->Draw();


  muMoPdgId_em->Update();
  //  muMoPdgId_em->SetNDC();
  TPaveStats *fake_rate_muIdstats1 = (TPaveStats*)(hmuMoPdgId_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muIdstats1 != 0 ) {
    fake_rate_muIdstats1->SetX1NDC(0.8);
    fake_rate_muIdstats1->SetY1NDC(0.75);
    fake_rate_muIdstats1->SetX2NDC(0.99);
    fake_rate_muIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_muMoPdgIdstats2 = (TPaveStats*)(hmuMoPdgId_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muMoPdgIdstats2 != 0 ) {
    fake_rate_muMoPdgIdstats2->SetX1NDC(0.8);
    fake_rate_muMoPdgIdstats2->SetY1NDC(0.5);
    fake_rate_muMoPdgIdstats2->SetX2NDC(0.99);
    fake_rate_muMoPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  muMoPdgId_em->Update();
   TLatex *   lable_refmuMoPdgId = new TLatex(0.55,0.75,"");
   lable_refmuMoPdgId->SetNDC();
   lable_refmuMoPdgId->SetTextSize(0.04);
   lable_refmuMoPdgId->SetTextColor(kRed);
   TLatex *   lable_newmuMoPdgId = new TLatex(0.55,0.7,"");
   lable_newmuMoPdgId->SetNDC();
   lable_newmuMoPdgId->SetTextSize(0.04);
   lable_newmuMoPdgId->SetTextColor(kBlue);
   lable_newmuMoPdgId->Draw();
   lable_refmuMoPdgId->Draw();


   //EEE222
   if(42 == 42) {
  // plot muPdgIdCat comparison
  muPdgIdCat_em->cd();
  //  muPdgIdCat_em->SetLogy();
  TH1F* hmuPdgIdCat_em_observed = (TH1F*) ((_file0->Get(sample+"muPdgIdCat_em"))->Clone("hmuPdgIdCat_em_observed"));
  //  hmuPdgIdCat_em_observed->Rebin(5);
  hmuPdgIdCat_em_observed->SetLineColor(kRed);
  hmuPdgIdCat_em_observed->SetFillColor(kWhite);
  hmuPdgIdCat_em_observed->SetLineWidth(2.);
  hmuPdgIdCat_em_observed->SetName(observed);
  hmuPdgIdCat_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuPdgIdCat_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuPdgIdCat_em_observed->GetYaxis()->SetTitle("Events");
  hmuPdgIdCat_em_observed->GetXaxis()->SetTitle("muon Pdg ID Category");
  hmuPdgIdCat_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuPdgIdCat_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuPdgIdCat_em_observed->SetMarkerStyle(20);
  hmuPdgIdCat_em_observed->SetMarkerColor(kRed);
  hmuPdgIdCat_em_observed->SetMarkerSize(1.1);
  //  hmuPdgIdCat_em_observed->Draw();

  TH1F* hmuPdgIdCat_em_predicted = (TH1F*) ((_file2->Get(sample+"muPdgIdCat_em"))->Clone("hmuPdgIdCat_em_predicted"));
  //  hmuPdgIdCat_em_predicted->Rebin(5);
  hmuPdgIdCat_em_predicted->SetLineColor(kBlue);
  hmuPdgIdCat_em_predicted->SetFillColor(kWhite);
  hmuPdgIdCat_em_predicted->SetLineWidth(2.);
  hmuPdgIdCat_em_predicted->SetName(predicted);
  hmuPdgIdCat_em_predicted->GetYaxis()->SetTitle("Events");
  hmuPdgIdCat_em_predicted->GetXaxis()->SetTitle("muon Pdg ID Category");
  hmuPdgIdCat_em_predicted->SetMarkerStyle(28);
  hmuPdgIdCat_em_predicted->SetMarkerColor(kBlue);
  hmuPdgIdCat_em_predicted->SetMarkerSize(1.2);
  //  hmuPdgIdCat_em_predicted->Draw("sames");
  if(hmuPdgIdCat_em_observed->GetMaximum() >= hmuPdgIdCat_em_predicted->GetMaximum() ) {
    hmuPdgIdCat_em_observed->Draw();
    hmuPdgIdCat_em_predicted->Draw("P sames");
  }
  else {
    hmuPdgIdCat_em_predicted->Draw("P");
    hmuPdgIdCat_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  muPdgIdCat_em->Update();
  letmuPdgIdCat = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmuPdgIdCat->SetLineColor(1);
  letmuPdgIdCat->SetLineStyle(1);
  letmuPdgIdCat->SetLineWidth(1);
  letmuPdgIdCat->SetFillColor(10);
  letmuPdgIdCat->SetBorderSize(1);
  //  letmuPdgIdCat->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmuPdgIdCat->AddEntry(hmuPdgIdCat_em_observed,       observed,"lpf");
  letmuPdgIdCat->AddEntry(hmuPdgIdCat_em_predicted,      predicted,"lpf");
  letmuPdgIdCat->Draw();


  muPdgIdCat_em->Update();
  //  muPdgIdCat_em->SetNDC();
  TPaveStats *fake_rate_muIdstats1 = (TPaveStats*)(hmuPdgIdCat_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muIdstats1 != 0 ) {
    fake_rate_muIdstats1->SetX1NDC(0.8);
    fake_rate_muIdstats1->SetY1NDC(0.75);
    fake_rate_muIdstats1->SetX2NDC(0.99);
    fake_rate_muIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_muPdgIdCatstats2 = (TPaveStats*)(hmuPdgIdCat_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_muPdgIdCatstats2 != 0 ) {
    fake_rate_muPdgIdCatstats2->SetX1NDC(0.8);
    fake_rate_muPdgIdCatstats2->SetY1NDC(0.5);
    fake_rate_muPdgIdCatstats2->SetX2NDC(0.99);
    fake_rate_muPdgIdCatstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  muPdgIdCat_em->Update();
   TLatex *   lable_refmuPdgIdCat = new TLatex(0.55,0.75,"");
   lable_refmuPdgIdCat->SetNDC();
   lable_refmuPdgIdCat->SetTextSize(0.04);
   lable_refmuPdgIdCat->SetTextColor(kRed);
   TLatex *   lable_newmuPdgIdCat = new TLatex(0.55,0.7,"");
   lable_newmuPdgIdCat->SetNDC();
   lable_newmuPdgIdCat->SetTextSize(0.04);
   lable_newmuPdgIdCat->SetTextColor(kBlue);
   lable_newmuPdgIdCat->Draw();
   lable_refmuPdgIdCat->Draw();
   }
   //IBLPasteB
  // plot muRelIso comparison
  cmuRelIso_em->cd();
  //  cmuRelIso_em->SetLogy();
  TH1F* hmuRelIso_em_observed = (TH1F*) ((_file0->Get(sample+"muRelIso_em"))->Clone("hmuRelIso_em_observed"));
  //  hmuRelIso_em_observed->Rebin(5);
  hmuRelIso_em_observed->SetLineColor(kRed);
  hmuRelIso_em_observed->SetFillColor(kWhite);
  hmuRelIso_em_observed->SetLineWidth(2.);
  hmuRelIso_em_observed->SetName(observed);
  hmuRelIso_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hmuRelIso_em_observed->GetYaxis()->SetTitleSize(0.04);
  hmuRelIso_em_observed->GetYaxis()->SetTitle("Events");
  hmuRelIso_em_observed->GetXaxis()->SetTitle("rel. Iso");
  hmuRelIso_em_observed->GetXaxis()->SetTitleOffset(1.2);
  hmuRelIso_em_observed->GetXaxis()->SetTitleSize(0.04);
  hmuRelIso_em_observed->SetMarkerStyle(20);
  hmuRelIso_em_observed->SetMarkerColor(kRed);
  hmuRelIso_em_observed->SetMarkerSize(1.1);
  //  hmuRelIso_em_observed->Draw();

  TH1F* hmuRelIso_em_predicted = (TH1F*) ((_file2->Get(sample+"muRelIso_em"))->Clone("hmuRelIso_em_predicted"));
  //  hmuRelIso_em_predicted->Rebin(5);
  hmuRelIso_em_predicted->SetLineColor(kBlue);
  hmuRelIso_em_predicted->SetFillColor(kWhite);
  hmuRelIso_em_predicted->SetLineWidth(2.);
  hmuRelIso_em_predicted->SetName(predicted);
  hmuRelIso_em_predicted->GetYaxis()->SetTitle("Events");
  hmuRelIso_em_predicted->GetXaxis()->SetTitle("rel. Iso");
  hmuRelIso_em_predicted->SetMarkerStyle(28);
  hmuRelIso_em_predicted->SetMarkerColor(kBlue);
  hmuRelIso_em_predicted->SetMarkerSize(1.2);
  //  hmuRelIso_em_predicted->Draw("sames");
  if(hmuRelIso_em_observed->GetMaximum() >= hmuRelIso_em_predicted->GetMaximum() ) {
    hmuRelIso_em_observed->Draw();
    hmuRelIso_em_predicted->Draw("hist P sames");
  }
  else {
    hmuRelIso_em_predicted->Draw("hist P");
    hmuRelIso_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  cmuRelIso_em->Update();
  letmuRelIso = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmuRelIso->SetLineColor(1);
  letmuRelIso->SetLineStyle(1);
  letmuRelIso->SetLineWidth(1);
  letmuRelIso->SetFillColor(10);
  letmuRelIso->SetBorderSize(1);
  //  letmuRelIso->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmuRelIso->AddEntry(hmuRelIso_em_observed,       observed,"lpf");
  letmuRelIso->AddEntry(hmuRelIso_em_predicted,      predicted,"lpf");
  letmuRelIso->Draw();


  cmuRelIso_em->Update();
  //  muRelIso_em->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(hmuRelIso_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(hmuRelIso_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  cmuRelIso_em->Update();
   TLatex *   lable_ref3i = new TLatex(0.55,0.75,"");
   lable_ref3i->SetNDC();
   lable_ref3i->SetTextSize(0.04);
   lable_ref3i->SetTextColor(kRed);
   TLatex *   lable_new3 = new TLatex(0.55,0.7,"");
   lable_new3->SetNDC();
   lable_new3->SetTextSize(0.04);
   lable_new3->SetTextColor(kBlue);
   lable_new3->Draw();
   lable_ref3i->Draw();

   //IBLPasteE

//    muopt_em->Print("muopt_em.pdf");
//    muoeta_em->Print("muoeta_em.pdf");
//    nJet_em->Print("nJet_em.pdf");
//    MET_em->Print("MET_em_em.pdf");

   muopt_em->Print(sample+"muopt_em.png");
   muoeta_em->Print(sample+"muoeta_em.png");
//    muoEtaPtBelow20_em->Print("muoEtaPtBelow20_em.png");
//    muoEtaPtAbove20_em->Print("muoEtaPtAbove20_em.png");
   nJet_em->Print(sample+"nJet_em.png");
   cmuRelIso_em->Print(sample+"muRelIso_em.png");
   muMoPdgId_em->Print(sample+"muMoPdgId_em.png");
   muPdgId_em->Print(sample+"muPdgId_em.png");
   muPdgIdCat_em->Print(sample+"muPdgIdCat_em.png");
   //   MET_em->Print("MET_em_em.png");
   

}
//   //IBLB
//   // Plot eta comparison
//   muoEtaPtBelow20_em->cd();
//   TH1F* hmuoEtaPtBelow20_em_observed = (TH1F*) ((_file0->Get(sample+"muEtaPtBelow20_em"))->Clone("hmuoEtaPtBelow20_em_observed"));
//   hmuoEtaPtBelow20_em_observed->Rebin(rebinvalue_eta);
//   hmuoEtaPtBelow20_em_observed->SetLineColor(kRed);
//   hmuoEtaPtBelow20_em_observed->SetFillColor(kWhite);
//   hmuoEtaPtBelow20_em_observed->SetLineWidth(2.);
//   hmuoEtaPtBelow20_em_observed->SetName(observed);
//   hmuoEtaPtBelow20_em_observed->GetYaxis()->SetTitle("Events");
//   hmuoEtaPtBelow20_em_observed->GetYaxis()->SetTitleOffset(1.2);
//   hmuoEtaPtBelow20_em_observed->GetYaxis()->SetTitleSize(0.04);
//   hmuoEtaPtBelow20_em_observed->GetXaxis()->SetTitle("#eta^{muon}");
//   hmuoEtaPtBelow20_em_observed->GetXaxis()->SetTitleOffset(1.2);
//   hmuoEtaPtBelow20_em_observed->GetXaxis()->SetTitleSize(0.04);
//   hmuoEtaPtBelow20_em_observed->SetMarkerStyle(20);
//   hmuoEtaPtBelow20_em_observed->SetMarkerColor(kRed);
//   hmuoEtaPtBelow20_em_observed->SetMarkerSize(1.1);
//   hmuoEtaPtBelow20_em_observed->Draw();

//   TH1F* hmuoEtaPtBelow20_em_predicted = (TH1F*) ((_file2->Get(sample+"muEtaPtBelow20_em"))->Clone("hmuoEtaPtBelow20_em_predicted"));
//   hmuoEtaPtBelow20_em_predicted->Rebin(rebinvalue_eta);
//   hmuoEtaPtBelow20_em_predicted->SetLineColor(kBlue);
//   hmuoEtaPtBelow20_em_predicted->SetFillColor(kWhite);
//   hmuoEtaPtBelow20_em_predicted->SetLineWidth(2.);
//   hmuoEtaPtBelow20_em_predicted->SetName(predicted);
//   hmuoEtaPtBelow20_em_predicted->SetMarkerStyle(28);
//   hmuoEtaPtBelow20_em_predicted->SetMarkerColor(kBlue);
//   hmuoEtaPtBelow20_em_predicted->SetMarkerSize(1.2);
//   hmuoEtaPtBelow20_em_predicted->Draw("sames");

//   muoEtaPtBelow20_em->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   letmuoEtaPtBelow20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   letmuoEtaPtBelow20->SetLineColor(1);
//   letmuoEtaPtBelow20->SetLineStyle(1);
//   letmuoEtaPtBelow20->SetLineWidth(1);
//   letmuoEtaPtBelow20->SetFillColor(10);
//   letmuoEtaPtBelow20->SetBorderSize(1);
//   //  letmuoEtaPtBelow20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   letmuoEtaPtBelow20->AddEntry(hmuoEtaPtBelow20_em_observed,       observed,"lpf");
//   letmuoEtaPtBelow20->AddEntry(hmuoEtaPtBelow20_em_predicted,      predicted,"lpf");
//   letmuoEtaPtBelow20->Draw();


//   TPaveStats *fake_rate_stats5b = (TPaveStats*)(hmuoEtaPtBelow20_em_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5b != 0 ) {
//     fake_rate_stats5b->SetX1NDC(0.8);
//     fake_rate_stats5b->SetY1NDC(0.75);
//     fake_rate_stats5b->SetX2NDC(0.99);
//     fake_rate_stats5b->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6b = (TPaveStats*)(hmuoEtaPtBelow20_em_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6b != 0 ) {
//     fake_rate_stats6b->SetX1NDC(0.8);
//     fake_rate_stats6b->SetY1NDC(0.5);
//     fake_rate_stats6b->SetX2NDC(0.99);
//     fake_rate_stats6b->SetY2NDC(0.74);
//   }
//   muoEtaPtBelow20_em->Update();
//   leg->Draw();
//   //IBLE
//   // Plot eta comparison
//   muoEtaPtAbove20_em->cd();
//   TH1F* hmuoEtaPtAbove20_em_observed = (TH1F*) ((_file0->Get(sample+"muEtaPtAbove20_em"))->Clone("hmuoEtaPtAbove20_em_observed"));
//   hmuoEtaPtAbove20_em_observed->Rebin(rebinvalue_eta);
//   hmuoEtaPtAbove20_em_observed->SetLineColor(kRed);
//   hmuoEtaPtAbove20_em_observed->SetFillColor(kWhite);
//   hmuoEtaPtAbove20_em_observed->SetLineWidth(2.);
//   hmuoEtaPtAbove20_em_observed->SetName(observed);
//   hmuoEtaPtAbove20_em_observed->GetYaxis()->SetTitle("Events");
//   hmuoEtaPtAbove20_em_observed->GetYaxis()->SetTitleOffset(1.2);
//   hmuoEtaPtAbove20_em_observed->GetYaxis()->SetTitleSize(0.04);
//   hmuoEtaPtAbove20_em_observed->GetXaxis()->SetTitle("#eta^{muon}");
//   hmuoEtaPtAbove20_em_observed->GetXaxis()->SetTitleOffset(1.2);
//   hmuoEtaPtAbove20_em_observed->GetXaxis()->SetTitleSize(0.04);
//   hmuoEtaPtAbove20_em_observed->SetMarkerStyle(20);
//   hmuoEtaPtAbove20_em_observed->SetMarkerColor(kRed);
//   hmuoEtaPtAbove20_em_observed->SetMarkerSize(1.1);
//   hmuoEtaPtAbove20_em_observed->Draw();

//   TH1F* hmuoEtaPtAbove20_em_predicted = (TH1F*) ((_file2->Get(sample+"muEtaPtAbove20_em"))->Clone("hmuoEtaPtAbove20_em_predicted"));
//   hmuoEtaPtAbove20_em_predicted->Rebin(rebinvalue_eta);
//   hmuoEtaPtAbove20_em_predicted->SetLineColor(kBlue);
//   hmuoEtaPtAbove20_em_predicted->SetFillColor(kWhite);
//   hmuoEtaPtAbove20_em_predicted->SetLineWidth(2.);
//   hmuoEtaPtAbove20_em_predicted->SetName(predicted);
//   hmuoEtaPtAbove20_em_predicted->SetMarkerStyle(28);
//   hmuoEtaPtAbove20_em_predicted->SetMarkerColor(kBlue);
//   hmuoEtaPtAbove20_em_predicted->SetMarkerSize(1.2);
//   hmuoEtaPtAbove20_em_predicted->Draw("sames");

//   muoEtaPtAbove20_em->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   letmuoEtaPtAbove20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   letmuoEtaPtAbove20->SetLineColor(1);
//   letmuoEtaPtAbove20->SetLineStyle(1);
//   letmuoEtaPtAbove20->SetLineWidth(1);
//   letmuoEtaPtAbove20->SetFillColor(10);
//   letmuoEtaPtAbove20->SetBorderSize(1);
//   //  letmuoEtaPtAbove20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   letmuoEtaPtAbove20->AddEntry(hmuoEtaPtAbove20_em_observed,       observed,"lpf");
//   letmuoEtaPtAbove20->AddEntry(hmuoEtaPtAbove20_em_predicted,      predicted,"lpf");
//   letmuoEtaPtAbove20->Draw();


//   TPaveStats *fake_rate_stats5a = (TPaveStats*)(hmuoEtaPtAbove20_em_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5a != 0 ) {
//     fake_rate_stats5a->SetX1NDC(0.8);
//     fake_rate_stats5a->SetY1NDC(0.75);
//     fake_rate_stats5a->SetX2NDC(0.99);
//     fake_rate_stats5a->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6a = (TPaveStats*)(hmuoEtaPtAbove20_em_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6a != 0 ) {
//     fake_rate_stats6a->SetX1NDC(0.8);
//     fake_rate_stats6a->SetY1NDC(0.5);
//     fake_rate_stats6a->SetX2NDC(0.99);
//     fake_rate_stats6a->SetY2NDC(0.74);
//   }
//   muoEtaPtAbove20_em->Update();
//   leg->Draw();

  //IBLE2

//   // plot MET comparison
//   MET_em->cd();
//   TH1F* hMET_em_observed = (TH1F*) ((_file0->Get(sample+"met_em"))->Clone("hMET_em_observed"));
//   //  hMET_em_observed->Rebin(5);
//   hMET_em_observed->SetLineColor(kRed);
//   hMET_em_observed->SetFillColor(kWhite);
//   hMET_em_observed->SetLineWidth(2.);
//   hMET_em_observed->SetName(observed);
//   hMET_em_observed->GetYaxis()->SetTitle("Events");
//   hMET_em_observed->GetYaxis()->SetTitleOffset(1.2);
//   hMET_em_observed->GetYaxis()->SetTitleSize(0.04);
//   hMET_em_observed->GetXaxis()->SetTitle("MET (GeV)");
//   hMET_em_observed->GetXaxis()->SetTitleOffset(1.2);
//   hMET_em_observed->GetXaxis()->SetTitleSize(0.04);
//   hMET_em_observed->SetMarkerStyle(20);
//   hMET_em_observed->SetMarkerColor(kRed);
//   hMET_em_observed->SetMarkerSize(1.1);
//   hMET_em_observed->Draw();

//   TH1F* hMET_em_predicted = (TH1F*) ((_file2->Get(sample+"met_em"))->Clone("hMET_em_predicted"));
//   //  hMET_em_predicted->Rebin(5);
//   hMET_em_predicted->SetLineColor(kBlue);
//   hMET_em_predicted->SetFillColor(kWhite);
//   hMET_em_predicted->SetLineWidth(2.);
//   hMET_em_predicted->SetName(predicted);
//   hMET_em_predicted->SetMarkerStyle(28);
//   hMET_em_predicted->SetMarkerColor(kBlue);
//   hMET_em_predicted->SetMarkerSize(1.2);
//   hMET_em_predicted->Draw("sames");
// //    lable_new->Draw();
// //    lable_ref->Draw();

//   letMET = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   letMET->SetLineColor(1);
//   letMET->SetLineStyle(1);
//   letMET->SetLineWidth(1);
//   letMET->SetFillColor(10);
//   letMET->SetBorderSize(1);
//   //  letMET->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   letMET->AddEntry(hnJet_em_observed,       observed,"lpf");
//   letMET->AddEntry(hnJet_em_predicted,      predicted,"lpf");
//   letMET->Draw();


//   MET_em->Update();
//   //  MET_em->SetNDC();
//   TPaveStats *fake_rate_stats2 = (TPaveStats*)(hMET_em_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats2 != 0 ) {
//     fake_rate_stats2->SetX1NDC(0.8);
//     fake_rate_stats2->SetY1NDC(0.75);
//     fake_rate_stats2->SetX2NDC(0.99);
//     fake_rate_stats2->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats3 = (TPaveStats*)(hMET_em_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats3 != 0 ) {
//     fake_rate_stats3->SetX1NDC(0.8);
//     fake_rate_stats3->SetY1NDC(0.5);
//     fake_rate_stats3->SetX2NDC(0.99);
//     fake_rate_stats3->SetY2NDC(0.74);
//   }
//   //  gPad->SetLeftMargin(0.17);

//   //  leg->Draw();
//   MET_em->Update();
//    TLatex *   lable_ref2 = new TLatex(0.55,0.75,"");
//    lable_ref2->SetNDC();
//    lable_ref2->SetTextSize(0.04);
//    lable_ref2->SetTextColor(kRed);
//    TLatex *   lable_new2 = new TLatex(0.55,0.7,"");
//    lable_new2->SetNDC();
//    lable_new2->SetTextSize(0.04);
//    lable_new2->SetTextColor(kBlue);
//    lable_new2->Draw();
//    lable_ref2->Draw();

