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
  bool doSS = false;

  //  TString sample="wjetsAlpgen_";
    TString sample="wjets_";
  //  TString sample="ttbar_";
  
  TString observed = "Observed (Numerator)";
  //  TString predicted = "FO not Numerator";
  TString predicted = "Predicted";
  
  int rebinvalue_pt = 1;
  int rebinvalue_eta = 1;

  if( sample.Contains("wjets") ) {
    if(doSS) {
      TFile *_file0 = TFile::Open("Wjets_SS_Numerator.root");
      TFile *_file2 = TFile::Open("Wjets_SS_Fakerate.root");
    }
    else {
      TFile *_file0 = TFile::Open("Wjets_Numerator.root");
      TFile *_file2 = TFile::Open("Wjets_Fakerate.root");
    }
  }
  else {
    if(doSS) {
      TFile *_file0 = TFile::Open("Ttbar_SS_Numerator.root");
      TFile *_file2 = TFile::Open("Ttbar_SS_Fakerate.root");
    }
    else {
      TFile *_file0 = TFile::Open("Ttbar_Numerator.root");
      TFile *_file2 = TFile::Open("Ttbar_Fakerate.root");
    }

  }
  //   TFile *_file0 = TFile::Open("Wjets_SS_FOs_Not_Numerator.root");
  //  TFile *_file2 = TFile::Open("Wjets_FOs_Not_Numerator.root");

//   TFile *_file0 = TFile::Open("Wjets_SS_Numerator.root");
//   TFile *_file2 = TFile::Open("Wjets_Numerator.root");

  TCanvas * elept_em = new TCanvas("elept_em","elept_em");
  TCanvas * eleeta_em = new TCanvas("eleeta_em","eleeta_em");
  TCanvas * eleEtaPtBelow20_em = new TCanvas("eleEtaPtBelow20_em","eleEtaPtBelow20_em");
  TCanvas * eleEtaPtAbove20_em = new TCanvas("eleEtaPtAbove20_em","eleEtaPtAbove20_em");
  TCanvas * nJet_em = new TCanvas("nJet_em","nJet_em");
  TCanvas * celeRelIso_em = new TCanvas("celeRelIsoN1_em","celeRelIsoN1_em");
  TCanvas * MET_em = new TCanvas("MET_em","MET_em");
  TCanvas * elPdgId_em = new TCanvas("elPdgId_em","elPdgId_em");
  TCanvas * elMoPdgId_em = new TCanvas("elMoPdgId_em","elMoPdgId_em");
  TCanvas * elPdgIdCat_em = new TCanvas("elPdgIdCat_em","elPdgIdCat_em");

  //  TCanvas * elept_em_true = new TCanvas("elept_em_true","elept_em_true");

  // plot Pt comparison
  
  elept_em->cd();
  elept_em->SetLogy();
  TH1F* helePt_em_observed = (TH1F*) ((_file0->Get(sample+"elPt_em"))->Clone("helePt_em_observed"));
  helePt_em_observed->Rebin(rebinvalue_pt);
  helePt_em_observed->SetLineColor(kRed);
  helePt_em_observed->SetFillColor(kWhite);
  helePt_em_observed->SetLineWidth(2.);
  helePt_em_observed->SetName(observed);
  helePt_em_observed->GetYaxis()->SetTitle("Events");
  helePt_em_observed->GetYaxis()->SetTitleOffset(1.2);
  helePt_em_observed->GetYaxis()->SetTitleSize(0.04);
  helePt_em_observed->GetXaxis()->SetTitle("p_{T}^{electron} (GeV)");
  helePt_em_observed->GetXaxis()->SetTitleOffset(1.2);
  helePt_em_observed->GetXaxis()->SetTitleSize(0.04);
  helePt_em_observed->SetMarkerStyle(20);
  helePt_em_observed->SetMarkerColor(kRed);
  helePt_em_observed->SetMarkerSize(1.1);


  TH1F* helePt_em_predicted = (TH1F*) ((_file2->Get(sample+"elPt_em"))->Clone("helePt_em_predicted"));
  helePt_em_predicted->Rebin(rebinvalue_pt);
  helePt_em_predicted->SetLineColor(kBlue);
  helePt_em_predicted->SetFillColor(kWhite);
  helePt_em_predicted->SetLineWidth(2.);
  helePt_em_predicted->SetName(predicted);
  helePt_em_predicted->SetMarkerStyle(28);
  helePt_em_predicted->SetMarkerColor(kBlue);
  helePt_em_predicted->SetMarkerSize(1.2);
  if(helePt_em_observed->GetMaximum() >= helePt_em_predicted->GetMaximum() ) {
    helePt_em_observed->Draw();
    helePt_em_predicted->Draw("sames");
  }
  else {
    helePt_em_predicted->Draw();
    helePt_em_observed->Draw("sames");
  }

  leg = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  //  leg->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leg->AddEntry(helePt_em_observed,       observed,"lpf");
  leg->AddEntry(helePt_em_predicted,      predicted,"lpf");
  leg->Draw();


  elept_em->Update();

  TPaveStats *fake_rate_stats = (TPaveStats*)(helePt_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats != 0 ) {
    fake_rate_stats->SetX1NDC(0.8);
    fake_rate_stats->SetY1NDC(0.75);
    fake_rate_stats->SetX2NDC(0.99);
    fake_rate_stats->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats1 = (TPaveStats*)(helePt_em_predicted->GetListOfFunctions()->FindObject("stats"));
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
  TH1F* heleEta_em_observed = (TH1F*) ((_file0->Get(sample+"elEta_em"))->Clone("heleEta_em_observed"));
  heleEta_em_observed->Rebin(rebinvalue_eta);
  heleEta_em_observed->SetLineColor(kRed);
  heleEta_em_observed->SetFillColor(kWhite);
  heleEta_em_observed->SetLineWidth(2.);
  heleEta_em_observed->SetName(observed);
  heleEta_em_observed->GetYaxis()->SetTitle("Events");
  heleEta_em_observed->GetYaxis()->SetTitleOffset(1.2);
  heleEta_em_observed->GetYaxis()->SetTitleSize(0.04);
  heleEta_em_observed->GetXaxis()->SetTitle("#eta^{electron}");
  heleEta_em_observed->GetXaxis()->SetTitleOffset(1.2);
  heleEta_em_observed->GetXaxis()->SetTitleSize(0.04);
  heleEta_em_observed->SetMarkerStyle(20);
  heleEta_em_observed->SetMarkerColor(kRed);
  heleEta_em_observed->SetMarkerSize(1.1);
  //  heleEta_em_observed->Draw();

  TH1F* heleEta_em_predicted = (TH1F*) ((_file2->Get(sample+"elEta_em"))->Clone("heleEta_em_predicted"));
  heleEta_em_predicted->Rebin(rebinvalue_eta);
  heleEta_em_predicted->SetLineColor(kBlue);
  heleEta_em_predicted->SetFillColor(kWhite);
  heleEta_em_predicted->SetLineWidth(2.);
  heleEta_em_predicted->SetName(predicted);
  heleEta_em_predicted->SetMarkerStyle(28);
  heleEta_em_predicted->SetMarkerColor(kBlue);
  heleEta_em_predicted->SetMarkerSize(1.2);
  //  heleEta_em_predicted->Draw("sames");
  if(heleEta_em_observed->GetMaximum() >= heleEta_em_predicted->GetMaximum() ) {
    heleEta_em_observed->Draw();
    heleEta_em_predicted->Draw("sames");
  }
  else {
    heleEta_em_predicted->Draw();
    heleEta_em_observed->Draw("sames");
  }

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
  leteleEta->AddEntry(heleEta_em_observed,       observed,"lpf");
  leteleEta->AddEntry(heleEta_em_predicted,      predicted,"lpf");
  leteleEta->Draw();


  TPaveStats *fake_rate_stats5 = (TPaveStats*)(heleEta_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats5 != 0 ) {
    fake_rate_stats5->SetX1NDC(0.8);
    fake_rate_stats5->SetY1NDC(0.75);
    fake_rate_stats5->SetX2NDC(0.99);
    fake_rate_stats5->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats6 = (TPaveStats*)(heleEta_em_predicted->GetListOfFunctions()->FindObject("stats"));
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
  TH1F* hnJet_em_observed = (TH1F*) ((_file0->Get(sample+"nJet_em"))->Clone("hnJet_em_observed"));
  //  hnJet_em_observed->Rebin(5);
  hnJet_em_observed->SetLineColor(kRed);
  hnJet_em_observed->SetFillColor(kWhite);
  hnJet_em_observed->SetLineWidth(2.);
  hnJet_em_observed->SetName(observed);
  hnJet_em_observed->GetYaxis()->SetTitle("Events");
  hnJet_em_observed->GetYaxis()->SetTitleOffset(1.2);
  hnJet_em_observed->GetYaxis()->SetTitleSize(0.04);
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
  hnJet_em_predicted->SetName(predicted);
  hnJet_em_predicted->SetMarkerStyle(28);
  hnJet_em_predicted->SetMarkerColor(kBlue);
  hnJet_em_predicted->SetMarkerSize(1.2);
  //  hnJet_em_predicted->Draw("sames");
  if(hnJet_em_observed->GetMaximum() >= hnJet_em_predicted->GetMaximum() ) {
    hnJet_em_observed->Draw();
    hnJet_em_predicted->Draw("sames");
  }
  else {
    hnJet_em_predicted->Draw();
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
  // plot elPdgId comparison
  elPdgId_em->cd();
  //  elPdgId_em->SetLogy();
  TH1F* helPdgId_em_observed = (TH1F*) ((_file0->Get(sample+"elPdgId_em"))->Clone("helPdgId_em_observed"));
  //  helPdgId_em_observed->Rebin(5);
  helPdgId_em_observed->SetLineColor(kRed);
  helPdgId_em_observed->SetFillColor(kWhite);
  helPdgId_em_observed->SetLineWidth(2.);
  helPdgId_em_observed->SetName(observed);
  helPdgId_em_observed->GetYaxis()->SetTitle("Events");
  helPdgId_em_observed->GetYaxis()->SetTitleOffset(1.2);
  helPdgId_em_observed->GetYaxis()->SetTitleSize(0.04);
  helPdgId_em_observed->GetXaxis()->SetTitle("electron Pdg ID");
  helPdgId_em_observed->GetXaxis()->SetTitleOffset(1.2);
  helPdgId_em_observed->GetXaxis()->SetTitleSize(0.04);
  helPdgId_em_observed->SetMarkerStyle(20);
  helPdgId_em_observed->SetMarkerColor(kRed);
  helPdgId_em_observed->SetMarkerSize(1.1);
  //  helPdgId_em_observed->Draw();

  TH1F* helPdgId_em_predicted = (TH1F*) ((_file2->Get(sample+"elPdgId_em"))->Clone("helPdgId_em_predicted"));
  //  helPdgId_em_predicted->Rebin(5);
  helPdgId_em_predicted->SetLineColor(kBlue);
  helPdgId_em_predicted->SetFillColor(kWhite);
  helPdgId_em_predicted->SetLineWidth(2.);
  helPdgId_em_predicted->SetName(predicted);
  helPdgId_em_predicted->SetMarkerStyle(28);
  helPdgId_em_predicted->SetMarkerColor(kBlue);
  helPdgId_em_predicted->SetMarkerSize(1.2);
  //  helPdgId_em_predicted->Draw("sames");
  if(helPdgId_em_observed->GetMaximum() >= helPdgId_em_predicted->GetMaximum() ) {
    helPdgId_em_observed->Draw();
    helPdgId_em_predicted->Draw("sames");
  }
  else {
    helPdgId_em_predicted->Draw();
    helPdgId_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elPdgId_em->Update();
  letelPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelPdgId->SetLineColor(1);
  letelPdgId->SetLineStyle(1);
  letelPdgId->SetLineWidth(1);
  letelPdgId->SetFillColor(10);
  letelPdgId->SetBorderSize(1);
  //  letelPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelPdgId->AddEntry(helPdgId_em_observed,       observed,"lpf");
  letelPdgId->AddEntry(helPdgId_em_predicted,      predicted,"lpf");
  letelPdgId->Draw();


  elPdgId_em->Update();
  //  elPdgId_em->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helPdgId_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elPdgIdstats2 = (TPaveStats*)(helPdgId_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elPdgIdstats2 != 0 ) {
    fake_rate_elPdgIdstats2->SetX1NDC(0.8);
    fake_rate_elPdgIdstats2->SetY1NDC(0.5);
    fake_rate_elPdgIdstats2->SetX2NDC(0.99);
    fake_rate_elPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elPdgId_em->Update();
   TLatex *   lable_refelPdgId = new TLatex(0.55,0.75,"");
   lable_refelPdgId->SetNDC();
   lable_refelPdgId->SetTextSize(0.04);
   lable_refelPdgId->SetTextColor(kRed);
   TLatex *   lable_newelPdgId = new TLatex(0.55,0.7,"");
   lable_newelPdgId->SetNDC();
   lable_newelPdgId->SetTextSize(0.04);
   lable_newelPdgId->SetTextColor(kBlue);
   lable_newelPdgId->Draw();
   lable_refelPdgId->Draw();

   //EEEE
  // plot elMoPdgId comparison
  elMoPdgId_em->cd();
  //  elMoPdgId_em->SetLogy();
  TH1F* helMoPdgId_em_observed = (TH1F*) ((_file0->Get(sample+"elMoPdgId_em"))->Clone("helMoPdgId_em_observed"));
  //  helMoPdgId_em_observed->Rebin(5);
  helMoPdgId_em_observed->SetLineColor(kRed);
  helMoPdgId_em_observed->SetFillColor(kWhite);
  helMoPdgId_em_observed->SetLineWidth(2.);
  helMoPdgId_em_observed->SetName(observed);
  helMoPdgId_em_observed->GetYaxis()->SetTitle("Events");
  helMoPdgId_em_observed->GetYaxis()->SetTitleOffset(1.2);
  helMoPdgId_em_observed->GetYaxis()->SetTitleSize(0.04);
  helMoPdgId_em_observed->GetXaxis()->SetTitle("electron Mother Pdg ID");
  helMoPdgId_em_observed->GetXaxis()->SetTitleOffset(1.2);
  helMoPdgId_em_observed->GetXaxis()->SetTitleSize(0.04);
  helMoPdgId_em_observed->SetMarkerStyle(20);
  helMoPdgId_em_observed->SetMarkerColor(kRed);
  helMoPdgId_em_observed->SetMarkerSize(1.1);
  //  helMoPdgId_em_observed->Draw();

  TH1F* helMoPdgId_em_predicted = (TH1F*) ((_file2->Get(sample+"elMoPdgId_em"))->Clone("helMoPdgId_em_predicted"));
  //  helMoPdgId_em_predicted->Rebin(5);
  helMoPdgId_em_predicted->SetLineColor(kBlue);
  helMoPdgId_em_predicted->SetFillColor(kWhite);
  helMoPdgId_em_predicted->SetLineWidth(2.);
  helMoPdgId_em_predicted->SetName(predicted);
  helMoPdgId_em_predicted->SetMarkerStyle(28);
  helMoPdgId_em_predicted->SetMarkerColor(kBlue);
  helMoPdgId_em_predicted->SetMarkerSize(1.2);
  //  helMoPdgId_em_predicted->Draw("sames");
  if(helMoPdgId_em_observed->GetMaximum() >= helMoPdgId_em_predicted->GetMaximum() ) {
    helMoPdgId_em_observed->Draw();
    helMoPdgId_em_predicted->Draw("sames");
  }
  else {
    helMoPdgId_em_predicted->Draw();
    helMoPdgId_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elMoPdgId_em->Update();
  letelMoPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelMoPdgId->SetLineColor(1);
  letelMoPdgId->SetLineStyle(1);
  letelMoPdgId->SetLineWidth(1);
  letelMoPdgId->SetFillColor(10);
  letelMoPdgId->SetBorderSize(1);
  //  letelMoPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelMoPdgId->AddEntry(helMoPdgId_em_observed,       observed,"lpf");
  letelMoPdgId->AddEntry(helMoPdgId_em_predicted,      predicted,"lpf");
  letelMoPdgId->Draw();


  elMoPdgId_em->Update();
  //  elMoPdgId_em->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helMoPdgId_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elMoPdgIdstats2 = (TPaveStats*)(helMoPdgId_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elMoPdgIdstats2 != 0 ) {
    fake_rate_elMoPdgIdstats2->SetX1NDC(0.8);
    fake_rate_elMoPdgIdstats2->SetY1NDC(0.5);
    fake_rate_elMoPdgIdstats2->SetX2NDC(0.99);
    fake_rate_elMoPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elMoPdgId_em->Update();
   TLatex *   lable_refelMoPdgId = new TLatex(0.55,0.75,"");
   lable_refelMoPdgId->SetNDC();
   lable_refelMoPdgId->SetTextSize(0.04);
   lable_refelMoPdgId->SetTextColor(kRed);
   TLatex *   lable_newelMoPdgId = new TLatex(0.55,0.7,"");
   lable_newelMoPdgId->SetNDC();
   lable_newelMoPdgId->SetTextSize(0.04);
   lable_newelMoPdgId->SetTextColor(kBlue);
   lable_newelMoPdgId->Draw();
   lable_refelMoPdgId->Draw();


   //EEE222
   if(42 == 42) {
  // plot elPdgIdCat comparison
  elPdgIdCat_em->cd();
  //  elPdgIdCat_em->SetLogy();
  TH1F* helPdgIdCat_em_observed = (TH1F*) ((_file0->Get(sample+"elPdgIdCat_em"))->Clone("helPdgIdCat_em_observed"));
  //  helPdgIdCat_em_observed->Rebin(5);
  helPdgIdCat_em_observed->SetLineColor(kRed);
  helPdgIdCat_em_observed->SetFillColor(kWhite);
  helPdgIdCat_em_observed->SetLineWidth(2.);
  helPdgIdCat_em_observed->SetName(observed);
  helPdgIdCat_em_observed->GetYaxis()->SetTitle("Events");
  helPdgIdCat_em_observed->GetYaxis()->SetTitleOffset(1.2);
  helPdgIdCat_em_observed->GetYaxis()->SetTitleSize(0.04);
  helPdgIdCat_em_observed->GetXaxis()->SetTitle("electron Pdg ID");
  helPdgIdCat_em_observed->GetXaxis()->SetTitleOffset(1.2);
  helPdgIdCat_em_observed->GetXaxis()->SetTitleSize(0.04);
  helPdgIdCat_em_observed->SetMarkerStyle(20);
  helPdgIdCat_em_observed->SetMarkerColor(kRed);
  helPdgIdCat_em_observed->SetMarkerSize(1.1);
  //  helPdgIdCat_em_observed->Draw();

  TH1F* helPdgIdCat_em_predicted = (TH1F*) ((_file2->Get(sample+"elPdgIdCat_em"))->Clone("helPdgIdCat_em_predicted"));
  //  helPdgIdCat_em_predicted->Rebin(5);
  helPdgIdCat_em_predicted->SetLineColor(kBlue);
  helPdgIdCat_em_predicted->SetFillColor(kWhite);
  helPdgIdCat_em_predicted->SetLineWidth(2.);
  helPdgIdCat_em_predicted->SetName(predicted);
  helPdgIdCat_em_predicted->SetMarkerStyle(28);
  helPdgIdCat_em_predicted->SetMarkerColor(kBlue);
  helPdgIdCat_em_predicted->SetMarkerSize(1.2);
  //  helPdgIdCat_em_predicted->Draw("sames");
  if(helPdgIdCat_em_observed->GetMaximum() >= helPdgIdCat_em_predicted->GetMaximum() ) {
    helPdgIdCat_em_observed->Draw();
    helPdgIdCat_em_predicted->Draw("sames");
  }
  else {
    helPdgIdCat_em_predicted->Draw();
    helPdgIdCat_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elPdgIdCat_em->Update();
  letelPdgIdCat = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelPdgIdCat->SetLineColor(1);
  letelPdgIdCat->SetLineStyle(1);
  letelPdgIdCat->SetLineWidth(1);
  letelPdgIdCat->SetFillColor(10);
  letelPdgIdCat->SetBorderSize(1);
  //  letelPdgIdCat->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelPdgIdCat->AddEntry(helPdgIdCat_em_observed,       observed,"lpf");
  letelPdgIdCat->AddEntry(helPdgIdCat_em_predicted,      predicted,"lpf");
  letelPdgIdCat->Draw();


  elPdgIdCat_em->Update();
  //  elPdgIdCat_em->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helPdgIdCat_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elPdgIdCatstats2 = (TPaveStats*)(helPdgIdCat_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elPdgIdCatstats2 != 0 ) {
    fake_rate_elPdgIdCatstats2->SetX1NDC(0.8);
    fake_rate_elPdgIdCatstats2->SetY1NDC(0.5);
    fake_rate_elPdgIdCatstats2->SetX2NDC(0.99);
    fake_rate_elPdgIdCatstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elPdgIdCat_em->Update();
   TLatex *   lable_refelPdgIdCat = new TLatex(0.55,0.75,"");
   lable_refelPdgIdCat->SetNDC();
   lable_refelPdgIdCat->SetTextSize(0.04);
   lable_refelPdgIdCat->SetTextColor(kRed);
   TLatex *   lable_newelPdgIdCat = new TLatex(0.55,0.7,"");
   lable_newelPdgIdCat->SetNDC();
   lable_newelPdgIdCat->SetTextSize(0.04);
   lable_newelPdgIdCat->SetTextColor(kBlue);
   lable_newelPdgIdCat->Draw();
   lable_refelPdgIdCat->Draw();
   }
   //IBLPasteB
  // plot eleRelIso comparison
  celeRelIso_em->cd();
  //  celeRelIso_em->SetLogy();
  TH1F* heleRelIso_em_observed = (TH1F*) ((_file0->Get(sample+"eleRelIso_em"))->Clone("heleRelIso_em_observed"));
  //  heleRelIso_em_observed->Rebin(5);
  heleRelIso_em_observed->SetLineColor(kRed);
  heleRelIso_em_observed->SetFillColor(kWhite);
  heleRelIso_em_observed->SetLineWidth(2.);
  heleRelIso_em_observed->SetName(observed);
  heleRelIso_em_observed->GetYaxis()->SetTitle("Events");
  heleRelIso_em_observed->GetYaxis()->SetTitleOffset(1.2);
  heleRelIso_em_observed->GetYaxis()->SetTitleSize(0.04);
  heleRelIso_em_observed->GetXaxis()->SetTitle("Number of Jets");
  heleRelIso_em_observed->GetXaxis()->SetTitleOffset(1.2);
  heleRelIso_em_observed->GetXaxis()->SetTitleSize(0.04);
  heleRelIso_em_observed->SetMarkerStyle(20);
  heleRelIso_em_observed->SetMarkerColor(kRed);
  heleRelIso_em_observed->SetMarkerSize(1.1);
  //  heleRelIso_em_observed->Draw();

  TH1F* heleRelIso_em_predicted = (TH1F*) ((_file2->Get(sample+"eleRelIso_em"))->Clone("heleRelIso_em_predicted"));
  //  heleRelIso_em_predicted->Rebin(5);
  heleRelIso_em_predicted->SetLineColor(kBlue);
  heleRelIso_em_predicted->SetFillColor(kWhite);
  heleRelIso_em_predicted->SetLineWidth(2.);
  heleRelIso_em_predicted->SetName(predicted);
  heleRelIso_em_predicted->SetMarkerStyle(28);
  heleRelIso_em_predicted->SetMarkerColor(kBlue);
  heleRelIso_em_predicted->SetMarkerSize(1.2);
  //  heleRelIso_em_predicted->Draw("sames");
  if(heleRelIso_em_observed->GetMaximum() >= heleRelIso_em_predicted->GetMaximum() ) {
    heleRelIso_em_observed->Draw();
    heleRelIso_em_predicted->Draw("sames");
  }
  else {
    heleRelIso_em_predicted->Draw();
    heleRelIso_em_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  celeRelIso_em->Update();
  leteleRelIso = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leteleRelIso->SetLineColor(1);
  leteleRelIso->SetLineStyle(1);
  leteleRelIso->SetLineWidth(1);
  leteleRelIso->SetFillColor(10);
  leteleRelIso->SetBorderSize(1);
  //  leteleRelIso->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leteleRelIso->AddEntry(heleRelIso_em_observed,       observed,"lpf");
  leteleRelIso->AddEntry(heleRelIso_em_predicted,      predicted,"lpf");
  leteleRelIso->Draw();


  celeRelIso_em->Update();
  //  eleRelIso_em->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(heleRelIso_em_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(heleRelIso_em_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  celeRelIso_em->Update();
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

//    elept_em->Print("elept_em.pdf");
//    eleeta_em->Print("eleeta_em.pdf");
//    nJet_em->Print("nJet_em.pdf");
//    MET_em->Print("MET_em_em.pdf");

   elept_em->Print(sample+"elept_em.png");
   eleeta_em->Print(sample+"eleeta_em.png");
//    eleEtaPtBelow20_em->Print("eleEtaPtBelow20_em.png");
//    eleEtaPtAbove20_em->Print("eleEtaPtAbove20_em.png");
   nJet_em->Print(sample+"nJet_em.png");
   celeRelIso_em->Print(sample+"eleRelIso_em.png");
   elMoPdgId_em->Print(sample+"elMoPdgId_em.png");
   elPdgId_em->Print(sample+"elPdgId_em.png");
   elPdgIdCat_em->Print(sample+"elPdgIdCat_em.png");
   //   MET_em->Print("MET_em_em.png");
   

}
//   //IBLB
//   // Plot eta comparison
//   eleEtaPtBelow20_em->cd();
//   TH1F* heleEtaPtBelow20_em_observed = (TH1F*) ((_file0->Get(sample+"elEtaPtBelow20_em"))->Clone("heleEtaPtBelow20_em_observed"));
//   heleEtaPtBelow20_em_observed->Rebin(rebinvalue_eta);
//   heleEtaPtBelow20_em_observed->SetLineColor(kRed);
//   heleEtaPtBelow20_em_observed->SetFillColor(kWhite);
//   heleEtaPtBelow20_em_observed->SetLineWidth(2.);
//   heleEtaPtBelow20_em_observed->SetName(observed);
//   heleEtaPtBelow20_em_observed->GetYaxis()->SetTitle("Events");
//   heleEtaPtBelow20_em_observed->GetYaxis()->SetTitleOffset(1.2);
//   heleEtaPtBelow20_em_observed->GetYaxis()->SetTitleSize(0.04);
//   heleEtaPtBelow20_em_observed->GetXaxis()->SetTitle("#eta^{electron}");
//   heleEtaPtBelow20_em_observed->GetXaxis()->SetTitleOffset(1.2);
//   heleEtaPtBelow20_em_observed->GetXaxis()->SetTitleSize(0.04);
//   heleEtaPtBelow20_em_observed->SetMarkerStyle(20);
//   heleEtaPtBelow20_em_observed->SetMarkerColor(kRed);
//   heleEtaPtBelow20_em_observed->SetMarkerSize(1.1);
//   heleEtaPtBelow20_em_observed->Draw();

//   TH1F* heleEtaPtBelow20_em_predicted = (TH1F*) ((_file2->Get(sample+"elEtaPtBelow20_em"))->Clone("heleEtaPtBelow20_em_predicted"));
//   heleEtaPtBelow20_em_predicted->Rebin(rebinvalue_eta);
//   heleEtaPtBelow20_em_predicted->SetLineColor(kBlue);
//   heleEtaPtBelow20_em_predicted->SetFillColor(kWhite);
//   heleEtaPtBelow20_em_predicted->SetLineWidth(2.);
//   heleEtaPtBelow20_em_predicted->SetName(predicted);
//   heleEtaPtBelow20_em_predicted->SetMarkerStyle(28);
//   heleEtaPtBelow20_em_predicted->SetMarkerColor(kBlue);
//   heleEtaPtBelow20_em_predicted->SetMarkerSize(1.2);
//   heleEtaPtBelow20_em_predicted->Draw("sames");

//   eleEtaPtBelow20_em->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   leteleEtaPtBelow20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   leteleEtaPtBelow20->SetLineColor(1);
//   leteleEtaPtBelow20->SetLineStyle(1);
//   leteleEtaPtBelow20->SetLineWidth(1);
//   leteleEtaPtBelow20->SetFillColor(10);
//   leteleEtaPtBelow20->SetBorderSize(1);
//   //  leteleEtaPtBelow20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   leteleEtaPtBelow20->AddEntry(heleEtaPtBelow20_em_observed,       observed,"lpf");
//   leteleEtaPtBelow20->AddEntry(heleEtaPtBelow20_em_predicted,      predicted,"lpf");
//   leteleEtaPtBelow20->Draw();


//   TPaveStats *fake_rate_stats5b = (TPaveStats*)(heleEtaPtBelow20_em_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5b != 0 ) {
//     fake_rate_stats5b->SetX1NDC(0.8);
//     fake_rate_stats5b->SetY1NDC(0.75);
//     fake_rate_stats5b->SetX2NDC(0.99);
//     fake_rate_stats5b->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6b = (TPaveStats*)(heleEtaPtBelow20_em_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6b != 0 ) {
//     fake_rate_stats6b->SetX1NDC(0.8);
//     fake_rate_stats6b->SetY1NDC(0.5);
//     fake_rate_stats6b->SetX2NDC(0.99);
//     fake_rate_stats6b->SetY2NDC(0.74);
//   }
//   eleEtaPtBelow20_em->Update();
//   leg->Draw();
//   //IBLE
//   // Plot eta comparison
//   eleEtaPtAbove20_em->cd();
//   TH1F* heleEtaPtAbove20_em_observed = (TH1F*) ((_file0->Get(sample+"elEtaPtAbove20_em"))->Clone("heleEtaPtAbove20_em_observed"));
//   heleEtaPtAbove20_em_observed->Rebin(rebinvalue_eta);
//   heleEtaPtAbove20_em_observed->SetLineColor(kRed);
//   heleEtaPtAbove20_em_observed->SetFillColor(kWhite);
//   heleEtaPtAbove20_em_observed->SetLineWidth(2.);
//   heleEtaPtAbove20_em_observed->SetName(observed);
//   heleEtaPtAbove20_em_observed->GetYaxis()->SetTitle("Events");
//   heleEtaPtAbove20_em_observed->GetYaxis()->SetTitleOffset(1.2);
//   heleEtaPtAbove20_em_observed->GetYaxis()->SetTitleSize(0.04);
//   heleEtaPtAbove20_em_observed->GetXaxis()->SetTitle("#eta^{electron}");
//   heleEtaPtAbove20_em_observed->GetXaxis()->SetTitleOffset(1.2);
//   heleEtaPtAbove20_em_observed->GetXaxis()->SetTitleSize(0.04);
//   heleEtaPtAbove20_em_observed->SetMarkerStyle(20);
//   heleEtaPtAbove20_em_observed->SetMarkerColor(kRed);
//   heleEtaPtAbove20_em_observed->SetMarkerSize(1.1);
//   heleEtaPtAbove20_em_observed->Draw();

//   TH1F* heleEtaPtAbove20_em_predicted = (TH1F*) ((_file2->Get(sample+"elEtaPtAbove20_em"))->Clone("heleEtaPtAbove20_em_predicted"));
//   heleEtaPtAbove20_em_predicted->Rebin(rebinvalue_eta);
//   heleEtaPtAbove20_em_predicted->SetLineColor(kBlue);
//   heleEtaPtAbove20_em_predicted->SetFillColor(kWhite);
//   heleEtaPtAbove20_em_predicted->SetLineWidth(2.);
//   heleEtaPtAbove20_em_predicted->SetName(predicted);
//   heleEtaPtAbove20_em_predicted->SetMarkerStyle(28);
//   heleEtaPtAbove20_em_predicted->SetMarkerColor(kBlue);
//   heleEtaPtAbove20_em_predicted->SetMarkerSize(1.2);
//   heleEtaPtAbove20_em_predicted->Draw("sames");

//   eleEtaPtAbove20_em->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   leteleEtaPtAbove20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   leteleEtaPtAbove20->SetLineColor(1);
//   leteleEtaPtAbove20->SetLineStyle(1);
//   leteleEtaPtAbove20->SetLineWidth(1);
//   leteleEtaPtAbove20->SetFillColor(10);
//   leteleEtaPtAbove20->SetBorderSize(1);
//   //  leteleEtaPtAbove20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   leteleEtaPtAbove20->AddEntry(heleEtaPtAbove20_em_observed,       observed,"lpf");
//   leteleEtaPtAbove20->AddEntry(heleEtaPtAbove20_em_predicted,      predicted,"lpf");
//   leteleEtaPtAbove20->Draw();


//   TPaveStats *fake_rate_stats5a = (TPaveStats*)(heleEtaPtAbove20_em_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5a != 0 ) {
//     fake_rate_stats5a->SetX1NDC(0.8);
//     fake_rate_stats5a->SetY1NDC(0.75);
//     fake_rate_stats5a->SetX2NDC(0.99);
//     fake_rate_stats5a->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6a = (TPaveStats*)(heleEtaPtAbove20_em_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6a != 0 ) {
//     fake_rate_stats6a->SetX1NDC(0.8);
//     fake_rate_stats6a->SetY1NDC(0.5);
//     fake_rate_stats6a->SetX2NDC(0.99);
//     fake_rate_stats6a->SetY2NDC(0.74);
//   }
//   eleEtaPtAbove20_em->Update();
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

