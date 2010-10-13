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
  //  TString sample="dyeeAlpgen_";
  //    TString sample="dy20ee_";
  //  TString sample="dy20tt_";
  //  TString sample="dyee_";
  
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
      TFile *_file0 = TFile::Open("Ttbar_SS_FlipObserved.root");
      TFile *_file2 = TFile::Open("Ttbar_SS_FlipPredicted.root");
    }
    else {
//       TFile *_file0 = TFile::Open("Ttbar_Numerator.root");
//       TFile *_file2 = TFile::Open("Ttbar_Fakerate.root");
    }

  }
  //   TFile *_file0 = TFile::Open("Wjets_SS_FOs_Not_Numerator.root");
  //  TFile *_file2 = TFile::Open("Wjets_FOs_Not_Numerator.root");

//   TFile *_file0 = TFile::Open("Wjets_SS_Numerator.root");
//   TFile *_file2 = TFile::Open("Wjets_Numerator.root");

  TCanvas * elept_all = new TCanvas("elept_all","elept_all");
  TCanvas * eleeta_all = new TCanvas("eleeta_all","eleeta_all");
  TCanvas * eleEtaPtBelow20_all = new TCanvas("eleEtaPtBelow20_all","eleEtaPtBelow20_all");
  TCanvas * eleEtaPtAbove20_all = new TCanvas("eleEtaPtAbove20_all","eleEtaPtAbove20_all");
  TCanvas * nJet_all = new TCanvas("nJet_all","nJet_all");
  TCanvas * cdilMass_all = new TCanvas("cdilMass_all","cdilMass_all");
  TCanvas * celeRelIso_all = new TCanvas("celeRelIsoN1_all","celeRelIsoN1_all");
  TCanvas * elPdgId_all = new TCanvas("elPdgId_all","elPdgId_all");
  TCanvas * elMoPdgId_all = new TCanvas("elMoPdgId_all","elMoPdgId_all");
  TCanvas * elPdgIdCat_all = new TCanvas("elPdgIdCat_all","elPdgIdCat_all");
  TCanvas * cnHyp_all = new TCanvas("cnHyp_all","cnHyp_all");
  TCanvas * cmet_all = new TCanvas("cmet_all","cmet_all");

  //  TCanvas * elept_all_true = new TCanvas("elept_all_true","elept_all_true");

  // plot Pt comparison
  
  elept_all->cd();
  //  elept_all->SetLogy();
  TH1F* helePt_all_observed = (TH1F*) ((_file0->Get(sample+"elPt_all"))->Clone("helePt_all_observed"));
  helePt_all_observed->Rebin(rebinvalue_pt);
  helePt_all_observed->SetLineColor(kRed);
  helePt_all_observed->SetFillColor(kWhite);
  helePt_all_observed->SetLineWidth(2.);
  helePt_all_observed->SetName(observed);
  helePt_all_observed->GetYaxis()->SetTitle("Events");
  helePt_all_observed->GetYaxis()->SetTitleOffset(1.2);
  helePt_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  helePt_all_observed->GetXaxis()->SetTitle("p_{T}^{electron} (GeV)");
  helePt_all_observed->GetXaxis()->SetTitleOffset(1.2);
  helePt_all_observed->GetXaxis()->SetTitleSize(0.04);
  helePt_all_observed->SetMarkerStyle(20);
  helePt_all_observed->SetMarkerColor(kRed);
  helePt_all_observed->SetMarkerSize(1.1);


  TH1F* helePt_all_predicted = (TH1F*) ((_file2->Get(sample+"elPt_all"))->Clone("helePt_all_predicted"));
  helePt_all_predicted->Rebin(rebinvalue_pt);
  helePt_all_predicted->SetLineColor(kBlue);
  helePt_all_predicted->SetFillColor(kWhite);
  helePt_all_predicted->SetLineWidth(2.);
  helePt_all_predicted->SetName(predicted);
  helePt_all_predicted->SetMarkerStyle(28);
  helePt_all_predicted->SetMarkerColor(kBlue);
  helePt_all_predicted->SetMarkerSize(1.2);
  if(helePt_all_observed->GetMaximum() >= helePt_all_predicted->GetMaximum() ) {
    helePt_all_observed->Draw();
    helePt_all_predicted->Draw("sames");
  }
  else {
    helePt_all_predicted->Draw();
    helePt_all_observed->Draw("sames");
  }

  leg = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetBorderSize(1);
  //  leg->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leg->AddEntry(helePt_all_observed,       observed,"lpf");
  leg->AddEntry(helePt_all_predicted,      predicted,"lpf");
  leg->Draw();


  elept_all->Update();

  TPaveStats *fake_rate_stats = (TPaveStats*)(helePt_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats != 0 ) {
    fake_rate_stats->SetX1NDC(0.8);
    fake_rate_stats->SetY1NDC(0.75);
    fake_rate_stats->SetX2NDC(0.99);
    fake_rate_stats->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats1 = (TPaveStats*)(helePt_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats1 != 0 ) {
    fake_rate_stats1->SetX1NDC(0.8);
    fake_rate_stats1->SetY1NDC(0.5);
    fake_rate_stats1->SetX2NDC(0.99);
    fake_rate_stats1->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);
  elept_all->Update();


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
  eleeta_all->cd();
  TH1F* heleEta_all_observed = (TH1F*) ((_file0->Get(sample+"elEta_all"))->Clone("heleEta_all_observed"));
  heleEta_all_observed->Rebin(rebinvalue_eta);
  heleEta_all_observed->SetLineColor(kRed);
  heleEta_all_observed->SetFillColor(kWhite);
  heleEta_all_observed->SetLineWidth(2.);
  heleEta_all_observed->SetName(observed);
  heleEta_all_observed->GetYaxis()->SetTitle("Events");
  heleEta_all_observed->GetYaxis()->SetTitleOffset(1.2);
  heleEta_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  heleEta_all_observed->GetXaxis()->SetTitle("#eta^{electron}");
  heleEta_all_observed->GetXaxis()->SetTitleOffset(1.2);
  heleEta_all_observed->GetXaxis()->SetTitleSize(0.04);
  heleEta_all_observed->SetMarkerStyle(20);
  heleEta_all_observed->SetMarkerColor(kRed);
  heleEta_all_observed->SetMarkerSize(1.1);
  //  heleEta_all_observed->Draw();

  TH1F* heleEta_all_predicted = (TH1F*) ((_file2->Get(sample+"elEta_all"))->Clone("heleEta_all_predicted"));
  heleEta_all_predicted->Rebin(rebinvalue_eta);
  heleEta_all_predicted->SetLineColor(kBlue);
  heleEta_all_predicted->SetFillColor(kWhite);
  heleEta_all_predicted->SetLineWidth(2.);
  heleEta_all_predicted->SetName(predicted);
  heleEta_all_predicted->SetMarkerStyle(28);
  heleEta_all_predicted->SetMarkerColor(kBlue);
  heleEta_all_predicted->SetMarkerSize(1.2);
  //  heleEta_all_predicted->Draw("sames");
  if(heleEta_all_observed->GetMaximum() >= heleEta_all_predicted->GetMaximum() ) {
    heleEta_all_observed->Draw();
    heleEta_all_predicted->Draw("sames");
  }
  else {
    heleEta_all_predicted->Draw();
    heleEta_all_observed->Draw("sames");
  }

  eleeta_all->Update();
//    lable_new->Draw();
//    lable_ref->Draw();
  leteleEta = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leteleEta->SetLineColor(1);
  leteleEta->SetLineStyle(1);
  leteleEta->SetLineWidth(1);
  leteleEta->SetFillColor(10);
  leteleEta->SetBorderSize(1);
  //  leteleEta->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leteleEta->AddEntry(heleEta_all_observed,       observed,"lpf");
  leteleEta->AddEntry(heleEta_all_predicted,      predicted,"lpf");
  leteleEta->Draw();


  TPaveStats *fake_rate_stats5 = (TPaveStats*)(heleEta_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats5 != 0 ) {
    fake_rate_stats5->SetX1NDC(0.8);
    fake_rate_stats5->SetY1NDC(0.75);
    fake_rate_stats5->SetX2NDC(0.99);
    fake_rate_stats5->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats6 = (TPaveStats*)(heleEta_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats6 != 0 ) {
    fake_rate_stats6->SetX1NDC(0.8);
    fake_rate_stats6->SetY1NDC(0.5);
    fake_rate_stats6->SetX2NDC(0.99);
    fake_rate_stats6->SetY2NDC(0.74);
  }
  eleeta_all->Update();
  leg->Draw();


  // plot nJet comparison
  nJet_all->cd();
  //  nJet_all->SetLogy();
  TH1F* hnJet_all_observed = (TH1F*) ((_file0->Get(sample+"nJet_all"))->Clone("hnJet_all_observed"));
  //  hnJet_all_observed->Rebin(5);
  hnJet_all_observed->SetLineColor(kRed);
  hnJet_all_observed->SetFillColor(kWhite);
  hnJet_all_observed->SetLineWidth(2.);
  hnJet_all_observed->SetName(observed);
  hnJet_all_observed->GetYaxis()->SetTitle("Events");
  hnJet_all_observed->GetYaxis()->SetTitleOffset(1.2);
  hnJet_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  hnJet_all_observed->GetXaxis()->SetTitle("Number of Jets");
  hnJet_all_observed->GetXaxis()->SetTitleOffset(1.2);
  hnJet_all_observed->GetXaxis()->SetTitleSize(0.04);
  hnJet_all_observed->SetMarkerStyle(20);
  hnJet_all_observed->SetMarkerColor(kRed);
  hnJet_all_observed->SetMarkerSize(1.1);
  //  hnJet_all_observed->Draw();

  TH1F* hnJet_all_predicted = (TH1F*) ((_file2->Get(sample+"nJet_all"))->Clone("hnJet_all_predicted"));
  //  hnJet_all_predicted->Rebin(5);
  hnJet_all_predicted->SetLineColor(kBlue);
  hnJet_all_predicted->SetFillColor(kWhite);
  hnJet_all_predicted->SetLineWidth(2.);
  hnJet_all_predicted->SetName(predicted);
  hnJet_all_predicted->SetMarkerStyle(28);
  hnJet_all_predicted->SetMarkerColor(kBlue);
  hnJet_all_predicted->SetMarkerSize(1.2);
  //  hnJet_all_predicted->Draw("sames");
  if(hnJet_all_observed->GetMaximum() >= hnJet_all_predicted->GetMaximum() ) {
    hnJet_all_observed->Draw();
    hnJet_all_predicted->Draw("sames");
  }
  else {
    hnJet_all_predicted->Draw();
    hnJet_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  nJet_all->Update();
  letnJet = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letnJet->SetLineColor(1);
  letnJet->SetLineStyle(1);
  letnJet->SetLineWidth(1);
  letnJet->SetFillColor(10);
  letnJet->SetBorderSize(1);
  //  letnJet->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letnJet->AddEntry(hnJet_all_observed,       observed,"lpf");
  letnJet->AddEntry(hnJet_all_predicted,      predicted,"lpf");
  letnJet->Draw();


  nJet_all->Update();
  //  nJet_all->SetNDC();
  TPaveStats *fake_rate_stats2 = (TPaveStats*)(hnJet_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2 != 0 ) {
    fake_rate_stats2->SetX1NDC(0.8);
    fake_rate_stats2->SetY1NDC(0.75);
    fake_rate_stats2->SetX2NDC(0.99);
    fake_rate_stats2->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3 = (TPaveStats*)(hnJet_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3 != 0 ) {
    fake_rate_stats3->SetX1NDC(0.8);
    fake_rate_stats3->SetY1NDC(0.5);
    fake_rate_stats3->SetX2NDC(0.99);
    fake_rate_stats3->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  nJet_all->Update();
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
  elPdgId_all->cd();
  //  elPdgId_all->SetLogy();
  TH1F* helPdgId_all_observed = (TH1F*) ((_file0->Get(sample+"elPdgId_all"))->Clone("helPdgId_all_observed"));
  //  helPdgId_all_observed->Rebin(5);
  helPdgId_all_observed->SetLineColor(kRed);
  helPdgId_all_observed->SetFillColor(kWhite);
  helPdgId_all_observed->SetLineWidth(2.);
  helPdgId_all_observed->SetName(observed);
  helPdgId_all_observed->GetYaxis()->SetTitle("Events");
  helPdgId_all_observed->GetYaxis()->SetTitleOffset(1.2);
  helPdgId_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  helPdgId_all_observed->GetXaxis()->SetTitle("electron Pdg ID");
  helPdgId_all_observed->GetXaxis()->SetTitleOffset(1.2);
  helPdgId_all_observed->GetXaxis()->SetTitleSize(0.04);
  helPdgId_all_observed->SetMarkerStyle(20);
  helPdgId_all_observed->SetMarkerColor(kRed);
  helPdgId_all_observed->SetMarkerSize(1.1);
  //  helPdgId_all_observed->Draw();

  TH1F* helPdgId_all_predicted = (TH1F*) ((_file2->Get(sample+"elPdgId_all"))->Clone("helPdgId_all_predicted"));
  //  helPdgId_all_predicted->Rebin(5);
  helPdgId_all_predicted->SetLineColor(kBlue);
  helPdgId_all_predicted->SetFillColor(kWhite);
  helPdgId_all_predicted->SetLineWidth(2.);
  helPdgId_all_predicted->SetName(predicted);
  helPdgId_all_predicted->SetMarkerStyle(28);
  helPdgId_all_predicted->SetMarkerColor(kBlue);
  helPdgId_all_predicted->SetMarkerSize(1.2);
  //  helPdgId_all_predicted->Draw("sames");
  if(helPdgId_all_observed->GetMaximum() >= helPdgId_all_predicted->GetMaximum() ) {
    helPdgId_all_observed->Draw();
    helPdgId_all_predicted->Draw("sames");
  }
  else {
    helPdgId_all_predicted->Draw();
    helPdgId_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elPdgId_all->Update();
  letelPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelPdgId->SetLineColor(1);
  letelPdgId->SetLineStyle(1);
  letelPdgId->SetLineWidth(1);
  letelPdgId->SetFillColor(10);
  letelPdgId->SetBorderSize(1);
  //  letelPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelPdgId->AddEntry(helPdgId_all_observed,       observed,"lpf");
  letelPdgId->AddEntry(helPdgId_all_predicted,      predicted,"lpf");
  letelPdgId->Draw();


  elPdgId_all->Update();
  //  elPdgId_all->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helPdgId_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elPdgIdstats2 = (TPaveStats*)(helPdgId_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elPdgIdstats2 != 0 ) {
    fake_rate_elPdgIdstats2->SetX1NDC(0.8);
    fake_rate_elPdgIdstats2->SetY1NDC(0.5);
    fake_rate_elPdgIdstats2->SetX2NDC(0.99);
    fake_rate_elPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elPdgId_all->Update();
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
  elMoPdgId_all->cd();
  //  elMoPdgId_all->SetLogy();
  TH1F* helMoPdgId_all_observed = (TH1F*) ((_file0->Get(sample+"elMoPdgId_all"))->Clone("helMoPdgId_all_observed"));
  //  helMoPdgId_all_observed->Rebin(5);
  helMoPdgId_all_observed->SetLineColor(kRed);
  helMoPdgId_all_observed->SetFillColor(kWhite);
  helMoPdgId_all_observed->SetLineWidth(2.);
  helMoPdgId_all_observed->SetName(observed);
  helMoPdgId_all_observed->GetYaxis()->SetTitle("Events");
  helMoPdgId_all_observed->GetYaxis()->SetTitleOffset(1.2);
  helMoPdgId_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  helMoPdgId_all_observed->GetXaxis()->SetTitle("electron Mother Pdg ID");
  helMoPdgId_all_observed->GetXaxis()->SetTitleOffset(1.2);
  helMoPdgId_all_observed->GetXaxis()->SetTitleSize(0.04);
  helMoPdgId_all_observed->SetMarkerStyle(20);
  helMoPdgId_all_observed->SetMarkerColor(kRed);
  helMoPdgId_all_observed->SetMarkerSize(1.1);
  //  helMoPdgId_all_observed->Draw();

  TH1F* helMoPdgId_all_predicted = (TH1F*) ((_file2->Get(sample+"elMoPdgId_all"))->Clone("helMoPdgId_all_predicted"));
  //  helMoPdgId_all_predicted->Rebin(5);
  helMoPdgId_all_predicted->SetLineColor(kBlue);
  helMoPdgId_all_predicted->SetFillColor(kWhite);
  helMoPdgId_all_predicted->SetLineWidth(2.);
  helMoPdgId_all_predicted->SetName(predicted);
  helMoPdgId_all_predicted->SetMarkerStyle(28);
  helMoPdgId_all_predicted->SetMarkerColor(kBlue);
  helMoPdgId_all_predicted->SetMarkerSize(1.2);
  //  helMoPdgId_all_predicted->Draw("sames");
  if(helMoPdgId_all_observed->GetMaximum() >= helMoPdgId_all_predicted->GetMaximum() ) {
    helMoPdgId_all_observed->Draw();
    helMoPdgId_all_predicted->Draw("sames");
  }
  else {
    helMoPdgId_all_predicted->Draw();
    helMoPdgId_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elMoPdgId_all->Update();
  letelMoPdgId = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelMoPdgId->SetLineColor(1);
  letelMoPdgId->SetLineStyle(1);
  letelMoPdgId->SetLineWidth(1);
  letelMoPdgId->SetFillColor(10);
  letelMoPdgId->SetBorderSize(1);
  //  letelMoPdgId->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelMoPdgId->AddEntry(helMoPdgId_all_observed,       observed,"lpf");
  letelMoPdgId->AddEntry(helMoPdgId_all_predicted,      predicted,"lpf");
  letelMoPdgId->Draw();


  elMoPdgId_all->Update();
  //  elMoPdgId_all->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helMoPdgId_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elMoPdgIdstats2 = (TPaveStats*)(helMoPdgId_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elMoPdgIdstats2 != 0 ) {
    fake_rate_elMoPdgIdstats2->SetX1NDC(0.8);
    fake_rate_elMoPdgIdstats2->SetY1NDC(0.5);
    fake_rate_elMoPdgIdstats2->SetX2NDC(0.99);
    fake_rate_elMoPdgIdstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elMoPdgId_all->Update();
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
  elPdgIdCat_all->cd();
  //  elPdgIdCat_all->SetLogy();
  TH1F* helPdgIdCat_all_observed = (TH1F*) ((_file0->Get(sample+"elPdgIdCat_all"))->Clone("helPdgIdCat_all_observed"));
  //  helPdgIdCat_all_observed->Rebin(5);
  helPdgIdCat_all_observed->SetLineColor(kRed);
  helPdgIdCat_all_observed->SetFillColor(kWhite);
  helPdgIdCat_all_observed->SetLineWidth(2.);
  helPdgIdCat_all_observed->SetName(observed);
  helPdgIdCat_all_observed->GetYaxis()->SetTitle("Events");
  helPdgIdCat_all_observed->GetYaxis()->SetTitleOffset(1.2);
  helPdgIdCat_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  helPdgIdCat_all_observed->GetXaxis()->SetTitle("electron Pdg ID");
  helPdgIdCat_all_observed->GetXaxis()->SetTitleOffset(1.2);
  helPdgIdCat_all_observed->GetXaxis()->SetTitleSize(0.04);
  helPdgIdCat_all_observed->SetMarkerStyle(20);
  helPdgIdCat_all_observed->SetMarkerColor(kRed);
  helPdgIdCat_all_observed->SetMarkerSize(1.1);
  //  helPdgIdCat_all_observed->Draw();

  TH1F* helPdgIdCat_all_predicted = (TH1F*) ((_file2->Get(sample+"elPdgIdCat_all"))->Clone("helPdgIdCat_all_predicted"));
  //  helPdgIdCat_all_predicted->Rebin(5);
  helPdgIdCat_all_predicted->SetLineColor(kBlue);
  helPdgIdCat_all_predicted->SetFillColor(kWhite);
  helPdgIdCat_all_predicted->SetLineWidth(2.);
  helPdgIdCat_all_predicted->SetName(predicted);
  helPdgIdCat_all_predicted->SetMarkerStyle(28);
  helPdgIdCat_all_predicted->SetMarkerColor(kBlue);
  helPdgIdCat_all_predicted->SetMarkerSize(1.2);
  //  helPdgIdCat_all_predicted->Draw("sames");
  if(helPdgIdCat_all_observed->GetMaximum() >= helPdgIdCat_all_predicted->GetMaximum() ) {
    helPdgIdCat_all_observed->Draw();
    helPdgIdCat_all_predicted->Draw("sames");
  }
  else {
    helPdgIdCat_all_predicted->Draw();
    helPdgIdCat_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  elPdgIdCat_all->Update();
  letelPdgIdCat = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letelPdgIdCat->SetLineColor(1);
  letelPdgIdCat->SetLineStyle(1);
  letelPdgIdCat->SetLineWidth(1);
  letelPdgIdCat->SetFillColor(10);
  letelPdgIdCat->SetBorderSize(1);
  //  letelPdgIdCat->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letelPdgIdCat->AddEntry(helPdgIdCat_all_observed,       observed,"lpf");
  letelPdgIdCat->AddEntry(helPdgIdCat_all_predicted,      predicted,"lpf");
  letelPdgIdCat->Draw();


  elPdgIdCat_all->Update();
  //  elPdgIdCat_all->SetNDC();
  TPaveStats *fake_rate_elIdstats1 = (TPaveStats*)(helPdgIdCat_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elIdstats1 != 0 ) {
    fake_rate_elIdstats1->SetX1NDC(0.8);
    fake_rate_elIdstats1->SetY1NDC(0.75);
    fake_rate_elIdstats1->SetX2NDC(0.99);
    fake_rate_elIdstats1->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_elPdgIdCatstats2 = (TPaveStats*)(helPdgIdCat_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_elPdgIdCatstats2 != 0 ) {
    fake_rate_elPdgIdCatstats2->SetX1NDC(0.8);
    fake_rate_elPdgIdCatstats2->SetY1NDC(0.5);
    fake_rate_elPdgIdCatstats2->SetX2NDC(0.99);
    fake_rate_elPdgIdCatstats2->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  elPdgIdCat_all->Update();
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
  celeRelIso_all->cd();
  //  celeRelIso_all->SetLogy();
  TH1F* heleRelIso_all_observed = (TH1F*) ((_file0->Get(sample+"eleRelIso_all"))->Clone("heleRelIso_all_observed"));
  //  heleRelIso_all_observed->Rebin(5);
  heleRelIso_all_observed->SetLineColor(kRed);
  heleRelIso_all_observed->SetFillColor(kWhite);
  heleRelIso_all_observed->SetLineWidth(2.);
  heleRelIso_all_observed->SetName(observed);
  heleRelIso_all_observed->GetYaxis()->SetTitle("Events");
  heleRelIso_all_observed->GetYaxis()->SetTitleOffset(1.2);
  heleRelIso_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  heleRelIso_all_observed->GetXaxis()->SetTitle("Rel Iso");
  heleRelIso_all_observed->GetXaxis()->SetTitleOffset(1.2);
  heleRelIso_all_observed->GetXaxis()->SetTitleSize(0.04);
  heleRelIso_all_observed->SetMarkerStyle(20);
  heleRelIso_all_observed->SetMarkerColor(kRed);
  heleRelIso_all_observed->SetMarkerSize(1.1);
  //  heleRelIso_all_observed->Draw();

  TH1F* heleRelIso_all_predicted = (TH1F*) ((_file2->Get(sample+"eleRelIso_all"))->Clone("heleRelIso_all_predicted"));
  //  heleRelIso_all_predicted->Rebin(5);
  heleRelIso_all_predicted->SetLineColor(kBlue);
  heleRelIso_all_predicted->SetFillColor(kWhite);
  heleRelIso_all_predicted->SetLineWidth(2.);
  heleRelIso_all_predicted->SetName(predicted);
  heleRelIso_all_predicted->SetMarkerStyle(28);
  heleRelIso_all_predicted->SetMarkerColor(kBlue);
  heleRelIso_all_predicted->SetMarkerSize(1.2);
  //  heleRelIso_all_predicted->Draw("sames");
  if(heleRelIso_all_observed->GetMaximum() >= heleRelIso_all_predicted->GetMaximum() ) {
    heleRelIso_all_observed->Draw();
    heleRelIso_all_predicted->Draw("sames");
  }
  else {
    heleRelIso_all_predicted->Draw();
    heleRelIso_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  celeRelIso_all->Update();
  leteleRelIso = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  leteleRelIso->SetLineColor(1);
  leteleRelIso->SetLineStyle(1);
  leteleRelIso->SetLineWidth(1);
  leteleRelIso->SetFillColor(10);
  leteleRelIso->SetBorderSize(1);
  //  leteleRelIso->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  leteleRelIso->AddEntry(heleRelIso_all_observed,       observed,"lpf");
  leteleRelIso->AddEntry(heleRelIso_all_predicted,      predicted,"lpf");
  leteleRelIso->Draw();


  celeRelIso_all->Update();
  //  eleRelIso_all->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(heleRelIso_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(heleRelIso_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  celeRelIso_all->Update();
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
  // plot dilMass comparison
  cdilMass_all->cd();
  //  cdilMass_all->SetLogy();
  TH1F* hdilMass_all_observed = (TH1F*) ((_file0->Get(sample+"dilMass_all"))->Clone("hdilMass_all_observed"));
  //  hdilMass_all_observed->Rebin(5);
  hdilMass_all_observed->SetLineColor(kRed);
  hdilMass_all_observed->SetFillColor(kWhite);
  hdilMass_all_observed->SetLineWidth(2.);
  hdilMass_all_observed->SetName(observed);
  hdilMass_all_observed->GetYaxis()->SetTitle("Events");
  hdilMass_all_observed->GetYaxis()->SetTitleOffset(1.2);
  hdilMass_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  hdilMass_all_observed->GetXaxis()->SetTitle("Number of Jets");
  hdilMass_all_observed->GetXaxis()->SetTitleOffset(1.2);
  hdilMass_all_observed->GetXaxis()->SetTitleSize(0.04);
  hdilMass_all_observed->SetMarkerStyle(20);
  hdilMass_all_observed->SetMarkerColor(kRed);
  hdilMass_all_observed->SetMarkerSize(1.1);
  //  hdilMass_all_observed->Draw();

  TH1F* hdilMass_all_predicted = (TH1F*) ((_file2->Get(sample+"dilMass_all"))->Clone("hdilMass_all_predicted"));
  //  hdilMass_all_predicted->Rebin(5);
  hdilMass_all_predicted->SetLineColor(kBlue);
  hdilMass_all_predicted->SetFillColor(kWhite);
  hdilMass_all_predicted->SetLineWidth(2.);
  hdilMass_all_predicted->SetName(predicted);
  hdilMass_all_predicted->SetMarkerStyle(28);
  hdilMass_all_predicted->SetMarkerColor(kBlue);
  hdilMass_all_predicted->SetMarkerSize(1.2);
  //  hdilMass_all_predicted->Draw("sames");
  if(hdilMass_all_observed->GetMaximum() >= hdilMass_all_predicted->GetMaximum() ) {
    hdilMass_all_observed->Draw();
    hdilMass_all_predicted->Draw("sames");
  }
  else {
    hdilMass_all_predicted->Draw();
    hdilMass_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  cdilMass_all->Update();
  letdilMass = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letdilMass->SetLineColor(1);
  letdilMass->SetLineStyle(1);
  letdilMass->SetLineWidth(1);
  letdilMass->SetFillColor(10);
  letdilMass->SetBorderSize(1);
  //  letdilMass->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letdilMass->AddEntry(hdilMass_all_observed,       observed,"lpf");
  letdilMass->AddEntry(hdilMass_all_predicted,      predicted,"lpf");
  letdilMass->Draw();


  cdilMass_all->Update();
  //  dilMass_all->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(hdilMass_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(hdilMass_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  cdilMass_all->Update();
   TLatex *   lable_refdilMass = new TLatex(0.55,0.75,"");
   lable_refdilMass->SetNDC();
   lable_refdilMass->SetTextSize(0.04);
   lable_refdilMass->SetTextColor(kRed);
   TLatex *   lable_newdilMass = new TLatex(0.55,0.7,"");
   lable_newdilMass->SetNDC();
   lable_newdilMass->SetTextSize(0.04);
   lable_newdilMass->SetTextColor(kBlue);
   lable_newdilMass->Draw();
   lable_refdilMass->Draw();

   //paste123
  cnHyp_all->cd();
  //  cnHyp_all->SetLogy();
  TH1F* hnHyp_all_observed = (TH1F*) ((_file0->Get(sample+"nHyp_all"))->Clone("hnHyp_all_observed"));
  //  hnHyp_all_observed->Rebin(5);
  hnHyp_all_observed->SetLineColor(kRed);
  hnHyp_all_observed->SetFillColor(kWhite);
  hnHyp_all_observed->SetLineWidth(2.);
  hnHyp_all_observed->SetName(observed);
  hnHyp_all_observed->GetYaxis()->SetTitle("Events");
  hnHyp_all_observed->GetYaxis()->SetTitleOffset(1.2);
  hnHyp_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  hnHyp_all_observed->GetXaxis()->SetTitle("Number of Jets");
  hnHyp_all_observed->GetXaxis()->SetTitleOffset(1.2);
  hnHyp_all_observed->GetXaxis()->SetTitleSize(0.04);
  hnHyp_all_observed->SetMarkerStyle(20);
  hnHyp_all_observed->SetMarkerColor(kRed);
  hnHyp_all_observed->SetMarkerSize(1.1);
  //  hnHyp_all_observed->Draw();

  TH1F* hnHyp_all_predicted = (TH1F*) ((_file2->Get(sample+"nHyp_all"))->Clone("hnHyp_all_predicted"));
  //  hnHyp_all_predicted->Rebin(5);
  hnHyp_all_predicted->SetLineColor(kBlue);
  hnHyp_all_predicted->SetFillColor(kWhite);
  hnHyp_all_predicted->SetLineWidth(2.);
  hnHyp_all_predicted->SetName(predicted);
  hnHyp_all_predicted->SetMarkerStyle(28);
  hnHyp_all_predicted->SetMarkerColor(kBlue);
  hnHyp_all_predicted->SetMarkerSize(1.2);
  //  hnHyp_all_predicted->Draw("sames");
  if(hnHyp_all_observed->GetMaximum() >= hnHyp_all_predicted->GetMaximum() ) {
    hnHyp_all_observed->Draw();
    hnHyp_all_predicted->Draw("sames");
  }
  else {
    hnHyp_all_predicted->Draw();
    hnHyp_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  cnHyp_all->Update();
  letnHyp = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letnHyp->SetLineColor(1);
  letnHyp->SetLineStyle(1);
  letnHyp->SetLineWidth(1);
  letnHyp->SetFillColor(10);
  letnHyp->SetBorderSize(1);
  //  letnHyp->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letnHyp->AddEntry(hnHyp_all_observed,       observed,"lpf");
  letnHyp->AddEntry(hnHyp_all_predicted,      predicted,"lpf");
  letnHyp->Draw();


  cnHyp_all->Update();
  //  nHyp_all->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(hnHyp_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(hnHyp_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  cnHyp_all->Update();
   TLatex *   lable_refnHyp = new TLatex(0.55,0.75,"");
   lable_refnHyp->SetNDC();
   lable_refnHyp->SetTextSize(0.04);
   lable_refnHyp->SetTextColor(kRed);
   TLatex *   lable_newnHyp = new TLatex(0.55,0.7,"");
   lable_newnHyp->SetNDC();
   lable_newnHyp->SetTextSize(0.04);
   lable_newnHyp->SetTextColor(kBlue);
   lable_newnHyp->Draw();
   lable_refnHyp->Draw();

  // plot met comparison
  cmet_all->cd();
  //  cmet_all->SetLogy();
  TH1F* hmet_all_observed = (TH1F*) ((_file0->Get(sample+"met_all"))->Clone("hmet_all_observed"));
  hmet_all_observed->Rebin(4);
  cout<<"Warning - rebinning MET plot"<<endl;
  hmet_all_observed->SetLineColor(kRed);
  hmet_all_observed->SetFillColor(kWhite);
  hmet_all_observed->SetLineWidth(2.);
  hmet_all_observed->SetName(observed);
  hmet_all_observed->GetYaxis()->SetTitle("Events");
  hmet_all_observed->GetYaxis()->SetTitleOffset(1.2);
  hmet_all_observed->GetYaxis()->SetTitleSize(0.04);
  //  hmet_all_observed->GetXaxis()->SetTitle("Number of Jets");
  hmet_all_observed->GetXaxis()->SetTitleOffset(1.2);
  hmet_all_observed->GetXaxis()->SetTitleSize(0.04);
  hmet_all_observed->SetMarkerStyle(20);
  hmet_all_observed->SetMarkerColor(kRed);
  hmet_all_observed->SetMarkerSize(1.1);
  //  hmet_all_observed->Draw();

  TH1F* hmet_all_predicted = (TH1F*) ((_file2->Get(sample+"met_all"))->Clone("hmet_all_predicted"));
  hmet_all_predicted->Rebin(4);
  hmet_all_predicted->SetLineColor(kBlue);
  hmet_all_predicted->SetFillColor(kWhite);
  hmet_all_predicted->SetLineWidth(2.);
  hmet_all_predicted->SetName(predicted);
  hmet_all_predicted->SetMarkerStyle(28);
  hmet_all_predicted->SetMarkerColor(kBlue);
  hmet_all_predicted->SetMarkerSize(1.2);
  //  hmet_all_predicted->Draw("sames");
  if(hmet_all_observed->GetMaximum() >= hmet_all_predicted->GetMaximum() ) {
    hmet_all_observed->Draw();
    hmet_all_predicted->Draw("sames");
  }
  else {
    hmet_all_predicted->Draw();
    hmet_all_observed->Draw("sames");
  }

//    lable_new->Draw();
//    lable_ref->Draw();
  cmet_all->Update();
  letmet = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
  letmet->SetLineColor(1);
  letmet->SetLineStyle(1);
  letmet->SetLineWidth(1);
  letmet->SetFillColor(10);
  letmet->SetBorderSize(1);
  //  letmet->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
  letmet->AddEntry(hmet_all_observed,       observed,"lpf");
  letmet->AddEntry(hmet_all_predicted,      predicted,"lpf");
  letmet->Draw();

  cmet_all->Update();
  //  met_all->SetNDC();
  TPaveStats *fake_rate_stats2i = (TPaveStats*)(hmet_all_observed->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats2i != 0 ) {
    fake_rate_stats2i->SetX1NDC(0.8);
    fake_rate_stats2i->SetY1NDC(0.75);
    fake_rate_stats2i->SetX2NDC(0.99);
    fake_rate_stats2i->SetY2NDC(0.99);
  }
  TPaveStats *fake_rate_stats3i = (TPaveStats*)(hmet_all_predicted->GetListOfFunctions()->FindObject("stats"));
  if ( fake_rate_stats3i != 0 ) {
    fake_rate_stats3i->SetX1NDC(0.8);
    fake_rate_stats3i->SetY1NDC(0.5);
    fake_rate_stats3i->SetX2NDC(0.99);
    fake_rate_stats3i->SetY2NDC(0.74);
  }
  //  gPad->SetLeftMargin(0.17);

  //  leg->Draw();
  cmet_all->Update();
   TLatex *   lable_refmet = new TLatex(0.55,0.75,"");
   lable_refmet->SetNDC();
   lable_refmet->SetTextSize(0.04);
   lable_refmet->SetTextColor(kRed);
   TLatex *   lable_newmet = new TLatex(0.55,0.7,"");
   lable_newmet->SetNDC();
   lable_newmet->SetTextSize(0.04);
   lable_newmet->SetTextColor(kBlue);
   lable_newmet->Draw();
   lable_refmet->Draw();


//    elept_all->Print("elept_all.pdf");
//    eleeta_all->Print("eleeta_all.pdf");
//    nJet_all->Print("nJet_all.pdf");
//    MET_all->Print("MET_all_all.pdf");

   elept_all->Print(sample+"elept_all.png");
   eleeta_all->Print(sample+"eleeta_all.png");
//    eleEtaPtBelow20_all->Print("eleEtaPtBelow20_all.png");
//    eleEtaPtAbove20_all->Print("eleEtaPtAbove20_all.png");
   nJet_all->Print(sample+"nJet_all.png");
   celeRelIso_all->Print(sample+"eleRelIso_all.png");
   elMoPdgId_all->Print(sample+"elMoPdgId_all.png");
   elPdgId_all->Print(sample+"elPdgId_all.png");
   elPdgIdCat_all->Print(sample+"elPdgIdCat_all.png");
   cdilMass_all->Print(sample+"dilMass_all.png");
   cmet_all->Print(sample+"met_all.png");
   cnHyp_all->Print(sample+"nHyp_all.png");
   //   MET_all->Print("MET_all_all.png");
   

}
//   //IBLB
//   // Plot eta comparison
//   eleEtaPtBelow20_all->cd();
//   TH1F* heleEtaPtBelow20_all_observed = (TH1F*) ((_file0->Get(sample+"elEtaPtBelow20_all"))->Clone("heleEtaPtBelow20_all_observed"));
//   heleEtaPtBelow20_all_observed->Rebin(rebinvalue_eta);
//   heleEtaPtBelow20_all_observed->SetLineColor(kRed);
//   heleEtaPtBelow20_all_observed->SetFillColor(kWhite);
//   heleEtaPtBelow20_all_observed->SetLineWidth(2.);
//   heleEtaPtBelow20_all_observed->SetName(observed);
//   heleEtaPtBelow20_all_observed->GetYaxis()->SetTitle("Events");
//   heleEtaPtBelow20_all_observed->GetYaxis()->SetTitleOffset(1.2);
//   heleEtaPtBelow20_all_observed->GetYaxis()->SetTitleSize(0.04);
//   heleEtaPtBelow20_all_observed->GetXaxis()->SetTitle("#eta^{electron}");
//   heleEtaPtBelow20_all_observed->GetXaxis()->SetTitleOffset(1.2);
//   heleEtaPtBelow20_all_observed->GetXaxis()->SetTitleSize(0.04);
//   heleEtaPtBelow20_all_observed->SetMarkerStyle(20);
//   heleEtaPtBelow20_all_observed->SetMarkerColor(kRed);
//   heleEtaPtBelow20_all_observed->SetMarkerSize(1.1);
//   heleEtaPtBelow20_all_observed->Draw();

//   TH1F* heleEtaPtBelow20_all_predicted = (TH1F*) ((_file2->Get(sample+"elEtaPtBelow20_all"))->Clone("heleEtaPtBelow20_all_predicted"));
//   heleEtaPtBelow20_all_predicted->Rebin(rebinvalue_eta);
//   heleEtaPtBelow20_all_predicted->SetLineColor(kBlue);
//   heleEtaPtBelow20_all_predicted->SetFillColor(kWhite);
//   heleEtaPtBelow20_all_predicted->SetLineWidth(2.);
//   heleEtaPtBelow20_all_predicted->SetName(predicted);
//   heleEtaPtBelow20_all_predicted->SetMarkerStyle(28);
//   heleEtaPtBelow20_all_predicted->SetMarkerColor(kBlue);
//   heleEtaPtBelow20_all_predicted->SetMarkerSize(1.2);
//   heleEtaPtBelow20_all_predicted->Draw("sames");

//   eleEtaPtBelow20_all->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   leteleEtaPtBelow20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   leteleEtaPtBelow20->SetLineColor(1);
//   leteleEtaPtBelow20->SetLineStyle(1);
//   leteleEtaPtBelow20->SetLineWidth(1);
//   leteleEtaPtBelow20->SetFillColor(10);
//   leteleEtaPtBelow20->SetBorderSize(1);
//   //  leteleEtaPtBelow20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   leteleEtaPtBelow20->AddEntry(heleEtaPtBelow20_all_observed,       observed,"lpf");
//   leteleEtaPtBelow20->AddEntry(heleEtaPtBelow20_all_predicted,      predicted,"lpf");
//   leteleEtaPtBelow20->Draw();


//   TPaveStats *fake_rate_stats5b = (TPaveStats*)(heleEtaPtBelow20_all_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5b != 0 ) {
//     fake_rate_stats5b->SetX1NDC(0.8);
//     fake_rate_stats5b->SetY1NDC(0.75);
//     fake_rate_stats5b->SetX2NDC(0.99);
//     fake_rate_stats5b->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6b = (TPaveStats*)(heleEtaPtBelow20_all_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6b != 0 ) {
//     fake_rate_stats6b->SetX1NDC(0.8);
//     fake_rate_stats6b->SetY1NDC(0.5);
//     fake_rate_stats6b->SetX2NDC(0.99);
//     fake_rate_stats6b->SetY2NDC(0.74);
//   }
//   eleEtaPtBelow20_all->Update();
//   leg->Draw();
//   //IBLE
//   // Plot eta comparison
//   eleEtaPtAbove20_all->cd();
//   TH1F* heleEtaPtAbove20_all_observed = (TH1F*) ((_file0->Get(sample+"elEtaPtAbove20_all"))->Clone("heleEtaPtAbove20_all_observed"));
//   heleEtaPtAbove20_all_observed->Rebin(rebinvalue_eta);
//   heleEtaPtAbove20_all_observed->SetLineColor(kRed);
//   heleEtaPtAbove20_all_observed->SetFillColor(kWhite);
//   heleEtaPtAbove20_all_observed->SetLineWidth(2.);
//   heleEtaPtAbove20_all_observed->SetName(observed);
//   heleEtaPtAbove20_all_observed->GetYaxis()->SetTitle("Events");
//   heleEtaPtAbove20_all_observed->GetYaxis()->SetTitleOffset(1.2);
//   heleEtaPtAbove20_all_observed->GetYaxis()->SetTitleSize(0.04);
//   heleEtaPtAbove20_all_observed->GetXaxis()->SetTitle("#eta^{electron}");
//   heleEtaPtAbove20_all_observed->GetXaxis()->SetTitleOffset(1.2);
//   heleEtaPtAbove20_all_observed->GetXaxis()->SetTitleSize(0.04);
//   heleEtaPtAbove20_all_observed->SetMarkerStyle(20);
//   heleEtaPtAbove20_all_observed->SetMarkerColor(kRed);
//   heleEtaPtAbove20_all_observed->SetMarkerSize(1.1);
//   heleEtaPtAbove20_all_observed->Draw();

//   TH1F* heleEtaPtAbove20_all_predicted = (TH1F*) ((_file2->Get(sample+"elEtaPtAbove20_all"))->Clone("heleEtaPtAbove20_all_predicted"));
//   heleEtaPtAbove20_all_predicted->Rebin(rebinvalue_eta);
//   heleEtaPtAbove20_all_predicted->SetLineColor(kBlue);
//   heleEtaPtAbove20_all_predicted->SetFillColor(kWhite);
//   heleEtaPtAbove20_all_predicted->SetLineWidth(2.);
//   heleEtaPtAbove20_all_predicted->SetName(predicted);
//   heleEtaPtAbove20_all_predicted->SetMarkerStyle(28);
//   heleEtaPtAbove20_all_predicted->SetMarkerColor(kBlue);
//   heleEtaPtAbove20_all_predicted->SetMarkerSize(1.2);
//   heleEtaPtAbove20_all_predicted->Draw("sames");

//   eleEtaPtAbove20_all->Update();
// //    lable_new->Draw();
// //    lable_ref->Draw();
//   leteleEtaPtAbove20 = new TLegend(0.5,0.75,0.79,0.99,NULL,"brNDC");
//   leteleEtaPtAbove20->SetLineColor(1);
//   leteleEtaPtAbove20->SetLineStyle(1);
//   leteleEtaPtAbove20->SetLineWidth(1);
//   leteleEtaPtAbove20->SetFillColor(10);
//   leteleEtaPtAbove20->SetBorderSize(1);
//   //  leteleEtaPtAbove20->SetHeader(                               "L1 CSC trigger efficiency - With ME42");
//   leteleEtaPtAbove20->AddEntry(heleEtaPtAbove20_all_observed,       observed,"lpf");
//   leteleEtaPtAbove20->AddEntry(heleEtaPtAbove20_all_predicted,      predicted,"lpf");
//   leteleEtaPtAbove20->Draw();


//   TPaveStats *fake_rate_stats5a = (TPaveStats*)(heleEtaPtAbove20_all_observed->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats5a != 0 ) {
//     fake_rate_stats5a->SetX1NDC(0.8);
//     fake_rate_stats5a->SetY1NDC(0.75);
//     fake_rate_stats5a->SetX2NDC(0.99);
//     fake_rate_stats5a->SetY2NDC(0.99);
//   }
//   TPaveStats *fake_rate_stats6a = (TPaveStats*)(heleEtaPtAbove20_all_predicted->GetListOfFunctions()->FindObject("stats"));
//   if ( fake_rate_stats6a != 0 ) {
//     fake_rate_stats6a->SetX1NDC(0.8);
//     fake_rate_stats6a->SetY1NDC(0.5);
//     fake_rate_stats6a->SetX2NDC(0.99);
//     fake_rate_stats6a->SetY2NDC(0.74);
//   }
//   eleEtaPtAbove20_all->Update();
//   leg->Draw();

  //IBLE2
