{

  TCanvas *canvas = new TCanvas;
  TPaveText *box = 0;
  double eff = 0.;

  TH1F *els_pt_recosim = (TH1F*)gFile->Get("els_pt_recosim");
  TH1F *els_eta_recosim = (TH1F*)gFile->Get("els_eta_recosim");

  TH1F *els_pt_reco = (TH1F*)gFile->Get("els_pt_reco");
  TH1F *els_eta_reco = (TH1F*)gFile->Get("els_eta_reco");

  TH1F *els_pt_recosim_incorCharge = (TH1F*)gFile->Get("els_pt_recosim_incorCharge");
  TH1F *els_eta_recosim_incorCharge = (TH1F*)gFile->Get("els_eta_recosim_incorCharge");

  TH1F *els_pt_recosim_corCharge = (TH1F*)gFile->Get("els_pt_recosim_corCharge");
  TH1F *els_eta_recosim_corCharge = (TH1F*)gFile->Get("els_eta_recosim_corCharge");

  TH1F *els_pt_reco_corCharge = (TH1F*)gFile->Get("els_pt_reco_corCharge");
  TH1F *els_eta_reco_corCharge = (TH1F*)gFile->Get("els_eta_reco_corCharge");

  TH1F *els_trkId = (TH1F*)gFile->Get("els_trkId");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  eff = els_pt_reco->Integral() / els_pt_sim->Integral() * 100;
  els_pt_reco->Divide(els_pt_sim);
  els_pt_reco->SetName("RecoEff");
  els_pt_reco->SetTitle("Reconstruction efficiency (reco/sim) in p_{T}");
  els_pt_reco->SetYTitle("Efficiency");
  els_pt_reco->SetMaximum(1.01);
  els_pt_reco->SetMinimum(0.);
  els_pt_reco->Draw();
  box = new TPaveText(40.,0.4,80.,0.5);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("reco eff: %2.1f %s",eff,"%"));
  box->Draw();
  canvas->Print("RecoEfficiencyPt.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  eff = els_eta_reco->Integral() / els_eta_sim->Integral() * 100;
  els_eta_reco->Divide(els_eta_sim);
  els_eta_reco->SetName("RecoEff");
  els_eta_reco->SetTitle("Reconstruction efficiency (reco/sim) in #eta");
  els_eta_reco->SetYTitle("Efficiency");
  els_eta_reco->SetMaximum(1.01);
  els_eta_reco->SetMinimum(0.);
  els_eta_reco->Draw();
  box = new TPaveText(-1,0.4,1,0.5);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("reco eff: %2.1f %s",eff,"%"));
  box->Draw();
  canvas->Print("RecoEfficiencyEta.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  eff = els_pt_recosim_incorCharge->Integral() / els_pt_recosim->Integral() * 100;
  els_pt_recosim_incorCharge->Divide(els_pt_recosim);
  els_pt_recosim_incorCharge->SetName("ChargeIDIneff");
  els_pt_recosim_incorCharge->SetTitle("Charge id inefficiency (recosim_incorCharge/recosim) in p_{T}");
  els_pt_recosim_incorCharge->SetYTitle("Inefficiency");
  els_pt_recosim_incorCharge->SetMaximum(0.1);
  els_pt_recosim_incorCharge->SetMinimum(0.);
  els_pt_recosim_incorCharge->Draw();
  box = new TPaveText(40,0.085,80,0.095);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("charge id ineff: %2.1f %s",eff,"%"));
  box->Draw();
  canvas->Print("ChargeIDInefficiencyPt.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  eff = els_eta_recosim_incorCharge->Integral() / els_eta_recosim->Integral() * 100;
  els_eta_recosim_incorCharge->Divide(els_eta_recosim);
  els_eta_recosim_incorCharge->SetName("ChargeIDIneff");
  els_eta_recosim_incorCharge->SetTitle("Charge id inefficiency (recosim_incorCharge/recosim) in #eta");
  els_eta_recosim_incorCharge->SetYTitle("Inefficiency");
  els_eta_recosim_incorCharge->SetMaximum(0.1);
  els_eta_recosim_incorCharge->SetMinimum(0.);
  els_eta_recosim_incorCharge->Draw();
  box = new TPaveText(-1,0.085,1,0.095);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("charge id ineff: %2.1f %s",eff,"%"));
  box->Draw();
  canvas->Print("ChargeIDInefficiencyEta.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  canvas->SetTopMargin(0.14);
  eff = els_trkId->Integral(0,10) / els_trkId->Integral(11,100) * 100;
  els_trkId->SetName("EleWOAssocTrk");
  els_trkId->SetTitle("Reconstructed electrons without associated track");
  els_trkId->SetYTitle("Electrons");
  els_trkId->Draw("HIST");
  box = new TPaveText(-800,600000,-200,800000);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("els wo assoc. trk: %2.1f %s",eff,"%"));
  box->Draw();
  canvas->Print("TrkWithoutAssocTrk.png");


//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_pt_recosim->Draw();
//   canvas->Print("els_pt_recosim.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_eta_recosim->Draw();
//   canvas->Print("els_eta_recosim.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_pt_recosim_incorCharge->Divide(els_pt_recosim);
//   els_pt_recosim_incorCharge->Draw();
//   els_pt_recosim_incorCharge->SetName("inCorChargeFractionPt");
//   els_pt_recosim_incorCharge->SetTitle("Incorrect Charge Fraction in p_{T}");
//   els_pt_recosim_incorCharge->SetYTitle("Fraction");
//   els_pt_recosim_incorCharge->SetMaximum(0.1);
//   canvas->Print("InCorChargeFractionPt.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_eta_recosim_incorCharge->Divide(els_eta_recosim);
//   els_eta_recosim_incorCharge->Draw();
//   els_eta_recosim_incorCharge->SetName("inCorChargeFraction");
//   els_eta_recosim_incorCharge->SetTitle("Incorrect Charge Fraction in #eta");
//   els_eta_recosim_incorCharge->SetYTitle("Fraction");
//   els_eta_recosim_incorCharge->SetMaximum(0.1);
//   canvas->Print("InCorChargeFractionEta.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_pt_recosim_corCharge->Divide(els_pt_recosim);
//   els_pt_recosim_corCharge->Draw();
//   els_pt_recosim_corCharge->SetName("ChargeIdEfficiency");
//   els_pt_recosim_corCharge->SetTitle("Charge Id Efficiency in p_{T}");
//   els_pt_recosim_corCharge->SetYTitle("Efficiency");
//   els_pt_recosim_corCharge->SetMaximum(1.02);
//   els_pt_recosim_corCharge->SetMinimum(0.9);
//   canvas->Print("ChargeIdEfficiencyPt.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_eta_recosim_corCharge->Divide(els_eta_recosim);
//   els_eta_recosim_corCharge->Draw();
//   els_eta_recosim_corCharge->SetName("ChargeIdEfficiency");
//   els_eta_recosim_corCharge->SetTitle("Charge Id Efficiency in #eta");
//   els_eta_recosim_corCharge->SetYTitle("Efficiency");
//   els_eta_recosim_corCharge->SetMaximum(1.02);
//   els_eta_recosim_corCharge->SetMinimum(0.9);
//   canvas->Print("ChargeIdEfficiencyEta.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_pt_reco_corCharge->Divide(els_pt_reco);
//   els_pt_reco_corCharge->Draw();
//   els_pt_reco_corCharge->SetName("ChargeIdPurity");
//   els_pt_reco_corCharge->SetTitle("Charge Id Purity in p_{T}");
//   els_pt_reco_corCharge->SetYTitle("Purity");
//   els_pt_reco_corCharge->SetMaximum(1.02);
//   els_pt_reco_corCharge->SetMinimum(0.9);
//   canvas->Print("ChargeIdPurityPt.png");

//   canvas->Clear();
//   canvas->SetLeftMargin(0.155);
//   canvas->SetBottomMargin(0.155);
//   els_eta_reco_corCharge->Divide(els_eta_reco);
//   els_eta_reco_corCharge->Draw();
//   els_eta_reco_corCharge->SetName("ChargeIdPurity");
//   els_eta_reco_corCharge->SetTitle("Charge Id Purity in p_{T}");
//   els_eta_reco_corCharge->SetYTitle("Purity");
//   els_eta_reco_corCharge->SetMaximum(1.02);
//   els_eta_reco_corCharge->SetMinimum(0.9);
//   canvas->Print("ChargeIdPurityEta.png");

}
