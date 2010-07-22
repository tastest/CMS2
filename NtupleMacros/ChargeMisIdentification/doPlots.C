{

  TCanvas *canvas = new TCanvas;
  TPaveText *box = 0;
  double eff = 0.;
  double barrel_num = 0.;
  double barrel_den = 0.;
  double forward_num = 0.;
  double forward_den = 0.;
  double barrel_eff = 0.;
  double forward_eff = 0.;
  double ineff = 0.;
  double barrel_ineff = 0.;
  double forward_ineff = 0.;
  double trk_eff = 0.;

  double f_b_eta_divider = 1.479;

  TH1F *els_pt_sim = (TH1F*)gFile->Get("els_pt_sim");
  TH1F *els_eta_sim = (TH1F*)gFile->Get("els_eta_sim");

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
  box->InsertText(Form("reco eff: %2.3f %s",eff,"%"));
  box->Draw();
  canvas->Print("RecoEfficiencyPt.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  barrel_num = els_eta_reco->Integral(els_eta_reco->FindBin(-1*f_b_eta_divider),els_eta_reco->FindBin(f_b_eta_divider));
  barrel_den = els_eta_sim->Integral(els_eta_sim->FindBin(-1*f_b_eta_divider),els_eta_sim->FindBin(f_b_eta_divider));
  barrel_eff = barrel_num / barrel_den * 100;
  forward_num = els_eta_reco->Integral(0,els_eta_reco->FindBin(-1*f_b_eta_divider)) + els_eta_reco->Integral(els_eta_reco->FindBin(f_b_eta_divider),els_eta_reco->GetNbinsX()+1);
  forward_den = els_eta_sim->Integral(0,els_eta_sim->FindBin(-1*f_b_eta_divider)) + els_eta_sim->Integral(els_eta_sim->FindBin(f_b_eta_divider),els_eta_sim->GetNbinsX()+1);
  forward_eff = forward_num / forward_den * 100;
  els_eta_reco->Divide(els_eta_sim);
  els_eta_reco->SetName("RecoEff");
  els_eta_reco->SetTitle("Reconstruction efficiency (reco/sim) in #eta");
  els_eta_reco->SetYTitle("Efficiency");
  els_eta_reco->SetMaximum(1.01);
  els_eta_reco->SetMinimum(0.);
  els_eta_reco->Draw();
  box = new TPaveText(-1.5,0.3,1.5,0.5);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("reco eff barrel (|#eta| #leq 1.479): %2.3f %s",barrel_eff,"%"));
  box->InsertText(Form("reco eff forward (|#eta| > 1.479) : %2.3f %s",forward_eff,"%"));
  box->Draw();
  canvas->Print("RecoEfficiencyEta.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  ineff = els_pt_recosim_incorCharge->Integral() / els_pt_recosim->Integral() * 100;
  els_pt_recosim_incorCharge->Divide(els_pt_recosim);
  els_pt_recosim_incorCharge->SetName("ChargeMisIDRate");
  els_pt_recosim_incorCharge->SetTitle("Charge MisId rate (recosim_incorCharge/recosim) in p_{T}");
  els_pt_recosim_incorCharge->SetYTitle("Charge MisId rate");
  els_pt_recosim_incorCharge->SetMaximum(0.1);
  els_pt_recosim_incorCharge->SetMinimum(0.);
  els_pt_recosim_incorCharge->Draw();
  box = new TPaveText(40,0.085,80,0.095);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("charge MisId rate: %2.3f %s",ineff,"%"));
  box->Draw();
  canvas->Print("ChargeMisIdRatePt.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  barrel_num = els_eta_recosim_incorCharge->Integral(els_eta_recosim_incorCharge->FindBin(-1*f_b_eta_divider),els_eta_recosim_incorCharge->FindBin(f_b_eta_divider));
  barrel_den = els_eta_recosim->Integral(els_eta_recosim->FindBin(-1*f_b_eta_divider),els_eta_recosim->FindBin(f_b_eta_divider));
  barrel_ineff = barrel_num / barrel_den * 100;
  forward_num = els_eta_recosim_incorCharge->Integral(0,els_eta_recosim_incorCharge->FindBin(-1*f_b_eta_divider)) + els_eta_recosim_incorCharge->Integral(els_eta_recosim_incorCharge->FindBin(f_b_eta_divider),els_eta_recosim_incorCharge->GetNbinsX()+1);
  forward_den = els_eta_recosim->Integral(0,els_eta_recosim->FindBin(-1*f_b_eta_divider)) + els_eta_recosim->Integral(els_eta_recosim->FindBin(f_b_eta_divider),els_eta_recosim->GetNbinsX()+1);
  forward_ineff = forward_num / forward_den * 100;
  els_eta_recosim_incorCharge->Divide(els_eta_recosim);
  els_eta_recosim_incorCharge->SetName("ChargeMisIdRate");
  els_eta_recosim_incorCharge->SetTitle("Charge MisId rate (recosim_incorCharge/recosim) in #eta");
  els_eta_recosim_incorCharge->SetYTitle("Charge MisId rate");
  els_eta_recosim_incorCharge->SetMaximum(0.1);
  els_eta_recosim_incorCharge->SetMinimum(0.);
  els_eta_recosim_incorCharge->Draw();
  box = new TPaveText(-1.5,0.075,1.5,0.095);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("charge MisId rate barrel (|#eta| #leq 1.479): %2.3f %s",barrel_ineff,"%"));
  box->InsertText(Form("charge MisId rate forward (|#eta| > 1.479) : %2.3f %s",forward_ineff,"%"));
  box->Draw();
  canvas->Print("ChargeMisIdRateEta.png");

  canvas->Clear();
  canvas->SetLeftMargin(0.155);
  canvas->SetBottomMargin(0.155);
  canvas->SetTopMargin(0.14);
  trk_eff = els_trkId->Integral(0,10) / els_trkId->Integral(11,100) * 100;
  els_trkId->SetName("EleWOAssocTrk");
  els_trkId->SetTitle("Reconstructed electrons without associated track");
  els_trkId->SetYTitle("Electrons");
  els_trkId->Draw("HIST");
  box = new TPaveText(-800,600000,-200,800000);
  box->SetBorderSize(1);
  box->SetFillColor(0);
  box->InsertText(Form("els wo assoc. trk: %2.3f %s",trk_eff,"%"));
  box->Draw();
  canvas->Print("TrkWithoutAssocTrk.png");

  std::cout << "| *comment* | *reco efficiency [%]* | *barrel reco efficiency [%]* | *forward reco efficiency [%]* | *charge misidentification rate [%]* | *barrel charge misidentification rate [%]* | *forward charge misidentification rate [%]* |" << std::endl;
  std::cout << "| | " << Form("%3.3f",eff) << " | " << Form("%3.3f",barrel_eff) << " | " << Form("%3.3f",forward_eff) << " | " << Form("%3.3f",ineff) << " | " << Form("%3.3f",barrel_ineff) << " | " << Form("%3.3f",forward_ineff) << " |" << std::endl;


}
