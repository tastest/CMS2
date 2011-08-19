{
  f = TFile::Open("processed_data_tag.root");
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  c1 = new TCanvas("c1","c1",500,500);
  s = new THStack("s","");
  
  ttbar_hdilMass_all->Rebin(10);
  tw_hdilMass_all = (TH1F*)f->Get("tw_hdilMass_all");
  tw_hdilMass_all->Rebin(10);
  ttbar_hdilMass_all->Add(tw_hdilMass_all);
  ttbar_hdilMass_all->SetFillColor(kMagenta);
  ttbar_hdilMass_all->SetLineColor(kBlack);
  ttbar_hdilMass_all->SetMarkerColor(kBlack);
  s->Add(ttbar_hdilMass_all);
  
  dytt_hdilMass_all->Rebin(10);
  dyee_hdilMass_all = (TH1F*)f->Get("dyee_hdilMass_all");
  dyee_hdilMass_all->Rebin(10);
  dytt_hdilMass_all->Add(dyee_hdilMass_all);
  dymm_hdilMass_all = (TH1F*)f->Get("dymm_hdilMass_all");
  dymm_hdilMass_all->Rebin(10);
  dytt_hdilMass_all->Add(dymm_hdilMass_all);
  dytt_hdilMass_all->SetFillColor(kBlack);
  dytt_hdilMass_all->SetLineColor(kBlack);
  dytt_hdilMass_all->SetMarkerColor(kBlack);
  s->Add(dytt_hdilMass_all);
  
  wz_hdilMass_all->Rebin(10);
  zz_hdilMass_all = (TH1F*)f->Get("zz_hdilMass_all");
  zz_hdilMass_all->Rebin(10);
  wz_hdilMass_all->Add(zz_hdilMass_all);
  wz_hdilMass_all->SetFillColor(kBlue);
  wz_hdilMass_all->SetLineColor(kBlack);
  wz_hdilMass_all->SetMarkerColor(kBlack);
  s->Add(wz_hdilMass_all);

  wjets_hdilMass_all->Rebin(10);
  wjets_hdilMass_all->SetFillColor(kGray);
  wjets_hdilMass_all->SetLineColor(kBlack);
  wjets_hdilMass_all->SetMarkerColor(kBlack);
  s->Add(wjets_hdilMass_all);
    
  ww_hdilMass_all->Rebin(10);
  ww_hdilMass_all->SetFillColor(kRed);
  ww_hdilMass_all->SetLineColor(kBlack);
  ww_hdilMass_all->SetMarkerColor(kBlack);
  s->Add(ww_hdilMass_all);
	
  s->Draw("hist");
  s->GetXaxis()->SetTitle("Mass [GeV/c^{2}]");
  s->GetYaxis()->SetTitle("Events/(30 GeV/c^{2})");
  /*
  // make it bigger and uglier for ARC
  s->GetXaxis()->SetLabelSize(0.06);
  s->GetXaxis()->SetTitleSize(0.06);
  s->GetYaxis()->SetLabelSize(0.06);
  s->GetYaxis()->SetLabelOffset(0.02);
  */
  s->SetTitle("CMS Preliminary");
  s->Draw("hist");
  
  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);
  leg->AddEntry(ww_hdilMass_all, "WW", "lpf");	
  leg->AddEntry(wjets_hdilMass_all, "Wjets", "lpf");	
  leg->AddEntry(wz_hdilMass_all, "ZZ/WZ", "lpf");	
  leg->AddEntry(dytt_hdilMass_all, "DY", "lpf");	
  leg->AddEntry(ttbar_hdilMass_all, "Top", "lpf");	
  leg->Draw();
  c1->Print("mass.eps");

  f = TFile::Open("mass.root","RECREATE");
  c1->Write();
  f->Close();
}
