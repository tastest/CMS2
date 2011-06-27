{
  f = TFile::Open("processed_data.root");
  gROOT->ProcessLine(".L tdrstyle.C");
  setTDRStyle();
  //  gROOT->SetStyle("Plain");
  c1 = new TCanvas("c1","c1",500,500);
  s = new THStack("s","");
  unsigned int nRebin = 1;
  
  ttbar_hmet_all->Rebin(nRebin);
  ttbar_hmet_all->SetFillColor(kMagenta);
  ttbar_hmet_all->SetLineColor(kBlack);
  ttbar_hmet_all->SetMarkerColor(kBlack);
  tw_hmet_all = (TH1F*)f->Get("tw_hmet_all");
  tw_hmet_all->Rebin(nRebin);
  ttbar_hmet_all->Add(tw_hmet_all);
  s->Add(ttbar_hmet_all);
  
  dytt_hmet_all->Rebin(nRebin);
  dyee_hmet_all = (TH1F*)f->Get("dyee_hmet_all");
  dyee_hmet_all->Rebin(nRebin);
  dytt_hmet_all->Add(dyee_hmet_all);
  dymm_hmet_all = (TH1F*)f->Get("dymm_hmet_all");
  dymm_hmet_all->Rebin(nRebin);
  dytt_hmet_all->Add(dymm_hmet_all);
  dytt_hmet_all->SetFillColor(kBlack);
  s->Add(dytt_hmet_all);
  
  wz_hmet_all->Rebin(nRebin);
  zz_hmet_all = (TH1F*)f->Get("zz_hmet_all");
  zz_hmet_all->Rebin(nRebin);
  wz_hmet_all->Add(zz_hmet_all);
  wz_hmet_all->SetFillColor(kBlue);
  wz_hmet_all->SetMarkerColor(kBlack);
  wz_hmet_all->SetLineColor(kBlack);
  s->Add(wz_hmet_all);

  wjets_hmet_all->Rebin(nRebin);
  wjets_hmet_all->SetFillColor(kGray);
  wjets_hmet_all->SetMarkerColor(kBlack);
  wjets_hmet_all->SetLineColor(kBlack);
  s->Add(wjets_hmet_all);
    
  ww_hmet_all->Rebin(nRebin);
  ww_hmet_all->SetFillColor(kRed);
  ww_hmet_all->SetLineColor(kBlack);
  s->Add(ww_hmet_all);
	
  s->Draw("hist");
  s->GetXaxis()->SetTitle("MET [GeV]");
  s->GetYaxis()->SetTitle("Events/(20 GeV)");
  s->SetTitle("CMS Preliminary");
  s->Draw("hist");

  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  leg->SetShadowColor(0);
  leg->SetFillColor(0);
  wwe = leg->AddEntry(ww_hmet_all, "WW", "lpf");	
  wje = leg->AddEntry(wjets_hmet_all, "Wjets", "lpf");	
  dbe = leg->AddEntry(wz_hmet_all, "ZZ/WZ", "lpf");	
  dye = leg->AddEntry(dytt_hmet_all, "DY", "lpf");	
  tope = leg->AddEntry(ttbar_hmet_all, "Top", "lpf");	
  tope->SetLineColor(616);
  leg->Draw();
  c1->Print("met.eps");

  f = TFile::Open("met.root","RECREATE");
  c1->Write();
  f->Close();
}
