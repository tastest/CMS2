{
  f = TFile::Open("processed_data_tag.root");
  gROOT->SetStyle("Plain");
  c1 = new TCanvas("c1","c1",500,500);
  s = new THStack("s","");
  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  
  ttbar_hmet_all->Rebin(10);
  tw_hmet_all = (TH1F*)f->Get("tw_hmet_all");
  tw_hmet_all->Rebin(10);
  ttbar_hmet_all->Add(tw_hmet_all);
  ttbar_hmet_all->SetFillColor(kYellow);
  s->Add(ttbar_hmet_all);
  leg->AddEntry(ttbar_hmet_all, "Top", "lpf");	
  
  dytt_hmet_all->Rebin(10);
  dyee_hmet_all = (TH1F*)f->Get("dyee_hmet_all");
  dyee_hmet_all->Rebin(10);
  dytt_hmet_all->Add(dyee_hmet_all);
  dymm_hmet_all = (TH1F*)f->Get("dymm_hmet_all");
  dymm_hmet_all->Rebin(10);
  dytt_hmet_all->Add(dymm_hmet_all);
  dytt_hmet_all->SetFillColor(kBlack);
  s->Add(dytt_hmet_all);
  leg->AddEntry(dytt_hmet_all, "DY", "lpf");	
  
  wz_hmet_all->Rebin(10);
  zz_hmet_all = (TH1F*)f->Get("zz_hmet_all");
  zz_hmet_all->Rebin(10);
  wz_hmet_all->Add(zz_hmet_all);
  wz_hmet_all->SetFillColor(kBlue);
  s->Add(wz_hmet_all);
  leg->AddEntry(wz_hmet_all, "ZZ/WZ", "lpf");	

  wjets_hmet_all->Rebin(10);
  wjets_hmet_all->SetFillColor(kGray);
  s->Add(wjets_hmet_all);
  leg->AddEntry(wjets_hmet_all, "Wjets", "lpf");	
    
  ww_hmet_all->Rebin(10);
  ww_hmet_all->SetFillColor(kRed);
  s->Add(ww_hmet_all);
  leg->AddEntry(ww_hmet_all, "WW", "lpf");	
	
  s->Draw("hist");
  s->GetXaxis()->SetTitle("MET, [GeV]");
  s->Draw("hist");
  leg->Draw();
}
