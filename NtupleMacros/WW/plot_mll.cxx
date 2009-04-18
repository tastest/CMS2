{
  f = TFile::Open("processed_data_tag.root");
  gROOT->SetStyle("Plain");
  c1 = new TCanvas("c1","c1",500,500);
  s = new THStack("s","");
  TLegend* leg = new TLegend(0.75, 0.75, 0.99, 0.99);
  
  ttbar_hdilMass_all->Rebin(10);
  tw_hdilMass_all = (TH1F*)f->Get("tw_hdilMass_all");
  tw_hdilMass_all->Rebin(10);
  ttbar_hdilMass_all->Add(tw_hdilMass_all);
  ttbar_hdilMass_all->SetFillColor(kYellow);
  s->Add(ttbar_hdilMass_all);
  leg->AddEntry(ttbar_hdilMass_all, "Top", "lpf");	
  
  dytt_hdilMass_all->Rebin(10);
  dyee_hdilMass_all = (TH1F*)f->Get("dyee_hdilMass_all");
  dyee_hdilMass_all->Rebin(10);
  dytt_hdilMass_all->Add(dyee_hdilMass_all);
  dymm_hdilMass_all = (TH1F*)f->Get("dymm_hdilMass_all");
  dymm_hdilMass_all->Rebin(10);
  dytt_hdilMass_all->Add(dymm_hdilMass_all);
  dytt_hdilMass_all->SetFillColor(kBlack);
  s->Add(dytt_hdilMass_all);
  leg->AddEntry(dytt_hdilMass_all, "DY", "lpf");	
  
  wz_hdilMass_all->Rebin(10);
  zz_hdilMass_all = (TH1F*)f->Get("zz_hdilMass_all");
  zz_hdilMass_all->Rebin(10);
  wz_hdilMass_all->Add(zz_hdilMass_all);
  wz_hdilMass_all->SetFillColor(kBlue);
  s->Add(wz_hdilMass_all);
  leg->AddEntry(wz_hdilMass_all, "ZZ/WZ", "lpf");	

  wjets_hdilMass_all->Rebin(10);
  wjets_hdilMass_all->SetFillColor(kGray);
  s->Add(wjets_hdilMass_all);
  leg->AddEntry(wjets_hdilMass_all, "Wjets", "lpf");	
    
  ww_hdilMass_all->Rebin(10);
  ww_hdilMass_all->SetFillColor(kRed);
  s->Add(ww_hdilMass_all);
  leg->AddEntry(ww_hdilMass_all, "WW", "lpf");	
	
  s->Draw("hist");
  s->GetXaxis()->SetTitle("Mass, [GeV]");
  s->Draw("hist");
  leg->Draw();
}
