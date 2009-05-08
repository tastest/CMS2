/*
.x setup.C
hist::loadHist("myHist_DYandTT_MG_1433600__OS_noDupWt_isoDil08_preDil08noIso_preMet08_hltMu9E15.root", "mgnozVeto", "*_hnJet_*");
hist::loadHist("myHist_DYandTT_PY_1433600__OS_noDupWt_isoDil08_preDil08noIso_preMet08_hltMu9E15.root", "pynozVeto", "*_hnJet_*");
hist::loadHist("myHist_DYandTT_MG_1957888__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15.root", "mgbase", "*_hnJet_*");
hist::loadHist("myHist_DYandTT_PY_1957888__OS_noDupWt_isoDil08_preDil08noIso_preMet08_outZ08_hltMu9E15.root", "pybase", "*_hnJet_*");
.L plotPYvsMG.C
plotPYvsMG("pybase", "mgbase", "sk_tdilgp032309_1957888");
plotPYvsMG("pynozVeto", "mgnozVeto", "sk_tdilgp032309_1433600");
//beware the long file names

*/

void plotPYvsMG(const char* pypfx, const char* mgpfx, const char* filePfx){
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  TCanvas* cDYee = new TCanvas("cDYee");
  TH1F* pydyeenjet = (TH1F*)gDirectory->Get(Form("%s_DYee_hnJet_ee", pypfx));
  TH1F* mgdyeenjet = (TH1F*)gDirectory->Get(Form("%s_DYee_hnJet_ee", mgpfx));
  pydyeenjet->SetStats(0);
  pydyeenjet->SetTitle("Z/#gamma^{*}#rightarrow e^{+}e^{-}");
  pydyeenjet->SetLineColor(2);
  pydyeenjet->SetFillColor(0);
  pydyeenjet->SetMarkerStyle(22);
  pydyeenjet->SetMarkerColor(2);
  mgdyeenjet->SetLineColor(1);
  mgdyeenjet->SetFillColor(0);
  mgdyeenjet->SetMarkerColor(1);
  cDYee->cd();
  pydyeenjet->Draw();
  mgdyeenjet->Draw("same");
  TLegend* leg = new TLegend(0.75,0.75, 0.95, 0.95);
  leg->AddEntry(pydyeenjet, "PYTHIA");
  leg->AddEntry(mgdyeenjet, "MadGraph");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->SaveAs(Form("%s_DYee_hnJet_ee_%s_v_%s.eps", filePfx, pypfx, mgpfx));

  TCanvas* cDYmm = new TCanvas("cDYmm");
  TH1F* pydymmnjet = (TH1F*)gDirectory->Get(Form("%s_DYmm_hnJet_mm", pypfx));
  TH1F* mgdymmnjet = (TH1F*)gDirectory->Get(Form("%s_DYmm_hnJet_mm", mgpfx));
  pydymmnjet->SetStats(0);
  pydymmnjet->SetTitle("Z/#gamma^{*}#rightarrow #mu^{+}#mu^{-}");
  pydymmnjet->SetLineColor(2);
  pydymmnjet->SetFillColor(0);
  pydymmnjet->SetMarkerStyle(22);
  pydymmnjet->SetMarkerColor(2);
  mgdymmnjet->SetLineColor(1);
  mgdymmnjet->SetFillColor(0);
  mgdymmnjet->SetMarkerColor(1);
  cDYmm->cd();
  pydymmnjet->Draw();
  mgdymmnjet->Draw("same");
  TLegend* leg = new TLegend(0.75,0.75, 0.95, 0.95);
  leg->AddEntry(pydymmnjet, "PYTHIA");
  leg->AddEntry(mgdymmnjet, "MadGraph");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->SaveAs(Form("%s_DYmm_hnJet_mm_%s_v_%s.eps", filePfx, pypfx, mgpfx));

  TCanvas* cTTee = new TCanvas("cTTee");
  TH1F* pytteenjet = (TH1F*)gDirectory->Get(Form("%s_ttdil_hnJet_ee", pypfx));
  TH1F* mgtteenjet = (TH1F*)gDirectory->Get(Form("%s_ttdil_hnJet_ee", mgpfx));
  pytteenjet->SetStats(0);
  pytteenjet->SetTitle("t#bar{t}#rightarrow e^{+}e^{-}");
  pytteenjet->SetLineColor(2);
  pytteenjet->SetFillColor(0);
  pytteenjet->SetMarkerStyle(22);
  pytteenjet->SetMarkerColor(2);
  mgtteenjet->SetLineColor(1);
  mgtteenjet->SetFillColor(0);
  mgtteenjet->SetMarkerColor(1);
  cTTee->cd();
  pytteenjet->Draw();
  mgtteenjet->Draw("same");
  TLegend* leg = new TLegend(0.75,0.75, 0.95, 0.95);
  leg->AddEntry(pytteenjet, "PYTHIA");
  leg->AddEntry(mgtteenjet, "MadGraph");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->SaveAs(Form("%s_ttdil_hnJet_ee_%s_v_%s.eps", filePfx, pypfx, mgpfx));

  TCanvas* cTTmm = new TCanvas("cTTmm");
  TH1F* pyttmmnjet = (TH1F*)gDirectory->Get(Form("%s_ttdil_hnJet_mm", pypfx));
  TH1F* mgttmmnjet = (TH1F*)gDirectory->Get(Form("%s_ttdil_hnJet_mm", mgpfx));
  pyttmmnjet->SetStats(0);
  pyttmmnjet->SetTitle("t#bar{t}#rightarrow #mu^{+}#mu^{-}");
  pyttmmnjet->SetLineColor(2);
  pyttmmnjet->SetFillColor(0);
  pyttmmnjet->SetMarkerStyle(22);
  pyttmmnjet->SetMarkerColor(2);
  mgttmmnjet->SetLineColor(1);
  mgttmmnjet->SetFillColor(0);
  mgttmmnjet->SetMarkerColor(1);
  cTTmm->cd();
  pyttmmnjet->Draw();pyttmmnjet->SetMaximum(1.2*pyttmmnjet->GetMaximum());
  mgttmmnjet->Draw("same");
  TLegend* leg = new TLegend(0.75,0.75, 0.95, 0.95);
  leg->AddEntry(pyttmmnjet, "PYTHIA");
  leg->AddEntry(mgttmmnjet, "MadGraph");
  leg->SetFillColor(0);
  leg->Draw();
  gPad->SaveAs(Form("%s_ttdil_hnJet_mm_%s_v_%s.eps", filePfx, pypfx, mgpfx));

}
