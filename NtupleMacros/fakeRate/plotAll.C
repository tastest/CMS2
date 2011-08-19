#include "plotFR.C"

void plotAll(void){
  
  gROOT->LoadMacro("eff.C");

  // Mu
  TChain *c_mu = new TChain("tree");
  c_mu->Add("august16/Mu_840nb.root");
  //c_mu->Add("Mu.root");

  // EG
  TChain *c_el = new TChain("tree");
  c_el->Add("august16/EG_840nb.root");
  //c_el->Add("EG.root");

  // Muons
  TCut mu_numCut    = "abs(id)==13 && ptj1>15 && mu9>1 && num";
  TCut mu_den04Cut  = "abs(id)==13 && ptj1>15 && mu9>1 && fo_04";
  plotFR( c_mu, mu_numCut, mu_den04Cut, "mu" );

  // V1 Electrons
  TCut elv1_numCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv1_denCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && v1";
  plotFR( c_el, elv1_numCut, elv1_denCut, "elv1" );

  // V2 Electrons
  TCut elv2_numCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv2_denCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && v2";
  plotFR( c_el, elv2_numCut, elv2_denCut, "elv2" );

  // V3 Electrons
  TCut elv3_numCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv3_denCut = "abs(id)==11 && ptj1>15 && el10_lw>1 && v3";
  plotFR( c_el, elv3_numCut, elv3_denCut, "elv3" );

  // Pt, Eta projections
  TH1F* mu_fr_pt  = eff( ((TH1F*)mu_den->ProjectionY()), ((TH1F*)mu_num->ProjectionY()), "mu_fr_pt" );
  TH1F* mu_fr_eta = eff( ((TH1F*)mu_den->ProjectionX()), ((TH1F*)mu_num->ProjectionX()), "mu_fr_eta" );
  TH1F* elv1_fr_pt  = eff( ((TH1F*)elv1_den->ProjectionY()), ((TH1F*)elv1_num->ProjectionY()), "elv1_fr_pt" );
  TH1F* elv1_fr_eta  = eff( ((TH1F*)elv1_den->ProjectionX()), ((TH1F*)elv1_num->ProjectionX()), "elv1_fr_eta" );
  TH1F* elv2_fr_pt  = eff( ((TH1F*)elv2_den->ProjectionY()), ((TH1F*)elv2_num->ProjectionY()), "elv2_fr_pt" );
  TH1F* elv2_fr_eta  = eff( ((TH1F*)elv2_den->ProjectionX()), ((TH1F*)elv2_num->ProjectionX()), "elv2_fr_eta" );
  TH1F* elv3_fr_pt  = eff( ((TH1F*)elv3_den->ProjectionY()), ((TH1F*)elv3_num->ProjectionY()), "elv3_fr_pt" );
  TH1F* elv3_fr_eta  = eff( ((TH1F*)elv3_den->ProjectionX()), ((TH1F*)elv3_num->ProjectionX()), "elv3_fr_eta" );

  // style
  gStyle->SetPaintTextFormat(".2f");
  mu_fr->SetMarkerSize(2.0);
  elv1_fr->SetMarkerSize(2.0);
  elv2_fr->SetMarkerSize(2.0);
  elv3_fr->SetMarkerSize(2.0);
  mu_fr_pt->SetLineColor(kRed);
  mu_fr_eta->SetLineColor(kRed);
  elv1_fr_pt->SetLineColor(kRed);
  elv1_fr_eta->SetLineColor(kRed);
  elv2_fr_pt->SetLineColor(kRed);
  elv2_fr_eta->SetLineColor(kRed);
  elv3_fr_pt->SetLineColor(kRed);
  elv3_fr_eta->SetLineColor(kRed);

  // plot muon FR
  TCanvas *c1 = new TCanvas();
  c1->SetWindowSize(1100,450);
  c1->Divide(3,1);
  c1->cd(1);
  mu_fr_pt->Draw();
  c1->cd(2);
  mu_fr_eta->Draw();
  c1->cd(3);
  mu_fr->Draw("TEXTE");

  // plot electron FR
  TCanvas *c2 = new TCanvas();
  c2->SetWindowSize(1100,450);
  c2->Divide(3,1);
  c2->cd(1);
  elv1_fr->Draw("TEXTE");
  c2->cd(2);
  elv2_fr->Draw("TEXTE");
  c2->cd(3);
  elv3_fr->Draw("TEXTE");

  TCanvas *c3 = new TCanvas();
  c3->SetWindowSize(1100,900);
  c3->Divide(3,2);
  c3->cd(1);
  elv1_fr_pt->Draw();
  c3->cd(2);
  elv2_fr_pt->Draw();
  c3->cd(3);
  elv3_fr_pt->Draw();
  c3->cd(4);
  elv1_fr_eta->Draw();
  c3->cd(5);
  elv2_fr_eta->Draw();
  c3->cd(6);
  elv3_fr_eta->Draw();


  
  

  return;
}
