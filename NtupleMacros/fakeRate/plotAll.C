#include "plotFR.C"

void plotAll(void){

  // Mu
  TChain *c_mu = new TChain("tree");
  c_mu->Add("august16/Mu_840nb.root");
  //c_mu->Add("Mu.root");

  // EG
  TChain *c_el = new TChain("tree");
  c_el->Add("august16/EG_840nb.root");
  //c_el->Add("EG.root");

  // Muons
  TCut mu_num   = "abs(id)==13 && ptj1>15 && mu9>1 && num";
  TCut mu_den04 = "abs(id)==13 && ptj1>15 && mu9>1 && fo_04";
  plotFR( c_mu, mu_num, mu_den04, "mu" );

  // V1 Electrons
  TCut elv1_num = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv1_den = "abs(id)==11 && ptj1>15 && el10_lw>1 && v1";
  plotFR( c_el, elv1_num, elv1_den, "eV1" );

  // V2 Electrons
  TCut elv2_num = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv2_den = "abs(id)==11 && ptj1>15 && el10_lw>1 && v2";
  plotFR( c_el, elv2_num, elv2_den, "eV2" );

  // V3 Electrons
  TCut elv3_num = "abs(id)==11 && ptj1>15 && el10_lw>1 && numAug9";
  TCut elv3_den = "abs(id)==11 && ptj1>15 && el10_lw>1 && v3";
  plotFR( c_el, elv3_num, elv3_den, "eV3" );

  // plot 2D
  gStyle->SetPaintTextFormat(".2f");
  mu_fr->SetMarkerSize(2.0);
  eV1_fr->SetMarkerSize(2.0);
  eV2_fr->SetMarkerSize(2.0);
  eV3_fr->SetMarkerSize(2.0);
  mu_fr->Draw("TEXTE");
  TCanvas *c1 = new TCanvas();
  c1->SetWindowSize(1100,450);
  c1->Divide(3,1);
  c1->cd(1);
  eV1_fr->Draw("TEXTE");
  c1->cd(2);
  eV2_fr->Draw("TEXTE");
  c1->cd(3);
  eV3_fr->Draw("TEXTE");

  // plot 1D projections

  return;
}
