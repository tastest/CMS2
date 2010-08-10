{
gROOT->LoadMacro("eff.C");
gStyle->SetOptStat("");

TChain *ch1 = new TChain("tree");
TChain *ch2 = new TChain("tree");
ch2->Add("EG.root");

// Book the histograms
//double xbin[4]={10.,15.,20.,40.};
//int nbinsx = 3;

double xbin[11]={10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.};
int nbinsx = 10;

double v1_xmax = 0.2;
double v2_xmax = 0.5;
double v3_xmax = 0.75;

//10u
TH1F* num10u = new TH1F("num10u","num10u",nbinsx,xbin);
TH1F* v110u  = new TH1F("v110u", "v110u", nbinsx,xbin);
TH1F* v210u  = new TH1F("v210u", "v210u", nbinsx,xbin);
TH1F* v310u  = new TH1F("v310u", "v310u", nbinsx,xbin);

//15u
TH1F* num15u = new TH1F("num15u","num15u",nbinsx,xbin);
TH1F* v115u  = new TH1F("v115u", "v115u", nbinsx,xbin);
TH1F* v215u  = new TH1F("v215u", "v215u", nbinsx,xbin);
TH1F* v315u  = new TH1F("v315u", "v315u", nbinsx,xbin);

//20u
TH1F* num20u = new TH1F("num20u","num20u",nbinsx,xbin);
TH1F* v120u  = new TH1F("v120u", "v120u", nbinsx,xbin);
TH1F* v220u  = new TH1F("v220u", "v220u", nbinsx,xbin);
TH1F* v320u  = new TH1F("v320u", "v320u", nbinsx,xbin);

// Fill the Histograms
TCanvas * ctemp = new TCanvas();

// 10u
ch2->Draw("min(pt,59.99)>>num10u","abs(id)==11 && el10>1 && ptj1>10 && tcmet<20 && num");
ch2->Draw("min(pt,59.99)>>v110u", "abs(id)==11 && el10>1 && ptj1>10 && tcmet<20 && v1");
ch2->Draw("min(pt,59.99)>>v210u", "abs(id)==11 && el10>1 && ptj1>10 && tcmet<20 && v2");
ch2->Draw("min(pt,59.99)>>v310u", "abs(id)==11 && el10>1 && ptj1>10 && tcmet<20 && v3");

// 15u
ch2->Draw("min(pt,59.99)>>num15u","abs(id)==11 && el10>1 && ptj1>15 && tcmet<20 && num");
ch2->Draw("min(pt,59.99)>>v115u", "abs(id)==11 && el10>1 && ptj1>15 && tcmet<20 && v1");
ch2->Draw("min(pt,59.99)>>v215u", "abs(id)==11 && el10>1 && ptj1>15 && tcmet<20 && v2");
ch2->Draw("min(pt,59.99)>>v315u", "abs(id)==11 && el10>1 && ptj1>15 && tcmet<20 && v3");

// 20u
ch2->Draw("min(pt,59.99)>>num20u","abs(id)==11 && el10>1 && ptj1>20 && tcmet<20 && num");
ch2->Draw("min(pt,59.99)>>v120u", "abs(id)==11 && el10>1 && ptj1>20 && tcmet<20 && v1");
ch2->Draw("min(pt,59.99)>>v220u", "abs(id)==11 && el10>1 && ptj1>20 && tcmet<20 && v2");
ch2->Draw("min(pt,59.99)>>v320u", "abs(id)==11 && el10>1 && ptj1>20 && tcmet<20 && v3");

delete ctemp;

// Get the fake rates

//10u
TH1F* effv110u = eff(v110u,num10u,"effv110u");
TH1F* effv210u = eff(v210u,num10u,"effv210u");
TH1F* effv310u = eff(v310u,num10u,"effv310u");

//15u
TH1F* effv115u = eff(v115u,num15u,"effv115u");
TH1F* effv215u = eff(v215u,num15u,"effv215u");
TH1F* effv315u = eff(v315u,num15u,"effv315u");

//20u
TH1F* effv120u = eff(v120u,num20u,"effv120u");
TH1F* effv220u = eff(v220u,num20u,"effv220u");
TH1F* effv320u = eff(v320u,num20u,"effv320u");

// 10u in blue
effv110u->SetMarkerStyle(21);
effv210u->SetMarkerStyle(21);
effv310u->SetMarkerStyle(21);
effv110u->SetLineColor(kRed);
effv210u->SetLineColor(kRed);
effv310u->SetLineColor(kRed);
effv110u->SetMarkerColor(kRed);
effv210u->SetMarkerColor(kRed);
effv310u->SetMarkerColor(kRed);

// 15U in blue
effv115u->SetMarkerStyle(21);
effv215u->SetMarkerStyle(21);
effv315u->SetMarkerStyle(21);
effv115u->SetLineColor(kBlue);
effv215u->SetLineColor(kBlue);
effv315u->SetLineColor(kBlue);
effv115u->SetMarkerColor(kBlue);
effv215u->SetMarkerColor(kBlue);
effv315u->SetMarkerColor(kBlue);

// 15U in blue
effv120u->SetMarkerStyle(21);
effv220u->SetMarkerStyle(21);
effv320u->SetMarkerStyle(21);
effv120u->SetLineColor(kBlack);
effv220u->SetLineColor(kBlack);
effv320u->SetLineColor(kBlack);
effv120u->SetMarkerColor(kBlack);
effv220u->SetMarkerColor(kBlack);
effv320u->SetMarkerColor(kBlack);

// Plot them
TCanvas *cnv = new TCanvas();
cnv->SetWindowSize(1100,450);
cnv->Divide(3,1);
// plot v1
cnv->cd(1);
effv115u->SetTitle("V1 Fake Rates");
effv115u->SetMaximum(v1_xmax);
effv115u->SetMinimum(0.0);
effv115u->Draw();
effv110u->Draw("same");
effv120u->Draw("same");
TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
leg1->AddEntry(effv110u,"10 GeV Jet");
leg1->AddEntry(effv115u,"15 GeV Jet");
leg1->AddEntry(effv120u,"20 GeV Jet");
leg1->SetTextSize(0.06);
leg1->Draw();
// plot v2
effv215u->SetTitle("V2 Fake Rates");
effv215u->SetMaximum(v2_xmax);
effv215u->SetMinimum(0.0);
cnv->cd(2);
effv215u->Draw();
effv210u->Draw("same");
effv220u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.75,0.99,0.99);
leg2->AddEntry(effv210u,"10 GeV Jet");
leg2->AddEntry(effv215u,"15 GeV Jet");
leg2->AddEntry(effv220u,"20 GeV Jet");
leg2->SetTextSize(0.06);
leg2->Draw();
// plot v3
effv315u->SetTitle("V3 Fake Rates");
effv315u->SetMaximum(v3_xmax);
effv315u->SetMinimum(0.0);
cnv->cd(3);
effv310u->Draw();
effv315u->Draw("same");
effv320u->Draw("same");
TLegend *leg3 = new TLegend(0.5,0.75,0.99,0.99);
leg3->AddEntry(effv310u,"10 GeV Jet");
leg3->AddEntry(effv315u,"15 GeV Jet");
leg3->AddEntry(effv320u,"20 GeV Jet");
leg3->SetTextSize(0.06);
leg3->Draw();

return;

}
