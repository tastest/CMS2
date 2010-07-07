{
gROOT->LoadMacro("eff.C");
gStyle->SetOptStat("");

TChain *ch1 = new TChain("tree");
TChain *ch2 = new TChain("tree");
ch1->Add("JMTMonitor.root");
ch2->Add("JMT.root");

// Book the histograms
double xbin[4]={10.,15.,20.,40.};
//double xbin[5]={0.,1., 1.479, 2., 2.5};
int nbinsx = 3;

double v1_xmax = 0.2;
double v2_xmax = 0.5;
double v3_xmax = 0.75;

// 6u
TH1F* num6u = new TH1F("num6u","num6u",nbinsx,xbin);
TH1F* v16u  = new TH1F("v16u", "v16u", nbinsx,xbin);
TH1F* v26u  = new TH1F("v26u", "v26u", nbinsx,xbin);
TH1F* v36u  = new TH1F("v36u", "v36u", nbinsx,xbin);
//10 u
TH1F* num10u = new TH1F("num10u","num10u",nbinsx,xbin);
TH1F* v110u  = new TH1F("v110u", "v110u", nbinsx,xbin);
TH1F* v210u  = new TH1F("v210u", "v210u", nbinsx,xbin);
TH1F* v310u  = new TH1F("v310u", "v310u", nbinsx,xbin);
//15 u
TH1F* num15u = new TH1F("num15u","num15u",nbinsx,xbin);
TH1F* v115u  = new TH1F("v115u", "v115u", nbinsx,xbin);
TH1F* v215u  = new TH1F("v215u", "v215u", nbinsx,xbin);
TH1F* v315u  = new TH1F("v315u", "v315u", nbinsx,xbin);
//30 u
TH1F* num30u = new TH1F("num30u","num30u",nbinsx,xbin);
TH1F* v130u  = new TH1F("v130u", "v130u", nbinsx,xbin);
TH1F* v230u  = new TH1F("v230u", "v230u", nbinsx,xbin);
TH1F* v330u  = new TH1F("v330u", "v330u", nbinsx,xbin);
//50 u
TH1F* num50u = new TH1F("num50u","num50u",nbinsx,xbin);
TH1F* v150u  = new TH1F("v150u", "v150u", nbinsx,xbin);
TH1F* v250u  = new TH1F("v250u", "v250u", nbinsx,xbin);
TH1F* v350u  = new TH1F("v350u", "v350u", nbinsx,xbin);

// Fill the Histograms
TCanvas * ctemp = new TCanvas();
// 6u
ch1->Draw("min(pt,39.99)>>num6u","num&&l16u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v16u","v1&&l16u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v26u","v2&&l16u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v36u","v3&&l16u>1&&abs(id)==11&&tcmet<20");
// 10u
ch1->Draw("min(pt,39.99)>>num10u","num&&l110u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v110u","v1&&l110u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v210u","v2&&l110u>1&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v310u","v3&&l110u>1&&abs(id)==11&&tcmet<20");
// 15u
ch2->Draw("min(pt,39.99)>>num15u","num&&hlt15u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v115u","v1&&hlt15u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v215u","v2&&hlt15u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v315u","v3&&hlt15u>1&&abs(id)==11&&tcmet<20");
// 30u
ch2->Draw("min(pt,39.99)>>num30u","num&&hlt30u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v130u","v1&&hlt30u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v230u","v2&&hlt30u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v330u","v3&&hlt30u>1&&abs(id)==11&&tcmet<20");
// 50u
ch2->Draw("min(pt,39.99)>>num50u","num&&hlt50u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v150u","v1&&hlt50u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v250u","v2&&hlt50u>1&&abs(id)==11&&tcmet<20");
ch2->Draw("min(pt,39.99)>>v350u","v3&&hlt50u>1&&abs(id)==11&&tcmet<20");


//ch1->Draw("eta>>num6u","num&&l16u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v16u","v1&&l16u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v26u","v2&&l16u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v36u","v3&&l16u>1&&abs(id)==11&&tcmet<20");
//// 10u
//ch1->Draw("eta>>num10u","num&&l110u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v110u","v1&&l110u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v210u","v2&&l110u>1&&abs(id)==11&&tcmet<20");
//ch1->Draw("eta>>v310u","v3&&l110u>1&&abs(id)==11&&tcmet<20");
//// 15u
//ch2->Draw("eta>>num15u","num&&hlt15u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v115u","v1&&hlt15u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v215u","v2&&hlt15u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v315u","v3&&hlt15u>1&&abs(id)==11&&tcmet<20");
//// 30u
//ch2->Draw("eta>>num30u","num&&hlt30u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v130u","v1&&hlt30u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v230u","v2&&hlt30u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v330u","v3&&hlt30u>1&&abs(id)==11&&tcmet<20");
//// 50u
//ch2->Draw("eta>>num50u","num&&hlt50u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v150u","v1&&hlt50u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v250u","v2&&hlt50u>1&&abs(id)==11&&tcmet<20");
//ch2->Draw("eta>>v350u","v3&&hlt50u>1&&abs(id)==11&&tcmet<20");


delete ctemp;

// Get the fake rates
// 6u
TH1F* effv16u = eff(v16u,num6u,"effv16u");
TH1F* effv26u = eff(v26u,num6u,"effv26u");
TH1F* effv36u = eff(v36u,num6u,"effv36u");
//10u
TH1F* effv110u = eff(v110u,num10u,"effv110u");
TH1F* effv210u = eff(v210u,num10u,"effv210u");
TH1F* effv310u = eff(v310u,num10u,"effv310u");
//15u
TH1F* effv115u = eff(v115u,num15u,"effv115u");
TH1F* effv215u = eff(v215u,num15u,"effv215u");
TH1F* effv315u = eff(v315u,num15u,"effv315u");
//30u
TH1F* effv130u = eff(v130u,num30u,"effv130u");
TH1F* effv230u = eff(v230u,num30u,"effv230u");
TH1F* effv330u = eff(v330u,num30u,"effv330u");
//50u
TH1F* effv150u = eff(v150u,num50u,"effv150u");
TH1F* effv250u = eff(v250u,num50u,"effv250u");
TH1F* effv350u = eff(v350u,num50u,"effv350u");

// 6U in black
effv16u->SetMarkerStyle(20);
effv26u->SetMarkerStyle(20);
effv36u->SetMarkerStyle(20);
effv16u->SetLineColor(1);
effv26u->SetLineColor(1);
effv36u->SetLineColor(1);
effv16u->SetMarkerColor(1);
effv26u->SetMarkerColor(1);
effv36u->SetMarkerColor(1);
// 10U in red
effv110u->SetMarkerStyle(23);
effv210u->SetMarkerStyle(23);
effv310u->SetMarkerStyle(23);
effv110u->SetLineColor(2);
effv210u->SetLineColor(2);
effv310u->SetLineColor(2);
effv110u->SetMarkerColor(2);
effv210u->SetMarkerColor(2);
effv310u->SetMarkerColor(2);
// 15U in blue
effv115u->SetMarkerStyle(21);
effv215u->SetMarkerStyle(21);
effv315u->SetMarkerStyle(21);
effv115u->SetLineColor(4);
effv215u->SetLineColor(4);
effv315u->SetLineColor(4);
effv115u->SetMarkerColor(4);
effv215u->SetMarkerColor(4);
effv315u->SetMarkerColor(4);
// 30U in blue
effv130u->SetMarkerStyle(21);
effv230u->SetMarkerStyle(21);
effv330u->SetMarkerStyle(21);
effv130u->SetLineColor(kGreen+2);
effv230u->SetLineColor(kGreen+2);
effv330u->SetLineColor(kGreen+2);
effv130u->SetMarkerColor(kGreen+2);
effv230u->SetMarkerColor(kGreen+2);
effv330u->SetMarkerColor(kGreen+2);
// 50U in blue
effv150u->SetMarkerStyle(21);
effv250u->SetMarkerStyle(21);
effv350u->SetMarkerStyle(21);
effv150u->SetLineColor(kCyan+1);
effv250u->SetLineColor(kCyan+1);
effv350u->SetLineColor(kCyan+1);
effv150u->SetMarkerColor(kCyan+1);
effv250u->SetMarkerColor(kCyan+1);
effv350u->SetMarkerColor(kCyan+1);

// Plot them
TCanvas *cnv = new TCanvas();
cnv->SetWindowSize(1100,450);
cnv->Divide(3,1);
// plot v1
cnv->cd(1);
effv110u->SetTitle("V1 Fake Rates");
effv110u->SetMaximum(v1_xmax);
effv110u->SetMinimum(0.0);
effv110u->Draw();
effv16u->Draw("same");
effv115u->Draw("same");
effv130u->Draw("same");
effv150u->Draw("same");
TLegend *leg1 = new TLegend(0.5,0.75,0.99,0.99);
leg1->AddEntry(effv16u,"HLT_L1Jet6U");
leg1->AddEntry(effv110u,"HLT_L1Jet10U");
leg1->AddEntry(effv115u,"HLT_Jet15U");
leg1->AddEntry(effv130u,"HLT_Jet30U");
leg1->AddEntry(effv150u,"HLT_Jet50U");
leg1->SetTextSize(0.06);
leg1->Draw();
// plot v2
effv210u->SetTitle("V2 Fake Rates");
effv210u->SetMaximum(v2_xmax);
effv210u->SetMinimum(0.0);
cnv->cd(2);
effv210u->Draw();
effv26u->Draw("same");
effv215u->Draw("same");
effv230u->Draw("same");
effv250u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.75,0.99,0.99);
leg2->AddEntry(effv26u,"HLT_L1Jet6U");
leg2->AddEntry(effv210u,"HLT_L1Jet10U");
leg2->AddEntry(effv215u,"HLT_Jet15U");
leg2->AddEntry(effv230u,"HLT_Jet30U");
leg2->AddEntry(effv250u,"HLT_Jet50U");
leg2->SetTextSize(0.06);
leg2->Draw();
// plot v3
effv310u->SetTitle("V3 Fake Rates");
effv310u->SetMaximum(v3_xmax);
effv310u->SetMinimum(0.0);
cnv->cd(3);
effv310u->Draw();
effv36u->Draw("same");
effv315u->Draw("same");
effv330u->Draw("same");
effv350u->Draw("same");
TLegend *leg3 = new TLegend(0.5,0.75,0.99,0.99);
leg3->AddEntry(effv36u,"HLT_L1Jet6U");
leg3->AddEntry(effv310u,"HLT_L1Jet10U");
leg3->AddEntry(effv315u,"HLT_Jet15U");
leg3->AddEntry(effv330u,"HLT_Jet30U");
leg3->AddEntry(effv350u,"HLT_Jet50U");
leg3->SetTextSize(0.06);
leg3->Draw();

return;

// num
num6u->SetMarkerStyle(20);
num6u->SetLineColor(1);
num6u->SetMarkerColor(1);

num10u->SetMarkerStyle(21);
num10u->SetLineColor(2);
num10u->SetMarkerColor(2);

num15u->SetMarkerStyle(22);
num15u->SetLineColor(4);
num15u->SetMarkerColor(4);

num30u->SetMarkerStyle(22);
num30u->SetLineColor(kGreen+2);
num30u->SetMarkerColor(kGreen+2);

num50u->SetMarkerStyle(22);
num50u->SetLineColor(kCyan+1);
num50u->SetMarkerColor(kCyan+1);

// 6U in black
v16u->SetMarkerStyle(20);
v26u->SetMarkerStyle(20);
v36u->SetMarkerStyle(20);
v16u->SetLineColor(1);
v26u->SetLineColor(1);
v36u->SetLineColor(1);
v16u->SetMarkerColor(1);
v26u->SetMarkerColor(1);
v36u->SetMarkerColor(1);
// 10U in red
v110u->SetMarkerStyle(23);
v210u->SetMarkerStyle(23);
v310u->SetMarkerStyle(23);
v110u->SetLineColor(2);
v210u->SetLineColor(2);
v310u->SetLineColor(2);
v110u->SetMarkerColor(2);
v210u->SetMarkerColor(2);
v310u->SetMarkerColor(2);
// 15U in blue
v115u->SetMarkerStyle(21);
v215u->SetMarkerStyle(21);
v315u->SetMarkerStyle(21);
v115u->SetLineColor(4);
v215u->SetLineColor(4);
v315u->SetLineColor(4);
v115u->SetMarkerColor(4);
v215u->SetMarkerColor(4);
v315u->SetMarkerColor(4);
// 30U in blue
v130u->SetMarkerStyle(21);
v230u->SetMarkerStyle(21);
v330u->SetMarkerStyle(21);
v130u->SetLineColor(kGreen+2);
v230u->SetLineColor(kGreen+2);
v330u->SetLineColor(kGreen+2);
v130u->SetMarkerColor(kGreen+2);
v230u->SetMarkerColor(kGreen+2);
v330u->SetMarkerColor(kGreen+2);
// 50U in blue
v150u->SetMarkerStyle(21);
v250u->SetMarkerStyle(21);
v350u->SetMarkerStyle(21);
v150u->SetLineColor(kCyan+1);
v250u->SetLineColor(kCyan+1);
v350u->SetLineColor(kCyan+1);
v150u->SetMarkerColor(kCyan+1);
v250u->SetMarkerColor(kCyan+1);
v350u->SetMarkerColor(kCyan+1);

// Plot them
TCanvas *cnv = new TCanvas();
cnv->SetWindowSize(1100,650);
cnv->Divide(3,2);

cnv->cd(2)->SetLogy();
num6u->SetTitle("Numerator");
num6u->SetMinimum(0.1);
num6u->Draw();
num10u->Draw("same");
num15u->Draw("same");
num30u->Draw("same");

// plot v1
cnv->cd(4)->SetLogy();
v110u->SetTitle("V1 Fake Rates");
v110u->SetMinimum(0.1);
v110u->Draw();
v16u->Draw("same");
v115u->Draw("same");
v130u->Draw("same");
//v150u->Draw("same");
TLegend *leg1 = new TLegend(0.5,0.15,0.99,0.5);
leg1->AddEntry(v16u,"HLT_L1Jet6U");
leg1->AddEntry(v110u,"HLT_L1Jet10U");
leg1->AddEntry(v115u,"HLT_Jet15U");
leg1->AddEntry(v130u,"HLT_Jet30U");
//leg1->AddEntry(v150u,"HLT_Jet50U");
leg1->SetTextSize(0.06);
leg1->Draw();
// plot v2
v210u->SetTitle("V2 Fake Rates");
v210u->SetMinimum(0.1);
cnv->cd(5)->SetLogy();
v210u->Draw();
v26u->Draw("same");
v215u->Draw("same");
v230u->Draw("same");
//v250u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.15,0.99,0.5);
leg2->AddEntry(v26u,"HLT_L1Jet6U");
leg2->AddEntry(v210u,"HLT_L1Jet10U");
leg2->AddEntry(v215u,"HLT_Jet15U");
leg2->AddEntry(v230u,"HLT_Jet30U");
//leg2->AddEntry(v250u,"HLT_Jet50U");
leg2->SetTextSize(0.06);
leg2->Draw();
// plot v3
v310u->SetTitle("V3 Fake Rates");
v310u->SetMinimum(0.1);
cnv->cd(6)->SetLogy();
v310u->Draw();
v36u->Draw("same");
v315u->Draw("same");
v330u->Draw("same");
//v350u->Draw("same");
TLegend *leg3 = new TLegend(0.5,0.15,0.99,0.5);
leg3->AddEntry(v36u,"HLT_L1Jet6U");
leg3->AddEntry(v310u,"HLT_L1Jet10U");
leg3->AddEntry(v315u,"HLT_Jet15U");
leg3->AddEntry(v330u,"HLT_Jet30U");
//leg3->AddEntry(v350u,"HLT_Jet50U");
leg3->SetTextSize(0.06);
leg3->Draw();



// Now egamma trigger study on 15U
// Photon 10 first
TH1F* num15uph10 = new TH1F("num15uph10","num15uph10",3,xbin);
TH1F* v115uph10  = new TH1F("v115uph10", "v115uph10", 3,xbin);
TH1F* v215uph10  = new TH1F("v215uph10", "v215uph10", 3,xbin);
TH1F* v315uph10  = new TH1F("v315uph10", "v315uph10", 3,xbin);

TCanvas * ctemp = new TCanvas();
ch1->Draw("min(pt,39.99)>>num15uph10","num&&hlt15u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v115uph10","v1&&hlt15u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v215uph10","v2&&hlt15u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v315uph10","v3&&hlt15u>1&&ph10==2&&abs(id)==11&&tcmet<20");
delete ctemp;

TH1F* effv115uph10 = eff(v115uph10,num15uph10,"effv115uph10");
TH1F* effv215uph10 = eff(v215uph10,num15uph10,"effv215uph10");
TH1F* effv315uph10 = eff(v315uph10,num15uph10,"effv315uph10");

effv115uph10->SetLineColor(8);
effv215uph10->SetLineColor(8);
effv315uph10->SetLineColor(8);
effv115uph10->SetMarkerStyle(29);
effv215uph10->SetMarkerStyle(29);
effv315uph10->SetMarkerStyle(29);
effv115uph10->SetMarkerColor(8);
effv215uph10->SetMarkerColor(8);
effv315uph10->SetMarkerColor(8);

TCanvas *cnvph10 = new TCanvas();
cnvph10->SetWindowSize(1100,450);
cnvph10->Divide(3,1);
cnvph10->cd(1);
effv115u->SetMaximum(0.1);
effv115u->Draw();
effv115uph10->Draw("same");
TLegend *leg1ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg1ph10->AddEntry(effv115u,"No EG trig required");
leg1ph10->AddEntry(effv115uph10,"ph10 required");
leg1ph10->SetTextSize(0.05);
leg1ph10->Draw();

cnvph10->cd(2);
effv215u->SetMaximum(0.3);
effv215u->Draw();
effv215uph10->Draw("same");
TLegend *leg2ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg2ph10->AddEntry(effv215u,"No EG trig required");
leg2ph10->AddEntry(effv215uph10,"ph10 required");
leg2ph10->SetTextSize(0.05);
leg2ph10->Draw();

cnvph10->cd(3);
effv315u->SetMaximum(0.6);
effv315u->Draw();
effv315uph10->Draw("same");
TLegend *leg3ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg3ph10->AddEntry(effv315u,"No EG trig required");
leg3ph10->AddEntry(effv315uph10,"ph10 required");
leg3ph10->SetTextSize(0.05);
leg3ph10->Draw();

return;

// Now egamma trigger study on 10U
// Photon 10 first
TH1F* num10uph10 = new TH1F("num10uph10","num10uph10",3,xbin);
TH1F* v110uph10  = new TH1F("v110uph10", "v110uph10", 3,xbin);
TH1F* v210uph10  = new TH1F("v210uph10", "v210uph10", 3,xbin);
TH1F* v310uph10  = new TH1F("v310uph10", "v310uph10", 3,xbin);

TCanvas * ctemp = new TCanvas();
ch1->Draw("min(pt,39.99)>>num10uph10","num&&l110u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v110uph10","v1&&l110u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v210uph10","v2&&l110u>1&&ph10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v310uph10","v3&&l110u>1&&ph10==2&&abs(id)==11&&tcmet<20");
delete ctemp;

TH1F* effv110uph10 = eff(v110uph10,num10uph10,"effv110uph10");
TH1F* effv210uph10 = eff(v210uph10,num10uph10,"effv210uph10");
TH1F* effv310uph10 = eff(v310uph10,num10uph10,"effv310uph10");

effv110uph10->SetLineColor(8);
effv210uph10->SetLineColor(8);
effv310uph10->SetLineColor(8);
effv110uph10->SetMarkerStyle(29);
effv210uph10->SetMarkerStyle(29);
effv310uph10->SetMarkerStyle(29);
effv110uph10->SetMarkerColor(8);
effv210uph10->SetMarkerColor(8);
effv310uph10->SetMarkerColor(8);

TCanvas *cnvph10 = new TCanvas();
cnvph10->Divide(3,1);
cnvph10->cd(1);
effv210u->SetMaximum(0.3);
effv110u->Draw();
effv110uph10->Draw("same");
TLegend *leg1ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg1ph10->AddEntry(effv110u,"No EG trig required");
leg1ph10->AddEntry(effv110uph10,"ph10 required");
leg1ph10->SetTextSize(0.05);
leg1ph10->Draw();

cnvph10->cd(2);
effv210u->SetMaximum(0.3);
effv210u->Draw();
effv210uph10->Draw("same");
TLegend *leg2ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg2ph10->AddEntry(effv210u,"No EG trig required");
leg2ph10->AddEntry(effv210uph10,"ph10 required");
leg2ph10->SetTextSize(0.05);
leg2ph10->Draw();

cnvph10->cd(3);
effv310u->SetMaximum(0.8);
effv310u->Draw();
effv310uph10->Draw("same");
TLegend *leg3ph10 = new TLegend(0.40,0.8,0.99,0.99);
leg3ph10->AddEntry(effv310u,"No EG trig required");
leg3ph10->AddEntry(effv310uph10,"ph10 required");
leg3ph10->SetTextSize(0.05);
leg3ph10->Draw();

return;

// Now el10
TH1F* num10uel10 = new TH1F("num10uel10","num10uel10",3,xbin);
TH1F* v110uel10  = new TH1F("v110uel10", "v110uel10", 3,xbin);
TH1F* v210uel10  = new TH1F("v210uel10", "v210uel10", 3,xbin);
TH1F* v310uel10  = new TH1F("v310uel10", "v310uel10", 3,xbin);

TCanvas * ctemp = new TCanvas();
ch1->Draw("min(pt,39.99)>>num10uel10","num&&l110u>1&&el10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v110uel10","v1&&l110u>1&&el10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v210uel10","v2&&l110u>1&&el10==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v310uel10","v3&&l110u>1&&el10==2&&abs(id)==11&&tcmet<20");
delete ctemp;

TH1F* effv110uel10 = eff(v110uel10,num10uel10,"effv110uel10&&abs(id)==11&&tcmet<20");
TH1F* effv210uel10 = eff(v210uel10,num10uel10,"effv210uel10&&abs(id)==11&&tcmet<20");
TH1F* effv310uel10 = eff(v310uel10,num10uel10,"effv310uel10&&abs(id)==11&&tcmet<20");

effv110uel10->SetLineColor(8);
effv210uel10->SetLineColor(8);
effv310uel10->SetLineColor(8);
effv110uel10->SetMarkerStyle(29);
effv210uel10->SetMarkerStyle(29);
effv310uel10->SetMarkerStyle(29);
effv110uel10->SetMarkerColor(8);
effv210uel10->SetMarkerColor(8);
effv310uel10->SetMarkerColor(8);

TCanvas *cnvel10 = new TCanvas();
cnvel10->Divide(3,1);
cnvel10->cd(1);
effv110u->Draw();
effv110uel10->Draw("same");
TLegend *leg1el10 = new TLegend(0.40,0.8,0.99,0.99);
leg1el10->AddEntry(effv110u,"No EG trig required");
leg1el10->AddEntry(effv110uel10,"el10 required");
leg1el10->SetTextSize(0.05);
leg1el10->Draw();

cnvel10->cd(2);
effv210u->Draw();
effv210uel10->Draw("same");
TLegend *leg2el10 = new TLegend(0.40,0.8,0.99,0.99);
leg2el10->AddEntry(effv210u,"No EG trig required");
leg2el10->AddEntry(effv210uel10,"el10 required");
leg2el10->SetTextSize(0.05);
leg2el10->Draw();

cnvel10->cd(3);
effv310u->Draw();
effv310uel10->Draw("same");
TLegend *leg3el10 = new TLegend(0.40,0.8,0.99,0.99);
leg3el10->AddEntry(effv310u,"No EG trig required");
leg3el10->AddEntry(effv310uel10,"el10 required");
leg3el10->SetTextSize(0.05);
leg3el10->Draw();

// Noe egamma 5
TH1F* num10ueg5 = new TH1F("num10ueg5","num10ueg5",3,xbin);
TH1F* v110ueg5  = new TH1F("v110ueg5", "v110ueg5", 3,xbin);
TH1F* v210ueg5  = new TH1F("v210ueg5", "v210ueg5", 3,xbin);
TH1F* v310ueg5  = new TH1F("v310ueg5", "v310ueg5", 3,xbin);

TCanvas * ctemp = new TCanvas();
ch1->Draw("min(pt,39.99)>>num10ueg5","num&&l110u>1&&eg5==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v110ueg5","v1&&l110u>1&&eg5==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v210ueg5","v2&&l110u>1&&eg5==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v310ueg5","v3&&l110u>1&&eg5==2&&abs(id)==11&&tcmet<20");
delete ctemp;

TH1F* effv110ueg5 = eff(v110ueg5,num10ueg5,"effv110ueg5");
TH1F* effv210ueg5 = eff(v210ueg5,num10ueg5,"effv210ueg5");
TH1F* effv310ueg5 = eff(v310ueg5,num10ueg5,"effv310ueg5");

effv110ueg5->SetLineColor(8);
effv210ueg5->SetLineColor(8);
effv310ueg5->SetLineColor(8);
effv110ueg5->SetMarkerStyle(29);
effv210ueg5->SetMarkerStyle(29);
effv310ueg5->SetMarkerStyle(29);
effv110ueg5->SetMarkerColor(8);
effv210ueg5->SetMarkerColor(8);
effv310ueg5->SetMarkerColor(8);

TCanvas *cnveg5 = new TCanvas();
cnveg5->Divide(3,1);
cnveg5->cd(1);
effv110u->Draw();
effv110ueg5->Draw("same");
TLegend *leg1eg5 = new TLegend(0.40,0.8,0.99,0.99);
leg1eg5->AddEntry(effv110u,"No EG trig required");
leg1eg5->AddEntry(effv110ueg5,"eg5 required");
leg1eg5->SetTextSize(0.05);
leg1eg5->Draw();

cnveg5->cd(2);
effv210u->Draw();
effv210ueg5->Draw("same");
TLegend *leg2eg5 = new TLegend(0.40,0.8,0.99,0.99);
leg2eg5->AddEntry(effv210u,"No EG trig required");
leg2eg5->AddEntry(effv210ueg5,"eg5 required");
leg2eg5->SetTextSize(0.05);
leg2eg5->Draw();

cnveg5->cd(3);
effv310u->Draw();
effv310ueg5->Draw("same");
TLegend *leg3eg5 = new TLegend(0.40,0.8,0.99,0.99);
leg3eg5->AddEntry(effv310u,"No EG trig required");
leg3eg5->AddEntry(effv310ueg5,"eg5 required");
leg3eg5->SetTextSize(0.05);
leg3eg5->Draw();

// Noe egamma 8
TH1F* num10ueg8 = new TH1F("num10ueg8","num10ueg8",3,xbin);
TH1F* v110ueg8  = new TH1F("v110ueg8", "v110ueg8", 3,xbin);
TH1F* v210ueg8  = new TH1F("v210ueg8", "v210ueg8", 3,xbin);
TH1F* v310ueg8  = new TH1F("v310ueg8", "v310ueg8", 3,xbin);

TCanvas * ctemp = new TCanvas();
ch1->Draw("min(pt,39.99)>>num10ueg8","num&&l110u>1&&eg8==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v110ueg8","v1&&l110u>1&&eg8==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v210ueg8","v2&&l110u>1&&eg8==2&&abs(id)==11&&tcmet<20");
ch1->Draw("min(pt,39.99)>>v310ueg8","v3&&l110u>1&&eg8==2&&abs(id)==11&&tcmet<20");
delete ctemp;

TH1F* effv110ueg8 = eff(v110ueg8,num10ueg8,"effv110ueg8");
TH1F* effv210ueg8 = eff(v210ueg8,num10ueg8,"effv210ueg8");
TH1F* effv310ueg8 = eff(v310ueg8,num10ueg8,"effv310ueg8");

effv110ueg8->SetLineColor(8);
effv210ueg8->SetLineColor(8);
effv310ueg8->SetLineColor(8);
effv110ueg8->SetMarkerStyle(29);
effv210ueg8->SetMarkerStyle(29);
effv310ueg8->SetMarkerStyle(29);
effv110ueg8->SetMarkerColor(8);
effv210ueg8->SetMarkerColor(8);
effv310ueg8->SetMarkerColor(8);

TCanvas *cnveg8 = new TCanvas();
cnveg8->Divide(3,1);
cnveg8->cd(1);
effv110u->Draw();
effv110ueg8->Draw("same");
TLegend *leg1eg8 = new TLegend(0.40,0.8,0.99,0.99);
leg1eg8->AddEntry(effv110u,"No EG trig required");
leg1eg8->AddEntry(effv110ueg8,"eg8 required");
leg1eg8->SetTextSize(0.05);
leg1eg8->Draw();

cnveg8->cd(2);
effv210u->Draw();
effv210ueg8->Draw("same");
TLegend *leg2eg8 = new TLegend(0.40,0.8,0.99,0.99);
leg2eg8->AddEntry(effv210u,"No EG trig required");
leg2eg8->AddEntry(effv210ueg8,"eg8 required");
leg2eg8->SetTextSize(0.05);
leg2eg8->Draw();

cnveg8->cd(3);
effv310u->Draw();
effv310ueg8->Draw("same");
TLegend *leg3eg8 = new TLegend(0.40,0.8,0.99,0.99);
leg3eg8->AddEntry(effv310u,"No EG trig required");
leg3eg8->AddEntry(effv310ueg8,"eg8 required");
leg3eg8->SetTextSize(0.05);
leg3eg8->Draw();


}
