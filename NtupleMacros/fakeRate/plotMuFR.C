{
gROOT->LoadMacro("eff.C");
gStyle->SetOptStat(0);

TChain *ch1 = new TChain("tree");
TChain *ch2 = new TChain("tree");
ch1->Add("JMTMonitor.root");
ch2->Add("JMT.root");

// Book the histograms
double xbin[3]={10.,20.,40.};
//double xbin[5]={0.,1., 1.479, 2., 2.5};
int nbinsx = 2;
double xmax = 0.5;
//double xmax = 0.9;

// 6u
TH1F* num6u = new TH1F("num6u","num6u",nbinsx,xbin);
TH1F* fo6u  = new TH1F("fo6u", "fo6u", nbinsx,xbin);
//10 u
TH1F* num10u = new TH1F("num10u","num10u",nbinsx,xbin);
TH1F* fo10u  = new TH1F("fo10u", "fo10u", nbinsx,xbin);
//15 u
TH1F* num15u = new TH1F("num15u","num15u",nbinsx,xbin);
TH1F* fo15u  = new TH1F("fo15u", "fo15u", nbinsx,xbin);
//30 u
TH1F* num30u = new TH1F("num30u","num30u",nbinsx,xbin);
TH1F* fo30u  = new TH1F("fo30u", "fo30u", nbinsx,xbin);
//50 u
TH1F* num50u = new TH1F("num50u","num50u",nbinsx,xbin);
TH1F* fo50u  = new TH1F("fo50u", "fo50u", nbinsx,xbin);
// 15u Mu5
TH1F* num15uMu5 = new TH1F("num15uMu5","num15uMu5",nbinsx,xbin);
TH1F* fo15uMu5  = new TH1F("fo15uMu5", "fo15uMu5", nbinsx,xbin);
// Mu9
TH1F* num15uMu9 = new TH1F("num15uMu9","num15uMu9",nbinsx,xbin);
TH1F* fo15uMu9  = new TH1F("fo15uMu9", "fo15uMu9", nbinsx,xbin);

//// Fill the Histograms
TCanvas * ctemp = new TCanvas();

// 6u
ch1->Draw("min(pt,39.99)>>num6u","num&&l16u>1&&abs(id)==13&&tcmet<20");
ch1->Draw("min(pt,39.99)>>fo6u",      "l16u>1&&abs(id)==13&&tcmet<20");
// 10u
ch1->Draw("min(pt,39.99)>>num10u","num&&l110u>1&&abs(id)==13&&tcmet<20");
ch1->Draw("min(pt,39.99)>>fo10u",      "l110u>1&&abs(id)==13&&tcmet<20");
// 15u
ch2->Draw("min(pt,39.99)>>num15u","num&&hlt15u>1&&abs(id)==13&&tcmet<20");
ch2->Draw("min(pt,39.99)>>fo15u",      "hlt15u>1&&abs(id)==13&&tcmet<20");
// 30u
ch2->Draw("min(pt,39.99)>>num30u","num&&hlt30u>1&&abs(id)==13&&tcmet<20");
ch2->Draw("min(pt,39.99)>>fo30u",      "hlt30u>1&&abs(id)==13&&tcmet<20");
// 50u
ch2->Draw("min(pt,39.99)>>num50u","num&&hlt50u>1&&abs(id)==13&&tcmet<20");
ch2->Draw("min(pt,39.99)>>fo50u",      "hlt50u>1&&abs(id)==13&&tcmet<20");
//
ch2->Draw( "min(pt,39.99)>>num15uMu5", "num&&hlt15u>1&&mu5==2&&abs(id)==13&&tcmet<20");
ch2->Draw( "min(pt,39.99)>>fo15uMu5", "hlt15u>1&&mu5==2&&abs(id)==13&&tcmet<20");
//
ch2->Draw( "min(pt,39.99)>>num15uMu9", "num&&hlt15u>1&&mu9==2&&abs(id)==13&&tcmet<20");
ch2->Draw( "min(pt,39.99)>>fo15uMu9", "hlt15u>1&&mu9==2&&abs(id)==13&&tcmet<20");

//// 6u
//ch1->Draw("eta>>num6u","num&&l16u>1&&abs(id)==13&&tcmet<20");
//ch1->Draw("eta>>fo6u",      "l16u>1&&abs(id)==13&&tcmet<20");
//// 10u
//ch1->Draw("eta>>num10u","num&&l110u>1&&abs(id)==13&&tcmet<20");
//ch1->Draw("eta>>fo10u",      "l110u>1&&abs(id)==13&&tcmet<20");
//// 15u
//ch2->Draw("eta>>num15u","num&&hlt15u>1&&abs(id)==13&&tcmet<20");
//ch2->Draw("eta>>fo15u",      "hlt15u>1&&abs(id)==13&&tcmet<20");
//// 30u
//ch2->Draw("eta>>num30u","num&&hlt30u>1&&abs(id)==13&&tcmet<20");
//ch2->Draw("eta>>fo30u",      "hlt30u>1&&abs(id)==13&&tcmet<20");
//// 50u
//ch2->Draw("eta>>num50u","num&&hlt50u>1&&abs(id)==13&&tcmet<20");
//ch2->Draw("eta>>fo50u",      "hlt50u>1&&abs(id)==13&&tcmet<20");
////
//ch2->Draw( "eta>>num15uMu5", "num&&hlt15u>1&&mu5==2&&abs(id)==13&&tcmet<20");
//ch2->Draw( "eta>>fo15uMu5", "hlt15u>1&&mu5==2&&abs(id)==13&&tcmet<20");
////
//ch2->Draw( "eta>>num15uMu9", "num&&hlt15u>1&&mu9==2&&abs(id)==13&&tcmet<20");
//ch2->Draw( "eta>>fo15uMu9", "hlt15u>1&&mu9==2&&abs(id)==13&&tcmet<20");


delete ctemp;


// Get the fake rates
TH1F* eff6u  = eff(fo6u, num6u, "eff6u");
TH1F* eff10u = eff(fo10u,num10u,"eff10u");
TH1F* eff15u = eff(fo15u,num15u,"eff15u");
TH1F* eff30u = eff(fo30u,num30u,"eff30u");
TH1F* eff50u = eff(fo50u,num50u,"eff50u");

// 6U in black
eff6u->SetMarkerStyle(20);
eff6u->SetLineColor(1);
eff6u->SetMarkerColor(1);

// 10U in black
eff10u->SetMarkerStyle(23);
eff10u->SetLineColor(2);
eff10u->SetMarkerColor(2);

// 15U in blue
eff15u->SetMarkerStyle(21);
eff15u->SetLineColor(4);
eff15u->SetMarkerColor(4);

// 30U in blue
eff30u->SetMarkerStyle(21);
eff30u->SetLineColor(kGreen+2);
eff30u->SetMarkerColor(kGreen+2);

// 50U in blue
eff50u->SetMarkerStyle(21);
eff50u->SetLineColor(kCyan+1);
eff50u->SetMarkerColor(kCyan+1);

//Plot it
eff6u->SetMinimum(0.);
eff6u->SetMaximum(xmax);
eff15u->SetTitle("Muon Fake Rate - HLT_Jet15U");
eff6u.Draw();
eff10u->Draw("same");
eff15u->Draw("same");
eff30u->Draw("same");
eff50u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.7,0.99,0.99);
leg2->AddEntry(eff6u,"HLT_L1Jet6U");
leg2->AddEntry(eff10u,"HLT_L1Jet10U");
leg2->AddEntry(eff15u,"HLT_Jet15U");
leg2->AddEntry(eff30u,"HLT_Jet30U");
leg2->AddEntry(eff50u,"HLT_Jet50U");
leg2->Draw();

return;

// 6U in black
num6u->SetMarkerStyle(20);
num6u->SetLineColor(1);
num6u->SetMarkerColor(1);
// 10U in black
num10u->SetMarkerStyle(23);
num10u->SetLineColor(2);
num10u->SetMarkerColor(2);
// 15U in blue
num15u->SetMarkerStyle(21);
num15u->SetLineColor(4);
num15u->SetMarkerColor(4);
// 30U in blue
num30u->SetMarkerStyle(21);
num30u->SetLineColor(kGreen+2);
num30u->SetMarkerColor(kGreen+2);
// 50U in blue
num50u->SetMarkerStyle(21);
num50u->SetLineColor(kCyan+1);
num50u->SetMarkerColor(kCyan+1);

//Plot it
TCanvas *c1 = new TCanvas();
c1->SetWindowSize(800,400);
c1->Divide(2,1);
c1->cd(1)->SetLogy();
num6u->SetMinimum(.1);
num6u->SetTitle("Muon Fake Rate");
num6u.Draw();
num10u->Draw("same");
num15u->Draw("same");
num30u->Draw("same");
//num50u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.7,0.99,0.99);
leg2->AddEntry(num6u,"HLT_L1Jet6U");
leg2->AddEntry(num10u,"HLT_L1Jet10U");
leg2->AddEntry(num15u,"HLT_Jet15U");
leg2->AddEntry(num30u,"HLT_Jet30U");
//leg2->AddEntry(num50u,"HLT_Jet50U");
leg2->Draw();


// 6U in black
fo6u->SetMarkerStyle(20);
fo6u->SetLineColor(1);
fo6u->SetMarkerColor(1);
// 10U in black
fo10u->SetMarkerStyle(23);
fo10u->SetLineColor(2);
fo10u->SetMarkerColor(2);
// 15U in blue
fo15u->SetMarkerStyle(21);
fo15u->SetLineColor(4);
fo15u->SetMarkerColor(4);
// 30U in blue
fo30u->SetMarkerStyle(21);
fo30u->SetLineColor(kGreen+2);
fo30u->SetMarkerColor(kGreen+2);
// 50U in blue
fo50u->SetMarkerStyle(21);
fo50u->SetLineColor(kCyan+1);
fo50u->SetMarkerColor(kCyan+1);

//Plot it
c1->cd(2)->SetLogy();
fo6u->SetMinimum(.1);
fo6u->SetTitle("Muon Fake Rate");
fo6u.Draw();
fo10u->Draw("same");
fo15u->Draw("same");
fo30u->Draw("same");
//fo50u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.7,0.99,0.99);
leg2->AddEntry(fo6u,"HLT_L1Jet6U");
leg2->AddEntry(fo10u,"HLT_L1Jet10U");
leg2->AddEntry(fo15u,"HLT_Jet15U");
leg2->AddEntry(fo30u,"HLT_Jet30U");
//leg2->AddEntry(fo50u,"HLT_Jet50U");
leg2->Draw();


// 15u Mu5
TH1F* eff15uMu5 = eff( fo15uMu5, num15uMu5, "eff15uMu5");

// Mu9
TH1F* eff15uMu9 = eff( fo15uMu9, num15uMu9, "eff15uMu9");

new TCanvas();
eff15u->Draw();
eff15uMu5->SetMarkerStyle(22);
eff15uMu5->SetMarkerColor(kRed);
eff15uMu5->SetLineColor(kRed);
eff15uMu5->Draw("same");
eff15uMu9->SetMarkerStyle(20);
//eff15uMu9->Draw("same");

TLegend *l1 = new TLegend(0.40,0.8,0.99,0.99);
l1->AddEntry(eff15u,"No HLT_Mu5 trigger required");
l1->AddEntry(eff15uMu5,"HLT_Mu5 trigger required");
//l1->AddEntry(eff15uMu9,"Mu9 trig required");
l1->Draw();
return;


}
