{
gROOT->LoadMacro("eff.C");
gStyle->SetOptStat(0);

TChain *qcd = new TChain("tree");
TChain *wmunu = new TChain("tree");

TChain *ch2 = new TChain("tree");

qcd->Add("inclMu.root");
wmunu->Add("Wmunu.root");
ch2->Add("mc/Mu_576nb.root");

// Book the histograms
double xbin[13]={10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.};
int nbinsx = 12;
double xmax = 1.0;

// 10u
TH1F* num10u = new TH1F("num10u","num10u",nbinsx,xbin);
TH1F* fo10u  = new TH1F("fo10u", "fo10u", nbinsx,xbin);
// 15u
TH1F* num15u = new TH1F("num15u","num15u",nbinsx,xbin);
TH1F* fo15u  = new TH1F("fo15u", "fo15u", nbinsx,xbin);
// 20u
TH1F* num20u = new TH1F("num20u","num20u",nbinsx,xbin);
TH1F* fo20u  = new TH1F("fo20u", "fo20u", nbinsx,xbin);

// 10u
TH1F* QCDnum10u = new TH1F("QCDnum10u","QCDnum10u",nbinsx,xbin);
TH1F* QCDfo10u  = new TH1F("QCDfo10u", "QCDfo10u", nbinsx,xbin);
// 15u
TH1F* QCDnum15u = new TH1F("QCDnum15u","QCDnum15u",nbinsx,xbin);
TH1F* QCDfo15u  = new TH1F("QCDfo15u", "QCDfo15u", nbinsx,xbin);
// 20u
TH1F* QCDnum20u = new TH1F("QCDnum20u","QCDnum20u",nbinsx,xbin);
TH1F* QCDfo20u  = new TH1F("QCDfo20u", "QCDfo20u", nbinsx,xbin);

// 10u
TH1F* Wnum10u = new TH1F("Wnum10u","Wnum10u",nbinsx,xbin);
TH1F* Wfo10u  = new TH1F("Wfo10u", "Wfo10u", nbinsx,xbin);
// 15u
TH1F* Wnum15u = new TH1F("Wnum15u","Wnum15u",nbinsx,xbin);
TH1F* Wfo15u  = new TH1F("Wfo15u", "Wfo15u", nbinsx,xbin);
// 20u
TH1F* Wnum20u = new TH1F("Wnum20u","Wnum20u",nbinsx,xbin);
TH1F* Wfo20u  = new TH1F("Wfo20u", "Wfo20u", nbinsx,xbin);

// 10u
TH1F* MCnum10u = new TH1F("MCnum10u","MCnum10u",nbinsx,xbin);
TH1F* MCfo10u  = new TH1F("MCfo10u", "MCfo10u", nbinsx,xbin);
// 15u
TH1F* MCnum15u = new TH1F("MCnum15u","MCnum15u",nbinsx,xbin);
TH1F* MCfo15u  = new TH1F("MCfo15u", "MCfo15u", nbinsx,xbin);
// 20u
TH1F* MCnum20u = new TH1F("MCnum20u","MCnum20u",nbinsx,xbin);
TH1F* MCfo20u  = new TH1F("MCfo20u", "MCfo20u", nbinsx,xbin);

//// Fill the Histograms
TCanvas * ctemp = new TCanvas();

// 10u
ch2->Draw("min(pt,69.99)>>num10u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10  && num");
ch2->Draw("min(pt,69.99)>>fo10u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10 ");
// 15u
ch2->Draw("min(pt,69.99)>>num15u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15  && num");
ch2->Draw("min(pt,69.99)>>fo15u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15 ");
// 20u
ch2->Draw("min(pt,69.99)>>num20u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20  && num");
ch2->Draw("min(pt,69.99)>>fo20u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20 ");

// 10u
qcd->Draw("min(pt,69.99)>>QCDnum10u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10  && num");
qcd->Draw("min(pt,69.99)>>QCDfo10u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10 ");
// 15u
qcd->Draw("min(pt,69.99)>>QCDnum15u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15  && num");
qcd->Draw("min(pt,69.99)>>QCDfo15u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15 ");
// 20u
qcd->Draw("min(pt,69.99)>>QCDnum20u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20  && num");
qcd->Draw("min(pt,69.99)>>QCDfo20u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20 ");

// 10u
wmunu->Draw("min(pt,69.99)>>Wnum10u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10  && num");
wmunu->Draw("min(pt,69.99)>>Wfo10u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>10 ");
// 15u
wmunu->Draw("min(pt,69.99)>>Wnum15u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15  && num");
wmunu->Draw("min(pt,69.99)>>Wfo15u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>15 ");
// 20u
wmunu->Draw("min(pt,69.99)>>Wnum20u",  "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20  && num");
wmunu->Draw("min(pt,69.99)>>Wfo20u",   "abs(id)==13 && tcmet<20 && mt<25 && ptj1_b2b>15 && mu5>1 && ptj1>20 ");


delete ctemp;

// Scale MC
QCDfo15u->Scale( (576*18.205289)/1000000 );
QCDnum15u->Scale( (576*18.205289)/1000000 );

Wfo15u->Scale( (576*3.7446570)/1000000 );
Wnum15u->Scale( (576*3.7446570)/1000000 );

// Add


MCfo10u = QCDfo15u;
MCfo10u->Add(Wfo15u);

MCfo15u = QCDfo15u;
MCfo15u->Add(Wfo15u);

MCfo20u = QCDfo15u;
MCfo20u->Add(Wfo15u);


MCnum10u = QCDnum15u;
MCnum10u->Add(Wnum15u);

MCnum15u = QCDnum15u;
MCnum15u->Add(Wnum15u);

MCnum20u = QCDnum15u;
MCnum20u->Add(Wnum15u);



// Get the fake rates
TH1F* eff10u = eff(fo10u,num10u,"eff10u");
TH1F* eff15u = eff(fo15u,num15u,"eff15u");
TH1F* eff20u = eff(fo20u,num20u,"eff20u");

TH1F* MCeff10u = eff( MCfo10u, MCnum10u,"MCeff10u");
TH1F* MCeff15u = eff( MCfo15u, MCnum15u,"MCeff15u");
TH1F* MCeff20u = eff( MCfo20u, MCnum20u,"MCeff20u");

// 10u
eff10u->SetMarkerStyle(23);
eff10u->SetLineColor(kRed);
eff10u->SetMarkerColor(kRed);
// 15u
eff15u->SetMarkerStyle(23);
eff15u->SetLineColor(kBlue);
eff15u->SetMarkerColor(kBlue);
// 20u
eff20u->SetMarkerStyle(23);
eff20u->SetLineColor(kBlack);
eff20u->SetMarkerColor(kBlack);

//Plot it
eff15u->SetMinimum(0.);
eff15u->SetMaximum(xmax);
eff15u->SetTitle("Muon Fake Rate - Mu5");
eff15u->Draw();
eff10u->Draw("same");
eff20u->Draw("same");
TLegend *leg2 = new TLegend(0.5,0.7,0.99,0.99);
leg2->AddEntry(eff10u,"10 GeV Jet");
leg2->AddEntry(eff15u,"15 GeV Jet");
leg2->AddEntry(eff20u,"20 GeV Jet");
leg2->Draw();

//
TCanvas *c = new TCanvas();
c->SetWindowSize(1100,600);
c->Divide(2,2);




fo15u->SetMarkerStyle(23);
num15u->SetMarkerStyle(23);
num15u->SetMarkerColor(kRed);
QCDnum15u->SetLineColor(kRed);
Wnum15u->SetLineColor(kRed);

c->cd(1)->SetLogy();
fo15u->SetMinimum(.1);
fo15u->Draw("P");
num15u->Draw("P same");

c->cd(2)->SetLogy();
num15u->SetMinimum(.1);
//num15u->Draw("P");
QCDnum15u->Draw();
Wnum15u->Draw("same");

c->cd(3);
eff15u->Draw();
MCeff15u->Draw("same");

c->cd(4)->SetLogy();
fo15u->SetMinimum(.1);
fo15u->Draw("P");
QCDfo15u->Draw("same");
Wfo15u->Draw("same");

return;

}
