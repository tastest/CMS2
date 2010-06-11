{
gROOT->LoadMacro("eff.C");
gROOT->LoadMacro("histio.cc");
gStyle->SetOptStat(0);

TChain *ch1 = new TChain("tree");
ch1->Add("JMTMonitor.root");


// Book the histograms
double xbin[14]={5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,16.,18.,20.,25.};


TH1F* ptnum    = new TH1F("ptnum","ptnum",13,xbin);
TH1F* ptnumtrg = new TH1F("ptnumtrg","ptnumtrg",13,xbin);
TH1F* scptnum    = new TH1F("scptnum","scptnum",13,xbin);
TH1F* scptnumtrg = new TH1F("scptnumtrg","scptnumtrg",13,xbin);


ch1->Draw("min(pt,24.99)>>ptnum","num&&abs(id)==11");
ch1->Draw("min(scet,24.99)>>scptnum","num&&abs(id)==11");
ch1->Draw("min(pt,24.99)>>ptnumtrg","num&&eg5==2&&abs(id)==11");
ch1->Draw("min(scet,24.99)>>scptnumtrg","num&&eg5==2&&abs(id)==11");

TH1F* effptnum = eff(ptnum,ptnumtrg,"effptnum");
TH1F* effscptnum = eff(scptnum,scptnumtrg,"effscptnum");



TH1F* ptFO    = new TH1F("ptFO","ptFO",13,xbin);
TH1F* ptFOtrg = new TH1F("ptFOtrg","ptFOtrg",13,xbin);
TH1F* scptFO    = new TH1F("scptFO","scptFO",13,xbin);
TH1F* scptFOtrg = new TH1F("scptFOtrg","scptFOtrg",13,xbin);


TCanvas* ctemp = new TCanvas();
ch1->Draw("min(pt,24.99)>>ptFO","(!num)&&abs(id)==11");
ch1->Draw("min(scet,24.99)>>scptFO","!num&&abs(id)==11");
ch1->Draw("min(pt,24.99)>>ptFOtrg","(!num)&&eg5==2&&abs(id)==11");
ch1->Draw("min(scet,24.99)>>scptFOtrg","(!num)&&eg5==2&&abs(id)==11");

TH1F* effptFO = eff(ptFO,ptFOtrg,"effptFO");
TH1F* effscptFO = eff(scptFO,scptFOtrg,"effscptFO");

delete ctemp;


gStyle->SetOptStat(0);
effptnum->SetLineColor(2);
effscptnum->SetLineColor(2);
effptnum->SetMaximum(1.1);
effscptnum->SetMaximum(1.1);

TCanvas* c11 = new TCanvas();
effptnum->Draw();
effptFO->Draw("same");
TCanvas* c12 = new TCanvas();
effscptnum->Draw();
effscptFO->Draw("same");

}
