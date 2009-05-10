
#include <vector>
#include <iostream>
#include "EffH1F.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"

using namespace std;

EffH1F::EffH1F( char* name, char* title, Int_t nbinsx, Double_t xlow, Double_t xup )
{
  numer = new TH1F( Form("%s_%s", name, "numer"), Form("%s_%s", name, "numer"), nbinsx, xlow, xup );
  denom = new TH1F( Form("%s_%s", name, "denom"), Form("%s_%s", name, "denom"), nbinsx, xlow, xup );
  eff = new TH1F( Form("%s_%s", name, "eff"), Form("%s_%s", name, "eff"), nbinsx, xlow, xup );
  gr_eff = new TGraphAsymmErrors(eff);

  numer->Sumw2();
  denom->Sumw2();
  eff->Sumw2();
}

EffH1F::~EffH1F() {
  delete numer;
  delete denom;
  delete eff;
}

void EffH1F::MakeEff( const Float_t n ) {
  eff->Divide( numer, denom );

  // Dave's code
  std::vector<Double_t> binEdges;
  //vector<double> binEdges;
  Int_t bin = denom->GetNbinsX() + 1;
  binEdges.push_back(denom->GetBinLowEdge(bin));
  Float_t nEntries = 0;
  //Float_t n = 500; //moved to arg.
  while (bin > 1) {
	bin --;
	nEntries = denom->GetBinContent(bin);
	while (nEntries < n && bin > 1) {
	  bin --;
	  nEntries += denom->GetBinContent(bin);
	}
	binEdges.push_back(denom->GetBinLowEdge(bin));
  }

  Double_t *array = new Double_t[binEdges.size()];

  Int_t j = 0;
  for (Int_t i = binEdges.size(); i > 0; --i) {
	array[j] = binEdges[i - 1];
	//std::cout << "array[" << j << "] = " << array[j] << std::endl;
	++j;
  }

  //std::cout << binEdges.size() << std::endl;
  //std::cout << h1_numerator->GetNbinsX() << std::endl;

  TH1F *h1_numer_rebin = (TH1F*)numer->Rebin(binEdges.size() - 1, "", array);
  TH1F *h1_denom_rebin = (TH1F*)denom->Rebin(binEdges.size() - 1, "", array);

  //gr_eff->SetMarkerSize(0.1);
  gr_eff->BayesDivide(h1_numer_rebin, h1_denom_rebin);

  //TFile outf("Results2.root","RECREATE");
  TFile outf("Results2.root","UPDATE");
  //gr_eff->Write();
  TCanvas *c1 = new TCanvas(eff->GetName(), eff->GetName());
  //c1->cd();
  gr_eff->Draw("AP");
  gr_eff->GetYaxis()->SetRangeUser(0.7, 1.05);
  c1->Write();
  c1->SaveAs((TString)eff->GetName()+".png");
}

