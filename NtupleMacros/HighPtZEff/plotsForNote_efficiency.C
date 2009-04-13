
#include <vector>
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TStyle.h"
#include "TRoot.h"

void plotsForNote_efficiency()
{
	
	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle();");
	gROOT->ForceStyle();
	gStyle->SetOptTitle(1);
	
	TFile file("Results_v4.root");
	TH1F *h1_numerator = (TH1F*)file.Get("dyee_nofilter_els_eff_p_iso_numer");
	//h1_numerator->Rebin(8);
	TH1F *h1_denominator = (TH1F*)file.Get("dyee_nofilter_Gen_Z_Momentum");
	//h1_denominator->Rebin(8);

	std::vector<Double_t> binEdges;
	Int_t bin = h1_denominator->GetNbinsX() + 1;
	binEdges.push_back(h1_denominator->GetBinLowEdge(bin));
	Float_t nEntries = 0;
	Float_t n = 100;
	while (bin > 1)
	{
		bin --;
		nEntries = h1_denominator->GetBinContent(bin);
		while (nEntries < n && bin > 1)
		{
			bin --;
			nEntries += h1_denominator->GetBinContent(bin);
		}
		binEdges.push_back(h1_denominator->GetBinLowEdge(bin));
	}	
	
	Double_t *array = new Double_t[binEdges.size()];
	
	Int_t j = 0;
	for (Int_t i = binEdges.size(); i > 0; --i)
	{
		array[j] = binEdges[i - 1];
		std::cout << "array[" << j << "] = " << array[j] << std::endl;
		++j;
	}
	
	std::cout << binEdges.size() << std::endl;
	std::cout << h1_numerator->GetNbinsX() << std::endl;
	
	
	TH1F *h1_numerator_rebin = (TH1F*)h1_numerator->Rebin(binEdges.size() - 1, "", array);
	TH1F *h1_denominator_rebin = (TH1F*)h1_denominator->Rebin(binEdges.size() - 1, "", array);
	
	TGraphAsymmErrors *gr_eff_zp = new TGraphAsymmErrors();
	gr_eff_zp->SetMarkerSize(0.1);
	gr_eff_zp->BayesDivide(h1_numerator_rebin, h1_denominator_rebin);
	
	TCanvas *c0 = new TCanvas();
	c0->Divide(1, 2);
	c0->cd(1);
	h1_numerator_rebin->Draw("HIST E1");
	c0->cd(2);
	h1_denominator_rebin->Draw("HIST E1");
	
	TCanvas *c1 = new TCanvas();
	c1->cd();
	gr_eff_zp->Draw("AP");
	gr_eff_zp->GetYaxis()->SetRangeUser(0.0, 1.2);	
	
}