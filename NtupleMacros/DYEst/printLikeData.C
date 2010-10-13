#include "processDYEstResults.C"

void printLikeData() {

        gROOT->ProcessLine(".L ~/tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(1);


	AnaHist a1;
	TFile f_45("DYEstResults_ForWW_MET45.root");:

	gROOT->cd();

	THStack *st_ee = a1.getMassStack(f_45, "0j", "ee");
	TH1F *h1_ee = a1.getLikeData(f_45, "mll", "0j", "ee");
	h1_ee->Rebin(10);
	Int_t integral = h1_ee->Integral();
	std::cout << integral << std::endl;

	TH1F *h1_pseudo_ee = new TH1F("h1_pseudo_ee", "h1_pseudo_ee", 200, 0, 200);
	h1_pseudo_ee->Rebin(10);
	h1_pseudo_ee->FillRandom(h1_ee, gRandom->Poisson(integral));

	TCanvas *c1 = new TCanvas();
	c1->Divide(1, 2);
	c1->cd(1);
	st_ee->Draw("HIST");
	c1->cd(2);
	h1_pseudo_ee->Draw("E1");


}


