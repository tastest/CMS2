
#include "ResultsHistograms.h"

ResultsHistograms::ResultsHistograms(Color_t bandColor, 
		TString title, TString titleX, TString titleY) 
{
	h1_true_        = new TH1F("h1_true_", "hist", 3, -0.5, 2.5);
        h1_true_errLow_ = new TH1F("h1_true_errLow_", "hist", 3, -0.5, 2.5);
        h1_true_errLow_ ->SetFillColor(10);
        h1_true_errHigh_= new TH1F("h1_true_errHigh_", title, 3, -0.5, 2.5);
	h1_true_errHigh_->GetXaxis()->SetTitle(titleX);
        h1_true_errHigh_->GetYaxis()->SetTitle(titleY);
        h1_true_errHigh_->SetFillColor(bandColor);
        h1_est_         = new TH1F("h1_est_", "hist", 3, -0.5, 2.5);
        h1_est_         ->SetMarkerStyle(22);

        lg_estimate_ = new TLegend(0.2, 0.6, 0.6, 0.9);
        lg_estimate_->SetFillColor(kWhite);
        lg_estimate_->SetLineColor(kWhite);
        //lg_estimate_->AddEntry(h1_true_, "True" + title, "l");
        lg_estimate_->AddEntry(h1_true_errHigh_, "True " + title + " #pm 1#sigma", "f");
        lg_estimate_->AddEntry(h1_est_, "Estimated " + title, "lp");	

}

void ResultsHistograms::add(JetBins_t bin, Float_t truth, Float_t errTruth,
				Float_t est, Float_t errEst) {
	h1_est_->SetBinContent(bin + 1, est);
        h1_est_->SetBinError(bin + 1, errEst);
        h1_true_->SetBinContent(bin + 1, truth);
        h1_true_errLow_->SetBinContent(bin + 1, truth - errTruth);
        h1_true_errHigh_->SetBinContent(bin + 1, truth + errTruth);
}

TCanvas* ResultsHistograms::results(Float_t min, Float_t max) {
	TCanvas *c1  = new TCanvas();
        c1->cd();
        h1_true_errHigh_->Draw("HIST");
        h1_true_errLow_->Draw("SAME HIST");
        //h1_true_->Draw("SAME");
        h1_est_->Draw("SAME E1");
	h1_true_errHigh_->GetYaxis()->SetRangeUser(min, max);
	lg_estimate_->Draw("SAME");
	c1->RedrawAxis();
        return c1;
}

