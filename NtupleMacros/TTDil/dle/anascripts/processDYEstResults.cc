
#include "processDYEstResults.h"

#include "TROOT.h"
#include "TStyle.h"

TH1F* getRHist(HistogramUtilities *hUtil, sources_t theSource, TString nJets, TString hyp_type)
{

	TH1F *h1_met_in = hUtil->getHistogram(theSource, "dyest_met_in", nJets, "_" + hyp_type);
    TH1F *h1_met_out = hUtil->getHistogram(theSource, "dyest_met_out", nJets, "_" + hyp_type);
	h1_met_in->Rebin(2);
	h1_met_out->Rebin(2);

	TH1F *h1_R = (TH1F*)h1_met_in->Clone();
	h1_R->Reset();

	for (Int_t i = 0; i < h1_met_in->GetNbinsX(); ++i)
	{
		Float_t nOut = h1_met_out->Integral(i, h1_met_out->GetNbinsX() + 1);
		Float_t nIn = h1_met_in->Integral(i, h1_met_in->GetNbinsX() + 1);
		Float_t err2_nIn = Utilities::getError2(h1_met_in, i, h1_met_in->GetNbinsX() + 1);
		Float_t err2_nOut = Utilities::getError2(h1_met_out, i, h1_met_in->GetNbinsX() + 1);

		Float_t R = 0;
		Float_t err = 0;

		if (nIn != 0 && nOut != 0) {
			R = nOut/nIn;
			err = sqrt( R*R*(err2_nOut/(nOut*nOut)) + (err2_nIn/(nIn*nIn)) );
			h1_R->SetBinContent(i, R);
			h1_R->SetBinError(i, err);
		}

	}

	h1_R->GetXaxis()->SetRangeUser(20, 100);
	delete h1_met_in;
	delete h1_met_out;
	return h1_R;

}

void get_R(HistogramUtilities *hUtil, TString hyp_type, TString nJets, Float_t &R, Float_t &err2)
{
	sources_t theSource = (1ll << H_DYEE);
	if (hyp_type == "mm") theSource = (1ll << H_DYMM);

	TH1F *h1_met_in = hUtil->getHistogram(theSource, "dyest_met_in", nJets, "_" + hyp_type);
	TH1F *h1_met_out = hUtil->getHistogram(theSource, "dyest_met_out", nJets, "_" + hyp_type);

	Float_t nOut = h1_met_out->Integral(metCutBin, h1_met_out->GetNbinsX() + 1);
	Float_t nIn = h1_met_in->Integral(metCutBin, h1_met_in->GetNbinsX() + 1);

	if (nIn != 0) R = nOut/nIn;
	else std::cout << "ERROR: nIn == 0" << std::endl;

	Float_t err2_nIn = Utilities::getError2(h1_met_in, metCutBin, h1_met_in->GetNbinsX() + 1);
	Float_t err2_nOut = Utilities::getError2(h1_met_out, metCutBin, h1_met_in->GetNbinsX() + 1);
	err2 = R*R;

	std::cout << "err2_R = " << err2 << std::endl;
	std::cout << "R = " << R << " $\\pm$ " << sqrt(err2) << " (" << nIn << ", " << nOut << ")" <<  std::endl;

	delete h1_met_in;
	delete h1_met_out;

}

void get_k(HistogramUtilities *hUtil, TString nJets, TString hyp_type, Float_t &k, Float_t &err2)
{
	Float_t dummyErr2;

	// k_eff needs to be found from the no met histogram - first do "DY ONLY"

	// get in and out (DY ONLY)
	TH1F *h1_met_in_mm = hUtil->getHistogram((1ll<<H_DYMM), "dyest_met_in", nJets, "_mm");
	TH1F *h1_met_in_ee = hUtil->getHistogram((1ll<<H_DYEE), "dyest_met_in", nJets, "_ee");

	Float_t nMC_mm = h1_met_in_mm->Integral(0, metCutBin);
	Float_t nMC_ee = h1_met_in_ee->Integral(0, metCutBin);
	Float_t k_effMC = sqrt(nMC_mm/nMC_ee);
	if (hyp_type == "ee") k_effMC = sqrt(nMC_ee/nMC_mm);
	Float_t err2_k_effMC = (1/4.0)*(1/(nMC_ee)) + ((nMC_mm)/(nMC_ee*nMC_ee));
	if (hyp_type == "ee") err2_k_effMC = (1/4.0)*(1/(nMC_mm)) + ((nMC_ee)/(nMC_mm*nMC_mm));

	// get in and out (FROM DATA)
	TH1F *h1_met_in_mm_data = hUtil->getHistogram(sources_all, "dyest_met_in", nJets, "_mm");
	TH1F *h1_met_in_ee_data = hUtil->getHistogram(sources_all, "dyest_met_in", nJets, "_ee");

	Float_t n_mm = h1_met_in_mm_data->Integral(0, metCutBin);
	Float_t n_ee = h1_met_in_ee_data->Integral(0, metCutBin);
	std::cout << n_mm << ", " << n_ee << std::endl;
	Float_t k_eff = sqrt(n_mm/n_ee);
	if (hyp_type == "ee") k_eff = sqrt(n_ee/n_mm);
	Float_t err2_k_eff = (1/4.0)*(1/(n_ee)) + ((n_mm)/(n_ee*n_ee));
	if (hyp_type == "ee") err2_k_eff = (1/4.0)*(1/(n_mm)) + ((n_ee)/(n_mm*n_mm));

	std::cout       << std::endl;
	std::cout       << "TABLE --- for comapring efficiency correction factors" << std::endl;
	std::cout       << "nJets \t& "
		<< "DY Only $\\varepsilon_{\\mu}/\\varepsilon_{e}$ \t& "
		<< "Realistic $\\varepsilon_{\\mu}/\\varepsilon_{e}$ \\\\ \\hline" << std::endl;
	std::cout       << " \\hline" << std::endl;
	std::cout       << nJets.Data() << " \t& ";
	printf("%.2f $\\pm$ %.4f \t&", k_effMC, sqrt(err2_k_effMC));
	printf("%.2f $\\pm$ %.4f \\\\ \\hline \n", k_eff, sqrt(err2_k_eff));

	delete h1_met_in_mm;
	delete h1_met_in_ee;
	delete h1_met_in_mm_data;
	delete h1_met_in_ee_data;

	// add the combinatorics factor
	k = 0.5 * k_eff;
	err2 = 0.5*err2_k_eff;

	printf("k including combinatorics = %.2f $\\pm$ %.4f \\\\ \\hline \n", k, sqrt(err2));

}


// input should be hyp_type and jet_bin

void estimate(HistogramUtilities *hUtil, JetBins_t bin, TString hyp_type, 
		ResultsHistograms &estResults, ResultsHistograms &nonpeakResults)
{

	TString nJets;
	if (bin == 0) nJets = "0j";
	if (bin == 1) nJets = "1j";
	if (bin == 2) nJets = "2j";

	//	hUtil->printInOutTruth(nJets, hyp_type);
	//	hUtil->printInOutTruth(nJets, "em");

	// values and errors
	Float_t n_in_data_em, 		err2_n_in_data_em;
	Float_t n_in_est_nonpeak_ll,	err2_n_in_est_nonpeak_ll;
	Float_t n_in_data_ll, 		err2_n_in_data_ll;
	Float_t n_in_est_peaking_ll, 	err2_n_in_est_peaking_ll;
	Float_t n_out_est_peaking_ll, 	err2_n_out_est_peaking_ll;
	Float_t n_out_peaking_ll, 	err2_n_out_peaking_ll;
	Float_t n_in_nonpeak_ll,	err2_n_in_nonpeak_ll;
	Float_t k, 			err2_k;
	Float_t R,			err2_R;

    TH1F *h1_met_in_data_em = hUtil->getHistogram(sources_all, "dyest_met_in", nJets, "_em");
    TH1F *h1_met_in_data_ll = hUtil->getHistogram(sources_all, "dyest_met_in", nJets, "_" + hyp_type);

	std::cout << "[estimate]" << std::endl;
	std::cout << "MET CUT is " << h1_met_in_data_em->GetBinLowEdge(metCutBin) << std::endl;

	n_in_data_em = h1_met_in_data_em->Integral(metCutBin, h1_met_in_data_em->GetNbinsX() + 1);
	err2_n_in_data_em = n_in_data_em;

	n_in_data_ll = h1_met_in_data_ll->Integral(metCutBin, h1_met_in_data_ll->GetNbinsX() + 1);
	err2_n_in_data_ll = n_in_data_ll;

	std::cout << n_in_data_em << ", " << n_in_data_ll << std::endl;


	get_R(hUtil, hyp_type, nJets, R, err2_R);
	get_k(hUtil, nJets, hyp_type, k, err2_k);

	// estimate the peaking background in the near-Z control region
	n_in_est_nonpeak_ll		= k*n_in_data_em;
	err2_n_in_est_nonpeak_ll        = (k*k*err2_n_in_data_em) + (n_in_data_em*n_in_data_em*err2_k);
	n_in_est_peaking_ll 		= n_in_data_ll - n_in_est_nonpeak_ll;
	err2_n_in_est_peaking_ll 	= err2_n_in_est_nonpeak_ll + err2_n_in_data_ll;

	// estimate the peaking background in the out region
	n_out_est_peaking_ll 		= n_in_est_peaking_ll * R;
	err2_n_out_est_peaking_ll 	= n_out_est_peaking_ll*n_out_est_peaking_ll*((err2_R/(R*R)) 
			+ (err2_n_in_est_peaking_ll/(n_in_est_peaking_ll*n_in_est_peaking_ll)));


	//
	//        err2 = Utilities::getError2(h_temp, zLow_, zHigh_);
	//

	// get the true peaking background in the out region


	TH1F *h1_met_out_peaking_ll = hUtil->getHistogram(sources_peaking, "dyest_met_out", nJets, "_" + hyp_type);
	n_out_peaking_ll = h1_met_out_peaking_ll->Integral(metCutBin, h1_met_out_peaking_ll->GetNbinsX() + 1);
	err2_n_out_peaking_ll = Utilities::getError2(h1_met_out_peaking_ll, metCutBin, h1_met_out_peaking_ll->GetNbinsX() + 1);

	TH1F *h1_met_in_nonpeak_ll = hUtil->getHistogram(sources_nonpeak, "dyest_met_in", nJets, "_" + hyp_type);
	n_in_nonpeak_ll = h1_met_in_nonpeak_ll->Integral(metCutBin, h1_met_in_nonpeak_ll->GetNbinsX() + 1);
	err2_n_in_nonpeak_ll = Utilities::getError2(h1_met_in_nonpeak_ll, metCutBin, h1_met_in_nonpeak_ll->GetNbinsX() + 1);


	// tables

	std::string hyp_tex_string = "ee";
	if (hyp_type == "em") hyp_tex_string = "e\\mu";
	else if (hyp_type == "mm") hyp_tex_string = "\\mu\\mu";

	//
	// TABLE comparing no-peak bg with truth
	std::cout       << std::endl;
	std::cout       << "TABLE --- for non peaking bg estimation" << std::endl;
	std::cout       << "nJets \t&"
		<< "$N_{"<<hyp_tex_string<<"}^{in}$ \t& "
		<< "Estimated $N_{" << hyp_tex_string << "}^{in~(no~peak)}$ \t& "
		<< "True $N_{"<<hyp_tex_string<<"}^{in~(no~peak)}$ \\\\ \\hline" << std::endl;
	std::cout       << " \\hline" << std::endl;
	std::cout       << nJets.Data() << " \t& ";
	printf("%.1f $\\pm$ %.1f \t&", n_in_data_ll, sqrt(err2_n_in_data_ll));
	printf("%.1f $\\pm$ %.1f \t&", n_in_est_nonpeak_ll, sqrt(err2_n_in_est_nonpeak_ll));
	printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_in_nonpeak_ll, sqrt(err2_n_in_nonpeak_ll));

	//
	// TABLE for WW NOTE (nopeak bg)
	std::cout       << std::endl;
	std::cout       << "TABLE --- for non peaking bg estimation WW NOTE" << std::endl;
	std::cout       << "$k$ \t & "
		<< "$N_{e\\mu}^{in}$ \t& "
		<< "Estimated $N_{" << hyp_tex_string << "}^{in~(no~peak)}$ \t& "
		<< "True $N_{"<<hyp_tex_string<<"}^{in~(no~peak)}$ \\\\ \\hline" << std::endl;
	std::cout       << " \\hline" << std::endl;
	printf("%.2f $\\pm$ %.4f \t&", k, sqrt(err2_k));
	printf("%.1f $\\pm$ %.1f \t&", n_in_data_em, sqrt(err2_n_in_data_em));
	printf("%.1f $\\pm$ %.1f \t&", n_in_est_nonpeak_ll, sqrt(err2_n_in_est_nonpeak_ll));
	printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_in_nonpeak_ll, sqrt(err2_n_in_nonpeak_ll));

	//
	// TABLE comparing estimates with truth
	std::cout       << std::endl;
	std::cout       << "TABLE --- for comparison of DY estimate and truth" << std::endl;
	std::cout       << "nJets \t&"
		<< "$N_{DY/ZZ}^{in~(est)}$ \t& "
		<< "$R_{out/in}$ \t& "
		<< "$N_{DY/ZZ}^{out~(est)}$ \t& "
		<< "True $N_{DY/ZZ}^{out}$ \\\\ \\hline" << std::endl;
	std::cout       << " \\hline" << std::endl;
	std::cout       << nJets.Data() << " \t& ";
	printf("%.1f $\\pm$ %.1f \t&", n_in_est_peaking_ll, sqrt(err2_n_in_est_peaking_ll));
	printf("%.2f $\\pm$ %.2f \t&", R, sqrt(err2_R));
	printf("%.1f $\\pm$ %.1f \t&", n_out_est_peaking_ll, sqrt(err2_n_out_est_peaking_ll));
	printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_out_peaking_ll, sqrt(n_out_peaking_ll));

	//
	// TABLE comparing estimates with truth for WW NOTE
	std::cout       << std::endl;
	std::cout       << "TABLE --- for comparison of DY estimate and truth for WW NOTE" << std::endl;
	std::cout 	<< "$R_{out/in}$ \t& "
		<< "$N_{ll}^{in}$ \t& "
		<< "$N_{DY/ZZ/WZ}^{in~(est)}$ \t& "
		<< "$N_{DY/ZZ/WZ}^{out~(est)}$ \t& "
		<< "True $N_{DY/ZZ/WZ}^{out}$ \\\\ \\hline" << std::endl;
	std::cout       << " \\hline" << std::endl;
	printf("%.2f $\\pm$ %.2f \t&", R, sqrt(err2_R));
	printf("%.1f $\\pm$ %.1f \t&", n_in_data_ll, sqrt(err2_n_in_data_ll));
	printf("%.1f $\\pm$ %.1f \t&", n_in_est_peaking_ll, sqrt(err2_n_in_est_peaking_ll));
	printf("%.1f $\\pm$ %.1f \t&", n_out_est_peaking_ll, sqrt(err2_n_out_est_peaking_ll));
	//        printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_out_peaking_ll, sqrt(n_out_peaking_ll));
	printf("%.1f $\\pm$ %.3f \\\\ \\hline \n", n_out_peaking_ll, err2_n_out_peaking_ll);

	// fill estimate results histograms
	estResults.add(bin, n_out_peaking_ll, sqrt(err2_n_out_peaking_ll), 
			n_out_est_peaking_ll, sqrt(err2_n_out_est_peaking_ll));

	// fill nonpeaking background estimate results histograms
	nonpeakResults.add(bin, n_in_nonpeak_ll, sqrt(err2_n_in_nonpeak_ll),
			n_in_est_nonpeak_ll, sqrt(err2_n_in_est_nonpeak_ll));

	delete h1_met_out_peaking_ll;
	delete h1_met_in_nonpeak_ll;
	delete h1_met_in_data_em;
	delete h1_met_in_data_ll;

}

void processDYEstResults(TString fileName)
{

	gROOT->ProcessLine(".L ../tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	gStyle->SetOptTitle(1);

	// filename, met cut, doWZ, lumiNorm
	std::vector<DataSource> sources;
	sources.push_back( fH_TTBAR() );
	sources.push_back( fH_DYMM() );
	sources.push_back( fH_DYEE() );
	sources.push_back( fH_ZZ() );
	HistogramUtilities *hUtil = new HistogramUtilities(fileName, sources, 1.0);


	gROOT->cd();

	ResultsHistograms estResults_ee(kYellow, "DY/ZZ/WZ (ee)", "n-Jets", "Number of Events");
	ResultsHistograms nonpeakResults_ee(kRed, "Non-peaking (ee)", "n-Jets", "Number of Events");

	ResultsHistograms estResults_mm(kYellow, "DY/ZZ/WZ (#mu#mu)", "n-Jets", "Number of Events");
	ResultsHistograms nonpeakResults_mm(kRed, "Non-peaking (#mu#mu)", "n-Jets", "Number of Events");


	estimate(hUtil, J0, "mm", estResults_mm, nonpeakResults_mm);
	estimate(hUtil, J1, "mm", estResults_mm, nonpeakResults_mm);
	estimate(hUtil, J2, "mm", estResults_mm, nonpeakResults_mm);

	estimate(hUtil, J0, "ee", estResults_ee, nonpeakResults_ee);
	estimate(hUtil, J1, "ee", estResults_ee, nonpeakResults_ee);
	estimate(hUtil, J2, "ee", estResults_ee, nonpeakResults_ee);

	TCanvas *c1 = (TCanvas*)estResults_mm.results(0.0, 25.0);
	c1->Draw();
    Utilities::saveCanvas(c1, "results/est_mm");

	TCanvas *c2 = (TCanvas*)nonpeakResults_mm.results(0.0, 10.0);
	c2->Draw();
    Utilities::saveCanvas(c2, "results/nopeak_mm");

	TCanvas *c3 = (TCanvas*)estResults_ee.results(0.0, 25.0);
	c3->Draw();
    Utilities::saveCanvas(c3, "results/est_ee");

	TCanvas *c4 = (TCanvas*)nonpeakResults_ee.results(0.0, 10.0);
	c4->Draw();
    Utilities::saveCanvas(c4, "results/nopeak_ee");

    TLegend *lg_R = new TLegend(0.2, 0.6, 0.6, 0.9);
    lg_R->SetFillColor(kWhite);
    lg_R->SetShadowColor(kWhite);
    lg_R->SetLineColor(kWhite);

    TH1F *h1_R_DY_ee = getRHist(hUtil, (1ll<<H_DYEE), "0j", "ee");
    h1_R_DY_ee->SetLineColor(kBlue);
	h1_R_DY_ee->SetMarkerColor(kBlue);
    h1_R_DY_ee->SetMarkerStyle(25);
    TH1F *h1_R_DY_mm = getRHist(hUtil, (1ll<<H_DYMM), "0j", "mm");
    h1_R_DY_mm->SetLineColor(kRed);
	h1_R_DY_mm->SetMarkerColor(kRed);
    h1_R_DY_mm->SetMarkerStyle(27);
	lg_R->AddEntry(h1_R_DY_ee, "DY->ee", "lp");
	lg_R->AddEntry(h1_R_DY_mm, "DY->mm", "lp");


    TCanvas *c5 = new TCanvas();
    c5->cd();
    h1_R_DY_ee->Draw("E1");
    h1_R_DY_mm->Draw("SAME E1");
	lg_R->Draw();
    h1_R_DY_ee->GetYaxis()->SetRangeUser(0, 1.0);
    Utilities::saveCanvas(c5, "results/R_DY");

	TH1F *h1_R_ZZ_ee = getRHist(hUtil, (1ll<<H_ZZ), "0j", "ee");
	h1_R_ZZ_ee->SetLineColor(kBlue);
	h1_R_ZZ_ee->SetMarkerColor(kBlue);
	h1_R_ZZ_ee->SetMarkerStyle(25);
	TH1F *h1_R_ZZ_mm = getRHist(hUtil, (1ll<<H_ZZ), "0j", "mm");
	h1_R_ZZ_mm->SetLineColor(kRed);
	h1_R_ZZ_mm->SetMarkerColor(kRed);
	h1_R_ZZ_mm->SetMarkerStyle(27);
	lg_R->Clear();
    lg_R->AddEntry(h1_R_ZZ_ee, "ZZ->ee", "lp");
    lg_R->AddEntry(h1_R_ZZ_mm, "ZZ->mm", "lp");

	TCanvas *c6 = new TCanvas();
	c6->cd();
	h1_R_ZZ_ee->Draw("E1");
	h1_R_ZZ_mm->Draw("SAME E1");
	lg_R->Draw();
	h1_R_ZZ_ee->GetYaxis()->SetRangeUser(0, 1.0);
    Utilities::saveCanvas(c6, "results/R_ZZ");

	delete hUtil;

}


