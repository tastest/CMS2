
#include "processDYEstResults.h"
#include "DataSource.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TH1F.h"

void get_R(HistogramUtilities *hUtil, TString hyp_type, TString nJets, Float_t &R, Float_t &err2)
{
        DataSource theSource = fH_DYEE();
        if (hyp_type == "mm") theSource = fH_DYMM();

        TH1F *h1_met_in = hUtil->getHistogram(theSource.getBit(), "met_in", nJets, hyp_type);
        TH1F *h1_met_out = hUtil->getHistogram(theSource.getBit(), "met_out", nJets, hyp_type);

        Float_t nOut = h1_met_out->Integral(hUtil->getMetCutBin(), h1_met_out->GetNbinsX() + 1);
        Float_t nIn = h1_met_in->Integral(hUtil->getMetCutBin(), h1_met_in->GetNbinsX() + 1);

        if (nIn != 0) R = nOut/nIn;
        else std::cout << "ERROR: nIn == 0" << std::endl;

        Float_t err2_nIn = Utilities::getError2(h1_met_in, hUtil->getMetCutBin(), h1_met_in->GetNbinsX() + 1);
        Float_t err2_nOut = Utilities::getError2(h1_met_out, hUtil->getMetCutBin(), h1_met_in->GetNbinsX() + 1);
//      if (nIn !=0 && nOut != 0) err2 = R*R*( (err2_nOut/(nOut*nOut)) + (err2_nIn/(nIn*nIn)) );
	err2 = R*R;

        std::cout << "err2_R = " << err2 << std::endl;
        std::cout << "R = " << R << " $\\pm$ " << sqrt(err2) << " (" << nIn << ", " << nOut << ")" <<  std::endl;

}

void get_k(TString nJets, TString hyp_type, Float_t &k, Float_t &err2)
{
	Float_t dummyErr2;

        // k_eff needs to be found from the no met histogram - first do "DY ONLY"
        HistogramUtilities *hUtil_kMC = new HistogramUtilities("DYEstResults_GetEMuEff_TRIG.root", 30.0, true, 1.0);
        HistogramUtilities *hUtil_k = new HistogramUtilities("DYEstResults_GetEMuEff_TRIG.root", 30.0, true, 0.1);

	// get in and out (DY ONLY)
        Float_t nMC_mm = hUtil_kMC->getN(true, sources_dy, nJets, "mm", dummyErr2);
        Float_t nMC_ee = hUtil_kMC->getN(true, sources_dy, nJets, "ee", dummyErr2);
        Float_t k_effMC = sqrt(nMC_mm/nMC_ee);
	if (hyp_type == "ee") k_effMC = sqrt(nMC_ee/nMC_mm);
        Float_t err2_k_effMC = (1/4.0)*(1/(nMC_ee)) + ((nMC_mm)/(nMC_ee*nMC_ee));
	if (hyp_type == "ee") err2_k_effMC = (1/4.0)*(1/(nMC_mm)) + ((nMC_ee)/(nMC_mm*nMC_mm));

        // get in and out (FROM DATA)
        Float_t n_mm = hUtil_k->getN(true, sources_all, nJets, "mm", dummyErr2);
        Float_t n_ee = hUtil_k->getN(true, sources_all, nJets, "ee", dummyErr2);
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

        delete hUtil_kMC;
	delete hUtil_k;

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

	hUtil->printInOutTruth(nJets, hyp_type);
        hUtil->printInOutTruth(nJets, "em");

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

	n_in_data_em = hUtil->getN(true, sources_all, nJets, "em", err2_n_in_data_em);
	err2_n_in_data_em = n_in_data_em;
	n_in_data_ll = hUtil->getN(true, sources_all, nJets, hyp_type, err2_n_in_data_ll);
	err2_n_in_data_ll = err2_n_in_data_ll;

	get_R(hUtil, hyp_type, nJets, R, err2_R);
	get_k(nJets, hyp_type, k, err2_k);

	// estimate the peaking background in the near-Z control region
	n_in_est_nonpeak_ll		= k*n_in_data_em;
        err2_n_in_est_nonpeak_ll        = (k*k*err2_n_in_data_em) + (n_in_data_em*n_in_data_em*err2_k);
	n_in_est_peaking_ll 		= n_in_data_ll - n_in_est_nonpeak_ll;
        err2_n_in_est_peaking_ll 	= err2_n_in_est_nonpeak_ll + err2_n_in_data_ll;

	// estimate the peaking background in the out region
	n_out_est_peaking_ll 		= n_in_est_peaking_ll * R;
        err2_n_out_est_peaking_ll 	= n_out_est_peaking_ll*n_out_est_peaking_ll*((err2_R/(R*R)) 
					+ (err2_n_in_est_peaking_ll/(n_in_est_peaking_ll*n_in_est_peaking_ll)));

	// get the true peaking background in the out region
	n_out_peaking_ll = hUtil->getN(false, sources_peaking, nJets, hyp_type, err2_n_out_peaking_ll);
	// expected error for 100pb
//	err2_n_out_peaking_ll = n_out_peaking_ll;

	// get the true non peaking background in the near-Z control region
        n_in_nonpeak_ll = hUtil->getN(true, sources_nonpeaking, nJets, hyp_type, err2_n_in_nonpeak_ll);
	// expected error for 100pb
	//err2_n_in_nonpeak_ll = n_in_nonpeak_ll;

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

}

void processDYEstResults(TString fileName)
{

        gROOT->ProcessLine(".L tdrStyle.C");
        gROOT->ProcessLine("setTDRStyle()");
        gStyle->SetOptTitle(1);

        // filename, met cut, doWZ, lumiNorm
        HistogramUtilities *hUtil = new HistogramUtilities(fileName, 20.0, true, 0.1);

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

        //TLegend *lg_estimate = new TLegend(0.2, 0.6, 0.6, 0.9);

	TCanvas *c1 = (TCanvas*)estResults_mm.results(0.0, 5.0);
	c1->Draw();

	TCanvas *c2 = (TCanvas*)nonpeakResults_mm.results(0.0, 30.0);
	c2->Draw();

        TCanvas *c3 = (TCanvas*)estResults_ee.results(0.0, 5.0);
        c3->Draw();

        TCanvas *c4 = (TCanvas*)nonpeakResults_ee.results(0.0, 30.0);
        c4->Draw();

	TH1F *h1_R_ZZ_ee = (TH1F*)hUtil->getRHist(fH_ZZ(), "0j", "ee");
		h1_R_ZZ_ee->SetLineColor(kBlue);
		h1_R_ZZ_ee->SetMarkerStyle(20);
        TH1F *h1_R_ZZ_mm = (TH1F*)hUtil->getRHist(fH_ZZ(), "0j", "mm");
		h1_R_ZZ_mm->SetLineColor(kRed);
		h1_R_ZZ_mm->SetMarkerStyle(22);
	
	TCanvas *c5 = new TCanvas();
	c5->cd();
	h1_R_ZZ_ee->Draw("E1");
	h1_R_ZZ_mm->Draw("SAME E1");

	delete hUtil;

}


