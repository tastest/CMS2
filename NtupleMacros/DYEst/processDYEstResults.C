
#include "TStyle.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"

#include <iostream>

#include "Utilities.h"

class AnaHist {

	public:
		AnaHist() {
                        sources[0] = "ww";
                        sources[1] = "ttbar";
			sources[2] = "dymm";
			sources[3] = "dyee";
			sources[4] = "dytt";
			sources[5] = "wjets";
			sources[6] = "tw";
			sources[7] = "zz";
			sources[8] = "wz";
			lumiNorm_ = 0.1;
			doWZ_ = true;
			doEfficiencyCorrection_ = true;
		}
		~AnaHist() {}

		Double_t getLumiNorm() { return lumiNorm_; }
		bool doEfficiencyCorrection() { return doEfficiencyCorrection_; }

		void setWZ(bool doWZ) {doWZ_ = doWZ;}
		bool doWZ() {return doWZ_;}

                void getInOutTruth(TFile &file, TString var, TString nJets, TString hyp_type) {
                        TString histNameSuffix = "_" + var + "_" + nJets + "-N-1_" + hyp_type;
			std::cout << std::endl;
			std::cout << "TRUTH for " << nJets.Data() << " " << hyp_type.Data() << std::endl;
			std::cout << "\t\t in \t\t out" << std::endl;
			Float_t totalIn = 0;
			Float_t totalOut = 0;
			Float_t totalIn_nopeak = 0;
			Float_t totalOut_nopeak = 0;
                        for (Int_t i = 0; i < 9; ++i)
                        {
				if (!doWZ_ && sources[i] == "wz") continue;
				TH1F *h_temp = ((TH1F*)file.Get(sources[i] + histNameSuffix));
        	                Int_t zLow = h_temp->FindBin(76.1);
	                        Int_t zHigh = h_temp->FindBin(105.9);
                                Float_t nIn = h_temp->Integral(zLow, zHigh);
				Float_t nOut = h_temp->Integral(0, h_temp->GetNbinsX() + 1) - nIn;
				totalIn += nIn;
				totalOut += nOut;
				if (sources[i] != "zz" && sources[i] != "dyee" && sources[i] != "dymm") {
					totalIn_nopeak += nIn;
					totalOut_nopeak += nOut;
				}
				std::cout << sources[i].Data() << "\t & \t ";
                                printf("%.2f \t & %.2f \t \\\\ \\hline \n", nIn*lumiNorm_, nOut*lumiNorm_);

                        } 
			
			std::cout << "TOTAL \t & \t";
                        printf("%.2f \t & %.2f \t \\\\ \\hline \n", totalIn*lumiNorm_, totalOut*lumiNorm_);
                        std::cout << "NOPEAK \t & \t";                         
			printf("%.2f \t & %.2f \t \\\\ \\hline \n", totalIn_nopeak*lumiNorm_, totalOut_nopeak*lumiNorm_);
			std::cout << std::endl;
                }
		
		Float_t getNOut(TFile &file, TString sourceName, TString nJets, TString hyp_type, Float_t &error2 = 0)
		{
              		TH1F *h_temp = (TH1F*)((TH1F*)file.Get(sourceName + "_mll_" + nJets + "-N-1_" + hyp_type))->Clone();
			h_temp->Scale(lumiNorm_);
                	Int_t zLow = h_temp->FindBin(76.1);
                	Int_t zHigh = h_temp->FindBin(105.9);
                	Float_t nIn = h_temp->Integral(zLow, zHigh);
                	Float_t nOut = h_temp->Integral(0, h_temp->GetNbinsX() + 1) - nIn;

			// error
			error2 = getError2(h_temp, 0, zLow) 
					+ getError2(h_temp, zHigh, h_temp->GetNbinsX() + 1);


			delete h_temp;
			return nOut;
		}

		TH1F *getLikeData(TFile &file, TString var, TString nJets, TString hyp_type) {
			TString histNameSuffix = "_" + var + "_" + nJets + "-N-1_" + hyp_type;
			TH1F *h1_data = (TH1F*)((file.Get(sources[0] + histNameSuffix))->Clone());			
			for (Int_t i = 1; i < 9; ++i)
			{
				//std::cout << sources[i] + histNameSuffix << std::endl;
				if (!doWZ_ && sources[i] == "wz") continue;
				h1_data->Add((TH1F*)file.Get(sources[i] + histNameSuffix));
			}

			h1_data->Scale(lumiNorm_);
			return h1_data;
		}

                TH1F *getNoPeakTruth(TFile &file, TString var, TString nJets, TString hyp_type) {
                        TString histNameSuffix = "_" + var + "_" + nJets + "-N-1_" + hyp_type;
                        TH1F *h1_data = (TH1F*)(file.Get(sources[0] + histNameSuffix)->Clone());
                        for (Int_t i = 1; i < 9; ++i)
                        {
				//if (!doWZ_ && sources[i] == "wz") continue;
                                if(sources[i] != "wz" &&
					sources[i] != "zz" && 
					sources[i] != "dyee" && 
					sources[i] != "dymm") 
					h1_data->Add((TH1F*)file.Get(sources[i] + histNameSuffix));
                        }

                        h1_data->Scale(lumiNorm_);
                        return h1_data;
                }


		TH1F *getMC(TFile &file, TString name) {
			return (TH1F*)(file.Get(name)->Clone());
		}

                TLegend *getLegend(TFile &file) {
                        TLegend *lg = new TLegend(0.7, 0.5, 0.9, 0.9);
			lg->SetFillColor(kWhite);
			lg->SetLineColor(kWhite);
                        TString histNameSuffix = "_mll_0j-N-1_mm";
                        for (Int_t i = 0; i < 9; ++i)
                        {
                                TH1F *h_temp = (TH1F*)file.Get(sources[i] + histNameSuffix)->Clone();
                                lg->AddEntry(h_temp, sources[i].Data(), "f");
                        }
                        return lg;
                }


                THStack *getMetStack(TFile &file, TString nJets, TString region, TString hyp_type) {
                        THStack *st_temp = new THStack("st_temp", "");
                        TString histNameSuffix = "_met_" + nJets + "_" + region + "-N-1_" + hyp_type;
                        //std::cout << histNameSuffix.Data() << std::endl;
                        for (Int_t i = 8; i >= 0; --i)
                        {
                                //TString temp = sources[i] + histNameSuffix;
                                //std::cout << "getting " << temp.Data() << std::endl;
                                TH1F *h1_temp = ((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5);
                                h1_temp->Scale(lumiNorm_);
                                st_temp->Add(h1_temp);
                                //st_temp->Add(((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5));
                                 
                        }
                        return st_temp;
                }

		THStack *getMassStack(TFile &file, TString nJets, TString hyp_type) {
			THStack *st_temp = new THStack("st_temp", "");
                        TString histNameSuffix = "_mll_" + nJets + "-N-1_" + hyp_type;
			//std::cout << histNameSuffix.Data() << std::endl;
                        for (Int_t i = 8; i >= 0; --i)
			{
				//TString temp = sources[i] + histNameSuffix;
				//std::cout << "getting " << temp.Data() << std::endl;
				TH1F *h1_temp = ((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5);
				h1_temp->Scale(lumiNorm_);
                                st_temp->Add(h1_temp);
                                //st_temp->Add(((TH1F*)(file.Get(sources[i] + histNameSuffix)->Clone()))->Rebin(5));
				
                        }
                        return st_temp;
		}

	private:
		TString sources[9];
		Double_t lumiNorm_;
		bool doWZ_;
		bool doEfficiencyCorrection_;
};

// do the actual estimate
//
void estimateDY(TFile &file, TString nJets, TString hyp_type, 
                Float_t &truth, Float_t &estimate, Float_t &error, Float_t &truth_error,
                Float_t &noPeak_estimate, Float_t &noPeak_error, Float_t &noPeak_truth, Float_t &noPeak_truth_error)
{

	// Utility class to get stuff
	AnaHist anaHist;
	anaHist.getInOutTruth(file, "mll", nJets, hyp_type);
        anaHist.getInOutTruth(file, "mll", nJets, "em");

	// ingrediants

	// ****** R_{out/in} for this jet bin and met cuts
	TString dyType = "dyee";
	if (hyp_type == "mm") dyType = "dymm";
        TH1F *h1_met_in = anaHist.getMC(file, dyType + "_met_" + nJets + "_in-N-1_" + hyp_type);
        TH1F *h1_met_out = anaHist.getMC(file, dyType + "_met_" + nJets + "_out-N-1_" + hyp_type);
        Int_t metCutBin = h1_met_out->FindBin(20.0);  
        Float_t nOut = h1_met_out->Integral(metCutBin, h1_met_out->GetNbinsX() + 1);
        Float_t nIn = h1_met_in->Integral(metCutBin, h1_met_in->GetNbinsX() + 1);


	Float_t R = 0;
	if (nIn != 0) R = nOut/nIn;
	Float_t err2_R = 0;
	Float_t err2_nIn = getError2(h1_met_in, metCutBin, h1_met_in->GetNbinsX() + 1);
        Float_t err2_nOut = getError2(h1_met_out, metCutBin, h1_met_in->GetNbinsX() + 1);
//	if (nIn !=0 && nOut != 0) err2_R = R*R*( (err2_nOut/(nOut*nOut)) + (err2_nIn/(nIn*nIn)) );
err2_R = R*R;

	else std::cout << "ERROR: nIn == 0" << std::endl;
	std::cout << "err2_R = " << err2_R << std::endl;
	std::cout << "R = " << R << " $\\pm$ " << sqrt(err2_R) << " (" << nIn << ", " << nOut << ")" <<  std::endl;

	// ******* n_in_ll_data
        TH1F *h1_mll_ll = anaHist.getMC(file, dyType + "_mll_" + nJets + "-N-1_" + hyp_type);
        Int_t zLow = h1_mll_ll->FindBin(76.1);
        Int_t zHigh = h1_mll_ll->FindBin(105.9);	
	TH1F *h1_mll_data = anaHist.getLikeData(file, "mll", nJets, hyp_type);	
	Float_t n_in_data_ll = h1_mll_data->Integral(zLow, zHigh);
	Float_t err2_n_in_data_ll = n_in_data_ll;
	std::cout << "n_in_data_ll = " << n_in_data_ll << " $\\pm$ " << sqrt(err2_n_in_data_ll) << std::endl;
	
        // ******* n_in_emu_data
        TH1F *h1_emu_data = anaHist.getLikeData(file, "mll", nJets, "em");
        Float_t n_in_data_emu = h1_emu_data->Integral(zLow, zHigh);
	Float_t err2_n_in_data_emu = n_in_data_emu;
        std::cout << "n_in_data_emu = " <<  n_in_data_emu << " $\\pm$ " << sqrt(err2_n_in_data_emu) << std::endl;

        // ******* k including correction for e/mu eff

	Float_t k_comb = 0.5;	// expect 0.5 * number of events in mm or ee as in em
	Float_t k = k_comb;	// need to also take into account differing eff for e and mu

	// k_eff needs to be found from the no met histogram - first do "DY ONLY"
	TFile file_noMetCuts("DYEstResults_GetEMuEff.root", "READ");
        TH1F *h1_mllNoMetMC_mm = anaHist.getMC(file_noMetCuts, "dymm_mll_" + nJets + "-N-1_mm");
        TH1F *h1_mllNoMetMC_ee = anaHist.getMC(file_noMetCuts, "dyee_mll_" + nJets + "-N-1_ee");
	Float_t nMC_mm = h1_mllNoMetMC_mm->Integral(zLow, zHigh);
	Float_t nMC_ee = h1_mllNoMetMC_ee->Integral(zLow, zHigh);
	// now do like "FROM DATA"
	TH1F *h1_mllNoMet_mm = anaHist.getLikeData(file_noMetCuts, "mll", nJets, "mm");
        TH1F *h1_mllNoMet_ee = anaHist.getLikeData(file_noMetCuts, "mll", nJets, "ee");
        Float_t n_mm = h1_mllNoMet_mm->Integral(zLow, zHigh);
        Float_t n_ee = h1_mllNoMet_ee->Integral(zLow, zHigh);
        file_noMetCuts.Close();

	std::cout << "computing k_eff (n_mm/n_ee) for 1fb MC " << nMC_mm << "/" << nMC_ee << std::endl;
	Float_t k_effMC = sqrt(nMC_mm/nMC_ee);
	Float_t err2_k_effMC = (1/(4*nMC_ee)) + ((4*nMC_mm)/(nMC_ee*nMC_ee));
	std::cout << "\t" << k_effMC << " $\\pm$ " << sqrt(err2_k_effMC) << std::endl;

	std::cout 	<< "computing k_eff from data from data... (lumiNorm_ = " << anaHist.getLumiNorm() << ") " 
			<< n_mm << "/" << n_ee << std::endl; 
        Float_t k_eff = sqrt(n_mm/n_ee);
	Float_t err2_k_eff = (1/(4*n_ee)) + ((4*n_mm)/(n_ee*n_ee));
	std::cout 	<< "\t" << k_eff << " $\\pm$ " << sqrt(err2_k_eff) << std::endl;

	// if the efficiency correction is to be applied then scale the 
	// k correction by this. If not then assign a systematic error 
	// to it that is consistent with the "meausred" value of the 
	// efficiency correction.
	bool doEfficiencyCorrection = anaHist.doEfficiencyCorrection();
	if (doEfficiencyCorrection && hyp_type == "mm") k *= k_eff;
	if (doEfficiencyCorrection && hyp_type == "ee") k /= k_eff;
	if (doEfficiencyCorrection) Float_t err2_k = k*k*err2_k_eff;
	else Float_t err2_k = k*k*k_eff*k_eff*(k_eff - 1)*(k_eff - 1);
	std::cout << "k = " << k << " $\\pm$ " << sqrt(err2_k) << std::endl;

        // ******* n_in_ll_dy_zz
	
	Float_t n_in_ll_dy_zz = n_in_data_ll - (k * n_in_data_emu);
	Float_t err2_n_in_ll_nopeak = (k*k*err2_n_in_data_emu) + (n_in_data_emu*n_in_data_emu*err2_k);
	Float_t err2_n_in_ll_dy_zz = err2_n_in_ll_nopeak + err2_n_in_data_ll;// / (n_in_ll_dy_zz*n_in_ll_dy_zz);
	std::cout << "n_in_ll_dy_zz " << n_in_ll_dy_zz << " $\\pm$ " << sqrt(err2_n_in_ll_dy_zz) << std::endl;	

        // ******* n_out_dy_zz
	
	Float_t n_out_ll_dy_zz = n_in_ll_dy_zz * R;
	//Float_t err2_n_out_ll_dy_zz = err2_R + err2_n_in_ll_dy_zz;
std::cout << nJets.Data() << ", " << R << ", " << n_in_ll_dy_zz << std::endl;
	Float_t err2_n_out_ll_dy_zz = n_out_ll_dy_zz*n_out_ll_dy_zz*((err2_R/(R*R)) + (err2_n_in_ll_dy_zz/(n_in_ll_dy_zz*n_in_ll_dy_zz)));
	std::cout << "n_out_ll_dy_zz = " << n_out_ll_dy_zz << " $\\pm$ " << sqrt(err2_n_out_ll_dy_zz) << std::endl;	

	// ******* compare to truth
	Float_t tmp_err2 = 0;
	Float_t err2_n_out_true_dy_zz = 0;
        Float_t n_out_true_dy_zz = anaHist.getNOut(file, "dy" + hyp_type, nJets, hyp_type, tmp_err2);
		//err2_n_out_true_dy_zz += tmp_err2;
        n_out_true_dy_zz += anaHist.getNOut(file, "zz", nJets, hyp_type, tmp_err2);
        if (anaHist.doWZ()) n_out_true_dy_zz += anaHist.getNOut(file, "wz", nJets, hyp_type, tmp_err2);

		//err2_n_out_true_dy_zz += tmp_err2;
	err2_n_out_true_dy_zz = sqrt(n_out_true_dy_zz);
        std::cout << "n_out_true_dy_zz = " << n_out_true_dy_zz << std::endl;

	// finish
	estimate = n_out_ll_dy_zz;
	error = sqrt(err2_n_out_ll_dy_zz);
	truth = n_out_true_dy_zz;
	truth_error = sqrt(err2_n_out_true_dy_zz);

	std::cout << "REPORT" << std::endl;
	// 
	// Do stuff for final report/tables
	TH1F *h1_noPeakTruth = anaHist.getNoPeakTruth(file, "mll", nJets, hyp_type);
        Float_t n_in_data_ll_nopeak = h1_noPeakTruth->Integral(zLow, zHigh);
        Float_t err2_n_in_data_ll_nopeak = n_in_data_ll_nopeak;
	std::string hyp_tex_string = "ee";
	if (hyp_type == "em") hyp_tex_string = "e\\mu";
	else if (hyp_type == "mm") hyp_tex_string = "\\mu\\mu";
	std::cout.setf(ios::fixed);
	std::cout.precision(4);

        noPeak_estimate = k*n_in_data_emu;
        noPeak_error = sqrt(err2_n_in_ll_nopeak);
        noPeak_truth = n_in_data_ll_nopeak;
        noPeak_truth_error = sqrt(n_in_data_ll_nopeak);

	//
	// TABLE comparing no-peak bg with truth
	std::cout 	<< std::endl;
	std::cout 	<< "TABLE --- for non peaking bg estimation" << std::endl;
	std::cout	<< "nJets \t&"
			<< "$N_{"<<hyp_tex_string<<"}^{in}$ \t& "    
			<< "Estimated $N_{\mu\mu}^{in~(no~peak)}$ \t& "
                        << "True $N_{"<<hyp_tex_string<<"}^{in~(no~peak)}$ \\\\ \\hline" << std::endl;
	std::cout 	<< " \\hline" << std::endl;
        std::cout       << nJets.Data() << " \t& ";
        printf("%.1f $\\pm$ %.1f \t&", n_in_data_ll, sqrt(err2_n_in_data_ll));
        printf("%.1f $\\pm$ %.1f \t&", k*n_in_data_emu, sqrt(err2_n_in_ll_nopeak));
        printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_in_data_ll_nopeak, sqrt(err2_n_in_data_ll_nopeak));
//	std::cout 	<< nJets.Data() << " \t& " 
//			<< n_in_data_ll << " $\\pm$ " << sqrt(err2_n_in_data_ll) << " \t& "
//                      << k*n_in_data_emu << " $\\pm$ " << sqrt(err2_n_in_ll_nopeak) << " \t& "
//			<< n_in_data_ll_nopeak << " $\\pm$ " << sqrt(err2_n_in_data_ll_nopeak) << " \\\\ \\hline " << std::endl;

	//
	// TABLE comparing estimates with truth
	std::cout 	<< std::endl;
	std::cout 	<< "TABLE --- for comparison of DY estimate and truth" << std::endl;
        std::cout       << "nJets \t&"
                        << "$N_{DY/ZZ}^{in~(est)}$ \t& "
                        << "$R_{out/in}$ \t& "
                        << "$N_{DY/ZZ}^{out~(est)}$ \t& "
			<< "True $N_{DY/ZZ}^{out}$ \\\\ \\hline" << std::endl;
        std::cout       << " \\hline" << std::endl;
	std::cout 	<< nJets.Data() << " \t& ";
        printf("%.1f $\\pm$ %.1f \t&", n_in_ll_dy_zz, sqrt(err2_n_in_ll_dy_zz));
        printf("%.2f $\\pm$ %.2f \t&", R, sqrt(err2_R));
        printf("%.1f $\\pm$ %.1f \t&", n_out_ll_dy_zz, sqrt(err2_n_out_ll_dy_zz));
        printf("%.1f $\\pm$ %.1f \\\\ \\hline \n", n_out_true_dy_zz, sqrt(err2_n_out_true_dy_zz));
//			<< n_in_ll_dy_zz << " $ \\pm $ " << sqrt(err2_n_in_ll_dy_zz) << " \t& "
//			<< R << " $\\pm$ " << sqrt(err2_R) << " \t& "
//			<< n_out_ll_dy_zz << " $\\pm$ " << sqrt(err2_n_out_ll_dy_zz) << " \t& "
//			<< n_out_true_dy_zz << " $\\pm$ " << sqrt(err2_n_out_true_dy_zz) << " \\\\ \\hline " << std::endl;
	//
	// TABLE comparing efficiency correction factors
	std::cout 	<< std::endl;
	std::cout	<< "TABLE --- for comapring efficiency correction factors" << std::endl;
	std::cout 	<< "nJets \t& "
			<< "DY Only $\\varepsilon_{\\mu}/\\varepsilon_{e}$ \t& "
			<< "Realistic $\\varepsilon_{\\mu}/\\varepsilon_{e}$ \\\\ \\hline" << std::endl;
	std::cout 	<< " \\hline" << std::endl;
	std::cout 	<< nJets.Data() << " \t& ";
        printf("%.2f $\\pm$ %.2f \t&", k_effMC, sqrt(err2_k_effMC));
        printf("%.2f $\\pm$ %.2f \\\\ \\hline \n", k_eff, sqrt(err2_k_eff));
//			<< k_effMC << " $\\pm$ " << sqrt(err2_k_effMC) << " \t& "
//			<< k_eff << " $\\pm$ " << sqrt(err2_k_eff) << " \\\\ \\hline" << std::endl;

/*
	delete h1_mllNoMet_mm;
	delete h1_mllNoMet_ee;
	delete h1_emu_data;
	delete h1_mll_data;
	delete h1_met_in;
	delete h1_met_out;
	delete h1_mll_ll;
*/

}

void printRForZZNoMET(TFile &f, TString hyp_type)
{
        TH1F *h1_met_in_mm = (TH1F*)(f.Get("zz_met_2j_in-N-1_" + hyp_type))->Clone();
        TH1F *h1_met_out_mm = (TH1F*)(f.Get("zz_met_2j_out-N-1_" + hyp_type))->Clone();

	Float_t weight = h1_met_in_mm->GetSumOfWeights() / h1_met_in_mm->GetEntries();
	std::cout << weight << std::endl;

	h1_met_in_mm->Scale(1/weight);
	h1_met_out_mm->Scale(1/weight);

	// no met
	Float_t nOut = h1_met_out_mm->Integral(0, h1_met_out_mm->GetNbinsX() + 1);
        Float_t nIn = h1_met_in_mm->Integral(0, h1_met_in_mm->GetNbinsX() + 1);
        Float_t R = 0;
        if (nOut != 0) R = nOut/nIn;
        Float_t err = 0;
        if (nOut != 0) err = sqrt( (1/nIn) + (1/nOut) );

        std::cout << "R = " << R << " $\\pm$ " << R*err << " (" << nIn << ", " << nOut << ")" <<  std::endl;

}

TH1F *getRHist(TFile &f, TString source, TString nJets, TString hyp_type, bool weighting)
{

	gROOT->cd();

        TH1F *h1_met_in = (TH1F*)(f.Get(source + "_met_" + nJets + "_in-N-1_" + hyp_type))->Clone();
        TH1F *h1_met_out = (TH1F*)(f.Get(source + "_met_" + nJets + "_out-N-1_" + hyp_type))->Clone();
        h1_met_in->Rebin(2);
        h1_met_out->Rebin(2);

	// un do the weighting - only works if all weights are the same
	// and is probably irrelevant now that I wrote the getError2 function...
	// FIXME FIXME FIXME
	if (!weighting)
	{
        	Float_t weight = h1_met_in->GetSumOfWeights() / h1_met_in->GetEntries();
        	h1_met_in->Scale(1/weight);
        	h1_met_out->Scale(1/weight);
	}

        TH1F h1_R("h1_R", "h1_R", h1_met_in->GetNbinsX(), 0, 200);
	Int_t minBin = h1_met_in->FindBin(20.0);	
        for (Int_t i = minBin; i < h1_met_in->GetNbinsX(); ++i)
        {
                Float_t nOut = h1_met_out->Integral(i, h1_met_out->GetNbinsX() + 1);
                Float_t nIn = h1_met_in->Integral(i, h1_met_in->GetNbinsX() + 1);
	        Float_t err2_nIn = getError2(h1_met_in, i, h1_met_in->GetNbinsX() + 1);
        	Float_t err2_nOut = getError2(h1_met_out, i, h1_met_in->GetNbinsX() + 1);

                Float_t R = 0;
                if (nIn != 0) R = nOut/nIn;
                Float_t err = 0;
                //if (nOut != 0) err = sqrt( (1/nIn) + (1/nOut) );
		if (nIn != 0 && nOut != 0) err = sqrt( R*R*(err2_nOut/(nOut*nOut)) + (err2_nIn/(nIn*nIn)) );
                h1_R.SetBinContent(i, R);
                h1_R.SetBinError(i, err);

        }

	h1_R.GetXaxis()->SetRangeUser(20, 100);	
	delete h1_met_in;
	delete h1_met_out;
	return (TH1F*)h1_R.Clone();

}

// main function
//
void processDYEstResults()
{

        gROOT->ProcessLine(".L ../Tools/histtools.cc");
	gROOT->ProcessLine(".L ~/tdrStyle.C");
	gROOT->ProcessLine("setTDRStyle()");
	gStyle->SetOptTitle(1);

        TFile f("DYEstResults_ForWW_MET45_INCL.root", "READ");
//	TFile f("DYEstResults_ForWW_MET45.root", "READ");
      //TFile f("DYEstResults_ForTTBAR.root", "READ");

	gROOT->cd();

	// now compare the estimates

        Float_t estimate = 0;
        Float_t error = 0; 
        Float_t truth = 0;
        Float_t truth_error = 0;

        Float_t noPeak_estimate = 0;
        Float_t noPeak_error = 0; 
        Float_t noPeak_truth = 0;
        Float_t noPeak_truth_error = 0;

        TH1F *h1_estimate_mm = new TH1F("h1_estimate_mm", "h1_estimate_mm", 3, -0.5, 2.5);
        TH1F *h1_estimate_ee = new TH1F("h1_estimate_ee", "h1_estimate_ee", 3, -0.5, 2.5);

        TH1F *h1_true_mm = new TH1F("h1_true_mm", "h1_true_mm", 3, -0.5, 2.5);
                h1_true_mm->SetFillColor(kWhite); 
                h1_true_mm->SetBinError(1, 0);
                h1_true_mm->SetBinError(2, 0);
                h1_true_mm->SetBinError(3, 0);
        TH1F *h1_true_mm_errDown = new TH1F("h1_true_mm_errDown", "h1_true_mm_errDown", 3, -0.5, 2.5);
                h1_true_mm_errDown->SetLineColor(10);
                h1_true_mm_errDown->SetFillColor(10);
                h1_true_mm_errDown->SetFillStyle(1001);
        TH1F *h1_true_mm_errUp = new TH1F("h1_true_mm_errUp", "h1_true_mm_errUp", 3, -0.5, 2.5);
                h1_true_mm_errUp->SetLineColor(10);
                h1_true_mm_errUp->SetFillColor(kYellow);
                h1_true_mm_errUp->GetXaxis()->SetTitle("n-Jets");
                h1_true_mm_errUp->GetYaxis()->SetTitle("Number of Events");

        TH1F *h1_true_ee = new TH1F("h1_true_ee", "h1_true_ee", 3, -0.5, 2.5);
                h1_true_ee->SetFillColor(kWhite);
                h1_true_ee->SetBinError(1, 0);
                h1_true_ee->SetBinError(2, 0);
                h1_true_ee->SetBinError(3, 0);
        TH1F *h1_true_ee_errDown = new TH1F("h1_true_ee_errDown", "h1_true_ee_errDown", 3, -0.5, 2.5);
                h1_true_ee_errDown->SetLineColor(10);
                h1_true_ee_errDown->SetFillColor(10);
                h1_true_ee_errDown->SetFillStyle(1001);
        TH1F *h1_true_ee_errUp = new TH1F("h1_true_ee_errUp", "h1_true_ee_errUp", 3, -0.5, 2.5);
                h1_true_ee_errUp->SetLineColor(10);
                h1_true_ee_errUp->SetFillColor(kYellow);
                h1_true_ee_errUp->GetXaxis()->SetTitle("n-Jets");
                h1_true_ee_errUp->GetYaxis()->SetTitle("Number of Events");

        TLegend *lg_estimate = new TLegend(0.2, 0.6, 0.6, 0.9);
        lg_estimate->SetFillColor(kWhite);
        lg_estimate->SetLineColor(kWhite);
        lg_estimate->AddEntry(h1_true_ee, "True DY/ZZ", "l");
        lg_estimate->AddEntry(h1_true_ee_errUp, "True DY/ZZ #pm 1#sigma", "f");
        lg_estimate->AddEntry(h1_estimate_ee, "Estimated DY/ZZ", "lp");

        TH1F *h1_estimateNoPeak_mm = new TH1F("h1_estimateNoPeak_mm", "h1_estimateNoPeak_mm", 3, -0.5, 2.5);
        TH1F *h1_estimateNoPeak_ee = new TH1F("h1_estimateNoPeak_ee", "h1_estimateNoPeak_ee", 3, -0.5, 2.5);

        TH1F *h1_trueNoPeak_mm = new TH1F("h1_trueNoPeak_mm", "h1_trueNoPeak_mm", 3, -0.5, 2.5);
                h1_trueNoPeak_mm->SetFillColor(kWhite);
                h1_trueNoPeak_mm->SetBinError(1, 0);
                h1_trueNoPeak_mm->SetBinError(2, 0);
                h1_trueNoPeak_mm->SetBinError(3, 0);
        TH1F *h1_trueNoPeak_mm_errDown = new TH1F("h1_trueNoPeak_mm_errDown", "h1_trueNoPeak_mm_errDown", 3, -0.5, 2.5);
                h1_trueNoPeak_mm_errDown->SetLineColor(10);
                h1_trueNoPeak_mm_errDown->SetFillColor(10);
                h1_trueNoPeak_mm_errDown->SetFillStyle(1001);
        TH1F *h1_trueNoPeak_mm_errUp = new TH1F("h1_trueNoPeak_mm_errUp", "h1_trueNoPeak_mm_errUp", 3, -0.5, 2.5);
                h1_trueNoPeak_mm_errUp->SetLineColor(10);
                h1_trueNoPeak_mm_errUp->SetFillColor(kRed);
                h1_trueNoPeak_mm_errUp->GetXaxis()->SetTitle("n-Jets");
                h1_trueNoPeak_mm_errUp->GetYaxis()->SetTitle("Number of Events");

        TH1F *h1_trueNoPeak_ee = new TH1F("h1_trueNoPeak_ee", "h1_trueNoPeak_ee", 3, -0.5, 2.5);
                h1_trueNoPeak_ee->SetFillColor(kWhite);
                h1_trueNoPeak_ee->SetBinError(1, 0);
                h1_trueNoPeak_ee->SetBinError(2, 0);
                h1_trueNoPeak_ee->SetBinError(3, 0);
        TH1F *h1_trueNoPeak_ee_errDown = new TH1F("h1_trueNoPeak_ee_errDown", "h1_trueNoPeak_ee_errDown", 3, -0.5, 2.5);
                h1_trueNoPeak_ee_errDown->SetLineColor(10);
                h1_trueNoPeak_ee_errDown->SetFillColor(10);
                h1_trueNoPeak_ee_errDown->SetFillStyle(1001);
        TH1F *h1_trueNoPeak_ee_errUp = new TH1F("h1_trueNoPeak_ee_errUp", "h1_trueNoPeak_ee_errUp", 3, -0.5, 2.5);
                h1_trueNoPeak_ee_errUp->SetLineColor(10);
                h1_trueNoPeak_ee_errUp->SetFillColor(kRed);
                h1_trueNoPeak_ee_errUp->GetXaxis()->SetTitle("n-Jets");
                h1_trueNoPeak_ee_errUp->GetYaxis()->SetTitle("Number of Events");

        TLegend *lg_estimateNoPeak = new TLegend(0.2, 0.6, 0.6, 0.9);
        lg_estimateNoPeak->SetFillColor(kWhite);
        lg_estimateNoPeak->SetLineColor(kWhite);
        lg_estimateNoPeak->AddEntry(h1_trueNoPeak_ee, "True Non-peaking BG", "l");
        lg_estimateNoPeak->AddEntry(h1_trueNoPeak_ee_errUp, "True Non-peaking BG #pm 1#sigma", "f");
        lg_estimateNoPeak->AddEntry(h1_estimateNoPeak_ee, "Estimated Non-peaking BG", "lp");

        TCanvas *c_ee = new TCanvas();
        TCanvas *c_mm = new TCanvas();
        TCanvas *cNoPeak_ee = new TCanvas();
        TCanvas *cNoPeak_mm = new TCanvas();

        //
        // ee
        //
        estimateDY(f, "0j", "ee", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_ee->SetBinContent(1, estimate);
        h1_estimate_ee->SetBinError(1, error);
        h1_true_ee->SetBinContent(1, truth);
        h1_true_ee_errDown->SetBinContent(1, truth - truth_error);
        h1_true_ee_errUp->SetBinContent(1, truth + truth_error);
        h1_estimateNoPeak_ee->SetBinContent(1, noPeak_estimate);
        h1_estimateNoPeak_ee->SetBinError(1, noPeak_error);
        h1_trueNoPeak_ee->SetBinContent(1, noPeak_truth);
        h1_trueNoPeak_ee_errDown->SetBinContent(1, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_ee_errUp->SetBinContent(1, noPeak_truth + noPeak_truth_error);


        estimateDY(f, "1j", "ee", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_ee->SetBinContent(2, estimate);
        h1_estimate_ee->SetBinError(2, error);
        h1_true_ee->SetBinContent(2, truth);
        h1_true_ee_errDown->SetBinContent(2, truth - truth_error);
        h1_true_ee_errUp->SetBinContent(2, truth + truth_error);
        h1_estimateNoPeak_ee->SetBinContent(2, noPeak_estimate);
        h1_estimateNoPeak_ee->SetBinError(2, noPeak_error);
        h1_trueNoPeak_ee->SetBinContent(2, noPeak_truth);
        h1_trueNoPeak_ee_errDown->SetBinContent(2, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_ee_errUp->SetBinContent(2, noPeak_truth + noPeak_truth_error);

        estimateDY(f, "2j", "ee", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_ee->SetBinContent(3, estimate);
        h1_estimate_ee->SetBinError(3, error);
        h1_true_ee->SetBinContent(3, truth);
        h1_true_ee_errDown->SetBinContent(3, truth - truth_error);
        h1_true_ee_errUp->SetBinContent(3, truth + truth_error);
        h1_estimateNoPeak_ee->SetBinContent(3, noPeak_estimate);
        h1_estimateNoPeak_ee->SetBinError(3, noPeak_error);
        h1_trueNoPeak_ee->SetBinContent(3, noPeak_truth);
        h1_trueNoPeak_ee_errDown->SetBinContent(3, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_ee_errUp->SetBinContent(3, noPeak_truth + noPeak_truth_error);

        c_ee->cd();
        h1_true_ee_errUp->Draw("");
        h1_true_ee_errDown->Draw("SAME HIST");
        h1_estimate_ee->Draw("SAME E1");
        h1_true_ee->Draw("SAME E1");
        lg_estimate->Draw();
        h1_true_ee_errUp->GetYaxis()->SetRangeUser(0, 50);
        c_ee->RedrawAxis();

        cNoPeak_ee->cd();
        h1_trueNoPeak_ee_errUp->Draw("");
        h1_trueNoPeak_ee_errDown->Draw("SAME HIST");
        h1_estimateNoPeak_ee->Draw("SAME E1");
        h1_trueNoPeak_ee->Draw("SAME E1");
        lg_estimateNoPeak->Draw();
        h1_trueNoPeak_ee_errUp->GetYaxis()->SetRangeUser(0, 50);
        cNoPeak_ee->RedrawAxis();


        //
        // mm
        //
        estimateDY(f, "0j", "mm", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_mm->SetBinContent(1, estimate);
        h1_estimate_mm->SetBinError(1, error);
        h1_true_mm->SetBinContent(1, truth);
        h1_true_mm_errDown->SetBinContent(1, truth - truth_error);
        h1_true_mm_errUp->SetBinContent(1, truth + truth_error);
        h1_estimateNoPeak_mm->SetBinContent(1, noPeak_estimate);
        h1_estimateNoPeak_mm->SetBinError(1, noPeak_error);
        h1_trueNoPeak_mm->SetBinContent(1, noPeak_truth);
        h1_trueNoPeak_mm_errDown->SetBinContent(1, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_mm_errUp->SetBinContent(1, noPeak_truth + noPeak_truth_error);

        estimateDY(f, "1j", "mm", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_mm->SetBinContent(2, estimate);
        h1_estimate_mm->SetBinError(2, error);
        h1_true_mm->SetBinContent(2, truth);
        h1_true_mm_errDown->SetBinContent(2, truth - truth_error);
        h1_true_mm_errUp->SetBinContent(2, truth + truth_error);
        h1_estimateNoPeak_mm->SetBinContent(2, noPeak_estimate);
        h1_estimateNoPeak_mm->SetBinError(2, noPeak_error);
        h1_trueNoPeak_mm->SetBinContent(2, noPeak_truth);
        h1_trueNoPeak_mm_errDown->SetBinContent(2, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_mm_errUp->SetBinContent(2, noPeak_truth + noPeak_truth_error);

        estimateDY(f, "2j", "mm", truth, estimate, error, truth_error,
                noPeak_estimate, noPeak_error, noPeak_truth, noPeak_truth_error);
        h1_estimate_mm->SetBinContent(3, estimate);
        h1_estimate_mm->SetBinError(3, error);
        h1_true_mm->SetBinContent(3, truth);
        h1_true_mm_errDown->SetBinContent(3, truth - truth_error);
        h1_true_mm_errUp->SetBinContent(3, truth + truth_error);
        h1_estimateNoPeak_mm->SetBinContent(3, noPeak_estimate);
        h1_estimateNoPeak_mm->SetBinError(3, noPeak_error);
        h1_trueNoPeak_mm->SetBinContent(3, noPeak_truth);
        h1_trueNoPeak_mm_errDown->SetBinContent(3, noPeak_truth - noPeak_truth_error);
        h1_trueNoPeak_mm_errUp->SetBinContent(3, noPeak_truth + noPeak_truth_error);

        c_mm->cd();
        h1_true_mm_errUp->Draw("");
        h1_true_mm_errDown->Draw("SAME HIST");
        h1_estimate_mm->Draw("SAME E1");
        h1_true_mm->Draw("SAME E1");
        lg_estimate->Draw();
        h1_true_mm_errUp->GetYaxis()->SetRangeUser(0, 50);
        c_mm->RedrawAxis();

        cNoPeak_mm->cd();
        h1_trueNoPeak_mm_errUp->Draw("");
        h1_trueNoPeak_mm_errDown->Draw("SAME HIST");
        h1_estimateNoPeak_mm->Draw("SAME E1");
        h1_trueNoPeak_mm->Draw("SAME E1");
        lg_estimateNoPeak->Draw();
        h1_trueNoPeak_mm_errUp->GetYaxis()->SetRangeUser(0, 50);
        cNoPeak_mm->RedrawAxis();


}


