
#include "HistogramUtilities.h"

HistogramUtilities::HistogramUtilities(TString fileName, Double_t metCut, 
				bool doWZ, Double_t lumiNorm) 
	: doWZ_(doWZ), lumiNorm_(lumiNorm)
{

	// define the data source names
	// input root files will contain
	// a histogram for each source...
	sources_.push_back(	fH_WW()		);
	sources_.push_back(	fH_TTBAR()	);
	sources_.push_back(	fH_DYMM()	);
	sources_.push_back(	fH_DYEE()	);
	sources_.push_back(	fH_DYTT()	);
	sources_.push_back(	fH_WJETS()	);
	sources_.push_back(	fH_TW()		);
	sources_.push_back(	fH_ZZ()		);
	sources_.push_back(	fH_WZ()		);

	// open root file
	file_ = new TFile(fileName, "READ");

	// set the bin position for the Z-veto
        TH1F *h_tempMll = (TH1F*)file_->Get("ww_mll_0j-N-1_ee");
        zLow_ = h_tempMll->FindBin(76.1);
        zHigh_ = h_tempMll->FindBin(105.9);
	std::cout << zLow_ << ", " << zHigh_ << std::endl;
	// set the met cut bin
        TH1F *h_tempMet = (TH1F*)file_->Get("ww_met_in_0j-N-1_ee");
        metCutBin_ = h_tempMet->FindBin(metCut + 0.01);

	delete h_tempMll;
	delete h_tempMet;

}

void HistogramUtilities::printInOutTruth(TString nJets, TString hyp_type) 
{

        Float_t totalIn = 0;
        Float_t totalOut = 0;
        Float_t totalIn_nopeak = 0;
        Float_t totalOut_nopeak = 0;

        //
        // TABLE comparing no-peak bg with truth
        std::cout       << "TABLE --- TRUTH " << nJets.Data() << std::endl;
	std::cout << " \\begin{table}[ht]" << std::endl;
	std::cout << " \\caption{in and out truth for " << hyp_type.Data() << "}" << std::endl;
	std::cout << " \\begin{center}" << std::endl;
	std::cout << " \\begin{tabular}{|c|c|c|}" << std::endl;
	std::cout << " \\hline" << std::endl;
        std::cout << "\t & \t in \t & \t out \\\\ \\hline \\hline" << std::endl;

        TString histNameSuffix = "_mll_" + nJets + "-N-1_" + hyp_type;
        for (unsigned int i = 0; i < sources_.size(); ++i)
        {

		// get the histogram for this source
                TH1F *h_temp = ((TH1F*)file_->Get(sources_[i].getName() + histNameSuffix));
		// get weight for this histogram
		Float_t weight = h_temp->Integral(0, h_temp->GetNbinsX() + 1)/h_temp->GetEntries();

                // compute nIn and nOut
                Float_t nIn = h_temp->Integral(zLow_, zHigh_);
                Float_t nOut = h_temp->Integral(0, h_temp->GetNbinsX() + 1) - nIn;
                totalIn += nIn;
                totalOut += nOut;

                // keep track of the "non-peaking" contribution
                if ((sources_nonpeaking & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
                	totalIn_nopeak += nIn;
                        totalOut_nopeak += nOut;
                }

		// print this row
                std::cout << sources_[i].getName().Data() << "\t & \t ";
                printf("%.2f (%.2f) \t & %.2f (%.2f) \t \\\\ \\hline \n", nIn*lumiNorm_, nIn/weight, nOut*lumiNorm_, nOut/weight);

	}

	// print totals
        std::cout << "TOTAL \t & \t";
        printf("%.2f \t & %.2f \t \\\\ \\hline \n", totalIn*lumiNorm_, totalOut*lumiNorm_);
        std::cout << "NOPEAK \t & \t";
        printf("%.2f \t & %.2f \t \\\\ \\hline \n", totalIn_nopeak*lumiNorm_, totalOut_nopeak*lumiNorm_);
	std::cout << "\\end{tabular}" << std::endl;
	std::cout << "\\end{center}" << std::endl;
	std::cout << "\\label{default}" << std::endl;
	std::cout << "\\end{table}" << std::endl;

}

Float_t HistogramUtilities::getN(bool in, const sources_t &theSources, 
			TString nJets, TString hyp_type, Float_t &err2)
{

	// get the relevant histogram
        TH1F *h_temp = getHistogram(theSources, "mll", nJets, hyp_type);

        // compute n ( in our out )
	Float_t n = 0;

	if (in) {
		// get integral and error "in"
       		n = h_temp->Integral(zLow_, zHigh_);
		err2 = Utilities::getError2(h_temp, zLow_, zHigh_);
	}
	else {
		// get integral and error "out"
		n = h_temp->Integral(0, h_temp->GetNbinsX() + 1) 
			- h_temp->Integral(zLow_, zHigh_);
		err2 = Utilities::getError2(h_temp, 0, zLow_)
			+ Utilities::getError2(h_temp, zHigh_, h_temp->GetNbinsX() + 1);
	}

	// tidy up and return
        delete h_temp;
	return n;
}

TH1F* HistogramUtilities::getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type) 
{
        TString histNameSuffix = "_" + var + "_" + nJets + "-N-1_" + hyp_type;
        TH1F *h1_data = (TH1F*)(file_->Get(sources_[0].getName() + histNameSuffix)->Clone());
	h1_data->Reset();
        for (unsigned int i = 0; i < sources_.size(); ++i)
        {
		if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) )
                	h1_data->Add((TH1F*)file_->Get(sources_[i].getName() + histNameSuffix));
        }       
        h1_data->Scale(lumiNorm_);
        return h1_data;
}

THStack* HistogramUtilities::getStack(sources_t theSources, TString var, TString nJets, TString hyp_type) 
{

	// create a new stack object
	THStack *st_temp = new THStack("st_temp", "");

	// get each constituent in turn and add to the stack
        TString histNameSuffix = "_" + var + "_" + nJets + "-N-1_" + hyp_type;
        for (int i = sources_.size() - 1; i >= 0; --i)
        {
                if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
	                TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone()));
			h1_temp->Rebin(5);
        	        h1_temp->Scale(lumiNorm_);
                	st_temp->Add(h1_temp);
		}
         }
         return st_temp;
}

TLegend* HistogramUtilities::getLegend(sources_t theSources) 
{

	// create a new legend object and make it look nice
	TLegend *lg = new TLegend(0.7, 0.5, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
	lg->SetShadowColor(kWhite);

	// get each constituent in turn and add to the legend
        TString histNameSuffix = "_mll_0j-N-1_mm";
        for (unsigned int i = 0; i < sources_.size(); ++i)
        {
                if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
		        TH1F *h_temp = (TH1F*)file_->Get(sources_[i].getName() + histNameSuffix)->Clone();
        		lg->AddEntry(h_temp, sources_[i].getName(), "f");
		}
        }
        return lg;
}

TH1F* HistogramUtilities::getRHist(DataSource theSource, TString nJets, TString hyp_type)
{

        TH1F *h1_met_in = (TH1F*)(file_->Get(theSource.getName() + "_met_in_" + nJets + "-N-1_" + hyp_type));
        TH1F *h1_met_out = (TH1F*)(file_->Get(theSource.getName() + "_met_out_" + nJets + "-N-1_" + hyp_type));
        h1_met_in->Rebin(5);
        h1_met_out->Rebin(5);

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


