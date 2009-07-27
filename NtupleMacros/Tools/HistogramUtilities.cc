
#include "Tools/HistogramUtilities.h"

HistogramUtilities::HistogramUtilities(TString fileName, Double_t lumiNorm) 
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

        sources_.push_back(     fH_WENU()         );
        sources_.push_back(     fH_EM30_80()         );
        sources_.push_back(     fH_BC30_80()         );

	// open root file
	file_ = new TFile(fileName, "READ");

	// leave the luminosity norm as in the root file
	lumiNorm_ = lumiNorm;

}

TH1F* HistogramUtilities::getHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin) 
{
        TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
        TH1F *h1_data = 0;
        for (unsigned int i = 0; i < sources_.size(); ++i)
        {
		if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
	                std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
			if (!h1_data) h1_data = (TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone());
                	else h1_data->Add((TH1F*)file_->Get(sources_[i].getName() + histNameSuffix));
		}
        }       
        h1_data->Scale(lumiNorm_);
	h1_data->Rebin(rebin);
        return h1_data;
}

THStack* HistogramUtilities::getStack(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin) 
{

	// create a new stack object
	THStack *st_temp = new THStack("st_temp", "");

	// get each constituent in turn and add to the stack
        TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
        for (int i = sources_.size() - 1; i >= 0; --i)
        {
                if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
			std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
	                TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone()));
			h1_temp->Rebin(rebin);
        	        h1_temp->Scale(lumiNorm_);
                	st_temp->Add(h1_temp);
		}
         }
         return st_temp;
}

TLegend* HistogramUtilities::getLegend(sources_t theSources, TString var, TString nJets, TString hyp_type) 
{

	// create a new legend object and make it look nice
	TLegend *lg = new TLegend(0.7, 0.5, 0.9, 0.9);
        lg->SetFillColor(kWhite);
        lg->SetLineColor(kWhite);
	lg->SetShadowColor(kWhite);

	// get each constituent in turn and add to the legend
        TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
        for (unsigned int i = 0; i < sources_.size(); ++i)
        {
                if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
		        TH1F *h_temp = (TH1F*)file_->Get(sources_[i].getName() + histNameSuffix)->Clone();
        		lg->AddEntry(h_temp, sources_[i].getName(), "f");
		}
        }
        return lg;
}

