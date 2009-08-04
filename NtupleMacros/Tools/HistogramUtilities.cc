
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

  sources_.push_back(     fH_WJET_ALP()   );	
  sources_.push_back(	fH_ZEEJET_ALP()	);
  sources_.push_back(     fH_ZMMJET_ALP() );
  sources_.push_back(     fH_ZTTJET_ALP() );

  sources_.push_back(	fH_QCD30()	);
  sources_.push_back(     fH_QCD80()      );

  // open root file
  file_ = new TFile(fileName, "READ");

  // leave the luminosity norm as in the root file
  lumiNorm_ = lumiNorm;

}

void HistogramUtilities::setOrder(std::vector<DataSource> potentialSources)
{
	sources_ = potentialSources;
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

TH2F* HistogramUtilities::get2dHistogram(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin) 
{
  //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
  TString histNameSuffix = "_" + var + nJets + "_" + hyp_type;

  TH2F *h_data = 0;
  for (unsigned int i = 0; i < sources_.size(); ++i)
	{
	  if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
		//std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
		if (!h_data) h_data = (TH2F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone());
		else h_data->Add((TH2F*)file_->Get(sources_[i].getName() + histNameSuffix));
	  }
	}       
  //h_data->Scale(lumiNorm_);
  //h_data->Rebin(rebin);
  return h_data;
}

//get a stack of single var, hyp for given sources
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
			if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
			h1_temp->Rebin(rebin);
        	        h1_temp->Scale(lumiNorm_);
                	st_temp->Add(h1_temp);
		}
         }
         return st_temp;
}

//for combining two hyps (not two vars, if you'd ever want to do that anyway)
THStack* HistogramUtilities::getSumStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, Int_t rebin) {

  // create a new stack object
  TString name = var + nJets + "_" + hyp1 + hyp2; //name includes both
  THStack *st_temp = new THStack(name, name);

  // get each constituent in turn and add to the stack
  //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
  TString histNameSuffix1 = "_" + var + nJets + "_" + hyp1;
  TString histNameSuffix2 = "_" + var + nJets + "_" + hyp2;
  
  for (int i = sources_.size() - 1; i >= 0; --i) {
	if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
	  //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
	  TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix1)->Clone()));
	  TH1F *h2_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix2)->Clone()));
	  if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
	  if (sources_[i].getColor() != 0) h2_temp->SetFillColor(sources_[i].getColor());
	  h1_temp->Rebin(rebin);
	  h1_temp->Scale(lumiNorm_);
	  h2_temp->Rebin(rebin);
	  h2_temp->Scale(lumiNorm_);

	  h1_temp->Add(h2_temp); //now h1 is h1+h2
	  st_temp->Add(h1_temp); //put the sum in the stack
	}
  }
  return st_temp;
}


//for combining three hyps, sum of the first two minus the third
THStack* HistogramUtilities::getSumDifStack(sources_t theSources, TString var, TString nJets, TString hyp1, TString hyp2, TString hyp3, Int_t rebin) {

  // create a new stack object
  TString name = var + nJets + "_" + hyp1 + hyp2 + "-" + hyp3; //name includes all
  THStack *st_temp = new THStack(name, name);

  // get each constituent in turn and add to the stack
  //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
  TString histNameSuffix1 = "_" + var + nJets + "_" + hyp1;
  TString histNameSuffix2 = "_" + var + nJets + "_" + hyp2;
  TString histNameSuffix3 = "_" + var + nJets + "_" + hyp3;
  
  for (int i = sources_.size() - 1; i >= 0; --i) {
	if( (theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
	  //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
	  TH1F *h1_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix1)->Clone()));
	  TH1F *h2_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix2)->Clone()));
	  TH1F *h3_temp = ((TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix3)->Clone()));
	  if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
	  if (sources_[i].getColor() != 0) h2_temp->SetFillColor(sources_[i].getColor());
	  if (sources_[i].getColor() != 0) h3_temp->SetFillColor(sources_[i].getColor());
	  h1_temp->Rebin(rebin);
	  h1_temp->Scale(lumiNorm_);
	  h2_temp->Rebin(rebin);
	  h2_temp->Scale(lumiNorm_);
	  h3_temp->Rebin(rebin);
	  h3_temp->Scale(lumiNorm_);

	  h1_temp->Add(h2_temp); //now h1 is h1+h2
	  h1_temp->Add(h3_temp, -1.0); //this adds -1*h3, ie, subtracts h3
	  st_temp->Add(h1_temp); //put the sum in the stack
	}
  }
  return st_temp;
}

TH1F* HistogramUtilities::getHistogramSum(sources_t theSources, TString var, TString nJets, TString hyp_type, Int_t rebin) 
{
  //create a new histogram object
  //TString histNameSuffix = "_" + var + "_" + nJets + hyp_type;
  TString histNameSuffix = "_" + var + nJets + "_" + hyp_type;

  TH1F *h_temp = (TH1F*)(file_->Get(sources_[0].getName() + histNameSuffix)->Clone());
  TH1F *st_temp = new TH1F(var, var, h_temp->GetNbinsX(), h_temp->GetXaxis()->GetXmin(), h_temp->GetXaxis()->GetXmax());

  // get each constituent in turn and add to the stack
  for (int i = sources_.size() - 1; i >= 0; --i) {
	if ((theSources & makeBit(sources_[i].getSource()) ) == makeBit(sources_[i].getSource()) ) {
	  //std::cout << "getting " << sources_[i].getName() + histNameSuffix << std::endl;
	  TH1F *h1_temp = (TH1F*)(file_->Get(sources_[i].getName() + histNameSuffix)->Clone());
	  if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
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
		        TH1F *h1_temp = (TH1F*)file_->Get(sources_[i].getName() + histNameSuffix)->Clone();
                        if (sources_[i].getColor() != 0) h1_temp->SetFillColor(sources_[i].getColor());
        		lg->AddEntry(h1_temp, sources_[i].getName(), "f");
		}
        }
        return lg;
}

