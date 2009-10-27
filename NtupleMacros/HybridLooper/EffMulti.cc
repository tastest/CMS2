
#include "HybridLooper/EffMulti.h"

#include <iostream>

EffMulti::EffMulti(bool lessThan, float thresholdEB, float thresholdEE, 
	std::string source, std::string var, std::string det) 
	: lessThan_(lessThan), thresholdEB_(thresholdEB), thresholdEE_(thresholdEE)
{
	e1_pt_ = new EffH1F(source + "_" + var + "_pt_" + det, var + ";p_{T} (GeV);", 100, 0, 200);
        e1_eta_ = new EffH1F(source + "_" + var + "_eta_" + det, var + ";#eta;", 24, -3, 3);
        e1_phi_ = new EffH1F(source + "_" + var + "_phi_" + det, var + ";#phi (radians);", 72, -3.14159, 3.14159);
}


EffMulti::~EffMulti()
{
	delete e1_pt_;
	delete e1_eta_;
	delete e1_phi_;
}

void EffMulti::Fill(float value, float pt, float eta, float phi, float weight)
{
	e1_pt_->denom->Fill(pt, weight);
	if (pt > 20.0) {
		e1_eta_->denom->Fill(eta, weight);
		e1_phi_->denom->Fill(phi, weight);
	}

	float threshold = thresholdEB_;
	if (fabs(eta) > 1.5) threshold = thresholdEE_;

	if (lessThan_ && (value < threshold)) {
	        e1_pt_->numer->Fill(pt, weight);
        	if (pt > 20.0) e1_eta_->numer->Fill(eta, weight);
	        e1_phi_->numer->Fill(phi, weight);
	}
	else if (!lessThan_ && (value > threshold)) {
                e1_pt_->numer->Fill(pt, weight);
                if (pt > 20.0) e1_eta_->numer->Fill(eta, weight);
                e1_phi_->numer->Fill(phi, weight);
	}
}


