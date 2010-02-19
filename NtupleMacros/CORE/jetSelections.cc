// $Id: jetSelections.cc,v 1.3 2010/02/19 21:01:14 jmuelmen Exp $

#include <algorithm>
#include <utility>
#include "Math/VectorUtil.h"
#include "jetSelections.h"

using std::vector;
using std::pair;

// define this type for speed: allows us to get a vector of selected
// jets, potentially with a correction factor, without having to make
// copies
typedef vector<pair<const LorentzVector *, double> > jets_with_corr_t;

// function to give us the indices of jets passing kinematic and cleaning cuts
static jets_with_corr_t getJets_fast (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
				      double deltaR, double min_pt, double max_eta)
{
     // JPT, PF or calo jets?  Introduce this variable so we only have to decide once
     const vector<LorentzVector> *jets = 0;
     switch (type) {
     case JETS_TYPE_JPT:
	  jets = &cms2.jpts_p4();
	  break;
     case JETS_TYPE_CALO_CORR: case JETS_TYPE_CALO_UNCORR:
	  jets = &cms2.jets_p4();
	  break;
     case JETS_TYPE_PF_CORR: case JETS_TYPE_PF_UNCORR:
	  jets = &cms2.pfjets_p4();
	  break;
     }
     jets_with_corr_t ret;
     ret.reserve(jets->size()); // reserve so we don't have to realloc later, which is slow
     for (unsigned int i = 0; i < jets->size(); ++i) {
	  //------------------------------------------------------------
	  // min pt cut
	  //------------------------------------------------------------
	  double corr = 1;
	  // CALO_CORR and PF_CORR need to be pt-corrected
	  switch (type) {
	  case JETS_TYPE_CALO_CORR:
	       corr = cms2.jets_cor().at(i);
	       break;
	  case JETS_TYPE_PF_CORR:
	       corr = cms2.pfjets_cor().at(i);
	       break;
	  case JETS_TYPE_JPT: 
	  case JETS_TYPE_CALO_UNCORR: 
	  case JETS_TYPE_PF_UNCORR:
	       break;
	  }
	  const double pt = jets->at(i).pt() * corr;
	  if (pt < min_pt) 
	       continue;
	  //------------------------------------------------------------
	  // max |eta| cut
	  //------------------------------------------------------------
	  if (fabs(jets->at(i).eta()) > max_eta)
	       continue;
	  //------------------------------------------------------------
	  // lepton cleaning
	  //------------------------------------------------------------
	  switch (cleaning) {
	  case JETS_CLEAN_NONE:
	       break;
	  case JETS_CLEAN_HYP_E_MU:
	  {
	       const LorentzVector &lt = cms2.hyp_lt_p4().at(i_hyp);
	       if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), lt) < deltaR)
		    continue;
	       const LorentzVector &ll = cms2.hyp_ll_p4().at(i_hyp);
	       if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), ll) < deltaR)
		    continue;
	       break;
	  }
	  case JETS_CLEAN_HYP_E:
	  {
	       const LorentzVector &lt = cms2.hyp_lt_p4().at(i_hyp);
	       const int lt_id = cms2.hyp_lt_id().at(i_hyp);
	       if (abs(lt_id) == 11 && ROOT::Math::VectorUtil::DeltaR(jets->at(i), lt) < deltaR)
		    continue;
	       const LorentzVector &ll = cms2.hyp_ll_p4().at(i_hyp);
	       const int ll_id = cms2.hyp_ll_id().at(i_hyp);
	       if (abs(ll_id) == 11 && ROOT::Math::VectorUtil::DeltaR(jets->at(i), ll) < deltaR)
		    continue;
	       break;
	  }
	  case JETS_CLEAN_SINGLE_E:
	  {
	       const LorentzVector &e = cms2.els_p4().at(i_hyp);
	       if (ROOT::Math::VectorUtil::DeltaR(jets->at(i), e) < deltaR)
		    continue;
	       break;
	  }
	  default:
	       assert(false);
	  }
	  //------------------------------------------------------------
	  // lepton passed all cuts
	  //------------------------------------------------------------
	  ret.push_back(pair<const LorentzVector *, double>(&jets->at(i), corr));
     }
     return ret;
}

struct jets_pt_gt {
     bool operator () (const LorentzVector &v1, const LorentzVector &v2) 
	  {
	       return v1.pt() > v2.pt();
	  }
};

// functions that we let other people use
vector<LorentzVector> getJets (unsigned int i_hyp, bool sort_, 
			       enum JetType type, enum CleaningType cleaning,
			       double deltaR, double min_pt, double max_eta)
{
     jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta);
     vector<LorentzVector> ret;
     ret.reserve(jets.size());
     for (unsigned int i = 0; i < jets.size(); ++i) {
	  // correct the jet momentum if a corrected jet type was requested
	  ret.push_back(*jets[i].first * jets[i].second);
     }
     if (sort_)
	  sort(ret.begin(), ret.end(), jets_pt_gt());
     return ret;
}

int nJets (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
	   double deltaR, double min_pt, double max_eta)
{
     jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta);
     return jets.size();
}

double sumPt (unsigned int i_hyp, enum JetType type, enum CleaningType cleaning,
	      double deltaR, double min_pt, double max_eta)
{
     jets_with_corr_t jets = getJets_fast(i_hyp, type, cleaning, deltaR, min_pt, max_eta);
     double ret = 0;
     for (unsigned int i = 0; i < jets.size(); ++i) {
	  // correct the jet momentum if a corrected jet type was requested
	  ret += jets[i].first->pt() * jets[i].second;
     }
     return ret;
}
