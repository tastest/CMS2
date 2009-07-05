#include <stdio.h>
#include <vector>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "signedImpact.h"
#include "CORE/CMS2.h"
#include <math.h>

using std::vector;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

bool passTrkjetCuts (int i_trk)
{

  return  fabs(cms2.trks_d0corr()[i_trk]) < 0.1 &&
       cms2.trks_chi2()[i_trk] / cms2.trks_ndof()[i_trk] < 5 &&
       cms2.trks_validHits()[i_trk]   > 10 &&
       cms2.trks_trk_p4()[i_trk].Pt() > 1.0 &&
    cms2.trks_trk_p4()[i_trk].Pt() < 200.0;
}

void calculateJetProb (const std::vector<std::pair<LorentzVector, std::vector<unsigned int> > > &trackJets,
		       vector<double> *jetprobs)
{
     jetprobs->clear();
     jetprobs->reserve(trackJets.size());
     for (unsigned int itrkjet = 0; itrkjet < trackJets.size(); ++itrkjet) {
	  double Pi = 1;
	  int n_jetprob_trks = 0;
	  for (unsigned int i = 0; i < trackJets[itrkjet].second.size(); ++i) {
	       const std::pair<double, double> sip = 
		    signedImpact(trackJets[itrkjet].first, trackJets[itrkjet].second[i]);
	       const double d0 = sip.first;
	       const double d0err = sip.second;
	       const double sipsig = d0 / d0err / 1.1;
	       if (sipsig > 0) {
		 Pi *= (1 - erf( sipsig / TMath::Sqrt(2) ));
		 n_jetprob_trks++;
	       }
	  }
	  double sum = 0;
	  int fact = 1;
	  for (unsigned int i = 0; i < trackJets[itrkjet].second.size(); fact *= ++i) {
	    if (i > 1)
	      //	      std::cout << "fact = " << fact << ", " << i << std::endl;
		 sum += ( ::pow(-log(Pi), i) / fact );
	  }
	  double jetprob = Pi * sum;
	  jetprobs->push_back(jetprob);	  
     }
}

