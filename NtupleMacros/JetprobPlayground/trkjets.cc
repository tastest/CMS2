#include <stdio.h>
#include <vector>
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
//#include "KtJet/KtEvent.h"
//#include "KtJet/KtLorentzVector.h"
#include "signedImpact.h"
#include "CORE/CMS2.h"
#include <math.h>

//using namespace KtJet;
using std::vector;

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

bool passTrkjetCuts (int i_trk)
{

  return  fabs(cms2.trks_d0corr()[i_trk]) < 0.1 &&
       cms2.trks_chi2()[i_trk] / cms2.trks_ndof()[i_trk] < 5 &&
       cms2.trks_validHits()[i_trk]   > 10 &&
       cms2.trks_trk_p4()[i_trk].Pt() > 1.0 &&
    cms2.trks_trk_p4()[i_trk].Pt() < 200.0;
    //    cms2.trks_d0Err()[i_trk] > 20e-4;
}

/*pair<vector<KtLorentzVector>, vector<unsigned int> > trkjetTracks (int i_hyp)
{
     pair<vector<KtLorentzVector>, vector<unsigned int> > tracks;
     for (unsigned int i = 0; i < cms2.trks_trk_p4().size(); ++i) {
	  // track quality cuts
	  if (!passTrkjetCuts(i))
	       continue;
	  // skip tracks close to lt or ll
	  if (ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp], 
					     cms2.trks_trk_p4()[i]) < 0.4)
	       continue;
	  if (ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp], 
					     cms2.trks_trk_p4()[i]) < 0.4)
	       continue;
	  // skip tracks with large delta z0
	  if (fabs(cms2.trks_z0()[i] - cms2.hyp_lt_z0()[i_hyp]) > 1)
	       continue;
	  // otherwise use for clustering
	  const LorentzVector &p4 = cms2.trks_trk_p4()[i];
	  KtLorentzVector v(p4.x(), p4.y(), p4.z(), p4.E());
	  tracks.first.push_back(v);
	  tracks.second.push_back(i);
     }
     return tracks;
}

vector<KtLorentzVector> findTrkJetsByEt (const pair<vector<KtLorentzVector>, 
					 vector<unsigned int> > &tracks)
{
//      for (unsigned int i = 0; i < cms2.jets_p4().size(); ++i) {
// 	  const LorentzVector &p4 = cms2.jets_p4()[i];
//  	  if (cms2.jets_emFrac()[i] > emf) {
// 	       KtLorentzVector v(p4.x(), p4.y(), p4.z(), p4.E());
// 	       tracks.push_back(scale * v);
// 	  }
//      }
     KtEvent ev(tracks.first, 4, 2, 4, 1);
     return ev.getJetsEt();
}*/

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
	       const double sipsig = d0 / d0err;
	       if (sipsig > 0) {
		 Pi *= exp(-0.5 * sipsig * sipsig / 1.1 / 1.1);
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

/*vector<KtLorentzVector> findTrkJetsByEt (int i_hyp)
{
     return findTrkJetsByEt(trkjetTracks(i_hyp));
}

std::pair<double, double> maxTrkJetEt (int i_hyp)
{
     vector<KtLorentzVector> jets = findTrkJetsByEt(i_hyp);
     printf("b quarks: ");
     for (unsigned int i = 0; i < cms2.genps_id().size(); ++i) {
	  if (abs(cms2.genps_id()[i]) != 5)
	       continue;
 	  printf("%5.1f %4.1f %4.1f\t", 
		 cms2.genps_p4()[i].Et(), cms2.genps_p4()[i].eta(), 
		 cms2.genps_p4()[i].phi());
     }
     printf("track jets: ");
     for (unsigned int i = 0; i < jets.size(); ++i) {
 	  printf("%5.1f %4.1f %4.1f (%2d)\t", jets[i].et(), jets[i].eta(), jets[i].phi(), 
		 jets[i].getNConstituents());
     }
     printf("other jets: ");
     for (unsigned int i = 0, i_jet = 0; i < cms2.hyp_other_jets_p4()[i_hyp].size(); ++i, ++i_jet) {
	  while (fabs(cms2.hyp_other_jets_p4()[i_hyp][i].Et() - cms2.jets_p4()[i_jet].Et()) < 
		 1e-5 * fabs(cms2.jets_p4()[i_jet].Et()))
	       i_jet++;
  	  printf("%5.1f %4.1f %4.1f (%4.2f)\t", cms2.hyp_other_jets_p4()[i_hyp][i].Et() * cms2.jets_tq_noCorrF()[i_jet], 
		 cms2.hyp_other_jets_p4()[i_hyp][i].eta(), cms2.hyp_other_jets_p4()[i_hyp][i].phi(),
		 cms2.hyp_other_jets_emFrac()[i_hyp][i]);
     }
//      printf("yet other jets: ");
//      for (unsigned int i = 0; i < cms2.jets_p4().size(); ++i) {
//  	  printf("%5.1f %4.1f %4.1f\t", cms2.jets_p4()[i].Et(), cms2.jets_p4()[i].eta(), cms2.jets_p4()[i].phi());
//      }
     printf("\n");
//      for (unsigned int i = 0; i < jets.size(); ++i) {
// 	  const LorentzVector j(jets[i].x(), jets[i].y(), jets[i].z(), jets[i].e());
// 	  if (ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp], j) < 0.4)
// 	       continue;
// 	  if (ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp], j) < 0.4)
// 	       continue;
// 	  return jets[i].et();
//      }
//      return 0;
#if 0
     if (jets.size() == 0)
	  return 0;
     LorentzVector trkjet(jets[0].x(), jets[0].y(), jets[0].z(), jets[0].e());
     LorentzVector trk_plus_calo_jet = trkjet;
     for (unsigned int i = 0, i_jet = 0; i < cms2.hyp_other_jets_p4()[i_hyp].size(); ++i, ++i_jet) {
	  while (fabs(cms2.hyp_other_jets_p4()[i_hyp][i].Et() - cms2.jets_p4()[i_jet].Et()) < 
		 1e-5 * fabs(cms2.jets_p4()[i_jet].Et()))
	       i_jet++;
	  LorentzVector calojet = cms2.hyp_other_jets_p4()[i_hyp][i] * 
	       cms2.jets_tq_noCorrF()[i_jet]; // use uncorrected Et
	  if (ROOT::Math::VectorUtil::DeltaR(calojet, trkjet) < 0.25) 
	       trk_plus_calo_jet += cms2.hyp_other_jets_emFrac()[i_hyp][i] * calojet;
     }
     for (unsigned int i = 0, i_jet = 0; i < cms2.hyp_jets_p4()[i_hyp].size(); ++i, ++i_jet) {
	  while (fabs(cms2.hyp_jets_p4()[i_hyp][i].Et() - cms2.jets_p4()[i_jet].Et()) < 
		 1e-5 * fabs(cms2.jets_p4()[i_jet].Et()))
	       i_jet++;
	  LorentzVector calojet = cms2.hyp_jets_p4()[i_hyp][i] * 
	       cms2.jets_tq_noCorrF()[i_jet]; // use uncorrected Et
	  if (ROOT::Math::VectorUtil::DeltaR(calojet, trkjet) < 0.25) 
	       trk_plus_calo_jet += cms2.hyp_jets_emFrac()[i_hyp][i] * calojet;
     }
     return trk_plus_calo_jet.Et();
#endif
     double trkjet_sumet = 0;
     for (unsigned int i = 0; i < jets.size(); ++i) {
	  if (jets[i].et() >= 15) {
	       trkjet_sumet += jets[i].et();
	  }
     }
     for (unsigned int i = 0, i_jet = 0; i < cms2.hyp_other_jets_p4()[i_hyp].size(); ++i, ++i_jet) {
	  while (fabs(cms2.hyp_other_jets_p4()[i_hyp][i].Et() - cms2.jets_p4()[i_jet].Et()) < 
		 1e-5 * fabs(cms2.jets_p4()[i_jet].Et()))
	       i_jet++;
	  LorentzVector calojet = cms2.hyp_other_jets_p4()[i_hyp][i] * 
	       cms2.jets_tq_noCorrF()[i_jet]; // use uncorrected Et
	  double em_et = cms2.hyp_other_jets_emFrac()[i_hyp][i] * calojet.Et();
	  if (em_et > 10)
	       trkjet_sumet += em_et;
     }
     double maxet = 0;
     if (jets.size() > 0)
	  maxet = jets[0].et();
     return std::pair<double, double>(maxet, trkjet_sumet);
}*/
