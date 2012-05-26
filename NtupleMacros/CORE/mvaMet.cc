//
// MVA MET 
#include <math.h>
#include <algorithm>
#include "TMath.h"
#include "TVector2.h"

#include "CMS2.h"
#include "trackSelections.h"
#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"

#include "mvaMet.h"

#include "../Tools/MVAMET/MetUtilities.h" 
#include "jetSelections.h"
#include "TLorentzVector.h"

//
// dZ of a pfcandidate to the primary vertex
//
double pfCandDz(const unsigned ipfcand, const unsigned ipv) {

	double lDz = -999;
	if(cms2.pfcands_trkidx()[ipfcand]>=0) 	lDz = fabs( trks_dz_pv( cms2.pfcands_trkidx()[ipfcand], ipv).first );
	/*  
		if( abs(cms2.pfcands_particleId()[ipfcand])==11   &&
		cms2.pfcands_pfelsidx()[ipfcand]>=0          &&
		cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx()[ipfcand])>=0) 
		{
		int gsfid = cms2.els_gsftrkidx().at(cms2.pfels_elsidx().at(cms2.pfcands_pfelsidx()[ipfcand]));
		lDz  = fabs(gsftrks_dz_pv_local( gsfid, ipv, false ).first);
		}  
	 */  
	return lDz;
}
		

//
// This function is a wrapper for GetMet(...) 
// Inputs to GetMet are prepared in the function 
//
float mvaMET( int hyp_index, double deltaZCut, MVAMet* mvaMet,  
			  FactorizedJetCorrector *jet_corrector )
{
	float mvamet_=-999;

	// Make Generic Objects
	std::vector<MVAMet::LorentzVector>                            lVisible;
	std::vector<std::pair<MVAMet::LorentzVector,double> >         lPFInfo;
	std::vector<MetUtilities::JetInfo>                            lJetInfo;
	std::vector<MVAMet::Vector>                                   lVtxInfo;

	//
	// Prepare PFCandidates and their dZ 
	//
	for(unsigned int i0 = 0; i0 < cms2.pfcands_p4().size(); ++i0)
	{
		double pDZ = -999;
		if(cms2.evt_nvtxs() > 0) pDZ  = pfCandDz(i0, 0); //If there is no vertex return negative number -999 for dZ
		std::pair<LorentzVector,double> pPFObject(cms2.pfcands_p4().at(i0),pDZ);
		lPFInfo.push_back(pPFObject);
	}

	//
	// Prepare jets     
	// 
	// This shouldn't be a double loop. I just copied from CMSSW.
	// Originally, the reason for double loop in CMSSW was to get JEC from the two jet collections.
	// Here in cms2, we don't need to do that way. I will fix this in the next update. 
	for(unsigned int i0   = 0; i0 < cms2.pfjets_p4().size(); i0++) {   // uncorrected jets collection                                                                                                        
		LorentzVector uncorjet = cms2.pfjets_p4().at(i0);
		for( unsigned int i1 = 0; i1 < cms2.pfjets_p4().size(); i1++) {   // corrected jets collection                                                                                                   
			// 
			// JEC 
			// 
			double jec = 1.0;
			jet_corrector->setRho(cms2.evt_ww_rho_vor());
			jet_corrector->setJetA(cms2.pfjets_area().at(i1));
			jet_corrector->setJetPt(cms2.pfjets_p4()[i1].pt());
			jet_corrector->setJetEta(cms2.pfjets_p4()[i1].eta());
			double corr = jet_corrector->getCorrection();
			jec *= corr; 
			LorentzVector corjet = jec*cms2.pfjets_p4().at(i1);

			if( cms2.pfjets_area().at(i0) != cms2.pfjets_area().at(i1) ) continue; 
			if( fabs( uncorjet.Eta() - corjet.Eta()) > 0.01 ) continue;
			if( corjet.Pt()  < 1 ) continue; 
			if( !passesPFJetID(i1) ) continue;
			double lJec = corjet.Pt()/uncorjet.Pt();

			// precomputed JetId MVA value
			double lMVA = cms2.pfjets_mvavalue().at(i1);

			double lNeuFrac = ( cms2.pfjets_neutralEmE().at(i1) + cms2.pfjets_neutralHadronE().at(i1) ) / cms2.pfjets_p4().at(i1).E();

			MetUtilities::JetInfo pJetObject;
			pJetObject.p4       = corjet;
			pJetObject.mva      = lMVA;
			pJetObject.neutFrac = lNeuFrac;
			lJetInfo.push_back(pJetObject);
			break;
		}
	}

	//
	// make Vertices 
	//
	for(unsigned int i0 = 0; i0 < cms2.vtxs_position().size(); i0++)
	{
		MVAMet::Vector pVec;
		pVec.SetCoordinates(cms2.vtxs_position().at(i0).px(),
				cms2.vtxs_position().at(i0).py(),
				cms2.vtxs_position().at(i0).pz() );
		lVtxInfo.push_back(pVec);
	}

	//
	// calculate MET 
	//

	// Dummy visible stuff 
	double lDummyPt0  = 10;
	double lDummyPhi0 = 0;
	double lDummyEta0 = 1;
	double lDummyPt1  = 10;
	double lDummyPhi1 = 2.;
	double lDummyEta1 = -1;
	TLorentzVector lD0; lD0.SetPtEtaPhiM(lDummyPt0,lDummyEta0,lDummyPhi0,0.);
	TLorentzVector lD1; lD1.SetPtEtaPhiM(lDummyPt1,lDummyEta1,lDummyPhi1,0.);

	LorentzVector lVis0; lVis0.SetCoordinates(lD0.Px(),lD0.Py(),lD0.Pz(),lD0.E());
	LorentzVector lVis1; lVis1.SetCoordinates(lD1.Px(),lD1.Py(),lD1.Pz(),lD1.E());
	lVisible.push_back(lVis0);
	lVisible.push_back(lVis1);

	//Calculate the MVA
	std::pair<LorentzVector,double> lMVAMetInfo = mvaMet->GetMet(lVisible,lJetInfo,lPFInfo,lVtxInfo,false);
	//std::cout << "Met---> " << lMVAMetInfo.first.pt() << " -- " << lMVAMetInfo.first.phi() << std::endl;
	std::cout << lMVAMetInfo.first.pt() << std::endl;

	mvamet_ = lMVAMetInfo.first.pt();

	return mvamet_;
}

