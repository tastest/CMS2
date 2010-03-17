#include "MT2.h"
#include "TMath.h"
#include "math.h"
#include "Math/VectorUtil.h"

//////////////////////////////////////////
// --- leptonic (WW) mt2 definition --- //
//////////////////////////////////////////

// MT2( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2 )
double MT2(
	const float met,
	const float metPhi,
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > v1,
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > v2
){

/* code here follows documentation in MT2Utility.cc */

	// initialize arrays for lepton 1, lepton 2, MET
	double pa[3]    = {0,0,0};
	double pb[3]    = {0,0,0};
	double pmiss[3] = {0,0,0};

//--- Set the masses -----------------------//

	//pa[0] = v1.M();
	//pb[0] = v2.M();
 	pa[0] = v1.M2() >= 0 ? v1.M() : 0.0;	
	pb[0] = v2.M2() >= 0 ? v2.M() : 0.0;

	/* Note on floating point precision & M^2 < 0:

		- Previously ( 2X ) we set pa[0] and pb[0] to their respective particle masses:

			pa[0] = v1.M()
			pb[0] = v2.M()
		
			This is what is recommended by the MT2Utility.cc authors

		- Now ( 3X ) floating point roundoff can lead to m^2 = E^2 - p^2 < 0, we choose to set M to zero for these cases 
		- Remember that many or all doubles were changed to floats at the NtupleMaker level from 2X to 3X to save space.
		- This presents with the following error from calling LorentzVector::M() :

			dlopen error: /home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C.so: undefined symbol: 
			_ZN4ROOT4Math9GenVector5ThrowEPKc
			Load Error: Failed to load Dynamic link library /home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C.so ...
			/home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C_ACLiC_dict.o(.gnu.linkonce.t._ZNK4ROOT4Math9PxPyPzE4DIfE1MEv+0x4e): In function `ROOT::Math::PxPyPzE4D<float>::M() const':
			: undefined reference to `ROOT::Math::GenVector::Throw(char const*)'

			JMu said it best:
			
			"The reason you see a linker error rather than an exception is a delicious morsel
			of ROOT rootiness.  The LorentzVector code is all included (because it's
			templated in all those funny ways), so it doesn't require linking to any
			library.  Except for the exceptions.  Those live in a library.  If you don't
			link it in, you get a linker error."
			
			- For the leptonic (WW) mt2 distribution, we could just set the masses to 0 since we're in the massless limit anyway (excluding taus)
			- The N jets (tt) mt2 distribution is significantly affected because the b quark mass is 4 GeV so the mass must be set properly
	*/

//-----------------------------------------------//

	// set the transverse momenta for the leptons & MET
	pa[1] = (double) v1.Px();
	pa[2] = (double) v1.Py();
	pb[1] = (double) v2.Px();
	pb[2] = (double) v2.Py();
	pmiss[0] = 0.0;	// not used
	pmiss[1] = (double) met*cos(metPhi);
	pmiss[2] = (double) met*sin(metPhi);

	// instantiate mt2 class, set momenta, and return mt2 
	mt2_bisect::mt2 mt2_event;
	mt2_event.set_momenta( pa, pb, pmiss );
	mt2_event.set_mn( 0.0 );
	return mt2_event.get_mt2();
}

//////////////////////////////////////////
// --- N Jets (tt) mt2 definition   --- //
//////////////////////////////////////////

// MT2J( MET_MAGNITUDE, MET_PHI, P4_LEPTON_1, P4_LEPTON_2, VECT_P4_Jets )
double MT2J(
	const float met,
	const float metPhi,
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_1,
	const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_2,
	const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vect_p4_jets
){
	if( vect_p4_jets.size() < 2 ){
		cout << "MT2.cc: error MT2J called with less than 2 jets... returning mt2 value of -1.0" << endl;
		return -1.0; 
	}
	double mt2_min = std::numeric_limits<double>::max();
	for(unsigned int j1=0; j1 < vect_p4_jets.size(); j1++){
	for(unsigned int j2=0; j2 < vect_p4_jets.size(); j2++){
		if(j1==j2) continue;
		double mt2_j1_j2 = MT2( met, metPhi, p4_lepton_1 + vect_p4_jets.at(j1), p4_lepton_2 + vect_p4_jets.at(j2) );
		if(mt2_j1_j2 < mt2_min){
			mt2_min = mt2_j1_j2;
		}
	}}
	return mt2_min;
}



Bool_t comparePt ( 
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lv1,
  ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lv2 
){
  return lv1.pt() > lv2.pt();
}
vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > CleanJets(
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vect_p4_jets,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_1,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > p4_lepton_2,
  const float jet_pt_threshold,
  const float jet_eta_threshold,
  const float jet_lepton_dR_veto_cone
){
  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > cleaned_vect_p4_jets;
  for (unsigned int ijet = 0; ijet < vect_p4_jets.size(); ijet++) {                                                   // loop on jets supplied by user
    if( vect_p4_jets.at(ijet).pt() < jet_pt_threshold ) continue;                                                     // apply pt threshold
    if( jet_eta_threshold > 0 ){                                                                                      // if the eta threshold is > 0
      if( fabs( vect_p4_jets.at(ijet).eta() ) > jet_eta_threshold ) continue;                                         // apply eta threshold
    }
    if( jet_lepton_dR_veto_cone > 0 ){
      if( ROOT::Math::VectorUtil::DeltaR( vect_p4_jets.at(ijet), p4_lepton_1 ) < jet_lepton_dR_veto_cone ) continue;  // apply lepton1-jet dR vetos
      if( ROOT::Math::VectorUtil::DeltaR( vect_p4_jets.at(ijet), p4_lepton_2 ) < jet_lepton_dR_veto_cone ) continue;  // apply lepton2-jet dR vetos
    }
    cleaned_vect_p4_jets.push_back( vect_p4_jets.at(ijet) );                                                          // keep jets within thresholds
  }
  // sort jets by pt
  sort( cleaned_vect_p4_jets.begin(), cleaned_vect_p4_jets.end(), comparePt );                                        // sort jets by pt
  return cleaned_vect_p4_jets;
}


