#include "MT2.h"

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
	if( v1.M2() >= 0 ){
		pa[0] = (double) v1.M();
	} else {
		pa[0] = (double) -1*sqrt( -1*v1.M2() );
	}
	//pb[0] = v2.M();
	if( v2.M2() >= 0 ){
		pb[0] = (double) v2.M();
	} else {
		pb[0] = (double) -1*sqrt( -1*v2.M2() );
	}

	/* Note on floating point precision & M^2 < 0:

		- Previously ( 2X ) we set pa[0] and pb[0] to their respective particle masses:

			pa[0] = v1.M()
			pb[0] = v2.M()
		
			This is what is recommended by the MT2Utility.cc authors

		- Now ( 3X ) when attempting to set pa[0] and pb[0] using LorentzVector::M() I get the following error:

			dlopen error: /home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C.so: undefined symbol: 
			_ZN4ROOT4Math9GenVector5ThrowEPKc
  		Load Error: Failed to load Dynamic link library /home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C.so
			.
			.
			.
			/home/users/dbarge/cvs/UserCode/JRibnik/CMS2/NtupleMacros/CORE/mt2/example_looper/./ScanChain_C_ACLiC_dict.o(.gnu.linkonce.t._ZNK4ROOT4Math9PxPyPzE4DIfE1MEv+0x4e): In function `ROOT::Math::PxPyPzE4D<float>::M() const':
			: undefined reference to `ROOT::Math::GenVector::Throw(char const*)'

			Remember that many or all doubles were changed to floats at the NtupleMaker level from 2X to 3X to save space.

			Jmu said it best:

			"As for your error, it's floating point roundoff.  You're asking for the mass.
			So first you calculate E * E - p * p.  E and p are both big numbers, so you can
			get a negative result if they're close enough together.  Then when you ask for
			m(), the LorentzVector class sees that m^2 < 0 and throws an exception.
			
			The reason you see a linker error rather than an exception is a delicious morsel
			of ROOT rootiness.  The LorentzVector code is all included (because it's
			templated in all those funny ways), so it doesn't require linking to any
			library.  Except for the exceptions.  Those live in a library.  If you don't
			link it in, you get a linker error.
			
			That probably doesn't interest you as much as how to solve the problem.
			Basically, ask for m2(), see if it's negative, and if it is, use
			-sqrt(-p4.m2())."

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
