//===========================================================
//
// Various matching tools are kept here
//
//============================================================
#include "Math/LorentzVector.h"
#include "TMath.h"
//#include <vector>
#include "TDatabasePDG.h"


double dRbetweenVectors(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec1, 
			ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > vec2 ){ 

  double dphi = TMath::Min(TMath::Abs(vec1.Phi() - vec2.Phi()), 2*TMath::Pi() - TMath::Abs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

int match4vector(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lvec, 
		 vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > vec, 
		 double cut=10.0 ){

  if( vec.size() == 0 ) return -1;
  //cout << "size of vec = " << vec.size() << endl;
  double dR = cut; 
  double x;
  int iret = -1;
  for ( unsigned int i=0; i < vec.size();++i) {
    x = dRbetweenVectors(lvec,vec[i]);
    if (x < dR ) {dR = x; iret = i;}
  }
  return iret;
}

