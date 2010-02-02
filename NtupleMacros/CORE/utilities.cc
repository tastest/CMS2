#include <algorithm>
#include <set>
#include "Math/VectorUtil.h"
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"
#include "TH1F.h"
#include "TFile.h"
#include "TMath.h"
#include "TList.h"
#include "TRegexp.h"
#include "TDirectory.h"
#include "TSystem.h"

#include "Math/LorentzVector.h"
#include "TMath.h"
#include <vector>
#include "TDatabasePDG.h"

// CMS2 includes
#include "CMS2.h"
#include "utilities.h"

using std::vector;
TH2D* rfhist = 0;

double dRbetweenVectors(const LorentzVector &vec1, 
			const LorentzVector &vec2 ){ 

  double dphi = std::min(::fabs(vec1.Phi() - vec2.Phi()), 2 * M_PI - fabs(vec1.Phi() - vec2.Phi()));
  double deta = vec1.Eta() - vec2.Eta();
  return sqrt(dphi*dphi + deta*deta);
}

int match4vector(const LorentzVector &lvec, 
		 const vector<LorentzVector> &vec, 
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

/* int num4vectorsmatched(ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lvec, 
		 vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > vec, 
		 double cut=0.4, double ptcut=0.0 ){

  if( vec.size() == 0 ) return 0;
  //cout << "size of vec = " << vec.size() << endl;
  double dR = cut; 
  int result = 0;
  for ( unsigned int i=0; i < vec.size();++i) {
    if ( vec[i].pt() > ptcut && dRbetweenVectors(lvec,vec[i]) < dR ) ++result;
  }
  return result;
}
*/

std::vector<LorentzVector> p4sInCone(const LorentzVector &refvec, 
				     const std::vector<LorentzVector> &invec, 
				     double coneSize=0.5 ) 
{
  vector<LorentzVector> result;
  if ( invec.size() == 0 ) return result;

  double dR = coneSize; 
  double x = 0.0;
  for ( unsigned int i=0; i < invec.size();++i) {
    x = dRbetweenVectors(refvec,invec[i]);
    if (x < dR ) {result.push_back(invec[i]);}
  }
  return result;
}

std::vector<unsigned int> idxInCone(const LorentzVector &refvec, 
				    const std::vector<LorentzVector> &invec, 
				    double coneSize=0.5 ) 
{
  vector<unsigned int > result;
  if ( invec.size() == 0 ) return result;

  double dR = coneSize; 
  for ( unsigned int i=0; i < invec.size();++i) {
    if (dRbetweenVectors(refvec,invec[i]) < dR ) {result.push_back(i);}
  }
  return result;
}

double trkIsolation(int trk_index) {
  //
  // calculate track isolation following electron isolation definition in ElectronMaker.cc
  //
  // sum up all track.pt around track if track fulfills:
  // dR < 0.3
  // dR > 0.01
  // d0 < 0.1
  // dZ < 0.5
  // pT >= 1.0
  // nHits > 7

  float dRConeMin = 0.01;
  float dRConeMax = 0.3;
  float vtxDiffZMax = 0.5;
  float tkVtxDMax = 0.1;
  float ptMin = 1.0;
  int nHits = 7;

  double isoResult = -10.;
  if( cms2.trks_trk_p4().size() == 0 ) {
    std::cout << "Configuration Error: track collection is not set!" <<std::endl;
    return isoResult;
  }

  double sumPt = 0;

  for ( unsigned int trk = 0;
	trk < cms2.trks_trk_p4().size();
	++trk ) {
    if ( cms2.trks_z0()[trk] > tkVtxDMax ) continue;
    double pt = cms2.trks_trk_p4()[trk].pt();
    if (  pt < ptMin ) continue;
    if ( cms2.trks_validHits()[trk] <= nHits ) continue;
    double dR = ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[trk_index], cms2.trks_trk_p4()[trk]);
    if (dR < dRConeMin) continue;
    if ( dR > dRConeMax ) continue;
    double dZ = TMath::Abs(cms2.trks_z0()[trk_index] - cms2.trks_z0()[trk]);
    if ( dZ >= vtxDiffZMax) continue;
    sumPt += pt;
  }

  return sumPt;
}
