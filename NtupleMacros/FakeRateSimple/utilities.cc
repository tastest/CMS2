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
#include "../CORE/utilities.h"

using std::vector;
TH2D* rfhist = 0;

//TVector3 correctMETforTracks()
//{
//  if ( ! rfhist ){
//    printf("correctMETforTracks(): loading RF histo\n");
//    TDirectory* currentDir = gDirectory;
//    TFile *metcorr_file = TFile::Open("$CMS2_LOCATION/NtupleMacros/data/metcorr.root", "read");
//    if ( metcorr_file == 0 ) {
//        std::cout << "$CMS2_LOCATION/NtupleMacros/data/metcorr.root could not be found!!" << std::endl;
//	std::cout << "Please make sure that $CMS2_LOCATION points to your CMS2 directory and that" << std::endl;
//	std::cout << "$CMS2_LOCATION/NtupleMacros/data/metcorr.root exists!" << std::endl;
//        gSystem->Exit(1); 
//    }
//    assert( metcorr_file );
//    rfhist = dynamic_cast<TH2D*>(metcorr_file->Get("rf_pt_mbin"));
//    assert( rfhist );
//    rfhist->SetDirectory(0);
//    metcorr_file->Close();
//    currentDir->cd();
//  }
//  // initialize temporary met variables
//  double met_x = 0;
//  double met_y = 0;
//
//  // set kinematic and quality cuts
//  const int    nhits_cut  = 7;
//  const double hoe_cut    = 0.1;
//  const double d0_cut     = 0.05;
//  const double nchisq_cut = 5.;
//  const double eta_cut    = 2.4;
//  const double low_pt_cut = 2.;
//  const double hi_pt_cut  = 100.;
//
//  for( unsigned int trkCount = 0; trkCount < cms2.trks_trk_p4().size(); ++trkCount ) {
//    // skip track if matched to a global muon
//    if( cms2.trk_musidx()[trkCount] != -999 &&  
//	ROOT::Math::VectorUtil::DeltaR(cms2.trks_trk_p4()[trkCount],
//				       cms2.mus_trk_p4()[cms2.trk_musidx()[trkCount]]) < 0.1) {
//      if( cms2.mus_trkidx()[cms2.trk_musidx()[trkCount]]  == int(trkCount) 
//	  // useless && cms2.mus_trkdr()[cms2.trk_musidx()[trkCount]] < 0.1 
//	  ) continue;
//    }
//
//    // skip track if matched to an "electron"
//    if( cms2.trks_elsidx()[trkCount] != -999 && cms2.trks_elsdr()[trkCount] < 0.1) {
//      if( cms2.els_hOverE()[ cms2.trks_elsidx()[trkCount] ] < hoe_cut ) {
//	if( cms2.els_trkidx()[cms2.trks_elsidx()[trkCount]]  == int(trkCount) && cms2.els_trkdr()[cms2.trks_elsidx()[trkCount]] < 0.1 ) continue;
//      }
//    }
//
//    // skip tracks at large eta or with large pt
//    if( TMath::Abs( cms2.trks_trk_p4()[trkCount].eta() ) > eta_cut || cms2.trks_trk_p4()[trkCount].pt() > hi_pt_cut ) continue;
//
//    // skip tracks that do no pass quality cuts
//    if( cms2.trks_validHits()[trkCount] < nhits_cut || ( cms2.trks_chi2()[trkCount] / cms2.trks_ndof()[trkCount] ) > nchisq_cut || TMath::Abs( cms2.trks_d0corr()[trkCount] ) > d0_cut ) continue;
//
//    // correct tracks w/ pt < 2 setting RF = 0
//    if( cms2.trks_trk_p4()[trkCount].pt() < low_pt_cut ) {
//      met_x -= cms2.trks_trk_p4()[trkCount].pt() * cos( cms2.trks_trk_p4()[trkCount].phi() );
//      met_y -= cms2.trks_trk_p4()[trkCount].pt() * sin( cms2.trks_trk_p4()[trkCount].phi() );
//      continue;
//    }
//
//    // skip any remaining tracks that don't have outerEta, outerPhi information
//    if( cms2.trks_outer_p4()[trkCount].x() == -999 || cms2.trks_outer_p4()[trkCount].y() == -999 ) continue;
//
//    // if we've made it this far, get the response from the histogram
//    int bin   = rfhist->FindBin( cms2.trks_trk_p4()[trkCount].eta(), cms2.trks_trk_p4()[trkCount].pt() );
//    double rf = rfhist->GetBinContent( bin );
//
//    // now, correct MET for track using RF
//    met_x +=  ( rf * cms2.trks_trk_p4()[trkCount].P() * ( 1 / cosh( cms2.trks_outer_p4()[trkCount].Eta() ) ) * cos( cms2.trks_outer_p4()[trkCount].Phi() ) - cms2.trks_trk_p4()[trkCount].pt() * cos( cms2.trks_trk_p4()[trkCount].phi() ) );
//    met_y +=  ( rf * cms2.trks_trk_p4()[trkCount].P() * ( 1 / cosh( cms2.trks_outer_p4()[trkCount].Eta() ) ) * sin( cms2.trks_outer_p4()[trkCount].phi() ) - cms2.trks_trk_p4()[trkCount].pt() * sin( cms2.trks_trk_p4()[trkCount].phi() ) );
//  }
//
//  // fill MET vector
//  TVector3 metvec( met_x, met_y, 0);
//  
//  // return corrected MET
//  return metvec;
//}

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
