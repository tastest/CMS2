#include <algorithm>
#include <set>
#include <utility>
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

// this is meant to be passed as the third argument, the predicate, of the
// standard library sort algorithm
bool sortByPt(const LorentzVector &vec1, 
			const LorentzVector &vec2 )
{
    return vec1.pt() > vec2.pt();
}

// this is a workaround for not having unique event id's in MC
class DorkyEvent
{
    public:
        DorkyEvent()
        {
            run_      = cms2.evt_run();
            lumi_     = cms2.evt_lumiBlock();
            event_    = cms2.evt_event();
            trks_d0_  = cms2.trks_d0().size()     ? cms2.trks_d0()[0]           : -9999.;
            trks_pt_  = cms2.trks_trk_p4().size() ? cms2.trks_trk_p4()[0].pt()  : -9999.;
            trks_eta_ = cms2.trks_trk_p4().size() ? cms2.trks_trk_p4()[0].eta() : -9999.;
            trks_phi_ = cms2.trks_trk_p4().size() ? cms2.trks_trk_p4()[0].phi() : -9999.;
        }
        ~DorkyEvent() {}

        bool operator < (const DorkyEvent &) const;
        bool operator == (const DorkyEvent &) const;

        unsigned int run      () const { return run_;      }
        unsigned int lumi     () const { return lumi_;     }
        unsigned int event    () const { return event_;    }
        float        trks_d0  () const { return trks_d0_;  }
        float        trks_pt  () const { return trks_pt_;  }
        float        trks_eta () const { return trks_eta_; }
        float        trks_phi () const { return trks_phi_; }

    private:
        unsigned int run_, lumi_, event_;
        float trks_d0_, trks_pt_, trks_eta_, trks_phi_;
};

class DorkyEventIdentifier
{
    public:
        DorkyEventIdentifier()
        {
            already_seen.clear();
            duplicates_total_n = 0;
            duplicates_total_weight = 0.;
        }
        ~DorkyEventIdentifier() {}

        bool is_duplicate(const DorkyEvent &id)
        {
            std::pair<std::set<DorkyEvent>::const_iterator, bool> ret =
                already_seen.insert(id);

            if (! ret.second) {
                duplicates_total_n++;
                duplicates_total_weight += cms2.evt_scale1fb();
                cout << "Duplicate event found. Run: " << ret.first->run() << ", Lumi: " << ret.first->lumi() << ", Event: " << ret.first->event() << endl;
                cout.precision(10);
                cout << "\td0:\t"  << ret.first->trks_d0()  << endl;
                cout << "\tpt:\t"  << ret.first->trks_pt()  << endl;
                cout << "\teta:\t" << ret.first->trks_eta() << endl;
                cout << "\tphi:\t" << ret.first->trks_phi() << endl;
            }

            return ! ret.second;
        }

    private:
        std::set<DorkyEvent> already_seen;
        int duplicates_total_n;
        double duplicates_total_weight;
};

bool DorkyEvent::operator < (const DorkyEvent &other) const
{
    if (run() != other.run())
        return run() < other.run();
    if (event() != other.event())
        return event() < other.event();
    // the floating point numbers are not easy, because we're
    // comapring ones that are truncated (because they were written
    // to file and read back in) with ones that are not truncated.
    if (fabs(trks_d0()  - other.trks_d0())  > 1e-6 * trks_d0())
        return trks_d0() < other.trks_d0();
    if (fabs(trks_pt()  - other.trks_pt())  > 1e-6 * trks_pt())
        return trks_pt() < other.trks_pt();
    if (fabs(trks_eta() - other.trks_eta()) > 1e-6 * trks_eta())
        return trks_eta() < other.trks_eta();
    if (fabs(trks_phi() - other.trks_phi()) > 1e-6 * trks_phi())
        return trks_phi() < other.trks_phi();
    // if the records are exactly the same, then r1 is not less than
    // r2.  Duh!
    return false;
}

bool DorkyEvent::operator == (const DorkyEvent &other) const
{
    if (run() != other.run())
        return false;
    if (event() != other.event())
        return false;
    // the floating point numbers are not easy, because we're
    // comapring ones that are truncated (because they were written
    // to file and read back in) with ones that are not truncated.
    if (fabs(trks_d0()  - other.trks_d0())  > 1e-6 * trks_d0())
        return false;
    if (fabs(trks_pt()  - other.trks_pt())  > 1e-6 * trks_pt())
        return false;
    if (fabs(trks_eta() - other.trks_eta()) > 1e-6 * trks_eta())
        return false;
    if (fabs(trks_phi() - other.trks_phi()) > 1e-6 * trks_phi())
        return false;
    return true;
}
