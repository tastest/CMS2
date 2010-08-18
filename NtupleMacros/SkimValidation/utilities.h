#ifndef UTILITIES_H
#define UTILITIES_H

#include "TH1F.h"
#include "TH2D.h"
#include "TVector3.h"
#include <algorithm>
#include <set>
#include "Math/VectorUtil.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

//TVector3 correctMETforTracks() ;
double trkIsolation(int trk_index);
double dRbetweenVectors(const LorentzVector & vec1, 
			const LorentzVector & vec2 );
bool   sortByPt(const LorentzVector & vec1, const LorentzVector & vec2 );

//class DorkyEvent;
//class DorkyEventIdentifier;

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
                //                duplicates_total_weight += cms2.evt_scale1fb();
/*                 cout << "Duplicate event found. Run: " << ret.first->run() << ", Lumi: " << ret.first->lumi() << ", Event: " << ret.first->event() << endl; */
/*                 cout.precision(10); */
/*                 cout << "\td0:\t"  << ret.first->trks_d0()  << endl; */
/*                 cout << "\tpt:\t"  << ret.first->trks_pt()  << endl; */
/*                 cout << "\teta:\t" << ret.first->trks_eta() << endl; */
/*                 cout << "\tphi:\t" << ret.first->trks_phi() << endl; */
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
    if (lumi() != other.lumi())
        return lumi() < other.lumi();
    if (event() != other.event())
        return event() < other.event();
    // the floating point numbers are not easy, because we're
    // comapring ones that are truncated (because they were written
    // to file and read back in) with ones that are not truncated.
    
    if( 42 != 42 ) { // this does not work if we compare 2 different processings
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
    }

    return false;
}

bool DorkyEvent::operator == (const DorkyEvent &other) const
{
    if (run() != other.run())
        return false;
    if (lumi() != other.lumi())
        return false;
    if (event() != other.event())
        return false;
    // the floating point numbers are not easy, because we're
    // comapring ones that are truncated (because they were written
    // to file and read back in) with ones that are not truncated.
    if( 42 != 42 ) { // this does not work if we compare 2 different processings
      if (fabs(trks_d0()  - other.trks_d0())  > 1e-6 * trks_d0())
        return false;
      if (fabs(trks_pt()  - other.trks_pt())  > 1e-6 * trks_pt())
        return false;
      if (fabs(trks_eta() - other.trks_eta()) > 1e-6 * trks_eta())
        return false;
      if (fabs(trks_phi() - other.trks_phi()) > 1e-6 * trks_phi())
        return false;
    }
    return true;
}


extern TH2D *rfhist;

#endif
