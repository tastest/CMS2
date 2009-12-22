#include <math.h>
#include "TCanvas.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
}

void Looper::BookHistos ()
{
     h_m_els = new TH1F("h_m_els", "e-e mass", 100, 0, 5);
     h_m_eltrks = new TH1F("h_m_eltrks", "etrk-etrk mass", 100, 0, 5);
     h_m_sctrks = new TH1F("h_m_sctrks", "sctrk-sctrk mass", 100, 0, 5);
}


bool Looper::FilterEvent()
{ 
  // duplicate filter, based on trk information and dilepton hyp
  //
  // comment in following lines
  // 

     if (cms2.trks_d0().size() == 0)
	  return true;
     DorkyEventIdentifier id(cms2);
     if (is_duplicate(id)) {
//     duplicates_total_n_++;
//     duplicates_total_weight_ += cms2.evt_scale1fb();
	  cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
	  return true;
     }
     return false; 
}

double Looper::DeltaR (int itrk, const LorentzVector &v)
{
     const double dphi = ROOT::Math::VectorUtil::Phi_mpi_pi(v.phi() - cms2.trks_outer_position()[itrk].phi());
     const double deta = fabs(v.eta() - cms2.trks_outer_position()[itrk].eta());
     return sqrt(deta * deta + dphi * dphi);
}

void Looper::FillEventHistos ()
{
//      for (int i = 0; i < cms2.scs_p4().size(); ++i) {
// 	  for (int j = i + 1; j < cms2.scs_p4().size(); ++j) {
// 	       h_m->Fill((cms2.scs_p4()[i] + cms2.scs_p4()[j]).M());
// 	  }
//      }
     for (int i = 0; i < cms2.els_p4().size(); ++i) {
// 	  if (abs(cms2.els_mc_id()[i]) != 11)
// 	       continue;
	  for (int j = i + 1; j < cms2.els_p4().size(); ++j) {
// 	       if (abs(cms2.els_mc_id()[j]) != 11)
// 		    continue;
	       h_m_els->Fill((cms2.els_p4()[i] + cms2.els_p4()[j]).M());
	       h_m_eltrks->Fill((cms2.els_trk_p4()[i] + cms2.els_trk_p4()[j]).M());
	  }
     }
     for (int i = 0; i < cms2.scs_p4().size(); ++i) {
	  LorentzVector trk1;
	  double mindr1 = 0.1;
	  for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
	       double dr = DeltaR(itrk, cms2.scs_p4()[i]);
	       if (dr < mindr1) {
		    trk1 = cms2.trks_trk_p4()[itrk];
		    mindr1 = dr;
	       }
	  }
 	  for (int j = i + 1; j < cms2.scs_p4().size(); ++j) {
	       LorentzVector trk2;
	       double mindr2 = 0.1;
	       for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
		    double dr = DeltaR(itrk, cms2.scs_p4()[j]);
		    if (dr < mindr2) {
			 trk2 = cms2.trks_trk_p4()[itrk];
			 mindr2 = dr;
		    }
	       }
	       if (mindr1 < 0.1 && mindr2 < 0.1) {
		    double m2 = (trk1 + trk2).M2();
		    if (m2 > 0)
			 h_m_sctrks->Fill(sqrt(m2));
	       }
	  }
     }
     
}

void Looper::End ()
{
     new TCanvas;
     h_m_els->Draw();
     new TCanvas;
     h_m_eltrks->Draw();
     new TCanvas;
     h_m_sctrks->Draw();
}
