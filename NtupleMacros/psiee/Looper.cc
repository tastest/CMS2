#include <math.h>
#include "TCanvas.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname),
       n_scs_el_matched(0),
       n_scs_el_not_matched(0)
{

}

void Looper::BookHistos ()
{
     h_m_els = new TH1F("h_m_els", "e-e mass", 100, 0, 5);
     h_m_eltrks = new TH1F("h_m_eltrks", "etrk-etrk mass", 100, 0, 5);
     h_m_sctrks = new TH1F("h_m_sctrks", "sctrk-sctrk mass", 100, 0, 5);
     h_m_esc = new TH1F("h_m_esc", "etrk-sc mass", 100, 0, 5);
     h_m_esctrks = new TH1F("h_m_esctrk", "etrk-sctrk mass", 100, 0, 5);
     h_m_esctrks_elmatched = new TH1F("h_m_esctrks_elmatched", "etrk-sctrk mass", 100, 0, 5);
     h_m_esctrks_elnotmatched = new TH1F("h_m_esctrks_elnotmatched", "etrk-sctrk mass", 100, 0, 5);
     h_scs_elsidx = new TH1F("h_scs_elsidx", "elsidx", 12, -1.5, 10.5);
     h_ncands = new TH1F("h_ncands", "N(#psi cands)", 11, -0.5, 10.5);
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
     // e-e mass
     for (int i = 0; i < cms2.els_p4().size(); ++i) {
// 	  if (abs(cms2.els_mc_id()[i]) != 11)
// 	       continue;
	  for (int j = i + 1; j < cms2.els_p4().size(); ++j) {
// 	       if (abs(cms2.els_mc_id()[j]) != 11)
// 		    continue;
	       if (cms2.els_charge()[i] * cms2.els_charge()[j] > 0)
		    continue;
	       h_m_els->Fill((cms2.els_p4()[i] + cms2.els_p4()[j]).M());
	       h_m_eltrks->Fill((cms2.els_trk_p4()[i] + cms2.els_trk_p4()[j]).M());
	  }
     }
     // sc-sc mass
     for (int i = 0; i < cms2.scs_p4().size(); ++i) {
	  LorentzVector trk1;
	  double mindr1 = 0.1;
	  int imin1 = -1;
	  for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
	       double dr = DeltaR(itrk, cms2.scs_p4()[i]);
	       if (dr < mindr1) {
		    trk1 = cms2.trks_trk_p4()[itrk];
		    mindr1 = dr;
		    imin1 = itrk;
	       }
	  }
 	  for (int j = i + 1; j < cms2.scs_p4().size(); ++j) {
	       LorentzVector trk2;
	       double mindr2 = 0.1;
	       int imin2 = -1;
	       for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
		    double dr = DeltaR(itrk, cms2.scs_p4()[j]);
		    if (dr < mindr2) {
			 trk2 = cms2.trks_trk_p4()[itrk];
			 mindr2 = dr;
			 imin2 = itrk;
		    }
	       }
	       if (mindr1 < 0.1 && mindr2 < 0.1) {
		    if (cms2.trks_charge()[imin1] * cms2.trks_charge()[imin2] > 0)
			 continue;
		    double m2 = (trk1 + trk2).M2();
		    if (m2 > 0)
			 h_m_sctrks->Fill(sqrt(m2));
	       }
	  }
     }
     // e-sc mass
     for (int i = 0; i < cms2.els_p4().size(); ++i) {
 	  for (int j = 0; j < cms2.scs_p4().size(); ++j) {
	       double m2 = (cms2.els_trk_p4()[i] + cms2.scs_p4()[j]).M2();
	       if (m2 > 0)
		    h_m_esc->Fill(sqrt(m2));
	  }
     }	       
     // e-sctrk mass
     // count psi candidates as well
     int ncands = 0;
     for (int i = 0; i < cms2.els_p4().size(); ++i) {
 	  for (int j = 0; j < cms2.scs_p4().size(); ++j) {
	       LorentzVector trk2;
	       double mindr2 = 0.1;
	       int imin2 = -1;
	       for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
		    double dr = DeltaR(itrk, cms2.scs_p4()[j]);
		    if (dr < mindr2) {
			 trk2 = cms2.trks_trk_p4()[itrk];
			 mindr2 = dr;
			 imin2 = itrk;
		    }
	       }
	       double m2 = (cms2.els_trk_p4()[i] + trk2).M2();
	       if (m2 > 0) {
		    const double m = sqrt(m2);
		    h_m_esctrks->Fill(m);
		    if (m > 2.5 && m < 5) {
			 h_scs_elsidx->Fill(cms2.scs_elsidx()[j]);
			 ncands++;
		    }
		    if (cms2.scs_elsidx()[j] != -9999) {
			 h_m_esctrks_elmatched->Fill(m);
		    } else {
			 h_m_esctrks_elnotmatched->Fill(m);
		    }
	       }
	  }
     }	       
     h_ncands->Fill(ncands);
}

void Looper::End ()
{
     TCanvas *c = new TCanvas;
     c->Divide(4, 2);
     int i = 0;
     c->cd(++i);
     h_m_els->Draw();
     c->cd(++i);
     h_m_eltrks->Draw();
     c->cd(++i);
     h_m_sctrks->Draw();
     c->cd(++i);
     h_m_esctrks->Draw();
     c->cd(++i);
     h_m_esctrks_elmatched->Draw();
     c->cd(++i);
     h_m_esctrks_elnotmatched->Draw();
     c->cd(++i);
     h_scs_elsidx->Draw();
     c->cd(++i);
     h_ncands->Draw();
}
