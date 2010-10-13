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

	h_m_els_pt_ee = new TH1F("h_m_els_pt_ee", "h_m_els_pt_ee", 100, 0, 20);
	h_m_els_pt_eb = new TH1F("h_m_els_pt_eb", "h_m_els_pt_eb", 100, 0, 20);
	h_m_els_et_ee = new TH1F("h_m_els_et_ee", "h_m_els_et_ee", 100, 0, 20);
	h_m_els_et_eb = new TH1F("h_m_els_et_eb", "h_m_els_et_eb", 100, 0, 20);

	h_m_elspf_pt_ee = new TH1F("h_m_elspf_pt_ee", "h_m_elspf_pt_ee", 100, 0, 20);
	h_m_elspf_pt_eb = new TH1F("h_m_elspf_pt_eb", "h_m_elspf_pt_eb", 100, 0, 20);
	h_m_elspf_et_ee = new TH1F("h_m_elspf_et_ee", "h_m_elspf_et_ee", 100, 0, 20);
	h_m_elspf_et_eb = new TH1F("h_m_elspf_et_eb", "h_m_elspf_et_eb", 100, 0, 20);



	h_m_els = new TH1F("h_m_els", "e-e mass", 100, 0, 5);
	h_m_eltrks = new TH1F("h_m_eltrks", "etrk-etrk mass", 100, 0, 5);
	h_m_sctrks = new TH1F("h_m_sctrks", "sctrk-sctrk mass", 100, 0, 5);

	h_m_elsctrk = new TH1F("h_m_elsctrk", "el-sctrk mass", 100, 0, 5);

	h1_mass_eltrkeltrk_ = new TH1F("h1_mass_eltrkeltrk_", "h1_mass_eltrkeltrk", 100, 0, 5);
        h1_mass_sctrkeltrk_ = new TH1F("h1_mass_sctrkeltrk_", "h1_mass_sctrkeltrk", 100, 0, 5);
        h1_mass_sctrksctrk_ = new TH1F("h1_mass_sctrksctrk_", "h1_mass_sctrksctrk", 100, 0, 5);

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

		int els_type = cms2.els_type()[i];
		if (fabs(cms2.els_etaSC()[i]) > 1.5) {
			h_m_els_pt_ee->Fill(cms2.els_trk_p4()[i].Pt());
			h_m_els_et_ee->Fill(cms2.els_p4()[i].Pt());
			if (els_type & (1<<3)) {
				h_m_elspf_pt_ee->Fill(cms2.els_trk_p4()[i].Pt());
				h_m_elspf_et_ee->Fill(cms2.els_p4()[i].Pt());
			}
		}
		else {
			h_m_els_pt_eb->Fill(cms2.els_trk_p4()[i].Pt());
			h_m_els_et_eb->Fill(cms2.els_p4()[i].Pt());
			if (els_type & (1<<3)) {
				h_m_elspf_pt_eb->Fill(cms2.els_trk_p4()[i].Pt());
				h_m_elspf_et_eb->Fill(cms2.els_p4()[i].Pt());
			}
		}

		// 	  if (abs(cms2.els_mc_id()[i]) != 11)
		// 	       continue;


		for (int j = i + 1; j < cms2.els_p4().size(); ++j) {


			// 	       if (abs(cms2.els_mc_id()[j]) != 11)
			// 		    continue;


			h_m_els->Fill((cms2.els_p4()[i] + cms2.els_p4()[j]).M());
			h_m_eltrks->Fill((cms2.els_trk_p4()[i] + cms2.els_trk_p4()[j]).M());


		}
	}


	// loop on tracks and make sc "candidates"
	std::vector<int> candFlags;
	std::vector<int> tkIndices;
	std::vector<int> elIndices;
	// bit (1) matches tk
	// bit (2) matches pflow electron
	for (int s = 0; s < cms2.scs_p4().size(); ++s) {
		
		int candFlag = 0;

		// check if this super cluster matches a track
		//
		int index = -1;
		double mindr = 999.9;
	        for (int t = 0; t < cms2.trks_trk_p4().size(); ++t) {
                        double dr = DeltaR(t, cms2.scs_p4()[s]);
                        if (dr < mindr) {
				index = t;
                                mindr = dr;
                        }
		}	

		if (mindr < 0.1) {
			candFlag |= (1<<1);
			tkIndices.push_back(index);
		}
		else tkIndices.push_back(-1);

		// check if this track matches a pflow electron
		//
                index = -1;
                mindr = 999.9;
                for (int i = 0; i < cms2.els_p4().size(); ++i) {

			// check it is a pflow electron
			if (!(cms2.els_type()[i] & (1<<3))) continue;
			double dEta = cms2.scs_eta()[s] - cms2.els_etaSC()[i];
			double dPhi = ROOT::Math::VectorUtil::Phi_mpi_pi(cms2.scs_phi()[s] - cms2.els_phiSC()[i]);
			double dr = sqrt(dEta*dEta + dPhi*dPhi);
			//double dr = ROOT::Math::VectorUtil::DeltaR(cms2.scs_p4()[s], cms2.els_p4()[i]);
                        if (dr < mindr) {
                                index = i;
                                mindr = dr;
                        }
                }

                if (mindr < 0.1) {
                        candFlag |= (1<<2);
                        elIndices.push_back(index);
                }
                else elIndices.push_back(-1);

		// set the flags for this cand
		//
		candFlags.push_back(candFlag);

	} // end loop on tracks


	// now do loop on sc cands
        for (int i = 0; i < cms2.scs_p4().size(); ++i) {
		for (int j = i + 1; j < cms2.scs_p4().size(); ++j) {

			// both are sc-trk
			if ((candFlags[i]&(1<<1)) == (1<<1) && (candFlags[j]&(1<<1)) == (1<<1)) {
				h1_mass_sctrksctrk_->Fill((
					cms2.trks_trk_p4()[tkIndices[i]] 
					+ cms2.trks_trk_p4()[tkIndices[j]]).M());
			}

			// both are el
                        if ((candFlags[i]&(1<<2)) == (1<<2)  && (candFlags[j]&(1<<2)) == (1<<2)) {
                                h1_mass_eltrkeltrk_->Fill((
                                        cms2.els_trk_p4()[elIndices[i]]
                                        + cms2.els_trk_p4()[elIndices[j]]).M());
                        }

			// one is sc-trk and other is el
                        if ((candFlags[i]&(1<<2)) == (1<<2) && (candFlags[j]&(1<<1)) == (1<<1)) {
                                h1_mass_sctrkeltrk_->Fill((
                                        cms2.trks_trk_p4()[tkIndices[i]]
                                        + cms2.els_trk_p4()[elIndices[j]]).M());
                        }
                        if ((candFlags[j]&(1<<2) == (1<<2)) && (candFlags[i]&(1<<1)) == (1<<1)) {
                                h1_mass_sctrkeltrk_->Fill((
                                        cms2.trks_trk_p4()[tkIndices[j]]
                                        + cms2.els_trk_p4()[elIndices[i]]).M());
                        }




		}
	}	




	// loop on superclusters
	for (int i = 0; i < cms2.scs_p4().size(); ++i) {
		float scEnergy1 = 0.0;
		LorentzVector trk1;
		double mindr1 = 0.1;
		for (int itrk = 0; itrk < cms2.trks_trk_p4().size(); ++itrk) {
			double dr = DeltaR(itrk, cms2.scs_p4()[i]);
			if (dr < mindr1) {
				trk1 = cms2.trks_trk_p4()[itrk];
				scEnergy1 = cms2.scs_energy()[i];
				mindr1 = dr;
			}
		}
		// loop on superclusters excluding the one already
		// considered in the previous loop
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

			// now consider pflow electron for a second leg
			// need to make sure it is not the same cluster as above

			for (int j = 0; j < cms2.els_p4().size(); ++j) {

				// is a pflow electron with a different cluster
				// to the selected sc-trk pair
				if (!(cms2.els_type()[j] & (1<<3))) continue;
				if (cms2.els_ecalEnergy()[j] == scEnergy1) {
					//std::cout << "found same" << std::endl;
					continue;
				}
				// fill sc-trk - electron mass hist
				h_m_elsctrk->Fill((cms2.els_trk_p4()[j] + trk1).M());
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
