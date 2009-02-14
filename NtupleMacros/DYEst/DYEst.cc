#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "DYEst/DYEst.h"

DYEst::DYEst (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters
     memset(cands_passing	, 0, sizeof(cands_passing       ));
     memset(cands_passing_w2	, 0, sizeof(cands_passing_w2    ));
     memset(cands_count		, 0, sizeof(cands_count         ));
}

void DYEst::BookHistos ()
{

	// for computing the R ratio / estimating the background
     hnm1_mll_0j_ = new NMinus1Hist(sample_, "mll_0j", 200, 0, 200, cuts_, jet_z_veto_cuts);
     hnm1_mll_1j_ = new NMinus1Hist(sample_, "mll_1j", 200, 0, 200, cuts_, jet_z_veto_cuts);
     hnm1_mll_2j_ = new NMinus1Hist(sample_, "mll_2j", 200, 0, 200, cuts_, jet_z_veto_cuts);

	// for studying the behavior of the R ratio as a function of MET
	// remove the straight met cut for this
     hnm1_met_0j_in_ = new NMinus1Hist(sample_, "met_0j_in", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);
     hnm1_met_1j_in_ = new NMinus1Hist(sample_, "met_1j_in", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);
     hnm1_met_2j_in_ = new NMinus1Hist(sample_, "met_2j_in", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);
     hnm1_met_0j_out_ = new NMinus1Hist(sample_, "met_0j_out", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);
     hnm1_met_1j_out_ = new NMinus1Hist(sample_, "met_1j_out", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);
     hnm1_met_2j_out_ = new NMinus1Hist(sample_, "met_2j_out", 200, 0, 200, cuts_, jet_z_veto_cuts | simple_met_cuts);

        n_WZ_ee_ = 0;
        n_ZZ_ee_ = 0;
        n_Total_ee_ = 0;
        
        n_WZ_mm_ = 0;
        n_ZZ_mm_ = 0;
        n_Total_mm_ = 0;

}

cuts_t DYEst::DilepSelect (int i_hyp)
{

     cuts_t ret = 0;
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

     // enough tracks?
     if (cms2.trks_trk_p4().size() > 2)
	  ret |= CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);
     // pt cuts
     if (cms2.hyp_lt_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LT_PT));
     if (cms2.hyp_ll_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LL_PT));
     // sign cuts
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) 
	  ret |= (CUT_BIT(CUT_OPP_SIGN));
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) 
	  ret |= (CUT_BIT(CUT_SAME_SIGN));
     // track corrected MET
     const TVector3 trkCorr(0, 0, 0); // well, don't correct until the tcmet in
			     // the 2_2 samples is figured out...
     // const TVector3 trkCorr = correctMETforTracks();
     if (pass4Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS4_TCMET));
     if (pass2Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS2_TCMET));
     // MET
     if (pass4Met(i_hyp, TVector3()))
	  ret |= (CUT_BIT(CUT_PASS4_MET));
     if (pass2Met(i_hyp, TVector3()))
	  ret |= (CUT_BIT(CUT_PASS2_MET));
     // muon quality
     int n_iso_mu = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
	  n_iso_mu++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
	  n_iso_mu++;
     }
     // electron quality
     int n_iso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_EL_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_EL_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_lt_index()[i_hyp], false)) {
	  ret |= CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_EL_ISO);
	  n_iso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_ll_index()[i_hyp], false)) {
	  ret |= CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_EL_ISO);
	  n_iso_el++;
     }     
     if (n_iso_mu + n_iso_el >= 1)
	  ret |= (CUT_BIT(CUT_ONE_ISO));
     if (n_iso_mu + n_iso_el >= 2)
	  ret |= (CUT_BIT(CUT_TWO_ISO));
     // electrons without d0
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_lt_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_ll_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
     // supertight cuts (only for electrons)
//      int n_supertight_el = 0;
//      if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
// 	  if (supertightElectron(cms2.hyp_lt_index()[i_hyp])) {
// // 	       ret |= (CUT_BIT(CUT_LT_SUPERTIGHT));
// 	       n_supertight_el++;
// 	  }
//      }
//      if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
// 	  if (supertightElectron(cms2.hyp_ll_index()[i_hyp])) {
// // 	       ret |= (CUT_BIT(CUT_LL_SUPERTIGHT));
// 	       n_supertight_el++;
// 	  }
//      }
//      if (n_supertight_el >= 1)
// 	  ret |= CUT_BIT(CUT_ONE_SUPERTIGHT);
//      if (n_supertight_el >= 2)
//    	  ret |= CUT_BIT(CUT_TWO_SUPERTIGHT);
     // supertight dphiin cut
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (deltaPhiInElectron(cms2.hyp_lt_index()[i_hyp])) {
	       ret |= CUT_BIT(CUT_LT_TIGHT_DPHIIN);
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (deltaPhiInElectron(cms2.hyp_ll_index()[i_hyp])) {
	       ret |= CUT_BIT(CUT_LL_TIGHT_DPHIIN);
	  }
     }
     // barrel
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (fabs(cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].eta()) < 1.479)
 	       ret |= (CUT_BIT(CUT_EL_BARREL));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (fabs(cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].eta()) < 1.479)
 	       ret |= (CUT_BIT(CUT_EL_BARREL));
     }
     // calo iso
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
     }
     int n_caloiso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
	  n_caloiso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
	  n_caloiso_el++;
     }     
     if (n_iso_mu + n_caloiso_el >= 1)
 	  ret |= (CUT_BIT(CUT_ONE_CALOISO));
     if (n_iso_mu + n_caloiso_el >= 2)
 	  ret |= (CUT_BIT(CUT_TWO_CALOISO));
     // jet veto
     if (cms2.hyp_njets()[i_hyp] == 0)
	  ret |= (CUT_BIT(CUT_PASS_JETVETO_CALO));
     // track jets
     if (passTrkJetVeto(i_hyp))
	  ret |= (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS));
     if (nJPTs(i_hyp, 20) == 0)
	  ret |= CUT_BIT(CUT_PASS_JETVETO_JPT20);
     if (nJPTs(i_hyp, 25) == 0)
	  ret |= CUT_BIT(CUT_PASS_JETVETO_JPT25);
     // muon b tag, with 20 GeV upper cut on the muon pt
     if (passMuonBVeto_1_6(i_hyp, true))
	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO));
     else ret |= (CUT_BIT(CUT_MUON_TAGGED));
     // muon b tag, with no upper cut on the muon pt
     if (numberOfExtraMuons(i_hyp, false) == 0)
	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT));
     else ret |= (CUT_BIT(CUT_MUON_TAGGED_WITHOUT_PTCUT));
     // Z veto
     if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2)
	  ret |= (CUT_BIT(CUT_PASS_ZVETO));
     else if (not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
	  ret |= (CUT_BIT(CUT_PASS_ZVETO));
     else ret |= (CUT_BIT(CUT_IN_Z_WINDOW));
     // Z veto using additional leptons in the event
     if (not additionalZveto())
	  ret |= (CUT_BIT(CUT_PASS_ADDZVETO));
     // any additional high-pt, isolated leptons?
     if (passTriLepVeto(i_hyp))
	  ret |= (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO));

//      if (myType == DILEPTON_MUMU) { // don't want to deal with electron overlap right now
// 	  if (passCaloTrkjetCombo())
// 	       ret |= CUT_BIT(CUT_PASS_JETVETO_CALOTRACKJETS_COMBO);
//      }

     //*****************************************************************
     // special handling for the fake rate cuts for now, because they
     // only work for emu
     //*****************************************************************
     if (myType != DILEPTON_EMU)
	  return ret;
     // in addition, for the muons, check that they pass tight+iso
     // (since the fake rate is electron only right now)
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  if ((ret & (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO))) != 
	      (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO)))
	       return ret;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  if ((ret & (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO))) != 
	      (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO)))
	       return ret;
     }

/*
     // now set the fake flags for the electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (isFakeable(cms2.hyp_lt_index()[i_hyp]))
	       ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
	  if (isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
	       ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
	  else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
     } else {
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	       if (isFakeable(cms2.hyp_ll_index()[i_hyp]))
		    ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
	       if (isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
		    ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
	       else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
	  }
     }
*/

	//
	// My stuff
	//

     if (pass5Met(i_hyp, trkCorr))
          ret |= (CUT_BIT(CUT_PASS5_MET));
     //if (metSimple(20.0, trkCorr))
     //     ret |= (CUT_BIT(CUT_MET_SIMPLE20));
     if (metSimple(35.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE35));
     if (metSimple(45.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE45));
     if (metBalance(i_hyp, trkCorr))
          ret |= (CUT_BIT(CUT_MET_BALLANCE));
     if (metProjected(i_hyp, trkCorr))
          ret |= (CUT_BIT(CUT_MET_PROJECTED));

     // Muon reco cleaning
     if (muonReconstructionCleaning(i_hyp, 0.10))
          ret |= (CUT_BIT(CUT_MUON_RECO_CLEANING));
     if (muonReconstructionCleaning(i_hyp, 0.20))
          ret |= (CUT_BIT(CUT_MUON_RECO_CLEANING20));


     return ret;

}

void DYEst::FillDilepHistos (int i_hyp)
{
     // everybody histogram needs to know what hypothesis he is 
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     // and what the event weight is 
     const double weight = Weight(i_hyp);

     // these are the cuts that the candidate passes:
     cuts_t cuts_passed = DilepSelect(i_hyp);

     // this is how to test that the candidate passes the cuts (which
     // we specified in the constructor when we made the looper)
     // (*note: the parentheses are important*):
     if ((cuts_passed & cuts_) == cuts_) {

	  // if the candidate passed, we count it
	  cands_passing[myType] += weight;
	  cands_passing_w2[myType] += weight * weight;
	  cands_count[myType]++;
	  cands_passing[DILEPTON_ALL] += weight;
	  cands_passing_w2[DILEPTON_ALL] += weight * weight;
	  cands_count[DILEPTON_ALL]++;
     }

     // for the NMinus1Hist, the histogram checks the cuts for us
     // now fill the N-1 histograms

     bool isInWindow = inZmassWindow(cms2.hyp_p4()[i_hyp].M());

	// get the tcMET
     const TVector3 trkCorr = correctMETforTracks();
     TVector3 hyp_met;
     hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
     hyp_met += trkCorr;

     // zero jet bin (including track jet veto)
     if (cms2.hyp_njets()[i_hyp] == 0 && nTrkJets(i_hyp) == 0) {
     	hnm1_mll_0j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

	//if (hyp_met.Pt() > 35.0) {
		//std::cout 	<< "hyp_type (" << myType << ") hyp_lt_mc_motherid, hyp_ll_mc_motherid: " 
		//		<< cms2.hyp_lt_mc_motherid()[i_hyp] << ", " << cms2.hyp_ll_mc_motherid()[i_hyp] << std::endl;

		//int eventCase = 0;
		//if ((cms2.hyp_lt_mc_motherid()[i_hyp] == 23 && abs(cms2.hyp_ll_mc_motherid()[i_hyp]) == 24)
		//     || (abs(cms2.hyp_lt_mc_motherid()[i_hyp]) == 24 && cms2.hyp_ll_mc_motherid()[i_hyp] == 23)) eventCase = 0;
		//else if (cms2.hyp_lt_mc_motherid()[i_hyp] == 23 && cms2.hyp_ll_mc_motherid()[i_hyp] == 23) eventCase = 1;

		//if (myType == 1) {
		//	if (eventCase == 0) n_WZ_mm_ ++;
		//	if (eventCase == 1) n_ZZ_mm_ ++;
		//}

                //if (myType == 3) {
                //        if (eventCase == 0) n_WZ_ee_ ++;
                //        if (eventCase == 1) n_ZZ_ee_ ++;
                //}

		//n_Total_mm_ ++;
		//n_Total_ee_ ++;

	//}

	if (isInWindow) hnm1_met_0j_in_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
	else hnm1_met_0j_out_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
     }

     // one jet bin
     if (cms2.hyp_njets()[i_hyp] == 1) {
        hnm1_mll_1j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

        if (isInWindow) hnm1_met_1j_in_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
        else hnm1_met_1j_out_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
     }

     // two or more jet bin
     if (cms2.hyp_njets()[i_hyp] >= 2) {
        hnm1_mll_2j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

        if (isInWindow) hnm1_met_2j_in_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
        else hnm1_met_2j_out_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
     }

}

void DYEst::End ()
{

	//std::cout << "mm" << std::endl;
	//std::cout << "n_WZ_mm_   " << n_WZ_mm_ << std::endl;
	//std::cout << "n_ZZ_mm_   " << n_ZZ_mm_ << std::endl;
	//std::cout << "n_Total_mm_" << n_Total_mm_ << std::endl;

        //std::cout << "ee" << std::endl;
        //std::cout << "n_WZ_ee_   " << n_WZ_ee_ << std::endl;
        //std::cout << "n_ZZ_ee_   " << n_ZZ_ee_ << std::endl;
        //std::cout << "n_Total_ee_" << n_Total_ee_ << std::endl;


     int ret = fprintf(logfile_, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample_.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");
}
