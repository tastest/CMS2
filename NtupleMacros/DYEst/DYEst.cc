#include <math.h>
#include <string>
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

     memset(cands_wz_all       , 0, sizeof(cands_wz_all       ));
     memset(cands_wz_zll       , 0, sizeof(cands_wz_zll       ));
     memset(cands_wz_all_w2    , 0, sizeof(cands_wz_all_w2    ));
     memset(cands_wz_zll_w2    , 0, sizeof(cands_wz_zll_w2    ));

}

void DYEst::BookHistos ()
{

	// for computing the R ratio / estimating the background
	// only remove the z veto in the zero jet bin
     hnm1_mll_0j_ = new NMinus1Hist(sample_, "mll_0j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)));
     hnm1_mll_1j_ = new NMinus1Hist(sample_, "mll_1j", 200, 0, 200, cuts_, jet_z_veto_cuts);
     hnm1_mll_2j_ = new NMinus1Hist(sample_, "mll_2j", 200, 0, 200, cuts_, jet_z_veto_cuts);

	// for studying the behavior of the R ratio as a function of MET
	// remove the straight met cut for this
	// note: remove the z veto from all because this is handled elsewhere
        // remove the jet veto from the 1 and 2 jet hists as this is handled elsewhere
     hnm1_met_in_0j_ = new NMinus1Hist(sample_, "met_in_0j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) | simple_met_cuts);
     hnm1_met_in_1j_ = new NMinus1Hist(sample_, "met_in_1j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) | (CUT_BIT(CUT_PASS_JETVETO_JPT20)) | simple_met_cuts);
     hnm1_met_in_2j_ = new NMinus1Hist(sample_, "met_in_2j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) | (CUT_BIT(CUT_PASS_JETVETO_JPT20)) | simple_met_cuts);
     hnm1_met_out_0j_ = new NMinus1Hist(sample_, "met_out_0j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) |simple_met_cuts);
     hnm1_met_out_1j_ = new NMinus1Hist(sample_, "met_out_1j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) | (CUT_BIT(CUT_PASS_JETVETO_JPT20)) | simple_met_cuts);
     hnm1_met_out_2j_ = new NMinus1Hist(sample_, "met_out_2j", 200, 0, 200, cuts_, (CUT_BIT(CUT_PASS_ZVETO)) | (CUT_BIT(CUT_PASS_JETVETO_JPT20)) | simple_met_cuts);

	std::string logFileString = std::string(fname_);
	std::string fileNameBase = logFileString.substr(0, logFileString.length() - 4);
	std::string outFileName = "LooperTree_" + fileNameBase + "_" + sample_.name + ".root";

//     std::string outFileName = "Looper_" + sample_.name + ".root";
     outFile_ = new TFile(outFileName.c_str(), "RECREATE");
     outFile_->cd();
     outTree_ = new TTree("T1", "Tree");

     // book the branches
     outTree_->Branch("sample_id", &sample_id_, "sample_id/I");
     outTree_->Branch("hyp_type", &hyp_type_, "hyp_type/I");
     outTree_->Branch("weight", &weight_, "weight/F");

     outTree_->Branch("evt_lumiblock", &evt_lumiblock_, "evt_lumiblock/I");
     outTree_->Branch("evt_run", &evt_run_, "evt_run/I");
     outTree_->Branch("evt_event", &evt_event_, "evt_event/I");

     outTree_->Branch("tcmet", &tcmet_, "tcmet/F");
     outTree_->Branch("phi_tcmet", &phi_tcmet_, "phi_tcmet/F");

     outTree_->Branch("mll", &mll_, "mll/F");
     outTree_->Branch("n_jptjets", &n_jptjets_, "n_jptjets/I");

     outTree_->Branch("pt_ll_trk", &pt_ll_trk_, "pt_ll_trk/F");
     outTree_->Branch("pt_ll_glb", &pt_ll_glb_, "pt_ll_glb/F");
     outTree_->Branch("phi_ll", &phi_ll_, "phi_ll/F");
     outTree_->Branch("pt_lt_trk", &pt_lt_trk_, "pt_lt_trk/F");
     outTree_->Branch("pt_lt_glb", &pt_lt_glb_, "pt_lt_glb/F");
     outTree_->Branch("phi_lt", &phi_lt_, "phi_lt/F");

}

cuts_t DYEst::DilepSelect (int i_hyp)
{

     cuts_t ret = 0;
     //const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

	// pass trigger?
	if (passTriggersMu9orLisoE15(cms2.hyp_type()[i_hyp]))
 		ret |= CUT_BIT(CUT_PASS_TRIGGER);

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
/*
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
*/

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
     if (metSimple(20.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE20));

     if (metSimple(25.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE25));

     if (metSimple(30.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE30));

     if (metSimple(35.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE35));
     if (metSimple(40.0, trkCorr))
          ret |= (CUT_BIT(CUT_MET_SIMPLE40));
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

          cands_wz_all[myType] += weight;
          cands_wz_all_w2[myType] += weight*weight;
          if (cms2.hyp_ll_mc_motherid()[i_hyp] == 23
              && cms2.hyp_lt_mc_motherid()[i_hyp] == 23) {
          	cands_wz_zll[myType] += weight;
		cands_wz_zll_w2[myType] += weight*weight;
          }
          
	  // if both leptons come from Z
	  if (myType == DILEPTON_MUMU) {
		std::cout << "MM: " << cms2.hyp_ll_mc_motherid()[i_hyp] << ", " << cms2.hyp_lt_mc_motherid()[i_hyp] << std::endl;
	  }
          // if both leptons come from Z
          if (myType == DILEPTON_EE) {
                std::cout << "EE: " << cms2.hyp_ll_mc_motherid()[i_hyp] << ", " << cms2.hyp_lt_mc_motherid()[i_hyp] << std::endl;
          }

     }


     if ((cuts_passed & baseline_cuts_nometsimple_nozveto) == baseline_cuts_nometsimple_nozveto) {

          // fill tree branches

	  evt_lumiblock_ = cms2.evt_lumiBlock();
	  evt_run_ = cms2.evt_run();
 	  evt_event_ = cms2.evt_event();

          sample_id_ = sample_.process;
          hyp_type_ = myType;
          weight_ = weight;
          n_jptjets_ = nJPTs(i_hyp);
          tcmet_ = cms2.evt_tcmet();
          phi_tcmet_ = cms2.evt_tcmetPhi();
          mll_ = cms2.hyp_p4()[i_hyp].M();

          if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) 
		pt_ll_glb_ = cms2.mus_trk_p4()[cms2.hyp_ll_index()[i_hyp]].Pt();
          else pt_ll_glb_ = cms2.hyp_ll_p4()[i_hyp].Pt();

          if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) 
		pt_lt_glb_ = cms2.mus_trk_p4()[cms2.hyp_lt_index()[i_hyp]].Pt();
          else pt_lt_glb_ = cms2.hyp_lt_p4()[i_hyp].Pt();
      
          pt_ll_trk_ = cms2.hyp_ll_p4()[i_hyp].Pt();
          phi_ll_ = cms2.hyp_ll_p4()[i_hyp].Phi();

          pt_lt_trk_ = cms2.hyp_lt_p4()[i_hyp].Pt();
          phi_lt_ = cms2.hyp_lt_p4()[i_hyp].Phi();


          // when everything is set
          // fill this entry in the tree
          outTree_->Fill();
     }

     // for the NMinus1Hist, the histogram checks the cuts for us
     // now fill the N-1 histograms

     bool isInWindow = inZmassWindow(cms2.hyp_p4()[i_hyp].M());

	// get the tcMET
     const TVector3 trkCorr(0, 0, 0); // = correctMETforTracks();
     TVector3 hyp_met;
     hyp_met.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
     hyp_met += trkCorr;

     // zero jet bin (including track jet veto)
////     if (cuts_passed & CUT_BIT(CUT_PASS_JETVETO_JPT20)) {
     	hnm1_mll_0j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

	if (isInWindow) hnm1_met_in_0j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
	else hnm1_met_out_0j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
////     }

     // one jet bin
     if (cms2.hyp_njets()[i_hyp] == 1) {
        hnm1_mll_1j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

        if (isInWindow) hnm1_met_in_1j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
        else hnm1_met_out_1j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
     }

     // two or more jet bin
     if (cms2.hyp_njets()[i_hyp] >= 2) {
        hnm1_mll_2j_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].M(), weight);

        if (isInWindow) hnm1_met_in_2j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
        else hnm1_met_out_2j_->Fill(cuts_passed, myType, hyp_met.Pt(), weight);
     }

}

void DYEst::End ()
{

	// tidy up tree stuff
        outFile_->cd();
        outTree_->Write();
        outFile_->Close();
        delete outFile_; 

     // mm
     std::cout << "... mumu " << std::endl;
     std::cout << "...... all " << cands_wz_all[DILEPTON_MUMU] << " $\\pm$ " << sqrt(cands_wz_all_w2[DILEPTON_MUMU]) << std::endl;
     std::cout << "...... zll " << cands_wz_zll[DILEPTON_MUMU] << " $\\pm$ " << sqrt(cands_wz_zll_w2[DILEPTON_MUMU]) << std::endl;

     // ee
     std::cout << "... ee " << std::endl;
     std::cout << "...... all " << cands_wz_all[DILEPTON_EE] << " $\\pm$ " << sqrt(cands_wz_all_w2[DILEPTON_EE]) << std::endl;
     std::cout << "...... zll " << cands_wz_zll[DILEPTON_EE] << " $\\pm$ " << sqrt(cands_wz_zll_w2[DILEPTON_EE]) << std::endl;


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
