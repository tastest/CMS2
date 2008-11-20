#include <math.h>
#include "TVector3.h"
#include "selections.h"
#include "utilities.h"
#include "CMS2.h"
#include "tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters
     memset(cands_passing	, 0, sizeof(cands_passing       ));
     memset(cands_passing_w2	, 0, sizeof(cands_passing_w2    ));
     memset(cands_count		, 0, sizeof(cands_count         ));
}

void Looper::BookHistos ()
{
     // book histograms the manual way:
     for (int i = 0; i < 4; ++i) {
	  helPt[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "elPt", dilepton_hypo_names[i]), ";el pt", 100, 0, 100);
	  hmuPt[i] = new TH1F(Form("%s_%s_%s", SampleName().c_str(), "muPt", dilepton_hypo_names[i]), ";mu pt", 100, 0, 100);
	  hCaloEtaPt[i] = new TH2F(Form("%s_%s_%s", SampleName().c_str(), "CaloEtaPt", dilepton_hypo_names[i]), ";pt;eta", 100, 0, 100, 10, -5, 5);
     }
     // or use the N - 1 technology (see NMinus1Hist.h)
     // arguments are as follows: sample, name, binning, required cuts, cuts that are relaxed for the N - 1 plot
     // for the lt N - 1 plot, we relax the lt pt requirement
     hltPt		= new NMinus1Hist(sample, "ltPt"   ,	 150, 0, 150, cuts, CUT_BIT(CUT_LT_PT));
     // same for ll
     hllPt		= new NMinus1Hist(sample, "llPt"   ,	 150, 0, 150, cuts, CUT_BIT(CUT_LL_PT));
     // for the dilepton mass plot, we relax any cut to do with the Z 
     hdilMass		= new NMinus1Hist(sample, "dilMass",	 100, 0, 300, cuts, 
					  CUT_BIT(CUT_PASS_ZVETO) | CUT_BIT(CUT_PASS_ADDZVETO) | CUT_BIT(CUT_IN_Z_WINDOW));
}

cuts_t Looper::DilepSelect (int i_hyp)
{
     // cuts are failed until proven otherwise
     cuts_t ret = 0;

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
     const TVector3 trkCorr = correctMETforTracks();
     if (pass4Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS4_METCORR));
     if (pass2Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS2_METCORR));
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
     int n_supertight_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (supertightElectron(cms2.hyp_lt_index()[i_hyp])) {
	       n_supertight_el++;
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (supertightElectron(cms2.hyp_ll_index()[i_hyp])) {
	       n_supertight_el++;
	  }
     }
     if (n_supertight_el >= 1)
	  ret |= CUT_BIT(CUT_ONE_SUPERTIGHT);
     if (n_supertight_el >= 2)
   	  ret |= CUT_BIT(CUT_TWO_SUPERTIGHT);
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
     // electron barrel
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
     // muon b tag, with 20 GeV upper cut on the muon pt
     if (passMuonBVeto(i_hyp, true))
	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO));
     else ret |= (CUT_BIT(CUT_MUON_TAGGED));
     // muon b tag, with no upper cut on the muon pt
     if (passMuonBVeto(i_hyp, false))
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

     // the return value gets cached, too
     return ret;
}

void Looper::FillDilepHistos (int i_hyp)
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
     if ((cuts_passed & cuts) == cuts) {

	  // if the candidate passed, we count it
	  cands_passing[myType] += weight;
	  cands_passing_w2[myType] += weight * weight;
	  cands_count[myType]++;
	  cands_passing[DILEPTON_ALL] += weight;
	  cands_passing_w2[DILEPTON_ALL] += weight * weight;
	  cands_count[DILEPTON_ALL]++;
     }

     // for TH1/TH2, we have to check explicitly whether the candidate passes
     if ((cuts_passed & cuts) == cuts) {
	  // and then fill
	  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	       helPt[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  } else {
	       hmuPt[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  }
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	       helPt[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  } else {
	       hmuPt[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  }
	  for (unsigned int i = 0; i < cms2.jets_p4().size(); ++i) {
	       // histogram the eta and pt of all jets
	       hCaloEtaPt[myType]->Fill(cms2.jets_p4()[i].pt() * cms2.jets_tq_noCorrF()[i],
					cms2.jets_p4()[i].eta(),
					weight);
	  }
     }

     // for the NMinus1Hist, the histogram checks the cuts for us
     hltPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     hllPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     hdilMass->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
}

void Looper::End ()
{
     int ret = fprintf(logfile, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");
}
