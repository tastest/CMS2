#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

static const double d0_bins[] = { 0, 0.025, 0.03, 0.055 };
static const double dphiin_bins[] = { -0.04, 0, 0.04, 0.045, 0.085 };
static const double iso_bins[] = { 0, 0.82, 0.9, 0.92, 1.0001};

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_	, 0, sizeof(cands_passing_       ));
     memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
     memset(cands_count_		, 0, sizeof(cands_count_         ));
     memset(count_cuts_		, 0, sizeof(count_cuts_         ));
     memset(count_correlation_	, 0, sizeof(count_correlation_         ));
}

void Looper::BookHistos ()
{
     hnJet		= new NMinus1Hist(sample_, "nJet"            ,	 6	, -0.5, 5.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
     hnCaloJet		= new NMinus1Hist(sample_, "nCaloJet"        ,	 6	, -0.5, 5.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS))	);
     hnTrackJet		= new NMinus1Hist(sample_, "nTrackJet"       ,	 6	, -0.5, 5.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
     hnJPTJet		= new NMinus1Hist(sample_, "nJPTJet"       ,	 6	, -0.5, 5.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS) | (CUT_BIT(CUT_PASS_JETVETO_JPT20))));
     hcaloJetPt		= new NMinus1Hist(sample_, "caloJetPt"       ,	 150	, 0, 150	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) 	);
     htrackJetPt	= new NMinus1Hist(sample_, "trackJetPt"      ,	 150	, 0, 150	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
     hJPTJetPt		= new NMinus1Hist(sample_, "JPTJetPt"      ,	 30	, 0, 150	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	| (CUT_BIT(CUT_PASS_JETVETO_JPT20)));
       hminLepPt	= new NMinus1Hist(sample_, "minLepPt"        ,	 150	, 0, 150	, cuts_, CUT_BIT(CUT_LL_PT)	);
       hmaxLepPt	= new NMinus1Hist(sample_, "maxLepPt"        ,	 150	, 0, 150	, cuts_, CUT_BIT(CUT_LL_PT)  	);
       hltPt		= new NMinus1Hist(sample_, "ltPt"            ,	 15	, 0, 150	, cuts_, (CUT_BIT(CUT_LT_PT))	);
       hllPt		= new NMinus1Hist(sample_, "llPt"            ,	 15	, 0, 150	, cuts_, (CUT_BIT(CUT_LL_PT))	);
       helPt		= new NMinus1Hist(sample_, "elPt"            ,	 16	, 0, 160	, cuts_, CUT_BIT(CUT_LL_PT)	);
       hmuPt		= new NMinus1Hist(sample_, "muPt"            ,	 16	, 0, 160	, cuts_, CUT_BIT(CUT_LL_PT)	);
       helEta		= new NMinus1Hist(sample_, "elEta"           ,	 12	, -3, 3		, cuts_, 0);
       hmuEta		= new NMinus1Hist(sample_, "muEta"           ,	 12	, -3, 3		, cuts_, 0);
       hdphiLep		= new NMinus1Hist(sample_, "dphiLep"         ,	 50	, 0, M_PI	, cuts_, 0	);
       hdilMass		= new NMinus1Hist(sample_, "dilMass"         ,	 30	, 0, 300	, cuts_, CUT_BIT(CUT_PASS_ZVETO) | CUT_BIT(CUT_PASS_ADDZVETO) | CUT_BIT(CUT_IN_Z_WINDOW));
       hdilPt		= new NMinus1Hist(sample_, "dilPt"           ,	 100	, 0, 300	, cuts_, 0	);
       hmet		= new NMinus1Hist(sample_, "met"             ,	 20	, 0, 100	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)	);
       hmetSpec		= new NMinus1Hist(sample_, "metSpec"         ,	 20	, 0, 100	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)  );
       hmetTrkCorr	= new NMinus1Hist(sample_, "metTrkCorr"      ,	 100	, 0, 200	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)	);
       hptJet1		= new NMinus1Hist(sample_, "ptJet1"          ,	 100	, 0, 300	, cuts_, 0		);
       hptJet2		= new NMinus1Hist(sample_, "ptJet2"          ,	 100	, 0, 300	, cuts_, 0       	);
       hptJet3		= new NMinus1Hist(sample_, "ptJet3"          ,	 100	, 0, 300	, cuts_, 0       	);
       hptJet4		= new NMinus1Hist(sample_, "ptJet4"          ,	 100	, 0, 300	, cuts_, 0       	);
       hetaJet1		= new NMinus1Hist(sample_, "etaJet1"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet2		= new NMinus1Hist(sample_, "etaJet2"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet3		= new NMinus1Hist(sample_, "etaJet3"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet4		= new NMinus1Hist(sample_, "etaJet4"         ,	 50	, -4, 4		, cuts_, 0       	);
       hnumTightLep	= new NMinus1Hist(sample_, "numTightLep"     ,	 6	, -0.5, 5.5	, cuts_, 0             	);
       heleRelIso	= new NMinus1Hist(sample_, "eleRelIso"       ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO)) | (CUT_BIT(CUT_LT_CALOISO)) | (CUT_BIT(CUT_LL_CALOISO)));
       heleRelIsoTrk	= new NMinus1Hist(sample_, "eleRelIsoTrk"    ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO)) | (CUT_BIT(CUT_LT_CALOISO)) | (CUT_BIT(CUT_LL_CALOISO)));
       hmuRelIso	= new NMinus1Hist(sample_, "muRelIso"        ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO)) | (CUT_BIT(CUT_LT_CALOISO)) | (CUT_BIT(CUT_LL_CALOISO))	);
       hminRelIso	= new NMinus1Hist(sample_, "minRelIso"       ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO))	);
       hminRelIso_withCalo = new NMinus1Hist	(sample_, "minRelIso_withCalo", 101, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO))	);
       htagMuPt		= new NMinus1Hist(sample_, "tagMuPt"	      ,	 100	, 0, 100	, cuts_ & ~((CUT_BIT(CUT_PASS_MUON_B_VETO)) | (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO))), (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)));
       htagMuRelIso	= new NMinus1Hist(sample_, "tagMuRelIso"     ,	 101	, 0, 1.01	, cuts_ & ~((CUT_BIT(CUT_PASS_MUON_B_VETO)) | (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO))), (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)));
       hmuPdgId		= new NMinus1Hist(sample_, "muPdgId", 		2301, -0.5, 2300.5, cuts_, 0);
       hmuMoPdgId	= new NMinus1Hist(sample_, "muMoPdgId", 	2301, -0.5, 2300.5, cuts_, 0);
       helPdgId		= new NMinus1Hist(sample_, "elPdgId", 		2301, -0.5, 2300.5, cuts_, 0);
       helMoPdgId	= new NMinus1Hist(sample_, "elMoPdgId", 	2301, -0.5, 2300.5, cuts_, 0);
       helEop      	= new NMinus1Hist(sample_, "elEop"	      ,  10	, 0, 10	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       held0    	= new NMinus1Hist(sample_, "eld0"	      ,  50	, 0, 0.1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helfbrem    	= new NMinus1Hist(sample_, "elfbrem"	      ,  11	, -0.1, 1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helHE       	= new NMinus1Hist(sample_, "elHE"	      ,  11	, -0.03, 0.3	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsee      	= new NMinus1Hist(sample_, "elsee"	      ,  50	, 0, 0.05	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsppEB  	= new NMinus1Hist(sample_, "elsppEB"	      ,  50	, 0, 0.0	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsppEE	  	= new NMinus1Hist(sample_, "elsppEC"	      ,  50	, 0, 0.0	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       heldphiin   	= new NMinus1Hist(sample_, "eldphiin"	      ,  10	, -0.1, 0.1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       heldetain   	= new NMinus1Hist(sample_, "eldetain"	      ,  10	, -0.02, 0.02	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helEseedopin	= new NMinus1Hist(sample_, "elEseedopin"     ,  10	, 0, 20	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helConvDeltaPhi_ss = new NMinus1Hist   	(sample_, "elConvDeltaPhi_ss"	    ,  10	, 0, 0.1, cuts_, 0);
       helConvDeltaPhi_os = new NMinus1Hist	(sample_, "elConvDeltaPhi_os"     ,  10	, 0, 0.1, cuts_, 0);
       held0vsRelIso 		= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "d0vsRelIso"), ";d0;rel iso", sizeof(d0_bins) / sizeof(double) - 1, d0_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       heldphiinvsRelIso 	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "dphiinvsRelIso"), ";dphiin;rel iso", sizeof(dphiin_bins) / sizeof(double) - 1, dphiin_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       held0vsRelIsoMCgamma 	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "d0vsRelIsoMCgamma"), ";d0;rel iso", sizeof(d0_bins) / sizeof(double) - 1, d0_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       heldphiinvsRelIsoMCgamma = new TH2F(Form("%s_%s_em", sample_.name.c_str(), "dphiinvsRelIsoMCgamma"), ";dphiin;rel iso", sizeof(dphiin_bins) / sizeof(double) - 1, dphiin_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       htrkCalodRvsPtSum	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "trkCalodRvsPtSum"), ";pt sum;#DeltaR", 10, 0, 30, 60, 0, M_PI);
       hCaloEtaPt		= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "CaloEtaPt"), ";pt;#eta", 10, 0, 30, 10, -5, 5);
}

bool Looper::FilterEvent()
{ 
    //
     // duplicate filter, based on trk information and dilepton hyp
     //
     if (cms2.trks_d0().size() == 0)
	  return true;
     DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
				 cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
     return is_duplicate(id); 
}

cuts_t Looper::EventSelect ()
{
     //------------------------------------------------------------
     // In an event-based analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::DilepSelect (int i_hyp)
{
     cuts_t ret = 0;
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

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
     const TVector3 trkCorr; // well, don't correct until the tcmet in
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
	  ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO)) | CUT_BIT(CUT_LT_CALOISO_1_6);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO)) | CUT_BIT(CUT_LL_CALOISO_1_6);
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
     // 1_6 calo iso
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation_1_6(cms2.hyp_lt_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LT_CALOISO_1_6);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation_1_6(cms2.hyp_ll_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LL_CALOISO_1_6);
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
     // this logic is here for historical curiosity (fake rate for emu
     // that assumes the mu is real and the e is fakeable) (if you
     // switch this logic on, bad things will probably happen in your
     // looper)
     // *****************************************************************
#if 0
     if (myType == DILEPTON_EMU) {
	  // in addition, for the muons, check that they pass tight+iso
	  // (since the fake rate is electron only right now)
#if 1
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
#else
	  // require MC truth mu from W instead
	  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && !trueMuonFromW_WJets(cms2.hyp_lt_index()[i_hyp]) )
	       return ret;
	  
	  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && !trueMuonFromW_WJets(cms2.hyp_ll_index()[i_hyp]) )
	       return ret;
#endif
	  
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
     }
#endif
     return ret;
}

double Looper::Weight (int)
{
     return cms2.evt_scale1fb() * 0.1;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
     //------------------------------------------------------------
     // In a trilepton analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp)
{
     //------------------------------------------------------------
     // In a quadlepton analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

void Looper::FillEventHistos ()
{
     //------------------------------------------------------------
     // In an event-based analysis, you would fill your histos here
     //------------------------------------------------------------
}

void Looper::CountCuts (cuts_t cuts_passed, double weight) 
{
     for (int i = 0; i < 64; ++i) {
	  if (cuts_passed & CUT_BIT(i)) {
	       count_cuts_[i] += weight;
	       for (int j = 0; j < 64; ++j) {
		    if (cuts_passed & CUT_BIT(j)) 
			 count_correlation_[i][j] += weight;
	       }
	  }
     }
}

void Looper::FillDilepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // Example dilepton histo filling; edit for your application
     //------------------------------------------------------------

//      if (cms2.vtxs_position().size() != 0)
// 	  printf("size %d\n", cms2.vtxs_position().size());

     // every histogram needs to know what hypothesis he is 
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     // and what the event weight is 
     const double weight = Weight(i_hyp);
     
     // these are the cuts that the candidate passes:
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     CountCuts(cuts_passed, weight);
     if ((cuts_passed & cuts_) == cuts_) {
	  cands_passing_[myType] += weight;
	  cands_passing_w2_[myType] += weight * weight;
	  cands_count_[myType]++;
	  cands_passing_[DILEPTON_ALL] += weight;
	  cands_passing_w2_[DILEPTON_ALL] += weight * weight;
	  cands_count_[DILEPTON_ALL]++;
// 	  printf("candidate passing: run %10lu, event %10lu, lumi block %10lu (lt pt = %5.1f, ll pt = %5.1f)\n",
// 		 cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(),
// 		 cms2.hyp_lt_p4()[i_hyp].pt(), cms2.hyp_ll_p4()[i_hyp].pt());
// 	  for (int i = 0; i < cms2.els_p4().size(); ++i) {
// 	       printf("E = (%5.1f %6.3f %6.3f) trk p = (%5.1f %6.3f %6.3f) eSC = %5.1f E/p = %5.3f\t", 
// 		      cms2.els_p4()[i].Et(), cms2.els_p4()[i].eta(), cms2.els_p4()[i].phi(),
// 		      cms2.els_trk_p4()[i].pt(), cms2.els_trk_p4()[i].eta(), cms2.els_trk_p4()[i].phi(),
// 		      cms2.els_eSC()[i_hyp],
// 		      cms2.els_eOverPIn()[i]);
// 	  }
// 	  printf("\n");
     }

     // jet count
     hnJet->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     hnCaloJet	->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     hnTrackJet	->Fill(cuts_passed, myType, nTrkJets(i_hyp), weight);
     hnJPTJet	->Fill(cuts_passed, myType, nJPTs(i_hyp, 20), weight);

     // jet pt's (of the max jet) -- so much work we only do it for the ones we care about
     std::vector<LorentzVector> allJpts = JPTs(i_hyp, 0);
     double maxPtJPT = -1;
     for (unsigned int i = 0; i < allJpts.size(); ++i) {
	  if (allJpts[i].pt() > maxPtJPT)
	       maxPtJPT = allJpts[i].pt();
     }
     hJPTJetPt->Fill(cuts_passed, myType, maxPtJPT, weight);

     // lepton pt's
     hminLepPt->Fill(cuts_passed, myType, std::min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt->Fill(cuts_passed, myType, std::max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hltPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     hllPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  helPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  helEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  hmuEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  helPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  helEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  hmuEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     }
    
     // dilepton mass
     hdilMass->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
    
     // delta phi btw leptons
     double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
     if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
     hdphiLep->Fill(cuts_passed, myType, dphi, weight);
    
     // Relative isolation... muons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  hmuRelIso->Fill(cuts_passed, myType, reliso_lt(i_hyp), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  hmuRelIso->Fill(cuts_passed, myType, reliso_ll(i_hyp), weight);
     }

     // Relative isolation... electrons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  heleRelIso->Fill(cuts_passed, myType, reliso_lt(i_hyp, true), weight);
	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_lt(i_hyp, false), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  heleRelIso->Fill(cuts_passed, myType, reliso_ll(i_hyp, true), weight);
	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_ll(i_hyp, false), weight);
     }

     // lower of the two isolations, regardless of species (used for Dumbo method)
     hminRelIso->Fill(cuts_passed, myType, std::min(reliso_lt(i_hyp), reliso_ll(i_hyp)), weight);
     // lower of the two isolations, regardless of species (used for
     // Dumbo method) --- using calo iso for electrons
     hminRelIso_withCalo->Fill(cuts_passed, myType, 
			      std::min(reliso_lt(i_hyp, true), reliso_ll(i_hyp, true)), 
			      weight);

     // dilepton pt
     hdilPt->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met and Met special
     hmet->Fill(cuts_passed, myType, cms2.hyp_met()[i_hyp], weight);      
     hmetSpec->Fill(cuts_passed, myType, 
		   MetSpecial(cms2.hyp_met()[i_hyp], cms2.hyp_metPhi()[i_hyp], i_hyp),
		   weight);
     hmetTrkCorr->Fill(cuts_passed, myType, cms2.evt_tcmet(), weight);

     // tag muon pt and iso
     htagMuPt->Fill(cuts_passed, myType, tagMuonPt(i_hyp), weight);
     htagMuRelIso->Fill(cuts_passed, myType, tagMuonRelIso(i_hyp), weight);

     if (myType == DILEPTON_EMU) {
	  int mu_idx = -1;
	  int el_idx = -1;
	  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	       mu_idx = cms2.hyp_lt_index()[i_hyp];
	       el_idx = cms2.hyp_ll_index()[i_hyp];
	  }
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	       mu_idx = cms2.hyp_ll_index()[i_hyp];
	       el_idx = cms2.hyp_lt_index()[i_hyp];
	  }
	  assert(el_idx != -1 && mu_idx != -1);
	  hmuPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_id()[mu_idx]), weight);
	  hmuMoPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_motherid()[mu_idx]), weight);
	  helPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_id()[el_idx]), weight);
	  helMoPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_motherid()[el_idx]), weight);
	  helEop      ->Fill(cuts_passed, myType, cms2.els_eOverPIn	()[el_idx], weight);
	  held0	      ->Fill(cuts_passed, myType, fabs(cms2.els_d0corr		()[el_idx]), weight);
	  if ((cuts_passed & cuts_) == cuts_) {
	       held0vsRelIso->Fill(fabs(cms2.els_d0corr		()[el_idx]), el_rel_iso(el_idx, true), weight);
	       heldphiinvsRelIso->Fill(cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], el_rel_iso(el_idx, true), weight);
	       if (abs(cms2.els_mc_id()[el_idx]) == 22) {
		    held0vsRelIsoMCgamma->Fill(fabs(cms2.els_d0corr		()[el_idx]), el_rel_iso(el_idx, true), weight);
		    heldphiinvsRelIsoMCgamma->Fill(cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], el_rel_iso(el_idx, true), weight);
	       }
	  }
	  helfbrem    ->Fill(cuts_passed, myType, cms2.els_fBrem		()[el_idx], weight);
	  helHE       ->Fill(cuts_passed, myType, cms2.els_hOverE	()[el_idx], weight);
	  helsee      ->Fill(cuts_passed, myType, cms2.els_sigmaEtaEta	()[el_idx], weight);
	  if (cuts_passed & CUT_BIT(CUT_EL_BARREL))
	       helsppEB      ->Fill(cuts_passed, myType, cms2.els_sigmaPhiPhi	()[el_idx], weight);
	  else helsppEE      ->Fill(cuts_passed, myType, cms2.els_sigmaPhiPhi	()[el_idx], weight);
	  heldphiin   ->Fill(cuts_passed, myType, cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], weight);
// 	  heldetain   ->Fill(cuts_passed, myType, cms2.els_dEtaIn	()[el_idx], weight);
	  helEseedopin->Fill(cuts_passed, myType, cms2.els_eSeedOverPOut	()[el_idx], weight);
// 	  if (cms2.els_mc_id()[el_idx] == 22 && (cuts_passed & cuts) == cuts)
// 	       printf("run %10u, event %10u: weight %f\n", cms2.evt_run(), cms2.evt_event(),
// 		      cms2.evt_scale1fb());
     }

     // coversions
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  const int conv = conversionPartner(cms2.hyp_lt_index()[i_hyp]);
	  if (conv != -1) {
	       const double dphi = conversionDeltaPhi(conv, cms2.hyp_lt_index()[i_hyp]);
	       if (cms2.hyp_lt_id()[i_hyp] * cms2.trks_charge()[conv] > 0)
		    helConvDeltaPhi_os->Fill(cuts_passed, myType, dphi, weight);
	       else helConvDeltaPhi_ss->Fill(cuts_passed, myType, dphi, weight);
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  const int conv = conversionPartner(cms2.hyp_ll_index()[i_hyp]);
	  if (conv != -1) {
	       const double dphi = conversionDeltaPhi(conv, cms2.hyp_ll_index()[i_hyp]);
	       if (cms2.hyp_ll_id()[i_hyp] * cms2.trks_charge()[conv] > 0)
		    helConvDeltaPhi_os->Fill(cuts_passed, myType, dphi, weight);
	       else helConvDeltaPhi_ss->Fill(cuts_passed, myType, dphi, weight);
	  }
     }
}

void Looper::FillTrilepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // In a trilepton analysis, you would fill your histos here
     //------------------------------------------------------------
}

void Looper::FillQuadlepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // In a quadlepton analysis, you would fill your histos here
     //------------------------------------------------------------
}

void Looper::End ()
{
     //------------------------------------------------------------
     //Example status message at the end of a looper; edit for your
     //application
     //------------------------------------------------------------

     int ret = fprintf(logfile_, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.3f +- %10.3f %10.3f +- %10.3f %10.3f +- %10.3f %10.3f +- %10.3f\n",   
		       sample_.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
     if (ret < 0)
	  perror("writing to log file");
     for (int i = 0; i < 64; ++i) {
	  ret = fprintf(logfile_, 
			"Sample %10s: cands passing cut %c%2d: %10.1f\n",
			sample_.name.c_str(), cuts_ & CUT_BIT(i) ? '*' : ' ',
			i, count_cuts_[i]);
	  if (ret < 0)
	       perror("writing to log file");
     }
     ret = fprintf(logfile_, 
		   "Sample %10s: correlations\nSample %10s:          ", 
		   sample_.name.c_str(), sample_.name.c_str());
     if (ret < 0)
	  perror("writing to log file");
     for (int i = 0; i < 64; ++i) {
	  if (not (cuts_ & CUT_BIT(i)))
	       continue;
	  ret = fprintf(logfile_, " %c%2d  ", cuts_ & CUT_BIT(i) ? '*' : ' ', i);
	  if (ret < 0)
	       perror("writing to log file");
     }
     ret = fprintf(logfile_, "\n");
     if (ret < 0)
	  perror("writing to log file");
     for (int i = 0; i < 64; ++i) {
	  if (not (cuts_ & CUT_BIT(i)))
	       continue;
	  ret = fprintf(logfile_, 
			"Sample %10s: %c%2d     ", 
			sample_.name.c_str(), cuts_ & CUT_BIT(i) ? '*' : ' ', i);
	  if (ret < 0)
	       perror("writing to log file");
	  const double den = count_cuts_[i];
	  for (int j = 0; j < 64; ++j) {
	       if (not (cuts_ & CUT_BIT(j)))
		    continue;
	       if (den != 0) {
		    ret = fprintf(logfile_, 
				  "%5.1f ", count_correlation_[i][j] / den);
	       } else {
		    ret = fprintf(logfile_, 
				  "----- ");
	       }
	       if (ret < 0)
		    perror("writing to log file");
	  }
	  ret = fprintf(logfile_, "\n");
	  if (ret < 0)
	       perror("writing to log file");
     }
}
