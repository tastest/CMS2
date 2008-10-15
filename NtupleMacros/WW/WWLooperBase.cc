#include <math.h>
#include "TChain.h"
#include "TVector3.h"
#include "selections.h"
#include "utilities.h"
#include "fakerates.h"
#include "CMS2.h"
#include "WWLooper.h"

WWLooperBase::WWLooperBase (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s),
       cuts(c),
       hnJet		(s, "nJet"            ,	 6	, -0.5, 5	, cuts, (1 << WW_PASS_JETVETO_CALO) | (1 << WW_PASS_JETVETO_TRACKJETS) 	),
       hnCaloJet	(s, "nCaloJet"        ,	 6	, -0.5, 5	, cuts, (1 << WW_PASS_JETVETO_CALO) | (1 << WW_PASS_JETVETO_TRACKJETS)	),
       hnTrackJet	(s, "nTrackJet"       ,	 6	, -0.5, 5	, cuts, (1 << WW_PASS_JETVETO_CALO) | (1 << WW_PASS_JETVETO_TRACKJETS) 	),
       hcaloJetPt	(s, "caloJetPt"       ,	 6	, -0.5, 5	, cuts, (1 << WW_PASS_JETVETO_CALO) 	),
       htrackJetPt	(s, "trackJetPt"      ,	 6	, -0.5, 5	, cuts, (1 << WW_PASS_JETVETO_TRACKJETS) 	),
       hminLepPt	(s, "minLepPt"        ,	 150	, 0, 150	, cuts, 1 << WW_LL_PT	),
       hmaxLepPt	(s, "maxLepPt"        ,	 150	, 0, 150	, cuts, 1 << WW_LL_PT  	),
       hltPt		(s, "ltPt"            ,	 150	, 0, 150	, cuts, (1 << WW_LT_PT)	),
       hllPt		(s, "llPt"            ,	 150	, 0, 150	, cuts, (1 << WW_LL_PT)	),
       helPt		(s, "elPt"            ,	 16	, 0, 160	, cuts, 1 << WW_LL_PT	),
       hmuPt		(s, "muPt"            ,	 16	, 0, 160	, cuts, 1 << WW_LL_PT	),
       helEta		(s, "elEta"           ,	 12	, -3, 3		, cuts, 0),
       hmuEta		(s, "muEta"           ,	 12	, -3, 3		, cuts, 0),
       helHE		(s, "elHE"	      ,  11	, -0.1, 1	, cuts, (1 << WW_LT_GOOD | 1 << WW_LL_GOOD | 1 << WW_LT_ISO | 1 << WW_LL_ISO | 1 << WW_LT_CALOISO | 1 << WW_LL_CALOISO)),
       hdphiLep		(s, "dphiLep"         ,	 50	, 0, M_PI	, cuts, 0	),
       hdilMass		(s, "dilMass"         ,	 100	, 0, 300	, cuts, (1 << WW_PASS_ZVETO) | (1 << WW_PASS_ADDZVETO) 	),
       hdilPt		(s, "dilPt"           ,	 100	, 0, 300	, cuts, 0	),
       hmet		(s, "met"             ,	 100	, 0, 200	, cuts, (1 << WW_PASS4_MET) | (1 << WW_PASS2_MET)	),
       hmetSpec		(s, "metSpec"         ,	 100	, 0, 200	, cuts, (1 << WW_PASS4_MET) | (1 << WW_PASS2_MET)      	),
       hmetTrkCorr	(s, "metTrkCorr"      ,	 100	, 0, 200	, cuts, (1 << WW_PASS4_MET) | (1 << WW_PASS2_MET)	),
       hptJet1		(s, "ptJet1"          ,	 100	, 0, 300	, cuts, 0		),
       hptJet2		(s, "ptJet2"          ,	 100	, 0, 300	, cuts, 0       	),
       hptJet3		(s, "ptJet3"          ,	 100	, 0, 300	, cuts, 0       	),
       hptJet4		(s, "ptJet4"          ,	 100	, 0, 300	, cuts, 0       	),
       hetaJet1		(s, "etaJet1"         ,	 50	, -4, 4		, cuts, 0       	),
       hetaJet2		(s, "etaJet2"         ,	 50	, -4, 4		, cuts, 0       	),
       hetaJet3		(s, "etaJet3"         ,	 50	, -4, 4		, cuts, 0       	),
       hetaJet4		(s, "etaJet4"         ,	 50	, -4, 4		, cuts, 0       	),
       hnumTightLep	(s, "numTightLep"     ,	 6	, -0.5, 5.5	, cuts, 0             	),
       heleRelIso	(s, "eleRelIso"       ,	 101	, 0, 1.01	, cuts, (1 << WW_LT_ISO) | (1 << WW_LL_ISO)	),
       hmuRelIso	(s, "muRelIso"        ,	 101	, 0, 1.01	, cuts, (1 << WW_LT_ISO) | (1 << WW_LL_ISO)	),
       hminRelIso	(s, "minRelIso"       ,	 101	, 0, 1.01	, cuts, (1 << WW_LT_ISO) | (1 << WW_LL_ISO)	),
       hminRelIso_withCalo	(s, "minRelIso_withCalo", 101, 0, 1.01	, cuts, (1 << WW_LT_ISO) | (1 << WW_LL_ISO)	),
       htagMuPt		(s, "tagMuPt"	      ,	 100	, 0, 100	, cuts & ~((1 << WW_PASS_MUON_B_VETO) | (1 << WW_PASS_MUON_B_VETO_WITHOUT_PTCUT) | (1 << WW_PASS_EXTRALEPTON_VETO)), (1 << WW_PASS_JETVETO_TRACKJETS)),
       htagMuRelIso	(s, "tagMuRelIso"     ,	 101	, 0, 1.01	, cuts & ~((1 << WW_PASS_MUON_B_VETO) | (1 << WW_PASS_MUON_B_VETO_WITHOUT_PTCUT) | (1 << WW_PASS_EXTRALEPTON_VETO)), (1 << WW_PASS_JETVETO_TRACKJETS)),
       hmuPdgId		(s, "muPdgId", 		2301, -0.5, 2300.5),
       hmuMoPdgId	(s, "muMoPdgId", 	2301, -0.5, 2300.5),
       helPdgId		(s, "elPdgId", 		2301, -0.5, 2300.5),
       helMoPdgId	(s, "elMoPdgId", 	2301, -0.5, 2300.5)
{
     memset(cands_passing	, 0, sizeof(cands_passing       ));
     memset(cands_passing_w2	, 0, sizeof(cands_passing_w2    ));
     memset(cands_count		, 0, sizeof(cands_count         ));
     if (fname != 0 && strlen(fname) != 0) {
	  logfile = fopen(fname, "a");
	  if (logfile == 0)
	       perror("opening log file");
     } else {
	  logfile = stdout;
     }
}

cuts_t WWLooperBase::DilepSelect (int i_hyp)
{
     cuts_t ret = 0;

     // pt cuts
     if (cms2.hyp_lt_p4()[i_hyp].pt() > 20.0) 
	  ret |= (1 << WW_LT_PT);
     if (cms2.hyp_ll_p4()[i_hyp].pt() > 20.0) 
	  ret |= (1 << WW_LL_PT);
     // sign cuts
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) 
	  ret |= (1 << WW_OPP_SIGN);
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) 
	  ret |= (1 << WW_SAME_SIGN);
     // MET
     if (pass4Met(i_hyp))
	  ret |= (1 << WW_PASS4_MET);
     if (pass2Met(i_hyp))
	  ret |= (1 << WW_PASS2_MET);
     // muon quality
     int n_iso_mu = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
	  ret |= (1 << WW_LT_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
	  ret |= (1 << WW_LL_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (1 << WW_LT_GOOD) | (1 << WW_LT_ISO);
	  n_iso_mu++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (1 << WW_LL_GOOD) | (1 << WW_LL_ISO);
	  n_iso_mu++;
     }
     // electron quality
     int n_iso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
	  ret |= (1 << WW_LT_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
	  ret |= (1 << WW_LL_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], false)) {
	  ret |= (1 << WW_LT_GOOD) | (1 << WW_LT_ISO);
	  n_iso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], false)) {
	  ret |= (1 << WW_LL_GOOD) | (1 << WW_LL_ISO);
	  n_iso_el++;
     }     
     if (n_iso_mu + n_iso_el >= 1)
	  ret |= (1 << WW_ONE_ISO);
     if (n_iso_mu + n_iso_el >= 2)
	  ret |= (1 << WW_TWO_ISO);
     // calo iso
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (1 << WW_LT_GOOD) | (1 << WW_LT_CALOISO);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (1 << WW_LL_GOOD) | (1 << WW_LL_CALOISO);
     }
     int n_caloiso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
	  ret |= (1 << WW_LT_GOOD) | (1 << WW_LT_CALOISO);
	  n_caloiso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
	  ret |= (1 << WW_LL_GOOD) | (1 << WW_LL_CALOISO);
	  n_caloiso_el++;
     }     
     if (n_iso_mu + n_caloiso_el >= 1)
	  ret |= (1 << WW_ONE_CALOISO);
     if (n_iso_mu + n_caloiso_el >= 2)
	  ret |= (1 << WW_TWO_CALOISO);
     // jet veto
     if (cms2.hyp_njets()[i_hyp] == 0)
	  ret |= (1 << WW_PASS_JETVETO_CALO);
     // track jets
     if (passTrkJetVeto(i_hyp))
	  ret |= (1 << WW_PASS_JETVETO_TRACKJETS);
     // muon b tag, with 20 GeV upper cut on the muon pt
     if (passMuonBVeto(i_hyp, true))
	  ret |= (1 << WW_PASS_MUON_B_VETO);
     else ret |= (1 << WW_MUON_TAGGED);
     // muon b tag, with no upper cut on the muon pt
     if (passMuonBVeto(i_hyp, false))
	  ret |= (1 << WW_PASS_MUON_B_VETO_WITHOUT_PTCUT);
     else ret |= (1 << WW_MUON_TAGGED_WITHOUT_PTCUT);
     // Z veto
     if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2 ||
	 not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
	  ret |= (1 << WW_PASS_ZVETO);
     // Z veto using additional leptons in the event
     if (not additionalZveto())
	  ret |= (1 << WW_PASS_ADDZVETO);
     // any additional high-pt, isolated leptons?
     if (passTriLepVeto(i_hyp))
	  ret |= (1 << WW_PASS_EXTRALEPTON_VETO);
     // the return value gets cached, too
     cuts_passed = ret;

     //*****************************************************************
     // special handling for the fake rate cuts for now, because they
     // only work for emu
     //*****************************************************************
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     if (myType != DILEPTON_EMU)
	  return ret;
     // in addition, for the muons, check that they pass tight+iso
     // (since the fake rate is electron only right now)
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  if ((cuts_passed & (1 << WW_LT_GOOD | 1 << WW_LT_ISO)) != 
	      (1 << WW_LT_GOOD | 1 << WW_LT_ISO))
	       return ret;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  if ((cuts_passed & (1 << WW_LL_GOOD | 1 << WW_LL_ISO)) != 
	      (1 << WW_LL_GOOD | 1 << WW_LL_ISO))
	       return ret;
     }
     // now set the fake flags for the electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (isFakeable_v2_2(cms2.hyp_lt_index()[i_hyp]))
	       ret |= 1 << WW_ELFAKE_FAKEABLE_OBJECT;
	  if (isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
	       ret |= 1 << WW_ELFAKE_NUMERATOR;
	  else ret |= 1 << WW_ELFAKE_NOT_NUMERATOR;
     } else {
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	       if (isFakeable_v2_2(cms2.hyp_ll_index()[i_hyp]))
		    ret |= 1 << WW_ELFAKE_FAKEABLE_OBJECT;
	       if (isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
		    ret |= 1 << WW_ELFAKE_NUMERATOR;
	       else ret |= 1 << WW_ELFAKE_NOT_NUMERATOR;
	  }
     }

     return cuts_passed = ret;
}

double WWLooperBase::Weight (int)
{
     return cms2.evt_scale1fb() * sample.kFactor;
}

void WWLooperBase::FillHists (int i_hyp)
{
     // common to all Fill() statements:
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     const double weight = Weight(i_hyp);
     
     if ((cuts_passed & cuts) == cuts) {
	  cands_passing[myType] += weight;
	  cands_passing_w2[myType] += weight * weight;
	  cands_count[myType]++;
	  cands_passing[DILEPTON_ALL] += weight;
	  cands_passing_w2[DILEPTON_ALL] += weight * weight;
	  cands_count[DILEPTON_ALL]++;
     }
//      printf("weight: %f\n", weight);
     // CLEANUP a lot of these histograms are much nicer as N - 1
     // plots!  Is it really worth keeping these after-cut plots?

     // jet count
     hnJet.Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp] + nTrkJets(i_hyp), weight);
     hnCaloJet	.Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     hnTrackJet	.Fill(cuts_passed, myType, nTrkJets(i_hyp), weight);

     // lepton pt's
     hminLepPt.Fill(cuts_passed, myType, std::min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt.Fill(cuts_passed, myType, std::max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hltPt.Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     hllPt.Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  helPt.Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  helEta.Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt.Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  hmuEta.Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  helPt.Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  helEta.Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt.Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  hmuEta.Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     }
    
     // dilepton mass
     hdilMass.Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
    
     // delta phi btw leptons
     double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
     if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
     hdphiLep.Fill(cuts_passed, myType, dphi, weight);
    
     // Relative isolation... muons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  hmuRelIso.Fill(cuts_passed, myType, reliso_lt(i_hyp), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  hmuRelIso.Fill(cuts_passed, myType, reliso_ll(i_hyp), weight);
     }

     // Relative isolation... electrons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  heleRelIso.Fill(cuts_passed, myType, reliso_lt(i_hyp), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  heleRelIso.Fill(cuts_passed, myType, reliso_ll(i_hyp), weight);
     }

     // lower of the two isolations, regardless of species (used for Dumbo method)
     hminRelIso.Fill(cuts_passed, myType, std::min(reliso_lt(i_hyp), reliso_ll(i_hyp)), weight);
     // lower of the two isolations, regardless of species (used for
     // Dumbo method) --- using calo iso for electrons
     hminRelIso_withCalo.Fill(cuts_passed, myType, 
			      std::min(reliso_lt(i_hyp, true), reliso_ll(i_hyp, true)), 
			      weight);

     // dilepton pt
     hdilPt.Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met and Met special
     hmet.Fill(cuts_passed, myType, cms2.hyp_met()[i_hyp], weight);      
     hmetSpec.Fill(cuts_passed, myType, 
		   MetSpecial(cms2.hyp_met()[i_hyp], cms2.hyp_metPhi()[i_hyp], i_hyp),
		   weight);
     // track correction to the met
     const TVector3 trkCorr = correctMETforTracks();
     TVector3 hyp_met;
     hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
     hyp_met += trkCorr;
     hmetTrkCorr.Fill(cuts_passed, myType, hyp_met.Perp(), weight);

     // tag muon pt and iso
     htagMuPt.Fill(cuts_passed, myType, tagMuonPt(i_hyp), weight);
     htagMuRelIso.Fill(cuts_passed, myType, tagMuonRelIso(i_hyp), weight);

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
	  hmuPdgId.Fill(myType, abs(cms2.mus_mc_id()[mu_idx]), weight);
	  hmuMoPdgId.Fill(myType, abs(cms2.mus_mc_motherid()[mu_idx]), weight);
	  helPdgId.Fill(myType, abs(cms2.els_mc_id()[el_idx]), weight);
	  helMoPdgId.Fill(myType, abs(cms2.els_mc_motherid()[el_idx]), weight);
	  helHE.Fill(cuts_passed, myType, cms2.els_hOverE()[el_idx], weight);
     }
}

void WWLooperBase::End ()
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
     if (logfile != stdout)
	  fclose(logfile);
     if (sample.chain != 0)
	  delete sample.chain;
}
