#include <sstream>
#include <iomanip>
#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"
#include "TDirectory.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
  : LooperBase(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));

// mybits.push_back(std::string("CUT_LT_PT"));
// mybits.push_back(std::string("CUT_LL_PT"));
// mybits.push_back(std::string("CUT_SAME_SIGN"));
// mybits.push_back(std::string("CUT_OPP_SIGN"));
// mybits.push_back(std::string("CUT_PASS2_MET"));
// mybits.push_back(std::string("CUT_PASS4_MET"));
// mybits.push_back(std::string("CUT_PASS2_TCMET"));
// mybits.push_back(std::string("CUT_PASS4_TCMET"));
// mybits.push_back(std::string("CUT_LT_GOOD"));
// mybits.push_back(std::string("CUT_LL_GOOD"));
// mybits.push_back(std::string("CUT_EL_GOOD"));
// mybits.push_back(std::string("CUT_MU_GOOD"));
// mybits.push_back(std::string("CUT_LT_ISO"));
// mybits.push_back(std::string("CUT_LL_ISO"));
// mybits.push_back(std::string("CUT_ONE_ISO"));
// mybits.push_back(std::string("CUT_TWO_ISO"));
// mybits.push_back(std::string("CUT_EL_ISO"));
// mybits.push_back(std::string("CUT_MU_ISO"));
// mybits.push_back(std::string("CUT_LT_CALOISO"));
// mybits.push_back(std::string("CUT_LL_CALOISO"));
// mybits.push_back(std::string("CUT_ONE_CALOISO"));
// mybits.push_back(std::string("CUT_TWO_CALOISO"));
// mybits.push_back(std::string("CUT_EL_CALOISO"));
// mybits.push_back(std::string("CUT_PASS_ZVETO"));
// mybits.push_back(std::string("CUT_IN_Z_WINDOW"));
// mybits.push_back(std::string("CUT_MUON_TAGGED"));
// mybits.push_back(std::string("CUT_ELFAKE_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_ELFAKE_NUMERATOR"));
// mybits.push_back(std::string("CUT_ELFAKE_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_MUFAKE_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_LT_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_MUFAKE_LT_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_LT_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_LL_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_MUFAKE_LL_NUMERATOR"));
// mybits.push_back(std::string("CUT_MUFAKE_LL_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_ELFAKE_LT_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_ELFAKE_LT_NUMERATOR"));
// mybits.push_back(std::string("CUT_ELFAKE_LT_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_ELFAKE_LL_FAKEABLE_OBJECT"));
// mybits.push_back(std::string("CUT_ELFAKE_LL_NUMERATOR"));
// mybits.push_back(std::string("CUT_ELFAKE_LL_NOT_NUMERATOR"));
// mybits.push_back(std::string("CUT_MORE_THAN_TWO_TRACKS"));
// mybits.push_back(std::string("CUT_TRUE_MU_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_MU_FROM_W_WJETS"));
// mybits.push_back(std::string("CUT_TRUE_LT_MU_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_LT_MU_FROM_W_WJETS"));
// mybits.push_back(std::string("CUT_TRUE_LL_MU_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_LL_MU_FROM_W_WJETS"));
// mybits.push_back(std::string("CUT_TRUE_EL_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_EL_FROM_W_WJETS"));
// mybits.push_back(std::string("CUT_TRUE_LT_EL_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_LT_EL_FROM_W_WJETS"));
// mybits.push_back(std::string("CUT_TRUE_LL_EL_FROM_W"));
// mybits.push_back(std::string("CUT_TRUE_LL_EL_FROM_W_WJETS"));
}

void Looper::BookHistos ()
{

  gDirectory = histo_directory;

  const unsigned int jetNBins = 5;
  float jetBins[jetNBins+1] = {0.,1.,2.,3.,4.,5.};

  const unsigned int ptNBins = 16;
  float ptBins[ptNBins+1];
  for ( unsigned int ptBin = 0;
	ptBin <= ptNBins;
	++ptBin) {
    ptBins[ptBin] = float(ptBin)*160./16.;
  }

  const unsigned int metNBins = 100;
  float metBins[metNBins+1];
  for ( unsigned int metBin = 0;
	metBin <= metNBins;
	++metBin) {
    metBins[metBin] = float(metBin)*200./100.;
  }

  const unsigned int etaNBins = 12;
  float etaBins[etaNBins+1];
  for ( unsigned int etaBin = 0;
	etaBin <= etaNBins;
	++etaBin) {
    etaBins[etaBin] = float(etaBin)/2.-3.;
  }

  for (unsigned int bucket = 0;
       bucket < 4;
       ++bucket ) {
    hnJet_[bucket] = book1DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"nJet",dilepton_hypo_names[bucket]),
				   Form("%s_%s_%s",sample_.name.c_str(),"nJet",dilepton_hypo_names[bucket]),
				   jetNBins,jetBins,
				   "n_{jet}","Events",sample_.histo_color);
    helPt_[bucket] = book1DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"elPt",dilepton_hypo_names[bucket]),
				   Form("%s_%s_%s",sample_.name.c_str(),"elPt",dilepton_hypo_names[bucket]),
				   ptNBins,ptBins,
				   "p_{T}^{e} [GeV]","Events",sample_.histo_color);
    helEta_[bucket] = book1DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"elEta",dilepton_hypo_names[bucket]),
				    Form("%s_%s_%s",sample_.name.c_str(),"elEta",dilepton_hypo_names[bucket]),
				    etaNBins,etaBins,
				    "#eta^{e} [GeV]","Events",sample_.histo_color);
    hmet_[bucket] = book1DVarHist(Form("%s_%s_%s",sample_.name.c_str(),"met",dilepton_hypo_names[bucket]),
				  Form("%s_%s_%s",sample_.name.c_str(),"met",dilepton_hypo_names[bucket]),
				  metNBins,metBins,
				  "MET [GeV]","Events",sample_.histo_color);
  }
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
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
    //    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_MU_GOOD);
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
    //    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_MU_GOOD);
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
    //    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
    n_iso_mu++;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
    //    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
    n_iso_mu++;
  }
  // electron quality
  int n_iso_el = 0;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_EL_GOOD);
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_EL_GOOD);
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_lt_index()[i_hyp], false)) {
    ret |= CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_EL_ISO);
    n_iso_el++;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_ll_index()[i_hyp], false)) {
    ret |= CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_EL_ISO);
    n_iso_el++;
  }     
  if (n_iso_mu + n_iso_el >= 1)
    ret |= (CUT_BIT(CUT_ONE_ISO));
  if (n_iso_mu + n_iso_el >= 2)
    ret |= (CUT_BIT(CUT_TWO_ISO));
//   // electrons without d0
//   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_lt_index()[i_hyp]) )
//     ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
//   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_ll_index()[i_hyp]) )
//     ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
  //      // supertight cuts (only for electrons)
  //      int n_supertight_el = 0;
  //      if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
  // 	  if (supertightElectron(cms2.hyp_lt_index()[i_hyp])) {
  // // 	       ret |= (CUT_BIT(CUT_LT_SUPERTIGHT));
  // 	       n_supertight_el++;
  // 	  }
  //      }
  //      if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
  // 	  if (supertightElectron(cms2.hyp_ll_index()[i_hyp])) {
  // // 	       ret |= (CUT_BIT(CUT_LL_SUPERTIGHT));
  // 	       n_supertight_el++;
  // 	  }
  //      }
  //      if (n_supertight_el >= 1)
  // 	  ret |= CUT_BIT(CUT_ONE_SUPERTIGHT);
  //      if (n_supertight_el >= 2)
  //    	  ret |= CUT_BIT(CUT_TWO_SUPERTIGHT);
//   // supertight dphiin cut
//   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
//     if (deltaPhiInElectron(cms2.hyp_lt_index()[i_hyp])) {
//       ret |= CUT_BIT(CUT_LT_TIGHT_DPHIIN);
//     }
//   }
//   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
//     if (deltaPhiInElectron(cms2.hyp_ll_index()[i_hyp])) {
//       ret |= CUT_BIT(CUT_LL_TIGHT_DPHIIN);
//     }
//   }
//   // barrel
//   if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
//     if (TMath::Abs(cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].eta()) < 1.479)
//       ret |= (CUT_BIT(CUT_EL_BARREL));
//   }
//   if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
//     if (TMath::Abs(cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].eta()) < 1.479)
//       ret |= (CUT_BIT(CUT_EL_BARREL));
//   }
  // calo iso
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
  }
  int n_caloiso_el = 0;
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
    n_caloiso_el++;
  }
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
    n_caloiso_el++;
  }     
  if (n_iso_mu + n_caloiso_el >= 1)
    ret |= (CUT_BIT(CUT_ONE_CALOISO));
  if (n_iso_mu + n_caloiso_el >= 2)
    ret |= (CUT_BIT(CUT_TWO_CALOISO));
//   // jet veto
//   if (cms2.hyp_njets()[i_hyp] == 0)
//     ret |= (CUT_BIT(CUT_PASS_JETVETO_CALO));
//   // track jets
//   if (passTrkJetVeto(i_hyp))
//     ret |= (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS));
  //      // muon b tag, with 20 GeV upper cut on the muon pt
  //      if (passMuonBVeto(i_hyp, true))
  // 	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO));
  //      else ret |= (CUT_BIT(CUT_MUON_TAGGED));
  //      // muon b tag, with no upper cut on the muon pt
  //      if (passMuonBVeto(i_hyp, false))
  // 	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT));
  //      else ret |= (CUT_BIT(CUT_MUON_TAGGED_WITHOUT_PTCUT));
  // Z veto
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2)
    ret |= (CUT_BIT(CUT_PASS_ZVETO));
  else if (not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
    ret |= (CUT_BIT(CUT_PASS_ZVETO));
  else ret |= (CUT_BIT(CUT_IN_Z_WINDOW));
//   // Z veto using additional leptons in the event
//   if (not additionalZveto())
//     ret |= (CUT_BIT(CUT_PASS_ADDZVETO));
//   // any additional high-pt, isolated leptons?
//   if (passTriLepVeto(i_hyp))
//     ret |= (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO));

  //      if (myType == DILEPTON_MUMU) { // don't want to deal with electron overlap right now
  // 	  if (passCaloTrkjetCombo())
  // 	       ret |= CUT_BIT(CUT_PASS_JETVETO_CALOTRACKJETS_COMBO);
  //      }

  // tag mu from w
  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && trueMuonFromW(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_MU_FROM_W) | CUT_BIT(CUT_TRUE_LT_MU_FROM_W);

  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && trueMuonFromW(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_MU_FROM_W) | CUT_BIT(CUT_TRUE_LL_MU_FROM_W);

  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && trueMuonFromW_WJets(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_MU_FROM_W_WJETS) | CUT_BIT(CUT_TRUE_LT_MU_FROM_W_WJETS);

  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && trueMuonFromW_WJets(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_MU_FROM_W_WJETS) | CUT_BIT(CUT_TRUE_LL_MU_FROM_W_WJETS);

  // tag el from w
  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && trueElectronFromW(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_EL_FROM_W) | CUT_BIT(CUT_TRUE_LT_EL_FROM_W);

  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && trueElectronFromW(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_EL_FROM_W) | CUT_BIT(CUT_TRUE_LL_EL_FROM_W);

  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && trueElectronFromW_WJets(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_EL_FROM_W_WJETS) | CUT_BIT(CUT_TRUE_LT_EL_FROM_W_WJETS);

  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && trueElectronFromW_WJets(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_TRUE_EL_FROM_W_WJETS) | CUT_BIT(CUT_TRUE_LL_EL_FROM_W_WJETS);

  if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11 && !trueGammaFromMuon(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_NOT_TRUE_LT_GAMMA_FROM_MUON) | CUT_BIT(CUT_NOT_TRUE_GAMMA_FROM_MUON);

  if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11 && !trueGammaFromMuon(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_NOT_TRUE_LL_GAMMA_FROM_MUON) | CUT_BIT(CUT_NOT_TRUE_GAMMA_FROM_MUON);

  // *****************************************************************
  // special handling for the fake rate cuts for now, because they
  // only work for emu
  // *****************************************************************
  if (myType != DILEPTON_EMU)
    return ret;
  // Here we try to find an event with a true Muon from W and a
  // fakable (denominator) electron
  // True muon is LT and fake electron is LL:
  if( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) { // just a protection
    if ((ret & ( CUT_BIT(CUT_TRUE_LT_MU_FROM_W_WJETS ) ) ) ==
        ( CUT_BIT(CUT_TRUE_LT_MU_FROM_W_WJETS))) {
      if (isFakeable(cms2.hyp_ll_index()[i_hyp]))
        ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(CUT_ELFAKE_LL_FAKEABLE_OBJECT);
      if (isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
        ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR) | CUT_BIT(CUT_ELFAKE_LL_NUMERATOR);
      else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR) | CUT_BIT(CUT_ELFAKE_LL_NOT_NUMERATOR);
    }
  }
  // True muon is LL and fake electron is LT:
  if( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) { // just a protection
    if ((ret & ( CUT_BIT(CUT_TRUE_LL_MU_FROM_W_WJETS) ) ) ==
        ( CUT_BIT(CUT_TRUE_LL_MU_FROM_W_WJETS) ) ) {
      if (isFakeable(cms2.hyp_lt_index()[i_hyp]))
        ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(CUT_ELFAKE_LT_FAKEABLE_OBJECT);
      if (isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
        ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR) | CUT_BIT(CUT_ELFAKE_LT_NUMERATOR);
      else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR) | CUT_BIT(CUT_ELFAKE_LT_NOT_NUMERATOR);
    }
  }
  return ret;
}

double Looper::Weight (int)
{
  return cms2.evt_scale1fb() * sample_.kFactor;
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

void Looper::FillDilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton histo filling; edit for your application
  //------------------------------------------------------------

  // these are the cuts that the candidate passes:
  cuts_t cuts_passed = DilepSelect(i_hyp);

//   ostringstream stream;

//   if ( cms2.evt_event() == 104 && cms2.evt_lumiBlock() == 11188 ) {
//     for ( unsigned int i = 0; i < 64; ++i ) {
//       if ((cuts_passed & ( CUT_BIT(i) ) ) == ( CUT_BIT(i))) stream << "fulfills: " << mybits[i] << endl;
//     }
//   }
    
  if ((cuts_passed & cuts_) == cuts_) {
    // every histogram needs to know what hypothesis he is 
    const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
    // and what the event weight is 
    const double weight = Weight(i_hyp);
     
    cands_passing_[myType] += weight;
    cands_passing_w2_[myType] += weight * weight;
    cands_count_[myType]++;
    cands_passing_[DILEPTON_ALL] += weight;
    cands_passing_w2_[DILEPTON_ALL] += weight * weight;
    cands_count_[DILEPTON_ALL]++;


//     stream << "Run: " << cms2.evt_run() << " event: " << cms2.evt_event() << " lumi: " << cms2.evt_lumiBlock() << endl;

    // jet count
    hnJet_[myType]->Fill(cms2.hyp_njets()[i_hyp], weight);
    hnJet_[DILEPTON_ALL]->Fill(cms2.hyp_njets()[i_hyp], weight);

    if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
      helPt_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
      helPt_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i_hyp].pt(), weight);
      helEta_[myType]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
      helEta_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i_hyp].eta(), weight);
    }
    if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
      helPt_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
      helPt_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i_hyp].pt(), weight);
      helEta_[myType]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
      helEta_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    }
    
    // Met and Met special
    hmet_[myType]->Fill(cms2.hyp_met()[i_hyp], weight);      
    hmet_[DILEPTON_ALL]->Fill(cms2.hyp_met()[i_hyp], weight);      
  }
//   cout << stream.str();

//   int ret = fprintf(logfile_, stream.str().c_str());
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

  ostringstream stream;
  
  stream << endl << "=========" << endl;
  stream << "Sample: " << SampleName().c_str() << endl;
  stream << "=========" << endl;
  stream << "Total candidate count ee: " << CandsCount(DILEPTON_EE) 
	 << " em: " << CandsCount(DILEPTON_EMU) 
	 << " mm: " << CandsCount(DILEPTON_MUMU)
	 << " all: " << CandsCount(DILEPTON_ALL) << endl;
  stream << "Total weight: ee: " << fixed << setprecision(1) << CandsPassing(DILEPTON_EE) << "+-" << RMS(DILEPTON_EE)
	 << " em: " << CandsPassing(DILEPTON_EMU) << "+-" << RMS(DILEPTON_EMU)
	 << " mm: " << CandsPassing(DILEPTON_MUMU) << "+-" << RMS(DILEPTON_MUMU)
	 << " all: " << CandsPassing(DILEPTON_ALL) << "+-" << RMS(DILEPTON_ALL) << endl;
  stream << "=========" << endl;
  
  cout << stream.str();

  int ret = fprintf(logfile_, stream.str().c_str());
  if (ret < 0)
    perror("writing to log file");
}
