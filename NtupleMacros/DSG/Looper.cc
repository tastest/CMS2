#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
  : LooperBase(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_ ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_ ));
  memset(cands_count_		, 0, sizeof(cands_count_ ));
}

void Looper::BookHistos ()
{

  // njet plots
  hnJet_ = new NMinus1Hist(sample_, "nJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );
  hnCaloJet_ = new NMinus1Hist(sample_, "nCaloJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );
  hnTrackJet_ = new NMinus1Hist(sample_, "nTrackJet;n_{jets};Dilepton Cand." , 6 , -0.5, 5 , cuts_, 0 );

  // general lepton plots

  std::vector<cuts_t> lt_cuts;
  lt_cuts.push_back((CUT_BIT(CUT_LT_PT)));
  lt_cuts.push_back((CUT_BIT(CUT_LT_GOOD))); 
  lt_cuts.push_back((CUT_BIT(CUT_LT_CALOISO)));

  std::vector<std::string> lt_cut_names;
  lt_cut_names.push_back("_Cut_On_PT");
  lt_cut_names.push_back("_Cut_On_GOOD");
  lt_cut_names.push_back("_Cut_On_CALOISO");

  hltPt_ = new NMinus1Hist(sample_, "ltPt;p_{T} [GeV];Dilepton Cand." , 150 , 0, 150 , cuts_, lt_cuts, lt_cut_names);
  hltEta_ = new NMinus1Hist(sample_, "ltEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lt_cuts, lt_cut_names);
  hltCaloIso_ = new NMinus1Hist(sample_, "ltCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1., cuts_, lt_cuts, lt_cut_names);

  std::vector<cuts_t> ll_cuts;
  ll_cuts.push_back((CUT_BIT(CUT_LL_PT)));
  ll_cuts.push_back((CUT_BIT(CUT_LL_GOOD))); 
  ll_cuts.push_back((CUT_BIT(CUT_LL_CALOISO)));

  std::vector<std::string> ll_cut_names;
  ll_cut_names.push_back("_Cut_On_PT");
  ll_cut_names.push_back("_Cut_On_GOOD");
  ll_cut_names.push_back("_Cut_On_CALOISO");

  hllPt_ = new NMinus1Hist(sample_, "llPt;p_{T} [GeV];Dilepton Cand." , 150 , 0, 150 , cuts_, ll_cuts, ll_cut_names);
  hllEta_ = new NMinus1Hist(sample_, "llEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, ll_cuts, ll_cut_names);
  hllCaloIso_ = new NMinus1Hist(sample_, "llCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1., cuts_, ll_cuts, ll_cut_names);

  std::vector<cuts_t> lepton_cuts;
  lepton_cuts.push_back((CUT_BIT(CUT_LT_PT)));
  lepton_cuts.push_back((CUT_BIT(CUT_LT_GOOD))); 
  lepton_cuts.push_back((CUT_BIT(CUT_LT_CALOISO)));
  lepton_cuts.push_back((CUT_BIT(CUT_LL_PT)));
  lepton_cuts.push_back((CUT_BIT(CUT_LL_GOOD))); 
  lepton_cuts.push_back((CUT_BIT(CUT_LL_CALOISO)));

  std::vector<std::string> lepton_cut_names;
  lepton_cut_names.push_back("_Cut_On_lt_PT");
  lepton_cut_names.push_back("_Cut_On_lt_GOOD");
  lepton_cut_names.push_back("_Cut_On_lt_CALOISO");
  lepton_cut_names.push_back("_Cut_On_ll_PT");
  lepton_cut_names.push_back("_Cut_On_ll_GOOD");
  lepton_cut_names.push_back("_Cut_On_ll_CALOISO");

  helPt_ = new NMinus1Hist(sample_, "elPt;p_{T} [GeV];Dilepton Cand." , 16 , 0, 160 , cuts_, lepton_cuts, lepton_cut_names);
  hmuPt_ = new NMinus1Hist(sample_, "muPt;p_{T} [GeV];Dilepton Cand." , 16 , 0, 160 , cuts_, lepton_cuts, lepton_cut_names);
  helEta_ = new NMinus1Hist(sample_, "elEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lepton_cuts, lepton_cut_names);
  hmuEta_ = new NMinus1Hist(sample_, "muEta;#eta;Dilepton Cand." , 12 , -3, 3 , cuts_, lepton_cuts, lepton_cut_names);
  helCaloIso_ = new NMinus1Hist(sample_, "elCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1. , cuts_, lepton_cuts, lepton_cut_names);
  hmuCaloIso_ = new NMinus1Hist(sample_, "muCaloIso;calo iso;Dilepton Cand." , 50 , 0., 1. , cuts_, lepton_cuts, lepton_cut_names);

  // track plots
  hnTrack_ = new NMinus1Hist(sample_, "nTrack;n_{track};Dilepton Cand." , 200 , 0., 200. , cuts_, (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS)) );

  // dilepton plots
  hdphiLep_ = new NMinus1Hist(sample_, "dphiLep;#Delta#Phi;Dilepton Cand." , 50 , 0, M_PI , cuts_, 0 );
  hdilMass_ = new NMinus1Hist(sample_, "dilMass;inv. mass [GeV];Dilepton Cand." , 100 , 0, 300 , cuts_, 0 );
  hdilPt_ = new NMinus1Hist(sample_, "dilPt;p_{T} [GeV];Dilepton Cand." , 100 , 0, 300 , cuts_, 0 );

  // MET plots
  hmet_ = new NMinus1Hist(sample_, "met;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, (CUT_BIT(CUT_PASS4_MET)) | (CUT_BIT(CUT_PASS2_MET)) );
  hmetSpec_ = new NMinus1Hist(sample_, "metSpec;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, (CUT_BIT(CUT_PASS4_MET)) | (CUT_BIT(CUT_PASS2_MET)) );
  hmetTrkCorr_ = new NMinus1Hist(sample_, "metTrkCorr;MET [GeV];Dilepton Cand." , 100 , 0, 200 , cuts_, (CUT_BIT(CUT_PASS4_MET)) | (CUT_BIT(CUT_PASS2_MET)) );

}


bool Looper::FilterEvent()
{ 

  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  // comment in following lines
  // 

  if (cms2.trks_d0().size() == 0)
    return true;
  DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
			      cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
  if (is_duplicate(id)) {
    duplicates_total_n_++;
    duplicates_total_weight_ += cms2.evt_scale1fb();
//     cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
    return true;
  }

  return false; 
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
  //------------------------------------------------------------
  // Example dilepton cuts; edit for your application
  //------------------------------------------------------------

  // cuts are failed until proven otherwise
  cuts_t ret = 0;
  const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

  // temporarily calculate ONE global uncorrected jetSumEt:
  // prepare uncorrected jet collection
  vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets_p4(cms2.jets_p4());
  double sumEt = 0.;

  for (unsigned int icalojet=0; icalojet<calo_jets_p4.size(); ++icalojet) {
    // remove electron jets:
    if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)||
	(abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],calo_jets_p4[icalojet]) < 0.4)
	) continue;
    TLorentzVector p(calo_jets_p4[icalojet].Px()*cms2.jets_tq_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].Py()*cms2.jets_tq_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].Pz()*cms2.jets_tq_noCorrF()[icalojet],
		     calo_jets_p4[icalojet].E()*cms2.jets_tq_noCorrF()[icalojet]);//p is now uncorrected jet energy
    //       if (p.Perp() < jetet) continue;
    //       if (fabs(p.Eta()) > jeteta) continue;
    //        if (p.Perp() < 28 ) cout << p.Perp() << endl;
    if (fabs(p.Eta()) > 3.0 ) continue;
    if (p.Perp() < 15.) continue;
    //                 calo_jets->push_back(p);
    sumEt += p.Perp();
  }

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

  // DSG grid cuts
  if (met10(i_hyp, TVector3())) 
    ret |= (CUT_BIT(CUT_PASS_MET_10));
  if (met1(i_hyp, TVector3()))
    ret |= (CUT_BIT(CUT_PASS_MET_1));

  if (sumEt10( sumEt )) 
    ret |= (CUT_BIT(CUT_PASS_SUMET_10));
  if (sumEt1( sumEt ))
    ret |= (CUT_BIT(CUT_PASS_SUMET_1));

  // muon quality
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
    ret |= CUT_BIT(CUT_LL_GOOD);

  // electron quality
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LL_GOOD);
  
  // calo iso
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && passMuonIsolation(cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LT_CALOISO));
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && passMuonIsolation(cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LL_CALOISO));
  }
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_lt_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LT_CALOISO);
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_ll_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LL_CALOISO);
  } 

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

  // track z veto
  // use one of CUT_Z_TRACK_VETO_HYP, CUT_Z_TRACK_VETO_TRACKS or CUT_Z_TRACK_VETO_HYP_OR_TRACKS
  
  // int modus = passTrackZVeto(i_hyp);
  int modus = 999;

  // modus 1: combine trk-isolated tracks with hyp ll or lt, set CUT_Z_TRACK_VETO_HYP
  if ( ! ( modus == 1 || modus == 3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_HYP);

  // modus 2: combine trk-isolated tracks with another track, set CUT_Z_TRACK_VETO_TRACKS
  if ( ! ( modus == 2 || modus == 3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_TRACKS);

  // modus 3: combine trk-isolated tracks with another track or hyp ll or lt , set CUT_Z_TRACK_VETO_HYP_OR_TRACKS
  if ( ! ( modus == 1 || modus == 2 || modus ==3 ) ) ret |= CUT_BIT(CUT_Z_TRACK_VETO_HYP_OR_TRACKS);

  // *****************************************************************
  // special handling for the fake rate cuts for now, because they
  // only work for emu
  // *****************************************************************
  if (myType != DILEPTON_EMU)
    return ret;
  // in addition, for the muons, check that they pass tight+iso
  // (since the fake rate is electron only right now)
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
    if ((ret & (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO))) != 
	(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO)))
      return ret;
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
    if ((ret & (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO))) != 
	(CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO)))
      return ret;
  }
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

  return ret;
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
//     cout << "run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
    // if the candidate passed, we count it
    cands_passing_[myType] += weight;
    cands_passing_w2_[myType] += weight * weight;
    cands_count_[myType]++;
    cands_passing_[DILEPTON_ALL] += weight;
    cands_passing_w2_[DILEPTON_ALL] += weight * weight;
    cands_count_[DILEPTON_ALL]++;
  }

  // jet plots
  hnJet_->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp] + nTrkJets(i_hyp), weight);
  hnCaloJet_->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
  hnTrackJet_->Fill(cuts_passed, myType, nTrkJets(i_hyp), weight);

  // lepton plots
  hltPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
  hllPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
  hltEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
  hllEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);

  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
    helPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
    helEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
    hltCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_lt_index()[i_hyp],true), weight);
    helCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_lt_index()[i_hyp],true), weight);
  } else {
    hmuPt_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
    hmuEta_->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
    hltCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_lt_index()[i_hyp]), weight);
    hmuCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_lt_index()[i_hyp]), weight);
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
    helPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
    helEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    hllCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_ll_index()[i_hyp],true), weight);
    helCaloIso_->Fill(cuts_passed, myType, el_rel_iso(cms2.hyp_ll_index()[i_hyp],true), weight);
  } else {
    hmuPt_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
    hmuEta_->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
    hllCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_ll_index()[i_hyp]), weight);
    hmuCaloIso_->Fill(cuts_passed, myType, mu_rel_iso(cms2.hyp_ll_index()[i_hyp]), weight);
  }

  // track plots
  hnTrack_->Fill(cuts_passed, myType, cms2.trks_trk_p4().size(), weight);

  // dilepton mass
  hdilMass_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
    
  // delta phi btw leptons
  double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
  if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  hdphiLep_->Fill(cuts_passed, myType, dphi, weight);

  // dilepton pt
  hdilPt_->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].pt(), weight);
    
  // MET plots
  hmet_->Fill(cuts_passed, myType, cms2.hyp_met()[i_hyp], weight);      
  hmetSpec_->Fill(cuts_passed, myType, 
MetSpecial(cms2.hyp_met()[i_hyp], cms2.hyp_metPhi()[i_hyp], i_hyp),
weight);
  // track correction to the met
  const TVector3 trkCorr = correctMETforTracks();
  TVector3 hyp_met;
  hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
  hyp_met += trkCorr;
  hmetTrkCorr_->Fill(cuts_passed, myType, hyp_met.Perp(), weight);

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
    " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n", 
    sample_.name.c_str(),
    CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
    CandsPassing(DILEPTON_EE) , RMS(DILEPTON_EE), 
    CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU), 
    CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
    CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
    perror("writing to log file");
}
