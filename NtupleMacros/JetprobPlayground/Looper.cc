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

bool passTrkjetCuts (int i_trk);
void calculateJetProb (const std::vector<std::pair<LorentzVector, std::vector<unsigned int> > > &trackJets,
		       vector<double> *jetprobs);

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
  hnJPTJet     = new NMinus1Hist(sample_, "nJPTJet"    , 	  6,	 -0.5, 	    5,	 cuts_, 0);
  hntracks     = new NMinus1Hist(sample_, "ntracks"    , 	101,	 -0.5, 	100.5,	 cuts_, 0);
  hJPTJetPt    = new NMinus1Hist(sample_, "JPTJetPt"   , 	150,	    0, 	  150,	 cuts_, 0);

  htrd0 	       	= new NMinus1Hist(sample_, "trd0"           ,	100,	-0.1,	 0.1,	 cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
  htrd0sig        	= new NMinus1Hist(sample_, "trd0sig"        ,	100,	 -10,	  10,	 cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
  htrd0Strange    	= new NMinus1Hist(sample_, "trd0Strange"    ,	100,	-0.1,	 0.1,	 cuts_, 0); // mother is a Ks or Lambda
  htrd0sigStrange 	= new NMinus1Hist(sample_, "trd0sigStrange" ,	100,	 -10,	  10,	 cuts_, 0); // mother is a Ks or Lambda
  hjetpt 	       	= new NMinus1Hist(sample_, "jetpt"          ,	100,	   0,	 100,	 cuts_, 0);
  hjetntr	       	= new NMinus1Hist(sample_, "jetntr"         ,	 21,	-0.5,	20.5,	 cuts_, 0);
  htrkprob	       	= new NMinus1Hist(sample_, "trkprob"        ,	100,	   0,	   1,	 cuts_, 0);
  hjetprob	       	= new NMinus1Hist(sample_, "jetprob"        ,	100,	   0,	   1,	 cuts_, 0);
  hlogjetprob     	= new NMinus1Hist(sample_, "logjetprob"     ,	 21,	 -20,	   1,	 cuts_, 0);
  hminjetprob     	= new NMinus1Hist(sample_, "minjetprob"     ,	120,	   0,	 1.2,	 cuts_, 0);
  hmaxptjetprob   	= new NMinus1Hist(sample_, "maxptjetprob"   ,	120,	   0,	 1.2,	 cuts_, 0);

  hmaxntrksPt10	     = new NMinus1Hist(sample_, "maxntrksPt10"      ,   10,  0,      10, cuts_, 0);
  hjetprobNtrks2Pt10 = new NMinus1Hist(sample_, "jetprobNtrks2Pt10" ,  100,  0,       1, cuts_, 0);

  for (int i = 0; i < 3; ++i) {
    htrkd0[i]        = new NMinus1Hist(sample_, Form("%s%d", "trkd0"       , i), 100,   -1,    1,  cuts_, 0);
    htrknchi2[i]     = new NMinus1Hist(sample_, Form("%s%d", "trknchi2"    , i), 100,    0,   10,  cuts_, 0);
    htrkvalidhits[i] = new NMinus1Hist(sample_, Form("%s%d", "trkvalidhits", i),  21, -0.5, 20.5,  cuts_, 0);
  }

  for (int i = 0; i < 4; ++i) {
    htrd0ByPt[i]	    = new NMinus1Hist(sample_, Form("%s%d", "trd0ByPt"          , i),	100, -0.1,   0.1,	cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    htrd0sigByPt[i]	    = new NMinus1Hist(sample_, Form("%s%d", "trd0sigByPt"       , i),	100,  -10,    10,	cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    htrd0ByNtrks[i]	    = new NMinus1Hist(sample_, Form("%s%d", "trd0ByNtrks"       , i),	100, -0.1,   0.1,	cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    htrd0sigByNtrks[i]	    = new NMinus1Hist(sample_, Form("%s%d", "trd0sigByNtrks"    , i),	100,  -10,    10,	cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    hmaxptjetprobByPt[i]    = new NMinus1Hist(sample_, Form("%s%d", "maxptjetprobByPt"  , i),	100,    0,     1,	cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    hjetprobByPt[i]         = new NMinus1Hist(sample_, Form("%s%d", "jetprobByPt"       , i),   100,    0,     1,       cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
    hjetprobByNtrks[i]      = new NMinus1Hist(sample_, Form("%s%d", "jetprobByNtrks"    , i),   100,    0,     1,       cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
  }

  for (int i = 0; i < 5; ++i) {
    hmaxptjetprobNtrks5[i]  = new NMinus1Hist(sample_, Form("%s%d", "maxptjetprobNtrks5", i),  20,  0,  20, cuts_, (CUT_BIT(CUT_PASS_JETVETO_JPT20)) );
  }

  hntrksvsjptpt      = new TH2D(Form("%s_ntrksvsjptpt_all"    , sample_.name.c_str()), Form("%s_ntrksvsjptpt_all"    , sample_.name.c_str()),  10,  0,  10, 102,    -2,  100);

  hntrksvsjptpt->SetLineColor(sample_.histo_color);

  hmud0sig      = new NMinus1Hist(sample_, "mud0sig"  ,	100,    -10,    10, cuts_, 0);
  hmutrkprob	= new NMinus1Hist(sample_, "mutrkprob",	101, -0.005, 1.005, cuts_, 0);
}

bool Looper::FilterEvent()
{ 
  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  if (cms2.trks_d0corr().size() == 0)
    return true;
  DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0corr()[0], 
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

  // muon quality
  int n_iso_mu = 0;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
    ret |= CUT_BIT(CUT_LL_GOOD);
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= CUT_BIT(CUT_LT_GOOD);
    n_iso_mu++;
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= CUT_BIT(CUT_LL_GOOD);
    n_iso_mu++;
  }

  // electron quality
  int n_iso_el = 0;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LT_GOOD);
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
    ret |= CUT_BIT(CUT_LL_GOOD);

  // calo iso
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
    ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
  }

  int n_caloiso_el = 0;
  if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO);
    n_caloiso_el++;
  }
  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
    ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO);
    n_caloiso_el++;
  }     

  // track jets
  if (nJPTs(i_hyp, 20) == 0)
    ret |= CUT_BIT(CUT_PASS_JETVETO_JPT20);

  // muon b tag, with no upper cut on the muon pt
  if (numberOfExtraMuons(i_hyp, false) == 0)
    ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT));

  // Z veto
  if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2)
    ret |= (CUT_BIT(CUT_PASS_ZVETO));
  else if (not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
    ret |= (CUT_BIT(CUT_PASS_ZVETO));

  //*****************************************************************
  // track jet veto with signed impact parameter 
  //*****************************************************************

  double max_pt    = 0;
  int    max_ntrks = 0;
  double jetprob   = 1.1;

  for (unsigned int itrkjet = 0; itrkjet < TrackJets().size(); ++itrkjet) { 
       
    jetprob = TrackJetProbs()[itrkjet];

    if( TrackJets()[itrkjet].second.size() > 2 && TrackJets()[itrkjet].first.pt() > max_pt && jetprob < 1e-2 ) {
      max_pt = TrackJets()[itrkjet].first.pt();

      if( TrackJets()[itrkjet].second.size() > max_ntrks )
	max_ntrks = TrackJets()[itrkjet].second.size();
    }
  }

  if( max_pt < 10 ) {
    ret |= (CUT_BIT(CUT_PASS_JETVETO_SIP));
  }

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

bool compareJetpt (const std::pair<LorentzVector, std::vector<unsigned int> > &j1, 
		   const std::pair<LorentzVector, std::vector<unsigned int> > &j2) 
{
  return j1.first.pt() > j2.first.pt();
}

void Looper::MakeTrackJets (int i_hyp)
{
  trackjets_.clear();
  for ( unsigned int ijptjet = 0; ijptjet < cms2.jpts_p4().size(); ++ijptjet) {

    const double etaMax = 3.0;
    const double vetoCone = 0.4;

    // get rid of jpts that are hypothesis leptons
    if ( TMath::Abs(cms2.jpts_p4()[ijptjet].eta()) > etaMax ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[ijptjet])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[ijptjet])) < vetoCone ) continue;

    TLorentzVector jptJet(cms2.jpts_p4()[ijptjet].Px(), cms2.jpts_p4()[ijptjet].Py(),
			  cms2.jpts_p4()[ijptjet].Pz(), cms2.jpts_p4()[ijptjet].E());

    vector<unsigned int> trkidx;
    // Loop over all tracks and match to the track jet
    for (unsigned int trkiter = 0; trkiter < cms2.trks_trk_p4().size(); ++trkiter)
      {

	if (passTrkjetCuts(trkiter)) {
	  TLorentzVector trk(cms2.trks_trk_p4()[trkiter].Px(), cms2.trks_trk_p4()[trkiter].Py(),
			     cms2.trks_trk_p4()[trkiter].Pz(), cms2.trks_trk_p4()[trkiter].E() );
	       
	  if (jptJet.DeltaR(trk) < 0.5) {
	    trkidx.push_back(trkiter);
	  }
	}
      }
	  
    std::pair<LorentzVector, std::vector<unsigned int> > jptjet(cms2.jpts_p4()[ijptjet],
								trkidx);
    trackjets_.push_back(jptjet);
  }
  // sort track jets by pt
  sort(trackjets_.begin(), trackjets_.end(), compareJetpt);
}

void Looper::MakeTrackJetProbs ()
{
  calculateJetProb(TrackJets(), &trackjet_jp_);
}

void Looper::FillDilepHistos (int i_hyp)
{
  // make trackjets (from the ntuple or (by virtual function call) by kt clustering)
  MakeTrackJets(i_hyp);

  // and cache jet probabilities
  MakeTrackJetProbs();

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
  }

  // use muon tracks to study track probability
  if( myType == DILEPTON_MUMU ) {

    double ltd0corr  = cms2.hyp_lt_d0corr()[i_hyp];
    double ltd0err   = cms2.hyp_lt_d0Err()[i_hyp];
    double ltsig     = ltd0corr / ltd0err;
    double lttrkprob = 1 - erf( fabs(ltsig) / TMath::Sqrt(2) );

    hmud0sig  ->Fill(cuts_passed, myType,     ltsig, weight);
    hmutrkprob->Fill(cuts_passed, myType, lttrkprob, weight); 

    double lld0corr = cms2.hyp_ll_d0corr()[i_hyp];
    double lld0err  = cms2.hyp_ll_d0Err()[i_hyp];
    double llsig    = lld0corr / lld0err;
    double lltrkprob = 1 - erf( fabs(llsig) / TMath::Sqrt(2) );

    hmud0sig  ->Fill(cuts_passed, myType,     llsig, weight);
    hmutrkprob->Fill(cuts_passed, myType, lltrkprob, weight); 

  }

  hntracks->Fill(cuts_passed, myType, cms2.trks_trk_p4().size(), weight);

  // signed impact parameter study (FG3 && J.Mü.)
  int max_ntrks            = 0;
  double tmp_max_ptjet     = 0;
  double min_jetprob       = 1.1;
  double jetprob_max_ptjet = 1.1;
  double max_ptjet[5]      = {0, 0, 0, 0, 0};
  const double sigma       = 1.15;

  // loop over jets
  for (unsigned int itrkjet = 0; itrkjet < TrackJets().size(); ++itrkjet) {

    if ((cuts_passed & cuts_) == cuts_)
      hntrksvsjptpt->Fill( TrackJets()[itrkjet].second.size(), TrackJets()[itrkjet].first.pt() ); // ntrks vs pt for jetprob jets passing all cuts

    // store pt and ntrks for each jetprob jet
    hjetpt ->Fill(cuts_passed, myType, TrackJets()[itrkjet].first.pt()   , weight);
    hjetntr->Fill(cuts_passed, myType, TrackJets()[itrkjet].second.size(), weight);

    int i_pt;

    if (TrackJets()[itrkjet].first.pt() > 20)
      i_pt = 0;
    else if (TrackJets()[itrkjet].first.pt() > 15)
      i_pt = 1;
    else if (TrackJets()[itrkjet].first.pt() > 10)
      i_pt = 2;
    else
      i_pt = 3;
	  
    int i_ntr;

    if (TrackJets()[itrkjet].second.size() > 5)
      i_ntr = 0;
    else if (TrackJets()[itrkjet].second.size() > 3)
      i_ntr = 1;
    else if (TrackJets()[itrkjet].second.size() > 1)
      i_ntr = 2;
    else i_ntr = 3;
	  
    int n_jetprob_trks   = 0;
    double jet_sip       = 0;
    double jet_sipsig    = 0;
    double jet_maxsipsig = -999;

    // loop over tracks
    for (unsigned int i = 0; i < TrackJets()[itrkjet].second.size(); ++i) {

      std::pair<double, double> sip = signedImpact(TrackJets()[itrkjet].first, TrackJets()[itrkjet].second[i]);

      double d0           = sip.first;
      double d0err        = sip.second;
      const double sipsig = d0 / d0err;

      if (sipsig > jet_maxsipsig)
	jet_maxsipsig = sipsig;

      // store d0 for tracks in jetprob jets
      htrd0	            ->Fill(cuts_passed, myType,     d0, weight);
      htrd0ByPt[i_pt]       ->Fill(cuts_passed, myType,     d0, weight);
      htrd0ByNtrks[i_ntr]   ->Fill(cuts_passed, myType,     d0, weight);
 
      // store d0 significance for tracks in jetprob jets
      htrd0sig	            ->Fill(cuts_passed, myType, sipsig, weight);
      htrd0sigByPt[i_pt]    ->Fill(cuts_passed, myType, sipsig, weight);
      htrd0sigByNtrks[i_ntr]->Fill(cuts_passed, myType, sipsig, weight);

      if (sipsig > 0) {

	n_jetprob_trks++;

	double trkprob = 1 - erf( sipsig / TMath::Sqrt(2) / 1.1 );
	
	htrkprob->Fill(cuts_passed, myType, trkprob, weight); // track probability for jetprob jet tracks with sipsig > 0
      }

      unsigned int trkidx = TrackJets()[itrkjet].second[i]; 
      int momcid          = cms2.trk_mc_motherid()[trkidx];

      // see if we need to put a cutoff for strange hadrons
      switch (abs(momcid)) {
      case 130: case 310: case 311: // Kl/Ks/K0
      case 3122: // Lambda
      case 3222: case 3112: // Sigma+, Sigma-
      case 3322: case 3312: // Cascade0, Cascade-
      case 3334: // Omega-
	htrd0Strange   ->Fill(cuts_passed, myType,     d0, weight);
	htrd0sigStrange->Fill(cuts_passed, myType, sipsig, weight);
	break;
      default:
	break;
      }
    }
  
    double jetprob = TrackJetProbs()[itrkjet];

    if( jetprob < min_jetprob )
      min_jetprob = jetprob;

    // store jetprob
    hjetprob              ->Fill(cuts_passed, myType,      jetprob, weight);
    hlogjetprob           ->Fill(cuts_passed, myType, log(jetprob), weight);
    hjetprobByPt[i_pt]    ->Fill(cuts_passed, myType,      jetprob, weight);
    hjetprobByNtrks[i_ntr]->Fill(cuts_passed, myType,      jetprob, weight);

    if( TrackJets()[itrkjet].first.pt() > 10 && TrackJets()[itrkjet].second.size() > 2 )
      hjetprobNtrks2Pt10->Fill(cuts_passed, myType, jetprob, weight); // jetprob for jets with pt > 10 and ntrks > 2
    
    for( int numtrks = 1; numtrks < 6; numtrks++ ) {

      if(TrackJets()[itrkjet].second.size() > numtrks && TrackJets()[itrkjet].first.pt() > max_ptjet[numtrks-1] && jetprob < 1e-2) {
	max_ptjet[numtrks-1] = TrackJets()[itrkjet].first.pt();
	jetprob_max_ptjet = jetprob;
      }
    }

    if(TrackJets()[itrkjet].first.pt() > 10 && jetprob < 1e-2 && TrackJets()[itrkjet].second.size() > max_ntrks)
	max_ntrks = TrackJets()[itrkjet].second.size();
  }

  // store maximum pt and and ntrks and minimum jet prob for each event
  hminjetprob  ->Fill(cuts_passed, myType,       min_jetprob, weight);
  hmaxptjetprob->Fill(cuts_passed, myType, jetprob_max_ptjet, weight);
  hmaxntrksPt10->Fill(cuts_passed, myType,         max_ntrks, weight);

  for( unsigned int i = 0; i < 5; i++ ) {
    hmaxptjetprobNtrks5[i]->Fill(cuts_passed, myType, max_ptjet[i], weight);  // store jetprob of jet with highest pt with ntrks > X
  }
     
  int i_pt;
  if (max_ptjet[1] > 15)
    i_pt = 0;
  else if (max_ptjet[1] > 10)
    i_pt = 1;
  else if (max_ptjet[1] > 5)
    i_pt = 2;
  else i_pt = 3;

  hmaxptjetprobByPt[i_pt]->Fill(cuts_passed, myType, jetprob_max_ptjet, weight);
     
  // track quality
  for (unsigned int i = 0; i < cms2.trks_trk_p4().size(); ++i) {
    // quality index: 
    // 0 for tracks with |ip sig| < 2; 
    // 1 for 2 < // |ip sig| < 5; 
    // 2 for |ip sig| > 5
    int i_qual = 0;

    if (fabs(cms2.trks_d0corr()[i] / cms2.trks_d0Err()[i]) > 2)
      i_qual = 1;
    if (fabs(cms2.trks_d0corr()[i] / cms2.trks_d0Err()[i]) > 5)
      i_qual = 2;

    enum { D0, D0ERR, NCHI2, HITS };

    unsigned int trkcuts = 0;

    if (fabs(cms2.trks_d0corr()[i]) < 0.1)
      trkcuts |= 1 << D0;

    if( cms2.trks_d0Err()[i] > 20e-4 )
      trkcuts |= 1 << D0ERR;

    if (cms2.trks_chi2()[i] / cms2.trks_ndof()[i] < 5)
      trkcuts |= 1 << NCHI2;

    if (cms2.trks_validHits()[i] > 10)
      trkcuts |= 1 << HITS;

    const unsigned int all = 1 << D0 | 1 << D0ERR | 1 << NCHI2 | 1 << HITS;

    if (((trkcuts | 1 << D0) & all) == all)
      htrkd0[i_qual]->Fill(cuts_passed, myType, cms2.trks_d0corr()[i], weight);
	  
    if (((trkcuts | 1 << NCHI2) & all) == all)
      htrknchi2[i_qual]->Fill(cuts_passed, myType, cms2.trks_chi2()[i] / cms2.trks_ndof()[i], weight);

    if (((trkcuts | 1 << HITS) & all) == all)
      htrkvalidhits[i_qual]->Fill(cuts_passed, myType, cms2.trks_validHits()[i], weight);
  }

  hnJPTJet	->Fill(cuts_passed, myType, nJPTs(i_hyp, 20), weight); // JPT jet count

  // jet pt's (of the max jet) -- so much work we only do it for the ones we care about
  std::vector<LorentzVector> allJpts = JPTs(i_hyp, 0);
  double maxPtJPT = -1;

  for (unsigned int i = 0; i < allJpts.size(); ++i) {
    if (allJpts[i].pt() > maxPtJPT)
      maxPtJPT = allJpts[i].pt();
  }

  hJPTJetPt->Fill(cuts_passed, myType, maxPtJPT, weight); // pt of highest momentum JPT jet
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
