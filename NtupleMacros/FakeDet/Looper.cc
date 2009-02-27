#include <math.h>
#include <sstream>
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
  events_            = 0;
  events_weighted_   = 0;
  events_passing_    = 0;
  events_passing_w2_ = 0;
  events_count_	     = 0;

  nDenominator_ = 0;
  nDenominator_wo_leading_ = 0;
  nDenominator_wo_second_leading_ = 0;
  nNumerator_ll_ = 0;
  nNumerator_ll_wo_leading_ = 0;
  nNumerator_ll_wo_second_leading_ = 0;
  nNumerator_lt_ = 0;
  nNumerator_lt_wo_leading_ = 0;
  nNumerator_lt_wo_second_leading_ = 0;

  nDenominator_weighted_ = 0;
  nDenominator_wo_leading_weighted_ = 0;
  nDenominator_wo_second_leading_weighted_ = 0;
  nNumerator_ll_weighted_ = 0;
  nNumerator_ll_wo_leading_weighted_ = 0;
  nNumerator_ll_wo_second_leading_weighted_ = 0;
  nNumerator_lt_weighted_ = 0;
  nNumerator_lt_wo_leading_weighted_ = 0;
  nNumerator_lt_wo_second_leading_weighted_ = 0;

  nDenominator_weighted_w2_ = 0;
  nDenominator_wo_leading_weighted_w2_ = 0;
  nDenominator_wo_second_leading_weighted_w2_ = 0;
  nNumerator_ll_weighted_w2_ = 0;
  nNumerator_ll_wo_leading_weighted_w2_ = 0;
  nNumerator_ll_wo_second_leading_weighted_w2_ = 0;
  nNumerator_lt_weighted_w2_ = 0;
  nNumerator_lt_wo_leading_weighted_w2_ = 0;
  nNumerator_lt_wo_second_leading_weighted_w2_ = 0;

}

void Looper::BookHistos ()
{
  //------------------------------------------------------------
  // Example histo booking; edit for your application
  //------------------------------------------------------------

  // sample name is the prefix for all histograms
  const char * prefix = SampleName().c_str();

  // sample color is used for histogram colors
  int color = sample_.histo_color;

  // coarse or fine fake rate binning?
  bool fineBinning = false;

  vector<float> binsPt;
  vector<float> binsEta;

  if(fineBinning) {

    // set Pt bins
    binsPt.push_back(0.);
    //	binsPt.push_back(10.);
    binsPt.push_back(20.);
    //	binsPt.push_back(30.);
    binsPt.push_back(40.);
    //	binsPt.push_back(50.);
    binsPt.push_back(65.);
    //	binsPt.push_back(80.);
    binsPt.push_back(95.);
    //	binsPt.push_back(115.);
    binsPt.push_back(150.);

    // set eta bins
    binsEta.push_back(-2.5);
    //	binsEta.push_back(-2.1);
    binsEta.push_back(-1.9);
    //	binsEta.push_back(-1.7);
    binsEta.push_back(-1.5);
    //	binsEta.push_back(-1.3);
    binsEta.push_back(-1.1);
    //	binsEta.push_back(-0.9);
    binsEta.push_back(-0.7);
    //	binsEta.push_back(-0.5);
    binsEta.push_back(-0.3);
    // 	binsEta.push_back(-0.1);
    // 	binsEta.push_back(0.1);
    binsEta.push_back(0.3);
    //	binsEta.push_back(0.5);
    binsEta.push_back(0.7);
    //	binsEta.push_back(0.9);
    binsEta.push_back(1.1);
    //	binsEta.push_back(1.3);
    binsEta.push_back(1.5);
    //	binsEta.push_back(1.7);
    binsEta.push_back(1.9);
    //	binsEta.push_back(2.1);
    binsEta.push_back(2.5);
  }
  else {

    // set Pt bins
//     binsPt.push_back(0.);
    binsPt.push_back(20.);
    binsPt.push_back(60.);
    binsPt.push_back(150.);

    // set eta bins
    binsEta.push_back(-2.5);
    binsEta.push_back(-1.479);
    binsEta.push_back(1.479);
    binsEta.push_back(2.5);
  }

  pt_num_ell_  = book1DVarHist(Form("%s_hpt_num_ell",prefix),"p_{T} Numerator",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator",color);
  pt_num_elt_  = book1DVarHist(Form("%s_hpt_num_elt",prefix),"p_{T} Numerator",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator",color);
  pt_den_ele_  = book1DVarHist(Form("%s_hpt_den_ele",prefix),"p_{T} Denominator",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator",color);

  eta_num_ell_ = book1DVarHist(Form("%s_heta_num_ell",prefix),"#eta Numerator",binsEta,"#eta","#eta Loose Electron Numerator",color);
  eta_num_elt_ = book1DVarHist(Form("%s_heta_num_elt",prefix),"#eta Numerator",binsEta,"#eta","#eta Tight Electron Numerator",color);
  eta_den_ele_ = book1DVarHist(Form("%s_heta_den_ele",prefix),"#eta Denominator",binsEta,"#eta","#eta Electron Denominator",color);

  num_ell_ = book2DVarHist(Form("%s_hnum_ell",prefix),"Numerator",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator",color);
  num_elt_ = book2DVarHist(Form("%s_hnum_elt",prefix),"Numerator",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator",color);
  den_ele_ = book2DVarHist(Form("%s_hden_ele",prefix),"Denominator",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator",color);

  pt_num_wo_leading_ell_  = book1DVarHist(Form("%s_hpt_num_wo_leading_ell",prefix),"p_{T} Numerator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator",color);
  pt_num_wo_leading_elt_  = book1DVarHist(Form("%s_hpt_num_wo_leading_elt",prefix),"p_{T} Numerator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator",color);
  pt_den_wo_leading_ele_  = book1DVarHist(Form("%s_hpt_den_wo_leading_ele",prefix),"p_{T} Denominator without leading Jet",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator",color);

  eta_num_wo_leading_ell_ = book1DVarHist(Form("%s_heta_num_wo_leading_ell",prefix),"#eta Numerator without leading Jet",binsEta,"#eta","#eta Loose Electron Numerator",color);
  eta_num_wo_leading_elt_ = book1DVarHist(Form("%s_heta_num_wo_leading_elt",prefix),"#eta Numerator without leading Jet",binsEta,"#eta","#eta Tight Electron Numerator",color);
  eta_den_wo_leading_ele_ = book1DVarHist(Form("%s_heta_den_wo_leading_ele",prefix),"#eta Denominator without leading Jet",binsEta,"#eta","#eta Electron Denominator",color);

  num_wo_leading_ell_ = book2DVarHist(Form("%s_hnum_wo_leading_ell",prefix),"Numerator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator",color);
  num_wo_leading_elt_ = book2DVarHist(Form("%s_hnum_wo_leading_elt",prefix),"Numerator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator",color);
  den_wo_leading_ele_ = book2DVarHist(Form("%s_hden_wo_leading_ele",prefix),"Denominator without leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator",color);

  pt_num_wo_second_leading_ell_  = book1DVarHist(Form("%s_hpt_num_wo_second_leading_ell",prefix),"p_{T} Numerator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Loose Electron Numerator",color);
  pt_num_wo_second_leading_elt_  = book1DVarHist(Form("%s_hpt_num_wo_second_leading_elt",prefix),"p_{T} Numerator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Tight Electron Numerator",color);
  pt_den_wo_second_leading_ele_  = book1DVarHist(Form("%s_hpt_den_wo_second_leading_ele",prefix),"p_{T} Denominator without second leading Jet",binsPt,"p_{T} [GeV]","p_{T} Electron Denominator",color);

  eta_num_wo_second_leading_ell_ = book1DVarHist(Form("%s_heta_num_wo_second_leading_ell",prefix),"#eta Numerator without second leading Jet",binsEta,"#eta","#eta Loose Electron Numerator",color);
  eta_num_wo_second_leading_elt_ = book1DVarHist(Form("%s_heta_num_wo_second_leading_elt",prefix),"#eta Numerator without second leading Jet",binsEta,"#eta","#eta Tight Electron Numerator",color);
  eta_den_wo_second_leading_ele_ = book1DVarHist(Form("%s_heta_den_wo_second_leading_ele",prefix),"#eta Denominator without second leading Jet",binsEta,"#eta","#eta Electron Denominator",color);

  num_wo_second_leading_ell_ = book2DVarHist(Form("%s_hnum_wo_second_leading_ell",prefix),"Numerator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Loose Electron Numerator",color);
  num_wo_second_leading_elt_ = book2DVarHist(Form("%s_hnum_wo_second_leading_elt",prefix),"Numerator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Tight Electron Numerator",color);
  den_wo_second_leading_ele_ = book2DVarHist(Form("%s_hden_wo_second_leading_ele",prefix),"Denominator without second leading Jet",binsEta,binsPt,"#eta","p_{T} [GeV]","Electron Denominator",color);
	
  ele_n_   = book1DHist(Form("%s_hele_n",prefix),"number of electrons per event",15,0.,15.,"n_{e}","Events",color);

  unsigned int nJetsBins = 50;
  float        nJetsLow  = 0.;
  float        nJetsHigh = 50.;

  njets_ = book1DHist(Form("%s_hnjets",prefix),"Number of Jets",nJetsBins,nJetsLow,nJetsHigh,"n_{Jets}","Events",color);

  unsigned int jetPtBins = 4000;
  float        jetPtLow  = 0.;
  float        jetPtHigh = 4000.;

  jetpt_ = book1DHist(Form("%s_hjetpt",prefix),"p_{T} of Jets",jetPtBins,jetPtLow,jetPtHigh,"p_{T} [GeV]","Events",color);

  unsigned int jetEtaBins = 30;
  float        jetEtaLow  = -3.;
  float        jetEtaHigh =  3.;

  jeteta_ = book1DHist(Form("%s_hjeteta",prefix),"#eta of Jets",jetEtaBins,jetEtaLow,jetEtaHigh,"#eta","Events",color);

  unsigned int deltaRBins = 500;
  float        deltaRLow  = 0.;
  float        deltaRHigh = 10.;

  deltaR_ = book1DHist(Form("%s_hdeltaR",prefix),"deltaR between jet and denominator",deltaRBins,deltaRLow,deltaRHigh,"#Delta R","Events",color);

  unsigned int tkIsoBins = 400;
  float        tkIsoLow  = 0.;
  float        tkIsoHigh = 200.;

  tkIso_ = book1DHist(Form("%s_htkIso",prefix),"track isolation",tkIsoBins,tkIsoLow,tkIsoHigh,"tk_{iso}","Events",color);
  tkIso_uncut_ = book1DHist(Form("%s_htkIso_uncut",prefix),"uncut track isolation",tkIsoBins,tkIsoLow,tkIsoHigh,"tk_{iso}","Events",color);

  unsigned int EOverpBins = 30;
  float        EOverpLow  = -1.5;
  float        EOverpHigh = 1.5;

  EOverp_ = book1DHist(Form("%s_hEOverp",prefix),"E/p",EOverpBins,EOverpLow,EOverpHigh,"E/p","Events",color);
  EOverp_uncut_ = book1DHist(Form("%s_hEOverp_uncut",prefix),"uncut E/p",EOverpBins,EOverpLow,EOverpHigh,"E/p","Events",color);

  unsigned int HOverEBins = 100;
  float        HOverELow  = 0.;
  float        HOverEHigh = 0.1;

  HOverE_ = book1DHist(Form("%s_hHOverE",prefix),"H/E",HOverEBins,HOverELow,HOverEHigh,"H/E","Events",color);
  HOverE_uncut_ = book1DHist(Form("%s_hHOverE_uncut",prefix),"uncut H/E",HOverEBins,HOverELow,HOverEHigh,"H/E","Events",color);

}


bool Looper::FilterEvent()
{ 

  //
  // duplicate filter, based on trk information and lepton candidate
  //
  // comment in following lines
  // 

  if (cms2.els_p4().size() > 0 ) {
    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.els_d0corr()[0], 
				cms2.els_p4()[0].pt(), cms2.els_p4()[0].eta(), cms2.els_p4()[0].phi() };
    if (is_duplicate(id)) {
      duplicates_total_n_++;
      duplicates_total_weight_ += cms2.evt_scale1fb();
      //cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
      return true;
    }
  } else if ( cms2.mus_p4().size() > 0 ) {
    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.mus_d0corr()[0], 
				cms2.mus_p4()[0].pt(), cms2.mus_p4()[0].eta(), cms2.mus_p4()[0].phi() };
    if (is_duplicate(id)) {
      duplicates_total_n_++;
      duplicates_total_weight_ += cms2.evt_scale1fb();
      //cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
      return true;
    }
  } else {
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

  // leading jet pt cut
  if ( cms2.jets_p4().size() > 0 ) {
    float HLT_jet_approx_uncorr 	= 30.0; // require leading jet to be larger than N GeV uncorrected v2_3
    float HLT_jet_approx 		= HLT_jet_approx_uncorr / cms2.jets_pat_noCorrF()[0]; // V00-05-01 Fake tuples have TQ corrected jets. Uncorreced cut is done here.
    if ( cms2.jets_p4()[0].Pt() >= HLT_jet_approx ) {
      ret |= (CUT_BIT(CUT_PT_LEADING_JET));
    }
  }

  ret |= (CUT_BIT(CUT_NO_CUT));

  if ( cms2.genps_pthat() < sample_.upper_pthat )
    ret |= (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));

  if ( events_ % 2 ) {
    ret |= (CUT_BIT(CUT_ODD));
  } else {
    ret |= (CUT_BIT(CUT_EVEN));
  }

  return ret;
}

cuts_t Looper::DilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton cuts; edit for your application
  //------------------------------------------------------------

  // cuts are failed until proven otherwise
  cuts_t ret = 0;

  // the return value gets cached, too
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

  // and what the event weight is 
  const double weight = Weight(0);

  ++events_;
  events_weighted_ += weight;

  // these are the cuts that the candidate passes:
  cuts_t cuts_passed = EventSelect();

  // this is how to test that the candidate passes the cuts (which
  // we specified in the constructor when we made the looper)
  // (*note: the parentheses are important*):
  if ((cuts_passed & cuts_) == cuts_) {
    // if the candidate passed, we count it
    events_passing_ += weight;
    events_passing_w2_ += weight * weight;
    ++events_count_;

    // fill general histograms
    ele_n_->Fill(cms2.evt_nels(),weight);

    njets_->Fill(cms2.evt_njets(),weight);
    for ( unsigned int jet_counter = 0;
	  jet_counter < (unsigned int)cms2.evt_njets();
	  ++jet_counter ) {
      jetpt_->Fill(cms2.jets_p4()[jet_counter].Pt(),weight);
      jeteta_->Fill(cms2.jets_p4()[jet_counter].Eta(),weight);
    }

    // cuts
    float deltaRCut = 0.2;
    bool  useAbsEta = true;

    // electron loop
    for ( unsigned int electron_counter = 0;
	  electron_counter < cms2.els_p4().size();
	  ++electron_counter ) {

      // general electron histograms
      tkIso_uncut_->Fill(cms2.els_tkIso()[electron_counter],weight);
      EOverp_uncut_->Fill(cms2.els_eOverPIn()[electron_counter],weight);
      HOverE_uncut_->Fill(cms2.els_hOverE()[electron_counter],weight);

      // fill every electron with pt > 150 GeV into last bin
      float pt = cms2.els_p4()[electron_counter].Pt();
      if (pt >= 150.) pt = 149.;

      if ( isFakeable(electron_counter)){
	  
	pt_den_ele_->Fill(pt);
	if ( useAbsEta ) {
	  eta_den_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  eta_den_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  den_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  den_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	} else {
	  eta_den_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	  den_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	}

	// hard cut to calculate zero-bin fake rate correctly
	if ( pt >= 20. ) {
	  ++nDenominator_;
	  nDenominator_weighted_ += weight;
	  nDenominator_weighted_w2_ += weight * weight;
	}

	for ( unsigned int jet_counter = 0;
	      jet_counter < (unsigned int)cms2.evt_njets();
	      ++jet_counter ) {
	  deltaR_->Fill(dRbetweenVectors(cms2.els_p4()[electron_counter],cms2.jets_p4()[jet_counter]),weight);
	}
	tkIso_->Fill(cms2.els_tkIso()[electron_counter],weight);
	EOverp_->Fill(cms2.els_eOverPIn()[electron_counter],weight);
	HOverE_->Fill(cms2.els_hOverE()[electron_counter],weight);
      }

      // loose electrons
      if (isNumeratorElectron(electron_counter,1)){ 
	// 1=loose, 2=tight
	if ( !isFakeable(electron_counter) ) cout << "Loose electron: num is fullfilled but den is not!" << endl;
	pt_num_ell_->Fill(pt,weight);
	if ( useAbsEta ) {
	  eta_num_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  eta_num_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  num_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  num_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	} else {
	  eta_num_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	  num_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	}

	// hard cut to calculate zero-bin fake rate correctly
	if ( pt >= 20. ) {
	  ++nNumerator_ll_;
	  nNumerator_ll_weighted_ += weight;
	  nNumerator_ll_weighted_w2_ += weight * weight;
	}
      }

      // tight electrons
      if ( isNumeratorElectron(electron_counter,2)){ 
	// 1=loose, 2=tight
	if ( !isFakeable(electron_counter) ) cout << "Tight electron: num is fullfilled but den is not!" << endl;
	pt_num_elt_->Fill(pt,weight);
	if ( useAbsEta ) {
	  eta_num_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  eta_num_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	  num_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  num_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	} else {
	  eta_num_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	  num_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	}
	// hard cut to calculate zero-bin fake rate correctly
	if ( pt >= 20. ) {
	  ++nNumerator_lt_;
	  nNumerator_lt_weighted_ += weight;
	  nNumerator_lt_weighted_w2_ += weight * weight;
	}
      }

      // exclude leading jet
      if ( dRbetweenVectors(cms2.els_p4()[electron_counter],cms2.jets_p4()[0]) >= deltaRCut ) {

	if ( isFakeable(electron_counter)){
	  pt_den_wo_leading_ele_->Fill(pt,weight);
	  if ( useAbsEta ) {
	    eta_den_wo_leading_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    eta_den_wo_leading_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    den_wo_leading_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    den_wo_leading_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  } else {
	    eta_den_wo_leading_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	    den_wo_leading_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	  }
	  // hard cut to calculate zero-bin fake rate correctly
	  if ( pt >= 20. ) {
	    ++nDenominator_wo_leading_;
	    nDenominator_wo_leading_weighted_ += weight;
	    nDenominator_wo_leading_weighted_w2_ += weight * weight;
	  }
	}

	// loose electrons
	if (isNumeratorElectron(electron_counter,1)){ 
	  // 1=loose, 2=tight
	  if ( !isFakeable(electron_counter) ) cout << "Loose electron wo leading: num is fullfilled but den is not!" << endl;
	  pt_num_wo_leading_ell_->Fill(pt,weight);
	  if ( useAbsEta ) {
	    eta_num_wo_leading_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    eta_num_wo_leading_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    num_wo_leading_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    num_wo_leading_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  } else {
	    eta_num_wo_leading_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	    num_wo_leading_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	  }
	  // hard cut to calculate zero-bin fake rate correctly
	  if ( pt >= 20. ) {
	    ++nNumerator_ll_wo_leading_;
	    nNumerator_ll_wo_leading_weighted_ += weight;
	    nNumerator_ll_wo_leading_weighted_w2_ += weight * weight;
	  }
	}

	// tight electrons
	if ( isNumeratorElectron(electron_counter,2)) { 
	  // 1=loose, 2=tight
	  if ( !isFakeable(electron_counter) ) cout << "Tight electron wo leading: num is fullfilled but den is not!" << endl;
	  pt_num_wo_leading_elt_->Fill(pt,weight);
	  if ( useAbsEta ) {
	    eta_num_wo_leading_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    eta_num_wo_leading_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	    num_wo_leading_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    num_wo_leading_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	  } else {
	    eta_num_wo_leading_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	    num_wo_leading_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	  }
	  // hard cut to calculate zero-bin fake rate correctly
	  if ( pt >= 20. ) {
	    ++nNumerator_lt_wo_leading_;
	    nNumerator_lt_wo_leading_weighted_ += weight;
	    nNumerator_lt_wo_leading_weighted_w2_ += weight * weight;
	  }
	}

      }

      // exclude leading and second leading jet
      if ( cms2.evt_njets() > 1 ) {

	if ( (dRbetweenVectors(cms2.els_p4()[electron_counter],cms2.jets_p4()[0]) >= deltaRCut) &&
	     (dRbetweenVectors(cms2.els_p4()[electron_counter],cms2.jets_p4()[1]) >= deltaRCut) ) {

	  if ( isFakeable(electron_counter)){
	    pt_den_wo_second_leading_ele_->Fill(pt,weight);
	    if ( useAbsEta ) {
	      eta_den_wo_second_leading_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      eta_den_wo_second_leading_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      den_wo_second_leading_ele_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	      den_wo_second_leading_ele_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    } else {
	      eta_den_wo_second_leading_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	      den_wo_second_leading_ele_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	    }
	    // hard cut to calculate zero-bin fake rate correctly
	    if ( pt >= 20. ) {
	      ++nDenominator_wo_second_leading_;
	      nDenominator_wo_second_leading_weighted_ += weight;
	      nDenominator_wo_second_leading_weighted_w2_ += weight * weight;
	    }
	  }

	  // loose electrons
	  if (isNumeratorElectron(electron_counter,1)){ 
	    // 1=loose, 2=tight
	    if ( !isFakeable(electron_counter) ) cout << "Loose electron wo 2nd leading: num is fullfilled but den is not!" << endl;
	    pt_num_wo_second_leading_ell_->Fill(pt,weight);
	    if ( useAbsEta ) {
	      eta_num_wo_second_leading_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      eta_num_wo_second_leading_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      num_wo_second_leading_ell_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	      num_wo_second_leading_ell_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    } else {
	      eta_num_wo_second_leading_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	      num_wo_second_leading_ell_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	    }
	    // hard cut to calculate zero-bin fake rate correctly
	    if ( pt >= 20. ) {
	      ++nNumerator_ll_wo_second_leading_;
	      nNumerator_ll_wo_second_leading_weighted_ += weight;
	      nNumerator_ll_wo_second_leading_weighted_w2_ += weight * weight;
	    }
	  }

	  // tight electrons
	  if ( isNumeratorElectron(electron_counter,2)){ 
	    // 1=loose, 2=tight
	    if ( !isFakeable(electron_counter) ) cout << "Tight electron wo 2nd leading: num is fullfilled but den is not!" << endl;
	    pt_num_wo_second_leading_elt_->Fill(pt,weight);
	    if ( useAbsEta ) {
	      eta_num_wo_second_leading_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      eta_num_wo_second_leading_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),weight);
	      num_wo_second_leading_elt_->Fill(TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	      num_wo_second_leading_elt_->Fill(-TMath::Abs(cms2.els_p4()[electron_counter].Eta()),pt,weight);
	    } else {
	      eta_num_wo_second_leading_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),weight);
	      num_wo_second_leading_elt_->Fill(cms2.els_p4()[electron_counter].Eta(),pt,weight);
	    }
	    // hard cut to calculate zero-bin fake rate correctly
	    if ( pt >= 20. ) {
	      ++nNumerator_lt_wo_second_leading_;
	      nNumerator_lt_wo_second_leading_weighted_ += weight;
	      nNumerator_lt_wo_second_leading_weighted_w2_ += weight * weight;
	    }
	  }
	}
      }
    }
  }
}

void Looper::FillDilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // Example dilepton histo filling; edit for your application
  //------------------------------------------------------------

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
  stream << events_count_ << " of " << events_ << " events passed event selection" << endl;
  stream << events_passing_ << " of " << events_weighted_ << " weighted events passed event selection" << endl;
  stream << duplicates_total_n_ << " duplicates corresponding to " << duplicates_total_weight_ << " weighted events excluded" << endl;      
  stream << "=========" << endl;

  stream << "Tight: num: " << nNumerator_lt_ << " den: " << nDenominator_ << " fake rate: " << (double)nNumerator_lt_/(double)nDenominator_ << "+-" <<  calculateFakeRateError(nNumerator_lt_,nDenominator_) << endl;
  stream << "Tight (weighted): num: " << nNumerator_lt_weighted_ << "+-" << sqrt(nNumerator_lt_weighted_w2_) << " den: " << nDenominator_weighted_ << "+-" << sqrt(nDenominator_weighted_w2_) << " fake rate: " << nNumerator_lt_weighted_/nDenominator_weighted_ << "+-" << calculateWeightedFakeRateError(nNumerator_lt_weighted_,nNumerator_lt_weighted_w2_,nDenominator_weighted_,nDenominator_weighted_w2_) << endl;

  stream << "Tight wo leading: num: " << nNumerator_lt_wo_leading_ << " den: " << nDenominator_wo_leading_ << " fake rate: " << (double)nNumerator_lt_wo_leading_/(double)nDenominator_wo_leading_ << "+-" << calculateFakeRateError(nNumerator_lt_wo_leading_,nDenominator_wo_leading_) << endl;
  stream << "Tight wo leading (weighted): num: " << nNumerator_lt_wo_leading_weighted_ << "+-" << sqrt(nNumerator_lt_wo_leading_weighted_w2_) << " den: " << nDenominator_wo_leading_weighted_ << "+-" << sqrt(nDenominator_wo_leading_weighted_w2_) << " fake rate: " << nNumerator_lt_wo_leading_weighted_/nDenominator_wo_leading_weighted_ << "+-" << calculateWeightedFakeRateError(nNumerator_lt_wo_leading_weighted_,nNumerator_lt_wo_leading_weighted_w2_,nDenominator_wo_leading_weighted_,nDenominator_wo_leading_weighted_w2_) << endl;

  stream << "Tight wo second leading: num: " << nNumerator_lt_wo_second_leading_ << " den: " << nDenominator_wo_second_leading_ << " fake rate: " << (double)nNumerator_lt_wo_second_leading_/(double)nDenominator_wo_second_leading_ << "+-" << calculateFakeRateError(nNumerator_lt_wo_second_leading_,nDenominator_wo_second_leading_) << endl;
  stream << "Tight wo second leading (weighted): num: " << nNumerator_lt_wo_second_leading_weighted_ << "+-" << sqrt(nNumerator_lt_wo_second_leading_weighted_w2_) << " den: " << nDenominator_wo_second_leading_weighted_ << "+-" << sqrt(nDenominator_wo_second_leading_weighted_w2_) << " fake rate: " << nNumerator_lt_wo_second_leading_weighted_/nDenominator_wo_second_leading_weighted_ << "+-" << calculateWeightedFakeRateError(nNumerator_lt_wo_second_leading_weighted_,nNumerator_lt_wo_second_leading_weighted_w2_,nDenominator_wo_second_leading_weighted_,nDenominator_wo_second_leading_weighted_w2_) << endl;
  stream << "=========" << endl << endl;

  cout << stream.str();

  int ret = fprintf(logfile_, 
		    stream.str().c_str());
  if (ret < 0)
    perror("writing to log file");
}

double Looper::calculateFakeRateError(unsigned int nNum, unsigned int nDen) {
  return sqrt( (double)nNum/((double)nDen*(double)nDen) + ((double)nNum*(double)nNum)/((double)nDen*(double)nDen*(double)nDen) );
}

double Looper::calculateWeightedFakeRateError(double nNum, double nNumErr2, double nDen, double nDenErr2) {
  return sqrt( (1/(nDen*nDen))*nNumErr2 + (nNum*nNum)/(nDen*nDen*nDen*nDen) * nDenErr2 );
}

