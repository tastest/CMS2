#include <math.h>
#include <sstream>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

MuonLooper::MuonLooper (Sample s, cuts_t c, const char *fname) 
  : Looper(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  events_            = 0;
  events_weighted_   = 0;
  events_passing_    = 0;
  events_passing_w2_ = 0;
  events_count_	     = 0;

  nDenominator_ = 0;
  nNumerator_ = 0;

  nDenominator_weighted_ = 0;
  nNumerator_weighted_ = 0;

  nDenominator_weighted_w2_ = 0;
  nNumerator_weighted_w2_ = 0;
}

void MuonLooper::BookHistos ()
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

  // set Pt bins
  binsPt.push_back(10.);
  binsPt.push_back(20.);
  binsPt.push_back(60.);
  binsPt.push_back(150.);

  // set eta bins
  binsEta.push_back(-2.4); //changed 2.5 -> 2.4 090722
  binsEta.push_back(-1.5);
  binsEta.push_back(1.5);
  binsEta.push_back(2.4);


  pt_num_mu_  	= book1DVarHist(Form("%s_hpt_num_mu",prefix),"p_{T} Numerator",binsPt,"p_{T} [GeV]","p_{T} Muon Numerator",color);
  eta_num_mu_ 	= book1DVarHist(Form("%s_heta_num_mu",prefix),"#eta Numerator",binsEta,"#eta","#eta Muon Numerator",color);
  num_mu_ 	= book2DVarHist(Form("%s_hnum_mu",prefix),"Numerator",binsEta,binsPt,"#eta","p_{T} [GeV]","Muon Numerator",color);

  pt_den_mu_  	= book1DVarHist(Form("%s_hpt_den_mu",prefix),"p_{T} Denominator",binsPt,"p_{T} [GeV]","p_{T} Muon Denominator",color);
  eta_den_mu_ 	= book1DVarHist(Form("%s_heta_den_mu",prefix),"#eta Denominator",binsEta,"#eta","#eta Muon Denominator",color);
  den_mu_ 	= book2DVarHist(Form("%s_hden_mu",prefix),"Denominator",binsEta,binsPt,"#eta","p_{T} [GeV]","Muon Denominator",color);
	
  mu_n_   	= book1DHist(Form("%s_hmu_n",prefix),"number of muons per event",15,0.,15.,"n_{e}","Events",color);
}

cuts_t MuonLooper::EventSelect ()
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

void MuonLooper::FillEventHistos ()
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
    mu_n_->Fill(cms2.mus_p4().size(),weight);

//     njets_->Fill(cms2.evt_njets(),weight);
//     for ( unsigned int jet_counter = 0;
// 	  jet_counter < (unsigned int)cms2.evt_njets();
// 	  ++jet_counter ) {
//       jetpt_->Fill(cms2.jets_p4()[jet_counter].Pt(),weight);
//       jeteta_->Fill(cms2.jets_p4()[jet_counter].Eta(),weight);
//     }

    // cuts
    float deltaRCut = 0.2;
    bool  useAbsEta = true;

    // muon loop
    for ( unsigned int muon_counter = 0;
	  muon_counter < cms2.mus_p4().size();
	  ++muon_counter ) {

      // general electron histograms
//       tkIso_uncut_->Fill(cms2.els_tkIso()[electron_counter],weight);
//       EOverp_uncut_->Fill(cms2.els_eOverPIn()[electron_counter],weight);
//       HOverE_uncut_->Fill(cms2.els_hOverE()[electron_counter],weight);

      // fill every electron with pt > 150 GeV into last bin
      float pt = cms2.mus_p4()[muon_counter].Pt();
      if (pt >= 150.) pt = 149.;

      if (pt < 10 || fabs(cms2.mus_p4()[muon_counter].Eta()) > 2.4)
	   continue;

      if ( isFakeableMuon(muon_counter)){
	  
	   pt_den_mu_->Fill(pt,weight);
	if ( useAbsEta ) {
	  eta_den_mu_->Fill(TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),weight);
	  eta_den_mu_->Fill(-TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),weight);
	  den_mu_->Fill(TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),pt,weight);
	  den_mu_->Fill(-TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),pt,weight);
	} else {
	  eta_den_mu_->Fill(cms2.mus_p4()[muon_counter].Eta(),weight);
	  den_mu_->Fill(cms2.mus_p4()[muon_counter].Eta(),pt,weight);
	}

	// hard cut to calculate zero-bin fake rate correctly
	if ( pt >= 20. ) {
	  ++nDenominator_;
	  nDenominator_weighted_ += weight;
	  nDenominator_weighted_w2_ += weight * weight;
	}

// 	for ( unsigned int jet_counter = 0;
// 	      jet_counter < (unsigned int)cms2.evt_njets();
// 	      ++jet_counter ) {
// 	  deltaR_->Fill(dRbetweenVectors(cms2.mus_p4()[muon_counter],cms2.jets_p4()[jet_counter]),weight);
// 	}
// 	tkIso_->Fill(cms2.mus_tkIso()[muon_counter],weight);
// 	EOverp_->Fill(cms2.mus_eOverPIn()[muon_counter],weight);
// 	HOverE_->Fill(cms2.mus_hOverE()[muon_counter],weight);
      }

      if (isNumeratorMuon(muon_counter) && pt > 10 && fabs(cms2.mus_p4()[muon_counter].Eta()) < 2.4){ 
	if ( !isFakeableMuon(muon_counter) ) cout << "Loose muon: num is fullfilled but den is not!" << endl;
	pt_num_mu_->Fill(pt,weight);
	if ( useAbsEta ) {
	  eta_num_mu_->Fill(TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),weight);
	  eta_num_mu_->Fill(-TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),weight);
	  num_mu_->Fill(TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),pt,weight);
	  num_mu_->Fill(-TMath::Abs(cms2.mus_p4()[muon_counter].Eta()),pt,weight);
	} else {
	  eta_num_mu_->Fill(cms2.mus_p4()[muon_counter].Eta(),weight);
	  num_mu_->Fill(cms2.mus_p4()[muon_counter].Eta(),pt,weight);
	}

	// hard cut to calculate zero-bin fake rate correctly
	if ( pt >= 20. ) {
	  ++nNumerator_;
	  nNumerator_weighted_ += weight;
	  nNumerator_weighted_w2_ += weight * weight;
	}
      }
    }
  }
}

void MuonLooper::End ()
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

  stream << "Tight: num: " << nNumerator_ << " den: " << nDenominator_ << " fake rate: " << (double)nNumerator_/(double)nDenominator_ << "+-" <<  calculateFakeRateError(nNumerator_,nDenominator_) << endl;
  stream << "Tight (weighted): num: " << nNumerator_weighted_ << "+-" << sqrt(nNumerator_weighted_w2_) << " den: " << nDenominator_weighted_ << "+-" << sqrt(nDenominator_weighted_w2_) << " fake rate: " << nNumerator_weighted_/nDenominator_weighted_ << "+-" << calculateWeightedFakeRateError(nNumerator_weighted_,nNumerator_weighted_w2_,nDenominator_weighted_,nDenominator_weighted_w2_) << endl;

  stream << "=========" << endl << endl;

  cout << stream.str();

  int ret = fprintf(logfile_, 
		    stream.str().c_str());
  if (ret < 0)
    perror("writing to log file");
}

double MuonLooper::calculateFakeRateError(unsigned int nNum, unsigned int nDen) {
  return sqrt( (double)nNum/((double)nDen*(double)nDen) + ((double)nNum*(double)nNum)/((double)nDen*(double)nDen*(double)nDen) );
}

double MuonLooper::calculateWeightedFakeRateError(double nNum, double nNumErr2, double nDen, double nDenErr2) {
  return sqrt( (1/(nDen*nDen))*nNumErr2 + (nNum*nNum)/(nDen*nDen*nDen*nDen) * nDenErr2 );
}

