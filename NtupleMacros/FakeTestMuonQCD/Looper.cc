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
  cands_passing_    = 0; 
  cands_passing_w2_ = 0; 
  cands_count_      = 0; 
  events_           = 0;

}

void Looper::BookHistos ()
{

  gDirectory = histo_directory;

  const unsigned int ptNBins = 16;
  float ptBins[ptNBins+1];
  for ( unsigned int ptBin = 0;
	ptBin <= ptNBins;
	++ptBin) {
    ptBins[ptBin] = float(ptBin)*160./16.;
  }

  const unsigned int etaNBins = 12;
  float etaBins[etaNBins+1];
  for ( unsigned int etaBin = 0;
	etaBin <= etaNBins;
	++etaBin) {
    etaBins[etaBin] = float(etaBin)/2.-3.;
  }

  hmuPt_ = book1DVarHist(Form("%s_%s",sample_.name.c_str(),"muPt"),
			 Form("%s_%s",sample_.name.c_str(),"muPt"),
			 ptNBins,ptBins,
			 "p_{T}^{e} [GeV]","Events",sample_.histo_color);
  hmuEta_ = book1DVarHist(Form("%s_%s",sample_.name.c_str(),"muEta"),
			  Form("%s_%s",sample_.name.c_str(),"muEta"),
			  etaNBins,etaBins,
			  "#eta^{e} [GeV]","Events",sample_.histo_color);
}

bool Looper::FilterEvent()
{ 
  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  if (cms2.els_p4().size() > 0 ) { 
    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.els_d0corr()[0],  
                                cms2.els_p4()[0].pt(), cms2.els_p4()[0].eta(), cms2.els_p4()[0].phi() }; 
    if (is_duplicate(id)) { 
      duplicates_total_n_++; 
      duplicates_total_weight_ += cms2.evt_scale1fb(); 
      return true; 
    } 
  } else if ( cms2.mus_p4().size() > 0 ) { 
    DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.mus_d0corr()[0],  
                                cms2.mus_p4()[0].pt(), cms2.mus_p4()[0].eta(), cms2.mus_p4()[0].phi() }; 
    if (is_duplicate(id)) { 
      duplicates_total_n_++; 
      duplicates_total_weight_ += cms2.evt_scale1fb(); 
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

  // check if there is at least one muon fulfilling
  // pt cut, ID cut, iso cut, denominator, numerator
  for ( unsigned int mu = 0;
	mu < cms2.mus_p4().size();
	++mu ) {
    if (isFakeableMuon(mu)) {
      ret |= CUT_BIT(CUT_MU_DENOMINATOR);
    }
    if (isNumeratorMuon(mu)) {
      ret |= CUT_BIT(CUT_MU_NUMERATOR);
    }
  }

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
  cuts_t ret = 0;

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

  ++events_;

  // these are the cuts that the candidate passes: 
  cuts_t cuts_passed = EventSelect(); 

  // this is how to test that the candidate passes the cuts (which 
  // we specified in the constructor when we made the looper) 
  // (*note: the parentheses are important*): 
  if ((cuts_passed & cuts_) == cuts_) { 
    // muon loop 
    for ( unsigned int mu = 0; 
          mu < cms2.mus_p4().size(); 
          ++mu ) { 

      bool fill = false;
      // decide whether in observation or prediction mode
      // fill if mu is numerator or denominator muon
      if ( (cuts_ & CUT_BIT(CUT_MU_DENOMINATOR)) == CUT_BIT(CUT_MU_DENOMINATOR) ) {
	// prediction mode
	if ( isFakeableMuon(mu) && !isNumeratorMuon(mu,2)) {
	  fill = true;
	}
      } else {
	// observation mode
	if ( isNumeratorMuon(mu,2)) {
	  fill = true;
	}
      }

      if ( fill == true ) {
	// and what the weight is muon dependent in the case of the denominator, numerator weight is the event weight
	const double weight = Weight(mu); 
	// if the candidate passed, we count it 
	cands_passing_ += weight; 
	cands_passing_w2_ += weight * weight; 
	cands_count_++;
	// fill histos
	hmuPt_->Fill(cms2.mus_p4()[mu].pt(), weight);
	hmuEta_->Fill(cms2.mus_p4()[mu].eta(), weight);
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
  stream << "Total candidate count : " << CandsCount() << endl;
  stream << "Total weight: " << fixed << setprecision(1) << CandsPassing() << "+-" << RMS() << endl;
  stream << "=========" << endl;
  
  cout << stream.str();

  int ret = fprintf(logfile_, stream.str().c_str());
  if (ret < 0)
    perror("writing to log file");
}
