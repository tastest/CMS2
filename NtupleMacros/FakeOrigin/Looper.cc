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
  cands_passing_    = 0; 
  cands_passing_w2_ = 0; 
  cands_count_      = 0; 
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

  mcId_ = book1DHist(Form("%s_mcId",prefix),"mc id",3000,0.,3000.,"mc id","Events",color); 
  mcMotherId_ = book1DHist(Form("%s_mcMotherId",prefix),"mc mother id",3000,0.,3000.,"mc mother id","Events",color); 
  mcIdDenominator_ = book1DHist(Form("%s_mcIdDenominator",prefix),"mc id denominator",3000,0.,3000.,"mc id","Events",color); 
  mcMotherIdDenominator_ = book1DHist(Form("%s_mcMotherIdDenominator",prefix),"mc mother id denominator",3000,0.,3000.,"mc mother id","Events",color); 
  mcIdNumerator_ = book1DHist(Form("%s_mcIdNumerator",prefix),"mc id numerator",3000,0.,3000.,"mc id","Events",color); 
  mcMotherIdNumerator_ = book1DHist(Form("%s_mcMotherIdNumerator",prefix),"mc mother id numerator",3000,0.,3000.,"mc mother id","Events",color); 

}


bool Looper::FilterEvent()
{ 

  //
  // duplicate filter, based on trk information and dilepton hyp
  //
  // comment in following lines
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

     // electron, with pT>15 GeV and |eta| < 1
     bool foundElectron = false;
     for ( unsigned int ele = 0;
	   ele < cms2.els_p4().size();
	   ++ele ) {
       if ( cms2.els_p4()[ele].pt() > 15. ) {
	 if ( TMath::Abs(cms2.els_p4()[ele].eta()) <= 2.4 ) {
	   foundElectron = true;
	   break;
	 }
       }
     }
     if ( foundElectron ) {
       ret |= (CUT_BIT(CUT_ELE));
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
 
  // these are the cuts that the candidate passes: 
  cuts_t cuts_passed = EventSelect(); 
 
  // this is how to test that the candidate passes the cuts (which 
  // we specified in the constructor when we made the looper) 
  // (*note: the parentheses are important*): 
  if ((cuts_passed & cuts_) == cuts_) { 
    // if the candidate passed, we count it 
    cands_passing_ += weight; 
    cands_passing_w2_ += weight * weight; 
    cands_count_++; 
    // electron loop 
    for ( unsigned int ele = 0; 
          ele < cms2.els_p4().size(); 
          ++ele ) { 
      if ( cms2.els_p4()[ele].pt() > 15. && 
           TMath::Abs(cms2.els_p4()[ele].eta()) < 2.4 ) { 
        mcId_->Fill(TMath::Abs(cms2.els_mc_id()[ele]),weight); 
        mcMotherId_->Fill(TMath::Abs(cms2.els_mc_motherid()[ele]),weight); 
      }
      if ( isFakeable(ele) ) {
        mcIdDenominator_->Fill(TMath::Abs(cms2.els_mc_id()[ele]),weight); 
        mcMotherIdDenominator_->Fill(TMath::Abs(cms2.els_mc_motherid()[ele]),weight); 
      }
      if ( isNumeratorElectron(ele) ) {
        mcIdNumerator_->Fill(TMath::Abs(cms2.els_mc_id()[ele]),weight); 
        mcMotherIdNumerator_->Fill(TMath::Abs(cms2.els_mc_motherid()[ele]),weight); 
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

  int ret = fprintf(logfile_,  
                    "Sample %10s: Total candidate count: %8u %8u %8u %8u." 
                    " Total weight %10.1f +- %10.1f\n",    
                    sample_.name.c_str(), 
                    CandsCount(), 
                    CandsPassing()  , RMS()); 
     if (ret < 0)
	  perror("writing to log file");
}
