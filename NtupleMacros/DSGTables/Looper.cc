#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

ClassImp(DSGTable);
ClassImp(DSGTables);

const int DSGTable::nZcat 	;
const int DSGTable::nMETcat 	;
const int DSGTable::nSumJetcat  ;
const int DSGTable::nJetcat 	;
const int DSGTable::nBuckets 	;

static bool addDirectoryOff_ ()
{
     TH1::AddDirectory(false);
     return true;
}

const bool addDirectoryOff = addDirectoryOff_();

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname),
       dsgTable(s)
{
     // zero out the candidate counters (don't comment this out)
//      memset(cands_passing_	, 0, sizeof(cands_passing_       ));
//      memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
//      memset(cands_count_		, 0, sizeof(cands_count_         ));
     // hah, I commented it out!
}

void Looper::BookHistos ()
{
     //------------------------------------------------------------
     // Example histo booking; edit for your application
     //------------------------------------------------------------
}

bool Looper::FilterEvent()
{ 
     return false;
     //
     // duplicate filter, based on trk information and dilepton hyp
     //
     if (cms2.trks_d0().size() == 0)
	  return true;
     DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0()[0], 
				 cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), 
				 cms2.hyp_lt_p4()[0].phi() };
     if (is_duplicate(id)) {
	  duplicates_total_n_++;
	  duplicates_total_weight_ += cms2.evt_scale1fb();
// 	  cout << "Filtered duplicate run: " << cms2.evt_run() << " event: " << cms2.evt_event() << endl;
	  return true;
     }
     return false; 
}

cuts_t Looper::DilepSelect (int i_hyp)
{
     cuts_t ret = 0;
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

     // pt cuts
     if (cms2.hyp_lt_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LT_PT));
     if (cms2.hyp_ll_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LL_PT));
     // muon quality
     int n_iso_mu = 0;
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
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
     }
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO);
     }     
     return ret;
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

// 	  // if the candidate passed, we count it
// 	  cands_passing_[myType] += weight;
// 	  cands_passing_w2_[myType] += weight * weight;
// 	  cands_count_[myType]++;
// 	  cands_passing_[DILEPTON_ALL] += weight;
// 	  cands_passing_w2_[DILEPTON_ALL] += weight * weight;
// 	  cands_count_[DILEPTON_ALL]++;
//      }

	  const int zcat 	= Zcat(i_hyp);
	  const int metcat	= METcat(i_hyp);
	  const int jetcat	= Jetcat(i_hyp);
	  const int sumjetcat	= SumJetcat(i_hyp);
	  const int bucket	= Bucket(i_hyp);
	  
	  for (int i = 0; i <= metcat; ++i) {
	    for (int j = 0; j <= sumjetcat; ++j) {
	      dsgTable.Increment(zcat, i, j, jetcat, bucket, weight);
	      dsgTable.hmet_[zcat][i][j][jetcat][bucket]->Fill(cms2.evt_tcmet(), weight);
	      dsgTable.hmll_[zcat][i][j][jetcat][bucket]->Fill(cms2.hyp_p4()[i_hyp].M(), weight);
	    }
	    
	  }
     }
}     

int Looper::Zcat (int i_hyp) const
{
  if (additionalZcounter() > 0) return 1;
  /*
     // hypo in Z mass window
     if (abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_ll_id()[i_hyp])
	 && cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 &&
	 inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
	  return 1;
  */
     return 0;
}

int Looper::METcat (int) const
{
     if (cms2.evt_tcmet() > 175)
       return 2; 
    if (cms2.evt_tcmet() > 100)
      return 1;
     if (cms2.evt_tcmet()  > 35)
       return 0;
     return -1;
}


int Looper::SumJetcat (int i_hyp) const
{
  float SumJet=0;
  for (int i = 0; i < cms2.hyp_jets_p4().at(i_hyp).size(); ++i) {
    SumJet += cms2.hyp_jets_p4().at(i_hyp).at(i).pt(); }

  if (SumJet > 500)
    return 2; 
  if (SumJet > 300)
    return 1;
  if (SumJet > 0)
    return 0;
  return -1;
}

int Looper::Jetcat (int i_hyp) const
{
     switch (cms2.hyp_jets_p4().size()) {
     case 0: case 1: 
	  return 0;
     case 2: case 3: case 4:
	  return 1;
     default:
	  return 2;
     }
     assert(false);
}

int Looper::Bucket (int i_hyp) const
{
  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == -11) // e+e+
	  return 0;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == 11) // e-e-
	  return 1;

 if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == -13) // m+m+
	  return 2;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == 13) // m-m-
	  return 3;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == -13) // e+m+
	  return 4;
  if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == -11) // m+e+
	  return 4;
 
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == 13) // e-m-
	  return 5;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == 11) // m-e-
	  return 5;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == 13) // e+m-
	  return 6;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == -13) // e-m+
	  return 7;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == 11) // e+e-
	  return 8;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == -11) // e-e+
	  return 8;

  if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == 13) // m+m-
	  return 9;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == -13) // m-m+
	  return 9;
     assert(false);
}
