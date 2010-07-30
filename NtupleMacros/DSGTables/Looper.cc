#include <math.h>
#include "TVector3.h"
#include "CORE/electronSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/muonSelections.h"
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

Sample rename (const Sample &s, const std::string &name) 
{
     Sample ret = s;
     ret.name += name;
     return ret;
}

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname),
       dsgTable(s),
       dsgTable_emu(rename(s, "_emu"))
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
     DSGSearchWindow *highPtZ = new DSGSearchWindow(sample_, "highPtZ");
     highPtZ->SetH(new TH1F(Form("highPtZ%s", sample_.name.c_str()), "highPtZ", 10, 0, 500));
     highPtZ->dsgTable_ = new DSGTable(highPtZ->s_);
     dsgTable.search_windows_.push_back(highPtZ); 
     DSGSearchWindow *nonIsoHighPtMu = new DSGSearchWindow(sample_, "nonIsoHighPtMu");
     nonIsoHighPtMu->SetH(new TH1F(Form("nonIsoHighPtMu%s", sample_.name.c_str()), "nonIsoHighPtMu", 10, 0, 500));
     dsgTable.search_windows_.push_back(nonIsoHighPtMu);
}

bool Looper::FilterEvent()
{ 
     return false;
     //
     // duplicate filter, based on trk information and dilepton hyp
     //
     if (cms2.trks_d0().size() == 0)
	  return true;
     DorkyEventIdentifier id (cms2);
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

     cuts_t ele_lt_cuts = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11)
	  ele_lt_cuts = electronSelection(cms2.hyp_lt_index()[i_hyp]);
     cuts_t ele_ll_cuts = 0;
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11)
	  ele_ll_cuts = electronSelection(cms2.hyp_ll_index()[i_hyp]);
     cuts_t mask = electronSelection_ttbarV2;
     cuts_t mask_iso = cuts_t(1) << ELEISO_REL015;
     cuts_t mask_no_iso = mask & ~mask_iso;
     
     // pt cuts
     if (cms2.hyp_lt_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LT_PT));
     if (cms2.hyp_ll_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LL_PT));
     // muon quality
     int n_iso_mu = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && muonIdNotIsolated(cms2.hyp_lt_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LT_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && muonIdNotIsolated(cms2.hyp_ll_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LL_GOOD);
     // electron quality
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && (ele_lt_cuts & mask_no_iso) == mask_no_iso)
	  ret |= CUT_BIT(CUT_LT_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && (ele_ll_cuts & mask_no_iso) == mask_no_iso)
	  ret |= CUT_BIT(CUT_LL_GOOD);
     // calo iso
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && muonId(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && muonId(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
     }
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && (ele_lt_cuts & mask_iso) == mask_iso) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && (ele_ll_cuts & mask_iso) == mask_iso) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO);
     }     
     return ret;
}

void Looper::FillDSGTable (DSGTable &dsgTable, int i_hyp)
{
     //------------------------------------------------------------
     // Example dilepton histo filling; edit for your application
     //------------------------------------------------------------
     
     // everybody histogram needs to know what hypothesis he is 
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     // and what the event weight is 
     const double weight = Weight(i_hyp);

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
     
     float jsumet = 0;
     for (int i = 0; i < cms2.hyp_jets_p4()[i_hyp].size(); ++i) {
	  jsumet += cms2.hyp_jets_p4()[i_hyp][i].pt();
     }
     for (int i = 0; i <= metcat; ++i) {
	  for (int j = 0; j <= sumjetcat; ++j) {
	       dsgTable.Increment(zcat, i, j, jetcat, bucket, weight);
	       dsgTable.hmet_[zcat][i][j][jetcat][bucket]->Fill(cms2.evt_tcmet(), weight);
	       dsgTable.hmll_[zcat][i][j][jetcat][bucket]->Fill(cms2.hyp_p4()[i_hyp].M(), weight);
	       dsgTable.hht_[zcat][i][j][jetcat][bucket]
		    ->Fill( (cms2.evt_tcmet() +
			     cms2.hyp_lt_p4()[i_hyp].pt() + cms2.hyp_ll_p4()[i_hyp].pt() +
			     jsumet),
			    weight);
	       dsgTable.hjsumet_[zcat][i][j][jetcat][bucket]->Fill(jsumet, weight);
	       if (cms2.hyp_jets_p4()[i_hyp].size() == 0) 
		    dsgTable.hmaxjetpt_[zcat][i][j][jetcat][bucket]->Fill(0.0, weight);
	       else
		    dsgTable.hmaxjetpt_[zcat][i][j][jetcat][bucket]->Fill(cms2.hyp_jets_p4()[i_hyp][0].pt(), weight);
	      
	       dsgTable.hmaxleppt_[zcat][i][j][jetcat][bucket]
		    ->Fill(std::max(cms2.hyp_lt_p4()[i_hyp].pt(),cms2.hyp_ll_p4()[i_hyp].pt()) , weight);
	       dsgTable.hlepdphi_[zcat][i][j][jetcat][bucket]
		    ->Fill(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4()[i_hyp],cms2.hyp_ll_p4()[i_hyp]), weight);


	       // also fill the Z or no Z row
	       dsgTable.Increment(DSGTable::nZcat, i, j, jetcat, bucket, weight);
	       dsgTable.hmet_[DSGTable::nZcat][i][j][jetcat][bucket]->Fill(cms2.evt_tcmet(), weight);
	       dsgTable.hmll_[DSGTable::nZcat][i][j][jetcat][bucket]->Fill(cms2.hyp_p4()[i_hyp].M(), weight);    
	       dsgTable.hht_[DSGTable::nZcat][i][j][jetcat][bucket]
		    ->Fill( (cms2.evt_tcmet() +
			     cms2.hyp_lt_p4()[i_hyp].pt() + cms2.hyp_ll_p4()[i_hyp].pt() +
			     jsumet),
			    weight);
	       dsgTable.hjsumet_[DSGTable::nZcat][i][j][jetcat][bucket]->Fill(jsumet, weight);
	       if (cms2.hyp_jets_p4()[i_hyp].size() == 0) 
		    dsgTable.hmaxjetpt_[DSGTable::nZcat][i][j][jetcat][bucket]->Fill(0.0, weight);
	       else
		    dsgTable.hmaxjetpt_[DSGTable::nZcat][i][j][jetcat][bucket]->Fill(cms2.hyp_jets_p4()[i_hyp][0].pt(), weight);
	      
	       dsgTable.hmaxleppt_[DSGTable::nZcat][i][j][jetcat][bucket]
		    ->Fill(std::max(cms2.hyp_lt_p4()[i_hyp].pt(),cms2.hyp_ll_p4()[i_hyp].pt()) , weight);
	       dsgTable.hlepdphi_[DSGTable::nZcat][i][j][jetcat][bucket]
		    ->Fill(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4()[i_hyp],cms2.hyp_ll_p4()[i_hyp]), weight);
	  }
	    
     }
}

void Looper::FillDilepHistos (int i_hyp)
{
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     const double weight = Weight(i_hyp);
     // these are the cuts that the candidate passes:
     cuts_t cuts_passed = DilepSelect(i_hyp);
     // this is how to test that the candidate passes the cuts (which
     // we specified in the constructor when we made the looper)
     // (*note: the parentheses are important*):
     if ((cuts_passed & cuts_) == cuts_) {
	  FillDSGTable(dsgTable, i_hyp);
     }
     const cuts_t cuts_no_iso = cuts_ & ~(CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO));
     if ((cuts_passed & cuts_no_iso) == cuts_no_iso) {
	  // even if we fail iso, we still look for a Z
	  if ((myType == DILEPTON_EE || myType == DILEPTON_MUMU) && // ee or mu mu 
	      cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 && // opposite charge
	      inZmassWindow(cms2.hyp_p4()[i_hyp].M()) && // in Z mass window
	      cms2.hyp_p4()[i_hyp].pt() > 100 &&
	      (cms2.hyp_lt_ecaliso()[i_hyp] + cms2.hyp_lt_trkiso()[i_hyp]) / cms2.hyp_lt_p4()[i_hyp].pt() < 0.1 &&
	      (cms2.hyp_ll_ecaliso()[i_hyp] + cms2.hyp_ll_trkiso()[i_hyp]) / cms2.hyp_ll_p4()[i_hyp].pt() < 0.1) {
	       DSGSearchWindow *highPtZ = dsgTable.search_windows_[0];
	       highPtZ->Increment(cms2.hyp_p4()[i_hyp].pt(), weight);
	       // also fill the DSGTable in the search window
	       FillDSGTable(*highPtZ->dsgTable_, i_hyp);
	  }
     }
}     

int Looper::Zcat (int i_hyp) const
{
     return 0;
#if 0
  if (additionalZcounter() > 0) return 1;
  /*
     // hypo in Z mass window
     if (abs(cms2.hyp_lt_id()[i_hyp]) == abs(cms2.hyp_ll_id()[i_hyp])
	 && cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 &&
	 inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
	  return 1;
  */
     return 0;
#endif
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
	  return EPEP;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == 11) // e-e-
	  return EMEM;

 if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == -13) // m+m+
	  return MPMP;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == 13) // m-m-
	  return MMMM;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == -13) // e+m+
	  return EPMP;
  if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == -11) // m+e+
	  return EPMP;
 
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == 13) // e-m-
	  return EMMM;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == 11) // m-e-
	  return EMMM;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == 13) // e+m-
	  return EPMM;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == -13) // e-m+
	  return EMMP;

  if (cms2.hyp_lt_id()[i_hyp] == -11 && cms2.hyp_ll_id()[i_hyp] == 11) // e+e-
	  return EPEM;
  if (cms2.hyp_lt_id()[i_hyp] ==  11 && cms2.hyp_ll_id()[i_hyp] == -11) // e-e+
	  return EPEM;

  if (cms2.hyp_lt_id()[i_hyp] == -13 && cms2.hyp_ll_id()[i_hyp] == 13) // m+m-
	  return MPMM;
  if (cms2.hyp_lt_id()[i_hyp] ==  13 && cms2.hyp_ll_id()[i_hyp] == -13) // m-m+
	  return MPMM;
  assert(false);
}

void overflow (TH1F* h)
{
  float overflowBin = h->GetBinContent(h->GetNbinsX()+1);
  float lastBin = h->GetBinContent(h->GetNbinsX());
  h->SetBinContent(h->GetNbinsX(),lastBin+overflowBin);
  // consider zeroing out the ovflow...  h->SetBinContent(h->GetNbinsX()+1,0.0);
}

void Looper::overflows ()
{
  for (int i = 0; i <= DSGTable::nZcat; ++i) {
    for (int j = 0; j < DSGTable::nMETcat; ++j) {
      for (int jj = 0; jj < DSGTable::nSumJetcat; ++jj) {
	for (int k = 0; k < DSGTable::nJetcat; ++k) {
	  for (int l = 0; l < DSGTable::nBuckets; ++l) {
	    overflow(dsgTable.hmet_[i][j][jj][k][l]);
	    overflow(dsgTable.hmll_[i][j][jj][k][l]);
	    overflow(dsgTable.hht_[i][j][jj][k][l]);
	    overflow(dsgTable.hjsumet_[i][j][jj][k][l]);
	    overflow(dsgTable.hmaxjetpt_[i][j][jj][k][l]);
	    overflow(dsgTable.hmaxleppt_[i][j][jj][k][l]);
	    overflow(dsgTable.hlepdphi_[i][j][jj][k][l]);
	  }
	}
      }
    }
  }
}

void Looper::emu_subtraction()
{
     static const int subtractions[][2] = { 
	  { EPMP, -1 },
	  { EMMM, -1 },
	  { EPMP, -1 },
	  { EMMM, -1 },
	  { -1, -1, },
	  { -1, -1, },
	  { -1, -1, },
	  { -1, -1, },
	  { EPMM, EMMP },
	  { EPMM, EMMP },
     };
     const static DSGTable::table_t DSGTable::* histos[] = {
	  &DSGTable::hmet_    ,
	  &DSGTable::hmll_    ,
	  &DSGTable::hht_     ,
	  &DSGTable::hjsumet_ ,
	  &DSGTable::hmaxjetpt_,
	  &DSGTable::hmaxleppt_,
	  &DSGTable::hlepdphi_,
     };
     static const int n_histos = sizeof(histos) / sizeof(DSGTable::table_t DSGTable::*);
     for (int i = 0; i <= DSGTable::nZcat; ++i) {
	  for (int j = 0; j < DSGTable::nMETcat; ++j) {
	       for (int jj = 0; jj < DSGTable::nSumJetcat; ++jj) {
		    for (int k = 0; k < DSGTable::nJetcat; ++k) {
			 for (int i_bucket = 0; i_bucket < DSGTable::nBuckets; ++i_bucket) {
			      int sub1 = subtractions[i_bucket][0];
			      int sub2 = subtractions[i_bucket][1];
			      double s1 = 0;
			      double s2 = 0;
			      double s1w2 = 0;
			      double s2w2 = 0;
			      if (sub1 >= 0) {
				   s1 = dsgTable.events_[i][j][jj][k][sub1];
				   s1w2 = dsgTable.w2s_[i][j][jj][k][sub1];
			      }
			      if (sub2 >= 0) {
				   s2 = dsgTable.events_[i][j][jj][k][sub2];
				   s2w2 = dsgTable.w2s_[i][j][jj][k][sub2];
			      }
			      dsgTable_emu.events_[i][j][jj][k][i_bucket] = dsgTable.events_[i][j][jj][k][i_bucket] - 0.5 * (s1 + s2);
			      dsgTable_emu.w2s_[i][j][jj][k][i_bucket] = dsgTable.w2s_[i][j][jj][k][i_bucket] + 0.25 * (s1w2 + s2w2);
			      for (int i_histo = 0; i_histo < n_histos; ++i_histo) {
				   TH1F *h_emu = (dsgTable_emu.*(histos[i_histo]))[i][j][jj][k][i_bucket];
				   TH1F *h = (dsgTable.*(histos[i_histo]))[i][j][jj][k][i_bucket];
				   *h_emu = *h;
				   if (sub1 >= 0) {
					TH1F *h_sub1 = (dsgTable.*(histos[i_histo]))[i][j][jj][k][sub1];
					h_emu->Add(h_sub1, -0.5);
				   }
				   if (sub2 >= 0) {
					TH1F *h_sub2 = (dsgTable.*(histos[i_histo]))[i][j][jj][k][sub2];
					h_emu->Add(h_sub2, -0.5);
				   }
			      }
			 }
		    }
	       }
	  }
     }
}

void Looper::End()
{
     overflows();
     emu_subtraction();
}

