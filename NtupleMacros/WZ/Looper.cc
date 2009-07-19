#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname) 
  : LooperBase(s, c, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
  //------------------------------------------------------------
  // Example histo booking; edit for your application
  //------------------------------------------------------------

  // or use the N - 1 technology (see NMinus1Hist.h)
  // arguments are as follows: sample, name, binning, required cuts, cuts that are relaxed for the N - 1 plot
  // for the lt N - 1 plot, we relax the lt pt requirement
  htcmet_		= new TrilepNMinus1Hist(sample_, "tcMET"   ,	 40, 0, 200, cuts_,CUT_BIT(CUT_TCMET_15));
  hptcmet_		= new TrilepNMinus1Hist(sample_, "tcMETP"   ,	 40, 0, 200, cuts_,CUT_BIT(CUT_TCMET_15));
  h_highest_lep_pt_	= new TrilepNMinus1Hist(sample_, "PtHighestLeptonPt",	 50, 0, 100, cuts_,  (CUT_BIT(CUT_HIGHEST_PT_LEP_PT20)) |
						(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_PT20)) |
						(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_PT20)));
  h_second_highest_lep_pt_	= new TrilepNMinus1Hist(sample_, "PtSecondHighestLeptonPt",	 50, 0, 100, cuts_,  (CUT_BIT(CUT_HIGHEST_PT_LEP_PT20)) |
							(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_PT20)) |
							(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_PT20)));
  h_third_highest_lep_pt_	= new TrilepNMinus1Hist(sample_, "PtThirdHighestLeptonPt",	 50, 0, 100, cuts_,  (CUT_BIT(CUT_HIGHEST_PT_LEP_PT20)) |
							(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_PT20)) |
							(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_PT20)));

  h_highest_lep_iso_	= new TrilepNMinus1Hist(sample_, "PtHighestLeptonIso",	 100,0.,1., cuts_,CUT_BIT(CUT_HIGHEST_PT_LEP_ISO));
  h_second_highest_lep_iso_	= new TrilepNMinus1Hist(sample_, "PtSecondHighestLeptonIso",	 100,0.,1., cuts_,CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO));
  h_third_highest_lep_iso_	= new TrilepNMinus1Hist(sample_, "PtThirdHighestLeptonIso",	 100,0.,1., cuts_,CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO));


  h_highest_iso_lep_iso_	= new TrilepNMinus1Hist(sample_, "IsoHighestLeptonIso",	 100,0.,1., cuts_, (CUT_BIT(CUT_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO)));
  h_second_highest_iso_lep_iso_	= new TrilepNMinus1Hist(sample_, "IsoSecondHighestLeptonIso",	 100,0.,1., cuts_, (CUT_BIT(CUT_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO)));
  h_third_highest_iso_lep_iso_	= new TrilepNMinus1Hist(sample_, "IsoThirdHighestLeptonIso",	 100,0.,1., cuts_, (CUT_BIT(CUT_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO)) |
							(CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO)));

  h_counter_electrons_ =  new TrilepNMinus1Hist(sample_, "counterElectrons",	 10,-0.5,9.5, cuts_, CUT_BIT(CUT_ADD_ELECTRONS_VETO_CUT));
  h_counter_muons_ =  new TrilepNMinus1Hist(sample_, "counterMuons",	 10,-0.5,9.5, cuts_, CUT_BIT(CUT_ADD_MUONS_VETO_CUT));

  h_njets_ = new TrilepNMinus1Hist(sample_, "nJets",	 10,-0.5,9.5, cuts_, 0);

  h_DeltaPhiMETNearestLepton_ = new TrilepNMinus1Hist(sample_, "DeltaPhiMETNearestLepton",	 20,0,TMath::Pi(), cuts_, CUT_BIT(CUT_TCMET_15));

  h_primZMass_ = new TrilepNMinus1Hist(sample_, "ZMassPrim",	 200,0,200, cuts_, CUT_BIT(CUT_PRIM_Z) | CUT_BIT(CUT_NO_SEC_Z));
  h_addZMass_ = new TrilepNMinus1Hist(sample_, "ZMassAdd",	 200,0,200, cuts_, CUT_BIT(CUT_NO_SEC_Z));
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

  cuts_t ret = 0;
  return ret;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would make your cuts here
  //------------------------------------------------------------

  // cuts are failed until proven otherwise
  cuts_t ret = 0;

  trileptonPt_.clear();
  trileptonIso_.clear();

  // check that trilepton types can only be 1 or 2
  assert(abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 || abs(cms2.hyp_trilep_first_type()[i_hyp]) == 2);
  assert(abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 || abs(cms2.hyp_trilep_second_type()[i_hyp]) == 2);
  assert(abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1 || abs(cms2.hyp_trilep_third_type()[i_hyp]) == 2);

  // first lepton has pt >= 20 GeV
  if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 ) {
    if ( cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_FIRSTLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt(),mu_rel_iso(cms2.hyp_trilep_first_index()[i_hyp])));
    trileptonIso_.insert(std::pair<float,float>(mu_rel_iso(cms2.hyp_trilep_first_index()[i_hyp]),cms2.mus_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt()));
  } else {
    if ( cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_FIRSTLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt(),el_rel_iso(cms2.hyp_trilep_first_index()[i_hyp],true)));
    trileptonIso_.insert(std::pair<float,float>(el_rel_iso(cms2.hyp_trilep_first_index()[i_hyp],true),cms2.els_p4()[cms2.hyp_trilep_first_index()[i_hyp]].Pt()));
  }

  // second lepton has pt >= 20 GeV
  if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
    if ( cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_SECONDLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt(),mu_rel_iso(cms2.hyp_trilep_second_index()[i_hyp])));
    trileptonIso_.insert(std::pair<float,float>(mu_rel_iso(cms2.hyp_trilep_second_index()[i_hyp]),cms2.mus_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt()));
  } else {
    if ( cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_SECONDLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt(),el_rel_iso(cms2.hyp_trilep_second_index()[i_hyp],true)));
    trileptonIso_.insert(std::pair<float,float>(el_rel_iso(cms2.hyp_trilep_second_index()[i_hyp],true),cms2.els_p4()[cms2.hyp_trilep_second_index()[i_hyp]].Pt()));
  }

  // third lepton has pt >= 20 GeV
  if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1 ) {
    if ( cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_THIRDLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt(),mu_rel_iso(cms2.hyp_trilep_third_index()[i_hyp])));
    trileptonIso_.insert(std::pair<float,float>(mu_rel_iso(cms2.hyp_trilep_third_index()[i_hyp]),cms2.mus_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt()));
  } else {
    if ( cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt() >= 20. )
      ret |= CUT_BIT(CUT_THIRDLEP_PT20);
    trileptonPt_.insert(std::pair<float,float>(cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt(),el_rel_iso(cms2.hyp_trilep_third_index()[i_hyp],true)));
    trileptonIso_.insert(std::pair<float,float>(el_rel_iso(cms2.hyp_trilep_third_index()[i_hyp],true),cms2.els_p4()[cms2.hyp_trilep_third_index()[i_hyp]].Pt()));
  }

  // all leptons >= 20 GeV
  if ( (CUT_BIT(CUT_FIRSTLEP_PT20) & ret) &&
       (CUT_BIT(CUT_SECONDLEP_PT20) & ret) &&
       (CUT_BIT(CUT_THIRDLEP_PT20) & ret ) ) {
    ret |= CUT_BIT(CUT_ALLLEP_PT20);
  }
     
  // muon quality
  if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 ) {
    if ( goodMuonWithoutIsolation(cms2.hyp_trilep_first_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_FIRSTLEP_GOOD);
    if ( passMuonIsolation(cms2.hyp_trilep_first_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_FIRSTLEP_ISO);
  }
  if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 ) {
    if ( goodMuonWithoutIsolation(cms2.hyp_trilep_second_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_SECONDLEP_GOOD);
    if ( passMuonIsolation(cms2.hyp_trilep_second_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_SECONDLEP_ISO);
  }
  if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1 ) {
    if ( goodMuonWithoutIsolation(cms2.hyp_trilep_third_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_THIRDLEP_GOOD);
    if ( passMuonIsolation(cms2.hyp_trilep_third_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_THIRDLEP_ISO);
  }
  // electron quality
  if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 2 ) {
    if ( goodElectronWithoutIsolation(cms2.hyp_trilep_first_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_FIRSTLEP_GOOD);
    if ( passElectronIsolation(cms2.hyp_trilep_first_index()[i_hyp],true) )
      ret |= CUT_BIT(CUT_FIRSTLEP_ISO);
  }
  if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 2 ) {
    if ( goodElectronWithoutIsolation(cms2.hyp_trilep_second_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_SECONDLEP_GOOD);
    if ( passElectronIsolation(cms2.hyp_trilep_second_index()[i_hyp],true) )
      ret |= CUT_BIT(CUT_SECONDLEP_ISO);
  }
  if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 2 ) {
    if ( goodElectronWithoutIsolation(cms2.hyp_trilep_third_index()[i_hyp]) )
      ret |= CUT_BIT(CUT_THIRDLEP_GOOD);
    if ( passElectronIsolation(cms2.hyp_trilep_third_index()[i_hyp],true) )
      ret |= CUT_BIT(CUT_THIRDLEP_ISO);
  }

  // all lepton quality
  if ( (CUT_BIT(CUT_FIRSTLEP_GOOD) & ret) &&
       (CUT_BIT(CUT_SECONDLEP_GOOD) & ret) &&
       (CUT_BIT(CUT_THIRDLEP_GOOD) & ret ) ) {
    ret |= CUT_BIT(CUT_ALLLEP_GOOD);
  }
  if ( (CUT_BIT(CUT_FIRSTLEP_ISO) & ret) &&
       (CUT_BIT(CUT_SECONDLEP_ISO) & ret) &&
       (CUT_BIT(CUT_THIRDLEP_ISO) & ret ) ) {
    ret |= CUT_BIT(CUT_ALLLEP_ISO);
  }

  // determine cut flags for leptons ordered in pt
  assert(trileptonPt_.size() == 3);
  std::multimap<float,float,std::greater<float> >::iterator entry = trileptonPt_.begin();  

  if ( entry->first >= 20.)
    ret |= CUT_BIT(CUT_HIGHEST_PT_LEP_PT20);
  if ( entry->second > 0.9 )
    ret |= CUT_BIT(CUT_HIGHEST_PT_LEP_ISO);
  ++entry;
  if ( entry->first >= 10.)
    ret |= CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_PT20);
  if ( entry->second > 0.9 )
    ret |= CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO);
  ++entry;
  if ( entry->first >= 10.)
    ret |= CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_PT20);
  if ( entry->second > 0.9 )
    ret |= CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO);
    
  // primary Z
  primZMass_ = 0.;
  int notUsedLepton = findPrimTrilepZ(i_hyp,primZMass_);
  assert(notUsedLepton != 900);

  if ( notUsedLepton <= 3 )
    ret |= CUT_BIT(CUT_PRIM_Z);

  // additional Z veto (lepton + high pt isolated track (pt > 20 GeV)
  addZMass_ = -1;
  if ( notUsedLepton != 999 )
    if ( !vetoAddZ(i_hyp,notUsedLepton,addZMass_) )
      ret |= CUT_BIT(CUT_NO_SEC_Z);     

  // Z veto using additional leptons in the event

  // MET >= 15
  if ( cms2.evt_tcmet() >= 15. ) 
    ret |= CUT_BIT(CUT_TCMET_15);

  // PMET >= 15
  if ( MetSpecialTrilep(cms2.evt_tcmet(), cms2.evt_tcmetPhi(), i_hyp) >= 15. ) 
    ret |= CUT_BIT(CUT_PTCMET_15);

  // fourth lepton veto
  addElectronsCounter_ = 0;
  for ( unsigned int electron = 0;
	electron < cms2.els_p4().size();
	++electron ) {
    if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 2 && cms2.hyp_trilep_first_index()[i_hyp] == electron ) continue;
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 2 && cms2.hyp_trilep_second_index()[i_hyp] == electron ) continue;
    if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 2 && cms2.hyp_trilep_third_index()[i_hyp] == electron ) continue;
    if ( goodElectronWithoutIsolation(electron) ) ++addElectronsCounter_;
  }

  // cut on more than 1 additional good electron
  if ( addElectronsCounter_ < 1 )
    ret |= CUT_BIT(CUT_ADD_ELECTRONS_VETO_CUT);

  addMuonsCounter_ = 0;
  for ( unsigned int muon = 0;
	muon < cms2.mus_p4().size();
	++muon ) {
    if ( abs(cms2.hyp_trilep_first_type()[i_hyp]) == 1 && cms2.hyp_trilep_first_index()[i_hyp] == muon ) continue;
    if ( abs(cms2.hyp_trilep_second_type()[i_hyp]) == 1 && cms2.hyp_trilep_second_index()[i_hyp] == muon ) continue;
    if ( abs(cms2.hyp_trilep_third_type()[i_hyp]) == 1 && cms2.hyp_trilep_third_index()[i_hyp] == muon ) continue;
    if ( goodMuonWithoutIsolation(muon) ) ++addMuonsCounter_;
  }

  // cut on more than 1 additional good muon
  if ( addMuonsCounter_ < 1 )
    ret |= CUT_BIT(CUT_ADD_MUONS_VETO_CUT);

  // the return value gets cached, too
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

}

void Looper::FillTrilepHistos (int i_hyp)
{
  //------------------------------------------------------------
  // In a trilepton analysis, you would fill your histos here
  //------------------------------------------------------------

  // event weight
  const double weight = Weight(i_hyp);

  // these are the cuts that the candidate passes:
  cuts_t cuts_passed = TrilepSelect(i_hyp);

  // this is how to test that the candidate passes the cuts (which
  // we specified in the constructor when we made the looper)
  // (*note: the parentheses are important*):
  if ((cuts_passed & cuts_) == cuts_) {

    // if the candidate passed, we count it
    cands_passing_[cms2.hyp_trilep_bucket()[i_hyp]] += weight;
    cands_passing_w2_[cms2.hyp_trilep_bucket()[i_hyp]] += weight * weight;
    cands_count_[cms2.hyp_trilep_bucket()[i_hyp]]++;
    cands_passing_[TRILEPTON_ALL] += weight;
    cands_passing_w2_[TRILEPTON_ALL] += weight * weight;
    cands_count_[TRILEPTON_ALL]++;
  }

  // for the NMinus1Hist, the histogram checks the cuts for us
  htcmet_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), cms2.evt_tcmet(), weight);
  hptcmet_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), MetSpecialTrilep(cms2.evt_tcmet(), cms2.evt_tcmetPhi(), i_hyp), weight);

  std::multimap<float,float,std::greater<float> >::iterator entry = trileptonPt_.begin();  
  h_highest_lep_pt_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->first, weight);
  h_highest_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->second, weight);
  ++entry;
  h_second_highest_lep_pt_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->first, weight);
  h_second_highest_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->second, weight);
  ++entry;
  h_third_highest_lep_pt_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->first, weight);
  h_third_highest_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry->second, weight);

  std::multimap<float,float,std::greater<float> >::iterator entry_iso = trileptonIso_.begin();  
  h_highest_iso_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry_iso->first, weight);
  ++entry_iso;
  h_second_highest_iso_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry_iso->first, weight);
  ++entry_iso;
  h_third_highest_iso_lep_iso_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), entry_iso->first, weight);

  h_counter_electrons_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), addElectronsCounter_, weight);
//   if ( (cuts_passed & cuts_) == cuts_ && sample_.name == "wz" && electron_counter > 0 )
//     std::cout << "Additional electron: run: " << cms2.evt_run() << " event: " << cms2.evt_event() << " lumi: " << cms2.evt_lumiBlock() << std::endl;
  h_counter_muons_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), addMuonsCounter_, weight);
//   if ( (cuts_passed & cuts_) == cuts_ && sample_.name == "wz" && muon_counter > 0 )
//     std::cout << "Additional muon: run: " << cms2.evt_run() << " event: " << cms2.evt_event() << " lumi: " << cms2.evt_lumiBlock() << std::endl;

  h_njets_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), nJPTsTrilep(i_hyp, 20.), weight);
  h_DeltaPhiMETNearestLepton_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), nearestDeltaPhiTrilep(cms2.evt_tcmetPhi(), i_hyp), weight);

//   if ( !inZmassWindow(primZMass_) )
//     std::cout << "mass: " << primZMass_ << std::endl;
  h_primZMass_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), primZMass_, weight);
  h_addZMass_->Fill(cuts_passed, (TrileptonHypType)(cms2.hyp_trilep_bucket()[i_hyp]), addZMass_, weight);
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
		    "Sample %10s: Total candidate count (all): %8u."
		    " Total weight %10.1f +- %10.1f",sample_.name.c_str(),CandsCount(TRILEPTON_ALL),
		    CandsPassing(TRILEPTON_ALL) , RMS(TRILEPTON_ALL));
  // 		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
  // 		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
  // 		       sample_.name.c_str(),
  // 		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
  // 		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
  // 		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
  // 		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
  // 		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
    perror("writing to log file");
}
