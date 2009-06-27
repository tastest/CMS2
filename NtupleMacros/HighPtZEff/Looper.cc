#include <math.h>
#include "TVector3.h"
#include "TCanvas.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

Looper::Looper (Sample s, cuts_t c, const char *fname, bool usew) 
     : LooperBase(s, c, fname)
{
  useweight = usew;
  
  for( int i=0;i<ncounts;i++ ) {
    count[i].total = 0;
    count[i].denom = 0;
    count[i].numer = 0;
	count[0].failnumer = 0;
	count[i].onelep = 0;
    count[i].geneta = 0;
	count[i].gencuts = 0;
    //count[i].opp_sign = 0;
    count[i].same_flv = 0;
    count[i].pt20  = 0;
    count[i].el_good = 0;
    count[i].bad_mom = 0;
    count[i].no_match_z = 0;
    //count[i].multihyp = 0;
	count[i].dupematch = 0;
  }
  
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
  double drmatchmax = 0.1;
  int drmatchbins = 20;
  hels_dr_recostat1[0] = new TH1F( Form("%s_edr_recostat1_1", SampleName().c_str()), Form("%s_edr_recostat1_1", SampleName().c_str()), drmatchbins, 0, drmatchmax );
  hels_dr_recostat1[1] = new TH1F( Form("%s_edr_recostat1_2", SampleName().c_str()), Form("%s_edr_recostat1_2", SampleName().c_str()), drmatchbins, 0, drmatchmax );
  double drmax = 3.5;
  int drbins = int(drmax/0.05); // bin width = 0.05, was 0.1
  eff_edr_reco = new EffH1F( Form("%s_e_drstat1_reco", SampleName().c_str()), Form("%s_e_drstat1_reco", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_id = new EffH1F( Form("%s_e_drstat1_id", SampleName().c_str()), Form("%s_e_drstat1_id", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_iso = new EffH1F( Form("%s_e_drstat1_iso", SampleName().c_str()), Form("%s_e_drstat1_iso", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_id_iso = new EffH1F( Form("%s_e_drstat1_id_iso", SampleName().c_str()), Form("%s_e_drstat1_id_iso", SampleName().c_str()), drbins, 0, drmax );
  e_dr = new TH1F( Form("%s_e_drstat1", SampleName().c_str()), Form("%s_e_drstat1", SampleName().c_str()), drbins, 0, drmax );
  e_dr->Sumw2();
	
  double epmax = 600;
  int epbins = 60;
  double eptmax = 250;
  int eptbins = 50;
  for( int i=0;i<2;i++ ) {
	eff_ep_reco[i]  = new EffH1F( Form("%s_e%i_p_reco", SampleName().c_str(), i+1), Form("%s_e%i_p_reco", SampleName().c_str(), i+1), epbins, 0, epmax );
	eff_ept_reco[i] = new EffH1F( Form("%s_e%i_pt_reco", SampleName().c_str(), i+1), Form("%s_e%i_pt_reco", SampleName().c_str(), i+1), eptbins, 0, eptmax );
	eff_ep_id[i]  = new EffH1F( Form("%s_e%i_p_id", SampleName().c_str(), i+1), Form("%s_e%i_p_id", SampleName().c_str(), i+1), epbins, 0, epmax );
	eff_ept_id[i] = new EffH1F( Form("%s_e%i_pt_id", SampleName().c_str(), i+1), Form("%s_e%i_pt_id", SampleName().c_str(), i+1), eptbins, 0, eptmax );
	eff_ep_iso[i]  = new EffH1F( Form("%s_e%i_p_iso", SampleName().c_str(), i+1), Form("%s_e%i_p_iso", SampleName().c_str(), i+1), epbins, 0, epmax );
	eff_ept_iso[i] = new EffH1F( Form("%s_e%i_pt_iso", SampleName().c_str(), i+1), Form("%s_e%i_pt_iso", SampleName().c_str(), i+1), eptbins, 0, eptmax );
	eff_ep_id_iso[i]  = new EffH1F( Form("%s_e%i_p_id_iso", SampleName().c_str(), i+1), Form("%s_e%i_p_id_iso", SampleName().c_str(), i+1), epbins, 0, epmax );
	eff_ept_id_iso[i] = new EffH1F( Form("%s_e%i_pt_id_iso", SampleName().c_str(), i+1), Form("%s_e%i_pt_id_iso", SampleName().c_str(), i+1), eptbins, 0, eptmax );
	ineff_ep_reco[i]  = new EffH1F( Form("%s_e%i_p_reco_in", SampleName().c_str(), i+1), Form("%s_e%i_p_reco_in", SampleName().c_str(), i+1), epbins, 0, epmax );
	ineff_ept_reco[i] = new EffH1F( Form("%s_e%i_pt_reco_in", SampleName().c_str(), i+1), Form("%s_e%i_pt_reco_in", SampleName().c_str(), i+1), eptbins, 0, eptmax );
  }
  eff_ep_sum_reco  = new EffH1F( Form("%s_e_p_sum_reco", SampleName().c_str()), Form("%s_e_p_sum_reco", SampleName().c_str()), epbins, 0, epmax );
  eff_ept_sum_reco = new EffH1F( Form("%s_e_pt_sum_reco", SampleName().c_str()), Form("%s_e_pt_sum_reco", SampleName().c_str()), eptbins, 0, eptmax );

  double etamax = 2.6;
  int etabins = int((2*etamax)/0.1); // bin width = 0.1: hist width / n bins = width per bin
  for( int i=0;i<2;i++ ) {
	eff_eeta_reco[i]  = new EffH1F( Form("%s_e%i_eta_reco", SampleName().c_str(), i+1), Form("%s_e%i_eta_reco", SampleName().c_str(), i+1), etabins, -etamax, etamax );
	eff_eeta_id[i]  = new EffH1F( Form("%s_e%i_eta_id", SampleName().c_str(), i+1), Form("%s_e%i_eta_id", SampleName().c_str(), i+1), etabins, -etamax, etamax );
	eff_eeta_iso[i]  = new EffH1F( Form("%s_e%i_eta_iso", SampleName().c_str(), i+1), Form("%s_e%i_eta_iso", SampleName().c_str(), i+1), etabins, -etamax, etamax );
	eff_eeta_id_iso[i]  = new EffH1F( Form("%s_e%i_eta_id_iso", SampleName().c_str(), i+1), Form("%s_e%i_eta_id_iso", SampleName().c_str(), i+1), etabins, -etamax, etamax );
	e_eta[i] = new TH1F( Form("%s_e%i_eta", SampleName().c_str(), i+1), Form("%s_e%i_eta", SampleName().c_str(), i+1), etabins, -etamax, etamax );
	e_eta[i]->Sumw2();
  }

  double zpmax = 1000;
  int zpbins = 50;
  double zptmax = 600;
  int zptbins = 30;
  eff_zp_ind_reco = new EffH1F( Form("%s_Zp_ind_reco", SampleName().c_str()), Form("%s_Zp_ind_reco", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_ind_reco = new EffH1F( Form("%s_Zpt_ind_reco", SampleName().c_str()), Form("%s_Zpt_ind_reco", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_ind_id = new EffH1F( Form("%s_Zp_ind_id", SampleName().c_str()), Form("%s_Zp_ind_id", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_ind_id = new EffH1F( Form("%s_Zpt_ind_id", SampleName().c_str()), Form("%s_Zpt_ind_id", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_ind_iso = new EffH1F( Form("%s_Zp_ind_iso", SampleName().c_str()), Form("%s_Zp_ind_iso", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_ind_iso = new EffH1F( Form("%s_Zpt_ind_iso", SampleName().c_str()), Form("%s_Zpt_ind_iso", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_ind_id_iso = new EffH1F( Form("%s_Zp_ind_id_iso", SampleName().c_str()), Form("%s_Zp_ind_id_iso", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_ind_id_iso = new EffH1F( Form("%s_Zpt_ind_id_iso", SampleName().c_str()), Form("%s_Zpt_ind_id_iso", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_pr_reco = new EffH1F( Form("%s_Zp_pr_reco", SampleName().c_str()), Form("%s_Zp_pr_reco", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_pr_reco = new EffH1F( Form("%s_Zpt_pr_reco", SampleName().c_str()), Form("%s_Zpt_pr_reco", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_pr_id = new EffH1F( Form("%s_Zp_pr_id", SampleName().c_str()), Form("%s_Zp_pr_id", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_pr_id = new EffH1F( Form("%s_Zpt_pr_id", SampleName().c_str()), Form("%s_Zpt_pr_id", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_pr_iso = new EffH1F( Form("%s_Zp_pr_iso", SampleName().c_str()), Form("%s_Zp_pr_iso", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_pr_iso = new EffH1F( Form("%s_Zpt_pr_iso", SampleName().c_str()), Form("%s_Zpt_pr_iso", SampleName().c_str()), zptbins, 0, zptmax );
  eff_zp_pr_id_iso = new EffH1F( Form("%s_Zp_pr_id_iso", SampleName().c_str()), Form("%s_Zp_pr_id_iso", SampleName().c_str()), zpbins, 0, zpmax );
  eff_zpt_pr_id_iso = new EffH1F( Form("%s_Zpt_pr_id_iso", SampleName().c_str()), Form("%s_Zpt_pr_id_iso", SampleName().c_str()), zptbins, 0, zptmax );

  e_dr_zpt = new TH2F( Form("%s_e_drstat1_zpt", SampleName().c_str()), Form("%s_e_drstat1_zpt", SampleName().c_str()), drbins, 0, drmax, zptbins, 0, zptmax);
}


bool Looper::FilterEvent() { 
  return false; 
}


cuts_t Looper::Stat1Select(vector<int> idx) {
  cuts_t ret = 0;
  
  if( idx.size() != 2 ) {
	cout << "Improper idx size in Stat1Select: " << idx.size() << endl;
	return ret;
  }

  if( cms2.genps_id_mother()[cms2.genps_lepdaughter_idx()[idx[0]]] != 23 ||
	  cms2.genps_id_mother()[cms2.genps_lepdaughter_idx()[idx[1]]] != 23) {
	//int mother = cms2.genps_id_mother()[cms2.genps_lepdaughter_idx()[i]];
	//if( mother != 23 ) {
	//if( cms2.evt_event() < 5 )
	cout << "mother mismatch: " << "" << "\n";
	return ret;
  }
	 
  if( abs(cms2.genps_lepdaughter_p4()[idx[0]].eta()) < 2.4 &&
	  abs(cms2.genps_lepdaughter_p4()[idx[1]].eta()) < 2.4 )
	ret |= CUT_BIT(CUT_ETA24);

  if( cms2.genps_lepdaughter_p4()[idx[0]].pt() > 20. &&
	  cms2.genps_lepdaughter_p4()[idx[1]].pt() > 20. )
	ret |= CUT_BIT(CUT_PT20);

  double mass = ( cms2.genps_lepdaughter_p4()[idx[0]] + cms2.genps_lepdaughter_p4()[idx[1]] ).M();
  if( mass > 76. && mass < 106. )
	ret |= CUT_BIT(CUT_IN_Z_WINDOW);

  return ret;
}
//end Looper::Stat1Select

 
//2nd arg: 0=el, 1=mu
cuts_t Looper::LepSelect(int i, int flv) {
  cuts_t ret = 0;
  if( i == 999 ) //got no index--made default (initialized) to 999
	return 0;

  if( flv == 0 ) { //electron by reco
	//doing my own dr matching now--done in FillEventHistos
	
	//if( abs(cms2.els_p4()[i].eta()) < 2.4 ) //checked at gen
	//  ret |= CUT_BIT(CUT_ETA24);
  
	//if( cms2.els_p4()[i].pt() > 10 ) //MAKE SURE
	//  ret |= (CUT_BIT(CUT_PT20));

	if( goodElectronWithoutIsolation(i) )
	  //if( cms2.els_tightId22XMinMatteo().at(i) == 1 )
	  //cms2.els_closestMuon().at(i) == -1 &&
	  //TMath::Abs(cms2.els_d0corr().at(i)) < 0.025 )
	  ret |= CUT_BIT(CUT_EL_GOOD);

	if( passElectronIsolation(i, true) )
	  //el_rel_iso > 0.92 -- trk, calo, hcal
	  ret |= CUT_BIT(CUT_EL_ISO);
  }
  else if( flv == 1 ) { //muon
	//if( abs(cms2.mus_p4()[i].eta()) < 2.4 ) //checked at gen
	//  ret |= CUT_BIT(CUT_ETA24);

	//if( cms2.mus_p4()[i].pt() > 10 ) //MAKE SURE
	//  ret |= CUT_BIT(CUT_PT20);

	if( cms2.mus_type()[i] & 0x2 )
	  ret |= CUT_BIT(CUT_MU_GLOBAL);
	
	if( goodMuonWithoutIsolation( i ) )
	  ret |= CUT_BIT(CUT_MU_GOOD);

	if( passMuonIsolation( i ) )
	  ret |= CUT_BIT(CUT_MU_ISO);
  }
  else
	cout << "Bad flv value in LepSelect\n";

  return ret;
}
//end Looper::LepSelect


cuts_t Looper::EventSelect () {
  cuts_t ret = 0;
  return ret;
}


void Looper::FillEventHistos ()
{
  //cout << "start FillEventHistos\n";
  if( useweight )
	weight = Weight(0);
  else
	weight = 1;

  denomitr = 0;

  for( int i=0;i<ncounts;i++ )
    count[i].total += weight;

  //get status 3 Z by looping on genps status 3 block
  int zidx = 0;
  for(unsigned int i=0; i<cms2.genps_id().size(); i++) {
	if( cms2.genps_id()[i] == 23 ) {
	  zidx = i;
	}
  }

  vector<int> idxlep1;
  //get status 1 leptons
  for( unsigned int i=0;i<cms2.genps_lepdaughter_id().size();i++ ) {
	if( abs(cms2.genps_lepdaughter_id()[i]) == 11 ||
		abs(cms2.genps_lepdaughter_id()[i]) == 13  ) 
	  idxlep1.push_back(i);	  
  }
  if( idxlep1.size() != 2 ) {
	cout << "\n\nIMPROPER idxlep1 SIZE";
	return;
  }
  //order by pt
  if( cms2.genps_lepdaughter_p4()[idxlep1[1]].pt() > cms2.genps_lepdaughter_p4()[idxlep1[0]].pt() ) {
	int tmp = idxlep1[1];
	idxlep1[1] = idxlep1[0];
	idxlep1[0] = tmp;
  }

  cuts_t stat1pass = Stat1Select(idxlep1);
  if( (stat1pass & stat1cuts) != stat1cuts ) {
	for( int i=0;i<ncounts;i++ ) count[i].gencuts += weight;
	return; //don't do anything if don't pass stat1cuts, just because i defined the denominator this way--should not bias results
  }
  count[0].denom += weight;
  double drstat1 = abs( ROOT::Math::VectorUtil::DeltaR(cms2.genps_lepdaughter_p4()[idxlep1[0]], cms2.genps_lepdaughter_p4()[idxlep1[1]]) );

  for( int i=0;i<2;i++ ) { //all denoms same--see comment below
	eff_ep_reco[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	eff_ept_reco[i]->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	eff_eeta_reco[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	eff_zp_ind_reco->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_ind_reco->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
	//id
	eff_ep_id[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	eff_ept_id[i]->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	eff_eeta_id[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	eff_zp_ind_id->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_ind_id->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
	//iso
	eff_ep_iso[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	eff_ept_iso[i]->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	eff_eeta_iso[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	eff_zp_ind_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_ind_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
	//id_iso
	eff_ep_id_iso[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	eff_ept_id_iso[i]->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	eff_eeta_id_iso[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	eff_zp_ind_id_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_ind_id_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
	//no repeat
	ineff_ep_reco[i] ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	ineff_ept_reco[i]->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
  }
  //for this one, need both to match to pass because if you only reconstruct one of two leptons, you have no dr.
  eff_edr_reco->denom->Fill( drstat1, weight );
  eff_edr_id->denom->Fill( drstat1, weight );
  eff_edr_iso->denom->Fill( drstat1, weight );
  eff_edr_id_iso->denom->Fill( drstat1, weight );
  e_dr_zpt->Fill( drstat1, cms2.genps_p4()[zidx].pt(), weight );
  e_dr->Fill( drstat1, weight );
  //Z effs for both els
  eff_zp_pr_reco->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
  eff_zpt_pr_reco->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
  eff_zp_pr_id->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
  eff_zpt_pr_id->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
  eff_zp_pr_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
  eff_zpt_pr_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
  eff_zp_pr_id_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
  eff_zpt_pr_id_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );

  
  // dr matching
  //cout << "dr matching\n";
  const double maxcone = 0.1;
  double min_dr[2] = {999,999};
  unsigned int idxlepreco[2] = {999,999};
  for( unsigned int i=0;i<idxlep1.size();i++ ) {
	for( unsigned int j=0;j<cms2.els_p4().size();j++ ) {
	  double dr = ROOT::Math::VectorUtil::DeltaR(cms2.genps_lepdaughter_p4()[idxlep1[i]], cms2.els_p4()[j]);
	  if( dr < maxcone && dr < min_dr[i] ) {
		if( i == 0 || (i>0 && idxlepreco[0] != j) ) { //make sure both stat1 match the same reco
		  min_dr[i] = dr;
		  idxlepreco[i] = j;
		}
		else
		  count[0].dupematch += weight;
	  }
	}
  }

  //cout << "check denom \n";
  // reco efficiency passes if dr < maxcone
  // change for new (6/17/09): do individual electrons, but with same denom.
  //  this means that the numerator histograms versus zp/pt can be filled 0, 1, or 2 times for the number of els which pass matching
  //  the denom is always filled twice for Z p/pt, so long as it passes
  //  for the eff vs ep/pt, order by p/pt for all cases, but only require that one per event pass--so fill for the one that passes regardless of other
  // notes (6/25/09): we want four plots per x-axis: dr, dr+id, dr+iso, dr+id+iso
  //  therefore, do in three stages (even though i will have 4 for loops): first check dr and put in cut bit
  //  then fill other two vars of cut bit using LepSelect, then one loop for each of the other two cuts
  //  and after each loop check pair
  cuts_t e_reco_cuts[2] = {0,0};
  //fill reco
  for( int i=0;i<2;i++ ) {
	if( min_dr[i] < maxcone ) {
	  eff_ep_reco[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	  eff_ept_reco[i]->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	  eff_eeta_reco[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	  eff_zp_ind_reco->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	  eff_zpt_ind_reco->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	  
	  e_eta[i]->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight ); //just here
	  hels_dr_recostat1[i]->Fill( min_dr[i], weight ); //just here
	  e_reco_cuts[i] |= CUT_BIT(CUT_EL_DR); //just here
	}
	else {
	  ineff_ep_reco[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	  ineff_ept_reco[i]->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	}
  }

  //counting + reco pair plot
  if( min_dr[0] < maxcone && min_dr[1] < maxcone ) {
	eff_edr_reco->numer->Fill( drstat1, weight );
	eff_zp_pr_reco->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_pr_reco->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	count[0].numer += weight;
  }
  else if( min_dr[0] < maxcone || min_dr[1] < maxcone ) 
	count[0].onelep += weight;
  else
	count[0].failnumer += weight;

  //cout << "check numer \n";
  //***don't do this here, because idxlep is 999 if they don't match...but in that case, LepSelect returns zero, so ok...
  e_reco_cuts[0] |= LepSelect(idxlepreco[0],0); //second arg is for flavor, el = 0
  e_reco_cuts[1] |= LepSelect(idxlepreco[1],0);
  //cout << "fill numer \n";

  //fill id
  cuts_t id_cut = ( CUT_BIT(CUT_EL_DR) | CUT_BIT(CUT_EL_GOOD) );
  for( int i=0;i<2;i++ ) {
	//if( min_dr[i] < maxcone && ((e_reco_cuts[i] & id_cut) == id_cut) ) {
	if( (e_reco_cuts[i] & id_cut) == id_cut ) {
	  eff_ep_id[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	  eff_ept_id[i]->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	  eff_eeta_id[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	  eff_zp_ind_id->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	  eff_zpt_ind_id->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	}
  }
  
  //if( min_dr[0] < maxcone && min_dr[1] < maxcone && ((e_reco_cuts[0] & e_reco_cuts[1] & id_cut) == id_cut) ) {
  if( (e_reco_cuts[0] & e_reco_cuts[1] & id_cut) == id_cut ) {
	eff_edr_id->numer->Fill( drstat1, weight );
	eff_zp_pr_id->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_pr_id->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }

  //fill iso
  cuts_t iso_cut = ( CUT_BIT(CUT_EL_DR) | CUT_BIT(CUT_EL_ISO) );
  for( int i=0;i<2;i++ ) {
	//if( min_dr[i] < maxcone && ((e_reco_cuts[i] & iso_cut) == iso_cut) ) {
	if( (e_reco_cuts[i] & iso_cut) == iso_cut ) {
	  eff_ep_iso[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	  eff_ept_iso[i]->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	  eff_eeta_iso[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	  eff_zp_ind_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	  eff_zpt_ind_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	}
  }
  
  //if( min_dr[0] < maxcone && min_dr[1] < maxcone && ((e_reco_cuts[0] & e_reco_cuts[1] & iso_cut) == iso_cut) ) {
  if( (e_reco_cuts[0] & e_reco_cuts[1] & iso_cut) == iso_cut ) {
	eff_edr_iso->numer->Fill( drstat1, weight );
	eff_zp_pr_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_pr_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }

  //fill id_iso
  cuts_t id_iso_cut = ( CUT_BIT(CUT_EL_DR) | CUT_BIT(CUT_EL_GOOD) | CUT_BIT(CUT_EL_ISO) );
  for( int i=0;i<2;i++ ) {
	//if( min_dr[i] < maxcone && ((e_reco_cuts[i] & id_iso_cut) == id_iso_cut) ) {
	if( (e_reco_cuts[i] & id_iso_cut) == id_iso_cut ) {
	  eff_ep_id_iso[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].P(), weight );
	  eff_ept_id_iso[i]->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].pt(), weight );
	  eff_eeta_id_iso[i] ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[i]].eta(), weight );
	  eff_zp_ind_id_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	  eff_zpt_ind_id_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	}
  }
  
  //if( min_dr[0] < maxcone && min_dr[1] < maxcone && ((e_reco_cuts[0] & e_reco_cuts[1] & id_iso_cut) == id_iso_cut) ) {
  if( (e_reco_cuts[0] & e_reco_cuts[1] & id_iso_cut) == id_iso_cut ) {
	eff_edr_id_iso->numer->Fill( drstat1, weight );
	eff_zp_pr_id_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_zpt_pr_id_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }


  /*if( (e_reco_cuts & els_iso_denom) == els_iso_denom ) {
	eff_p_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_pt_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }
  if( (e_reco_cuts & els_iso_numer) == els_iso_numer ) {
	eff_p_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_pt_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
	}*/
  
  //cout << "end FillEventHistos\n";
}
//end Looper::FillEventHistos()

void Looper::End ()
{

  double etamin = 0.3;
  double epmin = 0.5;
  double zpmin = 0.5;
  double zptmin = 0.5;

  for( int i=0;i<2;i++ ) {
	eff_ep_reco[i] ->MakeEff(epmin);
	eff_ept_reco[i]->MakeEff( );
	eff_eeta_reco[i]->MakeEff(etamin);
	eff_ep_id[i] ->MakeEff(epmin);
	eff_ept_id[i]->MakeEff( );
	eff_eeta_id[i]->MakeEff(etamin);
	eff_ep_iso[i] ->MakeEff(epmin);
	eff_ept_iso[i]->MakeEff( );
	eff_eeta_iso[i]->MakeEff(etamin);
	eff_ep_id_iso[i] ->MakeEff(epmin);
	eff_ept_id_iso[i]->MakeEff( );
	eff_eeta_id_iso[i]->MakeEff(etamin);
	//no repeat
	ineff_ep_reco[i] ->MakeEff( 0, 0.3 );
	ineff_ept_reco[i]->MakeEff( 0, 0.3 );
  }
  //eff_edr_reco->xmin = 0.0;
  eff_edr_reco->xmax = 1.0;
  eff_edr_id->xmax = 1.0;
  eff_edr_iso->xmax = 1.0;
  eff_edr_id_iso->xmax = 1.0;
  eff_edr_reco->MakeEff( );// 0.0, 1.0 );
  eff_edr_id->MakeEff( );// 0.0, 1.0 );
  eff_edr_iso->MakeEff( );// 0.0, 1.0 );
  eff_edr_id_iso->MakeEff( );// 0.0, 1.0 );
  //note: Rebin should only be called after MakeEff
  eff_edr_reco->Rebin( );
  eff_edr_id->Rebin( );
  eff_edr_iso->Rebin( );
  eff_edr_id_iso->Rebin( );
  
  //z indiv
  eff_zp_ind_reco ->MakeEff(zpmin);
  eff_zpt_ind_reco->MakeEff(zptmin);
  eff_zp_ind_id ->MakeEff(zpmin);
  eff_zpt_ind_id->MakeEff(zptmin);
  eff_zp_ind_iso ->MakeEff(zpmin);
  eff_zpt_ind_iso->MakeEff(zptmin);
  eff_zp_ind_id_iso ->MakeEff(zpmin);
  eff_zpt_ind_id_iso->MakeEff(zptmin);
  //z pair
  eff_zp_pr_reco ->MakeEff(zpmin);
  eff_zp_pr_id ->MakeEff(zpmin);
  eff_zp_pr_iso ->MakeEff(zpmin);
  eff_zp_pr_id_iso ->MakeEff(zpmin);
  //eff_zpt_pr_reco->MakeEff(zptmin);
  //eff_zpt_pr_id->MakeEff(zptmin);
  //eff_zpt_pr_iso->MakeEff(zptmin);
  //eff_zpt_pr_id_iso->MakeEff(zptmin);
  eff_zpt_pr_reco->MakeEff();
  eff_zpt_pr_id->MakeEff();
  eff_zpt_pr_iso->MakeEff();
  eff_zpt_pr_id_iso->MakeEff();

  //sum
  eff_ep_sum_reco ->numer->Add( eff_ep_reco[0]->numer, eff_ep_reco[1]->numer );
  eff_ep_sum_reco ->denom->Add( eff_ep_reco[0]->denom, eff_ep_reco[1]->denom );
  eff_ep_sum_reco->MakeEff();
  eff_ept_sum_reco->numer->Add( eff_ept_reco[0]->numer, eff_ept_reco[1]->numer );
  eff_ept_sum_reco->denom->Add( eff_ept_reco[0]->denom, eff_ept_reco[1]->denom );
  eff_ept_sum_reco->MakeEff();

  //eff_e1p_reco->MakeEff( binentries ); //was for rebinning, now off
  //eff_e2p_reco->MakeEff( binentries );
  //eff_e1pt_reco->MakeEff( binentries );
  //eff_e2pt_reco->MakeEff( binentries );


  //double ymin = 0.7;
  //double ymax = 1.;
  TCanvas *c1 = new TCanvas();// eff->GetName(), eff->GetName());
  //e_dr_zpt->GetYaxis()->SetRangeUser(ymin, ymax);
  e_dr_zpt->Draw("box");
  c1->SaveAs((TString)e_dr_zpt->GetName()+".png");
  //e_dr->GetYaxis()->SetRangeUser(ymin, ymax);
  e_dr->Draw();
  c1->SaveAs((TString)e_dr->GetName()+".png");

  //cout << "\n\nweight = " << weight << endl;

  //print struct count content 
  for( int i=0;i<ncounts;i++ ) {
    if( i == 0 )
	  cout << "\n\nReco Efficiencies:";
    else if ( i == 1 )
	  cout << "\n\nIso Efficiencies:";
    //else if ( i == 2 )
	//cout << "\n\nReco(3) Efficiencies:";
    cout << "\nTotal events run on: " << count[i].total
         << "\ntotal denominator events: " << count[i].denom
         << "\ntotal numerator events (both pass): " << count[i].numer
		 << "\nnum events with one lepton fail numer: " << count[i].onelep 
		 << "\nnum events with both lepton fail numer: " << count[i].failnumer
		 << "\nsum of above 3: " << count[i].numer + count[i].onelep + count[i].failnumer
	  //<< "\nnum events without two gen leptons in eta 2.4: " << count[i].geneta
		 << "\nnum events failing gen cuts: " << count[i].gencuts
	  //<< "\nnum events which fail opp sign: " << count[i].opp_sign
		 << "\nnum events with two stat1 match to same reco: " << count[i].dupematch
      //<< "\nnum events which fail same flv: " << count[i].same_flv 
	  //<< "\nnum events which fail pt > 20: " << count[i].pt20
         << "\nnum events which fail good lepton: " << count[i].el_good
      //<< "\nnum leptons with mc_motherid != 23: " << count[i].bad_mom 
	  //<< "\nnum events with (num leptons with mc_motherid == 23) != 2 : " << count[i].no_match_z
	  //<< "\nnum events with > 1 passing hypothesis: " << count[i].multihyp
         << endl;
  }
  cout << endl << endl;

  int ret = fprintf(logfile_, 
		       "Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
		       " Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
		       sample_.name.c_str(),
		       CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU),
					   CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
		       CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
		       CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
		       CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
		       CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
	perror("writing to log file");
}
