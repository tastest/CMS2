#include <math.h>
#include "TVector3.h"
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
    count[i].geneta = 0;
	count[i].gencuts = 0;
    count[i].opp_sign = 0;
    count[i].same_flv = 0;
    count[i].pt20  = 0;
    count[i].el_good = 0;
    count[i].bad_mom = 0;
    count[i].no_match_z = 0;
    count[i].multihyp = 0;
  }
  
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
  double drmax = 0.1;
  int drbins = 20;
  hels_dr_recostat1[0] = new TH1F( Form("%s_dr_recostat1_1", SampleName().c_str()), Form("%s_dr_recostat1_1", SampleName().c_str()), drbins, 0, drmax );
  hels_dr_recostat1[1] = new TH1F( Form("%s_dr_recostat1_2", SampleName().c_str()), Form("%s_dr_recostat1_2", SampleName().c_str()), drbins, 0, drmax );

  double epmax = 600;
  int epbins = 60;
  eff_e1p_reco = new EffH1F( Form("%s_e1_p_reco", SampleName().c_str()), Form("%s_e1_p_reco", SampleName().c_str()), epbins, 0, epmax );
  eff_e2p_reco = new EffH1F( Form("%s_e2_p_reco", SampleName().c_str()), Form("%s_e2_p_reco", SampleName().c_str()), epbins, 0, epmax );
  double eptmax = 250;
  int eptbins = 50;
  eff_e1pt_reco = new EffH1F( Form("%s_e1_pt_reco", SampleName().c_str()), Form("%s_e1_pt_reco", SampleName().c_str()), eptbins, 0, eptmax );
  eff_e2pt_reco = new EffH1F( Form("%s_e2_pt_reco", SampleName().c_str()), Form("%s_e2_pt_reco", SampleName().c_str()), eptbins, 0, eptmax );

  double zpmax = 800;
  int zpbins = 80;
  eff_p_reco = new EffH1F( Form("%s_p_reco", SampleName().c_str()), Form("%s_p_reco", SampleName().c_str()), zpbins, 0, zpmax );
  eff_p_iso = new EffH1F( Form("%s_p_iso", SampleName().c_str()), Form("%s_p_iso", SampleName().c_str()), zpbins, 0, zpmax );
  double zptmax = 500;
  int zptbins = 50;
  eff_pt_reco = new EffH1F( Form("%s_pt_reco", SampleName().c_str()), Form("%s_pt_reco", SampleName().c_str()), zptbins, 0, zptmax );
  eff_pt_iso = new EffH1F( Form("%s_pt_iso", SampleName().c_str()), Form("%s_pt_iso", SampleName().c_str()), zptbins, 0, zptmax );
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

  if( flv == 0 ) { //electron by reco
	//doing my own dr matching now
	//if( abs(cms2.els_mc_id()[i]) == 11 ) { //mc truth
	//  ret |= (CUT_BIT(CUT_MC_EL));
	//  if( cms2.els_mc_motherid()[i] != 23 ) 
	//	count[denomitr].bad_mom += weight;
	//}
	//
	//if( abs(cms2.els_mc3_id()[i]) == 11 ) { //mc truth
	//  ret |= (CUT_BIT(CUT_MC3_EL));
	//}
	
	//if( abs(cms2.els_p4()[i].eta()) < 2.4 ) //checked at gen
	//  ret |= CUT_BIT(CUT_ETA24);
  
	//if( cms2.els_p4()[i].pt() > 10 ) //MAKE SURE
	//  ret |= (CUT_BIT(CUT_PT20));

	//if( goodElectronWithoutIsolation(cms2.els_mcidx()[i]) ) //replaced below
	if( cms2.els_tightId22XMinMatteo().at(i) == 1 )
	  //cms2.els_closestMuon().at(i) == -1 &&
	  //TMath::Abs(cms2.els_d0corr().at(i)) < 0.025 )
	  ret |= CUT_BIT(CUT_EL_GOOD);

	if( passElectronIsolation(i, true) )
	  //el_rel_iso > 0.92 -- trk, calo, hcal
	  ret |= CUT_BIT(CUT_EL_ISO);
  }
  else if( flv == 1 ) { //muon
	if( abs(cms2.mus_mc_id()[i]) == 13 ) {
	  ret |= (CUT_BIT(CUT_MC_MU));
	  if( cms2.mus_mc_motherid()[i] != 23 ) 
		count[denomitr].bad_mom += weight;
	}

	if( abs(cms2.mus_mc3_id()[i]) == 13 ) {
	  ret |= (CUT_BIT(CUT_MC3_MU));
	}

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
  //order by pt
  if( cms2.genps_lepdaughter_p4()[idxlep1[1]].pt() > cms2.genps_lepdaughter_p4()[idxlep1[0]].pt() ) {
	int tmp = idxlep1[1];
	idxlep1[1] = idxlep1[0];
	idxlep1[0] = tmp;
  }

  cuts_t stat1pass = Stat1Select(idxlep1);
  if( (stat1pass & stat1cuts) != stat1cuts ) {
	for( int i=0;i<ncounts;i++ ) count[i].gencuts += weight;
	return; //don't do anything if don't pass stat1cuts
  }

  count[0].denom += weight;
  eff_e1p_reco ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[0]].P(), weight );
  eff_e2p_reco ->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[1]].P(), weight );
  eff_e1pt_reco->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[0]].pt(), weight );
  eff_e2pt_reco->denom->Fill( cms2.genps_lepdaughter_p4()[idxlep1[1]].pt(), weight );
  eff_p_reco->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
  eff_pt_reco->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );

  // dr matching
  double maxcone = 0.1;
  double min_dr[2] = {999,999};
  int idxlepreco[2] = {0,0};
  for( unsigned int i=0;i<idxlep1.size();i++ ) {
	for( unsigned int j=0;j<cms2.els_p4().size();j++ ) {
	  double dr = ROOT::Math::VectorUtil::DeltaR(cms2.genps_lepdaughter_p4()[idxlep1[i]], cms2.els_p4()[j]);
	  if( dr < maxcone && dr < min_dr[i] ) {
		min_dr[i] = dr;
		idxlepreco[i] = j;
	  }
	}
  }

  // reco efficiency passes if dr < maxcone
  if( min_dr[0] < maxcone && min_dr[1] < maxcone ) {
	count[0].numer += weight;
	hels_dr_recostat1[0]->Fill( min_dr[0], weight );
	hels_dr_recostat1[1]->Fill( min_dr[1], weight );
	eff_e1p_reco ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[0]].P(), weight );
	eff_e2p_reco ->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[1]].P(), weight );
	eff_e1pt_reco->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[0]].pt(), weight );
	eff_e2pt_reco->numer->Fill( cms2.genps_lepdaughter_p4()[idxlep1[1]].pt(), weight );
	eff_p_reco->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_pt_reco->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }
  else
	return;

  cuts_t reco_e = ( LepSelect(idxlepreco[0],0) & LepSelect(idxlepreco[1],0) );

  if( (reco_e & els_iso_denom) == els_iso_denom ) {
	eff_p_iso->denom->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_pt_iso->denom->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }
  if( (reco_e & els_iso_numer) == els_iso_numer ) {
	eff_p_iso->numer->Fill( cms2.genps_p4()[zidx].P(), weight );
	eff_pt_iso->numer->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }
  

}
//end Looper::FillEventHistos()

void Looper::End ()
{

  eff_e1p_reco->MakeEff( );
  eff_e2p_reco->MakeEff( );
  eff_e1pt_reco->MakeEff( );
  eff_e2pt_reco->MakeEff( );
  eff_p_reco->MakeEff( );
  eff_pt_reco->MakeEff( );
  eff_p_iso->MakeEff( 0.6 );
  eff_pt_iso->MakeEff( 0.6 );

  //eff_e1p_reco->MakeEff( binentries );
  //eff_e2p_reco->MakeEff( binentries );
  //eff_e1pt_reco->MakeEff( binentries );
  //eff_e2pt_reco->MakeEff( binentries );


  //print struct count content 
  for( int i=0;i<ncounts;i++ ) {
    if( i == 0 )
      cout << "\n\nIso Efficiencies:";
    else if ( i == 1 )
      cout << "\n\nReco(1) Efficiencies:";
    else if ( i == 2 )
      cout << "\n\nReco(3) Efficiencies:";
    cout << "\nTotal events run on: " << count[i].total
         << "\ntotal denominator events: " << count[i].denom
         << "\ntotal numerator events: " << count[i].numer
	  //<< "\nnum events without two gen leptons in eta 2.4: " << count[i].geneta
		 << "\nnum events failing gen cuts: " << count[i].gencuts
         << "\nnum events which fail opp sign: " << count[i].opp_sign
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
