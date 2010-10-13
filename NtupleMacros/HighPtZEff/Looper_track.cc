#include <math.h>
#include <algorithm>
#include "Math/LorentzVector.h"
//#include "DataFormats/Math/interface/LorentzVector.h" //for cmssw....
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "TCanvas.h"
#include "Tools/tools.h"
#include "Looper.h"
#include "/home/users/wandrews/macros/comparison.C"

//typedef math::XYZTLorentzVector LorentzVector; //for cmssw...
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

//Looper::Looper (Sample s, cuts_t c, const char *fname) 
//     : LooperBase(s, c, fname)
//{
Looper::Looper (Sample s, cuts_t c, const char *fname, bool usew) 
  : LooperBase(s, c, fname)
{
  useweight = usew;
  numprint = 0;

  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
  int isobins = 201; //int((1.0+0.005)/0.005); //201 now : smaller bins--.005 which is 20 per 0.1
  double isomax = 1.0 + 0.005; // width = nbins * bin_width
  e_hcal_iso = new TH1F( Form("%s_e_hcal_iso", SampleName().c_str()), Form("%s_e_hcal_iso", SampleName().c_str()), isobins, 0, isomax );
  e_hcal_iso_dr05_1 = new TH1F( Form("%s_e_hcal_iso_dr05_1", SampleName().c_str()), Form("%s_e_hcal_iso_dr05_1", SampleName().c_str()), isobins, 0, isomax );
  e_ecal_iso = new TH1F( Form("%s_e_ecal_iso", SampleName().c_str()), Form("%s_e_ecal_iso", SampleName().c_str()), isobins, 0, isomax );
  for(int i=0;i<2;i++){
	e_ecal_iso_bare[i] = new TH1F( Form("%s_e%i_ecal_iso_bare", SampleName().c_str(),i+1), Form("%s_e%i_ecal_iso_bare", SampleName().c_str(),i+1), 100, 0, 20 );
  }
  e_ecal_iso_dr05_1 = new TH1F( Form("%s_e_ecal_iso_dr05_1", SampleName().c_str()), Form("%s_e_ecal_iso_dr05_1", SampleName().c_str()), isobins, 0, isomax );
  e_trck_iso = new TH1F( Form("%s_e_trck_iso", SampleName().c_str()), Form("%s_e_trck_iso", SampleName().c_str()), isobins, 0, isomax );
  e_trck_iso_dr05_1 = new TH1F( Form("%s_e_trck_iso_dr05_1", SampleName().c_str()), Form("%s_e_trck_iso_dr05_1", SampleName().c_str()), isobins, 0, isomax );
  e_trck_iso_recalc = new TH1F( Form("%s_e_trck_iso_recalc", SampleName().c_str()), Form("%s_e_trck_iso_recalc", SampleName().c_str()), isobins, 0, isomax );
  e_trck_iso_recalc_dr05_1 = new TH1F( Form("%s_e_trck_iso_recalc_dr05_1", SampleName().c_str()), Form("%s_e_trck_iso_recalc_dr05_1", SampleName().c_str()), isobins, 0, isomax );
  e_trck_iso_affble = new TH1F( Form("%s_e_trck_iso_affable", SampleName().c_str()), Form("%s_e_trck_iso_affable", SampleName().c_str()), isobins, 0, isomax ); 
  e_trck_iso_affble_dr05_1 = new TH1F( Form("%s_e_trck_iso_affable_dr05_1", SampleName().c_str()), Form("%s_e_trck_iso_affable_dr05_1", SampleName().c_str()), isobins, 0, isomax );
  //e_trck_iso_match015 = new TH1F( Form("%s_e_trck_iso_match015", SampleName().c_str()), Form("%s_e_trck_iso_match015", SampleName().c_str()), isobins, 0, isomax );

  string etaregname[etaregions];
  etaregname[0] = "bb";
  etaregname[1] = "be";
  etaregname[2] = "ee";

  double drmax = 3.5;
  int drbins = int(drmax/0.05); // bin width = 0.05, was 0.1
  //INDIVIDUAL--SEPARATELY FOR LEADING AND SUBLEADING
  for( int i=0; i<2; i++ ) {
	eff_edr_hcal_iso[i] = new EffH1F( Form("%s_e%i_drstat1_hcal_iso", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_hcal_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_ecal_iso[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_ecal_iso_soft[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_soft", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso_soft", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_ecal_iso_subt[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_subt", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso_subt", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_trck_iso[i] = new EffH1F( Form("%s_e%i_drstat1_trck_iso", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_trck_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_trck_iso_recalc[i] = new EffH1F( Form("%s_e%i_drstat1_trck_iso_recalc", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_trck_iso_recalc", SampleName().c_str(), i+1), drbins, 0, drmax );
	eff_edr_trck_iso_affble[i] = new EffH1F( Form("%s_e%i_drstat1_trck_iso_affable", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_trck_iso_affable", SampleName().c_str(), i+1), drbins, 0, drmax );
	//eff_edr_ecal_iso_bb[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_bb", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	//eff_edr_ecal_iso_be[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_be", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	//eff_edr_ecal_iso_ee[i] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_ee", SampleName().c_str(), i+1), Form("%s_e%i_drstat1_ecal_iso", SampleName().c_str(), i+1), drbins, 0, drmax );
	for( int j=0; j<etaregions; j++ ) {
	  eff_edr_ecal_iso_reg[i][j] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_%s", SampleName().c_str(), i+1, etaregname[j].c_str()), Form("%s_e%i_drstat1_ecal_iso_%s", SampleName().c_str(), i+1, etaregname[j].c_str()), drbins, 0, drmax );
	}
  }
  //BOTH ELS HAVE TO BE ISOLATED FOR NUMERATOR (PAIR)
  eff_edr_hcal_iso_pair = new EffH1F( Form("%s_e_drstat1_hcal_iso_pair", SampleName().c_str()), Form("%s_e_drstat1_hcal_iso_pair", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_ecal_iso_pair = new EffH1F( Form("%s_e_drstat1_ecal_iso_pair", SampleName().c_str()), Form("%s_e_drstat1_ecal_iso_pair", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_trck_iso_pair = new EffH1F( Form("%s_e_drstat1_trck_iso_pair", SampleName().c_str()), Form("%s_e_drstat1_trck_iso_pair", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_trck_iso_recalc_pair = new EffH1F( Form("%s_e_drstat1_trck_iso_recalc_pair", SampleName().c_str()), Form("%s_e_drstat1_trck_iso_recalc_pair", SampleName().c_str()), drbins, 0, drmax );
  eff_edr_trck_iso_affble_pair = new EffH1F( Form("%s_e_drstat1_trck_iso_affable_pair", SampleName().c_str()), Form("%s_e_drstat1_trck_iso_affable_pair", SampleName().c_str()), drbins, 0, drmax );

  int njetdistbins = 9; //for plotting njets, not iso
  e_njets = new TH1F( Form("%s_njets", SampleName().c_str()), Form("%s_njets", SampleName().c_str()), njetdistbins, 0, njetdistbins );
  e_njets->Sumw2();
  e_njets_clean = new TH1F( Form("%s_njets_clean", SampleName().c_str()), Form("%s_njets_clean", SampleName().c_str()), njetdistbins, 0, njetdistbins );
  e_njets_clean->Sumw2();
  for( int i=0; i<2; i++ ) {
	for( int j=0; j<njetbins; j++ ) {
	  eff_edr_ecal_iso_njet[i][j] = new EffH1F( Form("%s_e%i_drstat1_ecal_iso_njet_%i", SampleName().c_str(), i+1, j), Form("%s_e%i_drstat1_ecal_iso_njet_%i", SampleName().c_str(), i+1, j), drbins, 0, drmax );
	}
  }


  //print header for printing of evt and iso info
  //else if( i == 1 && ftrcknumer && strcknumer && fecalnumer ) //first passed both ecal+trck, second passed trck, failed ecal
  //cout << SampleName().c_str() << "   " << cms2.evt_run() << "   "
  //<< cms2.evt_event() << "   " << cms2.evt_lumiBlock() << "   "
  //<< ecaliso[i] << endl;
  cout << "sample\t\t\tRun\tevent\tlumi\tecaliso2\n\n";
}


bool Looper::FilterEvent() { return false; }

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
  //if( cms2.genps_lepdaughter_p4()[idx[0]].pt() > 60. &&
  //cms2.genps_lepdaughter_p4()[idx[1]].pt() > 60. )
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

//reproduce pat track isolation
double track_iso(int els_idx) {
  //input is index of electron in els block
  //return is the sum of the pt of all tracks in the range 0.015 < dR < 0.3 around the els

  double isolation = 0;
  //for( unsigned int i=0; i<cms2.trks_trk_p4().size(); i++ ) {
  for( unsigned int i=0; i<cms2.trks_trk_p4().size(); i++ ) {
	//cuts on track quality
	if( cms2.trks_trk_p4()[i].pt() <= 1.0 )
	//if( cms2.trks_trk_p4()[i].pt() <= 1.0 || abs( cms2.trks_z0corr()[i] ) >= 0.2 )
	//if( cms2.trks_trk_p4()[i].pt() <= 1.0 || abs( cms2.els_vertex_p4()[i].z() - cms2.trks_vertex_p4()[i].z() ) >= 0.2 )
	//if( cms2.trks_trk_p4()[i].pt() <= 1.0 || abs( cms2.els_z0corr()[i] - cms2.trks_z0corr()[i] ) > 0.2 )
	//if( cms2.trks_trk_p4()[i].pt() <= 1.0 || abs( cms2.els_z0()[i] - cms2.trks_z0()[i] ) > 0.2 )
	  continue;
	//double dR = ROOT::Math::VectorUtil::DeltaR( cms2.els_trk_p4()[els_idx], cms2.trks_trk_p4()[i] );
	double dR = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4In()[els_idx], cms2.trks_trk_p4()[i] );
	if( dR > 0.015 && dR < 0.3 )
	//if( dR < 0.3 )
	  isolation += cms2.trks_trk_p4()[i].pt();
  }

  return isolation;
}

//improve upon pat track isolation--old algorithm uses index in els block
double track_iso_affable_old(int eidx_pri, int eidx_sec) {
  //input is index of both electrons in els block
  //return is the sum of the pt of all tracks in the range 0.015 < dR < 0.3 around the els, and this time, exclude cone around both els

  double isolation = 0;
  for( unsigned int i=0; i<cms2.trks_trk_p4().size(); i++ ) {
	//cuts on track quality
	if( cms2.trks_trk_p4()[i].pt() <= 1.0 )
	  continue;

	double dR1 = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4In()[eidx_pri], cms2.trks_trk_p4()[i] );
	double dR2 = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4In()[eidx_sec], cms2.trks_trk_p4()[i] );
	if( (unsigned int)eidx_sec > cms2.els_p4In().size() ) //keep redundant check--shouldn't appear
	  cout << "\teidx_sec out of bounds  " << eidx_sec << "  " << cms2.els_p4In().size() << "   " << dR2 << endl;
	if( dR1 > 0.015 && dR1 < 0.3 && dR2 > 0.015 )
	  isolation += cms2.trks_trk_p4()[i].pt();
  }

  return isolation;
}

//new affable algorithm corrects for all electrons passing certain cuts
//for now, these cuts are hard-coded to be same as what i did for ecal iso:
// pt > 15, H/E < 0.1
double track_iso_affable(int eidx_pri) {
  //input is index of both electrons in els block
  //return is the sum of the pt of all tracks in the range 0.015 < dR < 0.3 around the els, and this time, exclude cone around both els

  vector<LorentzVector> els_p4_cln;
  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	if( cms2.els_p4()[i].pt() > 15 && cms2.els_hOverE()[i] < 0.1 )
	  els_p4_cln.push_back( cms2.els_p4()[i] );
  }
  
  double isolation = 0;
  for( unsigned int i=0; i<cms2.trks_trk_p4().size(); i++ ) {
	//cuts on track quality
	if( cms2.trks_trk_p4()[i].pt() <= 1.0 )
	  continue;

	double dR1 = ROOT::Math::VectorUtil::DeltaR( cms2.els_p4In()[eidx_pri], cms2.trks_trk_p4()[i] );
	bool exclude = false;
	for( unsigned int j=0; j < els_p4_cln.size(); j++ ) {
	  double dR2 = ROOT::Math::VectorUtil::DeltaR( els_p4_cln[j], cms2.trks_trk_p4()[i] );

	  if( ((dR1 < 0.015 || dR2 < 0.015) && dR1 < 0.3) || dR1 > 0.3 ) {
		exclude = true;
		break;
	  }
	  //MAKE SURE THIS ALGO IS RIGHT
	}
	if( exclude ) continue;
	isolation += cms2.trks_trk_p4()[i].pt();
  }

  return isolation;
}


cuts_t Looper::EventSelect () {
     cuts_t ret = 0;
     return ret;
}


void Looper::FillEventHistos ()
{
  //cout << "start FillEventHistos\n";
  //if( useweight )
  weight = Weight(0);
  //else
  //weight = 1;

  //denomitr = 0;

  //for( int i=0;i<ncounts;i++ )
  //count[i].total += weight;

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
	//for( int i=0;i<ncounts;i++ ) count[i].gencuts += weight;
	return; //don't do anything if don't pass stat1cuts, just because i defined the denominator this way--should not bias results
  }
  //count[0].denom += weight;
  const double drstat1 = ROOT::Math::VectorUtil::DeltaR(cms2.genps_lepdaughter_p4()[idxlep1[0]], cms2.genps_lepdaughter_p4()[idxlep1[1]]) ;
  const double ptstat1[2] = { cms2.genps_lepdaughter_p4()[idxlep1[0]].pt(), cms2.genps_lepdaughter_p4()[idxlep1[1]].pt() };
  //bool drstat1_05_1 = false;
  e_njets->Fill( cms2.evt_njets(), weight );
  double eta0 = fabs( cms2.genps_lepdaughter_p4()[idxlep1[0]].eta() );
  double eta1 = fabs( cms2.genps_lepdaughter_p4()[idxlep1[1]].eta() );

  int etareg = -999; //bb = 0, be = 1, ee = 2
  if( eta0 < 1.479 && eta1 < 1.479 )
	etareg = 0;
  else if( (eta0 < 1.479 && eta1 > 1.479) || (eta0 > 1.479 && eta1 < 1.479 ) )
	etareg = 1;
  else if( eta0 > 1.479 && eta1 > 1.479 )
	etareg = 2;
  else
	cout << "\n\nBad etareg value. Check eta.\n\n";
  
  //fill denom histos here
  //individual--fill each
  for( int i=0;i<2;i++ ){
	eff_edr_hcal_iso[i]->denom->Fill( drstat1, weight );
	eff_edr_ecal_iso[i]->denom->Fill( drstat1, weight );
	eff_edr_ecal_iso_soft[i]->denom->Fill( drstat1, weight );
	eff_edr_ecal_iso_subt[i]->denom->Fill( drstat1, weight );
	eff_edr_trck_iso[i]->denom->Fill( drstat1, weight );
	eff_edr_trck_iso_recalc[i]->denom->Fill( drstat1, weight );
	eff_edr_trck_iso_affble[i]->denom->Fill( drstat1, weight );
	eff_edr_ecal_iso_reg[i][etareg]->denom->Fill( drstat1, weight ); 
  }
  eff_edr_hcal_iso_pair->denom->Fill( drstat1, weight );
  eff_edr_ecal_iso_pair->denom->Fill( drstat1, weight );
  eff_edr_trck_iso_pair->denom->Fill( drstat1, weight );
  eff_edr_trck_iso_recalc_pair->denom->Fill( drstat1, weight );
  eff_edr_trck_iso_affble_pair->denom->Fill( drstat1, weight );

  // dr matching
  //cout << "dr matching\n";
  //const double maxcone = 0.1;
  const double maxcone = 0.015;
  //const double maxcone_small = 0.015; //this is the inner cone for ele track iso in ww note
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
		//else
		  //count[0].dupematch += weight;
	  }
	}
  }

  //cout << "\nFill\n";
  //cuts_t e_reco_cuts[2] = {0,0};
  double hcaliso[2] = {0,0};
  double ecaliso[2] = {0,0};
  double ecaliso_subt[2] = {0,0};
  double trckiso[2] = {0,0};
  double trckiso_calc[2] = {0,0};
  double trckiso_affa[2] = {0,0};
  //fill iso histos only for those reco els which match
  for( int i=0;i<2;i++ ) {
	if( min_dr[i] < maxcone ) {
	  double pt  = cms2.els_p4().at(idxlepreco[i]).pt(); //reco pt
	  //return pt / (pt + sum + 1e-5); //from selections.cc
	  hcaliso[i] = pt/(pt + cms2.els_pat_hcalIso().at(idxlepreco[i]) );
	  ecaliso[i] = pt/(pt + cms2.els_pat_ecalIso().at(idxlepreco[i]) );
	  ecaliso_subt[i] = pt/(pt + max(cms2.els_pat_ecalIso().at(idxlepreco[i])-2.0, 0.0) );
	  trckiso[i] = pt/(pt + cms2.els_pat_trackIso().at(idxlepreco[i]) );
	  trckiso_calc[i] = pt/(pt + track_iso(idxlepreco[i]) );
	  trckiso_affa[i] = pt/(pt + track_iso_affable(idxlepreco[i]) ); //new one excludes all els (see above)
	  //This is for old algorithm only
	  //must check that the OTHER index is found, and in cone. If not, isolation is same as recalc
	  //if( idxlepreco[ i==0?1:0 ] > cms2.els_p4().size() || min_dr[ i==0?1:0 ] > maxcone ) {
	  //	//cout << "IDXLEPRECO OUT OF BOUNDS   " << idxlepreco[ i==0?1:0 ] << endl;
	  //	trckiso_affa[i] = trckiso_calc[i];
	  //}
	  //else { //index is ok
	  //	trckiso_affa[i] = pt/(pt + track_iso_affable_old(idxlepreco[i], idxlepreco[ i==0?1:0 ]) ); //old one just excludes dR matched
	  //}
	  
	  e_hcal_iso->Fill( hcaliso[i], weight );
	  e_ecal_iso->Fill( ecaliso[i], weight );
	  e_trck_iso->Fill( trckiso[i], weight );
	  e_trck_iso_recalc->Fill( trckiso_calc[i], weight );
	  e_trck_iso_affble->Fill( trckiso_affa[i], weight );
	  
	  if( drstat1 >= 0.5 && drstat1 <= 1.0 ) { //note greater/less equal
	  //drstat1_05_1 = true;
	  //if( drstat1_05_1 ) {
		e_hcal_iso_dr05_1->Fill( hcaliso[i], weight );
		e_ecal_iso_dr05_1->Fill( ecaliso[i], weight );
		e_trck_iso_dr05_1->Fill( trckiso[i], weight );
		e_trck_iso_recalc_dr05_1->Fill( trckiso_calc[i], weight );
		e_trck_iso_affble_dr05_1->Fill( trckiso_affa[i], weight );
		//e_reco_cuts[i] |= CUT_BIT(CUT_EL_DR);
	  }
	  //another cone size, but leave the first---all same cone now
	  //if( min_dr[i] < maxcone_small ) {
	  //e_trck_iso_match015->Fill( trckiso[i], weight );
	  //}
	}
  }

  //clean jets
  int jetsize = cms2.evt_njets();
  for( unsigned int i=0; i<2; i++ ) { //loop over 2 electrons
	if( min_dr[i] > maxcone ) continue; //only clean jets for dR matched els

	for( unsigned int j=0; j<cms2.jets_p4().size(); j++ ) { //using 0.4 to match b'c it's jet cone size
	  if( ROOT::Math::VectorUtil::DeltaR(cms2.els_p4()[idxlepreco[i]], cms2.jets_p4()[j]) < 0.4 ) { 
		jetsize--;
		break; //only allow each el to match to at most 1 jet
	  }
	}
  }
  e_njets_clean->Fill( jetsize, weight );

  //cuts_t iso_cut = ( CUT_BIT(CUT_EL_DR) | CUT_BIT(CUT_EL_ISO) ); 
  //if( (e_reco_cuts[0] & e_reco_cuts[1] & iso_cut) == iso_cut ) {
  //fill iso
  //define iso cuts: values taken from output of below
  double hcalcut = 0.985;
  double ecalcut = 0.95;
  double trckcut = 0.99;

  //njet bin efficiencies
  int jetbin = jetsize;
  if( jetsize >= (njetbins-1) ) jetbin = (njetbins-1); //this is maximum bin: njetbins = 5

  for( unsigned int i=0; i<2; i++ ) {
	eff_edr_ecal_iso_njet[i][jetbin]->denom->Fill( drstat1, weight );
	if( min_dr[i] < maxcone && ecaliso[i] > ecalcut ) { //numerator
	  eff_edr_ecal_iso_njet[i][jetbin]->numer->Fill( drstat1, weight );
	  eff_edr_ecal_iso_reg[i][etareg]->numer->Fill( drstat1, weight );
	}
  }

  //bools for selecting events from event display--just care about
  bool ftrcknumer = false; //first passed track numer
  bool fecalnumer = false; //first passed ecal numer
  bool strcknumer = false; //second passed track numer
  //bool secalmatch = false; //second passed dr match -- no need b'c already have if for it

  //INDIVIDUAL EFFICIENCY
  for( int i=0;i<2;i++ ) {
	if( min_dr[i] < maxcone ) {
	  if( hcaliso[i] > hcalcut ) {
		eff_edr_hcal_iso[i]->numer->Fill( drstat1, weight );
	  }

	  if( trckiso[i] > trckcut ) {
		eff_edr_trck_iso[i]->numer->Fill( drstat1, weight );
		//plot the ecal iso variable here (not as ratio)--for those that pass track to see if pedestal
		e_ecal_iso_bare[i]->Fill( cms2.els_pat_ecalIso().at(idxlepreco[i]), weight );
		if( i==0 ) ftrcknumer = true;
		else if( i == 1 ) strcknumer = true;
	  }

	  if( (ecaliso[i] > ecalcut && ptstat1[i] > 60) || (cms2.els_pat_ecalIso().at(idxlepreco[i]) < 3 && ptstat1[i] <= 60) ) {
		eff_edr_ecal_iso_soft[i]->numer->Fill( drstat1, weight );
	  }
	  
	  if( ecaliso_subt[i] > ecalcut ) {
		eff_edr_ecal_iso_subt[i]->numer->Fill( drstat1, weight );
	  }
	  
	  if( ecaliso[i] > ecalcut ) {
		eff_edr_ecal_iso[i]->numer->Fill( drstat1, weight );
		if( i==0 ) fecalnumer = true;
	  }
	  //first passed both ecal+trck, second passed trck, failed ecal
	  else if( i == 1 && ftrcknumer && strcknumer && fecalnumer ) {
		if( numprint < 200 ) {
		  //cout << SampleName().c_str() << "   " << cms2.evt_run() << "   "
		  //cout << cms2.evt_run() << "   "
		  //   << cms2.evt_event() << "   " << cms2.evt_lumiBlock() << endl;
		  //<< "   " << ecaliso[i] << endl;
		}
		numprint++;
	  }
	
	  if( trckiso_calc[i] > trckcut ) {
		eff_edr_trck_iso_recalc[i]->numer->Fill( drstat1, weight );
	  }

	  if( trckiso_affa[i] > trckcut ) {
		eff_edr_trck_iso_affble[i]->numer->Fill( drstat1, weight );
	  }
	}
  }

  //PAIR EFFICIENCY
  if( min_dr[0] < maxcone && min_dr[1] < maxcone ) {
	if( hcaliso[0] > hcalcut && hcaliso[1] > hcalcut ) {
	  eff_edr_hcal_iso_pair->numer->Fill( drstat1, weight );
	}

	if( ecaliso[0] > ecalcut && ecaliso[1] > ecalcut ) {
	  eff_edr_ecal_iso_pair->numer->Fill( drstat1, weight );
	}

	if( trckiso[0] > trckcut && trckiso[1] > trckcut ) {
	  eff_edr_trck_iso_pair->numer->Fill( drstat1, weight );
	}
	
	if( trckiso_calc[0] > trckcut && trckiso_calc[1] > trckcut ) {
	  eff_edr_trck_iso_recalc_pair->numer->Fill( drstat1, weight );
	}

	if( trckiso_affa[0] > trckcut && trckiso_affa[1] > trckcut ) {
	  eff_edr_trck_iso_affble_pair->numer->Fill( drstat1, weight );
	}
  }

}
// end void Looper::FillEventHistos ()


//function returns the BIN number of the bin which includes the 90% point
int get_90_bin(TH1F* hist ) {
  int j=0;
  double total = hist->GetBinContent(0); //underflow bin
  //arguments of Integral() includes overflow bin b'c '+1'                      
  double integral90pc = 0.1*( hist->Integral(0, hist->GetNbinsX()+1) );
  while( total < integral90pc ) {
    j++;
    total += hist->GetBinContent(j);
  }
  return j;
}
//end get_90_bin                                                                

double get_90_bin_highedge(TH1F* hist) {
  return hist->GetBinLowEdge( get_90_bin( hist ) + 1 ) ;
}


void Looper::End ()
{
  cout << "\n\nNumber lines passed print:  " << numprint << "\n\n";


  TCanvas *c1 = new TCanvas();// eff->GetName(), eff->GetName());

  //c1->SetLogy();
  gPad->SetLogy();
  e_hcal_iso->Draw(); c1->SaveAs((TString)e_hcal_iso->GetName()+".png");
  e_ecal_iso->Draw(); c1->SaveAs((TString)e_ecal_iso->GetName()+".png");
  e_trck_iso->Draw(); c1->SaveAs((TString)e_trck_iso->GetName()+".png");
  //e_trck_iso_match015->Draw(); c1->SaveAs((TString)e_trck_iso_match015->GetName()+".png");
  
  e_hcal_iso_dr05_1->Draw(); c1->SaveAs((TString)e_hcal_iso_dr05_1->GetName()+".png"); 
  e_ecal_iso_dr05_1->Draw(); c1->SaveAs((TString)e_ecal_iso_dr05_1->GetName()+".png");
  e_trck_iso_dr05_1->Draw(); c1->SaveAs((TString)e_trck_iso_dr05_1->GetName()+".png");

  e_trck_iso_recalc->Draw(); c1->SaveAs((TString)e_trck_iso_recalc->GetName()+".png");
  e_trck_iso_recalc_dr05_1->Draw(); c1->SaveAs((TString)e_trck_iso_recalc_dr05_1->GetName()+".png");
  e_trck_iso_affble->Draw(); c1->SaveAs((TString)e_trck_iso_affble->GetName()+".png");
  e_trck_iso_affble_dr05_1->Draw(); c1->SaveAs((TString)e_trck_iso_affble_dr05_1->GetName()+".png");

  e_njets->Draw(); c1->SaveAs((TString)e_njets->GetName()+".png");
  e_njets_clean->Draw(); c1->SaveAs((TString)e_njets_clean->GetName()+".png");

  e_ecal_iso_bare[0]->Draw(); c1->SaveAs((TString)e_ecal_iso_bare[0]->GetName()+".png");
  e_ecal_iso_bare[1]->Draw(); c1->SaveAs((TString)e_ecal_iso_bare[1]->GetName()+".png");

  //individual
  for( int i=0;i<2;i++ ) {
	eff_edr_hcal_iso[i]->MakeEff( );
	eff_edr_ecal_iso[i]->MakeEff( );
	eff_edr_ecal_iso_soft[i]->MakeEff();
	eff_edr_ecal_iso_subt[i]->MakeEff();
	eff_edr_trck_iso[i]->MakeEff( );
	eff_edr_trck_iso_recalc[i]->MakeEff( );
	eff_edr_trck_iso_affble[i]->MakeEff( );
  }
  //pair
  eff_edr_hcal_iso_pair->MakeEff( );
  eff_edr_ecal_iso_pair->MakeEff( );
  eff_edr_trck_iso_pair->MakeEff( );
  eff_edr_trck_iso_recalc_pair->MakeEff( );
  eff_edr_trck_iso_affble_pair->MakeEff( );

  //njet, etareg
  for( int i=0; i<2; i++ ) {
	for( int j=0; j<njetbins; j++ ) {
	  eff_edr_ecal_iso_njet[i][j]->MakeEff();
	}
	for( int j=0; j<etaregions; j++ ) {
	  eff_edr_ecal_iso_reg[i][j]->MakeEff();
	}
  }
	
  //saves overlayed with legend and stuff...see #included file
  over_save( e_ecal_iso_bare[0], e_ecal_iso_bare[1]);
  over_save( e_ecal_iso_bare[0], e_ecal_iso_bare[1], true, true);
  over_save( e_trck_iso, e_trck_iso_recalc, true ); 
  over_save( e_trck_iso_dr05_1, e_trck_iso_recalc_dr05_1, true );
  over_save( e_trck_iso_recalc, e_trck_iso_affble, true );
  over_save( e_trck_iso_recalc_dr05_1, e_trck_iso_affble_dr05_1, true );
  //over_save( e_trck_iso_affble_dr05_1, e_trck_iso_recalc_dr05_1, true ); //this one is duplicate of immediately above for check that it's right
  for( int i=0;i<2;i++ ) {
	over_save( eff_edr_trck_iso[i]->eff, eff_edr_trck_iso_affble[i]->eff );
	over_save( eff_edr_trck_iso_recalc[i]->eff, eff_edr_trck_iso_affble[i]->eff );
	over_save( eff_edr_ecal_iso[i]->eff, eff_edr_trck_iso[i]->eff );
	over_save( eff_edr_ecal_iso_soft[i]->eff, eff_edr_trck_iso[i]->eff );
	over_save( eff_edr_ecal_iso_subt[i]->eff, eff_edr_trck_iso[i]->eff );
	over_save( eff_edr_trck_iso[i]->eff, eff_edr_ecal_iso[i]->eff, eff_edr_hcal_iso[i]->eff, TString(Form("%s_e%i_drstat1_tri_iso_comp", SampleName().c_str(), i+1)) );
	over_save( eff_edr_ecal_iso_njet[i][1]->eff, eff_edr_ecal_iso_njet[i][2]->eff, eff_edr_ecal_iso_njet[i][3]->eff, eff_edr_ecal_iso_njet[i][4]->eff, TString(Form("%s_e%i_drstat1_ecal_iso_njet_all_comp", SampleName().c_str(), i+1)) );
	over_save( eff_edr_ecal_iso_reg[i][0]->eff, eff_edr_ecal_iso_reg[i][2]->eff );
	//over_save( eff_edr_trck_iso_dr05_1[i]->eff, eff_edr_trck_iso_affble_dr05_1[i]->eff);
	//over_save( eff_edr_trck_iso_recalc_dr05_1[i]->eff, eff_edr_trck_iso_affble_dr05_1[i]->eff);
  }
  over_save( eff_edr_ecal_iso_pair->eff, eff_edr_trck_iso_pair->eff );
  over_save( eff_edr_trck_iso_pair->eff, eff_edr_trck_iso_recalc_pair->eff );
  over_save( eff_edr_trck_iso_affble_pair->eff, eff_edr_trck_iso_pair->eff );
  over_save( eff_edr_trck_iso_recalc_pair->eff, eff_edr_trck_iso_affble_pair->eff );
  over_save( e_njets, e_njets_clean );

  //print 90% points--90% is below this point
  // use up (high) edge of bin which get_90_bin returns + 1 because i actually want up edge of this bin
  cout << "\n\n";
  cout << "hist\t\t\t90 bin low edge\t90 bin\n";
  cout << e_ecal_iso->GetName() << "   " << e_ecal_iso->GetBinLowEdge( get_90_bin( e_ecal_iso ) ) << "   " << get_90_bin( e_ecal_iso ) << endl; 
  cout << e_trck_iso->GetName() << "   " << e_trck_iso->GetBinLowEdge( get_90_bin( e_trck_iso ) ) << "   " << get_90_bin( e_trck_iso ) << endl;
  cout << e_hcal_iso->GetName() << "   " << e_hcal_iso->GetBinLowEdge( get_90_bin( e_hcal_iso ) ) << "   " << get_90_bin( e_hcal_iso ) << endl;
  cout << e_ecal_iso_dr05_1->GetName() << "   " << e_ecal_iso_dr05_1->GetBinLowEdge( get_90_bin(e_ecal_iso_dr05_1) ) << "   " << get_90_bin( e_ecal_iso_dr05_1 ) << endl;
  cout << e_trck_iso_dr05_1->GetName() << "   " << e_trck_iso_dr05_1->GetBinLowEdge( get_90_bin(e_trck_iso_dr05_1) ) << "   " << get_90_bin( e_trck_iso_dr05_1 ) << endl;
  cout << e_hcal_iso_dr05_1->GetName() << "   " << e_hcal_iso_dr05_1->GetBinLowEdge( get_90_bin(e_hcal_iso_dr05_1) ) << "   " << get_90_bin( e_hcal_iso_dr05_1 ) << endl;
  /*cout << "hist\t\t\t90 bin high edge\t90 bin\n";
  cout << e_ecal_iso->GetName() << "   " << get_90_bin_highedge( e_ecal_iso ) << "   " << get_90_bin( e_ecal_iso ) << endl; 
  cout << e_trck_iso->GetName() << "   " << get_90_bin_highedge( e_trck_iso ) << "   " << get_90_bin( e_trck_iso ) << endl;
  cout << e_hcal_iso->GetName() << "   " << get_90_bin_highedge( e_hcal_iso ) << "   " << get_90_bin( e_hcal_iso ) << endl;
  cout << e_ecal_iso_dr05_1->GetName() << "   " << get_90_bin_highedge( e_ecal_iso_dr05_1 ) << "   " << get_90_bin( e_ecal_iso_dr05_1 ) << endl;
  cout << e_trck_iso_dr05_1->GetName() << "   " << get_90_bin_highedge( e_trck_iso_dr05_1 ) << "   " << get_90_bin( e_trck_iso_dr05_1 ) << endl;
  cout << e_hcal_iso_dr05_1->GetName() << "   " << get_90_bin_highedge( e_hcal_iso_dr05_1 ) << "   " << get_90_bin( e_hcal_iso_dr05_1 ) << endl;*/
  cout << "\n\n";

  
  int ret = fprintf(logfile_, 
					"Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
					" Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
					sample_.name.c_str(),
					CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
					CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
					CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
					CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
					CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
	perror("writing to log file");
}


