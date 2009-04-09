#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"
//#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

using namespace std;


Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
  numprint = 0;
  count.total = 0;
  count.denom = 0;
  count.numer = 0;
  count.opp_sign = 0;
  count.same_flv = 0;
  count.pt20  = 0;
  count.el_good = 0;
  count.bad_mom = 0;
  count.no_match_z = 0;
  // zero out the candidate counters (don't comment this out)
  memset(cands_passing_	, 0, sizeof(cands_passing_       ));
  memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
  memset(cands_count_		, 0, sizeof(cands_count_         ));

  //initialize cuts vector for N-1 hist: THESE ARE NUMERATOR CUTS
  els_cuts.push_back( CUT_BIT(CUT_OPP_SIGN) );
  els_cuts.push_back( CUT_BIT(CUT_ETA25) );
  els_cuts.push_back( CUT_BIT(CUT_PT20) );
  els_cuts.push_back( CUT_BIT(CUT_IN_Z_WINDOW) );
  els_cuts.push_back( CUT_BIT(CUT_MC_EL) );
  els_cuts.push_back( CUT_BIT(CUT_EL_GOOD) );
  els_cuts.push_back( CUT_BIT(CUT_EL_ISO) );
  //els_cuts.push_back( CUT_BIT(CUT_SAME_FLV) );

  els_cuts_names.push_back( "_opp_sign" );
  els_cuts_names.push_back( "_eta_2.5" );
  els_cuts_names.push_back( "_pt20" );
  els_cuts_names.push_back( "_z_window" );
  els_cuts_names.push_back( "_mc_truth" );
  els_cuts_names.push_back( "_elgood" );
  els_cuts_names.push_back( "_eliso" );
}

void Looper::BookHistos ()
{
  hgen_z_mass = new TH1F( Form("%s_Gen_Z_Mass", SampleName().c_str()), "Gen_Z_Mass", 100, 0, 200 );
  hgen_z_p = new TH1F( Form("%s_Gen_Z_Momentum", SampleName().c_str()), "Gen_Z_Momentum", 80, 0, 800 );
  hgen_z_pt = new TH1F( Form("%s_Gen_Z_Pt", SampleName().c_str()), "Gen_Z_Pt", 80, 0, 800 );
  hgen_z_eta = new TH1F( Form("%s_Gen_Z_Eta", SampleName().c_str()), "Gen_Z_Eta", 50, -5, 5 );

  hgen1_lep_mass = new TH1F( Form("%s_Gen1_Lep_Z_Mass", SampleName().c_str()), "Gen1_Lep_Z_Mass", 100, 0, 200 );
  
  for( int i = 0; i < nlepplots; i++) {
	hgen1_lep_pt[i] = new TH1F( Form("%s_Gen1_Lep_Z_Pt%i", SampleName().c_str(),i), "Gen1_Lep_Z_Pt", 100, 0, 200 );
	hgen1_lep_eta[i] = new TH1F( Form("%s_Gen1_Lep_Z_Eta%i", SampleName().c_str(),i), "Gen1_Lep_Z_Eta", 50, -5, 5 );
	hgen1_lep_pt[i]->Sumw2();
	hgen1_lep_eta[i]->Sumw2();
  }
  
  hels_size = new TH1F( Form("%s_els_size", SampleName().c_str()), "els_size", 10, 0, 10 );
  hels_iso = new TH1F( Form("%s_els_iso", SampleName().c_str()), "els_iso", 50, 0, 1 );
  heff_iso = new TH1F( Form("%s_eff_iso", SampleName().c_str()), "iso_eff", 80, 0, 800 );
  //heff_mc1match = new TH1F( Form("%s_MC1match_Z_Efficiency", SampleName().c_str()), "MC1match_Z_Efficiency", 50, 0, 500 );
  //heff_mc3match = new TH1F( Form("%s_MC3match_Z_Efficiency", SampleName().c_str()), "MC3match_Z_Efficiency", 50, 0, 500 );

  //CHANGE CUT VECTOR (NEW ONE) FOR MUS
  //cut vectors for nminus1hist
  heff_p_N1 = new NMinus1Hist(sample_, "Z_Efficiency_vZp;", 80, 0, 800, el_base, els_cuts, els_cuts_names );
  //heff_pm_N1 = new NMinus1Hist(sample_, "Z_Efficiency_vZp_mus;", 80, 0, 800, el_base, els_cuts, els_cuts_names );
  heff_pt_N1 = new NMinus1Hist(sample_, "Z_Efficiency_vZpt;", 80, 0, 800, el_base, els_cuts, els_cuts_names );
  //heff_ptm_N1 = new NMinus1Hist(sample_, "Z_Efficiency_vZpt_mus;", 80, 0, 800, el_base, els_cuts, els_cuts_names );

  
  hgen_z_mass->Sumw2();
  hgen_z_p->Sumw2();
  hgen_z_pt->Sumw2();
  hgen_z_eta->Sumw2();
  hgen1_lep_mass->Sumw2();
  hels_size->Sumw2();
  hels_iso->Sumw2();
  heff_iso->Sumw2();
  //heff_mc1match->Sumw2();
  //heff_mc3match->Sumw2();
  
  //hgen_z_mass->SetFillColor(sample_.histo_color);
  //hgen_z_mass->SetLineColor(sample_.histo_color);
  //hgen_z_mass->SetMarkerColor(sample_.histo_color);
  //hgen_z_p->SetFillColor(sample_.histo_color);
  //hgen_z_p->SetLineColor(sample_.histo_color);
  //hgen_z_p->SetMarkerColor(sample_.histo_color);
  //hels_size->SetFillColor(sample_.histo_color);
  //hels_size->SetLineColor(sample_.histo_color);
  //hels_size->SetMarkerColor(sample_.histo_color);
  //hgen1lep_z_mass->SetFillColor(sample_.histo_color); // this is annoying...
  //hgen1lep_z_mass->SetLineColor(sample_.histo_color);
  //hgen1lep_z_mass->SetMarkerColor(sample_.histo_color);
  //heff_mc1match->SetFillColor(sample_.histo_color);
  //heff_mc1match->SetLineColor(sample_.histo_color);
  //heff_mc1match->SetMarkerColor(sample_.histo_color);
  //heff_mc3match->SetFillColor(sample_.histo_color);
  //heff_mc3match->SetLineColor(sample_.histo_color);
  //heff_mc3match->SetMarkerColor(sample_.histo_color);
}


bool Looper::FilterEvent()
{ 

  return false; 
}

cuts_t Looper::Stat1Select(vector<int> idx) {
  //this function should check all cuts in cuts_t stat1cuts (Looper.h)
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
	//if( numprint < 30 )
	cout << "mother mismatch: " << "" << "\n";
	return ret;
  }
	 
  if( cms2.genps_lepdaughter_p4()[idx[0]].eta() < 2.5 &&
	  cms2.genps_lepdaughter_p4()[idx[1]].eta() < 2.5 )
	ret |= CUT_BIT(CUT_ETA25);

  double mass = ( cms2.genps_lepdaughter_p4()[idx[0]] + cms2.genps_lepdaughter_p4()[idx[1]] ).M();
  if( mass > 76. && mass < 106. )
	ret |= CUT_BIT(CUT_IN_Z_WINDOW);

  return ret;
}
//end Looper::Stat1Select
 
//cuts on numerator (reco) leptons--just electrons now!
//2nd arg: 0=el, 1=mu
cuts_t Looper::LepSelect(int i, int flv) {
  cuts_t ret = 0;

  if( flv == 0 ) { //electron by reco
	if( abs(cms2.els_mc_id()[i]) == 11 ) { //mc truth
	  ret |= (CUT_BIT(CUT_MC_EL));
	  if( cms2.els_mc_motherid()[i] != 23 ) {
		count.bad_mom += weight;
	  }
	}
	
	if( cms2.els_p4()[i].eta() < 2.5 )
	  ret |= CUT_BIT(CUT_ETA25);
  
	if( cms2.els_p4()[i].pt() > 20 )
	  ret |= (CUT_BIT(CUT_PT20));
	else
	  count.pt20 += weight;

	//if( goodElectronWithoutIsolation(cms2.els_mcidx()[i]) ) //replaced below
	if( cms2.els_tightId22XMinMatteo().at(i) == 1 )
	  //cms2.els_closestMuon().at(i) == -1 &&
	  //TMath::Abs(cms2.els_d0corr().at(i)) < 0.025 )
	  ret |= CUT_BIT(CUT_EL_GOOD);
	else
	  count.el_good += weight;

	if( passElectronIsolation(i, true) )
	  ret |= CUT_BIT(CUT_EL_ISO);
  }
  else
	cout << "Bad flv value in LepSelect\n";

  return ret;
}
//end Looper::LepSelect


//cuts on pair of leptons--maybe should be DilepSelect
cuts_t Looper::PairSelect(cuts_t numer, cuts_t denom) {
  //NEVER DO cut1 & CUT_BIT(CUT) > 0 because the first bit is the sign!!!
  cuts_t ret = 0, lep1 = 0, lep2 = 0;
  int passdenom = 0, idx[2];

  for( unsigned int i=0; i<cms2.els_p4().size(); i++) {
	lep1 = LepSelect(i, 0);
	//cout << "done first lepselect\n";
	for( unsigned int j = i+1; j<cms2.els_p4().size(); j++) {
	  lep2 = LepSelect(j, 0);
	  //cout << "done second lepselect. size = " <<cms2.els_p4().size()<< "\n";
	  if( cms2.els_charge()[i] == -1*cms2.els_charge()[j] ) {//opp sign
		lep1 |= (CUT_BIT(CUT_OPP_SIGN)); //both just for convenience
		lep2 |= (CUT_BIT(CUT_OPP_SIGN));		
	  }

	  double mass = (cms2.els_p4()[i] + cms2.els_p4()[j]).M();
	  if( mass > 76. && mass < 106. ) {
		lep1 |= (CUT_BIT(CUT_IN_Z_WINDOW));
		lep2 |= (CUT_BIT(CUT_IN_Z_WINDOW));
	  }	  

	  if( lep1 & lep2 & denom ) {
		ret = lep1 & lep2;
		idx[0] = i;
		idx[1] = j;
		passdenom++;
	  }
	}
  }

  if( (ret & CUT_BIT(CUT_OPP_SIGN)) == 0 ) //NEED EXTRA PARENS FOR PRECENDNCE
	count.opp_sign += weight;

  if( (ret & denom) == denom ) {
	count.denom += weight;
	for( int i = 0; i < 2; i++ )
	  hels_iso->Fill( el_rel_iso(idx[i], true), weight );
  }

  if( (ret & numer) == numer )
	count.numer += weight;
  
  return ret;
}
//end Looper::PairSelect

// keep track of denominator cuts--todo???
cuts_t Looper::DenSelect(int i) {
  cuts_t ret = 0;

  return ret;
}

cuts_t Looper::EventSelect () {
     // In an event-based analysis, you would make your cuts here
     cuts_t ret = 0;
     return ret;
}

void Looper::FillStat1Histos( int zidx, vector<int> idx) {
  //moved z_pt and z_p to end of FillEventHistos b'c need to pass denom cuts
  hgen_z_eta->Fill( cms2.genps_p4()[zidx].eta(), weight );
  //size of idx already checked in Stat1Select
  hgen1_lep_mass->Fill( (cms2.genps_lepdaughter_p4()[idx[0]] + cms2.genps_lepdaughter_p4()[idx[1]]).M() , weight );
  double pt1 = cms2.genps_lepdaughter_p4()[idx[0]].pt();
  double pt2 = cms2.genps_lepdaughter_p4()[idx[1]].pt();
  if( pt1 > pt2 ) {
	hgen1_lep_pt[0]->Fill( pt1 );
	hgen1_lep_eta[0]->Fill(cms2.genps_lepdaughter_p4()[idx[0]].eta(), weight );
	hgen1_lep_pt[1]->Fill( pt2 );
	hgen1_lep_eta[1]->Fill(cms2.genps_lepdaughter_p4()[idx[1]].eta(), weight );
  }
  else {
	hgen1_lep_pt[0]->Fill( pt2 );
	hgen1_lep_eta[0]->Fill(cms2.genps_lepdaughter_p4()[idx[1]].eta(), weight );
	hgen1_lep_pt[1]->Fill( pt1 );
	hgen1_lep_eta[1]->Fill(cms2.genps_lepdaughter_p4()[idx[0]].eta(), weight );
  }
	
}
//end Looper::FillStat1Histos() 

void Looper::FillEventHistos ()
{
  //cout << "\nstart\n";
  //see 'notes' for definitions of efficiencies
  weight = Weight(0);

  cuts_t numer = el_base; //everything
  cuts_t denom = el_base & (~CUT_BIT(CUT_EL_ISO)); //switch off iso

  //way i did lorentz in susy:
  //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep1;
  //ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > lep2;
  //lep1.SetPx(0); lep1.SetPy(0); lep1.SetPz(0); lep1.SetE(0);
  //lep2.SetPx(0); lep2.SetPy(0); lep2.SetPz(0); lep2.SetE(0);
  TLorentzVector lep1(0,0,0,0);
  TLorentzVector lep2(0,0,0,0);
  //double zp = 0, zpt = 0, zeta = 0;
  int zidx = 0;
  int nstat3lep = 0;
  
  count.total += weight;
  //get status 3 Z, leps by looping on genps status 3 block
  for(unsigned int i=0; i<cms2.genps_id().size(); i++) {
	if( cms2.genps_id()[i] == 23 ) {
	  hgen_z_mass->Fill( cms2.genps_p4()[i].M(), weight );
	  //zpt = cms2.genps_p4()[i].pt();
	  //zp = cms2.genps_p4()[i].P();
	  //zeta = cms2.genps_p4()[i].eta();
	  zidx = i;
	}
	else if( abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13 ){ 
	  nstat3lep++;
	}
  }

  //if( zp == 0 || zpt == 0 ) {
  if( zidx == 0 ) {
	cout << "\nPROBLEM WITH GEN Z\n\n"; return;
  }
  else if( nstat3lep != 2 )
	cout << "number status 3 leptons != 2 : " << nstat3lep << endl;
  
  //get status 1 leps using lepdaughter--for denominator
  vector<int> idxlep;
  for(unsigned int i=0; i<cms2.genps_lepdaughter_id().size(); i++) {
	if( abs(cms2.genps_lepdaughter_id()[i]) == 11 ||
		abs(cms2.genps_lepdaughter_id()[i]) == 13  ) 
	  idxlep.push_back(i);	  
  }
  if( Stat1Select(idxlep) == stat1cuts ) {
	FillStat1Histos( zidx, idxlep );
  }
	
  //get reco els
  //cout << "reco\n";
  //bool have1 = false; bool have2 = false;
  //vector<int> idx1, idx3;
  //vector<cuts_t> cuts_passed;
  hels_size->Fill( cms2.els_p4().size() , weight );
  cuts_t pair = PairSelect(numer, denom);

  //cout << "middle\n";
  //check cuts, fill histos
  //if( cuts_passed1.size() == 2 && idx1.size() == 2 ) {
  if( (denom & pair) == denom ) {
	hgen_z_p->Fill( cms2.genps_p4()[zidx].P() , weight );
	hgen_z_pt->Fill( cms2.genps_p4()[zidx].pt(), weight );
  }
  if( (numer & pair) == numer ) {
	heff_iso->Fill( cms2.genps_p4()[zidx].P(), weight);

  }
  heff_p_N1->Fill(pair, DILEPTON_EE, cms2.genps_p4()[zidx].P(), weight);
  heff_pt_N1->Fill(pair, DILEPTON_EE, cms2.genps_p4()[zidx].pt(), weight);
  //heff_mc1match_N1->Fill(pair, DILEPTON_EE, zp, weight);

  //cout << "end\n";
  //else
  //count.no_match_z += weight;
  //if( cuts_passed3.size() == 2 ) 
  //heff_mc3match_N1->Fill(cuts_passed3[0]&cuts_passed3[1], (DileptonHypType)0, zp, weight);
}

void Looper::End ()
{
  //divide by the z momentum to get efficiency
  heff_iso->Divide( hgen_z_p );
  //heff_mc1match->Divide( hgen_z_p );
  //heff_mc3match->Divide( hgen_z_p );

  //heff_mc1match_N1->N()[DILEPTON_EE].Divide( hgen_z_p );
  //DileptonHypType:DILEPTON_ALL, DILEPTON_MUMU, DILEPTON_EMU, DILEPTON_EE
  for( unsigned int j=0; j<4; j++ ) {
	//heff_mc1match_N1->N()[(DileptonHypType)j].Divide( hgen_z_p );
	//heff_mc3match_N1->N()[(DileptonHypType)j].Divide( hgen_z_p );
	heff_p_N1->N()[(DileptonHypType)j].Divide( hgen_z_p );
	heff_pt_N1->N()[(DileptonHypType)j].Divide( hgen_z_pt );
	for( unsigned int i=0;i<els_cuts.size(); i++ ) {
	  //heff_mc1match_N1->NMinus1(i)[(DileptonHypType)j].Divide( hgen_z_p );
	  //heff_mc3match_N1->NMinus1(i)[(DileptonHypType)j].Divide( hgen_z_p );
	  heff_p_N1->NMinus1(i)[(DileptonHypType)j].Divide( hgen_z_p );
	  heff_pt_N1->NMinus1(i)[(DileptonHypType)j].Divide( hgen_z_pt );
	}
  }

  //print struct count content
  cout << "\nTotal events run on: " << count.total
	   << "\ntotal denominator events: " << count.denom
	   << "\ntotal numerator events: " << count.numer
	   << "\nnum events which fail opp sign: " << count.opp_sign
	   << "\nnum events which fail same flv: " << count.same_flv
	   << "\nnum events which fail pt > 20: " << count.pt20
	   << "\nnum events which fail good electron: " << count.el_good
	   << "\nnum electrons with mc_motherid != 23: " << count.bad_mom 
	   << "\nnum events with (num leptons with mc_motherid == 23) != 2 : " << count.no_match_z
	   << endl;

  //Example status message at the end of a looper; edit for your application
  
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




///////////////////////////////////////////////////////////////
// DEAD CODE

/*
  was in stat 3 loop for gen3 mass
	else if( abs(cms2.genps_id()[i]) == 11 || abs(cms2.genps_id()[i]) == 13 ){ 
	  nstat3lep++;
	  if( lep1.Px() == 0 ){
		lep1.SetPx( cms2.genps_p4()[i].Px() );
		lep1.SetPy( cms2.genps_p4()[i].Py() );
		lep1.SetPz( cms2.genps_p4()[i].Pz() );
		lep1.SetE(  cms2.genps_p4()[i].E() );
	  }
	  else if( lep2.Px() == 0 ) {
		lep2.SetPx( cms2.genps_p4()[i].Px() );
		lep2.SetPy( cms2.genps_p4()[i].Py() );
		lep2.SetPz( cms2.genps_p4()[i].Pz() );
		lep2.SetE(  cms2.genps_p4()[i].E() );
	  }
	}
  if( (lep1+lep2).M() > 0 )
	hgen3lep_z_mass->Fill( (lep1+lep2).M(), weight );
*/

//was before stat1 loop
  //if( cms2.evt_event() < 25 )
  //if( numprint < 30 ) {
	//cout << "evt:" << cms2.evt_event() << " lepdaughter size:" << cms2.genps_lepdaughter_id().size() << endl;
	//numprint++;
  //}

//was in stat1 loop
	//numlep++;
	//stopping point: fill the histos
	//if( !(abs(cms2.genps_lepdaughter_p4()[i].eta()) < 2.5) ) {
	  //if( cms2.evt_event() < 5 )
	  //if( numprint < 30 ) cout << "outside eta range\n";
	  //continue;
	//}

	//if( cms2.genps_lepdaughter_p4()[i].pt() < 20. )
	//continue;

	//matchlep.push_back(true);
/*
  if( numlep < 2 ) {
	if( numprint < 30 ) {
	  cout << "less than two leptons in event at status 1\n";
	  for(unsigned int i=0; i<10; i++) {
		cout << i << " " << cms2.genps_id()[i] << " " << endl;
	  }
	  for(unsigned int i=0; i<cms2.genps_lepdaughter_id().size(); i++) {
		cout << cms2.genps_lepdaughter_idx()[i] << " "
			 << cms2.genps_lepdaughter_id()[i] << " "
			 << cms2.genps_lepdaughter_p4()[i].pt() << endl;
	  }
	}
	numprint++;
  }
*/
/* from reco section
  for( unsigned int i=0; i<cms2.els_mc3idx().size(); i++) {
	if( cms2.els_mc3_motherid()[i] == 23 ) {
	  idx3.push_back(i);
	  //cuts_passed3.push_back( LepSelect(i, 0) );
	}
  }
*/
//also  from reco section
  //if( have1 && have2 )
  //numerator of mc1(3)match is reco matches two same flavor opp sign leptons
  //if( idx1.size() == 2 && cms2.els_mc_id()[idx1[0]] == -1*cms2.els_mc_id()[idx1[1]] )
  //	heff_mc1match->Fill( zp, weight );
  //if( idx3.size() == 2 && cms2.els_mc3_id()[idx3[0]] == -1*cms2.els_mc3_id()[idx3[1]] )
  // heff_mc3match->Fill( zp, weight );


//was in pair select
  //currently, cuts_passed.size is required to be two before this is called
//ret |= (cuts_passed[0] & cuts_passed[1]);

  //if( status == 1 ) { //useless?
  //if( cuts_passed[0] & cuts_passed[1] & (CUT_BIT(CUT_IS_EL)) ) { //same flv
  //}
  //else if( cuts_passed[0] & cuts_passed[1] & (CUT_BIT(CUT_IS_MU)) ) {
  //}
  //else
  //count.same_flv += weight;
  
  //if( cms2.els_mc_id()[i1] == -1* cms2.els_mc_id()[i2] ) //opp sign
