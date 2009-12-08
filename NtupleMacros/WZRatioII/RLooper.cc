#include <math.h>
#include "TVector3.h"
//#include "Math/VectorUtil.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "RLooper.h"


RLooper::RLooper (Sample s, cuts_t wcuts, const char *fname, cuts_t zcuts) 
     : LooperBase(s, wcuts, fname)
{
  // zero out the candidate counters (don't comment this out)
  memset(dcands_passing_	    , 0, sizeof(dcands_passing_       ));
  memset(dcands_passing_w2_	, 0, sizeof(dcands_passing_w2_    ));
  memset(dcands_count_		, 0, sizeof(dcands_count_         ));
  memset(scands_passing_	    , 0, sizeof(scands_passing_       ));
  memset(scands_passing_w2_	, 0, sizeof(scands_passing_w2_    ));
  memset(scands_count_		, 0, sizeof(scands_count_         ));

  //initialize data members
  transmass = 0;
  njets_20 = 0;
  njets_30 = 0;
  dil_njets_20 = 0;
  dil_njets_30 = 0;
  elidxs[0] = -1;
  elidxs[1] = -1;
  muidxs[0] = -1;
  muidxs[1] = -1;
  mu_tot = 0;
  mu_multi = 0; //more than 1 mu, but 1 selected
  mu_mulma = 0; //more than 1, and makes z mass with selected

  isssignal_ = false; //default is bkg for both
  isdsignal_ = false;

  wcuts_ = wcuts; //in this looper, these aren't used. prep for RLooper
  zcuts_ = zcuts; 
}


void RLooper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
}


//to add: mass, transverse mass, njets
void RLooper::BookHistos ()
{
  //double metmax = 80.;
  //int metbins = 80;
  //double d0max = 0.05;
  //double d0sigmax = 10.;

  // single lepton histograms (two + 1 types)
  for (unsigned int i = 0; i < 3; ++i) {
	std::string hyp = "e";
	if (i == 1) hyp = "m";
	else if (i == 2) hyp = "all";
        
	NewHist( hlep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
  }

  // di-lepton histograms (three + 1 types)
  for (unsigned int i = 0; i < 4; ++i) {
	NewHist( hdilep_0_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]), "dilep_0_pt", 100, 0.0, 100.0);
	NewHist( hdilep_1_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]), "dilep_1_pt", 100, 0.0, 100.0);
  }

}

cuts_t RLooper::DilepSelect(const enum DileptonHypType myType, int idx1, int idx2) {
  cuts_t ret = 0; 
  //float ptcut = 20.0; //since many bits, do manually
  LorentzVector lep1_p4, lep2_p4;	

  //if( elidxs[0] != -1 && elidxs[1] != -1 ) {
  if( myType == DILEPTON_EE ) {

	lep1_p4 = cms2.els_p4()[idx1];
	lep2_p4 = cms2.els_p4()[idx2];

	//20,20 -- nominal nominal
	if( cms2.els_p4()[idx1].pt() >= 20. && cms2.els_p4()[idx2].pt() >= 20. )
	  ret |= CUT_BIT(DILEP_PT);
	//10,20 -- loose nominal
	if( (cms2.els_p4()[idx1].pt() >= 10. && cms2.els_p4()[idx2].pt() >= 20.) ||
		(cms2.els_p4()[idx1].pt() >= 20. && cms2.els_p4()[idx2].pt() >= 10.) )
	  ret |= CUT_BIT(DILEP_PT_1020);
	//10,30 -- loose tightest
	if( (cms2.els_p4()[idx1].pt() >= 10. && cms2.els_p4()[idx2].pt() >= 30.) ||
		(cms2.els_p4()[idx1].pt() >= 30. && cms2.els_p4()[idx2].pt() >= 10.) )
  	  ret |= CUT_BIT(DILEP_PT_1030);
	//20,25 -- nominal tighter
	if( (cms2.els_p4()[idx1].pt() >= 25. && cms2.els_p4()[idx2].pt() >= 20.) ||
		(cms2.els_p4()[idx1].pt() >= 20. && cms2.els_p4()[idx2].pt() >= 25.) )
	  ret |= CUT_BIT(DILEP_PT_2025);
	//20,30 -- nominal tightest
	if( (cms2.els_p4()[idx1].pt() >= 30. && cms2.els_p4()[idx2].pt() >= 20.) ||
		(cms2.els_p4()[idx1].pt() >= 20. && cms2.els_p4()[idx2].pt() >= 30.) )
	  ret |= CUT_BIT(DILEP_PT_2030); 

	if( cms2.els_charge()[idx1] * cms2.els_charge()[idx2] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.els_p4()[idx1] + cms2.els_p4()[idx2]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);
	if( mass > 71. && mass < 111. )
	  ret |= CUT_BIT(DILEP_MASS_WIDE);

  }
  //else if( muidxs[0] != -1 && muidxs[1] != -1 ) {
  else if( myType == DILEPTON_MUMU ) {

	lep1_p4 = cms2.mus_p4()[idx1];
	lep2_p4 = cms2.mus_p4()[idx2];

	//20,20 -- nominal nominal
	if( cms2.mus_p4()[idx1].pt() >= 20. && cms2.mus_p4()[idx2].pt() >= 20.)
	  ret |= CUT_BIT(DILEP_PT);
	//10,20 -- loose nominal
	if( (cms2.mus_p4()[idx1].pt() >= 10. && cms2.mus_p4()[idx2].pt() >= 20.) ||
		(cms2.mus_p4()[idx1].pt() >= 20. && cms2.mus_p4()[idx2].pt() >= 10.) )
	  ret |= CUT_BIT(DILEP_PT_1020);
	//10,30 -- loose tightest
	if( (cms2.mus_p4()[idx1].pt() >= 10. && cms2.mus_p4()[idx2].pt() >= 30.) ||
		(cms2.mus_p4()[idx1].pt() >= 30. && cms2.mus_p4()[idx2].pt() >= 10.) )
  	  ret |= CUT_BIT(DILEP_PT_1030);
	//20,25 -- nominal tighter
	if( (cms2.mus_p4()[idx1].pt() >= 25. && cms2.mus_p4()[idx2].pt() >= 20.) ||
		(cms2.mus_p4()[idx1].pt() >= 20. && cms2.mus_p4()[idx2].pt() >= 25.) )
	  ret |= CUT_BIT(DILEP_PT_2025);
	//20,30 -- nominal tightest
	if( (cms2.mus_p4()[idx1].pt() >= 30. && cms2.mus_p4()[idx2].pt() >= 20.) ||
		(cms2.mus_p4()[idx1].pt() >= 20. && cms2.mus_p4()[idx2].pt() >= 30.) )
	  ret |= CUT_BIT(DILEP_PT_2030); 

	if( cms2.mus_charge()[idx1] * cms2.mus_charge()[idx2] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.mus_p4()[idx1] + cms2.mus_p4()[idx2]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);
	if( mass > 71. && mass < 111. )
	  ret |= CUT_BIT(DILEP_MASS_WIDE);

  }
  else
	cout << "BAD DILEPSELECT CALL" << endl;

  //jet vars--dilep case
  dil_njets_20 = 0;
  dil_njets_30 = 0;
  //this code ripped from selections.cc->getCaloJets
  for( unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
	if( ROOT::Math::VectorUtil::DeltaR(lep1_p4, cms2.jets_p4()[jj]) < 0.4 ||
		ROOT::Math::VectorUtil::DeltaR(lep2_p4, cms2.jets_p4()[jj]) < 0.4 )
	  continue;
	if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4)
	  continue;
	//count
	if (cms2.jets_p4()[jj].pt() > 20)
	  dil_njets_20++;
	if (cms2.jets_p4()[jj].pt() > 30)
	  dil_njets_30++;
  }
  if( dil_njets_20 == 0 )
	ret |= CUT_BIT(JET_VETO_20);
  if( dil_njets_30 == 0 )
	ret |= CUT_BIT(JET_VETO_30);

  return ret;
}

cuts_t RLooper::LepSelect(int lep_type, int i)
{
  cuts_t ret = 0;
  //float ptcut = 20.0; //do manually
  double vtx_d0 = -1.;
  LorentzVector lep_p4;
    
  // e
  if( lep_type == 0 ) {

	lep_p4 = cms2.els_p4()[i];
	//will hve to figure out which vtx this is--cut on delta z0
 	vtx_d0 = getCorrd0( cms2.vtxs_position()[0], cms2.els_p4()[i], cms2.els_vertex_p4()[i]);

	if( cms2.els_p4()[i].pt() >= 20. ) {
	  ret |= CUT_BIT(LEP_PT);    //i know i said don't use in RLooper, but i'll set it to be safe
	  ret |= CUT_BIT(LEP_PT_20);
	}
	if( cms2.els_p4()[i].pt() >= 10. )
	  ret |= CUT_BIT(LEP_PT_10);
	if( cms2.els_p4()[i].pt() >= 25. )
	  ret |= CUT_BIT(LEP_PT_25);
	if( cms2.els_p4()[i].pt() >= 30. )
	  ret |= CUT_BIT(LEP_PT_30);

	if( GoodSusyElectronWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyElectronIsolation(i, true) ) //bool is for use calo
	  ret |= CUT_BIT(LEP_ISO);

	//put in all cuts from GoodSusyElectronWithoutIsolation
	//if ( cms2.els_egamma_tightId().at(index)     !=  1) return false; 
	//if ( fabs(cms2.els_d0corr().at(index)) >= 0.02)   return false;
	//if ( cms2.els_closestMuon().at(index) != -1) return false; 
	//if ( TMath::Abs(cms2.els_p4()[index].eta()) > 2.4) return false;

	if ( cms2.els_egamma_tightId().at(i) ==  1
		 && cms2.els_closestMuon().at(i) == -1
		 && TMath::Abs(cms2.els_p4()[i].eta()) < 2.4 )
	  ret |= CUT_BIT(LEP_GOOD_NOD0);

	if ( fabs(cms2.els_d0corr().at(i)) <= 0.02)
	  ret |= CUT_BIT(LEP_D0);
  }

  // m
  else if( lep_type == 1 ) {

	lep_p4 = cms2.mus_p4()[i];
	//will hve to figure out which vtx this is--cut on delta z0
 	vtx_d0 = getCorrd0( cms2.vtxs_position()[0], cms2.mus_p4()[i], cms2.mus_vertex_p4()[i]);

	if( cms2.mus_p4()[i].pt() >= 20. ) {
	  ret |= CUT_BIT(LEP_PT);    //i know i said don't use in RLooper, but i'll set it to be safe
	  ret |= CUT_BIT(LEP_PT_20);
	}
	if( cms2.mus_p4()[i].pt() >= 10. )
	  ret |= CUT_BIT(LEP_PT_10);
	if( cms2.mus_p4()[i].pt() >= 25. )
	  ret |= CUT_BIT(LEP_PT_25);
	if( cms2.mus_p4()[i].pt() >= 30. )
	  ret |= CUT_BIT(LEP_PT_30);

	if( GoodSusyMuonWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyMuonIsolation(i) )
	  ret |= CUT_BIT(LEP_ISO);

	//put in all cuts from GoodSusyMuonWithoutIsolation
	//if (((cms2.mus_type().at(index)) & (1<<1)) == 0) return false; // global muon
	//if (((cms2.mus_type().at(index)) & (1<<2)) == 0) return false; // tracker muon
	//if (cms2.mus_validHits().at(index) < 11)    return false; 
	//if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; 
	//if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) return false;
	//if (cms2.mus_pat_ecalvetoDep().at(index) >= 4) return false; // ECalE < 4 
	//if (cms2.mus_pat_hcalvetoDep().at(index) >= 6) return false; // HCalE < 6 
	//if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) return false;
	
	if( (cms2.mus_type().at(i) & (1<<1)) == (1<<1) // global muon
		&& (cms2.mus_type().at(i) & (1<<2)) == (1<<2) // tracker muon
		&& cms2.mus_validHits().at(i) >= 11
		&& cms2.mus_gfit_chi2().at(i)/cms2.mus_gfit_ndof().at(i) < 10
		&& cms2.mus_pat_ecalvetoDep().at(i) < 4 // ECalE < 4 
		&& cms2.mus_pat_hcalvetoDep().at(i) < 6 // HCalE < 6 
		&& TMath::Abs(cms2.mus_p4()[i].eta()) < 2.4)
	  ret |= CUT_BIT(LEP_GOOD_NOD0);
	
	if (fabs(cms2.mus_d0corr().at(i)) <= 0.02) 
	  ret |= CUT_BIT(LEP_D0);
	
  }

  //jet vars--single lepton: can't use this for njet in dileptons b'c doesn't clean for two leptons
  //these are the globals, in eventhistos, copy globals on selection
  //change this so that when the lepton is selected, copy these into a
  //make new ones njets_lep_20 for the lepton, then copy that into njets_20 for leps which pass x cuts.
  njets_20 = 0;
  njets_30 = 0;
  //this code ripped from selections.cc->getCaloJets
  //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets;
  //calo_jets.clear();
  
  for( unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    //if( dRbetweenVectors(lep_p4, cms2.jets_p4()[jj]) < 0.4 )
	if(  ROOT::Math::VectorUtil::DeltaR(lep_p4, cms2.jets_p4()[jj]) < 0.4 )
		//(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
	  continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4)
	  continue;
    //if (cms2.jets_emFrac()[jj] < 0.1) continue;
	//count
    if (cms2.jets_p4()[jj].pt() > 20)
	  njets_20++;
    if (cms2.jets_p4()[jj].pt() > 30)
	  njets_30++;
    //calo_jets.push_back(cms2.jets_p4()[jj]);
  }

  if( njets_20 == 0 )
	ret |= CUT_BIT(JET_VETO_20); //if njets_20 > 0, this cut fails
  if( njets_30 == 0 )
	ret |= CUT_BIT(JET_VETO_30); //if njets_30 > 0, this cut fails

  //transverse mass
  double dphi = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector( cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()),
																 cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()),
																 0, cms2.evt_tcmet())
												  , lep_p4 );
  double masst = sqrt( 2 * cms2.evt_tcmet() * lep_p4.Et() * ( 1 - cos(dphi) ) );
  
  //check yanjun's code:
  TVector3 tcMET;
  tcMET.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  double massyj = sqrt( ( tcMET.Pt() + lep_p4.Et())*( tcMET.Pt() + lep_p4.Et())
					   - ( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )*( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )
					   - ( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() )*( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() ) );

  if( fabs( tcMET.Phi() - cms2.evt_tcmetPhi() ) > 0.01 )
	cout << "Phi error " << tcMET.Phi() << "   " << cms2.evt_tcmetPhi() << endl;

  if( fabs( massyj - masst ) > 0.1 && masst > 2 && massyj > 2 )
	cout << "Mass disagreement " << masst << "   " << massyj << endl;

  //set looper member
  transmass = masst;

  if( masst > 40 && masst < 100 )
	ret |= CUT_BIT(TMASS_100);


  return ret;
}

void RLooper::FillEventHistos ()
{

  // get the event weight
  if( sample_.kFactor != 1 ) cout << "kFactor non-unity " << sample_.kFactor << endl;
  double weight = RLooper::Weight();

  //reset indicies!!!!!!!!
  elidxs[0] = elidxs[1] = -1;
  muidxs[0] = muidxs[1] = -1;
  
  // need to determine if this is a di-lepton
  // or a single lepton event
  int nels = 0, nmus = 0;
  cuts_t pass_lep_cut = 0; //cuts of lepton which passes lepcuts--only good for single lepton
  cuts_t pass_nlep_cuts = 0; //cuts of all leptons which pass lepcuts: & of individual

  //only cuts we care about are end cuts, and lep select cuts
  cuts_t w_cuts = wcuts_;
  cuts_t z_cuts = zcuts_;
  //plan for baselepcuts: have a cut bit in Cuts.h which is all the W/Z specific cuts. Make baselepcuts be wcuts_ but with all those turned off.
  //cuts_t baselepcuts = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) ); //for base lep counting--OLD--SEE Cuts.h/cc
  cuts_t lepcuts = (z_cuts & ~non_lepcuts); // z_cuts MUST include loose LEP_PT_X cut
 
  //check tcMET here b'c no point doing for every lepton--don't use for Z
  cuts_t tcmetcut = 0;
  if( cms2.evt_tcmet() > 20 )
	tcmetcut |= CUT_BIT(TCMET_20);
  if( cms2.evt_tcmet() > 25 )
	tcmetcut |= CUT_BIT(TCMET_25);
  if( cms2.evt_tcmet() > 30 )
	tcmetcut |= CUT_BIT(TCMET_30);
  if( cms2.evt_tcmet() > 35 )
	tcmetcut |= CUT_BIT(TCMET_35);

  //select els
  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	transmass = 0;
	cuts_t elcut = LepSelect(0, i); //0 for els
	elcut |= tcmetcut;
	if( (elcut & lepcuts) == lepcuts ) {
	  nels++;
	  pass_lep_cut = elcut;
	  //njetels_20 = njets_20; //njets_20 is set in LepSelect
	  //njetels_30 = njets_30;
	  if( elidxs[0] == -1 ) {
		elidxs[0] = i;
		pass_nlep_cuts = pass_lep_cut;
	  }
	  else if( elidxs[1] == -1 ) {
		elidxs[1] = i;
		pass_nlep_cuts &= pass_lep_cut;
	  }
	  else { //if > 2 els, check that they're not higher pt than selected
		if( cms2.els_p4()[i].pt() > cms2.els_p4()[elidxs[0]].pt()
			|| cms2.els_p4()[i].pt() > cms2.els_p4()[elidxs[1]].pt() )
		  cout << "extra el which is > pt than selected" << endl;
	  }
	}
  }

  //select mus
  for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
	transmass = 0;
	cuts_t mucut = LepSelect(1, i); //1 for mus
	mucut |= tcmetcut;
	if( (mucut & lepcuts) == lepcuts ) { //all cuts
	  nmus++;
	  pass_lep_cut = mucut;
	  //njetmus_20 = njets_20;
	  //njetmus_30 = njets_30;
	  if( muidxs[0] == -1 ) {
		muidxs[0] = i;
		pass_nlep_cuts = pass_lep_cut;
	  }
	  else if( muidxs[1] == -1 ) {
		muidxs[1] = i;
		pass_nlep_cuts &= pass_lep_cut;
	  }
	  else { //if > 2 mus, check that they're not higher pt than selected
		if( cms2.mus_p4()[i].pt() > cms2.mus_p4()[muidxs[0]].pt()
			|| cms2.mus_p4()[i].pt() > cms2.mus_p4()[muidxs[1]].pt() )
		  cout << "extra mu which is > pt than selected" << endl;
	  }
	}
  }

  //Z--enforce exactly two leptons and SF requirement
  //remember, nels, nmus require NO MET, NOR MASS, NOR NJET
  if( (nels == 2 && nmus == 0) || (nmus == 2 && nels == 0) ) {

	const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);

	//const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i]);
	//DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE

	int idx1 = (nels == 2 ? elidxs[0] : muidxs[0]);
	int idx2 = (nels == 2 ? elidxs[1] : muidxs[1]);
	pass_nlep_cuts |= DilepSelect(myType, idx1, idx2); 

	if( (pass_nlep_cuts & z_cuts) == z_cuts ) {
	  double pt1 = (nels == 2 ? cms2.els_p4()[idx1].pt() : cms2.mus_p4()[idx1].pt());
	  double pt2 = (nels == 2 ? cms2.els_p4()[idx2].pt() : cms2.mus_p4()[idx2].pt());
	  hdilep_0_pt[myType]->Fill( pt1, weight );
	  hdilep_0_pt[DILEPTON_ALL]->Fill( pt1, weight );
	  hdilep_1_pt[myType]->Fill( pt2, weight );
	  hdilep_1_pt[DILEPTON_ALL]->Fill( pt2, weight );
	  ZEvent();
	}
  }
  else if( ((nels == 1 && nmus == 0) || (nmus == 1 && nels == 0)) 
		   && (pass_lep_cut & w_cuts) == w_cuts ) {
	unsigned int lep_type = (nels == 1 ? 0 : 1);
	unsigned int idx = (nels == 1 ? elidxs[0] : muidxs[0]);
	double pt = (nels == 1 ? cms2.els_p4()[idx].pt() : cms2.mus_p4()[idx].pt());

	hlep_pt[lep_type]->Fill( pt, weight );
	hlep_pt[2]->Fill( pt, weight );

	//muon checks -- only for events which pass full selection
	if( nmus == 1 ) { // && (pass_lep_cut & w_cuts) == w_cuts ) {
	  mu_tot += weight;
	  if( cms2.mus_p4().size() > 1 ) {
		mu_multi += weight; //more than 1 mu, but 1 selected
		for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
		  if( (int)i == muidxs[0] ) continue;
		  double mass = (cms2.mus_p4()[i] + cms2.mus_p4()[muidxs[0]]).mass();
		  if( mass > 71 && mass < 111 ) //use wide window
			mu_mulma += weight; //more than 1, and makes z mass with selected
		}
	  }
	}

	//for WEvent, need to pass tmass, met
	//if( (pass_lep_cut & w_cuts) == w_cuts )
	WEvent();
  }

}
//end FillEventHistos

void RLooper::WEvent() {
  
  double weight = RLooper::Weight();

  unsigned int lep_type = 0; //default el
  if( elidxs[0] == -1 )	lep_type = 1;

  scands_passing_[lep_type] += weight;
  scands_passing_w2_[lep_type] += weight * weight;
  scands_count_[lep_type]++;
  scands_passing_[2] += weight;
  scands_passing_w2_[2] += weight * weight;
  scands_count_[2]++;

}

void RLooper::ZEvent ()
{
  double weight = RLooper::Weight();
  const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);

  dcands_passing_[myType] += weight;
  dcands_passing_w2_[myType] += weight * weight;
  dcands_count_[myType]++;
  dcands_passing_[DILEPTON_ALL] += weight;
  dcands_passing_w2_[DILEPTON_ALL] += weight * weight;
  dcands_count_[DILEPTON_ALL]++;

}


void RLooper::End ()
{
  cout << "mu events " << mu_tot
	   << "  mu events w > 1 mu " <<  mu_multi //more than 1 mu, but 1 selected
	   << "  which make z mass " << mu_mulma << endl; //more than 1, and makes z mass with selected

  /*
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
  */
}

