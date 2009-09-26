#include <math.h>
#include "TVector3.h"
//#include "Math/VectorUtil.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(dcands_passing_	    , 0, sizeof(dcands_passing_       ));
     memset(dcands_passing_w2_	, 0, sizeof(dcands_passing_w2_    ));
     memset(dcands_count_		, 0, sizeof(dcands_count_         ));
     memset(scands_passing_	    , 0, sizeof(scands_passing_       ));
     memset(scands_passing_w2_	, 0, sizeof(scands_passing_w2_    ));
     memset(scands_count_		, 0, sizeof(scands_count_         ));

	 weight = 0;
	 //initialize indicies
	 elidxs[0] = -1;
	 elidxs[1] = -1;
	 muidxs[0] = -1;
	 muidxs[1] = -1;

}

void Looper::FormatHist(TH1* hist)
{
	hist->SetFillColor(sample_.histo_color);
}

void Looper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
}


//to add: mass, transverse mass, njets
void Looper::BookHistos ()
{

  // single lepton histograms (two + 1 types)
  for (unsigned int i = 0; i < 3; ++i) {
	std::string hyp = "e";
	if (i == 1) hyp = "m";
	else if (i == 2) hyp = "all";

	h1_lep_Highpt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_pt", hyp.c_str()), 
								  "Highlep_pt", 100, 0.0, 100.0);
	FormatHist(h1_lep_Highpt_[i]);
        
	h1_lep_HighptMet_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_Met", hyp.c_str()), 
									 "Highlep_Met", 100, 0.0, 100.0);
	FormatHist(h1_lep_HighptMet_[i]);
        
	h1_lep_HighptRelIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_RelIso", hyp.c_str()), 
										"Highlep_RelIso", 120, -0.1, 1.1);
	FormatHist(h1_lep_HighptRelIso_[i]);
        
	h1_lep_HighptRelIsoPtLg20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_RelIsoPtLg20", hyp.c_str()), 
											  "Highlep_RelIsoPtLg20", 120, -0.1, 1.1);
	FormatHist(h1_lep_HighptRelIsoPtLg20_[i]);

	////
	h1_lep_Lowpt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_pt", hyp.c_str()), 
								 "Lowlep_pt", 100, 0.0, 100.0);
	FormatHist(h1_lep_Lowpt_[i]);
        
	h1_lep_LowptMet_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_Met", hyp.c_str()), 
									"Lowlep_Met", 100, 0.0, 100.0);
	FormatHist(h1_lep_LowptMet_[i]);
        
	h1_lep_LowptRelIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_RelIso", hyp.c_str()), 
									   "Lowlep_RelIso", 120, -0.1, 1.1);
	FormatHist(h1_lep_LowptRelIso_[i]);
        
	h1_lep_LowptRelIsoPtLg20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_RelIsoPtLg20", hyp.c_str()), 
											 "Lowlep_RelIsoPtLg20", 120, -0.1, 1.1);
	FormatHist(h1_lep_LowptRelIsoPtLg20_[i]);

	h1_lep_LowptNLepGt10Lt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_NLepGt10Lt20", hyp.c_str()), 
											 "Lowlep_NLepGt10Lt20", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt10Lt20_[i]);

	h1_lep_LowptNLepGt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Lowlep_NLepGt20", hyp.c_str()), 
										 "Lowlep_NLepGt20", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20_[i]);

	NewHist( hlep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
	NewHist( hlep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "lep_transmass", hyp.c_str()), "lep_transmass", 200, 0.0, 200.0);
	NewHist( hlep_met[i], Form("%s_%s_%s", SampleName().c_str(), "lep_met", hyp.c_str()), "lep_met", 100, 0.0, 100.0);
	NewHist( hlep_met_dphi[i], Form("%s_%s_%s", SampleName().c_str(), "lep_met_dphi", hyp.c_str()), "lep_met_dphi", 100, 0, 2 * 3.14159);
	NewHist( hlep_trckIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_trckIso", hyp.c_str()), "lep_trckIso", 100, 0.0, 100.0);
	NewHist( hlep_ecalIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_ecalIso", hyp.c_str()), "lep_ecalIso", 100, 0.0, 100.0);

	//for nlep, fill before cutting on it--same for W,Z
	NewHist( hlep_nlep[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep", hyp.c_str()), "lep_nlep", 10, -0.5, 9.5);
	NewHist( hlep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet20", hyp.c_str()), "lep_njet20", 10, -0.5, 9.5);
	NewHist( hlep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet30", hyp.c_str()), "lep_njet30", 10, -0.5, 9.5);
	NewHist( hlep_conv[i], Form("%s_%s_%s", SampleName().c_str(), "lep_conversions", hyp.c_str()), "lep_conversions", 2, -0.5, 1.5);
  }

  // di-lepton histograms (three + 1 types)
  for (unsigned int i = 0; i < 4; ++i) {
	NewHist( hdilep_0_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]), "dilep_0_pt", 100, 0.0, 100.0);
	NewHist( hdilep_1_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]), "dilep_1_pt", 100, 0.0, 100.0);
	NewHist( hdilep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_pt", dilepton_hypo_names[i]), "dilep_pt", 100, 0.0, 100.0);
	NewHist( hdilep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_mass", dilepton_hypo_names[i]), "dilep_mass", 200, 0.0, 200.0);
	NewHist( hdilep_met[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_met", dilepton_hypo_names[i]), "dilep_met", 100, 0.0, 100.0);
	NewHist( hdilep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet20", dilepton_hypo_names[i]), "dilep_njet20", 10, -0.5, 9.5);
	NewHist( hdilep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet30", dilepton_hypo_names[i]), "dilep_njet30", 10, -0.5, 9.5);

  }
	
  // event level histograms
  //NewHist( hdilep_nhyp, Form("%s_%s_%s", SampleName().c_str(), "dilep_nhyp", "all"), "dilep_nhyp", 10, -0.5, 9.5);

}

cuts_t Looper::DilepSelect() //(int i_hyp), no hyp, just idxs
{
  cuts_t ret = 0;
  float ptcut = 20.0;
  //int idx1 = (elidxs[0] != -1 ? elidxs[0] : muidxs[0]);
  //int idx2 = (elidxs[1] != -1 ? elidxs[1] : muidxs[1]);
  //int idx1, idx2;

  if( elidxs[0] != -1 && elidxs[1] != -1 ) {
	 
	//if (cms2.hyp_lt_p4()[i_hyp].pt() > ptcut && cms2.hyp_ll_p4()[i_hyp].pt() > ptcut)
	if( cms2.els_p4()[elidxs[0]].pt() >= ptcut && cms2.els_p4()[elidxs[1]].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	//if( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 )
	if( cms2.els_charge()[elidxs[0]] * cms2.els_charge()[elidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  else if( muidxs[0] != -1 && muidxs[1] != -1 ) {

	if( cms2.mus_p4()[muidxs[0]].pt() >= ptcut && cms2.mus_p4()[muidxs[1]].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	if( cms2.mus_charge()[muidxs[0]] * cms2.mus_charge()[muidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  else
	cout << "BAD DILEPSELECT CALL" << endl;

  return ret;
}

cuts_t Looper::LepSelect(int lep_type, int i)
{
  cuts_t ret = 0;

  float ptcut = 20.0;

  // e
  if (lep_type == 0) {

	if (cms2.els_p4()[i].pt() >= ptcut)
	  ret |= CUT_BIT(LEP_PT);

	if( GoodSusyElectronWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyElectronIsolation(i, true) ) //bool is for use calo
	  ret |= CUT_BIT(LEP_ISO);
	
  }

  // m
  else if (lep_type == 1) {

	if (cms2.mus_p4()[i].pt() >= ptcut)
	  ret |= CUT_BIT(LEP_PT);

	if( GoodSusyMuonWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyMuonIsolation(i) )
	  ret |= CUT_BIT(LEP_ISO);
	
  }

  return ret;
}

void Looper::FillEventHistos ()
{

  // a first loop over the most energetic 
  // muon and electron to determine a good cut 
  // in pT and iso
  int   hiPtIdx    = -1;
  float hiPtmax    = -1.;
  
  // get the event weight
  if( sample_.kFactor != 1 ) cout << "kFactor non-unity " << sample_.kFactor << endl;
  weight = cms2.evt_scale1fb() * sample_.kFactor / 1000; //1pb
  
  // have a look at the highest pt electron
  for(int ele = 0; ele < (int)cms2.els_p4().size(); ++ele) {
    if(cms2.els_p4()[ele].pt() > hiPtmax && GoodSusyElectronWithoutIsolation(ele) ) hiPtIdx = ele;
  }

  // histogram indices are e, m, all (0, 1, 2)
  if(hiPtIdx != -1)  {
    h1_lep_Highpt_[0]                                                    ->Fill(cms2.els_p4()[hiPtIdx].pt(), weight);
    h1_lep_HighptMet_[0]                                                 ->Fill(cms2.evt_tcmet(), weight);
    h1_lep_HighptRelIso_[0]                                              ->Fill(inv_el_relsusy_iso(hiPtIdx, true), weight);
    if(cms2.els_p4()[hiPtIdx].pt() > 20. ) h1_lep_HighptRelIsoPtLg20_[0] ->Fill(inv_el_relsusy_iso(hiPtIdx, true), weight);
  }

 // have a look at all but the highest pt electron
  uint nEleGt10Lt20 = 0;
  uint nEleGt20     = 0;
  for(int ele = 0; ele <  (int)cms2.els_p4().size(); ++ele) {
    if(hiPtIdx != -1 && hiPtIdx != ele) {
      h1_lep_Lowpt_[0]                                                ->Fill(cms2.els_p4()[ele].pt(), weight);
      h1_lep_LowptMet_[0]                                             ->Fill(cms2.evt_tcmet(), weight);
      h1_lep_LowptRelIso_[0]                                          ->Fill(inv_el_relsusy_iso(ele, true), weight);
      if(cms2.els_p4()[ele].pt() > 20. ) h1_lep_LowptRelIsoPtLg20_[0] ->Fill(inv_el_relsusy_iso(hiPtIdx, true), weight);
      if(cms2.els_p4()[ele].pt() > 10. && (cms2.els_p4()[ele].pt() < 20. )) ++nEleGt10Lt20;
      if(cms2.els_p4()[ele].pt() > 20. )                                    ++nEleGt20;
    }
  }
  h1_lep_LowptNLepGt10Lt20_[0]->Fill(nEleGt10Lt20, weight);
  h1_lep_LowptNLepGt20_[0]    ->Fill(nEleGt20, weight);

  // have a look at the highest pt muon
  // reset highpt index
  hiPtIdx    = -1;
  hiPtmax    = -1.;
  for(int muo = 0; muo <  (int)cms2.mus_p4().size(); ++muo) {
    if(cms2.mus_p4()[muo].pt() > hiPtmax && GoodSusyMuonWithoutIsolation(muo) ) hiPtIdx = muo;
  }
  
  // histogram indices are e, m, all (0, 1, 2)
  if(hiPtIdx != -1)  {
    h1_lep_Highpt_[1]                                                    ->Fill(cms2.mus_p4()[hiPtIdx].pt(), weight);
    h1_lep_HighptMet_[1]                                                 ->Fill(cms2.evt_tcmet(), weight);
    h1_lep_HighptRelIso_[1]                                              ->Fill(inv_mu_relsusy_iso(hiPtIdx), weight);
    if(cms2.mus_p4()[hiPtIdx].pt() > 20. ) h1_lep_HighptRelIsoPtLg20_[1] ->Fill(inv_mu_relsusy_iso(hiPtIdx), weight);
  }

  // have a look at all but the highest pt muoctron
  uint nMuoGt10Lt20 = 0;
  uint nMuoGt20     = 0;
  for(int muo = 0; muo <  (int)cms2.mus_p4().size(); ++muo) {
    if(hiPtIdx != -1 && hiPtIdx != muo) {
      h1_lep_Lowpt_[1]                                                ->Fill(cms2.mus_p4()[muo].pt(), weight);
      h1_lep_LowptMet_[1]                                             ->Fill(cms2.evt_tcmet(), weight);
      h1_lep_LowptRelIso_[1]                                          ->Fill(inv_mu_relsusy_iso(muo), weight);
      if(cms2.mus_p4()[muo].pt() > 20. ) h1_lep_LowptRelIsoPtLg20_[1] ->Fill(inv_mu_relsusy_iso(muo), weight);
      if(cms2.mus_p4()[muo].pt() > 10. && (cms2.mus_p4()[muo].pt() < 20. )) ++nMuoGt10Lt20;
      if(cms2.mus_p4()[muo].pt() > 20. )                                    ++nMuoGt20;
    }
  }
  h1_lep_LowptNLepGt10Lt20_[1]->Fill(nMuoGt10Lt20, weight);
  h1_lep_LowptNLepGt20_[1]    ->Fill(nMuoGt20, weight);

  
  // need to determine if this is a di-lepton
  // or a single lepton event
  int nels = 0, nmus = 0;
  int nels_nopt = 0, nmus_nopt = 0;
  elidxs[0] = elidxs[1] = -1;
  muidxs[0] = muidxs[1] = -1;
  int elidxs_nopt[] = {-1, -1}; //this one can be local--still just 2
  int muidxs_nopt[] = {-1, -1};
  cuts_t lepcuts = (CUT_BIT(LEP_PT)
					| CUT_BIT(LEP_GOOD)
					| CUT_BIT(LEP_ISO)
					);
  cuts_t lepcuts_nopt = (CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO));
  cuts_t lepcuts_noiso = (CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD));

  //select els
  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	cuts_t elcut = LepSelect(0, i); //0 for els
	if( (elcut & lepcuts) == lepcuts ) { //all cuts
	  nels++;
	  if( elidxs[0] == -1 )
		elidxs[0] = i;
	  else if( elidxs[1] == -1 )
		elidxs[1] = i;
	  //if > 2 els, ignore the rest--should be sorted by pt already
	}
	if( (elcut & lepcuts_nopt) == lepcuts_nopt ) { //no pt cut
	  nels_nopt++;
	  if( elidxs_nopt[0] == -1 )
		elidxs_nopt[0] = i;
	  else if( elidxs_nopt[1] == -1 )
		elidxs_nopt[1] = i;
	}
	if( (elcut & lepcuts_noiso) == lepcuts_noiso ) {
	  hlep_ecalIso[0]->Fill( cms2.els_ecalIso()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.els_ecalIso()[i], weight );
	}
  }

  //select mus
  for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
	cuts_t mucut = LepSelect(1, i); //1 for mus
	if( (mucut & lepcuts) == lepcuts ) { //all cuts
	  nmus++;
	  if( muidxs[0] == -1 )
		muidxs[0] = i;
	  else if( muidxs[1] == -1 )
		muidxs[1] = i;
	}
	if( (mucut & lepcuts_nopt) == lepcuts_nopt ) { //no pt cut
	  nmus_nopt++;
	  if( muidxs_nopt[0] == -1 )
		muidxs_nopt[0] = i;
	  else if( muidxs_nopt[1] == -1 )
		muidxs_nopt[1] = i;
	}
	if( (mucut & lepcuts_noiso) == lepcuts_noiso ) {
	  hlep_ecalIso[1]->Fill( cms2.mus_iso03_emEt()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.mus_iso03_emEt()[i], weight );
	}
  }

  //fill nlep cut before checking nels, nmus
  //0=el, 1=mu, 2=all
  //if( nels > 0 )
  hlep_nlep[0]->Fill( nels, weight );
  //if( nmus > 0 )
  hlep_nlep[1]->Fill( nmus, weight );
  //if( nels > 0 || nmus > 0 )
  hlep_nlep[2]->Fill( nels+nmus, weight );

  //fill lep pt hists based on n nopt--single leps (each flavor independent of other)
  if( nels_nopt == 1 ) { //allow mus
	hlep_pt[0]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hlep_pt[2]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
  }
  if( nmus_nopt == 1 ) {
	hlep_pt[1]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hlep_pt[2]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
  }
  //now for N leps  //DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE
  if( nels_nopt > 1 ) { //allow mus
	hdilep_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_0_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_0_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_1_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_1_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
  }
  if( nmus_nopt > 1 ) {
	hdilep_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_0_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_0_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_1_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_1_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
  }
  
  //enforce exactly two leptons and SF requirement
  //if( cms2.evt_nels() == 2 || cms2.mus_p4().size() == 2 )
  //if( cms2.evt_nels() > 1 || cms2.mus_p4().size() > 1 )
  if( (nels == 2 && nmus == 0)
	  || (nmus == 2 && nels == 0) )
	ZEvent();
  //else if( (cms2.evt_nels() == 0 && cms2.mus_p4().size() == 1) ||
  //	   (cms2.evt_nels() == 1 && cms2.mus_p4().size() == 0))
  else if( (nels == 1 && nmus == 0)
		   || (nmus == 1 && nels == 0) )
	WEvent();

}
//end FillEventHistos

void Looper::WEvent() {
  
  // histogram indices are e, m, all (0, 1, 2)
  // get lep_type, lep_p4
  unsigned int lep_type = 0; //default el
  LorentzVector lep_p4;
  if( elidxs[0] == -1 ) { //no els
	lep_type = 1;
	lep_p4 = cms2.mus_p4()[muidxs[0]];
  }
  else 
	lep_p4 = cms2.els_p4()[elidxs[0]];
	
  hlep_met[lep_type]->Fill(cms2.evt_tcmet(), weight);
  hlep_met[2]->Fill(cms2.evt_tcmet(), weight);

  //put the event-ish cuts here
  if( cms2.evt_tcmet() <= 20 )
	return;

  //check yanjun's code:
  TVector3 tcMET;
  tcMET.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  double masst = sqrt( ( tcMET.Pt() + lep_p4.Et())*( tcMET.Pt() + lep_p4.Et())
						//double massyj = sqrt( ( tcMET.Pt() + lep_p4.E() )*( tcMET.Pt() + lep_p4.E() )
					   - ( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )*( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )
					   - ( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() )*( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() ) );

  //this is xyze
  //lep_p4.SetPz(0); //set z comp to zero so can use .mass because .mt isn't what we want
  //double masst = ( LorentzVector(cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()), cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()), 0, cms2.evt_tcmet())
  //			   + lep_p4 ).mass();
  //+ (lep_type == 0 ? cms2.els_p4()[elidxs[0]] : cms2.mus_p4()[muidxs[0]]) ).mt();

  hlep_mass[lep_type]->Fill( masst, weight ); //check all cuts but mass for mass plot (n-1)
  hlep_mass[2]->Fill( masst, weight );

  //masst = massyj; //Now, I get agreement w/ YJ if she uses E, but NOT if she uses Et, which she does. Her mt's are less, so fewer fail cut.

  if( fabs( tcMET.Phi() - cms2.evt_tcmetPhi() ) > 0.01 )
	cout << "Phi error " << tcMET.Phi() << "   " << cms2.evt_tcmetPhi() << endl;

  //if( fabs( massyj - masst ) > 0.1 )
  //cout << "Mass disagreement " << masst << "   " << massyj << endl;

  if( masst <= 40 || masst >= 100 )
	return;

  //conversion plot has ALL cuts applied
  if( lep_type == 0 ) {
	hlep_conv[lep_type]->Fill( conversionElectron(elidxs[0]), weight );
	hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );
  }
  
  //if (cms2.mus_p4().size() == 0) lep_type = 0;

  //for W, all checking is done
  // define the cuts to be used
  //cuts_t cuts = (CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) );
  // find out what cuts passed
  //cuts_t cuts_passed = LepSelect(lep_type, 0);

  //if ((cuts_passed & cuts) == cuts) {

  //jet vars--single lepton
  int njets_20 = 0;
  int njets_30 = 0;
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
  
  //if (calo_jets.size() > 1) {
  //sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  //}
  //return calo_jets;

  hlep_njet20[lep_type]->Fill(njets_20, weight);
  hlep_njet20[2]->Fill(njets_20, weight);
  hlep_njet30[lep_type]->Fill(njets_30, weight);
  hlep_njet30[2]->Fill(njets_30, weight);

  if( lep_type == 0 ) {
	//fill lep pt before/during lep selection
	//hlep_pt[lep_type]->Fill(cms2.els_p4()[elidxs[0]].pt(), weight);
	//hlep_pt[2]->Fill(cms2.els_p4()[elidxs[0]].pt(), weight);

	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[elidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);

	hlep_trckIso[lep_type]->Fill(cms2.els_tkIso()[elidxs[0]], weight);
	hlep_trckIso[2]->Fill(cms2.els_tkIso()[elidxs[0]], weight);
  }
  else if (lep_type == 1) {
	//hlep_pt[lep_type]->Fill(cms2.mus_p4()[muidxs[0]].pt(), weight);
	//hlep_pt[2]->Fill(cms2.mus_p4()[muidxs[0]].pt(), weight);
	
	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);
	
	hlep_trckIso[lep_type]->Fill(cms2.mus_iso03_sumPt()[muidxs[0]], weight);
	hlep_trckIso[2]->Fill(cms2.mus_iso03_sumPt()[muidxs[0]], weight);

  }

  scands_passing_[lep_type] += weight;
  scands_passing_w2_[lep_type] += weight * weight;
  scands_count_[lep_type]++;
  scands_passing_[2] += weight;
  scands_passing_w2_[2] += weight * weight;
  scands_count_[2]++;

  //} //end if passed cuts

}

void Looper::ZEvent ()
{
  //hdilep_nhyp_->Fill(cms2.hyp_p4().size(), weight);

  // define the cuts to be used
  cuts_t cuts = (CUT_BIT(DILEP_PT)
				 | CUT_BIT(DILEP_OS)
				 | CUT_BIT(DILEP_MASS)
				 //| CUT_BIT(LEP_GOOD) //these already checked
				 //| CUT_BIT(LEP_ISO)
				 );
  cuts_t cuts_nomass = (CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_OS));

  //get leptons from indicies, not hyps
  //for( unsigned int i = 0; i < cms2.hyp_p4().size(); ++i) {
  // get what type of di-lepton hypothesis this is
  //const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i]);
  //DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE

  const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);
  //int myType = (elidxs[0] != -1 ? 3 : 1);

  // does this hypothesis pass the required cuts?
  cuts_t cuts_passed = DilepSelect();

  //these already checked
  //require both hyps to pass lep select (1st arg:0=e,1=m)
  //cuts_passed |= LepSelect( abs(cms2.hyp_lt_id()[i]) == 11 ? 0 : 1, cms2.hyp_lt_index()[i] );
  //cuts_passed |= LepSelect( abs(cms2.hyp_ll_id()[i]) == 11 ? 0 : 1, cms2.hyp_ll_index()[i] );

  //lep vars
  LorentzVector lep1_p4, lep2_p4;	
  double pt1=0, pt2=0, mass=0;
  if( elidxs[0] != -1 && elidxs[1] != -1 ) {
	lep1_p4 = cms2.els_p4()[elidxs[0]];
	lep2_p4 = cms2.els_p4()[elidxs[1]];
	pt1 = cms2.els_p4()[elidxs[0]].pt();
	pt2 = cms2.els_p4()[elidxs[1]].pt();
	mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
  }
  else if( muidxs[0] != -1 && muidxs[1] != -1 ) {
	lep1_p4 = cms2.mus_p4()[muidxs[0]];
	lep2_p4 = cms2.mus_p4()[muidxs[1]];
	pt1 = cms2.mus_p4()[muidxs[0]].pt();
	pt2 = cms2.mus_p4()[muidxs[1]].pt();
	mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
  }
  else {
	cout << "BAD LEP INDXS IN ZEVENT\n" << endl;
	exit(1);
  }

  //fill pt hists before/during lep selection
  // fill histograms for the 0th lepton
  //hdilep_0_pt[myType]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
  //hdilep_0_pt[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
  //hdilep_0_pt[myType]->Fill(pt1, weight);
  //hdilep_0_pt[DILEPTON_ALL]->Fill(pt1, weight);
	
  // fill histograms for the 1th lepton
  //hdilep_1_pt[myType]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
  //hdilep_1_pt[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
  //hdilep_1_pt[myType]->Fill(pt2, weight);
  //hdilep_1_pt[DILEPTON_ALL]->Fill(pt2, weight);
	
  // fill hypothesis level histograms
  //hdilep_mass[myType]->Fill(cms2.hyp_p4()[i].mass(), weight);
  //hdilep_mass[DILEPTON_ALL]->Fill(cms2.hyp_p4()[i].mass(), weight);
  if( (cuts_passed & cuts_nomass) == cuts_nomass ) {
	hdilep_mass[myType]->Fill(mass, weight);
	hdilep_mass[DILEPTON_ALL]->Fill(mass, weight);
  }

  if( (cuts_passed & cuts) == cuts ) {

	//jet vars--dilep case
	int njets_20 = 0;
	int njets_30 = 0;
	//this code ripped from selections.cc->getCaloJets
	//vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets;
	//calo_jets.clear();
  
	for( unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
	  //if( dRbetweenVectors(lep1_p4, cms2.jets_p4()[jj]) < 0.4 ||
	  //  dRbetweenVectors(lep2_p4, cms2.jets_p4()[jj]) < 0.4 )
	  if(  ROOT::Math::VectorUtil::DeltaR(lep1_p4, cms2.jets_p4()[jj]) < 0.4 ||
		   ROOT::Math::VectorUtil::DeltaR(lep2_p4, cms2.jets_p4()[jj]) < 0.4 )
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
  
	//if (calo_jets.size() > 1) {
	//sort(calo_jets.begin(), calo_jets.end(),  comparePt);
	//}
	//return calo_jets;

	if( myType == DILEPTON_EE ) {
	  hlep_conv[0]->Fill( conversionElectron(elidxs[0]), weight );
	  hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );
	}
	
	hdilep_njet20[myType]->Fill(njets_20, weight);
	hdilep_njet20[DILEPTON_ALL]->Fill(njets_20, weight);
	hdilep_njet30[myType]->Fill(njets_30, weight);
	hdilep_njet30[DILEPTON_ALL]->Fill(njets_30, weight);

	//no met cut, so require all cuts for met plot
	hdilep_met[myType]->Fill(cms2.evt_tcmet(), weight);
	hdilep_met[DILEPTON_ALL]->Fill(cms2.evt_tcmet(), weight);

	//all other hists need to fill before checking cuts

	dcands_passing_[myType] += weight;
	dcands_passing_w2_[myType] += weight * weight;
	dcands_count_[myType]++;
	dcands_passing_[DILEPTON_ALL] += weight;
	dcands_passing_w2_[DILEPTON_ALL] += weight * weight;
	dcands_count_[DILEPTON_ALL]++;

  }
  //} //end loop on hyp

}


void Looper::End ()
{
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

