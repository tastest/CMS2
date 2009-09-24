# include <math.h>
#include "TVector3.h"
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
/*
void Looper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
}
*/

//to add: mass, transverse mass, njets
void Looper::BookHistos ()
{

  // single lepton histograms (two + 1 types)
  for (unsigned int i = 0; i < 3; ++i) {
	std::string hyp = "e";
	if (i == 1) hyp = "m";
	if (i == 2) hyp = "all";

	h1_lep_pt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), 
							  "lep_tkIso", 100, 0.0, 100.0);
	FormatHist(h1_lep_pt_[i]);

	h1_lep_met_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "lep_met", hyp.c_str()),
							   "lep_met", 100, 0.0, 100.0);
	FormatHist(h1_lep_met_[i]);

	h1_lep_met_dphi_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "lep_met_dphi", hyp.c_str()),
									"lep_met_dphi", 100, 0, 2 * 3.14159);
	FormatHist(h1_lep_met_dphi_[i]);

	h1_lep_tkIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "lep_tkIso", hyp.c_str()),
								 "lep_tkIso", 100, 0.0, 100.0);
	FormatHist(h1_lep_tkIso_[i]);

  }

  // di-lepton histograms (three + 1 types)
  for (unsigned int i = 0; i < 4; ++i) {
	h1_dilep_0_pt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]),
								  "dilep_0_pt", 100, 0.0, 100.0);
	FormatHist(h1_dilep_0_pt_[i]);

	h1_dilep_1_pt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]),
								  "dilep_1_pt", 100, 0.0, 100.0);
	FormatHist(h1_dilep_1_pt_[i]);

	h1_dilep_mass_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "dilep_mass", dilepton_hypo_names[i]),
								  "dilep_mass", 200, 0.0, 200.0);
	FormatHist(h1_dilep_mass_[i]);

	h1_dilep_met_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "dilep_met", dilepton_hypo_names[i]),
								 "lep_met", 100, 0.0, 100.0);
	FormatHist(h1_dilep_met_[i]);


  }
	
  // event level histograms
  h1_dilep_nhyp_ = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "dilep_nhyp", "all"),
							"dilep_nhyp", 10, -0.5, 9.5);
  FormatHist(h1_dilep_nhyp_);

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
	if( cms2.els_p4()[elidxs[0]].pt() > ptcut && cms2.els_p4()[elidxs[1]].pt() > ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	//if( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 )
	if( cms2.els_charge()[elidxs[0]] * cms2.els_charge()[elidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
	if( mass > 76 && mass < 106 )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  else if( muidxs[0] != -1 && muidxs[1] != -1 ) {

	if( cms2.mus_p4()[muidxs[0]].pt() > ptcut && cms2.mus_p4()[muidxs[1]].pt() > ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	if( cms2.mus_charge()[muidxs[0]] * cms2.mus_charge()[muidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
	if( mass > 76 && mass < 106 )
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

	if (cms2.els_p4()[i].pt() > ptcut)
	  ret |= CUT_BIT(LEP_PT);

	if( GoodSusyElectronWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyElectronIsolation(i, true) ) //bool is for use calo
	  ret |= CUT_BIT(LEP_ISO);
	
  }

  // m
  else if (lep_type == 1) {

	if (cms2.mus_p4()[i].pt() > ptcut)
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

  // need to determine if this is a di-lepton
  // or a single lepton event
  int nels = 0;
  int nmus = 0;
  elidxs[0] = elidxs[1] = -1;
  muidxs[0] = muidxs[1] = -1;
  cuts_t lepcuts = (CUT_BIT(LEP_PT)
					| CUT_BIT(LEP_GOOD)
					| CUT_BIT(LEP_ISO)
					);


  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	cuts_t elcut = LepSelect(0, i); //0 for els
	if( (elcut & lepcuts) == lepcuts ) {
	  nels++;
	  if( elidxs[0] == -1 )
		elidxs[0] = i;
	  else if( elidxs[1] == -1 )
		elidxs[1] = i;
	  //if > 2 els, ignore the rest--should be sorted by pt already
	}
  }

  for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
	cuts_t mucut = LepSelect(1, i); //1 for mus
	if( (mucut & lepcuts) == lepcuts ) {
	  nmus++;
	  if( muidxs[0] == -1 )
		muidxs[0] = i;
	  else if( muidxs[1] == -1 )
		muidxs[1] = i;
	}
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

void Looper::WEvent()
{
  //put the event-ish cuts here
  if( cms2.evt_tcmet() <= 20 )
	return;

  // histogram indices are e, m, all (0, 1, 2)
  unsigned int lep_type = 0; //default el
  if( elidxs[0] == -1 ) //no els
	lep_type = 1;

  //this is xyze
  double masst = ( LorentzVector(cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()), cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()), 0, cms2.evt_tcmet())
				   + (lep_type == 0 ? cms2.els_p4()[elidxs[0]] : cms2.mus_p4()[muidxs[0]]) ).mt();
  if( masst <= 40 || masst >= 100 )
	return;

  // get the event weight
  float weight = cms2.evt_scale1fb() * sample_.kFactor / 1000; //1pb

  //if (cms2.mus_p4().size() == 0) lep_type = 0;

  //for W, all checking is done
  // define the cuts to be used
  //cuts_t cuts = (CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) );
  // find out what cuts passed
  //cuts_t cuts_passed = LepSelect(lep_type, 0);

  //if ((cuts_passed & cuts) == cuts) {

  if( lep_type == 0 ) {	
	h1_lep_pt_[lep_type]->Fill(cms2.els_p4()[elidxs[0]].pt(), weight);
	h1_lep_pt_[2]->Fill(cms2.els_p4()[elidxs[0]].pt(), weight);
		
	h1_lep_met_[lep_type]->Fill(cms2.evt_tcmet(), weight);
	h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);

	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[elidxs[0]].Phi() ));
	h1_lep_met_dphi_[lep_type]->Fill(dphi, weight);
	h1_lep_met_dphi_[2]->Fill(dphi, weight);

	h1_lep_tkIso_[lep_type]->Fill(cms2.els_tkIso()[elidxs[0]], weight);
	h1_lep_tkIso_[2]->Fill(cms2.els_tkIso()[elidxs[0]], weight);
  }
  if (lep_type == 1) {
	h1_lep_pt_[lep_type]->Fill(cms2.mus_p4()[muidxs[0]].pt(), weight);
	h1_lep_pt_[2]->Fill(cms2.mus_p4()[muidxs[0]].pt(), weight);
                
	h1_lep_met_[lep_type]->Fill(cms2.evt_tcmet(), weight);
	h1_lep_met_[2]->Fill(cms2.evt_tcmet(), weight);
	
	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muidxs[0]].Phi() ));
	h1_lep_met_dphi_[lep_type]->Fill(dphi, weight);
	h1_lep_met_dphi_[2]->Fill(dphi, weight);
	
	h1_lep_tkIso_[lep_type]->Fill(cms2.mus_iso03_sumPt()[muidxs[0]], weight);
	h1_lep_tkIso_[2]->Fill(cms2.mus_iso03_sumPt()[muidxs[0]], weight);

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

  // get the event weight
  float weight = cms2.evt_scale1fb() * sample_.kFactor / 1000; //1pb
  h1_dilep_nhyp_->Fill(cms2.hyp_p4().size(), weight);

  // define the cuts to be used
  cuts_t cuts = (CUT_BIT(DILEP_PT)
				 | CUT_BIT(DILEP_OS)
				 | CUT_BIT(DILEP_MASS)
				 //| CUT_BIT(LEP_GOOD) //these already checked
				 //| CUT_BIT(LEP_ISO)
				 );

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

  if( (cuts_passed & cuts) == cuts ) {

	double pt1=0, pt2=0, mass=0;
	if( elidxs[0] != -1 ) {
	  pt1 = cms2.els_p4()[elidxs[0]].pt();
	  pt2 = cms2.els_p4()[elidxs[1]].pt();
	  mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
	}
	else {
	  pt1 = cms2.mus_p4()[muidxs[0]].pt();
	  pt2 = cms2.mus_p4()[muidxs[1]].pt();
	  mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
	}
	
	// fill histograms for the 0th lepton
	//h1_dilep_0_pt_[myType]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
	//h1_dilep_0_pt_[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
	h1_dilep_0_pt_[myType]->Fill(pt1, weight);
	h1_dilep_0_pt_[DILEPTON_ALL]->Fill(pt1, weight);
	
	// fill histograms for the 1th lepton
	//h1_dilep_1_pt_[myType]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
	//h1_dilep_1_pt_[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
	h1_dilep_1_pt_[myType]->Fill(pt2, weight);
	h1_dilep_1_pt_[DILEPTON_ALL]->Fill(pt2, weight);
	
	// fill hypothesis level histograms
	//h1_dilep_mass_[myType]->Fill(cms2.hyp_p4()[i].mass(), weight);
	//h1_dilep_mass_[DILEPTON_ALL]->Fill(cms2.hyp_p4()[i].mass(), weight);
	h1_dilep_mass_[myType]->Fill(mass, weight);
	h1_dilep_mass_[DILEPTON_ALL]->Fill(mass, weight);

	h1_dilep_met_[myType]->Fill(cms2.evt_tcmet(), weight);
	h1_dilep_met_[DILEPTON_ALL]->Fill(cms2.evt_tcmet(), weight);

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

