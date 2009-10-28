#include <math.h>
#include "TVector3.h"
//#include "Math/VectorUtil.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "ABCDLooper.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

ABCDLooper::ABCDLooper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
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

}


void ABCDLooper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
}


//to add: mass, transverse mass, njets
void ABCDLooper::BookHistos ()
{
  double metmax = 80.;
  int metbins = 80;

  // single lepton histograms (two + 1 types)
  for (unsigned int i = 0; i < 3; ++i) {
	std::string hyp = "e";
	if (i == 1) hyp = "m";
	else if (i == 2) hyp = "all";
	double d0max = 0.1;
        
	NewHist( hlep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
	NewHist( hlep_genpt[i], Form("%s_%s_%s", SampleName().c_str(), "lep_genpt", hyp.c_str()), "lep_genpt", 100, 0.0, 100.0);
	NewHist( hlep_genpt_mch[i], Form("%s_%s_%s", SampleName().c_str(), "lep_genpt_mch", hyp.c_str()), "lep_genpt_mch", 100, 0.0, 100.0);
	NewHist( hlep_pt_f[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt_f", hyp.c_str()), "lep_pt_f", 100, 0.0, 100.0);
	NewHist( hlep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "lep_transmass", hyp.c_str()), "lep_transmass", 200, 0.0, 200.0);
	NewHist( hlep_tcmet[i], Form("%s_%s_%s", SampleName().c_str(), "lep_tcmet", hyp.c_str()), "lep_met", 100, 0.0, 100.0);
	NewHist( hlep_clmumet[i], Form("%s_%s_%s", SampleName().c_str(), "lep_calomet_muon", hyp.c_str()), "lep_calomet_muon", 100, 0.0, 100.0);
	NewHist( hlep_met_dphi[i], Form("%s_%s_%s", SampleName().c_str(), "lep_met_dphi", hyp.c_str()), "lep_met_dphi", 100, 0, 2 * 3.14159);
	NewHist( hlep_trckIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_trckIso", hyp.c_str()), "lep_trckIso", 100, 0.0, 10.0);
	NewHist( hlep_ecalIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_ecalIso", hyp.c_str()), "lep_ecalIso", 100, 0.0, 10.0);
	NewHist( hlep_hcalIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_hcalIso", hyp.c_str()), "lep_hcalIso", 100, 0.0, 10.0);
	NewHist( hlep_relIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_relIso", hyp.c_str()), "lep_relIso", 100, 0.0, 1.0);
	NewHist( hlep_d0[i], Form("%s_%s_%s", SampleName().c_str(), "lep_d0", hyp.c_str()), "lep_d0", 100, 0.0, d0max);
	NewHist( hlep_d0Sig[i], Form("%s_%s_%s", SampleName().c_str(), "lep_d0Sig", hyp.c_str()), "lep_d0Sig", 100, 0.0, 100*d0max);

	//for nlep, fill before cutting on it--same for W,Z
	NewHist( hlep_nlep_nod0iso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep_nod0iso", hyp.c_str()), "lep_nlep_nod0iso", 10, -0.5, 9.5);
	NewHist( hlep_nlep_nometiso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep_nometiso", hyp.c_str()), "lep_nlep_nometiso", 10, -0.5, 9.5);
	NewHist( hlep_nlep[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep", hyp.c_str()), "lep_nlep", 10, -0.5, 9.5);
	NewHist( hlep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet20", hyp.c_str()), "lep_njet20", 10, -0.5, 9.5);
	NewHist( hlep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet30", hyp.c_str()), "lep_njet30", 10, -0.5, 9.5);
	NewHist( hlep_conv[i], Form("%s_%s_%s", SampleName().c_str(), "lep_conversions", hyp.c_str()), "lep_conversions", 2, -0.5, 1.5);

	//TH2F's for ABCD
	hlep_d0_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_trckIso", hyp.c_str()), "lep_d0_trckIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_ecalIso", hyp.c_str()), "lep_d0_ecalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_hcalIso", hyp.c_str()), "lep_d0_hcalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_relIso",  hyp.c_str()), "lep_d0_relIso", 100, 0.0, d0max, 100, 0.0, 1.0);
	hlep_met_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_trckIso", hyp.c_str()), "lep_tcMet_trckIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hlep_met_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_ecalIso", hyp.c_str()), "lep_tcMet_ecalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hlep_met_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_hcalIso", hyp.c_str()), "lep_tcMet_hcalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hlep_met_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_relIso",  hyp.c_str()), "lep_tcMet_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//only if neutrino is in eta 2.4
	hlep_accmet_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_acc_tcMet_relIso",  hyp.c_str()), "lep_acc_tcMet_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
  }

  // di-lepton histograms (three + 1 types)
  for (unsigned int i = 0; i < 4; ++i) {
	NewHist( hdilep_0_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]), "dilep_0_pt", 100, 0.0, 100.0);
	NewHist( hdilep_1_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]), "dilep_1_pt", 100, 0.0, 100.0);
	NewHist( hdilep_genpt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_genpt", dilepton_hypo_names[i]), "dilep_genpt", 100, 0.0, 100.0);
	//NewHist( hdilep_0_genpt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_genpt", dilepton_hypo_names[i]), "dilep_0_genpt", 100, 0.0, 100.0);
	//NewHist( hdilep_1_genpt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_genpt", dilepton_hypo_names[i]), "dilep_1_genpt", 100, 0.0, 100.0);
	NewHist( hdilep_0_genpt_mch[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_genpt_mch", dilepton_hypo_names[i]), "dilep_0_genpt_mch", 100, 0.0, 100.0);
	NewHist( hdilep_1_genpt_mch[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_genpt_mch", dilepton_hypo_names[i]), "dilep_1_genpt_mch", 100, 0.0, 100.0);
	NewHist( hdilep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_pt", dilepton_hypo_names[i]), "dilep_pt", 100, 0.0, 100.0);
	NewHist( hdilep_pt_mgen[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_pt_mgen", dilepton_hypo_names[i]), "dilep_pt_mgen", 200, -1.0, 1.0);
	NewHist( hdilep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_mass", dilepton_hypo_names[i]), "dilep_mass", 200, 0.0, 200.0);
	NewHist( hdilep_mass_ll20[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_mass_ll20", dilepton_hypo_names[i]), "dilep_mass_ll20", 200, 0.0, 200.0);
	NewHist( hdilep_tmass[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_tmass", dilepton_hypo_names[i]), "dilep_tmass", 200, 0.0, 200.0);
	NewHist( hdilep_tcmet[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_tcmet", dilepton_hypo_names[i]), "dilep_tcmet", 100, 0.0, 100.0);
	NewHist( hdilep_clmumet[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_calomet_muon", dilepton_hypo_names[i]), "dilep_calomet_muon", 100, 0.0, 100.0);
	NewHist( hdilep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet20", dilepton_hypo_names[i]), "dilep_njet20", 10, -0.5, 9.5);
	NewHist( hdilep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet30", dilepton_hypo_names[i]), "dilep_njet30", 10, -0.5, 9.5);
	NewHist( hdilep_reliso_ll[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_reliso_ll", dilepton_hypo_names[i]), "dilep_reliso_ll", 100, 0., 1.);

	//Z TH2F's for supp ABCD
	hdilep_ll_pt_eta[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_ll_pt_eta", dilepton_hypo_names[i]), "dilep_ll_pt_eta", metbins, 0.0, metmax, 48, -2.4, 2.4);
	hdilep_lt_pt_eta[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lt_pt_eta", dilepton_hypo_names[i]), "dilep_lt_pt_eta", metbins, 0.0, metmax, 48, -2.4, 2.4);
	hdilep_lt_ll20_pt_eta[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lt_ll20_pt_eta", dilepton_hypo_names[i]), "dilep_lt_ll20_pt_eta", metbins, 0.0, metmax, 48, -2.4, 2.4);
	//hdilep_lliso_pt_eta[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_pt_eta", dilepton_hypo_names[i]), "dilep_pt_eta", metbins, 0.0, metmax, 48, -2.4, 2.4);
	//hdilep_d0_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_d0_trckIso", dilepton_hypo_names[i]), "dilep_d0_trckIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	//hdilep_d0_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_d0_ecalIso", dilepton_hypo_names[i]), "dilep_d0_ecalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	//hdilep_d0_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_d0_hcalIso", dilepton_hypo_names[i]), "dilep_d0_hcalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	//hdilep_d0_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_d0_relIso",  dilepton_hypo_names[i]), "dilep_d0_relIso", 100, 0.0, d0max, 100, 0.0, 1.0);
	hdilep_lpt_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_trckIso", dilepton_hypo_names[i]), "dilep_lepPt_trckIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_ecalIso", dilepton_hypo_names[i]), "dilep_lepPt_ecalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_hcalIso", dilepton_hypo_names[i]), "dilep_lepPt_hcalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_relIso",  dilepton_hypo_names[i]), "dilep_lepPt_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//x-axis is met+lepPt
	hdilep_lepmet_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_trckIso", dilepton_hypo_names[i]), "dilep_lepMet_trckIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_ecalIso", dilepton_hypo_names[i]), "dilep_lepMet_ecalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_hcalIso", dilepton_hypo_names[i]), "dilep_lepMet_hcalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//same as above, only scale the pt by mw/mz
	hdilep_lpt_scl_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_Scl_trckIso", dilepton_hypo_names[i]), "dilep_lepPt_Scl_trckIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_scl_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_Scl_ecalIso", dilepton_hypo_names[i]), "dilep_lepPt_Scl_ecalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_scl_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_Scl_hcalIso", dilepton_hypo_names[i]), "dilep_lepPt_Scl_hcalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lpt_scl_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_Scl_relIso",  dilepton_hypo_names[i]), "dilep_lepPt_Scl_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//hdilep_lpt_rscl_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepPt_rScl_relIso",  dilepton_hypo_names[i]), "dilep_lepPt_rScl_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//x-axis is met+lepPt
	hdilep_lepmet_scl_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_trckIso", dilepton_hypo_names[i]), "dilep_lepMet_Scl_trckIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_scl_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_ecalIso", dilepton_hypo_names[i]), "dilep_lepMet_Scl_ecalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_scl_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_hcalIso", dilepton_hypo_names[i]), "dilep_lepMet_Scl_hcalIso", metbins, 0.0, metmax, 100, 0.0, 10.0);
	hdilep_lepmet_scl_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_Scl_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	hdilep_lepmet_rscl_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_rScl_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_rScl_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//require mc truth on scaled with met+lepPt
	hdilep_lepmet_scl_trth_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_trth_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_Scl_trth_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	hdilep_lepmet_rscl_trth_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_rScl_trth_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_rScl_trth_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//require transverse mass
	hdilep_lepmet_scl_tmas_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_tmas_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_Scl_tmas_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	hdilep_lepmet_rscl_tmas_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_rScl_tmas_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_rScl_tmas_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//require transverse mass + truth
	hdilep_lepmet_scl_tmast_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_tmast_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_Scl_tmast_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	hdilep_lepmet_rscl_tmast_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_rScl_tmast_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_rScl_tmast_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	//require mismeasurement be below some threshold, say 10% of gen
	hdilep_lepmet_scl_tmastmes_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_Scl_tmastmes_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_Scl_tmastmes_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
	hdilep_lepmet_rscl_tmastmes_relIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "dilep_lepMet_rScl_tmastmes_relIso",  dilepton_hypo_names[i]), "dilep_lepMet_rScl_tmastmes_relIso",  metbins, 0.0, metmax, 100, 0.0, 1.0);
  }
	
  // event level histograms
  //NewHist( hdilep_nhyp, Form("%s_%s_%s", SampleName().c_str(), "dilep_nhyp", "all"), "dilep_nhyp", 10, -0.5, 9.5);

}

cuts_t ABCDLooper::DilepSelect(const enum DileptonHypType myType, int idx1, int idx2) {
  cuts_t ret = 0; 
  float ptcut = 20.0;
  LorentzVector lep1_p4, lep2_p4;	

  //if( elidxs[0] != -1 && elidxs[1] != -1 ) {
  if( myType == DILEPTON_EE ) {

	lep1_p4 = cms2.els_p4()[idx1];
	lep2_p4 = cms2.els_p4()[idx2];
	 
	//if (cms2.hyp_lt_p4()[i_hyp].pt() > ptcut && cms2.hyp_ll_p4()[i_hyp].pt() > ptcut)
	if( cms2.els_p4()[idx1].pt() >= ptcut && cms2.els_p4()[idx2].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	//if( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 )
	if( cms2.els_charge()[idx1] * cms2.els_charge()[idx2] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.els_p4()[idx1] + cms2.els_p4()[idx2]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  //else if( muidxs[0] != -1 && muidxs[1] != -1 ) {
  else if( myType == DILEPTON_MUMU ) {

	lep1_p4 = cms2.mus_p4()[idx1];
	lep2_p4 = cms2.mus_p4()[idx2];

	if( cms2.mus_p4()[idx1].pt() >= ptcut && cms2.mus_p4()[idx2].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	if( cms2.mus_charge()[idx1] * cms2.mus_charge()[idx2] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.mus_p4()[idx1] + cms2.mus_p4()[idx2]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

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

cuts_t ABCDLooper::LepSelect(int lep_type, int i)
{
  cuts_t ret = 0;
  float ptcut = 20.0;
  LorentzVector lep_p4;
    
  // e
  if( lep_type == 0 ) {

	lep_p4 = cms2.els_p4()[i];

	if( cms2.els_p4()[i].pt() >= ptcut )
	  ret |= CUT_BIT(LEP_PT);

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

	if( cms2.mus_p4()[i].pt() >= ptcut )
	  ret |= CUT_BIT(LEP_PT);

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
  //if (calo_jets.size() > 1) {
  //sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  //}
  //return calo_jets;

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
	ret |= CUT_BIT(TMASS);


  return ret;
}

void ABCDLooper::FillEventHistos ()
{

  // get the event weight
  if( sample_.kFactor != 1 ) cout << "kFactor non-unity " << sample_.kFactor << endl;
  double weight = ABCDLooper::Weight();
  
  // need to determine if this is a di-lepton
  // or a single lepton event
  int nels = 0, nmus = 0;
  int nels_nopt = 0, nmus_nopt = 0;
  int nels_noiso = 0, nmus_noiso = 0;
  int nels_nod0iso = 0, nmus_nod0iso = 0;
  int nels_noptiso = 0, nmus_noptiso = 0;
  //int nels_nojets20 = 0, nmus_nojets20 = 0;
  int njetels_20 = 0, njetmus_20 = 0;
  //int nels_nojets30 = 0, nmus_nojets30 = 0; //just count 20 for now
  elidxs[0] = elidxs[1] = -1;
  muidxs[0] = muidxs[1] = -1;
  int elidxs_nopt[] = {-1, -1}; 
  int muidxs_nopt[] = {-1, -1};
  int elidxs_noiso[] = {-1, -1};
  int muidxs_noiso[] = {-1, -1};
  //int elidxs_nod0iso[] = {-1, -1};
  //int muidxs_nod0iso[] = {-1, -1};
  int elidxs_noptiso[] = {-1, -1};
  int muidxs_noptiso[] = {-1, -1};
  cuts_t pass_lep_cut = 0; //cuts of lepton which passes lepcuts--only good for single lepton
  cuts_t pass_nlep_cuts = 0; //cuts of all leptons which pass lepcuts: & of individual
  int pass_elcuts_noptiso[] = {0, 0}; //fill based on baselepcuts_noptiso for aviBCD
  int pass_mucuts_noptiso[] = {0, 0};
  int pass_elcuts_noiso[] = {0, 0}; //fill based on baselepcuts_noiso w abcd
  int pass_mucuts_noiso[] = {0, 0};
  //DEFINE CUTS
  //put a jet veto on everything
  // need to count leps before jets: if this cut gives 1 lep, use its cut bit already set in LepSelect for jet reqmnt. If 2 leps, don't fill W abcd
  cuts_t baselepcuts 			= ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) ); //FOR BASE LEP COUNTING
  cuts_t baselepcuts_noiso		= ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) ); //for w abcd lep counting
  cuts_t baselepcuts_nopt  		= ( CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) ); //cut for pt plot, also no mass, met cuts
  cuts_t baselepcuts_nod0iso    = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD_NOD0) );
  cuts_t baselepcuts_noptiso  	= ( CUT_BIT(LEP_GOOD) ); //for z abcd lep counting
  //N-1 plots for single leptons: require w cuts--for (d0,iso) plots, relax (d0,iso), apply everything else including tmass, tcmet
  cuts_t w_evt_sel      = ( CUT_BIT(TCMET) | CUT_BIT(TMASS) | CUT_BIT(JET_VETO_20) );
  cuts_t w_cuts_nod0    = ( w_evt_sel | CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD_NOD0) );
  cuts_t w_cuts_noiso   = ( w_evt_sel | CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) );
  cuts_t w_cuts_nod0iso = ( w_evt_sel | CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD_NOD0) );
  cuts_t w_cuts_notm    = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD) | CUT_BIT(TCMET) | CUT_BIT(JET_VETO_20) );
  cuts_t w_cuts_nomet   = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD) | CUT_BIT(TMASS) | CUT_BIT(JET_VETO_20) );
  cuts_t w_cuts_nojet20 = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD) | CUT_BIT(TCMET) | CUT_BIT(TMASS) ); 

  cuts_t w_cuts = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | w_evt_sel | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_20) ); //ALL W CUTS

  cuts_t z_basecuts     = ( CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_MASS) | CUT_BIT(DILEP_OS) );
  cuts_t z_cuts         = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | z_basecuts | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_20) ); //ALL Z CUTS
  cuts_t z_cuts_nomass  = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_OS) | CUT_BIT(LEP_ISO) | CUT_BIT(JET_VETO_20) ); //all but zmass
  cuts_t z_cuts_nojet20 = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | z_basecuts | CUT_BIT(LEP_ISO) ); //all Z cuts but njet

  //cuts for ABCD: CUT ON TMASS, MET
  cuts_t w_cuts_nometiso = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(TMASS) | CUT_BIT(JET_VETO_20) ); //for W only
  cuts_t z_cuts_nometiso = ( CUT_BIT(LEP_GOOD) | CUT_BIT(DILEP_MASS) | CUT_BIT(DILEP_OS) | CUT_BIT(JET_VETO_20) ); //for Z only : met means lep pt here

  //check tcMET here b'c no point doing for every lepton--don't use for Z
  cuts_t tcmetcut = 0;
  if( cms2.evt_tcmet() > 20 )
	tcmetcut |= CUT_BIT(TCMET);

  //select els
  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	transmass = 0;
	cuts_t elcut = LepSelect(0, i); //0 for els
	elcut |= tcmetcut;
	if( (elcut & baselepcuts) == baselepcuts ) {
	  nels++;
	  pass_lep_cut = elcut;
	  njetels_20 = njets_20; //njets_20 is set in LepSelect
	  if( elidxs[0] == -1 ) {
		elidxs[0] = i;
		pass_nlep_cuts = pass_lep_cut;
	  }
	  else if( elidxs[1] == -1 ) {
		elidxs[1] = i;
		pass_nlep_cuts &= pass_lep_cut;
	  }
	  //if > 2 els, ignore the rest--should be sorted by pt already
	}
	if( (elcut & baselepcuts_nopt) == baselepcuts_nopt ) { //no pt cut
	  nels_nopt++;
	  if( elidxs_nopt[0] == -1 )
		elidxs_nopt[0] = i;
	  else if( elidxs_nopt[1] == -1 )
		elidxs_nopt[1] = i;
	}
	if( (elcut & baselepcuts_noiso) == baselepcuts_noiso ) { //no iso cut
	  nels_noiso++;
	  if( elidxs_noiso[0] == -1 ) {
		elidxs_noiso[0] = i;
		pass_elcuts_noiso[0] = elcut;
	  }
	  else if( elidxs_noiso[1] == -1 ) {
		elidxs_noiso[1] = i;
		pass_elcuts_noiso[1] = elcut;
	  }
	}
	if( (elcut & baselepcuts_noptiso) == baselepcuts_noptiso ) { //no pt nor iso cuts
	  nels_noptiso++;
	  if( elidxs_noptiso[0] == -1 ) {
		elidxs_noptiso[0] = i;
		pass_elcuts_noptiso[0] = elcut;
	  }
	  else if( elidxs_noptiso[1] == -1 ) {
		elidxs_noptiso[1] = i;
		pass_elcuts_noptiso[1] = elcut;
	  }
	}
	if( (elcut & w_cuts_noiso) == w_cuts_noiso ) {
	  hlep_trckIso[0]->Fill( cms2.els_tkIso()[i], weight );
	  hlep_trckIso[2]->Fill( cms2.els_tkIso()[i], weight );
	  hlep_ecalIso[0]->Fill( cms2.els_pat_ecalIso()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.els_pat_ecalIso()[i], weight );
	  hlep_hcalIso[0]->Fill( cms2.els_pat_hcalIso()[i], weight );
	  hlep_hcalIso[2]->Fill( cms2.els_pat_hcalIso()[i], weight );
	  hlep_relIso[0]->Fill( inv_el_relsusy_iso(i, true), weight );
	  hlep_relIso[2]->Fill( inv_el_relsusy_iso(i, true), weight );
	}
	if( (elcut & w_cuts_nod0) == w_cuts_nod0 ) {
	  hlep_d0[0]->Fill( fabs(cms2.els_d0corr()[i]), weight );
	  hlep_d0[2]->Fill( fabs(cms2.els_d0corr()[i]), weight );
	  hlep_d0Sig[0]->Fill( fabs(cms2.els_d0corr()[i]/cms2.els_d0Err()[i]), weight );
	  hlep_d0Sig[2]->Fill( fabs(cms2.els_d0corr()[i]/cms2.els_d0Err()[i]), weight );
	}
	//for mass plots, check all but mass
	if( (elcut & w_cuts_notm) == w_cuts_notm ) {
	  hlep_mass[0]->Fill( transmass, weight ); //check all cuts but mass for mass plot (n-1)
	  hlep_mass[2]->Fill( transmass, weight );
	}
	//ABCD in d0, Iso
	if( (elcut & baselepcuts_nod0iso) == baselepcuts_nod0iso ) {
	  nels_nod0iso++;
	  hlep_d0_trckIso[0]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_tkIso()[i], weight );
	  hlep_d0_trckIso[2]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_tkIso()[i], weight );
	  hlep_d0_ecalIso[0]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_pat_ecalIso()[i], weight );
	  hlep_d0_ecalIso[2]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_pat_ecalIso()[i], weight );
	  hlep_d0_hcalIso[0]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_pat_hcalIso()[i], weight );
	  hlep_d0_hcalIso[2]->Fill( fabs(cms2.els_d0corr()[i]), cms2.els_pat_hcalIso()[i], weight );
	  hlep_d0_relIso[0] ->Fill( fabs(cms2.els_d0corr()[i]), inv_el_relsusy_iso(i, true), weight );
	  hlep_d0_relIso[2] ->Fill( fabs(cms2.els_d0corr()[i]), inv_el_relsusy_iso(i, true), weight );
	}
	//if( (elcut & w_cuts_nojet20) == w_cuts_nojets20 ) {
	//  nels_nojets20++;
	//  njetels_20 = njets_20;
	//}
  }

  //select mus
  for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
	transmass = 0;
	cuts_t mucut = LepSelect(1, i); //1 for mus
	mucut |= tcmetcut;
	if( (mucut & baselepcuts) == baselepcuts ) { //all cuts
	  nmus++;
	  pass_lep_cut = mucut;
	  njetmus_20 = njets_20;
	  if( muidxs[0] == -1 ) {
		muidxs[0] = i;
		pass_nlep_cuts = pass_lep_cut;
	  }
	  else if( muidxs[1] == -1 ) {
		muidxs[1] = i;
		pass_nlep_cuts &= pass_lep_cut;
	  }
	}
	if( (mucut & baselepcuts_nopt) == baselepcuts_nopt ) { //no pt cut
	  nmus_nopt++;
	  if( muidxs_nopt[0] == -1 )
		muidxs_nopt[0] = i;
	  else if( muidxs_nopt[1] == -1 )
		muidxs_nopt[1] = i;
	}
	if( (mucut & baselepcuts_noiso) == baselepcuts_noiso ) { //no iso cut
	  nmus_noiso++;
	  if( muidxs_noiso[0] == -1 ) {
		muidxs_noiso[0] = i;
		pass_mucuts_noiso[0] = mucut;
	  }
	  else if( muidxs_noiso[1] == -1 ) {
		muidxs_noiso[1] = i;
		pass_mucuts_noiso[1] = mucut;
	  }
	}
	if( (mucut & baselepcuts_noptiso) == baselepcuts_noptiso ) { //no iso,pt cut
	  nmus_noptiso++;
	  if( muidxs_noptiso[0] == -1 ) {
		muidxs_noptiso[0] = i;
		pass_mucuts_noptiso[0] = mucut;
	  }
	  else if( muidxs_noptiso[1] == -1 ) {
		muidxs_noptiso[1] = i;
		pass_mucuts_noptiso[1] = mucut;
	  }
	}
	if( (mucut & w_cuts_noiso) == w_cuts_noiso ) {
	  hlep_trckIso[1]->Fill( cms2.mus_iso03_sumPt()[i], weight);
	  hlep_trckIso[2]->Fill( cms2.mus_iso03_sumPt()[i], weight);
	  hlep_ecalIso[1]->Fill( cms2.mus_iso03_emEt()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.mus_iso03_emEt()[i], weight );
	  hlep_hcalIso[1]->Fill( cms2.mus_iso03_hadEt()[i], weight );
	  hlep_hcalIso[2]->Fill( cms2.mus_iso03_hadEt()[i], weight );
	  hlep_relIso[1]->Fill( inv_mu_relsusy_iso(i), weight );
	  hlep_relIso[2]->Fill( inv_mu_relsusy_iso(i), weight );
	}
	if( (mucut & w_cuts_nod0) == w_cuts_nod0 ) {
	  hlep_d0[1]->Fill( fabs(cms2.mus_d0corr()[i]), weight );
	  hlep_d0[2]->Fill( fabs(cms2.mus_d0corr()[i]), weight );
	  hlep_d0Sig[1]->Fill( fabs(cms2.mus_d0corr()[i]/cms2.mus_d0Err()[i]), weight );
	  hlep_d0Sig[2]->Fill( fabs(cms2.mus_d0corr()[i]/cms2.mus_d0Err()[i]), weight );
	  //correct formula for d0 sig:
	  //els_d0corr / sqrt(els_d0Err*els_d0Err + evt_bs_width *evt_bs_width)
	}
	//for mass plots, check all but mass
	if( (mucut & w_cuts_notm) == w_cuts_notm ) {
	  hlep_mass[1]->Fill( transmass, weight ); //check all cuts but mass for mass plot (n-1)
	  hlep_mass[2]->Fill( transmass, weight );
	}
	//ABCD in d0, Iso
	if( (mucut & w_cuts_nod0iso) == w_cuts_nod0iso ) {
	  nmus_nod0iso++;
	  hlep_d0_trckIso[1]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_sumPt()[i], weight );
	  hlep_d0_trckIso[2]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_sumPt()[i], weight );
	  hlep_d0_ecalIso[1]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_emEt()[i], weight );
	  hlep_d0_ecalIso[2]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_emEt()[i], weight );
	  hlep_d0_hcalIso[1]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_hadEt()[i], weight );
	  hlep_d0_hcalIso[2]->Fill( fabs(cms2.mus_d0corr()[i]), cms2.mus_iso03_hadEt()[i], weight );
	  hlep_d0_relIso[1] ->Fill( fabs(cms2.mus_d0corr()[i]), inv_mu_relsusy_iso(i), weight );
	  hlep_d0_relIso[2] ->Fill( fabs(cms2.mus_d0corr()[i]), inv_mu_relsusy_iso(i), weight );
	}
	//if( (mucut & w_cuts_nojets20) == w_cuts_nojets20 ) {
	//  nmus_nojets20++;
	//  njetmus_20 = njets_20;
	//}
  }

  //fill nlep cut before checking nels, nmus
  //0=el, 1=mu, 2=all
  hlep_nlep[0]->Fill( nels, weight );
  hlep_nlep[1]->Fill( nmus, weight );
  hlep_nlep[2]->Fill( nels+nmus, weight );

  hlep_nlep_nod0iso[0]->Fill( nels_nod0iso, weight );
  hlep_nlep_nod0iso[1]->Fill( nmus_nod0iso, weight );
  hlep_nlep_nod0iso[2]->Fill( nels_nod0iso + nmus_nod0iso, weight );
  //hlep_nlep_nometiso[0]->Fill( nels_nometiso, weight ); //these no longer necessary since now lep counting is prereq for this abcd
  //hlep_nlep_nometiso[1]->Fill( nmus_nometiso, weight );
  //hlep_nlep_nometiso[2]->Fill( nels_nometiso + nmus_nometiso, weight );

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

  //find tight leg of Z for avi abcd
  cuts_t zabcd = 0;
  if( (nels_noptiso == 2 && nmus_noptiso == 0) || (nmus_noptiso == 2 && nels_noptiso == 0) ) { //won't use event if has op flv
	const enum DileptonHypType myType = (nels_noptiso == 2 ? DILEPTON_EE : DILEPTON_MUMU);
	int idx1 = (nels_noptiso == 2 ? elidxs_noptiso[0] : muidxs_noptiso[0]);
	int idx2 = (nels_noptiso == 2 ? elidxs_noptiso[1] : muidxs_noptiso[1]);
	zabcd = DilepSelect(myType, idx1, idx2); 

	if( myType == DILEPTON_EE ) {
	  //if both pass tight, 50/50 chance to flip: lep idx 1 will be used as met
	  if( (pass_elcuts_noptiso[0] & pass_elcuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso ) {
	  	if( cms2.evt_event()%2 == 0 ) {
	  	  int tmp = elidxs_noptiso[0];
	  	  elidxs_noptiso[0] = elidxs_noptiso[1];
	  	  elidxs_noptiso[1] = tmp;
	  	}
	  }
	  else if( ((pass_elcuts_noptiso[0] & baselepcuts_noiso) == baselepcuts_noiso) &&
	  		   ((pass_elcuts_noptiso[1] & baselepcuts_noiso) != baselepcuts_noiso) ) { } //idx 1 is loose--ok
	  else if( ((pass_elcuts_noptiso[0] & baselepcuts_noiso) != baselepcuts_noiso) &&
			   ((pass_elcuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso) ) {
	  //REQUIRE loose leg to pass iso
	  //if( ((pass_elcuts_noptiso[0] & baselepcuts_noiso) == baselepcuts_noiso) && //idx 0 is tight: use for iso
	  //  ((pass_elcuts_noptiso[1] & baselepcuts_nopt) == baselepcuts_nopt) ) { } //idx 1 is loose--ok
	  //else if( ((pass_elcuts_noptiso[0] & baselepcuts_nopt) == baselepcuts_nopt) &&
	  //		   ((pass_elcuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso) ) {
	  	int tmp = elidxs_noptiso[0];
	  	elidxs_noptiso[0] = elidxs_noptiso[1];
	  	elidxs_noptiso[1] = tmp;
	  }
	  else //both fail
		zabcd = 0;
	  
	  zabcd |= ( pass_elcuts_noptiso[0] | pass_elcuts_noptiso[1] ); //this is ok b'c will fail mass cut from dilepselect
	}
	else {
	  //if both pass tight, 50/50 chance to flip: lep idx 1 will be used as met
	  if( (pass_mucuts_noptiso[0] & pass_mucuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso ) {
	  	if( cms2.evt_event()%2 == 0 ) {
	  	  int tmp = muidxs_noptiso[0];
	  	  muidxs_noptiso[0] = muidxs_noptiso[1];
	  	  muidxs_noptiso[1] = tmp;
	  	}
	  }
	  else if( ((pass_mucuts_noptiso[0] & baselepcuts_noiso) == baselepcuts_noiso) &&
	  		   ((pass_mucuts_noptiso[1] & baselepcuts_noiso) != baselepcuts_noiso) ) { } //idx 1 is loose--ok
	  else if( ((pass_mucuts_noptiso[0] & baselepcuts_noiso) != baselepcuts_noiso) &&
	  		   ((pass_mucuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso) ) {
	  //REQUIRE loose leg to pass iso
	  //if( ((pass_mucuts_noptiso[0] & baselepcuts_noiso) == baselepcuts_noiso) &&
	  //	  ((pass_mucuts_noptiso[1] & baselepcuts_nopt) == baselepcuts_nopt) ) { } 
	  //else if( ((pass_mucuts_noptiso[0] & baselepcuts_nopt) == baselepcuts_nopt) &&
	  //	   ((pass_mucuts_noptiso[1] & baselepcuts_noiso) == baselepcuts_noiso) ) {
		int tmp = muidxs_noptiso[0];
		muidxs_noptiso[0] = muidxs_noptiso[1];
		muidxs_noptiso[1] = tmp;
	  }
	  else //both fail
		zabcd = 0;
	  
	  zabcd |= ( pass_mucuts_noptiso[0] | pass_mucuts_noptiso[1] );
	}
  }


  //generator quantities
  double wmass = 80.4, zmass = 91.19;
  double scale = wmass/zmass; 
  double rscale = 1.031*scale; //refine: increase by ratio of means: 32.9/31.9 == 1.031
  //double wwin_low = 65., wwin_hgh = 95.;
  //double zwin_low = 76., zwin_hgh = 106.;
  double win_width = 15.;
  double wwin_low = wmass - win_width, wwin_hgh = wmass + win_width;
  double zwin_low = zmass - win_width, zwin_hgh = zmass + win_width;
  double gen_lt_pt=0, gen_ll_pt=0; //only for truth matched to leps selected above
  double gen_boson_mass=0;
  for( int i=0;i<(int)cms2.genps_p4().size();i++ ) {
	if( abs(cms2.genps_id()[i]) == 23 || abs(cms2.genps_id()[i]) == 24 ) //Z==23, W==24
	  gen_boson_mass = cms2.genps_p4()[i].mass();
	else if( abs(cms2.genps_id()[i]) == 11 ) { //e
	  if( elidxs_noptiso[0] != -1 && cms2.els_mc3idx()[elidxs_noptiso[0]] == i ) //truth match
		gen_lt_pt = cms2.genps_p4()[i].pt();
	  else if( elidxs_noptiso[1] != -1 && cms2.els_mc3idx()[elidxs_noptiso[1]] == i ) //truth match
		gen_ll_pt = cms2.genps_p4()[i].pt();
	  //for gen alone
	  if( abs(cms2.genps_id_mother()[i]) == 23 && gen_boson_mass > zwin_low && gen_boson_mass < zwin_hgh ) {
		hdilep_genpt[3]->Fill( scale*cms2.genps_p4()[i].pt(), weight ); //scale Zs
		hdilep_genpt[0]->Fill( scale*cms2.genps_p4()[i].pt(), weight ); //scale Zs
	  }
	  else if( abs(cms2.genps_id_mother()[i]) == 24 && gen_boson_mass > wwin_low && gen_boson_mass < wwin_hgh ) {
		hlep_genpt[0]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
		hlep_genpt[2]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
	  }
	}
	else if( abs(cms2.genps_id()[i]) == 13 ) { //mu
	  if( muidxs_noptiso[0] != -1 && cms2.mus_mc3idx()[muidxs_noptiso[0]] == i ) //truth match
		gen_lt_pt = cms2.genps_p4()[i].pt();
	  else if( muidxs_noptiso[1] != -1 && cms2.mus_mc3idx()[muidxs_noptiso[1]] == i ) //truth match
		gen_ll_pt = cms2.genps_p4()[i].pt();
	  //for gen alone
	  if( abs(cms2.genps_id_mother()[i]) == 23 && gen_boson_mass > zwin_low && gen_boson_mass < zwin_hgh ) {
		hdilep_genpt[1]->Fill( scale*cms2.genps_p4()[i].pt(), weight ); //scale Zs
		hdilep_genpt[0]->Fill( scale*cms2.genps_p4()[i].pt(), weight ); //scale Zs
	  }
	  else if( abs(cms2.genps_id_mother()[i]) == 24 && gen_boson_mass > wwin_low && gen_boson_mass < wwin_hgh ) {
		hlep_genpt[1]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
		hlep_genpt[2]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
	  }
	}
	else if( abs(cms2.genps_id()[i]) == 12 ) { //nu_e
	  if( abs(cms2.genps_id_mother()[i]) == 24 && gen_boson_mass > wwin_low && gen_boson_mass < wwin_hgh ) {
		hlep_genpt[0]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
		hlep_genpt[2]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
	  }
	}
	else if( abs(cms2.genps_id()[i]) == 14 ) { //nu_mu
	  if( abs(cms2.genps_id_mother()[i]) == 24 && gen_boson_mass > wwin_low && gen_boson_mass < wwin_hgh ) {
		hlep_genpt[1]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
		hlep_genpt[2]->Fill( cms2.genps_p4()[i].pt(), weight ); //don't scale Ws
	  }
	}
  }
  
  //ABCD for Zs
  //apply gen level mass cut
  if( (zabcd & z_cuts_nometiso) == z_cuts_nometiso ) {
	if( nels_noptiso == 2 && nmus_noptiso == 0 ) {
	  //met is from idx 1, iso is from idx 0 b'c this is the lep
	  double pt = cms2.els_p4()[elidxs_noptiso[1]].pt();
	  double lepmet = (LorentzVector( cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()),
									  cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()),
									  0, cms2.evt_tcmet()) + cms2.els_p4()[elidxs_noptiso[1]] ).pt();
	  hdilep_lpt_trckIso[3]				->Fill( pt, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lpt_trckIso[DILEPTON_ALL]	->Fill( pt, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lpt_ecalIso[3]				->Fill( pt, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_ecalIso[DILEPTON_ALL]	->Fill( pt, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_hcalIso[3]				->Fill( pt, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_hcalIso[DILEPTON_ALL]	->Fill( pt, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_relIso[3] 				->Fill( pt, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lpt_relIso[DILEPTON_ALL] 	->Fill( pt, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_trckIso[3]				->Fill( lepmet, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lepmet_trckIso[DILEPTON_ALL]	->Fill( lepmet, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lepmet_ecalIso[3]				->Fill( lepmet, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_ecalIso[DILEPTON_ALL]	->Fill( lepmet, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_hcalIso[3]				->Fill( lepmet, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_hcalIso[DILEPTON_ALL]	->Fill( lepmet, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_relIso[3] 				->Fill( lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_relIso[DILEPTON_ALL] 	->Fill( lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  //plots with scale factor
	  hdilep_lpt_scl_trckIso[3]				->Fill( scale*pt, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lpt_scl_trckIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lpt_scl_ecalIso[3]				->Fill( scale*pt, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_ecalIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_hcalIso[3]				->Fill( scale*pt, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_hcalIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_relIso[3] 				->Fill( scale*pt, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lpt_scl_relIso[DILEPTON_ALL] 	->Fill( scale*pt, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_scl_trckIso[3]				->Fill( scale*lepmet, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lepmet_scl_trckIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.els_tkIso()[elidxs_noptiso[0]], weight );		
	  hdilep_lepmet_scl_ecalIso[3]				->Fill( scale*lepmet, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_ecalIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.els_pat_ecalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_hcalIso[3]				->Fill( scale*lepmet, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_hcalIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.els_pat_hcalIso()[elidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_relIso[3] 				->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_scl_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_rscl_relIso[3] 				->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  hdilep_lepmet_rscl_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
	  //apply transverse mass cut on Z to match W selection
	  double dphi = ROOT::Math::VectorUtil::DeltaPhi( cms2.els_p4()[elidxs_noptiso[0]], cms2.els_p4()[elidxs_noptiso[1]] );
	  double masst = sqrt( 2 * scale * scale * pt * cms2.els_p4()[elidxs_noptiso[0]].pt() * ( 1 - cos(dphi) ) );
	  hdilep_tmass[3]->Fill( masst, weight );
	  hdilep_tmass[DILEPTON_ALL]->Fill( masst, weight );
	  if( masst > 40 && masst < 100 ) {
		hdilep_lepmet_scl_tmas_relIso[3] 				->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		hdilep_lepmet_scl_tmas_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		hdilep_lepmet_rscl_tmas_relIso[3] 				->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		hdilep_lepmet_rscl_tmas_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		hdilep_ll_pt_eta[3]           ->Fill( scale*pt, cms2.els_p4()[elidxs_noptiso[1]].eta(), weight );
		hdilep_ll_pt_eta[DILEPTON_ALL]->Fill( scale*pt, cms2.els_p4()[elidxs_noptiso[1]].eta(), weight );
		hdilep_lt_pt_eta[3]           ->Fill( scale*cms2.els_p4()[elidxs_noptiso[0]].pt(), cms2.els_p4()[elidxs_noptiso[0]].eta(), weight );
		hdilep_lt_pt_eta[DILEPTON_ALL]->Fill( scale*cms2.els_p4()[elidxs_noptiso[0]].pt(), cms2.els_p4()[elidxs_noptiso[0]].eta(), weight );
		hdilep_reliso_ll[3]           ->Fill( inv_el_relsusy_iso(elidxs_noptiso[1], true), weight );
		hdilep_reliso_ll[DILEPTON_ALL]->Fill( inv_el_relsusy_iso(elidxs_noptiso[1], true), weight );
		hdilep_0_genpt_mch[3]           ->Fill( scale*gen_lt_pt, weight);
		hdilep_0_genpt_mch[DILEPTON_ALL]->Fill( scale*gen_lt_pt, weight);
		hdilep_1_genpt_mch[3]           ->Fill( scale*gen_ll_pt, weight);
		hdilep_1_genpt_mch[DILEPTON_ALL]->Fill( scale*gen_ll_pt, weight);
		if( abs(cms2.els_mc3_id()[elidxs_noptiso[1]]) == 11 ) { //truth match met leg of z
		  hdilep_lepmet_scl_tmast_relIso[3] 			->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		  hdilep_lepmet_scl_tmast_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		  hdilep_lepmet_rscl_tmast_relIso[3] 			->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		  hdilep_lepmet_rscl_tmast_relIso[DILEPTON_ALL]	->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		  if( fabs(cms2.els_mc3_p4()[elidxs_noptiso[1]].pt() - pt)/cms2.els_mc3_p4()[elidxs_noptiso[1]].pt() < 0.1 ) { //pt measured to 10% of gen pt
			hdilep_lepmet_scl_tmastmes_relIso[3] 			->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
			hdilep_lepmet_scl_tmastmes_relIso[DILEPTON_ALL] ->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
			hdilep_lepmet_rscl_tmastmes_relIso[3] 			->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
			hdilep_lepmet_rscl_tmastmes_relIso[DILEPTON_ALL]->Fill( rscale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		  }
		}
	  }
	  if( abs(cms2.els_mc3_id()[elidxs_noptiso[1]]) == 11 ) { //truth match met leg of z
		hdilep_lepmet_scl_trth_relIso[3]		   ->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		hdilep_lepmet_scl_trth_relIso[DILEPTON_ALL]->Fill( scale*lepmet, inv_el_relsusy_iso(elidxs_noptiso[0], true), weight );
		if( scale*pt < 20 ) { //loose leg is in the low pt range
		  double ptgen = cms2.genps_p4()[cms2.els_mc3idx()[elidxs_noptiso[1]]].pt();
		  hdilep_pt_mgen[3]->Fill( (cms2.els_p4()[elidxs_noptiso[1]].pt()-ptgen)/ptgen, weight );
		  hdilep_pt_mgen[DILEPTON_ALL]->Fill( (cms2.els_p4()[elidxs_noptiso[1]].pt()-ptgen)/ptgen, weight );
		}
	  }
	  if( scale*pt < 20 ) { //loose leg is in the low pt range
		hdilep_lt_ll20_pt_eta[3]           ->Fill( scale*cms2.els_p4()[elidxs_noptiso[0]].pt(), cms2.els_p4()[elidxs_noptiso[0]].eta(), weight );
		hdilep_lt_ll20_pt_eta[DILEPTON_ALL]->Fill( scale*cms2.els_p4()[elidxs_noptiso[0]].pt(), cms2.els_p4()[elidxs_noptiso[0]].eta(), weight );
		double mass = (cms2.els_p4()[elidxs_noptiso[0]] + cms2.els_p4()[elidxs_noptiso[1]]).mass();
		hdilep_mass_ll20[3]->Fill( mass, weight );
		hdilep_mass_ll20[DILEPTON_ALL]->Fill( mass, weight );
	  }
	}
	else if( nmus_noptiso == 2 && nels_noptiso == 0 ) {
	  double pt = cms2.mus_p4()[muidxs_noptiso[1]].pt();
	  double lepmet = (LorentzVector( cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()),
									  cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()),
									  0, cms2.evt_tcmet()) + cms2.mus_p4()[muidxs_noptiso[1]] ).pt();
	  hdilep_lpt_trckIso[1]				->Fill( pt, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lpt_trckIso[DILEPTON_ALL]	->Fill( pt, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lpt_ecalIso[1]				->Fill( pt, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_ecalIso[DILEPTON_ALL]	->Fill( pt, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_hcalIso[1]				->Fill( pt, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_hcalIso[DILEPTON_ALL]	->Fill( pt, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_relIso[1] 				->Fill( pt, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lpt_relIso[DILEPTON_ALL] 	->Fill( pt, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_trckIso[1]				->Fill( lepmet, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lepmet_trckIso[DILEPTON_ALL]	->Fill( lepmet, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lepmet_ecalIso[1]				->Fill( lepmet, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_ecalIso[DILEPTON_ALL]	->Fill( lepmet, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_hcalIso[1]				->Fill( lepmet, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_hcalIso[DILEPTON_ALL]	->Fill( lepmet, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_relIso[1] 				->Fill( lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_relIso[DILEPTON_ALL] 	->Fill( lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  //scale plots
	  hdilep_lpt_scl_trckIso[1]				->Fill( scale*pt, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lpt_scl_trckIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lpt_scl_ecalIso[1]				->Fill( scale*pt, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_ecalIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_hcalIso[1]				->Fill( scale*pt, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_hcalIso[DILEPTON_ALL]	->Fill( scale*pt, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lpt_scl_relIso[1] 				->Fill( scale*pt, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lpt_scl_relIso[DILEPTON_ALL] 	->Fill( scale*pt, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_scl_trckIso[1]				->Fill( scale*lepmet, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lepmet_scl_trckIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.mus_iso03_sumPt()[muidxs_noptiso[0]], weight );		
	  hdilep_lepmet_scl_ecalIso[1]				->Fill( scale*lepmet, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_ecalIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.mus_iso03_emEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_hcalIso[1]				->Fill( scale*lepmet, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_hcalIso[DILEPTON_ALL]	->Fill( scale*lepmet, cms2.mus_iso03_hadEt()[muidxs_noptiso[0]], weight );	
	  hdilep_lepmet_scl_relIso[1] 				->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_scl_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_rscl_relIso[1] 				->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  hdilep_lepmet_rscl_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
	  //apply transverse mass cut on Z to match W selection
	  double dphi = ROOT::Math::VectorUtil::DeltaPhi( cms2.mus_p4()[muidxs_noptiso[0]], cms2.mus_p4()[muidxs_noptiso[1]] );
	  double masst = sqrt( 2 * scale * scale * pt * cms2.mus_p4()[muidxs_noptiso[0]].pt() * ( 1 - cos(dphi) ) );
	  hdilep_tmass[1]->Fill( masst, weight );
	  hdilep_tmass[DILEPTON_ALL]->Fill( masst, weight );
	  if( masst > 40 && masst < 100 ) {
		hdilep_lepmet_scl_tmas_relIso[1] 				->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_scl_tmas_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_rscl_tmas_relIso[1] 				->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_rscl_tmas_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_ll_pt_eta[1]           ->Fill( scale*pt, cms2.mus_p4()[muidxs_noptiso[1]].eta(), weight );
		hdilep_ll_pt_eta[DILEPTON_ALL]->Fill( scale*pt, cms2.mus_p4()[muidxs_noptiso[1]].eta(), weight );
		hdilep_lt_pt_eta[1]           ->Fill( scale*cms2.mus_p4()[muidxs_noptiso[0]].pt(), cms2.mus_p4()[muidxs_noptiso[0]].eta(), weight );
		hdilep_lt_pt_eta[DILEPTON_ALL]->Fill( scale*cms2.mus_p4()[muidxs_noptiso[0]].pt(), cms2.mus_p4()[muidxs_noptiso[0]].eta(), weight );
		hdilep_reliso_ll[1]           ->Fill( inv_mu_relsusy_iso(muidxs_noptiso[1]), weight );
		hdilep_reliso_ll[DILEPTON_ALL]->Fill( inv_mu_relsusy_iso(muidxs_noptiso[1]), weight );
		hdilep_0_genpt_mch[1]           ->Fill( gen_lt_pt, weight);
		hdilep_0_genpt_mch[DILEPTON_ALL]->Fill( gen_lt_pt, weight);
		hdilep_1_genpt_mch[1]           ->Fill( gen_ll_pt, weight);
		hdilep_1_genpt_mch[DILEPTON_ALL]->Fill( gen_ll_pt, weight);
		if( abs(cms2.mus_mc3_id()[muidxs_noptiso[1]]) == 13 ) { //truth match met leg
		  hdilep_lepmet_scl_tmast_relIso[1]				->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		  hdilep_lepmet_scl_tmast_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		  hdilep_lepmet_rscl_tmast_relIso[1]				->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		  hdilep_lepmet_rscl_tmast_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		  if( fabs(cms2.mus_mc3_p4()[muidxs_noptiso[1]].pt() - pt)/cms2.mus_mc3_p4()[muidxs_noptiso[1]].pt() < 0.1 ) { //pt measured to 10% of gen pt
			hdilep_lepmet_scl_tmastmes_relIso[3] 			->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
			hdilep_lepmet_scl_tmastmes_relIso[DILEPTON_ALL] ->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
			hdilep_lepmet_rscl_tmastmes_relIso[3] 			->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
			hdilep_lepmet_rscl_tmastmes_relIso[DILEPTON_ALL]->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		  }
		}
	  }
	  if( abs(cms2.mus_mc3_id()[muidxs_noptiso[1]]) == 13 ) { //truth match met leg
		hdilep_lepmet_scl_trth_relIso[1] 				->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_scl_trth_relIso[DILEPTON_ALL] 	->Fill( scale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_rscl_trth_relIso[1] 				->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		hdilep_lepmet_rscl_trth_relIso[DILEPTON_ALL] 	->Fill( rscale*lepmet, inv_mu_relsusy_iso(muidxs_noptiso[0]), weight );
		if( scale*pt < 20 ) { //loose leg is in the low pt range
		  double ptgen = cms2.genps_p4()[cms2.mus_mc3idx()[muidxs_noptiso[1]]].pt();
		  hdilep_pt_mgen[1]->Fill( (cms2.mus_p4()[muidxs_noptiso[1]].pt()-ptgen)/ptgen, weight );
		  hdilep_pt_mgen[DILEPTON_ALL]->Fill( (cms2.mus_p4()[muidxs_noptiso[1]].pt()-ptgen)/ptgen, weight );
		}
	  }
	  if( scale*pt < 20 ) { //loose leg is in the low pt range
		hdilep_lt_ll20_pt_eta[1]           ->Fill( scale*cms2.mus_p4()[muidxs_noptiso[0]].pt(), cms2.mus_p4()[muidxs_noptiso[0]].eta(), weight );
		hdilep_lt_ll20_pt_eta[DILEPTON_ALL]->Fill( scale*cms2.mus_p4()[muidxs_noptiso[0]].pt(), cms2.mus_p4()[muidxs_noptiso[0]].eta(), weight );
		double mass = (cms2.mus_p4()[muidxs_noptiso[0]] + cms2.mus_p4()[muidxs_noptiso[1]]).mass();
		hdilep_mass_ll20[1]->Fill( mass, weight );
		hdilep_mass_ll20[DILEPTON_ALL]->Fill( mass, weight );
	  }
	}
	else
	  cout << "ERROR IN Z ABCD SELECT\n\n";
  }

  //ABCD for Ws
  //NOTE: since the counting for zs is done for id only, and for ws for id+pt, can get event which satisfies both (one lep for z fails pt, other passes) But tmass, zmass, and met will help you
  cuts_t wabcd = 0;
  int lep_type = -1;
  if( nels_noiso == 1 && nmus_noiso == 0 ) {
	wabcd = pass_elcuts_noiso[0];
	lep_type = 0;
  }
  else if( nels_noiso == 0 && nmus_noiso == 1 ) {
	wabcd = pass_mucuts_noiso[0];
	lep_type = 1;
  }


  if( (w_cuts_nometiso & wabcd) == w_cuts_nometiso ) { // && gen_boson_mass > 65 && gen_boson_mass < 95
	//for plot with neutrino eta cut: find neutrino
	int nuidx = -1;
	//ABCD in tcmet, Iso
	if( lep_type == 0 ) {
	  //nels_nometiso++;
	  hlep_met_trckIso[0]->Fill( cms2.evt_tcmet(), cms2.els_tkIso()[elidxs_noiso[0]], weight );		
	  hlep_met_trckIso[2]->Fill( cms2.evt_tcmet(), cms2.els_tkIso()[elidxs_noiso[0]], weight );		
	  hlep_met_ecalIso[0]->Fill( cms2.evt_tcmet(), cms2.els_pat_ecalIso()[elidxs_noiso[0]], weight );	
	  hlep_met_ecalIso[2]->Fill( cms2.evt_tcmet(), cms2.els_pat_ecalIso()[elidxs_noiso[0]], weight );	
	  hlep_met_hcalIso[0]->Fill( cms2.evt_tcmet(), cms2.els_pat_hcalIso()[elidxs_noiso[0]], weight );	
	  hlep_met_hcalIso[2]->Fill( cms2.evt_tcmet(), cms2.els_pat_hcalIso()[elidxs_noiso[0]], weight );	
	  hlep_met_relIso[0] ->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(elidxs_noiso[0], true), weight );
	  hlep_met_relIso[2] ->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(elidxs_noiso[0], true), weight );
	  for( unsigned int i=0;i<cms2.genps_p4().size();i++ ) {
		if( abs( cms2.genps_id()[i] ) == 12 ) nuidx = i; //el nu only 
	  }
	  if( nuidx != -1 && fabs( cms2.genps_p4()[nuidx].eta() ) < 2.4 ) {
		hlep_accmet_relIso[0]->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(elidxs_noiso[0], true), weight );
		hlep_accmet_relIso[2]->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(elidxs_noiso[0], true), weight );
	  }
	  hlep_genpt_mch[0]->Fill( gen_lt_pt, weight);
	  hlep_genpt_mch[2]->Fill( gen_lt_pt, weight);
	}
	else if( lep_type == 1) {
	  //nmus_nometiso++;
	  hlep_met_trckIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_sumPt()[muidxs_noiso[0]], weight );		
	  hlep_met_trckIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_sumPt()[muidxs_noiso[0]], weight );		
	  hlep_met_ecalIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_emEt()[muidxs_noiso[0]], weight );	
	  hlep_met_ecalIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_emEt()[muidxs_noiso[0]], weight );	
	  hlep_met_hcalIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_hadEt()[muidxs_noiso[0]], weight );	
	  hlep_met_hcalIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_hadEt()[muidxs_noiso[0]], weight );	
	  hlep_met_relIso[1] ->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(muidxs_noiso[0]), weight );
	  hlep_met_relIso[2] ->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(muidxs_noiso[0]), weight );
	  for( unsigned int i=0;i<cms2.genps_p4().size();i++ ) {
		if( abs( cms2.genps_id()[i] ) == 14 ) nuidx = i;
	  }
	  if( nuidx != -1 && fabs( cms2.genps_p4()[nuidx].eta() ) < 2.4 ) {
		hlep_accmet_relIso[1]->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(muidxs_noiso[0]), weight );
		hlep_accmet_relIso[2]->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(muidxs_noiso[0]), weight );
	  }
	  hlep_genpt_mch[1]->Fill( gen_lt_pt, weight);
	  hlep_genpt_mch[2]->Fill( gen_lt_pt, weight);
	}
	else
	  cout << "ERROR IN SELECTING W ABCD\n"; //since cuts are checked, must be e or mu 
  }

  //Z--enforce exactly two leptons and SF requirement
  //remember, nels, nmus require NO MET, NOR MASS, NOR NJET
  if( (nels == 2 && nmus == 0) || (nmus == 2 && nels == 0) ) {

	const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);

	//const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i]);
	//DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE

	//cuts_t dilepcuts_passed = DilepSelect(); // does this hypothesis pass the required cuts?
	int idx1 = (elidxs[0] != -1 ? elidxs[0] : muidxs[0]);
	int idx2 = (elidxs[1] != -1 ? elidxs[1] : muidxs[1]);
	pass_nlep_cuts |= DilepSelect(myType, idx1, idx2); 

	double mass=0;
	if( elidxs[0] != -1 && elidxs[1] != -1 && myType == DILEPTON_EE ) {
	  mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
	}
	else if( muidxs[0] != -1 && muidxs[1] != -1 && myType == DILEPTON_MUMU ) {
	  mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
	}
	else {
	  cout << "BAD LEP INDXS IN Z SELECTION\n" << endl;
	  exit(1);
	}
	
	if( (pass_nlep_cuts & z_cuts_nomass) == z_cuts_nomass ) {
	  hdilep_mass[myType]->Fill(mass, weight);
	  hdilep_mass[DILEPTON_ALL]->Fill(mass, weight);
	}

	//don't check njet cut for njet histo
	int njet = (elidxs[0] != -1 ? njetels_20 : njetmus_20 );
	if( (pass_nlep_cuts & z_cuts_nojet20) == z_cuts_nojet20 ) {
	  hdilep_njet20[myType]->Fill(njet, weight);
	  hdilep_njet20[DILEPTON_ALL]->Fill(njet, weight);
	}

	if( (pass_nlep_cuts & z_cuts) == z_cuts )
	  ZEvent();
  }

  else if( (nels == 1 && nmus == 0) || (nmus == 1 && nels == 0) ) {
	unsigned int lep_type = (nels == 1 ? 0 : 1);

	if( (pass_lep_cut & w_cuts_nomet) == w_cuts_nomet ) {
	  hlep_tcmet[lep_type]->Fill(cms2.evt_tcmet(), weight);
	  hlep_tcmet[2]->Fill(cms2.evt_tcmet(), weight);
	  hlep_clmumet[lep_type]->Fill(cms2.evt_metMuonCorr(), weight);
	  hlep_clmumet[2]->Fill(cms2.evt_metMuonCorr(), weight);
	}

	//don't check njet cut for njet histo
	int njet = (elidxs[0] != -1 ? njetels_20 : njetmus_20 );
	if( (pass_lep_cut & w_cuts_nojet20) == w_cuts_nojet20 ) {
	  hlep_njet20[lep_type]->Fill(njet, weight);
	  hlep_njet20[2]->Fill(njet, weight);
	}

	//for WEvent, need to pass tmass, met
	if( (pass_lep_cut & w_cuts) == w_cuts ) //new, counting excludes iso, met, njet
	  WEvent();
  }

}
//end FillEventHistos

void ABCDLooper::WEvent() {
  
  double weight = ABCDLooper::Weight();

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

  if( lep_type == 0 ) {
	hlep_pt_f[lep_type]->Fill( cms2.els_p4()[elidxs[0]].pt(), weight );
	hlep_pt_f[2]->Fill( cms2.els_p4()[elidxs[0]].pt(), weight );

	//conversion plot has ALL cuts applied
	hlep_conv[lep_type]->Fill( conversionElectron(elidxs[0]), weight );
	hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );

	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[elidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);
  }
  else if (lep_type == 1) {
	hlep_pt_f[lep_type]->Fill( cms2.mus_p4()[muidxs[0]].pt(), weight );
	hlep_pt_f[2]->Fill( cms2.mus_p4()[muidxs[0]].pt(), weight );

	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);
  }

  scands_passing_[lep_type] += weight;
  scands_passing_w2_[lep_type] += weight * weight;
  scands_count_[lep_type]++;
  scands_passing_[2] += weight;
  scands_passing_w2_[2] += weight * weight;
  scands_count_[2]++;

}

void ABCDLooper::ZEvent ()
{
  double weight = ABCDLooper::Weight();
  const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);

  //if( (cuts_passed & cuts) == cuts ) { //all checking done already

  if( myType == DILEPTON_EE ) {
	hlep_conv[0]->Fill( conversionElectron(elidxs[0]), weight );
	hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );
	hlep_conv[0]->Fill( conversionElectron(elidxs[1]), weight );
	hlep_conv[2]->Fill( conversionElectron(elidxs[1]), weight );
  }

  //no met cut, so require all cuts for met plot
  hdilep_tcmet[myType]->Fill(cms2.evt_tcmet(), weight);
  hdilep_tcmet[DILEPTON_ALL]->Fill(cms2.evt_tcmet(), weight);
  hdilep_clmumet[myType]->Fill(cms2.evt_metMuonCorr(), weight);
  hdilep_clmumet[DILEPTON_ALL]->Fill(cms2.evt_metMuonCorr(), weight);

  //all other hists need to fill before checking cuts

  dcands_passing_[myType] += weight;
  dcands_passing_w2_[myType] += weight * weight;
  dcands_count_[myType]++;
  dcands_passing_[DILEPTON_ALL] += weight;
  dcands_passing_w2_[DILEPTON_ALL] += weight * weight;
  dcands_count_[DILEPTON_ALL]++;

  //}

}


void ABCDLooper::End ()
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

