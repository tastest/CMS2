// -*- C++ -*-

#ifndef ABCDLOOPER_H
#define ABCDLOOPER_H

#include "Tools/LooperBase.h"
#include "Cuts.h"
#include "CORE/CMS2.h"

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

class ABCDLooper : public LooperBase {

public:
  ABCDLooper (Sample, cuts_t cuts, const char *logfilename = 0);
  virtual ~ABCDLooper () { }

protected:
  // this is where we book our histograms
  virtual void	BookHistos ();

  virtual cuts_t DilepSelect(const enum DileptonHypType myType, int idx1, int idx2);
  virtual cuts_t LepSelect (int lep_type, int i);

  void WEvent();
  void ZEvent();

  virtual void	FillEventHistos ();
  virtual void	End		();

  // Weight for WZ analysis is pb
  double Weight() { return cms2.evt_scale1fb() * sample_.kFactor / 1000; }

  // do stuff with histogram
  virtual void NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max);

  void FormatHist(TH1* hist);



public:
  // these functions are called by the table-printing code
  virtual double	DCandsPassing (enum DileptonHypType i) const 		{ return dcands_passing_[i]; }
  virtual int      	DCandsCount (enum DileptonHypType i) const 			{ return dcands_count_[i]; }
  virtual double	DRMS (enum DileptonHypType i) const 				{ return sqrt(dcands_passing_w2_[i]); }
  virtual double	SCandsPassing (enum DileptonHypType i) const 		{ return scands_passing_[i]; }
  virtual int      	SCandsCount (enum DileptonHypType i) const 			{ return scands_count_[i]; }
  virtual double	SRMS (enum DileptonHypType i) const 				{ return sqrt(scands_passing_w2_[i]); }

protected:

  //W
  TH1F *hlep_pt[3];
  TH1F *hlep_genpt[3];
  TH1F *hlep_genpt_mch[3];
  TH1F *hlep_pt_f[3];
  TH1F *hlep_mass[3];
  TH1F *hlep_tcmet[3];
  TH1F *hlep_genmet[3];
  TH1F *hlep_accgenmet[3];
  TH1F *hlep_clmumet[3];
  TH1F *hlep_met_dphi[3];
  TH1F *hlep_trckIso[3];
  TH1F *hlep_ecalIso[3];
  TH1F *hlep_hcalIso[3];
  TH1F *hlep_relIso[3];
  TH1F *hlep_nlep[3];
  TH1F *hlep_nlep_nod0iso[3];
  TH1F *hlep_nlep_nometiso[3];
  TH1F *hlep_njet20[3];
  TH1F *hlep_njet30[3];
  TH1F *hlep_conv[3];
  TH1F *hlep_d0[3];
  TH1F *hlep_d0Sig[3];
  //TH2Fs for ABCD
  TH2F* hlep_d0_trckIso[3]; 
  TH2F* hlep_d0_ecalIso[3]; 
  TH2F* hlep_d0_hcalIso[3]; 
  TH2F* hlep_d0_relIso[3];  
  TH2F* hlep_met_trckIso[3];
  TH2F* hlep_met_ecalIso[3];
  TH2F* hlep_met_hcalIso[3];
  TH2F* hlep_met_relIso[3]; 
  TH2F* hlep_accmet_relIso[3]; //met with eta cut on neutrino

  //Z
  TH1F *hdilep_0_pt[4];
  TH1F *hdilep_1_pt[4];
  TH1F *hdilep_pt[4];
  TH1F *hdilep_pt_mgen[4];
  TH1F *hdilep_genpt[4];
  TH1F *hdilep_0_genpt[4];
  TH1F *hdilep_1_genpt[4];
  TH1F *hdilep_0_genpt_mch[4];
  TH1F *hdilep_1_genpt_mch[4];
  TH1F *hdilep_mass[4];
  TH1F *hdilep_mass_ll20[4];
  TH1F *hdilep_tmass[4];
  TH1F *hdilep_tcmet[4];
  TH1F *hdilep_clmumet[4];
  TH1F *hdilep_njet20[4];
  TH1F *hdilep_njet30[4];
  TH1F *hdilep_reliso_ll[4];  

  //Z TH2Fs for supplemental ABCD 
  TH2F* hdilep_ll_pt_eta[4];
  TH2F* hdilep_lt_pt_eta[4];
  TH2F* hdilep_lt_ll20_pt_eta[4];
  //TH2F* hdilep_lliso_pt_eta[4];
  //TH2F* hdilep_d0_trckIso[4]; 
  //TH2F* hdilep_d0_ecalIso[4]; 
  //TH2F* hdilep_d0_hcalIso[4]; 
  //TH2F* hdilep_d0_relIso[4];  
  TH2F* hdilep_lpt_trckIso[4]; //call it lpt because it's actually lepton pt that we call met
  TH2F* hdilep_lpt_ecalIso[4];
  TH2F* hdilep_lpt_hcalIso[4];
  TH2F* hdilep_lpt_relIso[4]; 
  TH2F* hdilep_lepmet_trckIso[4];
  TH2F* hdilep_lepmet_ecalIso[4];
  TH2F* hdilep_lepmet_hcalIso[4];
  TH2F* hdilep_lepmet_relIso[4]; 
  //same as above 8, but scaled by mw/mz
  TH2F* hdilep_lpt_scl_trckIso[4];
  TH2F* hdilep_lpt_scl_ecalIso[4];
  TH2F* hdilep_lpt_scl_hcalIso[4];
  TH2F* hdilep_lpt_scl_relIso[4]; 
  //TH2F* hdilep_lpt_rscl_relIso[4]; 
  TH2F* hdilep_lepmet_scl_trckIso[4];
  TH2F* hdilep_lepmet_scl_ecalIso[4];
  TH2F* hdilep_lepmet_scl_hcalIso[4];
  TH2F* hdilep_lepmet_scl_relIso[4]; 
  TH2F* hdilep_lepmet_rscl_relIso[4]; 
  TH2F* hdilep_lepmet_scl_trth_relIso[4];
  TH2F* hdilep_lepmet_rscl_trth_relIso[4];
  TH2F* hdilep_lepmet_scl_tmas_relIso[4];
  TH2F* hdilep_lepmet_rscl_tmas_relIso[4];
  TH2F* hdilep_lepmet_scl_tmast_relIso[4];
  TH2F* hdilep_lepmet_rscl_tmast_relIso[4];
  TH2F* hdilep_lepmet_scl_tmastmes_relIso[4];
  TH2F* hdilep_lepmet_rscl_tmastmes_relIso[4];


  /*
  //TH1F *hdilep_nhyp;
  //TH1F *hdilep_nlep;
  TH1F	*h1_lep_Highpt_[3];
  TH1F	*h1_lep_HighptMet_[3];
  TH1F	*h1_lep_HighptRelIso_[3];
  TH1F	*h1_lep_HighptRelIsoPtLg20_[3];

  TH1F	*h1_lep_Lowpt_[3];
  TH1F	*h1_lep_LowptMet_[3];
  TH1F	*h1_lep_LowptRelIso_[3];
  TH1F	*h1_lep_LowptRelIsoPtLg20_[3];
  TH1F	*h1_lep_LowptNLepGt10Lt20_[3];
  TH1F	*h1_lep_LowptNLepGt20_[3];
  TH1F	*h1_lep_LowptNLepGt20tightIDIso0_1_[3];
  TH1F	*h1_lep_LowptNLepGt20looseIDIso0_1_[3];
  TH1F	*h1_lep_LowptNLepGt20vlooseIDIso0_1_[3];
  TH1F	*h1_lep_LowptNLepGt20NOIDIso0_1_[3];

  TH1F* h1_lep_LowpthOverE_[3];
  TH1F* h1_lep_lowpteOverPIn_[3];
//            TH1F*      h1_lep_lowpteSeedOverPOut_[3];
//            TH1F*      h1_lep_lowpteSeedOverPIn_[3];
  TH1F* h1_lep_lowptfBrem_[3];
  TH1F* h1_lep_lowptdEtaIn_[3];
  TH1F* h1_lep_lowptdEtaOut_[3];
  TH1F* h1_lep_lowptdPhiIn_[3];
  TH1F* h1_lep_lowptdPhiInPhiOut_[3];
  TH1F* h1_lep_lowptdPhiOut_[3];

  TH1F* h1_lep_lowptsigmaPhiPhi_[3];
  TH1F* h1_lep_lowptsigmaIPhiIPhi_[3];
  TH1F* h1_lep_lowptsigmaEtaEta_[3];
  TH1F* h1_lep_lowptsigmaIEtaIEta_[3];

  TH1F* h1_lep_lowptegamma_robustLooseId_[3];
  TH1F* h1_lep_lowptegamma_robustTightId_[3];
  TH1F* h1_lep_lowptegamma_looseId_[3];
  TH1F* h1_lep_lowptegamma_tightId_[3];
  TH1F* h1_lep_lowptegamma_robustHighEnergy_[3];
  */
  //int njets;//need to do differently for W,Z anyway
  double transmass;
  int njets_20;
  int njets_30;
  int elidxs[2];
  int muidxs[2];
  int dil_njets_20;
  int dil_njets_30;

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  //these are for dilepton (Z)
  double		dcands_passing_[4];
  double		dcands_passing_w2_[4];
  unsigned int	dcands_count_[4];
  //these are for single lepton (W)
  double		scands_passing_[3];
  double		scands_passing_w2_[3];
  unsigned int	scands_count_[3];
};
#endif

