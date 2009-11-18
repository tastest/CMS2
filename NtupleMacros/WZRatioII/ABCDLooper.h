// -*- C++ -*-

#ifndef ABCDLOOPER_H
#define ABCDLOOPER_H

#include "Math/Point3D.h"
#include "Tools/LooperBase.h"
#include "Cuts.h"
#include "CORE/CMS2.h"

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point;

class ABCDLooper : public LooperBase {

public:
  ABCDLooper (Sample, cuts_t wcuts, const char *logfilename = 0, cuts_t zcuts = 0);
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

  
  double getCorrd0(const Point& myBeamSpot, const LorentzVector &p, const LorentzVector &v) const { 
	return ( -(v.x() - myBeamSpot.x())*p.y() + (v.y() - myBeamSpot.y())*p.x() ) / p.pt(); 
  }


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

  virtual bool isSsignal() const { return isssignal_; } //single lep signal
  virtual bool isDsignal() const { return isdsignal_; } //single dilep signal

  virtual void setSsignal(bool val) { isssignal_ = val; } //single lep signal
  virtual void setDsignal(bool val) { isdsignal_ = val; } //single dilep signal

protected:

  //W
  TH1F *hlep_pt[3];
  TH1F *hlep_genpt[3];
  TH1F *hlep_genpt_mch[3];
  TH1F *hlep_pt_f[3];
  TH1F *hlep_tmass[3];
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
  TH1F *hlep_d0err[3];
  TH1F *hlep_d0Sig[3];
  TH1F *hlep_d0Sigtest[3];
  TH1F *hlep_vtxd0[3];
  //TH1F *hlep_vtxd0err[3];
  TH1F *hlep_vtxd0Sig[3];
  //TH2Fs for ABCD
  TH2F* hlep_d0_d0err[3];  
  TH2F* hlep_d0_trckIso[3]; 
  TH2F* hlep_d0_ecalIso[3]; 
  TH2F* hlep_d0_hcalIso[3]; 
  TH2F* hlep_d0_relIso[3];  
  TH2F* hlep_d0sig_trckIso[3]; 
  TH2F* hlep_d0sig_ecalIso[3]; 
  TH2F* hlep_d0sig_hcalIso[3]; 
  TH2F* hlep_d0sig_relIso[3];  
  TH2F* hlep_vtxd0_d0err[3];  
  TH2F* hlep_vtxd0_trckIso[3]; 
  TH2F* hlep_vtxd0_ecalIso[3]; 
  TH2F* hlep_vtxd0_hcalIso[3]; 
  TH2F* hlep_vtxd0_relIso[3];  
  TH2F* hlep_vtxd0sig_trckIso[3]; 
  TH2F* hlep_vtxd0sig_ecalIso[3]; 
  TH2F* hlep_vtxd0sig_hcalIso[3]; 
  TH2F* hlep_vtxd0sig_relIso[3];  
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
  TH1F *hdilep_reliso_lt[4];  
  TH1F *hdilep_reliso_ll[4];  
  TH1F *hdilep_nopt_reliso_ll[4];  

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


  //TH1F *hdilep_nhyp;
  //TH1F *hdilep_nlep;
  double transmass;
  int njets_20;
  int njets_30;
  int elidxs[2];
  int muidxs[2];
  int dil_njets_20;
  int dil_njets_30;

  cuts_t wcuts_;
  cuts_t zcuts_;

protected:
  bool isssignal_; //single lep signal  
  bool isdsignal_; //single dilep signal

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

