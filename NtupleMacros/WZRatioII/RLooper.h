// -*- C++ -*-

#ifndef RLOOPER_H
#define RLOOPER_H

#include "Math/Point3D.h"
#include "Tools/LooperBase.h"
#include "Cuts.h"
#include "CORE/CMS2.h"

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> Point;

class RLooper : public LooperBase {

public:
  RLooper (Sample, cuts_t wcuts, const char *logfilename = 0, cuts_t zcuts = 0);
  virtual ~RLooper () { }

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

  //Z
  TH1F *hdilep_0_pt[4];
  TH1F *hdilep_1_pt[4];

  double transmass;
  int njets_20;
  int njets_30;
  int elidxs[2];
  int muidxs[2];
  int dil_njets_20;
  int dil_njets_30;

  cuts_t wcuts_;
  cuts_t zcuts_;

  //for muon counting
  double mu_tot;
  double mu_multi; //more than 1 mu, but 1 selected
  double mu_mulma; //more than 1, and makes z mass with selected

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

