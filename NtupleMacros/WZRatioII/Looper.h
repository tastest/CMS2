// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "Cuts.h"

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

class Looper : public LooperBase {

public:
  Looper (Sample, cuts_t cuts, const char *logfilename = 0);
  virtual ~Looper () { }

protected:
  // this is where we book our histograms
  virtual void	BookHistos ();

  virtual cuts_t DilepSelect(); //(int i_hyp); //now works off idxs
  virtual cuts_t LepSelect (int lep_type, int i);

  void WEvent();
  void ZEvent();

  virtual void	FillEventHistos ();
  virtual void	End		();

  // do stuff with histogram
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
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  TH1F	*h1_lep_pt_[3];
  TH1F	*h1_lep_met_[3];
  TH1F  *h1_lep_met_dphi_[3];
  TH1F  *h1_lep_tkIso_[3];

  TH1F	*h1_dilep_0_pt_[4];
  TH1F	*h1_dilep_1_pt_[4];
  TH1F 	*h1_dilep_mass_[4];
  TH1F	*h1_dilep_met_[4];
  TH1F	*h1_dilep_nhyp_;

  int elidxs[2];
  int muidxs[2];

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

