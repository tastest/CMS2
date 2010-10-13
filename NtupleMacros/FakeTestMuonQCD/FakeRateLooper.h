// -*- C++ -*-

#ifndef FAKERATELOOPER_H
#define FAKERATELOOPER_H

#include "Looper.h"
#include "TH3F.h"

class FakeRateLooper : public Looper {
public:
  FakeRateLooper (Sample s, cuts_t cuts, const char *fname = 0);
  virtual double	CandsPassingSystHi () const { return cands_passing_syst_hi_; }
  virtual double	CandsPassingSystLo () const { return cands_passing_syst_lo_; }
  virtual double	CandsPassingEventWeightOnly () const { return cands_passing_event_weight_only_; }
  virtual double	RMSEventWeightOnly () const { return sqrt(cands_passing_event_weight_only_w2_); }
  virtual double	FakeSyst () const;
  virtual void	End		();
  bool fillErrorInPrediction(TH1F* prediction,
			     TH3F* predictionError,
			     bool addStatisticalError = false);
  bool extractBins(TAxis *axis, unsigned int nBins, float* bins);
protected:
  virtual void	        BookHistos 	();
  virtual cuts_t	EventSelect 	();
  virtual void	        FillEventHistos ();
  virtual double	Weight		(int idx);
  virtual double	Weight		(int idx, int n_sig_syst);

protected:
  double cands_passing_event_weight_only_;
  double cands_passing_event_weight_only_w2_;
  double cands_passing_syst_hi_;
  double cands_passing_syst_lo_;
  TH2F	 *fake_syst;
  TH3F	 *hmuPt3D_;
  TH3F	 *hmuEta3D_;

};

#endif
