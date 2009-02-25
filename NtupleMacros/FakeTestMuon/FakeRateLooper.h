// -*- C++ -*-

#ifndef FAKERATELOOPER_H
#define FAKERATELOOPER_H

#include "Looper.h"
#include "TH3F.h"

// background estimate for W+jets from fake rates (only electrons in
// emu for now)
class FakeRateLooper : public Looper {
public:
  FakeRateLooper (Sample s, cuts_t cuts, const char *fname = 0);
  virtual double	CandsPassingSystHi (enum DileptonHypType i) const { return cands_passing_syst_hi[i]; }
  virtual double	CandsPassingSystLo (enum DileptonHypType i) const { return cands_passing_syst_lo[i]; }
  virtual double	CandsPassingEventWeightOnly (enum DileptonHypType i) const { return cands_passing_event_weight_only_[i]; }
  virtual double	RMSEventWeightOnly (enum DileptonHypType i) const { return sqrt(cands_passing_event_weight_only_w2_[i]); }
  virtual double	FakeSyst (enum DileptonHypType i) const;
  virtual void	End		();
  bool fillErrorInPrediction(TH1F* prediction,
			     TH3F* predictionError,
			     bool addStatisticalError = true);
  bool extractBins(TAxis *axis, unsigned int nBins, float* bins);
protected:
  virtual void	BookHistos 	();
  virtual cuts_t	DilepSelect 	(int idx);
  virtual void	FillDilepHistos (int idx);
  virtual double	Weight		(int idx);
  virtual double	Weight		(int idx, int n_sig_syst);

protected:
  double cands_passing_event_weight_only_[4];
  double cands_passing_event_weight_only_w2_[4];
  double cands_passing_syst_hi[4];
  double cands_passing_syst_lo[4];
  TH2F	 *fake_syst;
  TH3F	 *hnJet3D_[4];
  TH3F	 *helPt3D_[4];
  TH3F	 *helEta3D_[4];
  TH3F	 *hmet3D_[4];

};

#endif
