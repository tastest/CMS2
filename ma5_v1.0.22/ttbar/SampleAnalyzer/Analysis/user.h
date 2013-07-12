#ifndef analysis_user_h
#define analysis_user_h

#include "Core/AnalysisBase.h"

class user : public AnalysisBase
{
  INIT_ANALYSIS(user,"ttbar")

 public:
  virtual void Initialize();
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual void Execute(const SampleFormat& sample, const EventFormat& event);

 private:
  TH1F* myHisto;
  TH1F* myHisto_neg;
  TH1F* myHisto2;
  TH1F* myHisto2_neg;
  TH1F* myHisto3;
  TH1F* myHisto4;
};

#endif