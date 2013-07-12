#ifndef analysis_user_h
#define analysis_user_h

#include "Core/AnalysisBase.h"

class user : public AnalysisBase
{
  INIT_ANALYSIS(user,"dPhi")

 public:
  virtual void Initialize();
  virtual void Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files);
  virtual void Execute(const SampleFormat& sample, const EventFormat& event);

 private:
  TH1F* myHisto;
};

#endif