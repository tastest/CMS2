#ifndef BTAG_PAYLOAD_B_H
#define BTAG_PAYLOAD_B_H

#include <vector>

namespace beff
{

  static unsigned int bins;
  static std::vector <float> ptmin;
  static std::vector <float> ptmax;
  static std::vector <float> etamin;
  static std::vector <float> etamax;

  void setpt();
  void seteta();
  void setbins();
  void setCSVL_SFb_error();
  void setCSVM_SFb_error();
  void setCSVT_SFb_error();

  // Tagger: CSVL within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVL_SFb(float x);
  static std::vector <float> CSVL_SFb_error;

  // Tagger: CSVM within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVM_SFb(float x);
  static std::vector <float> CSVM_SFb_error;

  // Tagger: CSVT within 20 < pt < 800 GeV, abs(eta) < 2.4, x = pt
  float CSVT_SFb(float x);
  static std::vector <float> CSVT_SFb_error;
}			   
#endif
