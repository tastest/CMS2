#include "doAnalysis.h"
#include "fitWjets.h"

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

// #pragma link C++ function ScanChain(TChain*, enum Sample, double, bool);
#pragma link C++ enum Sample;
#pragma link C++ function estimate_background(int, RooFitResult*, float, float);
#pragma link C++ function fit_isolation(RooAbsData*, RooAbsData*, int, int, const char*, TH1F*);
// #pragma link C++ function ProcessSample(std::string, enum Sample, double, RooDataSet*, Color_t, bool);
// #pragma link C++ function ProcessSample(std::vector<std::string>, enum Sample, double, RooDataSet*, Color_t, bool);

#endif
