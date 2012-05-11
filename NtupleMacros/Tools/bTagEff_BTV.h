#include <string.h>
double getMisTagRate(double jet_pt, double jet_eta, string algo);
double getMisTagRate_Err(double jet_pt, double jet_eta, string algo);
double getMisTagSF(double jet_pt, double jet_eta, string algo);
double getMisTagSF_Err(double jet_pt, double jet_eta, string algo);
float  getBTagSF(const string algo, const float discriminator);
float  getBTagSF_Err(const string algo);
float  getBTagEff(const string algo, const float x, const bool isgen);
float  getBTagEff_Err(const string algo, const float x, const bool isgen);
float  getCTagEff(const string algo, const float x, const bool isgen);
