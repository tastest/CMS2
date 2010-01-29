#ifndef WW_fitWjets_h
#define WW_fitWjets_h

class TPaveText;
class RooFitResult;
class RooAbsData;
class TH1F;

// best for asymmetrical distributions with long tails
// fake type:
//   0 - electron fakes
//   1 - muon fakes
// 
// pdf type:
//   0 - 1st degree polynomial (make lead to negative background estimates
//   1 - 0 degree polynomial (flat background)
//   2 - power law pdf
//   3 - non-parametric QCD based background pdf


TPaveText* estimate_background(int pdf_type, 
			       RooFitResult* result, 
			       float xmin, float xsig);
RooFitResult* fit_isolation(RooAbsData* control_sample, 
			    RooAbsData* analysis_sample,
			    int fake_type,
			    int pdf_type,
			    const char* title,
			    TH1F* bkg_control_sample = 0);

RooFitResult* fit_isolation(RooAbsData* full_sample, 
			    int fake_type, 
			    int pdf_type, 
			    const char* title,
			    TH1F* bkg_control_sample=0);

#endif
