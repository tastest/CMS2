#ifndef ROO_ATGC_PDF
#define ROO_ATGC_PDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooSetProxy.h"
#include "RooAICRegistry.h"
#include <vector>

class RooRealVar;
class TH1;

// x-section is parametrised as:
// A+B*x+C*y+D*x*y+E*x^2+F*y^2
struct Measurement{
  double x;
  double y;
  double xs;
  double xserr;
  Measurement(double _x, double _y, double _xs, double _xserr):
    x(_x), y(_y), xs(_xs), xserr(_xserr){}
  Measurement():x(0), y(0), xs(0), xserr(0){}
};

struct Parabola{
  double A,B,C,D,E,F;
  static std::vector<Measurement>* m;
  Parabola( std::vector<Measurement>* _m );
  Parabola():A(0),B(0),C(0),D(0),E(0),F(0){}
  static void myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t);
  double value(double x, double y) const
  {
    return A+B*x+C*y+D*x*y+E*x*x+F*y*y;
  }
};

class RooATGCPdf : public RooAbsPdf {
public:
  RooATGCPdf():hist(0){} 
  RooATGCPdf(const char *name, const char *title, 
	     RooRealVar& obs, RooRealVar& nSM, RooRealVar& x, RooRealVar& y);
  // obs - observable (leading lepton pt)
  // nSM - expected number of events for Standard Model
  // x - first parameter
  // y - second parameter

  void addPoint(Measurement, const TH1*);
  bool build();
  bool isConfigured() const {return configured;}
  RooATGCPdf(const RooATGCPdf& other, const char* name=0);
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const ;
  virtual TObject* clone(const char* newname) const { return new RooATGCPdf(*this,newname); }
  virtual ExtendMode extendMode() const { return CanBeExtended ; }
  // virtual ExtendMode extendMode() const { return CanNotBeExtended ; }
  virtual Double_t expectedEvents(const RooArgSet* nset) const ;
  virtual Double_t expectedEvents(const RooArgSet& nset) const { return expectedEvents(&nset) ; }
protected:
  std::vector<Parabola> params;
  Parabola expectedYield;
  std::vector<std::pair<Measurement,const TH1*> > input; 
  Double_t evaluate() const;
  RooRealProxy obs;
  RooRealProxy nSM;
  RooRealProxy x;
  RooRealProxy y;
  TH1* hist;
  bool configured;
  ClassDef(RooATGCPdf,1) 
};

#endif
