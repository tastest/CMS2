#include "RooFit.h"
#include "Riostream.h"
#include "RooATGCPdf.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "TH1.h"
#include "TVirtualFitter.h"

ClassImp(RooATGCPdf);

RooATGCPdf::RooATGCPdf(const char *name, const char *title, 
		       RooRealVar& _obs, RooRealVar& _nSM, RooRealVar& _x, RooRealVar& _y):
  RooAbsPdf(name,title),
  obs("obs","observable", this, _obs),
  nSM("nSM","Expected number of events for Standard Model", this, _nSM),
  x("x","first parameter", this, _x),
  y("y","second parameter", this, _y),
  hist(0),
  configured(false)
{
}

RooATGCPdf::RooATGCPdf(const RooATGCPdf& other, const char* name) : 
  RooAbsPdf(other,name), params(other.params), expectedYield(other.expectedYield),
  obs("obs",this,other.obs), nSM("nSM",this,other.nSM), x("x",this,other.x), y("y",this,other.y),
  hist(other.hist), configured(other.configured)
{
}

void RooATGCPdf::addPoint(Measurement m, const TH1* h)
{
  if (!hist){
    hist = dynamic_cast<TH1*>(h->Clone("hist"));
    assert(hist);
    hist->SetDirectory(0);
  }
  else
    assert(h->GetNbinsX()==hist->GetNbinsX());
  input.push_back(std::pair<Measurement,const TH1*>(m,h));
}

bool RooATGCPdf::build(){
  assert(input.size()>5);
  assert(hist);
  Int_t size = hist->GetNbinsX();
  params.resize(size+2);
  // skip under and overflow.
  for ( Int_t i=1; i<=size; ++i ){
    std::vector<Measurement> mb;
    for ( std::vector<std::pair<Measurement,const TH1*> >::const_iterator itr=input.begin();
	    itr!=input.end(); ++itr )
      {
	Measurement m = itr->first;
	double scale = itr->second->GetBinContent(i)/itr->second->Integral(1,size);
	// ignore histogram errors for now
	m.xs *= scale;
	m.xserr *= scale;
	mb.push_back(m);
      }
    params.at(i) = Parabola( &mb );
  }
  // set expected yeild
  std::vector<Measurement> mb;
  for ( std::vector<std::pair<Measurement,const TH1*> >::const_iterator itr=input.begin();
	itr!=input.end(); ++itr )
    {
      mb.push_back(itr->first);
    }
  expectedYield = Parabola( &mb );
  configured = true;
  // here should be a check if parameterization crosses zero
  return true;
}

Double_t RooATGCPdf::evaluate() const
{
  assert(hist);
  assert(configured);
  Int_t bin = hist->FindBin(obs);
  return params.at(bin).value(x,y);
}
 
Int_t RooATGCPdf::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const 
{
  // Only analytical integrals over the full range are defined
  if (rangeName!=0) {
    return 0 ;
  }
  if (matchArgs(allVars,analVars,obs)) return 1 ;
  return 0 ;
}

Double_t RooATGCPdf::analyticalIntegral(Int_t code, const char* /*rangeName*/) const
{
  assert(code==1);
  assert(hist);
  Int_t size = hist->GetNbinsX();
  double integral(0);
  for (Int_t i = 1; i<= size; ++i )
    integral+=params.at(i).value(x,y)*hist->GetBinWidth(i);
  return integral;
//   static const Double_t root2 = sqrt(2.) ;
//   static const Double_t rootPiBy2 = sqrt(atan2(0.0,-1.0)/2.0);
//   Double_t xscale = root2*sigma;
//   Double_t ret = 0;
//   if(code==1){  
//     ret = rootPiBy2*sigma*(RooMath::erf((x.max(rangeName)-mean)/xscale)-RooMath::erf((x.min(rangeName)-mean)/xscale));
//     //cout << "Int_gauss_dx(mean=" << mean << ",sigma=" << sigma << ", xmin=" << x.min(rangeName) << ", xmax=" << x.max(rangeName) << ")=" << ret << endl ;
//   } else if(code==2) {
//     ret = rootPiBy2*sigma*(RooMath::erf((mean.max(rangeName)-mean)/xscale)-RooMath::erf((mean.min(rangeName)-mean)/xscale));
//   } else{
//     cout << "Error in RooGaussian::analyticalIntegral" << endl;
//   }
//   return ret ;
}

Parabola::Parabola( std::vector<Measurement>* _m )
{
  m = _m;
  TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
  TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 6);
  fitter->SetFCN(myfcn);
  
  fitter->SetParameter(0, "A",   0, 0.1, 0,0);
  fitter->SetParameter(1, "B",   0, 0.1, 0,0);
  fitter->SetParameter(2, "C",   0, 0.1, 0,0);
  fitter->SetParameter(3, "D",   0, 0.1, 0,0);
  fitter->SetParameter(4, "E",   0, 0.1, 0,0);
  fitter->SetParameter(5, "F",   0, 0.1, 0,0);
  
  Double_t arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);
  A = fitter->GetParameter(0);
  B = fitter->GetParameter(1);
  C = fitter->GetParameter(2);
  D = fitter->GetParameter(3);
  E = fitter->GetParameter(4);
  F = fitter->GetParameter(5);
  m = 0;
}

void 
Parabola::myfcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t)
{
  f = 0;
  assert(m);
  for (unsigned int i=0; i<m->size(); i++) {
    Double_t v = par[0]+(*m)[i].x*par[1]+(*m)[i].y*par[2]+(*m)[i].x*(*m)[i].y*par[3]+
      (*m)[i].x*(*m)[i].x*par[4]+(*m)[i].y*(*m)[i].y*par[5];
    Double_t d = (v-(*m)[i].xs)/(*m)[i].xserr;
    f += d*d;
  }
}

std::vector<Measurement>* Parabola::m = 0;

Double_t RooATGCPdf::expectedEvents(const RooArgSet* /*nset*/) const
{
  assert(configured);
  return expectedYield.value(x,y)/expectedYield.value(0,0)*nSM;
}
