#include "Math/Random.h"
#include "Math/GSLRndmEngines.h"
#include "TMath.h"
#include "TH1F.h"
#include <iostream>
void estimate_signal_significance(){
  using namespace std;
  float expectedNumberOfEvents = 47.9;
  float expectedBackgroundYield = 10.4;
  float expectedBackgroundYieldUncertainty = 4.6;
  long long int   nBackgroundTests = 1e8;
  long long int   nTests = 1e5;
  int   nBins = int(expectedNumberOfEvents*10);
  const double nSigma = 5.0;
  const double threshold = TMath::Erfc(nSigma/sqrt(2));
  
  // setup generator
  // http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
  // gsl_rng_mt19937
  ROOT::Math::Random<ROOT::Math::GSLRngMT> gen;

  TH1D* h = new TH1D("h","h",nBins,0,nBins);
  int nMoreThanNsigmaDeviationsIsNeededTotal(0);
  const int nCycles = 1;
  for ( int kk=0; kk < nCycles; ++kk ){
    h->Reset();
    cout << "Generating background distribution..." <<endl;
    for (long long int i=0; i<nBackgroundTests; ++i){
      float a = gen.Gaus(expectedBackgroundYield,expectedBackgroundYieldUncertainty); 
      if (a < 0) a = 0; // background estimate cannot be less than zero 
      h->Fill(gen.Poisson(a));
    }
    cout << "done" <<endl;
    cout << "Fluctuating total number of events." <<endl;
    long long int nMoreThanNsigmaDeviationsIsNeeded(0);
    double integral = h->Integral(0,nBins+1);
    for (long long int i=0; i<nTests; ++i){
      int nTotal = gen.Poisson(expectedNumberOfEvents);
      if(nTotal<nBins && h->Integral(nTotal+1,nBins+1)/integral < threshold ) ++nMoreThanNsigmaDeviationsIsNeeded;
    } 
    cout << "Probability to get " << nSigma << " sigma: " << nMoreThanNsigmaDeviationsIsNeeded/double(nTests) << endl;
    nMoreThanNsigmaDeviationsIsNeededTotal+=nMoreThanNsigmaDeviationsIsNeeded;
  }
  cout << "Total probability to get " << nSigma << " sigma: " << nMoreThanNsigmaDeviationsIsNeededTotal/double(nTests*nCycles) << endl;
}
