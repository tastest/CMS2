#include "TMath.h"

#include "TMinuit.h"

double probMultinom(double n0, double n1, double n2,
		    double m0, double m1, double m2){
  double p0 = m0/(m0+m1+m2); if (p0==0){ p0=1; n0=0;}
  double p1 = m1/(m0+m1+m2); if (p1==0){ p1=1; n1=0;}
  double p2 = m2/(m0+m1+m2); if (p2==0){ p2=1; n2=0;}
  double prob = TMath::Gamma(n0+n1+n2+1)/
    TMath::Gamma(n0+1)/TMath::Gamma(n1+1)/TMath::Gamma(n2+1);
  prob *= TMath::Power(p0, n0)*TMath::Power(p1, n1)*TMath::Power(p2, n2);
  return prob;
}
double probBinom(double n0, double n1,
		 double m0, double m1){
  double p0 = m0/(m0+m1);
  double p1 = m1/(m0+m1);
  double prob = TMath::Gamma(n0+n1+1)/
    TMath::Gamma(n0+1)/TMath::Gamma(n1+1);
  prob *= TMath::Power(p0, n0)*TMath::Power(p1, n1);
  return prob;
}


double sigNPL(double nS, double nB, double f2S, double f2B){
  double f2B_exp = 0.3;
  double f2B_err = 0.1;
  double nB_exp = 2.7032;
  double nB_err = 1.37978;

  double obs1 = 21;
  double obs2 = 30;
  double obs = obs1+obs2;

  double n1 = (1.-f2S)*nS + (1.-f2B)*nB;
  double n2 = f2S*nS + f2B*nB;

  double nlP = -1.*(
		   log (TMath::Poisson(nS+nB,obs)) - log (TMath::Poisson(obs, obs))
		   + log (probBinom(n1,n2,obs1,obs2)) - log (probBinom(obs1,obs2,obs1,obs2))
		   + log (TMath::Gaus(f2B, f2B_exp, f2B_err)) 
		   - log( 1. - 0.5*TMath::Prob(TMath::Power(f2B_exp/f2B_err,2),1))
		   + log (TMath::Gaus(nB, nB_exp, nB_err)) 
		   - log( 1. - 0.5*TMath::Prob(TMath::Power(nB_exp/nB_err,2),1))
		   );
  return nlP;
}


void fcn(int& npar, double* gin, double& f, double* par, int iflag){
  f= sigNPL(par[0], par[1], par[2], par[3]);
}


void minimize(){
  gMinuit = new TMinuit(4);
  gMinuit->SetFCN(fcn);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 0.5;

  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  static Double_t vstart[4]  = {40,3,0.5,0.5};
  Double_t step[4] = {0.001, 0.001, 0.0001,0.0001};
  gMinuit->mnparm(0, "nS", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "nB", vstart[1], step[1], 0,0,ierflg);
  gMinuit->mnparm(2, "fS", vstart[2], step[2], 0,0,ierflg);
  gMinuit->mnparm(3, "fB", vstart[3], step[3], 0,0,ierflg);
  arglist[0] = 500;
  arglist[1] = 0.01;
  gMinuit->mnexcm("MIGRAD", arglist ,4,ierflg);
  arglist[0] = 2;
  arglist[1] = 4;
  gMinuit->mnexcm("FIX", arglist ,4,ierflg);
  arglist[0] = 500;
  arglist[1] = 0.01;
  gMinuit->mnexcm("MIGRAD", arglist ,4,ierflg);
//   gMinuit->mnexcm("FIX", arglist ,4,ierflg);
//   gMinuit->mnexcm("FIX", arglist ,4,ierflg);
  gMinuit->mnexcm("MINOS", arglist ,4,ierflg);

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(4,amin);
}
