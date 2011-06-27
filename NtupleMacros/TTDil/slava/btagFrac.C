#include "TMath.h"

#include "TMinuit.h"

#include <iostream>

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

  //the above is for TCHEL
  bool isSSVHEM = true;
  if(isSSVHEM){
    f2B_exp = 0.51/2.29;//==0.22 //using MC bgds only here
    f2B_err = 0.1; //keep the same as in TCHEL: this time it's ~50% uncty
    nB_exp = 2.29;
    nB_err = nB_exp*1.37978/2.7032; //==1.2; leave fractionally the same uncty as in TCHEL
    
    obs1 = 20;
    obs2 = 15;
  }

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

double sigNPL3Bins(double nS, double sfB){

  double n0S_exp = 11.25;
  double n1S_exp = 26.02;
  double n2S_exp = 15.48;
  double nS_exp = n0S_exp + n1S_exp + n2S_exp;
  double f0S = 0;
  double f1S = 0;
  double nB_exp = 6.7;
  double nB_err = 2.0;

  double n0B_expMC = 4.55;
  double n1B_expMC = 1.78;
  double n2B_expMC = 0.51;
  double nB_expMC = n0B_expMC + n1B_expMC + n2B_expMC;
  
  double f0B_exp = n0B_expMC/nB_expMC;
  double f1B_exp = n1B_expMC/nB_expMC;
  double f2B_exp = 1. - f0B_exp - f1B_exp;

  double f0B_err = nB_err/nB_exp;
  double f1B_err = f1B_exp*0.5;
  double f2B_err = f2B_exp*0.5;
  
  double obs0 = 25;
  double obs1 = 20;
  double obs2 = 15;

  double obs = obs0 + obs1 + obs2;

  double nB = nB_exp;
  double f0B = f0B_exp;
  double f1B = f1B_exp;

  double f2S = 1. - f0S - f1S;
  double f2B = 1. - f0B - f1B;
  double n0 = nS*f0S + nB*f0B;
  double n1 = nS*f1S + nB*f1B;
  double n2 = nS*f2S + nB*f2B;

//   std::cout<<probMultinom(n0,n1,n2,obs0,obs1,obs2)
// 	   <<" "<<n0<<" "<<n1<<" "<<n2<<" "<<obs0<<" "<<obs1<<" "<<obs2
// 	   <<" ... "<<f2S<<" "<<f2B<<std::endl;
  double nlP = -1.*(
		    log (TMath::Poisson(nS+nB,obs)) - log (TMath::Poisson(obs, obs))
		    + log (probMultinom(n0,n1,n2,obs0,obs1,obs2)) - log (probMultinom(obs0,obs1,obs2,obs0,obs1,obs2))
		    + log (TMath::Gaus(f0B, f0B_exp, f0B_err)) 
		    - log( 1. - 0.5*TMath::Prob(TMath::Power(f0B_exp/f0B_err,2),1))
		    + log (TMath::Gaus(f1B, f1B_exp, f1B_err)) 
		    - log( 1. - 0.5*TMath::Prob(TMath::Power(f1B_exp/f1B_err,2),1))
		    + log (TMath::Gaus(f2B, f2B_exp, f2B_err)) 
		    - log( 1. - 0.5*TMath::Prob(TMath::Power(f2B_exp/f2B_err,2),1))
		    + log (TMath::Gaus(nB, nB_exp, nB_err)) 
		    - log( 1. - 0.5*TMath::Prob(TMath::Power(nB_exp/nB_err,2),1))
		    );
//   std::cout<<nS<<" "<<f0S<<" "<<f1S<<" "<<nlP<<std::endl;
  return nlP;
}


void fcn(int& npar, double* gin, double& f, double* par, int iflag){
  f= sigNPL(par[0], par[1], par[2], par[3]);
}
void fcn3(int& npar, double* gin, double& f, double* par, int iflag){
  f= sigNPL3Bins(par[0], par[1]);
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

void minimize3Bins(){
  gMinuit = new TMinuit(4);
  gMinuit->SetFCN(fcn3);
  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 0.5;

  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // Set starting values and step sizes for parameters
  static Double_t vstart[4]  = {50,1.0,0.5};
  Double_t step[4] = {0.01, 0.001, 0.0001,0.0001};
  gMinuit->mnparm(0, "nS", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "sfB", vstart[1], step[1], 0,0,ierflg);
  arglist[0] = 500;
  arglist[1] = 0.01;
  gMinuit->mnexcm("MIGRAD", arglist ,4,ierflg);
//   arglist[0] = 500;
//   arglist[1] = 0.01;
//   gMinuit->mnexcm("MIGRAD", arglist ,4,ierflg);
//   gMinuit->mnexcm("FIX", arglist ,4,ierflg);
//   gMinuit->mnexcm("FIX", arglist ,4,ierflg);
//   gMinuit->mnexcm("MINOS", arglist ,4,ierflg);

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  gMinuit->mnprin(4,amin);
}
