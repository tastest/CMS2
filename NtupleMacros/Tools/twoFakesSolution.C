//solutions are taken from AN-2010-261, the notation is preserved to some extent

//same flavor: returns N prompt-prompt passing numerator cuts
// inputs: n2 both numerators, n1 only one passes numerator cuts, n0 both leptons fail numerator
// inputs: eff -- efficiency of prompt/signal wrt denominator cuts; fr -- fake rate
double frSF_pp(double n2, double n1, double n0, double eff, double fr){
  double epsilon = fr/(1. - fr);
  double eta = (1. - eff)/eff;

  double frac = (1.- epsilon*eta);
  frac*= frac;
  frac = 1./frac;

  return (n2 - epsilon*n1 + epsilon*epsilon*n0)*frac;
}

//same flavor: returns N prompt-fake passing numerator cuts
// inputs: n2 both numerators, n1 only one passes numerator cuts, n0 both leptons fail numerator
// inputs: eff -- efficiency of prompt/signal wrt denominator cuts; fr -- fake rate
double frSF_pf(double n2, double n1, double n0, double eff, double fr){
  double epsilon = fr/(1. - fr);
  double eta = (1. - eff)/eff;

  double frac = (1.- epsilon*eta);
  frac*= frac;
  frac = 1./frac;

  return epsilon*(-2.*epsilon*n0 + (1.+epsilon*eta)*n1 - 2.*eta*n2)*frac;
}

// same flavor: returns N fake-fake passing numerator cuts
// inputs: n2 both numerators, n1 only one passes numerator cuts, n0 both leptons fail numerator
// inputs: eff -- efficiency of prompt/signal wrt denominator cuts; fr -- fake rate
double frSF_ff(double n2, double n1, double n0, double eff, double fr){
  double epsilon = fr/(1. - fr);
  double eta = (1. - eff)/eff;

  double frac = (1.- epsilon*eta);
  frac*= frac;
  frac = 1./frac;

  return epsilon*epsilon*(n0 - eta*n1 + eta*eta*n2)*frac;
}

// opposite flavor: returns N prompt-prompt passing numerator cuts
// inputs: nE1M1 -- both numerators; n10 e-num mu-nonnum; n01 e-nonum mu-num; n0 -- both nonnumerators
//         effE (effM) -- efficiency of prompt/signal wrt denominator cuts for electrons (muons)
//         frE (frM) -- fake rate for electrons (muons)
double frOF_EpMp(double nE1M1, double nE1M0, double nE0M1, double nE0M0, 
	       double effE, double effM, double frE, double frM){

  double epsilonE = frE/(1. - frE);
  double epsilonM = frM/(1. - frM);

  double etaE = (1. - effE)/effE;
  double etaM = (1. - effM)/effM;

  double frac = 1./( (1. - epsilonE*etaE)* (1. - epsilonM*etaM ) );

  return frac* ( nE1M1 - epsilonM*nE1M0 - epsilonE*nE0M1 + epsilonE*epsilonM*nE0M0);
}

// opposite flavor: returns N e-prompt mu-fake passing numerator cuts
// inputs: nE1M1 -- both numerators; n10 e-num mu-nonnum; n01 e-nonum mu-num; n0 -- both nonnumerators
//         effE (effM) -- efficiency of prompt/signal wrt denominator cuts for electrons (muons)
//         frE (frM) -- fake rate for electrons (muons)
double frOF_EpMf(double nE1M1, double nE1M0, double nE0M1, double nE0M0, 
	       double effE, double effM, double frE, double frM){

  double epsilonE = frE/(1. - frE);
  double epsilonM = frM/(1. - frM);

  double etaE = (1. - effE)/effE;
  double etaM = (1. - effM)/effM;

  double frac = 1./( (1. - epsilonE*etaE)* (1. - epsilonM*etaM ) );

  return frac*epsilonM * ( -epsilonE*nE0M0 + nE1M0 + epsilonE*etaM*nE0M1 - etaM*nE1M1);
}

// opposite flavor: returns N e-fake mu-prompt passing numerator cuts
// inputs: nE1M1 -- both numerators; n10 e-num mu-nonnum; n01 e-nonum mu-num; n0 -- both nonnumerators
//         effE (effM) -- efficiency of prompt/signal wrt denominator cuts for electrons (muons)
//         frE (frM) -- fake rate for electrons (muons)
double frOF_EfMp(double nE1M1, double nE1M0, double nE0M1, double nE0M0, 
	       double effE, double effM, double frE, double frM){

  double epsilonE = frE/(1. - frE);
  double epsilonM = frM/(1. - frM);

  double etaE = (1. - effE)/effE;
  double etaM = (1. - effM)/effM;

  double frac = 1./( (1. - epsilonE*etaE)* (1. - epsilonM*etaM ) );

  return frac* epsilonE * (  -epsilonM*nE0M0 + etaE*epsilonM*nE1M0 + nE0M1 - etaE*nE1M1);
}

// opposite flavor: returns e-fake mu-fake passing numerator cuts
// inputs: nE1M1 -- both numerators; n10 e-num mu-nonnum; n01 e-nonum mu-num; n0 -- both nonnumerators
//         effE (effM) -- efficiency of prompt/signal wrt denominator cuts for electrons (muons)
//         frE (frM) -- fake rate for electrons (muons)
double frOF_EfMf(double nE1M1, double nE1M0, double nE0M1, double nE0M0, 
	       double effE, double effM, double frE, double frM){

  double epsilonE = frE/(1. - frE);
  double epsilonM = frM/(1. - frM);

  double etaE = (1. - effE)/effE;
  double etaM = (1. - effM)/effM;

  double frac = 1./( (1. - epsilonE*etaE)* (1. - epsilonM*etaM ) );

  return frac* epsilonE * epsilonM * (  nE0M0 -  etaE*nE1M0 - etaM * nE0M1 + etaE*etaM*nE1M1);
}
