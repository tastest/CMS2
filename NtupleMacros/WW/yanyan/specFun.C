double parabolicCF(double nu, double z){
  double res = pow(2.,nu/2.)*sqrt(TMath::Pi())*exp(-z*z/4.)
    * (1./TMath::Gamma((1.-nu)/2.)*ROOT::Math::conf_hyperg(-nu/2., 1./2., z*z/2.) 
       - sqrt(2.)*z/TMath::Gamma(-nu/2.)*ROOT::Math::conf_hyperg((1.-nu)/2., 3./2., z*z/2.));

  return res;
}

double conf_hypergU_mod(double a, double b, double x){
  if (a < 7.) return ROOT::Math::conf_hypergU(a,b,x);
  else {
    return conf_hypergU_mod(a-2.,b,x)/(1.-a)/(a-b)
      + conf_hypergU_mod(a-1.,b,x)*(2.-2.*a + b - x)/(1.-a)/(a-b);
  }
}

double poisson_smeared_prob(double n, double x0, double sigma){

  double res= 0 ;
  if ((sigma - x0/sigma)>0){
    res = exp(-x0)
      *pow(sigma,n)
      *exp(x0*(1.-x0/2./sigma/sigma))
      /sqrt(TMath::Pi())
      /pow(2.,n/2.)
      /(1. + TMath::Erf(x0/sqrt(2.)/sigma))
      * conf_hypergU_mod((n+1.)/2., 1./2., 0.5*(sigma - x0/sigma)*(sigma - x0/sigma));
  } else if ((sigma - x0/sigma)<0) {
    double hypergVal = 0;
    if (fabs((sigma - x0/sigma)) <=1){
      res = exp(-x0)
	*pow(sigma,n)
	*exp(x0*(1.-x0/2./sigma/sigma))
	/sqrt(TMath::Pi())
	/pow(2.,n/2.)
	/(1. + TMath::Erf(x0/sqrt(2.)/sigma))
	* ( 2.*sqrt(TMath::Pi())/TMath::Gamma((2.+n)/2.)*ROOT::Math::conf_hyperg((n+1.)/2., 1./2.,  0.5*(sigma - x0/sigma)*(sigma - x0/sigma))
	    - conf_hypergU_mod((n+1.)/2., 1./2., 0.5*(sigma - x0/sigma)*(sigma - x0/sigma)));
    } else {
      res = exp(-0.5 * x0 - 0.5*sigma*fabs(sigma-x0/sigma))
	*pow(sigma,n)
	/sqrt(TMath::Pi())
	/pow(2.,n/2.)
	/(1. + TMath::Erf(x0/sqrt(2.)/sigma))
	* ( 2.*sqrt(TMath::Pi())/TMath::Gamma((2.+n)/2.)*ROOT::Math::conf_hyperg(-n/2., 1./2., - 0.5*(sigma - x0/sigma)*(sigma - x0/sigma) )
	    - exp(-0.5*(sigma - x0/sigma)*(sigma - x0/sigma))*conf_hypergU_mod((n+1.)/2., 1./2., 0.5*(sigma - x0/sigma)*(sigma - x0/sigma)));
    }
  } else {
    res = exp(-x0/2.)
      *pow(sigma,n)
      /pow(2.,n/2.)
      /(1. + TMath::Erf(sigma/sqrt(2.)))
      /TMath::Gamma((2.+n)/2.);
  }
  return res;
}
