// program to calculate Zn significance

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
long double sigcpunc( long double, long double, long double);
long double sigcp(long double, long double);
long double sumP( long double, long double);
long double poisson_pdf(long double, long double);
long double lgamma(long double);
long double gausin(long double);
long double freq(long double);


//***** program to calculate ScP *****

main(int argc, char **argv)
{
  long double Ns, Nb, Nsb, Nb2, Sbgr2, Dbgr, Dbgr0;
  long double Scl, S12, S121, S122, S123;

  if( argc < 3 ) {
    printf("\n  Usage: zn.exe Nsig Nbgr [Sbgr^2 [Dbgr]]\n");
    printf("  where  Nsig - number of signal events\n");
    printf("         Nbgr - number of background events\n");
    printf("         Sbgr - statistical background uncertainty\n");
    printf("                (usually Sbgr^2 = Nbgr)\n");
    printf("         Dbgr - systematic background uncertainty (in percent)\n\n");
    exit(1);
  }

  sscanf(argv[1],"%Lf",&Ns);
  sscanf(argv[2],"%Lf",&Nb);

  if( Ns <= 0. || Nb <= 0. ) {
    printf("\nIt's not an interesting case ");
    printf("from the physical point of view!\n\n");
    exit(1);
  }
  
  if( argc > 3 ) sscanf(argv[3],"%Lf",&Sbgr2);
  if( argc > 4 ) sscanf(argv[4],"%Lf",&Dbgr);
  else Dbgr = 0.;
  Dbgr0 = Nb * Dbgr / 100.0;
  Nb2 = Nb + Dbgr0;

// display input parameters

  printf("\nSignal Nsig = %Lg\nBackground Nbgr = %Lg\n",Ns,Nb);
  if( argc > 3 ) {
    printf("Statistical uncertainty Sbgr^2 = %Lg",Sbgr2);
    printf(" (Sbgr = %Lg)\n",sqrtl(Sbgr2));
  }
  if( argc > 4 )
    printf("Systematic uncertainty Dbgr = Nbgr * %Lg\% = %Lg\n",Dbgr,Dbgr0);
  
// counting significance Scl

  Nsb = Ns + Nb;
  Scl = sqrtl(2.0*(Nsb*logl(Nsb/Nb) - Ns));
  printf("\n Significance  Scl  =  sqrt(2*((Nsig+Nbgr)*");
  printf("ln(1+Nsig/Nbgr)-Nsig)) = %.3Lf\n",Scl);

// counting significance S12

  S12 = 2.0*(sqrtl(Nsb) - sqrtl(Nb));
  printf("\n Significance  S12    = 2 * (sqrt(Nsig + Nbgr)");
  printf(" - sqrt(Nbgr))   = %.3Lf\n",S12);

  if( argc > 3 ) {
    S121 = S12 * sqrtl(Nb) / sqrtl(Nb+Sbgr2);
    printf("+ statistical: S12'   = ");   
    printf("S12 * sqrt(Nbgr / (Nbgr + Sbgr2))      = %.3Lf\n",S121);   
  }

  if( argc > 4 ) {
    S122 = 2.0*(sqrtl(Nsb) - sqrtl(Nb2));
    printf("+ systematic:  S12''  = ");   
    printf("2 * (sqrt(Nsig+Nbgr) - sqrt(Nbgr+Dbgr))   = %.3Lf\n",S122);   
    S123 = S122*sqrtl(Nb2)/sqrtl(Nb2+Sbgr2);
    printf("+ stat + syst: S12''' = ");   
    printf("S12''*sqrt((Nbgr+Dbgr)/(Nbgr+Sbgr2+Dbgr)) = %.3Lf\n",S123);   
  }

// counting significance ScP

  if( S12 < 7.0 ) {

    printf("\n Significance  ScP = %.3Lf\n",sigcp(Ns+Nb,Nb));
    if( argc > 3 ) {
      printf("+ statistical: ScP = %.3Lf\n",sigcpunc(Ns+Nb,Nb,Sbgr2));
    }

    if( argc > 4 ) {
      printf("+ systematic:  ScP = %.3Lf\n",sigcp(Ns+Nb,Nb2));
      printf("+ stat + syst: ScP = %.3Lf\n",sigcpunc(Ns+Nb,Nb2,Sbgr2));
    }
  }
  else printf("\n Significance  ScP: use S12 instead!\n");
  printf("\n");
}

 
//***** significance ScP with uncertainties *****

long double sigcpunc( long double Nsb,
                      long double Nb1,
		      long double Sbgr2 )
{
  if( Sbgr2 <= 0.0 ) return sigcp(Nsb,Nb1);

  long double Sbgr = sqrtl(Sbgr2);
  
// to average background around Nbgr+Dbgr in range +/- 3 Sbgr

  long double sumW = 0.0;
  long double sumS = 0.0;

  int Nst = 60;
  long double dSbgr = 3.0*Sbgr/Nst;
  long double Sig2 = 2.0*Sbgr2;
  
  for( int i= -Nst; i<=Nst; i++ ) {
    long double x1 = dSbgr*i;  // x1: from -3 Sbgr to +3 Sbgr
    if( Nb1+x1 > 0. ) {
      long double weight = expl(-x1*x1/Sig2);
      sumW += weight;
      sumS += weight * sumP(Nsb, Nb1+x1);
    }
  } 
  sumS /= sumW;  // averaged value of Poisson tail
  
  return gausin(sumS);
}


//***** significance ScP *****

long double sigcp( long double Nsb, long double Nb )
{
//  return gausin(sumP(Nsb, Nb));
  return -gausin(1.0-sumP(Nsb, Nb));
}


//***** sum of Poisson tail *****

long double sumP( long double Nsb, long double Nb )
{
  long int N = (long int)Nsb;

// initial value takes into account the case Nsb is non-integer 
  long double sum = 0.5*(poisson_pdf((long double)N,Nb)+poisson_pdf(Nsb,Nb))*(Nsb-N);
  
// start from the least members of row to increase accuracy
  for ( long int i = N-1; i >=0; i-- ) sum += poisson_pdf(i, Nb);
  
  return sum;
}


//***** probability of Poisson distribution *****

long double poisson_pdf(long double x, long double mu) {

     return expl(x*logl(mu) - lgamma(x+1.) - mu);
}


//***** logarithm of gamma-function *****

long double lgamma(long double z) {
 
  // Computation of ln[gamma(z)] for all z>0.
  //
  // C.Lanczos, SIAM Journal of Numerical Analysis B1 (1964), 86.
  //
  // The accuracy of the result is better than 2e-10.
  //
  //--- Nve 14-nov-1998 UU-SAP Utrecht

    if (z<=0) return 0;
 
    // Coefficients for the series expansion
    long double c[7] = { 2.5066282746310005, 76.18009172947146,
                  -86.50532032941677,   24.01409824083091,
                   -1.231739572450155, 0.1208650973866179e-2,
                   -0.5395239384953e-5};
 
    long double x   = z;
    long double y   = x;
    long double tmp = x+5.5;
    tmp = (x+0.5)*logl(tmp)-tmp;
    long double ser = 1.000000000190015;
    for (int i=1; i<7; i++) {
      y   += 1.0;
      ser += c[i]/y;
    }
    long double v = tmp+logl(c[0]*ser/x);
    return v;
}


//***** inverse normal distribution *****

//     Computes a "Normal Deviate"
//     Based on G.W. Hill & A.W. Davis, Algorithm 442 Normal Deviate
//     Collected Algorithms from CACM

long double gausin( long double p )
{
  long double c = 2.50662827463100050e0;
  long double z1 = 1.0;
  long double half = z1/2;
  long double c1 = 3*z1/4;
  long double c2 = 7*z1/8;
  long double c3 = z1/3;

  long double x, z, sq;

  long double h = 0.;
  
  if( p <= 0. || p >= 1.0 ) {
    printf("GAUSIN: Argument P = %f is out of range!\n",p);
    exit(1);
  }
  
  if( p == half ) {
    h = 0.;
  }
  else {
    x = p;
    if( p > half ) x = 1. - p;
    x = sqrtl(-2.0*log(x));
    x -= ((7.47395*x+494.877)*x+1637.720)/(((x+117.9407)*x+908.401)*x+659.935);
    if( p < half ) x = -x;
    sq = x*x;
  
    z = c * (p - freq(x))*expl(half*sq);
    h = (((((c1*sq+c2)*z+x)*x+half)*c3*z+half*x)*z+1.0)*z+x;
  }
  return h;
}


//***** frequency calculation *****

long double freq( long double X )
{
  long double Z1 = 1.;
  long double HF = Z1/2;
  long double C1 = 0.56418958354775629e0;
  long double W2 = 1.41421356237309505e0;
  long double RW2 = 1./W2;

  long double P1[4] = {+2.4266795523053175e+2, +2.1979261618294152e+1,
                  +6.9963834886191355e+0, -3.5609843701815385e-2};
		  
  long double Q1[4] = {+2.1505887586986120e+2, +9.1164905404514901e+1,
                  +1.5082797630407787e+1, +1.0};

  long double P2[8] = {+3.00459261020161601e+2,+4.51918953711872942e+2,
                  +3.39320816734343687e+2,+1.52989285046940404e+2,
                  +4.31622272220567353e+1,+7.21175825088309366e+0,
                  +5.64195517478973971e-1,-1.36864857382716707e-7};

  long double Q2[8] = {+3.00459260956983293e+2,+7.90950925327898027e+2,
                  +9.31354094850609621e+2,+6.38980264465631167e+2,
                  +2.77585444743987643e+2,+7.70001529352294730e+1,
                  +1.27827273196294235e+1,+1.0};

  long double P3[5] = {-2.99610707703542174e-3, -4.94730910623250734e-2,
                  -2.26956593539686930e-1, -2.78661308609647788e-1,
                  -2.23192459734184686e-2};

  long double Q3[5] = {+1.06209230528467918e-2, +1.91308926107829841e-1,
                  +1.05167510706793207e+0, +1.98733201817135256e+0,
                  +1.0};

  long double V, Y, AP, AQ, H, HC, FREQ;
  int I; 

  V = RW2 * fabsl(X);
  if( V < HF ) {
    Y = V*V;
    AP = P1[3];
    AQ = Q1[3];
    for( I = 2; I>=0; I-- ) {
      AP = P1[I] + Y * AP;
      AQ = Q1[I] + Y * AQ;
    }
    H = V * AP / AQ;
    HC = 1. - H;
  }
  else {
    if( V < 4.0) {
      AP = P2[7];
      AQ = Q2[7];
      for( I = 6; I >= 0; I-- ) {
        AP = P2[I] + V * AP;
        AQ = Q2[I] + V * AQ;
      }
      HC = exp(-V*V) * AP / AQ;
      H = 1.0 - HC;
    }
    else {
      Y = 1.0 / (V*V);
      AP = P3[4];
      AQ = Q3[4];
      for( I = 3; I >= 0; I-- ) {
        AP = P3[I] + Y * AP;
        AQ = Q3[I] + Y * AQ;
      }
      HC = expl(-V*V) * (C1 + Y * AP / AQ) / V;
      H = 1.0 - HC;
    }
  }
  if( X > 0. ) FREQ = HF + HF * H;
  else FREQ = HF * HC;

  return FREQ;
  
}

// functions gausin and freq are translated from fortran by YuA
