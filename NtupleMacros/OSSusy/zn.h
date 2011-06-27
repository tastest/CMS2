// program to calculate Zn significance
//***** program to calculate ScP *****

#ifndef CMS2_ZN_H 
#define CMS2_ZN_H 

double sigcpunc( double, double, double);
double sigcp(double, double);
double sumP( double, double);
double poisson_pdf(double, double);
double lgamma(double);
double gausin(double);
double freq(double);
double getzn(double, double, double, double, bool printout = false);
 
#endif

