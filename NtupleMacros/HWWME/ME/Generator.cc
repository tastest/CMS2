#include "iostream"
#include "Generator.hh"


//______________________________________________________________________________
double GetRandom(TH1* his, double r1) 
{
   // return a random number distributed according the histogram bin contents.
   // This function checks if the bins integral exists. If not, the integral
   // is evaluated, normalized to one.
   // The integral is automatically recomputed if the number of entries
   // is not the same then when the integral was computed.
   // NB Only valid for 1-d histograms. Use GetRandom2 or 3 otherwise.

   
   int nbinsx = his->GetNbinsX();
   double *fIntegral = his->GetIntegral();
   if (!fIntegral) {
      int integral = (int) his->ComputeIntegral();
      if (integral == 0 || !fIntegral ){ 
         cout << "Error Integration of Histogram " << his->GetName() <<" with error code " << integral<<endl;
         return 0;
      }
   } 
   int ibin = TMath::BinarySearch(nbinsx,fIntegral,r1);
   double x = his->GetBinLowEdge(ibin+1);
   if (r1 > fIntegral[ibin]) x +=
      his->GetBinWidth(ibin+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
   return x;
}




void  GenMwModified(double r, double mminsq, double mmaxsq,
                    double* msq, double* wgt){

//#include "mcfm/../src/Inc/constants.F"
//c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt
//c---- such that mminsq<msq<mmaxsq
//c---- points are generated around resonance position rmass, but
//c---- breit-wigner should still be included in the matrix element
//c     wt is the jacobian between integration in msq and integration in x1
      
   double M1sqr =      74.0*74.0 ;
   double f     =      0.3 ;
   double rmass =   8.04000e+01 ;
   double rwidth=  2.78183e+00 ;

 if(mmaxsq<M1sqr){
    //double tmp  = TMath::Exp((mminsq-mmaxsq)/W1);
    //double A    = 1./(1.-tmp)/W1;
    //double val  = r/W1+A*tmp;
    //*msq = mmaxsq+W1*TMath::Log(val/A);
    //*wgt = 1./val;

    double B=3.61614e-06;
    double A=6.53757878e-08;
    
    double NormA = (B+A*(mmaxsq+mminsq))*(mmaxsq-mminsq)/2.;
     
    double tmp= (B+A*mminsq)*mminsq/2.;
    *msq = (-B+TMath::Sqrt(B*B+2*A*(NormA*r+tmp)))/A;
    *wgt = NormA/(B+A*(*msq));
 }
 else if(mminsq>=M1sqr){

   double almin=TMath::ATan((M1sqr -rmass*rmass)/rmass/rwidth);
   double almax=TMath::ATan((mmaxsq-rmass*rmass)/rmass/rwidth);
   double al=(almax-almin)*r+almin;
   double tanal=TMath::Tan(al);
   *msq= rmass*rmass+rmass*rwidth*tanal;
   *wgt=(almax-almin)*rmass*rwidth*(1.0+tanal*tanal);

 }
 else {

   double Norm = (1.0+f); 
  
   if(r*Norm<1.){
//     double tmp  = TMath::Exp((mminsq-M1sqr)/W1);
//     double A    = 1./(1.-tmp)/W1;
//     double val  = r/W1+A*tmp/Norm;
//     *msq = M1sqr+W1*TMath::Log(val*Norm/A);
//     *wgt = 1./val;

    double B=3.61614e-06;
    double A=6.53757878e-08;
    
    double NormA = (B+A*(M1sqr+mminsq))*(M1sqr-mminsq)/2.;
     
    double tmp= (B+A*mminsq)*mminsq/2.;
    *msq = (-B+TMath::Sqrt(B*B+2*A*(Norm*NormA*r+tmp)))/A;
    *wgt = Norm*NormA/(B+A*(*msq));
   }
   else
     {
     double almin=TMath::ATan((M1sqr -rmass*rmass)/rmass/rwidth);
     double almax=TMath::ATan((mmaxsq-rmass*rmass)/rmass/rwidth);
     r=(r*Norm-1.)/f;
     double al=(almax-almin)*r+almin;
     double tanal=TMath::Tan(al);
     *msq= rmass*rmass+rmass*rwidth*tanal;
     *wgt=(almax-almin)*rmass*rwidth*(1.0+tanal*tanal)*Norm/f;
     }
 }
   
}


