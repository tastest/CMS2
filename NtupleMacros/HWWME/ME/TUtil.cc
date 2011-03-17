#include "fstream"
#include "TMath.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TUtil.hh"
#include "TCanvas.h"

/*
 * Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
 * by D.G. Hook and P.R. McAree
 * from "Graphics Gems", Academic Press, 1990
 * */

#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;


void EtExponential(double alpha,double x0,double* Et,double* wt){
    double sign=1.;
    if(x0>0.5) { 
        sign=-1.;
	x0=(x0-0.5)*2.;
    }else{
        x0=x0*2.;
    }
    double A = 1.-exp(-3500./alpha);
    double B =-alpha*TMath::Log(1-x0*A);

    *wt = 2*alpha*A/exp(-B/alpha);
    *Et = B*sign;
}



void  breitw(double x1,double mminsq,double mmaxsq,double rmass,double rwidth,double *msq,double *wt){
//#include "mcfm/../src/Inc/constants.F"
//c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt
//c---- such that mminsq<msq<mmaxsq
//c---- points are generated around resonance position rmass, but
//c---- breit-wigner should still be included in the matrix element
//c     wt is the jacobian between integration in msq and integration in x1

      double almin,almax,al,tanal;

      int zerowidth=0;
      if (zerowidth){
          tanal=0;
          almax=+TMath::Pi()/2.;
          almin=-TMath::Pi()/2.;
      }else{
          almin=TMath::ATan((mminsq-rmass*rmass)/rmass/rwidth);
          almax=TMath::ATan((mmaxsq-rmass*rmass)/rmass/rwidth);
          al=(almax-almin)*x1+almin;
          tanal=TMath::Tan(al);
      }

      *msq= rmass*rmass+rmass*rwidth*tanal;
      *wt=(almax-almin)*rmass*rwidth*(1.0+tanal*tanal);
      
}



void My_choose(TVar::Process process){
 //Wp_1jet  Wm_1jet
 //Wp_gamma Wm_gamma
 //WW
 //WpZ WmZ 
 //ZZ
 //ZZ_4l
 //HWW
 //HZZ
 if(process==TVar::Wp_1jet ||
    process==TVar::Wm_1jet ){
   //nrpoc =11  '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+f(p5)'
   //nproc=16 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~(p4))+f(p5)'

   //W+ nwz_.nwz= 1; 
   if (process==TVar::Wp_1jet) nwz_.nwz=1;
   //W- nwz_.nwz=-1; 
   if (process==TVar::Wm_1jet) nwz_.nwz=-1;   
   npart_.npart=3;
   nqcdjets_.nqcdjets=1;
   bveg1_mcfm_.ndim=7;
   masses_mcfm_.mb=0;
   breit_.n2=0;
   breit_.n3=1;
   breit_.mass3 =masses_mcfm_.wmass;
   breit_.width3=masses_mcfm_.wwidth;

    
 } 
 else if(process==TVar::Wp_gamma ||
         process==TVar::Wm_gamma ){

   //18 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4))+gamma(p5)'
   //19 '  f(p1)+f(p2) --> W^-(-->e^-(p3)+nu~+(p4))+gamma(p5)'
   //W+ nwz_.nwz= 1; 
   if(process==TVar::Wp_gamma ) nwz_.nwz= 1;
   //W- nwz_.nwz=-1;
   if(process==TVar::Wm_gamma ) nwz_.nwz=-1;
      
   npart_.npart=3; 
   nqcdjets_.nqcdjets=0;
   //case='Wgamma'
   bveg1_mcfm_.ndim=7;
   masses_mcfm_.mb=0;
   breit_.n2=0;
   breit_.n3=1;
   breit_.mass3 =masses_mcfm_.wmass;
   breit_.width3=masses_mcfm_.wwidth;


   //sprintf(plabel[3-1],"nl"); 
   //sprintf(plabel[4-1],"ea");
   //sprintf(plabel[5-1],"ga");
   //sprintf(plabel[6-1],"pp");
 
 }
 else if(process==TVar::WW){ 
 //61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
   //maximum weight: 177475
   //call readcoup
   nwz_.nwz=1;
   npart_.npart=4;
   nqcdjets_.nqcdjets=0;
   bveg1_mcfm_.ndim=10;
   masses_mcfm_.mb=0;
   breit_.n2=1;
   breit_.n3=1;
   breit_.mass2 =masses_mcfm_.wmass;
   breit_.width2=masses_mcfm_.wwidth;
   breit_.mass3 =masses_mcfm_.wmass;
   breit_.width3=masses_mcfm_.wwidth;
   zcouple_.l1=1.;

   //case='WWqqbr'
   //sprintf(plabel[3-1],"nl");  
   //sprintf(plabel[4-1],"ea");
   //sprintf(plabel[5-1],"el");
   //sprintf(plabel[6-1],"na");
   //sprintf(plabel[7-1],"pp");
 
 } 
 else if(process==TVar::WZ || process==TVar::WpZ_lostW || process==TVar::WpZ_lostZ){ 
    //71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'
    //76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))
    //nwz_.nwz=-1;
    //nwz_.nwz=+1;
    //call readcoup
    nwz_.nwz=1;

    npart_.npart=4;
    nqcdjets_.nqcdjets=0;
    //case='WZbbar'
    //call checkminzmass(2)
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;

    zcouple_.q1=-1.;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    //sprintf(plabel[3-1],"nl");	    
    //sprintf(plabel[4-1],"ea");
    //sprintf(plabel[5-1],"ml");
    //sprintf(plabel[6-1],"ma");
    //sprintf(plabel[7-1],"pp");

 }
 else if(process==TVar::ZZ ){ 
 
   //82 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->e^-(p5)+e^+(p6))'
   //87 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->3*(nu(p3)+nu~(p4))) (NO GAMMA*)'
    //call readcoup
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;
    //case='ZZlept'
    //call checkminzmass(1)
    //call checkminzmass(2)
    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    //sprintf(plabel[3-1],"el"); 
    //sprintf(plabel[4-1],"ea");
    //sprintf(plabel[5-1],"nl");
    //sprintf(plabel[6-1],"na");
    //sprintf(plabel[7-1],"pp");

    zcouple_.q2=0.;
    zcouple_.l2=zcouple_.ln*sqrt(3.);
    zcouple_.r2=zcouple_.re*sqrt(3.);
 }  
 else if(process==TVar::ZZ_4l ){ 
 
    //81 '  f(p1)+f(p2) --> Z^0(-->mu^-(p3)+mu^+(p4)) + Z^0(-->e^-(p5)+e^+(p6))'
    //86 '  f(p1)+f(p2) --> Z^0(-->e^-(p5)+e^+(p6))+Z^0(-->mu^-(p3)+mu^+(p4)) (NO GAMMA*)'
    //call readcoup
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;
    //case='ZZlept'
    //call checkminzmass(1)
    //call checkminzmass(2)
    nwz_.nwz=0;
    bveg1_mcfm_.ndim=10;
    masses_mcfm_.mb=0;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2=masses_mcfm_.zmass;
    breit_.width2=masses_mcfm_.zwidth;
    breit_.mass3=masses_mcfm_.zmass;
    breit_.width3=masses_mcfm_.zwidth;

    zcouple_.q1=-1.;
    zcouple_.l1=zcouple_.le;
    zcouple_.r1=zcouple_.re;

    //sprintf(plabel[3-1],"el"); 
    //sprintf(plabel[4-1],"ea");
    //sprintf(plabel[5-1],"ml");
    //sprintf(plabel[6-1],"ma");
    //sprintf(plabel[7-1],"pp");

    zcouple_.q2=-1.;
    zcouple_.l2=zcouple_.le;
    zcouple_.r2=zcouple_.re;
 }  
 else if(process==TVar::HWW ){ 
    //112 '  f(p1)+f(p2) --> H(-->W^- W^+ (for total Xsect)'
    // spira gives us real Higgs width
    //if (spira) then
    //    call higgsp(br,gamgambr,wwbr,zzbr)
    //else
    //    call higgsw(br)
    //endif
    //write(6,99) hmass,hwidth,br
    //case='HWW_4l'
    //sprintf(plabel[3-1],"nl"); 
    //sprintf(plabel[4-1],"ea");
    //sprintf(plabel[5-1],"el");
    //sprintf(plabel[6-1],"na");
    //sprintf(plabel[7-1],"pp");
    npart_.npart=4;
    nqcdjets_.nqcdjets=0;
    nqcdjets_.nqcdstart=7;


    bveg1_mcfm_.ndim=10;
    breit_.n2=1;
    breit_.n3=1;
    breit_.mass2 =masses_mcfm_.wmass;
    breit_.width2=masses_mcfm_.wwidth;
    breit_.mass3 =masses_mcfm_.wmass;
    breit_.width3=masses_mcfm_.wwidth;

    //call branch(brwen,brzee,brtau,brtop)
    //BrnRat=brwen**2*wwbr

  
 }
 else if(process==TVar::HZZ ){ 
    //121 '  f(p1)+f(p2) --> H(-->Z^0(mu^-(p3)+mu^+(p4)) + Z^0(e^-(p5)+e^+(p6))'
    //case='HZZ_4l'
     npart_.npart=4;
     nqcdjets_.nqcdjets=0;
     nqcdjets_.nqcdstart=7;

     bveg1_mcfm_.ndim=10;
     breit_.n2=1;
     breit_.n3=1;
     breit_.mass2 =masses_mcfm_.zmass;
     breit_.width2=masses_mcfm_.zwidth;
     breit_.mass3 =masses_mcfm_.zmass;
     breit_.width3=masses_mcfm_.zwidth;

     //if (spira) then
     //    call higgsp(bbbr,gamgambr,wwbr,br)
     //else
     //    call higgsw(br)
     //endif
     //write(6,99) hmass,hwidth,br

     zcouple_.l1=zcouple_.le;
     zcouple_.r1=zcouple_.re;
     zcouple_.l2=zcouple_.le;
     zcouple_.r2=zcouple_.re;

     //sprintf(plabel[31],"el"); 
     //sprintf(plabel[41],"ea");
     //sprintf(plabel[51],"ml");
     //sprintf(plabel[61],"ma");
     //sprintf(plabel[71],"pp");

 }
 else if(process==TVar::Z_2l){
    //31 '  f(p1)+f(p2) --> Z^0(-->e^-(p3)+e^+(p4))'
    //         case='Z_only'      
             nqcdjets_.nqcdjets=0;
             npart_.npart=2; 
             nwz_.nwz=0;
             breit_.mass3 =masses_mcfm_.zmass ;
             breit_.width3=masses_mcfm_.zwidth;
             breit_.n3=1;
             bveg1_mcfm_.ndim=4;
             //BrnRat=1d0
             //plabel(3)='el'
             //plabel(4)='ea'
             //plabel(5)='pp'
             zcouple_.q1=-1;
             zcouple_.l1=zcouple_.le;
             zcouple_.r1=zcouple_.re;
 
 }
 else{
     std::cerr <<"[My_choose]: Can't identify Process: " << process <<endl;
 } 
}

//###########################
// event slection cuts to avoid divergence
//
//###########################

bool My_eventcuts(TVar::Process process, mcfm_event_type* mcfm_evt, cdf_event_type* cdf_evt){


  TLorentzVector* nu =  &mcfm_evt->p[2];
  TLorentzVector* ep =  &mcfm_evt->p[3];
  

  if(process==TVar::Wp_gamma ||
     process==TVar::Wm_gamma ){
    if(nwz_.nwz==-1){
       ep = &mcfm_evt->p[2]; 
       nu = &mcfm_evt->p[3];
    }
    TLorentzVector* gamma = &mcfm_evt->p[4];
   
    
    if(ep->Pt()           <8   ) return true;
    if(fabs(ep->Eta())    >10  ) return true;
    if(gamma->Pt()        <5   ) return true;
    if(fabs(gamma->Eta()) >10  ) return true;
    if(ep->DeltaR(*gamma) <0.35) return true;
  }
  else if(process==TVar::Wp_1jet ||
          process==TVar::Wm_1jet){

   if(nwz_.nwz==-1){
       ep = &mcfm_evt->p[2];
       nu = &mcfm_evt->p[3];
    }
    TLorentzVector* parton= &mcfm_evt->p[4];
    if(fabs(ep->Eta())     >10  ) return true;
    if(parton->Pt()        <8   ) return true;
    if(fabs(parton->Eta()) >3   ) return true;
    if(ep->DeltaR(*parton) <0.2 ) return true;

  }
  else if(process==TVar::Z_2l){
    double emPt = nu->Pt();
    double epPt = ep->Pt();
    if(emPt    <5 ) return true; //em
    if(epPt    <5 ) return true; //ep
  } 
  return false;
}




double GenGaus(double* rand,Double_t mean, Double_t sigma)
{
   double x = rand[0] * 6.28318530717958623;
   double y = rand[1];
   double z =TMath::Sqrt(-2*TMath::Log(y));
   double results = mean + sigma*TMath::Sin(x)*z;
   return results;
}
void EtGaus(double* rand, double mean, double sigma,double* results, double* wgt){
   double x = rand[0] * 6.28318530717958623;
   double y = rand[1];
   double z =TMath::Sqrt(-2*TMath::Log(y));
   double sinx = TMath::Sin(x);
   *results = mean + sigma*sinx*z;
   *wgt     = sigma*2.5066282746310002/TMath::Exp(-sinx*sinx*z*z/2);
}


bool My_masscuts(double s[][12],TVar::Process process){

 double minZmassSqr=15*15;

 if(process==TVar::WZ || process==TVar::WpZ_lostW || process==TVar::WpZ_lostZ){
   if(s[4][5]< minZmassSqr) return true;
   
 }
 else if(process==TVar::ZZ){
   if(s[2][3]< minZmassSqr) return true;
   if(s[4][5]< minZmassSqr) return true;
 }
 else if(process==TVar::Z_2l){
   if(s[2][3]< minZmassSqr) return true; 
 }
	 
 return false;	 

}
//c----reject event if any s(i,j) is too small
bool My_smalls(double s[][12],int npart){
//cutoff is defined in technical.Dat
	
      if ( 
       npart == 3 &&
       (
        (-s[5-1][1-1]< cutoff_.cutoff)  //gamma p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //gamma p2
     || (-s[4-1][1-1]< cutoff_.cutoff)  //e+    p1
     || (-s[4-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[3-1][1-1]< cutoff_.cutoff)  //nu    p1
     || (-s[3-1][2-1]< cutoff_.cutoff)  //nu    p2
     || (+s[5-1][4-1]< cutoff_.cutoff)  //gamma e+
     || (+s[5-1][3-1]< cutoff_.cutoff)  //gamma nu
     || (+s[4-1][3-1]< cutoff_.cutoff)  //e+    nu
	)	 
      ) 
        return true;
     
     else if (
       npart == 4 &&     
      (
        (-s[5-1][1-1]< cutoff_.cutoff)  //e-    p1
     || (-s[5-1][2-1]< cutoff_.cutoff)  //e-    p2
     || (-s[6-1][1-1]< cutoff_.cutoff)  //nb    p1
     || (-s[6-1][2-1]< cutoff_.cutoff)  //nb    p2
     || (+s[6-1][5-1]< cutoff_.cutoff)  //e-    nb
       )

     )
     {

      return true;
     
     }

     return false;
}




//Make sure
// 1. tot Energy Sum < 2EBEAM
// 2. PartonEnergy Fraction minimum<x0,x1<1
// 3. number of final state particle is defined
//
double SumMatrixElementPDF(TVar::Process process, mcfm_event_type* mcfm_event,double flavor_msq[][nmsq],double* flux){


           int NPart=npart_.npart+2;
	   double p4[4][12];
 	   double fx1[nmsq];
           double fx2[nmsq];
           double msq[nmsq][nmsq];


           //Parton Density Function is always evalualted at pT=0 frame
	   //Make sure parton Level Energy fraction is [0,1]
           //phase space function already makes sure the parton energy fraction between [min,1]
	   //  x0 EBeam =>   <= -x1 EBeam

	   double sysPz=mcfm_event->p[0].Pz()    +mcfm_event->p[1].Pz();
	   double sysE =mcfm_event->p[0].Energy()+mcfm_event->p[1].Energy();

           //Ignore the Pt doesn't make significant effect
	   //double sysPt_sqr=sysPx*sysPx+sysPy*sysPy;
	   //if(sysPt_sqr>=1.0E-10)  sysE=TMath::Sqrt(sysE*sysE-sysPt_sqr);
	   
           double xx[2]={(sysE+sysPz)/EBEAM/2,(sysE-sysPz)/EBEAM/2};
           if(xx[0] > 1.0 || xx[0]<=xmin_.xmin) return 0.0;
	   if(xx[1] > 1.0 || xx[1]<=xmin_.xmin) return 0.0;

           //Convert TLorentzVector into 4x12 Matrix
           //reverse sign of incident partons back
	   for(int ipar=0;ipar<2;ipar++){    
             if(mcfm_event->p[ipar].Energy()>0){
               p4[0][ipar] = -mcfm_event->p[ipar].Px();
               p4[1][ipar] = -mcfm_event->p[ipar].Py();
               p4[2][ipar] = -mcfm_event->p[ipar].Pz();
               p4[3][ipar] = -mcfm_event->p[ipar].Energy();
	     }
           }
           //initialize decayed particles
           for(int ipar=2;ipar<NPart;ipar++){

             p4[0][ipar] = mcfm_event->p[ipar].Px();
             p4[1][ipar] = mcfm_event->p[ipar].Py();
             p4[2][ipar] = mcfm_event->p[ipar].Pz();
             p4[3][ipar] = mcfm_event->p[ipar].Energy();

	   }

    //calculate invariant masses between partons/final state particles
    double s[12][12];
    for(int jdx=0;jdx< NPart ;jdx++){
      s[jdx][jdx]=0;
      for(int kdx=jdx+1;kdx<NPart;kdx++){
        s[jdx][kdx]=2*(p4[3][jdx]*p4[3][kdx]-p4[2][jdx]*p4[2][kdx]-p4[1][jdx]*p4[1][kdx]-p4[0][jdx]*p4[0][kdx]);
        s[kdx][jdx]=s[jdx][kdx];
      }
    }

    //remove events has small invariant mass
    if(My_masscuts(s,process)) return 0.0;
    if(My_smalls(s,npart_.npart)) return 0.0;

	   
           //Calculate Pdf
           //Always pass address through fortran function
           fdist_ (&density_.ih1, &xx[0], &scale_.scale, fx1); //P+=> W+->e+nu
           fdist_ (&density_.ih2, &xx[1], &scale_.scale, fx2); //P-=> W-->e-nu

	   //	   	   cout<<"energy="<<scale_.scale<<"\n";
		   
                if(process==TVar::WW)	      qqb_ww_  (p4[0],msq[0]);    
           else if(process==TVar::WZ ||
                   process==TVar::WpZ_lostW ||
                   process==TVar::WpZ_lostZ ) qqb_wz_  (p4[0],msq[0]);
           else if(process==TVar::ZZ  ||
	           process==TVar::ZZ_4l)      qqb_zz_  (p4[0],msq[0]);
	   else if(process==TVar::Wp_gamma ||
                   process==TVar::Wm_gamma)   qqb_wgam_(p4[0],msq[0]);
           else if(process==TVar::Wp_1jet || 
	           process==TVar::Wm_1jet )   qqb_w_g_ (p4[0],msq[0]);
           else if(process==TVar::HWW)	      qqb_hww_ (p4[0],msq[0]);
           else if(process==TVar::HZZ)	      qqb_hzz_ (p4[0],msq[0]);
           else if(process==TVar::Z_2l)       qqb_z_   (p4[0],msq[0]);     
                     
	   double msqjk=0;
	   for(int ii=0;ii<11;ii++){
	   for(int jj=0;jj<11;jj++){

	     //2-D matrix is reversed in fortran
             // msq[ parton2 ] [ parton1 ]
	     flavor_msq[jj][ii] = fx1[ii]*fx2[jj]*msq[jj][ii];
             msqjk+=flavor_msq[jj][ii];
	   }//ii
	   }//jj
	   
	   (*flux)=fbGeV2/(8*xx[0]*xx[1]*EBEAM*EBEAM);
	  
           if(msqjk != msqjk || flux!=flux ){
              cout << "SumMatrixPDF: "<< TVar::ProcessName(process) << " msqjk="  << msqjk << " flux="<< flux <<endl;
              msqjk=0;
              flux=0;
           }
	   return msqjk;

}


double HiggsWidth(double mass){

     double width=0.80720E-01; 


     if (mass==100) width=0.2598E-02;
  else 
    if (mass==110) width=0.29590E-02;
  else 
    if (mass==120) width=0.36050E-02;
  else 
    if (mass==130) width=0.49630E-02;
  else 
    if (mass==140) width=0.81090E-02;
  else 
    if (mass==150) width=0.16860E-01;
  else 
    if (mass==160) width=0.80720E-01; 
  else 
    if (mass==170) width=0.38660E+00;
  else 
    if (mass==180) width=0.63030E+00;
  else 
    if (mass==190) width=0.10400E+01;
  else 
    if (mass==200) width=0.14280E+01;
  else 
    if (mass==210) width=0.18400E+01;
  else 
    if (mass==220) width=2.29900E+00;
  else 
    if (mass==230) width=2.81500E+00;
  else 
    if (mass==240) width=3.39700E+00;
  else
    if (mass==250) width=4.04900E+00;
  else
    if (mass==260) width=4.77500E+00;
  else 
    if (mass==270) width=5.58300E+00;
  else 
    if (mass==280) width=6.47200E+00;
  else
    if (mass==290) width=7.44700E+00;
  else
    if (mass==300) width=8.51200E+00;

  //5 GeV steps
  else
    if (mass==115) width=0.3228E-02;
  else
    if (mass==125) width=0.4152E-02;
  else
    if (mass==135) width=0.6192E-02;
  else
    if (mass==145) width=0.1124E-01;
  else
    if (mass==155) width=0.2919E-01;
  else
    if (mass==165) width=0.2546E+00;
  else
    if (mass==175) width=0.5037E+00;
  else
    if (mass==185) width=0.8307E+00;
  else
    if (mass==195) width=1.2330E+00;

 return width;




}

double SetTGCParameter(TString name,double value){
     if(name=="delg1_z" ){ anomcoup_.delg1_z =value; printf("Set anomcoup_.delg1_z = %8.3f\n",anomcoup_.delg1_z );}
else if(name=="delg1_g" ){ anomcoup_.delg1_g =value; printf("Set anomcoup_.delg1_g = %8.3f\n",anomcoup_.delg1_g );}
else if(name=="lambda_g"){ anomcoup_.lambda_g=value; printf("Set anomcoup_.lambda_g= %8.3f\n",anomcoup_.lambda_g);}
else if(name=="lambda_z"){ anomcoup_.lambda_z=value; printf("Set anomcoup_.lambda_z= %8.3f\n",anomcoup_.lambda_z);}
else if(name=="delk_g"  ){ anomcoup_.delk_g  =value; printf("Set anomcoup_.delk_g  = %8.3f\n",anomcoup_.delk_g  );}
else if(name=="delk_z"  ){ anomcoup_.delk_z  =value; printf("Set anomcoup_.delk_z  = %8.3f\n",anomcoup_.delk_z  );}
else if(name=="tevscale"){ anomcoup_.tevscale=value; printf("Set anomcoup_.tevscale= %8.3f\n",anomcoup_.tevscale);}
else {

 cerr <<"No such parameter " << name <<endl;

}

     return 0.0;
}


double GetTGCParameter(int i){
if(i==0)return anomcoup_.delg1_z;
if(i==1)return anomcoup_.delg1_g;
if(i==2)return anomcoup_.lambda_g;
if(i==3)return anomcoup_.lambda_z;
if(i==4)return anomcoup_.delk_g;
if(i==5)return anomcoup_.delk_z;
if(i==6)return anomcoup_.tevscale;
return 999;
}


double getProbAcceptanceEfficiency(cdf_event_type cdf_event, EffHist effhist)
{
  double eff = 1.0; 
  
  for (int i=0;i<2;i++) {

    if ( TMath::Abs(cdf_event.PdgCode[i]) == 11)
      eff = eff * effhist.els_eff_mc->GetBinContent( effhist.els_eff_mc->GetXaxis()->FindBin(cdf_event.p[i].Eta()), effhist.els_eff_mc->GetYaxis()->FindBin(cdf_event.p[i].Pt()) );
    
    if ( TMath::Abs(cdf_event.PdgCode[i]) == 13)
      eff = eff * effhist.mus_eff_mc->GetBinContent( effhist.mus_eff_mc->GetXaxis()->FindBin(cdf_event.p[i].Eta()), effhist.mus_eff_mc->GetYaxis()->FindBin(cdf_event.p[i].Pt()) );
  }
  
  return eff;
}

void getProbFromHist(double x0, double* kX, double *wt, TH1F *hkx)
{
  double c = hkx->GetXaxis()->GetXmax() - hkx->GetXaxis()->GetXmin();
  *kX = hkx->GetXaxis()->GetXmin() + x0*c;
  *wt = hkx->GetBinContent(hkx->GetXaxis()->FindBin(*kX))/hkx->Integral(1, hkx->GetNbinsX())*c/hkx->GetBinWidth(1);
}

