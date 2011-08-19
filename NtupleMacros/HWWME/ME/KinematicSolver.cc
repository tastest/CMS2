
#include "TMCFM.hh"
#include "KinematicSolver.hh"
#include "PolySolver.hpp"
#include "iostream"

void WWL1L2Sol_MHiggsMw1( cdf_event_type* temp,
                    double qX, double qY , 
		    double MHiggs_sqr, double Mw1_sqr, double nuX, double nuY,
                    mcfm_event_type* sol){

 double epX=temp->p[0].Px() ;
 double epY=temp->p[0].Py() ;
 double epZ=temp->p[0].Pz() ;
 double epE=temp->p[0].Energy();

 double emX=temp->p[1].Px() ;
 double emY=temp->p[1].Py() ;
 double emZ=temp->p[1].Pz() ;
 double emE=temp->p[1].Energy();


 double nbX=qX-epX-emX-nuX; 
 double nbY=qY-epY-emY-nuY; 

			
 double limits[2]={-EBEAM,EBEAM};
 double roots_nuZ [2];
 double roots_nbZ [2][2];
 double jacobian[4];
 for(int i=0;i<4;i++) jacobian[i]=1;
			
 //Solve nuZ 
 double B=Mw1_sqr/2+epX*nuX+epY*nuY;
 
 double k_nuZ[10];
 k_nuZ[2]=(epZ-epE)*(epZ+epE);
 k_nuZ[1]= 2*epZ*B;
 k_nuZ[0]= B*B-epE*epE*(nuX*nuX+nuY*nuY);
 PolySolver::FindRoots(2,k_nuZ,limits,roots_nuZ);
/*
std::cout <<"H1Debug in  epZ " << epZ <<" epE " << epE  <<std::endl;
std::cout <<"H1Debug in  nuX " << nuX << " nuY " << nuY <<" B "<< B <<" Mw "<< Mw1_sqr<<endl;
std::cout <<"H1Debgu eqnnuZ " << k_nuZ[2]<<" "<<k_nuZ[1] <<" "<<k_nuZ[0] <<std::endl;
std::cout <<"H1Debug solnuZ " << roots_nuZ[0]<<" "<<roots_nuZ[1]<<std::endl; 
*/
 //Solve nbZ for each nuZ  
 for(int i=0;i<2;i++){
    double nuZ = roots_nuZ[i];

    double eval_nuZ = k_nuZ[0] + k_nuZ[1]*nuZ+ k_nuZ[2]*nuZ*nuZ;
    if(fabs(eval_nuZ/k_nuZ[2])>1){
      jacobian[2*i]   =-1;
      jacobian[2*i+1] =-1;
      roots_nbZ[i][0]=0;
      roots_nbZ[i][1]=0;      
      continue;
    }
    //Debug
   /* 
     TLorentzVector tempNu;
     tempNu.SetXYZM(nuX,nuY,nuZ,0);
     TLorentzVector WpP4 = tempNu + temp->p[0];
     cout <<"Debug: nuZ=" 
          << nuZ 
	  <<" eval/k[2]=" 
	  << eval_nuZ/k_nuZ[2] 
	  << " WpM=" << WpP4.M() <<" =?= " << TMath::Sqrt(Mw1_sqr) <<endl;
    */

    double nuE = TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
    double E3  = epE+emE+nuE;
    double Pz  = epZ+emZ+nuZ;
    double C = MHiggs_sqr+qX*qX+qY*qY-E3*E3+Pz*Pz-nbX*nbX-nbY*nbY;

    double k_nbZ[3];
    k_nbZ[2]=4*(E3+Pz)*(E3-Pz);
    k_nbZ[1]=-4*Pz*C;
    k_nbZ[0]=4*E3*E3*(nbX*nbX+nbY*nbY)-C*C;

    double roots[10];
    PolySolver::FindRoots(2,k_nbZ,limits,roots);

    for(int j=0;j<2;j++){
        double nbZ=roots[j];    
        double eval_nbZ = k_nbZ[0]+k_nbZ[1]*nbZ+k_nbZ[2]*nbZ*nbZ;
        if(fabs(eval_nbZ/k_nbZ[2]) <1){
            roots_nbZ[i][j]=nbZ;
	}else{
            jacobian[2*i+j]=-2;
	    roots_nbZ[i][j]=0;
	}
	    //Debug
    /* 	    
     TLorentzVector tempNb;
     tempNb.SetXYZM(nbX,nbY,nbZ,0);
     TLorentzVector HiggsP4 = tempNu +tempNb + temp->p[0] + temp->p[1];
     cout <<"Debug: nbZ=" 
          << nbZ 
	  <<" => "<<k_nbZ[0] <<"+"<<k_nbZ[1]<<"*x+"<<k_nbZ[2] <<" *X^2=0 => " 
	  <<" eval/k[2]=" 
	  << eval_nbZ/k_nbZ[2] 
	  << " HiggsM=" << HiggsP4.M() <<" =?= " << TMath::Sqrt(MHiggs_sqr) <<endl;
     */

    }//check nbZ solutions
 }
 
 
 double w1X,w1Y,w1Z,w1E,w1M2;
 double H2X,H2Y,H2Z,H2E,H2M2;
 double qZ, qE;
 

 for(int i=0;i<2;i++){
  for(int j=0;j<2;j++){ 
   int solIdx=2*i+j;      
   mcfm_event_type& solevent=sol[solIdx];
   if(  jacobian[solIdx]<0){
    solevent.pswt=-1;
    continue;
   }
   
   double nuZ=roots_nuZ[i];
   double nbZ=roots_nbZ[i][j];
   double nuE=TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
   double nbE=TMath::Sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ);

   w1X = epX+nuX;
   w1Y = epY+nuY;
   w1Z = epZ+nuZ;
   w1E = epE+nuE;
   w1M2= w1E*w1E-w1X*w1X-w1Y*w1Y-w1Z*w1Z;

   H2X = emX+nbX+w1X;
   H2Y = emY+nbY+w1Y;
   H2Z = emZ+nbZ+w1Z;
   H2E = emE+nbE+w1E;
   H2M2= H2E*H2E-H2X*H2X-H2Y*H2Y-H2Z*H2Z;


   if(fabs(w1M2-Mw1_sqr)/Mw1_sqr>1E-4 ||  fabs(H2M2-MHiggs_sqr)/MHiggs_sqr>1E-4 ){
     solevent.pswt=-1;
     continue;
   }
 

   qZ=H2Z;
   qE=H2E;

   if(qE >2*EBEAM) {
        solevent.pswt=-1;
        continue;
   }
   
   double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
   double x0 = (qEold+qZ)/2/EBEAM;
   double x1 = (qEold-qZ)/2/EBEAM;

   double scale0=x0/(x0+x1);
   double scale1=x1/(x0+x1);

   if(x0> 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     solevent.pswt=-1;
     continue;
   }

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);
   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.p[4].SetPxPyPzE   (emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   (nbX, nbY, nbZ, nbE);

   //TVector3 boostV(H2X/H2E,H2Y/H2E,0);
   //for(int ipt=0;ipt<6;ipt++)solevent.p[ipt].Boost(-boostV);
   
   //w1E=solevent.p[2].Energy()+solevent.p[3].Energy();
   //w1Z=solevent.p[2].Pz()    +solevent.p[3].Pz();
   //nuZ=solevent.p[2].Pz();
   //nuE=solevent.p[2].Energy();
   //H2E=solevent.p[4].Energy()+solevent.p[5].Energy();
   //H2Z=solevent.p[4].Pz()    +solevent.p[5].Pz();
   //nbZ=solevent.p[5].Pz();
   //nbE=solevent.p[5].Energy();
   
   double J3nuZ=2*w1E*nuZ/nuE-2*w1Z;
   double J4nbZ=2*H2E*nbZ/nbE-2*H2Z;

   solevent.pswt=TMath::Abs(J3nuZ*J4nbZ);

  }//loop for nbZ solution
 }//loop for nuZ solution
		    
}
//=======================================================
// qX, qY, MHiggs, nuX, nuY, nuZ 
//
void WWL1L2Sol_MHiggs( cdf_event_type* temp,
                    double qX, double qY , 
		    double MHiggs_sqr, double nuX, double nuY, double nuZ,
                    mcfm_event_type* sol){

 double epX=temp->p[0].Px() ;
 double epY=temp->p[0].Py() ;
 double epZ=temp->p[0].Pz() ;
 double epE=temp->p[0].Energy();

 double emX=temp->p[1].Px() ;
 double emY=temp->p[1].Py() ;
 double emZ=temp->p[1].Pz() ;
 double emE=temp->p[1].Energy();

 double nbX=qX-epX-emX-nuX; 
 double nbY=qY-epY-emY-nuY; 
 double nuE=TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);

 double sumE=epE+emE+nuE;
 double sumZ=epZ+emZ+nuZ;
 double nbT_sqr = nbX*nbX+nbY*nbY;

 double A=MHiggs_sqr
         +(sumZ+sumE)*(sumZ-sumE)
         -nbT_sqr
         +qX*qX+qY*qY;
	 
 double k[10],roots[10];
 k[2]=(sumZ+sumE)*(sumZ-sumE);
 k[1]=A*sumZ;
 k[0]=A*A/4-sumE*sumE*nbT_sqr;

 double B2_4AC=k[1]*k[1]-4*k[2]*k[0];
 if(B2_4AC<0){
    sol[0].pswt=-2;
    sol[1].pswt=-2; 

/*
 double limits[2]={-EBEAM,EBEAM};
 PolySolver::FindRoots(2,k,limits,roots);

   double eval =k[0]
               +k[1]*roots[0]
               +k[2]*roots[0]*roots[0];
   if(fabs(eval/k[2])<1.0){
     cout << " MHiggs_sqr "<< MHiggs_sqr<< " sumE "<< sumE <<" sumZ "<< sumZ <<" nbT_sqr "<< nbT_sqr <<" qX "
     << qX <<" qY "<< qY<<endl;
     cout << "  A " << A <<" a b c " << k[2]<<" "<<k[1] <<" "<<k[0]<<endl;
     cout << " criteria "<< (sumZ+sumE)*(sumZ-sumE)*nbT_sqr+A*A <<endl;
     cout << " PolySolver Solution " << roots[0] <<" "<<roots[1]<<endl;
     cout << " nu " << nuX <<" "<<nuY<<" "<<nuZ<<endl;
     cout << " nb " << nbX <<" "<<nbY<<" "<<roots[0]<<endl;
     cout << " nb " << nbX <<" "<<nbY<<" "<<roots[1]<<endl;
   }
*/
    return;
 }
 else{
   B2_4AC=TMath::Sqrt(B2_4AC);
   roots[0]=(-k[1]+B2_4AC)/2/k[2];
   roots[1]=(-k[1]-B2_4AC)/2/k[2];
 }
 for(int i=0;i<2;i++){
   mcfm_event_type& solevent=sol[i];
   double nbZ=roots[i];
   double nbE=TMath::Sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ);
   double qE = sumE+nbE;
   double qZ = sumZ+nbZ;
  
   double Msys_sqr = qE*qE-qX*qX-qY*qY-qZ*qZ;
   if(fabs(Msys_sqr-MHiggs_sqr)/MHiggs_sqr > 1E-6) {
     solevent.pswt=-1;
     continue;
   }

   double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
   double x0 = (qEold+qZ)/2/EBEAM;
   double x1 = (qEold-qZ)/2/EBEAM;
   if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     solevent.pswt=-3;
     continue;
   }

   double scale0=x0/(x0+x1);
   double scale1=x1/(x0+x1);

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);
   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.p[4].SetPxPyPzE   (emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   (nbX, nbY, nbZ, nbE);

   //double Higgs_nbE =solevent.p[5].Energy();
   //double Higgs_nbZ =solevent.p[5].Pz();
   double Higgs_E=solevent.p[0].Energy()+solevent.p[1].Energy();
   double Higgs_Z=solevent.p[0].Pz()    +solevent.p[1].Pz();
   //solevent.pswt=TMath::Abs(2*Higgs_E*Higgs_nbZ/Higgs_nbE-2*Higgs_Z);
   solevent.pswt=2*Higgs_E*TMath::Abs(nbZ/nbE-Higgs_Z/Higgs_E);
 } 
}



//=======================================================
// qX, qY, MHiggs, YHiggs, nuX, nuY
//
void WWL1L2Sol_MHiggsYHiggs( cdf_event_type* temp,
                    double qX, double qY , 
		    double MHiggs, double YHiggs, double nuX, double nuY,
                    mcfm_event_type* sol){

 double epX=temp->p[0].Px() ;
 double epY=temp->p[0].Py() ;
 double epZ=temp->p[0].Pz() ;
 double epE=temp->p[0].Energy();

 double emX=temp->p[1].Px() ;
 double emY=temp->p[1].Py() ;
 double emZ=temp->p[1].Pz() ;
 double emE=temp->p[1].Energy();

 double nbX=qX-epX-emX-nuX; 
 double nbY=qY-epY-emY-nuY; 
 double nbT_sqr = nbX*nbX+nbY*nbY;
 double nuT_sqr = nuX*nuX+nuY*nuY;

 double exppY=TMath::Exp(YHiggs);
 double expmY=1./exppY;
 double ZHiggs=MHiggs*(exppY-expmY)/2;
 double qEold =MHiggs*(exppY+expmY)/2;
 double EHiggs=TMath::Sqrt(qEold*qEold+qX*qX+qY*qY);
 

 double sumE=EHiggs-epE-emE;
 double sumZ=ZHiggs-epZ-emZ;
 double A=(sumE*sumE-sumZ*sumZ+nuT_sqr-nbT_sqr)*0.5;


 double k[10],roots[10];
 k[2]=(sumZ+sumE)*(sumZ-sumE);
 k[1]=2*sumZ*A;
 k[0]=A*A-sumE*sumE*nuT_sqr;

 //Solve two solutions for nuZ
 double limits[2]={-3500.0,3500.0};
 PolySolver::FindRoots(2,k,limits,roots);


 //Loop Two Solutions
 for(int i=0;i<2;i++){
   //nbZ could be calculated from  nuZ
   mcfm_event_type& solevent=sol[i];
   double nuZ= roots[i];
   double nbZ= sumZ-nuZ;
   double nuE= TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
   double nbE= TMath::Sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ);

   //Check Higgs mass to make sure solutions could sum up to input Higgs Mass
   //within 1E-6 precisions
   double testE=epE+emE+nuE+nbE;
   double testZ=epZ+emZ+nuZ+nbZ;
   double testM_sqr =  testE*testE-testZ*testZ-qX*qX-qY*qY;
   
   if(fabs(TMath::Sqrt(testM_sqr)-MHiggs)/MHiggs > 1E-6) {
     solevent.pswt=-1;
     continue;
   }

   double x0 = (qEold+ZHiggs)/2/EBEAM;
   double x1 = (qEold-ZHiggs)/2/EBEAM;
   if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     solevent.pswt=-3;
     continue;
   }

   double scale0=x0/(x0+x1);
   double scale1=x1/(x0+x1);

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*EHiggs);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*EHiggs);
   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.p[4].SetPxPyPzE   (emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   (nbX, nbY, nbZ, nbE);

  
   //double JH1=2*EHiggs*nuZ/nuE-2*ZHiggs;
   //double JH2=2*EHiggs*nbZ/nbE-2*ZHiggs;
   //double JY1=(nuZ/nuE-nuZ)/(EHiggs+ZHiggs);
   //double JY2=(nbZ/nbE-nbZ)/(EHiggs+ZHiggs);
   //solevent.pswt=TMath::Abs(JH1*JY2-JH2*JY1);
   solevent.pswt=TMath::Abs(nuZ/nuE-nbZ/nbE);
 }
}
//=============================================================
//
//  Solve Kinematics for W(+)W(-)
//
//==============================================================


void WWL1L2Sol_Mw1Mw2( cdf_event_type* temp,
                       double  qX, double  qY,
                       double Mw1_sqr, double Mw2_sqr, double nuZ, double nbZ,
		       mcfm_event_type* sol){


 double epX=temp->p[0].Px() ;
 double epY=temp->p[0].Py() ;
 double epZ=temp->p[0].Pz() ;
 double epE=temp->p[0].Energy();
 
 double emX=temp->p[1].Px() ;
 double emY=temp->p[1].Py() ;
 double emZ=temp->p[1].Pz() ;
 double emE=temp->p[1].Energy();

 double A=Mw1_sqr+2*epZ*nuZ;

 double p1=  epX*epX-epE*epE;
 double p2=2*epX*epY;
 double p3=  epY*epY-epE*epE;
 double p4=epX*A;
 double p5=epY*A;
 double p6=A*A/4-epE*epE*nuZ*nuZ;

/* 
 double A=Mw1_sqr+epX*epX+epY*epY+epZ*epZ-epE*epE+2*epZ*nuZ;
 double B=2*epX;
 double C=2*epY;


 double p1=B*B-4*epE*epE;
 double p2=2*B*C;
 double p3=C*C-4*epE*epE;
 double p4=2*A*B;
 double p5=2*A*C;
 double p6=A*A-4*epE*epE*nuZ*nuZ;
*/ 
 double Qx=qX-epX-emX;
 double Qy=qY-epY-emY;

 double D= Mw2_sqr/2+(emX*Qx+emY*Qy+emZ*nbZ);

 double q1=  emX*emX-emE*emE;
 double q2=2*emX*emY;
 double q3=  emY*emY-emE*emE;
 double q4=-2*emX*D + 2*Qx*emE*emE;
 double q5=-2*emY*D + 2*Qy*emE*emE;
 double q6=D*D -emE*emE*(Qx*Qx+Qy*Qy+nbZ*nbZ) ;

/* 
 double D=Mw2_sqr-emE*emE
         +epX*epX-Qx*Qx
	 +epY*epY-Qy*Qy
	 +emZ*emZ+2*emZ*nbZ;
	 
 double E=-2*emX;
 double F=-2*emY;

 double q1=E*E-4*emE*emE;
 double q2=2*E*F;
 double q3=F*F-4*emE*emE;
 double q4=2*D*E+8*Qx*emE*emE;
 double q5=2*D*F+8*Qy*emE*emE;
 double q6=D*D-4*emE*emE*nbZ*nbZ-4*emE*emE*(Qx*Qx+Qy*Qy);
*/ 
// p1=1;p2=0;p3=1;p4=-2;p5=-2;p6=2;
// q1=1;q2=0;q3=1;q4=-2;q5=-2;q6=2;

 double M= q3*p1-p3*q1;
 double N= q3*p4-p3*q4;
 double O= q3*p6-p3*q6;
 double P=-q3*p2+p3*q2;
 double Q=-q3*p5+p3*q5;
 
 double R=     P*p2+M*p3;
 double S=P*p5+Q*p2+N*p3;
 double T=Q*p5     +O*p3;
       
 double k[5];
 k[4]=  	       P*P*p1+M*R;
 k[3]=        2*P*Q*p1+P*P*p4+M*S+N*R;
 k[2]= Q*Q*p1+2*P*Q*p4+P*P*p6+M*T+N*S+O*R;
 k[1]= Q*Q*p4+2*P*Q*p6  	 +N*T+O*S;
 k[0]= Q*Q*p6			     +O*T;
 
 double limits[2]={-EBEAM,EBEAM};
 double roots [10];
 PolySolver::FindRoots(4,k,limits,roots);
 
 double nuX,nuY,nuE;
 double nbX,nbY,nbE;
 double w1X,w1Y,w1Z,w1E,w1M2; 
 double w2X,w2Y,w2Z,w2E,w2M2; 
 double qZ, qE;
 
 for(int i=0;i<4;i++){
   mcfm_event_type& solevent=sol[i];
   nuX=roots[i];
   if(nuX!=nuX){
     solevent.pswt=-1; 
     continue;
   }

   double eval =k[0]
               +k[1]*nuX
               +k[2]*nuX*nuX
               +k[3]*nuX*nuX*nuX
               +k[4]*nuX*nuX*nuX*nuX;

   if(fabs(eval/k[4])>1.0) {
      solevent.pswt=-1;      
      continue;
   }
   
   nuY=(M*nuX*nuX+N*nuX+O)/(P*nuX+Q);
   nuE=TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);

   nbX=qX-(epX+emX+nuX);
   nbY=qY-(epY+emY+nuY); 
   nbE=TMath::Sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ);

   w1X = epX+nuX;
   w1Y = epY+nuY;
   w1Z = epZ+nuZ;
   w1E = epE+nuE; 
   w1M2= w1E*w1E-w1X*w1X-w1Y*w1Y-w1Z*w1Z; 

   w2X = emX+nbX;
   w2Y = emY+nbY;
   w2Z = emZ+nbZ;
   w2E = emE+nbE;    
   w2M2= w2E*w2E-w2X*w2X-w2Y*w2Y-w2Z*w2Z;
			       

   if(fabs(w1M2-Mw1_sqr)/Mw1_sqr>1E-4 ||  fabs(w2M2-Mw2_sqr)/Mw2_sqr>1E-4 ){
     solevent.pswt=-1;
     continue;
   }
   
   qZ = w1Z+w2Z;
   qE = w1E+w2E;

   if(qE >2*EBEAM) {
        solevent.pswt=-1;
	continue;
   }

   double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
   double x0 = (qEold+qZ)/2/EBEAM;
   double x1 = (qEold-qZ)/2/EBEAM;     


   double scale0=x0/(x0+x1);
   double scale1=x1/(x0+x1);

   if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {

     solevent.pswt=-1;
     continue;
   }

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);

   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.p[4].SetPxPyPzE   (emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   (nbX, nbY, nbZ, nbE);

   //Jacobian
   
   double J3x=2*(epE*nuX/nuE-epX); 
   double J3y=2*(epE*nuY/nuE-epY); 
   double J4x=2*(emE*nbX/nbE-emX); 
   double J4y=2*(emE*nbY/nbE-emY); 
   solevent.pswt=TMath::Abs(J3y*J4x-J3x*J4y);

 /*
   //  
   //Lorentz Boost to PT=0 frame
   //
   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.p[4].SetPxPyPzE   (emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   (nbX, nbY, nbZ, nbE);

   if((qX*qX+qY*qY)>1.0E-8){
     TVector3 boostV(qX/qE,qY/qE,0);
     for(int ipt=2;ipt<6;ipt++)  solevent.p[ipt].Boost(-boostV);          
   }
   
   double qEold=0;
   for(int ipt=2;ipt<6;ipt++){    qEold+=solevent.p[ipt].Energy(); }
   
   double x0EBEAM = (qEold+qZ)/2;
   double x1EBEAM = (qEold-qZ)/2;     
   solevent.p[0].SetPxPyPzE(0.0,0.0, x0EBEAM,x0EBEAM);
   solevent.p[1].SetPxPyPzE(0.0,0.0,-x1EBEAM,x1EBEAM);
     
   double J3x=2*(solevent.p[3].Energy()*solevent.p[2].Px()/solevent.p[2].Energy()-solevent.p[3].Px()); 
   double J3y=2*(solevent.p[3].Energy()*solevent.p[2].Py()/solevent.p[2].Energy()-solevent.p[3].Py()); 
   double J4x=2*(solevent.p[4].Energy()*solevent.p[5].Px()/solevent.p[5].Energy()-solevent.p[4].Px()); 
   double J4y=2*(solevent.p[4].Energy()*solevent.p[5].Py()/solevent.p[5].Energy()-solevent.p[4].Py()); 
   solevent.pswt=TMath::Abs(J3y*J4x-J3x*J4y);
  */ 
 }

}

//=============================================================
//
//  Solve Kinematics for W(+/-) gamma/parton
//
//==============================================================
void WWL1L2Sol_Mw ( cdf_event_type* temp,
                    double  qX, double   qY, double Msqr,
                    mcfm_event_type* sol){
   
   
   double lepX = temp->p[0].Px();
   double lepY = temp->p[0].Py();
   double lepZ = temp->p[0].Pz();
   double lepE = temp->p[0].Energy();
   double gamX = temp->p[1].Px();
   double gamY = temp->p[1].Py();
   double gamZ = temp->p[1].Pz();
   double gamE = temp->p[1].Energy();
   if(nwz_.nwz==-1){
      	  lepX = temp->p[1].Px();
	  lepY = temp->p[1].Py();
	  lepZ = temp->p[1].Pz();
	  lepE = temp->p[1].Energy();
	  gamX = temp->p[0].Px();
	  gamY = temp->p[0].Py();
	  gamZ = temp->p[0].Pz();
	  gamE = temp->p[0].Energy();
   }
   double nuX  = qX-lepX-gamX;
   double nuY  = qY-lepY-gamY;
		
   double A = Msqr/2.+lepX*nuX+lepY*nuY;		
   double k[3];
   k[2]=(lepE+lepZ)*(lepE-lepZ);
   k[1]=-2.0*A*lepZ;
   k[0]=lepE*lepE*(nuX*nuX+nuY*nuY)-A*A;

//cout <<"Debug q " << qX <<" "<<qY <<" nu " << nuX <<" "<< nuY << " Msqr " << Msqr<<endl;
//cout <<"Debug lep1 "<< lepX <<" " << lepY <<" "<< lepZ <<" "<<lepE <<endl;
//cout <<"Debug gam1 "<< gamX <<" " << gamY <<" "<< gamZ <<" "<<gamE <<endl;

//cout <<"Debug eqn " << k[2] <<" "<<k[1] <<" "<< k[0] <<endl; 
   double limits[2]={-EBEAM,EBEAM};
   double roots [2];
   PolySolver::FindRoots(2,k,limits,roots);

   for(int i=0;i<2;i++){
     mcfm_event_type& solevent = sol[i];
     double nuZ= roots[i];

     double eval = k[0]
                  +k[1]*nuZ
                  +k[2]*nuZ*nuZ;

//cout <<"Debug MwSol " << nuZ <<" input " << Msqr <<" eval=" << eval <<endl;

     if(fabs(eval/k[2])>1.0){
       solevent.pswt=-1; 
       continue;
     }
     double nuE= TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
     double wX = lepX+nuX;
     double wY = lepY+nuY;
     double wZ = lepZ+nuZ;
     double wE = lepE+nuE;
     double wMsqr = wE*wE-wX*wX-wY*wY-wZ*wZ;

//cout <<"Debug MwSol " << wMsqr <<" input " << Msqr <<endl;
     if(fabs(wMsqr-Msqr)/Msqr<1E-4 >1){
       solevent.pswt=-1; 
       continue;

     }

     double qZ = wZ+gamZ;
     double qE = wE+gamE;

     if(qE >2*EBEAM) {
        solevent.pswt=-1;
	continue;
     }

     double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
     double x0 = (qEold+qZ)/2/EBEAM;
     double x1 = (qEold-qZ)/2/EBEAM;     


     double scale0=x0/(x0+x1);
     double scale1=x1/(x0+x1);

     if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     solevent.pswt=-1;
     continue;
   }

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);

   if(nwz_.nwz==1){
     solevent.p[2].SetPxPyPzE   ( nuX, nuY, nuZ, nuE);
     solevent.p[3].SetPxPyPzE   (lepX,lepY,lepZ,lepE);
     solevent.p[4].SetPxPyPzE   (gamX,gamY,gamZ,gamE);
   }
   else{  
     solevent.p[2].SetPxPyPzE   (lepX,lepY,lepZ,lepE);
     solevent.p[3].SetPxPyPzE   ( nuX, nuY, nuZ, nuE);
     solevent.p[4].SetPxPyPzE   (gamX,gamY,gamZ,gamE);
   }
   
     double J3x=2*(lepE*nuZ/nuE-lepZ); 
     solevent.pswt=TMath::Abs(J3x);
     
   
 }//two solutions

}

//=============================================================
//
//  Solve Kinematics for Z(l+l-) Z(nu nu~)
//
//==============================================================


void WWL1L2Sol_MzNu3( cdf_event_type* temp,
                    double  qX, double   qY,
                    double  Msqr, double  nuX, double nuY, double nuZ,
                    mcfm_event_type* sol){


 double epX=temp->p[0].Px() ;
 double epY=temp->p[0].Py() ;
 double epZ=temp->p[0].Pz() ;
 double epE=temp->p[0].Energy();
 
 double emX=temp->p[1].Px() ;
 double emY=temp->p[1].Py() ;
 double emZ=temp->p[1].Pz() ;
 double emE=temp->p[1].Energy();

		    
   double nbX = qX-epX-emX-nuX;
   double nbY = qY-epY-emY-nuY;
   
   double nuE = TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
   double A = Msqr/2.+nuX*nbX+nuY*nbY;
   
 double k[3];
 k[2]= (nuE+nuZ)*(nuE-nuZ);
 k[1]= -2*A*nuZ;
 k[0]= nuE*nuE*(nbX*nbX+nbY*nbY)-A*A;
 
 double limits[2]={-EBEAM,EBEAM};
 double roots [2];
 PolySolver::FindRoots(2,k,limits,roots);

 
 for(int i=0;i<2;i++){
   mcfm_event_type& solevent=sol[i];
   double nbZ=roots[i];
 
   double eval =k[0]
               +k[1]*nbZ
               +k[2]*nbZ*nbZ;

   if(fabs(eval/k[2])>1.0){
     solevent.pswt=-1;
     continue;
   }

   double nbE=TMath::Sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ);
   
   double zX = nuX+nbX;
   double zY = nuY+nbY;
   double zZ = nuZ+nbZ;
   double zE = nuE+nbE;
   
   double zMsqr = zE*zE-zX*zX-zY*zY-zZ*zZ;
   
   if(fabs(zMsqr-Msqr)/Msqr>1E-4){
     solevent.pswt=-1;
     continue;   
   }


     double qZ = epZ+emZ+zZ;
     double qE = epE+emE+zE;

     if(qE >2*EBEAM) {
        solevent.pswt=-1;
	continue;
     }

     double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
     double x0 = (qEold+qZ)/2/EBEAM;
     double x1 = (qEold-qZ)/2/EBEAM;     


     double scale0=x0/(x0+x1);
     double scale1=x1/(x0+x1);

     if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {

     solevent.pswt=-1;
     continue;
   }
   //82 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->e^-(p5)+e^+(p6))'

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);
  
   solevent.p[2].SetPxPyPzE   ( nuX, nuY, nuZ, nuE);
   solevent.p[3].SetPxPyPzE   ( nbX, nbY, nbZ, nbE);
   solevent.p[4].SetPxPyPzE   ( emX, emY, emZ, emE);
   solevent.p[5].SetPxPyPzE   ( epX, epY, epZ, epE);
     
   double J3x=2*(nuE*nbZ/nbE-nuZ); 
   solevent.pswt=TMath::Abs(J3x);
     
 }//solutions
  		    
}




//########################
void WWL1L2Sol_Mw1( event_type* temp
                    ,double Mw1, double nu_X, double nu_Y, double nu_Z
                    ,event_type* sol,double* jacobian){
		    
   double nb_X = -(temp->ep.Px()+temp->em.Px()+nu_X);		    
   double nb_Y = -(temp->ep.Py()+temp->em.Py()+nu_Y);		    

   double nu_E = TMath::Sqrt(nu_X*nu_X+nu_Y*nu_Y+nu_Z*nu_Z);
   double A = (nu_E-nu_Z)*(nu_Z+nu_Z);
   double nb_T = TMath::Sqrt(nb_X*nb_X+nb_Y*nb_Y);
   
 double k[3];
 k[2]= A;
 k[1]= -2*A*nu_Z;
 k[0]= (nu_E*nb_T-A)*(nu_E*nb_T+A);
 
 double limits[2]={-EBEAM,EBEAM};
 double roots [2];
 PolySolver::FindRoots(2,k,limits,roots);

 
 for(int i=0;i<2;i++){
   event_type& solevent=sol[i];
   double nb_Z=roots[i];
 
   solevent.nu.SetXYZM(nu_X, nu_Y, nu_Z,0);
   solevent.nb.SetXYZM(nb_X, nb_Y, nb_Z,0);

   double nu_E=solevent.nu.E();
   double nb_E=solevent.nb.E();
   
   TLorentzVector p1=solevent.nu+solevent.nb;


   double xx=roots[i];
   double eval =k[0]
        +k[1]*xx
        +k[2]*xx*xx;


   if(fabs(p1.M()-Mw1)/Mw1<1E-4 &&   fabs(eval/k[2])<1.0){


     solevent.ep.SetPxPyPzE(temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy());
     solevent.em.SetPxPyPzE(temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy());
     
     double J3x=2*(nu_E*nb_Z/nb_E-nb_Z); 
     jacobian[i]=TMath::Abs(J3x);
     
   }
   else{
     jacobian[i]=-1;
   }
 }
   		    
}

void WWL1L2Sol_Mw1Mw2( event_type* temp
                     , double Mw1, double Mw2, double nuZ, double nbZ
		     , event_type* sol, double* jacobian){

 double epX=temp->ep.Px();double  epY=temp->ep.Py();double  epZ=temp->ep.Pz();double  epE=temp->ep.Energy();
 double emX=temp->em.Px();double  emY=temp->em.Py();double  emZ=temp->em.Pz();double  emE=temp->em.Energy();    
 
 double A=Mw1*Mw1+epX*epX+epY*epY+epZ*epZ-epE*epE+2*epZ*nuZ;
 double B=2*epX;
 double C=2*epY;

 double p1=B*B-4*epE*epE;
 double p2=2*B*C;
 double p3=C*C-4*epE*epE;
 double p4=2*A*B;
 double p5=2*A*C;
 double p6=A*A-4*epE*epE*nuZ*nuZ;
 
 double Qx=-epX-emX;
 double Qy=-epY-emY;
 
 double D=Mw2*Mw2-emE*emE
         +epX*epX-Qx*Qx
	 +epY*epY-Qy*Qy
	 +emZ*emZ+2*emZ*nbZ;
	 
 double E=-2*emX;
 double F=-2*emY;

 double q1=E*E-4*emE*emE;
 double q2=2*E*F;
 double q3=F*F-4*emE*emE;
 double q4=2*D*E+8*Qx*emE*emE;
 double q5=2*D*F+8*Qy*emE*emE;
 double q6=D*D-4*emE*emE*nbZ*nbZ-4*emE*emE*(Qx*Qx+Qy*Qy);
 
// p1=1;p2=0;p3=1;p4=-2;p5=-2;p6=2;
// q1=1;q2=0;q3=1;q4=-2;q5=-2;q6=2;

 double M= q3*p1-p3*q1;
 double N= q3*p4-p3*q4;
 double O= q3*p6-p3*q6;
 double P=-q3*p2+p3*q2;
 double Q=-q3*p5+p3*q5;
 
 double R=     P*p2+M*p3;
 double S=P*p5+Q*p2+N*p3;
 double T=Q*p5     +O*p3;
       
 double k[5];
 k[4]=  	       P*P*p1+M*R;
 k[3]=        2*P*Q*p1+P*P*p4+M*S+N*R;
 k[2]= Q*Q*p1+2*P*Q*p4+P*P*p6+M*T+N*S+O*R;
 k[1]= Q*Q*p4+2*P*Q*p6  	 +N*T+O*S;
 k[0]= Q*Q*p6			     +O*T;
 
 double limits[2]={-EBEAM,EBEAM};
 double roots [10];
 PolySolver::FindRoots(4,k,limits,roots);
 
 double nuX,nuY,nuE;
 double nbX,nbY,nbE;
  
 
 for(int i=0;i<4;i++){
   event_type& solevent=sol[i];
   nuX=roots[i];
   nuY=(M*nuX*nuX+N*nuX+O)/(P*nuX+Q);
   nbX=-(epX+emX+nuX);
   nbY=-(epY+emY+nuY);
 
   solevent.nu.SetXYZM(nuX, nuY, nuZ,0);
   solevent.nb.SetXYZM(nbX, nbY, nbZ,0);

   nuE=solevent.nu.E();
   nbE=solevent.nb.E();
   
   TLorentzVector p1=temp->ep+solevent.nu;
   TLorentzVector p2=temp->em+solevent.nb;


   double xx=roots[i];
   double eval =k[0]
        +k[1]*xx
        +k[2]*xx*xx
        +k[3]*xx*xx*xx
        +k[4]*xx*xx*xx*xx;


   if(fabs(p1.M()-Mw1)/Mw1<1E-4 &&  fabs(p2.M()-Mw2)/Mw2<1E-4 && fabs(eval/k[4])<1.0){

     solevent.ep.SetPxPyPzE(temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy());
     solevent.em.SetPxPyPzE(temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy());
     
     double J3x=2*(epE*nuX/nuE-epX); 
     double J3y=2*(epE*nuY/nuE-epY); 
     double J4x=2*(emE*nbX/nbE-emX); 
     double J4y=2*(emE*nbY/nbE-emY); 
     jacobian[i]=TMath::Abs(J3y*J4x-J3x*J4y);
     
   }
   else{
     jacobian[i]=-1;
/*
     if(DEBUGFLAG && fabs(p1.M()-Mw1)/Mw1<1E-4 &&  fabs(p2.M()-Mw2)/Mw2<1E-4){
        cout <<"= Bad Solution Begin ====================="<<endl;
        printf("ep: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->ep.M(),temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy() );
	printf("em: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->em.M(),temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy() );
	cout << "Mwp=" << Mw1<<" Mwm="<< Mw2 <<" nuPz="<< nuZ<<" nbPz="<< nbZ<<endl;
      	cout << "sol nuX=" << xx <<" eval/k[4]=" << eval/k[4] << " Mw1=" << p1.M() << " Mw2=" << p2.M() <<endl;
	cout <<"= Bad Solution End ==================="<<endl;
     }
*/
   }
   
 }

}

void WWL1L2Sol_Mw1Mw2_exact( event_type* temp
                     , double Mw1, double Mw2, double nuZ, double nbZ
		     , event_type* sol, double* jacobian){

 double epX=temp->ep.Px();double  epY=temp->ep.Py();double  epZ=temp->ep.Pz();double  epE=temp->ep.Energy();
 double emX=temp->em.Px();double  emY=temp->em.Py();double  emZ=temp->em.Pz();double  emE=temp->em.Energy();    
 double scale=1.;
 Mw1=Mw1/scale;Mw2=Mw2/scale;
 nuZ=nuZ/scale;nbZ=nbZ/scale;
 epX=epX/scale;epY=epY/scale;epZ=epZ/scale;epE=epE/scale;
 emX=emX/scale;emY=emY/scale;emZ=emZ/scale;emE=emE/scale;
  
 double A=Mw1*Mw1+epX*epX+epY*epY+epZ*epZ-epE*epE+2*epZ*nuZ;
 double B=2*epX;
 double C=2*epY;

 double p1=B*B-4*epE*epE;
 double p2=2*B*C;
 double p3=C*C-4*epE*epE;
 double p4=2*A*B;
 double p5=2*A*C;
 double p6=A*A-4*epE*epE*nuZ*nuZ;
 
 double Qx=-epX-emX;
 double Qy=-epY-emY;
 
 double D=Mw2*Mw2-emE*emE
         +epX*epX-Qx*Qx
	 +epY*epY-Qy*Qy
	 +emZ*emZ+2*emZ*nbZ;
	 
 double E=-2*emX;
 double F=-2*emY;

 double q1=E*E-4*emE*emE;
 double q2=2*E*F;
 double q3=F*F-4*emE*emE;
 double q4=2*D*E+8*Qx*emE*emE;
 double q5=2*D*F+8*Qy*emE*emE;
 double q6=D*D-4*emE*emE*nbZ*nbZ-4*emE*emE*(Qx*Qx+Qy*Qy);
 
// p1=1;p2=0;p3=1;p4=-2;p5=-2;p6=2;
// q1=1;q2=0;q3=1;q4=-2;q5=-2;q6=2;

 double M= q3*p1-p3*q1;
 double N= q3*p4-p3*q4;
 double O= q3*p6-p3*q6;
 double P=-q3*p2+p3*q2;
 double Q=-q3*p5+p3*q5;
 
 double R=     P*p2+M*p3;
 double S=P*p5+Q*p2+N*p3;
 double T=Q*p5     +O*p3;
       
 double k[5];
 k[4]=  	       P*P*p1+M*R;
 k[3]=        2*P*Q*p1+P*P*p4+M*S+N*R;
 k[2]= Q*Q*p1+2*P*Q*p4+P*P*p6+M*T+N*S+O*R;
 k[1]= Q*Q*p4+2*P*Q*p6  	 +N*T+O*S;
 k[0]= Q*Q*p6			     +O*T;


double nuX1=temp->nu.Px()/scale;
double nuY1=temp->nu.Py()/scale;
//double nuE1=sqrt(temp->nu.Px()*temp->nu.Px()+temp->nu.Py()*temp->nu.Py()+temp->nu.Pz()*temp->nu.Pz());
double nuE1=sqrt(nuX1*nuX1+nuY1*nuY1+nuZ*nuZ);


	
double checkEq1=(epX+ nuX1)*(epX+ nuX1)
               +(epY+ nuY1)*(epY+ nuY1)
	       +(epZ+ nuZ)*(epZ+ nuZ)
	       -(epE+nuE1)*(epE+nuE1)
	       +Mw1*Mw1;
double checkEq2=A + B*nuX1 + C*nuY1 - 2*epE*nuE1;

double checkEq3=p1*nuX1*nuX1 
              + p2*nuX1*nuY1
	      + p3*nuY1*nuY1 
	      + p4*nuX1 
	      + p5*nuY1 
	      + p6;

double checkEq4= A*A + B*nuX1*B*nuX1 + C*nuY1*C*nuY1 +2*epE*nuE1*2*epE*nuE1
               +A*(B*nuX1+C*nuY1   -2*epE*nuE1)
	       +B*nuX1*(A+C*nuY1   -2*epE*nuE1)
	       +C*nuY1*(A+B*nuX1   -2*epE*nuE1)
                                   -2*epE*nuE1*(A+B*nuX1+C*nuY1);

double checkEq5=(A+B*nuX1+C*nuY1)*(A+B*nuX1+C*nuY1)-2*epE*nuE1*2*epE*nuE1;

double checkEq6=A*A+ B*nuX1*B*nuX1 + C*nuY1*C*nuY1
               +2*A*B*nuX1 + 2*B*C*nuX1*nuY1 + 2*A*C*nuY1
	       -2*epE*nuE1*2*epE*nuE1;

double checkEq7=(A + B*nuX1 + C*nuY1)*(A + B*nuX1 + C*nuY1)
               - 2* 2*epE*nuE1 *(A + B*nuX1 + C*nuY1)
	       + 2*epE*nuE1*2*epE*nuE1;

double checkEq8= A*A + B*B*nuX1*nuX1 + C*C*nuY1*nuY1
               + 2*A*B*nuX1 + 2*A*C*nuY1 + 2*B*C*nuX1*nuY1
	       - 4*epE*epE * nuX1*nuX1 
	       - 4*epE*epE * nuY1*nuY1
	       - 4*epE*epE * nuZ*nuZ;

double checkEq9= (B*B-4*epE*epE)*nuX1*nuX1
               +  2*B*C         *nuX1*nuY1      
               + (C*C-4*epE*epE)*nuY1*nuY1
	       +  2*A*B*nuX1 
	       +  2*A*C*nuY1
	       + A*A - 4*epE*epE * nuZ*nuZ;
	       
cout <<"Debug: Check1 =" << checkEq1 <<" check2=" << checkEq2 
          <<" checkEq3=" << checkEq3 <<" check4=" << checkEq4 
	  <<" checkEq5=" << checkEq5 <<" check6=" << checkEq6
	  <<" checkEq7=" << checkEq7 <<" check8=" << checkEq8
	  <<" checkEq9=" << checkEq9 
	  << endl;

double nbX1=-(epX+emY+nuX1);
double nbY1=-(epY+emY+nuY1);
double nbE1=sqrt(nbX1*nbX1+nbY1*nbY1+nbZ*nbZ);

nbX1 = temp->nb.Px();
nbY1 = temp->nb.Py();
nbX1 = -1*(epX+emX+nuX1);
nbY1 = -1*(epY+emY+nuY1);
nbE1=sqrt(nbX1*nbX1+nbY1*nbY1+nbZ*nbZ);

checkEq2=(emX +nbX1)*(emX +nbX1)
               +(emY +nbY1)*(emY +nbY1)
	       +(emZ +nbZ )*(emZ +nbZ )
	       -(emE +nbE1)*(emE +nbE1)
	       +Mw2*Mw2;

cout <<"Debug: PtX=" << epX+emX+nuX1+nbX1 << " Pty=" << epY+emY+nuY1+nbY1 <<endl;

cout <<"Debug: Check2 =" << checkEq2 <<endl;
double xx=temp->nu.Px()/scale;
double yy=temp->nu.Py()/scale;
double evalp=p1*xx*xx + p2*xx*yy + p3*yy*yy + p4*xx + p5*yy + p6;
double evalq=q1*xx*xx + q2*xx*yy + q3*yy*yy + q4*xx + q5*yy + q6;
cout << "Debug1:x="<< xx<<" p="<<p1 <<" "<<p2 <<" "<<p3<<" "<<p4 <<" " << p5 <<" "<<p6 <<" => " << evalp << endl;
cout << "Debug1:x="<< xx<<" q="<<q1 <<" "<<q2 <<" "<<q3<<" "<<q4 <<" " << q5 <<" "<<q6 <<" => " << evalq << endl;
/*
k[0]=k[0]/k[4];
k[1]=k[1]/k[4];
k[2]=k[2]/k[4];
k[3]=k[3]/k[4];
k[4]=k[4]/k[4];
*/
double eval=k[0]
	   +k[1]*xx
	   +k[2]*xx*xx
           +k[3]*xx*xx*xx
	   +k[4]*xx*xx*xx*xx;

double eval1=k[0]/k[4]
            +k[1]/k[4]*xx
            +k[2]/k[4]*xx*xx
            +k[3]/k[4]*xx*xx*xx
            +          xx*xx*xx*xx;
cout << "M=" << M << " N=" << N << " O=" << O << " P=" << P<< " Q=" << Q 
     << " R=" << R <<" S=" << S << " T=" << T <<endl;
//cout << "Debug: A="<<A <<" B="<<B <<" C="<<C <<endl;

cout << "Debug2: Sol x= " << xx 
     << " coef: " << k[0] <<" "<<k[1]<<" "<<k[2] <<" "<<k[3]<<" "<<k[4]	
     << " eqn=" <<eval<<endl;

cout << "Debug2: Sol x= " << xx    
     << " coef1: " << k[0]/k[4] <<" "<<k[1]/k[4]<<" "<<k[2]/k[4] <<" "<<k[3]/k[4]
     << " eqn1=" <<eval1<<endl;
     
 double limits[2]={-EBEAM,EBEAM};
 double roots [10];
 //PolySolver::FindRoots4thOrder(k,roots,flag);
 //PolySolver::FindRoots4thOrder(k,roots,flag);
 PolySolver::FindRoots(4,k,limits,roots); 
 double nuX,nuY,nuE;
 double nbX,nbY,nbE;
  
 
 for(int i=0;i<4;i++){
   xx=roots[i];	 
          eval =k[0]
               +k[1]*xx
	       +k[2]*xx*xx
	       +k[3]*xx*xx*xx
	       +k[4]*xx*xx*xx*xx;
	       
   cout << "Debug: Sol "<< i <<" x= " << xx <<" eqn=" <<eval<<endl;
          eval1=k[0]/k[4]
               +k[1]/k[4]*xx
	       +k[2]/k[4]*xx*xx
	       +k[3]/k[4]*xx*xx*xx
	       +          xx*xx*xx*xx;
   cout << "Debug: Sol "<< i <<" x= " << xx <<" eqn1=" <<eval1<<endl;
   
   event_type& solevent=sol[i];
   nuX=roots[i];
   nuY=(M*nuX*nuX+N*nuX+O)/(P*nuX+Q);
   nbX=-(epX+emX+nuX);
   nbY=-(epY+emY+nuY);
 
   solevent.nu.SetXYZM(nuX, nuY, nuZ,0);
   solevent.nb.SetXYZM(nbX, nbY, nbZ,0);

   nuE=solevent.nu.E();
   nbE=solevent.nb.E();
   
   TLorentzVector p1=temp->ep+solevent.nu;
   TLorentzVector p2=temp->em+solevent.nb;
   
   if(fabs(p1.M()-Mw1)/Mw1<1E-4 &&  fabs(p2.M()-Mw2)/Mw2<1E-4 && fabs(eval/k[4])<1E-1){
	   
     solevent.ep.SetPxPyPzE(temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy());
     solevent.em.SetPxPyPzE(temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy());
     
     double J3x=2*(epE*nuX/nuE-epX); 
     double J3y=2*(epE*nuY/nuE-epY); 
     double J4x=2*(emE*nbX/nbE-emX); 
     double J4y=2*(emE*nbY/nbE-emY); 
     jacobian[i]=TMath::Abs(J3y*J4x-J3x*J4y);
     
   }
   else{
     jacobian[i]=-1;
   }
   
 }

}


void WWL1L2Sol(double theta,double phi,event_type* temp, double* J){

    double costh = TMath::Cos(theta);
    double sinth = TMath::Sin(theta);
    double cosphi= TMath::Cos(phi);
    double sinphi= TMath::Sin(phi);

   TLorentzVector Pww= -1*(temp->p1+temp->p2);

   
   TLorentzVector PwwMinus= Pww - temp->ep - temp->em;
   double A = PwwMinus.P()*PwwMinus.P();
   double B =(PwwMinus.Px()*sinth*cosphi
            + PwwMinus.Py()*sinth*sinphi
  	    + PwwMinus.Pz()*costh)*2;
   double C = Pww.Energy() - temp->ep.Energy() - temp->em.Energy();
     
   double P4=(C*C-A)/(2*C-B);      
   *J = TMath::Abs(1-(B-2*P4)/2/TMath::Sqrt(A-B*P4+P4*P4));

//  double B = (PwwMinus.Px()*sinth*cosphi
//            + PwwMinus.Py()*sinth*sinphi
// 	    + PwwMinus.Pz()*costh);
    
    
//  double P4 = PwwMinus.M()*PwwMinus.M()/2/(PwwMinus.Energy()-B);
//  *J = 1- (B-P4)/TMath::Sqrt(PwwMinus.P()*PwwMinus.P()-2*B*P4+P4*P4);
  
   temp->nu.SetPxPyPzE(P4*sinth*cosphi,P4*sinth*sinphi,P4*costh,P4);
   temp->nb = Pww - temp->ep - temp->em - temp->nu;


/*
   cout << "SCHSU debug: "
        << " A B C:" << A <<" " <<B <<" "<< C
	<< " Beta Kppa B**2-2K:" << Beta<<" "<<Kppa <<" "<< (Beta*Beta-4*Kppa)	
	<< " P4:" <<P4
	<< " nu nb: " << temp->nu.Energy()<<" "<<temp->nb.Energy()
	<<endl;       
*/
}
void WWL1L2Sol_MHiggs( event_type* temp
                     ,double qX, double qY
                     ,double MHiggs, double nuX, double nuY, double nuZ
                     ,event_type* sol,double* jacobian){
		     
 double epX=temp->ep.Px();double  epY=temp->ep.Py();double  epZ=temp->ep.Pz();double  epE=temp->ep.Energy();
 double emX=temp->em.Px();double  emY=temp->em.Py();double  emZ=temp->em.Pz();double  emE=temp->em.Energy();    
 double nbX=qX-epX-emX-nuX; double  nbY=qY-epY-emY-nuY; 
 double nuE=TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
 
 double sumE=(epE+emE+nuE);
 double sumZ=(epZ+emZ+nuZ);
 double nbT = TMath::Sqrt(nbX*nbX+nbY*nbY);
 
 double A=MHiggs*MHiggs
         +sumZ*sumZ
	 -sumE*sumE
	 -nbT*nbT
         +qX*qX+qY*qY;

 double criteria = A*A+(sumZ+sumE)*(sumZ-sumE)*nbT*nbT;
 if(criteria<0){
   jacobian[0]=-2;
   jacobian[1]=-2;
   return;
 }

 double k[10];
 k[2]=(sumZ+sumE)*(sumZ-sumE);
 k[1]=A*sumZ;
 k[0]=(A/2-sumE*nbT)*(A/2+sumE*nbT);
 
 
 double limits[2]={-EBEAM,EBEAM};
 double roots [10];
 PolySolver::FindRoots(2,k,limits,roots);
 
 double nbZ,nbE;  
 
 for(int i=0;i<2;i++){
   event_type& solevent=sol[i];
   nbZ=roots[i];
 
   solevent.nu.SetXYZM(nuX, nuY, nuZ,0);
   solevent.nb.SetXYZM(nbX, nbY, nbZ,0);

   nuE=solevent.nu.E();
   nbE=solevent.nb.E();
   
   TLorentzVector p_Higgs=temp->ep+solevent.nu
                         +temp->em+solevent.nb;
  

   double eval =k[0]
               +k[1]*nbZ
               +k[2]*nbZ*nbZ;


   if(fabs(p_Higgs.M()-MHiggs)/MHiggs<1E-4 &&   fabs(eval/k[2])<1.0){

   double qE = sumE+nbE;
   double qZ = sumZ+nbZ;

   double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
   double x0 = (qEold+qZ)/2/EBEAM;
   double x1 = (qEold-qZ)/2/EBEAM;
   if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     jacobian[i]=-1;
     continue;
   }

//   double scale0=x0/(x0+x1);
//   double scale1=x1/(x0+x1);
//   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
//   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);
//   solevent.p[2].SetPxPyPzE   (nuX, nuY, nuZ, nuE);
   solevent.ep.SetPxPyPzE   (epX, epY, epZ, epE);
   solevent.em.SetPxPyPzE   (emX, emY, emZ, emE);

//   solevent..SetPxPyPzE   (nbX, nbY, nbZ, nbE);

//     solevent.ep.SetPxPyPzE(temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy());
//     solevent.em.SetPxPyPzE(temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy());

     jacobian[i]=TMath::Abs(2*p_Higgs.Energy()*nbZ/nbE-2*p_Higgs.Pz());
     
   }
   else{
     jacobian[i]=-2;
     /*
     double B2_4AC=k[1]*k[1]-4*k[2]*k[0];
     cout << " sumE "<< sumE <<" sumZ "<< sumZ <<" nbT " << nbT <<" MHiggs " << MHiggs <<" A "<< A <<endl;
     
     cout << " Error Sol:" << fabs(p_Higgs.M()-MHiggs)/MHiggs <<" "<<   fabs(eval/k[2])
          << " MHiggs " << MHiggs <<" sol "<< p_Higgs.M() <<" "<< nbZ
          << " B*B-4AC "<<k[1]*k[1]-4*k[2]*k[0] <<endl;
     if(B2_4AC>0){
       B2_4AC=TMath::Sqrt(B2_4AC);
       cout <<" criteria "<< A*A+(sumZ+sumE)*(sumZ-sumE)*nbT*nbT <<" should be >0 "<<endl;
       cout <<" sol+ "<< (-k[1]+B2_4AC)/2/k[2]
            <<" sol- "<< (-k[1]-B2_4AC)/2/k[2]
	    <<" eqn "<< k[2]<<" "<<k[1]<<" "<<k[0]
	    <<endl;
     }	  
     cout <<" ep "<< temp->ep.Px() <<" "<< temp->ep.Py()<<" "<< temp->ep.Pz()<<" "<<temp->ep.Energy()<<endl;
     cout <<" em "<< temp->em.Px() <<" "<< temp->em.Py()<<" "<< temp->em.Pz()<<" "<<temp->em.Energy()<<endl;
     cout <<" nu "<< solevent.nu.Px()<<" "<<solevent.nu.Py()<<" "<<solevent.nu.Pz()<<" "<<solevent.nu.Energy()<<endl;
     cout <<" nb "<< solevent.nb.Px()<<" "<<solevent.nb.Py()<<" "<<solevent.nb.Pz()<<" "<<solevent.nb.Energy()<<endl;
     cout <<" hmass "<< (temp->ep+temp->em+solevent.nu+solevent.nb).M()<<endl;
     */
/*
     if(DEBUGFLAG && fabs(p1.M()-Mw1)/Mw1<1E-4 &&  fabs(p2.M()-Mw2)/Mw2<1E-4){
        cout <<"= Bad Solution Begin ====================="<<endl;
        printf("ep: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->ep.M(),temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy() );
	printf("em: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->em.M(),temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy() );
	cout << "Mwp=" << Mw1<<" Mwm="<< Mw2 <<" nuPz="<< nuZ<<" nbPz="<< nbZ<<endl;
      	cout << "sol nuX=" << xx <<" eval/k[4]=" << eval/k[4] << " Mw1=" << p1.M() << " Mw2=" << p2.M() <<endl;
	cout <<"= Bad Solution End ==================="<<endl;
     }
*/
   }
   
 }//loop for Two solutions
		     
		     
}

void WWL1L2Sol_Mw1MHiggs( event_type* temp
                        ,double Mw1, double MHiggs, double nuX, double nuY
                        ,event_type* sol,double* jacobian){
			
 double limits[2]={-2*3500,2*3500};
 double roots_nuZ [2];
 double roots_nbZ [2][2];
			
 double epX=temp->ep.Px();double  epY=temp->ep.Py();double  epZ=temp->ep.Pz();double  epE=temp->ep.Energy();
 double emX=temp->em.Px();double  emY=temp->em.Py();double  emZ=temp->em.Pz();double  emE=temp->em.Energy();    
 double nbX=-epX-emX-nuX; double  nbY=-epY-emY-nuY; 
 
 double B=Mw1*Mw1/2+epX*nuX+epY*nuY;
 
 double k_nuZ[10];
 k_nuZ[2]=(epZ-epE)*(epZ+epE);
 k_nuZ[1]= 2*epZ*B;
 k_nuZ[0]= B*B-epE*epE*(nuX*nuX+nuY*nuY);
  
 PolySolver::FindRoots(2,k_nuZ,limits,roots_nuZ);
/*
cout <<"HDebug in  epZ " << epZ <<" epE " << epE  <<endl;
std::cout <<"HDebug in  nuX " << nuX << " nuY " << nuY <<" B "<< B <<" Mw "<< Mw1*Mw1<<endl;
cout <<"HDebug eqn " << k_nuZ[2] <<" "<<k_nuZ[1]<<" "<<k_nuZ[0]<<endl;
cout <<"HDebug sol " << roots_nuZ[0]<<" "<< roots_nuZ[1] <<endl; 
*/
 for(int i=0;i<4;i++) jacobian[i]=1;
  
 for(int i=0;i<2;i++){
    double nuZ = roots_nuZ[i];

    double eval_nuZ = k_nuZ[0] + k_nuZ[1]*nuZ+ k_nuZ[2]*nuZ*nuZ;
    if(fabs(eval_nuZ/k_nuZ[2])>1){
      jacobian[2*i]   =-1;
      jacobian[2*i+1] =-1;
      roots_nbZ[i][0]=0;
      roots_nbZ[i][1]=0;      
      continue;
    }
    //Debug
    /*
     TLorentzVector tempNu;
     tempNu.SetXYZM(nuX,nuY,nuZ,0);
     TLorentzVector WpP4 = tempNu + temp->ep;
     cout <<"Debug: nuZ=" 
          << nuZ 
	  <<" eval/k[2]=" 
	  << eval_nuZ/k_nuZ[2] 
	  << " WpM=" << WpP4.M() <<" =?= " << Mw1 <<endl;
    */
    
    double nuE = TMath::Sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
    double E3  = epE+emE+nuE;
    double Pz  = epZ+emZ+nuZ;
    double C = MHiggs*MHiggs-E3*E3+Pz*Pz-nbX*nbX-nbY*nbY;

    double k_nbZ[3];
    k_nbZ[2]=4*(E3+Pz)*(E3-Pz);
    k_nbZ[1]=-4*Pz*C;
    k_nbZ[0]=4*E3*E3*(nbX*nbX+nbY*nbY)-C*C;

    double roots[10];
    PolySolver::FindRoots(2,k_nbZ,limits,roots);

    for(int j=0;j<2;j++){
        double nbZ=roots[j];    
        double eval_nbZ = k_nbZ[0]+k_nbZ[1]*nbZ+k_nbZ[2]*nbZ*nbZ;
        if(fabs(eval_nbZ/k_nbZ[2]) <1){
            roots_nbZ[i][j]=nbZ;
	}else{
            jacobian[2*i+j]=-2;
	    roots_nbZ[i][j]=0;
	}
	    //Debug
     /*	    
     TLorentzVector tempNb;
     tempNb.SetXYZM(nbX,nbY,nbZ,0);
     TLorentzVector HiggsP4 = tempNu +tempNb + temp->ep + temp->em;
     cout <<"Debug: nbZ=" 
          << nbZ 
	  <<" => "<<k_nbZ[0] <<"+"<<k_nbZ[1]<<"*x+"<<k_nbZ[2] <<" *X^2=0 => " 
	  <<" eval/k[2]=" 
	  << eval_nbZ/k_nbZ[2] 
	  << " HiggsM=" << HiggsP4.M() <<" =?= " << MHiggs <<endl;
     */

    }//check nbZ solutions
 }
 
 
 

 for(int i=0;i<2;i++){
  for(int j=0;j<2;j++){ 
   int solIdx=2*i+j;
   event_type& solevent=sol[solIdx];
   double nuZ=roots_nuZ[i];
   double nbZ=roots_nbZ[i][j];
   
   solevent.nu.SetXYZM(nuX, nuY, nuZ,0);
   solevent.nb.SetXYZM(nbX, nbY, nbZ,0);

   TLorentzVector p_Wp   =temp->ep+solevent.nu;
   
   TLorentzVector p_Higgs=temp->ep+solevent.nu
                         +temp->em+solevent.nb;
  

   if(jacobian[solIdx]<0) continue;
   
   if(fabs(p_Higgs.M()-MHiggs)/MHiggs<1E-4 &&   fabs(p_Wp.M()-Mw1)/Mw1<1E-4){


     solevent.ep.SetPxPyPzE(temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy());
     solevent.em.SetPxPyPzE(temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy());
     double J3nuZ=2*p_Wp.Energy()   *solevent.nu.Pz()/solevent.nu.Energy()-2*p_Wp.Pz();
     double J4nbZ=2*p_Higgs.Energy()*solevent.nb.Pz()/solevent.nb.Energy()-2*p_Higgs.Pz();
     jacobian[solIdx]=TMath::Abs(J3nuZ*J4nbZ);
         
   }
   else{
     jacobian[solIdx]=-3;
/*
     if(DEBUGFLAG && fabs(p1.M()-Mw1)/Mw1<1E-4 &&  fabs(p2.M()-Mw2)/Mw2<1E-4){
        cout <<"= Bad Solution Begin ====================="<<endl;
        printf("ep: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->ep.M(),temp->ep.Px(),temp->ep.Py(),temp->ep.Pz(),temp->ep.Energy() );
	printf("em: %4.8f %4.8f %4.8f %4.8f %4.8f\n",temp->em.M(),temp->em.Px(),temp->em.Py(),temp->em.Pz(),temp->em.Energy() );
	cout << "Mwp=" << Mw1<<" Mwm="<< Mw2 <<" nuPz="<< nuZ<<" nbPz="<< nbZ<<endl;
      	cout << "sol nuX=" << xx <<" eval/k[4]=" << eval/k[4] << " Mw1=" << p1.M() << " Mw2=" << p2.M() <<endl;
	cout <<"= Bad Solution End ==================="<<endl;
     }
*/
   }
  }//loop for nbZ solution
 }//loop for nuZ solution
			
			
}


