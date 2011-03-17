#include "iostream"
#include "PhaseSpace.hh"


//============================================
// Higgs PhaseSpace
//============================================

void genMHiggsYHiggs(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

    double wt[4];
    TLorentzVector P4ll=cdf_event.p[0]+cdf_event.p[1];
   
    double Mll=P4ll.M();
    double mminsq=Mll*Mll;
    double mmaxsq=EBEAM*EBEAM;
    double msqHiggs=1.,YHiggs;

    double width=masses_mcfm_.hwidth;
    width=2;
    breitw(r[2],mminsq,mmaxsq,masses_mcfm_.hmass,width,&msqHiggs,&wt[2]);

    double Ell=P4ll.Energy();
          
    double MHiggs=TMath::Sqrt(msqHiggs);
    
    
    double YMax=TMath::Log(4*EBEAM/MHiggs);
    double YMin=TMath::Log(2*Ell/MHiggs);
    if(YMin<0.0 ) YMin=0;
    if(YMin>=YMax) {
       cout <<"MHYH PhaseSpace Error: YMin " << YMin 
            <<" YMax "<< YMax <<" MHiggs " <<MHiggs
            <<" Ell " << Ell
            <<endl;
    } 

    YHiggs=(YMax-YMin)*(2*r[3]-1);


    if(YHiggs>0) YHiggs+=YMin;
    if(YHiggs<0) YHiggs-=YMin;
    wt[3]=2*(YMax-YMin);

    double nu_X,nu_Y;
    
//========================
// Boost to CM frame 
//========================

 double exppY=TMath::Exp(YHiggs);
 double expmY=1./exppY;
 double qEold =MHiggs*(exppY+expmY)/2;
 double qX=0.,qY=0.;
 
  double nuPtMax=qEold-cdf_event.p[0].Energy()-cdf_event.p[1].Energy();

 if(nuPtMax<0){
  cout << "Error: M= "<< MHiggs <<" Y= "<<YHiggs <<" E="<< qEold 
       << " E0="<< cdf_event.p[0].Energy()<<" E1="<<cdf_event.p[0].Energy()
       << " nuPtMax " << nuPtMax <<endl;
  PSList[0].pswt=-1;
  PSList[1].pswt=-1;
  return;
 }

 double nuPt=r[0]*nuPtMax; wt[0]=nuPtMax*nuPt;
 double phi =r[1]*2*TMath::Pi(); wt[1]=2*TMath::Pi();
 nu_X = nuPt*TMath::Cos(phi);
 nu_Y = nuPt*TMath::Sin(phi);

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
    

    WWL1L2Sol_MHiggsYHiggs( &cdf_event,
                      qX, qY, MHiggs, YHiggs, nu_X, nu_Y, 
                      PSList);



    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
      
       if(temp.pswt>0){
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8
			     /2;
       }
    }

}



//============================================
// Higgs PhaseSpace
//============================================

void genMHiggs(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist) {

    double wt[6];
     
    double nu_X,nu_Y,nu_Z;
    EtExponential(20.,r[0],&nu_X,&wt[0]);  
    EtExponential(20.,r[1],&nu_Y,&wt[1]);
    EtExponential(80.,r[2],&nu_Z,&wt[2]);
    
    TLorentzVector P4ll=cdf_event.p[0]+cdf_event.p[1];
    double Mll=P4ll.M();

    double mminsq=Mll*Mll;
    double mmaxsq=EBEAM*EBEAM;
    double msqHiggs=1.;
    breitw(r[3],mminsq,mmaxsq,masses_mcfm_.hmass,masses_mcfm_.hwidth,&msqHiggs,&wt[3]);
    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
    
    double qX (0.0), qY (0.0);

    if(SmearLevel >= 1) {
      getProbFromHist(r[4], & qX, & wt[4], boosthist.kx);
      getProbFromHist(r[5], & qY, & wt[5], boosthist.ky);
      PSWeight = PSWeight*wt[4]*wt[5];
      // cout << "; " << wt[4] << "; " << wt[5] ;  
    }
    
    // cout << "\n" ;

    // cout << "PhaseSpace::genMHiggs = " << qX <<"; qY = " <<qY <<endl;
    WWL1L2Sol_MHiggs( &cdf_event,
                      qX, qY, msqHiggs, nu_X, nu_Y, nu_Z,
                      PSList);

    
    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
      
       if(temp.pswt>0){
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;

       }
    }
}


//============================================
// HiggsW* PhaseSpace
//============================================

void genMHiggsMw1(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList,  BoostHist boosthist) {

    double wt[6];
     
    double nu_X,nu_Y;
    EtExponential(20.,r[0],&nu_X,&wt[0]);  
    EtExponential(20.,r[1],&nu_Y,&wt[1]);

    TLorentzVector P4ll=cdf_event.p[0]+cdf_event.p[1];
    double Mll=P4ll.M();
    double mminsq=Mll*Mll;
    double mmaxsq=EBEAM*EBEAM;
    double msqHiggs=1.;
    breitw(r[2],mminsq,mmaxsq,masses_mcfm_.hmass,masses_mcfm_.hwidth,&msqHiggs,&wt[2]);

    double msqW=1.;
    double rmass =80.419;
    double rwidth=2.06;
    mminsq=1e-15;
    mmaxsq=TMath::Sqrt(msqHiggs)-Mll;
    mmaxsq=mmaxsq*mmaxsq;
    breitw(r[3],mminsq,mmaxsq,rmass,rwidth,&msqW,&wt[3]);    

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];

    //System Pt
    double qX (0.0), qY (0.0);

      if(SmearLevel >= 1) {
	getProbFromHist(r[4], & qX, & wt[4], boosthist.kx);
	getProbFromHist(r[5], & qY, & wt[5], boosthist.ky);
	PSWeight = PSWeight*wt[4]*wt[5];
      }


    WWL1L2Sol_MHiggsMw1( &cdf_event,
                      qX, qY, msqHiggs, msqW, nu_X, nu_Y,
                      PSList);

    for(int i=0;i<4;i++){
        
       mcfm_event_type& temp = PSList[i];

       if(temp.pswt>0){
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;

       }
    }
}



//============================================
// Phase Space for DY Integration
// Dgf: 2
// Et1, Et2
//============================================
void genDY(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

    double PSWeight=1;
    //System Pt
      double qX=cdf_event.p[0].Px()+cdf_event.p[1].Px();
      double qY=cdf_event.p[0].Py()+cdf_event.p[1].Py();

   
   mcfm_event_type& solevent = PSList[0];
   
   double qE = cdf_event.p[0].Energy()+ cdf_event.p[1].Energy();   
   double qZ = cdf_event.p[0].Pz()    + cdf_event.p[1].Pz();   
   double qEold= TMath::Sqrt(qE*qE-qX*qX-qY*qY);
   double x0 = (qEold+qZ)/2/EBEAM;
   double x1 = (qEold-qZ)/2/EBEAM;

   double scale0=x0/(x0+x1);
   double scale1=x1/(x0+x1);

   solevent.p[0].SetPxPyPzE   (scale0*qX, scale0*qY, x0*EBEAM, scale0*qE);
   solevent.p[1].SetPxPyPzE   (scale1*qX, scale1*qY,-x1*EBEAM, scale1*qE);
   solevent.p[2] = cdf_event.p[0];
   solevent.p[3] = cdf_event.p[1];
       
   solevent.pswt = PSWeight/solevent.p[2].Energy()
			   /solevent.p[3].Energy()
                           /four_2Pi_to_2;


   if(x0 > 1 || x1>1 || x0<=xmin_.xmin || x1 <=xmin_.xmin) {
     solevent.pswt=-1;
   }

}
//============================================
// Phase Space for WW Integration
//
// Two on-shell WW PhaseSpace
//
//============================================

void genMw1Mw2(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist) {
    double wt[6];
    double nu_Z,nb_Z;
    EtExponential(80.,r[0],&nu_Z,&wt[0]);  
    EtExponential(80.,r[1],&nb_Z,&wt[1]);

    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq2=1.;
    double msq3=1.;
    breitw(r[2],mminsq,mmaxsq,breit_.mass2,breit_.width2,&msq2,&wt[2]);
    breitw(r[3],mminsq,mmaxsq,breit_.mass3,breit_.width3,&msq3,&wt[3]);

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
   
    //System Pt
    //  test how pT shift affects the Event prob
    //  double qX=10.,qY=10.;
    double qX=0, qY=0;
    if(SmearLevel >= 1) {
      getProbFromHist(r[4], & qX, & wt[4], boosthist.kx);
      getProbFromHist(r[5], & qY, & wt[5], boosthist.ky);
      PSWeight = PSWeight*wt[4]*wt[5];
      // cout << "; " << wt[4] << "; " << wt[5] ;  
    }
    //    cout << "\n" ;
    

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
    WWL1L2Sol_Mw1Mw2( &cdf_event,
                      qX, qY, msq2, msq3, nu_Z, nb_Z,
                      PSList);
		      
    for(int i=0;i<4;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
//         temp.pswt = PSWeight/temp.pswt/sixteen_2Pi_to_8;
       
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;
/*
cout <<"DebugE "<<   temp.p[2].Energy()
	<<" "<<  temp.p[3].Energy()
	<<" "<<	  temp.p[4].Energy()
	<<" "<<	  temp.p[5].Energy()
	<<endl;
*/
       }else if(temp.pswt!=temp.pswt){
         cout <<"Error: WWL1L2Sol_Mw1Mw2 pswt= " << temp.pswt <<endl;
       }
       
    }
}

//============================================
// Phase Space for WZ Integration
//
// Two on-shell WZ PhaseSpace
//
//============================================

/*void genMw1Mz2(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

//    double sixteen_twopi_to_eighth=16*TMath::Power(2*TMath::Pi(),8);
    double wt[4];
     
    double nu_Z,nb_Z;
    EtExponential(80.,r[0],&nu_Z,&wt[0]);  
    EtExponential(80.,r[1],&nb_Z,&wt[1]);

    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq2=1.;
    double msq3=1.;
   // wmass
    breitw(r[2],mminsq,mmaxsq,breit_.mass3,breit_.width3,&msq2,&wt[2]);
   // zmass
    breitw(r[3],mminsq,mmaxsq,breit_.mass2,breit_.width2,&msq3,&wt[3]);

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
   

    //System Pt
      double qX=0.,qY=0.;

    
//    mcfm_event_type  sol[4];

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
// 61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))'
// W+Z -> L+ from W and L- from Z
// 71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'  =>lost e^+(p6)

//  if(nwz_.nwz==-1){
// W-Z -> L- from W and L+ from Z
// Swap definition of 
// => Swap 
// 76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))' =>lost e^-(p5)
// }
      
    WWL1L2Sol_Mw1Mw2( &cdf_event,
                      qX, qY, msq2, msq3, nu_Z, nb_Z,
                      PSList);
		      
    for(int i=0;i<4;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
//         temp.pswt = PSWeight/temp.pswt/sixteen_2Pi_to_8;
       
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;

       }else if(temp.pswt!=temp.pswt){
         cout <<"Error: WWL1L2Sol_Mw1Mw2 pswt= " << temp.pswt <<endl;
       }
       
    }
}

*/

//#################################
//
// Phase Space for Wgamma integration
//
//#################################

void genMw_Wgamma(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

//    double sixteen_twopi_to_eighth=16*TMath::Power(2*TMath::Pi(),8);
    double PSWeight;
     
    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq=1.;
    
    //breitw(r[0],mminsq,mmaxsq,breit_.mass2,breit_.width2,&msq,&PSWeight);
    GenMwModified(r[0], mminsq, mmaxsq, &msq, &PSWeight);
   

    //System Pt
      double qX=0.,qY=0.;

    
//    mcfm_event_type  sol[4];

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
    WWL1L2Sol_Mw( &cdf_event,
                   qX, qY, msq,
                   PSList);
		      
    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
         //temp.pswt = PSWeight/temp.pswt/eight_2Pi_to_5; 

         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
                             /eight_2Pi_to_5;

       }
//cout <<Form("DebugPhase %d %8.4E  %8.4E  %8.4E qT %8.4E  %8.4E msq %8.4E E %8.4E  %8.4E  %8.4E pswt  %8.4E ",i,r[0],r[2],r[4],qX,qY,msq,temp.p[2].Energy(),temp.p[3].Energy(),temp.p[4].Energy(),temp.pswt )<<endl;
    }
}
//####################################
//
// Phase Space for W+1jet
//
// Dgf: 5
// MW, qX, qY, Et1, Et2
//
//####################################
void genMw_W1jet(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

//    double sixteen_twopi_to_eighth=16*TMath::Power(2*TMath::Pi(),8);
    double PSWeight;
     
    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq=1.;
    
    breitw(r[0],mminsq,mmaxsq,breit_.mass2,breit_.width2,&msq,&PSWeight);
    //GenMwModified(r[0], mminsq, mmaxsq, &msq, &PSWeight);
   

    //System Pt
      double qX=0.,qY=0.;

    
//    mcfm_event_type  sol[4];

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
    WWL1L2Sol_Mw( &cdf_event,
                   qX, qY, msq,
                   PSList);
		      
    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
         //temp.pswt = PSWeight/temp.pswt/eight_2Pi_to_5;
 
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
                             /eight_2Pi_to_5;

       }
    }
}


//=============================================================
//
//  Phase Space for Z(l+l-)Z(nu nu~)
//
//==============================================================


void genMzNu3(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

    double wt[4];
     
    double nu_X,nu_Y,nu_Z;
    EtExponential(40.,r[0],&nu_X,&wt[0]);  
    EtExponential(40.,r[1],&nu_Y,&wt[1]);
    EtExponential(80.,r[2],&nu_Z,&wt[2]);

    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq=1.;
    breitw(r[3],mminsq,mmaxsq,breit_.mass3,breit_.width3,&msq,&wt[3]);

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
   

    //System Pt
      double qX=0.,qY=0.;
    
//    mcfm_event_type  sol[4];

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
    WWL1L2Sol_MzNu3( &cdf_event,
                      qX, qY, msq, nu_X, nu_Y, nu_Z,
                      PSList);
		      
    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
         //temp.pswt = PSWeight/temp.pswt/sixteen_2Pi_to_8;
 
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;

       }
    }
}

//=============================================================
//
//  Phase Space for WZ->(l+/- nu) (l+l-)
//  71 '  f(p1)+f(p2) --> W^+(-->nu(p3)+mu^+(p4))+Z^0(-->e^-(p5)+e^+(p6))'
//  76 '  f(p1)+f(p2) --> W^-(-->mu^-(p3)+nu~(p4))+Z^0(-->e^-(p5)+e^+(p6))'
//==============================================================


void genMwNu3(double* r,int SmearLevel,cdf_event_type cdf_event, mcfm_event_type* PSList){

    double wt[4];
     
    double nu_X,nu_Y,nu_Z;
    EtExponential(40.,r[0],&nu_X,&wt[0]);  
    EtExponential(40.,r[1],&nu_Y,&wt[1]);
    EtExponential(80.,r[2],&nu_Z,&wt[2]);

    double mminsq=1e-15;
    double mmaxsq=EBEAM*EBEAM;
    double msq=1.;
   // wmass
    breitw(r[3],mminsq,mmaxsq,breit_.mass3,breit_.width3,&msq,&wt[3]);

    double PSWeight=wt[0]*wt[1]*wt[2]*wt[3];
   

    //System Pt
      double qX=0.,qY=0.;

//    mcfm_event_type  sol[4];

//nu_Z=-0.62694; nb_Z=16.8335;
//Mw1=81.0633; Mw2=78.2776;
    WWL1L2Sol_MzNu3( &cdf_event,
                      qX, qY, msq, nu_X, nu_Y, nu_Z,
                      PSList);
//------------------------------------------
// The solution is ordered as
// 82 '  f(p1)+f(p2) --> Z^0(-->3*(nu(p3)+nu~(p4))) + Z^0(-->e^-(p5)+e^+(p6))'
//
		      
    for(int i=0;i<2;i++){
    
       mcfm_event_type& temp = PSList[i];
       if(temp.pswt>0){
         //temp.pswt = PSWeight/temp.pswt/sixteen_2Pi_to_8;
 
         temp.pswt = PSWeight/temp.pswt
                             /temp.p[2].Energy()
			     /temp.p[3].Energy()
			     /temp.p[4].Energy()
			     /temp.p[5].Energy()
                             /sixteen_2Pi_to_8;

       }
    }
}


//=============================================================
//
//  Generate Phase Space for 3 dimension
//
//==============================================================


void gen3(int* jbranch,double* r,double p[][mxpart],double* wt3){
//c----generate 3 dimensional phase space weight and vectors p(7,4)
//c----and x1 and x2 given seven random numbers
//c----p(5,i) and p(4,i) are set equal to zero


//cout <<"Debug gen3: Rnd: "; for(int i=0;i<7;i++) cout<<" " <<r[i]; cout <<endl;
      *wt3=0.0;
      double tau=TMath::Exp(TMath::Log(taumin_.taumin)*r[6-1]);
      double y=0.5*TMath::Log(tau)*(1.-2.*r[7-1]);
      double xjac=TMath::Log(taumin_.taumin)*tau*TMath::Log(tau);

//cout <<"Debug gen3 tau=" << tau <<" y=" << y <<" xjac=" << xjac <<" taumin=" << taumin_.taumin <<endl;
      double xx[2];
      xx[1-1]=TMath::Sqrt(tau)*TMath::Exp(+y);
      xx[2-1]=TMath::Sqrt(tau)*TMath::Exp(-y);

//cout <<"Debug gen3: xx " << xx[0] <<" "<< xx[1] <<endl;
//c---if x's out of normal range alternative return
      if   ((xx[1-1] > 1.0) ||
     	    (xx[2-1] > 1.0) ||
     	    (xx[1-1] < xmin_.xmin)||
     	    (xx[2-1] < xmin_.xmin)) return ;

      double p1[4],p2[4],p3[4],p4[4],p5[4],p6[4],p7[4];
      
      p1[4-1]=-xx[1-1]*energy_.sqrts*0.5;
      p1[1-1]=0.0;
      p1[2-1]=0.0;
      p1[3-1]=-xx[1-1]*energy_.sqrts*0.5;

      p2[4-1]=-xx[2-1]*energy_.sqrts*0.5;
      p2[1-1]=0.0;
      p2[2-1]=0.0;
      p2[3-1]=+xx[2-1]*energy_.sqrts*0.5;


      double pswt=0;
      phase3(jbranch,r,p1,p2,p3,p4,p5,p6,p7,&pswt);

      for(int nu=0;nu<4;nu++){
        p[nu][1-1]=p1[nu];
        p[nu][2-1]=p2[nu];
        p[nu][3-1]=p3[nu];
        p[nu][4-1]=p4[nu];
        p[nu][5-1]=p5[nu];
        p[nu][6-1]=p6[nu];
        p[nu][7-1]=p7[nu];
      }	 

       (*wt3)=xjac*pswt;



}
void phase3(int* jbranch,double* r,double* p1,double* p2,double* p3,double* p4,double* p5,double* p6,double* p7,double* wt){
//c----generate phase space for 2-->3 process
//c----r(mxdim),p1(4),p2(4) are inputs
//c----incoming p1 and p2 reversed in sign from physical values
//c----i.e. phase space for -p1-p2 --> p3+p4+p5
//c----with all 2 pi's (ie 1/(2*pi)^5)
//c----(p4,p5) are dummies
     double wt0=1./2/TMath::Pi();
     double m5=0.0;

     double p12[4];
     for(int j=0;j<4;j++){
       p12[j]=-p1[j]-p2[j];
       p6 [j]=0.0;
       p7 [j]=0.0;
     }
      double smin=0.0;

      double wt125;
      double p34[4];
//c---generate p5 and p34,
//c---smin is the minimum inv mass of 34 system
//c---m5 is the mass of p5
      phi1_2m(jbranch,m5,r[1-1],r[2-1],r[3-1],smin,p12,p5,p34,&wt125);

      double wt34;
//c---decay 34-system   
      phi3m0_(&r[4-1],&r[5-1],p34,p3,p4,&wt34);

      *wt=wt0*wt125*wt34;

}
void phi1_2m(int* jbranch, double m2,double x3,double xth,double xphi,double s3min,double* p1,double* p2,double* p3,double* wt){

//c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
//c     with invariant mass of particle three s3 integrated over.
//c     s3min is the minimum value of s3.
//c     Vectors returned p2 and p3 are in the same frame as p1 is supplied.
//c     Expression evaluated is
//c     ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
//c     delta(p2^2-m2) delta(p3^2-s3)
    double wt0=1./8./TMath::Pi();

      double xexp=1.0;
      *wt=0.0;
      double s1=p1[4-1]*p1[4-1]-p1[1-1]*p1[1-1]-p1[2-1]*p1[2-1]-p1[3-1]*p1[3-1];
      if (s1 < 0.0) return;
      double m1=TMath::Sqrt(s1);
      double s2=m2*m2;
      double s3max=(m2-m1)*(m2-m1);
      if (s3min > s3max) return;

      double s3,w3;
      //breitw_(&x3,&s3min,&s3max,&breit_.mass3,&breit_.width3,&s3,&w3);
      GenMwModified(x3, s3min, s3max, &s3, &w3);


     double   rtxth=TMath::Power(xth,xexp);
     double   xjac=1.0/(xexp*TMath::Power(xth,(xexp-1.0)));

      if (*jbranch == 1){
        *jbranch=2;
      }else if (*jbranch == 2) {
        *jbranch=1;
	rtxth=1.0-rtxth;
      }


     double costh=2.*rtxth-1.;
     double phi=2.*TMath::Pi()*xphi;
     double sinth=TMath::Sqrt(1.-costh*costh);
     double cphi=TMath::Cos(phi);
     double sphi=TMath::Sin(phi);
     double lambda=((s1-s2-s3)*(s1-s2-s3)-4.0*s2*s3);
      lambda=TMath::Sqrt(lambda);
//c      Eg=s1+s2-s3

      *wt =wt0*w3*lambda/s1/xjac;

      double p3cm[4];
      p3cm[4-1]=m1/2.*(s1+s3-s2)/s1;
      p3cm[1-1]=m1/2.*lambda/s1*sinth*sphi;
      p3cm[2-1]=m1/2.*lambda/s1*sinth*cphi;
      p3cm[3-1]=m1/2.*lambda/s1*costh;


      boost_mcfm_(&m1,p1,p3cm,p3);
      for(int j=0;j<4;j++)      p2[j]=p1[j]-p3[j];
      


      if (  p1[4-1] < 0.0 ||
            p2[4-1] < 0.0 ||
            p3[4-1] < 0.0 )       return;


}




