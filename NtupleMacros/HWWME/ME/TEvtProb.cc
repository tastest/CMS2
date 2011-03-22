//-----------------------------------------------------------------------------
//
// Class EventProb Module
//
//   EventProb Module
//
// March 21 2011
// S. Jindariani (sergo@fnal.gov)
// Y. Gao (ygao@fnal.gov)
// K. Burkett (burkett@fnal.gov)
//-----------------------------------------------------------------------------

#include "TEvtProb.hh"
#include "TVar.hh"
#include "TMatrixElement.hh"


ClassImp(TEvtProb)

using namespace std;

//-----------------------------------------------------------------------------
// Constructors and Destructor
//-----------------------------------------------------------------------------
TEvtProb::TEvtProb() {
  _matrixElement = TVar::MCFM;
  TMatrixElement::inst()->SetMatrixElement(_matrixElement);
  _hwwPhaseSpace=TVar::MH;
  mcfm_init_();
}

TEvtProb::~TEvtProb() {}


//-----------------------------------------------
// Global Variables
//-----------------------------------------------
TVar::Process Global_process;
TVar::HWWPhaseSpace Global_HWWPhaseSpace;
cdf_event_type Global_cdf_event;
mcfm_event_type Global_mcfm_event[4];
int Global_IntegrandType=0;
int Global_Ncalls;
int Global_SmearLevel;
int Global_IsApplyFake;
int Global_RandIdxMet;
int Global_RandIdxL1L2;
int Global_NSol;


//-----------------------------------------------
//
// Importance Sampleing Integration
//
//-----------------------------------------------

void   TEvtProb::NeutrinoIntegrate(TVar::Process proc,
                                         cdf_event_type cdf_event,
                                         double *Xsec,
                                         double *XsecErr){

    Global_process = proc;
    Global_HWWPhaseSpace = _hwwPhaseSpace;
    Global_cdf_event=cdf_event;
    Global_SmearLevel = _smearLevel;
    Global_Ncalls= _ncalls;
    Global_IsApplyFake=_isApplyFake;
    
    //Initialize Process
    SetProcess(Global_process);
    My_choose(Global_process);

    //delta(L1) delta(L2)
    int NDim = bveg1_mcfm_.ndim-6;
    //qx,qy
    if (Global_SmearLevel>=1) NDim+=2;
    //dE1,dE2
    if (Global_SmearLevel>=2) NDim+=2;

 
    cout <<" [NeutrinoIntegrate]: Evaluate " << TVar::ProcessName(proc)
     <<" Ncalls " << Global_Ncalls
     <<" npart._npart=" << npart_.npart
     <<" bveg1_mcfm_.ndim= " << bveg1_mcfm_.ndim
     <<" NDim= "<<NDim<<endl;


    //Phase space generation
    TRandom3 myRandom;
    myRandom.SetSeed(_seed);

    double r  [22] = { 0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
		       0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,
		       0.5,0.5};
    
    
    int count_PS=0;
    double sumW=0,sumW2=0;

    double probAcceptanceEfficiency = getProbAcceptanceEfficiency(cdf_event, _effhist);
    // cout << "probAcceptanceEfficency = " << probAcceptanceEfficiency << endl;
    // double probAcceptanceEfficiency = 1.0;
    if(probAcceptanceEfficiency == 0) return;
    if(probAcceptanceEfficiency<0) {
      cout <<"Error: " << probAcceptanceEfficiency <<endl;
      return;
    }

    for(int idx=0;idx<_ncalls;idx++){
      count_PS++;
      myRandom.RndmArray(NDim,r); // generate NDim random numbers and set the first NDim entries of r arrary
      double dXsec=0;
    
      // dXsec=Integrand_NeutrinoIntegration(r,NDim,0);
      dXsec=Integrand_NeutrinoIntegration(r,NDim,0, _boosthist)*probAcceptanceEfficiency;
  
      if (dXsec<=0) continue;
      // count_PS++;
      sumW  += dXsec;
      sumW2 += dXsec*dXsec;
    }
    
    *Xsec = sumW/count_PS;
    //    cout << "TEvtProb:: count_PS = " << count_PS << "; _ncalls = " << _ncalls << endl;
    
    *XsecErr = sumW2/count_PS-sumW/count_PS*sumW/count_PS;
    if(*XsecErr>0.0) *XsecErr = sqrt(*XsecErr/count_PS);
    else             *XsecErr = -1;  
 

}

//=======================================
// Integrand
//=======================================
double Integrand_NeutrinoIntegration(double * r, unsigned int NDim, void * param, BoostHist boosthist){

    //constants
    double sqrts = 2.*EBEAM;
    double W=sqrts*sqrts;


    Global_NSol = 4;


    //Weight calculation
    // double probAcceptanceEfficiency;
    double PSWeight=1.;
    double flux=1.;
    double dXsec=0.;
    double sumW =0.;

    int count_PS=0;
    
    count_PS++;
    
    PSWeight=1.;
  
    
    if (Global_process==TVar::WW)
      {
	Global_NSol=4;         
	genMw1Mw2(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event, boosthist);
      }

    if (Global_process==TVar::HWW)
      {
	if(Global_HWWPhaseSpace==TVar::MWMW)
	  { 
	    Global_NSol=4;
	    genMw1Mw2(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event, boosthist);
	  }
	else if (Global_HWWPhaseSpace==TVar::MHMW)
	  { 
	    Global_NSol=4;
	    genMHiggsMw1(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event, boosthist);
	  }
	else if (Global_HWWPhaseSpace==TVar::MHYH)
	  { 
	    Global_NSol=2;
	    genMHiggsYHiggs(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event, boosthist);
	  }
	else if (Global_HWWPhaseSpace==TVar::MH)
	  { 
	    Global_NSol=2;
	    genMHiggs(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event, boosthist);
	  }
      }

    else if (Global_process==TVar::Wp_gamma || Global_process==TVar::Wm_gamma )
      { 
	Global_NSol=2;  
	genMw_Wgamma(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event);
      }
    else if (Global_process==TVar::Wp_1jet || Global_process==TVar::Wm_1jet  )
      { 
	Global_NSol=2; 
	genMw_W1jet(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event);
      }
    else if (Global_process==TVar::ZZ)
      { 
	Global_NSol=2;
	genMzNu3     (r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event);
      }
    else if (Global_process==TVar::Z_2l)
      { 
	Global_NSol=1;    
	genDY(r,Global_SmearLevel,Global_cdf_event,Global_mcfm_event);
      } 
    
    //loop over solutions
    for(int jps=0;jps<Global_NSol;jps++){
      
      mcfm_event_type& mcfm_event=Global_mcfm_event[jps];
    
      // cout << "mcfm_event.pswt = " << mcfm_event.pswt << endl;
      if(mcfm_event.pswt<=0) continue;
      
      //Matrix Element evaluation in qX=qY=0 frame
      //Evaluate f(x1)f(x2)|M(q)|/x1/x2 
      // 
      double qX=mcfm_event.p[0].Px()+mcfm_event.p[1].Px();
      double qY=mcfm_event.p[0].Py()+mcfm_event.p[1].Py();
      // cout<< "TEvtProb:: qX = " << qX<< "qY = "<<qY<<"\n";
      
      if((qX*qX+qY*qY)>0){
	double qE = mcfm_event.p[0].Energy()+mcfm_event.p[1].Energy();
	TVector3 boostV(qX/qE,qY/qE,0);
	for(int ipt=0;ipt<6;ipt++) mcfm_event.p[ipt].Boost(-boostV);
      }
      
      //event selections in Lab Frame
      if(My_eventcuts(Global_process,&mcfm_event,&Global_cdf_event)){ mcfm_event.pswt=0; continue;}
      double flavor_msq[nmsq][nmsq];
      double msqjk = SumMatrixElementPDF(Global_process, &mcfm_event, flavor_msq, &flux);
      if(msqjk<=0){ mcfm_event.pswt=0; continue;}
      
      //npart_ means number of final state particles. 
      //Don't forget the first two are 2 initial state particles
      //for(int ipt=2;ipt<npart_.npart+2;ipt++) msqjk = msqjk/mcfm_event.p[ipt].Energy();
      
      flux=fbGeV2/(mcfm_event.p[0].Energy()*mcfm_event.p[1].Energy())	/(4*W);
      //dXsec=msqjk*flux*mcfm_event.pswt*PSWeight*probAcceptanceEfficiency;
      dXsec=msqjk*flux*mcfm_event.pswt*PSWeight;
      
      if (dXsec>0.0){
	sumW  += dXsec;
	mcfm_event.pswt=dXsec;
      }  
      else
	{
	  /*
	  cout <<" NeutrinoIntegrate Warning: dXsec " << dXsec
	       << " dXsec==dXsec " << (dXsec==dXsec)
	       << " dXsec>0.0 " << (dXsec>0.0)<<" "
	       <<TVar::ProcessName(Global_process)
	       <<" Msq="<<msqjk 
	       <<" flux="<<flux 
	       <<" wgt="<<mcfm_event.pswt
	       <<" PS=" <<PSWeight
	    //<<" eff="<<probAcceptanceEfficiency
	       <<endl;
	  */
	}
    }//loop solutions
    
    return sumW;
    
}


