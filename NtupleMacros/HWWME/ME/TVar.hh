#ifndef EvtProb_VAR
#define EvtProb_VAR

#include <TLorentzVector.h>
#include "TH2F.h"
#include "TH1F.h"

#define EBEAM 3500.0
#define fbGeV2 0.389379E12
#define sixteen_2Pi_to_8 3.88650230418250561e+07
#define   eight_2Pi_to_5 7.83410393050320417e+04
#define    four_2Pi_to_2 39.478417604357432
class TVar{
  public:

       enum MatrixElement{
          MCFM,
          MadGraph
       };


       enum Process{
          WW      = 0,
          Wp_1jet = 1,
          Wm_1jet = 2,
          Wp_gamma= 3,
          Wm_gamma= 4,
          ZZ      = 5,
          HWW120  = 6,
          HWW130  = 7,
          HWW140  = 8,
          HWW150  = 9,
          HWW160  =10,
          HWW170  =11,
          HWW180  =12,
          HWW190  =13,
          HWW200  =14,
          HWW210  =15,
          HWW220  =16,
          HWW230  =17,
          HWW250  =18,
          HWW300  =19,
          HWW     =20,
	  Z_2l    =21,
          WpZ_lostW=22,
          WpZ_lostZ=23,
          ZZ_4l    =24,
          Wgamma   =25,
	  Wegamma  =26,
	  Wmugamma =27,
	  Wtaugamma=28,
          ttbar    =29,
          Zee      =30,
          Zmumu    =31,
          Ztautau  =32,
          Data     =33,
	  Data_runFakes=34,
          Wenu     =35,
          Wmunu    =36,
          Wtaunu   =37,
          We1jet   =38,
          WZ       =39,
          HZZ      =40,
          Wmu1jet  =41,
	  WWj      =42,
	  HWWj     =43,
	  Null
       };


       enum HWWPhaseSpace{
          MHMW=0,
          MH  =1,
          MWMW=2,
          MHYH=3
       };


  //---------------------------------
  // Function
  //---------------------------------
  static TString ProcessName(int temp){

       if(temp==TVar::WW      ) return TString("WW");
  else if(temp==TVar::Wp_gamma) return TString("Wp_gamma");
  else if(temp==TVar::Wm_gamma) return TString("Wm_gamma");
  else if(temp==TVar::Wp_1jet ) return TString("Wp_1jet");
  else if(temp==TVar::Wm_1jet ) return TString("Wm_1jet");
  else if(temp==TVar::ZZ      ) return TString("ZZ");
  else if(temp==TVar::WZ      ) return TString("WZ");
  else if(temp==TVar::HWW     ) return TString("HWW");
  else if(temp==TVar::HZZ     ) return TString("HZZ");
  else if(temp==TVar::ZZ_4l   ) return TString("ZZ_4l");
  else if(temp==TVar::Z_2l    ) return TString("Z_2l");
  else if(temp==TVar::WpZ_lostW) return TString("WpZ_lostW");
  else if(temp==TVar::WpZ_lostZ) return TString("WpZ_lostZ");

  else if(temp==TVar::Wegamma ) return TString("Wegamma");
  else if(temp==TVar::Wmugamma) return TString("Wmugamma");
  else if(temp==TVar::Wtaugamma ) return TString("Wtaugamma");

  else if(temp==TVar::We1jet ) return TString("We1jet");
  else if(temp==TVar::Wmu1jet) return TString("Wmu1jet");

  else if(temp==TVar::Wenu   ) return TString("Wenu");
  else if(temp==TVar::Wmunu  ) return TString("Wmunu");
  else if(temp==TVar::Wtaunu ) return TString("Wtaunu");

  else if(temp==TVar::Data   ) return TString("Data");
  else if(temp==TVar::Data_runFakes) return TString("Data_runFakes");

  else if(temp==TVar::Zee     ) return TString("Zee");
  else if(temp==TVar::Zmumu   ) return TString("Zmumu");
  else if(temp==TVar::Ztautau ) return TString("Ztautau");
  else if(temp==TVar::ttbar   ) return TString("ttbar");

  else if(temp==TVar::HWW120  ) return TString("HWW120");
  else if(temp==TVar::HWW130  ) return TString("HWW130");
  else if(temp==TVar::HWW140  ) return TString("HWW140");
  else if(temp==TVar::HWW150  ) return TString("HWW150");
  else if(temp==TVar::HWW160  ) return TString("HWW160");
  else if(temp==TVar::HWW170  ) return TString("HWW170");
  else if(temp==TVar::HWW180  ) return TString("HWW180");
  else if(temp==TVar::HWW190  ) return TString("HWW190");
  else if(temp==TVar::HWW200  ) return TString("HWW200");

  else if(temp==TVar::HWW210  ) return TString("HWW210");
  else if(temp==TVar::HWW220  ) return TString("HWW220");
  else if(temp==TVar::HWW230  ) return TString("HWW230");
  else if(temp==TVar::HWW250  ) return TString("HWW250");
  else if(temp==TVar::HWW300  ) return TString("HWW300");

  else if(temp==TVar::WWj      ) return TString("WWj");
  else if(temp==TVar::HWWj     ) return TString("HWWj");

  return TString("UnKnown");
    
  };
 
  ClassDef(TVar,0)
};

struct branch_particle {
  int   PdgCode   ;
  int   Charge    ;
  double Px       ;
  double Py       ;
  double Pz       ;
  double E        ;
  double Eta      ;
  double Phi      ;

};
static const TString branch_format_particle =
 "PdgCode/I:"
 "Charge/I:"
 "Px/D:"
 "Py/D:"
 "Pz/D:"
 "E/D:"
 "Eta/D:"
 "Phi/D";

  
struct cdf_event_type{
  int PdgCode[2];
  TLorentzVector p[2];
  double MetX;
  double MetY;

  double Xsec   [60];
  double XsecErr[60];  
};
struct mcfm_event_type{
  int PdgCode[6];
  TLorentzVector p[6];
  double pswt;
};
struct event_type{
  TLorentzVector p1,p2,ep,em,nu,nb;
  double PSWeight;
};

struct array_event_type{
  int PdgCode[6];
  double p4[4][6];
  double PSWeight;
};

struct event_XSec{
  double dXSec, dXSecErr;
};

struct phasespace_Mw1Mw2{
	           double nuZ,nbZ,Mw1,Mw2,PSWeight;
};


struct phasespace_4D{
	   double x1,x2,costh_nu,phi_nu,PSWeight;
};


struct rand_type{
	   double r[10];
};

struct anomcoup{
	   double delg1_z, delg1_g, lambda_g, lambda_z, delk_g, delk_z_,tevscale;
};

struct EffHist{
  TH2F* els_eff_mc;
  TH2F* mus_eff_mc;
};

struct BoostHist{
  TH1F* kx;
  TH1F* ky;
};

#endif
