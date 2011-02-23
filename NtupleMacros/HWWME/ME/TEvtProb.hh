#ifndef _TEVTPROB_HH_
#define _TEVTPROB_HH_
//-----------------------------------------------------------------------------
// Description: Class TEvtProb: EvtProb base class
// ------------
//
//      Event Probability Density Calculation
//
// Feb 21 2011
// Sergo Jindariani
// Yanyan Gao
// Kevin Burkett
//-----------------------------------------------------------------------------
#include <sstream>
#include <string>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

#include "TObject.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1D.h"
#include "TFile.h"

#include "TMCFM.hh"
#include "TVar.hh"
#include "TUtil.hh"
#include "PhaseSpace.hh"
#include "foam/TFoamIntegrand.hh"

#include "TH2F.h"
#include "TH1F.h"
#include "assert.h"
#include "TROOT.h"



//----------------------------------------
// Class TEvtProb
//----------------------------------------
class TEvtProb : public TObject {

public:
  //--------------------
  // Variables
  //--------------------
  int _seed;
  int _ncalls;
  int _smearLevel;
  int _isApplyFake;

  TVar::Process _process;
  TVar::MatrixElement _matrixElement;
  TVar::HWWPhaseSpace _hwwPhaseSpace;
  EffHist _effhist;
  BoostHist _boosthist;
  
  //---------------------------------------------------------------------------
  // Constructors and Destructor
  //---------------------------------------------------------------------------
  TEvtProb();
 ~TEvtProb();

  //----------------------
  // Function
  //----------------------
  void SetSeed(int tmp) { _seed = tmp; }
  void SetNcalls(int tmp) { _ncalls= tmp; }
  void SetSmearLevel(int tmp) { _smearLevel=tmp; }
  void SetApplyFake (int tmp) { _isApplyFake=tmp; }
  void SetProcess(TVar::Process tmp) { _process = tmp; }
  void SetMatrixElement(TVar::MatrixElement tmp){ _matrixElement = tmp; }
  void SetHWWPhaseSpace(TVar::HWWPhaseSpace tmp){ _hwwPhaseSpace = tmp; }
  void SetMCHist(TVar::Process proc) {
    TString effFileName;
    if(TVar::ProcessName(proc) == "HWW") effFileName = "../ggH160_MCUtil.root";
    if(TVar::ProcessName(proc) == "WW")  effFileName = "../WW_MCUtil.root";
    if(TVar::ProcessName(proc) == "ZZ")  effFileName = "../ZZ_MCUtil.root";
    if(TVar::ProcessName(proc) == "WZ")  effFileName = "../ZZ_MCUtil.root";
    
    std::cout << "TEvtProb::SetMCHist: " << effFileName << std::endl;
    TFile *f = TFile::Open(effFileName, "READ");
    assert(f);
    gROOT->cd();
    _effhist.els_eff_mc = (TH2F*) f->Get("els_eff_mc")->Clone();
    _effhist.mus_eff_mc = (TH2F*) f->Get("mus_eff_mc")->Clone();
    _boosthist.kx = (TH1F*) f->Get("kx")->Clone();
    _boosthist.ky = (TH1F*) f->Get("ky")->Clone();  
    f->Close();
  }
  

  void LOXsec(double* Xsec,double* Err);
  void LOXsec_foam(double* Xsec,double* Err); 
 
  void NeutrinoIntegrate(TVar::Process proc,
                         cdf_event_type cdf_event,
                         double *Xsec,
                         double *XsecErr);

  void SetHiggsMass(double mass){
    masses_mcfm_.hmass=mass;
    masses_mcfm_.hwidth=HiggsWidth(mass);
}

 
  
  ClassDef(TEvtProb,0);
};


//----------------------------------------
// Global Function
//----------------------------------------
double Integrand_LOXsec(double * arg, unsigned int dim, void * param);

class  Integrand_LOXsec_foam : public TFoamIntegrand{
 public:
   Double_t Density(Int_t ndim, Double_t* );
};


double Integrand_NeutrinoIntegration(double * r, unsigned int NDim, void * param, BoostHist boosthist);

   
#endif
