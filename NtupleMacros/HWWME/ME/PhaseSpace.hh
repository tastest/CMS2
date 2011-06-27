#ifndef WW_PHASE
#define WW_PHASE
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <TFile.h>
#include <TF1.h>
#include "TRandom3.h"
#include <TVar.hh>
#include "TMCFM.hh"
#include "TUtil.hh"

#include "Generator.hh"
#include "KinematicSolver.hh"


using namespace std;


//-------------------------
void gen3(int* jbranch,double* r,double p[][12],double* wt3);
void phase3(int* jbranch, double* r,double*p1,double*p2,double*p3,double*p4,double*p5,double*p6,double*p7,double*wt);
void phi1_2m(int* jbranch, double m2,double x3,double xth,double xphi,double s3min,double* p1,double* p2,double* p3,double* wt);

//-------------------------
void genDY    (double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);

//-------------------------
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, double hmass);
void genMHiggsMw1(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMHiggs   (double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMHiggsYHiggs(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);

void genMw1Mw2(double* r,int SmearLevel,     event_type  cdf_event, event_type* sol, double* jacobian);
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type* cdf_event, array_event_type* sol);

//--------------------------
void genMw1Mw2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList, BoostHist boosthist);
void genMw1Mz2(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);

//---------------------------
void genMw_Wgamma(double* r,int SmearLevel, cdf_event_type cdf_event,mcfm_event_type* PSList);
void genMw_W1jet (double* r,int SmearLevel, cdf_event_type cdf_event,mcfm_event_type* PSList);

//---------------------------
void genMzNu3(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);
void genMwNu3(double* r,int SmearLevel, cdf_event_type cdf_event, mcfm_event_type* PSList);



#endif
