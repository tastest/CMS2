#ifndef WW_GENERATOR
#define WW_GENERATOR
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


using namespace std;

void testGenerator();

void  breitw(double x1,double mminsq,double mmaxsq,double rmass,double rwidth,double *msq,double *wt);

double GenGaus(double* rand,Double_t mean, Double_t sigma);

void Gen_Wmass_Wgamma(int nbinsx, double* fIntegral,double r1, double* mass, double* wgt);
void Set_Wmass_Wgamma(double* hmass);

void GenMwModified(double r, double mminsq, double mmaxsq, double* msq, double* wgt);	
double GetRandom(TH1* his, double r1);
#endif
