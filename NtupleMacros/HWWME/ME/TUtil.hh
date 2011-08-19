#ifndef WW_COMMON
#define WW_COMMON
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
#include "TVar.hh"
#include "PolySolver.hpp"
#include "TMCFM.hh"



using namespace std;

TString DbnEventLepSelName(int i);

double GetContentsFromHistogram(TH1F*, double);
double Prob_phi_WW (double x);
double Prob_phi_HWW(double x);
double pdf40(double);
double pdf60(double);
double pdf80(double);
double pdfCosTh_nu(double);

void  EtGaus(double* rand, double mean, double sigma,double* results, double* wgt);
void  EtExponential(double alpha,double rand0,double* Et,double* wt);
void  breitw(double x1,double mminsq,double mmaxsq,double rmass,double rwidth,double *msq,double *wt);

double GenGaus(double* rand,Double_t mean, Double_t sigma);

int kt_pdf(double r0, double* ktX);	

bool Acceptance(event_type);

void testPdf(double);

void My_choose(TVar::Process process);
bool My_eventcuts(TVar::Process process, mcfm_event_type*, cdf_event_type*);
bool My_masscuts(double s[][12],TVar::Process process);
bool My_smalls(double s[][12], int npart);

bool Reject(TLorentzVector);
bool TriggerLepton(TLorentzVector);
double ProbConversion(TLorentzVector);
double ProbFake(TLorentzVector);


double SumMatrixElementPDF(TVar::Process procees, mcfm_event_type* mcfm_event,double flavor_msq[][11],double* flux);

double HiggsWidth(double);
double SetTGCParameter(TString ,double par);
double GetTGCParameter(int i);
double getProbAcceptanceEfficiency(cdf_event_type cdf_event, EffHist effhist);
void getProbFromHist(double x0, double* kX, double *wt, TH1F *hkx);

#endif
