#ifndef WW_KINEMATIC
#define WW_KINEMATIC
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

void WWL1L2Sol_Mw1( event_type* temp
                    ,double Mw1, double nu_X, double nu_Y, double nu_Z
                    ,event_type* sol,double* jacobian);


void WWL1L2Sol_Mw1Mw2( event_type* temp
                     ,double Mw1, double Mw2, double nuZ, double nbZ
                     ,event_type* sol,double* jacobian);
void WWL1L2Sol_Mw1Mw2_exact( event_type* temp
                     ,double Mw1, double Mw2, double nuZ, double nbZ
                     ,event_type* sol,double* jacobian);

void WWL1L2Sol_MHiggs( event_type* temp
                     ,double qX, double qY
                     ,double MHiggs, double nu_X, double nu_Y, double nu_Z
                     ,event_type* sol,double* jacobian);
void WWL1L2Sol_Mw1MHiggs( event_type* temp
		        ,double Mw1, double MHiggs, double nu_X, double nu_Y
	                ,event_type* sol,double* jacobian);

void WWL1L2Sol(double theta,double phi,event_type* temp, double* J);




void WWL1L2Sol_MHiggsMw1( cdf_event_type* temp,
                    double qx, double qy, 
                    double MHiggs, double Mw1, double nuX, double nuY,
                    mcfm_event_type* sol);

void WWL1L2Sol_MHiggs( cdf_event_type* temp,
                    double qX, double qY , 
		    double MHiggs_sqr, double nuX, double nuY, double nuZ,
                    mcfm_event_type* sol);

void WWL1L2Sol_MHiggsYHiggs(cdf_event_type* temp,
                    double qX, double qY , 
		    double MHiggs, double YHiggs, double nuX, double nuY,
                    mcfm_event_type* sol);

void WWL1L2Sol_Mw1Mw2( cdf_event_type* temp,
                    double qx, double qy , double Mw1, double Mw2, double nuZ, double nbZ,
                    mcfm_event_type* sol);

void WWL1L2Sol_Mw( cdf_event_type* temp,
                   double  qX, double   qY, double Msqr,
                   mcfm_event_type* sol);

void WWL1L2Sol_MzNu3( cdf_event_type* temp,
                      double qx, double qy,
                      double Mz, double nu_X, double nu_Y, double nu_Z,
                      mcfm_event_type* sol);

#endif
