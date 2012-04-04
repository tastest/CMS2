#ifndef MSUGRAXSEC_H
#define MSUGRAXSEC_H

#include <fstream>
#include <iostream>
#include "TH1.h"

double getMsugraCrossSection( double my_m0 , double my_m12, double my_tanb , bool verbose = false );
void set_msugra_file ( const char* filename = "goodModelNames_tanbeta10.txt" , bool verbose = false );
ifstream ifile;

//SMS stuff
enum sms_process { gg, ss }; //gluino-gluino, squark-squark--that's all there is for now

TH1F* sms_gluino_xsec_hist;
TH1F* sms_squark_xsec_hist;

float getSMSCrossSection( const float mgluino , const sms_process type=gg );
void  set_sms_xsec_hist ( const char* filename, const sms_process type=gg );

#endif

