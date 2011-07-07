#include <fstream>

#ifndef MSUGRAXSEC_H
#define MSUGRAXSEC_H

double getMsugraCrossSection( double my_m0 , double my_m12, double my_tanb , bool verbose = false );
void set_msugra_file ( const char* filename = "goodModelNames_tanbeta10.txt" , bool verbose = false );
ifstream ifile;

#endif

