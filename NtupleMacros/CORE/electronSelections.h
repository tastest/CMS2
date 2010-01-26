
#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

//
// electron selections 
//

//
// candidate electron id function
//
bool electronId_cand01();

//
// class based electron id that we have used before
//
bool electronId_classBased();

//
// candidate electron isolation function
//
bool electronIsolation_cand01();

//
// functions for individual subdet electron isolations
//
bool electronTkIsolation(const unsigned int index, const float threshold);
bool electronEcalIsolation(const unsigned int index, const float threshold);
bool electronHcalIsolation(const unsigned int index, const float threshold);

#endif

