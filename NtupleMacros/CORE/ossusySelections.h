#ifndef OSSUSYSELECTIONS_H
#define OSSUSYSELECTIONS_H

#include <vector>
#include "Math/LorentzVector.h"


typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

/*****************************************************************************************/
//passes the OS SUSY trigger selections
/*****************************************************************************************/
bool passTrigger_OSV1( bool isData , int hypType );

/*****************************************************************************************/
//hypothesis disambiguation. Returns the hypothesis that has mass closent to MZ
/*****************************************************************************************/
unsigned int selectBestZHyp(const vector<unsigned int> &v_goodHyps);

#endif

