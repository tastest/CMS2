#ifndef WW_analysisTools_h
#define WW_analysisTools_h

#include "TH2D.h"
#include "TH1D.h"

#include "../../../Smurf/Core/SmurfTree.h"
#include "analysisEnums.h"

#include <vector>

#ifndef __CINT__
#include "../CORE/CMS2.h"
#endif

//
// Process utilities
//

bool isDYee();
bool isDYmm();
bool isDYtt();
bool isWW();
bool isWZ();
bool isZZ();
unsigned int getDrellYanType();
unsigned int getVVType();

// filter events by process
bool filterByProcess( enum SmurfTree::DataType sample );
bool isIdentified( enum SmurfTree::DataType sample );

//
// event utilities
//

class CMS2;
class EventIdentifier {
  unsigned long int run, event, lumi;
  float trks_d0;
  bool data;
 public:
  EventIdentifier(CMS2& cms2, bool isData);
  bool operator < (const EventIdentifier &) const;
  bool operator == (const EventIdentifier &) const;
};

static std::set<EventIdentifier> already_seen;
bool is_duplicate (const EventIdentifier &id);

//
// cut utilities
//

// N-1
// return true if the cuts to apply - the cuts to remove were passed in the cuts that "passed"
bool CheckCutsNM1(wwcuts_t apply, wwcuts_t remove, wwcuts_t passed);
// Simple check if the desired cuts to apply are set in the cuts that "passed"
bool CheckCuts(wwcuts_t apply, wwcuts_t passed);

//
// weights
//

float getHiggsPt();
float getHiggsPtWeight(float pt, int mH, TH1D *HiggsPtKFactor);
float getZPt();
float getZRapidity();
float getZMass();
float getDYNNLOWeight(float pt, float rap, float mass, vector<TH2D*> fDYNNLOKFactorHists);
double mt(double pt1, double pt2, double dphi);

#endif

