#ifndef MT2_H
#define MT2_H

#include "MT2Utility.cc"
#include "Math/LorentzVector.h"
#include "TMath.h"

enum enum_mt2_method { BISECT, GRID };

// MT2 Container
class TMt2 {
  private:
    double Mt2_;
    double Nu1Px_;
    double Nu1Py_;
    double Nu2Px_;
    double Nu2Py_;
  public:
    TMt2( double, double, double, double, double );
    inline double Mt2   (void) { return Mt2_;   }
    inline double Nu1Px (void) { return Nu1Px_; }
    inline double Nu1Py (void) { return Nu1Py_; }
    inline double Nu2Px (void) { return Nu2Px_; }
    inline double Nu2Py (void) { return Nu2Py_; }
};
TMt2::TMt2( double Mt2_arg, double Nu1Px_arg, double Nu1Py_arg, double Nu2Px_arg, double Nu2Py_arg ){
  Mt2_   = Mt2_arg;
  Nu1Px_ = Nu1Px_arg;
  Nu1Py_ = Nu1Py_arg;
  Nu2Px_ = Nu2Px_arg;
  Nu2Py_ = Nu2Py_arg;
}

// MT2 declaration
double MT2(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  float = 0.0,
  bool = false
);

TMt2* MT2_GRID(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  int = 500,
  int = 1,
  float = 0.0,
  bool = false
);

// MT2J declaration
double MT2J(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  float = 0.0,
  enum_mt2_method = BISECT,
  bool = false
);

/*
double MT2J_GRID(
  const float,
  const float,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  float = 0.0,
  bool = false
);
*/

vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > CleanJets(
  const vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >,
  const float,
  const float,
  const float
);

#endif

