// ROOT includes
#include "Math/LorentzVector.h"

// C++ includes
#include <map>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;

class Tauify{
  private:
    // data member to hold tau rest frame variables
    // pid, |P|, costheta
    // not as complicated as it looks... think: index <--> vector(int, float, float)
    map<int, pair<int, pair<double, double> > > tau_data;

    LorentzVector p4_lepton_lab;
    LorentzVector p4_tau_lepton_cm;
    LorentzVector p4_tau_lepton_lab;

  public:
    // constructor & destructor
    Tauify( const char*, bool = false );
    ~Tauify(){};

    // Set methods
    void SetLepton( LorentzVector, float, float, float );  // lepton p4, met, iso, d0

    // Get methods
    LorentzVector TauP4(void);
    float TauMET(void);
    float TauIso(void);
    float TauIP(void);

    // io test
    unsigned int TauSize(void);
    int   First(int);
    double Second(int);
    double Third(int);

};
