// C++ includes
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>

// ROOT includes
#include "Math/AxisAngle.h"
#include "Math/Boost.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"

//
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;

class Tauify{

  private:
    // data member to hold tau rest frame variables
    // pid, |P|, costheta
    // not as complicated as it looks... think: index <--> vector(int, double, double)
    map<int, pair<int, pair<double, double> > > tau_data;

    LorentzVector p4_lepton_lab;
    LorentzVector p4_lepton_cm;
    LorentzVector p4_tau_lepton_cm;
    LorentzVector p4_tau_lepton_lab;

    int    id;
    double p_cm;
    double costheta_cm;

    float tau_met;
    float tau_iso;
    float tau_ip;

  public:
    // constructor & destructor
    Tauify( const char*, bool = false );
    ~Tauify(){};

    // Set methods
    void SetLepton( LorentzVector, float, float, float );  // lepton p4, met, iso, d0

    // Get methods
    LorentzVector TauP4(void);

    int   ParticleId(void);
    float MomentumCM(void);
    float CosThetaCM(void);
    float TauMET(void);
    float TauIso(void);
    float TauIP(void);

    // io test
    unsigned int TauSize(void);
    int    ParticleId(int);
    double MomentumCM(int);
    double CosThetaCM(int);

};
