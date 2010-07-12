// ROOT includes
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class Tauify{
  private:
    // data member to hold tau rest frame variables
    // pid, |P|, costheta

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

};
