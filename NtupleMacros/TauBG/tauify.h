// ROOT includes
#include "Math/LorentzVector.h"
#include <map>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

class Tauify{
  private:
    // data member to hold tau rest frame variables
    // pid, |P|, costheta
    // not as complicated as it looks... think: index <--> vector(int, float, float)
    map<int, pair<int, pair<float, float> > > tau_data;

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
    float Second(int);
    float Third(int);

};
