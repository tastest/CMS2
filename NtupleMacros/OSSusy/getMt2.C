#include "CORE/MT2/MT2Utility.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

double getMt2(LorentzVector v1, LorentzVector v2, double met, double metPhi)
{
    double pa[3];
    double pb[3];
    double pmiss[3];
    pa[0] = v1.mass();
    pa[1] = v1.Px();
    pa[2] = v1.Py();
    pb[0] = v2.mass();
    pb[1] = v2.Px();
    pb[2] = v2.Py();
    pmiss[0] = 0.;
    pmiss[1] = met*cos(metPhi);
    pmiss[2] = met*sin(metPhi);
    mt2_bisect::mt2 mt2_event;
    mt2_event.set_momenta( pa, pb, pmiss );
    mt2_event.set_mn(0.0);
    double mt2_value = mt2_event.get_mt2(); 
    return mt2_value;
}
