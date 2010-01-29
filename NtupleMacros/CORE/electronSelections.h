
#ifndef ELECTRONSELECTIONS_H
#define ELECTRONSELECTIONS_H

//
// electron selections 
//

//
// candidate electron id function
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/ElectronID
//
bool electronId_cand01(const unsigned int index);
bool electronImpact_cand01(const unsigned int index);

//
// candidate electron isolation function
// see http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/EfficiencyMeasurement
//
bool electronIsolation_cand01(const unsigned int index);

//
// class based electron id that we have used before
//
bool electronId_classBasedLoose(const unsigned int index);
bool electronId_classBasedTight(const unsigned int index);

//
// electron isolation definitions that we have used before
//
float electronIsolation_rel(const unsigned int index, bool use_calo_iso);
float electronIsolation_relsusy(const unsigned int index, bool use_calo_iso);
float electronIsolation_relsusy_cand0(const unsigned int index, bool use_calo_iso);

//
// conversion rejection
//
bool isFromConversionHitPattern(const unsigned int index);
bool isFromConversionPartnerTrack(const unsigned int index);

#endif

