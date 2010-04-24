
#ifndef ELECTRONSELECTIONSCICPARAMETERS_H
#define ELECTRONSELECTIONSCICPARAMETERS_H

//
// Data required for V02 of 
// Branson/Sani electron ID
//
// Settings provided
//
//--------------------------------
//eidVeryLoose
//eidLoose
//eidMedium
//eidTight
//eidSuperTight
//eidHyperTight1
//eidHyperTight2
//eidHyperTight3
//eidHyperTight4
//--------------------------------

#include <vector>

enum cic_tightness {
    CIC_VERYLOOSE,
    CIC_LOOSE,
    CIC_MEDIUM,
    CIC_TIGHT,
    CIC_SUPERTIGHT,
    CIC_HYPERTIGHT1,
    CIC_HYPERTIGHT2,
    CIC_HYPERTIGHT3,
    CIC_HYPERTIGHT4,
};

void eidGet(const cic_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

void eidAssign(std::vector<double> &cutarr, double cutvals[]);

void eidVeryLoose(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee); 

void eidLoose(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee); 

void eidMedium(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

void eidTight(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

void eidSuperTight(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

void eidHyperTight1(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee); 

void eidHyperTight2(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee); 

void eidHyperTight3(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee); 

void eidHyperTight4(std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

#endif
 
