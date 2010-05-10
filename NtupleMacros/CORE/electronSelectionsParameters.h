
#ifndef ELECTRONSELECTIONSPARAMETERS_H
#define ELECTRONSELECTIONSPARAMETERS_H

#include <vector>

//
// Data required for V02 of 
// Branson/Sani electron ID
//
// Settings provided
//
//--------------------------------
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
//--------------------------------

//
// Data required for VBTF ID
// optimised for us by N. Rompotis/C. Seez
//
// Settings provided
//
//--------------------------------
enum vbtf_tightness {
    VBTF_TOP95,
    VBTF_TOP90,
    VBTF_TOP85,
    VBTF_TOP80,
    VBTF_TOP70,
    VBTF_TOP55,
};
//--------------------------------

//
// Data required for CAND ID
//
// Settings provided
//
//--------------------------------
enum cand_tightness {
    CAND_01,
    CAND_02,
};
//--------------------------------

void eidGetCand(const cand_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuthoe, std::vector<double> &cutslat);

void eidGetVBTF(const vbtf_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuthoe, std::vector<double> &cutsee, std::vector<double> &cutreliso);

void eidGetCIC(const cic_tightness tightness, std::vector<double> &cutdeta, std::vector<double> &cutdphi, std::vector<double> &cuteopin, std::vector<double> &cutet, std::vector<double> &cuthoe, std::vector<double> &cutip, std::vector<double> &cutisoecal, std::vector<double> &cutisohcal, std::vector<double> &cutisotk, std::vector<double> &cutmishits, std::vector<double> &cutsee);

void eidAssign(std::vector<double> &cutarr, double cutvals[], unsigned int size);

#endif

