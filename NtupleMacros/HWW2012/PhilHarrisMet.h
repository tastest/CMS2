#ifndef PHILHARRISMET_H
#define PHILHARRISMET_H

#include "TString.h"

namespace TMVA {
    class Reader;
}               

class PhilHarrisMet {

    public:     

        PhilHarrisMet(); 
        ~PhilHarrisMet(); 

        void initialize(TString weights);
        bool isInitialized() const { return isInitialized_; }
        float value(const unsigned int jet);

    protected:

        //
        // MVA
        //

        TMVA::Reader    *reader_;
        Bool_t          isInitialized_;
        TString         lJetName_;

        //
        // data used in MVA
        //

        int   lNPV_    ;
        float lJPt1_   ; float lJEta1_  ; float lJPhi1_  ;  float lJM1_ ;  float lJD01_ ; float lJDZ1_ ;   float lNPart1_ ; 
        float lLPt1_   ; float lLEta1_  ; float lLPhi1_  ; 
        float lSPt1_   ; float lSEta1_  ; float lSPhi1_  ; 
        float lNEPt1_  ; float lNEEta1_ ; float lNEPhi1_ ; 
        float lEMPt1_  ; float lEMEta1_ ; float lEMPhi1_ ; 
        float lChPt1_  ; float lChEta1_ ; float lChPhi1_ ;  float lLFr1_   ;
        float lDRlC1_  ; float lDRLS1_  ; float lDRM1_   ;  float lDRMNE1_ ; float lDREM1_ ; float lDRCH1_ ;

        //
        // util functions
        //

        int     leadCand         (unsigned int jet, int iPFType, bool i2nd);
        double  dRMean           (unsigned int jet, int iPFType);
        double  impactParameter  (unsigned int jet, bool iDZ);

};

#endif

