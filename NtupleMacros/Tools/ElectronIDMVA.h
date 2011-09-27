//--------------------------------------------------------------------------------------------------
// $Id $
//
// ElectronIDMVA
//
// Helper Class for Electron Identification MVA
//
// Authors: S.Xie
// Original based on MitPhysics/Utils/interface/ElectronIDMVA.h?view=markup
// Modified by DLE
//--------------------------------------------------------------------------------------------------

#ifndef ELECTRONIDMVA_H
#define ELECTRONIDMVA_H

#include "TString.h"

class TRandom3;
namespace TMVA {
    class Reader;
}

class ElectronIDMVA {
    public:
        ElectronIDMVA();
        ~ElectronIDMVA(); 

        void    Initialize(TString methodName,
                TString Subdet0Pt10To20Weights , 
                TString Subdet1Pt10To20Weights , 
                TString Subdet2Pt10To20Weights,
                TString Subdet0Pt20ToInfWeights, 
                TString Subdet1Pt20ToInfWeights, 
                TString Subdet2Pt20ToInfWeights);

        Bool_t   IsInitialized() const { return fIsInitialized; }
        Double_t MVAValue(const unsigned int ele, const unsigned int vertex);


    protected:      
        TMVA::Reader            *fTMVAReader[6];
        TString                  fMethodname;
        Bool_t                    fIsInitialized;

        Float_t                   fMVAVar_EleSigmaIEtaIEta; 
        Float_t                   fMVAVar_EleDEtaIn; 
        Float_t                   fMVAVar_EleDPhiIn; 
        Float_t                   fMVAVar_EleHoverE; 
        Float_t                   fMVAVar_EleD0; 
        Float_t                   fMVAVar_EleDZ; 
        Float_t                   fMVAVar_EleFBrem; 
        Float_t                   fMVAVar_EleEOverP; 
        Float_t                   fMVAVar_EleESeedClusterOverPout; 
        Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
        Float_t                   fMVAVar_EleNBrem; 
        Float_t                   fMVAVar_EleOneOverEMinusOneOverP; 
        Float_t                   fMVAVar_EleESeedClusterOverPIn; 
        Float_t                   fMVAVar_EleIP3d; 
        Float_t                   fMVAVar_EleIP3dSig; 
        Float_t                   fMVAVar_EleStandardLikelihood; 


};

#endif

