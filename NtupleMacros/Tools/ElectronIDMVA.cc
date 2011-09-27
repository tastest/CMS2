
#include "ElectronIDMVA.h"

#include "../CORE/CMS2.h"
#include "../CORE/electronSelections.h"

#include "TFile.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//--------------------------------------------------------------------------------------------------
ElectronIDMVA::ElectronIDMVA() :
    fMethodname("BDTG method"),
    fIsInitialized(kFALSE)
{
    // Constructor.
    for(UInt_t i=0; i<6; ++i) {
        fTMVAReader[i] = 0;
    }
}


//--------------------------------------------------------------------------------------------------
ElectronIDMVA::~ElectronIDMVA()
{
    for(UInt_t i=0; i<6; ++i) {
        if (fTMVAReader[i]) delete fTMVAReader[i];
    }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize( TString methodName,
        TString Subdet0Pt10To20Weights , 
        TString Subdet1Pt10To20Weights , 
        TString Subdet2Pt10To20Weights,
        TString Subdet0Pt20ToInfWeights,
        TString Subdet1Pt20ToInfWeights, 
        TString Subdet2Pt20ToInfWeights) {

    fIsInitialized = kTRUE;
    fMethodname = methodName;

    for(UInt_t i=0; i<6; ++i) {
        if (fTMVAReader[i]) delete fTMVAReader[i];

        fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
        fTMVAReader[i]->SetVerbose(kTRUE);
        fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
        fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
        fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
        fTMVAReader[i]->AddVariable( "HoverE",                &fMVAVar_EleHoverE                );
        fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                    );
        fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
        fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
        fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
        fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
        fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
        fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );
        fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
        fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                  );
        fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig               );
        fTMVAReader[i]->AddVariable( "StandardLikelihood",    &fMVAVar_EleStandardLikelihood    );

        if (i==0) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt10To20Weights );
        if (i==1) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt10To20Weights );
        if (i==2) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt10To20Weights );
        if (i==3) fTMVAReader[i]->BookMVA(fMethodname , Subdet0Pt20ToInfWeights );
        if (i==4) fTMVAReader[i]->BookMVA(fMethodname , Subdet1Pt20ToInfWeights );
        if (i==5) fTMVAReader[i]->BookMVA(fMethodname , Subdet2Pt20ToInfWeights );

    }

    std::cout << "Electron ID MVA Initialization\n";
    std::cout << "MethodName : " << fMethodname << std::endl;
    std::cout << "Load weights file : " << Subdet0Pt10To20Weights << std::endl;
    std::cout << "Load weights file : " << Subdet1Pt10To20Weights << std::endl;
    std::cout << "Load weights file : " << Subdet2Pt10To20Weights << std::endl;
    std::cout << "Load weights file : " << Subdet0Pt20ToInfWeights << std::endl;
    std::cout << "Load weights file : " << Subdet1Pt20ToInfWeights << std::endl;
    std::cout << "Load weights file : " << Subdet2Pt20ToInfWeights << std::endl;

}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const unsigned int ele, const unsigned int vertex) {

    if (!fIsInitialized) { 
        std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
        return -9999;
    }

    Int_t subdet = 0;
    if (fabs(cms2.els_etaSC()[ele]) < 1.0) subdet = 0;
    else if (fabs(cms2.els_etaSC()[ele]) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (cms2.els_p4()[ele].Pt() > 20.0) ptBin = 1;

    //set all input variables
    fMVAVar_EleSigmaIEtaIEta =          cms2.els_sigmaIEtaIEta()[ele]; 
    fMVAVar_EleDEtaIn =                 cms2.els_dEtaIn()[ele]; 
    fMVAVar_EleDPhiIn =                 cms2.els_dPhiIn()[ele]; 
    fMVAVar_EleHoverE =                 cms2.els_hOverE()[ele];
    fMVAVar_EleD0 =                     electron_d0PV_wwV1(ele);
    fMVAVar_EleDZ =                     electron_dzPV_wwV1(ele);
    fMVAVar_EleFBrem =                  cms2.els_fbrem()[ele]; 
    fMVAVar_EleEOverP =                 cms2.els_eOverPIn()[ele];
    fMVAVar_EleESeedClusterOverPout =   cms2.els_eSeedOverPOut()[ele];
    fMVAVar_EleSigmaIPhiIPhi =          cms2.els_sigmaIPhiIPhi()[ele];
    fMVAVar_EleNBrem =                  cms2.els_nSeed()[ele];
    TVector3 pIn(cms2.els_trk_p4()[ele].px(), cms2.els_trk_p4()[ele].py(), cms2.els_trk_p4()[ele].pz());
    fMVAVar_EleOneOverEMinusOneOverP =  1./(cms2.els_eOverPIn()[ele]*pIn.Mag()) - 1./pIn.Mag();
    fMVAVar_EleESeedClusterOverPIn =    cms2.els_eSeedOverPIn()[ele];
    fMVAVar_EleIP3d =                   cms2.els_ubIp3d()[ele]; 
    if (cms2.els_ubIp3derr()[ele] == 0.0) fMVAVar_EleIP3dSig = 0.0;
    else fMVAVar_EleIP3dSig =           cms2.els_ubIp3d()[ele] / cms2.els_ubIp3derr()[ele]; 
    fMVAVar_EleStandardLikelihood =     cms2.els_lh()[ele];

    Double_t mva = -9999;  
    TMVA::Reader *reader = 0;
    Int_t MVABin = -1;
    if (subdet == 0 && ptBin == 0) MVABin = 0;
    if (subdet == 1 && ptBin == 0) MVABin = 1;
    if (subdet == 2 && ptBin == 0) MVABin = 2;
    if (subdet == 0 && ptBin == 1) MVABin = 3;
    if (subdet == 1 && ptBin == 1) MVABin = 4;
    if (subdet == 2 && ptBin == 1) MVABin = 5;
    assert(MVABin >= 0 && MVABin <= 5);
    reader = fTMVAReader[MVABin];

    mva = reader->EvaluateMVA( fMethodname );
    return mva;

}

