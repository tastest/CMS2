#include "ElectronIDMVA.h"

#include "../CORE/CMS2.h"
#include "../CORE/electronSelections.h"
#include "../CORE/trackSelections.h"

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
void ElectronIDMVA::Initialize( TString methodName, unsigned int version,
        TString Subdet0Pt10To20Weights , 
        TString Subdet1Pt10To20Weights , 
        TString Subdet2Pt10To20Weights,
        TString Subdet0Pt20ToInfWeights,
        TString Subdet1Pt20ToInfWeights, 
        TString Subdet2Pt20ToInfWeights) {

    if (version != 1 && version != 2 && version != 3) {
        std::cout << "[ElectronIDMVA::Initialize] Version must be 1 or 2 or 3.  Aborting." << std::endl;
        return;
    }

    version_ = version;

    fIsInitialized = kTRUE;
    fMethodname = methodName;

    for(UInt_t i=0; i<6; ++i) {
        if (fTMVAReader[i]) delete fTMVAReader[i];

	fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
	fTMVAReader[i]->SetVerbose(kTRUE);
        // order matters!

        if (version == 1) {
            fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );  
            fTMVAReader[i]->SetVerbose(kTRUE);
            fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
            fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
            fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
            fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                 );
            fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                );
            fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout  );
            fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi         );
            fTMVAReader[i]->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                 );
            fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP );
            fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn   );
        } else if (version == 2) {
            fTMVAReader[i] = new TMVA::Reader( "!Color:!Silent:Error" );
            fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta         );
            fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                );
            fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                );
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
        } else if (version == 3) {
	    fTMVAReader[i]->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
	    fTMVAReader[i]->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
	    fTMVAReader[i]->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
	    //fTMVAReader[i]->AddVariable( "HoverE",                &fMVAVar_EleHoverE                   );
	    fTMVAReader[i]->AddVariable( "D0",                    &fMVAVar_EleD0                       );
	    fTMVAReader[i]->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
	    fTMVAReader[i]->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
	    fTMVAReader[i]->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
	    fTMVAReader[i]->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
	    fTMVAReader[i]->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
	    fTMVAReader[i]->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
	    fTMVAReader[i]->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
	    fTMVAReader[i]->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );
	    
	    fTMVAReader[i]->AddVariable( "GsfTrackChi2OverNdof",  &fMVAVar_EleGsfTrackChi2OverNdof     );
	    fTMVAReader[i]->AddVariable( "dEtaCalo",              &fMVAVar_EledEtaCalo                 );
	    fTMVAReader[i]->AddVariable( "dPhiCalo",              &fMVAVar_EledPhiCalo                 );
	    fTMVAReader[i]->AddVariable( "R9",                    &fMVAVar_EleR9                       );
	    fTMVAReader[i]->AddVariable( "SCEtaWidth",            &fMVAVar_EleSCEtaWidth               );
	    fTMVAReader[i]->AddVariable( "SCPhiWidth",            &fMVAVar_EleSCPhiWidth               );
	    fTMVAReader[i]->AddVariable( "CovIEtaIPhi",           &fMVAVar_EleCovIEtaIPhi              );
	    if (i == 2 || i == 5) {
	      fTMVAReader[i]->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
	    }
	    fTMVAReader[i]->AddVariable( "ChargedIso03",          &fMVAVar_EleChargedIso03OverPt       );
	    fTMVAReader[i]->AddVariable( "NeutralHadronIso03",    &fMVAVar_EleNeutralHadronIso03OverPt );
	    fTMVAReader[i]->AddVariable( "GammaIso03",            &fMVAVar_EleGammaIso03OverPt         );
	    fTMVAReader[i]->AddVariable( "ChargedIso04",          &fMVAVar_EleChargedIso04OverPt       );
	    fTMVAReader[i]->AddVariable( "NeutralHadronIso04",    &fMVAVar_EleNeutralHadronIso04OverPt );
	    fTMVAReader[i]->AddVariable( "GammaIso04",            &fMVAVar_EleGammaIso04OverPt         );
	}

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
    if (version_ == 3) fMVAVar_EleHoverE =                 cms2.els_hOverE()[ele]; //this is new
    fMVAVar_EleD0 =                     electron_d0PV_wwV1(ele);
    fMVAVar_EleDZ =                     electron_dzPV_wwV1(ele);
    fMVAVar_EleFBrem =                  cms2.els_fbrem()[ele]; 
    fMVAVar_EleEOverP =                 cms2.els_eOverPIn()[ele];
    fMVAVar_EleESeedClusterOverPout =   cms2.els_eSeedOverPOut()[ele];
    fMVAVar_EleSigmaIPhiIPhi =          cms2.els_sigmaIPhiIPhi()[ele];
    fMVAVar_EleNBrem =                  cms2.els_nSeed()[ele];
    TVector3 pIn(cms2.els_trk_p4()[ele].px(), cms2.els_trk_p4()[ele].py(), cms2.els_trk_p4()[ele].pz());
    fMVAVar_EleOneOverEMinusOneOverP =  1./cms2.els_eSC()[ele] - 1./pIn.Mag();
    fMVAVar_EleESeedClusterOverPIn =    cms2.els_eSeedOverPIn()[ele];
    const double gsfsign   = ( (gsftrks_d0_pv(cms2.els_gsftrkidx().at(ele),0).first)   >=0 ) ? 1. : -1.;
    fMVAVar_EleIP3d =                   cms2.els_ip3d()[ele]*gsfsign; 
    if (cms2.els_ip3derr()[ele] == 0.0) fMVAVar_EleIP3dSig = 0.0;
    else fMVAVar_EleIP3dSig =           cms2.els_ip3d()[ele]*gsfsign / cms2.els_ip3derr()[ele]; 

    if (version_ == 3) {    //these are new
      Double_t ElePt = cms2.els_p4()[ele].pt();
      Double_t EleEta = cms2.els_etaSC()[ele];
      Double_t Rho = cms2.evt_ww_rho_vor();
      fMVAVar_EleGsfTrackChi2OverNdof = cms2.els_chi2()[ele] / cms2.els_ndof()[ele];
      fMVAVar_EledEtaCalo = cms2.els_dEtaOut()[ele]; 
      fMVAVar_EledPhiCalo = cms2.els_dPhiOut()[ele];
      fMVAVar_EleR9 = cms2.els_e3x3()[ele] / cms2.els_eSCRaw()[ele];
      fMVAVar_EleSCEtaWidth = cms2.els_etaSCwidth()[ele];
      fMVAVar_EleSCPhiWidth = cms2.els_phiSCwidth()[ele];
      fMVAVar_EleCovIEtaIPhi = cms2.scs_sigmaIEtaIPhi()[cms2.els_scindex()[ele]]>0 ? pow(cms2.scs_sigmaIEtaIPhi()[cms2.els_scindex()[ele]],2) : -1.* pow(cms2.scs_sigmaIEtaIPhi()[cms2.els_scindex()[ele]],2);
      fMVAVar_ElePreShowerOverRaw = cms2.els_eSCPresh()[ele] / cms2.els_eSCRaw()[ele];
      fMVAVar_EleChargedIso03OverPt 
	= (cms2.els_iso03_pf_ch()[ele] //electronIsoValuePF(ele, vertex, 0.3, 99999., 0.1, 0.07, 0.025, -999., 0)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleChargedIso03, EleEta)) / ElePt;
      fMVAVar_EleNeutralHadronIso03OverPt 
	= (cms2.els_iso03_pf_nhad05()[ele] //electronIsoValuePF(ele, vertex, 0.3, 0.5, 0.1, 0.07, 0.025, -999., 130)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso03, EleEta) 
	   + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso007,EleEta)) / ElePt;
      fMVAVar_EleGammaIso03OverPt 
	= (cms2.els_iso03_pf_gamma05()[ele] //electronIsoValuePF(ele, vertex, 0.3, 0.5, 0.1, 0.07, 0.025, -999., 22)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIso03, EleEta) 
	   + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIsoVetoEtaStrip03,EleEta))/ElePt;      
      fMVAVar_EleChargedIso04OverPt 
	= (cms2.els_iso04_pf_ch()[ele] //electronIsoValuePF(ele, vertex, 0.4, 99999., 0.1, 0.07, 0.025, -999., 0)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleChargedIso04, EleEta)) / ElePt;
      fMVAVar_EleNeutralHadronIso04OverPt 
	= (cms2.els_iso04_pf_nhad05()[ele] //electronIsoValuePF(ele, vertex, 0.4, 0.5, 0.1, 0.07, 0.025, -999., 130)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso04, EleEta) 
	   + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso007,EleEta)) / ElePt;
      fMVAVar_EleGammaIso04OverPt 
	= (cms2.els_iso04_pf_gamma05()[ele] //electronIsoValuePF(ele, vertex, 0.4, 0.5, 0.1, 0.07, 0.025, -999., 22)*ElePt
	   - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIso04, EleEta) 
	   + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIsoVetoEtaStrip04,EleEta))/ElePt;
    }

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

    
    //DEBUG
    if (0) {
      cout << endl;
      cout << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << " " << cms2.evt_ww_rho_vor() << endl;
      std::cout << "Debug Electron MVA: "
		<< cms2.els_p4()[ele].pt() << " " << " " << cms2.els_p4()[ele].Eta() << " " << cms2.els_p4()[ele].Phi() << " : "
		<< cms2.els_p4()[ele].pt() << " " << cms2.els_etaSC()[ele] << " --> MVABin " << MVABin << " : "     
		<< fMVAVar_EleSigmaIEtaIEta << " " 
		<< fMVAVar_EleDEtaIn << " " 
		<< fMVAVar_EleDPhiIn << " " 
		<< fMVAVar_EleHoverE << " " 
		<< fMVAVar_EleD0 << " " 
		<< fMVAVar_EleDZ << " " 
		<< fMVAVar_EleFBrem << " " 
		<< fMVAVar_EleEOverP << " " 
		<< fMVAVar_EleESeedClusterOverPout << " " 
		<< fMVAVar_EleSigmaIPhiIPhi << " " 
		<< fMVAVar_EleNBrem << " " 
		<< fMVAVar_EleOneOverEMinusOneOverP << " " 
		<< fMVAVar_EleESeedClusterOverPIn << " " 
		<< fMVAVar_EleIP3d << " " 
		<< fMVAVar_EleIP3dSig << " " 
		<< fMVAVar_EleGsfTrackChi2OverNdof << " "
		<< fMVAVar_EledEtaCalo << " "
		<< fMVAVar_EledPhiCalo << " "
		<< fMVAVar_EleR9 << " "
		<< fMVAVar_EleSCEtaWidth << " "
		<< fMVAVar_EleSCPhiWidth << " "
		<< fMVAVar_EleCovIEtaIPhi << " "
		<< fMVAVar_ElePreShowerOverRaw << " "
		<< fMVAVar_EleChargedIso03OverPt  << " "
		<< fMVAVar_EleNeutralHadronIso03OverPt  << " "
		<< fMVAVar_EleGammaIso03OverPt  << " "
		<< fMVAVar_EleChargedIso04OverPt  << " "
		<< fMVAVar_EleNeutralHadronIso04OverPt  << " "
		<< fMVAVar_EleGammaIso04OverPt  << " "
		<< " === : === "
		<< mva << " "    
		<< std::endl;
    }
    return mva;

}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleSCEta,
        Double_t EleSigmaIEtaIEta,
        Double_t EleDEtaIn,
        Double_t EleDPhiIn,
        Double_t EleD0,
        Double_t EleDZ,
        Double_t EleFBrem,
        Double_t EleEOverP,
        Double_t EleESeedClusterOverPout,
        Double_t EleSigmaIPhiIPhi,
        Double_t EleNBrem,
        Double_t EleOneOverEMinusOneOverP,
        Double_t EleESeedClusterOverPIn,
        Double_t EleIP3d,
        Double_t EleIP3dSig
        ) {

    if (!fIsInitialized) { 
        std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
        return -9999;
    }

    Int_t subdet = 0;
    if (fabs(EleSCEta) < 1.0) subdet = 0;
    else if (fabs(EleSCEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (ElePt > 20.0) ptBin = 1;

    //set all input variables
    fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
    fMVAVar_EleDEtaIn = EleDEtaIn;
    fMVAVar_EleDPhiIn = EleDPhiIn;
    fMVAVar_EleD0 = EleD0;
    fMVAVar_EleDZ = EleDZ;
    fMVAVar_EleFBrem = EleFBrem;
    fMVAVar_EleEOverP = EleEOverP;
    fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
    fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
    fMVAVar_EleNBrem = EleNBrem;
    fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
    fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
    fMVAVar_EleIP3d = EleIP3d;
    fMVAVar_EleIP3dSig = EleIP3dSig;

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
//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleEta, Double_t PileupEnergyDensity,
                                 Double_t EleSigmaIEtaIEta,
                                 Double_t EleDEtaIn,
                                 Double_t EleDPhiIn,
                                 Double_t EleHoverE,
                                 Double_t EleD0,
                                 Double_t EleDZ,
                                 Double_t EleFBrem,
                                 Double_t EleEOverP,
                                 Double_t EleESeedClusterOverPout,
                                 Double_t EleSigmaIPhiIPhi,
                                 Double_t EleNBrem,
                                 Double_t EleOneOverEMinusOneOverP,
                                 Double_t EleESeedClusterOverPIn,
                                 Double_t EleIP3d,
                                 Double_t EleIP3dSig,
                                 Double_t EleGsfTrackChi2OverNdof,
                                 Double_t EledEtaCalo,
                                 Double_t EledPhiCalo,
                                 Double_t EleR9,
                                 Double_t EleSCEtaWidth,
                                 Double_t EleSCPhiWidth,
                                 Double_t EleCovIEtaIPhi,
                                 Double_t ElePreShowerOverRaw,
                                 Double_t EleChargedIso03,
                                 Double_t EleNeutralHadronIso03,
                                 Double_t EleGammaIso03,
                                 Double_t EleChargedIso04,
                                 Double_t EleNeutralHadronIso04,
                                 Double_t EleGammaIso04,
                                 Bool_t printDebug
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Double_t Rho = 0;
  if (!(TMath::IsNaN(PileupEnergyDensity) || isinf(PileupEnergyDensity))) Rho = PileupEnergyDensity;

  Int_t subdet = 0;
  if (fabs(EleEta) < 1.0) subdet = 0;
  else if (fabs(EleEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ElePt > 20.0) ptBin = 1;
  
  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleDZ = EleDZ;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;
  fMVAVar_EleGsfTrackChi2OverNdof = EleGsfTrackChi2OverNdof;
  fMVAVar_EledEtaCalo = EledEtaCalo;
  fMVAVar_EledPhiCalo = EledPhiCalo;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleSCEtaWidth = EleSCEtaWidth;
  fMVAVar_EleSCPhiWidth = EleSCPhiWidth;
  fMVAVar_EleCovIEtaIPhi = EleCovIEtaIPhi;
  fMVAVar_ElePreShowerOverRaw = ElePreShowerOverRaw;
  fMVAVar_EleChargedIso03OverPt 
    = (EleChargedIso03 
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleChargedIso03, EleEta)) / ElePt;
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (EleNeutralHadronIso03
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso03, EleEta) 
       + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso007,EleEta)) / ElePt;
  fMVAVar_EleGammaIso03OverPt 
    = (EleGammaIso03 
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIso03, EleEta) 
       + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIsoVetoEtaStrip03,EleEta))/ElePt;
  fMVAVar_EleChargedIso04OverPt 
    = (EleChargedIso04 
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleChargedIso04, EleEta))/ElePt;
  fMVAVar_EleNeutralHadronIso04OverPt
    = (EleNeutralHadronIso04 
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso04, EleEta) 
       + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleNeutralHadronIso007,EleEta))/ElePt;
  fMVAVar_EleGammaIso04OverPt 
    = (EleGammaIso04 
       - Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIso04, EleEta) 
       + Rho * ElectronEffectiveArea(ElectronIDMVA::kEleGammaIsoVetoEtaStrip04,EleEta))/ElePt;
  
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

  if (printDebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
	 << ElePt << " " << EleEta << " " << " --> MVABin " << MVABin << " : "     
	 << fMVAVar_EleSigmaIEtaIEta << " " 
	 << fMVAVar_EleDEtaIn << " " 
	 << fMVAVar_EleDPhiIn << " " 
	 << fMVAVar_EleHoverE << " " 
	 << fMVAVar_EleD0 << " " 
	 << fMVAVar_EleDZ << " " 
	 << fMVAVar_EleFBrem << " " 
	 << fMVAVar_EleEOverP << " " 
	 << fMVAVar_EleESeedClusterOverPout << " " 
	 << fMVAVar_EleSigmaIPhiIPhi << " " 
	 << fMVAVar_EleNBrem << " " 
	 << fMVAVar_EleOneOverEMinusOneOverP << " " 
	 << fMVAVar_EleESeedClusterOverPIn << " " 
	 << fMVAVar_EleIP3d << " " 
	 << fMVAVar_EleIP3dSig << " " 
	 << fMVAVar_EleGsfTrackChi2OverNdof << " "
	 << fMVAVar_EledEtaCalo << " "
	 << fMVAVar_EledPhiCalo << " "
	 << fMVAVar_EleR9 << " "
	 << fMVAVar_EleSCEtaWidth << " "
	 << fMVAVar_EleSCPhiWidth << " "
	 << fMVAVar_EleCovIEtaIPhi << " "
	 << fMVAVar_ElePreShowerOverRaw << " "
	 << fMVAVar_EleChargedIso03OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso03OverPt  << " "
	 << fMVAVar_EleGammaIso03OverPt  << " "
	 << fMVAVar_EleChargedIso04OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso04OverPt  << " "
	 << fMVAVar_EleGammaIso04OverPt  << " "
	 << " === : === "
	 << mva 
	 << std::endl;
  }

  return mva;
}

Double_t ElectronIDMVA::ElectronEffectiveArea(EElectronEffectiveAreaType type, Double_t Eta) {

  Double_t EffectiveArea = 0;

  if (fabs(Eta) < 1.0) {
    if (type == ElectronIDMVA::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso03) EffectiveArea = 0.017;
    if (type == ElectronIDMVA::kEleGammaIso03) EffectiveArea = 0.045;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.014;
    if (type == ElectronIDMVA::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso04) EffectiveArea = 0.034;
    if (type == ElectronIDMVA::kEleGammaIso04) EffectiveArea = 0.079;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.014;
    if (type == ElectronIDMVA::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleHoverE) EffectiveArea = 0.00016;
    if (type == ElectronIDMVA::kEleHcalDepth1OverEcal) EffectiveArea = 0.00016;
    if (type == ElectronIDMVA::kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == ElectronIDMVA::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso03) EffectiveArea = 0.025;
    if (type == ElectronIDMVA::kEleGammaIso03) EffectiveArea = 0.052;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.030;
    if (type == ElectronIDMVA::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso04) EffectiveArea = 0.050;
    if (type == ElectronIDMVA::kEleGammaIso04) EffectiveArea = 0.073;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.030;
    if (type == ElectronIDMVA::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleHoverE) EffectiveArea = 0.00022;
    if (type == ElectronIDMVA::kEleHcalDepth1OverEcal) EffectiveArea = 0.00022;
    if (type == ElectronIDMVA::kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == ElectronIDMVA::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso03) EffectiveArea = 0.030;
    if (type == ElectronIDMVA::kEleGammaIso03) EffectiveArea = 0.170;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.134;
    if (type == ElectronIDMVA::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso04) EffectiveArea = 0.060;
    if (type == ElectronIDMVA::kEleGammaIso04) EffectiveArea = 0.187;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.134;
    if (type == ElectronIDMVA::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleHoverE) EffectiveArea = 0.00030;
    if (type == ElectronIDMVA::kEleHcalDepth1OverEcal) EffectiveArea = 0.00026;
    if (type == ElectronIDMVA::kEleHcalDepth2OverEcal) EffectiveArea = 0.00002;        
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == ElectronIDMVA::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso03) EffectiveArea = 0.022;
    if (type == ElectronIDMVA::kEleGammaIso03) EffectiveArea = 0.623;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.516;
    if (type == ElectronIDMVA::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso04) EffectiveArea = 0.055;
    if (type == ElectronIDMVA::kEleGammaIso04) EffectiveArea = 0.659;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.517;
    if (type == ElectronIDMVA::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleHoverE) EffectiveArea = 0.00054;
    if (type == ElectronIDMVA::kEleHcalDepth1OverEcal) EffectiveArea = 0.00045;
    if (type == ElectronIDMVA::kEleHcalDepth2OverEcal) EffectiveArea = 0.00003;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.5 ) {
    if (type == ElectronIDMVA::kEleChargedIso03) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso03) EffectiveArea = 0.018;
    if (type == ElectronIDMVA::kEleGammaIso03) EffectiveArea = 1.198;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip03) EffectiveArea = 1.049;
    if (type == ElectronIDMVA::kEleChargedIso04) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleNeutralHadronIso04) EffectiveArea = 0.073;
    if (type == ElectronIDMVA::kEleGammaIso04) EffectiveArea = 1.258;
    if (type == ElectronIDMVA::kEleGammaIsoVetoEtaStrip04) EffectiveArea = 1.051;
    if (type == ElectronIDMVA::kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == ElectronIDMVA::kEleHoverE) EffectiveArea = 0.00082;
    if (type == ElectronIDMVA::kEleHcalDepth1OverEcal) EffectiveArea = 0.00066;
    if (type == ElectronIDMVA::kEleHcalDepth2OverEcal) EffectiveArea = 0.00004;
  }
    
  return EffectiveArea;  
}

