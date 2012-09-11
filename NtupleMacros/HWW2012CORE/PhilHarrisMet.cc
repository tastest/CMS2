
#include "PhilHarrisMet.h"

#include "../CORE/CMS2.h"
#include "analysisSelections.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

PhilHarrisMet::PhilHarrisMet() :
    isInitialized_(false)
{
    reader_ = 0;
}

PhilHarrisMet::~PhilHarrisMet()
{
    if (reader_) delete reader_;
}

void PhilHarrisMet::initialize(TString weights)
{

    TMVA::Tools::Instance();
    reader_ = new TMVA::Reader( "!Color:!Silent" );    
    lJetName_   = "JetID";
    isInitialized_ = true;

    lNPV_    = 0;
    lJPt1_   = 0; lJEta1_  = 0; lJPhi1_  = 0;  lJM1_ = 0;  lJD01_ = 0; lJDZ1_ = 0;  lNPart1_ = 0; 
    lLPt1_   = 0; lLEta1_  = 0; lLPhi1_  = 0; 
    lSPt1_   = 0; lSEta1_  = 0; lSPhi1_  = 0; 
    lNEPt1_  = 0; lNEEta1_ = 0; lNEPhi1_ = 0; 
    lEMPt1_  = 0; lEMEta1_ = 0; lEMPhi1_ = 0; 
    lChPt1_  = 0; lChEta1_ = 0; lChPhi1_ = 0;  lLFr1_   = 0;
    lDRlC1_  = 0; lDRLS1_  = 0; lDRM1_   = 0;  lDRMNE1_ = 0; lDREM1_ = 0; lDRCH1_ = 0;

    reader_->AddSpectator( "npv"   , &lNPV_    ); 
    reader_->AddVariable( "jspt_1"   , &lJPt1_   ); 
    reader_->AddVariable( "jseta_1"  , &lJEta1_  );
    reader_->AddVariable( "jsphi_1"  , &lJPhi1_  );                                                                                 
    reader_->AddVariable( "jd0_1"    , &lJD01_   );
    reader_->AddVariable( "jdZ_1"    , &lJDZ1_   );
    reader_->AddVariable( "jm_1"     , &lJM1_    );
    reader_->AddVariable( "npart_1"  , &lNPart1_ );
    reader_->AddVariable( "lpt_1"    , &lLPt1_   );                                                                            
    reader_->AddVariable( "leta_1"   , &lLEta1_  );
    reader_->AddVariable( "lphi_1"   , &lLPhi1_  );                                                                                 
    reader_->AddVariable( "spt_1"    , &lSPt1_   );                                                                                  
    reader_->AddVariable( "seta_1"   , &lSEta1_  );
    reader_->AddVariable( "sphi_1"   , &lSPhi1_  );                                                                               
    reader_->AddVariable( "lnept_1"  , &lNEPt1_  );                                                                              
    reader_->AddVariable( "lneeta_1" , &lNEEta1_ );
    reader_->AddVariable( "lnephi_1" , &lNEPhi1_ );                                                                            
    reader_->AddVariable( "lempt_1"  , &lEMPt1_  );                                                                                   
    reader_->AddVariable( "lemeta_1" , &lEMEta1_ );
    reader_->AddVariable( "lemphi_1" , &lEMPhi1_ );                                                                              
    reader_->AddVariable( "lchpt_1"  , &lChPt1_  );                                                                                   
    reader_->AddVariable( "lchphi_1" , &lChPhi1_ );                                                                               
    reader_->AddVariable( "lLfr_1"   , &lLFr1_   );
    reader_->AddVariable( "drlc_1"   , &lDRlC1_  );
    reader_->AddVariable( "drls_1"   , &lDRLS1_  );
    reader_->AddVariable( "drm_1"    , &lDRM1_   );
    reader_->AddVariable( "drmne_1"  , &lDRMNE1_ );
    reader_->AddVariable( "drem_1"   , &lDREM1_  );
    reader_->AddVariable( "drch_1"   , &lDRCH1_  );
    reader_->BookMVA(lJetName_, weights);

}

float PhilHarrisMet::value(const unsigned int jet)
{

/*
00039       X=0,     // undefined
00040       h,       // charged hadron
00041       e,       // electron 
00042       mu,      // muon 
00043       gamma,   // photon
00044       h0,      // neutral hadron
00045       h_HF,        // HF tower identified as a hadron
00046       egamma_HF    // HF tower identified as an EM particle
*/

    // number of vertices
    lNPV_    = nGoodVertex();
    
    // jet kinematics
    // FIXME - need to apply corrections
    lJPt1_  = cms2.pfjets_p4()[jet].Pt();
    lJEta1_ = cms2.pfjets_p4()[jet].Eta();
    lJPhi1_ = cms2.pfjets_p4()[jet].Phi();
    lJM1_ = cms2.pfjets_p4()[jet].M();
    lJD01_ = 99999;
    lJDZ1_ = 99999;
    lNPart1_ = cms2.pfjets_pfcandIndicies()[jet].size();
    
    // leading pt cand
    int leadIndex = leadCand(jet, -1, false);
    lLPt1_   = cms2.pfcands_p4()[leadIndex].Pt(); 
    lLEta1_  = cms2.pfcands_p4()[leadIndex].Eta(); 
    lLPhi1_  = cms2.pfcands_p4()[leadIndex].Phi();

    // sub leading pt cand
    int subIndex = leadCand(jet, -1, true);
    lSPt1_   = cms2.pfcands_p4()[subIndex].Pt();
    lSEta1_  = cms2.pfcands_p4()[subIndex].Eta();
    lSPhi1_  = cms2.pfcands_p4()[subIndex].Phi();

    // leading neutral hadron cand
    int leadh0 = leadCand(jet, 5, false);
    lNEPt1_   = cms2.pfcands_p4()[leadh0].Pt();
    lNEEta1_  = cms2.pfcands_p4()[leadh0].Eta();
    lNEPhi1_  = cms2.pfcands_p4()[leadh0].Phi();

    // leading em cand
    int leadem = leadCand(jet, 4, false);
    lEMPt1_   = cms2.pfcands_p4()[leadem].Pt();
    lEMEta1_  = cms2.pfcands_p4()[leadem].Eta();
    lEMPhi1_  = cms2.pfcands_p4()[leadem].Phi();

    // leading charged hadron cand
    int leadch = leadCand(jet, 1, false);
    lChPt1_   = cms2.pfcands_p4()[leadch].Pt();
    lChEta1_  = cms2.pfcands_p4()[leadch].Eta();
    lChPhi1_  = cms2.pfcands_p4()[leadch].Phi();

    // dr leading cand to jet axis
    lDRlC1_ = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4()[leadIndex], cms2.pfjets_p4()[jet]);

    // dr leading to second cand
    lDRLS1_ = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4()[leadIndex], cms2.pfcands_p4()[subIndex]);

    // mean delta r weighted moment of all candidates
    lDRM1_   = dRMean(jet, -1);  

    // mean delta r weighted moment of neutral hadrons
    lDRMNE1_ = dRMean(jet, 5); 

    // mean delta r weighted moment of em cands
    lDREM1_ = dRMean(jet, 4); 

    // mean delta r weighted moment of charged hadrons
    lDRCH1_ = dRMean(jet, 1);

    //
    // compure and return the mva output
    //

    float mva = reader_->EvaluateMVA(lJetName_);    
    return mva;

}

int PhilHarrisMet::leadCand(unsigned int jet, int iPFType, bool i2nd) 
{ 
    int lCount = 0;
    int index = -1;
    const std::vector<int> candIndices =  cms2.pfjets_pfcandIndicies()[jet];
    for (unsigned int i = 0; i < candIndices.size(); ++i) {
        index = candIndices[i];
        if (cms2.pfcands_particleId()[index] != -1 
                 && cms2.pfcands_particleId()[index] != iPFType) continue;
        if(lCount == 0 && !i2nd)    break;
        if(lCount >  0)             break;
        lCount++;
    }
    return index; 
}

double PhilHarrisMet::dRMean(unsigned int jet, int iPFType) 
{ 
    double lDRMean = 0;
    const std::vector<int> candIndices =  cms2.pfjets_pfcandIndicies()[jet];
    for (unsigned int i = 0; i < candIndices.size(); ++i) {
        int index = candIndices[i];
        if (cms2.pfcands_particleId()[index] != -1 && cms2.pfcands_particleId()[index] != iPFType) continue;
        double pDR = ROOT::Math::VectorUtil::DeltaR(cms2.pfcands_p4()[index], cms2.pfjets_p4()[jet]);
        lDRMean += (pDR * cms2.pfcands_p4()[index].Pt()) / cms2.pfjets_p4()[jet].Pt();
    }
  return lDRMean;
}

double PhilHarrisMet::impactParameter(unsigned int jet, bool iDZ) 
{ 

/*
    double lDZCorr = -1000;
    for(UInt_t i0 = 0; i0 < iJet->NPFCands(); i0++) { 
        const PFCandidate *pCand = iJet->PFCand(i0);
        if(pCand->BestTrk() == 0) continue;
        if(pCand->Pt() < 1.) continue;
        if(iDZ)  lDZCorr = pCand->BestTrk()->DzCorrected(*fVertex);
        if(!iDZ) lDZCorr = pCand->BestTrk()->D0Corrected(*fVertex);
        break;
  }
  return lDZCorr;
*/
    return 999.9;
}

