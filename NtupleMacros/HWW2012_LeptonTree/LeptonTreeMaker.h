#ifndef WW_LeptonTreeMaker_h
#define WW_LeptonTreeMaker_h

#include "Math/LorentzVector.h"
#include "Rtypes.h"
#include <vector>
#include <set>
#include "../HWW2012CORE/wwtypes.h"
#include "TChain.h"
#include <fstream>
#include <vector>

#include "../../../Smurf/Core/LeptonTree.h"
#include "../../../Smurf/Core/SmurfTree.h"

#include "../Tools/ElectronIDMVA.h"
#include "../CORE/jetcorr/FactorizedJetCorrector.h"
#include "../CORE/trackSelections.h"

#include "../HWW2012CORE/analysisEnums.h"
#include "SmurfDataTypes.h"


class TChain;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVectorD;

class LeptonTreeMaker {

    public:

        LeptonTreeMaker() {};
        LeptonTreeMaker(bool lockToCoreSelectors, bool useLHeleId, 
                bool useMVAeleId, bool doDYNNLOw, unsigned int prescale, bool realData);
        ~LeptonTreeMaker();

        void SetBatchMode(bool on) {
            batchMode_ = on;
        }

        void ScanChain(TString outfileid, 
                TChain* chain,
				SmurfTree::DataType sample,
                double integratedLumi,
                double xsec,
                int nProcessedEvents,
                bool identifyEvents = false,
                bool realData = false,
                TString cms2_json_file = "");

    private:

        bool loosefo(unsigned int index);

        //
        // efficiencies
        //

        void MakeElectronTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);
        void MakeMuonTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

        // this one is a little bit special
        void MakeElectronRecoTagAndProbeTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

        //
        // fake rates
        //

        void MakeElectronFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection);
        void MakeMuonFakeRateTree(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample, const unsigned int eventSelection);

        //
        // utilities
        //

        float GetAwayJetPt(LorentzVector lep1, LorentzVector lep2);

        //
        // common variables
        //

        void SetCommonTreeVariables(LeptonTree &leptonTree, const double &weight, SmurfTree::DataType sample);

        //
        // data members 
        //

        bool lockToCoreSelectors_;
        bool useLHeleId_;
        bool useMVAeleId_;
        bool doDYNNLOw_;
        unsigned int prescale_;
        bool realData_;
        bool batchMode_;

        // jet correction
        FactorizedJetCorrector *jet_corrector_pfL1FastJetL2L3_;

        // ID MVAas
        ElectronIDMVA *electronIdMVA_;

		// triggers
		// electron
		unsigned int HLT_Ele27_WP80_tag_;
		unsigned int HLT_Ele27_WP80_probe_;
		unsigned int HLT_Ele17_Ele8_tag_;				
		unsigned int HLT_Ele17_Ele8_probe_;
		unsigned int HLT_Ele17_Ele8_LeadingLeg_tag_; 	
		unsigned int HLT_Ele17_Ele8_LeadingLeg_probe_;
		unsigned int HLT_Ele17_Ele8_TrailingLeg_tag_;	
		unsigned int HLT_Ele17_Ele8_TrailingLeg_probe_;
		// muon	
		unsigned int HLT_IsoMu24_eta2p1_tag_;			
		unsigned int HLT_IsoMu24_eta2p1_probe_;
		unsigned int HLT_Mu17_Mu8_tag_;				
		unsigned int HLT_Mu17_Mu8_probe_;		
		unsigned int HLT_Mu17_Mu8_LeadingLeg_tag_;		
		unsigned int HLT_Mu17_Mu8_LeadingLeg_probe_;		
		unsigned int HLT_Mu17_Mu8_TrailingLeg_tag_;	
		unsigned int HLT_Mu17_Mu8_TrailingLeg_probe_;
		unsigned int HLT_Mu17_TkMu8_tag_;				
		unsigned int HLT_Mu17_TkMu8_probe_;	
		unsigned int HLT_Mu17_TkMu8_LeadingLeg_tag_;	
		unsigned int HLT_Mu17_TkMu8_LeadingLeg_probe_;	
		unsigned int HLT_Mu17_TkMu8_TrailingLeg_tag_;	
		unsigned int HLT_Mu17_TkMu8_TrailingLeg_probe_;
		// FR trigger	
		unsigned int HLT_Ele8_CaloIdL_TrkIdVL_;
		unsigned int HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_;
		unsigned int HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_;
		unsigned int HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_;
		unsigned int HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_;
		unsigned int HLT_Mu8_;
		unsigned int HLT_Mu17_;

		// test
		float mva_;
		float iso_;
		float fbrem_; 
		float kfchi2_;
		int   kfhits_;
		float gsfchi2_;
		float deta_;
		float dphi_;
		float detacalo_;
		float see_;
		float spp_;
		float etawidth_;
		float phiwidth_;
		float e1x5e5x5_;
		float R9_;
		float HoE_;
		float EoP_;
		float IoEmIoP_;
		float eleEoPout_;
		float PreShowerOverRaw_;
		float d0_;
		float ip3d_;
		float eta_;
		float pt_;
};

#endif
