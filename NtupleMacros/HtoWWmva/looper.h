#ifndef looper_h
#define looper_h

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "CORE/jetSelections.h"

class TChain;

class looper
{
    public:
        looper() {};
        ~looper() {
            delete babyFile_;
            delete eventTree_;
            //delete trackTree_;
        };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, const char*, bool isData, int nEvents = -1);
        void bookHistos();
        void fillUnderOverFlow(TH1F *h1, float value, float weight = 1);
        void printEvent(  ostream& ostr = cout );
        float deltaPhi( float phi1 , float phi2 );
        bool isGoodTrack( int index );
        bool isMuon( int index );
        bool isElectron( int index );
	bool isGoodVertex(size_t ivtx) ;
	unsigned int nGoodVertex() ;
	std::vector<unsigned int> getJetIdVector(unsigned int i_hyp, float etMin, float etMax, float etaMax, bool doBtag=false, float discCut=2.1, float vetoCone=0.3);
	unsigned int getHardestJetId(std::vector<unsigned int> jetIdVec);
	unsigned int get2ndHardJetId(std::vector<unsigned int> jetIdVec);
	float getCorrectedJetPt(unsigned int jetId);
	unsigned int getNExtraLeptons(unsigned int i_hyp, double minLepPt = 10.);
	unsigned int getNSoftMu(unsigned int i_hyp, std::vector<unsigned int> jetIdVec);

    private:
        
        //ntuple, file
        TFile *babyFile_;
        TTree *eventTree_;

	//jet correction tool
	FactorizedJetCorrector* jet_corrector_pf;

        //histos
        TH1F* h_dPhi_ll;
        TH1F* h_dEta_ll;
        TH1F* h_mass_ll;
        TH1F* h_minPt_ll;
        TH1F* h_maxPt_ll;

        // event stuff
        Int_t   sample_id_;

        Int_t   event_evt_;
        Int_t   event_run_;
        Int_t   event_lumi_;

        Float_t event_type_;

	Float_t event_xsec_incl_;
	Float_t event_xsec_excl_;
	Float_t event_kfactor_;
	Float_t event_scale1fb_;

	Float_t dil_mass_;
	Float_t dil_pt_;
	Float_t dil_eta_;
	Float_t dil_phi_;
	Float_t dil_dphi_;

	Int_t   lephard_q_;
	Int_t   lephard_id_;
        Float_t lephard_pt_;
        Float_t lephard_eta_;
        Float_t lephard_phi_;
        Float_t lephard_hOverE_;
        Float_t lephard_dEtaIn_;
        Float_t lephard_dPhiIn_;
        Float_t lephard_sigmaIEtaIEta_;
        Float_t lephard_e2x5Max_;
        Float_t lephard_e1x5_;
        Float_t lephard_e5x5_;
        Float_t lephard_eSC_;
        Float_t lephard_etaSC_;
        Float_t lephard_eOverPIn_;
        Float_t lephard_eOverPOut_;
        Float_t lephard_fbrem_;
	Int_t   lephard_genId_;
	Int_t   lephard_genMotherId_;
        Float_t lephard_mva_;
        Float_t lephard_newconv_dist_;
        Float_t lephard_newconv_dcot_;
        Float_t lephard_newconv_rad_;
        Float_t lephard_newconv_dmh_;

	Int_t   lepsoft_q_;
	Int_t   lepsoft_id_;
        Float_t lepsoft_pt_;
        Float_t lepsoft_eta_;
        Float_t lepsoft_phi_;
        Float_t lepsoft_hOverE_;
        Float_t lepsoft_dEtaIn_;
        Float_t lepsoft_dPhiIn_;
        Float_t lepsoft_sigmaIEtaIEta_;
        Float_t lepsoft_e2x5Max_;
        Float_t lepsoft_e1x5_;
        Float_t lepsoft_e5x5_;
        Float_t lepsoft_eSC_;
        Float_t lepsoft_etaSC_;
        Float_t lepsoft_eOverPIn_;
        Float_t lepsoft_eOverPOut_;
        Float_t lepsoft_fbrem_;
	Int_t   lepsoft_genId_;
	Int_t   lepsoft_genMotherId_;
	Int_t   lepsoft_passTighterId_;
        Float_t lepsoft_fr_;
        Float_t lepsoft_mva_;
        Float_t lepsoft_newconv_dist_;
        Float_t lepsoft_newconv_dcot_;
        Float_t lepsoft_newconv_rad_;
        Float_t lepsoft_newconv_dmh_;
        
	Float_t met_pt_;
	Float_t met_phi_;
	Float_t met_projpt_;

        Int_t   jets_num_;
        Int_t   lowptbtags_num_;
        Int_t   btags_num_;
        Int_t   extralep_num_;
        Int_t   softmu_num_;

	Float_t jethard_pt_;
	Float_t jethard_eta_;
	Float_t jethard_phi_;
	Float_t jethard_disc_;
	Float_t jethard2_pt_;
	Float_t jethard2_eta_;
	Float_t jethard2_phi_;
	Float_t jethard2_disc_;
	Float_t ctrjethard_pt_;
	Float_t ctrjethard_eta_;
	Float_t ctrjethard_phi_;
	Float_t ctrjethard_disc_;

	Float_t jets_discmax_;

	Float_t met_sumet_;
	Float_t dil_metdphi_;
	Float_t dil_deta_;
	Float_t dil_dr_;

	Float_t mt_lephardmet_;
	Float_t mt_lepsoftmet_;
	Float_t dphi_lephardmet_;
	Float_t dphi_lepsoftmet_;
	Float_t mt_dilmet_;

	Float_t llm_sumpt_;
	Float_t llmj_sumpt_;

	Float_t gen_lp_pt_;
	Float_t gen_lp_eta_;
	Float_t gen_lp_phi_;
	Int_t   gen_lp_id_;
	Float_t gen_lm_pt_;
	Float_t gen_lm_eta_;
	Float_t gen_lm_phi_;
	Int_t   gen_lm_id_;
	Float_t gen_np_pt_;
	Float_t gen_np_eta_;
	Float_t gen_np_phi_;
	Int_t   gen_np_id_;
	Float_t gen_nm_pt_;
	Float_t gen_nm_eta_;
	Float_t gen_nm_phi_;
	Int_t   gen_nm_id_;

        ofstream ofile;
};

#endif
