#ifndef topAFB_looper_h
#define topAFB_looper_h

#include <stdint.h>
#include "TFile.h"
#include "TTree.h"

#define NCHANNELS 1
#define NHISTS 1
#define JETPTCUT 30.0
// top mass                                                                                                                                                       
#include "../CORE/topmass/ttdilepsolve.cpp" 
#include "../CORE/topmass/getTopMassEstimate.icc" 
#include "../CORE/jetSelections.cc"
#include "../CORE/jetSmearingTools.h"
//class TChain;

class topAFB_looper
{
    public:
        topAFB_looper();
        ~topAFB_looper();
	enum FREnum   { e_qcd = 0, e_wjets };
        void MakeBabyNtuple (const char *);
        void InitBabyNtuple ();
        void FillBabyNtuple ();
        void CloseBabyNtuple ();
        void ScanChain (TChain*, vector<TString> , string, 
			bool, float, float, bool, FREnum frmode = e_wjets, double AMWTmass = 172.5);
	void bookHistos(const char* sample, int nchannels, int nhistsets);
	void FillHistograms(const unsigned int hypIdx, const vector<unsigned int> v_jets, const vector<unsigned int> v_jetsNoEtaCut,
					 const pair<float, float> p_met, const float weight, const string prefix) ;
	bool passbTagging(const unsigned int jet_idx, const string jetAlgo, const string bTagDiscriminator) ;
	double getFRWeight(const int hypIdx, SimpleFakeRate *mufr, SimpleFakeRate *elfr, FREnum frmode, bool isData);
	double getBFRWeight(const int hypIdx, 	vector<LorentzVector> & v_goodNonBtagJets_p4,vector<LorentzVector> & v_goodBtagJets_p4, bool isData);
	double TopPtWeight(double topPt);
	double LeptonPtWeight(double leptonPt);
	double JetPtWeight(double jetPt);
	double JERsf(double eta);
	void fillUnderOverFlow(TH1D *h1, float value, double weight = 1., int Nsolns=1);
	void fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight = 1., int Nsolns=1);
	//void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
	void fillOverFlow(TH1D *h1, float value, float weight = 1.);
	void fillOverFlow(TH2D *h2, float xvalue, float yvalue, float weight = 1.);
	void fillHistos(TH1D *h1[4][4],float value, double weight, int myType, int nJetsIdx, int Nsolns=1);
	void fillHistos(TH2D *h2[4][4],float xvalue, float yvalue, double weight, int myType, int nJetsIdx, int Nsolns=1);
	void fillHistos(TProfile *h2[4][4],float xvalue, float yvalue,  int myType, int nJetsIdx);
	
	int antimatch4vector(const LorentzVector &lvec, 
					    const vector<LorentzVector> &vec);
	double triggerEff(const int hypIdx, bool scaleTrigSFup = false, bool scaleTrigSFdown = false);
	
    private:
        //
        // BABY NTUPLE VARIABLES
        //
	JetSmearer* jet_smearer_;
        TFile *babyFile_;
        TTree *babyTree_;
	bool applyNoCuts;
	bool getVtxDistOnly;
	bool usePtGt2020;
	bool usePtGt2010;
	bool excludePtGt2020;
	bool applylepIDCuts;
	bool applyFOv1Cuts;
	bool applyFOv2Cuts;
	bool applyFOv3Cuts;
	bool applylepIsoCuts;
	bool applylepLooseIsoCuts;
	bool applyTriggers;
	bool vetoZmass;
	bool requireZmass;
	bool hypDisamb;
	bool useCorMET;
	bool usetcMET;
	bool usetcMET35X;
	bool usepfMET;
	bool vetoMET;
	bool vetoMET50;
	bool vetoProjectedMET;
	bool usejptJets;
	bool usecaloJets;
	bool usepfJets;
	bool veto1Jet;
	bool veto2Jets;
	bool requireEcalEls;
	bool useOS;
	bool useSS;
  bool applyTopPtWeighting;
  bool applyLeptonPtWeighting;
  bool applyJetPtWeighting;
  bool useReweightingUncorrelated;
  bool useReweightingLeadingObject;
	bool applyAlignmentCorrection;
	bool vetoHypMassLt10;
	bool vetoHypMassLt12;
	bool scaleJESMETUp; 
	bool scaleJESMETDown; 
	bool scaleJER; 
  bool scaleLeptonEnergyUp;
  bool scaleLeptonEnergyDown;
  bool scaleBTAGSFup;
  bool scaleBTAGSFdown;
  bool scaleTrigSFup;
  bool scaleTrigSFdown;
  bool noVertexReweighting;
  bool weighttaudecay;
	bool doJESmear; 
	bool estimateQCD;
	bool estimateWJets;
	bool requireBTag;
	bool require2BTag;
	bool sortJetCandidatesbyPt;
	bool sortJetCandidatesbyDR;
	bool applyLeptonJetInvMassCut450;
	bool applyTopSystEta;
	bool requireExact2BTag;
	bool applyHTCut;
	bool matchLeptonJetbyMaxDR;
	bool generalLeptonVeto;
	bool BTagAlgTCHE;
	bool createBabyNtuples;
	bool doBFR;
	bool applyMinMassLBCut;
	float globalJESRescale;
	FactorizedJetCorrector *jptL2L3Corr;
	FactorizedJetCorrector *pfL2L3Corr;
	FactorizedJetCorrector *jptL2L3ResidualCorr;
	FactorizedJetCorrector *pfL2L3ResidualCorr;
	ttdilepsolve *d_llsol;
	
        // event identification

        Int_t   run_;
        Int_t   ls_;
        Int_t   evt_;
	Int_t ndavtx_;
	float   t_mass_;
	double weight_;
	Int_t Nsolns_;
	float massltb_;
	float massllb_;
	float dr_ltjet_gen_;
	float dr_lljet_gen_;
	float tt_mass_ ;
	float ttRapidity_ ;
	float ttRapidity2_ ;
	float ttPt_ ;
	float lep_charge_asymmetry_ ;
	float lep_pseudorap_diff_ ;
	float lep_azimuthal_asymmetry_ ;
	float lep_azimuthal_asymmetry2_ ;
	float top_spin_correlation_ ;
	float top_costheta_cms_    ;
	float top_rapiditydiff_cms_    ;
	float top_rapiditydiff_Marco_    ;
	float top_pseudorapiditydiff_cms_    ;
	float lepPlus_costheta_cms_ ;
	float lepMinus_costheta_cms_ ;
	float tt_mass_gen_ ;
	float ttRapidity_gen_ ;
	float ttRapidity2_gen_ ;
	float ttPt_gen_ ;
	float lep_charge_asymmetry_gen_ ;
	float lep_azimuthal_asymmetry_gen_ ;
	float lep_azimuthal_asymmetry_gen2_ ;
	float top_spin_correlation_gen_ ;
	float  top_costheta_cms_gen_    ;
	float lepPlus_costheta_cms_gen_ ;
	float lepMinus_costheta_cms_gen_ ;
	float top_rapiditydiff_cms_gen_    ;
	float top_rapiditydiff_Marco_gen_    ;
	float top_pseudorapiditydiff_cms_gen_;  
	float massllbb_ ;
	float massllbb_gen_;
	//	float llbbRapidityQuark_gen_;
	//	float llbbRapidityGluon_gen_;
	TH1D* hnJet[4][4];                   // Njet distributions                                                                                                          
	TH1D* hnBtagJet[4][4];                   // NBTagjet distributions  
	TH1D* hnVtx[4][4];
	TH1D* hNsolns[4][4];
	TH1D* hmaxAMWTweight[4][4];
	TH1D* haveAMWTweight[4][4];
	TH1D* hsumAMWTweight[4][4];
	TH1D* hAMWTweightnojetsmear[4][4];          
	//Top Mass Plots
	TH1D *httRapidity[4][4];
	TH1D *httRapidity2[4][4];     
	TH1D *hllbbMass[4][4];
	TH1D *hllbbMass_gen[4][4];
	TH1D *httpT[4][4]; 
	TH1D *httMass[4][4]; 
	TH1D *httMass_pull[4][4]; 
	TH1D *httMass_diff[4][4];
	TH1D *htopMass_diff_plus[4][4];
	TH1D *htopMass_diff_minus[4][4];                                                                                                                                           
	TH1D *htopPCM_diff_plus[4][4];
	TH1D *htopPCM_diff_minus[4][4];    

	TH1D *htopMass[4][4];
	TH1D *httpT_nojetsmear[4][4]; 
	TH1D *httMass_nojetsmear[4][4]; 
	TH1D *httMass_nojetsmear_pull[4][4]; 
	TH1D *httMass_nojetsmear_diff[4][4];                                                                                                                                            
	TH1D *htopMass_nojetsmear_diff_plus[4][4];
	TH1D *htopMass_nojetsmear_diff_minus[4][4];   
	TH1D *htopMass_nojetsmear[4][4];
	TH1D *httpT_gen[4][4]; 
	TH1D *httMass_gen[4][4]; 
	TH1D *htopMass_plus_gen[4][4]; 
	TH1D *htopMass_minus_gen[4][4]; 
	TH1D *hmassllb[4][4];
	TH1D *hmassltb[4][4];
	TH1D *hmassllb1Dmasscut[4][4];
	TH1D *hmassltb1Dmasscut[4][4];
	TH1D *htheSumJetPt[4][4];
	TH1D *htheSumBtagJetPt[4][4];
	TH1D *hthefirstJetPt[4][4];
	TH1D *hthesecondJetPt[4][4];
	TH1D *htheleadinglepPt[4][4];
	TH1D *hthesecondlepPt[4][4];
	TH1D *hthesumlepPt[4][4];
        TH1D *hlepEta[4][4];
	TH1D *hjetPt[4][4];
	TH1D *hjetEta[4][4];
	TH1D *hMET[4][4];
	TH1D *htheSumLepPt[4][4];
	TH1D *htopCosTheta[4][4];
	TH1D *hpseudorapiditydiff[4][4];
	TH1D *hrapiditydiff[4][4];
	TH1D *hlepCosTheta[4][4];
	TH1D *hlepChargeAsym[4][4];
	TH1D *hlepAzimAsym[4][4];
	TH1D *hlepAzimAsym2[4][4];
	TH1D *htopSpinCorr[4][4];
	TH1D *htopCosTheta_gen[4][4];
	TH1D *hlepCosTheta_gen[4][4];
	TH1D *hlepChargeAsym_gen[4][4];
	TH1D *hlepAzimAsym_gen[4][4];
	TH1D *hlepAzimAsym2_gen[4][4];
	TH1D *htopSpinCorr_gen[4][4];
	
	TH1D *hlepRapDiff[4][4];
	TH1D *hlepAngleBetween[4][4];
	TH1D *hlepAngleBetweenCMS[4][4];
	TH1D *hpseudorapiditydiff2[4][4];
	TH1D *hrapiditydiff2[4][4];
	TH1D *hlepPlusCosTheta[4][4];
	TH1D *hlepMinusCosTheta[4][4];
	TH1D *hjetAzimAsym[4][4];
	TH1D *hjetRapDiff[4][4];
	TH1D *hjetAngleBetween[4][4];
	TH1D *hjetAngleBetweenCMS[4][4];
	TH1D *hlepPhi[4][4];
	TH1D *hlepPlusPhi[4][4];
	TH1D *hlepMinusPhi[4][4];
	TH1D *hjetPhi[4][4];
	TH1D *hlepPlusEta[4][4];
	TH1D *hlepMinusEta[4][4];
	TH1D *hlepPlusPt[4][4];
	TH1D *hlepMinusPt[4][4];
	
	TH1D *hlepPt[4][4];
	TH1D *hrapiditydiffMarco[4][4];
	
	TH1D *hlepPlusCosTheta_gen[4][4];
	TH1D *hlepMinusCosTheta_gen[4][4];
	TH1D *hpseudorapiditydiff_gen[4][4];
	TH1D *hrapiditydiff_gen[4][4];
	TH1D *hrapiditydiffMarco_gen[4][4];

	TH2D *hlepChargeAsym_gen2d[4][4];
	TH2D *hlepAzimAsym_gen2d[4][4];
	TH2D *hlepAzimAsym2_gen2d[4][4];
	TH2D *htopSpinCorr_gen2d[4][4];
	TH2D *htopCosTheta_gen2d[4][4];
	TH2D *hlepCosTheta_gen2d[4][4];
	TH2D *hlepPlusCosTheta_gen2d[4][4];
	TH2D *hlepMinusCosTheta_gen2d[4][4];
	TH2D *hpseudorapiditydiff_gen2d[4][4];
	TH2D *hrapiditydiff_gen2d[4][4];
	TH2D *hrapiditydiffMarco_gen2d[4][4];	
	
	TH2D *hlepChargeAsym_ttpT_gen2d[4][4];
	TH2D *hlepAzimAsym_ttpT_gen2d[4][4];
	TH2D *hlepAzimAsym2_ttpT_gen2d[4][4];
	TH2D *htopSpinCorr_ttpT_gen2d[4][4];
	TH2D *htopCosTheta_ttpT_gen2d[4][4];
	TH2D *hlepCosTheta_ttpT_gen2d[4][4];
	TH2D *hlepPlusCosTheta_ttpT_gen2d[4][4];
	TH2D *hlepMinusCosTheta_ttpT_gen2d[4][4];
	TH2D *hpseudorapiditydiff_ttpT_gen2d[4][4];
	TH2D *hrapiditydiff_ttpT_gen2d[4][4];
	TH2D *hrapiditydiffMarco_ttpT_gen2d[4][4];	


	TH2D *hlepChargeAsym_ttRapidity2_gen2d[4][4];
	TH2D *hlepAzimAsym_ttRapidity2_gen2d[4][4];
	TH2D *hlepAzimAsym2_ttRapidity2_gen2d[4][4];
	TH2D *htopSpinCorr_ttRapidity2_gen2d[4][4];
	TH2D *htopCosTheta_ttRapidity2_gen2d[4][4];
	TH2D *hlepCosTheta_ttRapidity2_gen2d[4][4];
	TH2D *hlepPlusCosTheta_ttRapidity2_gen2d[4][4];
	TH2D *hlepMinusCosTheta_ttRapidity2_gen2d[4][4];
	TH2D *hpseudorapiditydiff_ttRapidity2_gen2d[4][4];
	TH2D *hrapiditydiff_ttRapidity2_gen2d[4][4];
	TH2D *hrapiditydiffMarco_ttRapidity2_gen2d[4][4];
	
	
	TH1D *hlepChargeAsymGenDiff[4][4];
	TH1D *hlepAzimAsymGenDiff[4][4];
	TH1D *hlepAzimAsym2GenDiff[4][4];
	TH1D *htopSpinCorrGenDiff[4][4];
	TH1D *htopCosThetaGenDiff[4][4];
	TH1D *hlepCosThetaGenDiff[4][4];
	TH1D *hlepPlusCosThetaGenDiff[4][4];
	TH1D *hlepMinusCosThetaGenDiff[4][4];
	TH1D *hpseudorapiditydiffGenDiff[4][4];
	TH1D *hrapiditydiffGenDiff[4][4];
	TH1D *hrapiditydiffMarcoGenDiff[4][4];


	TH2D *htopCosTheta_2d[4][4];
	TH2D *hlepCosTheta_2d[4][4];
	TH2D *hlepChargeAsym_2d[4][4];
	TH2D *hlepAzimAsym_2d[4][4];
	TH2D *htopSpinCorr_2d[4][4];
	TH2D *httMass_2d[4][4];
	TH2D *httMass_nojetsmear_2d[4][4];
	TH2D *httpT_2d[4][4];
	TH2D *httpT_nojetsmear_2d[4][4];
	TH2D *htopP_2d[4][4];
	TH2D *htopP_nojetsmear_2d[4][4];
	
	TH2D *hmasslb_2d[4][4];
	TH2D *habcd_2d[4][4];

	TH2D *httmasssm_2d[4][4];
	TH2D *htopmasssm_2d[4][4];

	TH2D *htop1pCMsm_2d[4][4];
	TH2D *htop2pCMsm_2d[4][4];

	TH2D *htop1pTsm_2d[4][4];
	TH2D *htop1etasm_2d[4][4];
	TH2D *htop1phism_2d[4][4];
	TH2D *htop2pTsm_2d[4][4];
	TH2D *htop2etasm_2d[4][4];
	TH2D *htop2phism_2d[4][4];


	TH1D *hdRlbtruegenp[4][4];
	TH1D *hdRlbfalsegenp[4][4];
	TH1D *hdRlbratiogenp[4][4];
	TH1D *htopptgenp[4][4];
	
	TH1D *htopMassgenp[4][4];
	TH1D *htheSumBtagJetPtgenp[4][4];
	TH1D *hthefirstBtagJetPtgenp[4][4];
	TH1D *hthesecondBtagJetPtgenp[4][4];
	TH1D *htheleadinglepPtgenp[4][4];
	TH1D *hthesecondlepPtgenp[4][4];
	TH1D *htheSumLepPtgenp[4][4];
	TH1D *htheleadingNuPtgenp[4][4];
	TH1D *hthesecondNuPtgenp[4][4];
	TH1D *hMETgenp[4][4];
	TH1D *htheSumLBPtgenp[4][4];
	TH1D *htheSumLBNPtgenp[4][4];
	TH1D *httMassGluongenp[4][4];
	TH1D *httMassQuarkgenp[4][4];
	TH1D *httRapidityGluongenp[4][4];
	TH1D *httRapidityQuarkgenp[4][4];
	TH1D *hllbbRapidityGluongenp[4][4];
	TH1D *hllbbRapidityQuarkgenp[4][4];
	

	TH2D *hmasslbgenp_2d[4][4];
	TH2D *hmasslbgenmatch1_2d[4][4];
	TH2D *hmasslbgenmatch_2d[4][4];
	TH2D *htopptdrgenp_2d[4][4];
	
	TH1D *hlepPlusCosThetaTau_gen[4][4];
	TH1D *hlepMinusCosThetaTau_gen[4][4];
	
	TH1D *hlepPlusxTau_gen[4][4];
	TH1D *hlepMinusxTau_gen[4][4];

	// For the DY Estimation
	TH1D* hmetInDYEst[4][4];         // MET in Z window for the DY Estimation
	TH1D* hmetOutDYEst[4][4];        // MET outside Z window for the DY Estimation
	TH1D* hdilMassWithMetDYEst[4][4];// Dilepton mass with MET requirement for DY estimation
	TH1D* hdilMassNoMetDYEst[4][4];  // Dilepton mass without MET requirement for DY estimation

 private:
	Int_t   ngoodlep_;
        Int_t   ngoodel_;
        Int_t   ngoodmu_;



};



#endif
