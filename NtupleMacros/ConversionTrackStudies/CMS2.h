// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include <vector> 

#define PARANOIA

using namespace std; 
class CMS2 {
private: 
protected: 
	unsigned int index;
	TString evt_CMS2tag_;
	TBranch *evt_CMS2tag_branch;
	bool evt_CMS2tag_isLoaded;
	TString evt_dataset_;
	TBranch *evt_dataset_branch;
	bool evt_dataset_isLoaded;
	vector<TString> l1_trigNames_;
	TBranch *l1_trigNames_branch;
	bool l1_trigNames_isLoaded;
	float evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float hcalnoise_eventChargeFraction_;
	TBranch *hcalnoise_eventChargeFraction_branch;
	bool hcalnoise_eventChargeFraction_isLoaded;
	float hcalnoise_eventEMEnergy_;
	TBranch *hcalnoise_eventEMEnergy_branch;
	bool hcalnoise_eventEMEnergy_isLoaded;
	float hcalnoise_eventEMFraction_;
	TBranch *hcalnoise_eventEMFraction_branch;
	bool hcalnoise_eventEMFraction_isLoaded;
	float hcalnoise_eventHadEnergy_;
	TBranch *hcalnoise_eventHadEnergy_branch;
	bool hcalnoise_eventHadEnergy_isLoaded;
	float hcalnoise_eventTrackEnergy_;
	TBranch *hcalnoise_eventTrackEnergy_branch;
	bool hcalnoise_eventTrackEnergy_isLoaded;
	float hcalnoise_max10GeVHitTime_;
	TBranch *hcalnoise_max10GeVHitTime_branch;
	bool hcalnoise_max10GeVHitTime_isLoaded;
	float hcalnoise_max25GeVHitTime_;
	TBranch *hcalnoise_max25GeVHitTime_branch;
	bool hcalnoise_max25GeVHitTime_isLoaded;
	float hcalnoise_min10GeVHitTime_;
	TBranch *hcalnoise_min10GeVHitTime_branch;
	bool hcalnoise_min10GeVHitTime_isLoaded;
	float hcalnoise_min25GeVHitTime_;
	TBranch *hcalnoise_min25GeVHitTime_branch;
	bool hcalnoise_min25GeVHitTime_isLoaded;
	float hcalnoise_minE10TS_;
	TBranch *hcalnoise_minE10TS_branch;
	bool hcalnoise_minE10TS_isLoaded;
	float hcalnoise_minE2Over10TS_;
	TBranch *hcalnoise_minE2Over10TS_branch;
	bool hcalnoise_minE2Over10TS_isLoaded;
	float hcalnoise_minE2TS_;
	TBranch *hcalnoise_minE2TS_branch;
	bool hcalnoise_minE2TS_isLoaded;
	float hcalnoise_minHPDEMF_;
	TBranch *hcalnoise_minHPDEMF_branch;
	bool hcalnoise_minHPDEMF_isLoaded;
	float hcalnoise_minRBXEMF_;
	TBranch *hcalnoise_minRBXEMF_branch;
	bool hcalnoise_minRBXEMF_isLoaded;
	float hcalnoise_rms10GeVHitTime_;
	TBranch *hcalnoise_rms10GeVHitTime_branch;
	bool hcalnoise_rms10GeVHitTime_isLoaded;
	float hcalnoise_rms25GeVHitTime_;
	TBranch *hcalnoise_rms25GeVHitTime_branch;
	bool hcalnoise_rms25GeVHitTime_isLoaded;
	float l1_met_etTot_;
	TBranch *l1_met_etTot_branch;
	bool l1_met_etTot_isLoaded;
	float l1_met_met_;
	TBranch *l1_met_met_branch;
	bool l1_met_met_isLoaded;
	float l1_mht_htTot_;
	TBranch *l1_mht_htTot_branch;
	bool l1_mht_htTot_isLoaded;
	float l1_mht_mht_;
	TBranch *l1_mht_mht_branch;
	bool l1_mht_mht_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  l1_met_p4_;
	TBranch *l1_met_p4_branch;
	bool l1_met_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  l1_mht_p4_;
	TBranch *l1_mht_p4_branch;
	bool l1_mht_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_emiso_p4_;
	TBranch *l1_emiso_p4_branch;
	bool l1_emiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_emnoiso_p4_;
	TBranch *l1_emnoiso_p4_branch;
	bool l1_emnoiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_jetsc_p4_;
	TBranch *l1_jetsc_p4_branch;
	bool l1_jetsc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_jetsf_p4_;
	TBranch *l1_jetsf_p4_branch;
	bool l1_jetsf_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_jetst_p4_;
	TBranch *l1_jetst_p4_branch;
	bool l1_jetst_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > l1_mus_p4_;
	TBranch *l1_mus_p4_branch;
	bool l1_mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_p4_;
	TBranch *scs_p4_branch;
	bool scs_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_pos_p4_;
	TBranch *scs_pos_p4_branch;
	bool scs_pos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_vtx_p4_;
	TBranch *scs_vtx_p4_branch;
	bool scs_vtx_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trks_inner_position_;
	TBranch *trks_inner_position_branch;
	bool trks_inner_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trks_outer_p4_;
	TBranch *trks_outer_p4_branch;
	bool trks_outer_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trks_outer_position_;
	TBranch *trks_outer_position_branch;
	bool trks_outer_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trks_trk_p4_;
	TBranch *trks_trk_p4_branch;
	bool trks_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trks_vertex_p4_;
	TBranch *trks_vertex_p4_branch;
	bool trks_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	vector<float> twrs_ecalTime_;
	TBranch *twrs_ecalTime_branch;
	bool twrs_ecalTime_isLoaded;
	vector<float> twrs_emEnergy_;
	TBranch *twrs_emEnergy_branch;
	bool twrs_emEnergy_isLoaded;
	vector<float> twrs_emEt_;
	TBranch *twrs_emEt_branch;
	bool twrs_emEt_isLoaded;
	vector<float> twrs_emEtcorr_;
	TBranch *twrs_emEtcorr_branch;
	bool twrs_emEtcorr_isLoaded;
	vector<float> twrs_eta_;
	TBranch *twrs_eta_branch;
	bool twrs_eta_isLoaded;
	vector<float> twrs_etcorr_;
	TBranch *twrs_etcorr_branch;
	bool twrs_etcorr_isLoaded;
	vector<float> twrs_hadEnergy_;
	TBranch *twrs_hadEnergy_branch;
	bool twrs_hadEnergy_isLoaded;
	vector<float> twrs_hadEt_;
	TBranch *twrs_hadEt_branch;
	bool twrs_hadEt_isLoaded;
	vector<float> twrs_hadEtcorr_;
	TBranch *twrs_hadEtcorr_branch;
	bool twrs_hadEtcorr_isLoaded;
	vector<float> twrs_hcalTime_;
	TBranch *twrs_hcalTime_branch;
	bool twrs_hcalTime_isLoaded;
	vector<float> twrs_outerEnergy_;
	TBranch *twrs_outerEnergy_branch;
	bool twrs_outerEnergy_isLoaded;
	vector<float> twrs_outerEt_;
	TBranch *twrs_outerEt_branch;
	bool twrs_outerEt_isLoaded;
	vector<float> twrs_outerEtcorr_;
	TBranch *twrs_outerEtcorr_branch;
	bool twrs_outerEtcorr_isLoaded;
	vector<float> twrs_pcorr_;
	TBranch *twrs_pcorr_branch;
	bool twrs_pcorr_isLoaded;
	vector<float> twrs_phi_;
	TBranch *twrs_phi_branch;
	bool twrs_phi_isLoaded;
	vector<float> pseudo_dRClosestTower_;
	TBranch *pseudo_dRClosestTower_branch;
	bool pseudo_dRClosestTower_isLoaded;
	vector<float> pseudo_dRClosestTowerEmEt_;
	TBranch *pseudo_dRClosestTowerEmEt_branch;
	bool pseudo_dRClosestTowerEmEt_isLoaded;
	vector<float> pseudo_ecalEta_;
	TBranch *pseudo_ecalEta_branch;
	bool pseudo_ecalEta_isLoaded;
	vector<float> pseudo_ecalIso03_;
	TBranch *pseudo_ecalIso03_branch;
	bool pseudo_ecalIso03_isLoaded;
	vector<float> pseudo_ecalPhi_;
	TBranch *pseudo_ecalPhi_branch;
	bool pseudo_ecalPhi_isLoaded;
	vector<float> pseudo_eta_;
	TBranch *pseudo_eta_branch;
	bool pseudo_eta_isLoaded;
	vector<float> pseudo_hcalD1Iso03_;
	TBranch *pseudo_hcalD1Iso03_branch;
	bool pseudo_hcalD1Iso03_isLoaded;
	vector<float> pseudo_hcalD2Iso03_;
	TBranch *pseudo_hcalD2Iso03_branch;
	bool pseudo_hcalD2Iso03_isLoaded;
	vector<float> pseudo_phi_;
	TBranch *pseudo_phi_branch;
	bool pseudo_phi_isLoaded;
	vector<float> pseudo_tkIso03_;
	TBranch *pseudo_tkIso03_branch;
	bool pseudo_tkIso03_isLoaded;
	vector<float> pseudo_towerEmEt_;
	TBranch *pseudo_towerEmEt_branch;
	bool pseudo_towerEmEt_isLoaded;
	vector<float> pseudo_towerHadEt_;
	TBranch *pseudo_towerHadEt_branch;
	bool pseudo_towerHadEt_isLoaded;
	vector<float> scs_clustersSize_;
	TBranch *scs_clustersSize_branch;
	bool scs_clustersSize_isLoaded;
	vector<float> scs_crystalsSize_;
	TBranch *scs_crystalsSize_branch;
	bool scs_crystalsSize_isLoaded;
	vector<float> scs_e1x3_;
	TBranch *scs_e1x3_branch;
	bool scs_e1x3_isLoaded;
	vector<float> scs_e1x5_;
	TBranch *scs_e1x5_branch;
	bool scs_e1x5_isLoaded;
	vector<float> scs_e2nd_;
	TBranch *scs_e2nd_branch;
	bool scs_e2nd_isLoaded;
	vector<float> scs_e2x2_;
	TBranch *scs_e2x2_branch;
	bool scs_e2x2_isLoaded;
	vector<float> scs_e2x5Max_;
	TBranch *scs_e2x5Max_branch;
	bool scs_e2x5Max_isLoaded;
	vector<float> scs_e3x1_;
	TBranch *scs_e3x1_branch;
	bool scs_e3x1_isLoaded;
	vector<float> scs_e3x2_;
	TBranch *scs_e3x2_branch;
	bool scs_e3x2_isLoaded;
	vector<float> scs_e3x3_;
	TBranch *scs_e3x3_branch;
	bool scs_e3x3_isLoaded;
	vector<float> scs_e4x4_;
	TBranch *scs_e4x4_branch;
	bool scs_e4x4_isLoaded;
	vector<float> scs_e5x5_;
	TBranch *scs_e5x5_branch;
	bool scs_e5x5_isLoaded;
	vector<float> scs_eMax_;
	TBranch *scs_eMax_branch;
	bool scs_eMax_isLoaded;
	vector<float> scs_eSeed_;
	TBranch *scs_eSeed_branch;
	bool scs_eSeed_isLoaded;
	vector<float> scs_energy_;
	TBranch *scs_energy_branch;
	bool scs_energy_isLoaded;
	vector<float> scs_eta_;
	TBranch *scs_eta_branch;
	bool scs_eta_isLoaded;
	vector<float> scs_hoe_;
	TBranch *scs_hoe_branch;
	bool scs_hoe_isLoaded;
	vector<float> scs_phi_;
	TBranch *scs_phi_branch;
	bool scs_phi_isLoaded;
	vector<float> scs_preshowerEnergy_;
	TBranch *scs_preshowerEnergy_branch;
	bool scs_preshowerEnergy_isLoaded;
	vector<float> scs_rawEnergy_;
	TBranch *scs_rawEnergy_branch;
	bool scs_rawEnergy_isLoaded;
	vector<float> scs_sigmaEtaEta_;
	TBranch *scs_sigmaEtaEta_branch;
	bool scs_sigmaEtaEta_isLoaded;
	vector<float> scs_sigmaEtaPhi_;
	TBranch *scs_sigmaEtaPhi_branch;
	bool scs_sigmaEtaPhi_isLoaded;
	vector<float> scs_sigmaIEtaIEta_;
	TBranch *scs_sigmaIEtaIEta_branch;
	bool scs_sigmaIEtaIEta_isLoaded;
	vector<float> scs_sigmaIEtaIPhi_;
	TBranch *scs_sigmaIEtaIPhi_branch;
	bool scs_sigmaIEtaIPhi_isLoaded;
	vector<float> scs_sigmaIPhiIPhi_;
	TBranch *scs_sigmaIPhiIPhi_branch;
	bool scs_sigmaIPhiIPhi_isLoaded;
	vector<float> scs_sigmaPhiPhi_;
	TBranch *scs_sigmaPhiPhi_branch;
	bool scs_sigmaPhiPhi_isLoaded;
	vector<float> trks_chi2_;
	TBranch *trks_chi2_branch;
	bool trks_chi2_isLoaded;
	vector<float> trks_d0_;
	TBranch *trks_d0_branch;
	bool trks_d0_isLoaded;
	vector<float> trks_d0Err_;
	TBranch *trks_d0Err_branch;
	bool trks_d0Err_isLoaded;
	vector<float> trks_d0corr_;
	TBranch *trks_d0corr_branch;
	bool trks_d0corr_isLoaded;
	vector<float> trks_d0corrPhi_;
	TBranch *trks_d0corrPhi_branch;
	bool trks_d0corrPhi_isLoaded;
	vector<float> trks_etaErr_;
	TBranch *trks_etaErr_branch;
	bool trks_etaErr_isLoaded;
	vector<float> trks_layer1_charge_;
	TBranch *trks_layer1_charge_branch;
	bool trks_layer1_charge_isLoaded;
	vector<float> trks_ndof_;
	TBranch *trks_ndof_branch;
	bool trks_ndof_isLoaded;
	vector<float> trks_phiErr_;
	TBranch *trks_phiErr_branch;
	bool trks_phiErr_isLoaded;
	vector<float> trks_ptErr_;
	TBranch *trks_ptErr_branch;
	bool trks_ptErr_isLoaded;
	vector<float> trks_z0_;
	TBranch *trks_z0_branch;
	bool trks_z0_isLoaded;
	vector<float> trks_z0Err_;
	TBranch *trks_z0Err_branch;
	bool trks_z0Err_isLoaded;
	vector<float> trks_z0corr_;
	TBranch *trks_z0corr_branch;
	bool trks_z0corr_isLoaded;
	vector<float> vtxs_chi2_;
	TBranch *vtxs_chi2_branch;
	bool vtxs_chi2_isLoaded;
	vector<float> vtxs_ndof_;
	TBranch *vtxs_ndof_branch;
	bool vtxs_ndof_isLoaded;
	vector<float> vtxs_xError_;
	TBranch *vtxs_xError_branch;
	bool vtxs_xError_isLoaded;
	vector<float> vtxs_yError_;
	TBranch *vtxs_yError_branch;
	bool vtxs_yError_isLoaded;
	vector<float> vtxs_zError_;
	TBranch *vtxs_zError_branch;
	bool vtxs_zError_isLoaded;
	vector<vector<float> > pseudo_ecalIso03_recHitE_;
	TBranch *pseudo_ecalIso03_recHitE_branch;
	bool pseudo_ecalIso03_recHitE_isLoaded;
	vector<vector<float> > pseudo_ecalIso03_recHitEt_;
	TBranch *pseudo_ecalIso03_recHitEt_branch;
	bool pseudo_ecalIso03_recHitEt_isLoaded;
	vector<vector<float> > pseudo_srDIdx_;
	TBranch *pseudo_srDIdx_branch;
	bool pseudo_srDIdx_isLoaded;
	vector<vector<float> > vtxs_covMatrix_;
	TBranch *vtxs_covMatrix_branch;
	bool vtxs_covMatrix_isLoaded;
	int evt_bunchCrossing_;
	TBranch *evt_bunchCrossing_branch;
	bool evt_bunchCrossing_isLoaded;
	int evt_experimentType_;
	TBranch *evt_experimentType_branch;
	bool evt_experimentType_isLoaded;
	int evt_orbitNumber_;
	TBranch *evt_orbitNumber_branch;
	bool evt_orbitNumber_isLoaded;
	int evt_storeNumber_;
	TBranch *evt_storeNumber_branch;
	bool evt_storeNumber_isLoaded;
	int hcalnoise_maxHPDHits_;
	TBranch *hcalnoise_maxHPDHits_branch;
	bool hcalnoise_maxHPDHits_isLoaded;
	int hcalnoise_maxRBXHits_;
	TBranch *hcalnoise_maxRBXHits_branch;
	bool hcalnoise_maxRBXHits_isLoaded;
	int hcalnoise_maxZeros_;
	TBranch *hcalnoise_maxZeros_branch;
	bool hcalnoise_maxZeros_isLoaded;
	int hcalnoise_noiseFilterStatus_;
	TBranch *hcalnoise_noiseFilterStatus_branch;
	bool hcalnoise_noiseFilterStatus_isLoaded;
	int hcalnoise_noiseType_;
	TBranch *hcalnoise_noiseType_branch;
	bool hcalnoise_noiseType_isLoaded;
	int hcalnoise_num10GeVHits_;
	TBranch *hcalnoise_num10GeVHits_branch;
	bool hcalnoise_num10GeVHits_isLoaded;
	int hcalnoise_num25GeVHits_;
	TBranch *hcalnoise_num25GeVHits_branch;
	bool hcalnoise_num25GeVHits_isLoaded;
	int hcalnoise_numProblematicRBXs_;
	TBranch *hcalnoise_numProblematicRBXs_branch;
	bool hcalnoise_numProblematicRBXs_isLoaded;
	int hcalnoise_passHighLevelNoiseFilter_;
	TBranch *hcalnoise_passHighLevelNoiseFilter_branch;
	bool hcalnoise_passHighLevelNoiseFilter_isLoaded;
	int hcalnoise_passLooseNoiseFilter_;
	TBranch *hcalnoise_passLooseNoiseFilter_branch;
	bool hcalnoise_passLooseNoiseFilter_isLoaded;
	int hcalnoise_passTightNoiseFilter_;
	TBranch *hcalnoise_passTightNoiseFilter_branch;
	bool hcalnoise_passTightNoiseFilter_isLoaded;
	int l1_nemiso_;
	TBranch *l1_nemiso_branch;
	bool l1_nemiso_isLoaded;
	int l1_nemnoiso_;
	TBranch *l1_nemnoiso_branch;
	bool l1_nemnoiso_isLoaded;
	int l1_njetsc_;
	TBranch *l1_njetsc_branch;
	bool l1_njetsc_isLoaded;
	int l1_njetsf_;
	TBranch *l1_njetsf_branch;
	bool l1_njetsf_isLoaded;
	int l1_njetst_;
	TBranch *l1_njetst_branch;
	bool l1_njetst_isLoaded;
	int l1_nmus_;
	TBranch *l1_nmus_branch;
	bool l1_nmus_isLoaded;
	vector<int> l1_emiso_ieta_;
	TBranch *l1_emiso_ieta_branch;
	bool l1_emiso_ieta_isLoaded;
	vector<int> l1_emiso_iphi_;
	TBranch *l1_emiso_iphi_branch;
	bool l1_emiso_iphi_isLoaded;
	vector<int> l1_emiso_rawId_;
	TBranch *l1_emiso_rawId_branch;
	bool l1_emiso_rawId_isLoaded;
	vector<int> l1_emiso_type_;
	TBranch *l1_emiso_type_branch;
	bool l1_emiso_type_isLoaded;
	vector<int> l1_emnoiso_ieta_;
	TBranch *l1_emnoiso_ieta_branch;
	bool l1_emnoiso_ieta_isLoaded;
	vector<int> l1_emnoiso_iphi_;
	TBranch *l1_emnoiso_iphi_branch;
	bool l1_emnoiso_iphi_isLoaded;
	vector<int> l1_emnoiso_rawId_;
	TBranch *l1_emnoiso_rawId_branch;
	bool l1_emnoiso_rawId_isLoaded;
	vector<int> l1_emnoiso_type_;
	TBranch *l1_emnoiso_type_branch;
	bool l1_emnoiso_type_isLoaded;
	vector<int> l1_jetsc_ieta_;
	TBranch *l1_jetsc_ieta_branch;
	bool l1_jetsc_ieta_isLoaded;
	vector<int> l1_jetsc_iphi_;
	TBranch *l1_jetsc_iphi_branch;
	bool l1_jetsc_iphi_isLoaded;
	vector<int> l1_jetsc_rawId_;
	TBranch *l1_jetsc_rawId_branch;
	bool l1_jetsc_rawId_isLoaded;
	vector<int> l1_jetsc_type_;
	TBranch *l1_jetsc_type_branch;
	bool l1_jetsc_type_isLoaded;
	vector<int> l1_jetsf_ieta_;
	TBranch *l1_jetsf_ieta_branch;
	bool l1_jetsf_ieta_isLoaded;
	vector<int> l1_jetsf_iphi_;
	TBranch *l1_jetsf_iphi_branch;
	bool l1_jetsf_iphi_isLoaded;
	vector<int> l1_jetsf_rawId_;
	TBranch *l1_jetsf_rawId_branch;
	bool l1_jetsf_rawId_isLoaded;
	vector<int> l1_jetsf_type_;
	TBranch *l1_jetsf_type_branch;
	bool l1_jetsf_type_isLoaded;
	vector<int> l1_jetst_ieta_;
	TBranch *l1_jetst_ieta_branch;
	bool l1_jetst_ieta_isLoaded;
	vector<int> l1_jetst_iphi_;
	TBranch *l1_jetst_iphi_branch;
	bool l1_jetst_iphi_isLoaded;
	vector<int> l1_jetst_rawId_;
	TBranch *l1_jetst_rawId_branch;
	bool l1_jetst_rawId_isLoaded;
	vector<int> l1_jetst_type_;
	TBranch *l1_jetst_type_branch;
	bool l1_jetst_type_isLoaded;
	vector<int> l1_mus_flags_;
	TBranch *l1_mus_flags_branch;
	bool l1_mus_flags_isLoaded;
	vector<int> l1_mus_q_;
	TBranch *l1_mus_q_branch;
	bool l1_mus_q_isLoaded;
	vector<int> l1_mus_qual_;
	TBranch *l1_mus_qual_branch;
	bool l1_mus_qual_isLoaded;
	vector<int> l1_mus_qualFlags_;
	TBranch *l1_mus_qualFlags_branch;
	bool l1_mus_qualFlags_isLoaded;
	vector<int> pxl_ndigis_pxb_;
	TBranch *pxl_ndigis_pxb_branch;
	bool pxl_ndigis_pxb_isLoaded;
	vector<int> pxl_ndigis_pxf_;
	TBranch *pxl_ndigis_pxf_branch;
	bool pxl_ndigis_pxf_isLoaded;
	vector<int> scs_elsidx_;
	TBranch *scs_elsidx_branch;
	bool scs_elsidx_isLoaded;
	vector<int> scs_severitySeed_;
	TBranch *scs_severitySeed_branch;
	bool scs_severitySeed_isLoaded;
	vector<int> trks_algo_;
	TBranch *trks_algo_branch;
	bool trks_algo_isLoaded;
	vector<int> trks_charge_;
	TBranch *trks_charge_branch;
	bool trks_charge_isLoaded;
	vector<int> trks_exp_innerlayers_;
	TBranch *trks_exp_innerlayers_branch;
	bool trks_exp_innerlayers_isLoaded;
	vector<int> trks_exp_outerlayers_;
	TBranch *trks_exp_outerlayers_branch;
	bool trks_exp_outerlayers_isLoaded;
	vector<int> trks_layer1_det_;
	TBranch *trks_layer1_det_branch;
	bool trks_layer1_det_isLoaded;
	vector<int> trks_layer1_layer_;
	TBranch *trks_layer1_layer_branch;
	bool trks_layer1_layer_isLoaded;
	vector<int> trks_layer1_sizerphi_;
	TBranch *trks_layer1_sizerphi_branch;
	bool trks_layer1_sizerphi_isLoaded;
	vector<int> trks_layer1_sizerz_;
	TBranch *trks_layer1_sizerz_branch;
	bool trks_layer1_sizerz_isLoaded;
	vector<int> trks_lostHits_;
	TBranch *trks_lostHits_branch;
	bool trks_lostHits_isLoaded;
	vector<int> trks_lost_pixelhits_;
	TBranch *trks_lost_pixelhits_branch;
	bool trks_lost_pixelhits_isLoaded;
	vector<int> trks_nlayers_;
	TBranch *trks_nlayers_branch;
	bool trks_nlayers_isLoaded;
	vector<int> trks_nlayers3D_;
	TBranch *trks_nlayers3D_branch;
	bool trks_nlayers3D_isLoaded;
	vector<int> trks_nlayersLost_;
	TBranch *trks_nlayersLost_branch;
	bool trks_nlayersLost_isLoaded;
	vector<int> trks_qualityMask_;
	TBranch *trks_qualityMask_branch;
	bool trks_qualityMask_isLoaded;
	vector<int> trks_validHits_;
	TBranch *trks_validHits_branch;
	bool trks_validHits_isLoaded;
	vector<int> trks_valid_pixelhits_;
	TBranch *trks_valid_pixelhits_branch;
	bool trks_valid_pixelhits_isLoaded;
	vector<int> vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int> vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<int> vtxs_tracksSize_;
	TBranch *vtxs_tracksSize_branch;
	bool vtxs_tracksSize_isLoaded;
	vector<vector<int> > pseudo_srFlags_;
	TBranch *pseudo_srFlags_branch;
	bool pseudo_srFlags_isLoaded;
	unsigned int evt_ntwrs_;
	TBranch *evt_ntwrs_branch;
	bool evt_ntwrs_isLoaded;
	unsigned int evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int l1_bits1_;
	TBranch *l1_bits1_branch;
	bool l1_bits1_isLoaded;
	unsigned int l1_bits2_;
	TBranch *l1_bits2_branch;
	bool l1_bits2_isLoaded;
	unsigned int l1_bits3_;
	TBranch *l1_bits3_branch;
	bool l1_bits3_isLoaded;
	unsigned int l1_bits4_;
	TBranch *l1_bits4_branch;
	bool l1_bits4_isLoaded;
	unsigned int l1_techbits1_;
	TBranch *l1_techbits1_branch;
	bool l1_techbits1_isLoaded;
	unsigned int l1_techbits2_;
	TBranch *l1_techbits2_branch;
	bool l1_techbits2_isLoaded;
	unsigned int evt_nscs_;
	TBranch *evt_nscs_branch;
	bool evt_nscs_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	vector<unsigned int> twrs_numBadEcalCells_;
	TBranch *twrs_numBadEcalCells_branch;
	bool twrs_numBadEcalCells_isLoaded;
	vector<unsigned int> twrs_numBadHcalCells_;
	TBranch *twrs_numBadHcalCells_branch;
	bool twrs_numBadHcalCells_isLoaded;
	vector<unsigned int> twrs_numProblematicEcalCells_;
	TBranch *twrs_numProblematicEcalCells_branch;
	bool twrs_numProblematicEcalCells_isLoaded;
	vector<unsigned int> twrs_numProblematicHcalCells_;
	TBranch *twrs_numProblematicHcalCells_branch;
	bool twrs_numProblematicHcalCells_isLoaded;
	vector<unsigned int> twrs_numRecoveredEcalCells_;
	TBranch *twrs_numRecoveredEcalCells_branch;
	bool twrs_numRecoveredEcalCells_isLoaded;
	vector<unsigned int> twrs_numRecoveredHcalCells_;
	TBranch *twrs_numRecoveredHcalCells_branch;
	bool twrs_numRecoveredHcalCells_isLoaded;
public: 
int ScanChain(class TChain* chain, int nEvents=-1, std::string skimFilePrefix="");
void Init(TTree *tree) {
	l1_met_p4_branch = 0;
	if (tree->GetAlias("l1_met_p4") != 0) {
		l1_met_p4_branch = tree->GetBranch(tree->GetAlias("l1_met_p4"));
		l1_met_p4_branch->SetAddress(&l1_met_p4_);
	}
	if(l1_met_p4_branch == 0 ) {
	cout << "Branch l1_met_p4 does not exist." << endl;
	}
	l1_mht_p4_branch = 0;
	if (tree->GetAlias("l1_mht_p4") != 0) {
		l1_mht_p4_branch = tree->GetBranch(tree->GetAlias("l1_mht_p4"));
		l1_mht_p4_branch->SetAddress(&l1_mht_p4_);
	}
	if(l1_mht_p4_branch == 0 ) {
	cout << "Branch l1_mht_p4 does not exist." << endl;
	}
	l1_emiso_p4_branch = 0;
	if (tree->GetAlias("l1_emiso_p4") != 0) {
		l1_emiso_p4_branch = tree->GetBranch(tree->GetAlias("l1_emiso_p4"));
		l1_emiso_p4_branch->SetAddress(&l1_emiso_p4_);
	}
	if(l1_emiso_p4_branch == 0 ) {
	cout << "Branch l1_emiso_p4 does not exist." << endl;
	}
	l1_emnoiso_p4_branch = 0;
	if (tree->GetAlias("l1_emnoiso_p4") != 0) {
		l1_emnoiso_p4_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_p4"));
		l1_emnoiso_p4_branch->SetAddress(&l1_emnoiso_p4_);
	}
	if(l1_emnoiso_p4_branch == 0 ) {
	cout << "Branch l1_emnoiso_p4 does not exist." << endl;
	}
	l1_jetsc_p4_branch = 0;
	if (tree->GetAlias("l1_jetsc_p4") != 0) {
		l1_jetsc_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_p4"));
		l1_jetsc_p4_branch->SetAddress(&l1_jetsc_p4_);
	}
	if(l1_jetsc_p4_branch == 0 ) {
	cout << "Branch l1_jetsc_p4 does not exist." << endl;
	}
	l1_jetsf_p4_branch = 0;
	if (tree->GetAlias("l1_jetsf_p4") != 0) {
		l1_jetsf_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_p4"));
		l1_jetsf_p4_branch->SetAddress(&l1_jetsf_p4_);
	}
	if(l1_jetsf_p4_branch == 0 ) {
	cout << "Branch l1_jetsf_p4 does not exist." << endl;
	}
	l1_jetst_p4_branch = 0;
	if (tree->GetAlias("l1_jetst_p4") != 0) {
		l1_jetst_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetst_p4"));
		l1_jetst_p4_branch->SetAddress(&l1_jetst_p4_);
	}
	if(l1_jetst_p4_branch == 0 ) {
	cout << "Branch l1_jetst_p4 does not exist." << endl;
	}
	l1_mus_p4_branch = 0;
	if (tree->GetAlias("l1_mus_p4") != 0) {
		l1_mus_p4_branch = tree->GetBranch(tree->GetAlias("l1_mus_p4"));
		l1_mus_p4_branch->SetAddress(&l1_mus_p4_);
	}
	if(l1_mus_p4_branch == 0 ) {
	cout << "Branch l1_mus_p4 does not exist." << endl;
	}
	scs_p4_branch = 0;
	if (tree->GetAlias("scs_p4") != 0) {
		scs_p4_branch = tree->GetBranch(tree->GetAlias("scs_p4"));
		scs_p4_branch->SetAddress(&scs_p4_);
	}
	if(scs_p4_branch == 0 ) {
	cout << "Branch scs_p4 does not exist." << endl;
	}
	scs_pos_p4_branch = 0;
	if (tree->GetAlias("scs_pos_p4") != 0) {
		scs_pos_p4_branch = tree->GetBranch(tree->GetAlias("scs_pos_p4"));
		scs_pos_p4_branch->SetAddress(&scs_pos_p4_);
	}
	if(scs_pos_p4_branch == 0 ) {
	cout << "Branch scs_pos_p4 does not exist." << endl;
	}
	scs_vtx_p4_branch = 0;
	if (tree->GetAlias("scs_vtx_p4") != 0) {
		scs_vtx_p4_branch = tree->GetBranch(tree->GetAlias("scs_vtx_p4"));
		scs_vtx_p4_branch->SetAddress(&scs_vtx_p4_);
	}
	if(scs_vtx_p4_branch == 0 ) {
	cout << "Branch scs_vtx_p4 does not exist." << endl;
	}
	trks_inner_position_branch = 0;
	if (tree->GetAlias("trks_inner_position") != 0) {
		trks_inner_position_branch = tree->GetBranch(tree->GetAlias("trks_inner_position"));
		trks_inner_position_branch->SetAddress(&trks_inner_position_);
	}
	if(trks_inner_position_branch == 0 ) {
	cout << "Branch trks_inner_position does not exist." << endl;
	}
	trks_outer_p4_branch = 0;
	if (tree->GetAlias("trks_outer_p4") != 0) {
		trks_outer_p4_branch = tree->GetBranch(tree->GetAlias("trks_outer_p4"));
		trks_outer_p4_branch->SetAddress(&trks_outer_p4_);
	}
	if(trks_outer_p4_branch == 0 ) {
	cout << "Branch trks_outer_p4 does not exist." << endl;
	}
	trks_outer_position_branch = 0;
	if (tree->GetAlias("trks_outer_position") != 0) {
		trks_outer_position_branch = tree->GetBranch(tree->GetAlias("trks_outer_position"));
		trks_outer_position_branch->SetAddress(&trks_outer_position_);
	}
	if(trks_outer_position_branch == 0 ) {
	cout << "Branch trks_outer_position does not exist." << endl;
	}
	trks_trk_p4_branch = 0;
	if (tree->GetAlias("trks_trk_p4") != 0) {
		trks_trk_p4_branch = tree->GetBranch(tree->GetAlias("trks_trk_p4"));
		trks_trk_p4_branch->SetAddress(&trks_trk_p4_);
	}
	if(trks_trk_p4_branch == 0 ) {
	cout << "Branch trks_trk_p4 does not exist." << endl;
	}
	trks_vertex_p4_branch = 0;
	if (tree->GetAlias("trks_vertex_p4") != 0) {
		trks_vertex_p4_branch = tree->GetBranch(tree->GetAlias("trks_vertex_p4"));
		trks_vertex_p4_branch->SetAddress(&trks_vertex_p4_);
	}
	if(trks_vertex_p4_branch == 0 ) {
	cout << "Branch trks_vertex_p4 does not exist." << endl;
	}
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		vtxs_position_branch->SetAddress(&vtxs_position_);
	}
	if(vtxs_position_branch == 0 ) {
	cout << "Branch vtxs_position does not exist." << endl;
	}
  tree->SetMakeClass(1);
	evt_CMS2tag_branch = 0;
	if (tree->GetAlias("evt_CMS2tag") != 0) {
		evt_CMS2tag_branch = tree->GetBranch(tree->GetAlias("evt_CMS2tag"));
		evt_CMS2tag_branch->SetAddress(&evt_CMS2tag_);
	}
	if(evt_CMS2tag_branch == 0 ) {
	cout << "Branch evt_CMS2tag does not exist." << endl;
	}
	evt_dataset_branch = 0;
	if (tree->GetAlias("evt_dataset") != 0) {
		evt_dataset_branch = tree->GetBranch(tree->GetAlias("evt_dataset"));
		evt_dataset_branch->SetAddress(&evt_dataset_);
	}
	if(evt_dataset_branch == 0 ) {
	cout << "Branch evt_dataset does not exist." << endl;
	}
	l1_trigNames_branch = 0;
	if (tree->GetAlias("l1_trigNames") != 0) {
		l1_trigNames_branch = tree->GetBranch(tree->GetAlias("l1_trigNames"));
		l1_trigNames_branch->SetAddress(&l1_trigNames_);
	}
	if(l1_trigNames_branch == 0 ) {
	cout << "Branch l1_trigNames does not exist." << endl;
	}
	evt_bField_branch = 0;
	if (tree->GetAlias("evt_bField") != 0) {
		evt_bField_branch = tree->GetBranch(tree->GetAlias("evt_bField"));
		evt_bField_branch->SetAddress(&evt_bField_);
	}
	if(evt_bField_branch == 0 ) {
	cout << "Branch evt_bField does not exist." << endl;
	}
	hcalnoise_eventChargeFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventChargeFraction") != 0) {
		hcalnoise_eventChargeFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventChargeFraction"));
		hcalnoise_eventChargeFraction_branch->SetAddress(&hcalnoise_eventChargeFraction_);
	}
	if(hcalnoise_eventChargeFraction_branch == 0 ) {
	cout << "Branch hcalnoise_eventChargeFraction does not exist." << endl;
	}
	hcalnoise_eventEMEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMEnergy") != 0) {
		hcalnoise_eventEMEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMEnergy"));
		hcalnoise_eventEMEnergy_branch->SetAddress(&hcalnoise_eventEMEnergy_);
	}
	if(hcalnoise_eventEMEnergy_branch == 0 ) {
	cout << "Branch hcalnoise_eventEMEnergy does not exist." << endl;
	}
	hcalnoise_eventEMFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMFraction") != 0) {
		hcalnoise_eventEMFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMFraction"));
		hcalnoise_eventEMFraction_branch->SetAddress(&hcalnoise_eventEMFraction_);
	}
	if(hcalnoise_eventEMFraction_branch == 0 ) {
	cout << "Branch hcalnoise_eventEMFraction does not exist." << endl;
	}
	hcalnoise_eventHadEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventHadEnergy") != 0) {
		hcalnoise_eventHadEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventHadEnergy"));
		hcalnoise_eventHadEnergy_branch->SetAddress(&hcalnoise_eventHadEnergy_);
	}
	if(hcalnoise_eventHadEnergy_branch == 0 ) {
	cout << "Branch hcalnoise_eventHadEnergy does not exist." << endl;
	}
	hcalnoise_eventTrackEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventTrackEnergy") != 0) {
		hcalnoise_eventTrackEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventTrackEnergy"));
		hcalnoise_eventTrackEnergy_branch->SetAddress(&hcalnoise_eventTrackEnergy_);
	}
	if(hcalnoise_eventTrackEnergy_branch == 0 ) {
	cout << "Branch hcalnoise_eventTrackEnergy does not exist." << endl;
	}
	hcalnoise_max10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max10GeVHitTime") != 0) {
		hcalnoise_max10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max10GeVHitTime"));
		hcalnoise_max10GeVHitTime_branch->SetAddress(&hcalnoise_max10GeVHitTime_);
	}
	if(hcalnoise_max10GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_max10GeVHitTime does not exist." << endl;
	}
	hcalnoise_max25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max25GeVHitTime") != 0) {
		hcalnoise_max25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max25GeVHitTime"));
		hcalnoise_max25GeVHitTime_branch->SetAddress(&hcalnoise_max25GeVHitTime_);
	}
	if(hcalnoise_max25GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_max25GeVHitTime does not exist." << endl;
	}
	hcalnoise_min10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min10GeVHitTime") != 0) {
		hcalnoise_min10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min10GeVHitTime"));
		hcalnoise_min10GeVHitTime_branch->SetAddress(&hcalnoise_min10GeVHitTime_);
	}
	if(hcalnoise_min10GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_min10GeVHitTime does not exist." << endl;
	}
	hcalnoise_min25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min25GeVHitTime") != 0) {
		hcalnoise_min25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min25GeVHitTime"));
		hcalnoise_min25GeVHitTime_branch->SetAddress(&hcalnoise_min25GeVHitTime_);
	}
	if(hcalnoise_min25GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_min25GeVHitTime does not exist." << endl;
	}
	hcalnoise_minE10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE10TS") != 0) {
		hcalnoise_minE10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE10TS"));
		hcalnoise_minE10TS_branch->SetAddress(&hcalnoise_minE10TS_);
	}
	if(hcalnoise_minE10TS_branch == 0 ) {
	cout << "Branch hcalnoise_minE10TS does not exist." << endl;
	}
	hcalnoise_minE2Over10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2Over10TS") != 0) {
		hcalnoise_minE2Over10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2Over10TS"));
		hcalnoise_minE2Over10TS_branch->SetAddress(&hcalnoise_minE2Over10TS_);
	}
	if(hcalnoise_minE2Over10TS_branch == 0 ) {
	cout << "Branch hcalnoise_minE2Over10TS does not exist." << endl;
	}
	hcalnoise_minE2TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2TS") != 0) {
		hcalnoise_minE2TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2TS"));
		hcalnoise_minE2TS_branch->SetAddress(&hcalnoise_minE2TS_);
	}
	if(hcalnoise_minE2TS_branch == 0 ) {
	cout << "Branch hcalnoise_minE2TS does not exist." << endl;
	}
	hcalnoise_minHPDEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minHPDEMF") != 0) {
		hcalnoise_minHPDEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minHPDEMF"));
		hcalnoise_minHPDEMF_branch->SetAddress(&hcalnoise_minHPDEMF_);
	}
	if(hcalnoise_minHPDEMF_branch == 0 ) {
	cout << "Branch hcalnoise_minHPDEMF does not exist." << endl;
	}
	hcalnoise_minRBXEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minRBXEMF") != 0) {
		hcalnoise_minRBXEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minRBXEMF"));
		hcalnoise_minRBXEMF_branch->SetAddress(&hcalnoise_minRBXEMF_);
	}
	if(hcalnoise_minRBXEMF_branch == 0 ) {
	cout << "Branch hcalnoise_minRBXEMF does not exist." << endl;
	}
	hcalnoise_rms10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms10GeVHitTime") != 0) {
		hcalnoise_rms10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms10GeVHitTime"));
		hcalnoise_rms10GeVHitTime_branch->SetAddress(&hcalnoise_rms10GeVHitTime_);
	}
	if(hcalnoise_rms10GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_rms10GeVHitTime does not exist." << endl;
	}
	hcalnoise_rms25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms25GeVHitTime") != 0) {
		hcalnoise_rms25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms25GeVHitTime"));
		hcalnoise_rms25GeVHitTime_branch->SetAddress(&hcalnoise_rms25GeVHitTime_);
	}
	if(hcalnoise_rms25GeVHitTime_branch == 0 ) {
	cout << "Branch hcalnoise_rms25GeVHitTime does not exist." << endl;
	}
	l1_met_etTot_branch = 0;
	if (tree->GetAlias("l1_met_etTot") != 0) {
		l1_met_etTot_branch = tree->GetBranch(tree->GetAlias("l1_met_etTot"));
		l1_met_etTot_branch->SetAddress(&l1_met_etTot_);
	}
	if(l1_met_etTot_branch == 0 ) {
	cout << "Branch l1_met_etTot does not exist." << endl;
	}
	l1_met_met_branch = 0;
	if (tree->GetAlias("l1_met_met") != 0) {
		l1_met_met_branch = tree->GetBranch(tree->GetAlias("l1_met_met"));
		l1_met_met_branch->SetAddress(&l1_met_met_);
	}
	if(l1_met_met_branch == 0 ) {
	cout << "Branch l1_met_met does not exist." << endl;
	}
	l1_mht_htTot_branch = 0;
	if (tree->GetAlias("l1_mht_htTot") != 0) {
		l1_mht_htTot_branch = tree->GetBranch(tree->GetAlias("l1_mht_htTot"));
		l1_mht_htTot_branch->SetAddress(&l1_mht_htTot_);
	}
	if(l1_mht_htTot_branch == 0 ) {
	cout << "Branch l1_mht_htTot does not exist." << endl;
	}
	l1_mht_mht_branch = 0;
	if (tree->GetAlias("l1_mht_mht") != 0) {
		l1_mht_mht_branch = tree->GetBranch(tree->GetAlias("l1_mht_mht"));
		l1_mht_mht_branch->SetAddress(&l1_mht_mht_);
	}
	if(l1_mht_mht_branch == 0 ) {
	cout << "Branch l1_mht_mht does not exist." << endl;
	}
	twrs_ecalTime_branch = 0;
	if (tree->GetAlias("twrs_ecalTime") != 0) {
		twrs_ecalTime_branch = tree->GetBranch(tree->GetAlias("twrs_ecalTime"));
		twrs_ecalTime_branch->SetAddress(&twrs_ecalTime_);
	}
	if(twrs_ecalTime_branch == 0 ) {
	cout << "Branch twrs_ecalTime does not exist." << endl;
	}
	twrs_emEnergy_branch = 0;
	if (tree->GetAlias("twrs_emEnergy") != 0) {
		twrs_emEnergy_branch = tree->GetBranch(tree->GetAlias("twrs_emEnergy"));
		twrs_emEnergy_branch->SetAddress(&twrs_emEnergy_);
	}
	if(twrs_emEnergy_branch == 0 ) {
	cout << "Branch twrs_emEnergy does not exist." << endl;
	}
	twrs_emEt_branch = 0;
	if (tree->GetAlias("twrs_emEt") != 0) {
		twrs_emEt_branch = tree->GetBranch(tree->GetAlias("twrs_emEt"));
		twrs_emEt_branch->SetAddress(&twrs_emEt_);
	}
	if(twrs_emEt_branch == 0 ) {
	cout << "Branch twrs_emEt does not exist." << endl;
	}
	twrs_emEtcorr_branch = 0;
	if (tree->GetAlias("twrs_emEtcorr") != 0) {
		twrs_emEtcorr_branch = tree->GetBranch(tree->GetAlias("twrs_emEtcorr"));
		twrs_emEtcorr_branch->SetAddress(&twrs_emEtcorr_);
	}
	if(twrs_emEtcorr_branch == 0 ) {
	cout << "Branch twrs_emEtcorr does not exist." << endl;
	}
	twrs_eta_branch = 0;
	if (tree->GetAlias("twrs_eta") != 0) {
		twrs_eta_branch = tree->GetBranch(tree->GetAlias("twrs_eta"));
		twrs_eta_branch->SetAddress(&twrs_eta_);
	}
	if(twrs_eta_branch == 0 ) {
	cout << "Branch twrs_eta does not exist." << endl;
	}
	twrs_etcorr_branch = 0;
	if (tree->GetAlias("twrs_etcorr") != 0) {
		twrs_etcorr_branch = tree->GetBranch(tree->GetAlias("twrs_etcorr"));
		twrs_etcorr_branch->SetAddress(&twrs_etcorr_);
	}
	if(twrs_etcorr_branch == 0 ) {
	cout << "Branch twrs_etcorr does not exist." << endl;
	}
	twrs_hadEnergy_branch = 0;
	if (tree->GetAlias("twrs_hadEnergy") != 0) {
		twrs_hadEnergy_branch = tree->GetBranch(tree->GetAlias("twrs_hadEnergy"));
		twrs_hadEnergy_branch->SetAddress(&twrs_hadEnergy_);
	}
	if(twrs_hadEnergy_branch == 0 ) {
	cout << "Branch twrs_hadEnergy does not exist." << endl;
	}
	twrs_hadEt_branch = 0;
	if (tree->GetAlias("twrs_hadEt") != 0) {
		twrs_hadEt_branch = tree->GetBranch(tree->GetAlias("twrs_hadEt"));
		twrs_hadEt_branch->SetAddress(&twrs_hadEt_);
	}
	if(twrs_hadEt_branch == 0 ) {
	cout << "Branch twrs_hadEt does not exist." << endl;
	}
	twrs_hadEtcorr_branch = 0;
	if (tree->GetAlias("twrs_hadEtcorr") != 0) {
		twrs_hadEtcorr_branch = tree->GetBranch(tree->GetAlias("twrs_hadEtcorr"));
		twrs_hadEtcorr_branch->SetAddress(&twrs_hadEtcorr_);
	}
	if(twrs_hadEtcorr_branch == 0 ) {
	cout << "Branch twrs_hadEtcorr does not exist." << endl;
	}
	twrs_hcalTime_branch = 0;
	if (tree->GetAlias("twrs_hcalTime") != 0) {
		twrs_hcalTime_branch = tree->GetBranch(tree->GetAlias("twrs_hcalTime"));
		twrs_hcalTime_branch->SetAddress(&twrs_hcalTime_);
	}
	if(twrs_hcalTime_branch == 0 ) {
	cout << "Branch twrs_hcalTime does not exist." << endl;
	}
	twrs_outerEnergy_branch = 0;
	if (tree->GetAlias("twrs_outerEnergy") != 0) {
		twrs_outerEnergy_branch = tree->GetBranch(tree->GetAlias("twrs_outerEnergy"));
		twrs_outerEnergy_branch->SetAddress(&twrs_outerEnergy_);
	}
	if(twrs_outerEnergy_branch == 0 ) {
	cout << "Branch twrs_outerEnergy does not exist." << endl;
	}
	twrs_outerEt_branch = 0;
	if (tree->GetAlias("twrs_outerEt") != 0) {
		twrs_outerEt_branch = tree->GetBranch(tree->GetAlias("twrs_outerEt"));
		twrs_outerEt_branch->SetAddress(&twrs_outerEt_);
	}
	if(twrs_outerEt_branch == 0 ) {
	cout << "Branch twrs_outerEt does not exist." << endl;
	}
	twrs_outerEtcorr_branch = 0;
	if (tree->GetAlias("twrs_outerEtcorr") != 0) {
		twrs_outerEtcorr_branch = tree->GetBranch(tree->GetAlias("twrs_outerEtcorr"));
		twrs_outerEtcorr_branch->SetAddress(&twrs_outerEtcorr_);
	}
	if(twrs_outerEtcorr_branch == 0 ) {
	cout << "Branch twrs_outerEtcorr does not exist." << endl;
	}
	twrs_pcorr_branch = 0;
	if (tree->GetAlias("twrs_pcorr") != 0) {
		twrs_pcorr_branch = tree->GetBranch(tree->GetAlias("twrs_pcorr"));
		twrs_pcorr_branch->SetAddress(&twrs_pcorr_);
	}
	if(twrs_pcorr_branch == 0 ) {
	cout << "Branch twrs_pcorr does not exist." << endl;
	}
	twrs_phi_branch = 0;
	if (tree->GetAlias("twrs_phi") != 0) {
		twrs_phi_branch = tree->GetBranch(tree->GetAlias("twrs_phi"));
		twrs_phi_branch->SetAddress(&twrs_phi_);
	}
	if(twrs_phi_branch == 0 ) {
	cout << "Branch twrs_phi does not exist." << endl;
	}
	pseudo_dRClosestTower_branch = 0;
	if (tree->GetAlias("pseudo_dRClosestTower") != 0) {
		pseudo_dRClosestTower_branch = tree->GetBranch(tree->GetAlias("pseudo_dRClosestTower"));
		pseudo_dRClosestTower_branch->SetAddress(&pseudo_dRClosestTower_);
	}
	if(pseudo_dRClosestTower_branch == 0 ) {
	cout << "Branch pseudo_dRClosestTower does not exist." << endl;
	}
	pseudo_dRClosestTowerEmEt_branch = 0;
	if (tree->GetAlias("pseudo_dRClosestTowerEmEt") != 0) {
		pseudo_dRClosestTowerEmEt_branch = tree->GetBranch(tree->GetAlias("pseudo_dRClosestTowerEmEt"));
		pseudo_dRClosestTowerEmEt_branch->SetAddress(&pseudo_dRClosestTowerEmEt_);
	}
	if(pseudo_dRClosestTowerEmEt_branch == 0 ) {
	cout << "Branch pseudo_dRClosestTowerEmEt does not exist." << endl;
	}
	pseudo_ecalEta_branch = 0;
	if (tree->GetAlias("pseudo_ecalEta") != 0) {
		pseudo_ecalEta_branch = tree->GetBranch(tree->GetAlias("pseudo_ecalEta"));
		pseudo_ecalEta_branch->SetAddress(&pseudo_ecalEta_);
	}
	if(pseudo_ecalEta_branch == 0 ) {
	cout << "Branch pseudo_ecalEta does not exist." << endl;
	}
	pseudo_ecalIso03_branch = 0;
	if (tree->GetAlias("pseudo_ecalIso03") != 0) {
		pseudo_ecalIso03_branch = tree->GetBranch(tree->GetAlias("pseudo_ecalIso03"));
		pseudo_ecalIso03_branch->SetAddress(&pseudo_ecalIso03_);
	}
	if(pseudo_ecalIso03_branch == 0 ) {
	cout << "Branch pseudo_ecalIso03 does not exist." << endl;
	}
	pseudo_ecalPhi_branch = 0;
	if (tree->GetAlias("pseudo_ecalPhi") != 0) {
		pseudo_ecalPhi_branch = tree->GetBranch(tree->GetAlias("pseudo_ecalPhi"));
		pseudo_ecalPhi_branch->SetAddress(&pseudo_ecalPhi_);
	}
	if(pseudo_ecalPhi_branch == 0 ) {
	cout << "Branch pseudo_ecalPhi does not exist." << endl;
	}
	pseudo_eta_branch = 0;
	if (tree->GetAlias("pseudo_eta") != 0) {
		pseudo_eta_branch = tree->GetBranch(tree->GetAlias("pseudo_eta"));
		pseudo_eta_branch->SetAddress(&pseudo_eta_);
	}
	if(pseudo_eta_branch == 0 ) {
	cout << "Branch pseudo_eta does not exist." << endl;
	}
	pseudo_hcalD1Iso03_branch = 0;
	if (tree->GetAlias("pseudo_hcalD1Iso03") != 0) {
		pseudo_hcalD1Iso03_branch = tree->GetBranch(tree->GetAlias("pseudo_hcalD1Iso03"));
		pseudo_hcalD1Iso03_branch->SetAddress(&pseudo_hcalD1Iso03_);
	}
	if(pseudo_hcalD1Iso03_branch == 0 ) {
	cout << "Branch pseudo_hcalD1Iso03 does not exist." << endl;
	}
	pseudo_hcalD2Iso03_branch = 0;
	if (tree->GetAlias("pseudo_hcalD2Iso03") != 0) {
		pseudo_hcalD2Iso03_branch = tree->GetBranch(tree->GetAlias("pseudo_hcalD2Iso03"));
		pseudo_hcalD2Iso03_branch->SetAddress(&pseudo_hcalD2Iso03_);
	}
	if(pseudo_hcalD2Iso03_branch == 0 ) {
	cout << "Branch pseudo_hcalD2Iso03 does not exist." << endl;
	}
	pseudo_phi_branch = 0;
	if (tree->GetAlias("pseudo_phi") != 0) {
		pseudo_phi_branch = tree->GetBranch(tree->GetAlias("pseudo_phi"));
		pseudo_phi_branch->SetAddress(&pseudo_phi_);
	}
	if(pseudo_phi_branch == 0 ) {
	cout << "Branch pseudo_phi does not exist." << endl;
	}
	pseudo_tkIso03_branch = 0;
	if (tree->GetAlias("pseudo_tkIso03") != 0) {
		pseudo_tkIso03_branch = tree->GetBranch(tree->GetAlias("pseudo_tkIso03"));
		pseudo_tkIso03_branch->SetAddress(&pseudo_tkIso03_);
	}
	if(pseudo_tkIso03_branch == 0 ) {
	cout << "Branch pseudo_tkIso03 does not exist." << endl;
	}
	pseudo_towerEmEt_branch = 0;
	if (tree->GetAlias("pseudo_towerEmEt") != 0) {
		pseudo_towerEmEt_branch = tree->GetBranch(tree->GetAlias("pseudo_towerEmEt"));
		pseudo_towerEmEt_branch->SetAddress(&pseudo_towerEmEt_);
	}
	if(pseudo_towerEmEt_branch == 0 ) {
	cout << "Branch pseudo_towerEmEt does not exist." << endl;
	}
	pseudo_towerHadEt_branch = 0;
	if (tree->GetAlias("pseudo_towerHadEt") != 0) {
		pseudo_towerHadEt_branch = tree->GetBranch(tree->GetAlias("pseudo_towerHadEt"));
		pseudo_towerHadEt_branch->SetAddress(&pseudo_towerHadEt_);
	}
	if(pseudo_towerHadEt_branch == 0 ) {
	cout << "Branch pseudo_towerHadEt does not exist." << endl;
	}
	scs_clustersSize_branch = 0;
	if (tree->GetAlias("scs_clustersSize") != 0) {
		scs_clustersSize_branch = tree->GetBranch(tree->GetAlias("scs_clustersSize"));
		scs_clustersSize_branch->SetAddress(&scs_clustersSize_);
	}
	if(scs_clustersSize_branch == 0 ) {
	cout << "Branch scs_clustersSize does not exist." << endl;
	}
	scs_crystalsSize_branch = 0;
	if (tree->GetAlias("scs_crystalsSize") != 0) {
		scs_crystalsSize_branch = tree->GetBranch(tree->GetAlias("scs_crystalsSize"));
		scs_crystalsSize_branch->SetAddress(&scs_crystalsSize_);
	}
	if(scs_crystalsSize_branch == 0 ) {
	cout << "Branch scs_crystalsSize does not exist." << endl;
	}
	scs_e1x3_branch = 0;
	if (tree->GetAlias("scs_e1x3") != 0) {
		scs_e1x3_branch = tree->GetBranch(tree->GetAlias("scs_e1x3"));
		scs_e1x3_branch->SetAddress(&scs_e1x3_);
	}
	if(scs_e1x3_branch == 0 ) {
	cout << "Branch scs_e1x3 does not exist." << endl;
	}
	scs_e1x5_branch = 0;
	if (tree->GetAlias("scs_e1x5") != 0) {
		scs_e1x5_branch = tree->GetBranch(tree->GetAlias("scs_e1x5"));
		scs_e1x5_branch->SetAddress(&scs_e1x5_);
	}
	if(scs_e1x5_branch == 0 ) {
	cout << "Branch scs_e1x5 does not exist." << endl;
	}
	scs_e2nd_branch = 0;
	if (tree->GetAlias("scs_e2nd") != 0) {
		scs_e2nd_branch = tree->GetBranch(tree->GetAlias("scs_e2nd"));
		scs_e2nd_branch->SetAddress(&scs_e2nd_);
	}
	if(scs_e2nd_branch == 0 ) {
	cout << "Branch scs_e2nd does not exist." << endl;
	}
	scs_e2x2_branch = 0;
	if (tree->GetAlias("scs_e2x2") != 0) {
		scs_e2x2_branch = tree->GetBranch(tree->GetAlias("scs_e2x2"));
		scs_e2x2_branch->SetAddress(&scs_e2x2_);
	}
	if(scs_e2x2_branch == 0 ) {
	cout << "Branch scs_e2x2 does not exist." << endl;
	}
	scs_e2x5Max_branch = 0;
	if (tree->GetAlias("scs_e2x5Max") != 0) {
		scs_e2x5Max_branch = tree->GetBranch(tree->GetAlias("scs_e2x5Max"));
		scs_e2x5Max_branch->SetAddress(&scs_e2x5Max_);
	}
	if(scs_e2x5Max_branch == 0 ) {
	cout << "Branch scs_e2x5Max does not exist." << endl;
	}
	scs_e3x1_branch = 0;
	if (tree->GetAlias("scs_e3x1") != 0) {
		scs_e3x1_branch = tree->GetBranch(tree->GetAlias("scs_e3x1"));
		scs_e3x1_branch->SetAddress(&scs_e3x1_);
	}
	if(scs_e3x1_branch == 0 ) {
	cout << "Branch scs_e3x1 does not exist." << endl;
	}
	scs_e3x2_branch = 0;
	if (tree->GetAlias("scs_e3x2") != 0) {
		scs_e3x2_branch = tree->GetBranch(tree->GetAlias("scs_e3x2"));
		scs_e3x2_branch->SetAddress(&scs_e3x2_);
	}
	if(scs_e3x2_branch == 0 ) {
	cout << "Branch scs_e3x2 does not exist." << endl;
	}
	scs_e3x3_branch = 0;
	if (tree->GetAlias("scs_e3x3") != 0) {
		scs_e3x3_branch = tree->GetBranch(tree->GetAlias("scs_e3x3"));
		scs_e3x3_branch->SetAddress(&scs_e3x3_);
	}
	if(scs_e3x3_branch == 0 ) {
	cout << "Branch scs_e3x3 does not exist." << endl;
	}
	scs_e4x4_branch = 0;
	if (tree->GetAlias("scs_e4x4") != 0) {
		scs_e4x4_branch = tree->GetBranch(tree->GetAlias("scs_e4x4"));
		scs_e4x4_branch->SetAddress(&scs_e4x4_);
	}
	if(scs_e4x4_branch == 0 ) {
	cout << "Branch scs_e4x4 does not exist." << endl;
	}
	scs_e5x5_branch = 0;
	if (tree->GetAlias("scs_e5x5") != 0) {
		scs_e5x5_branch = tree->GetBranch(tree->GetAlias("scs_e5x5"));
		scs_e5x5_branch->SetAddress(&scs_e5x5_);
	}
	if(scs_e5x5_branch == 0 ) {
	cout << "Branch scs_e5x5 does not exist." << endl;
	}
	scs_eMax_branch = 0;
	if (tree->GetAlias("scs_eMax") != 0) {
		scs_eMax_branch = tree->GetBranch(tree->GetAlias("scs_eMax"));
		scs_eMax_branch->SetAddress(&scs_eMax_);
	}
	if(scs_eMax_branch == 0 ) {
	cout << "Branch scs_eMax does not exist." << endl;
	}
	scs_eSeed_branch = 0;
	if (tree->GetAlias("scs_eSeed") != 0) {
		scs_eSeed_branch = tree->GetBranch(tree->GetAlias("scs_eSeed"));
		scs_eSeed_branch->SetAddress(&scs_eSeed_);
	}
	if(scs_eSeed_branch == 0 ) {
	cout << "Branch scs_eSeed does not exist." << endl;
	}
	scs_energy_branch = 0;
	if (tree->GetAlias("scs_energy") != 0) {
		scs_energy_branch = tree->GetBranch(tree->GetAlias("scs_energy"));
		scs_energy_branch->SetAddress(&scs_energy_);
	}
	if(scs_energy_branch == 0 ) {
	cout << "Branch scs_energy does not exist." << endl;
	}
	scs_eta_branch = 0;
	if (tree->GetAlias("scs_eta") != 0) {
		scs_eta_branch = tree->GetBranch(tree->GetAlias("scs_eta"));
		scs_eta_branch->SetAddress(&scs_eta_);
	}
	if(scs_eta_branch == 0 ) {
	cout << "Branch scs_eta does not exist." << endl;
	}
	scs_hoe_branch = 0;
	if (tree->GetAlias("scs_hoe") != 0) {
		scs_hoe_branch = tree->GetBranch(tree->GetAlias("scs_hoe"));
		scs_hoe_branch->SetAddress(&scs_hoe_);
	}
	if(scs_hoe_branch == 0 ) {
	cout << "Branch scs_hoe does not exist." << endl;
	}
	scs_phi_branch = 0;
	if (tree->GetAlias("scs_phi") != 0) {
		scs_phi_branch = tree->GetBranch(tree->GetAlias("scs_phi"));
		scs_phi_branch->SetAddress(&scs_phi_);
	}
	if(scs_phi_branch == 0 ) {
	cout << "Branch scs_phi does not exist." << endl;
	}
	scs_preshowerEnergy_branch = 0;
	if (tree->GetAlias("scs_preshowerEnergy") != 0) {
		scs_preshowerEnergy_branch = tree->GetBranch(tree->GetAlias("scs_preshowerEnergy"));
		scs_preshowerEnergy_branch->SetAddress(&scs_preshowerEnergy_);
	}
	if(scs_preshowerEnergy_branch == 0 ) {
	cout << "Branch scs_preshowerEnergy does not exist." << endl;
	}
	scs_rawEnergy_branch = 0;
	if (tree->GetAlias("scs_rawEnergy") != 0) {
		scs_rawEnergy_branch = tree->GetBranch(tree->GetAlias("scs_rawEnergy"));
		scs_rawEnergy_branch->SetAddress(&scs_rawEnergy_);
	}
	if(scs_rawEnergy_branch == 0 ) {
	cout << "Branch scs_rawEnergy does not exist." << endl;
	}
	scs_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("scs_sigmaEtaEta") != 0) {
		scs_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("scs_sigmaEtaEta"));
		scs_sigmaEtaEta_branch->SetAddress(&scs_sigmaEtaEta_);
	}
	if(scs_sigmaEtaEta_branch == 0 ) {
	cout << "Branch scs_sigmaEtaEta does not exist." << endl;
	}
	scs_sigmaEtaPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaEtaPhi") != 0) {
		scs_sigmaEtaPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaEtaPhi"));
		scs_sigmaEtaPhi_branch->SetAddress(&scs_sigmaEtaPhi_);
	}
	if(scs_sigmaEtaPhi_branch == 0 ) {
	cout << "Branch scs_sigmaEtaPhi does not exist." << endl;
	}
	scs_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIEta") != 0) {
		scs_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIEta"));
		scs_sigmaIEtaIEta_branch->SetAddress(&scs_sigmaIEtaIEta_);
	}
	if(scs_sigmaIEtaIEta_branch == 0 ) {
	cout << "Branch scs_sigmaIEtaIEta does not exist." << endl;
	}
	scs_sigmaIEtaIPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIPhi") != 0) {
		scs_sigmaIEtaIPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIPhi"));
		scs_sigmaIEtaIPhi_branch->SetAddress(&scs_sigmaIEtaIPhi_);
	}
	if(scs_sigmaIEtaIPhi_branch == 0 ) {
	cout << "Branch scs_sigmaIEtaIPhi does not exist." << endl;
	}
	scs_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaIPhiIPhi") != 0) {
		scs_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIPhiIPhi"));
		scs_sigmaIPhiIPhi_branch->SetAddress(&scs_sigmaIPhiIPhi_);
	}
	if(scs_sigmaIPhiIPhi_branch == 0 ) {
	cout << "Branch scs_sigmaIPhiIPhi does not exist." << endl;
	}
	scs_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaPhiPhi") != 0) {
		scs_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaPhiPhi"));
		scs_sigmaPhiPhi_branch->SetAddress(&scs_sigmaPhiPhi_);
	}
	if(scs_sigmaPhiPhi_branch == 0 ) {
	cout << "Branch scs_sigmaPhiPhi does not exist." << endl;
	}
	trks_chi2_branch = 0;
	if (tree->GetAlias("trks_chi2") != 0) {
		trks_chi2_branch = tree->GetBranch(tree->GetAlias("trks_chi2"));
		trks_chi2_branch->SetAddress(&trks_chi2_);
	}
	if(trks_chi2_branch == 0 ) {
	cout << "Branch trks_chi2 does not exist." << endl;
	}
	trks_d0_branch = 0;
	if (tree->GetAlias("trks_d0") != 0) {
		trks_d0_branch = tree->GetBranch(tree->GetAlias("trks_d0"));
		trks_d0_branch->SetAddress(&trks_d0_);
	}
	if(trks_d0_branch == 0 ) {
	cout << "Branch trks_d0 does not exist." << endl;
	}
	trks_d0Err_branch = 0;
	if (tree->GetAlias("trks_d0Err") != 0) {
		trks_d0Err_branch = tree->GetBranch(tree->GetAlias("trks_d0Err"));
		trks_d0Err_branch->SetAddress(&trks_d0Err_);
	}
	if(trks_d0Err_branch == 0 ) {
	cout << "Branch trks_d0Err does not exist." << endl;
	}
	trks_d0corr_branch = 0;
	if (tree->GetAlias("trks_d0corr") != 0) {
		trks_d0corr_branch = tree->GetBranch(tree->GetAlias("trks_d0corr"));
		trks_d0corr_branch->SetAddress(&trks_d0corr_);
	}
	if(trks_d0corr_branch == 0 ) {
	cout << "Branch trks_d0corr does not exist." << endl;
	}
	trks_d0corrPhi_branch = 0;
	if (tree->GetAlias("trks_d0corrPhi") != 0) {
		trks_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("trks_d0corrPhi"));
		trks_d0corrPhi_branch->SetAddress(&trks_d0corrPhi_);
	}
	if(trks_d0corrPhi_branch == 0 ) {
	cout << "Branch trks_d0corrPhi does not exist." << endl;
	}
	trks_etaErr_branch = 0;
	if (tree->GetAlias("trks_etaErr") != 0) {
		trks_etaErr_branch = tree->GetBranch(tree->GetAlias("trks_etaErr"));
		trks_etaErr_branch->SetAddress(&trks_etaErr_);
	}
	if(trks_etaErr_branch == 0 ) {
	cout << "Branch trks_etaErr does not exist." << endl;
	}
	trks_layer1_charge_branch = 0;
	if (tree->GetAlias("trks_layer1_charge") != 0) {
		trks_layer1_charge_branch = tree->GetBranch(tree->GetAlias("trks_layer1_charge"));
		trks_layer1_charge_branch->SetAddress(&trks_layer1_charge_);
	}
	if(trks_layer1_charge_branch == 0 ) {
	cout << "Branch trks_layer1_charge does not exist." << endl;
	}
	trks_ndof_branch = 0;
	if (tree->GetAlias("trks_ndof") != 0) {
		trks_ndof_branch = tree->GetBranch(tree->GetAlias("trks_ndof"));
		trks_ndof_branch->SetAddress(&trks_ndof_);
	}
	if(trks_ndof_branch == 0 ) {
	cout << "Branch trks_ndof does not exist." << endl;
	}
	trks_phiErr_branch = 0;
	if (tree->GetAlias("trks_phiErr") != 0) {
		trks_phiErr_branch = tree->GetBranch(tree->GetAlias("trks_phiErr"));
		trks_phiErr_branch->SetAddress(&trks_phiErr_);
	}
	if(trks_phiErr_branch == 0 ) {
	cout << "Branch trks_phiErr does not exist." << endl;
	}
	trks_ptErr_branch = 0;
	if (tree->GetAlias("trks_ptErr") != 0) {
		trks_ptErr_branch = tree->GetBranch(tree->GetAlias("trks_ptErr"));
		trks_ptErr_branch->SetAddress(&trks_ptErr_);
	}
	if(trks_ptErr_branch == 0 ) {
	cout << "Branch trks_ptErr does not exist." << endl;
	}
	trks_z0_branch = 0;
	if (tree->GetAlias("trks_z0") != 0) {
		trks_z0_branch = tree->GetBranch(tree->GetAlias("trks_z0"));
		trks_z0_branch->SetAddress(&trks_z0_);
	}
	if(trks_z0_branch == 0 ) {
	cout << "Branch trks_z0 does not exist." << endl;
	}
	trks_z0Err_branch = 0;
	if (tree->GetAlias("trks_z0Err") != 0) {
		trks_z0Err_branch = tree->GetBranch(tree->GetAlias("trks_z0Err"));
		trks_z0Err_branch->SetAddress(&trks_z0Err_);
	}
	if(trks_z0Err_branch == 0 ) {
	cout << "Branch trks_z0Err does not exist." << endl;
	}
	trks_z0corr_branch = 0;
	if (tree->GetAlias("trks_z0corr") != 0) {
		trks_z0corr_branch = tree->GetBranch(tree->GetAlias("trks_z0corr"));
		trks_z0corr_branch->SetAddress(&trks_z0corr_);
	}
	if(trks_z0corr_branch == 0 ) {
	cout << "Branch trks_z0corr does not exist." << endl;
	}
	vtxs_chi2_branch = 0;
	if (tree->GetAlias("vtxs_chi2") != 0) {
		vtxs_chi2_branch = tree->GetBranch(tree->GetAlias("vtxs_chi2"));
		vtxs_chi2_branch->SetAddress(&vtxs_chi2_);
	}
	if(vtxs_chi2_branch == 0 ) {
	cout << "Branch vtxs_chi2 does not exist." << endl;
	}
	vtxs_ndof_branch = 0;
	if (tree->GetAlias("vtxs_ndof") != 0) {
		vtxs_ndof_branch = tree->GetBranch(tree->GetAlias("vtxs_ndof"));
		vtxs_ndof_branch->SetAddress(&vtxs_ndof_);
	}
	if(vtxs_ndof_branch == 0 ) {
	cout << "Branch vtxs_ndof does not exist." << endl;
	}
	vtxs_xError_branch = 0;
	if (tree->GetAlias("vtxs_xError") != 0) {
		vtxs_xError_branch = tree->GetBranch(tree->GetAlias("vtxs_xError"));
		vtxs_xError_branch->SetAddress(&vtxs_xError_);
	}
	if(vtxs_xError_branch == 0 ) {
	cout << "Branch vtxs_xError does not exist." << endl;
	}
	vtxs_yError_branch = 0;
	if (tree->GetAlias("vtxs_yError") != 0) {
		vtxs_yError_branch = tree->GetBranch(tree->GetAlias("vtxs_yError"));
		vtxs_yError_branch->SetAddress(&vtxs_yError_);
	}
	if(vtxs_yError_branch == 0 ) {
	cout << "Branch vtxs_yError does not exist." << endl;
	}
	vtxs_zError_branch = 0;
	if (tree->GetAlias("vtxs_zError") != 0) {
		vtxs_zError_branch = tree->GetBranch(tree->GetAlias("vtxs_zError"));
		vtxs_zError_branch->SetAddress(&vtxs_zError_);
	}
	if(vtxs_zError_branch == 0 ) {
	cout << "Branch vtxs_zError does not exist." << endl;
	}
	pseudo_ecalIso03_recHitE_branch = 0;
	if (tree->GetAlias("pseudo_ecalIso03_recHitE") != 0) {
		pseudo_ecalIso03_recHitE_branch = tree->GetBranch(tree->GetAlias("pseudo_ecalIso03_recHitE"));
		pseudo_ecalIso03_recHitE_branch->SetAddress(&pseudo_ecalIso03_recHitE_);
	}
	if(pseudo_ecalIso03_recHitE_branch == 0 ) {
	cout << "Branch pseudo_ecalIso03_recHitE does not exist." << endl;
	}
	pseudo_ecalIso03_recHitEt_branch = 0;
	if (tree->GetAlias("pseudo_ecalIso03_recHitEt") != 0) {
		pseudo_ecalIso03_recHitEt_branch = tree->GetBranch(tree->GetAlias("pseudo_ecalIso03_recHitEt"));
		pseudo_ecalIso03_recHitEt_branch->SetAddress(&pseudo_ecalIso03_recHitEt_);
	}
	if(pseudo_ecalIso03_recHitEt_branch == 0 ) {
	cout << "Branch pseudo_ecalIso03_recHitEt does not exist." << endl;
	}
	pseudo_srDIdx_branch = 0;
	if (tree->GetAlias("pseudo_srDIdx") != 0) {
		pseudo_srDIdx_branch = tree->GetBranch(tree->GetAlias("pseudo_srDIdx"));
		pseudo_srDIdx_branch->SetAddress(&pseudo_srDIdx_);
	}
	if(pseudo_srDIdx_branch == 0 ) {
	cout << "Branch pseudo_srDIdx does not exist." << endl;
	}
	vtxs_covMatrix_branch = 0;
	if (tree->GetAlias("vtxs_covMatrix") != 0) {
		vtxs_covMatrix_branch = tree->GetBranch(tree->GetAlias("vtxs_covMatrix"));
		vtxs_covMatrix_branch->SetAddress(&vtxs_covMatrix_);
	}
	if(vtxs_covMatrix_branch == 0 ) {
	cout << "Branch vtxs_covMatrix does not exist." << endl;
	}
	evt_bunchCrossing_branch = 0;
	if (tree->GetAlias("evt_bunchCrossing") != 0) {
		evt_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("evt_bunchCrossing"));
		evt_bunchCrossing_branch->SetAddress(&evt_bunchCrossing_);
	}
	if(evt_bunchCrossing_branch == 0 ) {
	cout << "Branch evt_bunchCrossing does not exist." << endl;
	}
	evt_experimentType_branch = 0;
	if (tree->GetAlias("evt_experimentType") != 0) {
		evt_experimentType_branch = tree->GetBranch(tree->GetAlias("evt_experimentType"));
		evt_experimentType_branch->SetAddress(&evt_experimentType_);
	}
	if(evt_experimentType_branch == 0 ) {
	cout << "Branch evt_experimentType does not exist." << endl;
	}
	evt_orbitNumber_branch = 0;
	if (tree->GetAlias("evt_orbitNumber") != 0) {
		evt_orbitNumber_branch = tree->GetBranch(tree->GetAlias("evt_orbitNumber"));
		evt_orbitNumber_branch->SetAddress(&evt_orbitNumber_);
	}
	if(evt_orbitNumber_branch == 0 ) {
	cout << "Branch evt_orbitNumber does not exist." << endl;
	}
	evt_storeNumber_branch = 0;
	if (tree->GetAlias("evt_storeNumber") != 0) {
		evt_storeNumber_branch = tree->GetBranch(tree->GetAlias("evt_storeNumber"));
		evt_storeNumber_branch->SetAddress(&evt_storeNumber_);
	}
	if(evt_storeNumber_branch == 0 ) {
	cout << "Branch evt_storeNumber does not exist." << endl;
	}
	hcalnoise_maxHPDHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxHPDHits") != 0) {
		hcalnoise_maxHPDHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxHPDHits"));
		hcalnoise_maxHPDHits_branch->SetAddress(&hcalnoise_maxHPDHits_);
	}
	if(hcalnoise_maxHPDHits_branch == 0 ) {
	cout << "Branch hcalnoise_maxHPDHits does not exist." << endl;
	}
	hcalnoise_maxRBXHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxRBXHits") != 0) {
		hcalnoise_maxRBXHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxRBXHits"));
		hcalnoise_maxRBXHits_branch->SetAddress(&hcalnoise_maxRBXHits_);
	}
	if(hcalnoise_maxRBXHits_branch == 0 ) {
	cout << "Branch hcalnoise_maxRBXHits does not exist." << endl;
	}
	hcalnoise_maxZeros_branch = 0;
	if (tree->GetAlias("hcalnoise_maxZeros") != 0) {
		hcalnoise_maxZeros_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxZeros"));
		hcalnoise_maxZeros_branch->SetAddress(&hcalnoise_maxZeros_);
	}
	if(hcalnoise_maxZeros_branch == 0 ) {
	cout << "Branch hcalnoise_maxZeros does not exist." << endl;
	}
	hcalnoise_noiseFilterStatus_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseFilterStatus") != 0) {
		hcalnoise_noiseFilterStatus_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseFilterStatus"));
		hcalnoise_noiseFilterStatus_branch->SetAddress(&hcalnoise_noiseFilterStatus_);
	}
	if(hcalnoise_noiseFilterStatus_branch == 0 ) {
	cout << "Branch hcalnoise_noiseFilterStatus does not exist." << endl;
	}
	hcalnoise_noiseType_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseType") != 0) {
		hcalnoise_noiseType_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseType"));
		hcalnoise_noiseType_branch->SetAddress(&hcalnoise_noiseType_);
	}
	if(hcalnoise_noiseType_branch == 0 ) {
	cout << "Branch hcalnoise_noiseType does not exist." << endl;
	}
	hcalnoise_num10GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num10GeVHits") != 0) {
		hcalnoise_num10GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num10GeVHits"));
		hcalnoise_num10GeVHits_branch->SetAddress(&hcalnoise_num10GeVHits_);
	}
	if(hcalnoise_num10GeVHits_branch == 0 ) {
	cout << "Branch hcalnoise_num10GeVHits does not exist." << endl;
	}
	hcalnoise_num25GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num25GeVHits") != 0) {
		hcalnoise_num25GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num25GeVHits"));
		hcalnoise_num25GeVHits_branch->SetAddress(&hcalnoise_num25GeVHits_);
	}
	if(hcalnoise_num25GeVHits_branch == 0 ) {
	cout << "Branch hcalnoise_num25GeVHits does not exist." << endl;
	}
	hcalnoise_numProblematicRBXs_branch = 0;
	if (tree->GetAlias("hcalnoise_numProblematicRBXs") != 0) {
		hcalnoise_numProblematicRBXs_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numProblematicRBXs"));
		hcalnoise_numProblematicRBXs_branch->SetAddress(&hcalnoise_numProblematicRBXs_);
	}
	if(hcalnoise_numProblematicRBXs_branch == 0 ) {
	cout << "Branch hcalnoise_numProblematicRBXs does not exist." << endl;
	}
	hcalnoise_passHighLevelNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passHighLevelNoiseFilter") != 0) {
		hcalnoise_passHighLevelNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passHighLevelNoiseFilter"));
		hcalnoise_passHighLevelNoiseFilter_branch->SetAddress(&hcalnoise_passHighLevelNoiseFilter_);
	}
	if(hcalnoise_passHighLevelNoiseFilter_branch == 0 ) {
	cout << "Branch hcalnoise_passHighLevelNoiseFilter does not exist." << endl;
	}
	hcalnoise_passLooseNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passLooseNoiseFilter") != 0) {
		hcalnoise_passLooseNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passLooseNoiseFilter"));
		hcalnoise_passLooseNoiseFilter_branch->SetAddress(&hcalnoise_passLooseNoiseFilter_);
	}
	if(hcalnoise_passLooseNoiseFilter_branch == 0 ) {
	cout << "Branch hcalnoise_passLooseNoiseFilter does not exist." << endl;
	}
	hcalnoise_passTightNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passTightNoiseFilter") != 0) {
		hcalnoise_passTightNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passTightNoiseFilter"));
		hcalnoise_passTightNoiseFilter_branch->SetAddress(&hcalnoise_passTightNoiseFilter_);
	}
	if(hcalnoise_passTightNoiseFilter_branch == 0 ) {
	cout << "Branch hcalnoise_passTightNoiseFilter does not exist." << endl;
	}
	l1_nemiso_branch = 0;
	if (tree->GetAlias("l1_nemiso") != 0) {
		l1_nemiso_branch = tree->GetBranch(tree->GetAlias("l1_nemiso"));
		l1_nemiso_branch->SetAddress(&l1_nemiso_);
	}
	if(l1_nemiso_branch == 0 ) {
	cout << "Branch l1_nemiso does not exist." << endl;
	}
	l1_nemnoiso_branch = 0;
	if (tree->GetAlias("l1_nemnoiso") != 0) {
		l1_nemnoiso_branch = tree->GetBranch(tree->GetAlias("l1_nemnoiso"));
		l1_nemnoiso_branch->SetAddress(&l1_nemnoiso_);
	}
	if(l1_nemnoiso_branch == 0 ) {
	cout << "Branch l1_nemnoiso does not exist." << endl;
	}
	l1_njetsc_branch = 0;
	if (tree->GetAlias("l1_njetsc") != 0) {
		l1_njetsc_branch = tree->GetBranch(tree->GetAlias("l1_njetsc"));
		l1_njetsc_branch->SetAddress(&l1_njetsc_);
	}
	if(l1_njetsc_branch == 0 ) {
	cout << "Branch l1_njetsc does not exist." << endl;
	}
	l1_njetsf_branch = 0;
	if (tree->GetAlias("l1_njetsf") != 0) {
		l1_njetsf_branch = tree->GetBranch(tree->GetAlias("l1_njetsf"));
		l1_njetsf_branch->SetAddress(&l1_njetsf_);
	}
	if(l1_njetsf_branch == 0 ) {
	cout << "Branch l1_njetsf does not exist." << endl;
	}
	l1_njetst_branch = 0;
	if (tree->GetAlias("l1_njetst") != 0) {
		l1_njetst_branch = tree->GetBranch(tree->GetAlias("l1_njetst"));
		l1_njetst_branch->SetAddress(&l1_njetst_);
	}
	if(l1_njetst_branch == 0 ) {
	cout << "Branch l1_njetst does not exist." << endl;
	}
	l1_nmus_branch = 0;
	if (tree->GetAlias("l1_nmus") != 0) {
		l1_nmus_branch = tree->GetBranch(tree->GetAlias("l1_nmus"));
		l1_nmus_branch->SetAddress(&l1_nmus_);
	}
	if(l1_nmus_branch == 0 ) {
	cout << "Branch l1_nmus does not exist." << endl;
	}
	l1_emiso_ieta_branch = 0;
	if (tree->GetAlias("l1_emiso_ieta") != 0) {
		l1_emiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1_emiso_ieta"));
		l1_emiso_ieta_branch->SetAddress(&l1_emiso_ieta_);
	}
	if(l1_emiso_ieta_branch == 0 ) {
	cout << "Branch l1_emiso_ieta does not exist." << endl;
	}
	l1_emiso_iphi_branch = 0;
	if (tree->GetAlias("l1_emiso_iphi") != 0) {
		l1_emiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1_emiso_iphi"));
		l1_emiso_iphi_branch->SetAddress(&l1_emiso_iphi_);
	}
	if(l1_emiso_iphi_branch == 0 ) {
	cout << "Branch l1_emiso_iphi does not exist." << endl;
	}
	l1_emiso_rawId_branch = 0;
	if (tree->GetAlias("l1_emiso_rawId") != 0) {
		l1_emiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1_emiso_rawId"));
		l1_emiso_rawId_branch->SetAddress(&l1_emiso_rawId_);
	}
	if(l1_emiso_rawId_branch == 0 ) {
	cout << "Branch l1_emiso_rawId does not exist." << endl;
	}
	l1_emiso_type_branch = 0;
	if (tree->GetAlias("l1_emiso_type") != 0) {
		l1_emiso_type_branch = tree->GetBranch(tree->GetAlias("l1_emiso_type"));
		l1_emiso_type_branch->SetAddress(&l1_emiso_type_);
	}
	if(l1_emiso_type_branch == 0 ) {
	cout << "Branch l1_emiso_type does not exist." << endl;
	}
	l1_emnoiso_ieta_branch = 0;
	if (tree->GetAlias("l1_emnoiso_ieta") != 0) {
		l1_emnoiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_ieta"));
		l1_emnoiso_ieta_branch->SetAddress(&l1_emnoiso_ieta_);
	}
	if(l1_emnoiso_ieta_branch == 0 ) {
	cout << "Branch l1_emnoiso_ieta does not exist." << endl;
	}
	l1_emnoiso_iphi_branch = 0;
	if (tree->GetAlias("l1_emnoiso_iphi") != 0) {
		l1_emnoiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_iphi"));
		l1_emnoiso_iphi_branch->SetAddress(&l1_emnoiso_iphi_);
	}
	if(l1_emnoiso_iphi_branch == 0 ) {
	cout << "Branch l1_emnoiso_iphi does not exist." << endl;
	}
	l1_emnoiso_rawId_branch = 0;
	if (tree->GetAlias("l1_emnoiso_rawId") != 0) {
		l1_emnoiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_rawId"));
		l1_emnoiso_rawId_branch->SetAddress(&l1_emnoiso_rawId_);
	}
	if(l1_emnoiso_rawId_branch == 0 ) {
	cout << "Branch l1_emnoiso_rawId does not exist." << endl;
	}
	l1_emnoiso_type_branch = 0;
	if (tree->GetAlias("l1_emnoiso_type") != 0) {
		l1_emnoiso_type_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_type"));
		l1_emnoiso_type_branch->SetAddress(&l1_emnoiso_type_);
	}
	if(l1_emnoiso_type_branch == 0 ) {
	cout << "Branch l1_emnoiso_type does not exist." << endl;
	}
	l1_jetsc_ieta_branch = 0;
	if (tree->GetAlias("l1_jetsc_ieta") != 0) {
		l1_jetsc_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_ieta"));
		l1_jetsc_ieta_branch->SetAddress(&l1_jetsc_ieta_);
	}
	if(l1_jetsc_ieta_branch == 0 ) {
	cout << "Branch l1_jetsc_ieta does not exist." << endl;
	}
	l1_jetsc_iphi_branch = 0;
	if (tree->GetAlias("l1_jetsc_iphi") != 0) {
		l1_jetsc_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_iphi"));
		l1_jetsc_iphi_branch->SetAddress(&l1_jetsc_iphi_);
	}
	if(l1_jetsc_iphi_branch == 0 ) {
	cout << "Branch l1_jetsc_iphi does not exist." << endl;
	}
	l1_jetsc_rawId_branch = 0;
	if (tree->GetAlias("l1_jetsc_rawId") != 0) {
		l1_jetsc_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_rawId"));
		l1_jetsc_rawId_branch->SetAddress(&l1_jetsc_rawId_);
	}
	if(l1_jetsc_rawId_branch == 0 ) {
	cout << "Branch l1_jetsc_rawId does not exist." << endl;
	}
	l1_jetsc_type_branch = 0;
	if (tree->GetAlias("l1_jetsc_type") != 0) {
		l1_jetsc_type_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_type"));
		l1_jetsc_type_branch->SetAddress(&l1_jetsc_type_);
	}
	if(l1_jetsc_type_branch == 0 ) {
	cout << "Branch l1_jetsc_type does not exist." << endl;
	}
	l1_jetsf_ieta_branch = 0;
	if (tree->GetAlias("l1_jetsf_ieta") != 0) {
		l1_jetsf_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_ieta"));
		l1_jetsf_ieta_branch->SetAddress(&l1_jetsf_ieta_);
	}
	if(l1_jetsf_ieta_branch == 0 ) {
	cout << "Branch l1_jetsf_ieta does not exist." << endl;
	}
	l1_jetsf_iphi_branch = 0;
	if (tree->GetAlias("l1_jetsf_iphi") != 0) {
		l1_jetsf_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_iphi"));
		l1_jetsf_iphi_branch->SetAddress(&l1_jetsf_iphi_);
	}
	if(l1_jetsf_iphi_branch == 0 ) {
	cout << "Branch l1_jetsf_iphi does not exist." << endl;
	}
	l1_jetsf_rawId_branch = 0;
	if (tree->GetAlias("l1_jetsf_rawId") != 0) {
		l1_jetsf_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_rawId"));
		l1_jetsf_rawId_branch->SetAddress(&l1_jetsf_rawId_);
	}
	if(l1_jetsf_rawId_branch == 0 ) {
	cout << "Branch l1_jetsf_rawId does not exist." << endl;
	}
	l1_jetsf_type_branch = 0;
	if (tree->GetAlias("l1_jetsf_type") != 0) {
		l1_jetsf_type_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_type"));
		l1_jetsf_type_branch->SetAddress(&l1_jetsf_type_);
	}
	if(l1_jetsf_type_branch == 0 ) {
	cout << "Branch l1_jetsf_type does not exist." << endl;
	}
	l1_jetst_ieta_branch = 0;
	if (tree->GetAlias("l1_jetst_ieta") != 0) {
		l1_jetst_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetst_ieta"));
		l1_jetst_ieta_branch->SetAddress(&l1_jetst_ieta_);
	}
	if(l1_jetst_ieta_branch == 0 ) {
	cout << "Branch l1_jetst_ieta does not exist." << endl;
	}
	l1_jetst_iphi_branch = 0;
	if (tree->GetAlias("l1_jetst_iphi") != 0) {
		l1_jetst_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetst_iphi"));
		l1_jetst_iphi_branch->SetAddress(&l1_jetst_iphi_);
	}
	if(l1_jetst_iphi_branch == 0 ) {
	cout << "Branch l1_jetst_iphi does not exist." << endl;
	}
	l1_jetst_rawId_branch = 0;
	if (tree->GetAlias("l1_jetst_rawId") != 0) {
		l1_jetst_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetst_rawId"));
		l1_jetst_rawId_branch->SetAddress(&l1_jetst_rawId_);
	}
	if(l1_jetst_rawId_branch == 0 ) {
	cout << "Branch l1_jetst_rawId does not exist." << endl;
	}
	l1_jetst_type_branch = 0;
	if (tree->GetAlias("l1_jetst_type") != 0) {
		l1_jetst_type_branch = tree->GetBranch(tree->GetAlias("l1_jetst_type"));
		l1_jetst_type_branch->SetAddress(&l1_jetst_type_);
	}
	if(l1_jetst_type_branch == 0 ) {
	cout << "Branch l1_jetst_type does not exist." << endl;
	}
	l1_mus_flags_branch = 0;
	if (tree->GetAlias("l1_mus_flags") != 0) {
		l1_mus_flags_branch = tree->GetBranch(tree->GetAlias("l1_mus_flags"));
		l1_mus_flags_branch->SetAddress(&l1_mus_flags_);
	}
	if(l1_mus_flags_branch == 0 ) {
	cout << "Branch l1_mus_flags does not exist." << endl;
	}
	l1_mus_q_branch = 0;
	if (tree->GetAlias("l1_mus_q") != 0) {
		l1_mus_q_branch = tree->GetBranch(tree->GetAlias("l1_mus_q"));
		l1_mus_q_branch->SetAddress(&l1_mus_q_);
	}
	if(l1_mus_q_branch == 0 ) {
	cout << "Branch l1_mus_q does not exist." << endl;
	}
	l1_mus_qual_branch = 0;
	if (tree->GetAlias("l1_mus_qual") != 0) {
		l1_mus_qual_branch = tree->GetBranch(tree->GetAlias("l1_mus_qual"));
		l1_mus_qual_branch->SetAddress(&l1_mus_qual_);
	}
	if(l1_mus_qual_branch == 0 ) {
	cout << "Branch l1_mus_qual does not exist." << endl;
	}
	l1_mus_qualFlags_branch = 0;
	if (tree->GetAlias("l1_mus_qualFlags") != 0) {
		l1_mus_qualFlags_branch = tree->GetBranch(tree->GetAlias("l1_mus_qualFlags"));
		l1_mus_qualFlags_branch->SetAddress(&l1_mus_qualFlags_);
	}
	if(l1_mus_qualFlags_branch == 0 ) {
	cout << "Branch l1_mus_qualFlags does not exist." << endl;
	}
	pxl_ndigis_pxb_branch = 0;
	if (tree->GetAlias("pxl_ndigis_pxb") != 0) {
		pxl_ndigis_pxb_branch = tree->GetBranch(tree->GetAlias("pxl_ndigis_pxb"));
		pxl_ndigis_pxb_branch->SetAddress(&pxl_ndigis_pxb_);
	}
	if(pxl_ndigis_pxb_branch == 0 ) {
	cout << "Branch pxl_ndigis_pxb does not exist." << endl;
	}
	pxl_ndigis_pxf_branch = 0;
	if (tree->GetAlias("pxl_ndigis_pxf") != 0) {
		pxl_ndigis_pxf_branch = tree->GetBranch(tree->GetAlias("pxl_ndigis_pxf"));
		pxl_ndigis_pxf_branch->SetAddress(&pxl_ndigis_pxf_);
	}
	if(pxl_ndigis_pxf_branch == 0 ) {
	cout << "Branch pxl_ndigis_pxf does not exist." << endl;
	}
	scs_elsidx_branch = 0;
	if (tree->GetAlias("scs_elsidx") != 0) {
		scs_elsidx_branch = tree->GetBranch(tree->GetAlias("scs_elsidx"));
		scs_elsidx_branch->SetAddress(&scs_elsidx_);
	}
	if(scs_elsidx_branch == 0 ) {
	cout << "Branch scs_elsidx does not exist." << endl;
	}
	scs_severitySeed_branch = 0;
	if (tree->GetAlias("scs_severitySeed") != 0) {
		scs_severitySeed_branch = tree->GetBranch(tree->GetAlias("scs_severitySeed"));
		scs_severitySeed_branch->SetAddress(&scs_severitySeed_);
	}
	if(scs_severitySeed_branch == 0 ) {
	cout << "Branch scs_severitySeed does not exist." << endl;
	}
	trks_algo_branch = 0;
	if (tree->GetAlias("trks_algo") != 0) {
		trks_algo_branch = tree->GetBranch(tree->GetAlias("trks_algo"));
		trks_algo_branch->SetAddress(&trks_algo_);
	}
	if(trks_algo_branch == 0 ) {
	cout << "Branch trks_algo does not exist." << endl;
	}
	trks_charge_branch = 0;
	if (tree->GetAlias("trks_charge") != 0) {
		trks_charge_branch = tree->GetBranch(tree->GetAlias("trks_charge"));
		trks_charge_branch->SetAddress(&trks_charge_);
	}
	if(trks_charge_branch == 0 ) {
	cout << "Branch trks_charge does not exist." << endl;
	}
	trks_exp_innerlayers_branch = 0;
	if (tree->GetAlias("trks_exp_innerlayers") != 0) {
		trks_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("trks_exp_innerlayers"));
		trks_exp_innerlayers_branch->SetAddress(&trks_exp_innerlayers_);
	}
	if(trks_exp_innerlayers_branch == 0 ) {
	cout << "Branch trks_exp_innerlayers does not exist." << endl;
	}
	trks_exp_outerlayers_branch = 0;
	if (tree->GetAlias("trks_exp_outerlayers") != 0) {
		trks_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("trks_exp_outerlayers"));
		trks_exp_outerlayers_branch->SetAddress(&trks_exp_outerlayers_);
	}
	if(trks_exp_outerlayers_branch == 0 ) {
	cout << "Branch trks_exp_outerlayers does not exist." << endl;
	}
	trks_layer1_det_branch = 0;
	if (tree->GetAlias("trks_layer1_det") != 0) {
		trks_layer1_det_branch = tree->GetBranch(tree->GetAlias("trks_layer1_det"));
		trks_layer1_det_branch->SetAddress(&trks_layer1_det_);
	}
	if(trks_layer1_det_branch == 0 ) {
	cout << "Branch trks_layer1_det does not exist." << endl;
	}
	trks_layer1_layer_branch = 0;
	if (tree->GetAlias("trks_layer1_layer") != 0) {
		trks_layer1_layer_branch = tree->GetBranch(tree->GetAlias("trks_layer1_layer"));
		trks_layer1_layer_branch->SetAddress(&trks_layer1_layer_);
	}
	if(trks_layer1_layer_branch == 0 ) {
	cout << "Branch trks_layer1_layer does not exist." << endl;
	}
	trks_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("trks_layer1_sizerphi") != 0) {
		trks_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("trks_layer1_sizerphi"));
		trks_layer1_sizerphi_branch->SetAddress(&trks_layer1_sizerphi_);
	}
	if(trks_layer1_sizerphi_branch == 0 ) {
	cout << "Branch trks_layer1_sizerphi does not exist." << endl;
	}
	trks_layer1_sizerz_branch = 0;
	if (tree->GetAlias("trks_layer1_sizerz") != 0) {
		trks_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("trks_layer1_sizerz"));
		trks_layer1_sizerz_branch->SetAddress(&trks_layer1_sizerz_);
	}
	if(trks_layer1_sizerz_branch == 0 ) {
	cout << "Branch trks_layer1_sizerz does not exist." << endl;
	}
	trks_lostHits_branch = 0;
	if (tree->GetAlias("trks_lostHits") != 0) {
		trks_lostHits_branch = tree->GetBranch(tree->GetAlias("trks_lostHits"));
		trks_lostHits_branch->SetAddress(&trks_lostHits_);
	}
	if(trks_lostHits_branch == 0 ) {
	cout << "Branch trks_lostHits does not exist." << endl;
	}
	trks_lost_pixelhits_branch = 0;
	if (tree->GetAlias("trks_lost_pixelhits") != 0) {
		trks_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("trks_lost_pixelhits"));
		trks_lost_pixelhits_branch->SetAddress(&trks_lost_pixelhits_);
	}
	if(trks_lost_pixelhits_branch == 0 ) {
	cout << "Branch trks_lost_pixelhits does not exist." << endl;
	}
	trks_nlayers_branch = 0;
	if (tree->GetAlias("trks_nlayers") != 0) {
		trks_nlayers_branch = tree->GetBranch(tree->GetAlias("trks_nlayers"));
		trks_nlayers_branch->SetAddress(&trks_nlayers_);
	}
	if(trks_nlayers_branch == 0 ) {
	cout << "Branch trks_nlayers does not exist." << endl;
	}
	trks_nlayers3D_branch = 0;
	if (tree->GetAlias("trks_nlayers3D") != 0) {
		trks_nlayers3D_branch = tree->GetBranch(tree->GetAlias("trks_nlayers3D"));
		trks_nlayers3D_branch->SetAddress(&trks_nlayers3D_);
	}
	if(trks_nlayers3D_branch == 0 ) {
	cout << "Branch trks_nlayers3D does not exist." << endl;
	}
	trks_nlayersLost_branch = 0;
	if (tree->GetAlias("trks_nlayersLost") != 0) {
		trks_nlayersLost_branch = tree->GetBranch(tree->GetAlias("trks_nlayersLost"));
		trks_nlayersLost_branch->SetAddress(&trks_nlayersLost_);
	}
	if(trks_nlayersLost_branch == 0 ) {
	cout << "Branch trks_nlayersLost does not exist." << endl;
	}
	trks_qualityMask_branch = 0;
	if (tree->GetAlias("trks_qualityMask") != 0) {
		trks_qualityMask_branch = tree->GetBranch(tree->GetAlias("trks_qualityMask"));
		trks_qualityMask_branch->SetAddress(&trks_qualityMask_);
	}
	if(trks_qualityMask_branch == 0 ) {
	cout << "Branch trks_qualityMask does not exist." << endl;
	}
	trks_validHits_branch = 0;
	if (tree->GetAlias("trks_validHits") != 0) {
		trks_validHits_branch = tree->GetBranch(tree->GetAlias("trks_validHits"));
		trks_validHits_branch->SetAddress(&trks_validHits_);
	}
	if(trks_validHits_branch == 0 ) {
	cout << "Branch trks_validHits does not exist." << endl;
	}
	trks_valid_pixelhits_branch = 0;
	if (tree->GetAlias("trks_valid_pixelhits") != 0) {
		trks_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("trks_valid_pixelhits"));
		trks_valid_pixelhits_branch->SetAddress(&trks_valid_pixelhits_);
	}
	if(trks_valid_pixelhits_branch == 0 ) {
	cout << "Branch trks_valid_pixelhits does not exist." << endl;
	}
	vtxs_isFake_branch = 0;
	if (tree->GetAlias("vtxs_isFake") != 0) {
		vtxs_isFake_branch = tree->GetBranch(tree->GetAlias("vtxs_isFake"));
		vtxs_isFake_branch->SetAddress(&vtxs_isFake_);
	}
	if(vtxs_isFake_branch == 0 ) {
	cout << "Branch vtxs_isFake does not exist." << endl;
	}
	vtxs_isValid_branch = 0;
	if (tree->GetAlias("vtxs_isValid") != 0) {
		vtxs_isValid_branch = tree->GetBranch(tree->GetAlias("vtxs_isValid"));
		vtxs_isValid_branch->SetAddress(&vtxs_isValid_);
	}
	if(vtxs_isValid_branch == 0 ) {
	cout << "Branch vtxs_isValid does not exist." << endl;
	}
	vtxs_tracksSize_branch = 0;
	if (tree->GetAlias("vtxs_tracksSize") != 0) {
		vtxs_tracksSize_branch = tree->GetBranch(tree->GetAlias("vtxs_tracksSize"));
		vtxs_tracksSize_branch->SetAddress(&vtxs_tracksSize_);
	}
	if(vtxs_tracksSize_branch == 0 ) {
	cout << "Branch vtxs_tracksSize does not exist." << endl;
	}
	pseudo_srFlags_branch = 0;
	if (tree->GetAlias("pseudo_srFlags") != 0) {
		pseudo_srFlags_branch = tree->GetBranch(tree->GetAlias("pseudo_srFlags"));
		pseudo_srFlags_branch->SetAddress(&pseudo_srFlags_);
	}
	if(pseudo_srFlags_branch == 0 ) {
	cout << "Branch pseudo_srFlags does not exist." << endl;
	}
	evt_ntwrs_branch = 0;
	if (tree->GetAlias("evt_ntwrs") != 0) {
		evt_ntwrs_branch = tree->GetBranch(tree->GetAlias("evt_ntwrs"));
		evt_ntwrs_branch->SetAddress(&evt_ntwrs_);
	}
	if(evt_ntwrs_branch == 0 ) {
	cout << "Branch evt_ntwrs does not exist." << endl;
	}
	evt_event_branch = 0;
	if (tree->GetAlias("evt_event") != 0) {
		evt_event_branch = tree->GetBranch(tree->GetAlias("evt_event"));
		evt_event_branch->SetAddress(&evt_event_);
	}
	if(evt_event_branch == 0 ) {
	cout << "Branch evt_event does not exist." << endl;
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetAlias("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch(tree->GetAlias("evt_lumiBlock"));
		evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);
	}
	if(evt_lumiBlock_branch == 0 ) {
	cout << "Branch evt_lumiBlock does not exist." << endl;
	}
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		evt_run_branch->SetAddress(&evt_run_);
	}
	if(evt_run_branch == 0 ) {
	cout << "Branch evt_run does not exist." << endl;
	}
	l1_bits1_branch = 0;
	if (tree->GetAlias("l1_bits1") != 0) {
		l1_bits1_branch = tree->GetBranch(tree->GetAlias("l1_bits1"));
		l1_bits1_branch->SetAddress(&l1_bits1_);
	}
	if(l1_bits1_branch == 0 ) {
	cout << "Branch l1_bits1 does not exist." << endl;
	}
	l1_bits2_branch = 0;
	if (tree->GetAlias("l1_bits2") != 0) {
		l1_bits2_branch = tree->GetBranch(tree->GetAlias("l1_bits2"));
		l1_bits2_branch->SetAddress(&l1_bits2_);
	}
	if(l1_bits2_branch == 0 ) {
	cout << "Branch l1_bits2 does not exist." << endl;
	}
	l1_bits3_branch = 0;
	if (tree->GetAlias("l1_bits3") != 0) {
		l1_bits3_branch = tree->GetBranch(tree->GetAlias("l1_bits3"));
		l1_bits3_branch->SetAddress(&l1_bits3_);
	}
	if(l1_bits3_branch == 0 ) {
	cout << "Branch l1_bits3 does not exist." << endl;
	}
	l1_bits4_branch = 0;
	if (tree->GetAlias("l1_bits4") != 0) {
		l1_bits4_branch = tree->GetBranch(tree->GetAlias("l1_bits4"));
		l1_bits4_branch->SetAddress(&l1_bits4_);
	}
	if(l1_bits4_branch == 0 ) {
	cout << "Branch l1_bits4 does not exist." << endl;
	}
	l1_techbits1_branch = 0;
	if (tree->GetAlias("l1_techbits1") != 0) {
		l1_techbits1_branch = tree->GetBranch(tree->GetAlias("l1_techbits1"));
		l1_techbits1_branch->SetAddress(&l1_techbits1_);
	}
	if(l1_techbits1_branch == 0 ) {
	cout << "Branch l1_techbits1 does not exist." << endl;
	}
	l1_techbits2_branch = 0;
	if (tree->GetAlias("l1_techbits2") != 0) {
		l1_techbits2_branch = tree->GetBranch(tree->GetAlias("l1_techbits2"));
		l1_techbits2_branch->SetAddress(&l1_techbits2_);
	}
	if(l1_techbits2_branch == 0 ) {
	cout << "Branch l1_techbits2 does not exist." << endl;
	}
	evt_nscs_branch = 0;
	if (tree->GetAlias("evt_nscs") != 0) {
		evt_nscs_branch = tree->GetBranch(tree->GetAlias("evt_nscs"));
		evt_nscs_branch->SetAddress(&evt_nscs_);
	}
	if(evt_nscs_branch == 0 ) {
	cout << "Branch evt_nscs does not exist." << endl;
	}
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	if(evt_nvtxs_branch == 0 ) {
	cout << "Branch evt_nvtxs does not exist." << endl;
	}
	twrs_numBadEcalCells_branch = 0;
	if (tree->GetAlias("twrs_numBadEcalCells") != 0) {
		twrs_numBadEcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numBadEcalCells"));
		twrs_numBadEcalCells_branch->SetAddress(&twrs_numBadEcalCells_);
	}
	if(twrs_numBadEcalCells_branch == 0 ) {
	cout << "Branch twrs_numBadEcalCells does not exist." << endl;
	}
	twrs_numBadHcalCells_branch = 0;
	if (tree->GetAlias("twrs_numBadHcalCells") != 0) {
		twrs_numBadHcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numBadHcalCells"));
		twrs_numBadHcalCells_branch->SetAddress(&twrs_numBadHcalCells_);
	}
	if(twrs_numBadHcalCells_branch == 0 ) {
	cout << "Branch twrs_numBadHcalCells does not exist." << endl;
	}
	twrs_numProblematicEcalCells_branch = 0;
	if (tree->GetAlias("twrs_numProblematicEcalCells") != 0) {
		twrs_numProblematicEcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numProblematicEcalCells"));
		twrs_numProblematicEcalCells_branch->SetAddress(&twrs_numProblematicEcalCells_);
	}
	if(twrs_numProblematicEcalCells_branch == 0 ) {
	cout << "Branch twrs_numProblematicEcalCells does not exist." << endl;
	}
	twrs_numProblematicHcalCells_branch = 0;
	if (tree->GetAlias("twrs_numProblematicHcalCells") != 0) {
		twrs_numProblematicHcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numProblematicHcalCells"));
		twrs_numProblematicHcalCells_branch->SetAddress(&twrs_numProblematicHcalCells_);
	}
	if(twrs_numProblematicHcalCells_branch == 0 ) {
	cout << "Branch twrs_numProblematicHcalCells does not exist." << endl;
	}
	twrs_numRecoveredEcalCells_branch = 0;
	if (tree->GetAlias("twrs_numRecoveredEcalCells") != 0) {
		twrs_numRecoveredEcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numRecoveredEcalCells"));
		twrs_numRecoveredEcalCells_branch->SetAddress(&twrs_numRecoveredEcalCells_);
	}
	if(twrs_numRecoveredEcalCells_branch == 0 ) {
	cout << "Branch twrs_numRecoveredEcalCells does not exist." << endl;
	}
	twrs_numRecoveredHcalCells_branch = 0;
	if (tree->GetAlias("twrs_numRecoveredHcalCells") != 0) {
		twrs_numRecoveredHcalCells_branch = tree->GetBranch(tree->GetAlias("twrs_numRecoveredHcalCells"));
		twrs_numRecoveredHcalCells_branch->SetAddress(&twrs_numRecoveredHcalCells_);
	}
	if(twrs_numRecoveredHcalCells_branch == 0 ) {
	cout << "Branch twrs_numRecoveredHcalCells does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		evt_CMS2tag_isLoaded = false;
		evt_dataset_isLoaded = false;
		l1_trigNames_isLoaded = false;
		evt_bField_isLoaded = false;
		hcalnoise_eventChargeFraction_isLoaded = false;
		hcalnoise_eventEMEnergy_isLoaded = false;
		hcalnoise_eventEMFraction_isLoaded = false;
		hcalnoise_eventHadEnergy_isLoaded = false;
		hcalnoise_eventTrackEnergy_isLoaded = false;
		hcalnoise_max10GeVHitTime_isLoaded = false;
		hcalnoise_max25GeVHitTime_isLoaded = false;
		hcalnoise_min10GeVHitTime_isLoaded = false;
		hcalnoise_min25GeVHitTime_isLoaded = false;
		hcalnoise_minE10TS_isLoaded = false;
		hcalnoise_minE2Over10TS_isLoaded = false;
		hcalnoise_minE2TS_isLoaded = false;
		hcalnoise_minHPDEMF_isLoaded = false;
		hcalnoise_minRBXEMF_isLoaded = false;
		hcalnoise_rms10GeVHitTime_isLoaded = false;
		hcalnoise_rms25GeVHitTime_isLoaded = false;
		l1_met_etTot_isLoaded = false;
		l1_met_met_isLoaded = false;
		l1_mht_htTot_isLoaded = false;
		l1_mht_mht_isLoaded = false;
		l1_met_p4_isLoaded = false;
		l1_mht_p4_isLoaded = false;
		l1_emiso_p4_isLoaded = false;
		l1_emnoiso_p4_isLoaded = false;
		l1_jetsc_p4_isLoaded = false;
		l1_jetsf_p4_isLoaded = false;
		l1_jetst_p4_isLoaded = false;
		l1_mus_p4_isLoaded = false;
		scs_p4_isLoaded = false;
		scs_pos_p4_isLoaded = false;
		scs_vtx_p4_isLoaded = false;
		trks_inner_position_isLoaded = false;
		trks_outer_p4_isLoaded = false;
		trks_outer_position_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		trks_vertex_p4_isLoaded = false;
		vtxs_position_isLoaded = false;
		twrs_ecalTime_isLoaded = false;
		twrs_emEnergy_isLoaded = false;
		twrs_emEt_isLoaded = false;
		twrs_emEtcorr_isLoaded = false;
		twrs_eta_isLoaded = false;
		twrs_etcorr_isLoaded = false;
		twrs_hadEnergy_isLoaded = false;
		twrs_hadEt_isLoaded = false;
		twrs_hadEtcorr_isLoaded = false;
		twrs_hcalTime_isLoaded = false;
		twrs_outerEnergy_isLoaded = false;
		twrs_outerEt_isLoaded = false;
		twrs_outerEtcorr_isLoaded = false;
		twrs_pcorr_isLoaded = false;
		twrs_phi_isLoaded = false;
		pseudo_dRClosestTower_isLoaded = false;
		pseudo_dRClosestTowerEmEt_isLoaded = false;
		pseudo_ecalEta_isLoaded = false;
		pseudo_ecalIso03_isLoaded = false;
		pseudo_ecalPhi_isLoaded = false;
		pseudo_eta_isLoaded = false;
		pseudo_hcalD1Iso03_isLoaded = false;
		pseudo_hcalD2Iso03_isLoaded = false;
		pseudo_phi_isLoaded = false;
		pseudo_tkIso03_isLoaded = false;
		pseudo_towerEmEt_isLoaded = false;
		pseudo_towerHadEt_isLoaded = false;
		scs_clustersSize_isLoaded = false;
		scs_crystalsSize_isLoaded = false;
		scs_e1x3_isLoaded = false;
		scs_e1x5_isLoaded = false;
		scs_e2nd_isLoaded = false;
		scs_e2x2_isLoaded = false;
		scs_e2x5Max_isLoaded = false;
		scs_e3x1_isLoaded = false;
		scs_e3x2_isLoaded = false;
		scs_e3x3_isLoaded = false;
		scs_e4x4_isLoaded = false;
		scs_e5x5_isLoaded = false;
		scs_eMax_isLoaded = false;
		scs_eSeed_isLoaded = false;
		scs_energy_isLoaded = false;
		scs_eta_isLoaded = false;
		scs_hoe_isLoaded = false;
		scs_phi_isLoaded = false;
		scs_preshowerEnergy_isLoaded = false;
		scs_rawEnergy_isLoaded = false;
		scs_sigmaEtaEta_isLoaded = false;
		scs_sigmaEtaPhi_isLoaded = false;
		scs_sigmaIEtaIEta_isLoaded = false;
		scs_sigmaIEtaIPhi_isLoaded = false;
		scs_sigmaIPhiIPhi_isLoaded = false;
		scs_sigmaPhiPhi_isLoaded = false;
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_d0corr_isLoaded = false;
		trks_d0corrPhi_isLoaded = false;
		trks_etaErr_isLoaded = false;
		trks_layer1_charge_isLoaded = false;
		trks_ndof_isLoaded = false;
		trks_phiErr_isLoaded = false;
		trks_ptErr_isLoaded = false;
		trks_z0_isLoaded = false;
		trks_z0Err_isLoaded = false;
		trks_z0corr_isLoaded = false;
		vtxs_chi2_isLoaded = false;
		vtxs_ndof_isLoaded = false;
		vtxs_xError_isLoaded = false;
		vtxs_yError_isLoaded = false;
		vtxs_zError_isLoaded = false;
		pseudo_ecalIso03_recHitE_isLoaded = false;
		pseudo_ecalIso03_recHitEt_isLoaded = false;
		pseudo_srDIdx_isLoaded = false;
		vtxs_covMatrix_isLoaded = false;
		evt_bunchCrossing_isLoaded = false;
		evt_experimentType_isLoaded = false;
		evt_orbitNumber_isLoaded = false;
		evt_storeNumber_isLoaded = false;
		hcalnoise_maxHPDHits_isLoaded = false;
		hcalnoise_maxRBXHits_isLoaded = false;
		hcalnoise_maxZeros_isLoaded = false;
		hcalnoise_noiseFilterStatus_isLoaded = false;
		hcalnoise_noiseType_isLoaded = false;
		hcalnoise_num10GeVHits_isLoaded = false;
		hcalnoise_num25GeVHits_isLoaded = false;
		hcalnoise_numProblematicRBXs_isLoaded = false;
		hcalnoise_passHighLevelNoiseFilter_isLoaded = false;
		hcalnoise_passLooseNoiseFilter_isLoaded = false;
		hcalnoise_passTightNoiseFilter_isLoaded = false;
		l1_nemiso_isLoaded = false;
		l1_nemnoiso_isLoaded = false;
		l1_njetsc_isLoaded = false;
		l1_njetsf_isLoaded = false;
		l1_njetst_isLoaded = false;
		l1_nmus_isLoaded = false;
		l1_emiso_ieta_isLoaded = false;
		l1_emiso_iphi_isLoaded = false;
		l1_emiso_rawId_isLoaded = false;
		l1_emiso_type_isLoaded = false;
		l1_emnoiso_ieta_isLoaded = false;
		l1_emnoiso_iphi_isLoaded = false;
		l1_emnoiso_rawId_isLoaded = false;
		l1_emnoiso_type_isLoaded = false;
		l1_jetsc_ieta_isLoaded = false;
		l1_jetsc_iphi_isLoaded = false;
		l1_jetsc_rawId_isLoaded = false;
		l1_jetsc_type_isLoaded = false;
		l1_jetsf_ieta_isLoaded = false;
		l1_jetsf_iphi_isLoaded = false;
		l1_jetsf_rawId_isLoaded = false;
		l1_jetsf_type_isLoaded = false;
		l1_jetst_ieta_isLoaded = false;
		l1_jetst_iphi_isLoaded = false;
		l1_jetst_rawId_isLoaded = false;
		l1_jetst_type_isLoaded = false;
		l1_mus_flags_isLoaded = false;
		l1_mus_q_isLoaded = false;
		l1_mus_qual_isLoaded = false;
		l1_mus_qualFlags_isLoaded = false;
		pxl_ndigis_pxb_isLoaded = false;
		pxl_ndigis_pxf_isLoaded = false;
		scs_elsidx_isLoaded = false;
		scs_severitySeed_isLoaded = false;
		trks_algo_isLoaded = false;
		trks_charge_isLoaded = false;
		trks_exp_innerlayers_isLoaded = false;
		trks_exp_outerlayers_isLoaded = false;
		trks_layer1_det_isLoaded = false;
		trks_layer1_layer_isLoaded = false;
		trks_layer1_sizerphi_isLoaded = false;
		trks_layer1_sizerz_isLoaded = false;
		trks_lostHits_isLoaded = false;
		trks_lost_pixelhits_isLoaded = false;
		trks_nlayers_isLoaded = false;
		trks_nlayers3D_isLoaded = false;
		trks_nlayersLost_isLoaded = false;
		trks_qualityMask_isLoaded = false;
		trks_validHits_isLoaded = false;
		trks_valid_pixelhits_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		vtxs_tracksSize_isLoaded = false;
		pseudo_srFlags_isLoaded = false;
		evt_ntwrs_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		l1_bits1_isLoaded = false;
		l1_bits2_isLoaded = false;
		l1_bits3_isLoaded = false;
		l1_bits4_isLoaded = false;
		l1_techbits1_isLoaded = false;
		l1_techbits2_isLoaded = false;
		evt_nscs_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		twrs_numBadEcalCells_isLoaded = false;
		twrs_numBadHcalCells_isLoaded = false;
		twrs_numProblematicEcalCells_isLoaded = false;
		twrs_numProblematicHcalCells_isLoaded = false;
		twrs_numRecoveredEcalCells_isLoaded = false;
		twrs_numRecoveredHcalCells_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (evt_CMS2tag_branch != 0) evt_CMS2tag();
	if (evt_dataset_branch != 0) evt_dataset();
	if (l1_trigNames_branch != 0) l1_trigNames();
	if (evt_bField_branch != 0) evt_bField();
	if (hcalnoise_eventChargeFraction_branch != 0) hcalnoise_eventChargeFraction();
	if (hcalnoise_eventEMEnergy_branch != 0) hcalnoise_eventEMEnergy();
	if (hcalnoise_eventEMFraction_branch != 0) hcalnoise_eventEMFraction();
	if (hcalnoise_eventHadEnergy_branch != 0) hcalnoise_eventHadEnergy();
	if (hcalnoise_eventTrackEnergy_branch != 0) hcalnoise_eventTrackEnergy();
	if (hcalnoise_max10GeVHitTime_branch != 0) hcalnoise_max10GeVHitTime();
	if (hcalnoise_max25GeVHitTime_branch != 0) hcalnoise_max25GeVHitTime();
	if (hcalnoise_min10GeVHitTime_branch != 0) hcalnoise_min10GeVHitTime();
	if (hcalnoise_min25GeVHitTime_branch != 0) hcalnoise_min25GeVHitTime();
	if (hcalnoise_minE10TS_branch != 0) hcalnoise_minE10TS();
	if (hcalnoise_minE2Over10TS_branch != 0) hcalnoise_minE2Over10TS();
	if (hcalnoise_minE2TS_branch != 0) hcalnoise_minE2TS();
	if (hcalnoise_minHPDEMF_branch != 0) hcalnoise_minHPDEMF();
	if (hcalnoise_minRBXEMF_branch != 0) hcalnoise_minRBXEMF();
	if (hcalnoise_rms10GeVHitTime_branch != 0) hcalnoise_rms10GeVHitTime();
	if (hcalnoise_rms25GeVHitTime_branch != 0) hcalnoise_rms25GeVHitTime();
	if (l1_met_etTot_branch != 0) l1_met_etTot();
	if (l1_met_met_branch != 0) l1_met_met();
	if (l1_mht_htTot_branch != 0) l1_mht_htTot();
	if (l1_mht_mht_branch != 0) l1_mht_mht();
	if (l1_met_p4_branch != 0) l1_met_p4();
	if (l1_mht_p4_branch != 0) l1_mht_p4();
	if (l1_emiso_p4_branch != 0) l1_emiso_p4();
	if (l1_emnoiso_p4_branch != 0) l1_emnoiso_p4();
	if (l1_jetsc_p4_branch != 0) l1_jetsc_p4();
	if (l1_jetsf_p4_branch != 0) l1_jetsf_p4();
	if (l1_jetst_p4_branch != 0) l1_jetst_p4();
	if (l1_mus_p4_branch != 0) l1_mus_p4();
	if (scs_p4_branch != 0) scs_p4();
	if (scs_pos_p4_branch != 0) scs_pos_p4();
	if (scs_vtx_p4_branch != 0) scs_vtx_p4();
	if (trks_inner_position_branch != 0) trks_inner_position();
	if (trks_outer_p4_branch != 0) trks_outer_p4();
	if (trks_outer_position_branch != 0) trks_outer_position();
	if (trks_trk_p4_branch != 0) trks_trk_p4();
	if (trks_vertex_p4_branch != 0) trks_vertex_p4();
	if (vtxs_position_branch != 0) vtxs_position();
	if (twrs_ecalTime_branch != 0) twrs_ecalTime();
	if (twrs_emEnergy_branch != 0) twrs_emEnergy();
	if (twrs_emEt_branch != 0) twrs_emEt();
	if (twrs_emEtcorr_branch != 0) twrs_emEtcorr();
	if (twrs_eta_branch != 0) twrs_eta();
	if (twrs_etcorr_branch != 0) twrs_etcorr();
	if (twrs_hadEnergy_branch != 0) twrs_hadEnergy();
	if (twrs_hadEt_branch != 0) twrs_hadEt();
	if (twrs_hadEtcorr_branch != 0) twrs_hadEtcorr();
	if (twrs_hcalTime_branch != 0) twrs_hcalTime();
	if (twrs_outerEnergy_branch != 0) twrs_outerEnergy();
	if (twrs_outerEt_branch != 0) twrs_outerEt();
	if (twrs_outerEtcorr_branch != 0) twrs_outerEtcorr();
	if (twrs_pcorr_branch != 0) twrs_pcorr();
	if (twrs_phi_branch != 0) twrs_phi();
	if (pseudo_dRClosestTower_branch != 0) pseudo_dRClosestTower();
	if (pseudo_dRClosestTowerEmEt_branch != 0) pseudo_dRClosestTowerEmEt();
	if (pseudo_ecalEta_branch != 0) pseudo_ecalEta();
	if (pseudo_ecalIso03_branch != 0) pseudo_ecalIso03();
	if (pseudo_ecalPhi_branch != 0) pseudo_ecalPhi();
	if (pseudo_eta_branch != 0) pseudo_eta();
	if (pseudo_hcalD1Iso03_branch != 0) pseudo_hcalD1Iso03();
	if (pseudo_hcalD2Iso03_branch != 0) pseudo_hcalD2Iso03();
	if (pseudo_phi_branch != 0) pseudo_phi();
	if (pseudo_tkIso03_branch != 0) pseudo_tkIso03();
	if (pseudo_towerEmEt_branch != 0) pseudo_towerEmEt();
	if (pseudo_towerHadEt_branch != 0) pseudo_towerHadEt();
	if (scs_clustersSize_branch != 0) scs_clustersSize();
	if (scs_crystalsSize_branch != 0) scs_crystalsSize();
	if (scs_e1x3_branch != 0) scs_e1x3();
	if (scs_e1x5_branch != 0) scs_e1x5();
	if (scs_e2nd_branch != 0) scs_e2nd();
	if (scs_e2x2_branch != 0) scs_e2x2();
	if (scs_e2x5Max_branch != 0) scs_e2x5Max();
	if (scs_e3x1_branch != 0) scs_e3x1();
	if (scs_e3x2_branch != 0) scs_e3x2();
	if (scs_e3x3_branch != 0) scs_e3x3();
	if (scs_e4x4_branch != 0) scs_e4x4();
	if (scs_e5x5_branch != 0) scs_e5x5();
	if (scs_eMax_branch != 0) scs_eMax();
	if (scs_eSeed_branch != 0) scs_eSeed();
	if (scs_energy_branch != 0) scs_energy();
	if (scs_eta_branch != 0) scs_eta();
	if (scs_hoe_branch != 0) scs_hoe();
	if (scs_phi_branch != 0) scs_phi();
	if (scs_preshowerEnergy_branch != 0) scs_preshowerEnergy();
	if (scs_rawEnergy_branch != 0) scs_rawEnergy();
	if (scs_sigmaEtaEta_branch != 0) scs_sigmaEtaEta();
	if (scs_sigmaEtaPhi_branch != 0) scs_sigmaEtaPhi();
	if (scs_sigmaIEtaIEta_branch != 0) scs_sigmaIEtaIEta();
	if (scs_sigmaIEtaIPhi_branch != 0) scs_sigmaIEtaIPhi();
	if (scs_sigmaIPhiIPhi_branch != 0) scs_sigmaIPhiIPhi();
	if (scs_sigmaPhiPhi_branch != 0) scs_sigmaPhiPhi();
	if (trks_chi2_branch != 0) trks_chi2();
	if (trks_d0_branch != 0) trks_d0();
	if (trks_d0Err_branch != 0) trks_d0Err();
	if (trks_d0corr_branch != 0) trks_d0corr();
	if (trks_d0corrPhi_branch != 0) trks_d0corrPhi();
	if (trks_etaErr_branch != 0) trks_etaErr();
	if (trks_layer1_charge_branch != 0) trks_layer1_charge();
	if (trks_ndof_branch != 0) trks_ndof();
	if (trks_phiErr_branch != 0) trks_phiErr();
	if (trks_ptErr_branch != 0) trks_ptErr();
	if (trks_z0_branch != 0) trks_z0();
	if (trks_z0Err_branch != 0) trks_z0Err();
	if (trks_z0corr_branch != 0) trks_z0corr();
	if (vtxs_chi2_branch != 0) vtxs_chi2();
	if (vtxs_ndof_branch != 0) vtxs_ndof();
	if (vtxs_xError_branch != 0) vtxs_xError();
	if (vtxs_yError_branch != 0) vtxs_yError();
	if (vtxs_zError_branch != 0) vtxs_zError();
	if (pseudo_ecalIso03_recHitE_branch != 0) pseudo_ecalIso03_recHitE();
	if (pseudo_ecalIso03_recHitEt_branch != 0) pseudo_ecalIso03_recHitEt();
	if (pseudo_srDIdx_branch != 0) pseudo_srDIdx();
	if (vtxs_covMatrix_branch != 0) vtxs_covMatrix();
	if (evt_bunchCrossing_branch != 0) evt_bunchCrossing();
	if (evt_experimentType_branch != 0) evt_experimentType();
	if (evt_orbitNumber_branch != 0) evt_orbitNumber();
	if (evt_storeNumber_branch != 0) evt_storeNumber();
	if (hcalnoise_maxHPDHits_branch != 0) hcalnoise_maxHPDHits();
	if (hcalnoise_maxRBXHits_branch != 0) hcalnoise_maxRBXHits();
	if (hcalnoise_maxZeros_branch != 0) hcalnoise_maxZeros();
	if (hcalnoise_noiseFilterStatus_branch != 0) hcalnoise_noiseFilterStatus();
	if (hcalnoise_noiseType_branch != 0) hcalnoise_noiseType();
	if (hcalnoise_num10GeVHits_branch != 0) hcalnoise_num10GeVHits();
	if (hcalnoise_num25GeVHits_branch != 0) hcalnoise_num25GeVHits();
	if (hcalnoise_numProblematicRBXs_branch != 0) hcalnoise_numProblematicRBXs();
	if (hcalnoise_passHighLevelNoiseFilter_branch != 0) hcalnoise_passHighLevelNoiseFilter();
	if (hcalnoise_passLooseNoiseFilter_branch != 0) hcalnoise_passLooseNoiseFilter();
	if (hcalnoise_passTightNoiseFilter_branch != 0) hcalnoise_passTightNoiseFilter();
	if (l1_nemiso_branch != 0) l1_nemiso();
	if (l1_nemnoiso_branch != 0) l1_nemnoiso();
	if (l1_njetsc_branch != 0) l1_njetsc();
	if (l1_njetsf_branch != 0) l1_njetsf();
	if (l1_njetst_branch != 0) l1_njetst();
	if (l1_nmus_branch != 0) l1_nmus();
	if (l1_emiso_ieta_branch != 0) l1_emiso_ieta();
	if (l1_emiso_iphi_branch != 0) l1_emiso_iphi();
	if (l1_emiso_rawId_branch != 0) l1_emiso_rawId();
	if (l1_emiso_type_branch != 0) l1_emiso_type();
	if (l1_emnoiso_ieta_branch != 0) l1_emnoiso_ieta();
	if (l1_emnoiso_iphi_branch != 0) l1_emnoiso_iphi();
	if (l1_emnoiso_rawId_branch != 0) l1_emnoiso_rawId();
	if (l1_emnoiso_type_branch != 0) l1_emnoiso_type();
	if (l1_jetsc_ieta_branch != 0) l1_jetsc_ieta();
	if (l1_jetsc_iphi_branch != 0) l1_jetsc_iphi();
	if (l1_jetsc_rawId_branch != 0) l1_jetsc_rawId();
	if (l1_jetsc_type_branch != 0) l1_jetsc_type();
	if (l1_jetsf_ieta_branch != 0) l1_jetsf_ieta();
	if (l1_jetsf_iphi_branch != 0) l1_jetsf_iphi();
	if (l1_jetsf_rawId_branch != 0) l1_jetsf_rawId();
	if (l1_jetsf_type_branch != 0) l1_jetsf_type();
	if (l1_jetst_ieta_branch != 0) l1_jetst_ieta();
	if (l1_jetst_iphi_branch != 0) l1_jetst_iphi();
	if (l1_jetst_rawId_branch != 0) l1_jetst_rawId();
	if (l1_jetst_type_branch != 0) l1_jetst_type();
	if (l1_mus_flags_branch != 0) l1_mus_flags();
	if (l1_mus_q_branch != 0) l1_mus_q();
	if (l1_mus_qual_branch != 0) l1_mus_qual();
	if (l1_mus_qualFlags_branch != 0) l1_mus_qualFlags();
	if (pxl_ndigis_pxb_branch != 0) pxl_ndigis_pxb();
	if (pxl_ndigis_pxf_branch != 0) pxl_ndigis_pxf();
	if (scs_elsidx_branch != 0) scs_elsidx();
	if (scs_severitySeed_branch != 0) scs_severitySeed();
	if (trks_algo_branch != 0) trks_algo();
	if (trks_charge_branch != 0) trks_charge();
	if (trks_exp_innerlayers_branch != 0) trks_exp_innerlayers();
	if (trks_exp_outerlayers_branch != 0) trks_exp_outerlayers();
	if (trks_layer1_det_branch != 0) trks_layer1_det();
	if (trks_layer1_layer_branch != 0) trks_layer1_layer();
	if (trks_layer1_sizerphi_branch != 0) trks_layer1_sizerphi();
	if (trks_layer1_sizerz_branch != 0) trks_layer1_sizerz();
	if (trks_lostHits_branch != 0) trks_lostHits();
	if (trks_lost_pixelhits_branch != 0) trks_lost_pixelhits();
	if (trks_nlayers_branch != 0) trks_nlayers();
	if (trks_nlayers3D_branch != 0) trks_nlayers3D();
	if (trks_nlayersLost_branch != 0) trks_nlayersLost();
	if (trks_qualityMask_branch != 0) trks_qualityMask();
	if (trks_validHits_branch != 0) trks_validHits();
	if (trks_valid_pixelhits_branch != 0) trks_valid_pixelhits();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (vtxs_tracksSize_branch != 0) vtxs_tracksSize();
	if (pseudo_srFlags_branch != 0) pseudo_srFlags();
	if (evt_ntwrs_branch != 0) evt_ntwrs();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (l1_bits1_branch != 0) l1_bits1();
	if (l1_bits2_branch != 0) l1_bits2();
	if (l1_bits3_branch != 0) l1_bits3();
	if (l1_bits4_branch != 0) l1_bits4();
	if (l1_techbits1_branch != 0) l1_techbits1();
	if (l1_techbits2_branch != 0) l1_techbits2();
	if (evt_nscs_branch != 0) evt_nscs();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (twrs_numBadEcalCells_branch != 0) twrs_numBadEcalCells();
	if (twrs_numBadHcalCells_branch != 0) twrs_numBadHcalCells();
	if (twrs_numProblematicEcalCells_branch != 0) twrs_numProblematicEcalCells();
	if (twrs_numProblematicHcalCells_branch != 0) twrs_numProblematicHcalCells();
	if (twrs_numRecoveredEcalCells_branch != 0) twrs_numRecoveredEcalCells();
	if (twrs_numRecoveredHcalCells_branch != 0) twrs_numRecoveredHcalCells();
}

	TString &evt_CMS2tag()
	{
		if (not evt_CMS2tag_isLoaded) {
			if (evt_CMS2tag_branch != 0) {
				evt_CMS2tag_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_CMS2tag_branch does not exist!\n");
				exit(1);
			}
			evt_CMS2tag_isLoaded = true;
		}
		return evt_CMS2tag_;
	}
	TString &evt_dataset()
	{
		if (not evt_dataset_isLoaded) {
			if (evt_dataset_branch != 0) {
				evt_dataset_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_dataset_branch does not exist!\n");
				exit(1);
			}
			evt_dataset_isLoaded = true;
		}
		return evt_dataset_;
	}
	vector<TString> &l1_trigNames()
	{
		if (not l1_trigNames_isLoaded) {
			if (l1_trigNames_branch != 0) {
				l1_trigNames_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_trigNames_branch does not exist!\n");
				exit(1);
			}
			l1_trigNames_isLoaded = true;
		}
		return l1_trigNames_;
	}
	float &evt_bField()
	{
		if (not evt_bField_isLoaded) {
			if (evt_bField_branch != 0) {
				evt_bField_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_bField_)) {
					printf("branch evt_bField_branch contains a bad float: %f\n", evt_bField_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bField_branch does not exist!\n");
				exit(1);
			}
			evt_bField_isLoaded = true;
		}
		return evt_bField_;
	}
	float &hcalnoise_eventChargeFraction()
	{
		if (not hcalnoise_eventChargeFraction_isLoaded) {
			if (hcalnoise_eventChargeFraction_branch != 0) {
				hcalnoise_eventChargeFraction_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_eventChargeFraction_)) {
					printf("branch hcalnoise_eventChargeFraction_branch contains a bad float: %f\n", hcalnoise_eventChargeFraction_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_eventChargeFraction_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_eventChargeFraction_isLoaded = true;
		}
		return hcalnoise_eventChargeFraction_;
	}
	float &hcalnoise_eventEMEnergy()
	{
		if (not hcalnoise_eventEMEnergy_isLoaded) {
			if (hcalnoise_eventEMEnergy_branch != 0) {
				hcalnoise_eventEMEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_eventEMEnergy_)) {
					printf("branch hcalnoise_eventEMEnergy_branch contains a bad float: %f\n", hcalnoise_eventEMEnergy_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_eventEMEnergy_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_eventEMEnergy_isLoaded = true;
		}
		return hcalnoise_eventEMEnergy_;
	}
	float &hcalnoise_eventEMFraction()
	{
		if (not hcalnoise_eventEMFraction_isLoaded) {
			if (hcalnoise_eventEMFraction_branch != 0) {
				hcalnoise_eventEMFraction_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_eventEMFraction_)) {
					printf("branch hcalnoise_eventEMFraction_branch contains a bad float: %f\n", hcalnoise_eventEMFraction_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_eventEMFraction_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_eventEMFraction_isLoaded = true;
		}
		return hcalnoise_eventEMFraction_;
	}
	float &hcalnoise_eventHadEnergy()
	{
		if (not hcalnoise_eventHadEnergy_isLoaded) {
			if (hcalnoise_eventHadEnergy_branch != 0) {
				hcalnoise_eventHadEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_eventHadEnergy_)) {
					printf("branch hcalnoise_eventHadEnergy_branch contains a bad float: %f\n", hcalnoise_eventHadEnergy_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_eventHadEnergy_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_eventHadEnergy_isLoaded = true;
		}
		return hcalnoise_eventHadEnergy_;
	}
	float &hcalnoise_eventTrackEnergy()
	{
		if (not hcalnoise_eventTrackEnergy_isLoaded) {
			if (hcalnoise_eventTrackEnergy_branch != 0) {
				hcalnoise_eventTrackEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_eventTrackEnergy_)) {
					printf("branch hcalnoise_eventTrackEnergy_branch contains a bad float: %f\n", hcalnoise_eventTrackEnergy_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_eventTrackEnergy_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_eventTrackEnergy_isLoaded = true;
		}
		return hcalnoise_eventTrackEnergy_;
	}
	float &hcalnoise_max10GeVHitTime()
	{
		if (not hcalnoise_max10GeVHitTime_isLoaded) {
			if (hcalnoise_max10GeVHitTime_branch != 0) {
				hcalnoise_max10GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_max10GeVHitTime_)) {
					printf("branch hcalnoise_max10GeVHitTime_branch contains a bad float: %f\n", hcalnoise_max10GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_max10GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_max10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_max10GeVHitTime_;
	}
	float &hcalnoise_max25GeVHitTime()
	{
		if (not hcalnoise_max25GeVHitTime_isLoaded) {
			if (hcalnoise_max25GeVHitTime_branch != 0) {
				hcalnoise_max25GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_max25GeVHitTime_)) {
					printf("branch hcalnoise_max25GeVHitTime_branch contains a bad float: %f\n", hcalnoise_max25GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_max25GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_max25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_max25GeVHitTime_;
	}
	float &hcalnoise_min10GeVHitTime()
	{
		if (not hcalnoise_min10GeVHitTime_isLoaded) {
			if (hcalnoise_min10GeVHitTime_branch != 0) {
				hcalnoise_min10GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_min10GeVHitTime_)) {
					printf("branch hcalnoise_min10GeVHitTime_branch contains a bad float: %f\n", hcalnoise_min10GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_min10GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_min10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_min10GeVHitTime_;
	}
	float &hcalnoise_min25GeVHitTime()
	{
		if (not hcalnoise_min25GeVHitTime_isLoaded) {
			if (hcalnoise_min25GeVHitTime_branch != 0) {
				hcalnoise_min25GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_min25GeVHitTime_)) {
					printf("branch hcalnoise_min25GeVHitTime_branch contains a bad float: %f\n", hcalnoise_min25GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_min25GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_min25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_min25GeVHitTime_;
	}
	float &hcalnoise_minE10TS()
	{
		if (not hcalnoise_minE10TS_isLoaded) {
			if (hcalnoise_minE10TS_branch != 0) {
				hcalnoise_minE10TS_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_minE10TS_)) {
					printf("branch hcalnoise_minE10TS_branch contains a bad float: %f\n", hcalnoise_minE10TS_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_minE10TS_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_minE10TS_isLoaded = true;
		}
		return hcalnoise_minE10TS_;
	}
	float &hcalnoise_minE2Over10TS()
	{
		if (not hcalnoise_minE2Over10TS_isLoaded) {
			if (hcalnoise_minE2Over10TS_branch != 0) {
				hcalnoise_minE2Over10TS_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_minE2Over10TS_)) {
					printf("branch hcalnoise_minE2Over10TS_branch contains a bad float: %f\n", hcalnoise_minE2Over10TS_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_minE2Over10TS_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_minE2Over10TS_isLoaded = true;
		}
		return hcalnoise_minE2Over10TS_;
	}
	float &hcalnoise_minE2TS()
	{
		if (not hcalnoise_minE2TS_isLoaded) {
			if (hcalnoise_minE2TS_branch != 0) {
				hcalnoise_minE2TS_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_minE2TS_)) {
					printf("branch hcalnoise_minE2TS_branch contains a bad float: %f\n", hcalnoise_minE2TS_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_minE2TS_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_minE2TS_isLoaded = true;
		}
		return hcalnoise_minE2TS_;
	}
	float &hcalnoise_minHPDEMF()
	{
		if (not hcalnoise_minHPDEMF_isLoaded) {
			if (hcalnoise_minHPDEMF_branch != 0) {
				hcalnoise_minHPDEMF_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_minHPDEMF_)) {
					printf("branch hcalnoise_minHPDEMF_branch contains a bad float: %f\n", hcalnoise_minHPDEMF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_minHPDEMF_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_minHPDEMF_isLoaded = true;
		}
		return hcalnoise_minHPDEMF_;
	}
	float &hcalnoise_minRBXEMF()
	{
		if (not hcalnoise_minRBXEMF_isLoaded) {
			if (hcalnoise_minRBXEMF_branch != 0) {
				hcalnoise_minRBXEMF_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_minRBXEMF_)) {
					printf("branch hcalnoise_minRBXEMF_branch contains a bad float: %f\n", hcalnoise_minRBXEMF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_minRBXEMF_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_minRBXEMF_isLoaded = true;
		}
		return hcalnoise_minRBXEMF_;
	}
	float &hcalnoise_rms10GeVHitTime()
	{
		if (not hcalnoise_rms10GeVHitTime_isLoaded) {
			if (hcalnoise_rms10GeVHitTime_branch != 0) {
				hcalnoise_rms10GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_rms10GeVHitTime_)) {
					printf("branch hcalnoise_rms10GeVHitTime_branch contains a bad float: %f\n", hcalnoise_rms10GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_rms10GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_rms10GeVHitTime_isLoaded = true;
		}
		return hcalnoise_rms10GeVHitTime_;
	}
	float &hcalnoise_rms25GeVHitTime()
	{
		if (not hcalnoise_rms25GeVHitTime_isLoaded) {
			if (hcalnoise_rms25GeVHitTime_branch != 0) {
				hcalnoise_rms25GeVHitTime_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(hcalnoise_rms25GeVHitTime_)) {
					printf("branch hcalnoise_rms25GeVHitTime_branch contains a bad float: %f\n", hcalnoise_rms25GeVHitTime_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_rms25GeVHitTime_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_rms25GeVHitTime_isLoaded = true;
		}
		return hcalnoise_rms25GeVHitTime_;
	}
	float &l1_met_etTot()
	{
		if (not l1_met_etTot_isLoaded) {
			if (l1_met_etTot_branch != 0) {
				l1_met_etTot_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1_met_etTot_)) {
					printf("branch l1_met_etTot_branch contains a bad float: %f\n", l1_met_etTot_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_met_etTot_branch does not exist!\n");
				exit(1);
			}
			l1_met_etTot_isLoaded = true;
		}
		return l1_met_etTot_;
	}
	float &l1_met_met()
	{
		if (not l1_met_met_isLoaded) {
			if (l1_met_met_branch != 0) {
				l1_met_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1_met_met_)) {
					printf("branch l1_met_met_branch contains a bad float: %f\n", l1_met_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_met_met_branch does not exist!\n");
				exit(1);
			}
			l1_met_met_isLoaded = true;
		}
		return l1_met_met_;
	}
	float &l1_mht_htTot()
	{
		if (not l1_mht_htTot_isLoaded) {
			if (l1_mht_htTot_branch != 0) {
				l1_mht_htTot_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1_mht_htTot_)) {
					printf("branch l1_mht_htTot_branch contains a bad float: %f\n", l1_mht_htTot_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mht_htTot_branch does not exist!\n");
				exit(1);
			}
			l1_mht_htTot_isLoaded = true;
		}
		return l1_mht_htTot_;
	}
	float &l1_mht_mht()
	{
		if (not l1_mht_mht_isLoaded) {
			if (l1_mht_mht_branch != 0) {
				l1_mht_mht_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1_mht_mht_)) {
					printf("branch l1_mht_mht_branch contains a bad float: %f\n", l1_mht_mht_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mht_mht_branch does not exist!\n");
				exit(1);
			}
			l1_mht_mht_isLoaded = true;
		}
		return l1_mht_mht_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_met_p4()
	{
		if (not l1_met_p4_isLoaded) {
			if (l1_met_p4_branch != 0) {
				l1_met_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				int e;
				frexp(l1_met_p4_.pt(), &e);
				if (not isfinite(l1_met_p4_.pt()) || e > 30) {
					printf("branch l1_met_p4_branch contains a bad float: %f\n", l1_met_p4_.pt());
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_met_p4_branch does not exist!\n");
				exit(1);
			}
			l1_met_p4_isLoaded = true;
		}
		return l1_met_p4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_mht_p4()
	{
		if (not l1_mht_p4_isLoaded) {
			if (l1_mht_p4_branch != 0) {
				l1_mht_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				int e;
				frexp(l1_mht_p4_.pt(), &e);
				if (not isfinite(l1_mht_p4_.pt()) || e > 30) {
					printf("branch l1_mht_p4_branch contains a bad float: %f\n", l1_mht_p4_.pt());
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mht_p4_branch does not exist!\n");
				exit(1);
			}
			l1_mht_p4_isLoaded = true;
		}
		return l1_mht_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emiso_p4()
	{
		if (not l1_emiso_p4_isLoaded) {
			if (l1_emiso_p4_branch != 0) {
				l1_emiso_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_emiso_p4_.begin(); i != l1_emiso_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_emiso_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emiso_p4_branch does not exist!\n");
				exit(1);
			}
			l1_emiso_p4_isLoaded = true;
		}
		return l1_emiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emnoiso_p4()
	{
		if (not l1_emnoiso_p4_isLoaded) {
			if (l1_emnoiso_p4_branch != 0) {
				l1_emnoiso_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_emnoiso_p4_.begin(); i != l1_emnoiso_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_emnoiso_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emnoiso_p4_branch does not exist!\n");
				exit(1);
			}
			l1_emnoiso_p4_isLoaded = true;
		}
		return l1_emnoiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsc_p4()
	{
		if (not l1_jetsc_p4_isLoaded) {
			if (l1_jetsc_p4_branch != 0) {
				l1_jetsc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_jetsc_p4_.begin(); i != l1_jetsc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_jetsc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsc_p4_branch does not exist!\n");
				exit(1);
			}
			l1_jetsc_p4_isLoaded = true;
		}
		return l1_jetsc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsf_p4()
	{
		if (not l1_jetsf_p4_isLoaded) {
			if (l1_jetsf_p4_branch != 0) {
				l1_jetsf_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_jetsf_p4_.begin(); i != l1_jetsf_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_jetsf_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsf_p4_branch does not exist!\n");
				exit(1);
			}
			l1_jetsf_p4_isLoaded = true;
		}
		return l1_jetsf_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetst_p4()
	{
		if (not l1_jetst_p4_isLoaded) {
			if (l1_jetst_p4_branch != 0) {
				l1_jetst_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_jetst_p4_.begin(); i != l1_jetst_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_jetst_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetst_p4_branch does not exist!\n");
				exit(1);
			}
			l1_jetst_p4_isLoaded = true;
		}
		return l1_jetst_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_mus_p4()
	{
		if (not l1_mus_p4_isLoaded) {
			if (l1_mus_p4_branch != 0) {
				l1_mus_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = l1_mus_p4_.begin(); i != l1_mus_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1_mus_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mus_p4_branch does not exist!\n");
				exit(1);
			}
			l1_mus_p4_isLoaded = true;
		}
		return l1_mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_p4()
	{
		if (not scs_p4_isLoaded) {
			if (scs_p4_branch != 0) {
				scs_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = scs_p4_.begin(); i != scs_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch scs_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_p4_branch does not exist!\n");
				exit(1);
			}
			scs_p4_isLoaded = true;
		}
		return scs_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_pos_p4()
	{
		if (not scs_pos_p4_isLoaded) {
			if (scs_pos_p4_branch != 0) {
				scs_pos_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = scs_pos_p4_.begin(); i != scs_pos_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch scs_pos_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_pos_p4_branch does not exist!\n");
				exit(1);
			}
			scs_pos_p4_isLoaded = true;
		}
		return scs_pos_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_vtx_p4()
	{
		if (not scs_vtx_p4_isLoaded) {
			if (scs_vtx_p4_branch != 0) {
				scs_vtx_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = scs_vtx_p4_.begin(); i != scs_vtx_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch scs_vtx_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_vtx_p4_branch does not exist!\n");
				exit(1);
			}
			scs_vtx_p4_isLoaded = true;
		}
		return scs_vtx_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_inner_position()
	{
		if (not trks_inner_position_isLoaded) {
			if (trks_inner_position_branch != 0) {
				trks_inner_position_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = trks_inner_position_.begin(); i != trks_inner_position_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trks_inner_position_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_inner_position_branch does not exist!\n");
				exit(1);
			}
			trks_inner_position_isLoaded = true;
		}
		return trks_inner_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_p4()
	{
		if (not trks_outer_p4_isLoaded) {
			if (trks_outer_p4_branch != 0) {
				trks_outer_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = trks_outer_p4_.begin(); i != trks_outer_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trks_outer_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_outer_p4_branch does not exist!\n");
				exit(1);
			}
			trks_outer_p4_isLoaded = true;
		}
		return trks_outer_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_position()
	{
		if (not trks_outer_position_isLoaded) {
			if (trks_outer_position_branch != 0) {
				trks_outer_position_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = trks_outer_position_.begin(); i != trks_outer_position_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trks_outer_position_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_outer_position_branch does not exist!\n");
				exit(1);
			}
			trks_outer_position_isLoaded = true;
		}
		return trks_outer_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_trk_p4()
	{
		if (not trks_trk_p4_isLoaded) {
			if (trks_trk_p4_branch != 0) {
				trks_trk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = trks_trk_p4_.begin(); i != trks_trk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trks_trk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_trk_p4_branch does not exist!\n");
				exit(1);
			}
			trks_trk_p4_isLoaded = true;
		}
		return trks_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_vertex_p4()
	{
		if (not trks_vertex_p4_isLoaded) {
			if (trks_vertex_p4_branch != 0) {
				trks_vertex_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = trks_vertex_p4_.begin(); i != trks_vertex_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trks_vertex_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			trks_vertex_p4_isLoaded = true;
		}
		return trks_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position()
	{
		if (not vtxs_position_isLoaded) {
			if (vtxs_position_branch != 0) {
				vtxs_position_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > >::const_iterator i = vtxs_position_.begin(); i != vtxs_position_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch vtxs_position_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_position_branch does not exist!\n");
				exit(1);
			}
			vtxs_position_isLoaded = true;
		}
		return vtxs_position_;
	}
	vector<float> &twrs_ecalTime()
	{
		if (not twrs_ecalTime_isLoaded) {
			if (twrs_ecalTime_branch != 0) {
				twrs_ecalTime_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_ecalTime_.begin(); i != twrs_ecalTime_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_ecalTime_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_ecalTime_branch does not exist!\n");
				exit(1);
			}
			twrs_ecalTime_isLoaded = true;
		}
		return twrs_ecalTime_;
	}
	vector<float> &twrs_emEnergy()
	{
		if (not twrs_emEnergy_isLoaded) {
			if (twrs_emEnergy_branch != 0) {
				twrs_emEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_emEnergy_.begin(); i != twrs_emEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_emEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_emEnergy_branch does not exist!\n");
				exit(1);
			}
			twrs_emEnergy_isLoaded = true;
		}
		return twrs_emEnergy_;
	}
	vector<float> &twrs_emEt()
	{
		if (not twrs_emEt_isLoaded) {
			if (twrs_emEt_branch != 0) {
				twrs_emEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_emEt_.begin(); i != twrs_emEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_emEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_emEt_branch does not exist!\n");
				exit(1);
			}
			twrs_emEt_isLoaded = true;
		}
		return twrs_emEt_;
	}
	vector<float> &twrs_emEtcorr()
	{
		if (not twrs_emEtcorr_isLoaded) {
			if (twrs_emEtcorr_branch != 0) {
				twrs_emEtcorr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_emEtcorr_.begin(); i != twrs_emEtcorr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_emEtcorr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_emEtcorr_branch does not exist!\n");
				exit(1);
			}
			twrs_emEtcorr_isLoaded = true;
		}
		return twrs_emEtcorr_;
	}
	vector<float> &twrs_eta()
	{
		if (not twrs_eta_isLoaded) {
			if (twrs_eta_branch != 0) {
				twrs_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_eta_.begin(); i != twrs_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_eta_branch does not exist!\n");
				exit(1);
			}
			twrs_eta_isLoaded = true;
		}
		return twrs_eta_;
	}
	vector<float> &twrs_etcorr()
	{
		if (not twrs_etcorr_isLoaded) {
			if (twrs_etcorr_branch != 0) {
				twrs_etcorr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_etcorr_.begin(); i != twrs_etcorr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_etcorr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_etcorr_branch does not exist!\n");
				exit(1);
			}
			twrs_etcorr_isLoaded = true;
		}
		return twrs_etcorr_;
	}
	vector<float> &twrs_hadEnergy()
	{
		if (not twrs_hadEnergy_isLoaded) {
			if (twrs_hadEnergy_branch != 0) {
				twrs_hadEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_hadEnergy_.begin(); i != twrs_hadEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_hadEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_hadEnergy_branch does not exist!\n");
				exit(1);
			}
			twrs_hadEnergy_isLoaded = true;
		}
		return twrs_hadEnergy_;
	}
	vector<float> &twrs_hadEt()
	{
		if (not twrs_hadEt_isLoaded) {
			if (twrs_hadEt_branch != 0) {
				twrs_hadEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_hadEt_.begin(); i != twrs_hadEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_hadEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_hadEt_branch does not exist!\n");
				exit(1);
			}
			twrs_hadEt_isLoaded = true;
		}
		return twrs_hadEt_;
	}
	vector<float> &twrs_hadEtcorr()
	{
		if (not twrs_hadEtcorr_isLoaded) {
			if (twrs_hadEtcorr_branch != 0) {
				twrs_hadEtcorr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_hadEtcorr_.begin(); i != twrs_hadEtcorr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_hadEtcorr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_hadEtcorr_branch does not exist!\n");
				exit(1);
			}
			twrs_hadEtcorr_isLoaded = true;
		}
		return twrs_hadEtcorr_;
	}
	vector<float> &twrs_hcalTime()
	{
		if (not twrs_hcalTime_isLoaded) {
			if (twrs_hcalTime_branch != 0) {
				twrs_hcalTime_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_hcalTime_.begin(); i != twrs_hcalTime_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_hcalTime_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_hcalTime_branch does not exist!\n");
				exit(1);
			}
			twrs_hcalTime_isLoaded = true;
		}
		return twrs_hcalTime_;
	}
	vector<float> &twrs_outerEnergy()
	{
		if (not twrs_outerEnergy_isLoaded) {
			if (twrs_outerEnergy_branch != 0) {
				twrs_outerEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_outerEnergy_.begin(); i != twrs_outerEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_outerEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_outerEnergy_branch does not exist!\n");
				exit(1);
			}
			twrs_outerEnergy_isLoaded = true;
		}
		return twrs_outerEnergy_;
	}
	vector<float> &twrs_outerEt()
	{
		if (not twrs_outerEt_isLoaded) {
			if (twrs_outerEt_branch != 0) {
				twrs_outerEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_outerEt_.begin(); i != twrs_outerEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_outerEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_outerEt_branch does not exist!\n");
				exit(1);
			}
			twrs_outerEt_isLoaded = true;
		}
		return twrs_outerEt_;
	}
	vector<float> &twrs_outerEtcorr()
	{
		if (not twrs_outerEtcorr_isLoaded) {
			if (twrs_outerEtcorr_branch != 0) {
				twrs_outerEtcorr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_outerEtcorr_.begin(); i != twrs_outerEtcorr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_outerEtcorr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_outerEtcorr_branch does not exist!\n");
				exit(1);
			}
			twrs_outerEtcorr_isLoaded = true;
		}
		return twrs_outerEtcorr_;
	}
	vector<float> &twrs_pcorr()
	{
		if (not twrs_pcorr_isLoaded) {
			if (twrs_pcorr_branch != 0) {
				twrs_pcorr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_pcorr_.begin(); i != twrs_pcorr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_pcorr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_pcorr_branch does not exist!\n");
				exit(1);
			}
			twrs_pcorr_isLoaded = true;
		}
		return twrs_pcorr_;
	}
	vector<float> &twrs_phi()
	{
		if (not twrs_phi_isLoaded) {
			if (twrs_phi_branch != 0) {
				twrs_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = twrs_phi_.begin(); i != twrs_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch twrs_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_phi_branch does not exist!\n");
				exit(1);
			}
			twrs_phi_isLoaded = true;
		}
		return twrs_phi_;
	}
	vector<float> &pseudo_dRClosestTower()
	{
		if (not pseudo_dRClosestTower_isLoaded) {
			if (pseudo_dRClosestTower_branch != 0) {
				pseudo_dRClosestTower_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_dRClosestTower_.begin(); i != pseudo_dRClosestTower_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_dRClosestTower_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_dRClosestTower_branch does not exist!\n");
				exit(1);
			}
			pseudo_dRClosestTower_isLoaded = true;
		}
		return pseudo_dRClosestTower_;
	}
	vector<float> &pseudo_dRClosestTowerEmEt()
	{
		if (not pseudo_dRClosestTowerEmEt_isLoaded) {
			if (pseudo_dRClosestTowerEmEt_branch != 0) {
				pseudo_dRClosestTowerEmEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_dRClosestTowerEmEt_.begin(); i != pseudo_dRClosestTowerEmEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_dRClosestTowerEmEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_dRClosestTowerEmEt_branch does not exist!\n");
				exit(1);
			}
			pseudo_dRClosestTowerEmEt_isLoaded = true;
		}
		return pseudo_dRClosestTowerEmEt_;
	}
	vector<float> &pseudo_ecalEta()
	{
		if (not pseudo_ecalEta_isLoaded) {
			if (pseudo_ecalEta_branch != 0) {
				pseudo_ecalEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_ecalEta_.begin(); i != pseudo_ecalEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_ecalEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_ecalEta_branch does not exist!\n");
				exit(1);
			}
			pseudo_ecalEta_isLoaded = true;
		}
		return pseudo_ecalEta_;
	}
	vector<float> &pseudo_ecalIso03()
	{
		if (not pseudo_ecalIso03_isLoaded) {
			if (pseudo_ecalIso03_branch != 0) {
				pseudo_ecalIso03_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_ecalIso03_.begin(); i != pseudo_ecalIso03_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_ecalIso03_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_ecalIso03_branch does not exist!\n");
				exit(1);
			}
			pseudo_ecalIso03_isLoaded = true;
		}
		return pseudo_ecalIso03_;
	}
	vector<float> &pseudo_ecalPhi()
	{
		if (not pseudo_ecalPhi_isLoaded) {
			if (pseudo_ecalPhi_branch != 0) {
				pseudo_ecalPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_ecalPhi_.begin(); i != pseudo_ecalPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_ecalPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_ecalPhi_branch does not exist!\n");
				exit(1);
			}
			pseudo_ecalPhi_isLoaded = true;
		}
		return pseudo_ecalPhi_;
	}
	vector<float> &pseudo_eta()
	{
		if (not pseudo_eta_isLoaded) {
			if (pseudo_eta_branch != 0) {
				pseudo_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_eta_.begin(); i != pseudo_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_eta_branch does not exist!\n");
				exit(1);
			}
			pseudo_eta_isLoaded = true;
		}
		return pseudo_eta_;
	}
	vector<float> &pseudo_hcalD1Iso03()
	{
		if (not pseudo_hcalD1Iso03_isLoaded) {
			if (pseudo_hcalD1Iso03_branch != 0) {
				pseudo_hcalD1Iso03_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_hcalD1Iso03_.begin(); i != pseudo_hcalD1Iso03_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_hcalD1Iso03_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_hcalD1Iso03_branch does not exist!\n");
				exit(1);
			}
			pseudo_hcalD1Iso03_isLoaded = true;
		}
		return pseudo_hcalD1Iso03_;
	}
	vector<float> &pseudo_hcalD2Iso03()
	{
		if (not pseudo_hcalD2Iso03_isLoaded) {
			if (pseudo_hcalD2Iso03_branch != 0) {
				pseudo_hcalD2Iso03_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_hcalD2Iso03_.begin(); i != pseudo_hcalD2Iso03_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_hcalD2Iso03_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_hcalD2Iso03_branch does not exist!\n");
				exit(1);
			}
			pseudo_hcalD2Iso03_isLoaded = true;
		}
		return pseudo_hcalD2Iso03_;
	}
	vector<float> &pseudo_phi()
	{
		if (not pseudo_phi_isLoaded) {
			if (pseudo_phi_branch != 0) {
				pseudo_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_phi_.begin(); i != pseudo_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_phi_branch does not exist!\n");
				exit(1);
			}
			pseudo_phi_isLoaded = true;
		}
		return pseudo_phi_;
	}
	vector<float> &pseudo_tkIso03()
	{
		if (not pseudo_tkIso03_isLoaded) {
			if (pseudo_tkIso03_branch != 0) {
				pseudo_tkIso03_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_tkIso03_.begin(); i != pseudo_tkIso03_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_tkIso03_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_tkIso03_branch does not exist!\n");
				exit(1);
			}
			pseudo_tkIso03_isLoaded = true;
		}
		return pseudo_tkIso03_;
	}
	vector<float> &pseudo_towerEmEt()
	{
		if (not pseudo_towerEmEt_isLoaded) {
			if (pseudo_towerEmEt_branch != 0) {
				pseudo_towerEmEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_towerEmEt_.begin(); i != pseudo_towerEmEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_towerEmEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_towerEmEt_branch does not exist!\n");
				exit(1);
			}
			pseudo_towerEmEt_isLoaded = true;
		}
		return pseudo_towerEmEt_;
	}
	vector<float> &pseudo_towerHadEt()
	{
		if (not pseudo_towerHadEt_isLoaded) {
			if (pseudo_towerHadEt_branch != 0) {
				pseudo_towerHadEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = pseudo_towerHadEt_.begin(); i != pseudo_towerHadEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch pseudo_towerHadEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_towerHadEt_branch does not exist!\n");
				exit(1);
			}
			pseudo_towerHadEt_isLoaded = true;
		}
		return pseudo_towerHadEt_;
	}
	vector<float> &scs_clustersSize()
	{
		if (not scs_clustersSize_isLoaded) {
			if (scs_clustersSize_branch != 0) {
				scs_clustersSize_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_clustersSize_.begin(); i != scs_clustersSize_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_clustersSize_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_clustersSize_branch does not exist!\n");
				exit(1);
			}
			scs_clustersSize_isLoaded = true;
		}
		return scs_clustersSize_;
	}
	vector<float> &scs_crystalsSize()
	{
		if (not scs_crystalsSize_isLoaded) {
			if (scs_crystalsSize_branch != 0) {
				scs_crystalsSize_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_crystalsSize_.begin(); i != scs_crystalsSize_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_crystalsSize_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_crystalsSize_branch does not exist!\n");
				exit(1);
			}
			scs_crystalsSize_isLoaded = true;
		}
		return scs_crystalsSize_;
	}
	vector<float> &scs_e1x3()
	{
		if (not scs_e1x3_isLoaded) {
			if (scs_e1x3_branch != 0) {
				scs_e1x3_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e1x3_.begin(); i != scs_e1x3_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e1x3_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e1x3_branch does not exist!\n");
				exit(1);
			}
			scs_e1x3_isLoaded = true;
		}
		return scs_e1x3_;
	}
	vector<float> &scs_e1x5()
	{
		if (not scs_e1x5_isLoaded) {
			if (scs_e1x5_branch != 0) {
				scs_e1x5_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e1x5_.begin(); i != scs_e1x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e1x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e1x5_branch does not exist!\n");
				exit(1);
			}
			scs_e1x5_isLoaded = true;
		}
		return scs_e1x5_;
	}
	vector<float> &scs_e2nd()
	{
		if (not scs_e2nd_isLoaded) {
			if (scs_e2nd_branch != 0) {
				scs_e2nd_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e2nd_.begin(); i != scs_e2nd_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e2nd_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e2nd_branch does not exist!\n");
				exit(1);
			}
			scs_e2nd_isLoaded = true;
		}
		return scs_e2nd_;
	}
	vector<float> &scs_e2x2()
	{
		if (not scs_e2x2_isLoaded) {
			if (scs_e2x2_branch != 0) {
				scs_e2x2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e2x2_.begin(); i != scs_e2x2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e2x2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e2x2_branch does not exist!\n");
				exit(1);
			}
			scs_e2x2_isLoaded = true;
		}
		return scs_e2x2_;
	}
	vector<float> &scs_e2x5Max()
	{
		if (not scs_e2x5Max_isLoaded) {
			if (scs_e2x5Max_branch != 0) {
				scs_e2x5Max_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e2x5Max_.begin(); i != scs_e2x5Max_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e2x5Max_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			scs_e2x5Max_isLoaded = true;
		}
		return scs_e2x5Max_;
	}
	vector<float> &scs_e3x1()
	{
		if (not scs_e3x1_isLoaded) {
			if (scs_e3x1_branch != 0) {
				scs_e3x1_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e3x1_.begin(); i != scs_e3x1_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e3x1_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e3x1_branch does not exist!\n");
				exit(1);
			}
			scs_e3x1_isLoaded = true;
		}
		return scs_e3x1_;
	}
	vector<float> &scs_e3x2()
	{
		if (not scs_e3x2_isLoaded) {
			if (scs_e3x2_branch != 0) {
				scs_e3x2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e3x2_.begin(); i != scs_e3x2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e3x2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e3x2_branch does not exist!\n");
				exit(1);
			}
			scs_e3x2_isLoaded = true;
		}
		return scs_e3x2_;
	}
	vector<float> &scs_e3x3()
	{
		if (not scs_e3x3_isLoaded) {
			if (scs_e3x3_branch != 0) {
				scs_e3x3_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e3x3_.begin(); i != scs_e3x3_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e3x3_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e3x3_branch does not exist!\n");
				exit(1);
			}
			scs_e3x3_isLoaded = true;
		}
		return scs_e3x3_;
	}
	vector<float> &scs_e4x4()
	{
		if (not scs_e4x4_isLoaded) {
			if (scs_e4x4_branch != 0) {
				scs_e4x4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e4x4_.begin(); i != scs_e4x4_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e4x4_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e4x4_branch does not exist!\n");
				exit(1);
			}
			scs_e4x4_isLoaded = true;
		}
		return scs_e4x4_;
	}
	vector<float> &scs_e5x5()
	{
		if (not scs_e5x5_isLoaded) {
			if (scs_e5x5_branch != 0) {
				scs_e5x5_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_e5x5_.begin(); i != scs_e5x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_e5x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_e5x5_branch does not exist!\n");
				exit(1);
			}
			scs_e5x5_isLoaded = true;
		}
		return scs_e5x5_;
	}
	vector<float> &scs_eMax()
	{
		if (not scs_eMax_isLoaded) {
			if (scs_eMax_branch != 0) {
				scs_eMax_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_eMax_.begin(); i != scs_eMax_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_eMax_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_eMax_branch does not exist!\n");
				exit(1);
			}
			scs_eMax_isLoaded = true;
		}
		return scs_eMax_;
	}
	vector<float> &scs_eSeed()
	{
		if (not scs_eSeed_isLoaded) {
			if (scs_eSeed_branch != 0) {
				scs_eSeed_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_eSeed_.begin(); i != scs_eSeed_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_eSeed_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_eSeed_branch does not exist!\n");
				exit(1);
			}
			scs_eSeed_isLoaded = true;
		}
		return scs_eSeed_;
	}
	vector<float> &scs_energy()
	{
		if (not scs_energy_isLoaded) {
			if (scs_energy_branch != 0) {
				scs_energy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_energy_.begin(); i != scs_energy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_energy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_energy_branch does not exist!\n");
				exit(1);
			}
			scs_energy_isLoaded = true;
		}
		return scs_energy_;
	}
	vector<float> &scs_eta()
	{
		if (not scs_eta_isLoaded) {
			if (scs_eta_branch != 0) {
				scs_eta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_eta_.begin(); i != scs_eta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_eta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_eta_branch does not exist!\n");
				exit(1);
			}
			scs_eta_isLoaded = true;
		}
		return scs_eta_;
	}
	vector<float> &scs_hoe()
	{
		if (not scs_hoe_isLoaded) {
			if (scs_hoe_branch != 0) {
				scs_hoe_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_hoe_.begin(); i != scs_hoe_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_hoe_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_hoe_branch does not exist!\n");
				exit(1);
			}
			scs_hoe_isLoaded = true;
		}
		return scs_hoe_;
	}
	vector<float> &scs_phi()
	{
		if (not scs_phi_isLoaded) {
			if (scs_phi_branch != 0) {
				scs_phi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_phi_.begin(); i != scs_phi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_phi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_phi_branch does not exist!\n");
				exit(1);
			}
			scs_phi_isLoaded = true;
		}
		return scs_phi_;
	}
	vector<float> &scs_preshowerEnergy()
	{
		if (not scs_preshowerEnergy_isLoaded) {
			if (scs_preshowerEnergy_branch != 0) {
				scs_preshowerEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_preshowerEnergy_.begin(); i != scs_preshowerEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_preshowerEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_preshowerEnergy_branch does not exist!\n");
				exit(1);
			}
			scs_preshowerEnergy_isLoaded = true;
		}
		return scs_preshowerEnergy_;
	}
	vector<float> &scs_rawEnergy()
	{
		if (not scs_rawEnergy_isLoaded) {
			if (scs_rawEnergy_branch != 0) {
				scs_rawEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_rawEnergy_.begin(); i != scs_rawEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_rawEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_rawEnergy_branch does not exist!\n");
				exit(1);
			}
			scs_rawEnergy_isLoaded = true;
		}
		return scs_rawEnergy_;
	}
	vector<float> &scs_sigmaEtaEta()
	{
		if (not scs_sigmaEtaEta_isLoaded) {
			if (scs_sigmaEtaEta_branch != 0) {
				scs_sigmaEtaEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaEtaEta_.begin(); i != scs_sigmaEtaEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaEtaEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaEtaEta_isLoaded = true;
		}
		return scs_sigmaEtaEta_;
	}
	vector<float> &scs_sigmaEtaPhi()
	{
		if (not scs_sigmaEtaPhi_isLoaded) {
			if (scs_sigmaEtaPhi_branch != 0) {
				scs_sigmaEtaPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaEtaPhi_.begin(); i != scs_sigmaEtaPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaEtaPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaEtaPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaEtaPhi_isLoaded = true;
		}
		return scs_sigmaEtaPhi_;
	}
	vector<float> &scs_sigmaIEtaIEta()
	{
		if (not scs_sigmaIEtaIEta_isLoaded) {
			if (scs_sigmaIEtaIEta_branch != 0) {
				scs_sigmaIEtaIEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaIEtaIEta_.begin(); i != scs_sigmaIEtaIEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaIEtaIEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIEta_isLoaded = true;
		}
		return scs_sigmaIEtaIEta_;
	}
	vector<float> &scs_sigmaIEtaIPhi()
	{
		if (not scs_sigmaIEtaIPhi_isLoaded) {
			if (scs_sigmaIEtaIPhi_branch != 0) {
				scs_sigmaIEtaIPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaIEtaIPhi_.begin(); i != scs_sigmaIEtaIPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaIEtaIPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaIEtaIPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIPhi_isLoaded = true;
		}
		return scs_sigmaIEtaIPhi_;
	}
	vector<float> &scs_sigmaIPhiIPhi()
	{
		if (not scs_sigmaIPhiIPhi_isLoaded) {
			if (scs_sigmaIPhiIPhi_branch != 0) {
				scs_sigmaIPhiIPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaIPhiIPhi_.begin(); i != scs_sigmaIPhiIPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaIPhiIPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaIPhiIPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIPhiIPhi_isLoaded = true;
		}
		return scs_sigmaIPhiIPhi_;
	}
	vector<float> &scs_sigmaPhiPhi()
	{
		if (not scs_sigmaPhiPhi_isLoaded) {
			if (scs_sigmaPhiPhi_branch != 0) {
				scs_sigmaPhiPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = scs_sigmaPhiPhi_.begin(); i != scs_sigmaPhiPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch scs_sigmaPhiPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_sigmaPhiPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaPhiPhi_isLoaded = true;
		}
		return scs_sigmaPhiPhi_;
	}
	vector<float> &trks_chi2()
	{
		if (not trks_chi2_isLoaded) {
			if (trks_chi2_branch != 0) {
				trks_chi2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_chi2_.begin(); i != trks_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_chi2_branch does not exist!\n");
				exit(1);
			}
			trks_chi2_isLoaded = true;
		}
		return trks_chi2_;
	}
	vector<float> &trks_d0()
	{
		if (not trks_d0_isLoaded) {
			if (trks_d0_branch != 0) {
				trks_d0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_d0_.begin(); i != trks_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_d0_branch does not exist!\n");
				exit(1);
			}
			trks_d0_isLoaded = true;
		}
		return trks_d0_;
	}
	vector<float> &trks_d0Err()
	{
		if (not trks_d0Err_isLoaded) {
			if (trks_d0Err_branch != 0) {
				trks_d0Err_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_d0Err_.begin(); i != trks_d0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_d0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_d0Err_branch does not exist!\n");
				exit(1);
			}
			trks_d0Err_isLoaded = true;
		}
		return trks_d0Err_;
	}
	vector<float> &trks_d0corr()
	{
		if (not trks_d0corr_isLoaded) {
			if (trks_d0corr_branch != 0) {
				trks_d0corr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_d0corr_.begin(); i != trks_d0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_d0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_d0corr_branch does not exist!\n");
				exit(1);
			}
			trks_d0corr_isLoaded = true;
		}
		return trks_d0corr_;
	}
	vector<float> &trks_d0corrPhi()
	{
		if (not trks_d0corrPhi_isLoaded) {
			if (trks_d0corrPhi_branch != 0) {
				trks_d0corrPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_d0corrPhi_.begin(); i != trks_d0corrPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_d0corrPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_d0corrPhi_branch does not exist!\n");
				exit(1);
			}
			trks_d0corrPhi_isLoaded = true;
		}
		return trks_d0corrPhi_;
	}
	vector<float> &trks_etaErr()
	{
		if (not trks_etaErr_isLoaded) {
			if (trks_etaErr_branch != 0) {
				trks_etaErr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_etaErr_.begin(); i != trks_etaErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_etaErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_etaErr_branch does not exist!\n");
				exit(1);
			}
			trks_etaErr_isLoaded = true;
		}
		return trks_etaErr_;
	}
	vector<float> &trks_layer1_charge()
	{
		if (not trks_layer1_charge_isLoaded) {
			if (trks_layer1_charge_branch != 0) {
				trks_layer1_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_layer1_charge_.begin(); i != trks_layer1_charge_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_layer1_charge_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_layer1_charge_branch does not exist!\n");
				exit(1);
			}
			trks_layer1_charge_isLoaded = true;
		}
		return trks_layer1_charge_;
	}
	vector<float> &trks_ndof()
	{
		if (not trks_ndof_isLoaded) {
			if (trks_ndof_branch != 0) {
				trks_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_ndof_.begin(); i != trks_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_ndof_branch does not exist!\n");
				exit(1);
			}
			trks_ndof_isLoaded = true;
		}
		return trks_ndof_;
	}
	vector<float> &trks_phiErr()
	{
		if (not trks_phiErr_isLoaded) {
			if (trks_phiErr_branch != 0) {
				trks_phiErr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_phiErr_.begin(); i != trks_phiErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_phiErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_phiErr_branch does not exist!\n");
				exit(1);
			}
			trks_phiErr_isLoaded = true;
		}
		return trks_phiErr_;
	}
	vector<float> &trks_ptErr()
	{
		if (not trks_ptErr_isLoaded) {
			if (trks_ptErr_branch != 0) {
				trks_ptErr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_ptErr_.begin(); i != trks_ptErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_ptErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_ptErr_branch does not exist!\n");
				exit(1);
			}
			trks_ptErr_isLoaded = true;
		}
		return trks_ptErr_;
	}
	vector<float> &trks_z0()
	{
		if (not trks_z0_isLoaded) {
			if (trks_z0_branch != 0) {
				trks_z0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_z0_.begin(); i != trks_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_z0_branch does not exist!\n");
				exit(1);
			}
			trks_z0_isLoaded = true;
		}
		return trks_z0_;
	}
	vector<float> &trks_z0Err()
	{
		if (not trks_z0Err_isLoaded) {
			if (trks_z0Err_branch != 0) {
				trks_z0Err_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_z0Err_.begin(); i != trks_z0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_z0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_z0Err_branch does not exist!\n");
				exit(1);
			}
			trks_z0Err_isLoaded = true;
		}
		return trks_z0Err_;
	}
	vector<float> &trks_z0corr()
	{
		if (not trks_z0corr_isLoaded) {
			if (trks_z0corr_branch != 0) {
				trks_z0corr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_z0corr_.begin(); i != trks_z0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_z0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_z0corr_branch does not exist!\n");
				exit(1);
			}
			trks_z0corr_isLoaded = true;
		}
		return trks_z0corr_;
	}
	vector<float> &vtxs_chi2()
	{
		if (not vtxs_chi2_isLoaded) {
			if (vtxs_chi2_branch != 0) {
				vtxs_chi2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = vtxs_chi2_.begin(); i != vtxs_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch vtxs_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_chi2_branch does not exist!\n");
				exit(1);
			}
			vtxs_chi2_isLoaded = true;
		}
		return vtxs_chi2_;
	}
	vector<float> &vtxs_ndof()
	{
		if (not vtxs_ndof_isLoaded) {
			if (vtxs_ndof_branch != 0) {
				vtxs_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = vtxs_ndof_.begin(); i != vtxs_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch vtxs_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_ndof_branch does not exist!\n");
				exit(1);
			}
			vtxs_ndof_isLoaded = true;
		}
		return vtxs_ndof_;
	}
	vector<float> &vtxs_xError()
	{
		if (not vtxs_xError_isLoaded) {
			if (vtxs_xError_branch != 0) {
				vtxs_xError_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = vtxs_xError_.begin(); i != vtxs_xError_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch vtxs_xError_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_xError_branch does not exist!\n");
				exit(1);
			}
			vtxs_xError_isLoaded = true;
		}
		return vtxs_xError_;
	}
	vector<float> &vtxs_yError()
	{
		if (not vtxs_yError_isLoaded) {
			if (vtxs_yError_branch != 0) {
				vtxs_yError_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = vtxs_yError_.begin(); i != vtxs_yError_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch vtxs_yError_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_yError_branch does not exist!\n");
				exit(1);
			}
			vtxs_yError_isLoaded = true;
		}
		return vtxs_yError_;
	}
	vector<float> &vtxs_zError()
	{
		if (not vtxs_zError_isLoaded) {
			if (vtxs_zError_branch != 0) {
				vtxs_zError_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = vtxs_zError_.begin(); i != vtxs_zError_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch vtxs_zError_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_zError_branch does not exist!\n");
				exit(1);
			}
			vtxs_zError_isLoaded = true;
		}
		return vtxs_zError_;
	}
	vector<vector<float> > &pseudo_ecalIso03_recHitE()
	{
		if (not pseudo_ecalIso03_recHitE_isLoaded) {
			if (pseudo_ecalIso03_recHitE_branch != 0) {
				pseudo_ecalIso03_recHitE_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = pseudo_ecalIso03_recHitE_.begin(); i != pseudo_ecalIso03_recHitE_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch pseudo_ecalIso03_recHitE_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_ecalIso03_recHitE_branch does not exist!\n");
				exit(1);
			}
			pseudo_ecalIso03_recHitE_isLoaded = true;
		}
		return pseudo_ecalIso03_recHitE_;
	}
	vector<vector<float> > &pseudo_ecalIso03_recHitEt()
	{
		if (not pseudo_ecalIso03_recHitEt_isLoaded) {
			if (pseudo_ecalIso03_recHitEt_branch != 0) {
				pseudo_ecalIso03_recHitEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = pseudo_ecalIso03_recHitEt_.begin(); i != pseudo_ecalIso03_recHitEt_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch pseudo_ecalIso03_recHitEt_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_ecalIso03_recHitEt_branch does not exist!\n");
				exit(1);
			}
			pseudo_ecalIso03_recHitEt_isLoaded = true;
		}
		return pseudo_ecalIso03_recHitEt_;
	}
	vector<vector<float> > &pseudo_srDIdx()
	{
		if (not pseudo_srDIdx_isLoaded) {
			if (pseudo_srDIdx_branch != 0) {
				pseudo_srDIdx_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = pseudo_srDIdx_.begin(); i != pseudo_srDIdx_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch pseudo_srDIdx_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_srDIdx_branch does not exist!\n");
				exit(1);
			}
			pseudo_srDIdx_isLoaded = true;
		}
		return pseudo_srDIdx_;
	}
	vector<vector<float> > &vtxs_covMatrix()
	{
		if (not vtxs_covMatrix_isLoaded) {
			if (vtxs_covMatrix_branch != 0) {
				vtxs_covMatrix_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = vtxs_covMatrix_.begin(); i != vtxs_covMatrix_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch vtxs_covMatrix_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_covMatrix_branch does not exist!\n");
				exit(1);
			}
			vtxs_covMatrix_isLoaded = true;
		}
		return vtxs_covMatrix_;
	}
	int &evt_bunchCrossing()
	{
		if (not evt_bunchCrossing_isLoaded) {
			if (evt_bunchCrossing_branch != 0) {
				evt_bunchCrossing_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bunchCrossing_branch does not exist!\n");
				exit(1);
			}
			evt_bunchCrossing_isLoaded = true;
		}
		return evt_bunchCrossing_;
	}
	int &evt_experimentType()
	{
		if (not evt_experimentType_isLoaded) {
			if (evt_experimentType_branch != 0) {
				evt_experimentType_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_experimentType_branch does not exist!\n");
				exit(1);
			}
			evt_experimentType_isLoaded = true;
		}
		return evt_experimentType_;
	}
	int &evt_orbitNumber()
	{
		if (not evt_orbitNumber_isLoaded) {
			if (evt_orbitNumber_branch != 0) {
				evt_orbitNumber_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_orbitNumber_branch does not exist!\n");
				exit(1);
			}
			evt_orbitNumber_isLoaded = true;
		}
		return evt_orbitNumber_;
	}
	int &evt_storeNumber()
	{
		if (not evt_storeNumber_isLoaded) {
			if (evt_storeNumber_branch != 0) {
				evt_storeNumber_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_storeNumber_branch does not exist!\n");
				exit(1);
			}
			evt_storeNumber_isLoaded = true;
		}
		return evt_storeNumber_;
	}
	int &hcalnoise_maxHPDHits()
	{
		if (not hcalnoise_maxHPDHits_isLoaded) {
			if (hcalnoise_maxHPDHits_branch != 0) {
				hcalnoise_maxHPDHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_maxHPDHits_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_maxHPDHits_isLoaded = true;
		}
		return hcalnoise_maxHPDHits_;
	}
	int &hcalnoise_maxRBXHits()
	{
		if (not hcalnoise_maxRBXHits_isLoaded) {
			if (hcalnoise_maxRBXHits_branch != 0) {
				hcalnoise_maxRBXHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_maxRBXHits_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_maxRBXHits_isLoaded = true;
		}
		return hcalnoise_maxRBXHits_;
	}
	int &hcalnoise_maxZeros()
	{
		if (not hcalnoise_maxZeros_isLoaded) {
			if (hcalnoise_maxZeros_branch != 0) {
				hcalnoise_maxZeros_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_maxZeros_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_maxZeros_isLoaded = true;
		}
		return hcalnoise_maxZeros_;
	}
	int &hcalnoise_noiseFilterStatus()
	{
		if (not hcalnoise_noiseFilterStatus_isLoaded) {
			if (hcalnoise_noiseFilterStatus_branch != 0) {
				hcalnoise_noiseFilterStatus_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_noiseFilterStatus_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_noiseFilterStatus_isLoaded = true;
		}
		return hcalnoise_noiseFilterStatus_;
	}
	int &hcalnoise_noiseType()
	{
		if (not hcalnoise_noiseType_isLoaded) {
			if (hcalnoise_noiseType_branch != 0) {
				hcalnoise_noiseType_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_noiseType_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_noiseType_isLoaded = true;
		}
		return hcalnoise_noiseType_;
	}
	int &hcalnoise_num10GeVHits()
	{
		if (not hcalnoise_num10GeVHits_isLoaded) {
			if (hcalnoise_num10GeVHits_branch != 0) {
				hcalnoise_num10GeVHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_num10GeVHits_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_num10GeVHits_isLoaded = true;
		}
		return hcalnoise_num10GeVHits_;
	}
	int &hcalnoise_num25GeVHits()
	{
		if (not hcalnoise_num25GeVHits_isLoaded) {
			if (hcalnoise_num25GeVHits_branch != 0) {
				hcalnoise_num25GeVHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_num25GeVHits_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_num25GeVHits_isLoaded = true;
		}
		return hcalnoise_num25GeVHits_;
	}
	int &hcalnoise_numProblematicRBXs()
	{
		if (not hcalnoise_numProblematicRBXs_isLoaded) {
			if (hcalnoise_numProblematicRBXs_branch != 0) {
				hcalnoise_numProblematicRBXs_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_numProblematicRBXs_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_numProblematicRBXs_isLoaded = true;
		}
		return hcalnoise_numProblematicRBXs_;
	}
	int &hcalnoise_passHighLevelNoiseFilter()
	{
		if (not hcalnoise_passHighLevelNoiseFilter_isLoaded) {
			if (hcalnoise_passHighLevelNoiseFilter_branch != 0) {
				hcalnoise_passHighLevelNoiseFilter_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_passHighLevelNoiseFilter_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_passHighLevelNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passHighLevelNoiseFilter_;
	}
	int &hcalnoise_passLooseNoiseFilter()
	{
		if (not hcalnoise_passLooseNoiseFilter_isLoaded) {
			if (hcalnoise_passLooseNoiseFilter_branch != 0) {
				hcalnoise_passLooseNoiseFilter_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_passLooseNoiseFilter_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_passLooseNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passLooseNoiseFilter_;
	}
	int &hcalnoise_passTightNoiseFilter()
	{
		if (not hcalnoise_passTightNoiseFilter_isLoaded) {
			if (hcalnoise_passTightNoiseFilter_branch != 0) {
				hcalnoise_passTightNoiseFilter_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hcalnoise_passTightNoiseFilter_branch does not exist!\n");
				exit(1);
			}
			hcalnoise_passTightNoiseFilter_isLoaded = true;
		}
		return hcalnoise_passTightNoiseFilter_;
	}
	int &l1_nemiso()
	{
		if (not l1_nemiso_isLoaded) {
			if (l1_nemiso_branch != 0) {
				l1_nemiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_nemiso_branch does not exist!\n");
				exit(1);
			}
			l1_nemiso_isLoaded = true;
		}
		return l1_nemiso_;
	}
	int &l1_nemnoiso()
	{
		if (not l1_nemnoiso_isLoaded) {
			if (l1_nemnoiso_branch != 0) {
				l1_nemnoiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_nemnoiso_branch does not exist!\n");
				exit(1);
			}
			l1_nemnoiso_isLoaded = true;
		}
		return l1_nemnoiso_;
	}
	int &l1_njetsc()
	{
		if (not l1_njetsc_isLoaded) {
			if (l1_njetsc_branch != 0) {
				l1_njetsc_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_njetsc_branch does not exist!\n");
				exit(1);
			}
			l1_njetsc_isLoaded = true;
		}
		return l1_njetsc_;
	}
	int &l1_njetsf()
	{
		if (not l1_njetsf_isLoaded) {
			if (l1_njetsf_branch != 0) {
				l1_njetsf_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_njetsf_branch does not exist!\n");
				exit(1);
			}
			l1_njetsf_isLoaded = true;
		}
		return l1_njetsf_;
	}
	int &l1_njetst()
	{
		if (not l1_njetst_isLoaded) {
			if (l1_njetst_branch != 0) {
				l1_njetst_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_njetst_branch does not exist!\n");
				exit(1);
			}
			l1_njetst_isLoaded = true;
		}
		return l1_njetst_;
	}
	int &l1_nmus()
	{
		if (not l1_nmus_isLoaded) {
			if (l1_nmus_branch != 0) {
				l1_nmus_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_nmus_branch does not exist!\n");
				exit(1);
			}
			l1_nmus_isLoaded = true;
		}
		return l1_nmus_;
	}
	vector<int> &l1_emiso_ieta()
	{
		if (not l1_emiso_ieta_isLoaded) {
			if (l1_emiso_ieta_branch != 0) {
				l1_emiso_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emiso_ieta_branch does not exist!\n");
				exit(1);
			}
			l1_emiso_ieta_isLoaded = true;
		}
		return l1_emiso_ieta_;
	}
	vector<int> &l1_emiso_iphi()
	{
		if (not l1_emiso_iphi_isLoaded) {
			if (l1_emiso_iphi_branch != 0) {
				l1_emiso_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emiso_iphi_branch does not exist!\n");
				exit(1);
			}
			l1_emiso_iphi_isLoaded = true;
		}
		return l1_emiso_iphi_;
	}
	vector<int> &l1_emiso_rawId()
	{
		if (not l1_emiso_rawId_isLoaded) {
			if (l1_emiso_rawId_branch != 0) {
				l1_emiso_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emiso_rawId_branch does not exist!\n");
				exit(1);
			}
			l1_emiso_rawId_isLoaded = true;
		}
		return l1_emiso_rawId_;
	}
	vector<int> &l1_emiso_type()
	{
		if (not l1_emiso_type_isLoaded) {
			if (l1_emiso_type_branch != 0) {
				l1_emiso_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emiso_type_branch does not exist!\n");
				exit(1);
			}
			l1_emiso_type_isLoaded = true;
		}
		return l1_emiso_type_;
	}
	vector<int> &l1_emnoiso_ieta()
	{
		if (not l1_emnoiso_ieta_isLoaded) {
			if (l1_emnoiso_ieta_branch != 0) {
				l1_emnoiso_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emnoiso_ieta_branch does not exist!\n");
				exit(1);
			}
			l1_emnoiso_ieta_isLoaded = true;
		}
		return l1_emnoiso_ieta_;
	}
	vector<int> &l1_emnoiso_iphi()
	{
		if (not l1_emnoiso_iphi_isLoaded) {
			if (l1_emnoiso_iphi_branch != 0) {
				l1_emnoiso_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emnoiso_iphi_branch does not exist!\n");
				exit(1);
			}
			l1_emnoiso_iphi_isLoaded = true;
		}
		return l1_emnoiso_iphi_;
	}
	vector<int> &l1_emnoiso_rawId()
	{
		if (not l1_emnoiso_rawId_isLoaded) {
			if (l1_emnoiso_rawId_branch != 0) {
				l1_emnoiso_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emnoiso_rawId_branch does not exist!\n");
				exit(1);
			}
			l1_emnoiso_rawId_isLoaded = true;
		}
		return l1_emnoiso_rawId_;
	}
	vector<int> &l1_emnoiso_type()
	{
		if (not l1_emnoiso_type_isLoaded) {
			if (l1_emnoiso_type_branch != 0) {
				l1_emnoiso_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_emnoiso_type_branch does not exist!\n");
				exit(1);
			}
			l1_emnoiso_type_isLoaded = true;
		}
		return l1_emnoiso_type_;
	}
	vector<int> &l1_jetsc_ieta()
	{
		if (not l1_jetsc_ieta_isLoaded) {
			if (l1_jetsc_ieta_branch != 0) {
				l1_jetsc_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsc_ieta_branch does not exist!\n");
				exit(1);
			}
			l1_jetsc_ieta_isLoaded = true;
		}
		return l1_jetsc_ieta_;
	}
	vector<int> &l1_jetsc_iphi()
	{
		if (not l1_jetsc_iphi_isLoaded) {
			if (l1_jetsc_iphi_branch != 0) {
				l1_jetsc_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsc_iphi_branch does not exist!\n");
				exit(1);
			}
			l1_jetsc_iphi_isLoaded = true;
		}
		return l1_jetsc_iphi_;
	}
	vector<int> &l1_jetsc_rawId()
	{
		if (not l1_jetsc_rawId_isLoaded) {
			if (l1_jetsc_rawId_branch != 0) {
				l1_jetsc_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsc_rawId_branch does not exist!\n");
				exit(1);
			}
			l1_jetsc_rawId_isLoaded = true;
		}
		return l1_jetsc_rawId_;
	}
	vector<int> &l1_jetsc_type()
	{
		if (not l1_jetsc_type_isLoaded) {
			if (l1_jetsc_type_branch != 0) {
				l1_jetsc_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsc_type_branch does not exist!\n");
				exit(1);
			}
			l1_jetsc_type_isLoaded = true;
		}
		return l1_jetsc_type_;
	}
	vector<int> &l1_jetsf_ieta()
	{
		if (not l1_jetsf_ieta_isLoaded) {
			if (l1_jetsf_ieta_branch != 0) {
				l1_jetsf_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsf_ieta_branch does not exist!\n");
				exit(1);
			}
			l1_jetsf_ieta_isLoaded = true;
		}
		return l1_jetsf_ieta_;
	}
	vector<int> &l1_jetsf_iphi()
	{
		if (not l1_jetsf_iphi_isLoaded) {
			if (l1_jetsf_iphi_branch != 0) {
				l1_jetsf_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsf_iphi_branch does not exist!\n");
				exit(1);
			}
			l1_jetsf_iphi_isLoaded = true;
		}
		return l1_jetsf_iphi_;
	}
	vector<int> &l1_jetsf_rawId()
	{
		if (not l1_jetsf_rawId_isLoaded) {
			if (l1_jetsf_rawId_branch != 0) {
				l1_jetsf_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsf_rawId_branch does not exist!\n");
				exit(1);
			}
			l1_jetsf_rawId_isLoaded = true;
		}
		return l1_jetsf_rawId_;
	}
	vector<int> &l1_jetsf_type()
	{
		if (not l1_jetsf_type_isLoaded) {
			if (l1_jetsf_type_branch != 0) {
				l1_jetsf_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetsf_type_branch does not exist!\n");
				exit(1);
			}
			l1_jetsf_type_isLoaded = true;
		}
		return l1_jetsf_type_;
	}
	vector<int> &l1_jetst_ieta()
	{
		if (not l1_jetst_ieta_isLoaded) {
			if (l1_jetst_ieta_branch != 0) {
				l1_jetst_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetst_ieta_branch does not exist!\n");
				exit(1);
			}
			l1_jetst_ieta_isLoaded = true;
		}
		return l1_jetst_ieta_;
	}
	vector<int> &l1_jetst_iphi()
	{
		if (not l1_jetst_iphi_isLoaded) {
			if (l1_jetst_iphi_branch != 0) {
				l1_jetst_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetst_iphi_branch does not exist!\n");
				exit(1);
			}
			l1_jetst_iphi_isLoaded = true;
		}
		return l1_jetst_iphi_;
	}
	vector<int> &l1_jetst_rawId()
	{
		if (not l1_jetst_rawId_isLoaded) {
			if (l1_jetst_rawId_branch != 0) {
				l1_jetst_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetst_rawId_branch does not exist!\n");
				exit(1);
			}
			l1_jetst_rawId_isLoaded = true;
		}
		return l1_jetst_rawId_;
	}
	vector<int> &l1_jetst_type()
	{
		if (not l1_jetst_type_isLoaded) {
			if (l1_jetst_type_branch != 0) {
				l1_jetst_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_jetst_type_branch does not exist!\n");
				exit(1);
			}
			l1_jetst_type_isLoaded = true;
		}
		return l1_jetst_type_;
	}
	vector<int> &l1_mus_flags()
	{
		if (not l1_mus_flags_isLoaded) {
			if (l1_mus_flags_branch != 0) {
				l1_mus_flags_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mus_flags_branch does not exist!\n");
				exit(1);
			}
			l1_mus_flags_isLoaded = true;
		}
		return l1_mus_flags_;
	}
	vector<int> &l1_mus_q()
	{
		if (not l1_mus_q_isLoaded) {
			if (l1_mus_q_branch != 0) {
				l1_mus_q_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mus_q_branch does not exist!\n");
				exit(1);
			}
			l1_mus_q_isLoaded = true;
		}
		return l1_mus_q_;
	}
	vector<int> &l1_mus_qual()
	{
		if (not l1_mus_qual_isLoaded) {
			if (l1_mus_qual_branch != 0) {
				l1_mus_qual_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mus_qual_branch does not exist!\n");
				exit(1);
			}
			l1_mus_qual_isLoaded = true;
		}
		return l1_mus_qual_;
	}
	vector<int> &l1_mus_qualFlags()
	{
		if (not l1_mus_qualFlags_isLoaded) {
			if (l1_mus_qualFlags_branch != 0) {
				l1_mus_qualFlags_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_mus_qualFlags_branch does not exist!\n");
				exit(1);
			}
			l1_mus_qualFlags_isLoaded = true;
		}
		return l1_mus_qualFlags_;
	}
	vector<int> &pxl_ndigis_pxb()
	{
		if (not pxl_ndigis_pxb_isLoaded) {
			if (pxl_ndigis_pxb_branch != 0) {
				pxl_ndigis_pxb_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pxl_ndigis_pxb_branch does not exist!\n");
				exit(1);
			}
			pxl_ndigis_pxb_isLoaded = true;
		}
		return pxl_ndigis_pxb_;
	}
	vector<int> &pxl_ndigis_pxf()
	{
		if (not pxl_ndigis_pxf_isLoaded) {
			if (pxl_ndigis_pxf_branch != 0) {
				pxl_ndigis_pxf_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pxl_ndigis_pxf_branch does not exist!\n");
				exit(1);
			}
			pxl_ndigis_pxf_isLoaded = true;
		}
		return pxl_ndigis_pxf_;
	}
	vector<int> &scs_elsidx()
	{
		if (not scs_elsidx_isLoaded) {
			if (scs_elsidx_branch != 0) {
				scs_elsidx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_elsidx_branch does not exist!\n");
				exit(1);
			}
			scs_elsidx_isLoaded = true;
		}
		return scs_elsidx_;
	}
	vector<int> &scs_severitySeed()
	{
		if (not scs_severitySeed_isLoaded) {
			if (scs_severitySeed_branch != 0) {
				scs_severitySeed_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_severitySeed_branch does not exist!\n");
				exit(1);
			}
			scs_severitySeed_isLoaded = true;
		}
		return scs_severitySeed_;
	}
	vector<int> &trks_algo()
	{
		if (not trks_algo_isLoaded) {
			if (trks_algo_branch != 0) {
				trks_algo_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_algo_branch does not exist!\n");
				exit(1);
			}
			trks_algo_isLoaded = true;
		}
		return trks_algo_;
	}
	vector<int> &trks_charge()
	{
		if (not trks_charge_isLoaded) {
			if (trks_charge_branch != 0) {
				trks_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_charge_branch does not exist!\n");
				exit(1);
			}
			trks_charge_isLoaded = true;
		}
		return trks_charge_;
	}
	vector<int> &trks_exp_innerlayers()
	{
		if (not trks_exp_innerlayers_isLoaded) {
			if (trks_exp_innerlayers_branch != 0) {
				trks_exp_innerlayers_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_exp_innerlayers_branch does not exist!\n");
				exit(1);
			}
			trks_exp_innerlayers_isLoaded = true;
		}
		return trks_exp_innerlayers_;
	}
	vector<int> &trks_exp_outerlayers()
	{
		if (not trks_exp_outerlayers_isLoaded) {
			if (trks_exp_outerlayers_branch != 0) {
				trks_exp_outerlayers_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_exp_outerlayers_branch does not exist!\n");
				exit(1);
			}
			trks_exp_outerlayers_isLoaded = true;
		}
		return trks_exp_outerlayers_;
	}
	vector<int> &trks_layer1_det()
	{
		if (not trks_layer1_det_isLoaded) {
			if (trks_layer1_det_branch != 0) {
				trks_layer1_det_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_layer1_det_branch does not exist!\n");
				exit(1);
			}
			trks_layer1_det_isLoaded = true;
		}
		return trks_layer1_det_;
	}
	vector<int> &trks_layer1_layer()
	{
		if (not trks_layer1_layer_isLoaded) {
			if (trks_layer1_layer_branch != 0) {
				trks_layer1_layer_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_layer1_layer_branch does not exist!\n");
				exit(1);
			}
			trks_layer1_layer_isLoaded = true;
		}
		return trks_layer1_layer_;
	}
	vector<int> &trks_layer1_sizerphi()
	{
		if (not trks_layer1_sizerphi_isLoaded) {
			if (trks_layer1_sizerphi_branch != 0) {
				trks_layer1_sizerphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_layer1_sizerphi_branch does not exist!\n");
				exit(1);
			}
			trks_layer1_sizerphi_isLoaded = true;
		}
		return trks_layer1_sizerphi_;
	}
	vector<int> &trks_layer1_sizerz()
	{
		if (not trks_layer1_sizerz_isLoaded) {
			if (trks_layer1_sizerz_branch != 0) {
				trks_layer1_sizerz_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_layer1_sizerz_branch does not exist!\n");
				exit(1);
			}
			trks_layer1_sizerz_isLoaded = true;
		}
		return trks_layer1_sizerz_;
	}
	vector<int> &trks_lostHits()
	{
		if (not trks_lostHits_isLoaded) {
			if (trks_lostHits_branch != 0) {
				trks_lostHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_lostHits_branch does not exist!\n");
				exit(1);
			}
			trks_lostHits_isLoaded = true;
		}
		return trks_lostHits_;
	}
	vector<int> &trks_lost_pixelhits()
	{
		if (not trks_lost_pixelhits_isLoaded) {
			if (trks_lost_pixelhits_branch != 0) {
				trks_lost_pixelhits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_lost_pixelhits_branch does not exist!\n");
				exit(1);
			}
			trks_lost_pixelhits_isLoaded = true;
		}
		return trks_lost_pixelhits_;
	}
	vector<int> &trks_nlayers()
	{
		if (not trks_nlayers_isLoaded) {
			if (trks_nlayers_branch != 0) {
				trks_nlayers_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_nlayers_branch does not exist!\n");
				exit(1);
			}
			trks_nlayers_isLoaded = true;
		}
		return trks_nlayers_;
	}
	vector<int> &trks_nlayers3D()
	{
		if (not trks_nlayers3D_isLoaded) {
			if (trks_nlayers3D_branch != 0) {
				trks_nlayers3D_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_nlayers3D_branch does not exist!\n");
				exit(1);
			}
			trks_nlayers3D_isLoaded = true;
		}
		return trks_nlayers3D_;
	}
	vector<int> &trks_nlayersLost()
	{
		if (not trks_nlayersLost_isLoaded) {
			if (trks_nlayersLost_branch != 0) {
				trks_nlayersLost_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_nlayersLost_branch does not exist!\n");
				exit(1);
			}
			trks_nlayersLost_isLoaded = true;
		}
		return trks_nlayersLost_;
	}
	vector<int> &trks_qualityMask()
	{
		if (not trks_qualityMask_isLoaded) {
			if (trks_qualityMask_branch != 0) {
				trks_qualityMask_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_qualityMask_branch does not exist!\n");
				exit(1);
			}
			trks_qualityMask_isLoaded = true;
		}
		return trks_qualityMask_;
	}
	vector<int> &trks_validHits()
	{
		if (not trks_validHits_isLoaded) {
			if (trks_validHits_branch != 0) {
				trks_validHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_validHits_branch does not exist!\n");
				exit(1);
			}
			trks_validHits_isLoaded = true;
		}
		return trks_validHits_;
	}
	vector<int> &trks_valid_pixelhits()
	{
		if (not trks_valid_pixelhits_isLoaded) {
			if (trks_valid_pixelhits_branch != 0) {
				trks_valid_pixelhits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_valid_pixelhits_branch does not exist!\n");
				exit(1);
			}
			trks_valid_pixelhits_isLoaded = true;
		}
		return trks_valid_pixelhits_;
	}
	vector<int> &vtxs_isFake()
	{
		if (not vtxs_isFake_isLoaded) {
			if (vtxs_isFake_branch != 0) {
				vtxs_isFake_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_isFake_branch does not exist!\n");
				exit(1);
			}
			vtxs_isFake_isLoaded = true;
		}
		return vtxs_isFake_;
	}
	vector<int> &vtxs_isValid()
	{
		if (not vtxs_isValid_isLoaded) {
			if (vtxs_isValid_branch != 0) {
				vtxs_isValid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_isValid_branch does not exist!\n");
				exit(1);
			}
			vtxs_isValid_isLoaded = true;
		}
		return vtxs_isValid_;
	}
	vector<int> &vtxs_tracksSize()
	{
		if (not vtxs_tracksSize_isLoaded) {
			if (vtxs_tracksSize_branch != 0) {
				vtxs_tracksSize_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch vtxs_tracksSize_branch does not exist!\n");
				exit(1);
			}
			vtxs_tracksSize_isLoaded = true;
		}
		return vtxs_tracksSize_;
	}
	vector<vector<int> > &pseudo_srFlags()
	{
		if (not pseudo_srFlags_isLoaded) {
			if (pseudo_srFlags_branch != 0) {
				pseudo_srFlags_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pseudo_srFlags_branch does not exist!\n");
				exit(1);
			}
			pseudo_srFlags_isLoaded = true;
		}
		return pseudo_srFlags_;
	}
	unsigned int &evt_ntwrs()
	{
		if (not evt_ntwrs_isLoaded) {
			if (evt_ntwrs_branch != 0) {
				evt_ntwrs_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_ntwrs_branch does not exist!\n");
				exit(1);
			}
			evt_ntwrs_isLoaded = true;
		}
		return evt_ntwrs_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) {
				evt_event_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_event_branch does not exist!\n");
				exit(1);
			}
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	unsigned int &evt_lumiBlock()
	{
		if (not evt_lumiBlock_isLoaded) {
			if (evt_lumiBlock_branch != 0) {
				evt_lumiBlock_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_lumiBlock_branch does not exist!\n");
				exit(1);
			}
			evt_lumiBlock_isLoaded = true;
		}
		return evt_lumiBlock_;
	}
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) {
				evt_run_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_run_branch does not exist!\n");
				exit(1);
			}
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &l1_bits1()
	{
		if (not l1_bits1_isLoaded) {
			if (l1_bits1_branch != 0) {
				l1_bits1_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_bits1_branch does not exist!\n");
				exit(1);
			}
			l1_bits1_isLoaded = true;
		}
		return l1_bits1_;
	}
	unsigned int &l1_bits2()
	{
		if (not l1_bits2_isLoaded) {
			if (l1_bits2_branch != 0) {
				l1_bits2_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_bits2_branch does not exist!\n");
				exit(1);
			}
			l1_bits2_isLoaded = true;
		}
		return l1_bits2_;
	}
	unsigned int &l1_bits3()
	{
		if (not l1_bits3_isLoaded) {
			if (l1_bits3_branch != 0) {
				l1_bits3_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_bits3_branch does not exist!\n");
				exit(1);
			}
			l1_bits3_isLoaded = true;
		}
		return l1_bits3_;
	}
	unsigned int &l1_bits4()
	{
		if (not l1_bits4_isLoaded) {
			if (l1_bits4_branch != 0) {
				l1_bits4_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_bits4_branch does not exist!\n");
				exit(1);
			}
			l1_bits4_isLoaded = true;
		}
		return l1_bits4_;
	}
	unsigned int &l1_techbits1()
	{
		if (not l1_techbits1_isLoaded) {
			if (l1_techbits1_branch != 0) {
				l1_techbits1_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_techbits1_branch does not exist!\n");
				exit(1);
			}
			l1_techbits1_isLoaded = true;
		}
		return l1_techbits1_;
	}
	unsigned int &l1_techbits2()
	{
		if (not l1_techbits2_isLoaded) {
			if (l1_techbits2_branch != 0) {
				l1_techbits2_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1_techbits2_branch does not exist!\n");
				exit(1);
			}
			l1_techbits2_isLoaded = true;
		}
		return l1_techbits2_;
	}
	unsigned int &evt_nscs()
	{
		if (not evt_nscs_isLoaded) {
			if (evt_nscs_branch != 0) {
				evt_nscs_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nscs_branch does not exist!\n");
				exit(1);
			}
			evt_nscs_isLoaded = true;
		}
		return evt_nscs_;
	}
	unsigned int &evt_nvtxs()
	{
		if (not evt_nvtxs_isLoaded) {
			if (evt_nvtxs_branch != 0) {
				evt_nvtxs_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nvtxs_branch does not exist!\n");
				exit(1);
			}
			evt_nvtxs_isLoaded = true;
		}
		return evt_nvtxs_;
	}
	vector<unsigned int> &twrs_numBadEcalCells()
	{
		if (not twrs_numBadEcalCells_isLoaded) {
			if (twrs_numBadEcalCells_branch != 0) {
				twrs_numBadEcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numBadEcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numBadEcalCells_isLoaded = true;
		}
		return twrs_numBadEcalCells_;
	}
	vector<unsigned int> &twrs_numBadHcalCells()
	{
		if (not twrs_numBadHcalCells_isLoaded) {
			if (twrs_numBadHcalCells_branch != 0) {
				twrs_numBadHcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numBadHcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numBadHcalCells_isLoaded = true;
		}
		return twrs_numBadHcalCells_;
	}
	vector<unsigned int> &twrs_numProblematicEcalCells()
	{
		if (not twrs_numProblematicEcalCells_isLoaded) {
			if (twrs_numProblematicEcalCells_branch != 0) {
				twrs_numProblematicEcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numProblematicEcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numProblematicEcalCells_isLoaded = true;
		}
		return twrs_numProblematicEcalCells_;
	}
	vector<unsigned int> &twrs_numProblematicHcalCells()
	{
		if (not twrs_numProblematicHcalCells_isLoaded) {
			if (twrs_numProblematicHcalCells_branch != 0) {
				twrs_numProblematicHcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numProblematicHcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numProblematicHcalCells_isLoaded = true;
		}
		return twrs_numProblematicHcalCells_;
	}
	vector<unsigned int> &twrs_numRecoveredEcalCells()
	{
		if (not twrs_numRecoveredEcalCells_isLoaded) {
			if (twrs_numRecoveredEcalCells_branch != 0) {
				twrs_numRecoveredEcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numRecoveredEcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numRecoveredEcalCells_isLoaded = true;
		}
		return twrs_numRecoveredEcalCells_;
	}
	vector<unsigned int> &twrs_numRecoveredHcalCells()
	{
		if (not twrs_numRecoveredHcalCells_isLoaded) {
			if (twrs_numRecoveredHcalCells_branch != 0) {
				twrs_numRecoveredHcalCells_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch twrs_numRecoveredHcalCells_branch does not exist!\n");
				exit(1);
			}
			twrs_numRecoveredHcalCells_isLoaded = true;
		}
		return twrs_numRecoveredHcalCells_;
	}
	bool passL1Trigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = l1_trigNames().begin();
		vector<TString>::const_iterator end_it = l1_trigNames().end();
		vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
		if(found_it != end_it)
			trigIndx = found_it - begin_it;
		else {
			cout << "Cannot find Trigger " << trigName << endl; 
			return 0;
		}

		if(trigIndx <= 31) {
			unsigned int bitmask = 1;
			bitmask <<= trigIndx;
			return l1_bits1() & bitmask;
		}
		if(trigIndx >= 32 && trigIndx <= 63) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 32); 
			return l1_bits2() & bitmask;
		}
		if(trigIndx >= 64 && trigIndx <= 95) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 64); 
			return l1_bits3() & bitmask;
		}
		if(trigIndx >= 96 && trigIndx <= 127) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 96); 
			return l1_bits4() & bitmask;
		}
	return 0;
	}
};

#ifndef __CINT__
extern CMS2 cms2;
#endif

namespace tas {
	TString &evt_CMS2tag() { return cms2.evt_CMS2tag(); }
	TString &evt_dataset() { return cms2.evt_dataset(); }
	vector<TString> &l1_trigNames() { return cms2.l1_trigNames(); }
	float &evt_bField() { return cms2.evt_bField(); }
	float &hcalnoise_eventChargeFraction() { return cms2.hcalnoise_eventChargeFraction(); }
	float &hcalnoise_eventEMEnergy() { return cms2.hcalnoise_eventEMEnergy(); }
	float &hcalnoise_eventEMFraction() { return cms2.hcalnoise_eventEMFraction(); }
	float &hcalnoise_eventHadEnergy() { return cms2.hcalnoise_eventHadEnergy(); }
	float &hcalnoise_eventTrackEnergy() { return cms2.hcalnoise_eventTrackEnergy(); }
	float &hcalnoise_max10GeVHitTime() { return cms2.hcalnoise_max10GeVHitTime(); }
	float &hcalnoise_max25GeVHitTime() { return cms2.hcalnoise_max25GeVHitTime(); }
	float &hcalnoise_min10GeVHitTime() { return cms2.hcalnoise_min10GeVHitTime(); }
	float &hcalnoise_min25GeVHitTime() { return cms2.hcalnoise_min25GeVHitTime(); }
	float &hcalnoise_minE10TS() { return cms2.hcalnoise_minE10TS(); }
	float &hcalnoise_minE2Over10TS() { return cms2.hcalnoise_minE2Over10TS(); }
	float &hcalnoise_minE2TS() { return cms2.hcalnoise_minE2TS(); }
	float &hcalnoise_minHPDEMF() { return cms2.hcalnoise_minHPDEMF(); }
	float &hcalnoise_minRBXEMF() { return cms2.hcalnoise_minRBXEMF(); }
	float &hcalnoise_rms10GeVHitTime() { return cms2.hcalnoise_rms10GeVHitTime(); }
	float &hcalnoise_rms25GeVHitTime() { return cms2.hcalnoise_rms25GeVHitTime(); }
	float &l1_met_etTot() { return cms2.l1_met_etTot(); }
	float &l1_met_met() { return cms2.l1_met_met(); }
	float &l1_mht_htTot() { return cms2.l1_mht_htTot(); }
	float &l1_mht_mht() { return cms2.l1_mht_mht(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_met_p4() { return cms2.l1_met_p4(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_mht_p4() { return cms2.l1_mht_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emiso_p4() { return cms2.l1_emiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emnoiso_p4() { return cms2.l1_emnoiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsc_p4() { return cms2.l1_jetsc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsf_p4() { return cms2.l1_jetsf_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetst_p4() { return cms2.l1_jetst_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_mus_p4() { return cms2.l1_mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_p4() { return cms2.scs_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_pos_p4() { return cms2.scs_pos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_vtx_p4() { return cms2.scs_vtx_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_inner_position() { return cms2.trks_inner_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_p4() { return cms2.trks_outer_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_position() { return cms2.trks_outer_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_trk_p4() { return cms2.trks_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_vertex_p4() { return cms2.trks_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position() { return cms2.vtxs_position(); }
	vector<float> &twrs_ecalTime() { return cms2.twrs_ecalTime(); }
	vector<float> &twrs_emEnergy() { return cms2.twrs_emEnergy(); }
	vector<float> &twrs_emEt() { return cms2.twrs_emEt(); }
	vector<float> &twrs_emEtcorr() { return cms2.twrs_emEtcorr(); }
	vector<float> &twrs_eta() { return cms2.twrs_eta(); }
	vector<float> &twrs_etcorr() { return cms2.twrs_etcorr(); }
	vector<float> &twrs_hadEnergy() { return cms2.twrs_hadEnergy(); }
	vector<float> &twrs_hadEt() { return cms2.twrs_hadEt(); }
	vector<float> &twrs_hadEtcorr() { return cms2.twrs_hadEtcorr(); }
	vector<float> &twrs_hcalTime() { return cms2.twrs_hcalTime(); }
	vector<float> &twrs_outerEnergy() { return cms2.twrs_outerEnergy(); }
	vector<float> &twrs_outerEt() { return cms2.twrs_outerEt(); }
	vector<float> &twrs_outerEtcorr() { return cms2.twrs_outerEtcorr(); }
	vector<float> &twrs_pcorr() { return cms2.twrs_pcorr(); }
	vector<float> &twrs_phi() { return cms2.twrs_phi(); }
	vector<float> &pseudo_dRClosestTower() { return cms2.pseudo_dRClosestTower(); }
	vector<float> &pseudo_dRClosestTowerEmEt() { return cms2.pseudo_dRClosestTowerEmEt(); }
	vector<float> &pseudo_ecalEta() { return cms2.pseudo_ecalEta(); }
	vector<float> &pseudo_ecalIso03() { return cms2.pseudo_ecalIso03(); }
	vector<float> &pseudo_ecalPhi() { return cms2.pseudo_ecalPhi(); }
	vector<float> &pseudo_eta() { return cms2.pseudo_eta(); }
	vector<float> &pseudo_hcalD1Iso03() { return cms2.pseudo_hcalD1Iso03(); }
	vector<float> &pseudo_hcalD2Iso03() { return cms2.pseudo_hcalD2Iso03(); }
	vector<float> &pseudo_phi() { return cms2.pseudo_phi(); }
	vector<float> &pseudo_tkIso03() { return cms2.pseudo_tkIso03(); }
	vector<float> &pseudo_towerEmEt() { return cms2.pseudo_towerEmEt(); }
	vector<float> &pseudo_towerHadEt() { return cms2.pseudo_towerHadEt(); }
	vector<float> &scs_clustersSize() { return cms2.scs_clustersSize(); }
	vector<float> &scs_crystalsSize() { return cms2.scs_crystalsSize(); }
	vector<float> &scs_e1x3() { return cms2.scs_e1x3(); }
	vector<float> &scs_e1x5() { return cms2.scs_e1x5(); }
	vector<float> &scs_e2nd() { return cms2.scs_e2nd(); }
	vector<float> &scs_e2x2() { return cms2.scs_e2x2(); }
	vector<float> &scs_e2x5Max() { return cms2.scs_e2x5Max(); }
	vector<float> &scs_e3x1() { return cms2.scs_e3x1(); }
	vector<float> &scs_e3x2() { return cms2.scs_e3x2(); }
	vector<float> &scs_e3x3() { return cms2.scs_e3x3(); }
	vector<float> &scs_e4x4() { return cms2.scs_e4x4(); }
	vector<float> &scs_e5x5() { return cms2.scs_e5x5(); }
	vector<float> &scs_eMax() { return cms2.scs_eMax(); }
	vector<float> &scs_eSeed() { return cms2.scs_eSeed(); }
	vector<float> &scs_energy() { return cms2.scs_energy(); }
	vector<float> &scs_eta() { return cms2.scs_eta(); }
	vector<float> &scs_hoe() { return cms2.scs_hoe(); }
	vector<float> &scs_phi() { return cms2.scs_phi(); }
	vector<float> &scs_preshowerEnergy() { return cms2.scs_preshowerEnergy(); }
	vector<float> &scs_rawEnergy() { return cms2.scs_rawEnergy(); }
	vector<float> &scs_sigmaEtaEta() { return cms2.scs_sigmaEtaEta(); }
	vector<float> &scs_sigmaEtaPhi() { return cms2.scs_sigmaEtaPhi(); }
	vector<float> &scs_sigmaIEtaIEta() { return cms2.scs_sigmaIEtaIEta(); }
	vector<float> &scs_sigmaIEtaIPhi() { return cms2.scs_sigmaIEtaIPhi(); }
	vector<float> &scs_sigmaIPhiIPhi() { return cms2.scs_sigmaIPhiIPhi(); }
	vector<float> &scs_sigmaPhiPhi() { return cms2.scs_sigmaPhiPhi(); }
	vector<float> &trks_chi2() { return cms2.trks_chi2(); }
	vector<float> &trks_d0() { return cms2.trks_d0(); }
	vector<float> &trks_d0Err() { return cms2.trks_d0Err(); }
	vector<float> &trks_d0corr() { return cms2.trks_d0corr(); }
	vector<float> &trks_d0corrPhi() { return cms2.trks_d0corrPhi(); }
	vector<float> &trks_etaErr() { return cms2.trks_etaErr(); }
	vector<float> &trks_layer1_charge() { return cms2.trks_layer1_charge(); }
	vector<float> &trks_ndof() { return cms2.trks_ndof(); }
	vector<float> &trks_phiErr() { return cms2.trks_phiErr(); }
	vector<float> &trks_ptErr() { return cms2.trks_ptErr(); }
	vector<float> &trks_z0() { return cms2.trks_z0(); }
	vector<float> &trks_z0Err() { return cms2.trks_z0Err(); }
	vector<float> &trks_z0corr() { return cms2.trks_z0corr(); }
	vector<float> &vtxs_chi2() { return cms2.vtxs_chi2(); }
	vector<float> &vtxs_ndof() { return cms2.vtxs_ndof(); }
	vector<float> &vtxs_xError() { return cms2.vtxs_xError(); }
	vector<float> &vtxs_yError() { return cms2.vtxs_yError(); }
	vector<float> &vtxs_zError() { return cms2.vtxs_zError(); }
	vector<vector<float> > &pseudo_ecalIso03_recHitE() { return cms2.pseudo_ecalIso03_recHitE(); }
	vector<vector<float> > &pseudo_ecalIso03_recHitEt() { return cms2.pseudo_ecalIso03_recHitEt(); }
	vector<vector<float> > &pseudo_srDIdx() { return cms2.pseudo_srDIdx(); }
	vector<vector<float> > &vtxs_covMatrix() { return cms2.vtxs_covMatrix(); }
	int &evt_bunchCrossing() { return cms2.evt_bunchCrossing(); }
	int &evt_experimentType() { return cms2.evt_experimentType(); }
	int &evt_orbitNumber() { return cms2.evt_orbitNumber(); }
	int &evt_storeNumber() { return cms2.evt_storeNumber(); }
	int &hcalnoise_maxHPDHits() { return cms2.hcalnoise_maxHPDHits(); }
	int &hcalnoise_maxRBXHits() { return cms2.hcalnoise_maxRBXHits(); }
	int &hcalnoise_maxZeros() { return cms2.hcalnoise_maxZeros(); }
	int &hcalnoise_noiseFilterStatus() { return cms2.hcalnoise_noiseFilterStatus(); }
	int &hcalnoise_noiseType() { return cms2.hcalnoise_noiseType(); }
	int &hcalnoise_num10GeVHits() { return cms2.hcalnoise_num10GeVHits(); }
	int &hcalnoise_num25GeVHits() { return cms2.hcalnoise_num25GeVHits(); }
	int &hcalnoise_numProblematicRBXs() { return cms2.hcalnoise_numProblematicRBXs(); }
	int &hcalnoise_passHighLevelNoiseFilter() { return cms2.hcalnoise_passHighLevelNoiseFilter(); }
	int &hcalnoise_passLooseNoiseFilter() { return cms2.hcalnoise_passLooseNoiseFilter(); }
	int &hcalnoise_passTightNoiseFilter() { return cms2.hcalnoise_passTightNoiseFilter(); }
	int &l1_nemiso() { return cms2.l1_nemiso(); }
	int &l1_nemnoiso() { return cms2.l1_nemnoiso(); }
	int &l1_njetsc() { return cms2.l1_njetsc(); }
	int &l1_njetsf() { return cms2.l1_njetsf(); }
	int &l1_njetst() { return cms2.l1_njetst(); }
	int &l1_nmus() { return cms2.l1_nmus(); }
	vector<int> &l1_emiso_ieta() { return cms2.l1_emiso_ieta(); }
	vector<int> &l1_emiso_iphi() { return cms2.l1_emiso_iphi(); }
	vector<int> &l1_emiso_rawId() { return cms2.l1_emiso_rawId(); }
	vector<int> &l1_emiso_type() { return cms2.l1_emiso_type(); }
	vector<int> &l1_emnoiso_ieta() { return cms2.l1_emnoiso_ieta(); }
	vector<int> &l1_emnoiso_iphi() { return cms2.l1_emnoiso_iphi(); }
	vector<int> &l1_emnoiso_rawId() { return cms2.l1_emnoiso_rawId(); }
	vector<int> &l1_emnoiso_type() { return cms2.l1_emnoiso_type(); }
	vector<int> &l1_jetsc_ieta() { return cms2.l1_jetsc_ieta(); }
	vector<int> &l1_jetsc_iphi() { return cms2.l1_jetsc_iphi(); }
	vector<int> &l1_jetsc_rawId() { return cms2.l1_jetsc_rawId(); }
	vector<int> &l1_jetsc_type() { return cms2.l1_jetsc_type(); }
	vector<int> &l1_jetsf_ieta() { return cms2.l1_jetsf_ieta(); }
	vector<int> &l1_jetsf_iphi() { return cms2.l1_jetsf_iphi(); }
	vector<int> &l1_jetsf_rawId() { return cms2.l1_jetsf_rawId(); }
	vector<int> &l1_jetsf_type() { return cms2.l1_jetsf_type(); }
	vector<int> &l1_jetst_ieta() { return cms2.l1_jetst_ieta(); }
	vector<int> &l1_jetst_iphi() { return cms2.l1_jetst_iphi(); }
	vector<int> &l1_jetst_rawId() { return cms2.l1_jetst_rawId(); }
	vector<int> &l1_jetst_type() { return cms2.l1_jetst_type(); }
	vector<int> &l1_mus_flags() { return cms2.l1_mus_flags(); }
	vector<int> &l1_mus_q() { return cms2.l1_mus_q(); }
	vector<int> &l1_mus_qual() { return cms2.l1_mus_qual(); }
	vector<int> &l1_mus_qualFlags() { return cms2.l1_mus_qualFlags(); }
	vector<int> &pxl_ndigis_pxb() { return cms2.pxl_ndigis_pxb(); }
	vector<int> &pxl_ndigis_pxf() { return cms2.pxl_ndigis_pxf(); }
	vector<int> &scs_elsidx() { return cms2.scs_elsidx(); }
	vector<int> &scs_severitySeed() { return cms2.scs_severitySeed(); }
	vector<int> &trks_algo() { return cms2.trks_algo(); }
	vector<int> &trks_charge() { return cms2.trks_charge(); }
	vector<int> &trks_exp_innerlayers() { return cms2.trks_exp_innerlayers(); }
	vector<int> &trks_exp_outerlayers() { return cms2.trks_exp_outerlayers(); }
	vector<int> &trks_layer1_det() { return cms2.trks_layer1_det(); }
	vector<int> &trks_layer1_layer() { return cms2.trks_layer1_layer(); }
	vector<int> &trks_layer1_sizerphi() { return cms2.trks_layer1_sizerphi(); }
	vector<int> &trks_layer1_sizerz() { return cms2.trks_layer1_sizerz(); }
	vector<int> &trks_lostHits() { return cms2.trks_lostHits(); }
	vector<int> &trks_lost_pixelhits() { return cms2.trks_lost_pixelhits(); }
	vector<int> &trks_nlayers() { return cms2.trks_nlayers(); }
	vector<int> &trks_nlayers3D() { return cms2.trks_nlayers3D(); }
	vector<int> &trks_nlayersLost() { return cms2.trks_nlayersLost(); }
	vector<int> &trks_qualityMask() { return cms2.trks_qualityMask(); }
	vector<int> &trks_validHits() { return cms2.trks_validHits(); }
	vector<int> &trks_valid_pixelhits() { return cms2.trks_valid_pixelhits(); }
	vector<int> &vtxs_isFake() { return cms2.vtxs_isFake(); }
	vector<int> &vtxs_isValid() { return cms2.vtxs_isValid(); }
	vector<int> &vtxs_tracksSize() { return cms2.vtxs_tracksSize(); }
	vector<vector<int> > &pseudo_srFlags() { return cms2.pseudo_srFlags(); }
	unsigned int &evt_ntwrs() { return cms2.evt_ntwrs(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	unsigned int &l1_bits1() { return cms2.l1_bits1(); }
	unsigned int &l1_bits2() { return cms2.l1_bits2(); }
	unsigned int &l1_bits3() { return cms2.l1_bits3(); }
	unsigned int &l1_bits4() { return cms2.l1_bits4(); }
	unsigned int &l1_techbits1() { return cms2.l1_techbits1(); }
	unsigned int &l1_techbits2() { return cms2.l1_techbits2(); }
	unsigned int &evt_nscs() { return cms2.evt_nscs(); }
	unsigned int &evt_nvtxs() { return cms2.evt_nvtxs(); }
	vector<unsigned int> &twrs_numBadEcalCells() { return cms2.twrs_numBadEcalCells(); }
	vector<unsigned int> &twrs_numBadHcalCells() { return cms2.twrs_numBadHcalCells(); }
	vector<unsigned int> &twrs_numProblematicEcalCells() { return cms2.twrs_numProblematicEcalCells(); }
	vector<unsigned int> &twrs_numProblematicHcalCells() { return cms2.twrs_numProblematicHcalCells(); }
	vector<unsigned int> &twrs_numRecoveredEcalCells() { return cms2.twrs_numRecoveredEcalCells(); }
	vector<unsigned int> &twrs_numRecoveredHcalCells() { return cms2.twrs_numRecoveredHcalCells(); }
	static bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }
}
#endif
