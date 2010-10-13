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
	vector<TString> evt_HLT_trigNames_;
	TBranch *evt_HLT_trigNames_branch;
	bool evt_HLT_trigNames_isLoaded;
	vector<TString> evt_L1_trigNames_;
	TBranch *evt_L1_trigNames_branch;
	bool evt_L1_trigNames_isLoaded;
	double genps_pthat_;
	TBranch *genps_pthat_branch;
	bool genps_pthat_isLoaded;
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  evt_bsp4_;
	TBranch *evt_bsp4_branch;
	bool evt_bsp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_mc3_motherp4_;
	TBranch *els_mc3_motherp4_branch;
	bool els_mc3_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_mc3_p4_;
	TBranch *els_mc3_p4_branch;
	bool els_mc3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_mc_motherp4_;
	TBranch *els_mc_motherp4_branch;
	bool els_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_mc_p4_;
	TBranch *els_mc_p4_branch;
	bool els_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_mc_gp_p4_;
	TBranch *jets_mc_gp_p4_branch;
	bool jets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_mc_p4_;
	TBranch *jets_mc_p4_branch;
	bool jets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_mc3_motherp4_;
	TBranch *mus_mc3_motherp4_branch;
	bool mus_mc3_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_mc3_p4_;
	TBranch *mus_mc3_p4_branch;
	bool mus_mc3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_mc_motherp4_;
	TBranch *mus_mc_motherp4_branch;
	bool mus_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_mc_p4_;
	TBranch *mus_mc_p4_branch;
	bool mus_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trk_mc3p4_;
	TBranch *trk_mc3p4_branch;
	bool trk_mc3p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trk_mcp4_;
	TBranch *trk_mcp4_branch;
	bool trk_mcp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_p4In_;
	TBranch *els_p4In_branch;
	bool els_p4In_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_p4Out_;
	TBranch *els_p4Out_branch;
	bool els_p4Out_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_trk_p4_;
	TBranch *els_trk_p4_branch;
	bool els_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_vertex_p4_;
	TBranch *els_vertex_p4_branch;
	bool els_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_lepdaughter_p4_;
	TBranch *genps_lepdaughter_p4_branch;
	bool genps_lepdaughter_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_prod_vtx_;
	TBranch *genps_prod_vtx_branch;
	bool genps_prod_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_gfit_outerPos_p4_;
	TBranch *mus_gfit_outerPos_p4_branch;
	bool mus_gfit_outerPos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_trk_p4_;
	TBranch *mus_trk_p4_branch;
	bool mus_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_vertex_p4_;
	TBranch *mus_vertex_p4_branch;
	bool mus_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trks_trk_p4_;
	TBranch *trks_trk_p4_branch;
	bool trks_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trks_vertex_p4_;
	TBranch *trks_vertex_p4_branch;
	bool trks_vertex_p4_isLoaded;
	vector<double> els_inner_positionx_;
	TBranch *els_inner_positionx_branch;
	bool els_inner_positionx_isLoaded;
	vector<double> els_inner_positiony_;
	TBranch *els_inner_positiony_branch;
	bool els_inner_positiony_isLoaded;
	vector<double> els_inner_positionz_;
	TBranch *els_inner_positionz_branch;
	bool els_inner_positionz_isLoaded;
	vector<double> trks_inner_positionx_;
	TBranch *trks_inner_positionx_branch;
	bool trks_inner_positionx_isLoaded;
	vector<double> trks_inner_positiony_;
	TBranch *trks_inner_positiony_branch;
	bool trks_inner_positiony_isLoaded;
	vector<double> trks_inner_positionz_;
	TBranch *trks_inner_positionz_branch;
	bool trks_inner_positionz_isLoaded;
	float evt_bs_dxdz_;
	TBranch *evt_bs_dxdz_branch;
	bool evt_bs_dxdz_isLoaded;
	float evt_bs_dxdzErr_;
	TBranch *evt_bs_dxdzErr_branch;
	bool evt_bs_dxdzErr_isLoaded;
	float evt_bs_dydz_;
	TBranch *evt_bs_dydz_branch;
	bool evt_bs_dydz_isLoaded;
	float evt_bs_dydzErr_;
	TBranch *evt_bs_dydzErr_branch;
	bool evt_bs_dydzErr_isLoaded;
	float evt_bs_sigmaZ_;
	TBranch *evt_bs_sigmaZ_branch;
	bool evt_bs_sigmaZ_isLoaded;
	float evt_bs_sigmaZErr_;
	TBranch *evt_bs_sigmaZErr_branch;
	bool evt_bs_sigmaZErr_isLoaded;
	float evt_bs_widthx_;
	TBranch *evt_bs_widthx_branch;
	bool evt_bs_widthx_isLoaded;
	float evt_bs_widthxErr_;
	TBranch *evt_bs_widthxErr_branch;
	bool evt_bs_widthxErr_isLoaded;
	float evt_bs_widthy_;
	TBranch *evt_bs_widthy_branch;
	bool evt_bs_widthy_isLoaded;
	float evt_bs_widthyErr_;
	TBranch *evt_bs_widthyErr_branch;
	bool evt_bs_widthyErr_isLoaded;
	float evt_bs_xErr_;
	TBranch *evt_bs_xErr_branch;
	bool evt_bs_xErr_isLoaded;
	float evt_bs_yErr_;
	TBranch *evt_bs_yErr_branch;
	bool evt_bs_yErr_isLoaded;
	float evt_bs_zErr_;
	TBranch *evt_bs_zErr_branch;
	bool evt_bs_zErr_isLoaded;
	float evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float evt_weight_;
	TBranch *evt_weight_branch;
	bool evt_weight_isLoaded;
	float evt_xsec_excl_;
	TBranch *evt_xsec_excl_branch;
	bool evt_xsec_excl_isLoaded;
	float evt_xsec_incl_;
	TBranch *evt_xsec_incl_branch;
	bool evt_xsec_incl_isLoaded;
	float gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float gen_metPhi_;
	TBranch *gen_metPhi_branch;
	bool gen_metPhi_isLoaded;
	vector<float> els_mc3dr_;
	TBranch *els_mc3dr_branch;
	bool els_mc3dr_isLoaded;
	vector<float> els_mcdr_;
	TBranch *els_mcdr_branch;
	bool els_mcdr_isLoaded;
	vector<float> jets_mcdr_;
	TBranch *jets_mcdr_branch;
	bool jets_mcdr_isLoaded;
	vector<float> jets_mc_emEnergy_;
	TBranch *jets_mc_emEnergy_branch;
	bool jets_mc_emEnergy_isLoaded;
	vector<float> jets_mc_gpdr_;
	TBranch *jets_mc_gpdr_branch;
	bool jets_mc_gpdr_isLoaded;
	vector<float> jets_mc_hadEnergy_;
	TBranch *jets_mc_hadEnergy_branch;
	bool jets_mc_hadEnergy_isLoaded;
	vector<float> jets_mc_invEnergy_;
	TBranch *jets_mc_invEnergy_branch;
	bool jets_mc_invEnergy_isLoaded;
	vector<float> jets_mc_otherEnergy_;
	TBranch *jets_mc_otherEnergy_branch;
	bool jets_mc_otherEnergy_isLoaded;
	vector<float> mus_mc3dr_;
	TBranch *mus_mc3dr_branch;
	bool mus_mc3dr_isLoaded;
	vector<float> mus_mcdr_;
	TBranch *mus_mcdr_branch;
	bool mus_mcdr_isLoaded;
	vector<float> trk_mc3dr_;
	TBranch *trk_mc3dr_branch;
	bool trk_mc3dr_isLoaded;
	vector<float> trk_mcdr_;
	TBranch *trk_mcdr_branch;
	bool trk_mcdr_isLoaded;
	vector<float> els_conv_dcot_;
	TBranch *els_conv_dcot_branch;
	bool els_conv_dcot_isLoaded;
	vector<float> els_conv_dist_;
	TBranch *els_conv_dist_branch;
	bool els_conv_dist_isLoaded;
	vector<float> trks_conv_dcot_;
	TBranch *trks_conv_dcot_branch;
	bool trks_conv_dcot_isLoaded;
	vector<float> trks_conv_dist_;
	TBranch *trks_conv_dist_branch;
	bool trks_conv_dist_isLoaded;
	vector<float> els_musdr_;
	TBranch *els_musdr_branch;
	bool els_musdr_isLoaded;
	vector<float> els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float> els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
	vector<float> els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float> els_d0_;
	TBranch *els_d0_branch;
	bool els_d0_isLoaded;
	vector<float> els_d0Err_;
	TBranch *els_d0Err_branch;
	bool els_d0Err_isLoaded;
	vector<float> els_d0corr_;
	TBranch *els_d0corr_branch;
	bool els_d0corr_isLoaded;
	vector<float> els_dEtaIn_;
	TBranch *els_dEtaIn_branch;
	bool els_dEtaIn_isLoaded;
	vector<float> els_dEtaOut_;
	TBranch *els_dEtaOut_branch;
	bool els_dEtaOut_isLoaded;
	vector<float> els_dPhiIn_;
	TBranch *els_dPhiIn_branch;
	bool els_dPhiIn_isLoaded;
	vector<float> els_dPhiInPhiOut_;
	TBranch *els_dPhiInPhiOut_branch;
	bool els_dPhiInPhiOut_isLoaded;
	vector<float> els_dPhiOut_;
	TBranch *els_dPhiOut_branch;
	bool els_dPhiOut_isLoaded;
	vector<float> els_e1x5_;
	TBranch *els_e1x5_branch;
	bool els_e1x5_isLoaded;
	vector<float> els_e2x5Max_;
	TBranch *els_e2x5Max_branch;
	bool els_e2x5Max_isLoaded;
	vector<float> els_e3x3_;
	TBranch *els_e3x3_branch;
	bool els_e3x3_isLoaded;
	vector<float> els_e5x5_;
	TBranch *els_e5x5_branch;
	bool els_e5x5_isLoaded;
	vector<float> els_eMax_;
	TBranch *els_eMax_branch;
	bool els_eMax_isLoaded;
	vector<float> els_eOverPIn_;
	TBranch *els_eOverPIn_branch;
	bool els_eOverPIn_isLoaded;
	vector<float> els_eSC_;
	TBranch *els_eSC_branch;
	bool els_eSC_isLoaded;
	vector<float> els_eSCPresh_;
	TBranch *els_eSCPresh_branch;
	bool els_eSCPresh_isLoaded;
	vector<float> els_eSCRaw_;
	TBranch *els_eSCRaw_branch;
	bool els_eSCRaw_isLoaded;
	vector<float> els_eSeed_;
	TBranch *els_eSeed_branch;
	bool els_eSeed_isLoaded;
	vector<float> els_eSeedOverPIn_;
	TBranch *els_eSeedOverPIn_branch;
	bool els_eSeedOverPIn_isLoaded;
	vector<float> els_eSeedOverPOut_;
	TBranch *els_eSeedOverPOut_branch;
	bool els_eSeedOverPOut_isLoaded;
	vector<float> els_ecalIso03_;
	TBranch *els_ecalIso03_branch;
	bool els_ecalIso03_isLoaded;
	vector<float> els_ecalIso04_;
	TBranch *els_ecalIso04_branch;
	bool els_ecalIso04_isLoaded;
	vector<float> els_egamma_looseId_;
	TBranch *els_egamma_looseId_branch;
	bool els_egamma_looseId_isLoaded;
	vector<float> els_egamma_robustHighEnergy_;
	TBranch *els_egamma_robustHighEnergy_branch;
	bool els_egamma_robustHighEnergy_isLoaded;
	vector<float> els_egamma_robustLooseId_;
	TBranch *els_egamma_robustLooseId_branch;
	bool els_egamma_robustLooseId_isLoaded;
	vector<float> els_egamma_robustTightId_;
	TBranch *els_egamma_robustTightId_branch;
	bool els_egamma_robustTightId_isLoaded;
	vector<float> els_egamma_tightId_;
	TBranch *els_egamma_tightId_branch;
	bool els_egamma_tightId_isLoaded;
	vector<float> els_etaErr_;
	TBranch *els_etaErr_branch;
	bool els_etaErr_isLoaded;
	vector<float> els_etaSC_;
	TBranch *els_etaSC_branch;
	bool els_etaSC_isLoaded;
	vector<float> els_fBrem_;
	TBranch *els_fBrem_branch;
	bool els_fBrem_isLoaded;
	vector<float> els_hOverE_;
	TBranch *els_hOverE_branch;
	bool els_hOverE_isLoaded;
	vector<float> els_hcalIso03_;
	TBranch *els_hcalIso03_branch;
	bool els_hcalIso03_isLoaded;
	vector<float> els_hcalIso04_;
	TBranch *els_hcalIso04_branch;
	bool els_hcalIso04_isLoaded;
	vector<float> els_layer1_charge_;
	TBranch *els_layer1_charge_branch;
	bool els_layer1_charge_isLoaded;
	vector<float> els_ndof_;
	TBranch *els_ndof_branch;
	bool els_ndof_isLoaded;
	vector<float> els_outerEta_;
	TBranch *els_outerEta_branch;
	bool els_outerEta_isLoaded;
	vector<float> els_outerPhi_;
	TBranch *els_outerPhi_branch;
	bool els_outerPhi_isLoaded;
	vector<float> els_phiErr_;
	TBranch *els_phiErr_branch;
	bool els_phiErr_isLoaded;
	vector<float> els_phiSC_;
	TBranch *els_phiSC_branch;
	bool els_phiSC_isLoaded;
	vector<float> els_pin_;
	TBranch *els_pin_branch;
	bool els_pin_isLoaded;
	vector<float> els_pout_;
	TBranch *els_pout_branch;
	bool els_pout_isLoaded;
	vector<float> els_ptErr_;
	TBranch *els_ptErr_branch;
	bool els_ptErr_isLoaded;
	vector<float> els_sigmaEtaEta_;
	TBranch *els_sigmaEtaEta_branch;
	bool els_sigmaEtaEta_isLoaded;
	vector<float> els_sigmaIEtaIEta_;
	TBranch *els_sigmaIEtaIEta_branch;
	bool els_sigmaIEtaIEta_isLoaded;
	vector<float> els_sigmaIPhiIPhi_;
	TBranch *els_sigmaIPhiIPhi_branch;
	bool els_sigmaIPhiIPhi_isLoaded;
	vector<float> els_sigmaPhiPhi_;
	TBranch *els_sigmaPhiPhi_branch;
	bool els_sigmaPhiPhi_isLoaded;
	vector<float> els_tkIso03_;
	TBranch *els_tkIso03_branch;
	bool els_tkIso03_isLoaded;
	vector<float> els_tkIso04_;
	TBranch *els_tkIso04_branch;
	bool els_tkIso04_isLoaded;
	vector<float> els_vertexphi_;
	TBranch *els_vertexphi_branch;
	bool els_vertexphi_isLoaded;
	vector<float> els_z0_;
	TBranch *els_z0_branch;
	bool els_z0_isLoaded;
	vector<float> els_z0Err_;
	TBranch *els_z0Err_branch;
	bool els_z0Err_isLoaded;
	vector<float> els_z0corr_;
	TBranch *els_z0corr_branch;
	bool els_z0corr_isLoaded;
	vector<float> mus_caloCompatibility_;
	TBranch *mus_caloCompatibility_branch;
	bool mus_caloCompatibility_isLoaded;
	vector<float> mus_chi2_;
	TBranch *mus_chi2_branch;
	bool mus_chi2_isLoaded;
	vector<float> mus_d0_;
	TBranch *mus_d0_branch;
	bool mus_d0_isLoaded;
	vector<float> mus_d0Err_;
	TBranch *mus_d0Err_branch;
	bool mus_d0Err_isLoaded;
	vector<float> mus_d0corr_;
	TBranch *mus_d0corr_branch;
	bool mus_d0corr_isLoaded;
	vector<float> mus_e_em_;
	TBranch *mus_e_em_branch;
	bool mus_e_em_isLoaded;
	vector<float> mus_e_emS9_;
	TBranch *mus_e_emS9_branch;
	bool mus_e_emS9_isLoaded;
	vector<float> mus_e_had_;
	TBranch *mus_e_had_branch;
	bool mus_e_had_isLoaded;
	vector<float> mus_e_hadS9_;
	TBranch *mus_e_hadS9_branch;
	bool mus_e_hadS9_isLoaded;
	vector<float> mus_e_ho_;
	TBranch *mus_e_ho_branch;
	bool mus_e_ho_isLoaded;
	vector<float> mus_e_hoS9_;
	TBranch *mus_e_hoS9_branch;
	bool mus_e_hoS9_isLoaded;
	vector<float> mus_etaErr_;
	TBranch *mus_etaErr_branch;
	bool mus_etaErr_isLoaded;
	vector<float> mus_gfit_chi2_;
	TBranch *mus_gfit_chi2_branch;
	bool mus_gfit_chi2_isLoaded;
	vector<float> mus_gfit_ndof_;
	TBranch *mus_gfit_ndof_branch;
	bool mus_gfit_ndof_isLoaded;
	vector<float> mus_iso_;
	TBranch *mus_iso_branch;
	bool mus_iso_isLoaded;
	vector<float> mus_iso03_emEt_;
	TBranch *mus_iso03_emEt_branch;
	bool mus_iso03_emEt_isLoaded;
	vector<float> mus_iso03_hadEt_;
	TBranch *mus_iso03_hadEt_branch;
	bool mus_iso03_hadEt_isLoaded;
	vector<float> mus_iso03_hoEt_;
	TBranch *mus_iso03_hoEt_branch;
	bool mus_iso03_hoEt_isLoaded;
	vector<float> mus_iso03_sumPt_;
	TBranch *mus_iso03_sumPt_branch;
	bool mus_iso03_sumPt_isLoaded;
	vector<float> mus_iso05_emEt_;
	TBranch *mus_iso05_emEt_branch;
	bool mus_iso05_emEt_isLoaded;
	vector<float> mus_iso05_hadEt_;
	TBranch *mus_iso05_hadEt_branch;
	bool mus_iso05_hadEt_isLoaded;
	vector<float> mus_iso05_hoEt_;
	TBranch *mus_iso05_hoEt_branch;
	bool mus_iso05_hoEt_isLoaded;
	vector<float> mus_iso05_sumPt_;
	TBranch *mus_iso05_sumPt_branch;
	bool mus_iso05_sumPt_isLoaded;
	vector<float> mus_ndof_;
	TBranch *mus_ndof_branch;
	bool mus_ndof_isLoaded;
	vector<float> mus_outerEta_;
	TBranch *mus_outerEta_branch;
	bool mus_outerEta_isLoaded;
	vector<float> mus_outerPhi_;
	TBranch *mus_outerPhi_branch;
	bool mus_outerPhi_isLoaded;
	vector<float> mus_phiErr_;
	TBranch *mus_phiErr_branch;
	bool mus_phiErr_isLoaded;
	vector<float> mus_ptErr_;
	TBranch *mus_ptErr_branch;
	bool mus_ptErr_isLoaded;
	vector<float> mus_qoverp_;
	TBranch *mus_qoverp_branch;
	bool mus_qoverp_isLoaded;
	vector<float> mus_qoverpError_;
	TBranch *mus_qoverpError_branch;
	bool mus_qoverpError_isLoaded;
	vector<float> mus_vertexphi_;
	TBranch *mus_vertexphi_branch;
	bool mus_vertexphi_isLoaded;
	vector<float> mus_z0_;
	TBranch *mus_z0_branch;
	bool mus_z0_isLoaded;
	vector<float> mus_z0Err_;
	TBranch *mus_z0Err_branch;
	bool mus_z0Err_isLoaded;
	vector<float> mus_z0corr_;
	TBranch *mus_z0corr_branch;
	bool mus_z0corr_isLoaded;
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
	vector<float> trks_outerEt_;
	TBranch *trks_outerEt_branch;
	bool trks_outerEt_isLoaded;
	vector<float> trks_outerEta_;
	TBranch *trks_outerEta_branch;
	bool trks_outerEta_isLoaded;
	vector<float> trks_outerPhi_;
	TBranch *trks_outerPhi_branch;
	bool trks_outerPhi_isLoaded;
	vector<float> trks_phiErr_;
	TBranch *trks_phiErr_branch;
	bool trks_phiErr_isLoaded;
	vector<float> trks_ptErr_;
	TBranch *trks_ptErr_branch;
	bool trks_ptErr_isLoaded;
	vector<float> trks_tkIso_;
	TBranch *trks_tkIso_branch;
	bool trks_tkIso_isLoaded;
	vector<float> trks_vertexphi_;
	TBranch *trks_vertexphi_branch;
	bool trks_vertexphi_isLoaded;
	vector<float> trks_z0_;
	TBranch *trks_z0_branch;
	bool trks_z0_isLoaded;
	vector<float> trks_z0Err_;
	TBranch *trks_z0Err_branch;
	bool trks_z0Err_isLoaded;
	vector<float> trks_z0corr_;
	TBranch *trks_z0corr_branch;
	bool trks_z0corr_isLoaded;
	vector<float> trks_elsdr_;
	TBranch *trks_elsdr_branch;
	bool trks_elsdr_isLoaded;
	vector<float> trks_elsshFrac_;
	TBranch *trks_elsshFrac_branch;
	bool trks_elsshFrac_isLoaded;
	vector<float> trk_musdr_;
	TBranch *trk_musdr_branch;
	bool trk_musdr_isLoaded;
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
	vector<vector<float> > trks_residualX_;
	TBranch *trks_residualX_branch;
	bool trks_residualX_isLoaded;
	vector<vector<float> > trks_residualY_;
	TBranch *trks_residualY_branch;
	bool trks_residualY_isLoaded;
	int evt_HLT1_;
	TBranch *evt_HLT1_branch;
	bool evt_HLT1_isLoaded;
	int evt_HLT2_;
	TBranch *evt_HLT2_branch;
	bool evt_HLT2_isLoaded;
	int evt_HLT3_;
	TBranch *evt_HLT3_branch;
	bool evt_HLT3_isLoaded;
	int evt_HLT4_;
	TBranch *evt_HLT4_branch;
	bool evt_HLT4_isLoaded;
	int evt_HLT5_;
	TBranch *evt_HLT5_branch;
	bool evt_HLT5_isLoaded;
	int evt_HLT6_;
	TBranch *evt_HLT6_branch;
	bool evt_HLT6_isLoaded;
	int evt_HLT7_;
	TBranch *evt_HLT7_branch;
	bool evt_HLT7_isLoaded;
	int evt_HLT8_;
	TBranch *evt_HLT8_branch;
	bool evt_HLT8_isLoaded;
	int evt_L1_1_;
	TBranch *evt_L1_1_branch;
	bool evt_L1_1_isLoaded;
	int evt_L1_2_;
	TBranch *evt_L1_2_branch;
	bool evt_L1_2_isLoaded;
	int evt_L1_3_;
	TBranch *evt_L1_3_branch;
	bool evt_L1_3_isLoaded;
	int evt_L1_4_;
	TBranch *evt_L1_4_branch;
	bool evt_L1_4_isLoaded;
	vector<int> els_mc3_id_;
	TBranch *els_mc3_id_branch;
	bool els_mc3_id_isLoaded;
	vector<int> els_mc3idx_;
	TBranch *els_mc3idx_branch;
	bool els_mc3idx_isLoaded;
	vector<int> els_mc3_motherid_;
	TBranch *els_mc3_motherid_branch;
	bool els_mc3_motherid_isLoaded;
	vector<int> els_mc_id_;
	TBranch *els_mc_id_branch;
	bool els_mc_id_isLoaded;
	vector<int> els_mcidx_;
	TBranch *els_mcidx_branch;
	bool els_mcidx_isLoaded;
	vector<int> els_mc_motherid_;
	TBranch *els_mc_motherid_branch;
	bool els_mc_motherid_isLoaded;
	vector<int> jets_mc_id_;
	TBranch *jets_mc_id_branch;
	bool jets_mc_id_isLoaded;
	vector<int> mus_mc3_id_;
	TBranch *mus_mc3_id_branch;
	bool mus_mc3_id_isLoaded;
	vector<int> mus_mc3idx_;
	TBranch *mus_mc3idx_branch;
	bool mus_mc3idx_isLoaded;
	vector<int> mus_mc3_motherid_;
	TBranch *mus_mc3_motherid_branch;
	bool mus_mc3_motherid_isLoaded;
	vector<int> mus_mc_id_;
	TBranch *mus_mc_id_branch;
	bool mus_mc_id_isLoaded;
	vector<int> mus_mcidx_;
	TBranch *mus_mcidx_branch;
	bool mus_mcidx_isLoaded;
	vector<int> mus_mc_motherid_;
	TBranch *mus_mc_motherid_branch;
	bool mus_mc_motherid_isLoaded;
	vector<int> trk_mc3_id_;
	TBranch *trk_mc3_id_branch;
	bool trk_mc3_id_isLoaded;
	vector<int> trk_mc3idx_;
	TBranch *trk_mc3idx_branch;
	bool trk_mc3idx_isLoaded;
	vector<int> trk_mc3_motherid_;
	TBranch *trk_mc3_motherid_branch;
	bool trk_mc3_motherid_isLoaded;
	vector<int> trk_mc_id_;
	TBranch *trk_mc_id_branch;
	bool trk_mc_id_isLoaded;
	vector<int> trk_mcidx_;
	TBranch *trk_mcidx_branch;
	bool trk_mcidx_isLoaded;
	vector<int> trk_mc_motherid_;
	TBranch *trk_mc_motherid_branch;
	bool trk_mc_motherid_isLoaded;
	vector<int> els_conv_tkidx_;
	TBranch *els_conv_tkidx_branch;
	bool els_conv_tkidx_isLoaded;
	vector<int> trks_conv_tkidx_;
	TBranch *trks_conv_tkidx_branch;
	bool trks_conv_tkidx_isLoaded;
	vector<int> els_closestMuon_;
	TBranch *els_closestMuon_branch;
	bool els_closestMuon_isLoaded;
	vector<int> els_trkidx_;
	TBranch *els_trkidx_branch;
	bool els_trkidx_isLoaded;
	vector<int> els_category_;
	TBranch *els_category_branch;
	bool els_category_isLoaded;
	vector<int> els_categoryold_;
	TBranch *els_categoryold_branch;
	bool els_categoryold_isLoaded;
	vector<int> els_charge_;
	TBranch *els_charge_branch;
	bool els_charge_isLoaded;
	vector<int> els_class_;
	TBranch *els_class_branch;
	bool els_class_isLoaded;
	vector<int> els_layer1_det_;
	TBranch *els_layer1_det_branch;
	bool els_layer1_det_isLoaded;
	vector<int> els_layer1_layer_;
	TBranch *els_layer1_layer_branch;
	bool els_layer1_layer_isLoaded;
	vector<int> els_layer1_sizerphi_;
	TBranch *els_layer1_sizerphi_branch;
	bool els_layer1_sizerphi_isLoaded;
	vector<int> els_layer1_sizerz_;
	TBranch *els_layer1_sizerz_branch;
	bool els_layer1_sizerz_isLoaded;
	vector<int> els_looseId_;
	TBranch *els_looseId_branch;
	bool els_looseId_isLoaded;
	vector<int> els_lostHits_;
	TBranch *els_lostHits_branch;
	bool els_lostHits_isLoaded;
	vector<int> els_lost_pixelhits_;
	TBranch *els_lost_pixelhits_branch;
	bool els_lost_pixelhits_isLoaded;
	vector<int> els_nSeed_;
	TBranch *els_nSeed_branch;
	bool els_nSeed_isLoaded;
	vector<int> els_n_inner_layers_;
	TBranch *els_n_inner_layers_branch;
	bool els_n_inner_layers_isLoaded;
	vector<int> els_n_outer_layers_;
	TBranch *els_n_outer_layers_branch;
	bool els_n_outer_layers_isLoaded;
	vector<int> els_pass3looseId_;
	TBranch *els_pass3looseId_branch;
	bool els_pass3looseId_isLoaded;
	vector<int> els_pass3simpleId_;
	TBranch *els_pass3simpleId_branch;
	bool els_pass3simpleId_isLoaded;
	vector<int> els_pass3tightId_;
	TBranch *els_pass3tightId_branch;
	bool els_pass3tightId_isLoaded;
	vector<int> els_robustId_;
	TBranch *els_robustId_branch;
	bool els_robustId_isLoaded;
	vector<int> els_simpleIdPlus_;
	TBranch *els_simpleIdPlus_branch;
	bool els_simpleIdPlus_isLoaded;
	vector<int> els_tightId_;
	TBranch *els_tightId_branch;
	bool els_tightId_isLoaded;
	vector<int> els_tightId22XMaxMatteo_;
	TBranch *els_tightId22XMaxMatteo_branch;
	bool els_tightId22XMaxMatteo_isLoaded;
	vector<int> els_tightId22XMinMatteo_;
	TBranch *els_tightId22XMinMatteo_branch;
	bool els_tightId22XMinMatteo_isLoaded;
	vector<int> els_validHits_;
	TBranch *els_validHits_branch;
	bool els_validHits_isLoaded;
	vector<int> els_valid_pixelhits_;
	TBranch *els_valid_pixelhits_branch;
	bool els_valid_pixelhits_isLoaded;
	vector<int> genps_id_;
	TBranch *genps_id_branch;
	bool genps_id_isLoaded;
	vector<int> genps_id_mother_;
	TBranch *genps_id_mother_branch;
	bool genps_id_mother_isLoaded;
	vector<int> genps_lepdaughter_id_;
	TBranch *genps_lepdaughter_id_branch;
	bool genps_lepdaughter_id_isLoaded;
	vector<int> genps_lepdaughter_idx_;
	TBranch *genps_lepdaughter_idx_branch;
	bool genps_lepdaughter_idx_isLoaded;
	vector<int> genps_status_;
	TBranch *genps_status_branch;
	bool genps_status_isLoaded;
	vector<int> mus_charge_;
	TBranch *mus_charge_branch;
	bool mus_charge_isLoaded;
	vector<int> mus_gfit_validHits_;
	TBranch *mus_gfit_validHits_branch;
	bool mus_gfit_validHits_isLoaded;
	vector<int> mus_goodmask_;
	TBranch *mus_goodmask_branch;
	bool mus_goodmask_isLoaded;
	vector<int> mus_iso03_ntrk_;
	TBranch *mus_iso03_ntrk_branch;
	bool mus_iso03_ntrk_isLoaded;
	vector<int> mus_iso05_ntrk_;
	TBranch *mus_iso05_ntrk_branch;
	bool mus_iso05_ntrk_isLoaded;
	vector<int> mus_lostHits_;
	TBranch *mus_lostHits_branch;
	bool mus_lostHits_isLoaded;
	vector<int> mus_nmatches_;
	TBranch *mus_nmatches_branch;
	bool mus_nmatches_isLoaded;
	vector<int> mus_pid_TM2DCompatibilityLoose_;
	TBranch *mus_pid_TM2DCompatibilityLoose_branch;
	bool mus_pid_TM2DCompatibilityLoose_isLoaded;
	vector<int> mus_pid_TM2DCompatibilityTight_;
	TBranch *mus_pid_TM2DCompatibilityTight_branch;
	bool mus_pid_TM2DCompatibilityTight_isLoaded;
	vector<int> mus_pid_TMLastStationLoose_;
	TBranch *mus_pid_TMLastStationLoose_branch;
	bool mus_pid_TMLastStationLoose_isLoaded;
	vector<int> mus_pid_TMLastStationTight_;
	TBranch *mus_pid_TMLastStationTight_branch;
	bool mus_pid_TMLastStationTight_isLoaded;
	vector<int> mus_trk_charge_;
	TBranch *mus_trk_charge_branch;
	bool mus_trk_charge_isLoaded;
	vector<int> mus_trkrefkey_;
	TBranch *mus_trkrefkey_branch;
	bool mus_trkrefkey_isLoaded;
	vector<int> mus_type_;
	TBranch *mus_type_branch;
	bool mus_type_isLoaded;
	vector<int> mus_validHits_;
	TBranch *mus_validHits_branch;
	bool mus_validHits_isLoaded;
	vector<int> trks_charge_;
	TBranch *trks_charge_branch;
	bool trks_charge_isLoaded;
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
	vector<int> trks_n_inner_layers_;
	TBranch *trks_n_inner_layers_branch;
	bool trks_n_inner_layers_isLoaded;
	vector<int> trks_n_outer_layers_;
	TBranch *trks_n_outer_layers_branch;
	bool trks_n_outer_layers_isLoaded;
	vector<int> trks_validHits_;
	TBranch *trks_validHits_branch;
	bool trks_validHits_isLoaded;
	vector<int> trks_valid_pixelhits_;
	TBranch *trks_valid_pixelhits_branch;
	bool trks_valid_pixelhits_isLoaded;
	vector<int> trks_elsidx_;
	TBranch *trks_elsidx_branch;
	bool trks_elsidx_isLoaded;
	vector<int> trk_musidx_;
	TBranch *trk_musidx_branch;
	bool trk_musidx_isLoaded;
	vector<int> vtxs_hasRefittedTracks_;
	TBranch *vtxs_hasRefittedTracks_branch;
	bool vtxs_hasRefittedTracks_isLoaded;
	vector<int> vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int> vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<int> vtxs_tracksSize_;
	TBranch *vtxs_tracksSize_branch;
	bool vtxs_tracksSize_isLoaded;
	vector<vector<int> > trks_hit_type_;
	TBranch *trks_hit_type_branch;
	bool trks_hit_type_isLoaded;
	unsigned int evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
public: 
int ScanChain(class TChain* chain, int nEvents=-1, std::string skimFilePrefix="");
void Init(TTree *tree) {
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		vtxs_position_branch->SetAddress(&vtxs_position_);
	}
	if(vtxs_position_branch == 0 ) {
	cout << "Branch vtxs_position does not exist." << endl;
	}
	evt_bsp4_branch = 0;
	if (tree->GetAlias("evt_bsp4") != 0) {
		evt_bsp4_branch = tree->GetBranch(tree->GetAlias("evt_bsp4"));
		evt_bsp4_branch->SetAddress(&evt_bsp4_);
	}
	if(evt_bsp4_branch == 0 ) {
	cout << "Branch evt_bsp4 does not exist." << endl;
	}
	els_mc3_motherp4_branch = 0;
	if (tree->GetAlias("els_mc3_motherp4") != 0) {
		els_mc3_motherp4_branch = tree->GetBranch(tree->GetAlias("els_mc3_motherp4"));
		els_mc3_motherp4_branch->SetAddress(&els_mc3_motherp4_);
	}
	if(els_mc3_motherp4_branch == 0 ) {
	cout << "Branch els_mc3_motherp4 does not exist." << endl;
	}
	els_mc3_p4_branch = 0;
	if (tree->GetAlias("els_mc3_p4") != 0) {
		els_mc3_p4_branch = tree->GetBranch(tree->GetAlias("els_mc3_p4"));
		els_mc3_p4_branch->SetAddress(&els_mc3_p4_);
	}
	if(els_mc3_p4_branch == 0 ) {
	cout << "Branch els_mc3_p4 does not exist." << endl;
	}
	els_mc_motherp4_branch = 0;
	if (tree->GetAlias("els_mc_motherp4") != 0) {
		els_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("els_mc_motherp4"));
		els_mc_motherp4_branch->SetAddress(&els_mc_motherp4_);
	}
	if(els_mc_motherp4_branch == 0 ) {
	cout << "Branch els_mc_motherp4 does not exist." << endl;
	}
	els_mc_p4_branch = 0;
	if (tree->GetAlias("els_mc_p4") != 0) {
		els_mc_p4_branch = tree->GetBranch(tree->GetAlias("els_mc_p4"));
		els_mc_p4_branch->SetAddress(&els_mc_p4_);
	}
	if(els_mc_p4_branch == 0 ) {
	cout << "Branch els_mc_p4 does not exist." << endl;
	}
	jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("jets_mc_gp_p4") != 0) {
		jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("jets_mc_gp_p4"));
		jets_mc_gp_p4_branch->SetAddress(&jets_mc_gp_p4_);
	}
	if(jets_mc_gp_p4_branch == 0 ) {
	cout << "Branch jets_mc_gp_p4 does not exist." << endl;
	}
	jets_mc_p4_branch = 0;
	if (tree->GetAlias("jets_mc_p4") != 0) {
		jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("jets_mc_p4"));
		jets_mc_p4_branch->SetAddress(&jets_mc_p4_);
	}
	if(jets_mc_p4_branch == 0 ) {
	cout << "Branch jets_mc_p4 does not exist." << endl;
	}
	mus_mc3_motherp4_branch = 0;
	if (tree->GetAlias("mus_mc3_motherp4") != 0) {
		mus_mc3_motherp4_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motherp4"));
		mus_mc3_motherp4_branch->SetAddress(&mus_mc3_motherp4_);
	}
	if(mus_mc3_motherp4_branch == 0 ) {
	cout << "Branch mus_mc3_motherp4 does not exist." << endl;
	}
	mus_mc3_p4_branch = 0;
	if (tree->GetAlias("mus_mc3_p4") != 0) {
		mus_mc3_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc3_p4"));
		mus_mc3_p4_branch->SetAddress(&mus_mc3_p4_);
	}
	if(mus_mc3_p4_branch == 0 ) {
	cout << "Branch mus_mc3_p4 does not exist." << endl;
	}
	mus_mc_motherp4_branch = 0;
	if (tree->GetAlias("mus_mc_motherp4") != 0) {
		mus_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherp4"));
		mus_mc_motherp4_branch->SetAddress(&mus_mc_motherp4_);
	}
	if(mus_mc_motherp4_branch == 0 ) {
	cout << "Branch mus_mc_motherp4 does not exist." << endl;
	}
	mus_mc_p4_branch = 0;
	if (tree->GetAlias("mus_mc_p4") != 0) {
		mus_mc_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc_p4"));
		mus_mc_p4_branch->SetAddress(&mus_mc_p4_);
	}
	if(mus_mc_p4_branch == 0 ) {
	cout << "Branch mus_mc_p4 does not exist." << endl;
	}
	trk_mc3p4_branch = 0;
	if (tree->GetAlias("trk_mc3p4") != 0) {
		trk_mc3p4_branch = tree->GetBranch(tree->GetAlias("trk_mc3p4"));
		trk_mc3p4_branch->SetAddress(&trk_mc3p4_);
	}
	if(trk_mc3p4_branch == 0 ) {
	cout << "Branch trk_mc3p4 does not exist." << endl;
	}
	trk_mcp4_branch = 0;
	if (tree->GetAlias("trk_mcp4") != 0) {
		trk_mcp4_branch = tree->GetBranch(tree->GetAlias("trk_mcp4"));
		trk_mcp4_branch->SetAddress(&trk_mcp4_);
	}
	if(trk_mcp4_branch == 0 ) {
	cout << "Branch trk_mcp4 does not exist." << endl;
	}
	els_p4_branch = 0;
	if (tree->GetAlias("els_p4") != 0) {
		els_p4_branch = tree->GetBranch(tree->GetAlias("els_p4"));
		els_p4_branch->SetAddress(&els_p4_);
	}
	if(els_p4_branch == 0 ) {
	cout << "Branch els_p4 does not exist." << endl;
	}
	els_p4In_branch = 0;
	if (tree->GetAlias("els_p4In") != 0) {
		els_p4In_branch = tree->GetBranch(tree->GetAlias("els_p4In"));
		els_p4In_branch->SetAddress(&els_p4In_);
	}
	if(els_p4In_branch == 0 ) {
	cout << "Branch els_p4In does not exist." << endl;
	}
	els_p4Out_branch = 0;
	if (tree->GetAlias("els_p4Out") != 0) {
		els_p4Out_branch = tree->GetBranch(tree->GetAlias("els_p4Out"));
		els_p4Out_branch->SetAddress(&els_p4Out_);
	}
	if(els_p4Out_branch == 0 ) {
	cout << "Branch els_p4Out does not exist." << endl;
	}
	els_trk_p4_branch = 0;
	if (tree->GetAlias("els_trk_p4") != 0) {
		els_trk_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_p4"));
		els_trk_p4_branch->SetAddress(&els_trk_p4_);
	}
	if(els_trk_p4_branch == 0 ) {
	cout << "Branch els_trk_p4 does not exist." << endl;
	}
	els_vertex_p4_branch = 0;
	if (tree->GetAlias("els_vertex_p4") != 0) {
		els_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_vertex_p4"));
		els_vertex_p4_branch->SetAddress(&els_vertex_p4_);
	}
	if(els_vertex_p4_branch == 0 ) {
	cout << "Branch els_vertex_p4 does not exist." << endl;
	}
	genps_lepdaughter_p4_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_p4") != 0) {
		genps_lepdaughter_p4_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_p4"));
		genps_lepdaughter_p4_branch->SetAddress(&genps_lepdaughter_p4_);
	}
	if(genps_lepdaughter_p4_branch == 0 ) {
	cout << "Branch genps_lepdaughter_p4 does not exist." << endl;
	}
	genps_p4_branch = 0;
	if (tree->GetAlias("genps_p4") != 0) {
		genps_p4_branch = tree->GetBranch(tree->GetAlias("genps_p4"));
		genps_p4_branch->SetAddress(&genps_p4_);
	}
	if(genps_p4_branch == 0 ) {
	cout << "Branch genps_p4 does not exist." << endl;
	}
	genps_prod_vtx_branch = 0;
	if (tree->GetAlias("genps_prod_vtx") != 0) {
		genps_prod_vtx_branch = tree->GetBranch(tree->GetAlias("genps_prod_vtx"));
		genps_prod_vtx_branch->SetAddress(&genps_prod_vtx_);
	}
	if(genps_prod_vtx_branch == 0 ) {
	cout << "Branch genps_prod_vtx does not exist." << endl;
	}
	mus_gfit_outerPos_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_outerPos_p4") != 0) {
		mus_gfit_outerPos_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_outerPos_p4"));
		mus_gfit_outerPos_p4_branch->SetAddress(&mus_gfit_outerPos_p4_);
	}
	if(mus_gfit_outerPos_p4_branch == 0 ) {
	cout << "Branch mus_gfit_outerPos_p4 does not exist." << endl;
	}
	mus_p4_branch = 0;
	if (tree->GetAlias("mus_p4") != 0) {
		mus_p4_branch = tree->GetBranch(tree->GetAlias("mus_p4"));
		mus_p4_branch->SetAddress(&mus_p4_);
	}
	if(mus_p4_branch == 0 ) {
	cout << "Branch mus_p4 does not exist." << endl;
	}
	mus_trk_p4_branch = 0;
	if (tree->GetAlias("mus_trk_p4") != 0) {
		mus_trk_p4_branch = tree->GetBranch(tree->GetAlias("mus_trk_p4"));
		mus_trk_p4_branch->SetAddress(&mus_trk_p4_);
	}
	if(mus_trk_p4_branch == 0 ) {
	cout << "Branch mus_trk_p4 does not exist." << endl;
	}
	mus_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_vertex_p4") != 0) {
		mus_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_vertex_p4"));
		mus_vertex_p4_branch->SetAddress(&mus_vertex_p4_);
	}
	if(mus_vertex_p4_branch == 0 ) {
	cout << "Branch mus_vertex_p4 does not exist." << endl;
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
	evt_HLT_trigNames_branch = 0;
	if (tree->GetAlias("evt_HLT_trigNames") != 0) {
		evt_HLT_trigNames_branch = tree->GetBranch(tree->GetAlias("evt_HLT_trigNames"));
		evt_HLT_trigNames_branch->SetAddress(&evt_HLT_trigNames_);
	}
	if(evt_HLT_trigNames_branch == 0 ) {
	cout << "Branch evt_HLT_trigNames does not exist." << endl;
	}
	evt_L1_trigNames_branch = 0;
	if (tree->GetAlias("evt_L1_trigNames") != 0) {
		evt_L1_trigNames_branch = tree->GetBranch(tree->GetAlias("evt_L1_trigNames"));
		evt_L1_trigNames_branch->SetAddress(&evt_L1_trigNames_);
	}
	if(evt_L1_trigNames_branch == 0 ) {
	cout << "Branch evt_L1_trigNames does not exist." << endl;
	}
	genps_pthat_branch = 0;
	if (tree->GetAlias("genps_pthat") != 0) {
		genps_pthat_branch = tree->GetBranch(tree->GetAlias("genps_pthat"));
		genps_pthat_branch->SetAddress(&genps_pthat_);
	}
	if(genps_pthat_branch == 0 ) {
	cout << "Branch genps_pthat does not exist." << endl;
	}
	els_inner_positionx_branch = 0;
	if (tree->GetAlias("els_inner_positionx") != 0) {
		els_inner_positionx_branch = tree->GetBranch(tree->GetAlias("els_inner_positionx"));
		els_inner_positionx_branch->SetAddress(&els_inner_positionx_);
	}
	if(els_inner_positionx_branch == 0 ) {
	cout << "Branch els_inner_positionx does not exist." << endl;
	}
	els_inner_positiony_branch = 0;
	if (tree->GetAlias("els_inner_positiony") != 0) {
		els_inner_positiony_branch = tree->GetBranch(tree->GetAlias("els_inner_positiony"));
		els_inner_positiony_branch->SetAddress(&els_inner_positiony_);
	}
	if(els_inner_positiony_branch == 0 ) {
	cout << "Branch els_inner_positiony does not exist." << endl;
	}
	els_inner_positionz_branch = 0;
	if (tree->GetAlias("els_inner_positionz") != 0) {
		els_inner_positionz_branch = tree->GetBranch(tree->GetAlias("els_inner_positionz"));
		els_inner_positionz_branch->SetAddress(&els_inner_positionz_);
	}
	if(els_inner_positionz_branch == 0 ) {
	cout << "Branch els_inner_positionz does not exist." << endl;
	}
	trks_inner_positionx_branch = 0;
	if (tree->GetAlias("trks_inner_positionx") != 0) {
		trks_inner_positionx_branch = tree->GetBranch(tree->GetAlias("trks_inner_positionx"));
		trks_inner_positionx_branch->SetAddress(&trks_inner_positionx_);
	}
	if(trks_inner_positionx_branch == 0 ) {
	cout << "Branch trks_inner_positionx does not exist." << endl;
	}
	trks_inner_positiony_branch = 0;
	if (tree->GetAlias("trks_inner_positiony") != 0) {
		trks_inner_positiony_branch = tree->GetBranch(tree->GetAlias("trks_inner_positiony"));
		trks_inner_positiony_branch->SetAddress(&trks_inner_positiony_);
	}
	if(trks_inner_positiony_branch == 0 ) {
	cout << "Branch trks_inner_positiony does not exist." << endl;
	}
	trks_inner_positionz_branch = 0;
	if (tree->GetAlias("trks_inner_positionz") != 0) {
		trks_inner_positionz_branch = tree->GetBranch(tree->GetAlias("trks_inner_positionz"));
		trks_inner_positionz_branch->SetAddress(&trks_inner_positionz_);
	}
	if(trks_inner_positionz_branch == 0 ) {
	cout << "Branch trks_inner_positionz does not exist." << endl;
	}
	evt_bs_dxdz_branch = 0;
	if (tree->GetAlias("evt_bs_dxdz") != 0) {
		evt_bs_dxdz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdz"));
		evt_bs_dxdz_branch->SetAddress(&evt_bs_dxdz_);
	}
	if(evt_bs_dxdz_branch == 0 ) {
	cout << "Branch evt_bs_dxdz does not exist." << endl;
	}
	evt_bs_dxdzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dxdzErr") != 0) {
		evt_bs_dxdzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdzErr"));
		evt_bs_dxdzErr_branch->SetAddress(&evt_bs_dxdzErr_);
	}
	if(evt_bs_dxdzErr_branch == 0 ) {
	cout << "Branch evt_bs_dxdzErr does not exist." << endl;
	}
	evt_bs_dydz_branch = 0;
	if (tree->GetAlias("evt_bs_dydz") != 0) {
		evt_bs_dydz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydz"));
		evt_bs_dydz_branch->SetAddress(&evt_bs_dydz_);
	}
	if(evt_bs_dydz_branch == 0 ) {
	cout << "Branch evt_bs_dydz does not exist." << endl;
	}
	evt_bs_dydzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dydzErr") != 0) {
		evt_bs_dydzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydzErr"));
		evt_bs_dydzErr_branch->SetAddress(&evt_bs_dydzErr_);
	}
	if(evt_bs_dydzErr_branch == 0 ) {
	cout << "Branch evt_bs_dydzErr does not exist." << endl;
	}
	evt_bs_sigmaZ_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZ") != 0) {
		evt_bs_sigmaZ_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZ"));
		evt_bs_sigmaZ_branch->SetAddress(&evt_bs_sigmaZ_);
	}
	if(evt_bs_sigmaZ_branch == 0 ) {
	cout << "Branch evt_bs_sigmaZ does not exist." << endl;
	}
	evt_bs_sigmaZErr_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZErr") != 0) {
		evt_bs_sigmaZErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZErr"));
		evt_bs_sigmaZErr_branch->SetAddress(&evt_bs_sigmaZErr_);
	}
	if(evt_bs_sigmaZErr_branch == 0 ) {
	cout << "Branch evt_bs_sigmaZErr does not exist." << endl;
	}
	evt_bs_widthx_branch = 0;
	if (tree->GetAlias("evt_bs_widthx") != 0) {
		evt_bs_widthx_branch = tree->GetBranch(tree->GetAlias("evt_bs_widthx"));
		evt_bs_widthx_branch->SetAddress(&evt_bs_widthx_);
	}
	if(evt_bs_widthx_branch == 0 ) {
	cout << "Branch evt_bs_widthx does not exist." << endl;
	}
	evt_bs_widthxErr_branch = 0;
	if (tree->GetAlias("evt_bs_widthxErr") != 0) {
		evt_bs_widthxErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_widthxErr"));
		evt_bs_widthxErr_branch->SetAddress(&evt_bs_widthxErr_);
	}
	if(evt_bs_widthxErr_branch == 0 ) {
	cout << "Branch evt_bs_widthxErr does not exist." << endl;
	}
	evt_bs_widthy_branch = 0;
	if (tree->GetAlias("evt_bs_widthy") != 0) {
		evt_bs_widthy_branch = tree->GetBranch(tree->GetAlias("evt_bs_widthy"));
		evt_bs_widthy_branch->SetAddress(&evt_bs_widthy_);
	}
	if(evt_bs_widthy_branch == 0 ) {
	cout << "Branch evt_bs_widthy does not exist." << endl;
	}
	evt_bs_widthyErr_branch = 0;
	if (tree->GetAlias("evt_bs_widthyErr") != 0) {
		evt_bs_widthyErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_widthyErr"));
		evt_bs_widthyErr_branch->SetAddress(&evt_bs_widthyErr_);
	}
	if(evt_bs_widthyErr_branch == 0 ) {
	cout << "Branch evt_bs_widthyErr does not exist." << endl;
	}
	evt_bs_xErr_branch = 0;
	if (tree->GetAlias("evt_bs_xErr") != 0) {
		evt_bs_xErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_xErr"));
		evt_bs_xErr_branch->SetAddress(&evt_bs_xErr_);
	}
	if(evt_bs_xErr_branch == 0 ) {
	cout << "Branch evt_bs_xErr does not exist." << endl;
	}
	evt_bs_yErr_branch = 0;
	if (tree->GetAlias("evt_bs_yErr") != 0) {
		evt_bs_yErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_yErr"));
		evt_bs_yErr_branch->SetAddress(&evt_bs_yErr_);
	}
	if(evt_bs_yErr_branch == 0 ) {
	cout << "Branch evt_bs_yErr does not exist." << endl;
	}
	evt_bs_zErr_branch = 0;
	if (tree->GetAlias("evt_bs_zErr") != 0) {
		evt_bs_zErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_zErr"));
		evt_bs_zErr_branch->SetAddress(&evt_bs_zErr_);
	}
	if(evt_bs_zErr_branch == 0 ) {
	cout << "Branch evt_bs_zErr does not exist." << endl;
	}
	evt_bField_branch = 0;
	if (tree->GetAlias("evt_bField") != 0) {
		evt_bField_branch = tree->GetBranch(tree->GetAlias("evt_bField"));
		evt_bField_branch->SetAddress(&evt_bField_);
	}
	if(evt_bField_branch == 0 ) {
	cout << "Branch evt_bField does not exist." << endl;
	}
	evt_kfactor_branch = 0;
	if (tree->GetAlias("evt_kfactor") != 0) {
		evt_kfactor_branch = tree->GetBranch(tree->GetAlias("evt_kfactor"));
		evt_kfactor_branch->SetAddress(&evt_kfactor_);
	}
	if(evt_kfactor_branch == 0 ) {
	cout << "Branch evt_kfactor does not exist." << endl;
	}
	evt_weight_branch = 0;
	if (tree->GetAlias("evt_weight") != 0) {
		evt_weight_branch = tree->GetBranch(tree->GetAlias("evt_weight"));
		evt_weight_branch->SetAddress(&evt_weight_);
	}
	if(evt_weight_branch == 0 ) {
	cout << "Branch evt_weight does not exist." << endl;
	}
	evt_xsec_excl_branch = 0;
	if (tree->GetAlias("evt_xsec_excl") != 0) {
		evt_xsec_excl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_excl"));
		evt_xsec_excl_branch->SetAddress(&evt_xsec_excl_);
	}
	if(evt_xsec_excl_branch == 0 ) {
	cout << "Branch evt_xsec_excl does not exist." << endl;
	}
	evt_xsec_incl_branch = 0;
	if (tree->GetAlias("evt_xsec_incl") != 0) {
		evt_xsec_incl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_incl"));
		evt_xsec_incl_branch->SetAddress(&evt_xsec_incl_);
	}
	if(evt_xsec_incl_branch == 0 ) {
	cout << "Branch evt_xsec_incl does not exist." << endl;
	}
	gen_met_branch = 0;
	if (tree->GetAlias("gen_met") != 0) {
		gen_met_branch = tree->GetBranch(tree->GetAlias("gen_met"));
		gen_met_branch->SetAddress(&gen_met_);
	}
	if(gen_met_branch == 0 ) {
	cout << "Branch gen_met does not exist." << endl;
	}
	gen_metPhi_branch = 0;
	if (tree->GetAlias("gen_metPhi") != 0) {
		gen_metPhi_branch = tree->GetBranch(tree->GetAlias("gen_metPhi"));
		gen_metPhi_branch->SetAddress(&gen_metPhi_);
	}
	if(gen_metPhi_branch == 0 ) {
	cout << "Branch gen_metPhi does not exist." << endl;
	}
	els_mc3dr_branch = 0;
	if (tree->GetAlias("els_mc3dr") != 0) {
		els_mc3dr_branch = tree->GetBranch(tree->GetAlias("els_mc3dr"));
		els_mc3dr_branch->SetAddress(&els_mc3dr_);
	}
	if(els_mc3dr_branch == 0 ) {
	cout << "Branch els_mc3dr does not exist." << endl;
	}
	els_mcdr_branch = 0;
	if (tree->GetAlias("els_mcdr") != 0) {
		els_mcdr_branch = tree->GetBranch(tree->GetAlias("els_mcdr"));
		els_mcdr_branch->SetAddress(&els_mcdr_);
	}
	if(els_mcdr_branch == 0 ) {
	cout << "Branch els_mcdr does not exist." << endl;
	}
	jets_mcdr_branch = 0;
	if (tree->GetAlias("jets_mcdr") != 0) {
		jets_mcdr_branch = tree->GetBranch(tree->GetAlias("jets_mcdr"));
		jets_mcdr_branch->SetAddress(&jets_mcdr_);
	}
	if(jets_mcdr_branch == 0 ) {
	cout << "Branch jets_mcdr does not exist." << endl;
	}
	jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_emEnergy") != 0) {
		jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_emEnergy"));
		jets_mc_emEnergy_branch->SetAddress(&jets_mc_emEnergy_);
	}
	if(jets_mc_emEnergy_branch == 0 ) {
	cout << "Branch jets_mc_emEnergy does not exist." << endl;
	}
	jets_mc_gpdr_branch = 0;
	if (tree->GetAlias("jets_mc_gpdr") != 0) {
		jets_mc_gpdr_branch = tree->GetBranch(tree->GetAlias("jets_mc_gpdr"));
		jets_mc_gpdr_branch->SetAddress(&jets_mc_gpdr_);
	}
	if(jets_mc_gpdr_branch == 0 ) {
	cout << "Branch jets_mc_gpdr does not exist." << endl;
	}
	jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_hadEnergy") != 0) {
		jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_hadEnergy"));
		jets_mc_hadEnergy_branch->SetAddress(&jets_mc_hadEnergy_);
	}
	if(jets_mc_hadEnergy_branch == 0 ) {
	cout << "Branch jets_mc_hadEnergy does not exist." << endl;
	}
	jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_invEnergy") != 0) {
		jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_invEnergy"));
		jets_mc_invEnergy_branch->SetAddress(&jets_mc_invEnergy_);
	}
	if(jets_mc_invEnergy_branch == 0 ) {
	cout << "Branch jets_mc_invEnergy does not exist." << endl;
	}
	jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_otherEnergy") != 0) {
		jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_otherEnergy"));
		jets_mc_otherEnergy_branch->SetAddress(&jets_mc_otherEnergy_);
	}
	if(jets_mc_otherEnergy_branch == 0 ) {
	cout << "Branch jets_mc_otherEnergy does not exist." << endl;
	}
	mus_mc3dr_branch = 0;
	if (tree->GetAlias("mus_mc3dr") != 0) {
		mus_mc3dr_branch = tree->GetBranch(tree->GetAlias("mus_mc3dr"));
		mus_mc3dr_branch->SetAddress(&mus_mc3dr_);
	}
	if(mus_mc3dr_branch == 0 ) {
	cout << "Branch mus_mc3dr does not exist." << endl;
	}
	mus_mcdr_branch = 0;
	if (tree->GetAlias("mus_mcdr") != 0) {
		mus_mcdr_branch = tree->GetBranch(tree->GetAlias("mus_mcdr"));
		mus_mcdr_branch->SetAddress(&mus_mcdr_);
	}
	if(mus_mcdr_branch == 0 ) {
	cout << "Branch mus_mcdr does not exist." << endl;
	}
	trk_mc3dr_branch = 0;
	if (tree->GetAlias("trk_mc3dr") != 0) {
		trk_mc3dr_branch = tree->GetBranch(tree->GetAlias("trk_mc3dr"));
		trk_mc3dr_branch->SetAddress(&trk_mc3dr_);
	}
	if(trk_mc3dr_branch == 0 ) {
	cout << "Branch trk_mc3dr does not exist." << endl;
	}
	trk_mcdr_branch = 0;
	if (tree->GetAlias("trk_mcdr") != 0) {
		trk_mcdr_branch = tree->GetBranch(tree->GetAlias("trk_mcdr"));
		trk_mcdr_branch->SetAddress(&trk_mcdr_);
	}
	if(trk_mcdr_branch == 0 ) {
	cout << "Branch trk_mcdr does not exist." << endl;
	}
	els_conv_dcot_branch = 0;
	if (tree->GetAlias("els_conv_dcot") != 0) {
		els_conv_dcot_branch = tree->GetBranch(tree->GetAlias("els_conv_dcot"));
		els_conv_dcot_branch->SetAddress(&els_conv_dcot_);
	}
	if(els_conv_dcot_branch == 0 ) {
	cout << "Branch els_conv_dcot does not exist." << endl;
	}
	els_conv_dist_branch = 0;
	if (tree->GetAlias("els_conv_dist") != 0) {
		els_conv_dist_branch = tree->GetBranch(tree->GetAlias("els_conv_dist"));
		els_conv_dist_branch->SetAddress(&els_conv_dist_);
	}
	if(els_conv_dist_branch == 0 ) {
	cout << "Branch els_conv_dist does not exist." << endl;
	}
	trks_conv_dcot_branch = 0;
	if (tree->GetAlias("trks_conv_dcot") != 0) {
		trks_conv_dcot_branch = tree->GetBranch(tree->GetAlias("trks_conv_dcot"));
		trks_conv_dcot_branch->SetAddress(&trks_conv_dcot_);
	}
	if(trks_conv_dcot_branch == 0 ) {
	cout << "Branch trks_conv_dcot does not exist." << endl;
	}
	trks_conv_dist_branch = 0;
	if (tree->GetAlias("trks_conv_dist") != 0) {
		trks_conv_dist_branch = tree->GetBranch(tree->GetAlias("trks_conv_dist"));
		trks_conv_dist_branch->SetAddress(&trks_conv_dist_);
	}
	if(trks_conv_dist_branch == 0 ) {
	cout << "Branch trks_conv_dist does not exist." << endl;
	}
	els_musdr_branch = 0;
	if (tree->GetAlias("els_musdr") != 0) {
		els_musdr_branch = tree->GetBranch(tree->GetAlias("els_musdr"));
		els_musdr_branch->SetAddress(&els_musdr_);
	}
	if(els_musdr_branch == 0 ) {
	cout << "Branch els_musdr does not exist." << endl;
	}
	els_trkdr_branch = 0;
	if (tree->GetAlias("els_trkdr") != 0) {
		els_trkdr_branch = tree->GetBranch(tree->GetAlias("els_trkdr"));
		els_trkdr_branch->SetAddress(&els_trkdr_);
	}
	if(els_trkdr_branch == 0 ) {
	cout << "Branch els_trkdr does not exist." << endl;
	}
	els_trkshFrac_branch = 0;
	if (tree->GetAlias("els_trkshFrac") != 0) {
		els_trkshFrac_branch = tree->GetBranch(tree->GetAlias("els_trkshFrac"));
		els_trkshFrac_branch->SetAddress(&els_trkshFrac_);
	}
	if(els_trkshFrac_branch == 0 ) {
	cout << "Branch els_trkshFrac does not exist." << endl;
	}
	els_chi2_branch = 0;
	if (tree->GetAlias("els_chi2") != 0) {
		els_chi2_branch = tree->GetBranch(tree->GetAlias("els_chi2"));
		els_chi2_branch->SetAddress(&els_chi2_);
	}
	if(els_chi2_branch == 0 ) {
	cout << "Branch els_chi2 does not exist." << endl;
	}
	els_d0_branch = 0;
	if (tree->GetAlias("els_d0") != 0) {
		els_d0_branch = tree->GetBranch(tree->GetAlias("els_d0"));
		els_d0_branch->SetAddress(&els_d0_);
	}
	if(els_d0_branch == 0 ) {
	cout << "Branch els_d0 does not exist." << endl;
	}
	els_d0Err_branch = 0;
	if (tree->GetAlias("els_d0Err") != 0) {
		els_d0Err_branch = tree->GetBranch(tree->GetAlias("els_d0Err"));
		els_d0Err_branch->SetAddress(&els_d0Err_);
	}
	if(els_d0Err_branch == 0 ) {
	cout << "Branch els_d0Err does not exist." << endl;
	}
	els_d0corr_branch = 0;
	if (tree->GetAlias("els_d0corr") != 0) {
		els_d0corr_branch = tree->GetBranch(tree->GetAlias("els_d0corr"));
		els_d0corr_branch->SetAddress(&els_d0corr_);
	}
	if(els_d0corr_branch == 0 ) {
	cout << "Branch els_d0corr does not exist." << endl;
	}
	els_dEtaIn_branch = 0;
	if (tree->GetAlias("els_dEtaIn") != 0) {
		els_dEtaIn_branch = tree->GetBranch(tree->GetAlias("els_dEtaIn"));
		els_dEtaIn_branch->SetAddress(&els_dEtaIn_);
	}
	if(els_dEtaIn_branch == 0 ) {
	cout << "Branch els_dEtaIn does not exist." << endl;
	}
	els_dEtaOut_branch = 0;
	if (tree->GetAlias("els_dEtaOut") != 0) {
		els_dEtaOut_branch = tree->GetBranch(tree->GetAlias("els_dEtaOut"));
		els_dEtaOut_branch->SetAddress(&els_dEtaOut_);
	}
	if(els_dEtaOut_branch == 0 ) {
	cout << "Branch els_dEtaOut does not exist." << endl;
	}
	els_dPhiIn_branch = 0;
	if (tree->GetAlias("els_dPhiIn") != 0) {
		els_dPhiIn_branch = tree->GetBranch(tree->GetAlias("els_dPhiIn"));
		els_dPhiIn_branch->SetAddress(&els_dPhiIn_);
	}
	if(els_dPhiIn_branch == 0 ) {
	cout << "Branch els_dPhiIn does not exist." << endl;
	}
	els_dPhiInPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiInPhiOut") != 0) {
		els_dPhiInPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiInPhiOut"));
		els_dPhiInPhiOut_branch->SetAddress(&els_dPhiInPhiOut_);
	}
	if(els_dPhiInPhiOut_branch == 0 ) {
	cout << "Branch els_dPhiInPhiOut does not exist." << endl;
	}
	els_dPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiOut") != 0) {
		els_dPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiOut"));
		els_dPhiOut_branch->SetAddress(&els_dPhiOut_);
	}
	if(els_dPhiOut_branch == 0 ) {
	cout << "Branch els_dPhiOut does not exist." << endl;
	}
	els_e1x5_branch = 0;
	if (tree->GetAlias("els_e1x5") != 0) {
		els_e1x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5"));
		els_e1x5_branch->SetAddress(&els_e1x5_);
	}
	if(els_e1x5_branch == 0 ) {
	cout << "Branch els_e1x5 does not exist." << endl;
	}
	els_e2x5Max_branch = 0;
	if (tree->GetAlias("els_e2x5Max") != 0) {
		els_e2x5Max_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max"));
		els_e2x5Max_branch->SetAddress(&els_e2x5Max_);
	}
	if(els_e2x5Max_branch == 0 ) {
	cout << "Branch els_e2x5Max does not exist." << endl;
	}
	els_e3x3_branch = 0;
	if (tree->GetAlias("els_e3x3") != 0) {
		els_e3x3_branch = tree->GetBranch(tree->GetAlias("els_e3x3"));
		els_e3x3_branch->SetAddress(&els_e3x3_);
	}
	if(els_e3x3_branch == 0 ) {
	cout << "Branch els_e3x3 does not exist." << endl;
	}
	els_e5x5_branch = 0;
	if (tree->GetAlias("els_e5x5") != 0) {
		els_e5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5"));
		els_e5x5_branch->SetAddress(&els_e5x5_);
	}
	if(els_e5x5_branch == 0 ) {
	cout << "Branch els_e5x5 does not exist." << endl;
	}
	els_eMax_branch = 0;
	if (tree->GetAlias("els_eMax") != 0) {
		els_eMax_branch = tree->GetBranch(tree->GetAlias("els_eMax"));
		els_eMax_branch->SetAddress(&els_eMax_);
	}
	if(els_eMax_branch == 0 ) {
	cout << "Branch els_eMax does not exist." << endl;
	}
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		els_eOverPIn_branch->SetAddress(&els_eOverPIn_);
	}
	if(els_eOverPIn_branch == 0 ) {
	cout << "Branch els_eOverPIn does not exist." << endl;
	}
	els_eSC_branch = 0;
	if (tree->GetAlias("els_eSC") != 0) {
		els_eSC_branch = tree->GetBranch(tree->GetAlias("els_eSC"));
		els_eSC_branch->SetAddress(&els_eSC_);
	}
	if(els_eSC_branch == 0 ) {
	cout << "Branch els_eSC does not exist." << endl;
	}
	els_eSCPresh_branch = 0;
	if (tree->GetAlias("els_eSCPresh") != 0) {
		els_eSCPresh_branch = tree->GetBranch(tree->GetAlias("els_eSCPresh"));
		els_eSCPresh_branch->SetAddress(&els_eSCPresh_);
	}
	if(els_eSCPresh_branch == 0 ) {
	cout << "Branch els_eSCPresh does not exist." << endl;
	}
	els_eSCRaw_branch = 0;
	if (tree->GetAlias("els_eSCRaw") != 0) {
		els_eSCRaw_branch = tree->GetBranch(tree->GetAlias("els_eSCRaw"));
		els_eSCRaw_branch->SetAddress(&els_eSCRaw_);
	}
	if(els_eSCRaw_branch == 0 ) {
	cout << "Branch els_eSCRaw does not exist." << endl;
	}
	els_eSeed_branch = 0;
	if (tree->GetAlias("els_eSeed") != 0) {
		els_eSeed_branch = tree->GetBranch(tree->GetAlias("els_eSeed"));
		els_eSeed_branch->SetAddress(&els_eSeed_);
	}
	if(els_eSeed_branch == 0 ) {
	cout << "Branch els_eSeed does not exist." << endl;
	}
	els_eSeedOverPIn_branch = 0;
	if (tree->GetAlias("els_eSeedOverPIn") != 0) {
		els_eSeedOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPIn"));
		els_eSeedOverPIn_branch->SetAddress(&els_eSeedOverPIn_);
	}
	if(els_eSeedOverPIn_branch == 0 ) {
	cout << "Branch els_eSeedOverPIn does not exist." << endl;
	}
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut_);
	}
	if(els_eSeedOverPOut_branch == 0 ) {
	cout << "Branch els_eSeedOverPOut does not exist." << endl;
	}
	els_ecalIso03_branch = 0;
	if (tree->GetAlias("els_ecalIso03") != 0) {
		els_ecalIso03_branch = tree->GetBranch(tree->GetAlias("els_ecalIso03"));
		els_ecalIso03_branch->SetAddress(&els_ecalIso03_);
	}
	if(els_ecalIso03_branch == 0 ) {
	cout << "Branch els_ecalIso03 does not exist." << endl;
	}
	els_ecalIso04_branch = 0;
	if (tree->GetAlias("els_ecalIso04") != 0) {
		els_ecalIso04_branch = tree->GetBranch(tree->GetAlias("els_ecalIso04"));
		els_ecalIso04_branch->SetAddress(&els_ecalIso04_);
	}
	if(els_ecalIso04_branch == 0 ) {
	cout << "Branch els_ecalIso04 does not exist." << endl;
	}
	els_egamma_looseId_branch = 0;
	if (tree->GetAlias("els_egamma_looseId") != 0) {
		els_egamma_looseId_branch = tree->GetBranch(tree->GetAlias("els_egamma_looseId"));
		els_egamma_looseId_branch->SetAddress(&els_egamma_looseId_);
	}
	if(els_egamma_looseId_branch == 0 ) {
	cout << "Branch els_egamma_looseId does not exist." << endl;
	}
	els_egamma_robustHighEnergy_branch = 0;
	if (tree->GetAlias("els_egamma_robustHighEnergy") != 0) {
		els_egamma_robustHighEnergy_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustHighEnergy"));
		els_egamma_robustHighEnergy_branch->SetAddress(&els_egamma_robustHighEnergy_);
	}
	if(els_egamma_robustHighEnergy_branch == 0 ) {
	cout << "Branch els_egamma_robustHighEnergy does not exist." << endl;
	}
	els_egamma_robustLooseId_branch = 0;
	if (tree->GetAlias("els_egamma_robustLooseId") != 0) {
		els_egamma_robustLooseId_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustLooseId"));
		els_egamma_robustLooseId_branch->SetAddress(&els_egamma_robustLooseId_);
	}
	if(els_egamma_robustLooseId_branch == 0 ) {
	cout << "Branch els_egamma_robustLooseId does not exist." << endl;
	}
	els_egamma_robustTightId_branch = 0;
	if (tree->GetAlias("els_egamma_robustTightId") != 0) {
		els_egamma_robustTightId_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustTightId"));
		els_egamma_robustTightId_branch->SetAddress(&els_egamma_robustTightId_);
	}
	if(els_egamma_robustTightId_branch == 0 ) {
	cout << "Branch els_egamma_robustTightId does not exist." << endl;
	}
	els_egamma_tightId_branch = 0;
	if (tree->GetAlias("els_egamma_tightId") != 0) {
		els_egamma_tightId_branch = tree->GetBranch(tree->GetAlias("els_egamma_tightId"));
		els_egamma_tightId_branch->SetAddress(&els_egamma_tightId_);
	}
	if(els_egamma_tightId_branch == 0 ) {
	cout << "Branch els_egamma_tightId does not exist." << endl;
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		els_etaErr_branch->SetAddress(&els_etaErr_);
	}
	if(els_etaErr_branch == 0 ) {
	cout << "Branch els_etaErr does not exist." << endl;
	}
	els_etaSC_branch = 0;
	if (tree->GetAlias("els_etaSC") != 0) {
		els_etaSC_branch = tree->GetBranch(tree->GetAlias("els_etaSC"));
		els_etaSC_branch->SetAddress(&els_etaSC_);
	}
	if(els_etaSC_branch == 0 ) {
	cout << "Branch els_etaSC does not exist." << endl;
	}
	els_fBrem_branch = 0;
	if (tree->GetAlias("els_fBrem") != 0) {
		els_fBrem_branch = tree->GetBranch(tree->GetAlias("els_fBrem"));
		els_fBrem_branch->SetAddress(&els_fBrem_);
	}
	if(els_fBrem_branch == 0 ) {
	cout << "Branch els_fBrem does not exist." << endl;
	}
	els_hOverE_branch = 0;
	if (tree->GetAlias("els_hOverE") != 0) {
		els_hOverE_branch = tree->GetBranch(tree->GetAlias("els_hOverE"));
		els_hOverE_branch->SetAddress(&els_hOverE_);
	}
	if(els_hOverE_branch == 0 ) {
	cout << "Branch els_hOverE does not exist." << endl;
	}
	els_hcalIso03_branch = 0;
	if (tree->GetAlias("els_hcalIso03") != 0) {
		els_hcalIso03_branch = tree->GetBranch(tree->GetAlias("els_hcalIso03"));
		els_hcalIso03_branch->SetAddress(&els_hcalIso03_);
	}
	if(els_hcalIso03_branch == 0 ) {
	cout << "Branch els_hcalIso03 does not exist." << endl;
	}
	els_hcalIso04_branch = 0;
	if (tree->GetAlias("els_hcalIso04") != 0) {
		els_hcalIso04_branch = tree->GetBranch(tree->GetAlias("els_hcalIso04"));
		els_hcalIso04_branch->SetAddress(&els_hcalIso04_);
	}
	if(els_hcalIso04_branch == 0 ) {
	cout << "Branch els_hcalIso04 does not exist." << endl;
	}
	els_layer1_charge_branch = 0;
	if (tree->GetAlias("els_layer1_charge") != 0) {
		els_layer1_charge_branch = tree->GetBranch(tree->GetAlias("els_layer1_charge"));
		els_layer1_charge_branch->SetAddress(&els_layer1_charge_);
	}
	if(els_layer1_charge_branch == 0 ) {
	cout << "Branch els_layer1_charge does not exist." << endl;
	}
	els_ndof_branch = 0;
	if (tree->GetAlias("els_ndof") != 0) {
		els_ndof_branch = tree->GetBranch(tree->GetAlias("els_ndof"));
		els_ndof_branch->SetAddress(&els_ndof_);
	}
	if(els_ndof_branch == 0 ) {
	cout << "Branch els_ndof does not exist." << endl;
	}
	els_outerEta_branch = 0;
	if (tree->GetAlias("els_outerEta") != 0) {
		els_outerEta_branch = tree->GetBranch(tree->GetAlias("els_outerEta"));
		els_outerEta_branch->SetAddress(&els_outerEta_);
	}
	if(els_outerEta_branch == 0 ) {
	cout << "Branch els_outerEta does not exist." << endl;
	}
	els_outerPhi_branch = 0;
	if (tree->GetAlias("els_outerPhi") != 0) {
		els_outerPhi_branch = tree->GetBranch(tree->GetAlias("els_outerPhi"));
		els_outerPhi_branch->SetAddress(&els_outerPhi_);
	}
	if(els_outerPhi_branch == 0 ) {
	cout << "Branch els_outerPhi does not exist." << endl;
	}
	els_phiErr_branch = 0;
	if (tree->GetAlias("els_phiErr") != 0) {
		els_phiErr_branch = tree->GetBranch(tree->GetAlias("els_phiErr"));
		els_phiErr_branch->SetAddress(&els_phiErr_);
	}
	if(els_phiErr_branch == 0 ) {
	cout << "Branch els_phiErr does not exist." << endl;
	}
	els_phiSC_branch = 0;
	if (tree->GetAlias("els_phiSC") != 0) {
		els_phiSC_branch = tree->GetBranch(tree->GetAlias("els_phiSC"));
		els_phiSC_branch->SetAddress(&els_phiSC_);
	}
	if(els_phiSC_branch == 0 ) {
	cout << "Branch els_phiSC does not exist." << endl;
	}
	els_pin_branch = 0;
	if (tree->GetAlias("els_pin") != 0) {
		els_pin_branch = tree->GetBranch(tree->GetAlias("els_pin"));
		els_pin_branch->SetAddress(&els_pin_);
	}
	if(els_pin_branch == 0 ) {
	cout << "Branch els_pin does not exist." << endl;
	}
	els_pout_branch = 0;
	if (tree->GetAlias("els_pout") != 0) {
		els_pout_branch = tree->GetBranch(tree->GetAlias("els_pout"));
		els_pout_branch->SetAddress(&els_pout_);
	}
	if(els_pout_branch == 0 ) {
	cout << "Branch els_pout does not exist." << endl;
	}
	els_ptErr_branch = 0;
	if (tree->GetAlias("els_ptErr") != 0) {
		els_ptErr_branch = tree->GetBranch(tree->GetAlias("els_ptErr"));
		els_ptErr_branch->SetAddress(&els_ptErr_);
	}
	if(els_ptErr_branch == 0 ) {
	cout << "Branch els_ptErr does not exist." << endl;
	}
	els_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta") != 0) {
		els_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta"));
		els_sigmaEtaEta_branch->SetAddress(&els_sigmaEtaEta_);
	}
	if(els_sigmaEtaEta_branch == 0 ) {
	cout << "Branch els_sigmaEtaEta does not exist." << endl;
	}
	els_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta") != 0) {
		els_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta"));
		els_sigmaIEtaIEta_branch->SetAddress(&els_sigmaIEtaIEta_);
	}
	if(els_sigmaIEtaIEta_branch == 0 ) {
	cout << "Branch els_sigmaIEtaIEta does not exist." << endl;
	}
	els_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi") != 0) {
		els_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi"));
		els_sigmaIPhiIPhi_branch->SetAddress(&els_sigmaIPhiIPhi_);
	}
	if(els_sigmaIPhiIPhi_branch == 0 ) {
	cout << "Branch els_sigmaIPhiIPhi does not exist." << endl;
	}
	els_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("els_sigmaPhiPhi") != 0) {
		els_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaPhiPhi"));
		els_sigmaPhiPhi_branch->SetAddress(&els_sigmaPhiPhi_);
	}
	if(els_sigmaPhiPhi_branch == 0 ) {
	cout << "Branch els_sigmaPhiPhi does not exist." << endl;
	}
	els_tkIso03_branch = 0;
	if (tree->GetAlias("els_tkIso03") != 0) {
		els_tkIso03_branch = tree->GetBranch(tree->GetAlias("els_tkIso03"));
		els_tkIso03_branch->SetAddress(&els_tkIso03_);
	}
	if(els_tkIso03_branch == 0 ) {
	cout << "Branch els_tkIso03 does not exist." << endl;
	}
	els_tkIso04_branch = 0;
	if (tree->GetAlias("els_tkIso04") != 0) {
		els_tkIso04_branch = tree->GetBranch(tree->GetAlias("els_tkIso04"));
		els_tkIso04_branch->SetAddress(&els_tkIso04_);
	}
	if(els_tkIso04_branch == 0 ) {
	cout << "Branch els_tkIso04 does not exist." << endl;
	}
	els_vertexphi_branch = 0;
	if (tree->GetAlias("els_vertexphi") != 0) {
		els_vertexphi_branch = tree->GetBranch(tree->GetAlias("els_vertexphi"));
		els_vertexphi_branch->SetAddress(&els_vertexphi_);
	}
	if(els_vertexphi_branch == 0 ) {
	cout << "Branch els_vertexphi does not exist." << endl;
	}
	els_z0_branch = 0;
	if (tree->GetAlias("els_z0") != 0) {
		els_z0_branch = tree->GetBranch(tree->GetAlias("els_z0"));
		els_z0_branch->SetAddress(&els_z0_);
	}
	if(els_z0_branch == 0 ) {
	cout << "Branch els_z0 does not exist." << endl;
	}
	els_z0Err_branch = 0;
	if (tree->GetAlias("els_z0Err") != 0) {
		els_z0Err_branch = tree->GetBranch(tree->GetAlias("els_z0Err"));
		els_z0Err_branch->SetAddress(&els_z0Err_);
	}
	if(els_z0Err_branch == 0 ) {
	cout << "Branch els_z0Err does not exist." << endl;
	}
	els_z0corr_branch = 0;
	if (tree->GetAlias("els_z0corr") != 0) {
		els_z0corr_branch = tree->GetBranch(tree->GetAlias("els_z0corr"));
		els_z0corr_branch->SetAddress(&els_z0corr_);
	}
	if(els_z0corr_branch == 0 ) {
	cout << "Branch els_z0corr does not exist." << endl;
	}
	mus_caloCompatibility_branch = 0;
	if (tree->GetAlias("mus_caloCompatibility") != 0) {
		mus_caloCompatibility_branch = tree->GetBranch(tree->GetAlias("mus_caloCompatibility"));
		mus_caloCompatibility_branch->SetAddress(&mus_caloCompatibility_);
	}
	if(mus_caloCompatibility_branch == 0 ) {
	cout << "Branch mus_caloCompatibility does not exist." << endl;
	}
	mus_chi2_branch = 0;
	if (tree->GetAlias("mus_chi2") != 0) {
		mus_chi2_branch = tree->GetBranch(tree->GetAlias("mus_chi2"));
		mus_chi2_branch->SetAddress(&mus_chi2_);
	}
	if(mus_chi2_branch == 0 ) {
	cout << "Branch mus_chi2 does not exist." << endl;
	}
	mus_d0_branch = 0;
	if (tree->GetAlias("mus_d0") != 0) {
		mus_d0_branch = tree->GetBranch(tree->GetAlias("mus_d0"));
		mus_d0_branch->SetAddress(&mus_d0_);
	}
	if(mus_d0_branch == 0 ) {
	cout << "Branch mus_d0 does not exist." << endl;
	}
	mus_d0Err_branch = 0;
	if (tree->GetAlias("mus_d0Err") != 0) {
		mus_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_d0Err"));
		mus_d0Err_branch->SetAddress(&mus_d0Err_);
	}
	if(mus_d0Err_branch == 0 ) {
	cout << "Branch mus_d0Err does not exist." << endl;
	}
	mus_d0corr_branch = 0;
	if (tree->GetAlias("mus_d0corr") != 0) {
		mus_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_d0corr"));
		mus_d0corr_branch->SetAddress(&mus_d0corr_);
	}
	if(mus_d0corr_branch == 0 ) {
	cout << "Branch mus_d0corr does not exist." << endl;
	}
	mus_e_em_branch = 0;
	if (tree->GetAlias("mus_e_em") != 0) {
		mus_e_em_branch = tree->GetBranch(tree->GetAlias("mus_e_em"));
		mus_e_em_branch->SetAddress(&mus_e_em_);
	}
	if(mus_e_em_branch == 0 ) {
	cout << "Branch mus_e_em does not exist." << endl;
	}
	mus_e_emS9_branch = 0;
	if (tree->GetAlias("mus_e_emS9") != 0) {
		mus_e_emS9_branch = tree->GetBranch(tree->GetAlias("mus_e_emS9"));
		mus_e_emS9_branch->SetAddress(&mus_e_emS9_);
	}
	if(mus_e_emS9_branch == 0 ) {
	cout << "Branch mus_e_emS9 does not exist." << endl;
	}
	mus_e_had_branch = 0;
	if (tree->GetAlias("mus_e_had") != 0) {
		mus_e_had_branch = tree->GetBranch(tree->GetAlias("mus_e_had"));
		mus_e_had_branch->SetAddress(&mus_e_had_);
	}
	if(mus_e_had_branch == 0 ) {
	cout << "Branch mus_e_had does not exist." << endl;
	}
	mus_e_hadS9_branch = 0;
	if (tree->GetAlias("mus_e_hadS9") != 0) {
		mus_e_hadS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hadS9"));
		mus_e_hadS9_branch->SetAddress(&mus_e_hadS9_);
	}
	if(mus_e_hadS9_branch == 0 ) {
	cout << "Branch mus_e_hadS9 does not exist." << endl;
	}
	mus_e_ho_branch = 0;
	if (tree->GetAlias("mus_e_ho") != 0) {
		mus_e_ho_branch = tree->GetBranch(tree->GetAlias("mus_e_ho"));
		mus_e_ho_branch->SetAddress(&mus_e_ho_);
	}
	if(mus_e_ho_branch == 0 ) {
	cout << "Branch mus_e_ho does not exist." << endl;
	}
	mus_e_hoS9_branch = 0;
	if (tree->GetAlias("mus_e_hoS9") != 0) {
		mus_e_hoS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hoS9"));
		mus_e_hoS9_branch->SetAddress(&mus_e_hoS9_);
	}
	if(mus_e_hoS9_branch == 0 ) {
	cout << "Branch mus_e_hoS9 does not exist." << endl;
	}
	mus_etaErr_branch = 0;
	if (tree->GetAlias("mus_etaErr") != 0) {
		mus_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_etaErr"));
		mus_etaErr_branch->SetAddress(&mus_etaErr_);
	}
	if(mus_etaErr_branch == 0 ) {
	cout << "Branch mus_etaErr does not exist." << endl;
	}
	mus_gfit_chi2_branch = 0;
	if (tree->GetAlias("mus_gfit_chi2") != 0) {
		mus_gfit_chi2_branch = tree->GetBranch(tree->GetAlias("mus_gfit_chi2"));
		mus_gfit_chi2_branch->SetAddress(&mus_gfit_chi2_);
	}
	if(mus_gfit_chi2_branch == 0 ) {
	cout << "Branch mus_gfit_chi2 does not exist." << endl;
	}
	mus_gfit_ndof_branch = 0;
	if (tree->GetAlias("mus_gfit_ndof") != 0) {
		mus_gfit_ndof_branch = tree->GetBranch(tree->GetAlias("mus_gfit_ndof"));
		mus_gfit_ndof_branch->SetAddress(&mus_gfit_ndof_);
	}
	if(mus_gfit_ndof_branch == 0 ) {
	cout << "Branch mus_gfit_ndof does not exist." << endl;
	}
	mus_iso_branch = 0;
	if (tree->GetAlias("mus_iso") != 0) {
		mus_iso_branch = tree->GetBranch(tree->GetAlias("mus_iso"));
		mus_iso_branch->SetAddress(&mus_iso_);
	}
	if(mus_iso_branch == 0 ) {
	cout << "Branch mus_iso does not exist." << endl;
	}
	mus_iso03_emEt_branch = 0;
	if (tree->GetAlias("mus_iso03_emEt") != 0) {
		mus_iso03_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_emEt"));
		mus_iso03_emEt_branch->SetAddress(&mus_iso03_emEt_);
	}
	if(mus_iso03_emEt_branch == 0 ) {
	cout << "Branch mus_iso03_emEt does not exist." << endl;
	}
	mus_iso03_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hadEt") != 0) {
		mus_iso03_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hadEt"));
		mus_iso03_hadEt_branch->SetAddress(&mus_iso03_hadEt_);
	}
	if(mus_iso03_hadEt_branch == 0 ) {
	cout << "Branch mus_iso03_hadEt does not exist." << endl;
	}
	mus_iso03_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hoEt") != 0) {
		mus_iso03_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hoEt"));
		mus_iso03_hoEt_branch->SetAddress(&mus_iso03_hoEt_);
	}
	if(mus_iso03_hoEt_branch == 0 ) {
	cout << "Branch mus_iso03_hoEt does not exist." << endl;
	}
	mus_iso03_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso03_sumPt") != 0) {
		mus_iso03_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_sumPt"));
		mus_iso03_sumPt_branch->SetAddress(&mus_iso03_sumPt_);
	}
	if(mus_iso03_sumPt_branch == 0 ) {
	cout << "Branch mus_iso03_sumPt does not exist." << endl;
	}
	mus_iso05_emEt_branch = 0;
	if (tree->GetAlias("mus_iso05_emEt") != 0) {
		mus_iso05_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_emEt"));
		mus_iso05_emEt_branch->SetAddress(&mus_iso05_emEt_);
	}
	if(mus_iso05_emEt_branch == 0 ) {
	cout << "Branch mus_iso05_emEt does not exist." << endl;
	}
	mus_iso05_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hadEt") != 0) {
		mus_iso05_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hadEt"));
		mus_iso05_hadEt_branch->SetAddress(&mus_iso05_hadEt_);
	}
	if(mus_iso05_hadEt_branch == 0 ) {
	cout << "Branch mus_iso05_hadEt does not exist." << endl;
	}
	mus_iso05_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hoEt") != 0) {
		mus_iso05_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hoEt"));
		mus_iso05_hoEt_branch->SetAddress(&mus_iso05_hoEt_);
	}
	if(mus_iso05_hoEt_branch == 0 ) {
	cout << "Branch mus_iso05_hoEt does not exist." << endl;
	}
	mus_iso05_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso05_sumPt") != 0) {
		mus_iso05_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_sumPt"));
		mus_iso05_sumPt_branch->SetAddress(&mus_iso05_sumPt_);
	}
	if(mus_iso05_sumPt_branch == 0 ) {
	cout << "Branch mus_iso05_sumPt does not exist." << endl;
	}
	mus_ndof_branch = 0;
	if (tree->GetAlias("mus_ndof") != 0) {
		mus_ndof_branch = tree->GetBranch(tree->GetAlias("mus_ndof"));
		mus_ndof_branch->SetAddress(&mus_ndof_);
	}
	if(mus_ndof_branch == 0 ) {
	cout << "Branch mus_ndof does not exist." << endl;
	}
	mus_outerEta_branch = 0;
	if (tree->GetAlias("mus_outerEta") != 0) {
		mus_outerEta_branch = tree->GetBranch(tree->GetAlias("mus_outerEta"));
		mus_outerEta_branch->SetAddress(&mus_outerEta_);
	}
	if(mus_outerEta_branch == 0 ) {
	cout << "Branch mus_outerEta does not exist." << endl;
	}
	mus_outerPhi_branch = 0;
	if (tree->GetAlias("mus_outerPhi") != 0) {
		mus_outerPhi_branch = tree->GetBranch(tree->GetAlias("mus_outerPhi"));
		mus_outerPhi_branch->SetAddress(&mus_outerPhi_);
	}
	if(mus_outerPhi_branch == 0 ) {
	cout << "Branch mus_outerPhi does not exist." << endl;
	}
	mus_phiErr_branch = 0;
	if (tree->GetAlias("mus_phiErr") != 0) {
		mus_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_phiErr"));
		mus_phiErr_branch->SetAddress(&mus_phiErr_);
	}
	if(mus_phiErr_branch == 0 ) {
	cout << "Branch mus_phiErr does not exist." << endl;
	}
	mus_ptErr_branch = 0;
	if (tree->GetAlias("mus_ptErr") != 0) {
		mus_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_ptErr"));
		mus_ptErr_branch->SetAddress(&mus_ptErr_);
	}
	if(mus_ptErr_branch == 0 ) {
	cout << "Branch mus_ptErr does not exist." << endl;
	}
	mus_qoverp_branch = 0;
	if (tree->GetAlias("mus_qoverp") != 0) {
		mus_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_qoverp"));
		mus_qoverp_branch->SetAddress(&mus_qoverp_);
	}
	if(mus_qoverp_branch == 0 ) {
	cout << "Branch mus_qoverp does not exist." << endl;
	}
	mus_qoverpError_branch = 0;
	if (tree->GetAlias("mus_qoverpError") != 0) {
		mus_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_qoverpError"));
		mus_qoverpError_branch->SetAddress(&mus_qoverpError_);
	}
	if(mus_qoverpError_branch == 0 ) {
	cout << "Branch mus_qoverpError does not exist." << endl;
	}
	mus_vertexphi_branch = 0;
	if (tree->GetAlias("mus_vertexphi") != 0) {
		mus_vertexphi_branch = tree->GetBranch(tree->GetAlias("mus_vertexphi"));
		mus_vertexphi_branch->SetAddress(&mus_vertexphi_);
	}
	if(mus_vertexphi_branch == 0 ) {
	cout << "Branch mus_vertexphi does not exist." << endl;
	}
	mus_z0_branch = 0;
	if (tree->GetAlias("mus_z0") != 0) {
		mus_z0_branch = tree->GetBranch(tree->GetAlias("mus_z0"));
		mus_z0_branch->SetAddress(&mus_z0_);
	}
	if(mus_z0_branch == 0 ) {
	cout << "Branch mus_z0 does not exist." << endl;
	}
	mus_z0Err_branch = 0;
	if (tree->GetAlias("mus_z0Err") != 0) {
		mus_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_z0Err"));
		mus_z0Err_branch->SetAddress(&mus_z0Err_);
	}
	if(mus_z0Err_branch == 0 ) {
	cout << "Branch mus_z0Err does not exist." << endl;
	}
	mus_z0corr_branch = 0;
	if (tree->GetAlias("mus_z0corr") != 0) {
		mus_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_z0corr"));
		mus_z0corr_branch->SetAddress(&mus_z0corr_);
	}
	if(mus_z0corr_branch == 0 ) {
	cout << "Branch mus_z0corr does not exist." << endl;
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
	trks_outerEt_branch = 0;
	if (tree->GetAlias("trks_outerEt") != 0) {
		trks_outerEt_branch = tree->GetBranch(tree->GetAlias("trks_outerEt"));
		trks_outerEt_branch->SetAddress(&trks_outerEt_);
	}
	if(trks_outerEt_branch == 0 ) {
	cout << "Branch trks_outerEt does not exist." << endl;
	}
	trks_outerEta_branch = 0;
	if (tree->GetAlias("trks_outerEta") != 0) {
		trks_outerEta_branch = tree->GetBranch(tree->GetAlias("trks_outerEta"));
		trks_outerEta_branch->SetAddress(&trks_outerEta_);
	}
	if(trks_outerEta_branch == 0 ) {
	cout << "Branch trks_outerEta does not exist." << endl;
	}
	trks_outerPhi_branch = 0;
	if (tree->GetAlias("trks_outerPhi") != 0) {
		trks_outerPhi_branch = tree->GetBranch(tree->GetAlias("trks_outerPhi"));
		trks_outerPhi_branch->SetAddress(&trks_outerPhi_);
	}
	if(trks_outerPhi_branch == 0 ) {
	cout << "Branch trks_outerPhi does not exist." << endl;
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
	trks_tkIso_branch = 0;
	if (tree->GetAlias("trks_tkIso") != 0) {
		trks_tkIso_branch = tree->GetBranch(tree->GetAlias("trks_tkIso"));
		trks_tkIso_branch->SetAddress(&trks_tkIso_);
	}
	if(trks_tkIso_branch == 0 ) {
	cout << "Branch trks_tkIso does not exist." << endl;
	}
	trks_vertexphi_branch = 0;
	if (tree->GetAlias("trks_vertexphi") != 0) {
		trks_vertexphi_branch = tree->GetBranch(tree->GetAlias("trks_vertexphi"));
		trks_vertexphi_branch->SetAddress(&trks_vertexphi_);
	}
	if(trks_vertexphi_branch == 0 ) {
	cout << "Branch trks_vertexphi does not exist." << endl;
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
	trks_elsdr_branch = 0;
	if (tree->GetAlias("trks_elsdr") != 0) {
		trks_elsdr_branch = tree->GetBranch(tree->GetAlias("trks_elsdr"));
		trks_elsdr_branch->SetAddress(&trks_elsdr_);
	}
	if(trks_elsdr_branch == 0 ) {
	cout << "Branch trks_elsdr does not exist." << endl;
	}
	trks_elsshFrac_branch = 0;
	if (tree->GetAlias("trks_elsshFrac") != 0) {
		trks_elsshFrac_branch = tree->GetBranch(tree->GetAlias("trks_elsshFrac"));
		trks_elsshFrac_branch->SetAddress(&trks_elsshFrac_);
	}
	if(trks_elsshFrac_branch == 0 ) {
	cout << "Branch trks_elsshFrac does not exist." << endl;
	}
	trk_musdr_branch = 0;
	if (tree->GetAlias("trk_musdr") != 0) {
		trk_musdr_branch = tree->GetBranch(tree->GetAlias("trk_musdr"));
		trk_musdr_branch->SetAddress(&trk_musdr_);
	}
	if(trk_musdr_branch == 0 ) {
	cout << "Branch trk_musdr does not exist." << endl;
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
	trks_residualX_branch = 0;
	if (tree->GetAlias("trks_residualX") != 0) {
		trks_residualX_branch = tree->GetBranch(tree->GetAlias("trks_residualX"));
		trks_residualX_branch->SetAddress(&trks_residualX_);
	}
	if(trks_residualX_branch == 0 ) {
	cout << "Branch trks_residualX does not exist." << endl;
	}
	trks_residualY_branch = 0;
	if (tree->GetAlias("trks_residualY") != 0) {
		trks_residualY_branch = tree->GetBranch(tree->GetAlias("trks_residualY"));
		trks_residualY_branch->SetAddress(&trks_residualY_);
	}
	if(trks_residualY_branch == 0 ) {
	cout << "Branch trks_residualY does not exist." << endl;
	}
	evt_HLT1_branch = 0;
	if (tree->GetAlias("evt_HLT1") != 0) {
		evt_HLT1_branch = tree->GetBranch(tree->GetAlias("evt_HLT1"));
		evt_HLT1_branch->SetAddress(&evt_HLT1_);
	}
	if(evt_HLT1_branch == 0 ) {
	cout << "Branch evt_HLT1 does not exist." << endl;
	}
	evt_HLT2_branch = 0;
	if (tree->GetAlias("evt_HLT2") != 0) {
		evt_HLT2_branch = tree->GetBranch(tree->GetAlias("evt_HLT2"));
		evt_HLT2_branch->SetAddress(&evt_HLT2_);
	}
	if(evt_HLT2_branch == 0 ) {
	cout << "Branch evt_HLT2 does not exist." << endl;
	}
	evt_HLT3_branch = 0;
	if (tree->GetAlias("evt_HLT3") != 0) {
		evt_HLT3_branch = tree->GetBranch(tree->GetAlias("evt_HLT3"));
		evt_HLT3_branch->SetAddress(&evt_HLT3_);
	}
	if(evt_HLT3_branch == 0 ) {
	cout << "Branch evt_HLT3 does not exist." << endl;
	}
	evt_HLT4_branch = 0;
	if (tree->GetAlias("evt_HLT4") != 0) {
		evt_HLT4_branch = tree->GetBranch(tree->GetAlias("evt_HLT4"));
		evt_HLT4_branch->SetAddress(&evt_HLT4_);
	}
	if(evt_HLT4_branch == 0 ) {
	cout << "Branch evt_HLT4 does not exist." << endl;
	}
	evt_HLT5_branch = 0;
	if (tree->GetAlias("evt_HLT5") != 0) {
		evt_HLT5_branch = tree->GetBranch(tree->GetAlias("evt_HLT5"));
		evt_HLT5_branch->SetAddress(&evt_HLT5_);
	}
	if(evt_HLT5_branch == 0 ) {
	cout << "Branch evt_HLT5 does not exist." << endl;
	}
	evt_HLT6_branch = 0;
	if (tree->GetAlias("evt_HLT6") != 0) {
		evt_HLT6_branch = tree->GetBranch(tree->GetAlias("evt_HLT6"));
		evt_HLT6_branch->SetAddress(&evt_HLT6_);
	}
	if(evt_HLT6_branch == 0 ) {
	cout << "Branch evt_HLT6 does not exist." << endl;
	}
	evt_HLT7_branch = 0;
	if (tree->GetAlias("evt_HLT7") != 0) {
		evt_HLT7_branch = tree->GetBranch(tree->GetAlias("evt_HLT7"));
		evt_HLT7_branch->SetAddress(&evt_HLT7_);
	}
	if(evt_HLT7_branch == 0 ) {
	cout << "Branch evt_HLT7 does not exist." << endl;
	}
	evt_HLT8_branch = 0;
	if (tree->GetAlias("evt_HLT8") != 0) {
		evt_HLT8_branch = tree->GetBranch(tree->GetAlias("evt_HLT8"));
		evt_HLT8_branch->SetAddress(&evt_HLT8_);
	}
	if(evt_HLT8_branch == 0 ) {
	cout << "Branch evt_HLT8 does not exist." << endl;
	}
	evt_L1_1_branch = 0;
	if (tree->GetAlias("evt_L1_1") != 0) {
		evt_L1_1_branch = tree->GetBranch(tree->GetAlias("evt_L1_1"));
		evt_L1_1_branch->SetAddress(&evt_L1_1_);
	}
	if(evt_L1_1_branch == 0 ) {
	cout << "Branch evt_L1_1 does not exist." << endl;
	}
	evt_L1_2_branch = 0;
	if (tree->GetAlias("evt_L1_2") != 0) {
		evt_L1_2_branch = tree->GetBranch(tree->GetAlias("evt_L1_2"));
		evt_L1_2_branch->SetAddress(&evt_L1_2_);
	}
	if(evt_L1_2_branch == 0 ) {
	cout << "Branch evt_L1_2 does not exist." << endl;
	}
	evt_L1_3_branch = 0;
	if (tree->GetAlias("evt_L1_3") != 0) {
		evt_L1_3_branch = tree->GetBranch(tree->GetAlias("evt_L1_3"));
		evt_L1_3_branch->SetAddress(&evt_L1_3_);
	}
	if(evt_L1_3_branch == 0 ) {
	cout << "Branch evt_L1_3 does not exist." << endl;
	}
	evt_L1_4_branch = 0;
	if (tree->GetAlias("evt_L1_4") != 0) {
		evt_L1_4_branch = tree->GetBranch(tree->GetAlias("evt_L1_4"));
		evt_L1_4_branch->SetAddress(&evt_L1_4_);
	}
	if(evt_L1_4_branch == 0 ) {
	cout << "Branch evt_L1_4 does not exist." << endl;
	}
	els_mc3_id_branch = 0;
	if (tree->GetAlias("els_mc3_id") != 0) {
		els_mc3_id_branch = tree->GetBranch(tree->GetAlias("els_mc3_id"));
		els_mc3_id_branch->SetAddress(&els_mc3_id_);
	}
	if(els_mc3_id_branch == 0 ) {
	cout << "Branch els_mc3_id does not exist." << endl;
	}
	els_mc3idx_branch = 0;
	if (tree->GetAlias("els_mc3idx") != 0) {
		els_mc3idx_branch = tree->GetBranch(tree->GetAlias("els_mc3idx"));
		els_mc3idx_branch->SetAddress(&els_mc3idx_);
	}
	if(els_mc3idx_branch == 0 ) {
	cout << "Branch els_mc3idx does not exist." << endl;
	}
	els_mc3_motherid_branch = 0;
	if (tree->GetAlias("els_mc3_motherid") != 0) {
		els_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc3_motherid"));
		els_mc3_motherid_branch->SetAddress(&els_mc3_motherid_);
	}
	if(els_mc3_motherid_branch == 0 ) {
	cout << "Branch els_mc3_motherid does not exist." << endl;
	}
	els_mc_id_branch = 0;
	if (tree->GetAlias("els_mc_id") != 0) {
		els_mc_id_branch = tree->GetBranch(tree->GetAlias("els_mc_id"));
		els_mc_id_branch->SetAddress(&els_mc_id_);
	}
	if(els_mc_id_branch == 0 ) {
	cout << "Branch els_mc_id does not exist." << endl;
	}
	els_mcidx_branch = 0;
	if (tree->GetAlias("els_mcidx") != 0) {
		els_mcidx_branch = tree->GetBranch(tree->GetAlias("els_mcidx"));
		els_mcidx_branch->SetAddress(&els_mcidx_);
	}
	if(els_mcidx_branch == 0 ) {
	cout << "Branch els_mcidx does not exist." << endl;
	}
	els_mc_motherid_branch = 0;
	if (tree->GetAlias("els_mc_motherid") != 0) {
		els_mc_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc_motherid"));
		els_mc_motherid_branch->SetAddress(&els_mc_motherid_);
	}
	if(els_mc_motherid_branch == 0 ) {
	cout << "Branch els_mc_motherid does not exist." << endl;
	}
	jets_mc_id_branch = 0;
	if (tree->GetAlias("jets_mc_id") != 0) {
		jets_mc_id_branch = tree->GetBranch(tree->GetAlias("jets_mc_id"));
		jets_mc_id_branch->SetAddress(&jets_mc_id_);
	}
	if(jets_mc_id_branch == 0 ) {
	cout << "Branch jets_mc_id does not exist." << endl;
	}
	mus_mc3_id_branch = 0;
	if (tree->GetAlias("mus_mc3_id") != 0) {
		mus_mc3_id_branch = tree->GetBranch(tree->GetAlias("mus_mc3_id"));
		mus_mc3_id_branch->SetAddress(&mus_mc3_id_);
	}
	if(mus_mc3_id_branch == 0 ) {
	cout << "Branch mus_mc3_id does not exist." << endl;
	}
	mus_mc3idx_branch = 0;
	if (tree->GetAlias("mus_mc3idx") != 0) {
		mus_mc3idx_branch = tree->GetBranch(tree->GetAlias("mus_mc3idx"));
		mus_mc3idx_branch->SetAddress(&mus_mc3idx_);
	}
	if(mus_mc3idx_branch == 0 ) {
	cout << "Branch mus_mc3idx does not exist." << endl;
	}
	mus_mc3_motherid_branch = 0;
	if (tree->GetAlias("mus_mc3_motherid") != 0) {
		mus_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motherid"));
		mus_mc3_motherid_branch->SetAddress(&mus_mc3_motherid_);
	}
	if(mus_mc3_motherid_branch == 0 ) {
	cout << "Branch mus_mc3_motherid does not exist." << endl;
	}
	mus_mc_id_branch = 0;
	if (tree->GetAlias("mus_mc_id") != 0) {
		mus_mc_id_branch = tree->GetBranch(tree->GetAlias("mus_mc_id"));
		mus_mc_id_branch->SetAddress(&mus_mc_id_);
	}
	if(mus_mc_id_branch == 0 ) {
	cout << "Branch mus_mc_id does not exist." << endl;
	}
	mus_mcidx_branch = 0;
	if (tree->GetAlias("mus_mcidx") != 0) {
		mus_mcidx_branch = tree->GetBranch(tree->GetAlias("mus_mcidx"));
		mus_mcidx_branch->SetAddress(&mus_mcidx_);
	}
	if(mus_mcidx_branch == 0 ) {
	cout << "Branch mus_mcidx does not exist." << endl;
	}
	mus_mc_motherid_branch = 0;
	if (tree->GetAlias("mus_mc_motherid") != 0) {
		mus_mc_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherid"));
		mus_mc_motherid_branch->SetAddress(&mus_mc_motherid_);
	}
	if(mus_mc_motherid_branch == 0 ) {
	cout << "Branch mus_mc_motherid does not exist." << endl;
	}
	trk_mc3_id_branch = 0;
	if (tree->GetAlias("trk_mc3_id") != 0) {
		trk_mc3_id_branch = tree->GetBranch(tree->GetAlias("trk_mc3_id"));
		trk_mc3_id_branch->SetAddress(&trk_mc3_id_);
	}
	if(trk_mc3_id_branch == 0 ) {
	cout << "Branch trk_mc3_id does not exist." << endl;
	}
	trk_mc3idx_branch = 0;
	if (tree->GetAlias("trk_mc3idx") != 0) {
		trk_mc3idx_branch = tree->GetBranch(tree->GetAlias("trk_mc3idx"));
		trk_mc3idx_branch->SetAddress(&trk_mc3idx_);
	}
	if(trk_mc3idx_branch == 0 ) {
	cout << "Branch trk_mc3idx does not exist." << endl;
	}
	trk_mc3_motherid_branch = 0;
	if (tree->GetAlias("trk_mc3_motherid") != 0) {
		trk_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("trk_mc3_motherid"));
		trk_mc3_motherid_branch->SetAddress(&trk_mc3_motherid_);
	}
	if(trk_mc3_motherid_branch == 0 ) {
	cout << "Branch trk_mc3_motherid does not exist." << endl;
	}
	trk_mc_id_branch = 0;
	if (tree->GetAlias("trk_mc_id") != 0) {
		trk_mc_id_branch = tree->GetBranch(tree->GetAlias("trk_mc_id"));
		trk_mc_id_branch->SetAddress(&trk_mc_id_);
	}
	if(trk_mc_id_branch == 0 ) {
	cout << "Branch trk_mc_id does not exist." << endl;
	}
	trk_mcidx_branch = 0;
	if (tree->GetAlias("trk_mcidx") != 0) {
		trk_mcidx_branch = tree->GetBranch(tree->GetAlias("trk_mcidx"));
		trk_mcidx_branch->SetAddress(&trk_mcidx_);
	}
	if(trk_mcidx_branch == 0 ) {
	cout << "Branch trk_mcidx does not exist." << endl;
	}
	trk_mc_motherid_branch = 0;
	if (tree->GetAlias("trk_mc_motherid") != 0) {
		trk_mc_motherid_branch = tree->GetBranch(tree->GetAlias("trk_mc_motherid"));
		trk_mc_motherid_branch->SetAddress(&trk_mc_motherid_);
	}
	if(trk_mc_motherid_branch == 0 ) {
	cout << "Branch trk_mc_motherid does not exist." << endl;
	}
	els_conv_tkidx_branch = 0;
	if (tree->GetAlias("els_conv_tkidx") != 0) {
		els_conv_tkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_tkidx"));
		els_conv_tkidx_branch->SetAddress(&els_conv_tkidx_);
	}
	if(els_conv_tkidx_branch == 0 ) {
	cout << "Branch els_conv_tkidx does not exist." << endl;
	}
	trks_conv_tkidx_branch = 0;
	if (tree->GetAlias("trks_conv_tkidx") != 0) {
		trks_conv_tkidx_branch = tree->GetBranch(tree->GetAlias("trks_conv_tkidx"));
		trks_conv_tkidx_branch->SetAddress(&trks_conv_tkidx_);
	}
	if(trks_conv_tkidx_branch == 0 ) {
	cout << "Branch trks_conv_tkidx does not exist." << endl;
	}
	els_closestMuon_branch = 0;
	if (tree->GetAlias("els_closestMuon") != 0) {
		els_closestMuon_branch = tree->GetBranch(tree->GetAlias("els_closestMuon"));
		els_closestMuon_branch->SetAddress(&els_closestMuon_);
	}
	if(els_closestMuon_branch == 0 ) {
	cout << "Branch els_closestMuon does not exist." << endl;
	}
	els_trkidx_branch = 0;
	if (tree->GetAlias("els_trkidx") != 0) {
		els_trkidx_branch = tree->GetBranch(tree->GetAlias("els_trkidx"));
		els_trkidx_branch->SetAddress(&els_trkidx_);
	}
	if(els_trkidx_branch == 0 ) {
	cout << "Branch els_trkidx does not exist." << endl;
	}
	els_category_branch = 0;
	if (tree->GetAlias("els_category") != 0) {
		els_category_branch = tree->GetBranch(tree->GetAlias("els_category"));
		els_category_branch->SetAddress(&els_category_);
	}
	if(els_category_branch == 0 ) {
	cout << "Branch els_category does not exist." << endl;
	}
	els_categoryold_branch = 0;
	if (tree->GetAlias("els_categoryold") != 0) {
		els_categoryold_branch = tree->GetBranch(tree->GetAlias("els_categoryold"));
		els_categoryold_branch->SetAddress(&els_categoryold_);
	}
	if(els_categoryold_branch == 0 ) {
	cout << "Branch els_categoryold does not exist." << endl;
	}
	els_charge_branch = 0;
	if (tree->GetAlias("els_charge") != 0) {
		els_charge_branch = tree->GetBranch(tree->GetAlias("els_charge"));
		els_charge_branch->SetAddress(&els_charge_);
	}
	if(els_charge_branch == 0 ) {
	cout << "Branch els_charge does not exist." << endl;
	}
	els_class_branch = 0;
	if (tree->GetAlias("els_class") != 0) {
		els_class_branch = tree->GetBranch(tree->GetAlias("els_class"));
		els_class_branch->SetAddress(&els_class_);
	}
	if(els_class_branch == 0 ) {
	cout << "Branch els_class does not exist." << endl;
	}
	els_layer1_det_branch = 0;
	if (tree->GetAlias("els_layer1_det") != 0) {
		els_layer1_det_branch = tree->GetBranch(tree->GetAlias("els_layer1_det"));
		els_layer1_det_branch->SetAddress(&els_layer1_det_);
	}
	if(els_layer1_det_branch == 0 ) {
	cout << "Branch els_layer1_det does not exist." << endl;
	}
	els_layer1_layer_branch = 0;
	if (tree->GetAlias("els_layer1_layer") != 0) {
		els_layer1_layer_branch = tree->GetBranch(tree->GetAlias("els_layer1_layer"));
		els_layer1_layer_branch->SetAddress(&els_layer1_layer_);
	}
	if(els_layer1_layer_branch == 0 ) {
	cout << "Branch els_layer1_layer does not exist." << endl;
	}
	els_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("els_layer1_sizerphi") != 0) {
		els_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerphi"));
		els_layer1_sizerphi_branch->SetAddress(&els_layer1_sizerphi_);
	}
	if(els_layer1_sizerphi_branch == 0 ) {
	cout << "Branch els_layer1_sizerphi does not exist." << endl;
	}
	els_layer1_sizerz_branch = 0;
	if (tree->GetAlias("els_layer1_sizerz") != 0) {
		els_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerz"));
		els_layer1_sizerz_branch->SetAddress(&els_layer1_sizerz_);
	}
	if(els_layer1_sizerz_branch == 0 ) {
	cout << "Branch els_layer1_sizerz does not exist." << endl;
	}
	els_looseId_branch = 0;
	if (tree->GetAlias("els_looseId") != 0) {
		els_looseId_branch = tree->GetBranch(tree->GetAlias("els_looseId"));
		els_looseId_branch->SetAddress(&els_looseId_);
	}
	if(els_looseId_branch == 0 ) {
	cout << "Branch els_looseId does not exist." << endl;
	}
	els_lostHits_branch = 0;
	if (tree->GetAlias("els_lostHits") != 0) {
		els_lostHits_branch = tree->GetBranch(tree->GetAlias("els_lostHits"));
		els_lostHits_branch->SetAddress(&els_lostHits_);
	}
	if(els_lostHits_branch == 0 ) {
	cout << "Branch els_lostHits does not exist." << endl;
	}
	els_lost_pixelhits_branch = 0;
	if (tree->GetAlias("els_lost_pixelhits") != 0) {
		els_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_lost_pixelhits"));
		els_lost_pixelhits_branch->SetAddress(&els_lost_pixelhits_);
	}
	if(els_lost_pixelhits_branch == 0 ) {
	cout << "Branch els_lost_pixelhits does not exist." << endl;
	}
	els_nSeed_branch = 0;
	if (tree->GetAlias("els_nSeed") != 0) {
		els_nSeed_branch = tree->GetBranch(tree->GetAlias("els_nSeed"));
		els_nSeed_branch->SetAddress(&els_nSeed_);
	}
	if(els_nSeed_branch == 0 ) {
	cout << "Branch els_nSeed does not exist." << endl;
	}
	els_n_inner_layers_branch = 0;
	if (tree->GetAlias("els_n_inner_layers") != 0) {
		els_n_inner_layers_branch = tree->GetBranch(tree->GetAlias("els_n_inner_layers"));
		els_n_inner_layers_branch->SetAddress(&els_n_inner_layers_);
	}
	if(els_n_inner_layers_branch == 0 ) {
	cout << "Branch els_n_inner_layers does not exist." << endl;
	}
	els_n_outer_layers_branch = 0;
	if (tree->GetAlias("els_n_outer_layers") != 0) {
		els_n_outer_layers_branch = tree->GetBranch(tree->GetAlias("els_n_outer_layers"));
		els_n_outer_layers_branch->SetAddress(&els_n_outer_layers_);
	}
	if(els_n_outer_layers_branch == 0 ) {
	cout << "Branch els_n_outer_layers does not exist." << endl;
	}
	els_pass3looseId_branch = 0;
	if (tree->GetAlias("els_pass3looseId") != 0) {
		els_pass3looseId_branch = tree->GetBranch(tree->GetAlias("els_pass3looseId"));
		els_pass3looseId_branch->SetAddress(&els_pass3looseId_);
	}
	if(els_pass3looseId_branch == 0 ) {
	cout << "Branch els_pass3looseId does not exist." << endl;
	}
	els_pass3simpleId_branch = 0;
	if (tree->GetAlias("els_pass3simpleId") != 0) {
		els_pass3simpleId_branch = tree->GetBranch(tree->GetAlias("els_pass3simpleId"));
		els_pass3simpleId_branch->SetAddress(&els_pass3simpleId_);
	}
	if(els_pass3simpleId_branch == 0 ) {
	cout << "Branch els_pass3simpleId does not exist." << endl;
	}
	els_pass3tightId_branch = 0;
	if (tree->GetAlias("els_pass3tightId") != 0) {
		els_pass3tightId_branch = tree->GetBranch(tree->GetAlias("els_pass3tightId"));
		els_pass3tightId_branch->SetAddress(&els_pass3tightId_);
	}
	if(els_pass3tightId_branch == 0 ) {
	cout << "Branch els_pass3tightId does not exist." << endl;
	}
	els_robustId_branch = 0;
	if (tree->GetAlias("els_robustId") != 0) {
		els_robustId_branch = tree->GetBranch(tree->GetAlias("els_robustId"));
		els_robustId_branch->SetAddress(&els_robustId_);
	}
	if(els_robustId_branch == 0 ) {
	cout << "Branch els_robustId does not exist." << endl;
	}
	els_simpleIdPlus_branch = 0;
	if (tree->GetAlias("els_simpleIdPlus") != 0) {
		els_simpleIdPlus_branch = tree->GetBranch(tree->GetAlias("els_simpleIdPlus"));
		els_simpleIdPlus_branch->SetAddress(&els_simpleIdPlus_);
	}
	if(els_simpleIdPlus_branch == 0 ) {
	cout << "Branch els_simpleIdPlus does not exist." << endl;
	}
	els_tightId_branch = 0;
	if (tree->GetAlias("els_tightId") != 0) {
		els_tightId_branch = tree->GetBranch(tree->GetAlias("els_tightId"));
		els_tightId_branch->SetAddress(&els_tightId_);
	}
	if(els_tightId_branch == 0 ) {
	cout << "Branch els_tightId does not exist." << endl;
	}
	els_tightId22XMaxMatteo_branch = 0;
	if (tree->GetAlias("els_tightId22XMaxMatteo") != 0) {
		els_tightId22XMaxMatteo_branch = tree->GetBranch(tree->GetAlias("els_tightId22XMaxMatteo"));
		els_tightId22XMaxMatteo_branch->SetAddress(&els_tightId22XMaxMatteo_);
	}
	if(els_tightId22XMaxMatteo_branch == 0 ) {
	cout << "Branch els_tightId22XMaxMatteo does not exist." << endl;
	}
	els_tightId22XMinMatteo_branch = 0;
	if (tree->GetAlias("els_tightId22XMinMatteo") != 0) {
		els_tightId22XMinMatteo_branch = tree->GetBranch(tree->GetAlias("els_tightId22XMinMatteo"));
		els_tightId22XMinMatteo_branch->SetAddress(&els_tightId22XMinMatteo_);
	}
	if(els_tightId22XMinMatteo_branch == 0 ) {
	cout << "Branch els_tightId22XMinMatteo does not exist." << endl;
	}
	els_validHits_branch = 0;
	if (tree->GetAlias("els_validHits") != 0) {
		els_validHits_branch = tree->GetBranch(tree->GetAlias("els_validHits"));
		els_validHits_branch->SetAddress(&els_validHits_);
	}
	if(els_validHits_branch == 0 ) {
	cout << "Branch els_validHits does not exist." << endl;
	}
	els_valid_pixelhits_branch = 0;
	if (tree->GetAlias("els_valid_pixelhits") != 0) {
		els_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_valid_pixelhits"));
		els_valid_pixelhits_branch->SetAddress(&els_valid_pixelhits_);
	}
	if(els_valid_pixelhits_branch == 0 ) {
	cout << "Branch els_valid_pixelhits does not exist." << endl;
	}
	genps_id_branch = 0;
	if (tree->GetAlias("genps_id") != 0) {
		genps_id_branch = tree->GetBranch(tree->GetAlias("genps_id"));
		genps_id_branch->SetAddress(&genps_id_);
	}
	if(genps_id_branch == 0 ) {
	cout << "Branch genps_id does not exist." << endl;
	}
	genps_id_mother_branch = 0;
	if (tree->GetAlias("genps_id_mother") != 0) {
		genps_id_mother_branch = tree->GetBranch(tree->GetAlias("genps_id_mother"));
		genps_id_mother_branch->SetAddress(&genps_id_mother_);
	}
	if(genps_id_mother_branch == 0 ) {
	cout << "Branch genps_id_mother does not exist." << endl;
	}
	genps_lepdaughter_id_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_id") != 0) {
		genps_lepdaughter_id_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_id"));
		genps_lepdaughter_id_branch->SetAddress(&genps_lepdaughter_id_);
	}
	if(genps_lepdaughter_id_branch == 0 ) {
	cout << "Branch genps_lepdaughter_id does not exist." << endl;
	}
	genps_lepdaughter_idx_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_idx") != 0) {
		genps_lepdaughter_idx_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_idx"));
		genps_lepdaughter_idx_branch->SetAddress(&genps_lepdaughter_idx_);
	}
	if(genps_lepdaughter_idx_branch == 0 ) {
	cout << "Branch genps_lepdaughter_idx does not exist." << endl;
	}
	genps_status_branch = 0;
	if (tree->GetAlias("genps_status") != 0) {
		genps_status_branch = tree->GetBranch(tree->GetAlias("genps_status"));
		genps_status_branch->SetAddress(&genps_status_);
	}
	if(genps_status_branch == 0 ) {
	cout << "Branch genps_status does not exist." << endl;
	}
	mus_charge_branch = 0;
	if (tree->GetAlias("mus_charge") != 0) {
		mus_charge_branch = tree->GetBranch(tree->GetAlias("mus_charge"));
		mus_charge_branch->SetAddress(&mus_charge_);
	}
	if(mus_charge_branch == 0 ) {
	cout << "Branch mus_charge does not exist." << endl;
	}
	mus_gfit_validHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validHits") != 0) {
		mus_gfit_validHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validHits"));
		mus_gfit_validHits_branch->SetAddress(&mus_gfit_validHits_);
	}
	if(mus_gfit_validHits_branch == 0 ) {
	cout << "Branch mus_gfit_validHits does not exist." << endl;
	}
	mus_goodmask_branch = 0;
	if (tree->GetAlias("mus_goodmask") != 0) {
		mus_goodmask_branch = tree->GetBranch(tree->GetAlias("mus_goodmask"));
		mus_goodmask_branch->SetAddress(&mus_goodmask_);
	}
	if(mus_goodmask_branch == 0 ) {
	cout << "Branch mus_goodmask does not exist." << endl;
	}
	mus_iso03_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso03_ntrk") != 0) {
		mus_iso03_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso03_ntrk"));
		mus_iso03_ntrk_branch->SetAddress(&mus_iso03_ntrk_);
	}
	if(mus_iso03_ntrk_branch == 0 ) {
	cout << "Branch mus_iso03_ntrk does not exist." << endl;
	}
	mus_iso05_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso05_ntrk") != 0) {
		mus_iso05_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso05_ntrk"));
		mus_iso05_ntrk_branch->SetAddress(&mus_iso05_ntrk_);
	}
	if(mus_iso05_ntrk_branch == 0 ) {
	cout << "Branch mus_iso05_ntrk does not exist." << endl;
	}
	mus_lostHits_branch = 0;
	if (tree->GetAlias("mus_lostHits") != 0) {
		mus_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_lostHits"));
		mus_lostHits_branch->SetAddress(&mus_lostHits_);
	}
	if(mus_lostHits_branch == 0 ) {
	cout << "Branch mus_lostHits does not exist." << endl;
	}
	mus_nmatches_branch = 0;
	if (tree->GetAlias("mus_nmatches") != 0) {
		mus_nmatches_branch = tree->GetBranch(tree->GetAlias("mus_nmatches"));
		mus_nmatches_branch->SetAddress(&mus_nmatches_);
	}
	if(mus_nmatches_branch == 0 ) {
	cout << "Branch mus_nmatches does not exist." << endl;
	}
	mus_pid_TM2DCompatibilityLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityLoose") != 0) {
		mus_pid_TM2DCompatibilityLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityLoose"));
		mus_pid_TM2DCompatibilityLoose_branch->SetAddress(&mus_pid_TM2DCompatibilityLoose_);
	}
	if(mus_pid_TM2DCompatibilityLoose_branch == 0 ) {
	cout << "Branch mus_pid_TM2DCompatibilityLoose does not exist." << endl;
	}
	mus_pid_TM2DCompatibilityTight_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityTight") != 0) {
		mus_pid_TM2DCompatibilityTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityTight"));
		mus_pid_TM2DCompatibilityTight_branch->SetAddress(&mus_pid_TM2DCompatibilityTight_);
	}
	if(mus_pid_TM2DCompatibilityTight_branch == 0 ) {
	cout << "Branch mus_pid_TM2DCompatibilityTight does not exist." << endl;
	}
	mus_pid_TMLastStationLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationLoose") != 0) {
		mus_pid_TMLastStationLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationLoose"));
		mus_pid_TMLastStationLoose_branch->SetAddress(&mus_pid_TMLastStationLoose_);
	}
	if(mus_pid_TMLastStationLoose_branch == 0 ) {
	cout << "Branch mus_pid_TMLastStationLoose does not exist." << endl;
	}
	mus_pid_TMLastStationTight_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationTight") != 0) {
		mus_pid_TMLastStationTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationTight"));
		mus_pid_TMLastStationTight_branch->SetAddress(&mus_pid_TMLastStationTight_);
	}
	if(mus_pid_TMLastStationTight_branch == 0 ) {
	cout << "Branch mus_pid_TMLastStationTight does not exist." << endl;
	}
	mus_trk_charge_branch = 0;
	if (tree->GetAlias("mus_trk_charge") != 0) {
		mus_trk_charge_branch = tree->GetBranch(tree->GetAlias("mus_trk_charge"));
		mus_trk_charge_branch->SetAddress(&mus_trk_charge_);
	}
	if(mus_trk_charge_branch == 0 ) {
	cout << "Branch mus_trk_charge does not exist." << endl;
	}
	mus_trkrefkey_branch = 0;
	if (tree->GetAlias("mus_trkrefkey") != 0) {
		mus_trkrefkey_branch = tree->GetBranch(tree->GetAlias("mus_trkrefkey"));
		mus_trkrefkey_branch->SetAddress(&mus_trkrefkey_);
	}
	if(mus_trkrefkey_branch == 0 ) {
	cout << "Branch mus_trkrefkey does not exist." << endl;
	}
	mus_type_branch = 0;
	if (tree->GetAlias("mus_type") != 0) {
		mus_type_branch = tree->GetBranch(tree->GetAlias("mus_type"));
		mus_type_branch->SetAddress(&mus_type_);
	}
	if(mus_type_branch == 0 ) {
	cout << "Branch mus_type does not exist." << endl;
	}
	mus_validHits_branch = 0;
	if (tree->GetAlias("mus_validHits") != 0) {
		mus_validHits_branch = tree->GetBranch(tree->GetAlias("mus_validHits"));
		mus_validHits_branch->SetAddress(&mus_validHits_);
	}
	if(mus_validHits_branch == 0 ) {
	cout << "Branch mus_validHits does not exist." << endl;
	}
	trks_charge_branch = 0;
	if (tree->GetAlias("trks_charge") != 0) {
		trks_charge_branch = tree->GetBranch(tree->GetAlias("trks_charge"));
		trks_charge_branch->SetAddress(&trks_charge_);
	}
	if(trks_charge_branch == 0 ) {
	cout << "Branch trks_charge does not exist." << endl;
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
	trks_n_inner_layers_branch = 0;
	if (tree->GetAlias("trks_n_inner_layers") != 0) {
		trks_n_inner_layers_branch = tree->GetBranch(tree->GetAlias("trks_n_inner_layers"));
		trks_n_inner_layers_branch->SetAddress(&trks_n_inner_layers_);
	}
	if(trks_n_inner_layers_branch == 0 ) {
	cout << "Branch trks_n_inner_layers does not exist." << endl;
	}
	trks_n_outer_layers_branch = 0;
	if (tree->GetAlias("trks_n_outer_layers") != 0) {
		trks_n_outer_layers_branch = tree->GetBranch(tree->GetAlias("trks_n_outer_layers"));
		trks_n_outer_layers_branch->SetAddress(&trks_n_outer_layers_);
	}
	if(trks_n_outer_layers_branch == 0 ) {
	cout << "Branch trks_n_outer_layers does not exist." << endl;
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
	trks_elsidx_branch = 0;
	if (tree->GetAlias("trks_elsidx") != 0) {
		trks_elsidx_branch = tree->GetBranch(tree->GetAlias("trks_elsidx"));
		trks_elsidx_branch->SetAddress(&trks_elsidx_);
	}
	if(trks_elsidx_branch == 0 ) {
	cout << "Branch trks_elsidx does not exist." << endl;
	}
	trk_musidx_branch = 0;
	if (tree->GetAlias("trk_musidx") != 0) {
		trk_musidx_branch = tree->GetBranch(tree->GetAlias("trk_musidx"));
		trk_musidx_branch->SetAddress(&trk_musidx_);
	}
	if(trk_musidx_branch == 0 ) {
	cout << "Branch trk_musidx does not exist." << endl;
	}
	vtxs_hasRefittedTracks_branch = 0;
	if (tree->GetAlias("vtxs_hasRefittedTracks") != 0) {
		vtxs_hasRefittedTracks_branch = tree->GetBranch(tree->GetAlias("vtxs_hasRefittedTracks"));
		vtxs_hasRefittedTracks_branch->SetAddress(&vtxs_hasRefittedTracks_);
	}
	if(vtxs_hasRefittedTracks_branch == 0 ) {
	cout << "Branch vtxs_hasRefittedTracks does not exist." << endl;
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
	trks_hit_type_branch = 0;
	if (tree->GetAlias("trks_hit_type") != 0) {
		trks_hit_type_branch = tree->GetBranch(tree->GetAlias("trks_hit_type"));
		trks_hit_type_branch->SetAddress(&trks_hit_type_);
	}
	if(trks_hit_type_branch == 0 ) {
	cout << "Branch trks_hit_type does not exist." << endl;
	}
	evt_nels_branch = 0;
	if (tree->GetAlias("evt_nels") != 0) {
		evt_nels_branch = tree->GetBranch(tree->GetAlias("evt_nels"));
		evt_nels_branch->SetAddress(&evt_nels_);
	}
	if(evt_nels_branch == 0 ) {
	cout << "Branch evt_nels does not exist." << endl;
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
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	if(evt_nvtxs_branch == 0 ) {
	cout << "Branch evt_nvtxs does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		evt_CMS2tag_isLoaded = false;
		evt_dataset_isLoaded = false;
		evt_HLT_trigNames_isLoaded = false;
		evt_L1_trigNames_isLoaded = false;
		genps_pthat_isLoaded = false;
		vtxs_position_isLoaded = false;
		evt_bsp4_isLoaded = false;
		els_mc3_motherp4_isLoaded = false;
		els_mc3_p4_isLoaded = false;
		els_mc_motherp4_isLoaded = false;
		els_mc_p4_isLoaded = false;
		jets_mc_gp_p4_isLoaded = false;
		jets_mc_p4_isLoaded = false;
		mus_mc3_motherp4_isLoaded = false;
		mus_mc3_p4_isLoaded = false;
		mus_mc_motherp4_isLoaded = false;
		mus_mc_p4_isLoaded = false;
		trk_mc3p4_isLoaded = false;
		trk_mcp4_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
		els_vertex_p4_isLoaded = false;
		genps_lepdaughter_p4_isLoaded = false;
		genps_p4_isLoaded = false;
		genps_prod_vtx_isLoaded = false;
		mus_gfit_outerPos_p4_isLoaded = false;
		mus_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		mus_vertex_p4_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		trks_vertex_p4_isLoaded = false;
		els_inner_positionx_isLoaded = false;
		els_inner_positiony_isLoaded = false;
		els_inner_positionz_isLoaded = false;
		trks_inner_positionx_isLoaded = false;
		trks_inner_positiony_isLoaded = false;
		trks_inner_positionz_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_widthx_isLoaded = false;
		evt_bs_widthxErr_isLoaded = false;
		evt_bs_widthy_isLoaded = false;
		evt_bs_widthyErr_isLoaded = false;
		evt_bs_xErr_isLoaded = false;
		evt_bs_yErr_isLoaded = false;
		evt_bs_zErr_isLoaded = false;
		evt_bField_isLoaded = false;
		evt_kfactor_isLoaded = false;
		evt_weight_isLoaded = false;
		evt_xsec_excl_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metPhi_isLoaded = false;
		els_mc3dr_isLoaded = false;
		els_mcdr_isLoaded = false;
		jets_mcdr_isLoaded = false;
		jets_mc_emEnergy_isLoaded = false;
		jets_mc_gpdr_isLoaded = false;
		jets_mc_hadEnergy_isLoaded = false;
		jets_mc_invEnergy_isLoaded = false;
		jets_mc_otherEnergy_isLoaded = false;
		mus_mc3dr_isLoaded = false;
		mus_mcdr_isLoaded = false;
		trk_mc3dr_isLoaded = false;
		trk_mcdr_isLoaded = false;
		els_conv_dcot_isLoaded = false;
		els_conv_dist_isLoaded = false;
		trks_conv_dcot_isLoaded = false;
		trks_conv_dist_isLoaded = false;
		els_musdr_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_trkshFrac_isLoaded = false;
		els_chi2_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_d0corr_isLoaded = false;
		els_dEtaIn_isLoaded = false;
		els_dEtaOut_isLoaded = false;
		els_dPhiIn_isLoaded = false;
		els_dPhiInPhiOut_isLoaded = false;
		els_dPhiOut_isLoaded = false;
		els_e1x5_isLoaded = false;
		els_e2x5Max_isLoaded = false;
		els_e3x3_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_eMax_isLoaded = false;
		els_eOverPIn_isLoaded = false;
		els_eSC_isLoaded = false;
		els_eSCPresh_isLoaded = false;
		els_eSCRaw_isLoaded = false;
		els_eSeed_isLoaded = false;
		els_eSeedOverPIn_isLoaded = false;
		els_eSeedOverPOut_isLoaded = false;
		els_ecalIso03_isLoaded = false;
		els_ecalIso04_isLoaded = false;
		els_egamma_looseId_isLoaded = false;
		els_egamma_robustHighEnergy_isLoaded = false;
		els_egamma_robustLooseId_isLoaded = false;
		els_egamma_robustTightId_isLoaded = false;
		els_egamma_tightId_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_etaSC_isLoaded = false;
		els_fBrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_hcalIso03_isLoaded = false;
		els_hcalIso04_isLoaded = false;
		els_layer1_charge_isLoaded = false;
		els_ndof_isLoaded = false;
		els_outerEta_isLoaded = false;
		els_outerPhi_isLoaded = false;
		els_phiErr_isLoaded = false;
		els_phiSC_isLoaded = false;
		els_pin_isLoaded = false;
		els_pout_isLoaded = false;
		els_ptErr_isLoaded = false;
		els_sigmaEtaEta_isLoaded = false;
		els_sigmaIEtaIEta_isLoaded = false;
		els_sigmaIPhiIPhi_isLoaded = false;
		els_sigmaPhiPhi_isLoaded = false;
		els_tkIso03_isLoaded = false;
		els_tkIso04_isLoaded = false;
		els_vertexphi_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		els_z0corr_isLoaded = false;
		mus_caloCompatibility_isLoaded = false;
		mus_chi2_isLoaded = false;
		mus_d0_isLoaded = false;
		mus_d0Err_isLoaded = false;
		mus_d0corr_isLoaded = false;
		mus_e_em_isLoaded = false;
		mus_e_emS9_isLoaded = false;
		mus_e_had_isLoaded = false;
		mus_e_hadS9_isLoaded = false;
		mus_e_ho_isLoaded = false;
		mus_e_hoS9_isLoaded = false;
		mus_etaErr_isLoaded = false;
		mus_gfit_chi2_isLoaded = false;
		mus_gfit_ndof_isLoaded = false;
		mus_iso_isLoaded = false;
		mus_iso03_emEt_isLoaded = false;
		mus_iso03_hadEt_isLoaded = false;
		mus_iso03_hoEt_isLoaded = false;
		mus_iso03_sumPt_isLoaded = false;
		mus_iso05_emEt_isLoaded = false;
		mus_iso05_hadEt_isLoaded = false;
		mus_iso05_hoEt_isLoaded = false;
		mus_iso05_sumPt_isLoaded = false;
		mus_ndof_isLoaded = false;
		mus_outerEta_isLoaded = false;
		mus_outerPhi_isLoaded = false;
		mus_phiErr_isLoaded = false;
		mus_ptErr_isLoaded = false;
		mus_qoverp_isLoaded = false;
		mus_qoverpError_isLoaded = false;
		mus_vertexphi_isLoaded = false;
		mus_z0_isLoaded = false;
		mus_z0Err_isLoaded = false;
		mus_z0corr_isLoaded = false;
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_d0corr_isLoaded = false;
		trks_d0corrPhi_isLoaded = false;
		trks_etaErr_isLoaded = false;
		trks_layer1_charge_isLoaded = false;
		trks_ndof_isLoaded = false;
		trks_outerEt_isLoaded = false;
		trks_outerEta_isLoaded = false;
		trks_outerPhi_isLoaded = false;
		trks_phiErr_isLoaded = false;
		trks_ptErr_isLoaded = false;
		trks_tkIso_isLoaded = false;
		trks_vertexphi_isLoaded = false;
		trks_z0_isLoaded = false;
		trks_z0Err_isLoaded = false;
		trks_z0corr_isLoaded = false;
		trks_elsdr_isLoaded = false;
		trks_elsshFrac_isLoaded = false;
		trk_musdr_isLoaded = false;
		vtxs_chi2_isLoaded = false;
		vtxs_ndof_isLoaded = false;
		vtxs_xError_isLoaded = false;
		vtxs_yError_isLoaded = false;
		vtxs_zError_isLoaded = false;
		trks_residualX_isLoaded = false;
		trks_residualY_isLoaded = false;
		evt_HLT1_isLoaded = false;
		evt_HLT2_isLoaded = false;
		evt_HLT3_isLoaded = false;
		evt_HLT4_isLoaded = false;
		evt_HLT5_isLoaded = false;
		evt_HLT6_isLoaded = false;
		evt_HLT7_isLoaded = false;
		evt_HLT8_isLoaded = false;
		evt_L1_1_isLoaded = false;
		evt_L1_2_isLoaded = false;
		evt_L1_3_isLoaded = false;
		evt_L1_4_isLoaded = false;
		els_mc3_id_isLoaded = false;
		els_mc3idx_isLoaded = false;
		els_mc3_motherid_isLoaded = false;
		els_mc_id_isLoaded = false;
		els_mcidx_isLoaded = false;
		els_mc_motherid_isLoaded = false;
		jets_mc_id_isLoaded = false;
		mus_mc3_id_isLoaded = false;
		mus_mc3idx_isLoaded = false;
		mus_mc3_motherid_isLoaded = false;
		mus_mc_id_isLoaded = false;
		mus_mcidx_isLoaded = false;
		mus_mc_motherid_isLoaded = false;
		trk_mc3_id_isLoaded = false;
		trk_mc3idx_isLoaded = false;
		trk_mc3_motherid_isLoaded = false;
		trk_mc_id_isLoaded = false;
		trk_mcidx_isLoaded = false;
		trk_mc_motherid_isLoaded = false;
		els_conv_tkidx_isLoaded = false;
		trks_conv_tkidx_isLoaded = false;
		els_closestMuon_isLoaded = false;
		els_trkidx_isLoaded = false;
		els_category_isLoaded = false;
		els_categoryold_isLoaded = false;
		els_charge_isLoaded = false;
		els_class_isLoaded = false;
		els_layer1_det_isLoaded = false;
		els_layer1_layer_isLoaded = false;
		els_layer1_sizerphi_isLoaded = false;
		els_layer1_sizerz_isLoaded = false;
		els_looseId_isLoaded = false;
		els_lostHits_isLoaded = false;
		els_lost_pixelhits_isLoaded = false;
		els_nSeed_isLoaded = false;
		els_n_inner_layers_isLoaded = false;
		els_n_outer_layers_isLoaded = false;
		els_pass3looseId_isLoaded = false;
		els_pass3simpleId_isLoaded = false;
		els_pass3tightId_isLoaded = false;
		els_robustId_isLoaded = false;
		els_simpleIdPlus_isLoaded = false;
		els_tightId_isLoaded = false;
		els_tightId22XMaxMatteo_isLoaded = false;
		els_tightId22XMinMatteo_isLoaded = false;
		els_validHits_isLoaded = false;
		els_valid_pixelhits_isLoaded = false;
		genps_id_isLoaded = false;
		genps_id_mother_isLoaded = false;
		genps_lepdaughter_id_isLoaded = false;
		genps_lepdaughter_idx_isLoaded = false;
		genps_status_isLoaded = false;
		mus_charge_isLoaded = false;
		mus_gfit_validHits_isLoaded = false;
		mus_goodmask_isLoaded = false;
		mus_iso03_ntrk_isLoaded = false;
		mus_iso05_ntrk_isLoaded = false;
		mus_lostHits_isLoaded = false;
		mus_nmatches_isLoaded = false;
		mus_pid_TM2DCompatibilityLoose_isLoaded = false;
		mus_pid_TM2DCompatibilityTight_isLoaded = false;
		mus_pid_TMLastStationLoose_isLoaded = false;
		mus_pid_TMLastStationTight_isLoaded = false;
		mus_trk_charge_isLoaded = false;
		mus_trkrefkey_isLoaded = false;
		mus_type_isLoaded = false;
		mus_validHits_isLoaded = false;
		trks_charge_isLoaded = false;
		trks_layer1_det_isLoaded = false;
		trks_layer1_layer_isLoaded = false;
		trks_layer1_sizerphi_isLoaded = false;
		trks_layer1_sizerz_isLoaded = false;
		trks_lostHits_isLoaded = false;
		trks_lost_pixelhits_isLoaded = false;
		trks_n_inner_layers_isLoaded = false;
		trks_n_outer_layers_isLoaded = false;
		trks_validHits_isLoaded = false;
		trks_valid_pixelhits_isLoaded = false;
		trks_elsidx_isLoaded = false;
		trk_musidx_isLoaded = false;
		vtxs_hasRefittedTracks_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		vtxs_tracksSize_isLoaded = false;
		trks_hit_type_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		evt_nvtxs_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (evt_CMS2tag_branch != 0) evt_CMS2tag();
	if (evt_dataset_branch != 0) evt_dataset();
	if (evt_HLT_trigNames_branch != 0) evt_HLT_trigNames();
	if (evt_L1_trigNames_branch != 0) evt_L1_trigNames();
	if (genps_pthat_branch != 0) genps_pthat();
	if (vtxs_position_branch != 0) vtxs_position();
	if (evt_bsp4_branch != 0) evt_bsp4();
	if (els_mc3_motherp4_branch != 0) els_mc3_motherp4();
	if (els_mc3_p4_branch != 0) els_mc3_p4();
	if (els_mc_motherp4_branch != 0) els_mc_motherp4();
	if (els_mc_p4_branch != 0) els_mc_p4();
	if (jets_mc_gp_p4_branch != 0) jets_mc_gp_p4();
	if (jets_mc_p4_branch != 0) jets_mc_p4();
	if (mus_mc3_motherp4_branch != 0) mus_mc3_motherp4();
	if (mus_mc3_p4_branch != 0) mus_mc3_p4();
	if (mus_mc_motherp4_branch != 0) mus_mc_motherp4();
	if (mus_mc_p4_branch != 0) mus_mc_p4();
	if (trk_mc3p4_branch != 0) trk_mc3p4();
	if (trk_mcp4_branch != 0) trk_mcp4();
	if (els_p4_branch != 0) els_p4();
	if (els_p4In_branch != 0) els_p4In();
	if (els_p4Out_branch != 0) els_p4Out();
	if (els_trk_p4_branch != 0) els_trk_p4();
	if (els_vertex_p4_branch != 0) els_vertex_p4();
	if (genps_lepdaughter_p4_branch != 0) genps_lepdaughter_p4();
	if (genps_p4_branch != 0) genps_p4();
	if (genps_prod_vtx_branch != 0) genps_prod_vtx();
	if (mus_gfit_outerPos_p4_branch != 0) mus_gfit_outerPos_p4();
	if (mus_p4_branch != 0) mus_p4();
	if (mus_trk_p4_branch != 0) mus_trk_p4();
	if (mus_vertex_p4_branch != 0) mus_vertex_p4();
	if (trks_trk_p4_branch != 0) trks_trk_p4();
	if (trks_vertex_p4_branch != 0) trks_vertex_p4();
	if (els_inner_positionx_branch != 0) els_inner_positionx();
	if (els_inner_positiony_branch != 0) els_inner_positiony();
	if (els_inner_positionz_branch != 0) els_inner_positionz();
	if (trks_inner_positionx_branch != 0) trks_inner_positionx();
	if (trks_inner_positiony_branch != 0) trks_inner_positiony();
	if (trks_inner_positionz_branch != 0) trks_inner_positionz();
	if (evt_bs_dxdz_branch != 0) evt_bs_dxdz();
	if (evt_bs_dxdzErr_branch != 0) evt_bs_dxdzErr();
	if (evt_bs_dydz_branch != 0) evt_bs_dydz();
	if (evt_bs_dydzErr_branch != 0) evt_bs_dydzErr();
	if (evt_bs_sigmaZ_branch != 0) evt_bs_sigmaZ();
	if (evt_bs_sigmaZErr_branch != 0) evt_bs_sigmaZErr();
	if (evt_bs_widthx_branch != 0) evt_bs_widthx();
	if (evt_bs_widthxErr_branch != 0) evt_bs_widthxErr();
	if (evt_bs_widthy_branch != 0) evt_bs_widthy();
	if (evt_bs_widthyErr_branch != 0) evt_bs_widthyErr();
	if (evt_bs_xErr_branch != 0) evt_bs_xErr();
	if (evt_bs_yErr_branch != 0) evt_bs_yErr();
	if (evt_bs_zErr_branch != 0) evt_bs_zErr();
	if (evt_bField_branch != 0) evt_bField();
	if (evt_kfactor_branch != 0) evt_kfactor();
	if (evt_weight_branch != 0) evt_weight();
	if (evt_xsec_excl_branch != 0) evt_xsec_excl();
	if (evt_xsec_incl_branch != 0) evt_xsec_incl();
	if (gen_met_branch != 0) gen_met();
	if (gen_metPhi_branch != 0) gen_metPhi();
	if (els_mc3dr_branch != 0) els_mc3dr();
	if (els_mcdr_branch != 0) els_mcdr();
	if (jets_mcdr_branch != 0) jets_mcdr();
	if (jets_mc_emEnergy_branch != 0) jets_mc_emEnergy();
	if (jets_mc_gpdr_branch != 0) jets_mc_gpdr();
	if (jets_mc_hadEnergy_branch != 0) jets_mc_hadEnergy();
	if (jets_mc_invEnergy_branch != 0) jets_mc_invEnergy();
	if (jets_mc_otherEnergy_branch != 0) jets_mc_otherEnergy();
	if (mus_mc3dr_branch != 0) mus_mc3dr();
	if (mus_mcdr_branch != 0) mus_mcdr();
	if (trk_mc3dr_branch != 0) trk_mc3dr();
	if (trk_mcdr_branch != 0) trk_mcdr();
	if (els_conv_dcot_branch != 0) els_conv_dcot();
	if (els_conv_dist_branch != 0) els_conv_dist();
	if (trks_conv_dcot_branch != 0) trks_conv_dcot();
	if (trks_conv_dist_branch != 0) trks_conv_dist();
	if (els_musdr_branch != 0) els_musdr();
	if (els_trkdr_branch != 0) els_trkdr();
	if (els_trkshFrac_branch != 0) els_trkshFrac();
	if (els_chi2_branch != 0) els_chi2();
	if (els_d0_branch != 0) els_d0();
	if (els_d0Err_branch != 0) els_d0Err();
	if (els_d0corr_branch != 0) els_d0corr();
	if (els_dEtaIn_branch != 0) els_dEtaIn();
	if (els_dEtaOut_branch != 0) els_dEtaOut();
	if (els_dPhiIn_branch != 0) els_dPhiIn();
	if (els_dPhiInPhiOut_branch != 0) els_dPhiInPhiOut();
	if (els_dPhiOut_branch != 0) els_dPhiOut();
	if (els_e1x5_branch != 0) els_e1x5();
	if (els_e2x5Max_branch != 0) els_e2x5Max();
	if (els_e3x3_branch != 0) els_e3x3();
	if (els_e5x5_branch != 0) els_e5x5();
	if (els_eMax_branch != 0) els_eMax();
	if (els_eOverPIn_branch != 0) els_eOverPIn();
	if (els_eSC_branch != 0) els_eSC();
	if (els_eSCPresh_branch != 0) els_eSCPresh();
	if (els_eSCRaw_branch != 0) els_eSCRaw();
	if (els_eSeed_branch != 0) els_eSeed();
	if (els_eSeedOverPIn_branch != 0) els_eSeedOverPIn();
	if (els_eSeedOverPOut_branch != 0) els_eSeedOverPOut();
	if (els_ecalIso03_branch != 0) els_ecalIso03();
	if (els_ecalIso04_branch != 0) els_ecalIso04();
	if (els_egamma_looseId_branch != 0) els_egamma_looseId();
	if (els_egamma_robustHighEnergy_branch != 0) els_egamma_robustHighEnergy();
	if (els_egamma_robustLooseId_branch != 0) els_egamma_robustLooseId();
	if (els_egamma_robustTightId_branch != 0) els_egamma_robustTightId();
	if (els_egamma_tightId_branch != 0) els_egamma_tightId();
	if (els_etaErr_branch != 0) els_etaErr();
	if (els_etaSC_branch != 0) els_etaSC();
	if (els_fBrem_branch != 0) els_fBrem();
	if (els_hOverE_branch != 0) els_hOverE();
	if (els_hcalIso03_branch != 0) els_hcalIso03();
	if (els_hcalIso04_branch != 0) els_hcalIso04();
	if (els_layer1_charge_branch != 0) els_layer1_charge();
	if (els_ndof_branch != 0) els_ndof();
	if (els_outerEta_branch != 0) els_outerEta();
	if (els_outerPhi_branch != 0) els_outerPhi();
	if (els_phiErr_branch != 0) els_phiErr();
	if (els_phiSC_branch != 0) els_phiSC();
	if (els_pin_branch != 0) els_pin();
	if (els_pout_branch != 0) els_pout();
	if (els_ptErr_branch != 0) els_ptErr();
	if (els_sigmaEtaEta_branch != 0) els_sigmaEtaEta();
	if (els_sigmaIEtaIEta_branch != 0) els_sigmaIEtaIEta();
	if (els_sigmaIPhiIPhi_branch != 0) els_sigmaIPhiIPhi();
	if (els_sigmaPhiPhi_branch != 0) els_sigmaPhiPhi();
	if (els_tkIso03_branch != 0) els_tkIso03();
	if (els_tkIso04_branch != 0) els_tkIso04();
	if (els_vertexphi_branch != 0) els_vertexphi();
	if (els_z0_branch != 0) els_z0();
	if (els_z0Err_branch != 0) els_z0Err();
	if (els_z0corr_branch != 0) els_z0corr();
	if (mus_caloCompatibility_branch != 0) mus_caloCompatibility();
	if (mus_chi2_branch != 0) mus_chi2();
	if (mus_d0_branch != 0) mus_d0();
	if (mus_d0Err_branch != 0) mus_d0Err();
	if (mus_d0corr_branch != 0) mus_d0corr();
	if (mus_e_em_branch != 0) mus_e_em();
	if (mus_e_emS9_branch != 0) mus_e_emS9();
	if (mus_e_had_branch != 0) mus_e_had();
	if (mus_e_hadS9_branch != 0) mus_e_hadS9();
	if (mus_e_ho_branch != 0) mus_e_ho();
	if (mus_e_hoS9_branch != 0) mus_e_hoS9();
	if (mus_etaErr_branch != 0) mus_etaErr();
	if (mus_gfit_chi2_branch != 0) mus_gfit_chi2();
	if (mus_gfit_ndof_branch != 0) mus_gfit_ndof();
	if (mus_iso_branch != 0) mus_iso();
	if (mus_iso03_emEt_branch != 0) mus_iso03_emEt();
	if (mus_iso03_hadEt_branch != 0) mus_iso03_hadEt();
	if (mus_iso03_hoEt_branch != 0) mus_iso03_hoEt();
	if (mus_iso03_sumPt_branch != 0) mus_iso03_sumPt();
	if (mus_iso05_emEt_branch != 0) mus_iso05_emEt();
	if (mus_iso05_hadEt_branch != 0) mus_iso05_hadEt();
	if (mus_iso05_hoEt_branch != 0) mus_iso05_hoEt();
	if (mus_iso05_sumPt_branch != 0) mus_iso05_sumPt();
	if (mus_ndof_branch != 0) mus_ndof();
	if (mus_outerEta_branch != 0) mus_outerEta();
	if (mus_outerPhi_branch != 0) mus_outerPhi();
	if (mus_phiErr_branch != 0) mus_phiErr();
	if (mus_ptErr_branch != 0) mus_ptErr();
	if (mus_qoverp_branch != 0) mus_qoverp();
	if (mus_qoverpError_branch != 0) mus_qoverpError();
	if (mus_vertexphi_branch != 0) mus_vertexphi();
	if (mus_z0_branch != 0) mus_z0();
	if (mus_z0Err_branch != 0) mus_z0Err();
	if (mus_z0corr_branch != 0) mus_z0corr();
	if (trks_chi2_branch != 0) trks_chi2();
	if (trks_d0_branch != 0) trks_d0();
	if (trks_d0Err_branch != 0) trks_d0Err();
	if (trks_d0corr_branch != 0) trks_d0corr();
	if (trks_d0corrPhi_branch != 0) trks_d0corrPhi();
	if (trks_etaErr_branch != 0) trks_etaErr();
	if (trks_layer1_charge_branch != 0) trks_layer1_charge();
	if (trks_ndof_branch != 0) trks_ndof();
	if (trks_outerEt_branch != 0) trks_outerEt();
	if (trks_outerEta_branch != 0) trks_outerEta();
	if (trks_outerPhi_branch != 0) trks_outerPhi();
	if (trks_phiErr_branch != 0) trks_phiErr();
	if (trks_ptErr_branch != 0) trks_ptErr();
	if (trks_tkIso_branch != 0) trks_tkIso();
	if (trks_vertexphi_branch != 0) trks_vertexphi();
	if (trks_z0_branch != 0) trks_z0();
	if (trks_z0Err_branch != 0) trks_z0Err();
	if (trks_z0corr_branch != 0) trks_z0corr();
	if (trks_elsdr_branch != 0) trks_elsdr();
	if (trks_elsshFrac_branch != 0) trks_elsshFrac();
	if (trk_musdr_branch != 0) trk_musdr();
	if (vtxs_chi2_branch != 0) vtxs_chi2();
	if (vtxs_ndof_branch != 0) vtxs_ndof();
	if (vtxs_xError_branch != 0) vtxs_xError();
	if (vtxs_yError_branch != 0) vtxs_yError();
	if (vtxs_zError_branch != 0) vtxs_zError();
	if (trks_residualX_branch != 0) trks_residualX();
	if (trks_residualY_branch != 0) trks_residualY();
	if (evt_HLT1_branch != 0) evt_HLT1();
	if (evt_HLT2_branch != 0) evt_HLT2();
	if (evt_HLT3_branch != 0) evt_HLT3();
	if (evt_HLT4_branch != 0) evt_HLT4();
	if (evt_HLT5_branch != 0) evt_HLT5();
	if (evt_HLT6_branch != 0) evt_HLT6();
	if (evt_HLT7_branch != 0) evt_HLT7();
	if (evt_HLT8_branch != 0) evt_HLT8();
	if (evt_L1_1_branch != 0) evt_L1_1();
	if (evt_L1_2_branch != 0) evt_L1_2();
	if (evt_L1_3_branch != 0) evt_L1_3();
	if (evt_L1_4_branch != 0) evt_L1_4();
	if (els_mc3_id_branch != 0) els_mc3_id();
	if (els_mc3idx_branch != 0) els_mc3idx();
	if (els_mc3_motherid_branch != 0) els_mc3_motherid();
	if (els_mc_id_branch != 0) els_mc_id();
	if (els_mcidx_branch != 0) els_mcidx();
	if (els_mc_motherid_branch != 0) els_mc_motherid();
	if (jets_mc_id_branch != 0) jets_mc_id();
	if (mus_mc3_id_branch != 0) mus_mc3_id();
	if (mus_mc3idx_branch != 0) mus_mc3idx();
	if (mus_mc3_motherid_branch != 0) mus_mc3_motherid();
	if (mus_mc_id_branch != 0) mus_mc_id();
	if (mus_mcidx_branch != 0) mus_mcidx();
	if (mus_mc_motherid_branch != 0) mus_mc_motherid();
	if (trk_mc3_id_branch != 0) trk_mc3_id();
	if (trk_mc3idx_branch != 0) trk_mc3idx();
	if (trk_mc3_motherid_branch != 0) trk_mc3_motherid();
	if (trk_mc_id_branch != 0) trk_mc_id();
	if (trk_mcidx_branch != 0) trk_mcidx();
	if (trk_mc_motherid_branch != 0) trk_mc_motherid();
	if (els_conv_tkidx_branch != 0) els_conv_tkidx();
	if (trks_conv_tkidx_branch != 0) trks_conv_tkidx();
	if (els_closestMuon_branch != 0) els_closestMuon();
	if (els_trkidx_branch != 0) els_trkidx();
	if (els_category_branch != 0) els_category();
	if (els_categoryold_branch != 0) els_categoryold();
	if (els_charge_branch != 0) els_charge();
	if (els_class_branch != 0) els_class();
	if (els_layer1_det_branch != 0) els_layer1_det();
	if (els_layer1_layer_branch != 0) els_layer1_layer();
	if (els_layer1_sizerphi_branch != 0) els_layer1_sizerphi();
	if (els_layer1_sizerz_branch != 0) els_layer1_sizerz();
	if (els_looseId_branch != 0) els_looseId();
	if (els_lostHits_branch != 0) els_lostHits();
	if (els_lost_pixelhits_branch != 0) els_lost_pixelhits();
	if (els_nSeed_branch != 0) els_nSeed();
	if (els_n_inner_layers_branch != 0) els_n_inner_layers();
	if (els_n_outer_layers_branch != 0) els_n_outer_layers();
	if (els_pass3looseId_branch != 0) els_pass3looseId();
	if (els_pass3simpleId_branch != 0) els_pass3simpleId();
	if (els_pass3tightId_branch != 0) els_pass3tightId();
	if (els_robustId_branch != 0) els_robustId();
	if (els_simpleIdPlus_branch != 0) els_simpleIdPlus();
	if (els_tightId_branch != 0) els_tightId();
	if (els_tightId22XMaxMatteo_branch != 0) els_tightId22XMaxMatteo();
	if (els_tightId22XMinMatteo_branch != 0) els_tightId22XMinMatteo();
	if (els_validHits_branch != 0) els_validHits();
	if (els_valid_pixelhits_branch != 0) els_valid_pixelhits();
	if (genps_id_branch != 0) genps_id();
	if (genps_id_mother_branch != 0) genps_id_mother();
	if (genps_lepdaughter_id_branch != 0) genps_lepdaughter_id();
	if (genps_lepdaughter_idx_branch != 0) genps_lepdaughter_idx();
	if (genps_status_branch != 0) genps_status();
	if (mus_charge_branch != 0) mus_charge();
	if (mus_gfit_validHits_branch != 0) mus_gfit_validHits();
	if (mus_goodmask_branch != 0) mus_goodmask();
	if (mus_iso03_ntrk_branch != 0) mus_iso03_ntrk();
	if (mus_iso05_ntrk_branch != 0) mus_iso05_ntrk();
	if (mus_lostHits_branch != 0) mus_lostHits();
	if (mus_nmatches_branch != 0) mus_nmatches();
	if (mus_pid_TM2DCompatibilityLoose_branch != 0) mus_pid_TM2DCompatibilityLoose();
	if (mus_pid_TM2DCompatibilityTight_branch != 0) mus_pid_TM2DCompatibilityTight();
	if (mus_pid_TMLastStationLoose_branch != 0) mus_pid_TMLastStationLoose();
	if (mus_pid_TMLastStationTight_branch != 0) mus_pid_TMLastStationTight();
	if (mus_trk_charge_branch != 0) mus_trk_charge();
	if (mus_trkrefkey_branch != 0) mus_trkrefkey();
	if (mus_type_branch != 0) mus_type();
	if (mus_validHits_branch != 0) mus_validHits();
	if (trks_charge_branch != 0) trks_charge();
	if (trks_layer1_det_branch != 0) trks_layer1_det();
	if (trks_layer1_layer_branch != 0) trks_layer1_layer();
	if (trks_layer1_sizerphi_branch != 0) trks_layer1_sizerphi();
	if (trks_layer1_sizerz_branch != 0) trks_layer1_sizerz();
	if (trks_lostHits_branch != 0) trks_lostHits();
	if (trks_lost_pixelhits_branch != 0) trks_lost_pixelhits();
	if (trks_n_inner_layers_branch != 0) trks_n_inner_layers();
	if (trks_n_outer_layers_branch != 0) trks_n_outer_layers();
	if (trks_validHits_branch != 0) trks_validHits();
	if (trks_valid_pixelhits_branch != 0) trks_valid_pixelhits();
	if (trks_elsidx_branch != 0) trks_elsidx();
	if (trk_musidx_branch != 0) trk_musidx();
	if (vtxs_hasRefittedTracks_branch != 0) vtxs_hasRefittedTracks();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (vtxs_tracksSize_branch != 0) vtxs_tracksSize();
	if (trks_hit_type_branch != 0) trks_hit_type();
	if (evt_nels_branch != 0) evt_nels();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
}

	TString &evt_CMS2tag()
	{
		if (not evt_CMS2tag_isLoaded) {
			if (evt_CMS2tag_branch != 0) {
				evt_CMS2tag_branch->GetEntry(index);
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
			} else { 
				printf("branch evt_dataset_branch does not exist!\n");
				exit(1);
			}
			evt_dataset_isLoaded = true;
		}
		return evt_dataset_;
	}
	vector<TString> &evt_HLT_trigNames()
	{
		if (not evt_HLT_trigNames_isLoaded) {
			if (evt_HLT_trigNames_branch != 0) {
				evt_HLT_trigNames_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT_trigNames_branch does not exist!\n");
				exit(1);
			}
			evt_HLT_trigNames_isLoaded = true;
		}
		return evt_HLT_trigNames_;
	}
	vector<TString> &evt_L1_trigNames()
	{
		if (not evt_L1_trigNames_isLoaded) {
			if (evt_L1_trigNames_branch != 0) {
				evt_L1_trigNames_branch->GetEntry(index);
			} else { 
				printf("branch evt_L1_trigNames_branch does not exist!\n");
				exit(1);
			}
			evt_L1_trigNames_isLoaded = true;
		}
		return evt_L1_trigNames_;
	}
	double &genps_pthat()
	{
		if (not genps_pthat_isLoaded) {
			if (genps_pthat_branch != 0) {
				genps_pthat_branch->GetEntry(index);
			} else { 
				printf("branch genps_pthat_branch does not exist!\n");
				exit(1);
			}
			genps_pthat_isLoaded = true;
		}
		return genps_pthat_;
	}
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &vtxs_position()
	{
		if (not vtxs_position_isLoaded) {
			if (vtxs_position_branch != 0) {
				vtxs_position_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_position_branch does not exist!\n");
				exit(1);
			}
			vtxs_position_isLoaded = true;
		}
		return vtxs_position_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &evt_bsp4()
	{
		if (not evt_bsp4_isLoaded) {
			if (evt_bsp4_branch != 0) {
				evt_bsp4_branch->GetEntry(index);
			} else { 
				printf("branch evt_bsp4_branch does not exist!\n");
				exit(1);
			}
			evt_bsp4_isLoaded = true;
		}
		return evt_bsp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_motherp4()
	{
		if (not els_mc3_motherp4_isLoaded) {
			if (els_mc3_motherp4_branch != 0) {
				els_mc3_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3_motherp4_branch does not exist!\n");
				exit(1);
			}
			els_mc3_motherp4_isLoaded = true;
		}
		return els_mc3_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_p4()
	{
		if (not els_mc3_p4_isLoaded) {
			if (els_mc3_p4_branch != 0) {
				els_mc3_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3_p4_branch does not exist!\n");
				exit(1);
			}
			els_mc3_p4_isLoaded = true;
		}
		return els_mc3_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_motherp4()
	{
		if (not els_mc_motherp4_isLoaded) {
			if (els_mc_motherp4_branch != 0) {
				els_mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			els_mc_motherp4_isLoaded = true;
		}
		return els_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_p4()
	{
		if (not els_mc_p4_isLoaded) {
			if (els_mc_p4_branch != 0) {
				els_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_p4_branch does not exist!\n");
				exit(1);
			}
			els_mc_p4_isLoaded = true;
		}
		return els_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_gp_p4()
	{
		if (not jets_mc_gp_p4_isLoaded) {
			if (jets_mc_gp_p4_branch != 0) {
				jets_mc_gp_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_gp_p4_branch does not exist!\n");
				exit(1);
			}
			jets_mc_gp_p4_isLoaded = true;
		}
		return jets_mc_gp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_p4()
	{
		if (not jets_mc_p4_isLoaded) {
			if (jets_mc_p4_branch != 0) {
				jets_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_p4_branch does not exist!\n");
				exit(1);
			}
			jets_mc_p4_isLoaded = true;
		}
		return jets_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc3_motherp4()
	{
		if (not mus_mc3_motherp4_isLoaded) {
			if (mus_mc3_motherp4_branch != 0) {
				mus_mc3_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3_motherp4_branch does not exist!\n");
				exit(1);
			}
			mus_mc3_motherp4_isLoaded = true;
		}
		return mus_mc3_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc3_p4()
	{
		if (not mus_mc3_p4_isLoaded) {
			if (mus_mc3_p4_branch != 0) {
				mus_mc3_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3_p4_branch does not exist!\n");
				exit(1);
			}
			mus_mc3_p4_isLoaded = true;
		}
		return mus_mc3_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc_motherp4()
	{
		if (not mus_mc_motherp4_isLoaded) {
			if (mus_mc_motherp4_branch != 0) {
				mus_mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			mus_mc_motherp4_isLoaded = true;
		}
		return mus_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc_p4()
	{
		if (not mus_mc_p4_isLoaded) {
			if (mus_mc_p4_branch != 0) {
				mus_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc_p4_branch does not exist!\n");
				exit(1);
			}
			mus_mc_p4_isLoaded = true;
		}
		return mus_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trk_mc3p4()
	{
		if (not trk_mc3p4_isLoaded) {
			if (trk_mc3p4_branch != 0) {
				trk_mc3p4_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3p4_branch does not exist!\n");
				exit(1);
			}
			trk_mc3p4_isLoaded = true;
		}
		return trk_mc3p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trk_mcp4()
	{
		if (not trk_mcp4_isLoaded) {
			if (trk_mcp4_branch != 0) {
				trk_mcp4_branch->GetEntry(index);
			} else { 
				printf("branch trk_mcp4_branch does not exist!\n");
				exit(1);
			}
			trk_mcp4_isLoaded = true;
		}
		return trk_mcp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) {
				els_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_p4_branch does not exist!\n");
				exit(1);
			}
			els_p4_isLoaded = true;
		}
		return els_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4In()
	{
		if (not els_p4In_isLoaded) {
			if (els_p4In_branch != 0) {
				els_p4In_branch->GetEntry(index);
			} else { 
				printf("branch els_p4In_branch does not exist!\n");
				exit(1);
			}
			els_p4In_isLoaded = true;
		}
		return els_p4In_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4Out()
	{
		if (not els_p4Out_isLoaded) {
			if (els_p4Out_branch != 0) {
				els_p4Out_branch->GetEntry(index);
			} else { 
				printf("branch els_p4Out_branch does not exist!\n");
				exit(1);
			}
			els_p4Out_isLoaded = true;
		}
		return els_p4Out_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_trk_p4()
	{
		if (not els_trk_p4_isLoaded) {
			if (els_trk_p4_branch != 0) {
				els_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_trk_p4_branch does not exist!\n");
				exit(1);
			}
			els_trk_p4_isLoaded = true;
		}
		return els_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_vertex_p4()
	{
		if (not els_vertex_p4_isLoaded) {
			if (els_vertex_p4_branch != 0) {
				els_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			els_vertex_p4_isLoaded = true;
		}
		return els_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_lepdaughter_p4()
	{
		if (not genps_lepdaughter_p4_isLoaded) {
			if (genps_lepdaughter_p4_branch != 0) {
				genps_lepdaughter_p4_branch->GetEntry(index);
			} else { 
				printf("branch genps_lepdaughter_p4_branch does not exist!\n");
				exit(1);
			}
			genps_lepdaughter_p4_isLoaded = true;
		}
		return genps_lepdaughter_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_p4()
	{
		if (not genps_p4_isLoaded) {
			if (genps_p4_branch != 0) {
				genps_p4_branch->GetEntry(index);
			} else { 
				printf("branch genps_p4_branch does not exist!\n");
				exit(1);
			}
			genps_p4_isLoaded = true;
		}
		return genps_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_prod_vtx()
	{
		if (not genps_prod_vtx_isLoaded) {
			if (genps_prod_vtx_branch != 0) {
				genps_prod_vtx_branch->GetEntry(index);
			} else { 
				printf("branch genps_prod_vtx_branch does not exist!\n");
				exit(1);
			}
			genps_prod_vtx_isLoaded = true;
		}
		return genps_prod_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_gfit_outerPos_p4()
	{
		if (not mus_gfit_outerPos_p4_isLoaded) {
			if (mus_gfit_outerPos_p4_branch != 0) {
				mus_gfit_outerPos_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_outerPos_p4_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_outerPos_p4_isLoaded = true;
		}
		return mus_gfit_outerPos_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4()
	{
		if (not mus_p4_isLoaded) {
			if (mus_p4_branch != 0) {
				mus_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_p4_branch does not exist!\n");
				exit(1);
			}
			mus_p4_isLoaded = true;
		}
		return mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4()
	{
		if (not mus_trk_p4_isLoaded) {
			if (mus_trk_p4_branch != 0) {
				mus_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_trk_p4_branch does not exist!\n");
				exit(1);
			}
			mus_trk_p4_isLoaded = true;
		}
		return mus_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_vertex_p4()
	{
		if (not mus_vertex_p4_isLoaded) {
			if (mus_vertex_p4_branch != 0) {
				mus_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			mus_vertex_p4_isLoaded = true;
		}
		return mus_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4()
	{
		if (not trks_trk_p4_isLoaded) {
			if (trks_trk_p4_branch != 0) {
				trks_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch trks_trk_p4_branch does not exist!\n");
				exit(1);
			}
			trks_trk_p4_isLoaded = true;
		}
		return trks_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_vertex_p4()
	{
		if (not trks_vertex_p4_isLoaded) {
			if (trks_vertex_p4_branch != 0) {
				trks_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch trks_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			trks_vertex_p4_isLoaded = true;
		}
		return trks_vertex_p4_;
	}
	vector<double> &els_inner_positionx()
	{
		if (not els_inner_positionx_isLoaded) {
			if (els_inner_positionx_branch != 0) {
				els_inner_positionx_branch->GetEntry(index);
			} else { 
				printf("branch els_inner_positionx_branch does not exist!\n");
				exit(1);
			}
			els_inner_positionx_isLoaded = true;
		}
		return els_inner_positionx_;
	}
	vector<double> &els_inner_positiony()
	{
		if (not els_inner_positiony_isLoaded) {
			if (els_inner_positiony_branch != 0) {
				els_inner_positiony_branch->GetEntry(index);
			} else { 
				printf("branch els_inner_positiony_branch does not exist!\n");
				exit(1);
			}
			els_inner_positiony_isLoaded = true;
		}
		return els_inner_positiony_;
	}
	vector<double> &els_inner_positionz()
	{
		if (not els_inner_positionz_isLoaded) {
			if (els_inner_positionz_branch != 0) {
				els_inner_positionz_branch->GetEntry(index);
			} else { 
				printf("branch els_inner_positionz_branch does not exist!\n");
				exit(1);
			}
			els_inner_positionz_isLoaded = true;
		}
		return els_inner_positionz_;
	}
	vector<double> &trks_inner_positionx()
	{
		if (not trks_inner_positionx_isLoaded) {
			if (trks_inner_positionx_branch != 0) {
				trks_inner_positionx_branch->GetEntry(index);
			} else { 
				printf("branch trks_inner_positionx_branch does not exist!\n");
				exit(1);
			}
			trks_inner_positionx_isLoaded = true;
		}
		return trks_inner_positionx_;
	}
	vector<double> &trks_inner_positiony()
	{
		if (not trks_inner_positiony_isLoaded) {
			if (trks_inner_positiony_branch != 0) {
				trks_inner_positiony_branch->GetEntry(index);
			} else { 
				printf("branch trks_inner_positiony_branch does not exist!\n");
				exit(1);
			}
			trks_inner_positiony_isLoaded = true;
		}
		return trks_inner_positiony_;
	}
	vector<double> &trks_inner_positionz()
	{
		if (not trks_inner_positionz_isLoaded) {
			if (trks_inner_positionz_branch != 0) {
				trks_inner_positionz_branch->GetEntry(index);
			} else { 
				printf("branch trks_inner_positionz_branch does not exist!\n");
				exit(1);
			}
			trks_inner_positionz_isLoaded = true;
		}
		return trks_inner_positionz_;
	}
	float &evt_bs_dxdz()
	{
		if (not evt_bs_dxdz_isLoaded) {
			if (evt_bs_dxdz_branch != 0) {
				evt_bs_dxdz_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dxdz_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dxdz_isLoaded = true;
		}
		return evt_bs_dxdz_;
	}
	float &evt_bs_dxdzErr()
	{
		if (not evt_bs_dxdzErr_isLoaded) {
			if (evt_bs_dxdzErr_branch != 0) {
				evt_bs_dxdzErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dxdzErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dxdzErr_isLoaded = true;
		}
		return evt_bs_dxdzErr_;
	}
	float &evt_bs_dydz()
	{
		if (not evt_bs_dydz_isLoaded) {
			if (evt_bs_dydz_branch != 0) {
				evt_bs_dydz_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dydz_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dydz_isLoaded = true;
		}
		return evt_bs_dydz_;
	}
	float &evt_bs_dydzErr()
	{
		if (not evt_bs_dydzErr_isLoaded) {
			if (evt_bs_dydzErr_branch != 0) {
				evt_bs_dydzErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_dydzErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_dydzErr_isLoaded = true;
		}
		return evt_bs_dydzErr_;
	}
	float &evt_bs_sigmaZ()
	{
		if (not evt_bs_sigmaZ_isLoaded) {
			if (evt_bs_sigmaZ_branch != 0) {
				evt_bs_sigmaZ_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_sigmaZ_branch does not exist!\n");
				exit(1);
			}
			evt_bs_sigmaZ_isLoaded = true;
		}
		return evt_bs_sigmaZ_;
	}
	float &evt_bs_sigmaZErr()
	{
		if (not evt_bs_sigmaZErr_isLoaded) {
			if (evt_bs_sigmaZErr_branch != 0) {
				evt_bs_sigmaZErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_sigmaZErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_sigmaZErr_isLoaded = true;
		}
		return evt_bs_sigmaZErr_;
	}
	float &evt_bs_widthx()
	{
		if (not evt_bs_widthx_isLoaded) {
			if (evt_bs_widthx_branch != 0) {
				evt_bs_widthx_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_widthx_branch does not exist!\n");
				exit(1);
			}
			evt_bs_widthx_isLoaded = true;
		}
		return evt_bs_widthx_;
	}
	float &evt_bs_widthxErr()
	{
		if (not evt_bs_widthxErr_isLoaded) {
			if (evt_bs_widthxErr_branch != 0) {
				evt_bs_widthxErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_widthxErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_widthxErr_isLoaded = true;
		}
		return evt_bs_widthxErr_;
	}
	float &evt_bs_widthy()
	{
		if (not evt_bs_widthy_isLoaded) {
			if (evt_bs_widthy_branch != 0) {
				evt_bs_widthy_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_widthy_branch does not exist!\n");
				exit(1);
			}
			evt_bs_widthy_isLoaded = true;
		}
		return evt_bs_widthy_;
	}
	float &evt_bs_widthyErr()
	{
		if (not evt_bs_widthyErr_isLoaded) {
			if (evt_bs_widthyErr_branch != 0) {
				evt_bs_widthyErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_widthyErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_widthyErr_isLoaded = true;
		}
		return evt_bs_widthyErr_;
	}
	float &evt_bs_xErr()
	{
		if (not evt_bs_xErr_isLoaded) {
			if (evt_bs_xErr_branch != 0) {
				evt_bs_xErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_xErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_xErr_isLoaded = true;
		}
		return evt_bs_xErr_;
	}
	float &evt_bs_yErr()
	{
		if (not evt_bs_yErr_isLoaded) {
			if (evt_bs_yErr_branch != 0) {
				evt_bs_yErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_yErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_yErr_isLoaded = true;
		}
		return evt_bs_yErr_;
	}
	float &evt_bs_zErr()
	{
		if (not evt_bs_zErr_isLoaded) {
			if (evt_bs_zErr_branch != 0) {
				evt_bs_zErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_zErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_zErr_isLoaded = true;
		}
		return evt_bs_zErr_;
	}
	float &evt_bField()
	{
		if (not evt_bField_isLoaded) {
			if (evt_bField_branch != 0) {
				evt_bField_branch->GetEntry(index);
			} else { 
				printf("branch evt_bField_branch does not exist!\n");
				exit(1);
			}
			evt_bField_isLoaded = true;
		}
		return evt_bField_;
	}
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) {
				evt_kfactor_branch->GetEntry(index);
			} else { 
				printf("branch evt_kfactor_branch does not exist!\n");
				exit(1);
			}
			evt_kfactor_isLoaded = true;
		}
		return evt_kfactor_;
	}
	float &evt_weight()
	{
		if (not evt_weight_isLoaded) {
			if (evt_weight_branch != 0) {
				evt_weight_branch->GetEntry(index);
			} else { 
				printf("branch evt_weight_branch does not exist!\n");
				exit(1);
			}
			evt_weight_isLoaded = true;
		}
		return evt_weight_;
	}
	float &evt_xsec_excl()
	{
		if (not evt_xsec_excl_isLoaded) {
			if (evt_xsec_excl_branch != 0) {
				evt_xsec_excl_branch->GetEntry(index);
			} else { 
				printf("branch evt_xsec_excl_branch does not exist!\n");
				exit(1);
			}
			evt_xsec_excl_isLoaded = true;
		}
		return evt_xsec_excl_;
	}
	float &evt_xsec_incl()
	{
		if (not evt_xsec_incl_isLoaded) {
			if (evt_xsec_incl_branch != 0) {
				evt_xsec_incl_branch->GetEntry(index);
			} else { 
				printf("branch evt_xsec_incl_branch does not exist!\n");
				exit(1);
			}
			evt_xsec_incl_isLoaded = true;
		}
		return evt_xsec_incl_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) {
				gen_met_branch->GetEntry(index);
			} else { 
				printf("branch gen_met_branch does not exist!\n");
				exit(1);
			}
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &gen_metPhi()
	{
		if (not gen_metPhi_isLoaded) {
			if (gen_metPhi_branch != 0) {
				gen_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch gen_metPhi_branch does not exist!\n");
				exit(1);
			}
			gen_metPhi_isLoaded = true;
		}
		return gen_metPhi_;
	}
	vector<float> &els_mc3dr()
	{
		if (not els_mc3dr_isLoaded) {
			if (els_mc3dr_branch != 0) {
				els_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3dr_branch does not exist!\n");
				exit(1);
			}
			els_mc3dr_isLoaded = true;
		}
		return els_mc3dr_;
	}
	vector<float> &els_mcdr()
	{
		if (not els_mcdr_isLoaded) {
			if (els_mcdr_branch != 0) {
				els_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch els_mcdr_branch does not exist!\n");
				exit(1);
			}
			els_mcdr_isLoaded = true;
		}
		return els_mcdr_;
	}
	vector<float> &jets_mcdr()
	{
		if (not jets_mcdr_isLoaded) {
			if (jets_mcdr_branch != 0) {
				jets_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch jets_mcdr_branch does not exist!\n");
				exit(1);
			}
			jets_mcdr_isLoaded = true;
		}
		return jets_mcdr_;
	}
	vector<float> &jets_mc_emEnergy()
	{
		if (not jets_mc_emEnergy_isLoaded) {
			if (jets_mc_emEnergy_branch != 0) {
				jets_mc_emEnergy_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_emEnergy_branch does not exist!\n");
				exit(1);
			}
			jets_mc_emEnergy_isLoaded = true;
		}
		return jets_mc_emEnergy_;
	}
	vector<float> &jets_mc_gpdr()
	{
		if (not jets_mc_gpdr_isLoaded) {
			if (jets_mc_gpdr_branch != 0) {
				jets_mc_gpdr_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_gpdr_branch does not exist!\n");
				exit(1);
			}
			jets_mc_gpdr_isLoaded = true;
		}
		return jets_mc_gpdr_;
	}
	vector<float> &jets_mc_hadEnergy()
	{
		if (not jets_mc_hadEnergy_isLoaded) {
			if (jets_mc_hadEnergy_branch != 0) {
				jets_mc_hadEnergy_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_hadEnergy_branch does not exist!\n");
				exit(1);
			}
			jets_mc_hadEnergy_isLoaded = true;
		}
		return jets_mc_hadEnergy_;
	}
	vector<float> &jets_mc_invEnergy()
	{
		if (not jets_mc_invEnergy_isLoaded) {
			if (jets_mc_invEnergy_branch != 0) {
				jets_mc_invEnergy_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_invEnergy_branch does not exist!\n");
				exit(1);
			}
			jets_mc_invEnergy_isLoaded = true;
		}
		return jets_mc_invEnergy_;
	}
	vector<float> &jets_mc_otherEnergy()
	{
		if (not jets_mc_otherEnergy_isLoaded) {
			if (jets_mc_otherEnergy_branch != 0) {
				jets_mc_otherEnergy_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_otherEnergy_branch does not exist!\n");
				exit(1);
			}
			jets_mc_otherEnergy_isLoaded = true;
		}
		return jets_mc_otherEnergy_;
	}
	vector<float> &mus_mc3dr()
	{
		if (not mus_mc3dr_isLoaded) {
			if (mus_mc3dr_branch != 0) {
				mus_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3dr_branch does not exist!\n");
				exit(1);
			}
			mus_mc3dr_isLoaded = true;
		}
		return mus_mc3dr_;
	}
	vector<float> &mus_mcdr()
	{
		if (not mus_mcdr_isLoaded) {
			if (mus_mcdr_branch != 0) {
				mus_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch mus_mcdr_branch does not exist!\n");
				exit(1);
			}
			mus_mcdr_isLoaded = true;
		}
		return mus_mcdr_;
	}
	vector<float> &trk_mc3dr()
	{
		if (not trk_mc3dr_isLoaded) {
			if (trk_mc3dr_branch != 0) {
				trk_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3dr_branch does not exist!\n");
				exit(1);
			}
			trk_mc3dr_isLoaded = true;
		}
		return trk_mc3dr_;
	}
	vector<float> &trk_mcdr()
	{
		if (not trk_mcdr_isLoaded) {
			if (trk_mcdr_branch != 0) {
				trk_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch trk_mcdr_branch does not exist!\n");
				exit(1);
			}
			trk_mcdr_isLoaded = true;
		}
		return trk_mcdr_;
	}
	vector<float> &els_conv_dcot()
	{
		if (not els_conv_dcot_isLoaded) {
			if (els_conv_dcot_branch != 0) {
				els_conv_dcot_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_dcot_branch does not exist!\n");
				exit(1);
			}
			els_conv_dcot_isLoaded = true;
		}
		return els_conv_dcot_;
	}
	vector<float> &els_conv_dist()
	{
		if (not els_conv_dist_isLoaded) {
			if (els_conv_dist_branch != 0) {
				els_conv_dist_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_dist_branch does not exist!\n");
				exit(1);
			}
			els_conv_dist_isLoaded = true;
		}
		return els_conv_dist_;
	}
	vector<float> &trks_conv_dcot()
	{
		if (not trks_conv_dcot_isLoaded) {
			if (trks_conv_dcot_branch != 0) {
				trks_conv_dcot_branch->GetEntry(index);
			} else { 
				printf("branch trks_conv_dcot_branch does not exist!\n");
				exit(1);
			}
			trks_conv_dcot_isLoaded = true;
		}
		return trks_conv_dcot_;
	}
	vector<float> &trks_conv_dist()
	{
		if (not trks_conv_dist_isLoaded) {
			if (trks_conv_dist_branch != 0) {
				trks_conv_dist_branch->GetEntry(index);
			} else { 
				printf("branch trks_conv_dist_branch does not exist!\n");
				exit(1);
			}
			trks_conv_dist_isLoaded = true;
		}
		return trks_conv_dist_;
	}
	vector<float> &els_musdr()
	{
		if (not els_musdr_isLoaded) {
			if (els_musdr_branch != 0) {
				els_musdr_branch->GetEntry(index);
			} else { 
				printf("branch els_musdr_branch does not exist!\n");
				exit(1);
			}
			els_musdr_isLoaded = true;
		}
		return els_musdr_;
	}
	vector<float> &els_trkdr()
	{
		if (not els_trkdr_isLoaded) {
			if (els_trkdr_branch != 0) {
				els_trkdr_branch->GetEntry(index);
			} else { 
				printf("branch els_trkdr_branch does not exist!\n");
				exit(1);
			}
			els_trkdr_isLoaded = true;
		}
		return els_trkdr_;
	}
	vector<float> &els_trkshFrac()
	{
		if (not els_trkshFrac_isLoaded) {
			if (els_trkshFrac_branch != 0) {
				els_trkshFrac_branch->GetEntry(index);
			} else { 
				printf("branch els_trkshFrac_branch does not exist!\n");
				exit(1);
			}
			els_trkshFrac_isLoaded = true;
		}
		return els_trkshFrac_;
	}
	vector<float> &els_chi2()
	{
		if (not els_chi2_isLoaded) {
			if (els_chi2_branch != 0) {
				els_chi2_branch->GetEntry(index);
			} else { 
				printf("branch els_chi2_branch does not exist!\n");
				exit(1);
			}
			els_chi2_isLoaded = true;
		}
		return els_chi2_;
	}
	vector<float> &els_d0()
	{
		if (not els_d0_isLoaded) {
			if (els_d0_branch != 0) {
				els_d0_branch->GetEntry(index);
			} else { 
				printf("branch els_d0_branch does not exist!\n");
				exit(1);
			}
			els_d0_isLoaded = true;
		}
		return els_d0_;
	}
	vector<float> &els_d0Err()
	{
		if (not els_d0Err_isLoaded) {
			if (els_d0Err_branch != 0) {
				els_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch els_d0Err_branch does not exist!\n");
				exit(1);
			}
			els_d0Err_isLoaded = true;
		}
		return els_d0Err_;
	}
	vector<float> &els_d0corr()
	{
		if (not els_d0corr_isLoaded) {
			if (els_d0corr_branch != 0) {
				els_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch els_d0corr_branch does not exist!\n");
				exit(1);
			}
			els_d0corr_isLoaded = true;
		}
		return els_d0corr_;
	}
	vector<float> &els_dEtaIn()
	{
		if (not els_dEtaIn_isLoaded) {
			if (els_dEtaIn_branch != 0) {
				els_dEtaIn_branch->GetEntry(index);
			} else { 
				printf("branch els_dEtaIn_branch does not exist!\n");
				exit(1);
			}
			els_dEtaIn_isLoaded = true;
		}
		return els_dEtaIn_;
	}
	vector<float> &els_dEtaOut()
	{
		if (not els_dEtaOut_isLoaded) {
			if (els_dEtaOut_branch != 0) {
				els_dEtaOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dEtaOut_branch does not exist!\n");
				exit(1);
			}
			els_dEtaOut_isLoaded = true;
		}
		return els_dEtaOut_;
	}
	vector<float> &els_dPhiIn()
	{
		if (not els_dPhiIn_isLoaded) {
			if (els_dPhiIn_branch != 0) {
				els_dPhiIn_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiIn_branch does not exist!\n");
				exit(1);
			}
			els_dPhiIn_isLoaded = true;
		}
		return els_dPhiIn_;
	}
	vector<float> &els_dPhiInPhiOut()
	{
		if (not els_dPhiInPhiOut_isLoaded) {
			if (els_dPhiInPhiOut_branch != 0) {
				els_dPhiInPhiOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiInPhiOut_branch does not exist!\n");
				exit(1);
			}
			els_dPhiInPhiOut_isLoaded = true;
		}
		return els_dPhiInPhiOut_;
	}
	vector<float> &els_dPhiOut()
	{
		if (not els_dPhiOut_isLoaded) {
			if (els_dPhiOut_branch != 0) {
				els_dPhiOut_branch->GetEntry(index);
			} else { 
				printf("branch els_dPhiOut_branch does not exist!\n");
				exit(1);
			}
			els_dPhiOut_isLoaded = true;
		}
		return els_dPhiOut_;
	}
	vector<float> &els_e1x5()
	{
		if (not els_e1x5_isLoaded) {
			if (els_e1x5_branch != 0) {
				els_e1x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e1x5_branch does not exist!\n");
				exit(1);
			}
			els_e1x5_isLoaded = true;
		}
		return els_e1x5_;
	}
	vector<float> &els_e2x5Max()
	{
		if (not els_e2x5Max_isLoaded) {
			if (els_e2x5Max_branch != 0) {
				els_e2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch els_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			els_e2x5Max_isLoaded = true;
		}
		return els_e2x5Max_;
	}
	vector<float> &els_e3x3()
	{
		if (not els_e3x3_isLoaded) {
			if (els_e3x3_branch != 0) {
				els_e3x3_branch->GetEntry(index);
			} else { 
				printf("branch els_e3x3_branch does not exist!\n");
				exit(1);
			}
			els_e3x3_isLoaded = true;
		}
		return els_e3x3_;
	}
	vector<float> &els_e5x5()
	{
		if (not els_e5x5_isLoaded) {
			if (els_e5x5_branch != 0) {
				els_e5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_e5x5_branch does not exist!\n");
				exit(1);
			}
			els_e5x5_isLoaded = true;
		}
		return els_e5x5_;
	}
	vector<float> &els_eMax()
	{
		if (not els_eMax_isLoaded) {
			if (els_eMax_branch != 0) {
				els_eMax_branch->GetEntry(index);
			} else { 
				printf("branch els_eMax_branch does not exist!\n");
				exit(1);
			}
			els_eMax_isLoaded = true;
		}
		return els_eMax_;
	}
	vector<float> &els_eOverPIn()
	{
		if (not els_eOverPIn_isLoaded) {
			if (els_eOverPIn_branch != 0) {
				els_eOverPIn_branch->GetEntry(index);
			} else { 
				printf("branch els_eOverPIn_branch does not exist!\n");
				exit(1);
			}
			els_eOverPIn_isLoaded = true;
		}
		return els_eOverPIn_;
	}
	vector<float> &els_eSC()
	{
		if (not els_eSC_isLoaded) {
			if (els_eSC_branch != 0) {
				els_eSC_branch->GetEntry(index);
			} else { 
				printf("branch els_eSC_branch does not exist!\n");
				exit(1);
			}
			els_eSC_isLoaded = true;
		}
		return els_eSC_;
	}
	vector<float> &els_eSCPresh()
	{
		if (not els_eSCPresh_isLoaded) {
			if (els_eSCPresh_branch != 0) {
				els_eSCPresh_branch->GetEntry(index);
			} else { 
				printf("branch els_eSCPresh_branch does not exist!\n");
				exit(1);
			}
			els_eSCPresh_isLoaded = true;
		}
		return els_eSCPresh_;
	}
	vector<float> &els_eSCRaw()
	{
		if (not els_eSCRaw_isLoaded) {
			if (els_eSCRaw_branch != 0) {
				els_eSCRaw_branch->GetEntry(index);
			} else { 
				printf("branch els_eSCRaw_branch does not exist!\n");
				exit(1);
			}
			els_eSCRaw_isLoaded = true;
		}
		return els_eSCRaw_;
	}
	vector<float> &els_eSeed()
	{
		if (not els_eSeed_isLoaded) {
			if (els_eSeed_branch != 0) {
				els_eSeed_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeed_branch does not exist!\n");
				exit(1);
			}
			els_eSeed_isLoaded = true;
		}
		return els_eSeed_;
	}
	vector<float> &els_eSeedOverPIn()
	{
		if (not els_eSeedOverPIn_isLoaded) {
			if (els_eSeedOverPIn_branch != 0) {
				els_eSeedOverPIn_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeedOverPIn_branch does not exist!\n");
				exit(1);
			}
			els_eSeedOverPIn_isLoaded = true;
		}
		return els_eSeedOverPIn_;
	}
	vector<float> &els_eSeedOverPOut()
	{
		if (not els_eSeedOverPOut_isLoaded) {
			if (els_eSeedOverPOut_branch != 0) {
				els_eSeedOverPOut_branch->GetEntry(index);
			} else { 
				printf("branch els_eSeedOverPOut_branch does not exist!\n");
				exit(1);
			}
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	vector<float> &els_ecalIso03()
	{
		if (not els_ecalIso03_isLoaded) {
			if (els_ecalIso03_branch != 0) {
				els_ecalIso03_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalIso03_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso03_isLoaded = true;
		}
		return els_ecalIso03_;
	}
	vector<float> &els_ecalIso04()
	{
		if (not els_ecalIso04_isLoaded) {
			if (els_ecalIso04_branch != 0) {
				els_ecalIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalIso04_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso04_isLoaded = true;
		}
		return els_ecalIso04_;
	}
	vector<float> &els_egamma_looseId()
	{
		if (not els_egamma_looseId_isLoaded) {
			if (els_egamma_looseId_branch != 0) {
				els_egamma_looseId_branch->GetEntry(index);
			} else { 
				printf("branch els_egamma_looseId_branch does not exist!\n");
				exit(1);
			}
			els_egamma_looseId_isLoaded = true;
		}
		return els_egamma_looseId_;
	}
	vector<float> &els_egamma_robustHighEnergy()
	{
		if (not els_egamma_robustHighEnergy_isLoaded) {
			if (els_egamma_robustHighEnergy_branch != 0) {
				els_egamma_robustHighEnergy_branch->GetEntry(index);
			} else { 
				printf("branch els_egamma_robustHighEnergy_branch does not exist!\n");
				exit(1);
			}
			els_egamma_robustHighEnergy_isLoaded = true;
		}
		return els_egamma_robustHighEnergy_;
	}
	vector<float> &els_egamma_robustLooseId()
	{
		if (not els_egamma_robustLooseId_isLoaded) {
			if (els_egamma_robustLooseId_branch != 0) {
				els_egamma_robustLooseId_branch->GetEntry(index);
			} else { 
				printf("branch els_egamma_robustLooseId_branch does not exist!\n");
				exit(1);
			}
			els_egamma_robustLooseId_isLoaded = true;
		}
		return els_egamma_robustLooseId_;
	}
	vector<float> &els_egamma_robustTightId()
	{
		if (not els_egamma_robustTightId_isLoaded) {
			if (els_egamma_robustTightId_branch != 0) {
				els_egamma_robustTightId_branch->GetEntry(index);
			} else { 
				printf("branch els_egamma_robustTightId_branch does not exist!\n");
				exit(1);
			}
			els_egamma_robustTightId_isLoaded = true;
		}
		return els_egamma_robustTightId_;
	}
	vector<float> &els_egamma_tightId()
	{
		if (not els_egamma_tightId_isLoaded) {
			if (els_egamma_tightId_branch != 0) {
				els_egamma_tightId_branch->GetEntry(index);
			} else { 
				printf("branch els_egamma_tightId_branch does not exist!\n");
				exit(1);
			}
			els_egamma_tightId_isLoaded = true;
		}
		return els_egamma_tightId_;
	}
	vector<float> &els_etaErr()
	{
		if (not els_etaErr_isLoaded) {
			if (els_etaErr_branch != 0) {
				els_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch els_etaErr_branch does not exist!\n");
				exit(1);
			}
			els_etaErr_isLoaded = true;
		}
		return els_etaErr_;
	}
	vector<float> &els_etaSC()
	{
		if (not els_etaSC_isLoaded) {
			if (els_etaSC_branch != 0) {
				els_etaSC_branch->GetEntry(index);
			} else { 
				printf("branch els_etaSC_branch does not exist!\n");
				exit(1);
			}
			els_etaSC_isLoaded = true;
		}
		return els_etaSC_;
	}
	vector<float> &els_fBrem()
	{
		if (not els_fBrem_isLoaded) {
			if (els_fBrem_branch != 0) {
				els_fBrem_branch->GetEntry(index);
			} else { 
				printf("branch els_fBrem_branch does not exist!\n");
				exit(1);
			}
			els_fBrem_isLoaded = true;
		}
		return els_fBrem_;
	}
	vector<float> &els_hOverE()
	{
		if (not els_hOverE_isLoaded) {
			if (els_hOverE_branch != 0) {
				els_hOverE_branch->GetEntry(index);
			} else { 
				printf("branch els_hOverE_branch does not exist!\n");
				exit(1);
			}
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	vector<float> &els_hcalIso03()
	{
		if (not els_hcalIso03_isLoaded) {
			if (els_hcalIso03_branch != 0) {
				els_hcalIso03_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalIso03_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso03_isLoaded = true;
		}
		return els_hcalIso03_;
	}
	vector<float> &els_hcalIso04()
	{
		if (not els_hcalIso04_isLoaded) {
			if (els_hcalIso04_branch != 0) {
				els_hcalIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalIso04_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso04_isLoaded = true;
		}
		return els_hcalIso04_;
	}
	vector<float> &els_layer1_charge()
	{
		if (not els_layer1_charge_isLoaded) {
			if (els_layer1_charge_branch != 0) {
				els_layer1_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_layer1_charge_branch does not exist!\n");
				exit(1);
			}
			els_layer1_charge_isLoaded = true;
		}
		return els_layer1_charge_;
	}
	vector<float> &els_ndof()
	{
		if (not els_ndof_isLoaded) {
			if (els_ndof_branch != 0) {
				els_ndof_branch->GetEntry(index);
			} else { 
				printf("branch els_ndof_branch does not exist!\n");
				exit(1);
			}
			els_ndof_isLoaded = true;
		}
		return els_ndof_;
	}
	vector<float> &els_outerEta()
	{
		if (not els_outerEta_isLoaded) {
			if (els_outerEta_branch != 0) {
				els_outerEta_branch->GetEntry(index);
			} else { 
				printf("branch els_outerEta_branch does not exist!\n");
				exit(1);
			}
			els_outerEta_isLoaded = true;
		}
		return els_outerEta_;
	}
	vector<float> &els_outerPhi()
	{
		if (not els_outerPhi_isLoaded) {
			if (els_outerPhi_branch != 0) {
				els_outerPhi_branch->GetEntry(index);
			} else { 
				printf("branch els_outerPhi_branch does not exist!\n");
				exit(1);
			}
			els_outerPhi_isLoaded = true;
		}
		return els_outerPhi_;
	}
	vector<float> &els_phiErr()
	{
		if (not els_phiErr_isLoaded) {
			if (els_phiErr_branch != 0) {
				els_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch els_phiErr_branch does not exist!\n");
				exit(1);
			}
			els_phiErr_isLoaded = true;
		}
		return els_phiErr_;
	}
	vector<float> &els_phiSC()
	{
		if (not els_phiSC_isLoaded) {
			if (els_phiSC_branch != 0) {
				els_phiSC_branch->GetEntry(index);
			} else { 
				printf("branch els_phiSC_branch does not exist!\n");
				exit(1);
			}
			els_phiSC_isLoaded = true;
		}
		return els_phiSC_;
	}
	vector<float> &els_pin()
	{
		if (not els_pin_isLoaded) {
			if (els_pin_branch != 0) {
				els_pin_branch->GetEntry(index);
			} else { 
				printf("branch els_pin_branch does not exist!\n");
				exit(1);
			}
			els_pin_isLoaded = true;
		}
		return els_pin_;
	}
	vector<float> &els_pout()
	{
		if (not els_pout_isLoaded) {
			if (els_pout_branch != 0) {
				els_pout_branch->GetEntry(index);
			} else { 
				printf("branch els_pout_branch does not exist!\n");
				exit(1);
			}
			els_pout_isLoaded = true;
		}
		return els_pout_;
	}
	vector<float> &els_ptErr()
	{
		if (not els_ptErr_isLoaded) {
			if (els_ptErr_branch != 0) {
				els_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch els_ptErr_branch does not exist!\n");
				exit(1);
			}
			els_ptErr_isLoaded = true;
		}
		return els_ptErr_;
	}
	vector<float> &els_sigmaEtaEta()
	{
		if (not els_sigmaEtaEta_isLoaded) {
			if (els_sigmaEtaEta_branch != 0) {
				els_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			els_sigmaEtaEta_isLoaded = true;
		}
		return els_sigmaEtaEta_;
	}
	vector<float> &els_sigmaIEtaIEta()
	{
		if (not els_sigmaIEtaIEta_isLoaded) {
			if (els_sigmaIEtaIEta_branch != 0) {
				els_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIEtaIEta_isLoaded = true;
		}
		return els_sigmaIEtaIEta_;
	}
	vector<float> &els_sigmaIPhiIPhi()
	{
		if (not els_sigmaIPhiIPhi_isLoaded) {
			if (els_sigmaIPhiIPhi_branch != 0) {
				els_sigmaIPhiIPhi_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIPhiIPhi_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIPhiIPhi_isLoaded = true;
		}
		return els_sigmaIPhiIPhi_;
	}
	vector<float> &els_sigmaPhiPhi()
	{
		if (not els_sigmaPhiPhi_isLoaded) {
			if (els_sigmaPhiPhi_branch != 0) {
				els_sigmaPhiPhi_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaPhiPhi_branch does not exist!\n");
				exit(1);
			}
			els_sigmaPhiPhi_isLoaded = true;
		}
		return els_sigmaPhiPhi_;
	}
	vector<float> &els_tkIso03()
	{
		if (not els_tkIso03_isLoaded) {
			if (els_tkIso03_branch != 0) {
				els_tkIso03_branch->GetEntry(index);
			} else { 
				printf("branch els_tkIso03_branch does not exist!\n");
				exit(1);
			}
			els_tkIso03_isLoaded = true;
		}
		return els_tkIso03_;
	}
	vector<float> &els_tkIso04()
	{
		if (not els_tkIso04_isLoaded) {
			if (els_tkIso04_branch != 0) {
				els_tkIso04_branch->GetEntry(index);
			} else { 
				printf("branch els_tkIso04_branch does not exist!\n");
				exit(1);
			}
			els_tkIso04_isLoaded = true;
		}
		return els_tkIso04_;
	}
	vector<float> &els_vertexphi()
	{
		if (not els_vertexphi_isLoaded) {
			if (els_vertexphi_branch != 0) {
				els_vertexphi_branch->GetEntry(index);
			} else { 
				printf("branch els_vertexphi_branch does not exist!\n");
				exit(1);
			}
			els_vertexphi_isLoaded = true;
		}
		return els_vertexphi_;
	}
	vector<float> &els_z0()
	{
		if (not els_z0_isLoaded) {
			if (els_z0_branch != 0) {
				els_z0_branch->GetEntry(index);
			} else { 
				printf("branch els_z0_branch does not exist!\n");
				exit(1);
			}
			els_z0_isLoaded = true;
		}
		return els_z0_;
	}
	vector<float> &els_z0Err()
	{
		if (not els_z0Err_isLoaded) {
			if (els_z0Err_branch != 0) {
				els_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch els_z0Err_branch does not exist!\n");
				exit(1);
			}
			els_z0Err_isLoaded = true;
		}
		return els_z0Err_;
	}
	vector<float> &els_z0corr()
	{
		if (not els_z0corr_isLoaded) {
			if (els_z0corr_branch != 0) {
				els_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch els_z0corr_branch does not exist!\n");
				exit(1);
			}
			els_z0corr_isLoaded = true;
		}
		return els_z0corr_;
	}
	vector<float> &mus_caloCompatibility()
	{
		if (not mus_caloCompatibility_isLoaded) {
			if (mus_caloCompatibility_branch != 0) {
				mus_caloCompatibility_branch->GetEntry(index);
			} else { 
				printf("branch mus_caloCompatibility_branch does not exist!\n");
				exit(1);
			}
			mus_caloCompatibility_isLoaded = true;
		}
		return mus_caloCompatibility_;
	}
	vector<float> &mus_chi2()
	{
		if (not mus_chi2_isLoaded) {
			if (mus_chi2_branch != 0) {
				mus_chi2_branch->GetEntry(index);
			} else { 
				printf("branch mus_chi2_branch does not exist!\n");
				exit(1);
			}
			mus_chi2_isLoaded = true;
		}
		return mus_chi2_;
	}
	vector<float> &mus_d0()
	{
		if (not mus_d0_isLoaded) {
			if (mus_d0_branch != 0) {
				mus_d0_branch->GetEntry(index);
			} else { 
				printf("branch mus_d0_branch does not exist!\n");
				exit(1);
			}
			mus_d0_isLoaded = true;
		}
		return mus_d0_;
	}
	vector<float> &mus_d0Err()
	{
		if (not mus_d0Err_isLoaded) {
			if (mus_d0Err_branch != 0) {
				mus_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_d0Err_branch does not exist!\n");
				exit(1);
			}
			mus_d0Err_isLoaded = true;
		}
		return mus_d0Err_;
	}
	vector<float> &mus_d0corr()
	{
		if (not mus_d0corr_isLoaded) {
			if (mus_d0corr_branch != 0) {
				mus_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_d0corr_branch does not exist!\n");
				exit(1);
			}
			mus_d0corr_isLoaded = true;
		}
		return mus_d0corr_;
	}
	vector<float> &mus_e_em()
	{
		if (not mus_e_em_isLoaded) {
			if (mus_e_em_branch != 0) {
				mus_e_em_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_em_branch does not exist!\n");
				exit(1);
			}
			mus_e_em_isLoaded = true;
		}
		return mus_e_em_;
	}
	vector<float> &mus_e_emS9()
	{
		if (not mus_e_emS9_isLoaded) {
			if (mus_e_emS9_branch != 0) {
				mus_e_emS9_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_emS9_branch does not exist!\n");
				exit(1);
			}
			mus_e_emS9_isLoaded = true;
		}
		return mus_e_emS9_;
	}
	vector<float> &mus_e_had()
	{
		if (not mus_e_had_isLoaded) {
			if (mus_e_had_branch != 0) {
				mus_e_had_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_had_branch does not exist!\n");
				exit(1);
			}
			mus_e_had_isLoaded = true;
		}
		return mus_e_had_;
	}
	vector<float> &mus_e_hadS9()
	{
		if (not mus_e_hadS9_isLoaded) {
			if (mus_e_hadS9_branch != 0) {
				mus_e_hadS9_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_hadS9_branch does not exist!\n");
				exit(1);
			}
			mus_e_hadS9_isLoaded = true;
		}
		return mus_e_hadS9_;
	}
	vector<float> &mus_e_ho()
	{
		if (not mus_e_ho_isLoaded) {
			if (mus_e_ho_branch != 0) {
				mus_e_ho_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_ho_branch does not exist!\n");
				exit(1);
			}
			mus_e_ho_isLoaded = true;
		}
		return mus_e_ho_;
	}
	vector<float> &mus_e_hoS9()
	{
		if (not mus_e_hoS9_isLoaded) {
			if (mus_e_hoS9_branch != 0) {
				mus_e_hoS9_branch->GetEntry(index);
			} else { 
				printf("branch mus_e_hoS9_branch does not exist!\n");
				exit(1);
			}
			mus_e_hoS9_isLoaded = true;
		}
		return mus_e_hoS9_;
	}
	vector<float> &mus_etaErr()
	{
		if (not mus_etaErr_isLoaded) {
			if (mus_etaErr_branch != 0) {
				mus_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch mus_etaErr_branch does not exist!\n");
				exit(1);
			}
			mus_etaErr_isLoaded = true;
		}
		return mus_etaErr_;
	}
	vector<float> &mus_gfit_chi2()
	{
		if (not mus_gfit_chi2_isLoaded) {
			if (mus_gfit_chi2_branch != 0) {
				mus_gfit_chi2_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_chi2_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_chi2_isLoaded = true;
		}
		return mus_gfit_chi2_;
	}
	vector<float> &mus_gfit_ndof()
	{
		if (not mus_gfit_ndof_isLoaded) {
			if (mus_gfit_ndof_branch != 0) {
				mus_gfit_ndof_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_ndof_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_ndof_isLoaded = true;
		}
		return mus_gfit_ndof_;
	}
	vector<float> &mus_iso()
	{
		if (not mus_iso_isLoaded) {
			if (mus_iso_branch != 0) {
				mus_iso_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso_branch does not exist!\n");
				exit(1);
			}
			mus_iso_isLoaded = true;
		}
		return mus_iso_;
	}
	vector<float> &mus_iso03_emEt()
	{
		if (not mus_iso03_emEt_isLoaded) {
			if (mus_iso03_emEt_branch != 0) {
				mus_iso03_emEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso03_emEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso03_emEt_isLoaded = true;
		}
		return mus_iso03_emEt_;
	}
	vector<float> &mus_iso03_hadEt()
	{
		if (not mus_iso03_hadEt_isLoaded) {
			if (mus_iso03_hadEt_branch != 0) {
				mus_iso03_hadEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso03_hadEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso03_hadEt_isLoaded = true;
		}
		return mus_iso03_hadEt_;
	}
	vector<float> &mus_iso03_hoEt()
	{
		if (not mus_iso03_hoEt_isLoaded) {
			if (mus_iso03_hoEt_branch != 0) {
				mus_iso03_hoEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso03_hoEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso03_hoEt_isLoaded = true;
		}
		return mus_iso03_hoEt_;
	}
	vector<float> &mus_iso03_sumPt()
	{
		if (not mus_iso03_sumPt_isLoaded) {
			if (mus_iso03_sumPt_branch != 0) {
				mus_iso03_sumPt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso03_sumPt_branch does not exist!\n");
				exit(1);
			}
			mus_iso03_sumPt_isLoaded = true;
		}
		return mus_iso03_sumPt_;
	}
	vector<float> &mus_iso05_emEt()
	{
		if (not mus_iso05_emEt_isLoaded) {
			if (mus_iso05_emEt_branch != 0) {
				mus_iso05_emEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso05_emEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso05_emEt_isLoaded = true;
		}
		return mus_iso05_emEt_;
	}
	vector<float> &mus_iso05_hadEt()
	{
		if (not mus_iso05_hadEt_isLoaded) {
			if (mus_iso05_hadEt_branch != 0) {
				mus_iso05_hadEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso05_hadEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso05_hadEt_isLoaded = true;
		}
		return mus_iso05_hadEt_;
	}
	vector<float> &mus_iso05_hoEt()
	{
		if (not mus_iso05_hoEt_isLoaded) {
			if (mus_iso05_hoEt_branch != 0) {
				mus_iso05_hoEt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso05_hoEt_branch does not exist!\n");
				exit(1);
			}
			mus_iso05_hoEt_isLoaded = true;
		}
		return mus_iso05_hoEt_;
	}
	vector<float> &mus_iso05_sumPt()
	{
		if (not mus_iso05_sumPt_isLoaded) {
			if (mus_iso05_sumPt_branch != 0) {
				mus_iso05_sumPt_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso05_sumPt_branch does not exist!\n");
				exit(1);
			}
			mus_iso05_sumPt_isLoaded = true;
		}
		return mus_iso05_sumPt_;
	}
	vector<float> &mus_ndof()
	{
		if (not mus_ndof_isLoaded) {
			if (mus_ndof_branch != 0) {
				mus_ndof_branch->GetEntry(index);
			} else { 
				printf("branch mus_ndof_branch does not exist!\n");
				exit(1);
			}
			mus_ndof_isLoaded = true;
		}
		return mus_ndof_;
	}
	vector<float> &mus_outerEta()
	{
		if (not mus_outerEta_isLoaded) {
			if (mus_outerEta_branch != 0) {
				mus_outerEta_branch->GetEntry(index);
			} else { 
				printf("branch mus_outerEta_branch does not exist!\n");
				exit(1);
			}
			mus_outerEta_isLoaded = true;
		}
		return mus_outerEta_;
	}
	vector<float> &mus_outerPhi()
	{
		if (not mus_outerPhi_isLoaded) {
			if (mus_outerPhi_branch != 0) {
				mus_outerPhi_branch->GetEntry(index);
			} else { 
				printf("branch mus_outerPhi_branch does not exist!\n");
				exit(1);
			}
			mus_outerPhi_isLoaded = true;
		}
		return mus_outerPhi_;
	}
	vector<float> &mus_phiErr()
	{
		if (not mus_phiErr_isLoaded) {
			if (mus_phiErr_branch != 0) {
				mus_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch mus_phiErr_branch does not exist!\n");
				exit(1);
			}
			mus_phiErr_isLoaded = true;
		}
		return mus_phiErr_;
	}
	vector<float> &mus_ptErr()
	{
		if (not mus_ptErr_isLoaded) {
			if (mus_ptErr_branch != 0) {
				mus_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch mus_ptErr_branch does not exist!\n");
				exit(1);
			}
			mus_ptErr_isLoaded = true;
		}
		return mus_ptErr_;
	}
	vector<float> &mus_qoverp()
	{
		if (not mus_qoverp_isLoaded) {
			if (mus_qoverp_branch != 0) {
				mus_qoverp_branch->GetEntry(index);
			} else { 
				printf("branch mus_qoverp_branch does not exist!\n");
				exit(1);
			}
			mus_qoverp_isLoaded = true;
		}
		return mus_qoverp_;
	}
	vector<float> &mus_qoverpError()
	{
		if (not mus_qoverpError_isLoaded) {
			if (mus_qoverpError_branch != 0) {
				mus_qoverpError_branch->GetEntry(index);
			} else { 
				printf("branch mus_qoverpError_branch does not exist!\n");
				exit(1);
			}
			mus_qoverpError_isLoaded = true;
		}
		return mus_qoverpError_;
	}
	vector<float> &mus_vertexphi()
	{
		if (not mus_vertexphi_isLoaded) {
			if (mus_vertexphi_branch != 0) {
				mus_vertexphi_branch->GetEntry(index);
			} else { 
				printf("branch mus_vertexphi_branch does not exist!\n");
				exit(1);
			}
			mus_vertexphi_isLoaded = true;
		}
		return mus_vertexphi_;
	}
	vector<float> &mus_z0()
	{
		if (not mus_z0_isLoaded) {
			if (mus_z0_branch != 0) {
				mus_z0_branch->GetEntry(index);
			} else { 
				printf("branch mus_z0_branch does not exist!\n");
				exit(1);
			}
			mus_z0_isLoaded = true;
		}
		return mus_z0_;
	}
	vector<float> &mus_z0Err()
	{
		if (not mus_z0Err_isLoaded) {
			if (mus_z0Err_branch != 0) {
				mus_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_z0Err_branch does not exist!\n");
				exit(1);
			}
			mus_z0Err_isLoaded = true;
		}
		return mus_z0Err_;
	}
	vector<float> &mus_z0corr()
	{
		if (not mus_z0corr_isLoaded) {
			if (mus_z0corr_branch != 0) {
				mus_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_z0corr_branch does not exist!\n");
				exit(1);
			}
			mus_z0corr_isLoaded = true;
		}
		return mus_z0corr_;
	}
	vector<float> &trks_chi2()
	{
		if (not trks_chi2_isLoaded) {
			if (trks_chi2_branch != 0) {
				trks_chi2_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_ndof_branch does not exist!\n");
				exit(1);
			}
			trks_ndof_isLoaded = true;
		}
		return trks_ndof_;
	}
	vector<float> &trks_outerEt()
	{
		if (not trks_outerEt_isLoaded) {
			if (trks_outerEt_branch != 0) {
				trks_outerEt_branch->GetEntry(index);
			} else { 
				printf("branch trks_outerEt_branch does not exist!\n");
				exit(1);
			}
			trks_outerEt_isLoaded = true;
		}
		return trks_outerEt_;
	}
	vector<float> &trks_outerEta()
	{
		if (not trks_outerEta_isLoaded) {
			if (trks_outerEta_branch != 0) {
				trks_outerEta_branch->GetEntry(index);
			} else { 
				printf("branch trks_outerEta_branch does not exist!\n");
				exit(1);
			}
			trks_outerEta_isLoaded = true;
		}
		return trks_outerEta_;
	}
	vector<float> &trks_outerPhi()
	{
		if (not trks_outerPhi_isLoaded) {
			if (trks_outerPhi_branch != 0) {
				trks_outerPhi_branch->GetEntry(index);
			} else { 
				printf("branch trks_outerPhi_branch does not exist!\n");
				exit(1);
			}
			trks_outerPhi_isLoaded = true;
		}
		return trks_outerPhi_;
	}
	vector<float> &trks_phiErr()
	{
		if (not trks_phiErr_isLoaded) {
			if (trks_phiErr_branch != 0) {
				trks_phiErr_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_ptErr_branch does not exist!\n");
				exit(1);
			}
			trks_ptErr_isLoaded = true;
		}
		return trks_ptErr_;
	}
	vector<float> &trks_tkIso()
	{
		if (not trks_tkIso_isLoaded) {
			if (trks_tkIso_branch != 0) {
				trks_tkIso_branch->GetEntry(index);
			} else { 
				printf("branch trks_tkIso_branch does not exist!\n");
				exit(1);
			}
			trks_tkIso_isLoaded = true;
		}
		return trks_tkIso_;
	}
	vector<float> &trks_vertexphi()
	{
		if (not trks_vertexphi_isLoaded) {
			if (trks_vertexphi_branch != 0) {
				trks_vertexphi_branch->GetEntry(index);
			} else { 
				printf("branch trks_vertexphi_branch does not exist!\n");
				exit(1);
			}
			trks_vertexphi_isLoaded = true;
		}
		return trks_vertexphi_;
	}
	vector<float> &trks_z0()
	{
		if (not trks_z0_isLoaded) {
			if (trks_z0_branch != 0) {
				trks_z0_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_z0corr_branch does not exist!\n");
				exit(1);
			}
			trks_z0corr_isLoaded = true;
		}
		return trks_z0corr_;
	}
	vector<float> &trks_elsdr()
	{
		if (not trks_elsdr_isLoaded) {
			if (trks_elsdr_branch != 0) {
				trks_elsdr_branch->GetEntry(index);
			} else { 
				printf("branch trks_elsdr_branch does not exist!\n");
				exit(1);
			}
			trks_elsdr_isLoaded = true;
		}
		return trks_elsdr_;
	}
	vector<float> &trks_elsshFrac()
	{
		if (not trks_elsshFrac_isLoaded) {
			if (trks_elsshFrac_branch != 0) {
				trks_elsshFrac_branch->GetEntry(index);
			} else { 
				printf("branch trks_elsshFrac_branch does not exist!\n");
				exit(1);
			}
			trks_elsshFrac_isLoaded = true;
		}
		return trks_elsshFrac_;
	}
	vector<float> &trk_musdr()
	{
		if (not trk_musdr_isLoaded) {
			if (trk_musdr_branch != 0) {
				trk_musdr_branch->GetEntry(index);
			} else { 
				printf("branch trk_musdr_branch does not exist!\n");
				exit(1);
			}
			trk_musdr_isLoaded = true;
		}
		return trk_musdr_;
	}
	vector<float> &vtxs_chi2()
	{
		if (not vtxs_chi2_isLoaded) {
			if (vtxs_chi2_branch != 0) {
				vtxs_chi2_branch->GetEntry(index);
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
			} else { 
				printf("branch vtxs_zError_branch does not exist!\n");
				exit(1);
			}
			vtxs_zError_isLoaded = true;
		}
		return vtxs_zError_;
	}
	vector<vector<float> > &trks_residualX()
	{
		if (not trks_residualX_isLoaded) {
			if (trks_residualX_branch != 0) {
				trks_residualX_branch->GetEntry(index);
			} else { 
				printf("branch trks_residualX_branch does not exist!\n");
				exit(1);
			}
			trks_residualX_isLoaded = true;
		}
		return trks_residualX_;
	}
	vector<vector<float> > &trks_residualY()
	{
		if (not trks_residualY_isLoaded) {
			if (trks_residualY_branch != 0) {
				trks_residualY_branch->GetEntry(index);
			} else { 
				printf("branch trks_residualY_branch does not exist!\n");
				exit(1);
			}
			trks_residualY_isLoaded = true;
		}
		return trks_residualY_;
	}
	int &evt_HLT1()
	{
		if (not evt_HLT1_isLoaded) {
			if (evt_HLT1_branch != 0) {
				evt_HLT1_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT1_branch does not exist!\n");
				exit(1);
			}
			evt_HLT1_isLoaded = true;
		}
		return evt_HLT1_;
	}
	int &evt_HLT2()
	{
		if (not evt_HLT2_isLoaded) {
			if (evt_HLT2_branch != 0) {
				evt_HLT2_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT2_branch does not exist!\n");
				exit(1);
			}
			evt_HLT2_isLoaded = true;
		}
		return evt_HLT2_;
	}
	int &evt_HLT3()
	{
		if (not evt_HLT3_isLoaded) {
			if (evt_HLT3_branch != 0) {
				evt_HLT3_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT3_branch does not exist!\n");
				exit(1);
			}
			evt_HLT3_isLoaded = true;
		}
		return evt_HLT3_;
	}
	int &evt_HLT4()
	{
		if (not evt_HLT4_isLoaded) {
			if (evt_HLT4_branch != 0) {
				evt_HLT4_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT4_branch does not exist!\n");
				exit(1);
			}
			evt_HLT4_isLoaded = true;
		}
		return evt_HLT4_;
	}
	int &evt_HLT5()
	{
		if (not evt_HLT5_isLoaded) {
			if (evt_HLT5_branch != 0) {
				evt_HLT5_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT5_branch does not exist!\n");
				exit(1);
			}
			evt_HLT5_isLoaded = true;
		}
		return evt_HLT5_;
	}
	int &evt_HLT6()
	{
		if (not evt_HLT6_isLoaded) {
			if (evt_HLT6_branch != 0) {
				evt_HLT6_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT6_branch does not exist!\n");
				exit(1);
			}
			evt_HLT6_isLoaded = true;
		}
		return evt_HLT6_;
	}
	int &evt_HLT7()
	{
		if (not evt_HLT7_isLoaded) {
			if (evt_HLT7_branch != 0) {
				evt_HLT7_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT7_branch does not exist!\n");
				exit(1);
			}
			evt_HLT7_isLoaded = true;
		}
		return evt_HLT7_;
	}
	int &evt_HLT8()
	{
		if (not evt_HLT8_isLoaded) {
			if (evt_HLT8_branch != 0) {
				evt_HLT8_branch->GetEntry(index);
			} else { 
				printf("branch evt_HLT8_branch does not exist!\n");
				exit(1);
			}
			evt_HLT8_isLoaded = true;
		}
		return evt_HLT8_;
	}
	int &evt_L1_1()
	{
		if (not evt_L1_1_isLoaded) {
			if (evt_L1_1_branch != 0) {
				evt_L1_1_branch->GetEntry(index);
			} else { 
				printf("branch evt_L1_1_branch does not exist!\n");
				exit(1);
			}
			evt_L1_1_isLoaded = true;
		}
		return evt_L1_1_;
	}
	int &evt_L1_2()
	{
		if (not evt_L1_2_isLoaded) {
			if (evt_L1_2_branch != 0) {
				evt_L1_2_branch->GetEntry(index);
			} else { 
				printf("branch evt_L1_2_branch does not exist!\n");
				exit(1);
			}
			evt_L1_2_isLoaded = true;
		}
		return evt_L1_2_;
	}
	int &evt_L1_3()
	{
		if (not evt_L1_3_isLoaded) {
			if (evt_L1_3_branch != 0) {
				evt_L1_3_branch->GetEntry(index);
			} else { 
				printf("branch evt_L1_3_branch does not exist!\n");
				exit(1);
			}
			evt_L1_3_isLoaded = true;
		}
		return evt_L1_3_;
	}
	int &evt_L1_4()
	{
		if (not evt_L1_4_isLoaded) {
			if (evt_L1_4_branch != 0) {
				evt_L1_4_branch->GetEntry(index);
			} else { 
				printf("branch evt_L1_4_branch does not exist!\n");
				exit(1);
			}
			evt_L1_4_isLoaded = true;
		}
		return evt_L1_4_;
	}
	vector<int> &els_mc3_id()
	{
		if (not els_mc3_id_isLoaded) {
			if (els_mc3_id_branch != 0) {
				els_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3_id_branch does not exist!\n");
				exit(1);
			}
			els_mc3_id_isLoaded = true;
		}
		return els_mc3_id_;
	}
	vector<int> &els_mc3idx()
	{
		if (not els_mc3idx_isLoaded) {
			if (els_mc3idx_branch != 0) {
				els_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3idx_branch does not exist!\n");
				exit(1);
			}
			els_mc3idx_isLoaded = true;
		}
		return els_mc3idx_;
	}
	vector<int> &els_mc3_motherid()
	{
		if (not els_mc3_motherid_isLoaded) {
			if (els_mc3_motherid_branch != 0) {
				els_mc3_motherid_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3_motherid_branch does not exist!\n");
				exit(1);
			}
			els_mc3_motherid_isLoaded = true;
		}
		return els_mc3_motherid_;
	}
	vector<int> &els_mc_id()
	{
		if (not els_mc_id_isLoaded) {
			if (els_mc_id_branch != 0) {
				els_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_id_branch does not exist!\n");
				exit(1);
			}
			els_mc_id_isLoaded = true;
		}
		return els_mc_id_;
	}
	vector<int> &els_mcidx()
	{
		if (not els_mcidx_isLoaded) {
			if (els_mcidx_branch != 0) {
				els_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch els_mcidx_branch does not exist!\n");
				exit(1);
			}
			els_mcidx_isLoaded = true;
		}
		return els_mcidx_;
	}
	vector<int> &els_mc_motherid()
	{
		if (not els_mc_motherid_isLoaded) {
			if (els_mc_motherid_branch != 0) {
				els_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch els_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			els_mc_motherid_isLoaded = true;
		}
		return els_mc_motherid_;
	}
	vector<int> &jets_mc_id()
	{
		if (not jets_mc_id_isLoaded) {
			if (jets_mc_id_branch != 0) {
				jets_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_id_branch does not exist!\n");
				exit(1);
			}
			jets_mc_id_isLoaded = true;
		}
		return jets_mc_id_;
	}
	vector<int> &mus_mc3_id()
	{
		if (not mus_mc3_id_isLoaded) {
			if (mus_mc3_id_branch != 0) {
				mus_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3_id_branch does not exist!\n");
				exit(1);
			}
			mus_mc3_id_isLoaded = true;
		}
		return mus_mc3_id_;
	}
	vector<int> &mus_mc3idx()
	{
		if (not mus_mc3idx_isLoaded) {
			if (mus_mc3idx_branch != 0) {
				mus_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3idx_branch does not exist!\n");
				exit(1);
			}
			mus_mc3idx_isLoaded = true;
		}
		return mus_mc3idx_;
	}
	vector<int> &mus_mc3_motherid()
	{
		if (not mus_mc3_motherid_isLoaded) {
			if (mus_mc3_motherid_branch != 0) {
				mus_mc3_motherid_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3_motherid_branch does not exist!\n");
				exit(1);
			}
			mus_mc3_motherid_isLoaded = true;
		}
		return mus_mc3_motherid_;
	}
	vector<int> &mus_mc_id()
	{
		if (not mus_mc_id_isLoaded) {
			if (mus_mc_id_branch != 0) {
				mus_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc_id_branch does not exist!\n");
				exit(1);
			}
			mus_mc_id_isLoaded = true;
		}
		return mus_mc_id_;
	}
	vector<int> &mus_mcidx()
	{
		if (not mus_mcidx_isLoaded) {
			if (mus_mcidx_branch != 0) {
				mus_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch mus_mcidx_branch does not exist!\n");
				exit(1);
			}
			mus_mcidx_isLoaded = true;
		}
		return mus_mcidx_;
	}
	vector<int> &mus_mc_motherid()
	{
		if (not mus_mc_motherid_isLoaded) {
			if (mus_mc_motherid_branch != 0) {
				mus_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			mus_mc_motherid_isLoaded = true;
		}
		return mus_mc_motherid_;
	}
	vector<int> &trk_mc3_id()
	{
		if (not trk_mc3_id_isLoaded) {
			if (trk_mc3_id_branch != 0) {
				trk_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3_id_branch does not exist!\n");
				exit(1);
			}
			trk_mc3_id_isLoaded = true;
		}
		return trk_mc3_id_;
	}
	vector<int> &trk_mc3idx()
	{
		if (not trk_mc3idx_isLoaded) {
			if (trk_mc3idx_branch != 0) {
				trk_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3idx_branch does not exist!\n");
				exit(1);
			}
			trk_mc3idx_isLoaded = true;
		}
		return trk_mc3idx_;
	}
	vector<int> &trk_mc3_motherid()
	{
		if (not trk_mc3_motherid_isLoaded) {
			if (trk_mc3_motherid_branch != 0) {
				trk_mc3_motherid_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3_motherid_branch does not exist!\n");
				exit(1);
			}
			trk_mc3_motherid_isLoaded = true;
		}
		return trk_mc3_motherid_;
	}
	vector<int> &trk_mc_id()
	{
		if (not trk_mc_id_isLoaded) {
			if (trk_mc_id_branch != 0) {
				trk_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc_id_branch does not exist!\n");
				exit(1);
			}
			trk_mc_id_isLoaded = true;
		}
		return trk_mc_id_;
	}
	vector<int> &trk_mcidx()
	{
		if (not trk_mcidx_isLoaded) {
			if (trk_mcidx_branch != 0) {
				trk_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch trk_mcidx_branch does not exist!\n");
				exit(1);
			}
			trk_mcidx_isLoaded = true;
		}
		return trk_mcidx_;
	}
	vector<int> &trk_mc_motherid()
	{
		if (not trk_mc_motherid_isLoaded) {
			if (trk_mc_motherid_branch != 0) {
				trk_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			trk_mc_motherid_isLoaded = true;
		}
		return trk_mc_motherid_;
	}
	vector<int> &els_conv_tkidx()
	{
		if (not els_conv_tkidx_isLoaded) {
			if (els_conv_tkidx_branch != 0) {
				els_conv_tkidx_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_tkidx_branch does not exist!\n");
				exit(1);
			}
			els_conv_tkidx_isLoaded = true;
		}
		return els_conv_tkidx_;
	}
	vector<int> &trks_conv_tkidx()
	{
		if (not trks_conv_tkidx_isLoaded) {
			if (trks_conv_tkidx_branch != 0) {
				trks_conv_tkidx_branch->GetEntry(index);
			} else { 
				printf("branch trks_conv_tkidx_branch does not exist!\n");
				exit(1);
			}
			trks_conv_tkidx_isLoaded = true;
		}
		return trks_conv_tkidx_;
	}
	vector<int> &els_closestMuon()
	{
		if (not els_closestMuon_isLoaded) {
			if (els_closestMuon_branch != 0) {
				els_closestMuon_branch->GetEntry(index);
			} else { 
				printf("branch els_closestMuon_branch does not exist!\n");
				exit(1);
			}
			els_closestMuon_isLoaded = true;
		}
		return els_closestMuon_;
	}
	vector<int> &els_trkidx()
	{
		if (not els_trkidx_isLoaded) {
			if (els_trkidx_branch != 0) {
				els_trkidx_branch->GetEntry(index);
			} else { 
				printf("branch els_trkidx_branch does not exist!\n");
				exit(1);
			}
			els_trkidx_isLoaded = true;
		}
		return els_trkidx_;
	}
	vector<int> &els_category()
	{
		if (not els_category_isLoaded) {
			if (els_category_branch != 0) {
				els_category_branch->GetEntry(index);
			} else { 
				printf("branch els_category_branch does not exist!\n");
				exit(1);
			}
			els_category_isLoaded = true;
		}
		return els_category_;
	}
	vector<int> &els_categoryold()
	{
		if (not els_categoryold_isLoaded) {
			if (els_categoryold_branch != 0) {
				els_categoryold_branch->GetEntry(index);
			} else { 
				printf("branch els_categoryold_branch does not exist!\n");
				exit(1);
			}
			els_categoryold_isLoaded = true;
		}
		return els_categoryold_;
	}
	vector<int> &els_charge()
	{
		if (not els_charge_isLoaded) {
			if (els_charge_branch != 0) {
				els_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_charge_branch does not exist!\n");
				exit(1);
			}
			els_charge_isLoaded = true;
		}
		return els_charge_;
	}
	vector<int> &els_class()
	{
		if (not els_class_isLoaded) {
			if (els_class_branch != 0) {
				els_class_branch->GetEntry(index);
			} else { 
				printf("branch els_class_branch does not exist!\n");
				exit(1);
			}
			els_class_isLoaded = true;
		}
		return els_class_;
	}
	vector<int> &els_layer1_det()
	{
		if (not els_layer1_det_isLoaded) {
			if (els_layer1_det_branch != 0) {
				els_layer1_det_branch->GetEntry(index);
			} else { 
				printf("branch els_layer1_det_branch does not exist!\n");
				exit(1);
			}
			els_layer1_det_isLoaded = true;
		}
		return els_layer1_det_;
	}
	vector<int> &els_layer1_layer()
	{
		if (not els_layer1_layer_isLoaded) {
			if (els_layer1_layer_branch != 0) {
				els_layer1_layer_branch->GetEntry(index);
			} else { 
				printf("branch els_layer1_layer_branch does not exist!\n");
				exit(1);
			}
			els_layer1_layer_isLoaded = true;
		}
		return els_layer1_layer_;
	}
	vector<int> &els_layer1_sizerphi()
	{
		if (not els_layer1_sizerphi_isLoaded) {
			if (els_layer1_sizerphi_branch != 0) {
				els_layer1_sizerphi_branch->GetEntry(index);
			} else { 
				printf("branch els_layer1_sizerphi_branch does not exist!\n");
				exit(1);
			}
			els_layer1_sizerphi_isLoaded = true;
		}
		return els_layer1_sizerphi_;
	}
	vector<int> &els_layer1_sizerz()
	{
		if (not els_layer1_sizerz_isLoaded) {
			if (els_layer1_sizerz_branch != 0) {
				els_layer1_sizerz_branch->GetEntry(index);
			} else { 
				printf("branch els_layer1_sizerz_branch does not exist!\n");
				exit(1);
			}
			els_layer1_sizerz_isLoaded = true;
		}
		return els_layer1_sizerz_;
	}
	vector<int> &els_looseId()
	{
		if (not els_looseId_isLoaded) {
			if (els_looseId_branch != 0) {
				els_looseId_branch->GetEntry(index);
			} else { 
				printf("branch els_looseId_branch does not exist!\n");
				exit(1);
			}
			els_looseId_isLoaded = true;
		}
		return els_looseId_;
	}
	vector<int> &els_lostHits()
	{
		if (not els_lostHits_isLoaded) {
			if (els_lostHits_branch != 0) {
				els_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch els_lostHits_branch does not exist!\n");
				exit(1);
			}
			els_lostHits_isLoaded = true;
		}
		return els_lostHits_;
	}
	vector<int> &els_lost_pixelhits()
	{
		if (not els_lost_pixelhits_isLoaded) {
			if (els_lost_pixelhits_branch != 0) {
				els_lost_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch els_lost_pixelhits_branch does not exist!\n");
				exit(1);
			}
			els_lost_pixelhits_isLoaded = true;
		}
		return els_lost_pixelhits_;
	}
	vector<int> &els_nSeed()
	{
		if (not els_nSeed_isLoaded) {
			if (els_nSeed_branch != 0) {
				els_nSeed_branch->GetEntry(index);
			} else { 
				printf("branch els_nSeed_branch does not exist!\n");
				exit(1);
			}
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	vector<int> &els_n_inner_layers()
	{
		if (not els_n_inner_layers_isLoaded) {
			if (els_n_inner_layers_branch != 0) {
				els_n_inner_layers_branch->GetEntry(index);
			} else { 
				printf("branch els_n_inner_layers_branch does not exist!\n");
				exit(1);
			}
			els_n_inner_layers_isLoaded = true;
		}
		return els_n_inner_layers_;
	}
	vector<int> &els_n_outer_layers()
	{
		if (not els_n_outer_layers_isLoaded) {
			if (els_n_outer_layers_branch != 0) {
				els_n_outer_layers_branch->GetEntry(index);
			} else { 
				printf("branch els_n_outer_layers_branch does not exist!\n");
				exit(1);
			}
			els_n_outer_layers_isLoaded = true;
		}
		return els_n_outer_layers_;
	}
	vector<int> &els_pass3looseId()
	{
		if (not els_pass3looseId_isLoaded) {
			if (els_pass3looseId_branch != 0) {
				els_pass3looseId_branch->GetEntry(index);
			} else { 
				printf("branch els_pass3looseId_branch does not exist!\n");
				exit(1);
			}
			els_pass3looseId_isLoaded = true;
		}
		return els_pass3looseId_;
	}
	vector<int> &els_pass3simpleId()
	{
		if (not els_pass3simpleId_isLoaded) {
			if (els_pass3simpleId_branch != 0) {
				els_pass3simpleId_branch->GetEntry(index);
			} else { 
				printf("branch els_pass3simpleId_branch does not exist!\n");
				exit(1);
			}
			els_pass3simpleId_isLoaded = true;
		}
		return els_pass3simpleId_;
	}
	vector<int> &els_pass3tightId()
	{
		if (not els_pass3tightId_isLoaded) {
			if (els_pass3tightId_branch != 0) {
				els_pass3tightId_branch->GetEntry(index);
			} else { 
				printf("branch els_pass3tightId_branch does not exist!\n");
				exit(1);
			}
			els_pass3tightId_isLoaded = true;
		}
		return els_pass3tightId_;
	}
	vector<int> &els_robustId()
	{
		if (not els_robustId_isLoaded) {
			if (els_robustId_branch != 0) {
				els_robustId_branch->GetEntry(index);
			} else { 
				printf("branch els_robustId_branch does not exist!\n");
				exit(1);
			}
			els_robustId_isLoaded = true;
		}
		return els_robustId_;
	}
	vector<int> &els_simpleIdPlus()
	{
		if (not els_simpleIdPlus_isLoaded) {
			if (els_simpleIdPlus_branch != 0) {
				els_simpleIdPlus_branch->GetEntry(index);
			} else { 
				printf("branch els_simpleIdPlus_branch does not exist!\n");
				exit(1);
			}
			els_simpleIdPlus_isLoaded = true;
		}
		return els_simpleIdPlus_;
	}
	vector<int> &els_tightId()
	{
		if (not els_tightId_isLoaded) {
			if (els_tightId_branch != 0) {
				els_tightId_branch->GetEntry(index);
			} else { 
				printf("branch els_tightId_branch does not exist!\n");
				exit(1);
			}
			els_tightId_isLoaded = true;
		}
		return els_tightId_;
	}
	vector<int> &els_tightId22XMaxMatteo()
	{
		if (not els_tightId22XMaxMatteo_isLoaded) {
			if (els_tightId22XMaxMatteo_branch != 0) {
				els_tightId22XMaxMatteo_branch->GetEntry(index);
			} else { 
				printf("branch els_tightId22XMaxMatteo_branch does not exist!\n");
				exit(1);
			}
			els_tightId22XMaxMatteo_isLoaded = true;
		}
		return els_tightId22XMaxMatteo_;
	}
	vector<int> &els_tightId22XMinMatteo()
	{
		if (not els_tightId22XMinMatteo_isLoaded) {
			if (els_tightId22XMinMatteo_branch != 0) {
				els_tightId22XMinMatteo_branch->GetEntry(index);
			} else { 
				printf("branch els_tightId22XMinMatteo_branch does not exist!\n");
				exit(1);
			}
			els_tightId22XMinMatteo_isLoaded = true;
		}
		return els_tightId22XMinMatteo_;
	}
	vector<int> &els_validHits()
	{
		if (not els_validHits_isLoaded) {
			if (els_validHits_branch != 0) {
				els_validHits_branch->GetEntry(index);
			} else { 
				printf("branch els_validHits_branch does not exist!\n");
				exit(1);
			}
			els_validHits_isLoaded = true;
		}
		return els_validHits_;
	}
	vector<int> &els_valid_pixelhits()
	{
		if (not els_valid_pixelhits_isLoaded) {
			if (els_valid_pixelhits_branch != 0) {
				els_valid_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch els_valid_pixelhits_branch does not exist!\n");
				exit(1);
			}
			els_valid_pixelhits_isLoaded = true;
		}
		return els_valid_pixelhits_;
	}
	vector<int> &genps_id()
	{
		if (not genps_id_isLoaded) {
			if (genps_id_branch != 0) {
				genps_id_branch->GetEntry(index);
			} else { 
				printf("branch genps_id_branch does not exist!\n");
				exit(1);
			}
			genps_id_isLoaded = true;
		}
		return genps_id_;
	}
	vector<int> &genps_id_mother()
	{
		if (not genps_id_mother_isLoaded) {
			if (genps_id_mother_branch != 0) {
				genps_id_mother_branch->GetEntry(index);
			} else { 
				printf("branch genps_id_mother_branch does not exist!\n");
				exit(1);
			}
			genps_id_mother_isLoaded = true;
		}
		return genps_id_mother_;
	}
	vector<int> &genps_lepdaughter_id()
	{
		if (not genps_lepdaughter_id_isLoaded) {
			if (genps_lepdaughter_id_branch != 0) {
				genps_lepdaughter_id_branch->GetEntry(index);
			} else { 
				printf("branch genps_lepdaughter_id_branch does not exist!\n");
				exit(1);
			}
			genps_lepdaughter_id_isLoaded = true;
		}
		return genps_lepdaughter_id_;
	}
	vector<int> &genps_lepdaughter_idx()
	{
		if (not genps_lepdaughter_idx_isLoaded) {
			if (genps_lepdaughter_idx_branch != 0) {
				genps_lepdaughter_idx_branch->GetEntry(index);
			} else { 
				printf("branch genps_lepdaughter_idx_branch does not exist!\n");
				exit(1);
			}
			genps_lepdaughter_idx_isLoaded = true;
		}
		return genps_lepdaughter_idx_;
	}
	vector<int> &genps_status()
	{
		if (not genps_status_isLoaded) {
			if (genps_status_branch != 0) {
				genps_status_branch->GetEntry(index);
			} else { 
				printf("branch genps_status_branch does not exist!\n");
				exit(1);
			}
			genps_status_isLoaded = true;
		}
		return genps_status_;
	}
	vector<int> &mus_charge()
	{
		if (not mus_charge_isLoaded) {
			if (mus_charge_branch != 0) {
				mus_charge_branch->GetEntry(index);
			} else { 
				printf("branch mus_charge_branch does not exist!\n");
				exit(1);
			}
			mus_charge_isLoaded = true;
		}
		return mus_charge_;
	}
	vector<int> &mus_gfit_validHits()
	{
		if (not mus_gfit_validHits_isLoaded) {
			if (mus_gfit_validHits_branch != 0) {
				mus_gfit_validHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_validHits_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_validHits_isLoaded = true;
		}
		return mus_gfit_validHits_;
	}
	vector<int> &mus_goodmask()
	{
		if (not mus_goodmask_isLoaded) {
			if (mus_goodmask_branch != 0) {
				mus_goodmask_branch->GetEntry(index);
			} else { 
				printf("branch mus_goodmask_branch does not exist!\n");
				exit(1);
			}
			mus_goodmask_isLoaded = true;
		}
		return mus_goodmask_;
	}
	vector<int> &mus_iso03_ntrk()
	{
		if (not mus_iso03_ntrk_isLoaded) {
			if (mus_iso03_ntrk_branch != 0) {
				mus_iso03_ntrk_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso03_ntrk_branch does not exist!\n");
				exit(1);
			}
			mus_iso03_ntrk_isLoaded = true;
		}
		return mus_iso03_ntrk_;
	}
	vector<int> &mus_iso05_ntrk()
	{
		if (not mus_iso05_ntrk_isLoaded) {
			if (mus_iso05_ntrk_branch != 0) {
				mus_iso05_ntrk_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso05_ntrk_branch does not exist!\n");
				exit(1);
			}
			mus_iso05_ntrk_isLoaded = true;
		}
		return mus_iso05_ntrk_;
	}
	vector<int> &mus_lostHits()
	{
		if (not mus_lostHits_isLoaded) {
			if (mus_lostHits_branch != 0) {
				mus_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_lostHits_branch does not exist!\n");
				exit(1);
			}
			mus_lostHits_isLoaded = true;
		}
		return mus_lostHits_;
	}
	vector<int> &mus_nmatches()
	{
		if (not mus_nmatches_isLoaded) {
			if (mus_nmatches_branch != 0) {
				mus_nmatches_branch->GetEntry(index);
			} else { 
				printf("branch mus_nmatches_branch does not exist!\n");
				exit(1);
			}
			mus_nmatches_isLoaded = true;
		}
		return mus_nmatches_;
	}
	vector<int> &mus_pid_TM2DCompatibilityLoose()
	{
		if (not mus_pid_TM2DCompatibilityLoose_isLoaded) {
			if (mus_pid_TM2DCompatibilityLoose_branch != 0) {
				mus_pid_TM2DCompatibilityLoose_branch->GetEntry(index);
			} else { 
				printf("branch mus_pid_TM2DCompatibilityLoose_branch does not exist!\n");
				exit(1);
			}
			mus_pid_TM2DCompatibilityLoose_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityLoose_;
	}
	vector<int> &mus_pid_TM2DCompatibilityTight()
	{
		if (not mus_pid_TM2DCompatibilityTight_isLoaded) {
			if (mus_pid_TM2DCompatibilityTight_branch != 0) {
				mus_pid_TM2DCompatibilityTight_branch->GetEntry(index);
			} else { 
				printf("branch mus_pid_TM2DCompatibilityTight_branch does not exist!\n");
				exit(1);
			}
			mus_pid_TM2DCompatibilityTight_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityTight_;
	}
	vector<int> &mus_pid_TMLastStationLoose()
	{
		if (not mus_pid_TMLastStationLoose_isLoaded) {
			if (mus_pid_TMLastStationLoose_branch != 0) {
				mus_pid_TMLastStationLoose_branch->GetEntry(index);
			} else { 
				printf("branch mus_pid_TMLastStationLoose_branch does not exist!\n");
				exit(1);
			}
			mus_pid_TMLastStationLoose_isLoaded = true;
		}
		return mus_pid_TMLastStationLoose_;
	}
	vector<int> &mus_pid_TMLastStationTight()
	{
		if (not mus_pid_TMLastStationTight_isLoaded) {
			if (mus_pid_TMLastStationTight_branch != 0) {
				mus_pid_TMLastStationTight_branch->GetEntry(index);
			} else { 
				printf("branch mus_pid_TMLastStationTight_branch does not exist!\n");
				exit(1);
			}
			mus_pid_TMLastStationTight_isLoaded = true;
		}
		return mus_pid_TMLastStationTight_;
	}
	vector<int> &mus_trk_charge()
	{
		if (not mus_trk_charge_isLoaded) {
			if (mus_trk_charge_branch != 0) {
				mus_trk_charge_branch->GetEntry(index);
			} else { 
				printf("branch mus_trk_charge_branch does not exist!\n");
				exit(1);
			}
			mus_trk_charge_isLoaded = true;
		}
		return mus_trk_charge_;
	}
	vector<int> &mus_trkrefkey()
	{
		if (not mus_trkrefkey_isLoaded) {
			if (mus_trkrefkey_branch != 0) {
				mus_trkrefkey_branch->GetEntry(index);
			} else { 
				printf("branch mus_trkrefkey_branch does not exist!\n");
				exit(1);
			}
			mus_trkrefkey_isLoaded = true;
		}
		return mus_trkrefkey_;
	}
	vector<int> &mus_type()
	{
		if (not mus_type_isLoaded) {
			if (mus_type_branch != 0) {
				mus_type_branch->GetEntry(index);
			} else { 
				printf("branch mus_type_branch does not exist!\n");
				exit(1);
			}
			mus_type_isLoaded = true;
		}
		return mus_type_;
	}
	vector<int> &mus_validHits()
	{
		if (not mus_validHits_isLoaded) {
			if (mus_validHits_branch != 0) {
				mus_validHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_validHits_branch does not exist!\n");
				exit(1);
			}
			mus_validHits_isLoaded = true;
		}
		return mus_validHits_;
	}
	vector<int> &trks_charge()
	{
		if (not trks_charge_isLoaded) {
			if (trks_charge_branch != 0) {
				trks_charge_branch->GetEntry(index);
			} else { 
				printf("branch trks_charge_branch does not exist!\n");
				exit(1);
			}
			trks_charge_isLoaded = true;
		}
		return trks_charge_;
	}
	vector<int> &trks_layer1_det()
	{
		if (not trks_layer1_det_isLoaded) {
			if (trks_layer1_det_branch != 0) {
				trks_layer1_det_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_lost_pixelhits_branch does not exist!\n");
				exit(1);
			}
			trks_lost_pixelhits_isLoaded = true;
		}
		return trks_lost_pixelhits_;
	}
	vector<int> &trks_n_inner_layers()
	{
		if (not trks_n_inner_layers_isLoaded) {
			if (trks_n_inner_layers_branch != 0) {
				trks_n_inner_layers_branch->GetEntry(index);
			} else { 
				printf("branch trks_n_inner_layers_branch does not exist!\n");
				exit(1);
			}
			trks_n_inner_layers_isLoaded = true;
		}
		return trks_n_inner_layers_;
	}
	vector<int> &trks_n_outer_layers()
	{
		if (not trks_n_outer_layers_isLoaded) {
			if (trks_n_outer_layers_branch != 0) {
				trks_n_outer_layers_branch->GetEntry(index);
			} else { 
				printf("branch trks_n_outer_layers_branch does not exist!\n");
				exit(1);
			}
			trks_n_outer_layers_isLoaded = true;
		}
		return trks_n_outer_layers_;
	}
	vector<int> &trks_validHits()
	{
		if (not trks_validHits_isLoaded) {
			if (trks_validHits_branch != 0) {
				trks_validHits_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_valid_pixelhits_branch does not exist!\n");
				exit(1);
			}
			trks_valid_pixelhits_isLoaded = true;
		}
		return trks_valid_pixelhits_;
	}
	vector<int> &trks_elsidx()
	{
		if (not trks_elsidx_isLoaded) {
			if (trks_elsidx_branch != 0) {
				trks_elsidx_branch->GetEntry(index);
			} else { 
				printf("branch trks_elsidx_branch does not exist!\n");
				exit(1);
			}
			trks_elsidx_isLoaded = true;
		}
		return trks_elsidx_;
	}
	vector<int> &trk_musidx()
	{
		if (not trk_musidx_isLoaded) {
			if (trk_musidx_branch != 0) {
				trk_musidx_branch->GetEntry(index);
			} else { 
				printf("branch trk_musidx_branch does not exist!\n");
				exit(1);
			}
			trk_musidx_isLoaded = true;
		}
		return trk_musidx_;
	}
	vector<int> &vtxs_hasRefittedTracks()
	{
		if (not vtxs_hasRefittedTracks_isLoaded) {
			if (vtxs_hasRefittedTracks_branch != 0) {
				vtxs_hasRefittedTracks_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_hasRefittedTracks_branch does not exist!\n");
				exit(1);
			}
			vtxs_hasRefittedTracks_isLoaded = true;
		}
		return vtxs_hasRefittedTracks_;
	}
	vector<int> &vtxs_isFake()
	{
		if (not vtxs_isFake_isLoaded) {
			if (vtxs_isFake_branch != 0) {
				vtxs_isFake_branch->GetEntry(index);
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
			} else { 
				printf("branch vtxs_tracksSize_branch does not exist!\n");
				exit(1);
			}
			vtxs_tracksSize_isLoaded = true;
		}
		return vtxs_tracksSize_;
	}
	vector<vector<int> > &trks_hit_type()
	{
		if (not trks_hit_type_isLoaded) {
			if (trks_hit_type_branch != 0) {
				trks_hit_type_branch->GetEntry(index);
			} else { 
				printf("branch trks_hit_type_branch does not exist!\n");
				exit(1);
			}
			trks_hit_type_isLoaded = true;
		}
		return trks_hit_type_;
	}
	unsigned int &evt_nels()
	{
		if (not evt_nels_isLoaded) {
			if (evt_nels_branch != 0) {
				evt_nels_branch->GetEntry(index);
			} else { 
				printf("branch evt_nels_branch does not exist!\n");
				exit(1);
			}
			evt_nels_isLoaded = true;
		}
		return evt_nels_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) {
				evt_event_branch->GetEntry(index);
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
			} else { 
				printf("branch evt_run_branch does not exist!\n");
				exit(1);
			}
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_nvtxs()
	{
		if (not evt_nvtxs_isLoaded) {
			if (evt_nvtxs_branch != 0) {
				evt_nvtxs_branch->GetEntry(index);
			} else { 
				printf("branch evt_nvtxs_branch does not exist!\n");
				exit(1);
			}
			evt_nvtxs_isLoaded = true;
		}
		return evt_nvtxs_;
	}
	bool passHLTTrigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = evt_HLT_trigNames().begin();
		vector<TString>::const_iterator end_it = evt_HLT_trigNames().end();
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
			return evt_HLT1() & bitmask;
		}
		if(trigIndx >= 32 && trigIndx <= 63) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 32); 
			return evt_HLT2() & bitmask;
		}
		if(trigIndx >= 64 && trigIndx <= 95) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 64); 
			return evt_HLT3() & bitmask;
		}
		if(trigIndx >= 96 && trigIndx <= 127) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 96); 
			return evt_HLT4() & bitmask;
		}
		if(trigIndx >= 128 && trigIndx <= 159) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 128); 
			return evt_HLT5() & bitmask;
		}
		if(trigIndx >= 160 && trigIndx <= 191) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 160); 
			return evt_HLT6() & bitmask;
		}
		if(trigIndx >= 192 && trigIndx <= 223) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 192); 
			return evt_HLT7() & bitmask;
		}
		if(trigIndx >= 224 && trigIndx <= 255) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 224); 
			return evt_HLT8() & bitmask;
		}
	return 0;
	}
	bool passL1Trigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = evt_L1_trigNames().begin();
		vector<TString>::const_iterator end_it = evt_L1_trigNames().end();
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
			return evt_L1_1() & bitmask;
		}
		if(trigIndx >= 32 && trigIndx <= 63) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 32); 
			return evt_L1_2() & bitmask;
		}
		if(trigIndx >= 64 && trigIndx <= 95) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 64); 
			return evt_L1_3() & bitmask;
		}
		if(trigIndx >= 96 && trigIndx <= 127) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 96); 
			return evt_L1_4() & bitmask;
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
	vector<TString> &evt_HLT_trigNames() { return cms2.evt_HLT_trigNames(); }
	vector<TString> &evt_L1_trigNames() { return cms2.evt_L1_trigNames(); }
	double &genps_pthat() { return cms2.genps_pthat(); }
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &vtxs_position() { return cms2.vtxs_position(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &evt_bsp4() { return cms2.evt_bsp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_motherp4() { return cms2.els_mc3_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_p4() { return cms2.els_mc3_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_motherp4() { return cms2.els_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_p4() { return cms2.els_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_gp_p4() { return cms2.jets_mc_gp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_p4() { return cms2.jets_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc3_motherp4() { return cms2.mus_mc3_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc3_p4() { return cms2.mus_mc3_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc_motherp4() { return cms2.mus_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc_p4() { return cms2.mus_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trk_mc3p4() { return cms2.trk_mc3p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trk_mcp4() { return cms2.trk_mcp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4() { return cms2.els_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4In() { return cms2.els_p4In(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4Out() { return cms2.els_p4Out(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_trk_p4() { return cms2.els_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_vertex_p4() { return cms2.els_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_lepdaughter_p4() { return cms2.genps_lepdaughter_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_p4() { return cms2.genps_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_prod_vtx() { return cms2.genps_prod_vtx(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_gfit_outerPos_p4() { return cms2.mus_gfit_outerPos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4() { return cms2.mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4() { return cms2.mus_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_vertex_p4() { return cms2.mus_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4() { return cms2.trks_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_vertex_p4() { return cms2.trks_vertex_p4(); }
	vector<double> &els_inner_positionx() { return cms2.els_inner_positionx(); }
	vector<double> &els_inner_positiony() { return cms2.els_inner_positiony(); }
	vector<double> &els_inner_positionz() { return cms2.els_inner_positionz(); }
	vector<double> &trks_inner_positionx() { return cms2.trks_inner_positionx(); }
	vector<double> &trks_inner_positiony() { return cms2.trks_inner_positiony(); }
	vector<double> &trks_inner_positionz() { return cms2.trks_inner_positionz(); }
	float &evt_bs_dxdz() { return cms2.evt_bs_dxdz(); }
	float &evt_bs_dxdzErr() { return cms2.evt_bs_dxdzErr(); }
	float &evt_bs_dydz() { return cms2.evt_bs_dydz(); }
	float &evt_bs_dydzErr() { return cms2.evt_bs_dydzErr(); }
	float &evt_bs_sigmaZ() { return cms2.evt_bs_sigmaZ(); }
	float &evt_bs_sigmaZErr() { return cms2.evt_bs_sigmaZErr(); }
	float &evt_bs_widthx() { return cms2.evt_bs_widthx(); }
	float &evt_bs_widthxErr() { return cms2.evt_bs_widthxErr(); }
	float &evt_bs_widthy() { return cms2.evt_bs_widthy(); }
	float &evt_bs_widthyErr() { return cms2.evt_bs_widthyErr(); }
	float &evt_bs_xErr() { return cms2.evt_bs_xErr(); }
	float &evt_bs_yErr() { return cms2.evt_bs_yErr(); }
	float &evt_bs_zErr() { return cms2.evt_bs_zErr(); }
	float &evt_bField() { return cms2.evt_bField(); }
	float &evt_kfactor() { return cms2.evt_kfactor(); }
	float &evt_weight() { return cms2.evt_weight(); }
	float &evt_xsec_excl() { return cms2.evt_xsec_excl(); }
	float &evt_xsec_incl() { return cms2.evt_xsec_incl(); }
	float &gen_met() { return cms2.gen_met(); }
	float &gen_metPhi() { return cms2.gen_metPhi(); }
	vector<float> &els_mc3dr() { return cms2.els_mc3dr(); }
	vector<float> &els_mcdr() { return cms2.els_mcdr(); }
	vector<float> &jets_mcdr() { return cms2.jets_mcdr(); }
	vector<float> &jets_mc_emEnergy() { return cms2.jets_mc_emEnergy(); }
	vector<float> &jets_mc_gpdr() { return cms2.jets_mc_gpdr(); }
	vector<float> &jets_mc_hadEnergy() { return cms2.jets_mc_hadEnergy(); }
	vector<float> &jets_mc_invEnergy() { return cms2.jets_mc_invEnergy(); }
	vector<float> &jets_mc_otherEnergy() { return cms2.jets_mc_otherEnergy(); }
	vector<float> &mus_mc3dr() { return cms2.mus_mc3dr(); }
	vector<float> &mus_mcdr() { return cms2.mus_mcdr(); }
	vector<float> &trk_mc3dr() { return cms2.trk_mc3dr(); }
	vector<float> &trk_mcdr() { return cms2.trk_mcdr(); }
	vector<float> &els_conv_dcot() { return cms2.els_conv_dcot(); }
	vector<float> &els_conv_dist() { return cms2.els_conv_dist(); }
	vector<float> &trks_conv_dcot() { return cms2.trks_conv_dcot(); }
	vector<float> &trks_conv_dist() { return cms2.trks_conv_dist(); }
	vector<float> &els_musdr() { return cms2.els_musdr(); }
	vector<float> &els_trkdr() { return cms2.els_trkdr(); }
	vector<float> &els_trkshFrac() { return cms2.els_trkshFrac(); }
	vector<float> &els_chi2() { return cms2.els_chi2(); }
	vector<float> &els_d0() { return cms2.els_d0(); }
	vector<float> &els_d0Err() { return cms2.els_d0Err(); }
	vector<float> &els_d0corr() { return cms2.els_d0corr(); }
	vector<float> &els_dEtaIn() { return cms2.els_dEtaIn(); }
	vector<float> &els_dEtaOut() { return cms2.els_dEtaOut(); }
	vector<float> &els_dPhiIn() { return cms2.els_dPhiIn(); }
	vector<float> &els_dPhiInPhiOut() { return cms2.els_dPhiInPhiOut(); }
	vector<float> &els_dPhiOut() { return cms2.els_dPhiOut(); }
	vector<float> &els_e1x5() { return cms2.els_e1x5(); }
	vector<float> &els_e2x5Max() { return cms2.els_e2x5Max(); }
	vector<float> &els_e3x3() { return cms2.els_e3x3(); }
	vector<float> &els_e5x5() { return cms2.els_e5x5(); }
	vector<float> &els_eMax() { return cms2.els_eMax(); }
	vector<float> &els_eOverPIn() { return cms2.els_eOverPIn(); }
	vector<float> &els_eSC() { return cms2.els_eSC(); }
	vector<float> &els_eSCPresh() { return cms2.els_eSCPresh(); }
	vector<float> &els_eSCRaw() { return cms2.els_eSCRaw(); }
	vector<float> &els_eSeed() { return cms2.els_eSeed(); }
	vector<float> &els_eSeedOverPIn() { return cms2.els_eSeedOverPIn(); }
	vector<float> &els_eSeedOverPOut() { return cms2.els_eSeedOverPOut(); }
	vector<float> &els_ecalIso03() { return cms2.els_ecalIso03(); }
	vector<float> &els_ecalIso04() { return cms2.els_ecalIso04(); }
	vector<float> &els_egamma_looseId() { return cms2.els_egamma_looseId(); }
	vector<float> &els_egamma_robustHighEnergy() { return cms2.els_egamma_robustHighEnergy(); }
	vector<float> &els_egamma_robustLooseId() { return cms2.els_egamma_robustLooseId(); }
	vector<float> &els_egamma_robustTightId() { return cms2.els_egamma_robustTightId(); }
	vector<float> &els_egamma_tightId() { return cms2.els_egamma_tightId(); }
	vector<float> &els_etaErr() { return cms2.els_etaErr(); }
	vector<float> &els_etaSC() { return cms2.els_etaSC(); }
	vector<float> &els_fBrem() { return cms2.els_fBrem(); }
	vector<float> &els_hOverE() { return cms2.els_hOverE(); }
	vector<float> &els_hcalIso03() { return cms2.els_hcalIso03(); }
	vector<float> &els_hcalIso04() { return cms2.els_hcalIso04(); }
	vector<float> &els_layer1_charge() { return cms2.els_layer1_charge(); }
	vector<float> &els_ndof() { return cms2.els_ndof(); }
	vector<float> &els_outerEta() { return cms2.els_outerEta(); }
	vector<float> &els_outerPhi() { return cms2.els_outerPhi(); }
	vector<float> &els_phiErr() { return cms2.els_phiErr(); }
	vector<float> &els_phiSC() { return cms2.els_phiSC(); }
	vector<float> &els_pin() { return cms2.els_pin(); }
	vector<float> &els_pout() { return cms2.els_pout(); }
	vector<float> &els_ptErr() { return cms2.els_ptErr(); }
	vector<float> &els_sigmaEtaEta() { return cms2.els_sigmaEtaEta(); }
	vector<float> &els_sigmaIEtaIEta() { return cms2.els_sigmaIEtaIEta(); }
	vector<float> &els_sigmaIPhiIPhi() { return cms2.els_sigmaIPhiIPhi(); }
	vector<float> &els_sigmaPhiPhi() { return cms2.els_sigmaPhiPhi(); }
	vector<float> &els_tkIso03() { return cms2.els_tkIso03(); }
	vector<float> &els_tkIso04() { return cms2.els_tkIso04(); }
	vector<float> &els_vertexphi() { return cms2.els_vertexphi(); }
	vector<float> &els_z0() { return cms2.els_z0(); }
	vector<float> &els_z0Err() { return cms2.els_z0Err(); }
	vector<float> &els_z0corr() { return cms2.els_z0corr(); }
	vector<float> &mus_caloCompatibility() { return cms2.mus_caloCompatibility(); }
	vector<float> &mus_chi2() { return cms2.mus_chi2(); }
	vector<float> &mus_d0() { return cms2.mus_d0(); }
	vector<float> &mus_d0Err() { return cms2.mus_d0Err(); }
	vector<float> &mus_d0corr() { return cms2.mus_d0corr(); }
	vector<float> &mus_e_em() { return cms2.mus_e_em(); }
	vector<float> &mus_e_emS9() { return cms2.mus_e_emS9(); }
	vector<float> &mus_e_had() { return cms2.mus_e_had(); }
	vector<float> &mus_e_hadS9() { return cms2.mus_e_hadS9(); }
	vector<float> &mus_e_ho() { return cms2.mus_e_ho(); }
	vector<float> &mus_e_hoS9() { return cms2.mus_e_hoS9(); }
	vector<float> &mus_etaErr() { return cms2.mus_etaErr(); }
	vector<float> &mus_gfit_chi2() { return cms2.mus_gfit_chi2(); }
	vector<float> &mus_gfit_ndof() { return cms2.mus_gfit_ndof(); }
	vector<float> &mus_iso() { return cms2.mus_iso(); }
	vector<float> &mus_iso03_emEt() { return cms2.mus_iso03_emEt(); }
	vector<float> &mus_iso03_hadEt() { return cms2.mus_iso03_hadEt(); }
	vector<float> &mus_iso03_hoEt() { return cms2.mus_iso03_hoEt(); }
	vector<float> &mus_iso03_sumPt() { return cms2.mus_iso03_sumPt(); }
	vector<float> &mus_iso05_emEt() { return cms2.mus_iso05_emEt(); }
	vector<float> &mus_iso05_hadEt() { return cms2.mus_iso05_hadEt(); }
	vector<float> &mus_iso05_hoEt() { return cms2.mus_iso05_hoEt(); }
	vector<float> &mus_iso05_sumPt() { return cms2.mus_iso05_sumPt(); }
	vector<float> &mus_ndof() { return cms2.mus_ndof(); }
	vector<float> &mus_outerEta() { return cms2.mus_outerEta(); }
	vector<float> &mus_outerPhi() { return cms2.mus_outerPhi(); }
	vector<float> &mus_phiErr() { return cms2.mus_phiErr(); }
	vector<float> &mus_ptErr() { return cms2.mus_ptErr(); }
	vector<float> &mus_qoverp() { return cms2.mus_qoverp(); }
	vector<float> &mus_qoverpError() { return cms2.mus_qoverpError(); }
	vector<float> &mus_vertexphi() { return cms2.mus_vertexphi(); }
	vector<float> &mus_z0() { return cms2.mus_z0(); }
	vector<float> &mus_z0Err() { return cms2.mus_z0Err(); }
	vector<float> &mus_z0corr() { return cms2.mus_z0corr(); }
	vector<float> &trks_chi2() { return cms2.trks_chi2(); }
	vector<float> &trks_d0() { return cms2.trks_d0(); }
	vector<float> &trks_d0Err() { return cms2.trks_d0Err(); }
	vector<float> &trks_d0corr() { return cms2.trks_d0corr(); }
	vector<float> &trks_d0corrPhi() { return cms2.trks_d0corrPhi(); }
	vector<float> &trks_etaErr() { return cms2.trks_etaErr(); }
	vector<float> &trks_layer1_charge() { return cms2.trks_layer1_charge(); }
	vector<float> &trks_ndof() { return cms2.trks_ndof(); }
	vector<float> &trks_outerEt() { return cms2.trks_outerEt(); }
	vector<float> &trks_outerEta() { return cms2.trks_outerEta(); }
	vector<float> &trks_outerPhi() { return cms2.trks_outerPhi(); }
	vector<float> &trks_phiErr() { return cms2.trks_phiErr(); }
	vector<float> &trks_ptErr() { return cms2.trks_ptErr(); }
	vector<float> &trks_tkIso() { return cms2.trks_tkIso(); }
	vector<float> &trks_vertexphi() { return cms2.trks_vertexphi(); }
	vector<float> &trks_z0() { return cms2.trks_z0(); }
	vector<float> &trks_z0Err() { return cms2.trks_z0Err(); }
	vector<float> &trks_z0corr() { return cms2.trks_z0corr(); }
	vector<float> &trks_elsdr() { return cms2.trks_elsdr(); }
	vector<float> &trks_elsshFrac() { return cms2.trks_elsshFrac(); }
	vector<float> &trk_musdr() { return cms2.trk_musdr(); }
	vector<float> &vtxs_chi2() { return cms2.vtxs_chi2(); }
	vector<float> &vtxs_ndof() { return cms2.vtxs_ndof(); }
	vector<float> &vtxs_xError() { return cms2.vtxs_xError(); }
	vector<float> &vtxs_yError() { return cms2.vtxs_yError(); }
	vector<float> &vtxs_zError() { return cms2.vtxs_zError(); }
	vector<vector<float> > &trks_residualX() { return cms2.trks_residualX(); }
	vector<vector<float> > &trks_residualY() { return cms2.trks_residualY(); }
	int &evt_HLT1() { return cms2.evt_HLT1(); }
	int &evt_HLT2() { return cms2.evt_HLT2(); }
	int &evt_HLT3() { return cms2.evt_HLT3(); }
	int &evt_HLT4() { return cms2.evt_HLT4(); }
	int &evt_HLT5() { return cms2.evt_HLT5(); }
	int &evt_HLT6() { return cms2.evt_HLT6(); }
	int &evt_HLT7() { return cms2.evt_HLT7(); }
	int &evt_HLT8() { return cms2.evt_HLT8(); }
	int &evt_L1_1() { return cms2.evt_L1_1(); }
	int &evt_L1_2() { return cms2.evt_L1_2(); }
	int &evt_L1_3() { return cms2.evt_L1_3(); }
	int &evt_L1_4() { return cms2.evt_L1_4(); }
	vector<int> &els_mc3_id() { return cms2.els_mc3_id(); }
	vector<int> &els_mc3idx() { return cms2.els_mc3idx(); }
	vector<int> &els_mc3_motherid() { return cms2.els_mc3_motherid(); }
	vector<int> &els_mc_id() { return cms2.els_mc_id(); }
	vector<int> &els_mcidx() { return cms2.els_mcidx(); }
	vector<int> &els_mc_motherid() { return cms2.els_mc_motherid(); }
	vector<int> &jets_mc_id() { return cms2.jets_mc_id(); }
	vector<int> &mus_mc3_id() { return cms2.mus_mc3_id(); }
	vector<int> &mus_mc3idx() { return cms2.mus_mc3idx(); }
	vector<int> &mus_mc3_motherid() { return cms2.mus_mc3_motherid(); }
	vector<int> &mus_mc_id() { return cms2.mus_mc_id(); }
	vector<int> &mus_mcidx() { return cms2.mus_mcidx(); }
	vector<int> &mus_mc_motherid() { return cms2.mus_mc_motherid(); }
	vector<int> &trk_mc3_id() { return cms2.trk_mc3_id(); }
	vector<int> &trk_mc3idx() { return cms2.trk_mc3idx(); }
	vector<int> &trk_mc3_motherid() { return cms2.trk_mc3_motherid(); }
	vector<int> &trk_mc_id() { return cms2.trk_mc_id(); }
	vector<int> &trk_mcidx() { return cms2.trk_mcidx(); }
	vector<int> &trk_mc_motherid() { return cms2.trk_mc_motherid(); }
	vector<int> &els_conv_tkidx() { return cms2.els_conv_tkidx(); }
	vector<int> &trks_conv_tkidx() { return cms2.trks_conv_tkidx(); }
	vector<int> &els_closestMuon() { return cms2.els_closestMuon(); }
	vector<int> &els_trkidx() { return cms2.els_trkidx(); }
	vector<int> &els_category() { return cms2.els_category(); }
	vector<int> &els_categoryold() { return cms2.els_categoryold(); }
	vector<int> &els_charge() { return cms2.els_charge(); }
	vector<int> &els_class() { return cms2.els_class(); }
	vector<int> &els_layer1_det() { return cms2.els_layer1_det(); }
	vector<int> &els_layer1_layer() { return cms2.els_layer1_layer(); }
	vector<int> &els_layer1_sizerphi() { return cms2.els_layer1_sizerphi(); }
	vector<int> &els_layer1_sizerz() { return cms2.els_layer1_sizerz(); }
	vector<int> &els_looseId() { return cms2.els_looseId(); }
	vector<int> &els_lostHits() { return cms2.els_lostHits(); }
	vector<int> &els_lost_pixelhits() { return cms2.els_lost_pixelhits(); }
	vector<int> &els_nSeed() { return cms2.els_nSeed(); }
	vector<int> &els_n_inner_layers() { return cms2.els_n_inner_layers(); }
	vector<int> &els_n_outer_layers() { return cms2.els_n_outer_layers(); }
	vector<int> &els_pass3looseId() { return cms2.els_pass3looseId(); }
	vector<int> &els_pass3simpleId() { return cms2.els_pass3simpleId(); }
	vector<int> &els_pass3tightId() { return cms2.els_pass3tightId(); }
	vector<int> &els_robustId() { return cms2.els_robustId(); }
	vector<int> &els_simpleIdPlus() { return cms2.els_simpleIdPlus(); }
	vector<int> &els_tightId() { return cms2.els_tightId(); }
	vector<int> &els_tightId22XMaxMatteo() { return cms2.els_tightId22XMaxMatteo(); }
	vector<int> &els_tightId22XMinMatteo() { return cms2.els_tightId22XMinMatteo(); }
	vector<int> &els_validHits() { return cms2.els_validHits(); }
	vector<int> &els_valid_pixelhits() { return cms2.els_valid_pixelhits(); }
	vector<int> &genps_id() { return cms2.genps_id(); }
	vector<int> &genps_id_mother() { return cms2.genps_id_mother(); }
	vector<int> &genps_lepdaughter_id() { return cms2.genps_lepdaughter_id(); }
	vector<int> &genps_lepdaughter_idx() { return cms2.genps_lepdaughter_idx(); }
	vector<int> &genps_status() { return cms2.genps_status(); }
	vector<int> &mus_charge() { return cms2.mus_charge(); }
	vector<int> &mus_gfit_validHits() { return cms2.mus_gfit_validHits(); }
	vector<int> &mus_goodmask() { return cms2.mus_goodmask(); }
	vector<int> &mus_iso03_ntrk() { return cms2.mus_iso03_ntrk(); }
	vector<int> &mus_iso05_ntrk() { return cms2.mus_iso05_ntrk(); }
	vector<int> &mus_lostHits() { return cms2.mus_lostHits(); }
	vector<int> &mus_nmatches() { return cms2.mus_nmatches(); }
	vector<int> &mus_pid_TM2DCompatibilityLoose() { return cms2.mus_pid_TM2DCompatibilityLoose(); }
	vector<int> &mus_pid_TM2DCompatibilityTight() { return cms2.mus_pid_TM2DCompatibilityTight(); }
	vector<int> &mus_pid_TMLastStationLoose() { return cms2.mus_pid_TMLastStationLoose(); }
	vector<int> &mus_pid_TMLastStationTight() { return cms2.mus_pid_TMLastStationTight(); }
	vector<int> &mus_trk_charge() { return cms2.mus_trk_charge(); }
	vector<int> &mus_trkrefkey() { return cms2.mus_trkrefkey(); }
	vector<int> &mus_type() { return cms2.mus_type(); }
	vector<int> &mus_validHits() { return cms2.mus_validHits(); }
	vector<int> &trks_charge() { return cms2.trks_charge(); }
	vector<int> &trks_layer1_det() { return cms2.trks_layer1_det(); }
	vector<int> &trks_layer1_layer() { return cms2.trks_layer1_layer(); }
	vector<int> &trks_layer1_sizerphi() { return cms2.trks_layer1_sizerphi(); }
	vector<int> &trks_layer1_sizerz() { return cms2.trks_layer1_sizerz(); }
	vector<int> &trks_lostHits() { return cms2.trks_lostHits(); }
	vector<int> &trks_lost_pixelhits() { return cms2.trks_lost_pixelhits(); }
	vector<int> &trks_n_inner_layers() { return cms2.trks_n_inner_layers(); }
	vector<int> &trks_n_outer_layers() { return cms2.trks_n_outer_layers(); }
	vector<int> &trks_validHits() { return cms2.trks_validHits(); }
	vector<int> &trks_valid_pixelhits() { return cms2.trks_valid_pixelhits(); }
	vector<int> &trks_elsidx() { return cms2.trks_elsidx(); }
	vector<int> &trk_musidx() { return cms2.trk_musidx(); }
	vector<int> &vtxs_hasRefittedTracks() { return cms2.vtxs_hasRefittedTracks(); }
	vector<int> &vtxs_isFake() { return cms2.vtxs_isFake(); }
	vector<int> &vtxs_isValid() { return cms2.vtxs_isValid(); }
	vector<int> &vtxs_tracksSize() { return cms2.vtxs_tracksSize(); }
	vector<vector<int> > &trks_hit_type() { return cms2.trks_hit_type(); }
	unsigned int &evt_nels() { return cms2.evt_nels(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	unsigned int &evt_nvtxs() { return cms2.evt_nvtxs(); }
	static bool passHLTTrigger(TString trigName) { return cms2.passHLTTrigger(trigName); }
	static bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }
}
#endif
