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
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

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
	vector<TString> hlt_trigNames_;
	TBranch *hlt_trigNames_branch;
	bool hlt_trigNames_isLoaded;
	vector<TString> l1_techtrigNames_;
	TBranch *l1_techtrigNames_branch;
	bool l1_techtrigNames_isLoaded;
	vector<TString> l1_trigNames_;
	TBranch *l1_trigNames_branch;
	bool l1_trigNames_isLoaded;
	vector<TString> evt_errCategory_;
	TBranch *evt_errCategory_branch;
	bool evt_errCategory_isLoaded;
	vector<TString> evt_errModule_;
	TBranch *evt_errModule_branch;
	bool evt_errModule_isLoaded;
	vector<TString> evt_errSeverity_;
	TBranch *evt_errSeverity_branch;
	bool evt_errSeverity_isLoaded;
	bool evt_eventHasHalo_;
	TBranch *evt_eventHasHalo_branch;
	bool evt_eventHasHalo_isLoaded;
	bool evt_hbheFilter_;
	TBranch *evt_hbheFilter_branch;
	bool evt_hbheFilter_isLoaded;
	vector<bool> mus_tightMatch_;
	TBranch *mus_tightMatch_branch;
	bool mus_tightMatch_isLoaded;
	vector<bool> mus_updatedSta_;
	TBranch *mus_updatedSta_branch;
	bool mus_updatedSta_isLoaded;
	vector<bool> photons_haspixelSeed_;
	TBranch *photons_haspixelSeed_branch;
	bool photons_haspixelSeed_isLoaded;
	vector<double> jets_closestElectron_DR_;
	TBranch *jets_closestElectron_DR_branch;
	bool jets_closestElectron_DR_isLoaded;
	vector<double> jets_closestMuon_DR_;
	TBranch *jets_closestMuon_DR_branch;
	bool jets_closestMuon_DR_isLoaded;
	float evt_bs_Xwidth_;
	TBranch *evt_bs_Xwidth_branch;
	bool evt_bs_Xwidth_isLoaded;
	float evt_bs_XwidthErr_;
	TBranch *evt_bs_XwidthErr_branch;
	bool evt_bs_XwidthErr_isLoaded;
	float evt_bs_Ywidth_;
	TBranch *evt_bs_Ywidth_branch;
	bool evt_bs_Ywidth_isLoaded;
	float evt_bs_YwidthErr_;
	TBranch *evt_bs_YwidthErr_branch;
	bool evt_bs_YwidthErr_isLoaded;
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
	float evt_bs_xErr_;
	TBranch *evt_bs_xErr_branch;
	bool evt_bs_xErr_isLoaded;
	float evt_bs_yErr_;
	TBranch *evt_bs_yErr_branch;
	bool evt_bs_yErr_isLoaded;
	float evt_bs_zErr_;
	TBranch *evt_bs_zErr_branch;
	bool evt_bs_zErr_isLoaded;
	float evthcal_dmetx_;
	TBranch *evthcal_dmetx_branch;
	bool evthcal_dmetx_isLoaded;
	float evthcal_dmety_;
	TBranch *evthcal_dmety_branch;
	bool evthcal_dmety_isLoaded;
	float evthcal_dsumet_;
	TBranch *evthcal_dsumet_branch;
	bool evthcal_dsumet_isLoaded;
	float evthf_dmetx_;
	TBranch *evthf_dmetx_branch;
	bool evthf_dmetx_isLoaded;
	float evthf_dmety_;
	TBranch *evthf_dmety_branch;
	bool evthf_dmety_isLoaded;
	float evthf_dsumet_;
	TBranch *evthf_dsumet_branch;
	bool evthf_dsumet_isLoaded;
	float evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float	evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float	evt_scale1fb_;
	TBranch *evt_scale1fb_branch;
	bool evt_scale1fb_isLoaded;
	float	evt_xsec_excl_;
	TBranch *evt_xsec_excl_branch;
	bool evt_xsec_excl_isLoaded;
	float	evt_xsec_incl_;
	TBranch *evt_xsec_incl_branch;
	bool evt_xsec_incl_isLoaded;
	float gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float gen_metPhi_;
	TBranch *gen_metPhi_branch;
	bool gen_metPhi_isLoaded;
	float genps_alphaQCD_;
	TBranch *genps_alphaQCD_branch;
	bool genps_alphaQCD_isLoaded;
	float genps_pthat_;
	TBranch *genps_pthat_branch;
	bool genps_pthat_isLoaded;
	float genps_qScale_;
	TBranch *genps_qScale_branch;
	bool genps_qScale_isLoaded;
	float genps_weight_;
	TBranch *genps_weight_branch;
	bool genps_weight_isLoaded;
	float gen_sumEt_;
	TBranch *gen_sumEt_branch;
	bool gen_sumEt_isLoaded;
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
	float evt_ecalendcapm_met_;
	TBranch *evt_ecalendcapm_met_branch;
	bool evt_ecalendcapm_met_isLoaded;
	float evt_ecalendcapm_metPhi_;
	TBranch *evt_ecalendcapm_metPhi_branch;
	bool evt_ecalendcapm_metPhi_isLoaded;
	float evt_ecalendcapp_met_;
	TBranch *evt_ecalendcapp_met_branch;
	bool evt_ecalendcapp_met_isLoaded;
	float evt_ecalendcapp_metPhi_;
	TBranch *evt_ecalendcapp_metPhi_branch;
	bool evt_ecalendcapp_metPhi_isLoaded;
	float evt_ecalmet_;
	TBranch *evt_ecalmet_branch;
	bool evt_ecalmet_isLoaded;
	float evt_ecalmetPhi_;
	TBranch *evt_ecalmetPhi_branch;
	bool evt_ecalmetPhi_isLoaded;
	float evt_endcapm_met_;
	TBranch *evt_endcapm_met_branch;
	bool evt_endcapm_met_isLoaded;
	float evt_endcapm_metPhi_;
	TBranch *evt_endcapm_metPhi_branch;
	bool evt_endcapm_metPhi_isLoaded;
	float evt_endcapp_met_;
	TBranch *evt_endcapp_met_branch;
	bool evt_endcapp_met_isLoaded;
	float evt_endcapp_metPhi_;
	TBranch *evt_endcapp_metPhi_branch;
	bool evt_endcapp_metPhi_isLoaded;
	float evt_hcalendcapm_met_;
	TBranch *evt_hcalendcapm_met_branch;
	bool evt_hcalendcapm_met_isLoaded;
	float evt_hcalendcapm_metPhi_;
	TBranch *evt_hcalendcapm_metPhi_branch;
	bool evt_hcalendcapm_metPhi_isLoaded;
	float evt_hcalendcapp_met_;
	TBranch *evt_hcalendcapp_met_branch;
	bool evt_hcalendcapp_met_isLoaded;
	float evt_hcalendcapp_metPhi_;
	TBranch *evt_hcalendcapp_metPhi_branch;
	bool evt_hcalendcapp_metPhi_isLoaded;
	float evt_hcalmet_;
	TBranch *evt_hcalmet_branch;
	bool evt_hcalmet_isLoaded;
	float evt_hcalmetPhi_;
	TBranch *evt_hcalmetPhi_branch;
	bool evt_hcalmetPhi_isLoaded;
	float evt_met_;
	TBranch *evt_met_branch;
	bool evt_met_isLoaded;
	float evt_metHO_;
	TBranch *evt_metHO_branch;
	bool evt_metHO_isLoaded;
	float evt_metHOPhi_;
	TBranch *evt_metHOPhi_branch;
	bool evt_metHOPhi_isLoaded;
	float evt_metHOSig_;
	TBranch *evt_metHOSig_branch;
	bool evt_metHOSig_isLoaded;
	float evt_metMuonCorr_;
	TBranch *evt_metMuonCorr_branch;
	bool evt_metMuonCorr_isLoaded;
	float evt_metMuonCorrPhi_;
	TBranch *evt_metMuonCorrPhi_branch;
	bool evt_metMuonCorrPhi_isLoaded;
	float evt_metMuonCorrSig_;
	TBranch *evt_metMuonCorrSig_branch;
	bool evt_metMuonCorrSig_isLoaded;
	float evt_metMuonJESCorr_;
	TBranch *evt_metMuonJESCorr_branch;
	bool evt_metMuonJESCorr_isLoaded;
	float evt_metMuonJESCorrPhi_;
	TBranch *evt_metMuonJESCorrPhi_branch;
	bool evt_metMuonJESCorrPhi_isLoaded;
	float evt_metMuonJESCorrSig_;
	TBranch *evt_metMuonJESCorrSig_branch;
	bool evt_metMuonJESCorrSig_isLoaded;
	float evt_metNoHF_;
	TBranch *evt_metNoHF_branch;
	bool evt_metNoHF_isLoaded;
	float evt_metNoHFHO_;
	TBranch *evt_metNoHFHO_branch;
	bool evt_metNoHFHO_isLoaded;
	float evt_metNoHFHOPhi_;
	TBranch *evt_metNoHFHOPhi_branch;
	bool evt_metNoHFHOPhi_isLoaded;
	float evt_metNoHFHOSig_;
	TBranch *evt_metNoHFHOSig_branch;
	bool evt_metNoHFHOSig_isLoaded;
	float evt_metNoHFPhi_;
	TBranch *evt_metNoHFPhi_branch;
	bool evt_metNoHFPhi_isLoaded;
	float evt_metNoHFSig_;
	TBranch *evt_metNoHFSig_branch;
	bool evt_metNoHFSig_isLoaded;
	float evt_metOpt_;
	TBranch *evt_metOpt_branch;
	bool evt_metOpt_isLoaded;
	float evt_metOptHO_;
	TBranch *evt_metOptHO_branch;
	bool evt_metOptHO_isLoaded;
	float evt_metOptHOPhi_;
	TBranch *evt_metOptHOPhi_branch;
	bool evt_metOptHOPhi_isLoaded;
	float evt_metOptHOSig_;
	TBranch *evt_metOptHOSig_branch;
	bool evt_metOptHOSig_isLoaded;
	float evt_metOptNoHF_;
	TBranch *evt_metOptNoHF_branch;
	bool evt_metOptNoHF_isLoaded;
	float evt_metOptNoHFHO_;
	TBranch *evt_metOptNoHFHO_branch;
	bool evt_metOptNoHFHO_isLoaded;
	float evt_metOptNoHFHOPhi_;
	TBranch *evt_metOptNoHFHOPhi_branch;
	bool evt_metOptNoHFHOPhi_isLoaded;
	float evt_metOptNoHFHOSig_;
	TBranch *evt_metOptNoHFHOSig_branch;
	bool evt_metOptNoHFHOSig_isLoaded;
	float evt_metOptNoHFPhi_;
	TBranch *evt_metOptNoHFPhi_branch;
	bool evt_metOptNoHFPhi_isLoaded;
	float evt_metOptNoHFSig_;
	TBranch *evt_metOptNoHFSig_branch;
	bool evt_metOptNoHFSig_isLoaded;
	float evt_metOptPhi_;
	TBranch *evt_metOptPhi_branch;
	bool evt_metOptPhi_isLoaded;
	float evt_metOptSig_;
	TBranch *evt_metOptSig_branch;
	bool evt_metOptSig_isLoaded;
	float evt_metPhi_;
	TBranch *evt_metPhi_branch;
	bool evt_metPhi_isLoaded;
	float evt_metSig_;
	TBranch *evt_metSig_branch;
	bool evt_metSig_isLoaded;
	float evt_sumet_;
	TBranch *evt_sumet_branch;
	bool evt_sumet_isLoaded;
	float evt_sumetHO_;
	TBranch *evt_sumetHO_branch;
	bool evt_sumetHO_isLoaded;
	float evt_sumetMuonCorr_;
	TBranch *evt_sumetMuonCorr_branch;
	bool evt_sumetMuonCorr_isLoaded;
	float evt_sumetNoHF_;
	TBranch *evt_sumetNoHF_branch;
	bool evt_sumetNoHF_isLoaded;
	float evt_sumetNoHFHO_;
	TBranch *evt_sumetNoHFHO_branch;
	bool evt_sumetNoHFHO_isLoaded;
	float evt_sumetOpt_;
	TBranch *evt_sumetOpt_branch;
	bool evt_sumetOpt_isLoaded;
	float evt_sumetOptHO_;
	TBranch *evt_sumetOptHO_branch;
	bool evt_sumetOptHO_isLoaded;
	float evt_sumetOptNoHF_;
	TBranch *evt_sumetOptNoHF_branch;
	bool evt_sumetOptNoHF_isLoaded;
	float evt_sumetOptNoHFHO_;
	TBranch *evt_sumetOptNoHFHO_branch;
	bool evt_sumetOptNoHFHO_isLoaded;
	float met_pat_metCor_;
	TBranch *met_pat_metCor_branch;
	bool met_pat_metCor_isLoaded;
	float met_pat_metPhiCor_;
	TBranch *met_pat_metPhiCor_branch;
	bool met_pat_metPhiCor_isLoaded;
	float met_pat_metPhiUncor_;
	TBranch *met_pat_metPhiUncor_branch;
	bool met_pat_metPhiUncor_isLoaded;
	float met_pat_metPhiUncorJES_;
	TBranch *met_pat_metPhiUncorJES_branch;
	bool met_pat_metPhiUncorJES_isLoaded;
	float met_pat_metPhiUncorMuon_;
	TBranch *met_pat_metPhiUncorMuon_branch;
	bool met_pat_metPhiUncorMuon_isLoaded;
	float met_pat_metUncor_;
	TBranch *met_pat_metUncor_branch;
	bool met_pat_metUncor_isLoaded;
	float met_pat_metUncorJES_;
	TBranch *met_pat_metUncorJES_branch;
	bool met_pat_metUncorJES_isLoaded;
	float met_pat_metUncorMuon_;
	TBranch *met_pat_metUncorMuon_branch;
	bool met_pat_metUncorMuon_isLoaded;
	float pdfinfo_scale_;
	TBranch *pdfinfo_scale_branch;
	bool pdfinfo_scale_isLoaded;
	float pdfinfo_x1_;
	TBranch *pdfinfo_x1_branch;
	bool pdfinfo_x1_isLoaded;
	float pdfinfo_x2_;
	TBranch *pdfinfo_x2_branch;
	bool pdfinfo_x2_isLoaded;
	float evt_pfmet_;
	TBranch *evt_pfmet_branch;
	bool evt_pfmet_isLoaded;
	float evt_pfmetPhi_;
	TBranch *evt_pfmetPhi_branch;
	bool evt_pfmetPhi_isLoaded;
	float evt_pfmetSig_;
	TBranch *evt_pfmetSig_branch;
	bool evt_pfmetSig_isLoaded;
	float evt_pfsumet_;
	TBranch *evt_pfsumet_branch;
	bool evt_pfsumet_isLoaded;
	float evt_tcmet_;
	TBranch *evt_tcmet_branch;
	bool evt_tcmet_isLoaded;
	float evt_tcmetPhi_;
	TBranch *evt_tcmetPhi_branch;
	bool evt_tcmetPhi_isLoaded;
	float evt_tcmetSig_;
	TBranch *evt_tcmetSig_branch;
	bool evt_tcmetSig_isLoaded;
	float evt_tcsumet_;
	TBranch *evt_tcsumet_branch;
	bool evt_tcsumet_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  evt_bsp4_;
	TBranch *evt_bsp4_branch;
	bool evt_bsp4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  l1_met_p4_;
	TBranch *l1_met_p4_branch;
	bool l1_met_p4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  l1_mht_p4_;
	TBranch *l1_mht_p4_branch;
	bool l1_mht_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_motherp4_;
	TBranch *els_mc_motherp4_branch;
	bool els_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_mc_p4_;
	TBranch *els_mc_p4_branch;
	bool els_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_mc_gp_p4_;
	TBranch *jets_mc_gp_p4_branch;
	bool jets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_mc_motherp4_;
	TBranch *jets_mc_motherp4_branch;
	bool jets_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_mc_p4_;
	TBranch *jets_mc_p4_branch;
	bool jets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_mc_motherp4_;
	TBranch *mus_mc_motherp4_branch;
	bool mus_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_mc_p4_;
	TBranch *mus_mc_p4_branch;
	bool mus_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_gp_p4_;
	TBranch *pfjets_mc_gp_p4_branch;
	bool pfjets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_motherp4_;
	TBranch *pfjets_mc_motherp4_branch;
	bool pfjets_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_mc_p4_;
	TBranch *pfjets_mc_p4_branch;
	bool pfjets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_mc_motherp4_;
	TBranch *photons_mc_motherp4_branch;
	bool photons_mc_motherp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_mc_p4_;
	TBranch *photons_mc_p4_branch;
	bool photons_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trk_mcp4_;
	TBranch *trk_mcp4_branch;
	bool trk_mcp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_conv_pos_p4_;
	TBranch *els_conv_pos_p4_branch;
	bool els_conv_pos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_inner_position_;
	TBranch *els_inner_position_branch;
	bool els_inner_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_outer_position_;
	TBranch *els_outer_position_branch;
	bool els_outer_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4In_;
	TBranch *els_p4In_branch;
	bool els_p4In_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_p4Out_;
	TBranch *els_p4Out_branch;
	bool els_p4Out_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_trk_p4_;
	TBranch *els_trk_p4_branch;
	bool els_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_vertex_p4_;
	TBranch *els_vertex_p4_branch;
	bool els_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genjets_p4_;
	TBranch *genjets_p4_branch;
	bool genjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > genps_prod_vtx_;
	TBranch *genps_prod_vtx_branch;
	bool genps_prod_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gsftrks_inner_position_;
	TBranch *gsftrks_inner_position_branch;
	bool gsftrks_inner_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gsftrks_outer_p4_;
	TBranch *gsftrks_outer_p4_branch;
	bool gsftrks_outer_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gsftrks_outer_position_;
	TBranch *gsftrks_outer_position_branch;
	bool gsftrks_outer_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gsftrks_p4_;
	TBranch *gsftrks_p4_branch;
	bool gsftrks_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gsftrks_vertex_p4_;
	TBranch *gsftrks_vertex_p4_branch;
	bool gsftrks_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_ll_p4_;
	TBranch *hyp_ll_p4_branch;
	bool hyp_ll_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_ll_trk_p4_;
	TBranch *hyp_ll_trk_p4_branch;
	bool hyp_ll_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_lt_p4_;
	TBranch *hyp_lt_p4_branch;
	bool hyp_lt_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_lt_trk_p4_;
	TBranch *hyp_lt_trk_p4_branch;
	bool hyp_lt_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_p4_;
	TBranch *hyp_p4_branch;
	bool hyp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_FVFit_p4_;
	TBranch *hyp_FVFit_p4_branch;
	bool hyp_FVFit_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_FVFit_v4_;
	TBranch *hyp_FVFit_v4_branch;
	bool hyp_FVFit_v4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_ll_mc_p4_;
	TBranch *hyp_ll_mc_p4_branch;
	bool hyp_ll_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > hyp_lt_mc_p4_;
	TBranch *hyp_lt_mc_p4_branch;
	bool hyp_lt_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_p4_;
	TBranch *jets_p4_branch;
	bool jets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_vertex_p4_;
	TBranch *jets_vertex_p4_branch;
	bool jets_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jpts_p4_;
	TBranch *jpts_p4_branch;
	bool jpts_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_ecalpos_p4_;
	TBranch *mus_ecalpos_p4_branch;
	bool mus_ecalpos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_fitdefault_p4_;
	TBranch *mus_fitdefault_p4_branch;
	bool mus_fitdefault_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_fitfirsthit_p4_;
	TBranch *mus_fitfirsthit_p4_branch;
	bool mus_fitfirsthit_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_fitpicky_p4_;
	TBranch *mus_fitpicky_p4_branch;
	bool mus_fitpicky_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_fittev_p4_;
	TBranch *mus_fittev_p4_branch;
	bool mus_fittev_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_gfit_outerPos_p4_;
	TBranch *mus_gfit_outerPos_p4_branch;
	bool mus_gfit_outerPos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_gfit_p4_;
	TBranch *mus_gfit_p4_branch;
	bool mus_gfit_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_gfit_vertex_p4_;
	TBranch *mus_gfit_vertex_p4_branch;
	bool mus_gfit_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_sta_p4_;
	TBranch *mus_sta_p4_branch;
	bool mus_sta_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_sta_vertex_p4_;
	TBranch *mus_sta_vertex_p4_branch;
	bool mus_sta_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_trk_p4_;
	TBranch *mus_trk_p4_branch;
	bool mus_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_vertex_p4_;
	TBranch *mus_vertex_p4_branch;
	bool mus_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_pat_genMotherP4_;
	TBranch *els_pat_genMotherP4_branch;
	bool els_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_pat_genP4_;
	TBranch *els_pat_genP4_branch;
	bool els_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > els_pat_p4_;
	TBranch *els_pat_p4_branch;
	bool els_pat_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_pat_genJet_p4_;
	TBranch *jets_pat_genJet_p4_branch;
	bool jets_pat_genJet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_pat_genPartonMother_p4_;
	TBranch *jets_pat_genPartonMother_p4_branch;
	bool jets_pat_genPartonMother_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_pat_genParton_p4_;
	TBranch *jets_pat_genParton_p4_branch;
	bool jets_pat_genParton_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_pat_jet_p4_;
	TBranch *jets_pat_jet_p4_branch;
	bool jets_pat_jet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jets_pat_jet_uncorp4_;
	TBranch *jets_pat_jet_uncorp4_branch;
	bool jets_pat_jet_uncorp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_pat_genMotherP4_;
	TBranch *mus_pat_genMotherP4_branch;
	bool mus_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_pat_genP4_;
	TBranch *mus_pat_genP4_branch;
	bool mus_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > mus_pat_p4_;
	TBranch *mus_pat_p4_branch;
	bool mus_pat_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfels_p4_;
	TBranch *pfels_p4_branch;
	bool pfels_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfels_posAtEcal_p4_;
	TBranch *pfels_posAtEcal_p4_branch;
	bool pfels_posAtEcal_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfjets_p4_;
	TBranch *pfjets_p4_branch;
	bool pfjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfmus_p4_;
	TBranch *pfmus_p4_branch;
	bool pfmus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > pfmus_posAtEcal_p4_;
	TBranch *pfmus_posAtEcal_p4_branch;
	bool pfmus_posAtEcal_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > photons_p4_;
	TBranch *photons_p4_branch;
	bool photons_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_p4_;
	TBranch *scs_p4_branch;
	bool scs_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_pos_p4_;
	TBranch *scs_pos_p4_branch;
	bool scs_pos_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > scs_vtx_p4_;
	TBranch *scs_vtx_p4_branch;
	bool scs_vtx_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > svs_flight_;
	TBranch *svs_flight_branch;
	bool svs_flight_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > svs_mc3_p4_;
	TBranch *svs_mc3_p4_branch;
	bool svs_mc3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > svs_p4_;
	TBranch *svs_p4_branch;
	bool svs_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > svs_position_;
	TBranch *svs_position_branch;
	bool svs_position_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > svs_refitp4_;
	TBranch *svs_refitp4_branch;
	bool svs_refitp4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > trkjets_p4_;
	TBranch *trkjets_p4_branch;
	bool trkjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > genps_lepdaughter_p4_;
	TBranch *genps_lepdaughter_p4_branch;
	bool genps_lepdaughter_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > hlt_trigObjs_p4_;
	TBranch *hlt_trigObjs_p4_branch;
	bool hlt_trigObjs_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > hyp_jets_p4_;
	TBranch *hyp_jets_p4_branch;
	bool hyp_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > hyp_other_jets_p4_;
	TBranch *hyp_other_jets_p4_branch;
	bool hyp_other_jets_p4_isLoaded;
	vector<float> jpts_combinedSecondaryVertexBJetTag_;
	TBranch *jpts_combinedSecondaryVertexBJetTag_branch;
	bool jpts_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> jpts_combinedSecondaryVertexMVABJetTag_;
	TBranch *jpts_combinedSecondaryVertexMVABJetTag_branch;
	bool jpts_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> jpts_jetBProbabilityBJetTag_;
	TBranch *jpts_jetBProbabilityBJetTag_branch;
	bool jpts_jetBProbabilityBJetTag_isLoaded;
	vector<float> jpts_jetProbabilityBJetTag_;
	TBranch *jpts_jetProbabilityBJetTag_branch;
	bool jpts_jetProbabilityBJetTag_isLoaded;
	vector<float> jpts_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *jpts_simpleSecondaryVertexHighEffBJetTag_branch;
	bool jpts_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> jpts_simpleSecondaryVertexHighPurBJetTags_;
	TBranch *jpts_simpleSecondaryVertexHighPurBJetTags_branch;
	bool jpts_simpleSecondaryVertexHighPurBJetTags_isLoaded;
	vector<float> jpts_softElectronByIP3dBJetTag_;
	TBranch *jpts_softElectronByIP3dBJetTag_branch;
	bool jpts_softElectronByIP3dBJetTag_isLoaded;
	vector<float> jpts_softElectronByPtBJetTag_;
	TBranch *jpts_softElectronByPtBJetTag_branch;
	bool jpts_softElectronByPtBJetTag_isLoaded;
	vector<float> jpts_softMuonBJetTag_;
	TBranch *jpts_softMuonBJetTag_branch;
	bool jpts_softMuonBJetTag_isLoaded;
	vector<float> jpts_softMuonByIP3dBJetTag_;
	TBranch *jpts_softMuonByIP3dBJetTag_branch;
	bool jpts_softMuonByIP3dBJetTag_isLoaded;
	vector<float> jpts_softMuonByPtBJetTag_;
	TBranch *jpts_softMuonByPtBJetTag_branch;
	bool jpts_softMuonByPtBJetTag_isLoaded;
	vector<float> jpts_trackCountingHighEffBJetTag_;
	TBranch *jpts_trackCountingHighEffBJetTag_branch;
	bool jpts_trackCountingHighEffBJetTag_isLoaded;
	vector<float> jpts_trackCountingHighPurBJetTag_;
	TBranch *jpts_trackCountingHighPurBJetTag_branch;
	bool jpts_trackCountingHighPurBJetTag_isLoaded;
	vector<float> jets_combinedSecondaryVertexBJetTag_;
	TBranch *jets_combinedSecondaryVertexBJetTag_branch;
	bool jets_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_combinedSecondaryVertexMVABJetTag_;
	TBranch *jets_combinedSecondaryVertexMVABJetTag_branch;
	bool jets_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> jets_jetBProbabilityBJetTag_;
	TBranch *jets_jetBProbabilityBJetTag_branch;
	bool jets_jetBProbabilityBJetTag_isLoaded;
	vector<float> jets_jetProbabilityBJetTag_;
	TBranch *jets_jetProbabilityBJetTag_branch;
	bool jets_jetProbabilityBJetTag_isLoaded;
	vector<float> jets_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *jets_simpleSecondaryVertexHighEffBJetTag_branch;
	bool jets_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> jets_simpleSecondaryVertexHighPurBJetTags_;
	TBranch *jets_simpleSecondaryVertexHighPurBJetTags_branch;
	bool jets_simpleSecondaryVertexHighPurBJetTags_isLoaded;
	vector<float> jets_softElectronByIP3dBJetTag_;
	TBranch *jets_softElectronByIP3dBJetTag_branch;
	bool jets_softElectronByIP3dBJetTag_isLoaded;
	vector<float> jets_softElectronByPtBJetTag_;
	TBranch *jets_softElectronByPtBJetTag_branch;
	bool jets_softElectronByPtBJetTag_isLoaded;
	vector<float> jets_softMuonBJetTag_;
	TBranch *jets_softMuonBJetTag_branch;
	bool jets_softMuonBJetTag_isLoaded;
	vector<float> jets_softMuonByIP3dBJetTag_;
	TBranch *jets_softMuonByIP3dBJetTag_branch;
	bool jets_softMuonByIP3dBJetTag_isLoaded;
	vector<float> jets_softMuonByPtBJetTag_;
	TBranch *jets_softMuonByPtBJetTag_branch;
	bool jets_softMuonByPtBJetTag_isLoaded;
	vector<float> jets_trackCountingHighEffBJetTag_;
	TBranch *jets_trackCountingHighEffBJetTag_branch;
	bool jets_trackCountingHighEffBJetTag_isLoaded;
	vector<float> jets_trackCountingHighPurBJetTag_;
	TBranch *jets_trackCountingHighPurBJetTag_branch;
	bool jets_trackCountingHighPurBJetTag_isLoaded;
	vector<float> pfjets_combinedSecondaryVertexBJetTag_;
	TBranch *pfjets_combinedSecondaryVertexBJetTag_branch;
	bool pfjets_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> pfjets_combinedSecondaryVertexMVABJetTag_;
	TBranch *pfjets_combinedSecondaryVertexMVABJetTag_branch;
	bool pfjets_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> pfjets_jetBProbabilityBJetTag_;
	TBranch *pfjets_jetBProbabilityBJetTag_branch;
	bool pfjets_jetBProbabilityBJetTag_isLoaded;
	vector<float> pfjets_jetProbabilityBJetTag_;
	TBranch *pfjets_jetProbabilityBJetTag_branch;
	bool pfjets_jetProbabilityBJetTag_isLoaded;
	vector<float> pfjets_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *pfjets_simpleSecondaryVertexHighEffBJetTag_branch;
	bool pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> pfjets_simpleSecondaryVertexHighPurBJetTags_;
	TBranch *pfjets_simpleSecondaryVertexHighPurBJetTags_branch;
	bool pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded;
	vector<float> pfjets_softElectronByIP3dBJetTag_;
	TBranch *pfjets_softElectronByIP3dBJetTag_branch;
	bool pfjets_softElectronByIP3dBJetTag_isLoaded;
	vector<float> pfjets_softElectronByPtBJetTag_;
	TBranch *pfjets_softElectronByPtBJetTag_branch;
	bool pfjets_softElectronByPtBJetTag_isLoaded;
	vector<float> pfjets_softMuonBJetTag_;
	TBranch *pfjets_softMuonBJetTag_branch;
	bool pfjets_softMuonBJetTag_isLoaded;
	vector<float> pfjets_softMuonByIP3dBJetTag_;
	TBranch *pfjets_softMuonByIP3dBJetTag_branch;
	bool pfjets_softMuonByIP3dBJetTag_isLoaded;
	vector<float> pfjets_softMuonByPtBJetTag_;
	TBranch *pfjets_softMuonByPtBJetTag_branch;
	bool pfjets_softMuonByPtBJetTag_isLoaded;
	vector<float> pfjets_trackCountingHighEffBJetTag_;
	TBranch *pfjets_trackCountingHighEffBJetTag_branch;
	bool pfjets_trackCountingHighEffBJetTag_isLoaded;
	vector<float> pfjets_trackCountingHighPurBJetTag_;
	TBranch *pfjets_trackCountingHighPurBJetTag_branch;
	bool pfjets_trackCountingHighPurBJetTag_isLoaded;
	vector<float> trkjets_combinedSecondaryVertexBJetTag_;
	TBranch *trkjets_combinedSecondaryVertexBJetTag_branch;
	bool trkjets_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> trkjets_combinedSecondaryVertexMVABJetTag_;
	TBranch *trkjets_combinedSecondaryVertexMVABJetTag_branch;
	bool trkjets_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> trkjets_jetBProbabilityBJetTag_;
	TBranch *trkjets_jetBProbabilityBJetTag_branch;
	bool trkjets_jetBProbabilityBJetTag_isLoaded;
	vector<float> trkjets_jetProbabilityBJetTag_;
	TBranch *trkjets_jetProbabilityBJetTag_branch;
	bool trkjets_jetProbabilityBJetTag_isLoaded;
	vector<float> trkjets_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *trkjets_simpleSecondaryVertexHighEffBJetTag_branch;
	bool trkjets_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> trkjets_simpleSecondaryVertexHighPurBJetTags_;
	TBranch *trkjets_simpleSecondaryVertexHighPurBJetTags_branch;
	bool trkjets_simpleSecondaryVertexHighPurBJetTags_isLoaded;
	vector<float> trkjets_softElectronByIP3dBJetTag_;
	TBranch *trkjets_softElectronByIP3dBJetTag_branch;
	bool trkjets_softElectronByIP3dBJetTag_isLoaded;
	vector<float> trkjets_softElectronByPtBJetTag_;
	TBranch *trkjets_softElectronByPtBJetTag_branch;
	bool trkjets_softElectronByPtBJetTag_isLoaded;
	vector<float> trkjets_softMuonBJetTag_;
	TBranch *trkjets_softMuonBJetTag_branch;
	bool trkjets_softMuonBJetTag_isLoaded;
	vector<float> trkjets_softMuonByIP3dBJetTag_;
	TBranch *trkjets_softMuonByIP3dBJetTag_branch;
	bool trkjets_softMuonByIP3dBJetTag_isLoaded;
	vector<float> trkjets_softMuonByPtBJetTag_;
	TBranch *trkjets_softMuonByPtBJetTag_branch;
	bool trkjets_softMuonByPtBJetTag_isLoaded;
	vector<float> trkjets_trackCountingHighEffBJetTag_;
	TBranch *trkjets_trackCountingHighEffBJetTag_branch;
	bool trkjets_trackCountingHighEffBJetTag_isLoaded;
	vector<float> trkjets_trackCountingHighPurBJetTag_;
	TBranch *trkjets_trackCountingHighPurBJetTag_branch;
	bool trkjets_trackCountingHighPurBJetTag_isLoaded;
	vector<float> evt_bs_covMatrix_;
	TBranch *evt_bs_covMatrix_branch;
	bool evt_bs_covMatrix_isLoaded;
	vector<float> els_mc3dr_;
	TBranch *els_mc3dr_branch;
	bool els_mc3dr_isLoaded;
	vector<float> els_mcdr_;
	TBranch *els_mcdr_branch;
	bool els_mcdr_isLoaded;
	vector<float> jets_mc3dr_;
	TBranch *jets_mc3dr_branch;
	bool jets_mc3dr_isLoaded;
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
	vector<float> pfjets_mc3dr_;
	TBranch *pfjets_mc3dr_branch;
	bool pfjets_mc3dr_isLoaded;
	vector<float> pfjets_mcdr_;
	TBranch *pfjets_mcdr_branch;
	bool pfjets_mcdr_isLoaded;
	vector<float> pfjets_mc_emEnergy_;
	TBranch *pfjets_mc_emEnergy_branch;
	bool pfjets_mc_emEnergy_isLoaded;
	vector<float> pfjets_mc_gpdr_;
	TBranch *pfjets_mc_gpdr_branch;
	bool pfjets_mc_gpdr_isLoaded;
	vector<float> pfjets_mc_hadEnergy_;
	TBranch *pfjets_mc_hadEnergy_branch;
	bool pfjets_mc_hadEnergy_isLoaded;
	vector<float> pfjets_mc_invEnergy_;
	TBranch *pfjets_mc_invEnergy_branch;
	bool pfjets_mc_invEnergy_isLoaded;
	vector<float> pfjets_mc_otherEnergy_;
	TBranch *pfjets_mc_otherEnergy_branch;
	bool pfjets_mc_otherEnergy_isLoaded;
	vector<float> photons_mc3dr_;
	TBranch *photons_mc3dr_branch;
	bool photons_mc3dr_isLoaded;
	vector<float> photons_mcdr_;
	TBranch *photons_mcdr_branch;
	bool photons_mcdr_isLoaded;
	vector<float> trk_mc3dr_;
	TBranch *trk_mc3dr_branch;
	bool trk_mc3dr_isLoaded;
	vector<float> trk_mcdr_;
	TBranch *trk_mcdr_branch;
	bool trk_mcdr_isLoaded;
	vector<float> trks_conv_dcot_;
	TBranch *trks_conv_dcot_branch;
	bool trks_conv_dcot_isLoaded;
	vector<float> trks_conv_dist_;
	TBranch *trks_conv_dist_branch;
	bool trks_conv_dist_isLoaded;
	vector<float> els_ecalJuraIso_;
	TBranch *els_ecalJuraIso_branch;
	bool els_ecalJuraIso_isLoaded;
	vector<float> els_ecalJuraTowerIso_;
	TBranch *els_ecalJuraTowerIso_branch;
	bool els_ecalJuraTowerIso_isLoaded;
	vector<float> els_hcalConeIso_;
	TBranch *els_hcalConeIso_branch;
	bool els_hcalConeIso_isLoaded;
	vector<float> els_tkJuraIso_;
	TBranch *els_tkJuraIso_branch;
	bool els_tkJuraIso_isLoaded;
	vector<float> els_jetdr_;
	TBranch *els_jetdr_branch;
	bool els_jetdr_isLoaded;
	vector<float> els_musdr_;
	TBranch *els_musdr_branch;
	bool els_musdr_isLoaded;
	vector<float> els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float> els_conv_dcot_;
	TBranch *els_conv_dcot_branch;
	bool els_conv_dcot_isLoaded;
	vector<float> els_conv_dist_;
	TBranch *els_conv_dist_branch;
	bool els_conv_dist_isLoaded;
	vector<float> els_conv_radius_;
	TBranch *els_conv_radius_branch;
	bool els_conv_radius_isLoaded;
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
	vector<float> els_deltaEtaEleClusterTrackAtCalo_;
	TBranch *els_deltaEtaEleClusterTrackAtCalo_branch;
	bool els_deltaEtaEleClusterTrackAtCalo_isLoaded;
	vector<float> els_deltaPhiEleClusterTrackAtCalo_;
	TBranch *els_deltaPhiEleClusterTrackAtCalo_branch;
	bool els_deltaPhiEleClusterTrackAtCalo_isLoaded;
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
	vector<float> els_eOverPOut_;
	TBranch *els_eOverPOut_branch;
	bool els_eOverPOut_isLoaded;
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
	vector<float> els_ecalEnergy_;
	TBranch *els_ecalEnergy_branch;
	bool els_ecalEnergy_isLoaded;
	vector<float> els_ecalEnergyError_;
	TBranch *els_ecalEnergyError_branch;
	bool els_ecalEnergyError_isLoaded;
	vector<float> els_ecalIso_;
	TBranch *els_ecalIso_branch;
	bool els_ecalIso_isLoaded;
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
	vector<float> els_electronMomentumError_;
	TBranch *els_electronMomentumError_branch;
	bool els_electronMomentumError_isLoaded;
	vector<float> els_etaErr_;
	TBranch *els_etaErr_branch;
	bool els_etaErr_isLoaded;
	vector<float> els_etaSC_;
	TBranch *els_etaSC_branch;
	bool els_etaSC_isLoaded;
	vector<float> els_fbrem_;
	TBranch *els_fbrem_branch;
	bool els_fbrem_isLoaded;
	vector<float> els_hOverE_;
	TBranch *els_hOverE_branch;
	bool els_hOverE_isLoaded;
	vector<float> els_hcalDepth1OverEcal_;
	TBranch *els_hcalDepth1OverEcal_branch;
	bool els_hcalDepth1OverEcal_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt_;
	TBranch *els_hcalDepth1TowerSumEt_branch;
	bool els_hcalDepth1TowerSumEt_isLoaded;
	vector<float> els_hcalDepth1TowerSumEt04_;
	TBranch *els_hcalDepth1TowerSumEt04_branch;
	bool els_hcalDepth1TowerSumEt04_isLoaded;
	vector<float> els_hcalDepth2OverEcal_;
	TBranch *els_hcalDepth2OverEcal_branch;
	bool els_hcalDepth2OverEcal_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt_;
	TBranch *els_hcalDepth2TowerSumEt_branch;
	bool els_hcalDepth2TowerSumEt_isLoaded;
	vector<float> els_hcalDepth2TowerSumEt04_;
	TBranch *els_hcalDepth2TowerSumEt04_branch;
	bool els_hcalDepth2TowerSumEt04_isLoaded;
	vector<float> els_hcalIso_;
	TBranch *els_hcalIso_branch;
	bool els_hcalIso_isLoaded;
	vector<float> els_hcalIso04_;
	TBranch *els_hcalIso04_branch;
	bool els_hcalIso04_isLoaded;
	vector<float> els_layer1_charge_;
	TBranch *els_layer1_charge_branch;
	bool els_layer1_charge_isLoaded;
	vector<float> els_mva_;
	TBranch *els_mva_branch;
	bool els_mva_isLoaded;
	vector<float> els_ndof_;
	TBranch *els_ndof_branch;
	bool els_ndof_isLoaded;
	vector<float> els_phiErr_;
	TBranch *els_phiErr_branch;
	bool els_phiErr_isLoaded;
	vector<float> els_phiSC_;
	TBranch *els_phiSC_branch;
	bool els_phiSC_isLoaded;
	vector<float> els_ptErr_;
	TBranch *els_ptErr_branch;
	bool els_ptErr_isLoaded;
	vector<float> els_sigmaEtaEta_;
	TBranch *els_sigmaEtaEta_branch;
	bool els_sigmaEtaEta_isLoaded;
	vector<float> els_sigmaIEtaIEta_;
	TBranch *els_sigmaIEtaIEta_branch;
	bool els_sigmaIEtaIEta_isLoaded;
	vector<float> els_sigmaIEtaIEtaSC_;
	TBranch *els_sigmaIEtaIEtaSC_branch;
	bool els_sigmaIEtaIEtaSC_isLoaded;
	vector<float> els_sigmaIPhiIPhi_;
	TBranch *els_sigmaIPhiIPhi_branch;
	bool els_sigmaIPhiIPhi_isLoaded;
	vector<float> els_sigmaIPhiIPhiSC_;
	TBranch *els_sigmaIPhiIPhiSC_branch;
	bool els_sigmaIPhiIPhiSC_isLoaded;
	vector<float> els_sigmaPhiPhi_;
	TBranch *els_sigmaPhiPhi_branch;
	bool els_sigmaPhiPhi_isLoaded;
	vector<float> els_tkIso_;
	TBranch *els_tkIso_branch;
	bool els_tkIso_isLoaded;
	vector<float> els_tkIso04_;
	TBranch *els_tkIso04_branch;
	bool els_tkIso04_isLoaded;
	vector<float> els_trackMomentumError_;
	TBranch *els_trackMomentumError_branch;
	bool els_trackMomentumError_isLoaded;
	vector<float> els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float> els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
	vector<float> els_z0_;
	TBranch *els_z0_branch;
	bool els_z0_isLoaded;
	vector<float> els_z0Err_;
	TBranch *els_z0Err_branch;
	bool els_z0Err_isLoaded;
	vector<float> els_z0corr_;
	TBranch *els_z0corr_branch;
	bool els_z0corr_isLoaded;
	vector<float> gsftrks_chi2_;
	TBranch *gsftrks_chi2_branch;
	bool gsftrks_chi2_isLoaded;
	vector<float> gsftrks_d0_;
	TBranch *gsftrks_d0_branch;
	bool gsftrks_d0_isLoaded;
	vector<float> gsftrks_d0Err_;
	TBranch *gsftrks_d0Err_branch;
	bool gsftrks_d0Err_isLoaded;
	vector<float> gsftrks_d0corr_;
	TBranch *gsftrks_d0corr_branch;
	bool gsftrks_d0corr_isLoaded;
	vector<float> gsftrks_d0corrPhi_;
	TBranch *gsftrks_d0corrPhi_branch;
	bool gsftrks_d0corrPhi_isLoaded;
	vector<float> gsftrks_d0phiCov_;
	TBranch *gsftrks_d0phiCov_branch;
	bool gsftrks_d0phiCov_isLoaded;
	vector<float> gsftrks_etaErr_;
	TBranch *gsftrks_etaErr_branch;
	bool gsftrks_etaErr_isLoaded;
	vector<float> gsftrks_layer1_charge_;
	TBranch *gsftrks_layer1_charge_branch;
	bool gsftrks_layer1_charge_isLoaded;
	vector<float> gsftrks_ndof_;
	TBranch *gsftrks_ndof_branch;
	bool gsftrks_ndof_isLoaded;
	vector<float> gsftrks_phiErr_;
	TBranch *gsftrks_phiErr_branch;
	bool gsftrks_phiErr_isLoaded;
	vector<float> gsftrks_ptErr_;
	TBranch *gsftrks_ptErr_branch;
	bool gsftrks_ptErr_isLoaded;
	vector<float> gsftrks_z0_;
	TBranch *gsftrks_z0_branch;
	bool gsftrks_z0_isLoaded;
	vector<float> gsftrks_z0Err_;
	TBranch *gsftrks_z0Err_branch;
	bool gsftrks_z0Err_isLoaded;
	vector<float> gsftrks_z0corr_;
	TBranch *gsftrks_z0corr_branch;
	bool gsftrks_z0corr_isLoaded;
	vector<float> hyp_Ht_;
	TBranch *hyp_Ht_branch;
	bool hyp_Ht_isLoaded;
	vector<float> hyp_dPhi_nJet_metMuonJESCorr_;
	TBranch *hyp_dPhi_nJet_metMuonJESCorr_branch;
	bool hyp_dPhi_nJet_metMuonJESCorr_isLoaded;
	vector<float> hyp_dPhi_nJet_muCorrMet_;
	TBranch *hyp_dPhi_nJet_muCorrMet_branch;
	bool hyp_dPhi_nJet_muCorrMet_isLoaded;
	vector<float> hyp_dPhi_nJet_tcMet_;
	TBranch *hyp_dPhi_nJet_tcMet_branch;
	bool hyp_dPhi_nJet_tcMet_isLoaded;
	vector<float> hyp_dPhi_nJet_unCorrMet_;
	TBranch *hyp_dPhi_nJet_unCorrMet_branch;
	bool hyp_dPhi_nJet_unCorrMet_isLoaded;
	vector<float> hyp_ll_chi2_;
	TBranch *hyp_ll_chi2_branch;
	bool hyp_ll_chi2_isLoaded;
	vector<float> hyp_ll_d0_;
	TBranch *hyp_ll_d0_branch;
	bool hyp_ll_d0_isLoaded;
	vector<float> hyp_ll_d0Err_;
	TBranch *hyp_ll_d0Err_branch;
	bool hyp_ll_d0Err_isLoaded;
	vector<float> hyp_ll_d0corr_;
	TBranch *hyp_ll_d0corr_branch;
	bool hyp_ll_d0corr_isLoaded;
	vector<float> hyp_ll_dPhi_metMuonJESCorr_;
	TBranch *hyp_ll_dPhi_metMuonJESCorr_branch;
	bool hyp_ll_dPhi_metMuonJESCorr_isLoaded;
	vector<float> hyp_ll_dPhi_muCorrMet_;
	TBranch *hyp_ll_dPhi_muCorrMet_branch;
	bool hyp_ll_dPhi_muCorrMet_isLoaded;
	vector<float> hyp_ll_dPhi_tcMet_;
	TBranch *hyp_ll_dPhi_tcMet_branch;
	bool hyp_ll_dPhi_tcMet_isLoaded;
	vector<float> hyp_ll_dPhi_unCorrMet_;
	TBranch *hyp_ll_dPhi_unCorrMet_branch;
	bool hyp_ll_dPhi_unCorrMet_isLoaded;
	vector<float> hyp_ll_etaErr_;
	TBranch *hyp_ll_etaErr_branch;
	bool hyp_ll_etaErr_isLoaded;
	vector<float> hyp_ll_ndof_;
	TBranch *hyp_ll_ndof_branch;
	bool hyp_ll_ndof_isLoaded;
	vector<float> hyp_ll_phiErr_;
	TBranch *hyp_ll_phiErr_branch;
	bool hyp_ll_phiErr_isLoaded;
	vector<float> hyp_ll_ptErr_;
	TBranch *hyp_ll_ptErr_branch;
	bool hyp_ll_ptErr_isLoaded;
	vector<float> hyp_ll_z0_;
	TBranch *hyp_ll_z0_branch;
	bool hyp_ll_z0_isLoaded;
	vector<float> hyp_ll_z0Err_;
	TBranch *hyp_ll_z0Err_branch;
	bool hyp_ll_z0Err_isLoaded;
	vector<float> hyp_ll_z0corr_;
	TBranch *hyp_ll_z0corr_branch;
	bool hyp_ll_z0corr_isLoaded;
	vector<float> hyp_lt_chi2_;
	TBranch *hyp_lt_chi2_branch;
	bool hyp_lt_chi2_isLoaded;
	vector<float> hyp_lt_d0_;
	TBranch *hyp_lt_d0_branch;
	bool hyp_lt_d0_isLoaded;
	vector<float> hyp_lt_d0Err_;
	TBranch *hyp_lt_d0Err_branch;
	bool hyp_lt_d0Err_isLoaded;
	vector<float> hyp_lt_d0corr_;
	TBranch *hyp_lt_d0corr_branch;
	bool hyp_lt_d0corr_isLoaded;
	vector<float> hyp_lt_dPhi_metMuonJESCorr_;
	TBranch *hyp_lt_dPhi_metMuonJESCorr_branch;
	bool hyp_lt_dPhi_metMuonJESCorr_isLoaded;
	vector<float> hyp_lt_dPhi_muCorrMet_;
	TBranch *hyp_lt_dPhi_muCorrMet_branch;
	bool hyp_lt_dPhi_muCorrMet_isLoaded;
	vector<float> hyp_lt_dPhi_tcMet_;
	TBranch *hyp_lt_dPhi_tcMet_branch;
	bool hyp_lt_dPhi_tcMet_isLoaded;
	vector<float> hyp_lt_dPhi_unCorrMet_;
	TBranch *hyp_lt_dPhi_unCorrMet_branch;
	bool hyp_lt_dPhi_unCorrMet_isLoaded;
	vector<float> hyp_lt_etaErr_;
	TBranch *hyp_lt_etaErr_branch;
	bool hyp_lt_etaErr_isLoaded;
	vector<float> hyp_lt_ndof_;
	TBranch *hyp_lt_ndof_branch;
	bool hyp_lt_ndof_isLoaded;
	vector<float> hyp_lt_phiErr_;
	TBranch *hyp_lt_phiErr_branch;
	bool hyp_lt_phiErr_isLoaded;
	vector<float> hyp_lt_ptErr_;
	TBranch *hyp_lt_ptErr_branch;
	bool hyp_lt_ptErr_isLoaded;
	vector<float> hyp_lt_z0_;
	TBranch *hyp_lt_z0_branch;
	bool hyp_lt_z0_isLoaded;
	vector<float> hyp_lt_z0Err_;
	TBranch *hyp_lt_z0Err_branch;
	bool hyp_lt_z0Err_isLoaded;
	vector<float> hyp_lt_z0corr_;
	TBranch *hyp_lt_z0corr_branch;
	bool hyp_lt_z0corr_isLoaded;
	vector<float> hyp_mt2_metMuonJESCorr_;
	TBranch *hyp_mt2_metMuonJESCorr_branch;
	bool hyp_mt2_metMuonJESCorr_isLoaded;
	vector<float> hyp_mt2_muCorrMet_;
	TBranch *hyp_mt2_muCorrMet_branch;
	bool hyp_mt2_muCorrMet_isLoaded;
	vector<float> hyp_mt2_tcMet_;
	TBranch *hyp_mt2_tcMet_branch;
	bool hyp_mt2_tcMet_isLoaded;
	vector<float> hyp_sumJetPt_;
	TBranch *hyp_sumJetPt_branch;
	bool hyp_sumJetPt_isLoaded;
	vector<float> hyp_FVFit_chi2ndf_;
	TBranch *hyp_FVFit_chi2ndf_branch;
	bool hyp_FVFit_chi2ndf_isLoaded;
	vector<float> hyp_FVFit_prob_;
	TBranch *hyp_FVFit_prob_branch;
	bool hyp_FVFit_prob_isLoaded;
	vector<float> hyp_FVFit_v4cxx_;
	TBranch *hyp_FVFit_v4cxx_branch;
	bool hyp_FVFit_v4cxx_isLoaded;
	vector<float> hyp_FVFit_v4cxy_;
	TBranch *hyp_FVFit_v4cxy_branch;
	bool hyp_FVFit_v4cxy_isLoaded;
	vector<float> hyp_FVFit_v4cyy_;
	TBranch *hyp_FVFit_v4cyy_branch;
	bool hyp_FVFit_v4cyy_isLoaded;
	vector<float> hyp_FVFit_v4czx_;
	TBranch *hyp_FVFit_v4czx_branch;
	bool hyp_FVFit_v4czx_isLoaded;
	vector<float> hyp_FVFit_v4czy_;
	TBranch *hyp_FVFit_v4czy_branch;
	bool hyp_FVFit_v4czy_isLoaded;
	vector<float> hyp_FVFit_v4czz_;
	TBranch *hyp_FVFit_v4czz_branch;
	bool hyp_FVFit_v4czz_isLoaded;
	vector<float> hyp_ll_ecaliso_;
	TBranch *hyp_ll_ecaliso_branch;
	bool hyp_ll_ecaliso_isLoaded;
	vector<float> hyp_ll_trkiso_;
	TBranch *hyp_ll_trkiso_branch;
	bool hyp_ll_trkiso_isLoaded;
	vector<float> hyp_lt_ecaliso_;
	TBranch *hyp_lt_ecaliso_branch;
	bool hyp_lt_ecaliso_isLoaded;
	vector<float> hyp_lt_trkiso_;
	TBranch *hyp_lt_trkiso_branch;
	bool hyp_lt_trkiso_isLoaded;
	vector<float> jets_approximatefHPD_;
	TBranch *jets_approximatefHPD_branch;
	bool jets_approximatefHPD_isLoaded;
	vector<float> jets_approximatefRBX_;
	TBranch *jets_approximatefRBX_branch;
	bool jets_approximatefRBX_isLoaded;
	vector<float> jets_cor_;
	TBranch *jets_cor_branch;
	bool jets_cor_isLoaded;
	vector<float> jets_emFrac_;
	TBranch *jets_emFrac_branch;
	bool jets_emFrac_isLoaded;
	vector<float> jets_fHPD_;
	TBranch *jets_fHPD_branch;
	bool jets_fHPD_isLoaded;
	vector<float> jets_fRBX_;
	TBranch *jets_fRBX_branch;
	bool jets_fRBX_isLoaded;
	vector<float> jets_fSubDetector1_;
	TBranch *jets_fSubDetector1_branch;
	bool jets_fSubDetector1_isLoaded;
	vector<float> jets_fSubDetector2_;
	TBranch *jets_fSubDetector2_branch;
	bool jets_fSubDetector2_isLoaded;
	vector<float> jets_fSubDetector3_;
	TBranch *jets_fSubDetector3_branch;
	bool jets_fSubDetector3_isLoaded;
	vector<float> jets_fSubDetector4_;
	TBranch *jets_fSubDetector4_branch;
	bool jets_fSubDetector4_isLoaded;
	vector<float> jets_hitsInN90_;
	TBranch *jets_hitsInN90_branch;
	bool jets_hitsInN90_isLoaded;
	vector<float> jets_n90Hits_;
	TBranch *jets_n90Hits_branch;
	bool jets_n90Hits_isLoaded;
	vector<float> jets_nECALTowers_;
	TBranch *jets_nECALTowers_branch;
	bool jets_nECALTowers_isLoaded;
	vector<float> jets_nHCALTowers_;
	TBranch *jets_nHCALTowers_branch;
	bool jets_nHCALTowers_isLoaded;
	vector<float> jets_restrictedEMF_;
	TBranch *jets_restrictedEMF_branch;
	bool jets_restrictedEMF_isLoaded;
	vector<float> jpts_cor_;
	TBranch *jpts_cor_branch;
	bool jpts_cor_isLoaded;
	vector<float> jpts_emFrac_;
	TBranch *jpts_emFrac_branch;
	bool jpts_emFrac_isLoaded;
	vector<float> mus_met_deltax_;
	TBranch *mus_met_deltax_branch;
	bool mus_met_deltax_isLoaded;
	vector<float> mus_met_deltay_;
	TBranch *mus_met_deltay_branch;
	bool mus_met_deltay_isLoaded;
	vector<float> mus_eledr_;
	TBranch *mus_eledr_branch;
	bool mus_eledr_isLoaded;
	vector<float> mus_jetdr_;
	TBranch *mus_jetdr_branch;
	bool mus_jetdr_isLoaded;
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
	vector<float> mus_gfit_d0_;
	TBranch *mus_gfit_d0_branch;
	bool mus_gfit_d0_isLoaded;
	vector<float> mus_gfit_d0Err_;
	TBranch *mus_gfit_d0Err_branch;
	bool mus_gfit_d0Err_isLoaded;
	vector<float> mus_gfit_d0corr_;
	TBranch *mus_gfit_d0corr_branch;
	bool mus_gfit_d0corr_isLoaded;
	vector<float> mus_gfit_ndof_;
	TBranch *mus_gfit_ndof_branch;
	bool mus_gfit_ndof_isLoaded;
	vector<float> mus_gfit_qoverp_;
	TBranch *mus_gfit_qoverp_branch;
	bool mus_gfit_qoverp_isLoaded;
	vector<float> mus_gfit_qoverpError_;
	TBranch *mus_gfit_qoverpError_branch;
	bool mus_gfit_qoverpError_isLoaded;
	vector<float> mus_gfit_z0_;
	TBranch *mus_gfit_z0_branch;
	bool mus_gfit_z0_isLoaded;
	vector<float> mus_gfit_z0Err_;
	TBranch *mus_gfit_z0Err_branch;
	bool mus_gfit_z0Err_isLoaded;
	vector<float> mus_gfit_z0corr_;
	TBranch *mus_gfit_z0corr_branch;
	bool mus_gfit_z0corr_isLoaded;
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
	vector<float> mus_iso_ecalvetoDep_;
	TBranch *mus_iso_ecalvetoDep_branch;
	bool mus_iso_ecalvetoDep_isLoaded;
	vector<float> mus_iso_hcalvetoDep_;
	TBranch *mus_iso_hcalvetoDep_branch;
	bool mus_iso_hcalvetoDep_isLoaded;
	vector<float> mus_iso_hovetoDep_;
	TBranch *mus_iso_hovetoDep_branch;
	bool mus_iso_hovetoDep_isLoaded;
	vector<float> mus_iso_trckvetoDep_;
	TBranch *mus_iso_trckvetoDep_branch;
	bool mus_iso_trckvetoDep_isLoaded;
	vector<float> mus_ndof_;
	TBranch *mus_ndof_branch;
	bool mus_ndof_isLoaded;
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
	vector<float> mus_sta_chi2_;
	TBranch *mus_sta_chi2_branch;
	bool mus_sta_chi2_isLoaded;
	vector<float> mus_sta_d0_;
	TBranch *mus_sta_d0_branch;
	bool mus_sta_d0_isLoaded;
	vector<float> mus_sta_d0Err_;
	TBranch *mus_sta_d0Err_branch;
	bool mus_sta_d0Err_isLoaded;
	vector<float> mus_sta_d0corr_;
	TBranch *mus_sta_d0corr_branch;
	bool mus_sta_d0corr_isLoaded;
	vector<float> mus_sta_ndof_;
	TBranch *mus_sta_ndof_branch;
	bool mus_sta_ndof_isLoaded;
	vector<float> mus_sta_qoverp_;
	TBranch *mus_sta_qoverp_branch;
	bool mus_sta_qoverp_isLoaded;
	vector<float> mus_sta_qoverpError_;
	TBranch *mus_sta_qoverpError_branch;
	bool mus_sta_qoverpError_isLoaded;
	vector<float> mus_sta_z0_;
	TBranch *mus_sta_z0_branch;
	bool mus_sta_z0_isLoaded;
	vector<float> mus_sta_z0Err_;
	TBranch *mus_sta_z0Err_branch;
	bool mus_sta_z0Err_isLoaded;
	vector<float> mus_sta_z0corr_;
	TBranch *mus_sta_z0corr_branch;
	bool mus_sta_z0corr_isLoaded;
	vector<float> mus_timeAtIpInOut_;
	TBranch *mus_timeAtIpInOut_branch;
	bool mus_timeAtIpInOut_isLoaded;
	vector<float> mus_timeAtIpInOutErr_;
	TBranch *mus_timeAtIpInOutErr_branch;
	bool mus_timeAtIpInOutErr_isLoaded;
	vector<float> mus_timeAtIpOutIn_;
	TBranch *mus_timeAtIpOutIn_branch;
	bool mus_timeAtIpOutIn_isLoaded;
	vector<float> mus_timeAtIpOutInErr_;
	TBranch *mus_timeAtIpOutInErr_branch;
	bool mus_timeAtIpOutInErr_isLoaded;
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
	vector<float> els_pat_caloIso_;
	TBranch *els_pat_caloIso_branch;
	bool els_pat_caloIso_isLoaded;
	vector<float> els_pat_ecalIso_;
	TBranch *els_pat_ecalIso_branch;
	bool els_pat_ecalIso_isLoaded;
	vector<float> els_pat_hcalIso_;
	TBranch *els_pat_hcalIso_branch;
	bool els_pat_hcalIso_isLoaded;
	vector<float> els_pat_looseId_;
	TBranch *els_pat_looseId_branch;
	bool els_pat_looseId_isLoaded;
	vector<float> els_pat_robustHighEnergy_;
	TBranch *els_pat_robustHighEnergy_branch;
	bool els_pat_robustHighEnergy_isLoaded;
	vector<float> els_pat_robustLooseId_;
	TBranch *els_pat_robustLooseId_branch;
	bool els_pat_robustLooseId_isLoaded;
	vector<float> els_pat_robustTightId_;
	TBranch *els_pat_robustTightId_branch;
	bool els_pat_robustTightId_isLoaded;
	vector<float> els_pat_scE1x5_;
	TBranch *els_pat_scE1x5_branch;
	bool els_pat_scE1x5_isLoaded;
	vector<float> els_pat_scE2x5Max_;
	TBranch *els_pat_scE2x5Max_branch;
	bool els_pat_scE2x5Max_isLoaded;
	vector<float> els_pat_scE5x5_;
	TBranch *els_pat_scE5x5_branch;
	bool els_pat_scE5x5_isLoaded;
	vector<float> els_pat_sigmaEtaEta_;
	TBranch *els_pat_sigmaEtaEta_branch;
	bool els_pat_sigmaEtaEta_isLoaded;
	vector<float> els_pat_sigmaIEtaIEta_;
	TBranch *els_pat_sigmaIEtaIEta_branch;
	bool els_pat_sigmaIEtaIEta_isLoaded;
	vector<float> els_pat_tightId_;
	TBranch *els_pat_tightId_branch;
	bool els_pat_tightId_isLoaded;
	vector<float> els_pat_trackIso_;
	TBranch *els_pat_trackIso_branch;
	bool els_pat_trackIso_isLoaded;
	vector<float> jets_pat_combinedSecondaryVertexBJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexBJetTag_branch;
	bool jets_pat_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_pat_combinedSecondaryVertexMVABJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexMVABJetTag_branch;
	bool jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> jets_pat_coneIsolationTauJetTag_;
	TBranch *jets_pat_coneIsolationTauJetTag_branch;
	bool jets_pat_coneIsolationTauJetTag_isLoaded;
	vector<float> jets_pat_impactParameterMVABJetTag_;
	TBranch *jets_pat_impactParameterMVABJetTag_branch;
	bool jets_pat_impactParameterMVABJetTag_isLoaded;
	vector<float> jets_pat_jetBProbabilityBJetTag_;
	TBranch *jets_pat_jetBProbabilityBJetTag_branch;
	bool jets_pat_jetBProbabilityBJetTag_isLoaded;
	vector<float> jets_pat_jetCharge_;
	TBranch *jets_pat_jetCharge_branch;
	bool jets_pat_jetCharge_isLoaded;
	vector<float> jets_pat_jetProbabilityBJetTag_;
	TBranch *jets_pat_jetProbabilityBJetTag_branch;
	bool jets_pat_jetProbabilityBJetTag_isLoaded;
	vector<float> jets_pat_noCorrF_;
	TBranch *jets_pat_noCorrF_branch;
	bool jets_pat_noCorrF_isLoaded;
	vector<float> jets_pat_simpleSecondaryVertexHighEffBJetTag_;
	TBranch *jets_pat_simpleSecondaryVertexHighEffBJetTag_branch;
	bool jets_pat_simpleSecondaryVertexHighEffBJetTag_isLoaded;
	vector<float> jets_pat_simpleSecondaryVertexHighPurBJetTag_;
	TBranch *jets_pat_simpleSecondaryVertexHighPurBJetTag_branch;
	bool jets_pat_simpleSecondaryVertexHighPurBJetTag_isLoaded;
	vector<float> jets_pat_softElectronByIP3dBJetTag_;
	TBranch *jets_pat_softElectronByIP3dBJetTag_branch;
	bool jets_pat_softElectronByIP3dBJetTag_isLoaded;
	vector<float> jets_pat_softElectronByPtBJetTag_;
	TBranch *jets_pat_softElectronByPtBJetTag_branch;
	bool jets_pat_softElectronByPtBJetTag_isLoaded;
	vector<float> jets_pat_softMuonBJetTag_;
	TBranch *jets_pat_softMuonBJetTag_branch;
	bool jets_pat_softMuonBJetTag_isLoaded;
	vector<float> jets_pat_softMuonByIP3dBJetTag_;
	TBranch *jets_pat_softMuonByIP3dBJetTag_branch;
	bool jets_pat_softMuonByIP3dBJetTag_isLoaded;
	vector<float> jets_pat_softMuonByPtBJetTag_;
	TBranch *jets_pat_softMuonByPtBJetTag_branch;
	bool jets_pat_softMuonByPtBJetTag_isLoaded;
	vector<float> jets_pat_trackCountingHighEffBJetTag_;
	TBranch *jets_pat_trackCountingHighEffBJetTag_branch;
	bool jets_pat_trackCountingHighEffBJetTag_isLoaded;
	vector<float> jets_pat_trackCountingHighPurBJetTag_;
	TBranch *jets_pat_trackCountingHighPurBJetTag_branch;
	bool jets_pat_trackCountingHighPurBJetTag_isLoaded;
	vector<float> mus_pat_caloIso_;
	TBranch *mus_pat_caloIso_branch;
	bool mus_pat_caloIso_isLoaded;
	vector<float> mus_pat_calovetoDep_;
	TBranch *mus_pat_calovetoDep_branch;
	bool mus_pat_calovetoDep_isLoaded;
	vector<float> mus_pat_ecalIso_;
	TBranch *mus_pat_ecalIso_branch;
	bool mus_pat_ecalIso_isLoaded;
	vector<float> mus_pat_ecalvetoDep_;
	TBranch *mus_pat_ecalvetoDep_branch;
	bool mus_pat_ecalvetoDep_isLoaded;
	vector<float> mus_pat_hcalIso_;
	TBranch *mus_pat_hcalIso_branch;
	bool mus_pat_hcalIso_isLoaded;
	vector<float> mus_pat_hcalvetoDep_;
	TBranch *mus_pat_hcalvetoDep_branch;
	bool mus_pat_hcalvetoDep_isLoaded;
	vector<float> mus_pat_trackIso_;
	TBranch *mus_pat_trackIso_branch;
	bool mus_pat_trackIso_isLoaded;
	vector<float> mus_pat_trckvetoDep_;
	TBranch *mus_pat_trckvetoDep_branch;
	bool mus_pat_trckvetoDep_isLoaded;
	vector<float> pfels_deltaP_;
	TBranch *pfels_deltaP_branch;
	bool pfels_deltaP_isLoaded;
	vector<float> pfels_ecalE_;
	TBranch *pfels_ecalE_branch;
	bool pfels_ecalE_isLoaded;
	vector<float> pfels_hcalE_;
	TBranch *pfels_hcalE_branch;
	bool pfels_hcalE_isLoaded;
	vector<float> pfels_isoChargedHadrons_;
	TBranch *pfels_isoChargedHadrons_branch;
	bool pfels_isoChargedHadrons_isLoaded;
	vector<float> pfels_isoNeutralHadrons_;
	TBranch *pfels_isoNeutralHadrons_branch;
	bool pfels_isoNeutralHadrons_isLoaded;
	vector<float> pfels_isoPhotons_;
	TBranch *pfels_isoPhotons_branch;
	bool pfels_isoPhotons_isLoaded;
	vector<float> pfels_mva_emu_;
	TBranch *pfels_mva_emu_branch;
	bool pfels_mva_emu_isLoaded;
	vector<float> pfels_mva_epi_;
	TBranch *pfels_mva_epi_branch;
	bool pfels_mva_epi_isLoaded;
	vector<float> pfels_mva_nothing_gamma_;
	TBranch *pfels_mva_nothing_gamma_branch;
	bool pfels_mva_nothing_gamma_isLoaded;
	vector<float> pfels_mva_nothing_nh_;
	TBranch *pfels_mva_nothing_nh_branch;
	bool pfels_mva_nothing_nh_isLoaded;
	vector<float> pfels_mva_pimu_;
	TBranch *pfels_mva_pimu_branch;
	bool pfels_mva_pimu_isLoaded;
	vector<float> pfels_pS1E_;
	TBranch *pfels_pS1E_branch;
	bool pfels_pS1E_isLoaded;
	vector<float> pfels_pS2E_;
	TBranch *pfels_pS2E_branch;
	bool pfels_pS2E_isLoaded;
	vector<float> pfels_rawEcalE_;
	TBranch *pfels_rawEcalE_branch;
	bool pfels_rawEcalE_isLoaded;
	vector<float> pfels_rawHcalE_;
	TBranch *pfels_rawHcalE_branch;
	bool pfels_rawHcalE_isLoaded;
	vector<float> pfjets_chargedEmE_;
	TBranch *pfjets_chargedEmE_branch;
	bool pfjets_chargedEmE_isLoaded;
	vector<float> pfjets_chargedHadronE_;
	TBranch *pfjets_chargedHadronE_branch;
	bool pfjets_chargedHadronE_isLoaded;
	vector<float> pfjets_cor_;
	TBranch *pfjets_cor_branch;
	bool pfjets_cor_isLoaded;
	vector<float> pfjets_neutralEmE_;
	TBranch *pfjets_neutralEmE_branch;
	bool pfjets_neutralEmE_isLoaded;
	vector<float> pfjets_neutralHadronE_;
	TBranch *pfjets_neutralHadronE_branch;
	bool pfjets_neutralHadronE_isLoaded;
	vector<float> pfmus_deltaP_;
	TBranch *pfmus_deltaP_branch;
	bool pfmus_deltaP_isLoaded;
	vector<float> pfmus_ecalE_;
	TBranch *pfmus_ecalE_branch;
	bool pfmus_ecalE_isLoaded;
	vector<float> pfmus_hcalE_;
	TBranch *pfmus_hcalE_branch;
	bool pfmus_hcalE_isLoaded;
	vector<float> pfmus_isoChargedHadrons_;
	TBranch *pfmus_isoChargedHadrons_branch;
	bool pfmus_isoChargedHadrons_isLoaded;
	vector<float> pfmus_isoNeutralHadrons_;
	TBranch *pfmus_isoNeutralHadrons_branch;
	bool pfmus_isoNeutralHadrons_isLoaded;
	vector<float> pfmus_isoPhotons_;
	TBranch *pfmus_isoPhotons_branch;
	bool pfmus_isoPhotons_isLoaded;
	vector<float> pfmus_mva_emu_;
	TBranch *pfmus_mva_emu_branch;
	bool pfmus_mva_emu_isLoaded;
	vector<float> pfmus_mva_epi_;
	TBranch *pfmus_mva_epi_branch;
	bool pfmus_mva_epi_isLoaded;
	vector<float> pfmus_mva_nothing_gamma_;
	TBranch *pfmus_mva_nothing_gamma_branch;
	bool pfmus_mva_nothing_gamma_isLoaded;
	vector<float> pfmus_mva_nothing_nh_;
	TBranch *pfmus_mva_nothing_nh_branch;
	bool pfmus_mva_nothing_nh_isLoaded;
	vector<float> pfmus_mva_pimu_;
	TBranch *pfmus_mva_pimu_branch;
	bool pfmus_mva_pimu_isLoaded;
	vector<float> pfmus_pS1E_;
	TBranch *pfmus_pS1E_branch;
	bool pfmus_pS1E_isLoaded;
	vector<float> pfmus_pS2E_;
	TBranch *pfmus_pS2E_branch;
	bool pfmus_pS2E_isLoaded;
	vector<float> pfmus_rawEcalE_;
	TBranch *pfmus_rawEcalE_branch;
	bool pfmus_rawEcalE_isLoaded;
	vector<float> pfmus_rawHcalE_;
	TBranch *pfmus_rawHcalE_branch;
	bool pfmus_rawHcalE_isLoaded;
	vector<float> photons_e1x5_;
	TBranch *photons_e1x5_branch;
	bool photons_e1x5_isLoaded;
	vector<float> photons_e2x5Max_;
	TBranch *photons_e2x5Max_branch;
	bool photons_e2x5Max_isLoaded;
	vector<float> photons_e3x3_;
	TBranch *photons_e3x3_branch;
	bool photons_e3x3_isLoaded;
	vector<float> photons_e5x5_;
	TBranch *photons_e5x5_branch;
	bool photons_e5x5_isLoaded;
	vector<float> photons_ecalIso03_;
	TBranch *photons_ecalIso03_branch;
	bool photons_ecalIso03_isLoaded;
	vector<float> photons_ecalIso04_;
	TBranch *photons_ecalIso04_branch;
	bool photons_ecalIso04_isLoaded;
	vector<float> photons_hOverE_;
	TBranch *photons_hOverE_branch;
	bool photons_hOverE_isLoaded;
	vector<float> photons_hcalIso03_;
	TBranch *photons_hcalIso03_branch;
	bool photons_hcalIso03_isLoaded;
	vector<float> photons_hcalIso04_;
	TBranch *photons_hcalIso04_branch;
	bool photons_hcalIso04_isLoaded;
	vector<float> photons_ntkIsoHollow03_;
	TBranch *photons_ntkIsoHollow03_branch;
	bool photons_ntkIsoHollow03_isLoaded;
	vector<float> photons_ntkIsoHollow04_;
	TBranch *photons_ntkIsoHollow04_branch;
	bool photons_ntkIsoHollow04_isLoaded;
	vector<float> photons_ntkIsoSolid03_;
	TBranch *photons_ntkIsoSolid03_branch;
	bool photons_ntkIsoSolid03_isLoaded;
	vector<float> photons_ntkIsoSolid04_;
	TBranch *photons_ntkIsoSolid04_branch;
	bool photons_ntkIsoSolid04_isLoaded;
	vector<float> photons_sigmaEtaEta_;
	TBranch *photons_sigmaEtaEta_branch;
	bool photons_sigmaEtaEta_isLoaded;
	vector<float> photons_sigmaIEtaIEta_;
	TBranch *photons_sigmaIEtaIEta_branch;
	bool photons_sigmaIEtaIEta_isLoaded;
	vector<float> photons_swissSeed_;
	TBranch *photons_swissSeed_branch;
	bool photons_swissSeed_isLoaded;
	vector<float> photons_tkIsoHollow03_;
	TBranch *photons_tkIsoHollow03_branch;
	bool photons_tkIsoHollow03_isLoaded;
	vector<float> photons_tkIsoHollow04_;
	TBranch *photons_tkIsoHollow04_branch;
	bool photons_tkIsoHollow04_isLoaded;
	vector<float> photons_tkIsoSolid03_;
	TBranch *photons_tkIsoSolid03_branch;
	bool photons_tkIsoSolid03_isLoaded;
	vector<float> photons_tkIsoSolid04_;
	TBranch *photons_tkIsoSolid04_branch;
	bool photons_tkIsoSolid04_isLoaded;
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
	vector<float> scs_sigmaIEtaIEtaSC_;
	TBranch *scs_sigmaIEtaIEtaSC_branch;
	bool scs_sigmaIEtaIEtaSC_isLoaded;
	vector<float> scs_sigmaIEtaIPhi_;
	TBranch *scs_sigmaIEtaIPhi_branch;
	bool scs_sigmaIEtaIPhi_isLoaded;
	vector<float> scs_sigmaIEtaIPhiSC_;
	TBranch *scs_sigmaIEtaIPhiSC_branch;
	bool scs_sigmaIEtaIPhiSC_isLoaded;
	vector<float> scs_sigmaIPhiIPhi_;
	TBranch *scs_sigmaIPhiIPhi_branch;
	bool scs_sigmaIPhiIPhi_isLoaded;
	vector<float> scs_sigmaIPhiIPhiSC_;
	TBranch *scs_sigmaIPhiIPhiSC_branch;
	bool scs_sigmaIPhiIPhiSC_isLoaded;
	vector<float> scs_sigmaPhiPhi_;
	TBranch *scs_sigmaPhiPhi_branch;
	bool scs_sigmaPhiPhi_isLoaded;
	vector<float> scs_timeSeed_;
	TBranch *scs_timeSeed_branch;
	bool scs_timeSeed_isLoaded;
	vector<float> svs_anglePV_;
	TBranch *svs_anglePV_branch;
	bool svs_anglePV_isLoaded;
	vector<float> svs_chi2_;
	TBranch *svs_chi2_branch;
	bool svs_chi2_isLoaded;
	vector<float> svs_dist3Dsig_;
	TBranch *svs_dist3Dsig_branch;
	bool svs_dist3Dsig_isLoaded;
	vector<float> svs_dist3Dval_;
	TBranch *svs_dist3Dval_branch;
	bool svs_dist3Dval_isLoaded;
	vector<float> svs_distXYsig_;
	TBranch *svs_distXYsig_branch;
	bool svs_distXYsig_isLoaded;
	vector<float> svs_distXYval_;
	TBranch *svs_distXYval_branch;
	bool svs_distXYval_isLoaded;
	vector<float> svs_ndof_;
	TBranch *svs_ndof_branch;
	bool svs_ndof_isLoaded;
	vector<float> svs_prob_;
	TBranch *svs_prob_branch;
	bool svs_prob_isLoaded;
	vector<float> svs_xError_;
	TBranch *svs_xError_branch;
	bool svs_xError_isLoaded;
	vector<float> svs_yError_;
	TBranch *svs_yError_branch;
	bool svs_yError_isLoaded;
	vector<float> svs_zError_;
	TBranch *svs_zError_branch;
	bool svs_zError_isLoaded;
	vector<float> mus_tcmet_deltax_;
	TBranch *mus_tcmet_deltax_branch;
	bool mus_tcmet_deltax_isLoaded;
	vector<float> mus_tcmet_deltay_;
	TBranch *mus_tcmet_deltay_branch;
	bool mus_tcmet_deltay_isLoaded;
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
	vector<float> trks_d0phiCov_;
	TBranch *trks_d0phiCov_branch;
	bool trks_d0phiCov_isLoaded;
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
	vector<float> trkjets_cor_;
	TBranch *trkjets_cor_branch;
	bool trkjets_cor_isLoaded;
	vector<float> trks_d0Errvtx_;
	TBranch *trks_d0Errvtx_branch;
	bool trks_d0Errvtx_isLoaded;
	vector<float> trks_d0vtx_;
	TBranch *trks_d0vtx_branch;
	bool trks_d0vtx_isLoaded;
	vector<float> vtxs_chi2_;
	TBranch *vtxs_chi2_branch;
	bool vtxs_chi2_isLoaded;
	vector<float> vtxs_ndof_;
	TBranch *vtxs_ndof_branch;
	bool vtxs_ndof_isLoaded;
	vector<float> vtxs_sumpt_;
	TBranch *vtxs_sumpt_branch;
	bool vtxs_sumpt_isLoaded;
	vector<float> vtxs_xError_;
	TBranch *vtxs_xError_branch;
	bool vtxs_xError_isLoaded;
	vector<float> vtxs_yError_;
	TBranch *vtxs_yError_branch;
	bool vtxs_yError_isLoaded;
	vector<float> vtxs_zError_;
	TBranch *vtxs_zError_branch;
	bool vtxs_zError_isLoaded;
	vector<vector<float> > vtxs_covMatrix_;
	TBranch *vtxs_covMatrix_branch;
	bool vtxs_covMatrix_isLoaded;
	int evt_cscLooseHaloId_;
	TBranch *evt_cscLooseHaloId_branch;
	bool evt_cscLooseHaloId_isLoaded;
	int evt_cscTightHaloId_;
	TBranch *evt_cscTightHaloId_branch;
	bool evt_cscTightHaloId_isLoaded;
	int evt_ecalLooseHaloId_;
	TBranch *evt_ecalLooseHaloId_branch;
	bool evt_ecalLooseHaloId_isLoaded;
	int evt_ecalTightHaloId_;
	TBranch *evt_ecalTightHaloId_branch;
	bool evt_ecalTightHaloId_isLoaded;
	int evt_extremeTightHaloId_;
	TBranch *evt_extremeTightHaloId_branch;
	bool evt_extremeTightHaloId_isLoaded;
	int evt_globalLooseHaloId_;
	TBranch *evt_globalLooseHaloId_branch;
	bool evt_globalLooseHaloId_isLoaded;
	int evt_globalTightHaloId_;
	TBranch *evt_globalTightHaloId_branch;
	bool evt_globalTightHaloId_isLoaded;
	int evt_hcalLooseHaloId_;
	TBranch *evt_hcalLooseHaloId_branch;
	bool evt_hcalLooseHaloId_isLoaded;
	int evt_hcalTightHaloId_;
	TBranch *evt_hcalTightHaloId_branch;
	bool evt_hcalTightHaloId_isLoaded;
	int evt_looseHaloId_;
	TBranch *evt_looseHaloId_branch;
	bool evt_looseHaloId_isLoaded;
	int evt_nHaloLikeTracks_;
	TBranch *evt_nHaloLikeTracks_branch;
	bool evt_nHaloLikeTracks_isLoaded;
	int evt_nHaloTriggerCandidates_;
	TBranch *evt_nHaloTriggerCandidates_branch;
	bool evt_nHaloTriggerCandidates_isLoaded;
	int evt_tightHaloId_;
	TBranch *evt_tightHaloId_branch;
	bool evt_tightHaloId_isLoaded;
	int evt_bsType_;
	TBranch *evt_bsType_branch;
	bool evt_bsType_isLoaded;
	int evt_bunchCrossing_;
	TBranch *evt_bunchCrossing_branch;
	bool evt_bunchCrossing_isLoaded;
	int evt_experimentType_;
	TBranch *evt_experimentType_branch;
	bool evt_experimentType_isLoaded;
	int evt_isRealData_;
	TBranch *evt_isRealData_branch;
	bool evt_isRealData_isLoaded;
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
	int pdfinfo_id1_;
	TBranch *pdfinfo_id1_branch;
	bool pdfinfo_id1_isLoaded;
	int pdfinfo_id2_;
	TBranch *pdfinfo_id2_branch;
	bool pdfinfo_id2_isLoaded;
	vector<int> evt_ecaliPhiSuspects_;
	TBranch *evt_ecaliPhiSuspects_branch;
	bool evt_ecaliPhiSuspects_isLoaded;
	vector<int> evt_globaliPhiSuspects_;
	TBranch *evt_globaliPhiSuspects_branch;
	bool evt_globaliPhiSuspects_isLoaded;
	vector<int> evt_hcaliPhiSuspects_;
	TBranch *evt_hcaliPhiSuspects_branch;
	bool evt_hcaliPhiSuspects_isLoaded;
	vector<int> els_mc3_id_;
	TBranch *els_mc3_id_branch;
	bool els_mc3_id_isLoaded;
	vector<int> els_mc3idx_;
	TBranch *els_mc3idx_branch;
	bool els_mc3idx_isLoaded;
	vector<int> els_mc3_motherid_;
	TBranch *els_mc3_motherid_branch;
	bool els_mc3_motherid_isLoaded;
	vector<int> els_mc3_motheridx_;
	TBranch *els_mc3_motheridx_branch;
	bool els_mc3_motheridx_isLoaded;
	vector<int> els_mc_id_;
	TBranch *els_mc_id_branch;
	bool els_mc_id_isLoaded;
	vector<int> els_mcidx_;
	TBranch *els_mcidx_branch;
	bool els_mcidx_isLoaded;
	vector<int> els_mc_motherid_;
	TBranch *els_mc_motherid_branch;
	bool els_mc_motherid_isLoaded;
	vector<int> jets_mc3_id_;
	TBranch *jets_mc3_id_branch;
	bool jets_mc3_id_isLoaded;
	vector<int> jets_mc3idx_;
	TBranch *jets_mc3idx_branch;
	bool jets_mc3idx_isLoaded;
	vector<int> jets_mc_gpidx_;
	TBranch *jets_mc_gpidx_branch;
	bool jets_mc_gpidx_isLoaded;
	vector<int> jets_mc_id_;
	TBranch *jets_mc_id_branch;
	bool jets_mc_id_isLoaded;
	vector<int> jets_mcidx_;
	TBranch *jets_mcidx_branch;
	bool jets_mcidx_isLoaded;
	vector<int> jets_mc_motherid_;
	TBranch *jets_mc_motherid_branch;
	bool jets_mc_motherid_isLoaded;
	vector<int> mus_mc3_id_;
	TBranch *mus_mc3_id_branch;
	bool mus_mc3_id_isLoaded;
	vector<int> mus_mc3idx_;
	TBranch *mus_mc3idx_branch;
	bool mus_mc3idx_isLoaded;
	vector<int> mus_mc3_motherid_;
	TBranch *mus_mc3_motherid_branch;
	bool mus_mc3_motherid_isLoaded;
	vector<int> mus_mc3_motheridx_;
	TBranch *mus_mc3_motheridx_branch;
	bool mus_mc3_motheridx_isLoaded;
	vector<int> mus_mc_id_;
	TBranch *mus_mc_id_branch;
	bool mus_mc_id_isLoaded;
	vector<int> mus_mcidx_;
	TBranch *mus_mcidx_branch;
	bool mus_mcidx_isLoaded;
	vector<int> mus_mc_motherid_;
	TBranch *mus_mc_motherid_branch;
	bool mus_mc_motherid_isLoaded;
	vector<int> pfjets_mc3_id_;
	TBranch *pfjets_mc3_id_branch;
	bool pfjets_mc3_id_isLoaded;
	vector<int> pfjets_mc3idx_;
	TBranch *pfjets_mc3idx_branch;
	bool pfjets_mc3idx_isLoaded;
	vector<int> pfjets_mc_gpidx_;
	TBranch *pfjets_mc_gpidx_branch;
	bool pfjets_mc_gpidx_isLoaded;
	vector<int> pfjets_mc_id_;
	TBranch *pfjets_mc_id_branch;
	bool pfjets_mc_id_isLoaded;
	vector<int> pfjets_mcidx_;
	TBranch *pfjets_mcidx_branch;
	bool pfjets_mcidx_isLoaded;
	vector<int> pfjets_mc_motherid_;
	TBranch *pfjets_mc_motherid_branch;
	bool pfjets_mc_motherid_isLoaded;
	vector<int> photons_mc3_id_;
	TBranch *photons_mc3_id_branch;
	bool photons_mc3_id_isLoaded;
	vector<int> photons_mc3idx_;
	TBranch *photons_mc3idx_branch;
	bool photons_mc3idx_isLoaded;
	vector<int> photons_mc3_motherid_;
	TBranch *photons_mc3_motherid_branch;
	bool photons_mc3_motherid_isLoaded;
	vector<int> photons_mc3_motheridx_;
	TBranch *photons_mc3_motheridx_branch;
	bool photons_mc3_motheridx_isLoaded;
	vector<int> photons_mc_id_;
	TBranch *photons_mc_id_branch;
	bool photons_mc_id_isLoaded;
	vector<int> photons_mcidx_;
	TBranch *photons_mcidx_branch;
	bool photons_mcidx_isLoaded;
	vector<int> photons_mc_motherid_;
	TBranch *photons_mc_motherid_branch;
	bool photons_mc_motherid_isLoaded;
	vector<int> trk_mc3_id_;
	TBranch *trk_mc3_id_branch;
	bool trk_mc3_id_isLoaded;
	vector<int> trk_mc3idx_;
	TBranch *trk_mc3idx_branch;
	bool trk_mc3idx_isLoaded;
	vector<int> trk_mc3_motherid_;
	TBranch *trk_mc3_motherid_branch;
	bool trk_mc3_motherid_isLoaded;
	vector<int> trk_mc3_motheridx_;
	TBranch *trk_mc3_motheridx_branch;
	bool trk_mc3_motheridx_isLoaded;
	vector<int> trk_mc_id_;
	TBranch *trk_mc_id_branch;
	bool trk_mc_id_isLoaded;
	vector<int> trk_mcidx_;
	TBranch *trk_mcidx_branch;
	bool trk_mcidx_isLoaded;
	vector<int> trk_mc_motherid_;
	TBranch *trk_mc_motherid_branch;
	bool trk_mc_motherid_isLoaded;
	vector<int> trks_conv_tkidx_;
	TBranch *trks_conv_tkidx_branch;
	bool trks_conv_tkidx_isLoaded;
	vector<int> els_exp_innerlayers39X_;
	TBranch *els_exp_innerlayers39X_branch;
	bool els_exp_innerlayers39X_isLoaded;
	vector<int> els_closestJet_;
	TBranch *els_closestJet_branch;
	bool els_closestJet_isLoaded;
	vector<int> els_closestMuon_;
	TBranch *els_closestMuon_branch;
	bool els_closestMuon_isLoaded;
	vector<int> els_pfelsidx_;
	TBranch *els_pfelsidx_branch;
	bool els_pfelsidx_isLoaded;
	vector<int> els_category_;
	TBranch *els_category_branch;
	bool els_category_isLoaded;
	vector<int> els_charge_;
	TBranch *els_charge_branch;
	bool els_charge_isLoaded;
	vector<int> els_class_;
	TBranch *els_class_branch;
	bool els_class_isLoaded;
	vector<int> els_conv_tkidx_;
	TBranch *els_conv_tkidx_branch;
	bool els_conv_tkidx_isLoaded;
	vector<int> els_exp_innerlayers_;
	TBranch *els_exp_innerlayers_branch;
	bool els_exp_innerlayers_isLoaded;
	vector<int> els_exp_outerlayers_;
	TBranch *els_exp_outerlayers_branch;
	bool els_exp_outerlayers_isLoaded;
	vector<int> els_fiduciality_;
	TBranch *els_fiduciality_branch;
	bool els_fiduciality_isLoaded;
	vector<int> els_gsftrkidx_;
	TBranch *els_gsftrkidx_branch;
	bool els_gsftrkidx_isLoaded;
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
	vector<int> els_lostHits_;
	TBranch *els_lostHits_branch;
	bool els_lostHits_isLoaded;
	vector<int> els_lost_pixelhits_;
	TBranch *els_lost_pixelhits_branch;
	bool els_lost_pixelhits_isLoaded;
	vector<int> els_nSeed_;
	TBranch *els_nSeed_branch;
	bool els_nSeed_isLoaded;
	vector<int> els_sccharge_;
	TBranch *els_sccharge_branch;
	bool els_sccharge_isLoaded;
	vector<int> els_scindex_;
	TBranch *els_scindex_branch;
	bool els_scindex_isLoaded;
	vector<int> els_trk_charge_;
	TBranch *els_trk_charge_branch;
	bool els_trk_charge_isLoaded;
	vector<int> els_trkidx_;
	TBranch *els_trkidx_branch;
	bool els_trkidx_isLoaded;
	vector<int> els_type_;
	TBranch *els_type_branch;
	bool els_type_isLoaded;
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
	vector<int> genps_status_;
	TBranch *genps_status_branch;
	bool genps_status_isLoaded;
	vector<int> gsftrks_charge_;
	TBranch *gsftrks_charge_branch;
	bool gsftrks_charge_isLoaded;
	vector<int> gsftrks_exp_innerlayers_;
	TBranch *gsftrks_exp_innerlayers_branch;
	bool gsftrks_exp_innerlayers_isLoaded;
	vector<int> gsftrks_exp_outerlayers_;
	TBranch *gsftrks_exp_outerlayers_branch;
	bool gsftrks_exp_outerlayers_isLoaded;
	vector<int> gsftrks_layer1_det_;
	TBranch *gsftrks_layer1_det_branch;
	bool gsftrks_layer1_det_isLoaded;
	vector<int> gsftrks_layer1_layer_;
	TBranch *gsftrks_layer1_layer_branch;
	bool gsftrks_layer1_layer_isLoaded;
	vector<int> gsftrks_layer1_sizerphi_;
	TBranch *gsftrks_layer1_sizerphi_branch;
	bool gsftrks_layer1_sizerphi_isLoaded;
	vector<int> gsftrks_layer1_sizerz_;
	TBranch *gsftrks_layer1_sizerz_branch;
	bool gsftrks_layer1_sizerz_isLoaded;
	vector<int> gsftrks_lostHits_;
	TBranch *gsftrks_lostHits_branch;
	bool gsftrks_lostHits_isLoaded;
	vector<int> gsftrks_lost_pixelhits_;
	TBranch *gsftrks_lost_pixelhits_branch;
	bool gsftrks_lost_pixelhits_isLoaded;
	vector<int> gsftrks_nlayers_;
	TBranch *gsftrks_nlayers_branch;
	bool gsftrks_nlayers_isLoaded;
	vector<int> gsftrks_nlayers3D_;
	TBranch *gsftrks_nlayers3D_branch;
	bool gsftrks_nlayers3D_isLoaded;
	vector<int> gsftrks_nlayersLost_;
	TBranch *gsftrks_nlayersLost_branch;
	bool gsftrks_nlayersLost_isLoaded;
	vector<int> gsftrks_validHits_;
	TBranch *gsftrks_validHits_branch;
	bool gsftrks_validHits_isLoaded;
	vector<int> gsftrks_valid_pixelhits_;
	TBranch *gsftrks_valid_pixelhits_branch;
	bool gsftrks_valid_pixelhits_isLoaded;
	vector<int> hyp_ll_charge_;
	TBranch *hyp_ll_charge_branch;
	bool hyp_ll_charge_isLoaded;
	vector<int> hyp_ll_id_;
	TBranch *hyp_ll_id_branch;
	bool hyp_ll_id_isLoaded;
	vector<int> hyp_ll_index_;
	TBranch *hyp_ll_index_branch;
	bool hyp_ll_index_isLoaded;
	vector<int> hyp_ll_lostHits_;
	TBranch *hyp_ll_lostHits_branch;
	bool hyp_ll_lostHits_isLoaded;
	vector<int> hyp_ll_validHits_;
	TBranch *hyp_ll_validHits_branch;
	bool hyp_ll_validHits_isLoaded;
	vector<int> hyp_lt_charge_;
	TBranch *hyp_lt_charge_branch;
	bool hyp_lt_charge_isLoaded;
	vector<int> hyp_lt_id_;
	TBranch *hyp_lt_id_branch;
	bool hyp_lt_id_isLoaded;
	vector<int> hyp_lt_index_;
	TBranch *hyp_lt_index_branch;
	bool hyp_lt_index_isLoaded;
	vector<int> hyp_lt_lostHits_;
	TBranch *hyp_lt_lostHits_branch;
	bool hyp_lt_lostHits_isLoaded;
	vector<int> hyp_lt_validHits_;
	TBranch *hyp_lt_validHits_branch;
	bool hyp_lt_validHits_isLoaded;
	vector<int> hyp_njets_;
	TBranch *hyp_njets_branch;
	bool hyp_njets_isLoaded;
	vector<int> hyp_nojets_;
	TBranch *hyp_nojets_branch;
	bool hyp_nojets_isLoaded;
	vector<int> hyp_type_;
	TBranch *hyp_type_branch;
	bool hyp_type_isLoaded;
	vector<int> hyp_FVFit_ndf_;
	TBranch *hyp_FVFit_ndf_branch;
	bool hyp_FVFit_ndf_isLoaded;
	vector<int> hyp_FVFit_status_;
	TBranch *hyp_FVFit_status_branch;
	bool hyp_FVFit_status_isLoaded;
	vector<int> hyp_ll_mc_id_;
	TBranch *hyp_ll_mc_id_branch;
	bool hyp_ll_mc_id_isLoaded;
	vector<int> hyp_ll_mc_motherid_;
	TBranch *hyp_ll_mc_motherid_branch;
	bool hyp_ll_mc_motherid_isLoaded;
	vector<int> hyp_lt_mc_id_;
	TBranch *hyp_lt_mc_id_branch;
	bool hyp_lt_mc_id_isLoaded;
	vector<int> hyp_lt_mc_motherid_;
	TBranch *hyp_lt_mc_motherid_branch;
	bool hyp_lt_mc_motherid_isLoaded;
	vector<int> hyp_quadlep_first_type_;
	TBranch *hyp_quadlep_first_type_branch;
	bool hyp_quadlep_first_type_isLoaded;
	vector<int> hyp_quadlep_fourth_type_;
	TBranch *hyp_quadlep_fourth_type_branch;
	bool hyp_quadlep_fourth_type_isLoaded;
	vector<int> hyp_quadlep_second_type_;
	TBranch *hyp_quadlep_second_type_branch;
	bool hyp_quadlep_second_type_isLoaded;
	vector<int> hyp_quadlep_third_type_;
	TBranch *hyp_quadlep_third_type_branch;
	bool hyp_quadlep_third_type_isLoaded;
	vector<int> hyp_trilep_first_type_;
	TBranch *hyp_trilep_first_type_branch;
	bool hyp_trilep_first_type_isLoaded;
	vector<int> hyp_trilep_second_type_;
	TBranch *hyp_trilep_second_type_branch;
	bool hyp_trilep_second_type_isLoaded;
	vector<int> hyp_trilep_third_type_;
	TBranch *hyp_trilep_third_type_branch;
	bool hyp_trilep_third_type_isLoaded;
	vector<int> jets_closestElectron_;
	TBranch *jets_closestElectron_branch;
	bool jets_closestElectron_isLoaded;
	vector<int> jets_closestMuon_;
	TBranch *jets_closestMuon_branch;
	bool jets_closestMuon_isLoaded;
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
	vector<int> mus_met_flag_;
	TBranch *mus_met_flag_branch;
	bool mus_met_flag_isLoaded;
	vector<int> mus_closestEle_;
	TBranch *mus_closestEle_branch;
	bool mus_closestEle_isLoaded;
	vector<int> mus_closestJet_;
	TBranch *mus_closestJet_branch;
	bool mus_closestJet_isLoaded;
	vector<int> mus_pfmusidx_;
	TBranch *mus_pfmusidx_branch;
	bool mus_pfmusidx_isLoaded;
	vector<int> mus_charge_;
	TBranch *mus_charge_branch;
	bool mus_charge_isLoaded;
	vector<int> mus_chi2LocalMomentum_;
	TBranch *mus_chi2LocalMomentum_branch;
	bool mus_chi2LocalMomentum_isLoaded;
	vector<int> mus_chi2LocalPosition_;
	TBranch *mus_chi2LocalPosition_branch;
	bool mus_chi2LocalPosition_isLoaded;
	vector<int> mus_gfit_validHits_;
	TBranch *mus_gfit_validHits_branch;
	bool mus_gfit_validHits_isLoaded;
	vector<int> mus_gfit_validSTAHits_;
	TBranch *mus_gfit_validSTAHits_branch;
	bool mus_gfit_validSTAHits_isLoaded;
	vector<int> mus_gfit_validSiHits_;
	TBranch *mus_gfit_validSiHits_branch;
	bool mus_gfit_validSiHits_isLoaded;
	vector<int> mus_glbKink_;
	TBranch *mus_glbKink_branch;
	bool mus_glbKink_isLoaded;
	vector<int> mus_glbTrackProbability_;
	TBranch *mus_glbTrackProbability_branch;
	bool mus_glbTrackProbability_isLoaded;
	vector<int> mus_globalDeltaEtaPhi_;
	TBranch *mus_globalDeltaEtaPhi_branch;
	bool mus_globalDeltaEtaPhi_isLoaded;
	vector<int> mus_goodmask_;
	TBranch *mus_goodmask_branch;
	bool mus_goodmask_isLoaded;
	vector<int> mus_iso03_ntrk_;
	TBranch *mus_iso03_ntrk_branch;
	bool mus_iso03_ntrk_isLoaded;
	vector<int> mus_iso05_ntrk_;
	TBranch *mus_iso05_ntrk_branch;
	bool mus_iso05_ntrk_isLoaded;
	vector<int> mus_localDistance_;
	TBranch *mus_localDistance_branch;
	bool mus_localDistance_isLoaded;
	vector<int> mus_lostHits_;
	TBranch *mus_lostHits_branch;
	bool mus_lostHits_isLoaded;
	vector<int> mus_nOverlaps_;
	TBranch *mus_nOverlaps_branch;
	bool mus_nOverlaps_isLoaded;
	vector<int> mus_nmatches_;
	TBranch *mus_nmatches_branch;
	bool mus_nmatches_isLoaded;
	vector<int> mus_overlap0_;
	TBranch *mus_overlap0_branch;
	bool mus_overlap0_isLoaded;
	vector<int> mus_overlap1_;
	TBranch *mus_overlap1_branch;
	bool mus_overlap1_isLoaded;
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
	vector<int> mus_staRelChi2_;
	TBranch *mus_staRelChi2_branch;
	bool mus_staRelChi2_isLoaded;
	vector<int> mus_sta_validHits_;
	TBranch *mus_sta_validHits_branch;
	bool mus_sta_validHits_isLoaded;
	vector<int> mus_timeDirection_;
	TBranch *mus_timeDirection_branch;
	bool mus_timeDirection_isLoaded;
	vector<int> mus_timeNumStationsUsed_;
	TBranch *mus_timeNumStationsUsed_branch;
	bool mus_timeNumStationsUsed_isLoaded;
	vector<int> mus_trkKink_;
	TBranch *mus_trkKink_branch;
	bool mus_trkKink_isLoaded;
	vector<int> mus_trkRelChi2_;
	TBranch *mus_trkRelChi2_branch;
	bool mus_trkRelChi2_isLoaded;
	vector<int> mus_trk_charge_;
	TBranch *mus_trk_charge_branch;
	bool mus_trk_charge_isLoaded;
	vector<int> mus_trkidx_;
	TBranch *mus_trkidx_branch;
	bool mus_trkidx_isLoaded;
	vector<int> mus_type_;
	TBranch *mus_type_branch;
	bool mus_type_isLoaded;
	vector<int> mus_validHits_;
	TBranch *mus_validHits_branch;
	bool mus_validHits_isLoaded;
	vector<int> els_pat_genID_;
	TBranch *els_pat_genID_branch;
	bool els_pat_genID_isLoaded;
	vector<int> els_pat_genMotherID_;
	TBranch *els_pat_genMotherID_branch;
	bool els_pat_genMotherID_isLoaded;
	vector<int> jets_pat_genPartonMother_id_;
	TBranch *jets_pat_genPartonMother_id_branch;
	bool jets_pat_genPartonMother_id_isLoaded;
	vector<int> jets_pat_genParton_id_;
	TBranch *jets_pat_genParton_id_branch;
	bool jets_pat_genParton_id_isLoaded;
	vector<int> jets_pat_jetIDLoose_;
	TBranch *jets_pat_jetIDLoose_branch;
	bool jets_pat_jetIDLoose_isLoaded;
	vector<int> jets_pat_jetIDLooseAOD_;
	TBranch *jets_pat_jetIDLooseAOD_branch;
	bool jets_pat_jetIDLooseAOD_isLoaded;
	vector<int> jets_pat_jetIDMinimal_;
	TBranch *jets_pat_jetIDMinimal_branch;
	bool jets_pat_jetIDMinimal_isLoaded;
	vector<int> jets_pat_jetIDTight_;
	TBranch *jets_pat_jetIDTight_branch;
	bool jets_pat_jetIDTight_isLoaded;
	vector<int> jets_pat_partonFlavour_;
	TBranch *jets_pat_partonFlavour_branch;
	bool jets_pat_partonFlavour_isLoaded;
	vector<int> mus_pat_genID_;
	TBranch *mus_pat_genID_branch;
	bool mus_pat_genID_isLoaded;
	vector<int> mus_pat_genMotherID_;
	TBranch *mus_pat_genMotherID_branch;
	bool mus_pat_genMotherID_isLoaded;
	vector<int> pfels_elsidx_;
	TBranch *pfels_elsidx_branch;
	bool pfels_elsidx_isLoaded;
	vector<int> pfels_charge_;
	TBranch *pfels_charge_branch;
	bool pfels_charge_isLoaded;
	vector<int> pfels_flag_;
	TBranch *pfels_flag_branch;
	bool pfels_flag_isLoaded;
	vector<int> pfels_particleId_;
	TBranch *pfels_particleId_branch;
	bool pfels_particleId_isLoaded;
	vector<int> pfjets_chargedMultiplicity_;
	TBranch *pfjets_chargedMultiplicity_branch;
	bool pfjets_chargedMultiplicity_isLoaded;
	vector<int> pfjets_muonMultiplicity_;
	TBranch *pfjets_muonMultiplicity_branch;
	bool pfjets_muonMultiplicity_isLoaded;
	vector<int> pfjets_neutralMultiplicity_;
	TBranch *pfjets_neutralMultiplicity_branch;
	bool pfjets_neutralMultiplicity_isLoaded;
	vector<int> pfmus_musidx_;
	TBranch *pfmus_musidx_branch;
	bool pfmus_musidx_isLoaded;
	vector<int> pfmus_charge_;
	TBranch *pfmus_charge_branch;
	bool pfmus_charge_isLoaded;
	vector<int> pfmus_flag_;
	TBranch *pfmus_flag_branch;
	bool pfmus_flag_isLoaded;
	vector<int> pfmus_particleId_;
	TBranch *pfmus_particleId_branch;
	bool pfmus_particleId_isLoaded;
	vector<int> photons_fiduciality_;
	TBranch *photons_fiduciality_branch;
	bool photons_fiduciality_isLoaded;
	vector<int> photons_scindex_;
	TBranch *photons_scindex_branch;
	bool photons_scindex_isLoaded;
	vector<int> scs_detIdSeed_;
	TBranch *scs_detIdSeed_branch;
	bool scs_detIdSeed_isLoaded;
	vector<int> scs_elsidx_;
	TBranch *scs_elsidx_branch;
	bool scs_elsidx_isLoaded;
	vector<int> scs_severitySeed_;
	TBranch *scs_severitySeed_branch;
	bool scs_severitySeed_isLoaded;
	vector<int> svs_isKs_;
	TBranch *svs_isKs_branch;
	bool svs_isKs_isLoaded;
	vector<int> svs_isLambda_;
	TBranch *svs_isLambda_branch;
	bool svs_isLambda_isLoaded;
	vector<int> svs_mc3_id_;
	TBranch *svs_mc3_id_branch;
	bool svs_mc3_id_isLoaded;
	vector<int> svs_nTrks_;
	TBranch *svs_nTrks_branch;
	bool svs_nTrks_isLoaded;
	vector<int> mus_tcmet_flag_;
	TBranch *mus_tcmet_flag_branch;
	bool mus_tcmet_flag_isLoaded;
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
	vector<int> trks_elsidx_;
	TBranch *trks_elsidx_branch;
	bool trks_elsidx_isLoaded;
	vector<int> trk_musidx_;
	TBranch *trk_musidx_branch;
	bool trk_musidx_isLoaded;
	vector<int> trkjets_ntrks_;
	TBranch *trkjets_ntrks_branch;
	bool trkjets_ntrks_isLoaded;
	vector<int> trkjets_vtxs_idx_;
	TBranch *trkjets_vtxs_idx_branch;
	bool trkjets_vtxs_idx_isLoaded;
	vector<int> vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int> vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<int> vtxs_tracksSize_;
	TBranch *vtxs_tracksSize_branch;
	bool vtxs_tracksSize_isLoaded;
	vector<vector<int> > genps_lepdaughter_id_;
	TBranch *genps_lepdaughter_id_branch;
	bool genps_lepdaughter_id_isLoaded;
	vector<vector<int> > genps_lepdaughter_idx_;
	TBranch *genps_lepdaughter_idx_branch;
	bool genps_lepdaughter_idx_isLoaded;
	vector<vector<int> > hlt_trigObjs_id_;
	TBranch *hlt_trigObjs_id_branch;
	bool hlt_trigObjs_id_isLoaded;
	vector<vector<int> > hyp_jets_idx_;
	TBranch *hyp_jets_idx_branch;
	bool hyp_jets_idx_isLoaded;
	vector<vector<int> > hyp_other_jets_idx_;
	TBranch *hyp_other_jets_idx_branch;
	bool hyp_other_jets_idx_isLoaded;
	unsigned int evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int evt_detectorStatus_;
	TBranch *evt_detectorStatus_branch;
	bool evt_detectorStatus_isLoaded;
	unsigned int evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int genps_flavorHistoryFilterResult_;
	TBranch *genps_flavorHistoryFilterResult_branch;
	bool genps_flavorHistoryFilterResult_isLoaded;
	unsigned int evt_ngenjets_;
	TBranch *evt_ngenjets_branch;
	bool evt_ngenjets_isLoaded;
	unsigned int genps_signalProcessID_;
	TBranch *genps_signalProcessID_branch;
	bool genps_signalProcessID_isLoaded;
	unsigned int hlt_bits1_;
	TBranch *hlt_bits1_branch;
	bool hlt_bits1_isLoaded;
	unsigned int hlt_bits2_;
	TBranch *hlt_bits2_branch;
	bool hlt_bits2_isLoaded;
	unsigned int hlt_bits3_;
	TBranch *hlt_bits3_branch;
	bool hlt_bits3_isLoaded;
	unsigned int hlt_bits4_;
	TBranch *hlt_bits4_branch;
	bool hlt_bits4_isLoaded;
	unsigned int hlt_bits5_;
	TBranch *hlt_bits5_branch;
	bool hlt_bits5_isLoaded;
	unsigned int hlt_bits6_;
	TBranch *hlt_bits6_branch;
	bool hlt_bits6_isLoaded;
	unsigned int hlt_bits7_;
	TBranch *hlt_bits7_branch;
	bool hlt_bits7_isLoaded;
	unsigned int hlt_bits8_;
	TBranch *hlt_bits8_branch;
	bool hlt_bits8_isLoaded;
	unsigned int evt_njets_;
	TBranch *evt_njets_branch;
	bool evt_njets_isLoaded;
	unsigned int evt_njpts_;
	TBranch *evt_njpts_branch;
	bool evt_njpts_isLoaded;
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
	unsigned int evt_nphotons_;
	TBranch *evt_nphotons_branch;
	bool evt_nphotons_isLoaded;
	unsigned int evt_ecalRecoStatus_;
	TBranch *evt_ecalRecoStatus_branch;
	bool evt_ecalRecoStatus_isLoaded;
	unsigned int evt_nscs_;
	TBranch *evt_nscs_branch;
	bool evt_nscs_isLoaded;
	unsigned int evt_ntrkjets_;
	TBranch *evt_ntrkjets_branch;
	bool evt_ntrkjets_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	vector<unsigned int> hlt_prescales_;
	TBranch *hlt_prescales_branch;
	bool hlt_prescales_isLoaded;
	vector<unsigned int> hyp_quadlep_bucket_;
	TBranch *hyp_quadlep_bucket_branch;
	bool hyp_quadlep_bucket_isLoaded;
	vector<unsigned int> hyp_quadlep_first_index_;
	TBranch *hyp_quadlep_first_index_branch;
	bool hyp_quadlep_first_index_isLoaded;
	vector<unsigned int> hyp_quadlep_fourth_index_;
	TBranch *hyp_quadlep_fourth_index_branch;
	bool hyp_quadlep_fourth_index_isLoaded;
	vector<unsigned int> hyp_quadlep_second_index_;
	TBranch *hyp_quadlep_second_index_branch;
	bool hyp_quadlep_second_index_isLoaded;
	vector<unsigned int> hyp_quadlep_third_index_;
	TBranch *hyp_quadlep_third_index_branch;
	bool hyp_quadlep_third_index_isLoaded;
	vector<unsigned int> hyp_trilep_bucket_;
	TBranch *hyp_trilep_bucket_branch;
	bool hyp_trilep_bucket_isLoaded;
	vector<unsigned int> hyp_trilep_first_index_;
	TBranch *hyp_trilep_first_index_branch;
	bool hyp_trilep_first_index_isLoaded;
	vector<unsigned int> hyp_trilep_second_index_;
	TBranch *hyp_trilep_second_index_branch;
	bool hyp_trilep_second_index_isLoaded;
	vector<unsigned int> hyp_trilep_third_index_;
	TBranch *hyp_trilep_third_index_branch;
	bool hyp_trilep_third_index_isLoaded;
	vector<unsigned int> l1_prescales_;
	TBranch *l1_prescales_branch;
	bool l1_prescales_isLoaded;
	vector<unsigned int> l1_techtrigprescales_;
	TBranch *l1_techtrigprescales_branch;
	bool l1_techtrigprescales_isLoaded;
	vector<unsigned int> els_pat_flag_;
	TBranch *els_pat_flag_branch;
	bool els_pat_flag_isLoaded;
	vector<unsigned int> jets_pat_flag_;
	TBranch *jets_pat_flag_branch;
	bool jets_pat_flag_isLoaded;
	vector<unsigned int> mus_pat_flag_;
	TBranch *mus_pat_flag_branch;
	bool mus_pat_flag_isLoaded;
	int	evt_nEvts_;
	TBranch *evt_nEvts_branch;
	bool evt_nEvts_isLoaded;
	float	evt_filt_eff_;
	TBranch *evt_filt_eff_branch;
	bool evt_filt_eff_isLoaded;
public: 
void Init(TTree *tree) {
	evt_bsp4_branch = 0;
	if (tree->GetAlias("evt_bsp4") != 0) {
		evt_bsp4_branch = tree->GetBranch(tree->GetAlias("evt_bsp4"));
		evt_bsp4_branch->SetAddress(&evt_bsp4_);
	}
	l1_met_p4_branch = 0;
	if (tree->GetAlias("l1_met_p4") != 0) {
		l1_met_p4_branch = tree->GetBranch(tree->GetAlias("l1_met_p4"));
		l1_met_p4_branch->SetAddress(&l1_met_p4_);
	}
	l1_mht_p4_branch = 0;
	if (tree->GetAlias("l1_mht_p4") != 0) {
		l1_mht_p4_branch = tree->GetBranch(tree->GetAlias("l1_mht_p4"));
		l1_mht_p4_branch->SetAddress(&l1_mht_p4_);
	}
	els_mc_motherp4_branch = 0;
	if (tree->GetAlias("els_mc_motherp4") != 0) {
		els_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("els_mc_motherp4"));
		els_mc_motherp4_branch->SetAddress(&els_mc_motherp4_);
	}
	els_mc_p4_branch = 0;
	if (tree->GetAlias("els_mc_p4") != 0) {
		els_mc_p4_branch = tree->GetBranch(tree->GetAlias("els_mc_p4"));
		els_mc_p4_branch->SetAddress(&els_mc_p4_);
	}
	jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("jets_mc_gp_p4") != 0) {
		jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("jets_mc_gp_p4"));
		jets_mc_gp_p4_branch->SetAddress(&jets_mc_gp_p4_);
	}
	jets_mc_motherp4_branch = 0;
	if (tree->GetAlias("jets_mc_motherp4") != 0) {
		jets_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("jets_mc_motherp4"));
		jets_mc_motherp4_branch->SetAddress(&jets_mc_motherp4_);
	}
	jets_mc_p4_branch = 0;
	if (tree->GetAlias("jets_mc_p4") != 0) {
		jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("jets_mc_p4"));
		jets_mc_p4_branch->SetAddress(&jets_mc_p4_);
	}
	mus_mc_motherp4_branch = 0;
	if (tree->GetAlias("mus_mc_motherp4") != 0) {
		mus_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherp4"));
		mus_mc_motherp4_branch->SetAddress(&mus_mc_motherp4_);
	}
	mus_mc_p4_branch = 0;
	if (tree->GetAlias("mus_mc_p4") != 0) {
		mus_mc_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc_p4"));
		mus_mc_p4_branch->SetAddress(&mus_mc_p4_);
	}
	pfjets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("pfjets_mc_gp_p4") != 0) {
		pfjets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gp_p4"));
		pfjets_mc_gp_p4_branch->SetAddress(&pfjets_mc_gp_p4_);
	}
	pfjets_mc_motherp4_branch = 0;
	if (tree->GetAlias("pfjets_mc_motherp4") != 0) {
		pfjets_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_motherp4"));
		pfjets_mc_motherp4_branch->SetAddress(&pfjets_mc_motherp4_);
	}
	pfjets_mc_p4_branch = 0;
	if (tree->GetAlias("pfjets_mc_p4") != 0) {
		pfjets_mc_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_p4"));
		pfjets_mc_p4_branch->SetAddress(&pfjets_mc_p4_);
	}
	photons_mc_motherp4_branch = 0;
	if (tree->GetAlias("photons_mc_motherp4") != 0) {
		photons_mc_motherp4_branch = tree->GetBranch(tree->GetAlias("photons_mc_motherp4"));
		photons_mc_motherp4_branch->SetAddress(&photons_mc_motherp4_);
	}
	photons_mc_p4_branch = 0;
	if (tree->GetAlias("photons_mc_p4") != 0) {
		photons_mc_p4_branch = tree->GetBranch(tree->GetAlias("photons_mc_p4"));
		photons_mc_p4_branch->SetAddress(&photons_mc_p4_);
	}
	trk_mcp4_branch = 0;
	if (tree->GetAlias("trk_mcp4") != 0) {
		trk_mcp4_branch = tree->GetBranch(tree->GetAlias("trk_mcp4"));
		trk_mcp4_branch->SetAddress(&trk_mcp4_);
	}
	els_conv_pos_p4_branch = 0;
	if (tree->GetAlias("els_conv_pos_p4") != 0) {
		els_conv_pos_p4_branch = tree->GetBranch(tree->GetAlias("els_conv_pos_p4"));
		els_conv_pos_p4_branch->SetAddress(&els_conv_pos_p4_);
	}
	els_inner_position_branch = 0;
	if (tree->GetAlias("els_inner_position") != 0) {
		els_inner_position_branch = tree->GetBranch(tree->GetAlias("els_inner_position"));
		els_inner_position_branch->SetAddress(&els_inner_position_);
	}
	els_outer_position_branch = 0;
	if (tree->GetAlias("els_outer_position") != 0) {
		els_outer_position_branch = tree->GetBranch(tree->GetAlias("els_outer_position"));
		els_outer_position_branch->SetAddress(&els_outer_position_);
	}
	els_p4_branch = 0;
	if (tree->GetAlias("els_p4") != 0) {
		els_p4_branch = tree->GetBranch(tree->GetAlias("els_p4"));
		els_p4_branch->SetAddress(&els_p4_);
	}
	els_p4In_branch = 0;
	if (tree->GetAlias("els_p4In") != 0) {
		els_p4In_branch = tree->GetBranch(tree->GetAlias("els_p4In"));
		els_p4In_branch->SetAddress(&els_p4In_);
	}
	els_p4Out_branch = 0;
	if (tree->GetAlias("els_p4Out") != 0) {
		els_p4Out_branch = tree->GetBranch(tree->GetAlias("els_p4Out"));
		els_p4Out_branch->SetAddress(&els_p4Out_);
	}
	els_trk_p4_branch = 0;
	if (tree->GetAlias("els_trk_p4") != 0) {
		els_trk_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_p4"));
		els_trk_p4_branch->SetAddress(&els_trk_p4_);
	}
	els_vertex_p4_branch = 0;
	if (tree->GetAlias("els_vertex_p4") != 0) {
		els_vertex_p4_branch = tree->GetBranch(tree->GetAlias("els_vertex_p4"));
		els_vertex_p4_branch->SetAddress(&els_vertex_p4_);
	}
	genjets_p4_branch = 0;
	if (tree->GetAlias("genjets_p4") != 0) {
		genjets_p4_branch = tree->GetBranch(tree->GetAlias("genjets_p4"));
		genjets_p4_branch->SetAddress(&genjets_p4_);
	}
	genps_p4_branch = 0;
	if (tree->GetAlias("genps_p4") != 0) {
		genps_p4_branch = tree->GetBranch(tree->GetAlias("genps_p4"));
		genps_p4_branch->SetAddress(&genps_p4_);
	}
	genps_prod_vtx_branch = 0;
	if (tree->GetAlias("genps_prod_vtx") != 0) {
		genps_prod_vtx_branch = tree->GetBranch(tree->GetAlias("genps_prod_vtx"));
		genps_prod_vtx_branch->SetAddress(&genps_prod_vtx_);
	}
	gsftrks_inner_position_branch = 0;
	if (tree->GetAlias("gsftrks_inner_position") != 0) {
		gsftrks_inner_position_branch = tree->GetBranch(tree->GetAlias("gsftrks_inner_position"));
		gsftrks_inner_position_branch->SetAddress(&gsftrks_inner_position_);
	}
	gsftrks_outer_p4_branch = 0;
	if (tree->GetAlias("gsftrks_outer_p4") != 0) {
		gsftrks_outer_p4_branch = tree->GetBranch(tree->GetAlias("gsftrks_outer_p4"));
		gsftrks_outer_p4_branch->SetAddress(&gsftrks_outer_p4_);
	}
	gsftrks_outer_position_branch = 0;
	if (tree->GetAlias("gsftrks_outer_position") != 0) {
		gsftrks_outer_position_branch = tree->GetBranch(tree->GetAlias("gsftrks_outer_position"));
		gsftrks_outer_position_branch->SetAddress(&gsftrks_outer_position_);
	}
	gsftrks_p4_branch = 0;
	if (tree->GetAlias("gsftrks_p4") != 0) {
		gsftrks_p4_branch = tree->GetBranch(tree->GetAlias("gsftrks_p4"));
		gsftrks_p4_branch->SetAddress(&gsftrks_p4_);
	}
	gsftrks_vertex_p4_branch = 0;
	if (tree->GetAlias("gsftrks_vertex_p4") != 0) {
		gsftrks_vertex_p4_branch = tree->GetBranch(tree->GetAlias("gsftrks_vertex_p4"));
		gsftrks_vertex_p4_branch->SetAddress(&gsftrks_vertex_p4_);
	}
	hyp_ll_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_p4") != 0) {
		hyp_ll_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_p4"));
		hyp_ll_p4_branch->SetAddress(&hyp_ll_p4_);
	}
	hyp_ll_trk_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_trk_p4") != 0) {
		hyp_ll_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_trk_p4"));
		hyp_ll_trk_p4_branch->SetAddress(&hyp_ll_trk_p4_);
	}
	hyp_lt_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_p4") != 0) {
		hyp_lt_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_p4"));
		hyp_lt_p4_branch->SetAddress(&hyp_lt_p4_);
	}
	hyp_lt_trk_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_trk_p4") != 0) {
		hyp_lt_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_trk_p4"));
		hyp_lt_trk_p4_branch->SetAddress(&hyp_lt_trk_p4_);
	}
	hyp_p4_branch = 0;
	if (tree->GetAlias("hyp_p4") != 0) {
		hyp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_p4"));
		hyp_p4_branch->SetAddress(&hyp_p4_);
	}
	hyp_FVFit_p4_branch = 0;
	if (tree->GetAlias("hyp_FVFit_p4") != 0) {
		hyp_FVFit_p4_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_p4"));
		hyp_FVFit_p4_branch->SetAddress(&hyp_FVFit_p4_);
	}
	hyp_FVFit_v4_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4") != 0) {
		hyp_FVFit_v4_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4"));
		hyp_FVFit_v4_branch->SetAddress(&hyp_FVFit_v4_);
	}
	hyp_ll_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_p4") != 0) {
		hyp_ll_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_p4"));
		hyp_ll_mc_p4_branch->SetAddress(&hyp_ll_mc_p4_);
	}
	hyp_lt_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_p4") != 0) {
		hyp_lt_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_p4"));
		hyp_lt_mc_p4_branch->SetAddress(&hyp_lt_mc_p4_);
	}
	jets_p4_branch = 0;
	if (tree->GetAlias("jets_p4") != 0) {
		jets_p4_branch = tree->GetBranch(tree->GetAlias("jets_p4"));
		jets_p4_branch->SetAddress(&jets_p4_);
	}
	jets_vertex_p4_branch = 0;
	if (tree->GetAlias("jets_vertex_p4") != 0) {
		jets_vertex_p4_branch = tree->GetBranch(tree->GetAlias("jets_vertex_p4"));
		jets_vertex_p4_branch->SetAddress(&jets_vertex_p4_);
	}
	jpts_p4_branch = 0;
	if (tree->GetAlias("jpts_p4") != 0) {
		jpts_p4_branch = tree->GetBranch(tree->GetAlias("jpts_p4"));
		jpts_p4_branch->SetAddress(&jpts_p4_);
	}
	l1_emiso_p4_branch = 0;
	if (tree->GetAlias("l1_emiso_p4") != 0) {
		l1_emiso_p4_branch = tree->GetBranch(tree->GetAlias("l1_emiso_p4"));
		l1_emiso_p4_branch->SetAddress(&l1_emiso_p4_);
	}
	l1_emnoiso_p4_branch = 0;
	if (tree->GetAlias("l1_emnoiso_p4") != 0) {
		l1_emnoiso_p4_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_p4"));
		l1_emnoiso_p4_branch->SetAddress(&l1_emnoiso_p4_);
	}
	l1_jetsc_p4_branch = 0;
	if (tree->GetAlias("l1_jetsc_p4") != 0) {
		l1_jetsc_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_p4"));
		l1_jetsc_p4_branch->SetAddress(&l1_jetsc_p4_);
	}
	l1_jetsf_p4_branch = 0;
	if (tree->GetAlias("l1_jetsf_p4") != 0) {
		l1_jetsf_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_p4"));
		l1_jetsf_p4_branch->SetAddress(&l1_jetsf_p4_);
	}
	l1_jetst_p4_branch = 0;
	if (tree->GetAlias("l1_jetst_p4") != 0) {
		l1_jetst_p4_branch = tree->GetBranch(tree->GetAlias("l1_jetst_p4"));
		l1_jetst_p4_branch->SetAddress(&l1_jetst_p4_);
	}
	l1_mus_p4_branch = 0;
	if (tree->GetAlias("l1_mus_p4") != 0) {
		l1_mus_p4_branch = tree->GetBranch(tree->GetAlias("l1_mus_p4"));
		l1_mus_p4_branch->SetAddress(&l1_mus_p4_);
	}
	mus_ecalpos_p4_branch = 0;
	if (tree->GetAlias("mus_ecalpos_p4") != 0) {
		mus_ecalpos_p4_branch = tree->GetBranch(tree->GetAlias("mus_ecalpos_p4"));
		mus_ecalpos_p4_branch->SetAddress(&mus_ecalpos_p4_);
	}
	mus_fitdefault_p4_branch = 0;
	if (tree->GetAlias("mus_fitdefault_p4") != 0) {
		mus_fitdefault_p4_branch = tree->GetBranch(tree->GetAlias("mus_fitdefault_p4"));
		mus_fitdefault_p4_branch->SetAddress(&mus_fitdefault_p4_);
	}
	mus_fitfirsthit_p4_branch = 0;
	if (tree->GetAlias("mus_fitfirsthit_p4") != 0) {
		mus_fitfirsthit_p4_branch = tree->GetBranch(tree->GetAlias("mus_fitfirsthit_p4"));
		mus_fitfirsthit_p4_branch->SetAddress(&mus_fitfirsthit_p4_);
	}
	mus_fitpicky_p4_branch = 0;
	if (tree->GetAlias("mus_fitpicky_p4") != 0) {
		mus_fitpicky_p4_branch = tree->GetBranch(tree->GetAlias("mus_fitpicky_p4"));
		mus_fitpicky_p4_branch->SetAddress(&mus_fitpicky_p4_);
	}
	mus_fittev_p4_branch = 0;
	if (tree->GetAlias("mus_fittev_p4") != 0) {
		mus_fittev_p4_branch = tree->GetBranch(tree->GetAlias("mus_fittev_p4"));
		mus_fittev_p4_branch->SetAddress(&mus_fittev_p4_);
	}
	mus_gfit_outerPos_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_outerPos_p4") != 0) {
		mus_gfit_outerPos_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_outerPos_p4"));
		mus_gfit_outerPos_p4_branch->SetAddress(&mus_gfit_outerPos_p4_);
	}
	mus_gfit_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_p4") != 0) {
		mus_gfit_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_p4"));
		mus_gfit_p4_branch->SetAddress(&mus_gfit_p4_);
	}
	mus_gfit_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_gfit_vertex_p4") != 0) {
		mus_gfit_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_gfit_vertex_p4"));
		mus_gfit_vertex_p4_branch->SetAddress(&mus_gfit_vertex_p4_);
	}
	mus_p4_branch = 0;
	if (tree->GetAlias("mus_p4") != 0) {
		mus_p4_branch = tree->GetBranch(tree->GetAlias("mus_p4"));
		mus_p4_branch->SetAddress(&mus_p4_);
	}
	mus_sta_p4_branch = 0;
	if (tree->GetAlias("mus_sta_p4") != 0) {
		mus_sta_p4_branch = tree->GetBranch(tree->GetAlias("mus_sta_p4"));
		mus_sta_p4_branch->SetAddress(&mus_sta_p4_);
	}
	mus_sta_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_sta_vertex_p4") != 0) {
		mus_sta_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_sta_vertex_p4"));
		mus_sta_vertex_p4_branch->SetAddress(&mus_sta_vertex_p4_);
	}
	mus_trk_p4_branch = 0;
	if (tree->GetAlias("mus_trk_p4") != 0) {
		mus_trk_p4_branch = tree->GetBranch(tree->GetAlias("mus_trk_p4"));
		mus_trk_p4_branch->SetAddress(&mus_trk_p4_);
	}
	mus_vertex_p4_branch = 0;
	if (tree->GetAlias("mus_vertex_p4") != 0) {
		mus_vertex_p4_branch = tree->GetBranch(tree->GetAlias("mus_vertex_p4"));
		mus_vertex_p4_branch->SetAddress(&mus_vertex_p4_);
	}
	els_pat_genMotherP4_branch = 0;
	if (tree->GetAlias("els_pat_genMotherP4") != 0) {
		els_pat_genMotherP4_branch = tree->GetBranch(tree->GetAlias("els_pat_genMotherP4"));
		els_pat_genMotherP4_branch->SetAddress(&els_pat_genMotherP4_);
	}
	els_pat_genP4_branch = 0;
	if (tree->GetAlias("els_pat_genP4") != 0) {
		els_pat_genP4_branch = tree->GetBranch(tree->GetAlias("els_pat_genP4"));
		els_pat_genP4_branch->SetAddress(&els_pat_genP4_);
	}
	els_pat_p4_branch = 0;
	if (tree->GetAlias("els_pat_p4") != 0) {
		els_pat_p4_branch = tree->GetBranch(tree->GetAlias("els_pat_p4"));
		els_pat_p4_branch->SetAddress(&els_pat_p4_);
	}
	jets_pat_genJet_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genJet_p4") != 0) {
		jets_pat_genJet_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genJet_p4"));
		jets_pat_genJet_p4_branch->SetAddress(&jets_pat_genJet_p4_);
	}
	jets_pat_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genPartonMother_p4") != 0) {
		jets_pat_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genPartonMother_p4"));
		jets_pat_genPartonMother_p4_branch->SetAddress(&jets_pat_genPartonMother_p4_);
	}
	jets_pat_genParton_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genParton_p4") != 0) {
		jets_pat_genParton_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genParton_p4"));
		jets_pat_genParton_p4_branch->SetAddress(&jets_pat_genParton_p4_);
	}
	jets_pat_jet_p4_branch = 0;
	if (tree->GetAlias("jets_pat_jet_p4") != 0) {
		jets_pat_jet_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_jet_p4"));
		jets_pat_jet_p4_branch->SetAddress(&jets_pat_jet_p4_);
	}
	jets_pat_jet_uncorp4_branch = 0;
	if (tree->GetAlias("jets_pat_jet_uncorp4") != 0) {
		jets_pat_jet_uncorp4_branch = tree->GetBranch(tree->GetAlias("jets_pat_jet_uncorp4"));
		jets_pat_jet_uncorp4_branch->SetAddress(&jets_pat_jet_uncorp4_);
	}
	mus_pat_genMotherP4_branch = 0;
	if (tree->GetAlias("mus_pat_genMotherP4") != 0) {
		mus_pat_genMotherP4_branch = tree->GetBranch(tree->GetAlias("mus_pat_genMotherP4"));
		mus_pat_genMotherP4_branch->SetAddress(&mus_pat_genMotherP4_);
	}
	mus_pat_genP4_branch = 0;
	if (tree->GetAlias("mus_pat_genP4") != 0) {
		mus_pat_genP4_branch = tree->GetBranch(tree->GetAlias("mus_pat_genP4"));
		mus_pat_genP4_branch->SetAddress(&mus_pat_genP4_);
	}
	mus_pat_p4_branch = 0;
	if (tree->GetAlias("mus_pat_p4") != 0) {
		mus_pat_p4_branch = tree->GetBranch(tree->GetAlias("mus_pat_p4"));
		mus_pat_p4_branch->SetAddress(&mus_pat_p4_);
	}
	pfels_p4_branch = 0;
	if (tree->GetAlias("pfels_p4") != 0) {
		pfels_p4_branch = tree->GetBranch(tree->GetAlias("pfels_p4"));
		pfels_p4_branch->SetAddress(&pfels_p4_);
	}
	pfels_posAtEcal_p4_branch = 0;
	if (tree->GetAlias("pfels_posAtEcal_p4") != 0) {
		pfels_posAtEcal_p4_branch = tree->GetBranch(tree->GetAlias("pfels_posAtEcal_p4"));
		pfels_posAtEcal_p4_branch->SetAddress(&pfels_posAtEcal_p4_);
	}
	pfjets_p4_branch = 0;
	if (tree->GetAlias("pfjets_p4") != 0) {
		pfjets_p4_branch = tree->GetBranch(tree->GetAlias("pfjets_p4"));
		pfjets_p4_branch->SetAddress(&pfjets_p4_);
	}
	pfmus_p4_branch = 0;
	if (tree->GetAlias("pfmus_p4") != 0) {
		pfmus_p4_branch = tree->GetBranch(tree->GetAlias("pfmus_p4"));
		pfmus_p4_branch->SetAddress(&pfmus_p4_);
	}
	pfmus_posAtEcal_p4_branch = 0;
	if (tree->GetAlias("pfmus_posAtEcal_p4") != 0) {
		pfmus_posAtEcal_p4_branch = tree->GetBranch(tree->GetAlias("pfmus_posAtEcal_p4"));
		pfmus_posAtEcal_p4_branch->SetAddress(&pfmus_posAtEcal_p4_);
	}
	photons_p4_branch = 0;
	if (tree->GetAlias("photons_p4") != 0) {
		photons_p4_branch = tree->GetBranch(tree->GetAlias("photons_p4"));
		photons_p4_branch->SetAddress(&photons_p4_);
	}
	scs_p4_branch = 0;
	if (tree->GetAlias("scs_p4") != 0) {
		scs_p4_branch = tree->GetBranch(tree->GetAlias("scs_p4"));
		scs_p4_branch->SetAddress(&scs_p4_);
	}
	scs_pos_p4_branch = 0;
	if (tree->GetAlias("scs_pos_p4") != 0) {
		scs_pos_p4_branch = tree->GetBranch(tree->GetAlias("scs_pos_p4"));
		scs_pos_p4_branch->SetAddress(&scs_pos_p4_);
	}
	scs_vtx_p4_branch = 0;
	if (tree->GetAlias("scs_vtx_p4") != 0) {
		scs_vtx_p4_branch = tree->GetBranch(tree->GetAlias("scs_vtx_p4"));
		scs_vtx_p4_branch->SetAddress(&scs_vtx_p4_);
	}
	svs_flight_branch = 0;
	if (tree->GetAlias("svs_flight") != 0) {
		svs_flight_branch = tree->GetBranch(tree->GetAlias("svs_flight"));
		svs_flight_branch->SetAddress(&svs_flight_);
	}
	svs_mc3_p4_branch = 0;
	if (tree->GetAlias("svs_mc3_p4") != 0) {
		svs_mc3_p4_branch = tree->GetBranch(tree->GetAlias("svs_mc3_p4"));
		svs_mc3_p4_branch->SetAddress(&svs_mc3_p4_);
	}
	svs_p4_branch = 0;
	if (tree->GetAlias("svs_p4") != 0) {
		svs_p4_branch = tree->GetBranch(tree->GetAlias("svs_p4"));
		svs_p4_branch->SetAddress(&svs_p4_);
	}
	svs_position_branch = 0;
	if (tree->GetAlias("svs_position") != 0) {
		svs_position_branch = tree->GetBranch(tree->GetAlias("svs_position"));
		svs_position_branch->SetAddress(&svs_position_);
	}
	svs_refitp4_branch = 0;
	if (tree->GetAlias("svs_refitp4") != 0) {
		svs_refitp4_branch = tree->GetBranch(tree->GetAlias("svs_refitp4"));
		svs_refitp4_branch->SetAddress(&svs_refitp4_);
	}
	trks_inner_position_branch = 0;
	if (tree->GetAlias("trks_inner_position") != 0) {
		trks_inner_position_branch = tree->GetBranch(tree->GetAlias("trks_inner_position"));
		trks_inner_position_branch->SetAddress(&trks_inner_position_);
	}
	trks_outer_p4_branch = 0;
	if (tree->GetAlias("trks_outer_p4") != 0) {
		trks_outer_p4_branch = tree->GetBranch(tree->GetAlias("trks_outer_p4"));
		trks_outer_p4_branch->SetAddress(&trks_outer_p4_);
	}
	trks_outer_position_branch = 0;
	if (tree->GetAlias("trks_outer_position") != 0) {
		trks_outer_position_branch = tree->GetBranch(tree->GetAlias("trks_outer_position"));
		trks_outer_position_branch->SetAddress(&trks_outer_position_);
	}
	trks_trk_p4_branch = 0;
	if (tree->GetAlias("trks_trk_p4") != 0) {
		trks_trk_p4_branch = tree->GetBranch(tree->GetAlias("trks_trk_p4"));
		trks_trk_p4_branch->SetAddress(&trks_trk_p4_);
	}
	trks_vertex_p4_branch = 0;
	if (tree->GetAlias("trks_vertex_p4") != 0) {
		trks_vertex_p4_branch = tree->GetBranch(tree->GetAlias("trks_vertex_p4"));
		trks_vertex_p4_branch->SetAddress(&trks_vertex_p4_);
	}
	trkjets_p4_branch = 0;
	if (tree->GetAlias("trkjets_p4") != 0) {
		trkjets_p4_branch = tree->GetBranch(tree->GetAlias("trkjets_p4"));
		trkjets_p4_branch->SetAddress(&trkjets_p4_);
	}
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		vtxs_position_branch->SetAddress(&vtxs_position_);
	}
  tree->SetMakeClass(1);
	evt_CMS2tag_branch = 0;
	if (tree->GetAlias("evt_CMS2tag") != 0) {
		evt_CMS2tag_branch = tree->GetBranch(tree->GetAlias("evt_CMS2tag"));
		evt_CMS2tag_branch->SetAddress(&evt_CMS2tag_);
	}
	evt_dataset_branch = 0;
	if (tree->GetAlias("evt_dataset") != 0) {
		evt_dataset_branch = tree->GetBranch(tree->GetAlias("evt_dataset"));
		evt_dataset_branch->SetAddress(&evt_dataset_);
	}
	hlt_trigNames_branch = 0;
	if (tree->GetAlias("hlt_trigNames") != 0) {
		hlt_trigNames_branch = tree->GetBranch(tree->GetAlias("hlt_trigNames"));
		hlt_trigNames_branch->SetAddress(&hlt_trigNames_);
	}
	l1_techtrigNames_branch = 0;
	if (tree->GetAlias("l1_techtrigNames") != 0) {
		l1_techtrigNames_branch = tree->GetBranch(tree->GetAlias("l1_techtrigNames"));
		l1_techtrigNames_branch->SetAddress(&l1_techtrigNames_);
	}
	l1_trigNames_branch = 0;
	if (tree->GetAlias("l1_trigNames") != 0) {
		l1_trigNames_branch = tree->GetBranch(tree->GetAlias("l1_trigNames"));
		l1_trigNames_branch->SetAddress(&l1_trigNames_);
	}
	evt_errCategory_branch = 0;
	if (tree->GetAlias("evt_errCategory") != 0) {
		evt_errCategory_branch = tree->GetBranch(tree->GetAlias("evt_errCategory"));
		evt_errCategory_branch->SetAddress(&evt_errCategory_);
	}
	evt_errModule_branch = 0;
	if (tree->GetAlias("evt_errModule") != 0) {
		evt_errModule_branch = tree->GetBranch(tree->GetAlias("evt_errModule"));
		evt_errModule_branch->SetAddress(&evt_errModule_);
	}
	evt_errSeverity_branch = 0;
	if (tree->GetAlias("evt_errSeverity") != 0) {
		evt_errSeverity_branch = tree->GetBranch(tree->GetAlias("evt_errSeverity"));
		evt_errSeverity_branch->SetAddress(&evt_errSeverity_);
	}
	evt_eventHasHalo_branch = 0;
	if (tree->GetAlias("evt_eventHasHalo") != 0) {
		evt_eventHasHalo_branch = tree->GetBranch(tree->GetAlias("evt_eventHasHalo"));
		evt_eventHasHalo_branch->SetAddress(&evt_eventHasHalo_);
	}
	evt_hbheFilter_branch = 0;
	if (tree->GetAlias("evt_hbheFilter") != 0) {
		evt_hbheFilter_branch = tree->GetBranch(tree->GetAlias("evt_hbheFilter"));
		evt_hbheFilter_branch->SetAddress(&evt_hbheFilter_);
	}
	mus_tightMatch_branch = 0;
	if (tree->GetAlias("mus_tightMatch") != 0) {
		mus_tightMatch_branch = tree->GetBranch(tree->GetAlias("mus_tightMatch"));
		mus_tightMatch_branch->SetAddress(&mus_tightMatch_);
	}
	mus_updatedSta_branch = 0;
	if (tree->GetAlias("mus_updatedSta") != 0) {
		mus_updatedSta_branch = tree->GetBranch(tree->GetAlias("mus_updatedSta"));
		mus_updatedSta_branch->SetAddress(&mus_updatedSta_);
	}
	photons_haspixelSeed_branch = 0;
	if (tree->GetAlias("photons_haspixelSeed") != 0) {
		photons_haspixelSeed_branch = tree->GetBranch(tree->GetAlias("photons_haspixelSeed"));
		photons_haspixelSeed_branch->SetAddress(&photons_haspixelSeed_);
	}
	jets_closestElectron_DR_branch = 0;
	if (tree->GetAlias("jets_closestElectron_DR") != 0) {
		jets_closestElectron_DR_branch = tree->GetBranch(tree->GetAlias("jets_closestElectron_DR"));
		jets_closestElectron_DR_branch->SetAddress(&jets_closestElectron_DR_);
	}
	jets_closestMuon_DR_branch = 0;
	if (tree->GetAlias("jets_closestMuon_DR") != 0) {
		jets_closestMuon_DR_branch = tree->GetBranch(tree->GetAlias("jets_closestMuon_DR"));
		jets_closestMuon_DR_branch->SetAddress(&jets_closestMuon_DR_);
	}
	evt_bs_Xwidth_branch = 0;
	if (tree->GetAlias("evt_bs_Xwidth") != 0) {
		evt_bs_Xwidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Xwidth"));
		evt_bs_Xwidth_branch->SetAddress(&evt_bs_Xwidth_);
	}
	evt_bs_XwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_XwidthErr") != 0) {
		evt_bs_XwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_XwidthErr"));
		evt_bs_XwidthErr_branch->SetAddress(&evt_bs_XwidthErr_);
	}
	evt_bs_Ywidth_branch = 0;
	if (tree->GetAlias("evt_bs_Ywidth") != 0) {
		evt_bs_Ywidth_branch = tree->GetBranch(tree->GetAlias("evt_bs_Ywidth"));
		evt_bs_Ywidth_branch->SetAddress(&evt_bs_Ywidth_);
	}
	evt_bs_YwidthErr_branch = 0;
	if (tree->GetAlias("evt_bs_YwidthErr") != 0) {
		evt_bs_YwidthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_YwidthErr"));
		evt_bs_YwidthErr_branch->SetAddress(&evt_bs_YwidthErr_);
	}
	evt_bs_dxdz_branch = 0;
	if (tree->GetAlias("evt_bs_dxdz") != 0) {
		evt_bs_dxdz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdz"));
		evt_bs_dxdz_branch->SetAddress(&evt_bs_dxdz_);
	}
	evt_bs_dxdzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dxdzErr") != 0) {
		evt_bs_dxdzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dxdzErr"));
		evt_bs_dxdzErr_branch->SetAddress(&evt_bs_dxdzErr_);
	}
	evt_bs_dydz_branch = 0;
	if (tree->GetAlias("evt_bs_dydz") != 0) {
		evt_bs_dydz_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydz"));
		evt_bs_dydz_branch->SetAddress(&evt_bs_dydz_);
	}
	evt_bs_dydzErr_branch = 0;
	if (tree->GetAlias("evt_bs_dydzErr") != 0) {
		evt_bs_dydzErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_dydzErr"));
		evt_bs_dydzErr_branch->SetAddress(&evt_bs_dydzErr_);
	}
	evt_bs_sigmaZ_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZ") != 0) {
		evt_bs_sigmaZ_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZ"));
		evt_bs_sigmaZ_branch->SetAddress(&evt_bs_sigmaZ_);
	}
	evt_bs_sigmaZErr_branch = 0;
	if (tree->GetAlias("evt_bs_sigmaZErr") != 0) {
		evt_bs_sigmaZErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_sigmaZErr"));
		evt_bs_sigmaZErr_branch->SetAddress(&evt_bs_sigmaZErr_);
	}
	evt_bs_xErr_branch = 0;
	if (tree->GetAlias("evt_bs_xErr") != 0) {
		evt_bs_xErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_xErr"));
		evt_bs_xErr_branch->SetAddress(&evt_bs_xErr_);
	}
	evt_bs_yErr_branch = 0;
	if (tree->GetAlias("evt_bs_yErr") != 0) {
		evt_bs_yErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_yErr"));
		evt_bs_yErr_branch->SetAddress(&evt_bs_yErr_);
	}
	evt_bs_zErr_branch = 0;
	if (tree->GetAlias("evt_bs_zErr") != 0) {
		evt_bs_zErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_zErr"));
		evt_bs_zErr_branch->SetAddress(&evt_bs_zErr_);
	}
	evthcal_dmetx_branch = 0;
	if (tree->GetAlias("evthcal_dmetx") != 0) {
		evthcal_dmetx_branch = tree->GetBranch(tree->GetAlias("evthcal_dmetx"));
		evthcal_dmetx_branch->SetAddress(&evthcal_dmetx_);
	}
	evthcal_dmety_branch = 0;
	if (tree->GetAlias("evthcal_dmety") != 0) {
		evthcal_dmety_branch = tree->GetBranch(tree->GetAlias("evthcal_dmety"));
		evthcal_dmety_branch->SetAddress(&evthcal_dmety_);
	}
	evthcal_dsumet_branch = 0;
	if (tree->GetAlias("evthcal_dsumet") != 0) {
		evthcal_dsumet_branch = tree->GetBranch(tree->GetAlias("evthcal_dsumet"));
		evthcal_dsumet_branch->SetAddress(&evthcal_dsumet_);
	}
	evthf_dmetx_branch = 0;
	if (tree->GetAlias("evthf_dmetx") != 0) {
		evthf_dmetx_branch = tree->GetBranch(tree->GetAlias("evthf_dmetx"));
		evthf_dmetx_branch->SetAddress(&evthf_dmetx_);
	}
	evthf_dmety_branch = 0;
	if (tree->GetAlias("evthf_dmety") != 0) {
		evthf_dmety_branch = tree->GetBranch(tree->GetAlias("evthf_dmety"));
		evthf_dmety_branch->SetAddress(&evthf_dmety_);
	}
	evthf_dsumet_branch = 0;
	if (tree->GetAlias("evthf_dsumet") != 0) {
		evthf_dsumet_branch = tree->GetBranch(tree->GetAlias("evthf_dsumet"));
		evthf_dsumet_branch->SetAddress(&evthf_dsumet_);
	}
	evt_bField_branch = 0;
	if (tree->GetAlias("evt_bField") != 0) {
		evt_bField_branch = tree->GetBranch(tree->GetAlias("evt_bField"));
		evt_bField_branch->SetAddress(&evt_bField_);
	}
	evt_kfactor_branch = 0;
	if (tree->GetAlias("evt_kfactor") != 0) {
		evt_kfactor_branch = tree->GetBranch(tree->GetAlias("evt_kfactor"));
		evt_kfactor_branch->SetAddress(&evt_kfactor_);
	}
	evt_scale1fb_branch = 0;
	if (tree->GetAlias("evt_scale1fb") != 0) {
		evt_scale1fb_branch = tree->GetBranch(tree->GetAlias("evt_scale1fb"));
		evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
	}
	evt_xsec_excl_branch = 0;
	if (tree->GetAlias("evt_xsec_excl") != 0) {
		evt_xsec_excl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_excl"));
		evt_xsec_excl_branch->SetAddress(&evt_xsec_excl_);
	}
	evt_xsec_incl_branch = 0;
	if (tree->GetAlias("evt_xsec_incl") != 0) {
		evt_xsec_incl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_incl"));
		evt_xsec_incl_branch->SetAddress(&evt_xsec_incl_);
	}
	gen_met_branch = 0;
	if (tree->GetAlias("gen_met") != 0) {
		gen_met_branch = tree->GetBranch(tree->GetAlias("gen_met"));
		gen_met_branch->SetAddress(&gen_met_);
	}
	gen_metPhi_branch = 0;
	if (tree->GetAlias("gen_metPhi") != 0) {
		gen_metPhi_branch = tree->GetBranch(tree->GetAlias("gen_metPhi"));
		gen_metPhi_branch->SetAddress(&gen_metPhi_);
	}
	genps_alphaQCD_branch = 0;
	if (tree->GetAlias("genps_alphaQCD") != 0) {
		genps_alphaQCD_branch = tree->GetBranch(tree->GetAlias("genps_alphaQCD"));
		genps_alphaQCD_branch->SetAddress(&genps_alphaQCD_);
	}
	genps_pthat_branch = 0;
	if (tree->GetAlias("genps_pthat") != 0) {
		genps_pthat_branch = tree->GetBranch(tree->GetAlias("genps_pthat"));
		genps_pthat_branch->SetAddress(&genps_pthat_);
	}
	genps_qScale_branch = 0;
	if (tree->GetAlias("genps_qScale") != 0) {
		genps_qScale_branch = tree->GetBranch(tree->GetAlias("genps_qScale"));
		genps_qScale_branch->SetAddress(&genps_qScale_);
	}
	genps_weight_branch = 0;
	if (tree->GetAlias("genps_weight") != 0) {
		genps_weight_branch = tree->GetBranch(tree->GetAlias("genps_weight"));
		genps_weight_branch->SetAddress(&genps_weight_);
	}
	gen_sumEt_branch = 0;
	if (tree->GetAlias("gen_sumEt") != 0) {
		gen_sumEt_branch = tree->GetBranch(tree->GetAlias("gen_sumEt"));
		gen_sumEt_branch->SetAddress(&gen_sumEt_);
	}
	hcalnoise_eventChargeFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventChargeFraction") != 0) {
		hcalnoise_eventChargeFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventChargeFraction"));
		hcalnoise_eventChargeFraction_branch->SetAddress(&hcalnoise_eventChargeFraction_);
	}
	hcalnoise_eventEMEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMEnergy") != 0) {
		hcalnoise_eventEMEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMEnergy"));
		hcalnoise_eventEMEnergy_branch->SetAddress(&hcalnoise_eventEMEnergy_);
	}
	hcalnoise_eventEMFraction_branch = 0;
	if (tree->GetAlias("hcalnoise_eventEMFraction") != 0) {
		hcalnoise_eventEMFraction_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventEMFraction"));
		hcalnoise_eventEMFraction_branch->SetAddress(&hcalnoise_eventEMFraction_);
	}
	hcalnoise_eventHadEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventHadEnergy") != 0) {
		hcalnoise_eventHadEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventHadEnergy"));
		hcalnoise_eventHadEnergy_branch->SetAddress(&hcalnoise_eventHadEnergy_);
	}
	hcalnoise_eventTrackEnergy_branch = 0;
	if (tree->GetAlias("hcalnoise_eventTrackEnergy") != 0) {
		hcalnoise_eventTrackEnergy_branch = tree->GetBranch(tree->GetAlias("hcalnoise_eventTrackEnergy"));
		hcalnoise_eventTrackEnergy_branch->SetAddress(&hcalnoise_eventTrackEnergy_);
	}
	hcalnoise_max10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max10GeVHitTime") != 0) {
		hcalnoise_max10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max10GeVHitTime"));
		hcalnoise_max10GeVHitTime_branch->SetAddress(&hcalnoise_max10GeVHitTime_);
	}
	hcalnoise_max25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_max25GeVHitTime") != 0) {
		hcalnoise_max25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_max25GeVHitTime"));
		hcalnoise_max25GeVHitTime_branch->SetAddress(&hcalnoise_max25GeVHitTime_);
	}
	hcalnoise_min10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min10GeVHitTime") != 0) {
		hcalnoise_min10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min10GeVHitTime"));
		hcalnoise_min10GeVHitTime_branch->SetAddress(&hcalnoise_min10GeVHitTime_);
	}
	hcalnoise_min25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_min25GeVHitTime") != 0) {
		hcalnoise_min25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_min25GeVHitTime"));
		hcalnoise_min25GeVHitTime_branch->SetAddress(&hcalnoise_min25GeVHitTime_);
	}
	hcalnoise_minE10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE10TS") != 0) {
		hcalnoise_minE10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE10TS"));
		hcalnoise_minE10TS_branch->SetAddress(&hcalnoise_minE10TS_);
	}
	hcalnoise_minE2Over10TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2Over10TS") != 0) {
		hcalnoise_minE2Over10TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2Over10TS"));
		hcalnoise_minE2Over10TS_branch->SetAddress(&hcalnoise_minE2Over10TS_);
	}
	hcalnoise_minE2TS_branch = 0;
	if (tree->GetAlias("hcalnoise_minE2TS") != 0) {
		hcalnoise_minE2TS_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minE2TS"));
		hcalnoise_minE2TS_branch->SetAddress(&hcalnoise_minE2TS_);
	}
	hcalnoise_minHPDEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minHPDEMF") != 0) {
		hcalnoise_minHPDEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minHPDEMF"));
		hcalnoise_minHPDEMF_branch->SetAddress(&hcalnoise_minHPDEMF_);
	}
	hcalnoise_minRBXEMF_branch = 0;
	if (tree->GetAlias("hcalnoise_minRBXEMF") != 0) {
		hcalnoise_minRBXEMF_branch = tree->GetBranch(tree->GetAlias("hcalnoise_minRBXEMF"));
		hcalnoise_minRBXEMF_branch->SetAddress(&hcalnoise_minRBXEMF_);
	}
	hcalnoise_rms10GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms10GeVHitTime") != 0) {
		hcalnoise_rms10GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms10GeVHitTime"));
		hcalnoise_rms10GeVHitTime_branch->SetAddress(&hcalnoise_rms10GeVHitTime_);
	}
	hcalnoise_rms25GeVHitTime_branch = 0;
	if (tree->GetAlias("hcalnoise_rms25GeVHitTime") != 0) {
		hcalnoise_rms25GeVHitTime_branch = tree->GetBranch(tree->GetAlias("hcalnoise_rms25GeVHitTime"));
		hcalnoise_rms25GeVHitTime_branch->SetAddress(&hcalnoise_rms25GeVHitTime_);
	}
	l1_met_etTot_branch = 0;
	if (tree->GetAlias("l1_met_etTot") != 0) {
		l1_met_etTot_branch = tree->GetBranch(tree->GetAlias("l1_met_etTot"));
		l1_met_etTot_branch->SetAddress(&l1_met_etTot_);
	}
	l1_met_met_branch = 0;
	if (tree->GetAlias("l1_met_met") != 0) {
		l1_met_met_branch = tree->GetBranch(tree->GetAlias("l1_met_met"));
		l1_met_met_branch->SetAddress(&l1_met_met_);
	}
	l1_mht_htTot_branch = 0;
	if (tree->GetAlias("l1_mht_htTot") != 0) {
		l1_mht_htTot_branch = tree->GetBranch(tree->GetAlias("l1_mht_htTot"));
		l1_mht_htTot_branch->SetAddress(&l1_mht_htTot_);
	}
	l1_mht_mht_branch = 0;
	if (tree->GetAlias("l1_mht_mht") != 0) {
		l1_mht_mht_branch = tree->GetBranch(tree->GetAlias("l1_mht_mht"));
		l1_mht_mht_branch->SetAddress(&l1_mht_mht_);
	}
	evt_ecalendcapm_met_branch = 0;
	if (tree->GetAlias("evt_ecalendcapm_met") != 0) {
		evt_ecalendcapm_met_branch = tree->GetBranch(tree->GetAlias("evt_ecalendcapm_met"));
		evt_ecalendcapm_met_branch->SetAddress(&evt_ecalendcapm_met_);
	}
	evt_ecalendcapm_metPhi_branch = 0;
	if (tree->GetAlias("evt_ecalendcapm_metPhi") != 0) {
		evt_ecalendcapm_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_ecalendcapm_metPhi"));
		evt_ecalendcapm_metPhi_branch->SetAddress(&evt_ecalendcapm_metPhi_);
	}
	evt_ecalendcapp_met_branch = 0;
	if (tree->GetAlias("evt_ecalendcapp_met") != 0) {
		evt_ecalendcapp_met_branch = tree->GetBranch(tree->GetAlias("evt_ecalendcapp_met"));
		evt_ecalendcapp_met_branch->SetAddress(&evt_ecalendcapp_met_);
	}
	evt_ecalendcapp_metPhi_branch = 0;
	if (tree->GetAlias("evt_ecalendcapp_metPhi") != 0) {
		evt_ecalendcapp_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_ecalendcapp_metPhi"));
		evt_ecalendcapp_metPhi_branch->SetAddress(&evt_ecalendcapp_metPhi_);
	}
	evt_ecalmet_branch = 0;
	if (tree->GetAlias("evt_ecalmet") != 0) {
		evt_ecalmet_branch = tree->GetBranch(tree->GetAlias("evt_ecalmet"));
		evt_ecalmet_branch->SetAddress(&evt_ecalmet_);
	}
	evt_ecalmetPhi_branch = 0;
	if (tree->GetAlias("evt_ecalmetPhi") != 0) {
		evt_ecalmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_ecalmetPhi"));
		evt_ecalmetPhi_branch->SetAddress(&evt_ecalmetPhi_);
	}
	evt_endcapm_met_branch = 0;
	if (tree->GetAlias("evt_endcapm_met") != 0) {
		evt_endcapm_met_branch = tree->GetBranch(tree->GetAlias("evt_endcapm_met"));
		evt_endcapm_met_branch->SetAddress(&evt_endcapm_met_);
	}
	evt_endcapm_metPhi_branch = 0;
	if (tree->GetAlias("evt_endcapm_metPhi") != 0) {
		evt_endcapm_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_endcapm_metPhi"));
		evt_endcapm_metPhi_branch->SetAddress(&evt_endcapm_metPhi_);
	}
	evt_endcapp_met_branch = 0;
	if (tree->GetAlias("evt_endcapp_met") != 0) {
		evt_endcapp_met_branch = tree->GetBranch(tree->GetAlias("evt_endcapp_met"));
		evt_endcapp_met_branch->SetAddress(&evt_endcapp_met_);
	}
	evt_endcapp_metPhi_branch = 0;
	if (tree->GetAlias("evt_endcapp_metPhi") != 0) {
		evt_endcapp_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_endcapp_metPhi"));
		evt_endcapp_metPhi_branch->SetAddress(&evt_endcapp_metPhi_);
	}
	evt_hcalendcapm_met_branch = 0;
	if (tree->GetAlias("evt_hcalendcapm_met") != 0) {
		evt_hcalendcapm_met_branch = tree->GetBranch(tree->GetAlias("evt_hcalendcapm_met"));
		evt_hcalendcapm_met_branch->SetAddress(&evt_hcalendcapm_met_);
	}
	evt_hcalendcapm_metPhi_branch = 0;
	if (tree->GetAlias("evt_hcalendcapm_metPhi") != 0) {
		evt_hcalendcapm_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_hcalendcapm_metPhi"));
		evt_hcalendcapm_metPhi_branch->SetAddress(&evt_hcalendcapm_metPhi_);
	}
	evt_hcalendcapp_met_branch = 0;
	if (tree->GetAlias("evt_hcalendcapp_met") != 0) {
		evt_hcalendcapp_met_branch = tree->GetBranch(tree->GetAlias("evt_hcalendcapp_met"));
		evt_hcalendcapp_met_branch->SetAddress(&evt_hcalendcapp_met_);
	}
	evt_hcalendcapp_metPhi_branch = 0;
	if (tree->GetAlias("evt_hcalendcapp_metPhi") != 0) {
		evt_hcalendcapp_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_hcalendcapp_metPhi"));
		evt_hcalendcapp_metPhi_branch->SetAddress(&evt_hcalendcapp_metPhi_);
	}
	evt_hcalmet_branch = 0;
	if (tree->GetAlias("evt_hcalmet") != 0) {
		evt_hcalmet_branch = tree->GetBranch(tree->GetAlias("evt_hcalmet"));
		evt_hcalmet_branch->SetAddress(&evt_hcalmet_);
	}
	evt_hcalmetPhi_branch = 0;
	if (tree->GetAlias("evt_hcalmetPhi") != 0) {
		evt_hcalmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_hcalmetPhi"));
		evt_hcalmetPhi_branch->SetAddress(&evt_hcalmetPhi_);
	}
	evt_met_branch = 0;
	if (tree->GetAlias("evt_met") != 0) {
		evt_met_branch = tree->GetBranch(tree->GetAlias("evt_met"));
		evt_met_branch->SetAddress(&evt_met_);
	}
	evt_metHO_branch = 0;
	if (tree->GetAlias("evt_metHO") != 0) {
		evt_metHO_branch = tree->GetBranch(tree->GetAlias("evt_metHO"));
		evt_metHO_branch->SetAddress(&evt_metHO_);
	}
	evt_metHOPhi_branch = 0;
	if (tree->GetAlias("evt_metHOPhi") != 0) {
		evt_metHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metHOPhi"));
		evt_metHOPhi_branch->SetAddress(&evt_metHOPhi_);
	}
	evt_metHOSig_branch = 0;
	if (tree->GetAlias("evt_metHOSig") != 0) {
		evt_metHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metHOSig"));
		evt_metHOSig_branch->SetAddress(&evt_metHOSig_);
	}
	evt_metMuonCorr_branch = 0;
	if (tree->GetAlias("evt_metMuonCorr") != 0) {
		evt_metMuonCorr_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorr"));
		evt_metMuonCorr_branch->SetAddress(&evt_metMuonCorr_);
	}
	evt_metMuonCorrPhi_branch = 0;
	if (tree->GetAlias("evt_metMuonCorrPhi") != 0) {
		evt_metMuonCorrPhi_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorrPhi"));
		evt_metMuonCorrPhi_branch->SetAddress(&evt_metMuonCorrPhi_);
	}
	evt_metMuonCorrSig_branch = 0;
	if (tree->GetAlias("evt_metMuonCorrSig") != 0) {
		evt_metMuonCorrSig_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorrSig"));
		evt_metMuonCorrSig_branch->SetAddress(&evt_metMuonCorrSig_);
	}
	evt_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorr") != 0) {
		evt_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorr"));
		evt_metMuonJESCorr_branch->SetAddress(&evt_metMuonJESCorr_);
	}
	evt_metMuonJESCorrPhi_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorrPhi") != 0) {
		evt_metMuonJESCorrPhi_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorrPhi"));
		evt_metMuonJESCorrPhi_branch->SetAddress(&evt_metMuonJESCorrPhi_);
	}
	evt_metMuonJESCorrSig_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorrSig") != 0) {
		evt_metMuonJESCorrSig_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorrSig"));
		evt_metMuonJESCorrSig_branch->SetAddress(&evt_metMuonJESCorrSig_);
	}
	evt_metNoHF_branch = 0;
	if (tree->GetAlias("evt_metNoHF") != 0) {
		evt_metNoHF_branch = tree->GetBranch(tree->GetAlias("evt_metNoHF"));
		evt_metNoHF_branch->SetAddress(&evt_metNoHF_);
	}
	evt_metNoHFHO_branch = 0;
	if (tree->GetAlias("evt_metNoHFHO") != 0) {
		evt_metNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHO"));
		evt_metNoHFHO_branch->SetAddress(&evt_metNoHFHO_);
	}
	evt_metNoHFHOPhi_branch = 0;
	if (tree->GetAlias("evt_metNoHFHOPhi") != 0) {
		evt_metNoHFHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHOPhi"));
		evt_metNoHFHOPhi_branch->SetAddress(&evt_metNoHFHOPhi_);
	}
	evt_metNoHFHOSig_branch = 0;
	if (tree->GetAlias("evt_metNoHFHOSig") != 0) {
		evt_metNoHFHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHOSig"));
		evt_metNoHFHOSig_branch->SetAddress(&evt_metNoHFHOSig_);
	}
	evt_metNoHFPhi_branch = 0;
	if (tree->GetAlias("evt_metNoHFPhi") != 0) {
		evt_metNoHFPhi_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFPhi"));
		evt_metNoHFPhi_branch->SetAddress(&evt_metNoHFPhi_);
	}
	evt_metNoHFSig_branch = 0;
	if (tree->GetAlias("evt_metNoHFSig") != 0) {
		evt_metNoHFSig_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFSig"));
		evt_metNoHFSig_branch->SetAddress(&evt_metNoHFSig_);
	}
	evt_metOpt_branch = 0;
	if (tree->GetAlias("evt_metOpt") != 0) {
		evt_metOpt_branch = tree->GetBranch(tree->GetAlias("evt_metOpt"));
		evt_metOpt_branch->SetAddress(&evt_metOpt_);
	}
	evt_metOptHO_branch = 0;
	if (tree->GetAlias("evt_metOptHO") != 0) {
		evt_metOptHO_branch = tree->GetBranch(tree->GetAlias("evt_metOptHO"));
		evt_metOptHO_branch->SetAddress(&evt_metOptHO_);
	}
	evt_metOptHOPhi_branch = 0;
	if (tree->GetAlias("evt_metOptHOPhi") != 0) {
		evt_metOptHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptHOPhi"));
		evt_metOptHOPhi_branch->SetAddress(&evt_metOptHOPhi_);
	}
	evt_metOptHOSig_branch = 0;
	if (tree->GetAlias("evt_metOptHOSig") != 0) {
		evt_metOptHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptHOSig"));
		evt_metOptHOSig_branch->SetAddress(&evt_metOptHOSig_);
	}
	evt_metOptNoHF_branch = 0;
	if (tree->GetAlias("evt_metOptNoHF") != 0) {
		evt_metOptNoHF_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHF"));
		evt_metOptNoHF_branch->SetAddress(&evt_metOptNoHF_);
	}
	evt_metOptNoHFHO_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHO") != 0) {
		evt_metOptNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHO"));
		evt_metOptNoHFHO_branch->SetAddress(&evt_metOptNoHFHO_);
	}
	evt_metOptNoHFHOPhi_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHOPhi") != 0) {
		evt_metOptNoHFHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHOPhi"));
		evt_metOptNoHFHOPhi_branch->SetAddress(&evt_metOptNoHFHOPhi_);
	}
	evt_metOptNoHFHOSig_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHOSig") != 0) {
		evt_metOptNoHFHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHOSig"));
		evt_metOptNoHFHOSig_branch->SetAddress(&evt_metOptNoHFHOSig_);
	}
	evt_metOptNoHFPhi_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFPhi") != 0) {
		evt_metOptNoHFPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFPhi"));
		evt_metOptNoHFPhi_branch->SetAddress(&evt_metOptNoHFPhi_);
	}
	evt_metOptNoHFSig_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFSig") != 0) {
		evt_metOptNoHFSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFSig"));
		evt_metOptNoHFSig_branch->SetAddress(&evt_metOptNoHFSig_);
	}
	evt_metOptPhi_branch = 0;
	if (tree->GetAlias("evt_metOptPhi") != 0) {
		evt_metOptPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptPhi"));
		evt_metOptPhi_branch->SetAddress(&evt_metOptPhi_);
	}
	evt_metOptSig_branch = 0;
	if (tree->GetAlias("evt_metOptSig") != 0) {
		evt_metOptSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptSig"));
		evt_metOptSig_branch->SetAddress(&evt_metOptSig_);
	}
	evt_metPhi_branch = 0;
	if (tree->GetAlias("evt_metPhi") != 0) {
		evt_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_metPhi"));
		evt_metPhi_branch->SetAddress(&evt_metPhi_);
	}
	evt_metSig_branch = 0;
	if (tree->GetAlias("evt_metSig") != 0) {
		evt_metSig_branch = tree->GetBranch(tree->GetAlias("evt_metSig"));
		evt_metSig_branch->SetAddress(&evt_metSig_);
	}
	evt_sumet_branch = 0;
	if (tree->GetAlias("evt_sumet") != 0) {
		evt_sumet_branch = tree->GetBranch(tree->GetAlias("evt_sumet"));
		evt_sumet_branch->SetAddress(&evt_sumet_);
	}
	evt_sumetHO_branch = 0;
	if (tree->GetAlias("evt_sumetHO") != 0) {
		evt_sumetHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetHO"));
		evt_sumetHO_branch->SetAddress(&evt_sumetHO_);
	}
	evt_sumetMuonCorr_branch = 0;
	if (tree->GetAlias("evt_sumetMuonCorr") != 0) {
		evt_sumetMuonCorr_branch = tree->GetBranch(tree->GetAlias("evt_sumetMuonCorr"));
		evt_sumetMuonCorr_branch->SetAddress(&evt_sumetMuonCorr_);
	}
	evt_sumetNoHF_branch = 0;
	if (tree->GetAlias("evt_sumetNoHF") != 0) {
		evt_sumetNoHF_branch = tree->GetBranch(tree->GetAlias("evt_sumetNoHF"));
		evt_sumetNoHF_branch->SetAddress(&evt_sumetNoHF_);
	}
	evt_sumetNoHFHO_branch = 0;
	if (tree->GetAlias("evt_sumetNoHFHO") != 0) {
		evt_sumetNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetNoHFHO"));
		evt_sumetNoHFHO_branch->SetAddress(&evt_sumetNoHFHO_);
	}
	evt_sumetOpt_branch = 0;
	if (tree->GetAlias("evt_sumetOpt") != 0) {
		evt_sumetOpt_branch = tree->GetBranch(tree->GetAlias("evt_sumetOpt"));
		evt_sumetOpt_branch->SetAddress(&evt_sumetOpt_);
	}
	evt_sumetOptHO_branch = 0;
	if (tree->GetAlias("evt_sumetOptHO") != 0) {
		evt_sumetOptHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptHO"));
		evt_sumetOptHO_branch->SetAddress(&evt_sumetOptHO_);
	}
	evt_sumetOptNoHF_branch = 0;
	if (tree->GetAlias("evt_sumetOptNoHF") != 0) {
		evt_sumetOptNoHF_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptNoHF"));
		evt_sumetOptNoHF_branch->SetAddress(&evt_sumetOptNoHF_);
	}
	evt_sumetOptNoHFHO_branch = 0;
	if (tree->GetAlias("evt_sumetOptNoHFHO") != 0) {
		evt_sumetOptNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptNoHFHO"));
		evt_sumetOptNoHFHO_branch->SetAddress(&evt_sumetOptNoHFHO_);
	}
	met_pat_metCor_branch = 0;
	if (tree->GetAlias("met_pat_metCor") != 0) {
		met_pat_metCor_branch = tree->GetBranch(tree->GetAlias("met_pat_metCor"));
		met_pat_metCor_branch->SetAddress(&met_pat_metCor_);
	}
	met_pat_metPhiCor_branch = 0;
	if (tree->GetAlias("met_pat_metPhiCor") != 0) {
		met_pat_metPhiCor_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiCor"));
		met_pat_metPhiCor_branch->SetAddress(&met_pat_metPhiCor_);
	}
	met_pat_metPhiUncor_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncor") != 0) {
		met_pat_metPhiUncor_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncor"));
		met_pat_metPhiUncor_branch->SetAddress(&met_pat_metPhiUncor_);
	}
	met_pat_metPhiUncorJES_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncorJES") != 0) {
		met_pat_metPhiUncorJES_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncorJES"));
		met_pat_metPhiUncorJES_branch->SetAddress(&met_pat_metPhiUncorJES_);
	}
	met_pat_metPhiUncorMuon_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncorMuon") != 0) {
		met_pat_metPhiUncorMuon_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncorMuon"));
		met_pat_metPhiUncorMuon_branch->SetAddress(&met_pat_metPhiUncorMuon_);
	}
	met_pat_metUncor_branch = 0;
	if (tree->GetAlias("met_pat_metUncor") != 0) {
		met_pat_metUncor_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncor"));
		met_pat_metUncor_branch->SetAddress(&met_pat_metUncor_);
	}
	met_pat_metUncorJES_branch = 0;
	if (tree->GetAlias("met_pat_metUncorJES") != 0) {
		met_pat_metUncorJES_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncorJES"));
		met_pat_metUncorJES_branch->SetAddress(&met_pat_metUncorJES_);
	}
	met_pat_metUncorMuon_branch = 0;
	if (tree->GetAlias("met_pat_metUncorMuon") != 0) {
		met_pat_metUncorMuon_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncorMuon"));
		met_pat_metUncorMuon_branch->SetAddress(&met_pat_metUncorMuon_);
	}
	pdfinfo_scale_branch = 0;
	if (tree->GetAlias("pdfinfo_scale") != 0) {
		pdfinfo_scale_branch = tree->GetBranch(tree->GetAlias("pdfinfo_scale"));
		pdfinfo_scale_branch->SetAddress(&pdfinfo_scale_);
	}
	pdfinfo_x1_branch = 0;
	if (tree->GetAlias("pdfinfo_x1") != 0) {
		pdfinfo_x1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x1"));
		pdfinfo_x1_branch->SetAddress(&pdfinfo_x1_);
	}
	pdfinfo_x2_branch = 0;
	if (tree->GetAlias("pdfinfo_x2") != 0) {
		pdfinfo_x2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x2"));
		pdfinfo_x2_branch->SetAddress(&pdfinfo_x2_);
	}
	evt_pfmet_branch = 0;
	if (tree->GetAlias("evt_pfmet") != 0) {
		evt_pfmet_branch = tree->GetBranch(tree->GetAlias("evt_pfmet"));
		evt_pfmet_branch->SetAddress(&evt_pfmet_);
	}
	evt_pfmetPhi_branch = 0;
	if (tree->GetAlias("evt_pfmetPhi") != 0) {
		evt_pfmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_pfmetPhi"));
		evt_pfmetPhi_branch->SetAddress(&evt_pfmetPhi_);
	}
	evt_pfmetSig_branch = 0;
	if (tree->GetAlias("evt_pfmetSig") != 0) {
		evt_pfmetSig_branch = tree->GetBranch(tree->GetAlias("evt_pfmetSig"));
		evt_pfmetSig_branch->SetAddress(&evt_pfmetSig_);
	}
	evt_pfsumet_branch = 0;
	if (tree->GetAlias("evt_pfsumet") != 0) {
		evt_pfsumet_branch = tree->GetBranch(tree->GetAlias("evt_pfsumet"));
		evt_pfsumet_branch->SetAddress(&evt_pfsumet_);
	}
	evt_tcmet_branch = 0;
	if (tree->GetAlias("evt_tcmet") != 0) {
		evt_tcmet_branch = tree->GetBranch(tree->GetAlias("evt_tcmet"));
		evt_tcmet_branch->SetAddress(&evt_tcmet_);
	}
	evt_tcmetPhi_branch = 0;
	if (tree->GetAlias("evt_tcmetPhi") != 0) {
		evt_tcmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_tcmetPhi"));
		evt_tcmetPhi_branch->SetAddress(&evt_tcmetPhi_);
	}
	evt_tcmetSig_branch = 0;
	if (tree->GetAlias("evt_tcmetSig") != 0) {
		evt_tcmetSig_branch = tree->GetBranch(tree->GetAlias("evt_tcmetSig"));
		evt_tcmetSig_branch->SetAddress(&evt_tcmetSig_);
	}
	evt_tcsumet_branch = 0;
	if (tree->GetAlias("evt_tcsumet") != 0) {
		evt_tcsumet_branch = tree->GetBranch(tree->GetAlias("evt_tcsumet"));
		evt_tcsumet_branch->SetAddress(&evt_tcsumet_);
	}
	genps_lepdaughter_p4_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_p4") != 0) {
		genps_lepdaughter_p4_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_p4"));
		genps_lepdaughter_p4_branch->SetAddress(&genps_lepdaughter_p4_);
	}
	hlt_trigObjs_p4_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_p4") != 0) {
		hlt_trigObjs_p4_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_p4"));
		hlt_trigObjs_p4_branch->SetAddress(&hlt_trigObjs_p4_);
	}
	hyp_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_p4") != 0) {
		hyp_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_p4"));
		hyp_jets_p4_branch->SetAddress(&hyp_jets_p4_);
	}
	hyp_other_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_p4") != 0) {
		hyp_other_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_p4"));
		hyp_other_jets_p4_branch->SetAddress(&hyp_other_jets_p4_);
	}
	jpts_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jpts_combinedSecondaryVertexBJetTag") != 0) {
		jpts_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_combinedSecondaryVertexBJetTag"));
		jpts_combinedSecondaryVertexBJetTag_branch->SetAddress(&jpts_combinedSecondaryVertexBJetTag_);
	}
	jpts_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("jpts_combinedSecondaryVertexMVABJetTag") != 0) {
		jpts_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_combinedSecondaryVertexMVABJetTag"));
		jpts_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&jpts_combinedSecondaryVertexMVABJetTag_);
	}
	jpts_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jpts_jetBProbabilityBJetTag") != 0) {
		jpts_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_jetBProbabilityBJetTag"));
		jpts_jetBProbabilityBJetTag_branch->SetAddress(&jpts_jetBProbabilityBJetTag_);
	}
	jpts_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jpts_jetProbabilityBJetTag") != 0) {
		jpts_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_jetProbabilityBJetTag"));
		jpts_jetProbabilityBJetTag_branch->SetAddress(&jpts_jetProbabilityBJetTag_);
	}
	jpts_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jpts_simpleSecondaryVertexHighEffBJetTag") != 0) {
		jpts_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_simpleSecondaryVertexHighEffBJetTag"));
		jpts_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&jpts_simpleSecondaryVertexHighEffBJetTag_);
	}
	jpts_simpleSecondaryVertexHighPurBJetTags_branch = 0;
	if (tree->GetAlias("jpts_simpleSecondaryVertexHighPurBJetTags") != 0) {
		jpts_simpleSecondaryVertexHighPurBJetTags_branch = tree->GetBranch(tree->GetAlias("jpts_simpleSecondaryVertexHighPurBJetTags"));
		jpts_simpleSecondaryVertexHighPurBJetTags_branch->SetAddress(&jpts_simpleSecondaryVertexHighPurBJetTags_);
	}
	jpts_softElectronByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jpts_softElectronByIP3dBJetTag") != 0) {
		jpts_softElectronByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_softElectronByIP3dBJetTag"));
		jpts_softElectronByIP3dBJetTag_branch->SetAddress(&jpts_softElectronByIP3dBJetTag_);
	}
	jpts_softElectronByPtBJetTag_branch = 0;
	if (tree->GetAlias("jpts_softElectronByPtBJetTag") != 0) {
		jpts_softElectronByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_softElectronByPtBJetTag"));
		jpts_softElectronByPtBJetTag_branch->SetAddress(&jpts_softElectronByPtBJetTag_);
	}
	jpts_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("jpts_softMuonBJetTag") != 0) {
		jpts_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_softMuonBJetTag"));
		jpts_softMuonBJetTag_branch->SetAddress(&jpts_softMuonBJetTag_);
	}
	jpts_softMuonByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jpts_softMuonByIP3dBJetTag") != 0) {
		jpts_softMuonByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_softMuonByIP3dBJetTag"));
		jpts_softMuonByIP3dBJetTag_branch->SetAddress(&jpts_softMuonByIP3dBJetTag_);
	}
	jpts_softMuonByPtBJetTag_branch = 0;
	if (tree->GetAlias("jpts_softMuonByPtBJetTag") != 0) {
		jpts_softMuonByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_softMuonByPtBJetTag"));
		jpts_softMuonByPtBJetTag_branch->SetAddress(&jpts_softMuonByPtBJetTag_);
	}
	jpts_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jpts_trackCountingHighEffBJetTag") != 0) {
		jpts_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_trackCountingHighEffBJetTag"));
		jpts_trackCountingHighEffBJetTag_branch->SetAddress(&jpts_trackCountingHighEffBJetTag_);
	}
	jpts_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jpts_trackCountingHighPurBJetTag") != 0) {
		jpts_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jpts_trackCountingHighPurBJetTag"));
		jpts_trackCountingHighPurBJetTag_branch->SetAddress(&jpts_trackCountingHighPurBJetTag_);
	}
	jets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_combinedSecondaryVertexBJetTag") != 0) {
		jets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_combinedSecondaryVertexBJetTag"));
		jets_combinedSecondaryVertexBJetTag_branch->SetAddress(&jets_combinedSecondaryVertexBJetTag_);
	}
	jets_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_combinedSecondaryVertexMVABJetTag") != 0) {
		jets_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_combinedSecondaryVertexMVABJetTag"));
		jets_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&jets_combinedSecondaryVertexMVABJetTag_);
	}
	jets_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_jetBProbabilityBJetTag") != 0) {
		jets_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_jetBProbabilityBJetTag"));
		jets_jetBProbabilityBJetTag_branch->SetAddress(&jets_jetBProbabilityBJetTag_);
	}
	jets_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_jetProbabilityBJetTag") != 0) {
		jets_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_jetProbabilityBJetTag"));
		jets_jetProbabilityBJetTag_branch->SetAddress(&jets_jetProbabilityBJetTag_);
	}
	jets_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_simpleSecondaryVertexHighEffBJetTag") != 0) {
		jets_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_simpleSecondaryVertexHighEffBJetTag"));
		jets_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&jets_simpleSecondaryVertexHighEffBJetTag_);
	}
	jets_simpleSecondaryVertexHighPurBJetTags_branch = 0;
	if (tree->GetAlias("jets_simpleSecondaryVertexHighPurBJetTags") != 0) {
		jets_simpleSecondaryVertexHighPurBJetTags_branch = tree->GetBranch(tree->GetAlias("jets_simpleSecondaryVertexHighPurBJetTags"));
		jets_simpleSecondaryVertexHighPurBJetTags_branch->SetAddress(&jets_simpleSecondaryVertexHighPurBJetTags_);
	}
	jets_softElectronByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jets_softElectronByIP3dBJetTag") != 0) {
		jets_softElectronByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softElectronByIP3dBJetTag"));
		jets_softElectronByIP3dBJetTag_branch->SetAddress(&jets_softElectronByIP3dBJetTag_);
	}
	jets_softElectronByPtBJetTag_branch = 0;
	if (tree->GetAlias("jets_softElectronByPtBJetTag") != 0) {
		jets_softElectronByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softElectronByPtBJetTag"));
		jets_softElectronByPtBJetTag_branch->SetAddress(&jets_softElectronByPtBJetTag_);
	}
	jets_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("jets_softMuonBJetTag") != 0) {
		jets_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softMuonBJetTag"));
		jets_softMuonBJetTag_branch->SetAddress(&jets_softMuonBJetTag_);
	}
	jets_softMuonByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jets_softMuonByIP3dBJetTag") != 0) {
		jets_softMuonByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softMuonByIP3dBJetTag"));
		jets_softMuonByIP3dBJetTag_branch->SetAddress(&jets_softMuonByIP3dBJetTag_);
	}
	jets_softMuonByPtBJetTag_branch = 0;
	if (tree->GetAlias("jets_softMuonByPtBJetTag") != 0) {
		jets_softMuonByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softMuonByPtBJetTag"));
		jets_softMuonByPtBJetTag_branch->SetAddress(&jets_softMuonByPtBJetTag_);
	}
	jets_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_trackCountingHighEffBJetTag") != 0) {
		jets_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_trackCountingHighEffBJetTag"));
		jets_trackCountingHighEffBJetTag_branch->SetAddress(&jets_trackCountingHighEffBJetTag_);
	}
	jets_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jets_trackCountingHighPurBJetTag") != 0) {
		jets_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_trackCountingHighPurBJetTag"));
		jets_trackCountingHighPurBJetTag_branch->SetAddress(&jets_trackCountingHighPurBJetTag_);
	}
	pfjets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_combinedSecondaryVertexBJetTag") != 0) {
		pfjets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_combinedSecondaryVertexBJetTag"));
		pfjets_combinedSecondaryVertexBJetTag_branch->SetAddress(&pfjets_combinedSecondaryVertexBJetTag_);
	}
	pfjets_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("pfjets_combinedSecondaryVertexMVABJetTag") != 0) {
		pfjets_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_combinedSecondaryVertexMVABJetTag"));
		pfjets_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&pfjets_combinedSecondaryVertexMVABJetTag_);
	}
	pfjets_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_jetBProbabilityBJetTag") != 0) {
		pfjets_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_jetBProbabilityBJetTag"));
		pfjets_jetBProbabilityBJetTag_branch->SetAddress(&pfjets_jetBProbabilityBJetTag_);
	}
	pfjets_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_jetProbabilityBJetTag") != 0) {
		pfjets_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_jetProbabilityBJetTag"));
		pfjets_jetProbabilityBJetTag_branch->SetAddress(&pfjets_jetProbabilityBJetTag_);
	}
	pfjets_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_simpleSecondaryVertexHighEffBJetTag") != 0) {
		pfjets_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_simpleSecondaryVertexHighEffBJetTag"));
		pfjets_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&pfjets_simpleSecondaryVertexHighEffBJetTag_);
	}
	pfjets_simpleSecondaryVertexHighPurBJetTags_branch = 0;
	if (tree->GetAlias("pfjets_simpleSecondaryVertexHighPurBJetTags") != 0) {
		pfjets_simpleSecondaryVertexHighPurBJetTags_branch = tree->GetBranch(tree->GetAlias("pfjets_simpleSecondaryVertexHighPurBJetTags"));
		pfjets_simpleSecondaryVertexHighPurBJetTags_branch->SetAddress(&pfjets_simpleSecondaryVertexHighPurBJetTags_);
	}
	pfjets_softElectronByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_softElectronByIP3dBJetTag") != 0) {
		pfjets_softElectronByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_softElectronByIP3dBJetTag"));
		pfjets_softElectronByIP3dBJetTag_branch->SetAddress(&pfjets_softElectronByIP3dBJetTag_);
	}
	pfjets_softElectronByPtBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_softElectronByPtBJetTag") != 0) {
		pfjets_softElectronByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_softElectronByPtBJetTag"));
		pfjets_softElectronByPtBJetTag_branch->SetAddress(&pfjets_softElectronByPtBJetTag_);
	}
	pfjets_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_softMuonBJetTag") != 0) {
		pfjets_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_softMuonBJetTag"));
		pfjets_softMuonBJetTag_branch->SetAddress(&pfjets_softMuonBJetTag_);
	}
	pfjets_softMuonByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_softMuonByIP3dBJetTag") != 0) {
		pfjets_softMuonByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_softMuonByIP3dBJetTag"));
		pfjets_softMuonByIP3dBJetTag_branch->SetAddress(&pfjets_softMuonByIP3dBJetTag_);
	}
	pfjets_softMuonByPtBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_softMuonByPtBJetTag") != 0) {
		pfjets_softMuonByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_softMuonByPtBJetTag"));
		pfjets_softMuonByPtBJetTag_branch->SetAddress(&pfjets_softMuonByPtBJetTag_);
	}
	pfjets_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_trackCountingHighEffBJetTag") != 0) {
		pfjets_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_trackCountingHighEffBJetTag"));
		pfjets_trackCountingHighEffBJetTag_branch->SetAddress(&pfjets_trackCountingHighEffBJetTag_);
	}
	pfjets_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("pfjets_trackCountingHighPurBJetTag") != 0) {
		pfjets_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("pfjets_trackCountingHighPurBJetTag"));
		pfjets_trackCountingHighPurBJetTag_branch->SetAddress(&pfjets_trackCountingHighPurBJetTag_);
	}
	trkjets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_combinedSecondaryVertexBJetTag") != 0) {
		trkjets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_combinedSecondaryVertexBJetTag"));
		trkjets_combinedSecondaryVertexBJetTag_branch->SetAddress(&trkjets_combinedSecondaryVertexBJetTag_);
	}
	trkjets_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("trkjets_combinedSecondaryVertexMVABJetTag") != 0) {
		trkjets_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_combinedSecondaryVertexMVABJetTag"));
		trkjets_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&trkjets_combinedSecondaryVertexMVABJetTag_);
	}
	trkjets_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_jetBProbabilityBJetTag") != 0) {
		trkjets_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_jetBProbabilityBJetTag"));
		trkjets_jetBProbabilityBJetTag_branch->SetAddress(&trkjets_jetBProbabilityBJetTag_);
	}
	trkjets_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_jetProbabilityBJetTag") != 0) {
		trkjets_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_jetProbabilityBJetTag"));
		trkjets_jetProbabilityBJetTag_branch->SetAddress(&trkjets_jetProbabilityBJetTag_);
	}
	trkjets_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_simpleSecondaryVertexHighEffBJetTag") != 0) {
		trkjets_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_simpleSecondaryVertexHighEffBJetTag"));
		trkjets_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&trkjets_simpleSecondaryVertexHighEffBJetTag_);
	}
	trkjets_simpleSecondaryVertexHighPurBJetTags_branch = 0;
	if (tree->GetAlias("trkjets_simpleSecondaryVertexHighPurBJetTags") != 0) {
		trkjets_simpleSecondaryVertexHighPurBJetTags_branch = tree->GetBranch(tree->GetAlias("trkjets_simpleSecondaryVertexHighPurBJetTags"));
		trkjets_simpleSecondaryVertexHighPurBJetTags_branch->SetAddress(&trkjets_simpleSecondaryVertexHighPurBJetTags_);
	}
	trkjets_softElectronByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_softElectronByIP3dBJetTag") != 0) {
		trkjets_softElectronByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_softElectronByIP3dBJetTag"));
		trkjets_softElectronByIP3dBJetTag_branch->SetAddress(&trkjets_softElectronByIP3dBJetTag_);
	}
	trkjets_softElectronByPtBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_softElectronByPtBJetTag") != 0) {
		trkjets_softElectronByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_softElectronByPtBJetTag"));
		trkjets_softElectronByPtBJetTag_branch->SetAddress(&trkjets_softElectronByPtBJetTag_);
	}
	trkjets_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_softMuonBJetTag") != 0) {
		trkjets_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_softMuonBJetTag"));
		trkjets_softMuonBJetTag_branch->SetAddress(&trkjets_softMuonBJetTag_);
	}
	trkjets_softMuonByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_softMuonByIP3dBJetTag") != 0) {
		trkjets_softMuonByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_softMuonByIP3dBJetTag"));
		trkjets_softMuonByIP3dBJetTag_branch->SetAddress(&trkjets_softMuonByIP3dBJetTag_);
	}
	trkjets_softMuonByPtBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_softMuonByPtBJetTag") != 0) {
		trkjets_softMuonByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_softMuonByPtBJetTag"));
		trkjets_softMuonByPtBJetTag_branch->SetAddress(&trkjets_softMuonByPtBJetTag_);
	}
	trkjets_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_trackCountingHighEffBJetTag") != 0) {
		trkjets_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_trackCountingHighEffBJetTag"));
		trkjets_trackCountingHighEffBJetTag_branch->SetAddress(&trkjets_trackCountingHighEffBJetTag_);
	}
	trkjets_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("trkjets_trackCountingHighPurBJetTag") != 0) {
		trkjets_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("trkjets_trackCountingHighPurBJetTag"));
		trkjets_trackCountingHighPurBJetTag_branch->SetAddress(&trkjets_trackCountingHighPurBJetTag_);
	}
	evt_bs_covMatrix_branch = 0;
	if (tree->GetAlias("evt_bs_covMatrix") != 0) {
		evt_bs_covMatrix_branch = tree->GetBranch(tree->GetAlias("evt_bs_covMatrix"));
		evt_bs_covMatrix_branch->SetAddress(&evt_bs_covMatrix_);
	}
	els_mc3dr_branch = 0;
	if (tree->GetAlias("els_mc3dr") != 0) {
		els_mc3dr_branch = tree->GetBranch(tree->GetAlias("els_mc3dr"));
		els_mc3dr_branch->SetAddress(&els_mc3dr_);
	}
	els_mcdr_branch = 0;
	if (tree->GetAlias("els_mcdr") != 0) {
		els_mcdr_branch = tree->GetBranch(tree->GetAlias("els_mcdr"));
		els_mcdr_branch->SetAddress(&els_mcdr_);
	}
	jets_mc3dr_branch = 0;
	if (tree->GetAlias("jets_mc3dr") != 0) {
		jets_mc3dr_branch = tree->GetBranch(tree->GetAlias("jets_mc3dr"));
		jets_mc3dr_branch->SetAddress(&jets_mc3dr_);
	}
	jets_mcdr_branch = 0;
	if (tree->GetAlias("jets_mcdr") != 0) {
		jets_mcdr_branch = tree->GetBranch(tree->GetAlias("jets_mcdr"));
		jets_mcdr_branch->SetAddress(&jets_mcdr_);
	}
	jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_emEnergy") != 0) {
		jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_emEnergy"));
		jets_mc_emEnergy_branch->SetAddress(&jets_mc_emEnergy_);
	}
	jets_mc_gpdr_branch = 0;
	if (tree->GetAlias("jets_mc_gpdr") != 0) {
		jets_mc_gpdr_branch = tree->GetBranch(tree->GetAlias("jets_mc_gpdr"));
		jets_mc_gpdr_branch->SetAddress(&jets_mc_gpdr_);
	}
	jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_hadEnergy") != 0) {
		jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_hadEnergy"));
		jets_mc_hadEnergy_branch->SetAddress(&jets_mc_hadEnergy_);
	}
	jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_invEnergy") != 0) {
		jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_invEnergy"));
		jets_mc_invEnergy_branch->SetAddress(&jets_mc_invEnergy_);
	}
	jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("jets_mc_otherEnergy") != 0) {
		jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("jets_mc_otherEnergy"));
		jets_mc_otherEnergy_branch->SetAddress(&jets_mc_otherEnergy_);
	}
	mus_mc3dr_branch = 0;
	if (tree->GetAlias("mus_mc3dr") != 0) {
		mus_mc3dr_branch = tree->GetBranch(tree->GetAlias("mus_mc3dr"));
		mus_mc3dr_branch->SetAddress(&mus_mc3dr_);
	}
	mus_mcdr_branch = 0;
	if (tree->GetAlias("mus_mcdr") != 0) {
		mus_mcdr_branch = tree->GetBranch(tree->GetAlias("mus_mcdr"));
		mus_mcdr_branch->SetAddress(&mus_mcdr_);
	}
	pfjets_mc3dr_branch = 0;
	if (tree->GetAlias("pfjets_mc3dr") != 0) {
		pfjets_mc3dr_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3dr"));
		pfjets_mc3dr_branch->SetAddress(&pfjets_mc3dr_);
	}
	pfjets_mcdr_branch = 0;
	if (tree->GetAlias("pfjets_mcdr") != 0) {
		pfjets_mcdr_branch = tree->GetBranch(tree->GetAlias("pfjets_mcdr"));
		pfjets_mcdr_branch->SetAddress(&pfjets_mcdr_);
	}
	pfjets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_emEnergy") != 0) {
		pfjets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_emEnergy"));
		pfjets_mc_emEnergy_branch->SetAddress(&pfjets_mc_emEnergy_);
	}
	pfjets_mc_gpdr_branch = 0;
	if (tree->GetAlias("pfjets_mc_gpdr") != 0) {
		pfjets_mc_gpdr_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gpdr"));
		pfjets_mc_gpdr_branch->SetAddress(&pfjets_mc_gpdr_);
	}
	pfjets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_hadEnergy") != 0) {
		pfjets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_hadEnergy"));
		pfjets_mc_hadEnergy_branch->SetAddress(&pfjets_mc_hadEnergy_);
	}
	pfjets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_invEnergy") != 0) {
		pfjets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_invEnergy"));
		pfjets_mc_invEnergy_branch->SetAddress(&pfjets_mc_invEnergy_);
	}
	pfjets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("pfjets_mc_otherEnergy") != 0) {
		pfjets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_otherEnergy"));
		pfjets_mc_otherEnergy_branch->SetAddress(&pfjets_mc_otherEnergy_);
	}
	photons_mc3dr_branch = 0;
	if (tree->GetAlias("photons_mc3dr") != 0) {
		photons_mc3dr_branch = tree->GetBranch(tree->GetAlias("photons_mc3dr"));
		photons_mc3dr_branch->SetAddress(&photons_mc3dr_);
	}
	photons_mcdr_branch = 0;
	if (tree->GetAlias("photons_mcdr") != 0) {
		photons_mcdr_branch = tree->GetBranch(tree->GetAlias("photons_mcdr"));
		photons_mcdr_branch->SetAddress(&photons_mcdr_);
	}
	trk_mc3dr_branch = 0;
	if (tree->GetAlias("trk_mc3dr") != 0) {
		trk_mc3dr_branch = tree->GetBranch(tree->GetAlias("trk_mc3dr"));
		trk_mc3dr_branch->SetAddress(&trk_mc3dr_);
	}
	trk_mcdr_branch = 0;
	if (tree->GetAlias("trk_mcdr") != 0) {
		trk_mcdr_branch = tree->GetBranch(tree->GetAlias("trk_mcdr"));
		trk_mcdr_branch->SetAddress(&trk_mcdr_);
	}
	trks_conv_dcot_branch = 0;
	if (tree->GetAlias("trks_conv_dcot") != 0) {
		trks_conv_dcot_branch = tree->GetBranch(tree->GetAlias("trks_conv_dcot"));
		trks_conv_dcot_branch->SetAddress(&trks_conv_dcot_);
	}
	trks_conv_dist_branch = 0;
	if (tree->GetAlias("trks_conv_dist") != 0) {
		trks_conv_dist_branch = tree->GetBranch(tree->GetAlias("trks_conv_dist"));
		trks_conv_dist_branch->SetAddress(&trks_conv_dist_);
	}
	els_ecalJuraIso_branch = 0;
	if (tree->GetAlias("els_ecalJuraIso") != 0) {
		els_ecalJuraIso_branch = tree->GetBranch(tree->GetAlias("els_ecalJuraIso"));
		els_ecalJuraIso_branch->SetAddress(&els_ecalJuraIso_);
	}
	els_ecalJuraTowerIso_branch = 0;
	if (tree->GetAlias("els_ecalJuraTowerIso") != 0) {
		els_ecalJuraTowerIso_branch = tree->GetBranch(tree->GetAlias("els_ecalJuraTowerIso"));
		els_ecalJuraTowerIso_branch->SetAddress(&els_ecalJuraTowerIso_);
	}
	els_hcalConeIso_branch = 0;
	if (tree->GetAlias("els_hcalConeIso") != 0) {
		els_hcalConeIso_branch = tree->GetBranch(tree->GetAlias("els_hcalConeIso"));
		els_hcalConeIso_branch->SetAddress(&els_hcalConeIso_);
	}
	els_tkJuraIso_branch = 0;
	if (tree->GetAlias("els_tkJuraIso") != 0) {
		els_tkJuraIso_branch = tree->GetBranch(tree->GetAlias("els_tkJuraIso"));
		els_tkJuraIso_branch->SetAddress(&els_tkJuraIso_);
	}
	els_jetdr_branch = 0;
	if (tree->GetAlias("els_jetdr") != 0) {
		els_jetdr_branch = tree->GetBranch(tree->GetAlias("els_jetdr"));
		els_jetdr_branch->SetAddress(&els_jetdr_);
	}
	els_musdr_branch = 0;
	if (tree->GetAlias("els_musdr") != 0) {
		els_musdr_branch = tree->GetBranch(tree->GetAlias("els_musdr"));
		els_musdr_branch->SetAddress(&els_musdr_);
	}
	els_chi2_branch = 0;
	if (tree->GetAlias("els_chi2") != 0) {
		els_chi2_branch = tree->GetBranch(tree->GetAlias("els_chi2"));
		els_chi2_branch->SetAddress(&els_chi2_);
	}
	els_conv_dcot_branch = 0;
	if (tree->GetAlias("els_conv_dcot") != 0) {
		els_conv_dcot_branch = tree->GetBranch(tree->GetAlias("els_conv_dcot"));
		els_conv_dcot_branch->SetAddress(&els_conv_dcot_);
	}
	els_conv_dist_branch = 0;
	if (tree->GetAlias("els_conv_dist") != 0) {
		els_conv_dist_branch = tree->GetBranch(tree->GetAlias("els_conv_dist"));
		els_conv_dist_branch->SetAddress(&els_conv_dist_);
	}
	els_conv_radius_branch = 0;
	if (tree->GetAlias("els_conv_radius") != 0) {
		els_conv_radius_branch = tree->GetBranch(tree->GetAlias("els_conv_radius"));
		els_conv_radius_branch->SetAddress(&els_conv_radius_);
	}
	els_d0_branch = 0;
	if (tree->GetAlias("els_d0") != 0) {
		els_d0_branch = tree->GetBranch(tree->GetAlias("els_d0"));
		els_d0_branch->SetAddress(&els_d0_);
	}
	els_d0Err_branch = 0;
	if (tree->GetAlias("els_d0Err") != 0) {
		els_d0Err_branch = tree->GetBranch(tree->GetAlias("els_d0Err"));
		els_d0Err_branch->SetAddress(&els_d0Err_);
	}
	els_d0corr_branch = 0;
	if (tree->GetAlias("els_d0corr") != 0) {
		els_d0corr_branch = tree->GetBranch(tree->GetAlias("els_d0corr"));
		els_d0corr_branch->SetAddress(&els_d0corr_);
	}
	els_dEtaIn_branch = 0;
	if (tree->GetAlias("els_dEtaIn") != 0) {
		els_dEtaIn_branch = tree->GetBranch(tree->GetAlias("els_dEtaIn"));
		els_dEtaIn_branch->SetAddress(&els_dEtaIn_);
	}
	els_dEtaOut_branch = 0;
	if (tree->GetAlias("els_dEtaOut") != 0) {
		els_dEtaOut_branch = tree->GetBranch(tree->GetAlias("els_dEtaOut"));
		els_dEtaOut_branch->SetAddress(&els_dEtaOut_);
	}
	els_dPhiIn_branch = 0;
	if (tree->GetAlias("els_dPhiIn") != 0) {
		els_dPhiIn_branch = tree->GetBranch(tree->GetAlias("els_dPhiIn"));
		els_dPhiIn_branch->SetAddress(&els_dPhiIn_);
	}
	els_dPhiInPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiInPhiOut") != 0) {
		els_dPhiInPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiInPhiOut"));
		els_dPhiInPhiOut_branch->SetAddress(&els_dPhiInPhiOut_);
	}
	els_dPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiOut") != 0) {
		els_dPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiOut"));
		els_dPhiOut_branch->SetAddress(&els_dPhiOut_);
	}
	els_deltaEtaEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo") != 0) {
		els_deltaEtaEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaEtaEleClusterTrackAtCalo"));
		els_deltaEtaEleClusterTrackAtCalo_branch->SetAddress(&els_deltaEtaEleClusterTrackAtCalo_);
	}
	els_deltaPhiEleClusterTrackAtCalo_branch = 0;
	if (tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo") != 0) {
		els_deltaPhiEleClusterTrackAtCalo_branch = tree->GetBranch(tree->GetAlias("els_deltaPhiEleClusterTrackAtCalo"));
		els_deltaPhiEleClusterTrackAtCalo_branch->SetAddress(&els_deltaPhiEleClusterTrackAtCalo_);
	}
	els_e1x5_branch = 0;
	if (tree->GetAlias("els_e1x5") != 0) {
		els_e1x5_branch = tree->GetBranch(tree->GetAlias("els_e1x5"));
		els_e1x5_branch->SetAddress(&els_e1x5_);
	}
	els_e2x5Max_branch = 0;
	if (tree->GetAlias("els_e2x5Max") != 0) {
		els_e2x5Max_branch = tree->GetBranch(tree->GetAlias("els_e2x5Max"));
		els_e2x5Max_branch->SetAddress(&els_e2x5Max_);
	}
	els_e3x3_branch = 0;
	if (tree->GetAlias("els_e3x3") != 0) {
		els_e3x3_branch = tree->GetBranch(tree->GetAlias("els_e3x3"));
		els_e3x3_branch->SetAddress(&els_e3x3_);
	}
	els_e5x5_branch = 0;
	if (tree->GetAlias("els_e5x5") != 0) {
		els_e5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5"));
		els_e5x5_branch->SetAddress(&els_e5x5_);
	}
	els_eMax_branch = 0;
	if (tree->GetAlias("els_eMax") != 0) {
		els_eMax_branch = tree->GetBranch(tree->GetAlias("els_eMax"));
		els_eMax_branch->SetAddress(&els_eMax_);
	}
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		els_eOverPIn_branch->SetAddress(&els_eOverPIn_);
	}
	els_eOverPOut_branch = 0;
	if (tree->GetAlias("els_eOverPOut") != 0) {
		els_eOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eOverPOut"));
		els_eOverPOut_branch->SetAddress(&els_eOverPOut_);
	}
	els_eSC_branch = 0;
	if (tree->GetAlias("els_eSC") != 0) {
		els_eSC_branch = tree->GetBranch(tree->GetAlias("els_eSC"));
		els_eSC_branch->SetAddress(&els_eSC_);
	}
	els_eSCPresh_branch = 0;
	if (tree->GetAlias("els_eSCPresh") != 0) {
		els_eSCPresh_branch = tree->GetBranch(tree->GetAlias("els_eSCPresh"));
		els_eSCPresh_branch->SetAddress(&els_eSCPresh_);
	}
	els_eSCRaw_branch = 0;
	if (tree->GetAlias("els_eSCRaw") != 0) {
		els_eSCRaw_branch = tree->GetBranch(tree->GetAlias("els_eSCRaw"));
		els_eSCRaw_branch->SetAddress(&els_eSCRaw_);
	}
	els_eSeed_branch = 0;
	if (tree->GetAlias("els_eSeed") != 0) {
		els_eSeed_branch = tree->GetBranch(tree->GetAlias("els_eSeed"));
		els_eSeed_branch->SetAddress(&els_eSeed_);
	}
	els_eSeedOverPIn_branch = 0;
	if (tree->GetAlias("els_eSeedOverPIn") != 0) {
		els_eSeedOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPIn"));
		els_eSeedOverPIn_branch->SetAddress(&els_eSeedOverPIn_);
	}
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut_);
	}
	els_ecalEnergy_branch = 0;
	if (tree->GetAlias("els_ecalEnergy") != 0) {
		els_ecalEnergy_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergy"));
		els_ecalEnergy_branch->SetAddress(&els_ecalEnergy_);
	}
	els_ecalEnergyError_branch = 0;
	if (tree->GetAlias("els_ecalEnergyError") != 0) {
		els_ecalEnergyError_branch = tree->GetBranch(tree->GetAlias("els_ecalEnergyError"));
		els_ecalEnergyError_branch->SetAddress(&els_ecalEnergyError_);
	}
	els_ecalIso_branch = 0;
	if (tree->GetAlias("els_ecalIso") != 0) {
		els_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_ecalIso"));
		els_ecalIso_branch->SetAddress(&els_ecalIso_);
	}
	els_ecalIso04_branch = 0;
	if (tree->GetAlias("els_ecalIso04") != 0) {
		els_ecalIso04_branch = tree->GetBranch(tree->GetAlias("els_ecalIso04"));
		els_ecalIso04_branch->SetAddress(&els_ecalIso04_);
	}
	els_egamma_looseId_branch = 0;
	if (tree->GetAlias("els_egamma_looseId") != 0) {
		els_egamma_looseId_branch = tree->GetBranch(tree->GetAlias("els_egamma_looseId"));
		els_egamma_looseId_branch->SetAddress(&els_egamma_looseId_);
	}
	els_egamma_robustHighEnergy_branch = 0;
	if (tree->GetAlias("els_egamma_robustHighEnergy") != 0) {
		els_egamma_robustHighEnergy_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustHighEnergy"));
		els_egamma_robustHighEnergy_branch->SetAddress(&els_egamma_robustHighEnergy_);
	}
	els_egamma_robustLooseId_branch = 0;
	if (tree->GetAlias("els_egamma_robustLooseId") != 0) {
		els_egamma_robustLooseId_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustLooseId"));
		els_egamma_robustLooseId_branch->SetAddress(&els_egamma_robustLooseId_);
	}
	els_egamma_robustTightId_branch = 0;
	if (tree->GetAlias("els_egamma_robustTightId") != 0) {
		els_egamma_robustTightId_branch = tree->GetBranch(tree->GetAlias("els_egamma_robustTightId"));
		els_egamma_robustTightId_branch->SetAddress(&els_egamma_robustTightId_);
	}
	els_egamma_tightId_branch = 0;
	if (tree->GetAlias("els_egamma_tightId") != 0) {
		els_egamma_tightId_branch = tree->GetBranch(tree->GetAlias("els_egamma_tightId"));
		els_egamma_tightId_branch->SetAddress(&els_egamma_tightId_);
	}
	els_electronMomentumError_branch = 0;
	if (tree->GetAlias("els_electronMomentumError") != 0) {
		els_electronMomentumError_branch = tree->GetBranch(tree->GetAlias("els_electronMomentumError"));
		els_electronMomentumError_branch->SetAddress(&els_electronMomentumError_);
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		els_etaErr_branch->SetAddress(&els_etaErr_);
	}
	els_etaSC_branch = 0;
	if (tree->GetAlias("els_etaSC") != 0) {
		els_etaSC_branch = tree->GetBranch(tree->GetAlias("els_etaSC"));
		els_etaSC_branch->SetAddress(&els_etaSC_);
	}
	els_fbrem_branch = 0;
	if (tree->GetAlias("els_fbrem") != 0) {
		els_fbrem_branch = tree->GetBranch(tree->GetAlias("els_fbrem"));
		els_fbrem_branch->SetAddress(&els_fbrem_);
	}
	els_hOverE_branch = 0;
	if (tree->GetAlias("els_hOverE") != 0) {
		els_hOverE_branch = tree->GetBranch(tree->GetAlias("els_hOverE"));
		els_hOverE_branch->SetAddress(&els_hOverE_);
	}
	els_hcalDepth1OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth1OverEcal") != 0) {
		els_hcalDepth1OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1OverEcal"));
		els_hcalDepth1OverEcal_branch->SetAddress(&els_hcalDepth1OverEcal_);
	}
	els_hcalDepth1TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt") != 0) {
		els_hcalDepth1TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt"));
		els_hcalDepth1TowerSumEt_branch->SetAddress(&els_hcalDepth1TowerSumEt_);
	}
	els_hcalDepth1TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth1TowerSumEt04") != 0) {
		els_hcalDepth1TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth1TowerSumEt04"));
		els_hcalDepth1TowerSumEt04_branch->SetAddress(&els_hcalDepth1TowerSumEt04_);
	}
	els_hcalDepth2OverEcal_branch = 0;
	if (tree->GetAlias("els_hcalDepth2OverEcal") != 0) {
		els_hcalDepth2OverEcal_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2OverEcal"));
		els_hcalDepth2OverEcal_branch->SetAddress(&els_hcalDepth2OverEcal_);
	}
	els_hcalDepth2TowerSumEt_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt") != 0) {
		els_hcalDepth2TowerSumEt_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt"));
		els_hcalDepth2TowerSumEt_branch->SetAddress(&els_hcalDepth2TowerSumEt_);
	}
	els_hcalDepth2TowerSumEt04_branch = 0;
	if (tree->GetAlias("els_hcalDepth2TowerSumEt04") != 0) {
		els_hcalDepth2TowerSumEt04_branch = tree->GetBranch(tree->GetAlias("els_hcalDepth2TowerSumEt04"));
		els_hcalDepth2TowerSumEt04_branch->SetAddress(&els_hcalDepth2TowerSumEt04_);
	}
	els_hcalIso_branch = 0;
	if (tree->GetAlias("els_hcalIso") != 0) {
		els_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_hcalIso"));
		els_hcalIso_branch->SetAddress(&els_hcalIso_);
	}
	els_hcalIso04_branch = 0;
	if (tree->GetAlias("els_hcalIso04") != 0) {
		els_hcalIso04_branch = tree->GetBranch(tree->GetAlias("els_hcalIso04"));
		els_hcalIso04_branch->SetAddress(&els_hcalIso04_);
	}
	els_layer1_charge_branch = 0;
	if (tree->GetAlias("els_layer1_charge") != 0) {
		els_layer1_charge_branch = tree->GetBranch(tree->GetAlias("els_layer1_charge"));
		els_layer1_charge_branch->SetAddress(&els_layer1_charge_);
	}
	els_mva_branch = 0;
	if (tree->GetAlias("els_mva") != 0) {
		els_mva_branch = tree->GetBranch(tree->GetAlias("els_mva"));
		els_mva_branch->SetAddress(&els_mva_);
	}
	els_ndof_branch = 0;
	if (tree->GetAlias("els_ndof") != 0) {
		els_ndof_branch = tree->GetBranch(tree->GetAlias("els_ndof"));
		els_ndof_branch->SetAddress(&els_ndof_);
	}
	els_phiErr_branch = 0;
	if (tree->GetAlias("els_phiErr") != 0) {
		els_phiErr_branch = tree->GetBranch(tree->GetAlias("els_phiErr"));
		els_phiErr_branch->SetAddress(&els_phiErr_);
	}
	els_phiSC_branch = 0;
	if (tree->GetAlias("els_phiSC") != 0) {
		els_phiSC_branch = tree->GetBranch(tree->GetAlias("els_phiSC"));
		els_phiSC_branch->SetAddress(&els_phiSC_);
	}
	els_ptErr_branch = 0;
	if (tree->GetAlias("els_ptErr") != 0) {
		els_ptErr_branch = tree->GetBranch(tree->GetAlias("els_ptErr"));
		els_ptErr_branch->SetAddress(&els_ptErr_);
	}
	els_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta") != 0) {
		els_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta"));
		els_sigmaEtaEta_branch->SetAddress(&els_sigmaEtaEta_);
	}
	els_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEta") != 0) {
		els_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEta"));
		els_sigmaIEtaIEta_branch->SetAddress(&els_sigmaIEtaIEta_);
	}
	els_sigmaIEtaIEtaSC_branch = 0;
	if (tree->GetAlias("els_sigmaIEtaIEtaSC") != 0) {
		els_sigmaIEtaIEtaSC_branch = tree->GetBranch(tree->GetAlias("els_sigmaIEtaIEtaSC"));
		els_sigmaIEtaIEtaSC_branch->SetAddress(&els_sigmaIEtaIEtaSC_);
	}
	els_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhi") != 0) {
		els_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhi"));
		els_sigmaIPhiIPhi_branch->SetAddress(&els_sigmaIPhiIPhi_);
	}
	els_sigmaIPhiIPhiSC_branch = 0;
	if (tree->GetAlias("els_sigmaIPhiIPhiSC") != 0) {
		els_sigmaIPhiIPhiSC_branch = tree->GetBranch(tree->GetAlias("els_sigmaIPhiIPhiSC"));
		els_sigmaIPhiIPhiSC_branch->SetAddress(&els_sigmaIPhiIPhiSC_);
	}
	els_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("els_sigmaPhiPhi") != 0) {
		els_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaPhiPhi"));
		els_sigmaPhiPhi_branch->SetAddress(&els_sigmaPhiPhi_);
	}
	els_tkIso_branch = 0;
	if (tree->GetAlias("els_tkIso") != 0) {
		els_tkIso_branch = tree->GetBranch(tree->GetAlias("els_tkIso"));
		els_tkIso_branch->SetAddress(&els_tkIso_);
	}
	els_tkIso04_branch = 0;
	if (tree->GetAlias("els_tkIso04") != 0) {
		els_tkIso04_branch = tree->GetBranch(tree->GetAlias("els_tkIso04"));
		els_tkIso04_branch->SetAddress(&els_tkIso04_);
	}
	els_trackMomentumError_branch = 0;
	if (tree->GetAlias("els_trackMomentumError") != 0) {
		els_trackMomentumError_branch = tree->GetBranch(tree->GetAlias("els_trackMomentumError"));
		els_trackMomentumError_branch->SetAddress(&els_trackMomentumError_);
	}
	els_trkdr_branch = 0;
	if (tree->GetAlias("els_trkdr") != 0) {
		els_trkdr_branch = tree->GetBranch(tree->GetAlias("els_trkdr"));
		els_trkdr_branch->SetAddress(&els_trkdr_);
	}
	els_trkshFrac_branch = 0;
	if (tree->GetAlias("els_trkshFrac") != 0) {
		els_trkshFrac_branch = tree->GetBranch(tree->GetAlias("els_trkshFrac"));
		els_trkshFrac_branch->SetAddress(&els_trkshFrac_);
	}
	els_z0_branch = 0;
	if (tree->GetAlias("els_z0") != 0) {
		els_z0_branch = tree->GetBranch(tree->GetAlias("els_z0"));
		els_z0_branch->SetAddress(&els_z0_);
	}
	els_z0Err_branch = 0;
	if (tree->GetAlias("els_z0Err") != 0) {
		els_z0Err_branch = tree->GetBranch(tree->GetAlias("els_z0Err"));
		els_z0Err_branch->SetAddress(&els_z0Err_);
	}
	els_z0corr_branch = 0;
	if (tree->GetAlias("els_z0corr") != 0) {
		els_z0corr_branch = tree->GetBranch(tree->GetAlias("els_z0corr"));
		els_z0corr_branch->SetAddress(&els_z0corr_);
	}
	gsftrks_chi2_branch = 0;
	if (tree->GetAlias("gsftrks_chi2") != 0) {
		gsftrks_chi2_branch = tree->GetBranch(tree->GetAlias("gsftrks_chi2"));
		gsftrks_chi2_branch->SetAddress(&gsftrks_chi2_);
	}
	gsftrks_d0_branch = 0;
	if (tree->GetAlias("gsftrks_d0") != 0) {
		gsftrks_d0_branch = tree->GetBranch(tree->GetAlias("gsftrks_d0"));
		gsftrks_d0_branch->SetAddress(&gsftrks_d0_);
	}
	gsftrks_d0Err_branch = 0;
	if (tree->GetAlias("gsftrks_d0Err") != 0) {
		gsftrks_d0Err_branch = tree->GetBranch(tree->GetAlias("gsftrks_d0Err"));
		gsftrks_d0Err_branch->SetAddress(&gsftrks_d0Err_);
	}
	gsftrks_d0corr_branch = 0;
	if (tree->GetAlias("gsftrks_d0corr") != 0) {
		gsftrks_d0corr_branch = tree->GetBranch(tree->GetAlias("gsftrks_d0corr"));
		gsftrks_d0corr_branch->SetAddress(&gsftrks_d0corr_);
	}
	gsftrks_d0corrPhi_branch = 0;
	if (tree->GetAlias("gsftrks_d0corrPhi") != 0) {
		gsftrks_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("gsftrks_d0corrPhi"));
		gsftrks_d0corrPhi_branch->SetAddress(&gsftrks_d0corrPhi_);
	}
	gsftrks_d0phiCov_branch = 0;
	if (tree->GetAlias("gsftrks_d0phiCov") != 0) {
		gsftrks_d0phiCov_branch = tree->GetBranch(tree->GetAlias("gsftrks_d0phiCov"));
		gsftrks_d0phiCov_branch->SetAddress(&gsftrks_d0phiCov_);
	}
	gsftrks_etaErr_branch = 0;
	if (tree->GetAlias("gsftrks_etaErr") != 0) {
		gsftrks_etaErr_branch = tree->GetBranch(tree->GetAlias("gsftrks_etaErr"));
		gsftrks_etaErr_branch->SetAddress(&gsftrks_etaErr_);
	}
	gsftrks_layer1_charge_branch = 0;
	if (tree->GetAlias("gsftrks_layer1_charge") != 0) {
		gsftrks_layer1_charge_branch = tree->GetBranch(tree->GetAlias("gsftrks_layer1_charge"));
		gsftrks_layer1_charge_branch->SetAddress(&gsftrks_layer1_charge_);
	}
	gsftrks_ndof_branch = 0;
	if (tree->GetAlias("gsftrks_ndof") != 0) {
		gsftrks_ndof_branch = tree->GetBranch(tree->GetAlias("gsftrks_ndof"));
		gsftrks_ndof_branch->SetAddress(&gsftrks_ndof_);
	}
	gsftrks_phiErr_branch = 0;
	if (tree->GetAlias("gsftrks_phiErr") != 0) {
		gsftrks_phiErr_branch = tree->GetBranch(tree->GetAlias("gsftrks_phiErr"));
		gsftrks_phiErr_branch->SetAddress(&gsftrks_phiErr_);
	}
	gsftrks_ptErr_branch = 0;
	if (tree->GetAlias("gsftrks_ptErr") != 0) {
		gsftrks_ptErr_branch = tree->GetBranch(tree->GetAlias("gsftrks_ptErr"));
		gsftrks_ptErr_branch->SetAddress(&gsftrks_ptErr_);
	}
	gsftrks_z0_branch = 0;
	if (tree->GetAlias("gsftrks_z0") != 0) {
		gsftrks_z0_branch = tree->GetBranch(tree->GetAlias("gsftrks_z0"));
		gsftrks_z0_branch->SetAddress(&gsftrks_z0_);
	}
	gsftrks_z0Err_branch = 0;
	if (tree->GetAlias("gsftrks_z0Err") != 0) {
		gsftrks_z0Err_branch = tree->GetBranch(tree->GetAlias("gsftrks_z0Err"));
		gsftrks_z0Err_branch->SetAddress(&gsftrks_z0Err_);
	}
	gsftrks_z0corr_branch = 0;
	if (tree->GetAlias("gsftrks_z0corr") != 0) {
		gsftrks_z0corr_branch = tree->GetBranch(tree->GetAlias("gsftrks_z0corr"));
		gsftrks_z0corr_branch->SetAddress(&gsftrks_z0corr_);
	}
	hyp_Ht_branch = 0;
	if (tree->GetAlias("hyp_Ht") != 0) {
		hyp_Ht_branch = tree->GetBranch(tree->GetAlias("hyp_Ht"));
		hyp_Ht_branch->SetAddress(&hyp_Ht_);
	}
	hyp_dPhi_nJet_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_metMuonJESCorr") != 0) {
		hyp_dPhi_nJet_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_metMuonJESCorr"));
		hyp_dPhi_nJet_metMuonJESCorr_branch->SetAddress(&hyp_dPhi_nJet_metMuonJESCorr_);
	}
	hyp_dPhi_nJet_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_muCorrMet") != 0) {
		hyp_dPhi_nJet_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_muCorrMet"));
		hyp_dPhi_nJet_muCorrMet_branch->SetAddress(&hyp_dPhi_nJet_muCorrMet_);
	}
	hyp_dPhi_nJet_tcMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_tcMet") != 0) {
		hyp_dPhi_nJet_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_tcMet"));
		hyp_dPhi_nJet_tcMet_branch->SetAddress(&hyp_dPhi_nJet_tcMet_);
	}
	hyp_dPhi_nJet_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_unCorrMet") != 0) {
		hyp_dPhi_nJet_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_unCorrMet"));
		hyp_dPhi_nJet_unCorrMet_branch->SetAddress(&hyp_dPhi_nJet_unCorrMet_);
	}
	hyp_ll_chi2_branch = 0;
	if (tree->GetAlias("hyp_ll_chi2") != 0) {
		hyp_ll_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_ll_chi2"));
		hyp_ll_chi2_branch->SetAddress(&hyp_ll_chi2_);
	}
	hyp_ll_d0_branch = 0;
	if (tree->GetAlias("hyp_ll_d0") != 0) {
		hyp_ll_d0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0"));
		hyp_ll_d0_branch->SetAddress(&hyp_ll_d0_);
	}
	hyp_ll_d0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_d0Err") != 0) {
		hyp_ll_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0Err"));
		hyp_ll_d0Err_branch->SetAddress(&hyp_ll_d0Err_);
	}
	hyp_ll_d0corr_branch = 0;
	if (tree->GetAlias("hyp_ll_d0corr") != 0) {
		hyp_ll_d0corr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0corr"));
		hyp_ll_d0corr_branch->SetAddress(&hyp_ll_d0corr_);
	}
	hyp_ll_dPhi_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_metMuonJESCorr") != 0) {
		hyp_ll_dPhi_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_metMuonJESCorr"));
		hyp_ll_dPhi_metMuonJESCorr_branch->SetAddress(&hyp_ll_dPhi_metMuonJESCorr_);
	}
	hyp_ll_dPhi_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_muCorrMet") != 0) {
		hyp_ll_dPhi_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_muCorrMet"));
		hyp_ll_dPhi_muCorrMet_branch->SetAddress(&hyp_ll_dPhi_muCorrMet_);
	}
	hyp_ll_dPhi_tcMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_tcMet") != 0) {
		hyp_ll_dPhi_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_tcMet"));
		hyp_ll_dPhi_tcMet_branch->SetAddress(&hyp_ll_dPhi_tcMet_);
	}
	hyp_ll_dPhi_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_unCorrMet") != 0) {
		hyp_ll_dPhi_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_unCorrMet"));
		hyp_ll_dPhi_unCorrMet_branch->SetAddress(&hyp_ll_dPhi_unCorrMet_);
	}
	hyp_ll_etaErr_branch = 0;
	if (tree->GetAlias("hyp_ll_etaErr") != 0) {
		hyp_ll_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_etaErr"));
		hyp_ll_etaErr_branch->SetAddress(&hyp_ll_etaErr_);
	}
	hyp_ll_ndof_branch = 0;
	if (tree->GetAlias("hyp_ll_ndof") != 0) {
		hyp_ll_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ndof"));
		hyp_ll_ndof_branch->SetAddress(&hyp_ll_ndof_);
	}
	hyp_ll_phiErr_branch = 0;
	if (tree->GetAlias("hyp_ll_phiErr") != 0) {
		hyp_ll_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_phiErr"));
		hyp_ll_phiErr_branch->SetAddress(&hyp_ll_phiErr_);
	}
	hyp_ll_ptErr_branch = 0;
	if (tree->GetAlias("hyp_ll_ptErr") != 0) {
		hyp_ll_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ptErr"));
		hyp_ll_ptErr_branch->SetAddress(&hyp_ll_ptErr_);
	}
	hyp_ll_z0_branch = 0;
	if (tree->GetAlias("hyp_ll_z0") != 0) {
		hyp_ll_z0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0"));
		hyp_ll_z0_branch->SetAddress(&hyp_ll_z0_);
	}
	hyp_ll_z0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_z0Err") != 0) {
		hyp_ll_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0Err"));
		hyp_ll_z0Err_branch->SetAddress(&hyp_ll_z0Err_);
	}
	hyp_ll_z0corr_branch = 0;
	if (tree->GetAlias("hyp_ll_z0corr") != 0) {
		hyp_ll_z0corr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0corr"));
		hyp_ll_z0corr_branch->SetAddress(&hyp_ll_z0corr_);
	}
	hyp_lt_chi2_branch = 0;
	if (tree->GetAlias("hyp_lt_chi2") != 0) {
		hyp_lt_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_lt_chi2"));
		hyp_lt_chi2_branch->SetAddress(&hyp_lt_chi2_);
	}
	hyp_lt_d0_branch = 0;
	if (tree->GetAlias("hyp_lt_d0") != 0) {
		hyp_lt_d0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0"));
		hyp_lt_d0_branch->SetAddress(&hyp_lt_d0_);
	}
	hyp_lt_d0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_d0Err") != 0) {
		hyp_lt_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0Err"));
		hyp_lt_d0Err_branch->SetAddress(&hyp_lt_d0Err_);
	}
	hyp_lt_d0corr_branch = 0;
	if (tree->GetAlias("hyp_lt_d0corr") != 0) {
		hyp_lt_d0corr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0corr"));
		hyp_lt_d0corr_branch->SetAddress(&hyp_lt_d0corr_);
	}
	hyp_lt_dPhi_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_metMuonJESCorr") != 0) {
		hyp_lt_dPhi_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_metMuonJESCorr"));
		hyp_lt_dPhi_metMuonJESCorr_branch->SetAddress(&hyp_lt_dPhi_metMuonJESCorr_);
	}
	hyp_lt_dPhi_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_muCorrMet") != 0) {
		hyp_lt_dPhi_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_muCorrMet"));
		hyp_lt_dPhi_muCorrMet_branch->SetAddress(&hyp_lt_dPhi_muCorrMet_);
	}
	hyp_lt_dPhi_tcMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_tcMet") != 0) {
		hyp_lt_dPhi_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_tcMet"));
		hyp_lt_dPhi_tcMet_branch->SetAddress(&hyp_lt_dPhi_tcMet_);
	}
	hyp_lt_dPhi_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_unCorrMet") != 0) {
		hyp_lt_dPhi_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_unCorrMet"));
		hyp_lt_dPhi_unCorrMet_branch->SetAddress(&hyp_lt_dPhi_unCorrMet_);
	}
	hyp_lt_etaErr_branch = 0;
	if (tree->GetAlias("hyp_lt_etaErr") != 0) {
		hyp_lt_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_etaErr"));
		hyp_lt_etaErr_branch->SetAddress(&hyp_lt_etaErr_);
	}
	hyp_lt_ndof_branch = 0;
	if (tree->GetAlias("hyp_lt_ndof") != 0) {
		hyp_lt_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ndof"));
		hyp_lt_ndof_branch->SetAddress(&hyp_lt_ndof_);
	}
	hyp_lt_phiErr_branch = 0;
	if (tree->GetAlias("hyp_lt_phiErr") != 0) {
		hyp_lt_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_phiErr"));
		hyp_lt_phiErr_branch->SetAddress(&hyp_lt_phiErr_);
	}
	hyp_lt_ptErr_branch = 0;
	if (tree->GetAlias("hyp_lt_ptErr") != 0) {
		hyp_lt_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ptErr"));
		hyp_lt_ptErr_branch->SetAddress(&hyp_lt_ptErr_);
	}
	hyp_lt_z0_branch = 0;
	if (tree->GetAlias("hyp_lt_z0") != 0) {
		hyp_lt_z0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0"));
		hyp_lt_z0_branch->SetAddress(&hyp_lt_z0_);
	}
	hyp_lt_z0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_z0Err") != 0) {
		hyp_lt_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0Err"));
		hyp_lt_z0Err_branch->SetAddress(&hyp_lt_z0Err_);
	}
	hyp_lt_z0corr_branch = 0;
	if (tree->GetAlias("hyp_lt_z0corr") != 0) {
		hyp_lt_z0corr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0corr"));
		hyp_lt_z0corr_branch->SetAddress(&hyp_lt_z0corr_);
	}
	hyp_mt2_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("hyp_mt2_metMuonJESCorr") != 0) {
		hyp_mt2_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("hyp_mt2_metMuonJESCorr"));
		hyp_mt2_metMuonJESCorr_branch->SetAddress(&hyp_mt2_metMuonJESCorr_);
	}
	hyp_mt2_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_mt2_muCorrMet") != 0) {
		hyp_mt2_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_mt2_muCorrMet"));
		hyp_mt2_muCorrMet_branch->SetAddress(&hyp_mt2_muCorrMet_);
	}
	hyp_mt2_tcMet_branch = 0;
	if (tree->GetAlias("hyp_mt2_tcMet") != 0) {
		hyp_mt2_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_mt2_tcMet"));
		hyp_mt2_tcMet_branch->SetAddress(&hyp_mt2_tcMet_);
	}
	hyp_sumJetPt_branch = 0;
	if (tree->GetAlias("hyp_sumJetPt") != 0) {
		hyp_sumJetPt_branch = tree->GetBranch(tree->GetAlias("hyp_sumJetPt"));
		hyp_sumJetPt_branch->SetAddress(&hyp_sumJetPt_);
	}
	hyp_FVFit_chi2ndf_branch = 0;
	if (tree->GetAlias("hyp_FVFit_chi2ndf") != 0) {
		hyp_FVFit_chi2ndf_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_chi2ndf"));
		hyp_FVFit_chi2ndf_branch->SetAddress(&hyp_FVFit_chi2ndf_);
	}
	hyp_FVFit_prob_branch = 0;
	if (tree->GetAlias("hyp_FVFit_prob") != 0) {
		hyp_FVFit_prob_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_prob"));
		hyp_FVFit_prob_branch->SetAddress(&hyp_FVFit_prob_);
	}
	hyp_FVFit_v4cxx_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4cxx") != 0) {
		hyp_FVFit_v4cxx_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4cxx"));
		hyp_FVFit_v4cxx_branch->SetAddress(&hyp_FVFit_v4cxx_);
	}
	hyp_FVFit_v4cxy_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4cxy") != 0) {
		hyp_FVFit_v4cxy_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4cxy"));
		hyp_FVFit_v4cxy_branch->SetAddress(&hyp_FVFit_v4cxy_);
	}
	hyp_FVFit_v4cyy_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4cyy") != 0) {
		hyp_FVFit_v4cyy_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4cyy"));
		hyp_FVFit_v4cyy_branch->SetAddress(&hyp_FVFit_v4cyy_);
	}
	hyp_FVFit_v4czx_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4czx") != 0) {
		hyp_FVFit_v4czx_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4czx"));
		hyp_FVFit_v4czx_branch->SetAddress(&hyp_FVFit_v4czx_);
	}
	hyp_FVFit_v4czy_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4czy") != 0) {
		hyp_FVFit_v4czy_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4czy"));
		hyp_FVFit_v4czy_branch->SetAddress(&hyp_FVFit_v4czy_);
	}
	hyp_FVFit_v4czz_branch = 0;
	if (tree->GetAlias("hyp_FVFit_v4czz") != 0) {
		hyp_FVFit_v4czz_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_v4czz"));
		hyp_FVFit_v4czz_branch->SetAddress(&hyp_FVFit_v4czz_);
	}
	hyp_ll_ecaliso_branch = 0;
	if (tree->GetAlias("hyp_ll_ecaliso") != 0) {
		hyp_ll_ecaliso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ecaliso"));
		hyp_ll_ecaliso_branch->SetAddress(&hyp_ll_ecaliso_);
	}
	hyp_ll_trkiso_branch = 0;
	if (tree->GetAlias("hyp_ll_trkiso") != 0) {
		hyp_ll_trkiso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_trkiso"));
		hyp_ll_trkiso_branch->SetAddress(&hyp_ll_trkiso_);
	}
	hyp_lt_ecaliso_branch = 0;
	if (tree->GetAlias("hyp_lt_ecaliso") != 0) {
		hyp_lt_ecaliso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ecaliso"));
		hyp_lt_ecaliso_branch->SetAddress(&hyp_lt_ecaliso_);
	}
	hyp_lt_trkiso_branch = 0;
	if (tree->GetAlias("hyp_lt_trkiso") != 0) {
		hyp_lt_trkiso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_trkiso"));
		hyp_lt_trkiso_branch->SetAddress(&hyp_lt_trkiso_);
	}
	jets_approximatefHPD_branch = 0;
	if (tree->GetAlias("jets_approximatefHPD") != 0) {
		jets_approximatefHPD_branch = tree->GetBranch(tree->GetAlias("jets_approximatefHPD"));
		jets_approximatefHPD_branch->SetAddress(&jets_approximatefHPD_);
	}
	jets_approximatefRBX_branch = 0;
	if (tree->GetAlias("jets_approximatefRBX") != 0) {
		jets_approximatefRBX_branch = tree->GetBranch(tree->GetAlias("jets_approximatefRBX"));
		jets_approximatefRBX_branch->SetAddress(&jets_approximatefRBX_);
	}
	jets_cor_branch = 0;
	if (tree->GetAlias("jets_cor") != 0) {
		jets_cor_branch = tree->GetBranch(tree->GetAlias("jets_cor"));
		jets_cor_branch->SetAddress(&jets_cor_);
	}
	jets_emFrac_branch = 0;
	if (tree->GetAlias("jets_emFrac") != 0) {
		jets_emFrac_branch = tree->GetBranch(tree->GetAlias("jets_emFrac"));
		jets_emFrac_branch->SetAddress(&jets_emFrac_);
	}
	jets_fHPD_branch = 0;
	if (tree->GetAlias("jets_fHPD") != 0) {
		jets_fHPD_branch = tree->GetBranch(tree->GetAlias("jets_fHPD"));
		jets_fHPD_branch->SetAddress(&jets_fHPD_);
	}
	jets_fRBX_branch = 0;
	if (tree->GetAlias("jets_fRBX") != 0) {
		jets_fRBX_branch = tree->GetBranch(tree->GetAlias("jets_fRBX"));
		jets_fRBX_branch->SetAddress(&jets_fRBX_);
	}
	jets_fSubDetector1_branch = 0;
	if (tree->GetAlias("jets_fSubDetector1") != 0) {
		jets_fSubDetector1_branch = tree->GetBranch(tree->GetAlias("jets_fSubDetector1"));
		jets_fSubDetector1_branch->SetAddress(&jets_fSubDetector1_);
	}
	jets_fSubDetector2_branch = 0;
	if (tree->GetAlias("jets_fSubDetector2") != 0) {
		jets_fSubDetector2_branch = tree->GetBranch(tree->GetAlias("jets_fSubDetector2"));
		jets_fSubDetector2_branch->SetAddress(&jets_fSubDetector2_);
	}
	jets_fSubDetector3_branch = 0;
	if (tree->GetAlias("jets_fSubDetector3") != 0) {
		jets_fSubDetector3_branch = tree->GetBranch(tree->GetAlias("jets_fSubDetector3"));
		jets_fSubDetector3_branch->SetAddress(&jets_fSubDetector3_);
	}
	jets_fSubDetector4_branch = 0;
	if (tree->GetAlias("jets_fSubDetector4") != 0) {
		jets_fSubDetector4_branch = tree->GetBranch(tree->GetAlias("jets_fSubDetector4"));
		jets_fSubDetector4_branch->SetAddress(&jets_fSubDetector4_);
	}
	jets_hitsInN90_branch = 0;
	if (tree->GetAlias("jets_hitsInN90") != 0) {
		jets_hitsInN90_branch = tree->GetBranch(tree->GetAlias("jets_hitsInN90"));
		jets_hitsInN90_branch->SetAddress(&jets_hitsInN90_);
	}
	jets_n90Hits_branch = 0;
	if (tree->GetAlias("jets_n90Hits") != 0) {
		jets_n90Hits_branch = tree->GetBranch(tree->GetAlias("jets_n90Hits"));
		jets_n90Hits_branch->SetAddress(&jets_n90Hits_);
	}
	jets_nECALTowers_branch = 0;
	if (tree->GetAlias("jets_nECALTowers") != 0) {
		jets_nECALTowers_branch = tree->GetBranch(tree->GetAlias("jets_nECALTowers"));
		jets_nECALTowers_branch->SetAddress(&jets_nECALTowers_);
	}
	jets_nHCALTowers_branch = 0;
	if (tree->GetAlias("jets_nHCALTowers") != 0) {
		jets_nHCALTowers_branch = tree->GetBranch(tree->GetAlias("jets_nHCALTowers"));
		jets_nHCALTowers_branch->SetAddress(&jets_nHCALTowers_);
	}
	jets_restrictedEMF_branch = 0;
	if (tree->GetAlias("jets_restrictedEMF") != 0) {
		jets_restrictedEMF_branch = tree->GetBranch(tree->GetAlias("jets_restrictedEMF"));
		jets_restrictedEMF_branch->SetAddress(&jets_restrictedEMF_);
	}
	jpts_cor_branch = 0;
	if (tree->GetAlias("jpts_cor") != 0) {
		jpts_cor_branch = tree->GetBranch(tree->GetAlias("jpts_cor"));
		jpts_cor_branch->SetAddress(&jpts_cor_);
	}
	jpts_emFrac_branch = 0;
	if (tree->GetAlias("jpts_emFrac") != 0) {
		jpts_emFrac_branch = tree->GetBranch(tree->GetAlias("jpts_emFrac"));
		jpts_emFrac_branch->SetAddress(&jpts_emFrac_);
	}
	mus_met_deltax_branch = 0;
	if (tree->GetAlias("mus_met_deltax") != 0) {
		mus_met_deltax_branch = tree->GetBranch(tree->GetAlias("mus_met_deltax"));
		mus_met_deltax_branch->SetAddress(&mus_met_deltax_);
	}
	mus_met_deltay_branch = 0;
	if (tree->GetAlias("mus_met_deltay") != 0) {
		mus_met_deltay_branch = tree->GetBranch(tree->GetAlias("mus_met_deltay"));
		mus_met_deltay_branch->SetAddress(&mus_met_deltay_);
	}
	mus_eledr_branch = 0;
	if (tree->GetAlias("mus_eledr") != 0) {
		mus_eledr_branch = tree->GetBranch(tree->GetAlias("mus_eledr"));
		mus_eledr_branch->SetAddress(&mus_eledr_);
	}
	mus_jetdr_branch = 0;
	if (tree->GetAlias("mus_jetdr") != 0) {
		mus_jetdr_branch = tree->GetBranch(tree->GetAlias("mus_jetdr"));
		mus_jetdr_branch->SetAddress(&mus_jetdr_);
	}
	mus_caloCompatibility_branch = 0;
	if (tree->GetAlias("mus_caloCompatibility") != 0) {
		mus_caloCompatibility_branch = tree->GetBranch(tree->GetAlias("mus_caloCompatibility"));
		mus_caloCompatibility_branch->SetAddress(&mus_caloCompatibility_);
	}
	mus_chi2_branch = 0;
	if (tree->GetAlias("mus_chi2") != 0) {
		mus_chi2_branch = tree->GetBranch(tree->GetAlias("mus_chi2"));
		mus_chi2_branch->SetAddress(&mus_chi2_);
	}
	mus_d0_branch = 0;
	if (tree->GetAlias("mus_d0") != 0) {
		mus_d0_branch = tree->GetBranch(tree->GetAlias("mus_d0"));
		mus_d0_branch->SetAddress(&mus_d0_);
	}
	mus_d0Err_branch = 0;
	if (tree->GetAlias("mus_d0Err") != 0) {
		mus_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_d0Err"));
		mus_d0Err_branch->SetAddress(&mus_d0Err_);
	}
	mus_d0corr_branch = 0;
	if (tree->GetAlias("mus_d0corr") != 0) {
		mus_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_d0corr"));
		mus_d0corr_branch->SetAddress(&mus_d0corr_);
	}
	mus_e_em_branch = 0;
	if (tree->GetAlias("mus_e_em") != 0) {
		mus_e_em_branch = tree->GetBranch(tree->GetAlias("mus_e_em"));
		mus_e_em_branch->SetAddress(&mus_e_em_);
	}
	mus_e_emS9_branch = 0;
	if (tree->GetAlias("mus_e_emS9") != 0) {
		mus_e_emS9_branch = tree->GetBranch(tree->GetAlias("mus_e_emS9"));
		mus_e_emS9_branch->SetAddress(&mus_e_emS9_);
	}
	mus_e_had_branch = 0;
	if (tree->GetAlias("mus_e_had") != 0) {
		mus_e_had_branch = tree->GetBranch(tree->GetAlias("mus_e_had"));
		mus_e_had_branch->SetAddress(&mus_e_had_);
	}
	mus_e_hadS9_branch = 0;
	if (tree->GetAlias("mus_e_hadS9") != 0) {
		mus_e_hadS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hadS9"));
		mus_e_hadS9_branch->SetAddress(&mus_e_hadS9_);
	}
	mus_e_ho_branch = 0;
	if (tree->GetAlias("mus_e_ho") != 0) {
		mus_e_ho_branch = tree->GetBranch(tree->GetAlias("mus_e_ho"));
		mus_e_ho_branch->SetAddress(&mus_e_ho_);
	}
	mus_e_hoS9_branch = 0;
	if (tree->GetAlias("mus_e_hoS9") != 0) {
		mus_e_hoS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hoS9"));
		mus_e_hoS9_branch->SetAddress(&mus_e_hoS9_);
	}
	mus_etaErr_branch = 0;
	if (tree->GetAlias("mus_etaErr") != 0) {
		mus_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_etaErr"));
		mus_etaErr_branch->SetAddress(&mus_etaErr_);
	}
	mus_gfit_chi2_branch = 0;
	if (tree->GetAlias("mus_gfit_chi2") != 0) {
		mus_gfit_chi2_branch = tree->GetBranch(tree->GetAlias("mus_gfit_chi2"));
		mus_gfit_chi2_branch->SetAddress(&mus_gfit_chi2_);
	}
	mus_gfit_d0_branch = 0;
	if (tree->GetAlias("mus_gfit_d0") != 0) {
		mus_gfit_d0_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0"));
		mus_gfit_d0_branch->SetAddress(&mus_gfit_d0_);
	}
	mus_gfit_d0Err_branch = 0;
	if (tree->GetAlias("mus_gfit_d0Err") != 0) {
		mus_gfit_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0Err"));
		mus_gfit_d0Err_branch->SetAddress(&mus_gfit_d0Err_);
	}
	mus_gfit_d0corr_branch = 0;
	if (tree->GetAlias("mus_gfit_d0corr") != 0) {
		mus_gfit_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_d0corr"));
		mus_gfit_d0corr_branch->SetAddress(&mus_gfit_d0corr_);
	}
	mus_gfit_ndof_branch = 0;
	if (tree->GetAlias("mus_gfit_ndof") != 0) {
		mus_gfit_ndof_branch = tree->GetBranch(tree->GetAlias("mus_gfit_ndof"));
		mus_gfit_ndof_branch->SetAddress(&mus_gfit_ndof_);
	}
	mus_gfit_qoverp_branch = 0;
	if (tree->GetAlias("mus_gfit_qoverp") != 0) {
		mus_gfit_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_gfit_qoverp"));
		mus_gfit_qoverp_branch->SetAddress(&mus_gfit_qoverp_);
	}
	mus_gfit_qoverpError_branch = 0;
	if (tree->GetAlias("mus_gfit_qoverpError") != 0) {
		mus_gfit_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_gfit_qoverpError"));
		mus_gfit_qoverpError_branch->SetAddress(&mus_gfit_qoverpError_);
	}
	mus_gfit_z0_branch = 0;
	if (tree->GetAlias("mus_gfit_z0") != 0) {
		mus_gfit_z0_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0"));
		mus_gfit_z0_branch->SetAddress(&mus_gfit_z0_);
	}
	mus_gfit_z0Err_branch = 0;
	if (tree->GetAlias("mus_gfit_z0Err") != 0) {
		mus_gfit_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0Err"));
		mus_gfit_z0Err_branch->SetAddress(&mus_gfit_z0Err_);
	}
	mus_gfit_z0corr_branch = 0;
	if (tree->GetAlias("mus_gfit_z0corr") != 0) {
		mus_gfit_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_gfit_z0corr"));
		mus_gfit_z0corr_branch->SetAddress(&mus_gfit_z0corr_);
	}
	mus_iso03_emEt_branch = 0;
	if (tree->GetAlias("mus_iso03_emEt") != 0) {
		mus_iso03_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_emEt"));
		mus_iso03_emEt_branch->SetAddress(&mus_iso03_emEt_);
	}
	mus_iso03_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hadEt") != 0) {
		mus_iso03_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hadEt"));
		mus_iso03_hadEt_branch->SetAddress(&mus_iso03_hadEt_);
	}
	mus_iso03_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hoEt") != 0) {
		mus_iso03_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hoEt"));
		mus_iso03_hoEt_branch->SetAddress(&mus_iso03_hoEt_);
	}
	mus_iso03_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso03_sumPt") != 0) {
		mus_iso03_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_sumPt"));
		mus_iso03_sumPt_branch->SetAddress(&mus_iso03_sumPt_);
	}
	mus_iso05_emEt_branch = 0;
	if (tree->GetAlias("mus_iso05_emEt") != 0) {
		mus_iso05_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_emEt"));
		mus_iso05_emEt_branch->SetAddress(&mus_iso05_emEt_);
	}
	mus_iso05_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hadEt") != 0) {
		mus_iso05_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hadEt"));
		mus_iso05_hadEt_branch->SetAddress(&mus_iso05_hadEt_);
	}
	mus_iso05_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hoEt") != 0) {
		mus_iso05_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hoEt"));
		mus_iso05_hoEt_branch->SetAddress(&mus_iso05_hoEt_);
	}
	mus_iso05_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso05_sumPt") != 0) {
		mus_iso05_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_sumPt"));
		mus_iso05_sumPt_branch->SetAddress(&mus_iso05_sumPt_);
	}
	mus_iso_ecalvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_ecalvetoDep") != 0) {
		mus_iso_ecalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_ecalvetoDep"));
		mus_iso_ecalvetoDep_branch->SetAddress(&mus_iso_ecalvetoDep_);
	}
	mus_iso_hcalvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_hcalvetoDep") != 0) {
		mus_iso_hcalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_hcalvetoDep"));
		mus_iso_hcalvetoDep_branch->SetAddress(&mus_iso_hcalvetoDep_);
	}
	mus_iso_hovetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_hovetoDep") != 0) {
		mus_iso_hovetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_hovetoDep"));
		mus_iso_hovetoDep_branch->SetAddress(&mus_iso_hovetoDep_);
	}
	mus_iso_trckvetoDep_branch = 0;
	if (tree->GetAlias("mus_iso_trckvetoDep") != 0) {
		mus_iso_trckvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_iso_trckvetoDep"));
		mus_iso_trckvetoDep_branch->SetAddress(&mus_iso_trckvetoDep_);
	}
	mus_ndof_branch = 0;
	if (tree->GetAlias("mus_ndof") != 0) {
		mus_ndof_branch = tree->GetBranch(tree->GetAlias("mus_ndof"));
		mus_ndof_branch->SetAddress(&mus_ndof_);
	}
	mus_phiErr_branch = 0;
	if (tree->GetAlias("mus_phiErr") != 0) {
		mus_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_phiErr"));
		mus_phiErr_branch->SetAddress(&mus_phiErr_);
	}
	mus_ptErr_branch = 0;
	if (tree->GetAlias("mus_ptErr") != 0) {
		mus_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_ptErr"));
		mus_ptErr_branch->SetAddress(&mus_ptErr_);
	}
	mus_qoverp_branch = 0;
	if (tree->GetAlias("mus_qoverp") != 0) {
		mus_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_qoverp"));
		mus_qoverp_branch->SetAddress(&mus_qoverp_);
	}
	mus_qoverpError_branch = 0;
	if (tree->GetAlias("mus_qoverpError") != 0) {
		mus_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_qoverpError"));
		mus_qoverpError_branch->SetAddress(&mus_qoverpError_);
	}
	mus_sta_chi2_branch = 0;
	if (tree->GetAlias("mus_sta_chi2") != 0) {
		mus_sta_chi2_branch = tree->GetBranch(tree->GetAlias("mus_sta_chi2"));
		mus_sta_chi2_branch->SetAddress(&mus_sta_chi2_);
	}
	mus_sta_d0_branch = 0;
	if (tree->GetAlias("mus_sta_d0") != 0) {
		mus_sta_d0_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0"));
		mus_sta_d0_branch->SetAddress(&mus_sta_d0_);
	}
	mus_sta_d0Err_branch = 0;
	if (tree->GetAlias("mus_sta_d0Err") != 0) {
		mus_sta_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0Err"));
		mus_sta_d0Err_branch->SetAddress(&mus_sta_d0Err_);
	}
	mus_sta_d0corr_branch = 0;
	if (tree->GetAlias("mus_sta_d0corr") != 0) {
		mus_sta_d0corr_branch = tree->GetBranch(tree->GetAlias("mus_sta_d0corr"));
		mus_sta_d0corr_branch->SetAddress(&mus_sta_d0corr_);
	}
	mus_sta_ndof_branch = 0;
	if (tree->GetAlias("mus_sta_ndof") != 0) {
		mus_sta_ndof_branch = tree->GetBranch(tree->GetAlias("mus_sta_ndof"));
		mus_sta_ndof_branch->SetAddress(&mus_sta_ndof_);
	}
	mus_sta_qoverp_branch = 0;
	if (tree->GetAlias("mus_sta_qoverp") != 0) {
		mus_sta_qoverp_branch = tree->GetBranch(tree->GetAlias("mus_sta_qoverp"));
		mus_sta_qoverp_branch->SetAddress(&mus_sta_qoverp_);
	}
	mus_sta_qoverpError_branch = 0;
	if (tree->GetAlias("mus_sta_qoverpError") != 0) {
		mus_sta_qoverpError_branch = tree->GetBranch(tree->GetAlias("mus_sta_qoverpError"));
		mus_sta_qoverpError_branch->SetAddress(&mus_sta_qoverpError_);
	}
	mus_sta_z0_branch = 0;
	if (tree->GetAlias("mus_sta_z0") != 0) {
		mus_sta_z0_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0"));
		mus_sta_z0_branch->SetAddress(&mus_sta_z0_);
	}
	mus_sta_z0Err_branch = 0;
	if (tree->GetAlias("mus_sta_z0Err") != 0) {
		mus_sta_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0Err"));
		mus_sta_z0Err_branch->SetAddress(&mus_sta_z0Err_);
	}
	mus_sta_z0corr_branch = 0;
	if (tree->GetAlias("mus_sta_z0corr") != 0) {
		mus_sta_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_sta_z0corr"));
		mus_sta_z0corr_branch->SetAddress(&mus_sta_z0corr_);
	}
	mus_timeAtIpInOut_branch = 0;
	if (tree->GetAlias("mus_timeAtIpInOut") != 0) {
		mus_timeAtIpInOut_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpInOut"));
		mus_timeAtIpInOut_branch->SetAddress(&mus_timeAtIpInOut_);
	}
	mus_timeAtIpInOutErr_branch = 0;
	if (tree->GetAlias("mus_timeAtIpInOutErr") != 0) {
		mus_timeAtIpInOutErr_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpInOutErr"));
		mus_timeAtIpInOutErr_branch->SetAddress(&mus_timeAtIpInOutErr_);
	}
	mus_timeAtIpOutIn_branch = 0;
	if (tree->GetAlias("mus_timeAtIpOutIn") != 0) {
		mus_timeAtIpOutIn_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpOutIn"));
		mus_timeAtIpOutIn_branch->SetAddress(&mus_timeAtIpOutIn_);
	}
	mus_timeAtIpOutInErr_branch = 0;
	if (tree->GetAlias("mus_timeAtIpOutInErr") != 0) {
		mus_timeAtIpOutInErr_branch = tree->GetBranch(tree->GetAlias("mus_timeAtIpOutInErr"));
		mus_timeAtIpOutInErr_branch->SetAddress(&mus_timeAtIpOutInErr_);
	}
	mus_vertexphi_branch = 0;
	if (tree->GetAlias("mus_vertexphi") != 0) {
		mus_vertexphi_branch = tree->GetBranch(tree->GetAlias("mus_vertexphi"));
		mus_vertexphi_branch->SetAddress(&mus_vertexphi_);
	}
	mus_z0_branch = 0;
	if (tree->GetAlias("mus_z0") != 0) {
		mus_z0_branch = tree->GetBranch(tree->GetAlias("mus_z0"));
		mus_z0_branch->SetAddress(&mus_z0_);
	}
	mus_z0Err_branch = 0;
	if (tree->GetAlias("mus_z0Err") != 0) {
		mus_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_z0Err"));
		mus_z0Err_branch->SetAddress(&mus_z0Err_);
	}
	mus_z0corr_branch = 0;
	if (tree->GetAlias("mus_z0corr") != 0) {
		mus_z0corr_branch = tree->GetBranch(tree->GetAlias("mus_z0corr"));
		mus_z0corr_branch->SetAddress(&mus_z0corr_);
	}
	els_pat_caloIso_branch = 0;
	if (tree->GetAlias("els_pat_caloIso") != 0) {
		els_pat_caloIso_branch = tree->GetBranch(tree->GetAlias("els_pat_caloIso"));
		els_pat_caloIso_branch->SetAddress(&els_pat_caloIso_);
	}
	els_pat_ecalIso_branch = 0;
	if (tree->GetAlias("els_pat_ecalIso") != 0) {
		els_pat_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_pat_ecalIso"));
		els_pat_ecalIso_branch->SetAddress(&els_pat_ecalIso_);
	}
	els_pat_hcalIso_branch = 0;
	if (tree->GetAlias("els_pat_hcalIso") != 0) {
		els_pat_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_pat_hcalIso"));
		els_pat_hcalIso_branch->SetAddress(&els_pat_hcalIso_);
	}
	els_pat_looseId_branch = 0;
	if (tree->GetAlias("els_pat_looseId") != 0) {
		els_pat_looseId_branch = tree->GetBranch(tree->GetAlias("els_pat_looseId"));
		els_pat_looseId_branch->SetAddress(&els_pat_looseId_);
	}
	els_pat_robustHighEnergy_branch = 0;
	if (tree->GetAlias("els_pat_robustHighEnergy") != 0) {
		els_pat_robustHighEnergy_branch = tree->GetBranch(tree->GetAlias("els_pat_robustHighEnergy"));
		els_pat_robustHighEnergy_branch->SetAddress(&els_pat_robustHighEnergy_);
	}
	els_pat_robustLooseId_branch = 0;
	if (tree->GetAlias("els_pat_robustLooseId") != 0) {
		els_pat_robustLooseId_branch = tree->GetBranch(tree->GetAlias("els_pat_robustLooseId"));
		els_pat_robustLooseId_branch->SetAddress(&els_pat_robustLooseId_);
	}
	els_pat_robustTightId_branch = 0;
	if (tree->GetAlias("els_pat_robustTightId") != 0) {
		els_pat_robustTightId_branch = tree->GetBranch(tree->GetAlias("els_pat_robustTightId"));
		els_pat_robustTightId_branch->SetAddress(&els_pat_robustTightId_);
	}
	els_pat_scE1x5_branch = 0;
	if (tree->GetAlias("els_pat_scE1x5") != 0) {
		els_pat_scE1x5_branch = tree->GetBranch(tree->GetAlias("els_pat_scE1x5"));
		els_pat_scE1x5_branch->SetAddress(&els_pat_scE1x5_);
	}
	els_pat_scE2x5Max_branch = 0;
	if (tree->GetAlias("els_pat_scE2x5Max") != 0) {
		els_pat_scE2x5Max_branch = tree->GetBranch(tree->GetAlias("els_pat_scE2x5Max"));
		els_pat_scE2x5Max_branch->SetAddress(&els_pat_scE2x5Max_);
	}
	els_pat_scE5x5_branch = 0;
	if (tree->GetAlias("els_pat_scE5x5") != 0) {
		els_pat_scE5x5_branch = tree->GetBranch(tree->GetAlias("els_pat_scE5x5"));
		els_pat_scE5x5_branch->SetAddress(&els_pat_scE5x5_);
	}
	els_pat_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_pat_sigmaEtaEta") != 0) {
		els_pat_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_pat_sigmaEtaEta"));
		els_pat_sigmaEtaEta_branch->SetAddress(&els_pat_sigmaEtaEta_);
	}
	els_pat_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_pat_sigmaIEtaIEta") != 0) {
		els_pat_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_pat_sigmaIEtaIEta"));
		els_pat_sigmaIEtaIEta_branch->SetAddress(&els_pat_sigmaIEtaIEta_);
	}
	els_pat_tightId_branch = 0;
	if (tree->GetAlias("els_pat_tightId") != 0) {
		els_pat_tightId_branch = tree->GetBranch(tree->GetAlias("els_pat_tightId"));
		els_pat_tightId_branch->SetAddress(&els_pat_tightId_);
	}
	els_pat_trackIso_branch = 0;
	if (tree->GetAlias("els_pat_trackIso") != 0) {
		els_pat_trackIso_branch = tree->GetBranch(tree->GetAlias("els_pat_trackIso"));
		els_pat_trackIso_branch->SetAddress(&els_pat_trackIso_);
	}
	jets_pat_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_combinedSecondaryVertexBJetTag") != 0) {
		jets_pat_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_combinedSecondaryVertexBJetTag"));
		jets_pat_combinedSecondaryVertexBJetTag_branch->SetAddress(&jets_pat_combinedSecondaryVertexBJetTag_);
	}
	jets_pat_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_combinedSecondaryVertexMVABJetTag") != 0) {
		jets_pat_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_combinedSecondaryVertexMVABJetTag"));
		jets_pat_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&jets_pat_combinedSecondaryVertexMVABJetTag_);
	}
	jets_pat_coneIsolationTauJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_coneIsolationTauJetTag") != 0) {
		jets_pat_coneIsolationTauJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_coneIsolationTauJetTag"));
		jets_pat_coneIsolationTauJetTag_branch->SetAddress(&jets_pat_coneIsolationTauJetTag_);
	}
	jets_pat_impactParameterMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_impactParameterMVABJetTag") != 0) {
		jets_pat_impactParameterMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_impactParameterMVABJetTag"));
		jets_pat_impactParameterMVABJetTag_branch->SetAddress(&jets_pat_impactParameterMVABJetTag_);
	}
	jets_pat_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_jetBProbabilityBJetTag") != 0) {
		jets_pat_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetBProbabilityBJetTag"));
		jets_pat_jetBProbabilityBJetTag_branch->SetAddress(&jets_pat_jetBProbabilityBJetTag_);
	}
	jets_pat_jetCharge_branch = 0;
	if (tree->GetAlias("jets_pat_jetCharge") != 0) {
		jets_pat_jetCharge_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetCharge"));
		jets_pat_jetCharge_branch->SetAddress(&jets_pat_jetCharge_);
	}
	jets_pat_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_jetProbabilityBJetTag") != 0) {
		jets_pat_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetProbabilityBJetTag"));
		jets_pat_jetProbabilityBJetTag_branch->SetAddress(&jets_pat_jetProbabilityBJetTag_);
	}
	jets_pat_noCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_noCorrF") != 0) {
		jets_pat_noCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_noCorrF"));
		jets_pat_noCorrF_branch->SetAddress(&jets_pat_noCorrF_);
	}
	jets_pat_simpleSecondaryVertexHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_simpleSecondaryVertexHighEffBJetTag") != 0) {
		jets_pat_simpleSecondaryVertexHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_simpleSecondaryVertexHighEffBJetTag"));
		jets_pat_simpleSecondaryVertexHighEffBJetTag_branch->SetAddress(&jets_pat_simpleSecondaryVertexHighEffBJetTag_);
	}
	jets_pat_simpleSecondaryVertexHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_simpleSecondaryVertexHighPurBJetTag") != 0) {
		jets_pat_simpleSecondaryVertexHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_simpleSecondaryVertexHighPurBJetTag"));
		jets_pat_simpleSecondaryVertexHighPurBJetTag_branch->SetAddress(&jets_pat_simpleSecondaryVertexHighPurBJetTag_);
	}
	jets_pat_softElectronByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softElectronByIP3dBJetTag") != 0) {
		jets_pat_softElectronByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softElectronByIP3dBJetTag"));
		jets_pat_softElectronByIP3dBJetTag_branch->SetAddress(&jets_pat_softElectronByIP3dBJetTag_);
	}
	jets_pat_softElectronByPtBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softElectronByPtBJetTag") != 0) {
		jets_pat_softElectronByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softElectronByPtBJetTag"));
		jets_pat_softElectronByPtBJetTag_branch->SetAddress(&jets_pat_softElectronByPtBJetTag_);
	}
	jets_pat_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softMuonBJetTag") != 0) {
		jets_pat_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softMuonBJetTag"));
		jets_pat_softMuonBJetTag_branch->SetAddress(&jets_pat_softMuonBJetTag_);
	}
	jets_pat_softMuonByIP3dBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softMuonByIP3dBJetTag") != 0) {
		jets_pat_softMuonByIP3dBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softMuonByIP3dBJetTag"));
		jets_pat_softMuonByIP3dBJetTag_branch->SetAddress(&jets_pat_softMuonByIP3dBJetTag_);
	}
	jets_pat_softMuonByPtBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softMuonByPtBJetTag") != 0) {
		jets_pat_softMuonByPtBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softMuonByPtBJetTag"));
		jets_pat_softMuonByPtBJetTag_branch->SetAddress(&jets_pat_softMuonByPtBJetTag_);
	}
	jets_pat_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_trackCountingHighEffBJetTag") != 0) {
		jets_pat_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_trackCountingHighEffBJetTag"));
		jets_pat_trackCountingHighEffBJetTag_branch->SetAddress(&jets_pat_trackCountingHighEffBJetTag_);
	}
	jets_pat_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_trackCountingHighPurBJetTag") != 0) {
		jets_pat_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_trackCountingHighPurBJetTag"));
		jets_pat_trackCountingHighPurBJetTag_branch->SetAddress(&jets_pat_trackCountingHighPurBJetTag_);
	}
	mus_pat_caloIso_branch = 0;
	if (tree->GetAlias("mus_pat_caloIso") != 0) {
		mus_pat_caloIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_caloIso"));
		mus_pat_caloIso_branch->SetAddress(&mus_pat_caloIso_);
	}
	mus_pat_calovetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_calovetoDep") != 0) {
		mus_pat_calovetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_calovetoDep"));
		mus_pat_calovetoDep_branch->SetAddress(&mus_pat_calovetoDep_);
	}
	mus_pat_ecalIso_branch = 0;
	if (tree->GetAlias("mus_pat_ecalIso") != 0) {
		mus_pat_ecalIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_ecalIso"));
		mus_pat_ecalIso_branch->SetAddress(&mus_pat_ecalIso_);
	}
	mus_pat_ecalvetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_ecalvetoDep") != 0) {
		mus_pat_ecalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_ecalvetoDep"));
		mus_pat_ecalvetoDep_branch->SetAddress(&mus_pat_ecalvetoDep_);
	}
	mus_pat_hcalIso_branch = 0;
	if (tree->GetAlias("mus_pat_hcalIso") != 0) {
		mus_pat_hcalIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_hcalIso"));
		mus_pat_hcalIso_branch->SetAddress(&mus_pat_hcalIso_);
	}
	mus_pat_hcalvetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_hcalvetoDep") != 0) {
		mus_pat_hcalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_hcalvetoDep"));
		mus_pat_hcalvetoDep_branch->SetAddress(&mus_pat_hcalvetoDep_);
	}
	mus_pat_trackIso_branch = 0;
	if (tree->GetAlias("mus_pat_trackIso") != 0) {
		mus_pat_trackIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_trackIso"));
		mus_pat_trackIso_branch->SetAddress(&mus_pat_trackIso_);
	}
	mus_pat_trckvetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_trckvetoDep") != 0) {
		mus_pat_trckvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_trckvetoDep"));
		mus_pat_trckvetoDep_branch->SetAddress(&mus_pat_trckvetoDep_);
	}
	pfels_deltaP_branch = 0;
	if (tree->GetAlias("pfels_deltaP") != 0) {
		pfels_deltaP_branch = tree->GetBranch(tree->GetAlias("pfels_deltaP"));
		pfels_deltaP_branch->SetAddress(&pfels_deltaP_);
	}
	pfels_ecalE_branch = 0;
	if (tree->GetAlias("pfels_ecalE") != 0) {
		pfels_ecalE_branch = tree->GetBranch(tree->GetAlias("pfels_ecalE"));
		pfels_ecalE_branch->SetAddress(&pfels_ecalE_);
	}
	pfels_hcalE_branch = 0;
	if (tree->GetAlias("pfels_hcalE") != 0) {
		pfels_hcalE_branch = tree->GetBranch(tree->GetAlias("pfels_hcalE"));
		pfels_hcalE_branch->SetAddress(&pfels_hcalE_);
	}
	pfels_isoChargedHadrons_branch = 0;
	if (tree->GetAlias("pfels_isoChargedHadrons") != 0) {
		pfels_isoChargedHadrons_branch = tree->GetBranch(tree->GetAlias("pfels_isoChargedHadrons"));
		pfels_isoChargedHadrons_branch->SetAddress(&pfels_isoChargedHadrons_);
	}
	pfels_isoNeutralHadrons_branch = 0;
	if (tree->GetAlias("pfels_isoNeutralHadrons") != 0) {
		pfels_isoNeutralHadrons_branch = tree->GetBranch(tree->GetAlias("pfels_isoNeutralHadrons"));
		pfels_isoNeutralHadrons_branch->SetAddress(&pfels_isoNeutralHadrons_);
	}
	pfels_isoPhotons_branch = 0;
	if (tree->GetAlias("pfels_isoPhotons") != 0) {
		pfels_isoPhotons_branch = tree->GetBranch(tree->GetAlias("pfels_isoPhotons"));
		pfels_isoPhotons_branch->SetAddress(&pfels_isoPhotons_);
	}
	pfels_mva_emu_branch = 0;
	if (tree->GetAlias("pfels_mva_emu") != 0) {
		pfels_mva_emu_branch = tree->GetBranch(tree->GetAlias("pfels_mva_emu"));
		pfels_mva_emu_branch->SetAddress(&pfels_mva_emu_);
	}
	pfels_mva_epi_branch = 0;
	if (tree->GetAlias("pfels_mva_epi") != 0) {
		pfels_mva_epi_branch = tree->GetBranch(tree->GetAlias("pfels_mva_epi"));
		pfels_mva_epi_branch->SetAddress(&pfels_mva_epi_);
	}
	pfels_mva_nothing_gamma_branch = 0;
	if (tree->GetAlias("pfels_mva_nothing_gamma") != 0) {
		pfels_mva_nothing_gamma_branch = tree->GetBranch(tree->GetAlias("pfels_mva_nothing_gamma"));
		pfels_mva_nothing_gamma_branch->SetAddress(&pfels_mva_nothing_gamma_);
	}
	pfels_mva_nothing_nh_branch = 0;
	if (tree->GetAlias("pfels_mva_nothing_nh") != 0) {
		pfels_mva_nothing_nh_branch = tree->GetBranch(tree->GetAlias("pfels_mva_nothing_nh"));
		pfels_mva_nothing_nh_branch->SetAddress(&pfels_mva_nothing_nh_);
	}
	pfels_mva_pimu_branch = 0;
	if (tree->GetAlias("pfels_mva_pimu") != 0) {
		pfels_mva_pimu_branch = tree->GetBranch(tree->GetAlias("pfels_mva_pimu"));
		pfels_mva_pimu_branch->SetAddress(&pfels_mva_pimu_);
	}
	pfels_pS1E_branch = 0;
	if (tree->GetAlias("pfels_pS1E") != 0) {
		pfels_pS1E_branch = tree->GetBranch(tree->GetAlias("pfels_pS1E"));
		pfels_pS1E_branch->SetAddress(&pfels_pS1E_);
	}
	pfels_pS2E_branch = 0;
	if (tree->GetAlias("pfels_pS2E") != 0) {
		pfels_pS2E_branch = tree->GetBranch(tree->GetAlias("pfels_pS2E"));
		pfels_pS2E_branch->SetAddress(&pfels_pS2E_);
	}
	pfels_rawEcalE_branch = 0;
	if (tree->GetAlias("pfels_rawEcalE") != 0) {
		pfels_rawEcalE_branch = tree->GetBranch(tree->GetAlias("pfels_rawEcalE"));
		pfels_rawEcalE_branch->SetAddress(&pfels_rawEcalE_);
	}
	pfels_rawHcalE_branch = 0;
	if (tree->GetAlias("pfels_rawHcalE") != 0) {
		pfels_rawHcalE_branch = tree->GetBranch(tree->GetAlias("pfels_rawHcalE"));
		pfels_rawHcalE_branch->SetAddress(&pfels_rawHcalE_);
	}
	pfjets_chargedEmE_branch = 0;
	if (tree->GetAlias("pfjets_chargedEmE") != 0) {
		pfjets_chargedEmE_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedEmE"));
		pfjets_chargedEmE_branch->SetAddress(&pfjets_chargedEmE_);
	}
	pfjets_chargedHadronE_branch = 0;
	if (tree->GetAlias("pfjets_chargedHadronE") != 0) {
		pfjets_chargedHadronE_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedHadronE"));
		pfjets_chargedHadronE_branch->SetAddress(&pfjets_chargedHadronE_);
	}
	pfjets_cor_branch = 0;
	if (tree->GetAlias("pfjets_cor") != 0) {
		pfjets_cor_branch = tree->GetBranch(tree->GetAlias("pfjets_cor"));
		pfjets_cor_branch->SetAddress(&pfjets_cor_);
	}
	pfjets_neutralEmE_branch = 0;
	if (tree->GetAlias("pfjets_neutralEmE") != 0) {
		pfjets_neutralEmE_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralEmE"));
		pfjets_neutralEmE_branch->SetAddress(&pfjets_neutralEmE_);
	}
	pfjets_neutralHadronE_branch = 0;
	if (tree->GetAlias("pfjets_neutralHadronE") != 0) {
		pfjets_neutralHadronE_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralHadronE"));
		pfjets_neutralHadronE_branch->SetAddress(&pfjets_neutralHadronE_);
	}
	pfmus_deltaP_branch = 0;
	if (tree->GetAlias("pfmus_deltaP") != 0) {
		pfmus_deltaP_branch = tree->GetBranch(tree->GetAlias("pfmus_deltaP"));
		pfmus_deltaP_branch->SetAddress(&pfmus_deltaP_);
	}
	pfmus_ecalE_branch = 0;
	if (tree->GetAlias("pfmus_ecalE") != 0) {
		pfmus_ecalE_branch = tree->GetBranch(tree->GetAlias("pfmus_ecalE"));
		pfmus_ecalE_branch->SetAddress(&pfmus_ecalE_);
	}
	pfmus_hcalE_branch = 0;
	if (tree->GetAlias("pfmus_hcalE") != 0) {
		pfmus_hcalE_branch = tree->GetBranch(tree->GetAlias("pfmus_hcalE"));
		pfmus_hcalE_branch->SetAddress(&pfmus_hcalE_);
	}
	pfmus_isoChargedHadrons_branch = 0;
	if (tree->GetAlias("pfmus_isoChargedHadrons") != 0) {
		pfmus_isoChargedHadrons_branch = tree->GetBranch(tree->GetAlias("pfmus_isoChargedHadrons"));
		pfmus_isoChargedHadrons_branch->SetAddress(&pfmus_isoChargedHadrons_);
	}
	pfmus_isoNeutralHadrons_branch = 0;
	if (tree->GetAlias("pfmus_isoNeutralHadrons") != 0) {
		pfmus_isoNeutralHadrons_branch = tree->GetBranch(tree->GetAlias("pfmus_isoNeutralHadrons"));
		pfmus_isoNeutralHadrons_branch->SetAddress(&pfmus_isoNeutralHadrons_);
	}
	pfmus_isoPhotons_branch = 0;
	if (tree->GetAlias("pfmus_isoPhotons") != 0) {
		pfmus_isoPhotons_branch = tree->GetBranch(tree->GetAlias("pfmus_isoPhotons"));
		pfmus_isoPhotons_branch->SetAddress(&pfmus_isoPhotons_);
	}
	pfmus_mva_emu_branch = 0;
	if (tree->GetAlias("pfmus_mva_emu") != 0) {
		pfmus_mva_emu_branch = tree->GetBranch(tree->GetAlias("pfmus_mva_emu"));
		pfmus_mva_emu_branch->SetAddress(&pfmus_mva_emu_);
	}
	pfmus_mva_epi_branch = 0;
	if (tree->GetAlias("pfmus_mva_epi") != 0) {
		pfmus_mva_epi_branch = tree->GetBranch(tree->GetAlias("pfmus_mva_epi"));
		pfmus_mva_epi_branch->SetAddress(&pfmus_mva_epi_);
	}
	pfmus_mva_nothing_gamma_branch = 0;
	if (tree->GetAlias("pfmus_mva_nothing_gamma") != 0) {
		pfmus_mva_nothing_gamma_branch = tree->GetBranch(tree->GetAlias("pfmus_mva_nothing_gamma"));
		pfmus_mva_nothing_gamma_branch->SetAddress(&pfmus_mva_nothing_gamma_);
	}
	pfmus_mva_nothing_nh_branch = 0;
	if (tree->GetAlias("pfmus_mva_nothing_nh") != 0) {
		pfmus_mva_nothing_nh_branch = tree->GetBranch(tree->GetAlias("pfmus_mva_nothing_nh"));
		pfmus_mva_nothing_nh_branch->SetAddress(&pfmus_mva_nothing_nh_);
	}
	pfmus_mva_pimu_branch = 0;
	if (tree->GetAlias("pfmus_mva_pimu") != 0) {
		pfmus_mva_pimu_branch = tree->GetBranch(tree->GetAlias("pfmus_mva_pimu"));
		pfmus_mva_pimu_branch->SetAddress(&pfmus_mva_pimu_);
	}
	pfmus_pS1E_branch = 0;
	if (tree->GetAlias("pfmus_pS1E") != 0) {
		pfmus_pS1E_branch = tree->GetBranch(tree->GetAlias("pfmus_pS1E"));
		pfmus_pS1E_branch->SetAddress(&pfmus_pS1E_);
	}
	pfmus_pS2E_branch = 0;
	if (tree->GetAlias("pfmus_pS2E") != 0) {
		pfmus_pS2E_branch = tree->GetBranch(tree->GetAlias("pfmus_pS2E"));
		pfmus_pS2E_branch->SetAddress(&pfmus_pS2E_);
	}
	pfmus_rawEcalE_branch = 0;
	if (tree->GetAlias("pfmus_rawEcalE") != 0) {
		pfmus_rawEcalE_branch = tree->GetBranch(tree->GetAlias("pfmus_rawEcalE"));
		pfmus_rawEcalE_branch->SetAddress(&pfmus_rawEcalE_);
	}
	pfmus_rawHcalE_branch = 0;
	if (tree->GetAlias("pfmus_rawHcalE") != 0) {
		pfmus_rawHcalE_branch = tree->GetBranch(tree->GetAlias("pfmus_rawHcalE"));
		pfmus_rawHcalE_branch->SetAddress(&pfmus_rawHcalE_);
	}
	photons_e1x5_branch = 0;
	if (tree->GetAlias("photons_e1x5") != 0) {
		photons_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_e1x5"));
		photons_e1x5_branch->SetAddress(&photons_e1x5_);
	}
	photons_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_e2x5Max") != 0) {
		photons_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_e2x5Max"));
		photons_e2x5Max_branch->SetAddress(&photons_e2x5Max_);
	}
	photons_e3x3_branch = 0;
	if (tree->GetAlias("photons_e3x3") != 0) {
		photons_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_e3x3"));
		photons_e3x3_branch->SetAddress(&photons_e3x3_);
	}
	photons_e5x5_branch = 0;
	if (tree->GetAlias("photons_e5x5") != 0) {
		photons_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_e5x5"));
		photons_e5x5_branch->SetAddress(&photons_e5x5_);
	}
	photons_ecalIso03_branch = 0;
	if (tree->GetAlias("photons_ecalIso03") != 0) {
		photons_ecalIso03_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso03"));
		photons_ecalIso03_branch->SetAddress(&photons_ecalIso03_);
	}
	photons_ecalIso04_branch = 0;
	if (tree->GetAlias("photons_ecalIso04") != 0) {
		photons_ecalIso04_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso04"));
		photons_ecalIso04_branch->SetAddress(&photons_ecalIso04_);
	}
	photons_hOverE_branch = 0;
	if (tree->GetAlias("photons_hOverE") != 0) {
		photons_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_hOverE"));
		photons_hOverE_branch->SetAddress(&photons_hOverE_);
	}
	photons_hcalIso03_branch = 0;
	if (tree->GetAlias("photons_hcalIso03") != 0) {
		photons_hcalIso03_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso03"));
		photons_hcalIso03_branch->SetAddress(&photons_hcalIso03_);
	}
	photons_hcalIso04_branch = 0;
	if (tree->GetAlias("photons_hcalIso04") != 0) {
		photons_hcalIso04_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso04"));
		photons_hcalIso04_branch->SetAddress(&photons_hcalIso04_);
	}
	photons_ntkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow03") != 0) {
		photons_ntkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow03"));
		photons_ntkIsoHollow03_branch->SetAddress(&photons_ntkIsoHollow03_);
	}
	photons_ntkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoHollow04") != 0) {
		photons_ntkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoHollow04"));
		photons_ntkIsoHollow04_branch->SetAddress(&photons_ntkIsoHollow04_);
	}
	photons_ntkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid03") != 0) {
		photons_ntkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid03"));
		photons_ntkIsoSolid03_branch->SetAddress(&photons_ntkIsoSolid03_);
	}
	photons_ntkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_ntkIsoSolid04") != 0) {
		photons_ntkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_ntkIsoSolid04"));
		photons_ntkIsoSolid04_branch->SetAddress(&photons_ntkIsoSolid04_);
	}
	photons_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_sigmaEtaEta") != 0) {
		photons_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaEtaEta"));
		photons_sigmaEtaEta_branch->SetAddress(&photons_sigmaEtaEta_);
	}
	photons_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_sigmaIEtaIEta") != 0) {
		photons_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaIEtaIEta"));
		photons_sigmaIEtaIEta_branch->SetAddress(&photons_sigmaIEtaIEta_);
	}
	photons_swissSeed_branch = 0;
	if (tree->GetAlias("photons_swissSeed") != 0) {
		photons_swissSeed_branch = tree->GetBranch(tree->GetAlias("photons_swissSeed"));
		photons_swissSeed_branch->SetAddress(&photons_swissSeed_);
	}
	photons_tkIsoHollow03_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow03") != 0) {
		photons_tkIsoHollow03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow03"));
		photons_tkIsoHollow03_branch->SetAddress(&photons_tkIsoHollow03_);
	}
	photons_tkIsoHollow04_branch = 0;
	if (tree->GetAlias("photons_tkIsoHollow04") != 0) {
		photons_tkIsoHollow04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoHollow04"));
		photons_tkIsoHollow04_branch->SetAddress(&photons_tkIsoHollow04_);
	}
	photons_tkIsoSolid03_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid03") != 0) {
		photons_tkIsoSolid03_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid03"));
		photons_tkIsoSolid03_branch->SetAddress(&photons_tkIsoSolid03_);
	}
	photons_tkIsoSolid04_branch = 0;
	if (tree->GetAlias("photons_tkIsoSolid04") != 0) {
		photons_tkIsoSolid04_branch = tree->GetBranch(tree->GetAlias("photons_tkIsoSolid04"));
		photons_tkIsoSolid04_branch->SetAddress(&photons_tkIsoSolid04_);
	}
	scs_clustersSize_branch = 0;
	if (tree->GetAlias("scs_clustersSize") != 0) {
		scs_clustersSize_branch = tree->GetBranch(tree->GetAlias("scs_clustersSize"));
		scs_clustersSize_branch->SetAddress(&scs_clustersSize_);
	}
	scs_crystalsSize_branch = 0;
	if (tree->GetAlias("scs_crystalsSize") != 0) {
		scs_crystalsSize_branch = tree->GetBranch(tree->GetAlias("scs_crystalsSize"));
		scs_crystalsSize_branch->SetAddress(&scs_crystalsSize_);
	}
	scs_e1x3_branch = 0;
	if (tree->GetAlias("scs_e1x3") != 0) {
		scs_e1x3_branch = tree->GetBranch(tree->GetAlias("scs_e1x3"));
		scs_e1x3_branch->SetAddress(&scs_e1x3_);
	}
	scs_e1x5_branch = 0;
	if (tree->GetAlias("scs_e1x5") != 0) {
		scs_e1x5_branch = tree->GetBranch(tree->GetAlias("scs_e1x5"));
		scs_e1x5_branch->SetAddress(&scs_e1x5_);
	}
	scs_e2nd_branch = 0;
	if (tree->GetAlias("scs_e2nd") != 0) {
		scs_e2nd_branch = tree->GetBranch(tree->GetAlias("scs_e2nd"));
		scs_e2nd_branch->SetAddress(&scs_e2nd_);
	}
	scs_e2x2_branch = 0;
	if (tree->GetAlias("scs_e2x2") != 0) {
		scs_e2x2_branch = tree->GetBranch(tree->GetAlias("scs_e2x2"));
		scs_e2x2_branch->SetAddress(&scs_e2x2_);
	}
	scs_e2x5Max_branch = 0;
	if (tree->GetAlias("scs_e2x5Max") != 0) {
		scs_e2x5Max_branch = tree->GetBranch(tree->GetAlias("scs_e2x5Max"));
		scs_e2x5Max_branch->SetAddress(&scs_e2x5Max_);
	}
	scs_e3x1_branch = 0;
	if (tree->GetAlias("scs_e3x1") != 0) {
		scs_e3x1_branch = tree->GetBranch(tree->GetAlias("scs_e3x1"));
		scs_e3x1_branch->SetAddress(&scs_e3x1_);
	}
	scs_e3x2_branch = 0;
	if (tree->GetAlias("scs_e3x2") != 0) {
		scs_e3x2_branch = tree->GetBranch(tree->GetAlias("scs_e3x2"));
		scs_e3x2_branch->SetAddress(&scs_e3x2_);
	}
	scs_e3x3_branch = 0;
	if (tree->GetAlias("scs_e3x3") != 0) {
		scs_e3x3_branch = tree->GetBranch(tree->GetAlias("scs_e3x3"));
		scs_e3x3_branch->SetAddress(&scs_e3x3_);
	}
	scs_e4x4_branch = 0;
	if (tree->GetAlias("scs_e4x4") != 0) {
		scs_e4x4_branch = tree->GetBranch(tree->GetAlias("scs_e4x4"));
		scs_e4x4_branch->SetAddress(&scs_e4x4_);
	}
	scs_e5x5_branch = 0;
	if (tree->GetAlias("scs_e5x5") != 0) {
		scs_e5x5_branch = tree->GetBranch(tree->GetAlias("scs_e5x5"));
		scs_e5x5_branch->SetAddress(&scs_e5x5_);
	}
	scs_eMax_branch = 0;
	if (tree->GetAlias("scs_eMax") != 0) {
		scs_eMax_branch = tree->GetBranch(tree->GetAlias("scs_eMax"));
		scs_eMax_branch->SetAddress(&scs_eMax_);
	}
	scs_eSeed_branch = 0;
	if (tree->GetAlias("scs_eSeed") != 0) {
		scs_eSeed_branch = tree->GetBranch(tree->GetAlias("scs_eSeed"));
		scs_eSeed_branch->SetAddress(&scs_eSeed_);
	}
	scs_energy_branch = 0;
	if (tree->GetAlias("scs_energy") != 0) {
		scs_energy_branch = tree->GetBranch(tree->GetAlias("scs_energy"));
		scs_energy_branch->SetAddress(&scs_energy_);
	}
	scs_eta_branch = 0;
	if (tree->GetAlias("scs_eta") != 0) {
		scs_eta_branch = tree->GetBranch(tree->GetAlias("scs_eta"));
		scs_eta_branch->SetAddress(&scs_eta_);
	}
	scs_hoe_branch = 0;
	if (tree->GetAlias("scs_hoe") != 0) {
		scs_hoe_branch = tree->GetBranch(tree->GetAlias("scs_hoe"));
		scs_hoe_branch->SetAddress(&scs_hoe_);
	}
	scs_phi_branch = 0;
	if (tree->GetAlias("scs_phi") != 0) {
		scs_phi_branch = tree->GetBranch(tree->GetAlias("scs_phi"));
		scs_phi_branch->SetAddress(&scs_phi_);
	}
	scs_preshowerEnergy_branch = 0;
	if (tree->GetAlias("scs_preshowerEnergy") != 0) {
		scs_preshowerEnergy_branch = tree->GetBranch(tree->GetAlias("scs_preshowerEnergy"));
		scs_preshowerEnergy_branch->SetAddress(&scs_preshowerEnergy_);
	}
	scs_rawEnergy_branch = 0;
	if (tree->GetAlias("scs_rawEnergy") != 0) {
		scs_rawEnergy_branch = tree->GetBranch(tree->GetAlias("scs_rawEnergy"));
		scs_rawEnergy_branch->SetAddress(&scs_rawEnergy_);
	}
	scs_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("scs_sigmaEtaEta") != 0) {
		scs_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("scs_sigmaEtaEta"));
		scs_sigmaEtaEta_branch->SetAddress(&scs_sigmaEtaEta_);
	}
	scs_sigmaEtaPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaEtaPhi") != 0) {
		scs_sigmaEtaPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaEtaPhi"));
		scs_sigmaEtaPhi_branch->SetAddress(&scs_sigmaEtaPhi_);
	}
	scs_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIEta") != 0) {
		scs_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIEta"));
		scs_sigmaIEtaIEta_branch->SetAddress(&scs_sigmaIEtaIEta_);
	}
	scs_sigmaIEtaIEtaSC_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIEtaSC") != 0) {
		scs_sigmaIEtaIEtaSC_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIEtaSC"));
		scs_sigmaIEtaIEtaSC_branch->SetAddress(&scs_sigmaIEtaIEtaSC_);
	}
	scs_sigmaIEtaIPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIPhi") != 0) {
		scs_sigmaIEtaIPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIPhi"));
		scs_sigmaIEtaIPhi_branch->SetAddress(&scs_sigmaIEtaIPhi_);
	}
	scs_sigmaIEtaIPhiSC_branch = 0;
	if (tree->GetAlias("scs_sigmaIEtaIPhiSC") != 0) {
		scs_sigmaIEtaIPhiSC_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIEtaIPhiSC"));
		scs_sigmaIEtaIPhiSC_branch->SetAddress(&scs_sigmaIEtaIPhiSC_);
	}
	scs_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaIPhiIPhi") != 0) {
		scs_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIPhiIPhi"));
		scs_sigmaIPhiIPhi_branch->SetAddress(&scs_sigmaIPhiIPhi_);
	}
	scs_sigmaIPhiIPhiSC_branch = 0;
	if (tree->GetAlias("scs_sigmaIPhiIPhiSC") != 0) {
		scs_sigmaIPhiIPhiSC_branch = tree->GetBranch(tree->GetAlias("scs_sigmaIPhiIPhiSC"));
		scs_sigmaIPhiIPhiSC_branch->SetAddress(&scs_sigmaIPhiIPhiSC_);
	}
	scs_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("scs_sigmaPhiPhi") != 0) {
		scs_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("scs_sigmaPhiPhi"));
		scs_sigmaPhiPhi_branch->SetAddress(&scs_sigmaPhiPhi_);
	}
	scs_timeSeed_branch = 0;
	if (tree->GetAlias("scs_timeSeed") != 0) {
		scs_timeSeed_branch = tree->GetBranch(tree->GetAlias("scs_timeSeed"));
		scs_timeSeed_branch->SetAddress(&scs_timeSeed_);
	}
	svs_anglePV_branch = 0;
	if (tree->GetAlias("svs_anglePV") != 0) {
		svs_anglePV_branch = tree->GetBranch(tree->GetAlias("svs_anglePV"));
		svs_anglePV_branch->SetAddress(&svs_anglePV_);
	}
	svs_chi2_branch = 0;
	if (tree->GetAlias("svs_chi2") != 0) {
		svs_chi2_branch = tree->GetBranch(tree->GetAlias("svs_chi2"));
		svs_chi2_branch->SetAddress(&svs_chi2_);
	}
	svs_dist3Dsig_branch = 0;
	if (tree->GetAlias("svs_dist3Dsig") != 0) {
		svs_dist3Dsig_branch = tree->GetBranch(tree->GetAlias("svs_dist3Dsig"));
		svs_dist3Dsig_branch->SetAddress(&svs_dist3Dsig_);
	}
	svs_dist3Dval_branch = 0;
	if (tree->GetAlias("svs_dist3Dval") != 0) {
		svs_dist3Dval_branch = tree->GetBranch(tree->GetAlias("svs_dist3Dval"));
		svs_dist3Dval_branch->SetAddress(&svs_dist3Dval_);
	}
	svs_distXYsig_branch = 0;
	if (tree->GetAlias("svs_distXYsig") != 0) {
		svs_distXYsig_branch = tree->GetBranch(tree->GetAlias("svs_distXYsig"));
		svs_distXYsig_branch->SetAddress(&svs_distXYsig_);
	}
	svs_distXYval_branch = 0;
	if (tree->GetAlias("svs_distXYval") != 0) {
		svs_distXYval_branch = tree->GetBranch(tree->GetAlias("svs_distXYval"));
		svs_distXYval_branch->SetAddress(&svs_distXYval_);
	}
	svs_ndof_branch = 0;
	if (tree->GetAlias("svs_ndof") != 0) {
		svs_ndof_branch = tree->GetBranch(tree->GetAlias("svs_ndof"));
		svs_ndof_branch->SetAddress(&svs_ndof_);
	}
	svs_prob_branch = 0;
	if (tree->GetAlias("svs_prob") != 0) {
		svs_prob_branch = tree->GetBranch(tree->GetAlias("svs_prob"));
		svs_prob_branch->SetAddress(&svs_prob_);
	}
	svs_xError_branch = 0;
	if (tree->GetAlias("svs_xError") != 0) {
		svs_xError_branch = tree->GetBranch(tree->GetAlias("svs_xError"));
		svs_xError_branch->SetAddress(&svs_xError_);
	}
	svs_yError_branch = 0;
	if (tree->GetAlias("svs_yError") != 0) {
		svs_yError_branch = tree->GetBranch(tree->GetAlias("svs_yError"));
		svs_yError_branch->SetAddress(&svs_yError_);
	}
	svs_zError_branch = 0;
	if (tree->GetAlias("svs_zError") != 0) {
		svs_zError_branch = tree->GetBranch(tree->GetAlias("svs_zError"));
		svs_zError_branch->SetAddress(&svs_zError_);
	}
	mus_tcmet_deltax_branch = 0;
	if (tree->GetAlias("mus_tcmet_deltax") != 0) {
		mus_tcmet_deltax_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_deltax"));
		mus_tcmet_deltax_branch->SetAddress(&mus_tcmet_deltax_);
	}
	mus_tcmet_deltay_branch = 0;
	if (tree->GetAlias("mus_tcmet_deltay") != 0) {
		mus_tcmet_deltay_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_deltay"));
		mus_tcmet_deltay_branch->SetAddress(&mus_tcmet_deltay_);
	}
	trks_chi2_branch = 0;
	if (tree->GetAlias("trks_chi2") != 0) {
		trks_chi2_branch = tree->GetBranch(tree->GetAlias("trks_chi2"));
		trks_chi2_branch->SetAddress(&trks_chi2_);
	}
	trks_d0_branch = 0;
	if (tree->GetAlias("trks_d0") != 0) {
		trks_d0_branch = tree->GetBranch(tree->GetAlias("trks_d0"));
		trks_d0_branch->SetAddress(&trks_d0_);
	}
	trks_d0Err_branch = 0;
	if (tree->GetAlias("trks_d0Err") != 0) {
		trks_d0Err_branch = tree->GetBranch(tree->GetAlias("trks_d0Err"));
		trks_d0Err_branch->SetAddress(&trks_d0Err_);
	}
	trks_d0corr_branch = 0;
	if (tree->GetAlias("trks_d0corr") != 0) {
		trks_d0corr_branch = tree->GetBranch(tree->GetAlias("trks_d0corr"));
		trks_d0corr_branch->SetAddress(&trks_d0corr_);
	}
	trks_d0corrPhi_branch = 0;
	if (tree->GetAlias("trks_d0corrPhi") != 0) {
		trks_d0corrPhi_branch = tree->GetBranch(tree->GetAlias("trks_d0corrPhi"));
		trks_d0corrPhi_branch->SetAddress(&trks_d0corrPhi_);
	}
	trks_d0phiCov_branch = 0;
	if (tree->GetAlias("trks_d0phiCov") != 0) {
		trks_d0phiCov_branch = tree->GetBranch(tree->GetAlias("trks_d0phiCov"));
		trks_d0phiCov_branch->SetAddress(&trks_d0phiCov_);
	}
	trks_etaErr_branch = 0;
	if (tree->GetAlias("trks_etaErr") != 0) {
		trks_etaErr_branch = tree->GetBranch(tree->GetAlias("trks_etaErr"));
		trks_etaErr_branch->SetAddress(&trks_etaErr_);
	}
	trks_layer1_charge_branch = 0;
	if (tree->GetAlias("trks_layer1_charge") != 0) {
		trks_layer1_charge_branch = tree->GetBranch(tree->GetAlias("trks_layer1_charge"));
		trks_layer1_charge_branch->SetAddress(&trks_layer1_charge_);
	}
	trks_ndof_branch = 0;
	if (tree->GetAlias("trks_ndof") != 0) {
		trks_ndof_branch = tree->GetBranch(tree->GetAlias("trks_ndof"));
		trks_ndof_branch->SetAddress(&trks_ndof_);
	}
	trks_phiErr_branch = 0;
	if (tree->GetAlias("trks_phiErr") != 0) {
		trks_phiErr_branch = tree->GetBranch(tree->GetAlias("trks_phiErr"));
		trks_phiErr_branch->SetAddress(&trks_phiErr_);
	}
	trks_ptErr_branch = 0;
	if (tree->GetAlias("trks_ptErr") != 0) {
		trks_ptErr_branch = tree->GetBranch(tree->GetAlias("trks_ptErr"));
		trks_ptErr_branch->SetAddress(&trks_ptErr_);
	}
	trks_z0_branch = 0;
	if (tree->GetAlias("trks_z0") != 0) {
		trks_z0_branch = tree->GetBranch(tree->GetAlias("trks_z0"));
		trks_z0_branch->SetAddress(&trks_z0_);
	}
	trks_z0Err_branch = 0;
	if (tree->GetAlias("trks_z0Err") != 0) {
		trks_z0Err_branch = tree->GetBranch(tree->GetAlias("trks_z0Err"));
		trks_z0Err_branch->SetAddress(&trks_z0Err_);
	}
	trks_z0corr_branch = 0;
	if (tree->GetAlias("trks_z0corr") != 0) {
		trks_z0corr_branch = tree->GetBranch(tree->GetAlias("trks_z0corr"));
		trks_z0corr_branch->SetAddress(&trks_z0corr_);
	}
	trkjets_cor_branch = 0;
	if (tree->GetAlias("trkjets_cor") != 0) {
		trkjets_cor_branch = tree->GetBranch(tree->GetAlias("trkjets_cor"));
		trkjets_cor_branch->SetAddress(&trkjets_cor_);
	}
	trks_d0Errvtx_branch = 0;
	if (tree->GetAlias("trks_d0Errvtx") != 0) {
		trks_d0Errvtx_branch = tree->GetBranch(tree->GetAlias("trks_d0Errvtx"));
		trks_d0Errvtx_branch->SetAddress(&trks_d0Errvtx_);
	}
	trks_d0vtx_branch = 0;
	if (tree->GetAlias("trks_d0vtx") != 0) {
		trks_d0vtx_branch = tree->GetBranch(tree->GetAlias("trks_d0vtx"));
		trks_d0vtx_branch->SetAddress(&trks_d0vtx_);
	}
	vtxs_chi2_branch = 0;
	if (tree->GetAlias("vtxs_chi2") != 0) {
		vtxs_chi2_branch = tree->GetBranch(tree->GetAlias("vtxs_chi2"));
		vtxs_chi2_branch->SetAddress(&vtxs_chi2_);
	}
	vtxs_ndof_branch = 0;
	if (tree->GetAlias("vtxs_ndof") != 0) {
		vtxs_ndof_branch = tree->GetBranch(tree->GetAlias("vtxs_ndof"));
		vtxs_ndof_branch->SetAddress(&vtxs_ndof_);
	}
	vtxs_sumpt_branch = 0;
	if (tree->GetAlias("vtxs_sumpt") != 0) {
		vtxs_sumpt_branch = tree->GetBranch(tree->GetAlias("vtxs_sumpt"));
		vtxs_sumpt_branch->SetAddress(&vtxs_sumpt_);
	}
	vtxs_xError_branch = 0;
	if (tree->GetAlias("vtxs_xError") != 0) {
		vtxs_xError_branch = tree->GetBranch(tree->GetAlias("vtxs_xError"));
		vtxs_xError_branch->SetAddress(&vtxs_xError_);
	}
	vtxs_yError_branch = 0;
	if (tree->GetAlias("vtxs_yError") != 0) {
		vtxs_yError_branch = tree->GetBranch(tree->GetAlias("vtxs_yError"));
		vtxs_yError_branch->SetAddress(&vtxs_yError_);
	}
	vtxs_zError_branch = 0;
	if (tree->GetAlias("vtxs_zError") != 0) {
		vtxs_zError_branch = tree->GetBranch(tree->GetAlias("vtxs_zError"));
		vtxs_zError_branch->SetAddress(&vtxs_zError_);
	}
	vtxs_covMatrix_branch = 0;
	if (tree->GetAlias("vtxs_covMatrix") != 0) {
		vtxs_covMatrix_branch = tree->GetBranch(tree->GetAlias("vtxs_covMatrix"));
		vtxs_covMatrix_branch->SetAddress(&vtxs_covMatrix_);
	}
	evt_cscLooseHaloId_branch = 0;
	if (tree->GetAlias("evt_cscLooseHaloId") != 0) {
		evt_cscLooseHaloId_branch = tree->GetBranch(tree->GetAlias("evt_cscLooseHaloId"));
		evt_cscLooseHaloId_branch->SetAddress(&evt_cscLooseHaloId_);
	}
	evt_cscTightHaloId_branch = 0;
	if (tree->GetAlias("evt_cscTightHaloId") != 0) {
		evt_cscTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_cscTightHaloId"));
		evt_cscTightHaloId_branch->SetAddress(&evt_cscTightHaloId_);
	}
	evt_ecalLooseHaloId_branch = 0;
	if (tree->GetAlias("evt_ecalLooseHaloId") != 0) {
		evt_ecalLooseHaloId_branch = tree->GetBranch(tree->GetAlias("evt_ecalLooseHaloId"));
		evt_ecalLooseHaloId_branch->SetAddress(&evt_ecalLooseHaloId_);
	}
	evt_ecalTightHaloId_branch = 0;
	if (tree->GetAlias("evt_ecalTightHaloId") != 0) {
		evt_ecalTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_ecalTightHaloId"));
		evt_ecalTightHaloId_branch->SetAddress(&evt_ecalTightHaloId_);
	}
	evt_extremeTightHaloId_branch = 0;
	if (tree->GetAlias("evt_extremeTightHaloId") != 0) {
		evt_extremeTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_extremeTightHaloId"));
		evt_extremeTightHaloId_branch->SetAddress(&evt_extremeTightHaloId_);
	}
	evt_globalLooseHaloId_branch = 0;
	if (tree->GetAlias("evt_globalLooseHaloId") != 0) {
		evt_globalLooseHaloId_branch = tree->GetBranch(tree->GetAlias("evt_globalLooseHaloId"));
		evt_globalLooseHaloId_branch->SetAddress(&evt_globalLooseHaloId_);
	}
	evt_globalTightHaloId_branch = 0;
	if (tree->GetAlias("evt_globalTightHaloId") != 0) {
		evt_globalTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_globalTightHaloId"));
		evt_globalTightHaloId_branch->SetAddress(&evt_globalTightHaloId_);
	}
	evt_hcalLooseHaloId_branch = 0;
	if (tree->GetAlias("evt_hcalLooseHaloId") != 0) {
		evt_hcalLooseHaloId_branch = tree->GetBranch(tree->GetAlias("evt_hcalLooseHaloId"));
		evt_hcalLooseHaloId_branch->SetAddress(&evt_hcalLooseHaloId_);
	}
	evt_hcalTightHaloId_branch = 0;
	if (tree->GetAlias("evt_hcalTightHaloId") != 0) {
		evt_hcalTightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_hcalTightHaloId"));
		evt_hcalTightHaloId_branch->SetAddress(&evt_hcalTightHaloId_);
	}
	evt_looseHaloId_branch = 0;
	if (tree->GetAlias("evt_looseHaloId") != 0) {
		evt_looseHaloId_branch = tree->GetBranch(tree->GetAlias("evt_looseHaloId"));
		evt_looseHaloId_branch->SetAddress(&evt_looseHaloId_);
	}
	evt_nHaloLikeTracks_branch = 0;
	if (tree->GetAlias("evt_nHaloLikeTracks") != 0) {
		evt_nHaloLikeTracks_branch = tree->GetBranch(tree->GetAlias("evt_nHaloLikeTracks"));
		evt_nHaloLikeTracks_branch->SetAddress(&evt_nHaloLikeTracks_);
	}
	evt_nHaloTriggerCandidates_branch = 0;
	if (tree->GetAlias("evt_nHaloTriggerCandidates") != 0) {
		evt_nHaloTriggerCandidates_branch = tree->GetBranch(tree->GetAlias("evt_nHaloTriggerCandidates"));
		evt_nHaloTriggerCandidates_branch->SetAddress(&evt_nHaloTriggerCandidates_);
	}
	evt_tightHaloId_branch = 0;
	if (tree->GetAlias("evt_tightHaloId") != 0) {
		evt_tightHaloId_branch = tree->GetBranch(tree->GetAlias("evt_tightHaloId"));
		evt_tightHaloId_branch->SetAddress(&evt_tightHaloId_);
	}
	evt_bsType_branch = 0;
	if (tree->GetAlias("evt_bsType") != 0) {
		evt_bsType_branch = tree->GetBranch(tree->GetAlias("evt_bsType"));
		evt_bsType_branch->SetAddress(&evt_bsType_);
	}
	evt_bunchCrossing_branch = 0;
	if (tree->GetAlias("evt_bunchCrossing") != 0) {
		evt_bunchCrossing_branch = tree->GetBranch(tree->GetAlias("evt_bunchCrossing"));
		evt_bunchCrossing_branch->SetAddress(&evt_bunchCrossing_);
	}
	evt_experimentType_branch = 0;
	if (tree->GetAlias("evt_experimentType") != 0) {
		evt_experimentType_branch = tree->GetBranch(tree->GetAlias("evt_experimentType"));
		evt_experimentType_branch->SetAddress(&evt_experimentType_);
	}
	evt_isRealData_branch = 0;
	if (tree->GetAlias("evt_isRealData") != 0) {
		evt_isRealData_branch = tree->GetBranch(tree->GetAlias("evt_isRealData"));
		evt_isRealData_branch->SetAddress(&evt_isRealData_);
	}
	evt_orbitNumber_branch = 0;
	if (tree->GetAlias("evt_orbitNumber") != 0) {
		evt_orbitNumber_branch = tree->GetBranch(tree->GetAlias("evt_orbitNumber"));
		evt_orbitNumber_branch->SetAddress(&evt_orbitNumber_);
	}
	evt_storeNumber_branch = 0;
	if (tree->GetAlias("evt_storeNumber") != 0) {
		evt_storeNumber_branch = tree->GetBranch(tree->GetAlias("evt_storeNumber"));
		evt_storeNumber_branch->SetAddress(&evt_storeNumber_);
	}
	hcalnoise_maxHPDHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxHPDHits") != 0) {
		hcalnoise_maxHPDHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxHPDHits"));
		hcalnoise_maxHPDHits_branch->SetAddress(&hcalnoise_maxHPDHits_);
	}
	hcalnoise_maxRBXHits_branch = 0;
	if (tree->GetAlias("hcalnoise_maxRBXHits") != 0) {
		hcalnoise_maxRBXHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxRBXHits"));
		hcalnoise_maxRBXHits_branch->SetAddress(&hcalnoise_maxRBXHits_);
	}
	hcalnoise_maxZeros_branch = 0;
	if (tree->GetAlias("hcalnoise_maxZeros") != 0) {
		hcalnoise_maxZeros_branch = tree->GetBranch(tree->GetAlias("hcalnoise_maxZeros"));
		hcalnoise_maxZeros_branch->SetAddress(&hcalnoise_maxZeros_);
	}
	hcalnoise_noiseFilterStatus_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseFilterStatus") != 0) {
		hcalnoise_noiseFilterStatus_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseFilterStatus"));
		hcalnoise_noiseFilterStatus_branch->SetAddress(&hcalnoise_noiseFilterStatus_);
	}
	hcalnoise_noiseType_branch = 0;
	if (tree->GetAlias("hcalnoise_noiseType") != 0) {
		hcalnoise_noiseType_branch = tree->GetBranch(tree->GetAlias("hcalnoise_noiseType"));
		hcalnoise_noiseType_branch->SetAddress(&hcalnoise_noiseType_);
	}
	hcalnoise_num10GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num10GeVHits") != 0) {
		hcalnoise_num10GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num10GeVHits"));
		hcalnoise_num10GeVHits_branch->SetAddress(&hcalnoise_num10GeVHits_);
	}
	hcalnoise_num25GeVHits_branch = 0;
	if (tree->GetAlias("hcalnoise_num25GeVHits") != 0) {
		hcalnoise_num25GeVHits_branch = tree->GetBranch(tree->GetAlias("hcalnoise_num25GeVHits"));
		hcalnoise_num25GeVHits_branch->SetAddress(&hcalnoise_num25GeVHits_);
	}
	hcalnoise_numProblematicRBXs_branch = 0;
	if (tree->GetAlias("hcalnoise_numProblematicRBXs") != 0) {
		hcalnoise_numProblematicRBXs_branch = tree->GetBranch(tree->GetAlias("hcalnoise_numProblematicRBXs"));
		hcalnoise_numProblematicRBXs_branch->SetAddress(&hcalnoise_numProblematicRBXs_);
	}
	hcalnoise_passHighLevelNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passHighLevelNoiseFilter") != 0) {
		hcalnoise_passHighLevelNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passHighLevelNoiseFilter"));
		hcalnoise_passHighLevelNoiseFilter_branch->SetAddress(&hcalnoise_passHighLevelNoiseFilter_);
	}
	hcalnoise_passLooseNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passLooseNoiseFilter") != 0) {
		hcalnoise_passLooseNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passLooseNoiseFilter"));
		hcalnoise_passLooseNoiseFilter_branch->SetAddress(&hcalnoise_passLooseNoiseFilter_);
	}
	hcalnoise_passTightNoiseFilter_branch = 0;
	if (tree->GetAlias("hcalnoise_passTightNoiseFilter") != 0) {
		hcalnoise_passTightNoiseFilter_branch = tree->GetBranch(tree->GetAlias("hcalnoise_passTightNoiseFilter"));
		hcalnoise_passTightNoiseFilter_branch->SetAddress(&hcalnoise_passTightNoiseFilter_);
	}
	l1_nemiso_branch = 0;
	if (tree->GetAlias("l1_nemiso") != 0) {
		l1_nemiso_branch = tree->GetBranch(tree->GetAlias("l1_nemiso"));
		l1_nemiso_branch->SetAddress(&l1_nemiso_);
	}
	l1_nemnoiso_branch = 0;
	if (tree->GetAlias("l1_nemnoiso") != 0) {
		l1_nemnoiso_branch = tree->GetBranch(tree->GetAlias("l1_nemnoiso"));
		l1_nemnoiso_branch->SetAddress(&l1_nemnoiso_);
	}
	l1_njetsc_branch = 0;
	if (tree->GetAlias("l1_njetsc") != 0) {
		l1_njetsc_branch = tree->GetBranch(tree->GetAlias("l1_njetsc"));
		l1_njetsc_branch->SetAddress(&l1_njetsc_);
	}
	l1_njetsf_branch = 0;
	if (tree->GetAlias("l1_njetsf") != 0) {
		l1_njetsf_branch = tree->GetBranch(tree->GetAlias("l1_njetsf"));
		l1_njetsf_branch->SetAddress(&l1_njetsf_);
	}
	l1_njetst_branch = 0;
	if (tree->GetAlias("l1_njetst") != 0) {
		l1_njetst_branch = tree->GetBranch(tree->GetAlias("l1_njetst"));
		l1_njetst_branch->SetAddress(&l1_njetst_);
	}
	l1_nmus_branch = 0;
	if (tree->GetAlias("l1_nmus") != 0) {
		l1_nmus_branch = tree->GetBranch(tree->GetAlias("l1_nmus"));
		l1_nmus_branch->SetAddress(&l1_nmus_);
	}
	pdfinfo_id1_branch = 0;
	if (tree->GetAlias("pdfinfo_id1") != 0) {
		pdfinfo_id1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id1"));
		pdfinfo_id1_branch->SetAddress(&pdfinfo_id1_);
	}
	pdfinfo_id2_branch = 0;
	if (tree->GetAlias("pdfinfo_id2") != 0) {
		pdfinfo_id2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id2"));
		pdfinfo_id2_branch->SetAddress(&pdfinfo_id2_);
	}
	evt_ecaliPhiSuspects_branch = 0;
	if (tree->GetAlias("evt_ecaliPhiSuspects") != 0) {
		evt_ecaliPhiSuspects_branch = tree->GetBranch(tree->GetAlias("evt_ecaliPhiSuspects"));
		evt_ecaliPhiSuspects_branch->SetAddress(&evt_ecaliPhiSuspects_);
	}
	evt_globaliPhiSuspects_branch = 0;
	if (tree->GetAlias("evt_globaliPhiSuspects") != 0) {
		evt_globaliPhiSuspects_branch = tree->GetBranch(tree->GetAlias("evt_globaliPhiSuspects"));
		evt_globaliPhiSuspects_branch->SetAddress(&evt_globaliPhiSuspects_);
	}
	evt_hcaliPhiSuspects_branch = 0;
	if (tree->GetAlias("evt_hcaliPhiSuspects") != 0) {
		evt_hcaliPhiSuspects_branch = tree->GetBranch(tree->GetAlias("evt_hcaliPhiSuspects"));
		evt_hcaliPhiSuspects_branch->SetAddress(&evt_hcaliPhiSuspects_);
	}
	els_mc3_id_branch = 0;
	if (tree->GetAlias("els_mc3_id") != 0) {
		els_mc3_id_branch = tree->GetBranch(tree->GetAlias("els_mc3_id"));
		els_mc3_id_branch->SetAddress(&els_mc3_id_);
	}
	els_mc3idx_branch = 0;
	if (tree->GetAlias("els_mc3idx") != 0) {
		els_mc3idx_branch = tree->GetBranch(tree->GetAlias("els_mc3idx"));
		els_mc3idx_branch->SetAddress(&els_mc3idx_);
	}
	els_mc3_motherid_branch = 0;
	if (tree->GetAlias("els_mc3_motherid") != 0) {
		els_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc3_motherid"));
		els_mc3_motherid_branch->SetAddress(&els_mc3_motherid_);
	}
	els_mc3_motheridx_branch = 0;
	if (tree->GetAlias("els_mc3_motheridx") != 0) {
		els_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("els_mc3_motheridx"));
		els_mc3_motheridx_branch->SetAddress(&els_mc3_motheridx_);
	}
	els_mc_id_branch = 0;
	if (tree->GetAlias("els_mc_id") != 0) {
		els_mc_id_branch = tree->GetBranch(tree->GetAlias("els_mc_id"));
		els_mc_id_branch->SetAddress(&els_mc_id_);
	}
	els_mcidx_branch = 0;
	if (tree->GetAlias("els_mcidx") != 0) {
		els_mcidx_branch = tree->GetBranch(tree->GetAlias("els_mcidx"));
		els_mcidx_branch->SetAddress(&els_mcidx_);
	}
	els_mc_motherid_branch = 0;
	if (tree->GetAlias("els_mc_motherid") != 0) {
		els_mc_motherid_branch = tree->GetBranch(tree->GetAlias("els_mc_motherid"));
		els_mc_motherid_branch->SetAddress(&els_mc_motherid_);
	}
	jets_mc3_id_branch = 0;
	if (tree->GetAlias("jets_mc3_id") != 0) {
		jets_mc3_id_branch = tree->GetBranch(tree->GetAlias("jets_mc3_id"));
		jets_mc3_id_branch->SetAddress(&jets_mc3_id_);
	}
	jets_mc3idx_branch = 0;
	if (tree->GetAlias("jets_mc3idx") != 0) {
		jets_mc3idx_branch = tree->GetBranch(tree->GetAlias("jets_mc3idx"));
		jets_mc3idx_branch->SetAddress(&jets_mc3idx_);
	}
	jets_mc_gpidx_branch = 0;
	if (tree->GetAlias("jets_mc_gpidx") != 0) {
		jets_mc_gpidx_branch = tree->GetBranch(tree->GetAlias("jets_mc_gpidx"));
		jets_mc_gpidx_branch->SetAddress(&jets_mc_gpidx_);
	}
	jets_mc_id_branch = 0;
	if (tree->GetAlias("jets_mc_id") != 0) {
		jets_mc_id_branch = tree->GetBranch(tree->GetAlias("jets_mc_id"));
		jets_mc_id_branch->SetAddress(&jets_mc_id_);
	}
	jets_mcidx_branch = 0;
	if (tree->GetAlias("jets_mcidx") != 0) {
		jets_mcidx_branch = tree->GetBranch(tree->GetAlias("jets_mcidx"));
		jets_mcidx_branch->SetAddress(&jets_mcidx_);
	}
	jets_mc_motherid_branch = 0;
	if (tree->GetAlias("jets_mc_motherid") != 0) {
		jets_mc_motherid_branch = tree->GetBranch(tree->GetAlias("jets_mc_motherid"));
		jets_mc_motherid_branch->SetAddress(&jets_mc_motherid_);
	}
	mus_mc3_id_branch = 0;
	if (tree->GetAlias("mus_mc3_id") != 0) {
		mus_mc3_id_branch = tree->GetBranch(tree->GetAlias("mus_mc3_id"));
		mus_mc3_id_branch->SetAddress(&mus_mc3_id_);
	}
	mus_mc3idx_branch = 0;
	if (tree->GetAlias("mus_mc3idx") != 0) {
		mus_mc3idx_branch = tree->GetBranch(tree->GetAlias("mus_mc3idx"));
		mus_mc3idx_branch->SetAddress(&mus_mc3idx_);
	}
	mus_mc3_motherid_branch = 0;
	if (tree->GetAlias("mus_mc3_motherid") != 0) {
		mus_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motherid"));
		mus_mc3_motherid_branch->SetAddress(&mus_mc3_motherid_);
	}
	mus_mc3_motheridx_branch = 0;
	if (tree->GetAlias("mus_mc3_motheridx") != 0) {
		mus_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("mus_mc3_motheridx"));
		mus_mc3_motheridx_branch->SetAddress(&mus_mc3_motheridx_);
	}
	mus_mc_id_branch = 0;
	if (tree->GetAlias("mus_mc_id") != 0) {
		mus_mc_id_branch = tree->GetBranch(tree->GetAlias("mus_mc_id"));
		mus_mc_id_branch->SetAddress(&mus_mc_id_);
	}
	mus_mcidx_branch = 0;
	if (tree->GetAlias("mus_mcidx") != 0) {
		mus_mcidx_branch = tree->GetBranch(tree->GetAlias("mus_mcidx"));
		mus_mcidx_branch->SetAddress(&mus_mcidx_);
	}
	mus_mc_motherid_branch = 0;
	if (tree->GetAlias("mus_mc_motherid") != 0) {
		mus_mc_motherid_branch = tree->GetBranch(tree->GetAlias("mus_mc_motherid"));
		mus_mc_motherid_branch->SetAddress(&mus_mc_motherid_);
	}
	pfjets_mc3_id_branch = 0;
	if (tree->GetAlias("pfjets_mc3_id") != 0) {
		pfjets_mc3_id_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3_id"));
		pfjets_mc3_id_branch->SetAddress(&pfjets_mc3_id_);
	}
	pfjets_mc3idx_branch = 0;
	if (tree->GetAlias("pfjets_mc3idx") != 0) {
		pfjets_mc3idx_branch = tree->GetBranch(tree->GetAlias("pfjets_mc3idx"));
		pfjets_mc3idx_branch->SetAddress(&pfjets_mc3idx_);
	}
	pfjets_mc_gpidx_branch = 0;
	if (tree->GetAlias("pfjets_mc_gpidx") != 0) {
		pfjets_mc_gpidx_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_gpidx"));
		pfjets_mc_gpidx_branch->SetAddress(&pfjets_mc_gpidx_);
	}
	pfjets_mc_id_branch = 0;
	if (tree->GetAlias("pfjets_mc_id") != 0) {
		pfjets_mc_id_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_id"));
		pfjets_mc_id_branch->SetAddress(&pfjets_mc_id_);
	}
	pfjets_mcidx_branch = 0;
	if (tree->GetAlias("pfjets_mcidx") != 0) {
		pfjets_mcidx_branch = tree->GetBranch(tree->GetAlias("pfjets_mcidx"));
		pfjets_mcidx_branch->SetAddress(&pfjets_mcidx_);
	}
	pfjets_mc_motherid_branch = 0;
	if (tree->GetAlias("pfjets_mc_motherid") != 0) {
		pfjets_mc_motherid_branch = tree->GetBranch(tree->GetAlias("pfjets_mc_motherid"));
		pfjets_mc_motherid_branch->SetAddress(&pfjets_mc_motherid_);
	}
	photons_mc3_id_branch = 0;
	if (tree->GetAlias("photons_mc3_id") != 0) {
		photons_mc3_id_branch = tree->GetBranch(tree->GetAlias("photons_mc3_id"));
		photons_mc3_id_branch->SetAddress(&photons_mc3_id_);
	}
	photons_mc3idx_branch = 0;
	if (tree->GetAlias("photons_mc3idx") != 0) {
		photons_mc3idx_branch = tree->GetBranch(tree->GetAlias("photons_mc3idx"));
		photons_mc3idx_branch->SetAddress(&photons_mc3idx_);
	}
	photons_mc3_motherid_branch = 0;
	if (tree->GetAlias("photons_mc3_motherid") != 0) {
		photons_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("photons_mc3_motherid"));
		photons_mc3_motherid_branch->SetAddress(&photons_mc3_motherid_);
	}
	photons_mc3_motheridx_branch = 0;
	if (tree->GetAlias("photons_mc3_motheridx") != 0) {
		photons_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("photons_mc3_motheridx"));
		photons_mc3_motheridx_branch->SetAddress(&photons_mc3_motheridx_);
	}
	photons_mc_id_branch = 0;
	if (tree->GetAlias("photons_mc_id") != 0) {
		photons_mc_id_branch = tree->GetBranch(tree->GetAlias("photons_mc_id"));
		photons_mc_id_branch->SetAddress(&photons_mc_id_);
	}
	photons_mcidx_branch = 0;
	if (tree->GetAlias("photons_mcidx") != 0) {
		photons_mcidx_branch = tree->GetBranch(tree->GetAlias("photons_mcidx"));
		photons_mcidx_branch->SetAddress(&photons_mcidx_);
	}
	photons_mc_motherid_branch = 0;
	if (tree->GetAlias("photons_mc_motherid") != 0) {
		photons_mc_motherid_branch = tree->GetBranch(tree->GetAlias("photons_mc_motherid"));
		photons_mc_motherid_branch->SetAddress(&photons_mc_motherid_);
	}
	trk_mc3_id_branch = 0;
	if (tree->GetAlias("trk_mc3_id") != 0) {
		trk_mc3_id_branch = tree->GetBranch(tree->GetAlias("trk_mc3_id"));
		trk_mc3_id_branch->SetAddress(&trk_mc3_id_);
	}
	trk_mc3idx_branch = 0;
	if (tree->GetAlias("trk_mc3idx") != 0) {
		trk_mc3idx_branch = tree->GetBranch(tree->GetAlias("trk_mc3idx"));
		trk_mc3idx_branch->SetAddress(&trk_mc3idx_);
	}
	trk_mc3_motherid_branch = 0;
	if (tree->GetAlias("trk_mc3_motherid") != 0) {
		trk_mc3_motherid_branch = tree->GetBranch(tree->GetAlias("trk_mc3_motherid"));
		trk_mc3_motherid_branch->SetAddress(&trk_mc3_motherid_);
	}
	trk_mc3_motheridx_branch = 0;
	if (tree->GetAlias("trk_mc3_motheridx") != 0) {
		trk_mc3_motheridx_branch = tree->GetBranch(tree->GetAlias("trk_mc3_motheridx"));
		trk_mc3_motheridx_branch->SetAddress(&trk_mc3_motheridx_);
	}
	trk_mc_id_branch = 0;
	if (tree->GetAlias("trk_mc_id") != 0) {
		trk_mc_id_branch = tree->GetBranch(tree->GetAlias("trk_mc_id"));
		trk_mc_id_branch->SetAddress(&trk_mc_id_);
	}
	trk_mcidx_branch = 0;
	if (tree->GetAlias("trk_mcidx") != 0) {
		trk_mcidx_branch = tree->GetBranch(tree->GetAlias("trk_mcidx"));
		trk_mcidx_branch->SetAddress(&trk_mcidx_);
	}
	trk_mc_motherid_branch = 0;
	if (tree->GetAlias("trk_mc_motherid") != 0) {
		trk_mc_motherid_branch = tree->GetBranch(tree->GetAlias("trk_mc_motherid"));
		trk_mc_motherid_branch->SetAddress(&trk_mc_motherid_);
	}
	trks_conv_tkidx_branch = 0;
	if (tree->GetAlias("trks_conv_tkidx") != 0) {
		trks_conv_tkidx_branch = tree->GetBranch(tree->GetAlias("trks_conv_tkidx"));
		trks_conv_tkidx_branch->SetAddress(&trks_conv_tkidx_);
	}
	els_exp_innerlayers39X_branch = 0;
	if (tree->GetAlias("els_exp_innerlayers39X") != 0) {
		els_exp_innerlayers39X_branch = tree->GetBranch(tree->GetAlias("els_exp_innerlayers39X"));
		els_exp_innerlayers39X_branch->SetAddress(&els_exp_innerlayers39X_);
	}
	els_closestJet_branch = 0;
	if (tree->GetAlias("els_closestJet") != 0) {
		els_closestJet_branch = tree->GetBranch(tree->GetAlias("els_closestJet"));
		els_closestJet_branch->SetAddress(&els_closestJet_);
	}
	els_closestMuon_branch = 0;
	if (tree->GetAlias("els_closestMuon") != 0) {
		els_closestMuon_branch = tree->GetBranch(tree->GetAlias("els_closestMuon"));
		els_closestMuon_branch->SetAddress(&els_closestMuon_);
	}
	els_pfelsidx_branch = 0;
	if (tree->GetAlias("els_pfelsidx") != 0) {
		els_pfelsidx_branch = tree->GetBranch(tree->GetAlias("els_pfelsidx"));
		els_pfelsidx_branch->SetAddress(&els_pfelsidx_);
	}
	els_category_branch = 0;
	if (tree->GetAlias("els_category") != 0) {
		els_category_branch = tree->GetBranch(tree->GetAlias("els_category"));
		els_category_branch->SetAddress(&els_category_);
	}
	els_charge_branch = 0;
	if (tree->GetAlias("els_charge") != 0) {
		els_charge_branch = tree->GetBranch(tree->GetAlias("els_charge"));
		els_charge_branch->SetAddress(&els_charge_);
	}
	els_class_branch = 0;
	if (tree->GetAlias("els_class") != 0) {
		els_class_branch = tree->GetBranch(tree->GetAlias("els_class"));
		els_class_branch->SetAddress(&els_class_);
	}
	els_conv_tkidx_branch = 0;
	if (tree->GetAlias("els_conv_tkidx") != 0) {
		els_conv_tkidx_branch = tree->GetBranch(tree->GetAlias("els_conv_tkidx"));
		els_conv_tkidx_branch->SetAddress(&els_conv_tkidx_);
	}
	els_exp_innerlayers_branch = 0;
	if (tree->GetAlias("els_exp_innerlayers") != 0) {
		els_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_innerlayers"));
		els_exp_innerlayers_branch->SetAddress(&els_exp_innerlayers_);
	}
	els_exp_outerlayers_branch = 0;
	if (tree->GetAlias("els_exp_outerlayers") != 0) {
		els_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("els_exp_outerlayers"));
		els_exp_outerlayers_branch->SetAddress(&els_exp_outerlayers_);
	}
	els_fiduciality_branch = 0;
	if (tree->GetAlias("els_fiduciality") != 0) {
		els_fiduciality_branch = tree->GetBranch(tree->GetAlias("els_fiduciality"));
		els_fiduciality_branch->SetAddress(&els_fiduciality_);
	}
	els_gsftrkidx_branch = 0;
	if (tree->GetAlias("els_gsftrkidx") != 0) {
		els_gsftrkidx_branch = tree->GetBranch(tree->GetAlias("els_gsftrkidx"));
		els_gsftrkidx_branch->SetAddress(&els_gsftrkidx_);
	}
	els_layer1_det_branch = 0;
	if (tree->GetAlias("els_layer1_det") != 0) {
		els_layer1_det_branch = tree->GetBranch(tree->GetAlias("els_layer1_det"));
		els_layer1_det_branch->SetAddress(&els_layer1_det_);
	}
	els_layer1_layer_branch = 0;
	if (tree->GetAlias("els_layer1_layer") != 0) {
		els_layer1_layer_branch = tree->GetBranch(tree->GetAlias("els_layer1_layer"));
		els_layer1_layer_branch->SetAddress(&els_layer1_layer_);
	}
	els_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("els_layer1_sizerphi") != 0) {
		els_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerphi"));
		els_layer1_sizerphi_branch->SetAddress(&els_layer1_sizerphi_);
	}
	els_layer1_sizerz_branch = 0;
	if (tree->GetAlias("els_layer1_sizerz") != 0) {
		els_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("els_layer1_sizerz"));
		els_layer1_sizerz_branch->SetAddress(&els_layer1_sizerz_);
	}
	els_lostHits_branch = 0;
	if (tree->GetAlias("els_lostHits") != 0) {
		els_lostHits_branch = tree->GetBranch(tree->GetAlias("els_lostHits"));
		els_lostHits_branch->SetAddress(&els_lostHits_);
	}
	els_lost_pixelhits_branch = 0;
	if (tree->GetAlias("els_lost_pixelhits") != 0) {
		els_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_lost_pixelhits"));
		els_lost_pixelhits_branch->SetAddress(&els_lost_pixelhits_);
	}
	els_nSeed_branch = 0;
	if (tree->GetAlias("els_nSeed") != 0) {
		els_nSeed_branch = tree->GetBranch(tree->GetAlias("els_nSeed"));
		els_nSeed_branch->SetAddress(&els_nSeed_);
	}
	els_sccharge_branch = 0;
	if (tree->GetAlias("els_sccharge") != 0) {
		els_sccharge_branch = tree->GetBranch(tree->GetAlias("els_sccharge"));
		els_sccharge_branch->SetAddress(&els_sccharge_);
	}
	els_scindex_branch = 0;
	if (tree->GetAlias("els_scindex") != 0) {
		els_scindex_branch = tree->GetBranch(tree->GetAlias("els_scindex"));
		els_scindex_branch->SetAddress(&els_scindex_);
	}
	els_trk_charge_branch = 0;
	if (tree->GetAlias("els_trk_charge") != 0) {
		els_trk_charge_branch = tree->GetBranch(tree->GetAlias("els_trk_charge"));
		els_trk_charge_branch->SetAddress(&els_trk_charge_);
	}
	els_trkidx_branch = 0;
	if (tree->GetAlias("els_trkidx") != 0) {
		els_trkidx_branch = tree->GetBranch(tree->GetAlias("els_trkidx"));
		els_trkidx_branch->SetAddress(&els_trkidx_);
	}
	els_type_branch = 0;
	if (tree->GetAlias("els_type") != 0) {
		els_type_branch = tree->GetBranch(tree->GetAlias("els_type"));
		els_type_branch->SetAddress(&els_type_);
	}
	els_validHits_branch = 0;
	if (tree->GetAlias("els_validHits") != 0) {
		els_validHits_branch = tree->GetBranch(tree->GetAlias("els_validHits"));
		els_validHits_branch->SetAddress(&els_validHits_);
	}
	els_valid_pixelhits_branch = 0;
	if (tree->GetAlias("els_valid_pixelhits") != 0) {
		els_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("els_valid_pixelhits"));
		els_valid_pixelhits_branch->SetAddress(&els_valid_pixelhits_);
	}
	genps_id_branch = 0;
	if (tree->GetAlias("genps_id") != 0) {
		genps_id_branch = tree->GetBranch(tree->GetAlias("genps_id"));
		genps_id_branch->SetAddress(&genps_id_);
	}
	genps_id_mother_branch = 0;
	if (tree->GetAlias("genps_id_mother") != 0) {
		genps_id_mother_branch = tree->GetBranch(tree->GetAlias("genps_id_mother"));
		genps_id_mother_branch->SetAddress(&genps_id_mother_);
	}
	genps_status_branch = 0;
	if (tree->GetAlias("genps_status") != 0) {
		genps_status_branch = tree->GetBranch(tree->GetAlias("genps_status"));
		genps_status_branch->SetAddress(&genps_status_);
	}
	gsftrks_charge_branch = 0;
	if (tree->GetAlias("gsftrks_charge") != 0) {
		gsftrks_charge_branch = tree->GetBranch(tree->GetAlias("gsftrks_charge"));
		gsftrks_charge_branch->SetAddress(&gsftrks_charge_);
	}
	gsftrks_exp_innerlayers_branch = 0;
	if (tree->GetAlias("gsftrks_exp_innerlayers") != 0) {
		gsftrks_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("gsftrks_exp_innerlayers"));
		gsftrks_exp_innerlayers_branch->SetAddress(&gsftrks_exp_innerlayers_);
	}
	gsftrks_exp_outerlayers_branch = 0;
	if (tree->GetAlias("gsftrks_exp_outerlayers") != 0) {
		gsftrks_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("gsftrks_exp_outerlayers"));
		gsftrks_exp_outerlayers_branch->SetAddress(&gsftrks_exp_outerlayers_);
	}
	gsftrks_layer1_det_branch = 0;
	if (tree->GetAlias("gsftrks_layer1_det") != 0) {
		gsftrks_layer1_det_branch = tree->GetBranch(tree->GetAlias("gsftrks_layer1_det"));
		gsftrks_layer1_det_branch->SetAddress(&gsftrks_layer1_det_);
	}
	gsftrks_layer1_layer_branch = 0;
	if (tree->GetAlias("gsftrks_layer1_layer") != 0) {
		gsftrks_layer1_layer_branch = tree->GetBranch(tree->GetAlias("gsftrks_layer1_layer"));
		gsftrks_layer1_layer_branch->SetAddress(&gsftrks_layer1_layer_);
	}
	gsftrks_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("gsftrks_layer1_sizerphi") != 0) {
		gsftrks_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("gsftrks_layer1_sizerphi"));
		gsftrks_layer1_sizerphi_branch->SetAddress(&gsftrks_layer1_sizerphi_);
	}
	gsftrks_layer1_sizerz_branch = 0;
	if (tree->GetAlias("gsftrks_layer1_sizerz") != 0) {
		gsftrks_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("gsftrks_layer1_sizerz"));
		gsftrks_layer1_sizerz_branch->SetAddress(&gsftrks_layer1_sizerz_);
	}
	gsftrks_lostHits_branch = 0;
	if (tree->GetAlias("gsftrks_lostHits") != 0) {
		gsftrks_lostHits_branch = tree->GetBranch(tree->GetAlias("gsftrks_lostHits"));
		gsftrks_lostHits_branch->SetAddress(&gsftrks_lostHits_);
	}
	gsftrks_lost_pixelhits_branch = 0;
	if (tree->GetAlias("gsftrks_lost_pixelhits") != 0) {
		gsftrks_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("gsftrks_lost_pixelhits"));
		gsftrks_lost_pixelhits_branch->SetAddress(&gsftrks_lost_pixelhits_);
	}
	gsftrks_nlayers_branch = 0;
	if (tree->GetAlias("gsftrks_nlayers") != 0) {
		gsftrks_nlayers_branch = tree->GetBranch(tree->GetAlias("gsftrks_nlayers"));
		gsftrks_nlayers_branch->SetAddress(&gsftrks_nlayers_);
	}
	gsftrks_nlayers3D_branch = 0;
	if (tree->GetAlias("gsftrks_nlayers3D") != 0) {
		gsftrks_nlayers3D_branch = tree->GetBranch(tree->GetAlias("gsftrks_nlayers3D"));
		gsftrks_nlayers3D_branch->SetAddress(&gsftrks_nlayers3D_);
	}
	gsftrks_nlayersLost_branch = 0;
	if (tree->GetAlias("gsftrks_nlayersLost") != 0) {
		gsftrks_nlayersLost_branch = tree->GetBranch(tree->GetAlias("gsftrks_nlayersLost"));
		gsftrks_nlayersLost_branch->SetAddress(&gsftrks_nlayersLost_);
	}
	gsftrks_validHits_branch = 0;
	if (tree->GetAlias("gsftrks_validHits") != 0) {
		gsftrks_validHits_branch = tree->GetBranch(tree->GetAlias("gsftrks_validHits"));
		gsftrks_validHits_branch->SetAddress(&gsftrks_validHits_);
	}
	gsftrks_valid_pixelhits_branch = 0;
	if (tree->GetAlias("gsftrks_valid_pixelhits") != 0) {
		gsftrks_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("gsftrks_valid_pixelhits"));
		gsftrks_valid_pixelhits_branch->SetAddress(&gsftrks_valid_pixelhits_);
	}
	hyp_ll_charge_branch = 0;
	if (tree->GetAlias("hyp_ll_charge") != 0) {
		hyp_ll_charge_branch = tree->GetBranch(tree->GetAlias("hyp_ll_charge"));
		hyp_ll_charge_branch->SetAddress(&hyp_ll_charge_);
	}
	hyp_ll_id_branch = 0;
	if (tree->GetAlias("hyp_ll_id") != 0) {
		hyp_ll_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_id"));
		hyp_ll_id_branch->SetAddress(&hyp_ll_id_);
	}
	hyp_ll_index_branch = 0;
	if (tree->GetAlias("hyp_ll_index") != 0) {
		hyp_ll_index_branch = tree->GetBranch(tree->GetAlias("hyp_ll_index"));
		hyp_ll_index_branch->SetAddress(&hyp_ll_index_);
	}
	hyp_ll_lostHits_branch = 0;
	if (tree->GetAlias("hyp_ll_lostHits") != 0) {
		hyp_ll_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_lostHits"));
		hyp_ll_lostHits_branch->SetAddress(&hyp_ll_lostHits_);
	}
	hyp_ll_validHits_branch = 0;
	if (tree->GetAlias("hyp_ll_validHits") != 0) {
		hyp_ll_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_validHits"));
		hyp_ll_validHits_branch->SetAddress(&hyp_ll_validHits_);
	}
	hyp_lt_charge_branch = 0;
	if (tree->GetAlias("hyp_lt_charge") != 0) {
		hyp_lt_charge_branch = tree->GetBranch(tree->GetAlias("hyp_lt_charge"));
		hyp_lt_charge_branch->SetAddress(&hyp_lt_charge_);
	}
	hyp_lt_id_branch = 0;
	if (tree->GetAlias("hyp_lt_id") != 0) {
		hyp_lt_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_id"));
		hyp_lt_id_branch->SetAddress(&hyp_lt_id_);
	}
	hyp_lt_index_branch = 0;
	if (tree->GetAlias("hyp_lt_index") != 0) {
		hyp_lt_index_branch = tree->GetBranch(tree->GetAlias("hyp_lt_index"));
		hyp_lt_index_branch->SetAddress(&hyp_lt_index_);
	}
	hyp_lt_lostHits_branch = 0;
	if (tree->GetAlias("hyp_lt_lostHits") != 0) {
		hyp_lt_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_lostHits"));
		hyp_lt_lostHits_branch->SetAddress(&hyp_lt_lostHits_);
	}
	hyp_lt_validHits_branch = 0;
	if (tree->GetAlias("hyp_lt_validHits") != 0) {
		hyp_lt_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_validHits"));
		hyp_lt_validHits_branch->SetAddress(&hyp_lt_validHits_);
	}
	hyp_njets_branch = 0;
	if (tree->GetAlias("hyp_njets") != 0) {
		hyp_njets_branch = tree->GetBranch(tree->GetAlias("hyp_njets"));
		hyp_njets_branch->SetAddress(&hyp_njets_);
	}
	hyp_nojets_branch = 0;
	if (tree->GetAlias("hyp_nojets") != 0) {
		hyp_nojets_branch = tree->GetBranch(tree->GetAlias("hyp_nojets"));
		hyp_nojets_branch->SetAddress(&hyp_nojets_);
	}
	hyp_type_branch = 0;
	if (tree->GetAlias("hyp_type") != 0) {
		hyp_type_branch = tree->GetBranch(tree->GetAlias("hyp_type"));
		hyp_type_branch->SetAddress(&hyp_type_);
	}
	hyp_FVFit_ndf_branch = 0;
	if (tree->GetAlias("hyp_FVFit_ndf") != 0) {
		hyp_FVFit_ndf_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_ndf"));
		hyp_FVFit_ndf_branch->SetAddress(&hyp_FVFit_ndf_);
	}
	hyp_FVFit_status_branch = 0;
	if (tree->GetAlias("hyp_FVFit_status") != 0) {
		hyp_FVFit_status_branch = tree->GetBranch(tree->GetAlias("hyp_FVFit_status"));
		hyp_FVFit_status_branch->SetAddress(&hyp_FVFit_status_);
	}
	hyp_ll_mc_id_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_id") != 0) {
		hyp_ll_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_id"));
		hyp_ll_mc_id_branch->SetAddress(&hyp_ll_mc_id_);
	}
	hyp_ll_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_motherid") != 0) {
		hyp_ll_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_motherid"));
		hyp_ll_mc_motherid_branch->SetAddress(&hyp_ll_mc_motherid_);
	}
	hyp_lt_mc_id_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_id") != 0) {
		hyp_lt_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_id"));
		hyp_lt_mc_id_branch->SetAddress(&hyp_lt_mc_id_);
	}
	hyp_lt_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_motherid") != 0) {
		hyp_lt_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_motherid"));
		hyp_lt_mc_motherid_branch->SetAddress(&hyp_lt_mc_motherid_);
	}
	hyp_quadlep_first_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_type") != 0) {
		hyp_quadlep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_type"));
		hyp_quadlep_first_type_branch->SetAddress(&hyp_quadlep_first_type_);
	}
	hyp_quadlep_fourth_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_type") != 0) {
		hyp_quadlep_fourth_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_type"));
		hyp_quadlep_fourth_type_branch->SetAddress(&hyp_quadlep_fourth_type_);
	}
	hyp_quadlep_second_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_type") != 0) {
		hyp_quadlep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_type"));
		hyp_quadlep_second_type_branch->SetAddress(&hyp_quadlep_second_type_);
	}
	hyp_quadlep_third_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_type") != 0) {
		hyp_quadlep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_type"));
		hyp_quadlep_third_type_branch->SetAddress(&hyp_quadlep_third_type_);
	}
	hyp_trilep_first_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_type") != 0) {
		hyp_trilep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_type"));
		hyp_trilep_first_type_branch->SetAddress(&hyp_trilep_first_type_);
	}
	hyp_trilep_second_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_type") != 0) {
		hyp_trilep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_type"));
		hyp_trilep_second_type_branch->SetAddress(&hyp_trilep_second_type_);
	}
	hyp_trilep_third_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_type") != 0) {
		hyp_trilep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_type"));
		hyp_trilep_third_type_branch->SetAddress(&hyp_trilep_third_type_);
	}
	jets_closestElectron_branch = 0;
	if (tree->GetAlias("jets_closestElectron") != 0) {
		jets_closestElectron_branch = tree->GetBranch(tree->GetAlias("jets_closestElectron"));
		jets_closestElectron_branch->SetAddress(&jets_closestElectron_);
	}
	jets_closestMuon_branch = 0;
	if (tree->GetAlias("jets_closestMuon") != 0) {
		jets_closestMuon_branch = tree->GetBranch(tree->GetAlias("jets_closestMuon"));
		jets_closestMuon_branch->SetAddress(&jets_closestMuon_);
	}
	l1_emiso_ieta_branch = 0;
	if (tree->GetAlias("l1_emiso_ieta") != 0) {
		l1_emiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1_emiso_ieta"));
		l1_emiso_ieta_branch->SetAddress(&l1_emiso_ieta_);
	}
	l1_emiso_iphi_branch = 0;
	if (tree->GetAlias("l1_emiso_iphi") != 0) {
		l1_emiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1_emiso_iphi"));
		l1_emiso_iphi_branch->SetAddress(&l1_emiso_iphi_);
	}
	l1_emiso_rawId_branch = 0;
	if (tree->GetAlias("l1_emiso_rawId") != 0) {
		l1_emiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1_emiso_rawId"));
		l1_emiso_rawId_branch->SetAddress(&l1_emiso_rawId_);
	}
	l1_emiso_type_branch = 0;
	if (tree->GetAlias("l1_emiso_type") != 0) {
		l1_emiso_type_branch = tree->GetBranch(tree->GetAlias("l1_emiso_type"));
		l1_emiso_type_branch->SetAddress(&l1_emiso_type_);
	}
	l1_emnoiso_ieta_branch = 0;
	if (tree->GetAlias("l1_emnoiso_ieta") != 0) {
		l1_emnoiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_ieta"));
		l1_emnoiso_ieta_branch->SetAddress(&l1_emnoiso_ieta_);
	}
	l1_emnoiso_iphi_branch = 0;
	if (tree->GetAlias("l1_emnoiso_iphi") != 0) {
		l1_emnoiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_iphi"));
		l1_emnoiso_iphi_branch->SetAddress(&l1_emnoiso_iphi_);
	}
	l1_emnoiso_rawId_branch = 0;
	if (tree->GetAlias("l1_emnoiso_rawId") != 0) {
		l1_emnoiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_rawId"));
		l1_emnoiso_rawId_branch->SetAddress(&l1_emnoiso_rawId_);
	}
	l1_emnoiso_type_branch = 0;
	if (tree->GetAlias("l1_emnoiso_type") != 0) {
		l1_emnoiso_type_branch = tree->GetBranch(tree->GetAlias("l1_emnoiso_type"));
		l1_emnoiso_type_branch->SetAddress(&l1_emnoiso_type_);
	}
	l1_jetsc_ieta_branch = 0;
	if (tree->GetAlias("l1_jetsc_ieta") != 0) {
		l1_jetsc_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_ieta"));
		l1_jetsc_ieta_branch->SetAddress(&l1_jetsc_ieta_);
	}
	l1_jetsc_iphi_branch = 0;
	if (tree->GetAlias("l1_jetsc_iphi") != 0) {
		l1_jetsc_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_iphi"));
		l1_jetsc_iphi_branch->SetAddress(&l1_jetsc_iphi_);
	}
	l1_jetsc_rawId_branch = 0;
	if (tree->GetAlias("l1_jetsc_rawId") != 0) {
		l1_jetsc_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_rawId"));
		l1_jetsc_rawId_branch->SetAddress(&l1_jetsc_rawId_);
	}
	l1_jetsc_type_branch = 0;
	if (tree->GetAlias("l1_jetsc_type") != 0) {
		l1_jetsc_type_branch = tree->GetBranch(tree->GetAlias("l1_jetsc_type"));
		l1_jetsc_type_branch->SetAddress(&l1_jetsc_type_);
	}
	l1_jetsf_ieta_branch = 0;
	if (tree->GetAlias("l1_jetsf_ieta") != 0) {
		l1_jetsf_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_ieta"));
		l1_jetsf_ieta_branch->SetAddress(&l1_jetsf_ieta_);
	}
	l1_jetsf_iphi_branch = 0;
	if (tree->GetAlias("l1_jetsf_iphi") != 0) {
		l1_jetsf_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_iphi"));
		l1_jetsf_iphi_branch->SetAddress(&l1_jetsf_iphi_);
	}
	l1_jetsf_rawId_branch = 0;
	if (tree->GetAlias("l1_jetsf_rawId") != 0) {
		l1_jetsf_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_rawId"));
		l1_jetsf_rawId_branch->SetAddress(&l1_jetsf_rawId_);
	}
	l1_jetsf_type_branch = 0;
	if (tree->GetAlias("l1_jetsf_type") != 0) {
		l1_jetsf_type_branch = tree->GetBranch(tree->GetAlias("l1_jetsf_type"));
		l1_jetsf_type_branch->SetAddress(&l1_jetsf_type_);
	}
	l1_jetst_ieta_branch = 0;
	if (tree->GetAlias("l1_jetst_ieta") != 0) {
		l1_jetst_ieta_branch = tree->GetBranch(tree->GetAlias("l1_jetst_ieta"));
		l1_jetst_ieta_branch->SetAddress(&l1_jetst_ieta_);
	}
	l1_jetst_iphi_branch = 0;
	if (tree->GetAlias("l1_jetst_iphi") != 0) {
		l1_jetst_iphi_branch = tree->GetBranch(tree->GetAlias("l1_jetst_iphi"));
		l1_jetst_iphi_branch->SetAddress(&l1_jetst_iphi_);
	}
	l1_jetst_rawId_branch = 0;
	if (tree->GetAlias("l1_jetst_rawId") != 0) {
		l1_jetst_rawId_branch = tree->GetBranch(tree->GetAlias("l1_jetst_rawId"));
		l1_jetst_rawId_branch->SetAddress(&l1_jetst_rawId_);
	}
	l1_jetst_type_branch = 0;
	if (tree->GetAlias("l1_jetst_type") != 0) {
		l1_jetst_type_branch = tree->GetBranch(tree->GetAlias("l1_jetst_type"));
		l1_jetst_type_branch->SetAddress(&l1_jetst_type_);
	}
	l1_mus_flags_branch = 0;
	if (tree->GetAlias("l1_mus_flags") != 0) {
		l1_mus_flags_branch = tree->GetBranch(tree->GetAlias("l1_mus_flags"));
		l1_mus_flags_branch->SetAddress(&l1_mus_flags_);
	}
	l1_mus_q_branch = 0;
	if (tree->GetAlias("l1_mus_q") != 0) {
		l1_mus_q_branch = tree->GetBranch(tree->GetAlias("l1_mus_q"));
		l1_mus_q_branch->SetAddress(&l1_mus_q_);
	}
	l1_mus_qual_branch = 0;
	if (tree->GetAlias("l1_mus_qual") != 0) {
		l1_mus_qual_branch = tree->GetBranch(tree->GetAlias("l1_mus_qual"));
		l1_mus_qual_branch->SetAddress(&l1_mus_qual_);
	}
	l1_mus_qualFlags_branch = 0;
	if (tree->GetAlias("l1_mus_qualFlags") != 0) {
		l1_mus_qualFlags_branch = tree->GetBranch(tree->GetAlias("l1_mus_qualFlags"));
		l1_mus_qualFlags_branch->SetAddress(&l1_mus_qualFlags_);
	}
	mus_met_flag_branch = 0;
	if (tree->GetAlias("mus_met_flag") != 0) {
		mus_met_flag_branch = tree->GetBranch(tree->GetAlias("mus_met_flag"));
		mus_met_flag_branch->SetAddress(&mus_met_flag_);
	}
	mus_closestEle_branch = 0;
	if (tree->GetAlias("mus_closestEle") != 0) {
		mus_closestEle_branch = tree->GetBranch(tree->GetAlias("mus_closestEle"));
		mus_closestEle_branch->SetAddress(&mus_closestEle_);
	}
	mus_closestJet_branch = 0;
	if (tree->GetAlias("mus_closestJet") != 0) {
		mus_closestJet_branch = tree->GetBranch(tree->GetAlias("mus_closestJet"));
		mus_closestJet_branch->SetAddress(&mus_closestJet_);
	}
	mus_pfmusidx_branch = 0;
	if (tree->GetAlias("mus_pfmusidx") != 0) {
		mus_pfmusidx_branch = tree->GetBranch(tree->GetAlias("mus_pfmusidx"));
		mus_pfmusidx_branch->SetAddress(&mus_pfmusidx_);
	}
	mus_charge_branch = 0;
	if (tree->GetAlias("mus_charge") != 0) {
		mus_charge_branch = tree->GetBranch(tree->GetAlias("mus_charge"));
		mus_charge_branch->SetAddress(&mus_charge_);
	}
	mus_chi2LocalMomentum_branch = 0;
	if (tree->GetAlias("mus_chi2LocalMomentum") != 0) {
		mus_chi2LocalMomentum_branch = tree->GetBranch(tree->GetAlias("mus_chi2LocalMomentum"));
		mus_chi2LocalMomentum_branch->SetAddress(&mus_chi2LocalMomentum_);
	}
	mus_chi2LocalPosition_branch = 0;
	if (tree->GetAlias("mus_chi2LocalPosition") != 0) {
		mus_chi2LocalPosition_branch = tree->GetBranch(tree->GetAlias("mus_chi2LocalPosition"));
		mus_chi2LocalPosition_branch->SetAddress(&mus_chi2LocalPosition_);
	}
	mus_gfit_validHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validHits") != 0) {
		mus_gfit_validHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validHits"));
		mus_gfit_validHits_branch->SetAddress(&mus_gfit_validHits_);
	}
	mus_gfit_validSTAHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validSTAHits") != 0) {
		mus_gfit_validSTAHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validSTAHits"));
		mus_gfit_validSTAHits_branch->SetAddress(&mus_gfit_validSTAHits_);
	}
	mus_gfit_validSiHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validSiHits") != 0) {
		mus_gfit_validSiHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validSiHits"));
		mus_gfit_validSiHits_branch->SetAddress(&mus_gfit_validSiHits_);
	}
	mus_glbKink_branch = 0;
	if (tree->GetAlias("mus_glbKink") != 0) {
		mus_glbKink_branch = tree->GetBranch(tree->GetAlias("mus_glbKink"));
		mus_glbKink_branch->SetAddress(&mus_glbKink_);
	}
	mus_glbTrackProbability_branch = 0;
	if (tree->GetAlias("mus_glbTrackProbability") != 0) {
		mus_glbTrackProbability_branch = tree->GetBranch(tree->GetAlias("mus_glbTrackProbability"));
		mus_glbTrackProbability_branch->SetAddress(&mus_glbTrackProbability_);
	}
	mus_globalDeltaEtaPhi_branch = 0;
	if (tree->GetAlias("mus_globalDeltaEtaPhi") != 0) {
		mus_globalDeltaEtaPhi_branch = tree->GetBranch(tree->GetAlias("mus_globalDeltaEtaPhi"));
		mus_globalDeltaEtaPhi_branch->SetAddress(&mus_globalDeltaEtaPhi_);
	}
	mus_goodmask_branch = 0;
	if (tree->GetAlias("mus_goodmask") != 0) {
		mus_goodmask_branch = tree->GetBranch(tree->GetAlias("mus_goodmask"));
		mus_goodmask_branch->SetAddress(&mus_goodmask_);
	}
	mus_iso03_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso03_ntrk") != 0) {
		mus_iso03_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso03_ntrk"));
		mus_iso03_ntrk_branch->SetAddress(&mus_iso03_ntrk_);
	}
	mus_iso05_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso05_ntrk") != 0) {
		mus_iso05_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso05_ntrk"));
		mus_iso05_ntrk_branch->SetAddress(&mus_iso05_ntrk_);
	}
	mus_localDistance_branch = 0;
	if (tree->GetAlias("mus_localDistance") != 0) {
		mus_localDistance_branch = tree->GetBranch(tree->GetAlias("mus_localDistance"));
		mus_localDistance_branch->SetAddress(&mus_localDistance_);
	}
	mus_lostHits_branch = 0;
	if (tree->GetAlias("mus_lostHits") != 0) {
		mus_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_lostHits"));
		mus_lostHits_branch->SetAddress(&mus_lostHits_);
	}
	mus_nOverlaps_branch = 0;
	if (tree->GetAlias("mus_nOverlaps") != 0) {
		mus_nOverlaps_branch = tree->GetBranch(tree->GetAlias("mus_nOverlaps"));
		mus_nOverlaps_branch->SetAddress(&mus_nOverlaps_);
	}
	mus_nmatches_branch = 0;
	if (tree->GetAlias("mus_nmatches") != 0) {
		mus_nmatches_branch = tree->GetBranch(tree->GetAlias("mus_nmatches"));
		mus_nmatches_branch->SetAddress(&mus_nmatches_);
	}
	mus_overlap0_branch = 0;
	if (tree->GetAlias("mus_overlap0") != 0) {
		mus_overlap0_branch = tree->GetBranch(tree->GetAlias("mus_overlap0"));
		mus_overlap0_branch->SetAddress(&mus_overlap0_);
	}
	mus_overlap1_branch = 0;
	if (tree->GetAlias("mus_overlap1") != 0) {
		mus_overlap1_branch = tree->GetBranch(tree->GetAlias("mus_overlap1"));
		mus_overlap1_branch->SetAddress(&mus_overlap1_);
	}
	mus_pid_TM2DCompatibilityLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityLoose") != 0) {
		mus_pid_TM2DCompatibilityLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityLoose"));
		mus_pid_TM2DCompatibilityLoose_branch->SetAddress(&mus_pid_TM2DCompatibilityLoose_);
	}
	mus_pid_TM2DCompatibilityTight_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityTight") != 0) {
		mus_pid_TM2DCompatibilityTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityTight"));
		mus_pid_TM2DCompatibilityTight_branch->SetAddress(&mus_pid_TM2DCompatibilityTight_);
	}
	mus_pid_TMLastStationLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationLoose") != 0) {
		mus_pid_TMLastStationLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationLoose"));
		mus_pid_TMLastStationLoose_branch->SetAddress(&mus_pid_TMLastStationLoose_);
	}
	mus_pid_TMLastStationTight_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationTight") != 0) {
		mus_pid_TMLastStationTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationTight"));
		mus_pid_TMLastStationTight_branch->SetAddress(&mus_pid_TMLastStationTight_);
	}
	mus_staRelChi2_branch = 0;
	if (tree->GetAlias("mus_staRelChi2") != 0) {
		mus_staRelChi2_branch = tree->GetBranch(tree->GetAlias("mus_staRelChi2"));
		mus_staRelChi2_branch->SetAddress(&mus_staRelChi2_);
	}
	mus_sta_validHits_branch = 0;
	if (tree->GetAlias("mus_sta_validHits") != 0) {
		mus_sta_validHits_branch = tree->GetBranch(tree->GetAlias("mus_sta_validHits"));
		mus_sta_validHits_branch->SetAddress(&mus_sta_validHits_);
	}
	mus_timeDirection_branch = 0;
	if (tree->GetAlias("mus_timeDirection") != 0) {
		mus_timeDirection_branch = tree->GetBranch(tree->GetAlias("mus_timeDirection"));
		mus_timeDirection_branch->SetAddress(&mus_timeDirection_);
	}
	mus_timeNumStationsUsed_branch = 0;
	if (tree->GetAlias("mus_timeNumStationsUsed") != 0) {
		mus_timeNumStationsUsed_branch = tree->GetBranch(tree->GetAlias("mus_timeNumStationsUsed"));
		mus_timeNumStationsUsed_branch->SetAddress(&mus_timeNumStationsUsed_);
	}
	mus_trkKink_branch = 0;
	if (tree->GetAlias("mus_trkKink") != 0) {
		mus_trkKink_branch = tree->GetBranch(tree->GetAlias("mus_trkKink"));
		mus_trkKink_branch->SetAddress(&mus_trkKink_);
	}
	mus_trkRelChi2_branch = 0;
	if (tree->GetAlias("mus_trkRelChi2") != 0) {
		mus_trkRelChi2_branch = tree->GetBranch(tree->GetAlias("mus_trkRelChi2"));
		mus_trkRelChi2_branch->SetAddress(&mus_trkRelChi2_);
	}
	mus_trk_charge_branch = 0;
	if (tree->GetAlias("mus_trk_charge") != 0) {
		mus_trk_charge_branch = tree->GetBranch(tree->GetAlias("mus_trk_charge"));
		mus_trk_charge_branch->SetAddress(&mus_trk_charge_);
	}
	mus_trkidx_branch = 0;
	if (tree->GetAlias("mus_trkidx") != 0) {
		mus_trkidx_branch = tree->GetBranch(tree->GetAlias("mus_trkidx"));
		mus_trkidx_branch->SetAddress(&mus_trkidx_);
	}
	mus_type_branch = 0;
	if (tree->GetAlias("mus_type") != 0) {
		mus_type_branch = tree->GetBranch(tree->GetAlias("mus_type"));
		mus_type_branch->SetAddress(&mus_type_);
	}
	mus_validHits_branch = 0;
	if (tree->GetAlias("mus_validHits") != 0) {
		mus_validHits_branch = tree->GetBranch(tree->GetAlias("mus_validHits"));
		mus_validHits_branch->SetAddress(&mus_validHits_);
	}
	els_pat_genID_branch = 0;
	if (tree->GetAlias("els_pat_genID") != 0) {
		els_pat_genID_branch = tree->GetBranch(tree->GetAlias("els_pat_genID"));
		els_pat_genID_branch->SetAddress(&els_pat_genID_);
	}
	els_pat_genMotherID_branch = 0;
	if (tree->GetAlias("els_pat_genMotherID") != 0) {
		els_pat_genMotherID_branch = tree->GetBranch(tree->GetAlias("els_pat_genMotherID"));
		els_pat_genMotherID_branch->SetAddress(&els_pat_genMotherID_);
	}
	jets_pat_genPartonMother_id_branch = 0;
	if (tree->GetAlias("jets_pat_genPartonMother_id") != 0) {
		jets_pat_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("jets_pat_genPartonMother_id"));
		jets_pat_genPartonMother_id_branch->SetAddress(&jets_pat_genPartonMother_id_);
	}
	jets_pat_genParton_id_branch = 0;
	if (tree->GetAlias("jets_pat_genParton_id") != 0) {
		jets_pat_genParton_id_branch = tree->GetBranch(tree->GetAlias("jets_pat_genParton_id"));
		jets_pat_genParton_id_branch->SetAddress(&jets_pat_genParton_id_);
	}
	jets_pat_jetIDLoose_branch = 0;
	if (tree->GetAlias("jets_pat_jetIDLoose") != 0) {
		jets_pat_jetIDLoose_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetIDLoose"));
		jets_pat_jetIDLoose_branch->SetAddress(&jets_pat_jetIDLoose_);
	}
	jets_pat_jetIDLooseAOD_branch = 0;
	if (tree->GetAlias("jets_pat_jetIDLooseAOD") != 0) {
		jets_pat_jetIDLooseAOD_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetIDLooseAOD"));
		jets_pat_jetIDLooseAOD_branch->SetAddress(&jets_pat_jetIDLooseAOD_);
	}
	jets_pat_jetIDMinimal_branch = 0;
	if (tree->GetAlias("jets_pat_jetIDMinimal") != 0) {
		jets_pat_jetIDMinimal_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetIDMinimal"));
		jets_pat_jetIDMinimal_branch->SetAddress(&jets_pat_jetIDMinimal_);
	}
	jets_pat_jetIDTight_branch = 0;
	if (tree->GetAlias("jets_pat_jetIDTight") != 0) {
		jets_pat_jetIDTight_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetIDTight"));
		jets_pat_jetIDTight_branch->SetAddress(&jets_pat_jetIDTight_);
	}
	jets_pat_partonFlavour_branch = 0;
	if (tree->GetAlias("jets_pat_partonFlavour") != 0) {
		jets_pat_partonFlavour_branch = tree->GetBranch(tree->GetAlias("jets_pat_partonFlavour"));
		jets_pat_partonFlavour_branch->SetAddress(&jets_pat_partonFlavour_);
	}
	mus_pat_genID_branch = 0;
	if (tree->GetAlias("mus_pat_genID") != 0) {
		mus_pat_genID_branch = tree->GetBranch(tree->GetAlias("mus_pat_genID"));
		mus_pat_genID_branch->SetAddress(&mus_pat_genID_);
	}
	mus_pat_genMotherID_branch = 0;
	if (tree->GetAlias("mus_pat_genMotherID") != 0) {
		mus_pat_genMotherID_branch = tree->GetBranch(tree->GetAlias("mus_pat_genMotherID"));
		mus_pat_genMotherID_branch->SetAddress(&mus_pat_genMotherID_);
	}
	pfels_elsidx_branch = 0;
	if (tree->GetAlias("pfels_elsidx") != 0) {
		pfels_elsidx_branch = tree->GetBranch(tree->GetAlias("pfels_elsidx"));
		pfels_elsidx_branch->SetAddress(&pfels_elsidx_);
	}
	pfels_charge_branch = 0;
	if (tree->GetAlias("pfels_charge") != 0) {
		pfels_charge_branch = tree->GetBranch(tree->GetAlias("pfels_charge"));
		pfels_charge_branch->SetAddress(&pfels_charge_);
	}
	pfels_flag_branch = 0;
	if (tree->GetAlias("pfels_flag") != 0) {
		pfels_flag_branch = tree->GetBranch(tree->GetAlias("pfels_flag"));
		pfels_flag_branch->SetAddress(&pfels_flag_);
	}
	pfels_particleId_branch = 0;
	if (tree->GetAlias("pfels_particleId") != 0) {
		pfels_particleId_branch = tree->GetBranch(tree->GetAlias("pfels_particleId"));
		pfels_particleId_branch->SetAddress(&pfels_particleId_);
	}
	pfjets_chargedMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_chargedMultiplicity") != 0) {
		pfjets_chargedMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_chargedMultiplicity"));
		pfjets_chargedMultiplicity_branch->SetAddress(&pfjets_chargedMultiplicity_);
	}
	pfjets_muonMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_muonMultiplicity") != 0) {
		pfjets_muonMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_muonMultiplicity"));
		pfjets_muonMultiplicity_branch->SetAddress(&pfjets_muonMultiplicity_);
	}
	pfjets_neutralMultiplicity_branch = 0;
	if (tree->GetAlias("pfjets_neutralMultiplicity") != 0) {
		pfjets_neutralMultiplicity_branch = tree->GetBranch(tree->GetAlias("pfjets_neutralMultiplicity"));
		pfjets_neutralMultiplicity_branch->SetAddress(&pfjets_neutralMultiplicity_);
	}
	pfmus_musidx_branch = 0;
	if (tree->GetAlias("pfmus_musidx") != 0) {
		pfmus_musidx_branch = tree->GetBranch(tree->GetAlias("pfmus_musidx"));
		pfmus_musidx_branch->SetAddress(&pfmus_musidx_);
	}
	pfmus_charge_branch = 0;
	if (tree->GetAlias("pfmus_charge") != 0) {
		pfmus_charge_branch = tree->GetBranch(tree->GetAlias("pfmus_charge"));
		pfmus_charge_branch->SetAddress(&pfmus_charge_);
	}
	pfmus_flag_branch = 0;
	if (tree->GetAlias("pfmus_flag") != 0) {
		pfmus_flag_branch = tree->GetBranch(tree->GetAlias("pfmus_flag"));
		pfmus_flag_branch->SetAddress(&pfmus_flag_);
	}
	pfmus_particleId_branch = 0;
	if (tree->GetAlias("pfmus_particleId") != 0) {
		pfmus_particleId_branch = tree->GetBranch(tree->GetAlias("pfmus_particleId"));
		pfmus_particleId_branch->SetAddress(&pfmus_particleId_);
	}
	photons_fiduciality_branch = 0;
	if (tree->GetAlias("photons_fiduciality") != 0) {
		photons_fiduciality_branch = tree->GetBranch(tree->GetAlias("photons_fiduciality"));
		photons_fiduciality_branch->SetAddress(&photons_fiduciality_);
	}
	photons_scindex_branch = 0;
	if (tree->GetAlias("photons_scindex") != 0) {
		photons_scindex_branch = tree->GetBranch(tree->GetAlias("photons_scindex"));
		photons_scindex_branch->SetAddress(&photons_scindex_);
	}
	scs_detIdSeed_branch = 0;
	if (tree->GetAlias("scs_detIdSeed") != 0) {
		scs_detIdSeed_branch = tree->GetBranch(tree->GetAlias("scs_detIdSeed"));
		scs_detIdSeed_branch->SetAddress(&scs_detIdSeed_);
	}
	scs_elsidx_branch = 0;
	if (tree->GetAlias("scs_elsidx") != 0) {
		scs_elsidx_branch = tree->GetBranch(tree->GetAlias("scs_elsidx"));
		scs_elsidx_branch->SetAddress(&scs_elsidx_);
	}
	scs_severitySeed_branch = 0;
	if (tree->GetAlias("scs_severitySeed") != 0) {
		scs_severitySeed_branch = tree->GetBranch(tree->GetAlias("scs_severitySeed"));
		scs_severitySeed_branch->SetAddress(&scs_severitySeed_);
	}
	svs_isKs_branch = 0;
	if (tree->GetAlias("svs_isKs") != 0) {
		svs_isKs_branch = tree->GetBranch(tree->GetAlias("svs_isKs"));
		svs_isKs_branch->SetAddress(&svs_isKs_);
	}
	svs_isLambda_branch = 0;
	if (tree->GetAlias("svs_isLambda") != 0) {
		svs_isLambda_branch = tree->GetBranch(tree->GetAlias("svs_isLambda"));
		svs_isLambda_branch->SetAddress(&svs_isLambda_);
	}
	svs_mc3_id_branch = 0;
	if (tree->GetAlias("svs_mc3_id") != 0) {
		svs_mc3_id_branch = tree->GetBranch(tree->GetAlias("svs_mc3_id"));
		svs_mc3_id_branch->SetAddress(&svs_mc3_id_);
	}
	svs_nTrks_branch = 0;
	if (tree->GetAlias("svs_nTrks") != 0) {
		svs_nTrks_branch = tree->GetBranch(tree->GetAlias("svs_nTrks"));
		svs_nTrks_branch->SetAddress(&svs_nTrks_);
	}
	mus_tcmet_flag_branch = 0;
	if (tree->GetAlias("mus_tcmet_flag") != 0) {
		mus_tcmet_flag_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_flag"));
		mus_tcmet_flag_branch->SetAddress(&mus_tcmet_flag_);
	}
	trks_algo_branch = 0;
	if (tree->GetAlias("trks_algo") != 0) {
		trks_algo_branch = tree->GetBranch(tree->GetAlias("trks_algo"));
		trks_algo_branch->SetAddress(&trks_algo_);
	}
	trks_charge_branch = 0;
	if (tree->GetAlias("trks_charge") != 0) {
		trks_charge_branch = tree->GetBranch(tree->GetAlias("trks_charge"));
		trks_charge_branch->SetAddress(&trks_charge_);
	}
	trks_exp_innerlayers_branch = 0;
	if (tree->GetAlias("trks_exp_innerlayers") != 0) {
		trks_exp_innerlayers_branch = tree->GetBranch(tree->GetAlias("trks_exp_innerlayers"));
		trks_exp_innerlayers_branch->SetAddress(&trks_exp_innerlayers_);
	}
	trks_exp_outerlayers_branch = 0;
	if (tree->GetAlias("trks_exp_outerlayers") != 0) {
		trks_exp_outerlayers_branch = tree->GetBranch(tree->GetAlias("trks_exp_outerlayers"));
		trks_exp_outerlayers_branch->SetAddress(&trks_exp_outerlayers_);
	}
	trks_layer1_det_branch = 0;
	if (tree->GetAlias("trks_layer1_det") != 0) {
		trks_layer1_det_branch = tree->GetBranch(tree->GetAlias("trks_layer1_det"));
		trks_layer1_det_branch->SetAddress(&trks_layer1_det_);
	}
	trks_layer1_layer_branch = 0;
	if (tree->GetAlias("trks_layer1_layer") != 0) {
		trks_layer1_layer_branch = tree->GetBranch(tree->GetAlias("trks_layer1_layer"));
		trks_layer1_layer_branch->SetAddress(&trks_layer1_layer_);
	}
	trks_layer1_sizerphi_branch = 0;
	if (tree->GetAlias("trks_layer1_sizerphi") != 0) {
		trks_layer1_sizerphi_branch = tree->GetBranch(tree->GetAlias("trks_layer1_sizerphi"));
		trks_layer1_sizerphi_branch->SetAddress(&trks_layer1_sizerphi_);
	}
	trks_layer1_sizerz_branch = 0;
	if (tree->GetAlias("trks_layer1_sizerz") != 0) {
		trks_layer1_sizerz_branch = tree->GetBranch(tree->GetAlias("trks_layer1_sizerz"));
		trks_layer1_sizerz_branch->SetAddress(&trks_layer1_sizerz_);
	}
	trks_lostHits_branch = 0;
	if (tree->GetAlias("trks_lostHits") != 0) {
		trks_lostHits_branch = tree->GetBranch(tree->GetAlias("trks_lostHits"));
		trks_lostHits_branch->SetAddress(&trks_lostHits_);
	}
	trks_lost_pixelhits_branch = 0;
	if (tree->GetAlias("trks_lost_pixelhits") != 0) {
		trks_lost_pixelhits_branch = tree->GetBranch(tree->GetAlias("trks_lost_pixelhits"));
		trks_lost_pixelhits_branch->SetAddress(&trks_lost_pixelhits_);
	}
	trks_nlayers_branch = 0;
	if (tree->GetAlias("trks_nlayers") != 0) {
		trks_nlayers_branch = tree->GetBranch(tree->GetAlias("trks_nlayers"));
		trks_nlayers_branch->SetAddress(&trks_nlayers_);
	}
	trks_nlayers3D_branch = 0;
	if (tree->GetAlias("trks_nlayers3D") != 0) {
		trks_nlayers3D_branch = tree->GetBranch(tree->GetAlias("trks_nlayers3D"));
		trks_nlayers3D_branch->SetAddress(&trks_nlayers3D_);
	}
	trks_nlayersLost_branch = 0;
	if (tree->GetAlias("trks_nlayersLost") != 0) {
		trks_nlayersLost_branch = tree->GetBranch(tree->GetAlias("trks_nlayersLost"));
		trks_nlayersLost_branch->SetAddress(&trks_nlayersLost_);
	}
	trks_qualityMask_branch = 0;
	if (tree->GetAlias("trks_qualityMask") != 0) {
		trks_qualityMask_branch = tree->GetBranch(tree->GetAlias("trks_qualityMask"));
		trks_qualityMask_branch->SetAddress(&trks_qualityMask_);
	}
	trks_validHits_branch = 0;
	if (tree->GetAlias("trks_validHits") != 0) {
		trks_validHits_branch = tree->GetBranch(tree->GetAlias("trks_validHits"));
		trks_validHits_branch->SetAddress(&trks_validHits_);
	}
	trks_valid_pixelhits_branch = 0;
	if (tree->GetAlias("trks_valid_pixelhits") != 0) {
		trks_valid_pixelhits_branch = tree->GetBranch(tree->GetAlias("trks_valid_pixelhits"));
		trks_valid_pixelhits_branch->SetAddress(&trks_valid_pixelhits_);
	}
	trks_elsidx_branch = 0;
	if (tree->GetAlias("trks_elsidx") != 0) {
		trks_elsidx_branch = tree->GetBranch(tree->GetAlias("trks_elsidx"));
		trks_elsidx_branch->SetAddress(&trks_elsidx_);
	}
	trk_musidx_branch = 0;
	if (tree->GetAlias("trk_musidx") != 0) {
		trk_musidx_branch = tree->GetBranch(tree->GetAlias("trk_musidx"));
		trk_musidx_branch->SetAddress(&trk_musidx_);
	}
	trkjets_ntrks_branch = 0;
	if (tree->GetAlias("trkjets_ntrks") != 0) {
		trkjets_ntrks_branch = tree->GetBranch(tree->GetAlias("trkjets_ntrks"));
		trkjets_ntrks_branch->SetAddress(&trkjets_ntrks_);
	}
	trkjets_vtxs_idx_branch = 0;
	if (tree->GetAlias("trkjets_vtxs_idx") != 0) {
		trkjets_vtxs_idx_branch = tree->GetBranch(tree->GetAlias("trkjets_vtxs_idx"));
		trkjets_vtxs_idx_branch->SetAddress(&trkjets_vtxs_idx_);
	}
	vtxs_isFake_branch = 0;
	if (tree->GetAlias("vtxs_isFake") != 0) {
		vtxs_isFake_branch = tree->GetBranch(tree->GetAlias("vtxs_isFake"));
		vtxs_isFake_branch->SetAddress(&vtxs_isFake_);
	}
	vtxs_isValid_branch = 0;
	if (tree->GetAlias("vtxs_isValid") != 0) {
		vtxs_isValid_branch = tree->GetBranch(tree->GetAlias("vtxs_isValid"));
		vtxs_isValid_branch->SetAddress(&vtxs_isValid_);
	}
	vtxs_tracksSize_branch = 0;
	if (tree->GetAlias("vtxs_tracksSize") != 0) {
		vtxs_tracksSize_branch = tree->GetBranch(tree->GetAlias("vtxs_tracksSize"));
		vtxs_tracksSize_branch->SetAddress(&vtxs_tracksSize_);
	}
	genps_lepdaughter_id_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_id") != 0) {
		genps_lepdaughter_id_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_id"));
		genps_lepdaughter_id_branch->SetAddress(&genps_lepdaughter_id_);
	}
	genps_lepdaughter_idx_branch = 0;
	if (tree->GetAlias("genps_lepdaughter_idx") != 0) {
		genps_lepdaughter_idx_branch = tree->GetBranch(tree->GetAlias("genps_lepdaughter_idx"));
		genps_lepdaughter_idx_branch->SetAddress(&genps_lepdaughter_idx_);
	}
	hlt_trigObjs_id_branch = 0;
	if (tree->GetAlias("hlt_trigObjs_id") != 0) {
		hlt_trigObjs_id_branch = tree->GetBranch(tree->GetAlias("hlt_trigObjs_id"));
		hlt_trigObjs_id_branch->SetAddress(&hlt_trigObjs_id_);
	}
	hyp_jets_idx_branch = 0;
	if (tree->GetAlias("hyp_jets_idx") != 0) {
		hyp_jets_idx_branch = tree->GetBranch(tree->GetAlias("hyp_jets_idx"));
		hyp_jets_idx_branch->SetAddress(&hyp_jets_idx_);
	}
	hyp_other_jets_idx_branch = 0;
	if (tree->GetAlias("hyp_other_jets_idx") != 0) {
		hyp_other_jets_idx_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_idx"));
		hyp_other_jets_idx_branch->SetAddress(&hyp_other_jets_idx_);
	}
	evt_nels_branch = 0;
	if (tree->GetAlias("evt_nels") != 0) {
		evt_nels_branch = tree->GetBranch(tree->GetAlias("evt_nels"));
		evt_nels_branch->SetAddress(&evt_nels_);
	}
	evt_detectorStatus_branch = 0;
	if (tree->GetAlias("evt_detectorStatus") != 0) {
		evt_detectorStatus_branch = tree->GetBranch(tree->GetAlias("evt_detectorStatus"));
		evt_detectorStatus_branch->SetAddress(&evt_detectorStatus_);
	}
	evt_event_branch = 0;
	if (tree->GetAlias("evt_event") != 0) {
		evt_event_branch = tree->GetBranch(tree->GetAlias("evt_event"));
		evt_event_branch->SetAddress(&evt_event_);
	}
	evt_lumiBlock_branch = 0;
	if (tree->GetAlias("evt_lumiBlock") != 0) {
		evt_lumiBlock_branch = tree->GetBranch(tree->GetAlias("evt_lumiBlock"));
		evt_lumiBlock_branch->SetAddress(&evt_lumiBlock_);
	}
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		evt_run_branch->SetAddress(&evt_run_);
	}
	genps_flavorHistoryFilterResult_branch = 0;
	if (tree->GetAlias("genps_flavorHistoryFilterResult") != 0) {
		genps_flavorHistoryFilterResult_branch = tree->GetBranch(tree->GetAlias("genps_flavorHistoryFilterResult"));
		genps_flavorHistoryFilterResult_branch->SetAddress(&genps_flavorHistoryFilterResult_);
	}
	evt_ngenjets_branch = 0;
	if (tree->GetAlias("evt_ngenjets") != 0) {
		evt_ngenjets_branch = tree->GetBranch(tree->GetAlias("evt_ngenjets"));
		evt_ngenjets_branch->SetAddress(&evt_ngenjets_);
	}
	genps_signalProcessID_branch = 0;
	if (tree->GetAlias("genps_signalProcessID") != 0) {
		genps_signalProcessID_branch = tree->GetBranch(tree->GetAlias("genps_signalProcessID"));
		genps_signalProcessID_branch->SetAddress(&genps_signalProcessID_);
	}
	hlt_bits1_branch = 0;
	if (tree->GetAlias("hlt_bits1") != 0) {
		hlt_bits1_branch = tree->GetBranch(tree->GetAlias("hlt_bits1"));
		hlt_bits1_branch->SetAddress(&hlt_bits1_);
	}
	hlt_bits2_branch = 0;
	if (tree->GetAlias("hlt_bits2") != 0) {
		hlt_bits2_branch = tree->GetBranch(tree->GetAlias("hlt_bits2"));
		hlt_bits2_branch->SetAddress(&hlt_bits2_);
	}
	hlt_bits3_branch = 0;
	if (tree->GetAlias("hlt_bits3") != 0) {
		hlt_bits3_branch = tree->GetBranch(tree->GetAlias("hlt_bits3"));
		hlt_bits3_branch->SetAddress(&hlt_bits3_);
	}
	hlt_bits4_branch = 0;
	if (tree->GetAlias("hlt_bits4") != 0) {
		hlt_bits4_branch = tree->GetBranch(tree->GetAlias("hlt_bits4"));
		hlt_bits4_branch->SetAddress(&hlt_bits4_);
	}
	hlt_bits5_branch = 0;
	if (tree->GetAlias("hlt_bits5") != 0) {
		hlt_bits5_branch = tree->GetBranch(tree->GetAlias("hlt_bits5"));
		hlt_bits5_branch->SetAddress(&hlt_bits5_);
	}
	hlt_bits6_branch = 0;
	if (tree->GetAlias("hlt_bits6") != 0) {
		hlt_bits6_branch = tree->GetBranch(tree->GetAlias("hlt_bits6"));
		hlt_bits6_branch->SetAddress(&hlt_bits6_);
	}
	hlt_bits7_branch = 0;
	if (tree->GetAlias("hlt_bits7") != 0) {
		hlt_bits7_branch = tree->GetBranch(tree->GetAlias("hlt_bits7"));
		hlt_bits7_branch->SetAddress(&hlt_bits7_);
	}
	hlt_bits8_branch = 0;
	if (tree->GetAlias("hlt_bits8") != 0) {
		hlt_bits8_branch = tree->GetBranch(tree->GetAlias("hlt_bits8"));
		hlt_bits8_branch->SetAddress(&hlt_bits8_);
	}
	evt_njets_branch = 0;
	if (tree->GetAlias("evt_njets") != 0) {
		evt_njets_branch = tree->GetBranch(tree->GetAlias("evt_njets"));
		evt_njets_branch->SetAddress(&evt_njets_);
	}
	evt_njpts_branch = 0;
	if (tree->GetAlias("evt_njpts") != 0) {
		evt_njpts_branch = tree->GetBranch(tree->GetAlias("evt_njpts"));
		evt_njpts_branch->SetAddress(&evt_njpts_);
	}
	l1_bits1_branch = 0;
	if (tree->GetAlias("l1_bits1") != 0) {
		l1_bits1_branch = tree->GetBranch(tree->GetAlias("l1_bits1"));
		l1_bits1_branch->SetAddress(&l1_bits1_);
	}
	l1_bits2_branch = 0;
	if (tree->GetAlias("l1_bits2") != 0) {
		l1_bits2_branch = tree->GetBranch(tree->GetAlias("l1_bits2"));
		l1_bits2_branch->SetAddress(&l1_bits2_);
	}
	l1_bits3_branch = 0;
	if (tree->GetAlias("l1_bits3") != 0) {
		l1_bits3_branch = tree->GetBranch(tree->GetAlias("l1_bits3"));
		l1_bits3_branch->SetAddress(&l1_bits3_);
	}
	l1_bits4_branch = 0;
	if (tree->GetAlias("l1_bits4") != 0) {
		l1_bits4_branch = tree->GetBranch(tree->GetAlias("l1_bits4"));
		l1_bits4_branch->SetAddress(&l1_bits4_);
	}
	l1_techbits1_branch = 0;
	if (tree->GetAlias("l1_techbits1") != 0) {
		l1_techbits1_branch = tree->GetBranch(tree->GetAlias("l1_techbits1"));
		l1_techbits1_branch->SetAddress(&l1_techbits1_);
	}
	l1_techbits2_branch = 0;
	if (tree->GetAlias("l1_techbits2") != 0) {
		l1_techbits2_branch = tree->GetBranch(tree->GetAlias("l1_techbits2"));
		l1_techbits2_branch->SetAddress(&l1_techbits2_);
	}
	evt_nphotons_branch = 0;
	if (tree->GetAlias("evt_nphotons") != 0) {
		evt_nphotons_branch = tree->GetBranch(tree->GetAlias("evt_nphotons"));
		evt_nphotons_branch->SetAddress(&evt_nphotons_);
	}
	evt_ecalRecoStatus_branch = 0;
	if (tree->GetAlias("evt_ecalRecoStatus") != 0) {
		evt_ecalRecoStatus_branch = tree->GetBranch(tree->GetAlias("evt_ecalRecoStatus"));
		evt_ecalRecoStatus_branch->SetAddress(&evt_ecalRecoStatus_);
	}
	evt_nscs_branch = 0;
	if (tree->GetAlias("evt_nscs") != 0) {
		evt_nscs_branch = tree->GetBranch(tree->GetAlias("evt_nscs"));
		evt_nscs_branch->SetAddress(&evt_nscs_);
	}
	evt_ntrkjets_branch = 0;
	if (tree->GetAlias("evt_ntrkjets") != 0) {
		evt_ntrkjets_branch = tree->GetBranch(tree->GetAlias("evt_ntrkjets"));
		evt_ntrkjets_branch->SetAddress(&evt_ntrkjets_);
	}
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	hlt_prescales_branch = 0;
	if (tree->GetAlias("hlt_prescales") != 0) {
		hlt_prescales_branch = tree->GetBranch(tree->GetAlias("hlt_prescales"));
		hlt_prescales_branch->SetAddress(&hlt_prescales_);
	}
	hyp_quadlep_bucket_branch = 0;
	if (tree->GetAlias("hyp_quadlep_bucket") != 0) {
		hyp_quadlep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_bucket"));
		hyp_quadlep_bucket_branch->SetAddress(&hyp_quadlep_bucket_);
	}
	hyp_quadlep_first_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_index") != 0) {
		hyp_quadlep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_index"));
		hyp_quadlep_first_index_branch->SetAddress(&hyp_quadlep_first_index_);
	}
	hyp_quadlep_fourth_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_index") != 0) {
		hyp_quadlep_fourth_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_index"));
		hyp_quadlep_fourth_index_branch->SetAddress(&hyp_quadlep_fourth_index_);
	}
	hyp_quadlep_second_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_index") != 0) {
		hyp_quadlep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_index"));
		hyp_quadlep_second_index_branch->SetAddress(&hyp_quadlep_second_index_);
	}
	hyp_quadlep_third_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_index") != 0) {
		hyp_quadlep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_index"));
		hyp_quadlep_third_index_branch->SetAddress(&hyp_quadlep_third_index_);
	}
	hyp_trilep_bucket_branch = 0;
	if (tree->GetAlias("hyp_trilep_bucket") != 0) {
		hyp_trilep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_bucket"));
		hyp_trilep_bucket_branch->SetAddress(&hyp_trilep_bucket_);
	}
	hyp_trilep_first_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_index") != 0) {
		hyp_trilep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_index"));
		hyp_trilep_first_index_branch->SetAddress(&hyp_trilep_first_index_);
	}
	hyp_trilep_second_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_index") != 0) {
		hyp_trilep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_index"));
		hyp_trilep_second_index_branch->SetAddress(&hyp_trilep_second_index_);
	}
	hyp_trilep_third_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_index") != 0) {
		hyp_trilep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_index"));
		hyp_trilep_third_index_branch->SetAddress(&hyp_trilep_third_index_);
	}
	l1_prescales_branch = 0;
	if (tree->GetAlias("l1_prescales") != 0) {
		l1_prescales_branch = tree->GetBranch(tree->GetAlias("l1_prescales"));
		l1_prescales_branch->SetAddress(&l1_prescales_);
	}
	l1_techtrigprescales_branch = 0;
	if (tree->GetAlias("l1_techtrigprescales") != 0) {
		l1_techtrigprescales_branch = tree->GetBranch(tree->GetAlias("l1_techtrigprescales"));
		l1_techtrigprescales_branch->SetAddress(&l1_techtrigprescales_);
	}
	els_pat_flag_branch = 0;
	if (tree->GetAlias("els_pat_flag") != 0) {
		els_pat_flag_branch = tree->GetBranch(tree->GetAlias("els_pat_flag"));
		els_pat_flag_branch->SetAddress(&els_pat_flag_);
	}
	jets_pat_flag_branch = 0;
	if (tree->GetAlias("jets_pat_flag") != 0) {
		jets_pat_flag_branch = tree->GetBranch(tree->GetAlias("jets_pat_flag"));
		jets_pat_flag_branch->SetAddress(&jets_pat_flag_);
	}
	mus_pat_flag_branch = 0;
	if (tree->GetAlias("mus_pat_flag") != 0) {
		mus_pat_flag_branch = tree->GetBranch(tree->GetAlias("mus_pat_flag"));
		mus_pat_flag_branch->SetAddress(&mus_pat_flag_);
	}
	evt_nEvts_branch = 0;
	if (tree->GetAlias("evt_nEvts") != 0) {
		evt_nEvts_branch = tree->GetBranch(tree->GetAlias("evt_nEvts"));
		evt_nEvts_branch->SetAddress(&evt_nEvts_);
	}
	evt_filt_eff_branch = 0;
	if (tree->GetAlias("evt_filt_eff") != 0) {
		evt_filt_eff_branch = tree->GetBranch(tree->GetAlias("evt_filt_eff"));
		evt_filt_eff_branch->SetAddress(&evt_filt_eff_);
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		evt_CMS2tag_isLoaded = false;
		evt_dataset_isLoaded = false;
		hlt_trigNames_isLoaded = false;
		l1_techtrigNames_isLoaded = false;
		l1_trigNames_isLoaded = false;
		evt_errCategory_isLoaded = false;
		evt_errModule_isLoaded = false;
		evt_errSeverity_isLoaded = false;
		evt_eventHasHalo_isLoaded = false;
		evt_hbheFilter_isLoaded = false;
		mus_tightMatch_isLoaded = false;
		mus_updatedSta_isLoaded = false;
		photons_haspixelSeed_isLoaded = false;
		jets_closestElectron_DR_isLoaded = false;
		jets_closestMuon_DR_isLoaded = false;
		evt_bs_Xwidth_isLoaded = false;
		evt_bs_XwidthErr_isLoaded = false;
		evt_bs_Ywidth_isLoaded = false;
		evt_bs_YwidthErr_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_xErr_isLoaded = false;
		evt_bs_yErr_isLoaded = false;
		evt_bs_zErr_isLoaded = false;
		evthcal_dmetx_isLoaded = false;
		evthcal_dmety_isLoaded = false;
		evthcal_dsumet_isLoaded = false;
		evthf_dmetx_isLoaded = false;
		evthf_dmety_isLoaded = false;
		evthf_dsumet_isLoaded = false;
		evt_bField_isLoaded = false;
		evt_kfactor_isLoaded = false;
		evt_scale1fb_isLoaded = false;
		evt_xsec_excl_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metPhi_isLoaded = false;
		genps_alphaQCD_isLoaded = false;
		genps_pthat_isLoaded = false;
		genps_qScale_isLoaded = false;
		genps_weight_isLoaded = false;
		gen_sumEt_isLoaded = false;
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
		evt_ecalendcapm_met_isLoaded = false;
		evt_ecalendcapm_metPhi_isLoaded = false;
		evt_ecalendcapp_met_isLoaded = false;
		evt_ecalendcapp_metPhi_isLoaded = false;
		evt_ecalmet_isLoaded = false;
		evt_ecalmetPhi_isLoaded = false;
		evt_endcapm_met_isLoaded = false;
		evt_endcapm_metPhi_isLoaded = false;
		evt_endcapp_met_isLoaded = false;
		evt_endcapp_metPhi_isLoaded = false;
		evt_hcalendcapm_met_isLoaded = false;
		evt_hcalendcapm_metPhi_isLoaded = false;
		evt_hcalendcapp_met_isLoaded = false;
		evt_hcalendcapp_metPhi_isLoaded = false;
		evt_hcalmet_isLoaded = false;
		evt_hcalmetPhi_isLoaded = false;
		evt_met_isLoaded = false;
		evt_metHO_isLoaded = false;
		evt_metHOPhi_isLoaded = false;
		evt_metHOSig_isLoaded = false;
		evt_metMuonCorr_isLoaded = false;
		evt_metMuonCorrPhi_isLoaded = false;
		evt_metMuonCorrSig_isLoaded = false;
		evt_metMuonJESCorr_isLoaded = false;
		evt_metMuonJESCorrPhi_isLoaded = false;
		evt_metMuonJESCorrSig_isLoaded = false;
		evt_metNoHF_isLoaded = false;
		evt_metNoHFHO_isLoaded = false;
		evt_metNoHFHOPhi_isLoaded = false;
		evt_metNoHFHOSig_isLoaded = false;
		evt_metNoHFPhi_isLoaded = false;
		evt_metNoHFSig_isLoaded = false;
		evt_metOpt_isLoaded = false;
		evt_metOptHO_isLoaded = false;
		evt_metOptHOPhi_isLoaded = false;
		evt_metOptHOSig_isLoaded = false;
		evt_metOptNoHF_isLoaded = false;
		evt_metOptNoHFHO_isLoaded = false;
		evt_metOptNoHFHOPhi_isLoaded = false;
		evt_metOptNoHFHOSig_isLoaded = false;
		evt_metOptNoHFPhi_isLoaded = false;
		evt_metOptNoHFSig_isLoaded = false;
		evt_metOptPhi_isLoaded = false;
		evt_metOptSig_isLoaded = false;
		evt_metPhi_isLoaded = false;
		evt_metSig_isLoaded = false;
		evt_sumet_isLoaded = false;
		evt_sumetHO_isLoaded = false;
		evt_sumetMuonCorr_isLoaded = false;
		evt_sumetNoHF_isLoaded = false;
		evt_sumetNoHFHO_isLoaded = false;
		evt_sumetOpt_isLoaded = false;
		evt_sumetOptHO_isLoaded = false;
		evt_sumetOptNoHF_isLoaded = false;
		evt_sumetOptNoHFHO_isLoaded = false;
		met_pat_metCor_isLoaded = false;
		met_pat_metPhiCor_isLoaded = false;
		met_pat_metPhiUncor_isLoaded = false;
		met_pat_metPhiUncorJES_isLoaded = false;
		met_pat_metPhiUncorMuon_isLoaded = false;
		met_pat_metUncor_isLoaded = false;
		met_pat_metUncorJES_isLoaded = false;
		met_pat_metUncorMuon_isLoaded = false;
		pdfinfo_scale_isLoaded = false;
		pdfinfo_x1_isLoaded = false;
		pdfinfo_x2_isLoaded = false;
		evt_pfmet_isLoaded = false;
		evt_pfmetPhi_isLoaded = false;
		evt_pfmetSig_isLoaded = false;
		evt_pfsumet_isLoaded = false;
		evt_tcmet_isLoaded = false;
		evt_tcmetPhi_isLoaded = false;
		evt_tcmetSig_isLoaded = false;
		evt_tcsumet_isLoaded = false;
		evt_bsp4_isLoaded = false;
		l1_met_p4_isLoaded = false;
		l1_mht_p4_isLoaded = false;
		els_mc_motherp4_isLoaded = false;
		els_mc_p4_isLoaded = false;
		jets_mc_gp_p4_isLoaded = false;
		jets_mc_motherp4_isLoaded = false;
		jets_mc_p4_isLoaded = false;
		mus_mc_motherp4_isLoaded = false;
		mus_mc_p4_isLoaded = false;
		pfjets_mc_gp_p4_isLoaded = false;
		pfjets_mc_motherp4_isLoaded = false;
		pfjets_mc_p4_isLoaded = false;
		photons_mc_motherp4_isLoaded = false;
		photons_mc_p4_isLoaded = false;
		trk_mcp4_isLoaded = false;
		els_conv_pos_p4_isLoaded = false;
		els_inner_position_isLoaded = false;
		els_outer_position_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
		els_vertex_p4_isLoaded = false;
		genjets_p4_isLoaded = false;
		genps_p4_isLoaded = false;
		genps_prod_vtx_isLoaded = false;
		gsftrks_inner_position_isLoaded = false;
		gsftrks_outer_p4_isLoaded = false;
		gsftrks_outer_position_isLoaded = false;
		gsftrks_p4_isLoaded = false;
		gsftrks_vertex_p4_isLoaded = false;
		hyp_ll_p4_isLoaded = false;
		hyp_ll_trk_p4_isLoaded = false;
		hyp_lt_p4_isLoaded = false;
		hyp_lt_trk_p4_isLoaded = false;
		hyp_p4_isLoaded = false;
		hyp_FVFit_p4_isLoaded = false;
		hyp_FVFit_v4_isLoaded = false;
		hyp_ll_mc_p4_isLoaded = false;
		hyp_lt_mc_p4_isLoaded = false;
		jets_p4_isLoaded = false;
		jets_vertex_p4_isLoaded = false;
		jpts_p4_isLoaded = false;
		l1_emiso_p4_isLoaded = false;
		l1_emnoiso_p4_isLoaded = false;
		l1_jetsc_p4_isLoaded = false;
		l1_jetsf_p4_isLoaded = false;
		l1_jetst_p4_isLoaded = false;
		l1_mus_p4_isLoaded = false;
		mus_ecalpos_p4_isLoaded = false;
		mus_fitdefault_p4_isLoaded = false;
		mus_fitfirsthit_p4_isLoaded = false;
		mus_fitpicky_p4_isLoaded = false;
		mus_fittev_p4_isLoaded = false;
		mus_gfit_outerPos_p4_isLoaded = false;
		mus_gfit_p4_isLoaded = false;
		mus_gfit_vertex_p4_isLoaded = false;
		mus_p4_isLoaded = false;
		mus_sta_p4_isLoaded = false;
		mus_sta_vertex_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		mus_vertex_p4_isLoaded = false;
		els_pat_genMotherP4_isLoaded = false;
		els_pat_genP4_isLoaded = false;
		els_pat_p4_isLoaded = false;
		jets_pat_genJet_p4_isLoaded = false;
		jets_pat_genPartonMother_p4_isLoaded = false;
		jets_pat_genParton_p4_isLoaded = false;
		jets_pat_jet_p4_isLoaded = false;
		jets_pat_jet_uncorp4_isLoaded = false;
		mus_pat_genMotherP4_isLoaded = false;
		mus_pat_genP4_isLoaded = false;
		mus_pat_p4_isLoaded = false;
		pfels_p4_isLoaded = false;
		pfels_posAtEcal_p4_isLoaded = false;
		pfjets_p4_isLoaded = false;
		pfmus_p4_isLoaded = false;
		pfmus_posAtEcal_p4_isLoaded = false;
		photons_p4_isLoaded = false;
		scs_p4_isLoaded = false;
		scs_pos_p4_isLoaded = false;
		scs_vtx_p4_isLoaded = false;
		svs_flight_isLoaded = false;
		svs_mc3_p4_isLoaded = false;
		svs_p4_isLoaded = false;
		svs_position_isLoaded = false;
		svs_refitp4_isLoaded = false;
		trks_inner_position_isLoaded = false;
		trks_outer_p4_isLoaded = false;
		trks_outer_position_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		trks_vertex_p4_isLoaded = false;
		trkjets_p4_isLoaded = false;
		vtxs_position_isLoaded = false;
		genps_lepdaughter_p4_isLoaded = false;
		hlt_trigObjs_p4_isLoaded = false;
		hyp_jets_p4_isLoaded = false;
		hyp_other_jets_p4_isLoaded = false;
		jpts_combinedSecondaryVertexBJetTag_isLoaded = false;
		jpts_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		jpts_jetBProbabilityBJetTag_isLoaded = false;
		jpts_jetProbabilityBJetTag_isLoaded = false;
		jpts_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		jpts_simpleSecondaryVertexHighPurBJetTags_isLoaded = false;
		jpts_softElectronByIP3dBJetTag_isLoaded = false;
		jpts_softElectronByPtBJetTag_isLoaded = false;
		jpts_softMuonBJetTag_isLoaded = false;
		jpts_softMuonByIP3dBJetTag_isLoaded = false;
		jpts_softMuonByPtBJetTag_isLoaded = false;
		jpts_trackCountingHighEffBJetTag_isLoaded = false;
		jpts_trackCountingHighPurBJetTag_isLoaded = false;
		jets_combinedSecondaryVertexBJetTag_isLoaded = false;
		jets_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		jets_jetBProbabilityBJetTag_isLoaded = false;
		jets_jetProbabilityBJetTag_isLoaded = false;
		jets_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		jets_simpleSecondaryVertexHighPurBJetTags_isLoaded = false;
		jets_softElectronByIP3dBJetTag_isLoaded = false;
		jets_softElectronByPtBJetTag_isLoaded = false;
		jets_softMuonBJetTag_isLoaded = false;
		jets_softMuonByIP3dBJetTag_isLoaded = false;
		jets_softMuonByPtBJetTag_isLoaded = false;
		jets_trackCountingHighEffBJetTag_isLoaded = false;
		jets_trackCountingHighPurBJetTag_isLoaded = false;
		pfjets_combinedSecondaryVertexBJetTag_isLoaded = false;
		pfjets_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		pfjets_jetBProbabilityBJetTag_isLoaded = false;
		pfjets_jetProbabilityBJetTag_isLoaded = false;
		pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = false;
		pfjets_softElectronByIP3dBJetTag_isLoaded = false;
		pfjets_softElectronByPtBJetTag_isLoaded = false;
		pfjets_softMuonBJetTag_isLoaded = false;
		pfjets_softMuonByIP3dBJetTag_isLoaded = false;
		pfjets_softMuonByPtBJetTag_isLoaded = false;
		pfjets_trackCountingHighEffBJetTag_isLoaded = false;
		pfjets_trackCountingHighPurBJetTag_isLoaded = false;
		trkjets_combinedSecondaryVertexBJetTag_isLoaded = false;
		trkjets_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		trkjets_jetBProbabilityBJetTag_isLoaded = false;
		trkjets_jetProbabilityBJetTag_isLoaded = false;
		trkjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		trkjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = false;
		trkjets_softElectronByIP3dBJetTag_isLoaded = false;
		trkjets_softElectronByPtBJetTag_isLoaded = false;
		trkjets_softMuonBJetTag_isLoaded = false;
		trkjets_softMuonByIP3dBJetTag_isLoaded = false;
		trkjets_softMuonByPtBJetTag_isLoaded = false;
		trkjets_trackCountingHighEffBJetTag_isLoaded = false;
		trkjets_trackCountingHighPurBJetTag_isLoaded = false;
		evt_bs_covMatrix_isLoaded = false;
		els_mc3dr_isLoaded = false;
		els_mcdr_isLoaded = false;
		jets_mc3dr_isLoaded = false;
		jets_mcdr_isLoaded = false;
		jets_mc_emEnergy_isLoaded = false;
		jets_mc_gpdr_isLoaded = false;
		jets_mc_hadEnergy_isLoaded = false;
		jets_mc_invEnergy_isLoaded = false;
		jets_mc_otherEnergy_isLoaded = false;
		mus_mc3dr_isLoaded = false;
		mus_mcdr_isLoaded = false;
		pfjets_mc3dr_isLoaded = false;
		pfjets_mcdr_isLoaded = false;
		pfjets_mc_emEnergy_isLoaded = false;
		pfjets_mc_gpdr_isLoaded = false;
		pfjets_mc_hadEnergy_isLoaded = false;
		pfjets_mc_invEnergy_isLoaded = false;
		pfjets_mc_otherEnergy_isLoaded = false;
		photons_mc3dr_isLoaded = false;
		photons_mcdr_isLoaded = false;
		trk_mc3dr_isLoaded = false;
		trk_mcdr_isLoaded = false;
		trks_conv_dcot_isLoaded = false;
		trks_conv_dist_isLoaded = false;
		els_ecalJuraIso_isLoaded = false;
		els_ecalJuraTowerIso_isLoaded = false;
		els_hcalConeIso_isLoaded = false;
		els_tkJuraIso_isLoaded = false;
		els_jetdr_isLoaded = false;
		els_musdr_isLoaded = false;
		els_chi2_isLoaded = false;
		els_conv_dcot_isLoaded = false;
		els_conv_dist_isLoaded = false;
		els_conv_radius_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_d0corr_isLoaded = false;
		els_dEtaIn_isLoaded = false;
		els_dEtaOut_isLoaded = false;
		els_dPhiIn_isLoaded = false;
		els_dPhiInPhiOut_isLoaded = false;
		els_dPhiOut_isLoaded = false;
		els_deltaEtaEleClusterTrackAtCalo_isLoaded = false;
		els_deltaPhiEleClusterTrackAtCalo_isLoaded = false;
		els_e1x5_isLoaded = false;
		els_e2x5Max_isLoaded = false;
		els_e3x3_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_eMax_isLoaded = false;
		els_eOverPIn_isLoaded = false;
		els_eOverPOut_isLoaded = false;
		els_eSC_isLoaded = false;
		els_eSCPresh_isLoaded = false;
		els_eSCRaw_isLoaded = false;
		els_eSeed_isLoaded = false;
		els_eSeedOverPIn_isLoaded = false;
		els_eSeedOverPOut_isLoaded = false;
		els_ecalEnergy_isLoaded = false;
		els_ecalEnergyError_isLoaded = false;
		els_ecalIso_isLoaded = false;
		els_ecalIso04_isLoaded = false;
		els_egamma_looseId_isLoaded = false;
		els_egamma_robustHighEnergy_isLoaded = false;
		els_egamma_robustLooseId_isLoaded = false;
		els_egamma_robustTightId_isLoaded = false;
		els_egamma_tightId_isLoaded = false;
		els_electronMomentumError_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_etaSC_isLoaded = false;
		els_fbrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_hcalDepth1OverEcal_isLoaded = false;
		els_hcalDepth1TowerSumEt_isLoaded = false;
		els_hcalDepth1TowerSumEt04_isLoaded = false;
		els_hcalDepth2OverEcal_isLoaded = false;
		els_hcalDepth2TowerSumEt_isLoaded = false;
		els_hcalDepth2TowerSumEt04_isLoaded = false;
		els_hcalIso_isLoaded = false;
		els_hcalIso04_isLoaded = false;
		els_layer1_charge_isLoaded = false;
		els_mva_isLoaded = false;
		els_ndof_isLoaded = false;
		els_phiErr_isLoaded = false;
		els_phiSC_isLoaded = false;
		els_ptErr_isLoaded = false;
		els_sigmaEtaEta_isLoaded = false;
		els_sigmaIEtaIEta_isLoaded = false;
		els_sigmaIEtaIEtaSC_isLoaded = false;
		els_sigmaIPhiIPhi_isLoaded = false;
		els_sigmaIPhiIPhiSC_isLoaded = false;
		els_sigmaPhiPhi_isLoaded = false;
		els_tkIso_isLoaded = false;
		els_tkIso04_isLoaded = false;
		els_trackMomentumError_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_trkshFrac_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		els_z0corr_isLoaded = false;
		gsftrks_chi2_isLoaded = false;
		gsftrks_d0_isLoaded = false;
		gsftrks_d0Err_isLoaded = false;
		gsftrks_d0corr_isLoaded = false;
		gsftrks_d0corrPhi_isLoaded = false;
		gsftrks_d0phiCov_isLoaded = false;
		gsftrks_etaErr_isLoaded = false;
		gsftrks_layer1_charge_isLoaded = false;
		gsftrks_ndof_isLoaded = false;
		gsftrks_phiErr_isLoaded = false;
		gsftrks_ptErr_isLoaded = false;
		gsftrks_z0_isLoaded = false;
		gsftrks_z0Err_isLoaded = false;
		gsftrks_z0corr_isLoaded = false;
		hyp_Ht_isLoaded = false;
		hyp_dPhi_nJet_metMuonJESCorr_isLoaded = false;
		hyp_dPhi_nJet_muCorrMet_isLoaded = false;
		hyp_dPhi_nJet_tcMet_isLoaded = false;
		hyp_dPhi_nJet_unCorrMet_isLoaded = false;
		hyp_ll_chi2_isLoaded = false;
		hyp_ll_d0_isLoaded = false;
		hyp_ll_d0Err_isLoaded = false;
		hyp_ll_d0corr_isLoaded = false;
		hyp_ll_dPhi_metMuonJESCorr_isLoaded = false;
		hyp_ll_dPhi_muCorrMet_isLoaded = false;
		hyp_ll_dPhi_tcMet_isLoaded = false;
		hyp_ll_dPhi_unCorrMet_isLoaded = false;
		hyp_ll_etaErr_isLoaded = false;
		hyp_ll_ndof_isLoaded = false;
		hyp_ll_phiErr_isLoaded = false;
		hyp_ll_ptErr_isLoaded = false;
		hyp_ll_z0_isLoaded = false;
		hyp_ll_z0Err_isLoaded = false;
		hyp_ll_z0corr_isLoaded = false;
		hyp_lt_chi2_isLoaded = false;
		hyp_lt_d0_isLoaded = false;
		hyp_lt_d0Err_isLoaded = false;
		hyp_lt_d0corr_isLoaded = false;
		hyp_lt_dPhi_metMuonJESCorr_isLoaded = false;
		hyp_lt_dPhi_muCorrMet_isLoaded = false;
		hyp_lt_dPhi_tcMet_isLoaded = false;
		hyp_lt_dPhi_unCorrMet_isLoaded = false;
		hyp_lt_etaErr_isLoaded = false;
		hyp_lt_ndof_isLoaded = false;
		hyp_lt_phiErr_isLoaded = false;
		hyp_lt_ptErr_isLoaded = false;
		hyp_lt_z0_isLoaded = false;
		hyp_lt_z0Err_isLoaded = false;
		hyp_lt_z0corr_isLoaded = false;
		hyp_mt2_metMuonJESCorr_isLoaded = false;
		hyp_mt2_muCorrMet_isLoaded = false;
		hyp_mt2_tcMet_isLoaded = false;
		hyp_sumJetPt_isLoaded = false;
		hyp_FVFit_chi2ndf_isLoaded = false;
		hyp_FVFit_prob_isLoaded = false;
		hyp_FVFit_v4cxx_isLoaded = false;
		hyp_FVFit_v4cxy_isLoaded = false;
		hyp_FVFit_v4cyy_isLoaded = false;
		hyp_FVFit_v4czx_isLoaded = false;
		hyp_FVFit_v4czy_isLoaded = false;
		hyp_FVFit_v4czz_isLoaded = false;
		hyp_ll_ecaliso_isLoaded = false;
		hyp_ll_trkiso_isLoaded = false;
		hyp_lt_ecaliso_isLoaded = false;
		hyp_lt_trkiso_isLoaded = false;
		jets_approximatefHPD_isLoaded = false;
		jets_approximatefRBX_isLoaded = false;
		jets_cor_isLoaded = false;
		jets_emFrac_isLoaded = false;
		jets_fHPD_isLoaded = false;
		jets_fRBX_isLoaded = false;
		jets_fSubDetector1_isLoaded = false;
		jets_fSubDetector2_isLoaded = false;
		jets_fSubDetector3_isLoaded = false;
		jets_fSubDetector4_isLoaded = false;
		jets_hitsInN90_isLoaded = false;
		jets_n90Hits_isLoaded = false;
		jets_nECALTowers_isLoaded = false;
		jets_nHCALTowers_isLoaded = false;
		jets_restrictedEMF_isLoaded = false;
		jpts_cor_isLoaded = false;
		jpts_emFrac_isLoaded = false;
		mus_met_deltax_isLoaded = false;
		mus_met_deltay_isLoaded = false;
		mus_eledr_isLoaded = false;
		mus_jetdr_isLoaded = false;
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
		mus_gfit_d0_isLoaded = false;
		mus_gfit_d0Err_isLoaded = false;
		mus_gfit_d0corr_isLoaded = false;
		mus_gfit_ndof_isLoaded = false;
		mus_gfit_qoverp_isLoaded = false;
		mus_gfit_qoverpError_isLoaded = false;
		mus_gfit_z0_isLoaded = false;
		mus_gfit_z0Err_isLoaded = false;
		mus_gfit_z0corr_isLoaded = false;
		mus_iso03_emEt_isLoaded = false;
		mus_iso03_hadEt_isLoaded = false;
		mus_iso03_hoEt_isLoaded = false;
		mus_iso03_sumPt_isLoaded = false;
		mus_iso05_emEt_isLoaded = false;
		mus_iso05_hadEt_isLoaded = false;
		mus_iso05_hoEt_isLoaded = false;
		mus_iso05_sumPt_isLoaded = false;
		mus_iso_ecalvetoDep_isLoaded = false;
		mus_iso_hcalvetoDep_isLoaded = false;
		mus_iso_hovetoDep_isLoaded = false;
		mus_iso_trckvetoDep_isLoaded = false;
		mus_ndof_isLoaded = false;
		mus_phiErr_isLoaded = false;
		mus_ptErr_isLoaded = false;
		mus_qoverp_isLoaded = false;
		mus_qoverpError_isLoaded = false;
		mus_sta_chi2_isLoaded = false;
		mus_sta_d0_isLoaded = false;
		mus_sta_d0Err_isLoaded = false;
		mus_sta_d0corr_isLoaded = false;
		mus_sta_ndof_isLoaded = false;
		mus_sta_qoverp_isLoaded = false;
		mus_sta_qoverpError_isLoaded = false;
		mus_sta_z0_isLoaded = false;
		mus_sta_z0Err_isLoaded = false;
		mus_sta_z0corr_isLoaded = false;
		mus_timeAtIpInOut_isLoaded = false;
		mus_timeAtIpInOutErr_isLoaded = false;
		mus_timeAtIpOutIn_isLoaded = false;
		mus_timeAtIpOutInErr_isLoaded = false;
		mus_vertexphi_isLoaded = false;
		mus_z0_isLoaded = false;
		mus_z0Err_isLoaded = false;
		mus_z0corr_isLoaded = false;
		els_pat_caloIso_isLoaded = false;
		els_pat_ecalIso_isLoaded = false;
		els_pat_hcalIso_isLoaded = false;
		els_pat_looseId_isLoaded = false;
		els_pat_robustHighEnergy_isLoaded = false;
		els_pat_robustLooseId_isLoaded = false;
		els_pat_robustTightId_isLoaded = false;
		els_pat_scE1x5_isLoaded = false;
		els_pat_scE2x5Max_isLoaded = false;
		els_pat_scE5x5_isLoaded = false;
		els_pat_sigmaEtaEta_isLoaded = false;
		els_pat_sigmaIEtaIEta_isLoaded = false;
		els_pat_tightId_isLoaded = false;
		els_pat_trackIso_isLoaded = false;
		jets_pat_combinedSecondaryVertexBJetTag_isLoaded = false;
		jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		jets_pat_coneIsolationTauJetTag_isLoaded = false;
		jets_pat_impactParameterMVABJetTag_isLoaded = false;
		jets_pat_jetBProbabilityBJetTag_isLoaded = false;
		jets_pat_jetCharge_isLoaded = false;
		jets_pat_jetProbabilityBJetTag_isLoaded = false;
		jets_pat_noCorrF_isLoaded = false;
		jets_pat_simpleSecondaryVertexHighEffBJetTag_isLoaded = false;
		jets_pat_simpleSecondaryVertexHighPurBJetTag_isLoaded = false;
		jets_pat_softElectronByIP3dBJetTag_isLoaded = false;
		jets_pat_softElectronByPtBJetTag_isLoaded = false;
		jets_pat_softMuonBJetTag_isLoaded = false;
		jets_pat_softMuonByIP3dBJetTag_isLoaded = false;
		jets_pat_softMuonByPtBJetTag_isLoaded = false;
		jets_pat_trackCountingHighEffBJetTag_isLoaded = false;
		jets_pat_trackCountingHighPurBJetTag_isLoaded = false;
		mus_pat_caloIso_isLoaded = false;
		mus_pat_calovetoDep_isLoaded = false;
		mus_pat_ecalIso_isLoaded = false;
		mus_pat_ecalvetoDep_isLoaded = false;
		mus_pat_hcalIso_isLoaded = false;
		mus_pat_hcalvetoDep_isLoaded = false;
		mus_pat_trackIso_isLoaded = false;
		mus_pat_trckvetoDep_isLoaded = false;
		pfels_deltaP_isLoaded = false;
		pfels_ecalE_isLoaded = false;
		pfels_hcalE_isLoaded = false;
		pfels_isoChargedHadrons_isLoaded = false;
		pfels_isoNeutralHadrons_isLoaded = false;
		pfels_isoPhotons_isLoaded = false;
		pfels_mva_emu_isLoaded = false;
		pfels_mva_epi_isLoaded = false;
		pfels_mva_nothing_gamma_isLoaded = false;
		pfels_mva_nothing_nh_isLoaded = false;
		pfels_mva_pimu_isLoaded = false;
		pfels_pS1E_isLoaded = false;
		pfels_pS2E_isLoaded = false;
		pfels_rawEcalE_isLoaded = false;
		pfels_rawHcalE_isLoaded = false;
		pfjets_chargedEmE_isLoaded = false;
		pfjets_chargedHadronE_isLoaded = false;
		pfjets_cor_isLoaded = false;
		pfjets_neutralEmE_isLoaded = false;
		pfjets_neutralHadronE_isLoaded = false;
		pfmus_deltaP_isLoaded = false;
		pfmus_ecalE_isLoaded = false;
		pfmus_hcalE_isLoaded = false;
		pfmus_isoChargedHadrons_isLoaded = false;
		pfmus_isoNeutralHadrons_isLoaded = false;
		pfmus_isoPhotons_isLoaded = false;
		pfmus_mva_emu_isLoaded = false;
		pfmus_mva_epi_isLoaded = false;
		pfmus_mva_nothing_gamma_isLoaded = false;
		pfmus_mva_nothing_nh_isLoaded = false;
		pfmus_mva_pimu_isLoaded = false;
		pfmus_pS1E_isLoaded = false;
		pfmus_pS2E_isLoaded = false;
		pfmus_rawEcalE_isLoaded = false;
		pfmus_rawHcalE_isLoaded = false;
		photons_e1x5_isLoaded = false;
		photons_e2x5Max_isLoaded = false;
		photons_e3x3_isLoaded = false;
		photons_e5x5_isLoaded = false;
		photons_ecalIso03_isLoaded = false;
		photons_ecalIso04_isLoaded = false;
		photons_hOverE_isLoaded = false;
		photons_hcalIso03_isLoaded = false;
		photons_hcalIso04_isLoaded = false;
		photons_ntkIsoHollow03_isLoaded = false;
		photons_ntkIsoHollow04_isLoaded = false;
		photons_ntkIsoSolid03_isLoaded = false;
		photons_ntkIsoSolid04_isLoaded = false;
		photons_sigmaEtaEta_isLoaded = false;
		photons_sigmaIEtaIEta_isLoaded = false;
		photons_swissSeed_isLoaded = false;
		photons_tkIsoHollow03_isLoaded = false;
		photons_tkIsoHollow04_isLoaded = false;
		photons_tkIsoSolid03_isLoaded = false;
		photons_tkIsoSolid04_isLoaded = false;
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
		scs_sigmaIEtaIEtaSC_isLoaded = false;
		scs_sigmaIEtaIPhi_isLoaded = false;
		scs_sigmaIEtaIPhiSC_isLoaded = false;
		scs_sigmaIPhiIPhi_isLoaded = false;
		scs_sigmaIPhiIPhiSC_isLoaded = false;
		scs_sigmaPhiPhi_isLoaded = false;
		scs_timeSeed_isLoaded = false;
		svs_anglePV_isLoaded = false;
		svs_chi2_isLoaded = false;
		svs_dist3Dsig_isLoaded = false;
		svs_dist3Dval_isLoaded = false;
		svs_distXYsig_isLoaded = false;
		svs_distXYval_isLoaded = false;
		svs_ndof_isLoaded = false;
		svs_prob_isLoaded = false;
		svs_xError_isLoaded = false;
		svs_yError_isLoaded = false;
		svs_zError_isLoaded = false;
		mus_tcmet_deltax_isLoaded = false;
		mus_tcmet_deltay_isLoaded = false;
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_d0corr_isLoaded = false;
		trks_d0corrPhi_isLoaded = false;
		trks_d0phiCov_isLoaded = false;
		trks_etaErr_isLoaded = false;
		trks_layer1_charge_isLoaded = false;
		trks_ndof_isLoaded = false;
		trks_phiErr_isLoaded = false;
		trks_ptErr_isLoaded = false;
		trks_z0_isLoaded = false;
		trks_z0Err_isLoaded = false;
		trks_z0corr_isLoaded = false;
		trkjets_cor_isLoaded = false;
		trks_d0Errvtx_isLoaded = false;
		trks_d0vtx_isLoaded = false;
		vtxs_chi2_isLoaded = false;
		vtxs_ndof_isLoaded = false;
		vtxs_sumpt_isLoaded = false;
		vtxs_xError_isLoaded = false;
		vtxs_yError_isLoaded = false;
		vtxs_zError_isLoaded = false;
		vtxs_covMatrix_isLoaded = false;
		evt_cscLooseHaloId_isLoaded = false;
		evt_cscTightHaloId_isLoaded = false;
		evt_ecalLooseHaloId_isLoaded = false;
		evt_ecalTightHaloId_isLoaded = false;
		evt_extremeTightHaloId_isLoaded = false;
		evt_globalLooseHaloId_isLoaded = false;
		evt_globalTightHaloId_isLoaded = false;
		evt_hcalLooseHaloId_isLoaded = false;
		evt_hcalTightHaloId_isLoaded = false;
		evt_looseHaloId_isLoaded = false;
		evt_nHaloLikeTracks_isLoaded = false;
		evt_nHaloTriggerCandidates_isLoaded = false;
		evt_tightHaloId_isLoaded = false;
		evt_bsType_isLoaded = false;
		evt_bunchCrossing_isLoaded = false;
		evt_experimentType_isLoaded = false;
		evt_isRealData_isLoaded = false;
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
		pdfinfo_id1_isLoaded = false;
		pdfinfo_id2_isLoaded = false;
		evt_ecaliPhiSuspects_isLoaded = false;
		evt_globaliPhiSuspects_isLoaded = false;
		evt_hcaliPhiSuspects_isLoaded = false;
		els_mc3_id_isLoaded = false;
		els_mc3idx_isLoaded = false;
		els_mc3_motherid_isLoaded = false;
		els_mc3_motheridx_isLoaded = false;
		els_mc_id_isLoaded = false;
		els_mcidx_isLoaded = false;
		els_mc_motherid_isLoaded = false;
		jets_mc3_id_isLoaded = false;
		jets_mc3idx_isLoaded = false;
		jets_mc_gpidx_isLoaded = false;
		jets_mc_id_isLoaded = false;
		jets_mcidx_isLoaded = false;
		jets_mc_motherid_isLoaded = false;
		mus_mc3_id_isLoaded = false;
		mus_mc3idx_isLoaded = false;
		mus_mc3_motherid_isLoaded = false;
		mus_mc3_motheridx_isLoaded = false;
		mus_mc_id_isLoaded = false;
		mus_mcidx_isLoaded = false;
		mus_mc_motherid_isLoaded = false;
		pfjets_mc3_id_isLoaded = false;
		pfjets_mc3idx_isLoaded = false;
		pfjets_mc_gpidx_isLoaded = false;
		pfjets_mc_id_isLoaded = false;
		pfjets_mcidx_isLoaded = false;
		pfjets_mc_motherid_isLoaded = false;
		photons_mc3_id_isLoaded = false;
		photons_mc3idx_isLoaded = false;
		photons_mc3_motherid_isLoaded = false;
		photons_mc3_motheridx_isLoaded = false;
		photons_mc_id_isLoaded = false;
		photons_mcidx_isLoaded = false;
		photons_mc_motherid_isLoaded = false;
		trk_mc3_id_isLoaded = false;
		trk_mc3idx_isLoaded = false;
		trk_mc3_motherid_isLoaded = false;
		trk_mc3_motheridx_isLoaded = false;
		trk_mc_id_isLoaded = false;
		trk_mcidx_isLoaded = false;
		trk_mc_motherid_isLoaded = false;
		trks_conv_tkidx_isLoaded = false;
		els_exp_innerlayers39X_isLoaded = false;
		els_closestJet_isLoaded = false;
		els_closestMuon_isLoaded = false;
		els_pfelsidx_isLoaded = false;
		els_category_isLoaded = false;
		els_charge_isLoaded = false;
		els_class_isLoaded = false;
		els_conv_tkidx_isLoaded = false;
		els_exp_innerlayers_isLoaded = false;
		els_exp_outerlayers_isLoaded = false;
		els_fiduciality_isLoaded = false;
		els_gsftrkidx_isLoaded = false;
		els_layer1_det_isLoaded = false;
		els_layer1_layer_isLoaded = false;
		els_layer1_sizerphi_isLoaded = false;
		els_layer1_sizerz_isLoaded = false;
		els_lostHits_isLoaded = false;
		els_lost_pixelhits_isLoaded = false;
		els_nSeed_isLoaded = false;
		els_sccharge_isLoaded = false;
		els_scindex_isLoaded = false;
		els_trk_charge_isLoaded = false;
		els_trkidx_isLoaded = false;
		els_type_isLoaded = false;
		els_validHits_isLoaded = false;
		els_valid_pixelhits_isLoaded = false;
		genps_id_isLoaded = false;
		genps_id_mother_isLoaded = false;
		genps_status_isLoaded = false;
		gsftrks_charge_isLoaded = false;
		gsftrks_exp_innerlayers_isLoaded = false;
		gsftrks_exp_outerlayers_isLoaded = false;
		gsftrks_layer1_det_isLoaded = false;
		gsftrks_layer1_layer_isLoaded = false;
		gsftrks_layer1_sizerphi_isLoaded = false;
		gsftrks_layer1_sizerz_isLoaded = false;
		gsftrks_lostHits_isLoaded = false;
		gsftrks_lost_pixelhits_isLoaded = false;
		gsftrks_nlayers_isLoaded = false;
		gsftrks_nlayers3D_isLoaded = false;
		gsftrks_nlayersLost_isLoaded = false;
		gsftrks_validHits_isLoaded = false;
		gsftrks_valid_pixelhits_isLoaded = false;
		hyp_ll_charge_isLoaded = false;
		hyp_ll_id_isLoaded = false;
		hyp_ll_index_isLoaded = false;
		hyp_ll_lostHits_isLoaded = false;
		hyp_ll_validHits_isLoaded = false;
		hyp_lt_charge_isLoaded = false;
		hyp_lt_id_isLoaded = false;
		hyp_lt_index_isLoaded = false;
		hyp_lt_lostHits_isLoaded = false;
		hyp_lt_validHits_isLoaded = false;
		hyp_njets_isLoaded = false;
		hyp_nojets_isLoaded = false;
		hyp_type_isLoaded = false;
		hyp_FVFit_ndf_isLoaded = false;
		hyp_FVFit_status_isLoaded = false;
		hyp_ll_mc_id_isLoaded = false;
		hyp_ll_mc_motherid_isLoaded = false;
		hyp_lt_mc_id_isLoaded = false;
		hyp_lt_mc_motherid_isLoaded = false;
		hyp_quadlep_first_type_isLoaded = false;
		hyp_quadlep_fourth_type_isLoaded = false;
		hyp_quadlep_second_type_isLoaded = false;
		hyp_quadlep_third_type_isLoaded = false;
		hyp_trilep_first_type_isLoaded = false;
		hyp_trilep_second_type_isLoaded = false;
		hyp_trilep_third_type_isLoaded = false;
		jets_closestElectron_isLoaded = false;
		jets_closestMuon_isLoaded = false;
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
		mus_met_flag_isLoaded = false;
		mus_closestEle_isLoaded = false;
		mus_closestJet_isLoaded = false;
		mus_pfmusidx_isLoaded = false;
		mus_charge_isLoaded = false;
		mus_chi2LocalMomentum_isLoaded = false;
		mus_chi2LocalPosition_isLoaded = false;
		mus_gfit_validHits_isLoaded = false;
		mus_gfit_validSTAHits_isLoaded = false;
		mus_gfit_validSiHits_isLoaded = false;
		mus_glbKink_isLoaded = false;
		mus_glbTrackProbability_isLoaded = false;
		mus_globalDeltaEtaPhi_isLoaded = false;
		mus_goodmask_isLoaded = false;
		mus_iso03_ntrk_isLoaded = false;
		mus_iso05_ntrk_isLoaded = false;
		mus_localDistance_isLoaded = false;
		mus_lostHits_isLoaded = false;
		mus_nOverlaps_isLoaded = false;
		mus_nmatches_isLoaded = false;
		mus_overlap0_isLoaded = false;
		mus_overlap1_isLoaded = false;
		mus_pid_TM2DCompatibilityLoose_isLoaded = false;
		mus_pid_TM2DCompatibilityTight_isLoaded = false;
		mus_pid_TMLastStationLoose_isLoaded = false;
		mus_pid_TMLastStationTight_isLoaded = false;
		mus_staRelChi2_isLoaded = false;
		mus_sta_validHits_isLoaded = false;
		mus_timeDirection_isLoaded = false;
		mus_timeNumStationsUsed_isLoaded = false;
		mus_trkKink_isLoaded = false;
		mus_trkRelChi2_isLoaded = false;
		mus_trk_charge_isLoaded = false;
		mus_trkidx_isLoaded = false;
		mus_type_isLoaded = false;
		mus_validHits_isLoaded = false;
		els_pat_genID_isLoaded = false;
		els_pat_genMotherID_isLoaded = false;
		jets_pat_genPartonMother_id_isLoaded = false;
		jets_pat_genParton_id_isLoaded = false;
		jets_pat_jetIDLoose_isLoaded = false;
		jets_pat_jetIDLooseAOD_isLoaded = false;
		jets_pat_jetIDMinimal_isLoaded = false;
		jets_pat_jetIDTight_isLoaded = false;
		jets_pat_partonFlavour_isLoaded = false;
		mus_pat_genID_isLoaded = false;
		mus_pat_genMotherID_isLoaded = false;
		pfels_elsidx_isLoaded = false;
		pfels_charge_isLoaded = false;
		pfels_flag_isLoaded = false;
		pfels_particleId_isLoaded = false;
		pfjets_chargedMultiplicity_isLoaded = false;
		pfjets_muonMultiplicity_isLoaded = false;
		pfjets_neutralMultiplicity_isLoaded = false;
		pfmus_musidx_isLoaded = false;
		pfmus_charge_isLoaded = false;
		pfmus_flag_isLoaded = false;
		pfmus_particleId_isLoaded = false;
		photons_fiduciality_isLoaded = false;
		photons_scindex_isLoaded = false;
		scs_detIdSeed_isLoaded = false;
		scs_elsidx_isLoaded = false;
		scs_severitySeed_isLoaded = false;
		svs_isKs_isLoaded = false;
		svs_isLambda_isLoaded = false;
		svs_mc3_id_isLoaded = false;
		svs_nTrks_isLoaded = false;
		mus_tcmet_flag_isLoaded = false;
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
		trks_elsidx_isLoaded = false;
		trk_musidx_isLoaded = false;
		trkjets_ntrks_isLoaded = false;
		trkjets_vtxs_idx_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		vtxs_tracksSize_isLoaded = false;
		genps_lepdaughter_id_isLoaded = false;
		genps_lepdaughter_idx_isLoaded = false;
		hlt_trigObjs_id_isLoaded = false;
		hyp_jets_idx_isLoaded = false;
		hyp_other_jets_idx_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_detectorStatus_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		genps_flavorHistoryFilterResult_isLoaded = false;
		evt_ngenjets_isLoaded = false;
		genps_signalProcessID_isLoaded = false;
		hlt_bits1_isLoaded = false;
		hlt_bits2_isLoaded = false;
		hlt_bits3_isLoaded = false;
		hlt_bits4_isLoaded = false;
		hlt_bits5_isLoaded = false;
		hlt_bits6_isLoaded = false;
		hlt_bits7_isLoaded = false;
		hlt_bits8_isLoaded = false;
		evt_njets_isLoaded = false;
		evt_njpts_isLoaded = false;
		l1_bits1_isLoaded = false;
		l1_bits2_isLoaded = false;
		l1_bits3_isLoaded = false;
		l1_bits4_isLoaded = false;
		l1_techbits1_isLoaded = false;
		l1_techbits2_isLoaded = false;
		evt_nphotons_isLoaded = false;
		evt_ecalRecoStatus_isLoaded = false;
		evt_nscs_isLoaded = false;
		evt_ntrkjets_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		hlt_prescales_isLoaded = false;
		hyp_quadlep_bucket_isLoaded = false;
		hyp_quadlep_first_index_isLoaded = false;
		hyp_quadlep_fourth_index_isLoaded = false;
		hyp_quadlep_second_index_isLoaded = false;
		hyp_quadlep_third_index_isLoaded = false;
		hyp_trilep_bucket_isLoaded = false;
		hyp_trilep_first_index_isLoaded = false;
		hyp_trilep_second_index_isLoaded = false;
		hyp_trilep_third_index_isLoaded = false;
		l1_prescales_isLoaded = false;
		l1_techtrigprescales_isLoaded = false;
		els_pat_flag_isLoaded = false;
		jets_pat_flag_isLoaded = false;
		mus_pat_flag_isLoaded = false;
		evt_nEvts_isLoaded = false;
		evt_filt_eff_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (evt_CMS2tag_branch != 0) evt_CMS2tag();
	if (evt_dataset_branch != 0) evt_dataset();
	if (hlt_trigNames_branch != 0) hlt_trigNames();
	if (l1_techtrigNames_branch != 0) l1_techtrigNames();
	if (l1_trigNames_branch != 0) l1_trigNames();
	if (evt_errCategory_branch != 0) evt_errCategory();
	if (evt_errModule_branch != 0) evt_errModule();
	if (evt_errSeverity_branch != 0) evt_errSeverity();
	if (evt_eventHasHalo_branch != 0) evt_eventHasHalo();
	if (evt_hbheFilter_branch != 0) evt_hbheFilter();
	if (mus_tightMatch_branch != 0) mus_tightMatch();
	if (mus_updatedSta_branch != 0) mus_updatedSta();
	if (photons_haspixelSeed_branch != 0) photons_haspixelSeed();
	if (jets_closestElectron_DR_branch != 0) jets_closestElectron_DR();
	if (jets_closestMuon_DR_branch != 0) jets_closestMuon_DR();
	if (evt_bs_Xwidth_branch != 0) evt_bs_Xwidth();
	if (evt_bs_XwidthErr_branch != 0) evt_bs_XwidthErr();
	if (evt_bs_Ywidth_branch != 0) evt_bs_Ywidth();
	if (evt_bs_YwidthErr_branch != 0) evt_bs_YwidthErr();
	if (evt_bs_dxdz_branch != 0) evt_bs_dxdz();
	if (evt_bs_dxdzErr_branch != 0) evt_bs_dxdzErr();
	if (evt_bs_dydz_branch != 0) evt_bs_dydz();
	if (evt_bs_dydzErr_branch != 0) evt_bs_dydzErr();
	if (evt_bs_sigmaZ_branch != 0) evt_bs_sigmaZ();
	if (evt_bs_sigmaZErr_branch != 0) evt_bs_sigmaZErr();
	if (evt_bs_xErr_branch != 0) evt_bs_xErr();
	if (evt_bs_yErr_branch != 0) evt_bs_yErr();
	if (evt_bs_zErr_branch != 0) evt_bs_zErr();
	if (evthcal_dmetx_branch != 0) evthcal_dmetx();
	if (evthcal_dmety_branch != 0) evthcal_dmety();
	if (evthcal_dsumet_branch != 0) evthcal_dsumet();
	if (evthf_dmetx_branch != 0) evthf_dmetx();
	if (evthf_dmety_branch != 0) evthf_dmety();
	if (evthf_dsumet_branch != 0) evthf_dsumet();
	if (evt_bField_branch != 0) evt_bField();
	if (evt_kfactor_branch != 0) evt_kfactor();
	if (evt_scale1fb_branch != 0) evt_scale1fb();
	if (evt_xsec_excl_branch != 0) evt_xsec_excl();
	if (evt_xsec_incl_branch != 0) evt_xsec_incl();
	if (gen_met_branch != 0) gen_met();
	if (gen_metPhi_branch != 0) gen_metPhi();
	if (genps_alphaQCD_branch != 0) genps_alphaQCD();
	if (genps_pthat_branch != 0) genps_pthat();
	if (genps_qScale_branch != 0) genps_qScale();
	if (genps_weight_branch != 0) genps_weight();
	if (gen_sumEt_branch != 0) gen_sumEt();
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
	if (evt_ecalendcapm_met_branch != 0) evt_ecalendcapm_met();
	if (evt_ecalendcapm_metPhi_branch != 0) evt_ecalendcapm_metPhi();
	if (evt_ecalendcapp_met_branch != 0) evt_ecalendcapp_met();
	if (evt_ecalendcapp_metPhi_branch != 0) evt_ecalendcapp_metPhi();
	if (evt_ecalmet_branch != 0) evt_ecalmet();
	if (evt_ecalmetPhi_branch != 0) evt_ecalmetPhi();
	if (evt_endcapm_met_branch != 0) evt_endcapm_met();
	if (evt_endcapm_metPhi_branch != 0) evt_endcapm_metPhi();
	if (evt_endcapp_met_branch != 0) evt_endcapp_met();
	if (evt_endcapp_metPhi_branch != 0) evt_endcapp_metPhi();
	if (evt_hcalendcapm_met_branch != 0) evt_hcalendcapm_met();
	if (evt_hcalendcapm_metPhi_branch != 0) evt_hcalendcapm_metPhi();
	if (evt_hcalendcapp_met_branch != 0) evt_hcalendcapp_met();
	if (evt_hcalendcapp_metPhi_branch != 0) evt_hcalendcapp_metPhi();
	if (evt_hcalmet_branch != 0) evt_hcalmet();
	if (evt_hcalmetPhi_branch != 0) evt_hcalmetPhi();
	if (evt_met_branch != 0) evt_met();
	if (evt_metHO_branch != 0) evt_metHO();
	if (evt_metHOPhi_branch != 0) evt_metHOPhi();
	if (evt_metHOSig_branch != 0) evt_metHOSig();
	if (evt_metMuonCorr_branch != 0) evt_metMuonCorr();
	if (evt_metMuonCorrPhi_branch != 0) evt_metMuonCorrPhi();
	if (evt_metMuonCorrSig_branch != 0) evt_metMuonCorrSig();
	if (evt_metMuonJESCorr_branch != 0) evt_metMuonJESCorr();
	if (evt_metMuonJESCorrPhi_branch != 0) evt_metMuonJESCorrPhi();
	if (evt_metMuonJESCorrSig_branch != 0) evt_metMuonJESCorrSig();
	if (evt_metNoHF_branch != 0) evt_metNoHF();
	if (evt_metNoHFHO_branch != 0) evt_metNoHFHO();
	if (evt_metNoHFHOPhi_branch != 0) evt_metNoHFHOPhi();
	if (evt_metNoHFHOSig_branch != 0) evt_metNoHFHOSig();
	if (evt_metNoHFPhi_branch != 0) evt_metNoHFPhi();
	if (evt_metNoHFSig_branch != 0) evt_metNoHFSig();
	if (evt_metOpt_branch != 0) evt_metOpt();
	if (evt_metOptHO_branch != 0) evt_metOptHO();
	if (evt_metOptHOPhi_branch != 0) evt_metOptHOPhi();
	if (evt_metOptHOSig_branch != 0) evt_metOptHOSig();
	if (evt_metOptNoHF_branch != 0) evt_metOptNoHF();
	if (evt_metOptNoHFHO_branch != 0) evt_metOptNoHFHO();
	if (evt_metOptNoHFHOPhi_branch != 0) evt_metOptNoHFHOPhi();
	if (evt_metOptNoHFHOSig_branch != 0) evt_metOptNoHFHOSig();
	if (evt_metOptNoHFPhi_branch != 0) evt_metOptNoHFPhi();
	if (evt_metOptNoHFSig_branch != 0) evt_metOptNoHFSig();
	if (evt_metOptPhi_branch != 0) evt_metOptPhi();
	if (evt_metOptSig_branch != 0) evt_metOptSig();
	if (evt_metPhi_branch != 0) evt_metPhi();
	if (evt_metSig_branch != 0) evt_metSig();
	if (evt_sumet_branch != 0) evt_sumet();
	if (evt_sumetHO_branch != 0) evt_sumetHO();
	if (evt_sumetMuonCorr_branch != 0) evt_sumetMuonCorr();
	if (evt_sumetNoHF_branch != 0) evt_sumetNoHF();
	if (evt_sumetNoHFHO_branch != 0) evt_sumetNoHFHO();
	if (evt_sumetOpt_branch != 0) evt_sumetOpt();
	if (evt_sumetOptHO_branch != 0) evt_sumetOptHO();
	if (evt_sumetOptNoHF_branch != 0) evt_sumetOptNoHF();
	if (evt_sumetOptNoHFHO_branch != 0) evt_sumetOptNoHFHO();
	if (met_pat_metCor_branch != 0) met_pat_metCor();
	if (met_pat_metPhiCor_branch != 0) met_pat_metPhiCor();
	if (met_pat_metPhiUncor_branch != 0) met_pat_metPhiUncor();
	if (met_pat_metPhiUncorJES_branch != 0) met_pat_metPhiUncorJES();
	if (met_pat_metPhiUncorMuon_branch != 0) met_pat_metPhiUncorMuon();
	if (met_pat_metUncor_branch != 0) met_pat_metUncor();
	if (met_pat_metUncorJES_branch != 0) met_pat_metUncorJES();
	if (met_pat_metUncorMuon_branch != 0) met_pat_metUncorMuon();
	if (pdfinfo_scale_branch != 0) pdfinfo_scale();
	if (pdfinfo_x1_branch != 0) pdfinfo_x1();
	if (pdfinfo_x2_branch != 0) pdfinfo_x2();
	if (evt_pfmet_branch != 0) evt_pfmet();
	if (evt_pfmetPhi_branch != 0) evt_pfmetPhi();
	if (evt_pfmetSig_branch != 0) evt_pfmetSig();
	if (evt_pfsumet_branch != 0) evt_pfsumet();
	if (evt_tcmet_branch != 0) evt_tcmet();
	if (evt_tcmetPhi_branch != 0) evt_tcmetPhi();
	if (evt_tcmetSig_branch != 0) evt_tcmetSig();
	if (evt_tcsumet_branch != 0) evt_tcsumet();
	if (evt_bsp4_branch != 0) evt_bsp4();
	if (l1_met_p4_branch != 0) l1_met_p4();
	if (l1_mht_p4_branch != 0) l1_mht_p4();
	if (els_mc_motherp4_branch != 0) els_mc_motherp4();
	if (els_mc_p4_branch != 0) els_mc_p4();
	if (jets_mc_gp_p4_branch != 0) jets_mc_gp_p4();
	if (jets_mc_motherp4_branch != 0) jets_mc_motherp4();
	if (jets_mc_p4_branch != 0) jets_mc_p4();
	if (mus_mc_motherp4_branch != 0) mus_mc_motherp4();
	if (mus_mc_p4_branch != 0) mus_mc_p4();
	if (pfjets_mc_gp_p4_branch != 0) pfjets_mc_gp_p4();
	if (pfjets_mc_motherp4_branch != 0) pfjets_mc_motherp4();
	if (pfjets_mc_p4_branch != 0) pfjets_mc_p4();
	if (photons_mc_motherp4_branch != 0) photons_mc_motherp4();
	if (photons_mc_p4_branch != 0) photons_mc_p4();
	if (trk_mcp4_branch != 0) trk_mcp4();
	if (els_conv_pos_p4_branch != 0) els_conv_pos_p4();
	if (els_inner_position_branch != 0) els_inner_position();
	if (els_outer_position_branch != 0) els_outer_position();
	if (els_p4_branch != 0) els_p4();
	if (els_p4In_branch != 0) els_p4In();
	if (els_p4Out_branch != 0) els_p4Out();
	if (els_trk_p4_branch != 0) els_trk_p4();
	if (els_vertex_p4_branch != 0) els_vertex_p4();
	if (genjets_p4_branch != 0) genjets_p4();
	if (genps_p4_branch != 0) genps_p4();
	if (genps_prod_vtx_branch != 0) genps_prod_vtx();
	if (gsftrks_inner_position_branch != 0) gsftrks_inner_position();
	if (gsftrks_outer_p4_branch != 0) gsftrks_outer_p4();
	if (gsftrks_outer_position_branch != 0) gsftrks_outer_position();
	if (gsftrks_p4_branch != 0) gsftrks_p4();
	if (gsftrks_vertex_p4_branch != 0) gsftrks_vertex_p4();
	if (hyp_ll_p4_branch != 0) hyp_ll_p4();
	if (hyp_ll_trk_p4_branch != 0) hyp_ll_trk_p4();
	if (hyp_lt_p4_branch != 0) hyp_lt_p4();
	if (hyp_lt_trk_p4_branch != 0) hyp_lt_trk_p4();
	if (hyp_p4_branch != 0) hyp_p4();
	if (hyp_FVFit_p4_branch != 0) hyp_FVFit_p4();
	if (hyp_FVFit_v4_branch != 0) hyp_FVFit_v4();
	if (hyp_ll_mc_p4_branch != 0) hyp_ll_mc_p4();
	if (hyp_lt_mc_p4_branch != 0) hyp_lt_mc_p4();
	if (jets_p4_branch != 0) jets_p4();
	if (jets_vertex_p4_branch != 0) jets_vertex_p4();
	if (jpts_p4_branch != 0) jpts_p4();
	if (l1_emiso_p4_branch != 0) l1_emiso_p4();
	if (l1_emnoiso_p4_branch != 0) l1_emnoiso_p4();
	if (l1_jetsc_p4_branch != 0) l1_jetsc_p4();
	if (l1_jetsf_p4_branch != 0) l1_jetsf_p4();
	if (l1_jetst_p4_branch != 0) l1_jetst_p4();
	if (l1_mus_p4_branch != 0) l1_mus_p4();
	if (mus_ecalpos_p4_branch != 0) mus_ecalpos_p4();
	if (mus_fitdefault_p4_branch != 0) mus_fitdefault_p4();
	if (mus_fitfirsthit_p4_branch != 0) mus_fitfirsthit_p4();
	if (mus_fitpicky_p4_branch != 0) mus_fitpicky_p4();
	if (mus_fittev_p4_branch != 0) mus_fittev_p4();
	if (mus_gfit_outerPos_p4_branch != 0) mus_gfit_outerPos_p4();
	if (mus_gfit_p4_branch != 0) mus_gfit_p4();
	if (mus_gfit_vertex_p4_branch != 0) mus_gfit_vertex_p4();
	if (mus_p4_branch != 0) mus_p4();
	if (mus_sta_p4_branch != 0) mus_sta_p4();
	if (mus_sta_vertex_p4_branch != 0) mus_sta_vertex_p4();
	if (mus_trk_p4_branch != 0) mus_trk_p4();
	if (mus_vertex_p4_branch != 0) mus_vertex_p4();
	if (els_pat_genMotherP4_branch != 0) els_pat_genMotherP4();
	if (els_pat_genP4_branch != 0) els_pat_genP4();
	if (els_pat_p4_branch != 0) els_pat_p4();
	if (jets_pat_genJet_p4_branch != 0) jets_pat_genJet_p4();
	if (jets_pat_genPartonMother_p4_branch != 0) jets_pat_genPartonMother_p4();
	if (jets_pat_genParton_p4_branch != 0) jets_pat_genParton_p4();
	if (jets_pat_jet_p4_branch != 0) jets_pat_jet_p4();
	if (jets_pat_jet_uncorp4_branch != 0) jets_pat_jet_uncorp4();
	if (mus_pat_genMotherP4_branch != 0) mus_pat_genMotherP4();
	if (mus_pat_genP4_branch != 0) mus_pat_genP4();
	if (mus_pat_p4_branch != 0) mus_pat_p4();
	if (pfels_p4_branch != 0) pfels_p4();
	if (pfels_posAtEcal_p4_branch != 0) pfels_posAtEcal_p4();
	if (pfjets_p4_branch != 0) pfjets_p4();
	if (pfmus_p4_branch != 0) pfmus_p4();
	if (pfmus_posAtEcal_p4_branch != 0) pfmus_posAtEcal_p4();
	if (photons_p4_branch != 0) photons_p4();
	if (scs_p4_branch != 0) scs_p4();
	if (scs_pos_p4_branch != 0) scs_pos_p4();
	if (scs_vtx_p4_branch != 0) scs_vtx_p4();
	if (svs_flight_branch != 0) svs_flight();
	if (svs_mc3_p4_branch != 0) svs_mc3_p4();
	if (svs_p4_branch != 0) svs_p4();
	if (svs_position_branch != 0) svs_position();
	if (svs_refitp4_branch != 0) svs_refitp4();
	if (trks_inner_position_branch != 0) trks_inner_position();
	if (trks_outer_p4_branch != 0) trks_outer_p4();
	if (trks_outer_position_branch != 0) trks_outer_position();
	if (trks_trk_p4_branch != 0) trks_trk_p4();
	if (trks_vertex_p4_branch != 0) trks_vertex_p4();
	if (trkjets_p4_branch != 0) trkjets_p4();
	if (vtxs_position_branch != 0) vtxs_position();
	if (genps_lepdaughter_p4_branch != 0) genps_lepdaughter_p4();
	if (hlt_trigObjs_p4_branch != 0) hlt_trigObjs_p4();
	if (hyp_jets_p4_branch != 0) hyp_jets_p4();
	if (hyp_other_jets_p4_branch != 0) hyp_other_jets_p4();
	if (jpts_combinedSecondaryVertexBJetTag_branch != 0) jpts_combinedSecondaryVertexBJetTag();
	if (jpts_combinedSecondaryVertexMVABJetTag_branch != 0) jpts_combinedSecondaryVertexMVABJetTag();
	if (jpts_jetBProbabilityBJetTag_branch != 0) jpts_jetBProbabilityBJetTag();
	if (jpts_jetProbabilityBJetTag_branch != 0) jpts_jetProbabilityBJetTag();
	if (jpts_simpleSecondaryVertexHighEffBJetTag_branch != 0) jpts_simpleSecondaryVertexHighEffBJetTag();
	if (jpts_simpleSecondaryVertexHighPurBJetTags_branch != 0) jpts_simpleSecondaryVertexHighPurBJetTags();
	if (jpts_softElectronByIP3dBJetTag_branch != 0) jpts_softElectronByIP3dBJetTag();
	if (jpts_softElectronByPtBJetTag_branch != 0) jpts_softElectronByPtBJetTag();
	if (jpts_softMuonBJetTag_branch != 0) jpts_softMuonBJetTag();
	if (jpts_softMuonByIP3dBJetTag_branch != 0) jpts_softMuonByIP3dBJetTag();
	if (jpts_softMuonByPtBJetTag_branch != 0) jpts_softMuonByPtBJetTag();
	if (jpts_trackCountingHighEffBJetTag_branch != 0) jpts_trackCountingHighEffBJetTag();
	if (jpts_trackCountingHighPurBJetTag_branch != 0) jpts_trackCountingHighPurBJetTag();
	if (jets_combinedSecondaryVertexBJetTag_branch != 0) jets_combinedSecondaryVertexBJetTag();
	if (jets_combinedSecondaryVertexMVABJetTag_branch != 0) jets_combinedSecondaryVertexMVABJetTag();
	if (jets_jetBProbabilityBJetTag_branch != 0) jets_jetBProbabilityBJetTag();
	if (jets_jetProbabilityBJetTag_branch != 0) jets_jetProbabilityBJetTag();
	if (jets_simpleSecondaryVertexHighEffBJetTag_branch != 0) jets_simpleSecondaryVertexHighEffBJetTag();
	if (jets_simpleSecondaryVertexHighPurBJetTags_branch != 0) jets_simpleSecondaryVertexHighPurBJetTags();
	if (jets_softElectronByIP3dBJetTag_branch != 0) jets_softElectronByIP3dBJetTag();
	if (jets_softElectronByPtBJetTag_branch != 0) jets_softElectronByPtBJetTag();
	if (jets_softMuonBJetTag_branch != 0) jets_softMuonBJetTag();
	if (jets_softMuonByIP3dBJetTag_branch != 0) jets_softMuonByIP3dBJetTag();
	if (jets_softMuonByPtBJetTag_branch != 0) jets_softMuonByPtBJetTag();
	if (jets_trackCountingHighEffBJetTag_branch != 0) jets_trackCountingHighEffBJetTag();
	if (jets_trackCountingHighPurBJetTag_branch != 0) jets_trackCountingHighPurBJetTag();
	if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) pfjets_combinedSecondaryVertexBJetTag();
	if (pfjets_combinedSecondaryVertexMVABJetTag_branch != 0) pfjets_combinedSecondaryVertexMVABJetTag();
	if (pfjets_jetBProbabilityBJetTag_branch != 0) pfjets_jetBProbabilityBJetTag();
	if (pfjets_jetProbabilityBJetTag_branch != 0) pfjets_jetProbabilityBJetTag();
	if (pfjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) pfjets_simpleSecondaryVertexHighEffBJetTag();
	if (pfjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) pfjets_simpleSecondaryVertexHighPurBJetTags();
	if (pfjets_softElectronByIP3dBJetTag_branch != 0) pfjets_softElectronByIP3dBJetTag();
	if (pfjets_softElectronByPtBJetTag_branch != 0) pfjets_softElectronByPtBJetTag();
	if (pfjets_softMuonBJetTag_branch != 0) pfjets_softMuonBJetTag();
	if (pfjets_softMuonByIP3dBJetTag_branch != 0) pfjets_softMuonByIP3dBJetTag();
	if (pfjets_softMuonByPtBJetTag_branch != 0) pfjets_softMuonByPtBJetTag();
	if (pfjets_trackCountingHighEffBJetTag_branch != 0) pfjets_trackCountingHighEffBJetTag();
	if (pfjets_trackCountingHighPurBJetTag_branch != 0) pfjets_trackCountingHighPurBJetTag();
	if (trkjets_combinedSecondaryVertexBJetTag_branch != 0) trkjets_combinedSecondaryVertexBJetTag();
	if (trkjets_combinedSecondaryVertexMVABJetTag_branch != 0) trkjets_combinedSecondaryVertexMVABJetTag();
	if (trkjets_jetBProbabilityBJetTag_branch != 0) trkjets_jetBProbabilityBJetTag();
	if (trkjets_jetProbabilityBJetTag_branch != 0) trkjets_jetProbabilityBJetTag();
	if (trkjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) trkjets_simpleSecondaryVertexHighEffBJetTag();
	if (trkjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) trkjets_simpleSecondaryVertexHighPurBJetTags();
	if (trkjets_softElectronByIP3dBJetTag_branch != 0) trkjets_softElectronByIP3dBJetTag();
	if (trkjets_softElectronByPtBJetTag_branch != 0) trkjets_softElectronByPtBJetTag();
	if (trkjets_softMuonBJetTag_branch != 0) trkjets_softMuonBJetTag();
	if (trkjets_softMuonByIP3dBJetTag_branch != 0) trkjets_softMuonByIP3dBJetTag();
	if (trkjets_softMuonByPtBJetTag_branch != 0) trkjets_softMuonByPtBJetTag();
	if (trkjets_trackCountingHighEffBJetTag_branch != 0) trkjets_trackCountingHighEffBJetTag();
	if (trkjets_trackCountingHighPurBJetTag_branch != 0) trkjets_trackCountingHighPurBJetTag();
	if (evt_bs_covMatrix_branch != 0) evt_bs_covMatrix();
	if (els_mc3dr_branch != 0) els_mc3dr();
	if (els_mcdr_branch != 0) els_mcdr();
	if (jets_mc3dr_branch != 0) jets_mc3dr();
	if (jets_mcdr_branch != 0) jets_mcdr();
	if (jets_mc_emEnergy_branch != 0) jets_mc_emEnergy();
	if (jets_mc_gpdr_branch != 0) jets_mc_gpdr();
	if (jets_mc_hadEnergy_branch != 0) jets_mc_hadEnergy();
	if (jets_mc_invEnergy_branch != 0) jets_mc_invEnergy();
	if (jets_mc_otherEnergy_branch != 0) jets_mc_otherEnergy();
	if (mus_mc3dr_branch != 0) mus_mc3dr();
	if (mus_mcdr_branch != 0) mus_mcdr();
	if (pfjets_mc3dr_branch != 0) pfjets_mc3dr();
	if (pfjets_mcdr_branch != 0) pfjets_mcdr();
	if (pfjets_mc_emEnergy_branch != 0) pfjets_mc_emEnergy();
	if (pfjets_mc_gpdr_branch != 0) pfjets_mc_gpdr();
	if (pfjets_mc_hadEnergy_branch != 0) pfjets_mc_hadEnergy();
	if (pfjets_mc_invEnergy_branch != 0) pfjets_mc_invEnergy();
	if (pfjets_mc_otherEnergy_branch != 0) pfjets_mc_otherEnergy();
	if (photons_mc3dr_branch != 0) photons_mc3dr();
	if (photons_mcdr_branch != 0) photons_mcdr();
	if (trk_mc3dr_branch != 0) trk_mc3dr();
	if (trk_mcdr_branch != 0) trk_mcdr();
	if (trks_conv_dcot_branch != 0) trks_conv_dcot();
	if (trks_conv_dist_branch != 0) trks_conv_dist();
	if (els_ecalJuraIso_branch != 0) els_ecalJuraIso();
	if (els_ecalJuraTowerIso_branch != 0) els_ecalJuraTowerIso();
	if (els_hcalConeIso_branch != 0) els_hcalConeIso();
	if (els_tkJuraIso_branch != 0) els_tkJuraIso();
	if (els_jetdr_branch != 0) els_jetdr();
	if (els_musdr_branch != 0) els_musdr();
	if (els_chi2_branch != 0) els_chi2();
	if (els_conv_dcot_branch != 0) els_conv_dcot();
	if (els_conv_dist_branch != 0) els_conv_dist();
	if (els_conv_radius_branch != 0) els_conv_radius();
	if (els_d0_branch != 0) els_d0();
	if (els_d0Err_branch != 0) els_d0Err();
	if (els_d0corr_branch != 0) els_d0corr();
	if (els_dEtaIn_branch != 0) els_dEtaIn();
	if (els_dEtaOut_branch != 0) els_dEtaOut();
	if (els_dPhiIn_branch != 0) els_dPhiIn();
	if (els_dPhiInPhiOut_branch != 0) els_dPhiInPhiOut();
	if (els_dPhiOut_branch != 0) els_dPhiOut();
	if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) els_deltaEtaEleClusterTrackAtCalo();
	if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) els_deltaPhiEleClusterTrackAtCalo();
	if (els_e1x5_branch != 0) els_e1x5();
	if (els_e2x5Max_branch != 0) els_e2x5Max();
	if (els_e3x3_branch != 0) els_e3x3();
	if (els_e5x5_branch != 0) els_e5x5();
	if (els_eMax_branch != 0) els_eMax();
	if (els_eOverPIn_branch != 0) els_eOverPIn();
	if (els_eOverPOut_branch != 0) els_eOverPOut();
	if (els_eSC_branch != 0) els_eSC();
	if (els_eSCPresh_branch != 0) els_eSCPresh();
	if (els_eSCRaw_branch != 0) els_eSCRaw();
	if (els_eSeed_branch != 0) els_eSeed();
	if (els_eSeedOverPIn_branch != 0) els_eSeedOverPIn();
	if (els_eSeedOverPOut_branch != 0) els_eSeedOverPOut();
	if (els_ecalEnergy_branch != 0) els_ecalEnergy();
	if (els_ecalEnergyError_branch != 0) els_ecalEnergyError();
	if (els_ecalIso_branch != 0) els_ecalIso();
	if (els_ecalIso04_branch != 0) els_ecalIso04();
	if (els_egamma_looseId_branch != 0) els_egamma_looseId();
	if (els_egamma_robustHighEnergy_branch != 0) els_egamma_robustHighEnergy();
	if (els_egamma_robustLooseId_branch != 0) els_egamma_robustLooseId();
	if (els_egamma_robustTightId_branch != 0) els_egamma_robustTightId();
	if (els_egamma_tightId_branch != 0) els_egamma_tightId();
	if (els_electronMomentumError_branch != 0) els_electronMomentumError();
	if (els_etaErr_branch != 0) els_etaErr();
	if (els_etaSC_branch != 0) els_etaSC();
	if (els_fbrem_branch != 0) els_fbrem();
	if (els_hOverE_branch != 0) els_hOverE();
	if (els_hcalDepth1OverEcal_branch != 0) els_hcalDepth1OverEcal();
	if (els_hcalDepth1TowerSumEt_branch != 0) els_hcalDepth1TowerSumEt();
	if (els_hcalDepth1TowerSumEt04_branch != 0) els_hcalDepth1TowerSumEt04();
	if (els_hcalDepth2OverEcal_branch != 0) els_hcalDepth2OverEcal();
	if (els_hcalDepth2TowerSumEt_branch != 0) els_hcalDepth2TowerSumEt();
	if (els_hcalDepth2TowerSumEt04_branch != 0) els_hcalDepth2TowerSumEt04();
	if (els_hcalIso_branch != 0) els_hcalIso();
	if (els_hcalIso04_branch != 0) els_hcalIso04();
	if (els_layer1_charge_branch != 0) els_layer1_charge();
	if (els_mva_branch != 0) els_mva();
	if (els_ndof_branch != 0) els_ndof();
	if (els_phiErr_branch != 0) els_phiErr();
	if (els_phiSC_branch != 0) els_phiSC();
	if (els_ptErr_branch != 0) els_ptErr();
	if (els_sigmaEtaEta_branch != 0) els_sigmaEtaEta();
	if (els_sigmaIEtaIEta_branch != 0) els_sigmaIEtaIEta();
	if (els_sigmaIEtaIEtaSC_branch != 0) els_sigmaIEtaIEtaSC();
	if (els_sigmaIPhiIPhi_branch != 0) els_sigmaIPhiIPhi();
	if (els_sigmaIPhiIPhiSC_branch != 0) els_sigmaIPhiIPhiSC();
	if (els_sigmaPhiPhi_branch != 0) els_sigmaPhiPhi();
	if (els_tkIso_branch != 0) els_tkIso();
	if (els_tkIso04_branch != 0) els_tkIso04();
	if (els_trackMomentumError_branch != 0) els_trackMomentumError();
	if (els_trkdr_branch != 0) els_trkdr();
	if (els_trkshFrac_branch != 0) els_trkshFrac();
	if (els_z0_branch != 0) els_z0();
	if (els_z0Err_branch != 0) els_z0Err();
	if (els_z0corr_branch != 0) els_z0corr();
	if (gsftrks_chi2_branch != 0) gsftrks_chi2();
	if (gsftrks_d0_branch != 0) gsftrks_d0();
	if (gsftrks_d0Err_branch != 0) gsftrks_d0Err();
	if (gsftrks_d0corr_branch != 0) gsftrks_d0corr();
	if (gsftrks_d0corrPhi_branch != 0) gsftrks_d0corrPhi();
	if (gsftrks_d0phiCov_branch != 0) gsftrks_d0phiCov();
	if (gsftrks_etaErr_branch != 0) gsftrks_etaErr();
	if (gsftrks_layer1_charge_branch != 0) gsftrks_layer1_charge();
	if (gsftrks_ndof_branch != 0) gsftrks_ndof();
	if (gsftrks_phiErr_branch != 0) gsftrks_phiErr();
	if (gsftrks_ptErr_branch != 0) gsftrks_ptErr();
	if (gsftrks_z0_branch != 0) gsftrks_z0();
	if (gsftrks_z0Err_branch != 0) gsftrks_z0Err();
	if (gsftrks_z0corr_branch != 0) gsftrks_z0corr();
	if (hyp_Ht_branch != 0) hyp_Ht();
	if (hyp_dPhi_nJet_metMuonJESCorr_branch != 0) hyp_dPhi_nJet_metMuonJESCorr();
	if (hyp_dPhi_nJet_muCorrMet_branch != 0) hyp_dPhi_nJet_muCorrMet();
	if (hyp_dPhi_nJet_tcMet_branch != 0) hyp_dPhi_nJet_tcMet();
	if (hyp_dPhi_nJet_unCorrMet_branch != 0) hyp_dPhi_nJet_unCorrMet();
	if (hyp_ll_chi2_branch != 0) hyp_ll_chi2();
	if (hyp_ll_d0_branch != 0) hyp_ll_d0();
	if (hyp_ll_d0Err_branch != 0) hyp_ll_d0Err();
	if (hyp_ll_d0corr_branch != 0) hyp_ll_d0corr();
	if (hyp_ll_dPhi_metMuonJESCorr_branch != 0) hyp_ll_dPhi_metMuonJESCorr();
	if (hyp_ll_dPhi_muCorrMet_branch != 0) hyp_ll_dPhi_muCorrMet();
	if (hyp_ll_dPhi_tcMet_branch != 0) hyp_ll_dPhi_tcMet();
	if (hyp_ll_dPhi_unCorrMet_branch != 0) hyp_ll_dPhi_unCorrMet();
	if (hyp_ll_etaErr_branch != 0) hyp_ll_etaErr();
	if (hyp_ll_ndof_branch != 0) hyp_ll_ndof();
	if (hyp_ll_phiErr_branch != 0) hyp_ll_phiErr();
	if (hyp_ll_ptErr_branch != 0) hyp_ll_ptErr();
	if (hyp_ll_z0_branch != 0) hyp_ll_z0();
	if (hyp_ll_z0Err_branch != 0) hyp_ll_z0Err();
	if (hyp_ll_z0corr_branch != 0) hyp_ll_z0corr();
	if (hyp_lt_chi2_branch != 0) hyp_lt_chi2();
	if (hyp_lt_d0_branch != 0) hyp_lt_d0();
	if (hyp_lt_d0Err_branch != 0) hyp_lt_d0Err();
	if (hyp_lt_d0corr_branch != 0) hyp_lt_d0corr();
	if (hyp_lt_dPhi_metMuonJESCorr_branch != 0) hyp_lt_dPhi_metMuonJESCorr();
	if (hyp_lt_dPhi_muCorrMet_branch != 0) hyp_lt_dPhi_muCorrMet();
	if (hyp_lt_dPhi_tcMet_branch != 0) hyp_lt_dPhi_tcMet();
	if (hyp_lt_dPhi_unCorrMet_branch != 0) hyp_lt_dPhi_unCorrMet();
	if (hyp_lt_etaErr_branch != 0) hyp_lt_etaErr();
	if (hyp_lt_ndof_branch != 0) hyp_lt_ndof();
	if (hyp_lt_phiErr_branch != 0) hyp_lt_phiErr();
	if (hyp_lt_ptErr_branch != 0) hyp_lt_ptErr();
	if (hyp_lt_z0_branch != 0) hyp_lt_z0();
	if (hyp_lt_z0Err_branch != 0) hyp_lt_z0Err();
	if (hyp_lt_z0corr_branch != 0) hyp_lt_z0corr();
	if (hyp_mt2_metMuonJESCorr_branch != 0) hyp_mt2_metMuonJESCorr();
	if (hyp_mt2_muCorrMet_branch != 0) hyp_mt2_muCorrMet();
	if (hyp_mt2_tcMet_branch != 0) hyp_mt2_tcMet();
	if (hyp_sumJetPt_branch != 0) hyp_sumJetPt();
	if (hyp_FVFit_chi2ndf_branch != 0) hyp_FVFit_chi2ndf();
	if (hyp_FVFit_prob_branch != 0) hyp_FVFit_prob();
	if (hyp_FVFit_v4cxx_branch != 0) hyp_FVFit_v4cxx();
	if (hyp_FVFit_v4cxy_branch != 0) hyp_FVFit_v4cxy();
	if (hyp_FVFit_v4cyy_branch != 0) hyp_FVFit_v4cyy();
	if (hyp_FVFit_v4czx_branch != 0) hyp_FVFit_v4czx();
	if (hyp_FVFit_v4czy_branch != 0) hyp_FVFit_v4czy();
	if (hyp_FVFit_v4czz_branch != 0) hyp_FVFit_v4czz();
	if (hyp_ll_ecaliso_branch != 0) hyp_ll_ecaliso();
	if (hyp_ll_trkiso_branch != 0) hyp_ll_trkiso();
	if (hyp_lt_ecaliso_branch != 0) hyp_lt_ecaliso();
	if (hyp_lt_trkiso_branch != 0) hyp_lt_trkiso();
	if (jets_approximatefHPD_branch != 0) jets_approximatefHPD();
	if (jets_approximatefRBX_branch != 0) jets_approximatefRBX();
	if (jets_cor_branch != 0) jets_cor();
	if (jets_emFrac_branch != 0) jets_emFrac();
	if (jets_fHPD_branch != 0) jets_fHPD();
	if (jets_fRBX_branch != 0) jets_fRBX();
	if (jets_fSubDetector1_branch != 0) jets_fSubDetector1();
	if (jets_fSubDetector2_branch != 0) jets_fSubDetector2();
	if (jets_fSubDetector3_branch != 0) jets_fSubDetector3();
	if (jets_fSubDetector4_branch != 0) jets_fSubDetector4();
	if (jets_hitsInN90_branch != 0) jets_hitsInN90();
	if (jets_n90Hits_branch != 0) jets_n90Hits();
	if (jets_nECALTowers_branch != 0) jets_nECALTowers();
	if (jets_nHCALTowers_branch != 0) jets_nHCALTowers();
	if (jets_restrictedEMF_branch != 0) jets_restrictedEMF();
	if (jpts_cor_branch != 0) jpts_cor();
	if (jpts_emFrac_branch != 0) jpts_emFrac();
	if (mus_met_deltax_branch != 0) mus_met_deltax();
	if (mus_met_deltay_branch != 0) mus_met_deltay();
	if (mus_eledr_branch != 0) mus_eledr();
	if (mus_jetdr_branch != 0) mus_jetdr();
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
	if (mus_gfit_d0_branch != 0) mus_gfit_d0();
	if (mus_gfit_d0Err_branch != 0) mus_gfit_d0Err();
	if (mus_gfit_d0corr_branch != 0) mus_gfit_d0corr();
	if (mus_gfit_ndof_branch != 0) mus_gfit_ndof();
	if (mus_gfit_qoverp_branch != 0) mus_gfit_qoverp();
	if (mus_gfit_qoverpError_branch != 0) mus_gfit_qoverpError();
	if (mus_gfit_z0_branch != 0) mus_gfit_z0();
	if (mus_gfit_z0Err_branch != 0) mus_gfit_z0Err();
	if (mus_gfit_z0corr_branch != 0) mus_gfit_z0corr();
	if (mus_iso03_emEt_branch != 0) mus_iso03_emEt();
	if (mus_iso03_hadEt_branch != 0) mus_iso03_hadEt();
	if (mus_iso03_hoEt_branch != 0) mus_iso03_hoEt();
	if (mus_iso03_sumPt_branch != 0) mus_iso03_sumPt();
	if (mus_iso05_emEt_branch != 0) mus_iso05_emEt();
	if (mus_iso05_hadEt_branch != 0) mus_iso05_hadEt();
	if (mus_iso05_hoEt_branch != 0) mus_iso05_hoEt();
	if (mus_iso05_sumPt_branch != 0) mus_iso05_sumPt();
	if (mus_iso_ecalvetoDep_branch != 0) mus_iso_ecalvetoDep();
	if (mus_iso_hcalvetoDep_branch != 0) mus_iso_hcalvetoDep();
	if (mus_iso_hovetoDep_branch != 0) mus_iso_hovetoDep();
	if (mus_iso_trckvetoDep_branch != 0) mus_iso_trckvetoDep();
	if (mus_ndof_branch != 0) mus_ndof();
	if (mus_phiErr_branch != 0) mus_phiErr();
	if (mus_ptErr_branch != 0) mus_ptErr();
	if (mus_qoverp_branch != 0) mus_qoverp();
	if (mus_qoverpError_branch != 0) mus_qoverpError();
	if (mus_sta_chi2_branch != 0) mus_sta_chi2();
	if (mus_sta_d0_branch != 0) mus_sta_d0();
	if (mus_sta_d0Err_branch != 0) mus_sta_d0Err();
	if (mus_sta_d0corr_branch != 0) mus_sta_d0corr();
	if (mus_sta_ndof_branch != 0) mus_sta_ndof();
	if (mus_sta_qoverp_branch != 0) mus_sta_qoverp();
	if (mus_sta_qoverpError_branch != 0) mus_sta_qoverpError();
	if (mus_sta_z0_branch != 0) mus_sta_z0();
	if (mus_sta_z0Err_branch != 0) mus_sta_z0Err();
	if (mus_sta_z0corr_branch != 0) mus_sta_z0corr();
	if (mus_timeAtIpInOut_branch != 0) mus_timeAtIpInOut();
	if (mus_timeAtIpInOutErr_branch != 0) mus_timeAtIpInOutErr();
	if (mus_timeAtIpOutIn_branch != 0) mus_timeAtIpOutIn();
	if (mus_timeAtIpOutInErr_branch != 0) mus_timeAtIpOutInErr();
	if (mus_vertexphi_branch != 0) mus_vertexphi();
	if (mus_z0_branch != 0) mus_z0();
	if (mus_z0Err_branch != 0) mus_z0Err();
	if (mus_z0corr_branch != 0) mus_z0corr();
	if (els_pat_caloIso_branch != 0) els_pat_caloIso();
	if (els_pat_ecalIso_branch != 0) els_pat_ecalIso();
	if (els_pat_hcalIso_branch != 0) els_pat_hcalIso();
	if (els_pat_looseId_branch != 0) els_pat_looseId();
	if (els_pat_robustHighEnergy_branch != 0) els_pat_robustHighEnergy();
	if (els_pat_robustLooseId_branch != 0) els_pat_robustLooseId();
	if (els_pat_robustTightId_branch != 0) els_pat_robustTightId();
	if (els_pat_scE1x5_branch != 0) els_pat_scE1x5();
	if (els_pat_scE2x5Max_branch != 0) els_pat_scE2x5Max();
	if (els_pat_scE5x5_branch != 0) els_pat_scE5x5();
	if (els_pat_sigmaEtaEta_branch != 0) els_pat_sigmaEtaEta();
	if (els_pat_sigmaIEtaIEta_branch != 0) els_pat_sigmaIEtaIEta();
	if (els_pat_tightId_branch != 0) els_pat_tightId();
	if (els_pat_trackIso_branch != 0) els_pat_trackIso();
	if (jets_pat_combinedSecondaryVertexBJetTag_branch != 0) jets_pat_combinedSecondaryVertexBJetTag();
	if (jets_pat_combinedSecondaryVertexMVABJetTag_branch != 0) jets_pat_combinedSecondaryVertexMVABJetTag();
	if (jets_pat_coneIsolationTauJetTag_branch != 0) jets_pat_coneIsolationTauJetTag();
	if (jets_pat_impactParameterMVABJetTag_branch != 0) jets_pat_impactParameterMVABJetTag();
	if (jets_pat_jetBProbabilityBJetTag_branch != 0) jets_pat_jetBProbabilityBJetTag();
	if (jets_pat_jetCharge_branch != 0) jets_pat_jetCharge();
	if (jets_pat_jetProbabilityBJetTag_branch != 0) jets_pat_jetProbabilityBJetTag();
	if (jets_pat_noCorrF_branch != 0) jets_pat_noCorrF();
	if (jets_pat_simpleSecondaryVertexHighEffBJetTag_branch != 0) jets_pat_simpleSecondaryVertexHighEffBJetTag();
	if (jets_pat_simpleSecondaryVertexHighPurBJetTag_branch != 0) jets_pat_simpleSecondaryVertexHighPurBJetTag();
	if (jets_pat_softElectronByIP3dBJetTag_branch != 0) jets_pat_softElectronByIP3dBJetTag();
	if (jets_pat_softElectronByPtBJetTag_branch != 0) jets_pat_softElectronByPtBJetTag();
	if (jets_pat_softMuonBJetTag_branch != 0) jets_pat_softMuonBJetTag();
	if (jets_pat_softMuonByIP3dBJetTag_branch != 0) jets_pat_softMuonByIP3dBJetTag();
	if (jets_pat_softMuonByPtBJetTag_branch != 0) jets_pat_softMuonByPtBJetTag();
	if (jets_pat_trackCountingHighEffBJetTag_branch != 0) jets_pat_trackCountingHighEffBJetTag();
	if (jets_pat_trackCountingHighPurBJetTag_branch != 0) jets_pat_trackCountingHighPurBJetTag();
	if (mus_pat_caloIso_branch != 0) mus_pat_caloIso();
	if (mus_pat_calovetoDep_branch != 0) mus_pat_calovetoDep();
	if (mus_pat_ecalIso_branch != 0) mus_pat_ecalIso();
	if (mus_pat_ecalvetoDep_branch != 0) mus_pat_ecalvetoDep();
	if (mus_pat_hcalIso_branch != 0) mus_pat_hcalIso();
	if (mus_pat_hcalvetoDep_branch != 0) mus_pat_hcalvetoDep();
	if (mus_pat_trackIso_branch != 0) mus_pat_trackIso();
	if (mus_pat_trckvetoDep_branch != 0) mus_pat_trckvetoDep();
	if (pfels_deltaP_branch != 0) pfels_deltaP();
	if (pfels_ecalE_branch != 0) pfels_ecalE();
	if (pfels_hcalE_branch != 0) pfels_hcalE();
	if (pfels_isoChargedHadrons_branch != 0) pfels_isoChargedHadrons();
	if (pfels_isoNeutralHadrons_branch != 0) pfels_isoNeutralHadrons();
	if (pfels_isoPhotons_branch != 0) pfels_isoPhotons();
	if (pfels_mva_emu_branch != 0) pfels_mva_emu();
	if (pfels_mva_epi_branch != 0) pfels_mva_epi();
	if (pfels_mva_nothing_gamma_branch != 0) pfels_mva_nothing_gamma();
	if (pfels_mva_nothing_nh_branch != 0) pfels_mva_nothing_nh();
	if (pfels_mva_pimu_branch != 0) pfels_mva_pimu();
	if (pfels_pS1E_branch != 0) pfels_pS1E();
	if (pfels_pS2E_branch != 0) pfels_pS2E();
	if (pfels_rawEcalE_branch != 0) pfels_rawEcalE();
	if (pfels_rawHcalE_branch != 0) pfels_rawHcalE();
	if (pfjets_chargedEmE_branch != 0) pfjets_chargedEmE();
	if (pfjets_chargedHadronE_branch != 0) pfjets_chargedHadronE();
	if (pfjets_cor_branch != 0) pfjets_cor();
	if (pfjets_neutralEmE_branch != 0) pfjets_neutralEmE();
	if (pfjets_neutralHadronE_branch != 0) pfjets_neutralHadronE();
	if (pfmus_deltaP_branch != 0) pfmus_deltaP();
	if (pfmus_ecalE_branch != 0) pfmus_ecalE();
	if (pfmus_hcalE_branch != 0) pfmus_hcalE();
	if (pfmus_isoChargedHadrons_branch != 0) pfmus_isoChargedHadrons();
	if (pfmus_isoNeutralHadrons_branch != 0) pfmus_isoNeutralHadrons();
	if (pfmus_isoPhotons_branch != 0) pfmus_isoPhotons();
	if (pfmus_mva_emu_branch != 0) pfmus_mva_emu();
	if (pfmus_mva_epi_branch != 0) pfmus_mva_epi();
	if (pfmus_mva_nothing_gamma_branch != 0) pfmus_mva_nothing_gamma();
	if (pfmus_mva_nothing_nh_branch != 0) pfmus_mva_nothing_nh();
	if (pfmus_mva_pimu_branch != 0) pfmus_mva_pimu();
	if (pfmus_pS1E_branch != 0) pfmus_pS1E();
	if (pfmus_pS2E_branch != 0) pfmus_pS2E();
	if (pfmus_rawEcalE_branch != 0) pfmus_rawEcalE();
	if (pfmus_rawHcalE_branch != 0) pfmus_rawHcalE();
	if (photons_e1x5_branch != 0) photons_e1x5();
	if (photons_e2x5Max_branch != 0) photons_e2x5Max();
	if (photons_e3x3_branch != 0) photons_e3x3();
	if (photons_e5x5_branch != 0) photons_e5x5();
	if (photons_ecalIso03_branch != 0) photons_ecalIso03();
	if (photons_ecalIso04_branch != 0) photons_ecalIso04();
	if (photons_hOverE_branch != 0) photons_hOverE();
	if (photons_hcalIso03_branch != 0) photons_hcalIso03();
	if (photons_hcalIso04_branch != 0) photons_hcalIso04();
	if (photons_ntkIsoHollow03_branch != 0) photons_ntkIsoHollow03();
	if (photons_ntkIsoHollow04_branch != 0) photons_ntkIsoHollow04();
	if (photons_ntkIsoSolid03_branch != 0) photons_ntkIsoSolid03();
	if (photons_ntkIsoSolid04_branch != 0) photons_ntkIsoSolid04();
	if (photons_sigmaEtaEta_branch != 0) photons_sigmaEtaEta();
	if (photons_sigmaIEtaIEta_branch != 0) photons_sigmaIEtaIEta();
	if (photons_swissSeed_branch != 0) photons_swissSeed();
	if (photons_tkIsoHollow03_branch != 0) photons_tkIsoHollow03();
	if (photons_tkIsoHollow04_branch != 0) photons_tkIsoHollow04();
	if (photons_tkIsoSolid03_branch != 0) photons_tkIsoSolid03();
	if (photons_tkIsoSolid04_branch != 0) photons_tkIsoSolid04();
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
	if (scs_sigmaIEtaIEtaSC_branch != 0) scs_sigmaIEtaIEtaSC();
	if (scs_sigmaIEtaIPhi_branch != 0) scs_sigmaIEtaIPhi();
	if (scs_sigmaIEtaIPhiSC_branch != 0) scs_sigmaIEtaIPhiSC();
	if (scs_sigmaIPhiIPhi_branch != 0) scs_sigmaIPhiIPhi();
	if (scs_sigmaIPhiIPhiSC_branch != 0) scs_sigmaIPhiIPhiSC();
	if (scs_sigmaPhiPhi_branch != 0) scs_sigmaPhiPhi();
	if (scs_timeSeed_branch != 0) scs_timeSeed();
	if (svs_anglePV_branch != 0) svs_anglePV();
	if (svs_chi2_branch != 0) svs_chi2();
	if (svs_dist3Dsig_branch != 0) svs_dist3Dsig();
	if (svs_dist3Dval_branch != 0) svs_dist3Dval();
	if (svs_distXYsig_branch != 0) svs_distXYsig();
	if (svs_distXYval_branch != 0) svs_distXYval();
	if (svs_ndof_branch != 0) svs_ndof();
	if (svs_prob_branch != 0) svs_prob();
	if (svs_xError_branch != 0) svs_xError();
	if (svs_yError_branch != 0) svs_yError();
	if (svs_zError_branch != 0) svs_zError();
	if (mus_tcmet_deltax_branch != 0) mus_tcmet_deltax();
	if (mus_tcmet_deltay_branch != 0) mus_tcmet_deltay();
	if (trks_chi2_branch != 0) trks_chi2();
	if (trks_d0_branch != 0) trks_d0();
	if (trks_d0Err_branch != 0) trks_d0Err();
	if (trks_d0corr_branch != 0) trks_d0corr();
	if (trks_d0corrPhi_branch != 0) trks_d0corrPhi();
	if (trks_d0phiCov_branch != 0) trks_d0phiCov();
	if (trks_etaErr_branch != 0) trks_etaErr();
	if (trks_layer1_charge_branch != 0) trks_layer1_charge();
	if (trks_ndof_branch != 0) trks_ndof();
	if (trks_phiErr_branch != 0) trks_phiErr();
	if (trks_ptErr_branch != 0) trks_ptErr();
	if (trks_z0_branch != 0) trks_z0();
	if (trks_z0Err_branch != 0) trks_z0Err();
	if (trks_z0corr_branch != 0) trks_z0corr();
	if (trkjets_cor_branch != 0) trkjets_cor();
	if (trks_d0Errvtx_branch != 0) trks_d0Errvtx();
	if (trks_d0vtx_branch != 0) trks_d0vtx();
	if (vtxs_chi2_branch != 0) vtxs_chi2();
	if (vtxs_ndof_branch != 0) vtxs_ndof();
	if (vtxs_sumpt_branch != 0) vtxs_sumpt();
	if (vtxs_xError_branch != 0) vtxs_xError();
	if (vtxs_yError_branch != 0) vtxs_yError();
	if (vtxs_zError_branch != 0) vtxs_zError();
	if (vtxs_covMatrix_branch != 0) vtxs_covMatrix();
	if (evt_cscLooseHaloId_branch != 0) evt_cscLooseHaloId();
	if (evt_cscTightHaloId_branch != 0) evt_cscTightHaloId();
	if (evt_ecalLooseHaloId_branch != 0) evt_ecalLooseHaloId();
	if (evt_ecalTightHaloId_branch != 0) evt_ecalTightHaloId();
	if (evt_extremeTightHaloId_branch != 0) evt_extremeTightHaloId();
	if (evt_globalLooseHaloId_branch != 0) evt_globalLooseHaloId();
	if (evt_globalTightHaloId_branch != 0) evt_globalTightHaloId();
	if (evt_hcalLooseHaloId_branch != 0) evt_hcalLooseHaloId();
	if (evt_hcalTightHaloId_branch != 0) evt_hcalTightHaloId();
	if (evt_looseHaloId_branch != 0) evt_looseHaloId();
	if (evt_nHaloLikeTracks_branch != 0) evt_nHaloLikeTracks();
	if (evt_nHaloTriggerCandidates_branch != 0) evt_nHaloTriggerCandidates();
	if (evt_tightHaloId_branch != 0) evt_tightHaloId();
	if (evt_bsType_branch != 0) evt_bsType();
	if (evt_bunchCrossing_branch != 0) evt_bunchCrossing();
	if (evt_experimentType_branch != 0) evt_experimentType();
	if (evt_isRealData_branch != 0) evt_isRealData();
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
	if (pdfinfo_id1_branch != 0) pdfinfo_id1();
	if (pdfinfo_id2_branch != 0) pdfinfo_id2();
	if (evt_ecaliPhiSuspects_branch != 0) evt_ecaliPhiSuspects();
	if (evt_globaliPhiSuspects_branch != 0) evt_globaliPhiSuspects();
	if (evt_hcaliPhiSuspects_branch != 0) evt_hcaliPhiSuspects();
	if (els_mc3_id_branch != 0) els_mc3_id();
	if (els_mc3idx_branch != 0) els_mc3idx();
	if (els_mc3_motherid_branch != 0) els_mc3_motherid();
	if (els_mc3_motheridx_branch != 0) els_mc3_motheridx();
	if (els_mc_id_branch != 0) els_mc_id();
	if (els_mcidx_branch != 0) els_mcidx();
	if (els_mc_motherid_branch != 0) els_mc_motherid();
	if (jets_mc3_id_branch != 0) jets_mc3_id();
	if (jets_mc3idx_branch != 0) jets_mc3idx();
	if (jets_mc_gpidx_branch != 0) jets_mc_gpidx();
	if (jets_mc_id_branch != 0) jets_mc_id();
	if (jets_mcidx_branch != 0) jets_mcidx();
	if (jets_mc_motherid_branch != 0) jets_mc_motherid();
	if (mus_mc3_id_branch != 0) mus_mc3_id();
	if (mus_mc3idx_branch != 0) mus_mc3idx();
	if (mus_mc3_motherid_branch != 0) mus_mc3_motherid();
	if (mus_mc3_motheridx_branch != 0) mus_mc3_motheridx();
	if (mus_mc_id_branch != 0) mus_mc_id();
	if (mus_mcidx_branch != 0) mus_mcidx();
	if (mus_mc_motherid_branch != 0) mus_mc_motherid();
	if (pfjets_mc3_id_branch != 0) pfjets_mc3_id();
	if (pfjets_mc3idx_branch != 0) pfjets_mc3idx();
	if (pfjets_mc_gpidx_branch != 0) pfjets_mc_gpidx();
	if (pfjets_mc_id_branch != 0) pfjets_mc_id();
	if (pfjets_mcidx_branch != 0) pfjets_mcidx();
	if (pfjets_mc_motherid_branch != 0) pfjets_mc_motherid();
	if (photons_mc3_id_branch != 0) photons_mc3_id();
	if (photons_mc3idx_branch != 0) photons_mc3idx();
	if (photons_mc3_motherid_branch != 0) photons_mc3_motherid();
	if (photons_mc3_motheridx_branch != 0) photons_mc3_motheridx();
	if (photons_mc_id_branch != 0) photons_mc_id();
	if (photons_mcidx_branch != 0) photons_mcidx();
	if (photons_mc_motherid_branch != 0) photons_mc_motherid();
	if (trk_mc3_id_branch != 0) trk_mc3_id();
	if (trk_mc3idx_branch != 0) trk_mc3idx();
	if (trk_mc3_motherid_branch != 0) trk_mc3_motherid();
	if (trk_mc3_motheridx_branch != 0) trk_mc3_motheridx();
	if (trk_mc_id_branch != 0) trk_mc_id();
	if (trk_mcidx_branch != 0) trk_mcidx();
	if (trk_mc_motherid_branch != 0) trk_mc_motherid();
	if (trks_conv_tkidx_branch != 0) trks_conv_tkidx();
	if (els_exp_innerlayers39X_branch != 0) els_exp_innerlayers39X();
	if (els_closestJet_branch != 0) els_closestJet();
	if (els_closestMuon_branch != 0) els_closestMuon();
	if (els_pfelsidx_branch != 0) els_pfelsidx();
	if (els_category_branch != 0) els_category();
	if (els_charge_branch != 0) els_charge();
	if (els_class_branch != 0) els_class();
	if (els_conv_tkidx_branch != 0) els_conv_tkidx();
	if (els_exp_innerlayers_branch != 0) els_exp_innerlayers();
	if (els_exp_outerlayers_branch != 0) els_exp_outerlayers();
	if (els_fiduciality_branch != 0) els_fiduciality();
	if (els_gsftrkidx_branch != 0) els_gsftrkidx();
	if (els_layer1_det_branch != 0) els_layer1_det();
	if (els_layer1_layer_branch != 0) els_layer1_layer();
	if (els_layer1_sizerphi_branch != 0) els_layer1_sizerphi();
	if (els_layer1_sizerz_branch != 0) els_layer1_sizerz();
	if (els_lostHits_branch != 0) els_lostHits();
	if (els_lost_pixelhits_branch != 0) els_lost_pixelhits();
	if (els_nSeed_branch != 0) els_nSeed();
	if (els_sccharge_branch != 0) els_sccharge();
	if (els_scindex_branch != 0) els_scindex();
	if (els_trk_charge_branch != 0) els_trk_charge();
	if (els_trkidx_branch != 0) els_trkidx();
	if (els_type_branch != 0) els_type();
	if (els_validHits_branch != 0) els_validHits();
	if (els_valid_pixelhits_branch != 0) els_valid_pixelhits();
	if (genps_id_branch != 0) genps_id();
	if (genps_id_mother_branch != 0) genps_id_mother();
	if (genps_status_branch != 0) genps_status();
	if (gsftrks_charge_branch != 0) gsftrks_charge();
	if (gsftrks_exp_innerlayers_branch != 0) gsftrks_exp_innerlayers();
	if (gsftrks_exp_outerlayers_branch != 0) gsftrks_exp_outerlayers();
	if (gsftrks_layer1_det_branch != 0) gsftrks_layer1_det();
	if (gsftrks_layer1_layer_branch != 0) gsftrks_layer1_layer();
	if (gsftrks_layer1_sizerphi_branch != 0) gsftrks_layer1_sizerphi();
	if (gsftrks_layer1_sizerz_branch != 0) gsftrks_layer1_sizerz();
	if (gsftrks_lostHits_branch != 0) gsftrks_lostHits();
	if (gsftrks_lost_pixelhits_branch != 0) gsftrks_lost_pixelhits();
	if (gsftrks_nlayers_branch != 0) gsftrks_nlayers();
	if (gsftrks_nlayers3D_branch != 0) gsftrks_nlayers3D();
	if (gsftrks_nlayersLost_branch != 0) gsftrks_nlayersLost();
	if (gsftrks_validHits_branch != 0) gsftrks_validHits();
	if (gsftrks_valid_pixelhits_branch != 0) gsftrks_valid_pixelhits();
	if (hyp_ll_charge_branch != 0) hyp_ll_charge();
	if (hyp_ll_id_branch != 0) hyp_ll_id();
	if (hyp_ll_index_branch != 0) hyp_ll_index();
	if (hyp_ll_lostHits_branch != 0) hyp_ll_lostHits();
	if (hyp_ll_validHits_branch != 0) hyp_ll_validHits();
	if (hyp_lt_charge_branch != 0) hyp_lt_charge();
	if (hyp_lt_id_branch != 0) hyp_lt_id();
	if (hyp_lt_index_branch != 0) hyp_lt_index();
	if (hyp_lt_lostHits_branch != 0) hyp_lt_lostHits();
	if (hyp_lt_validHits_branch != 0) hyp_lt_validHits();
	if (hyp_njets_branch != 0) hyp_njets();
	if (hyp_nojets_branch != 0) hyp_nojets();
	if (hyp_type_branch != 0) hyp_type();
	if (hyp_FVFit_ndf_branch != 0) hyp_FVFit_ndf();
	if (hyp_FVFit_status_branch != 0) hyp_FVFit_status();
	if (hyp_ll_mc_id_branch != 0) hyp_ll_mc_id();
	if (hyp_ll_mc_motherid_branch != 0) hyp_ll_mc_motherid();
	if (hyp_lt_mc_id_branch != 0) hyp_lt_mc_id();
	if (hyp_lt_mc_motherid_branch != 0) hyp_lt_mc_motherid();
	if (hyp_quadlep_first_type_branch != 0) hyp_quadlep_first_type();
	if (hyp_quadlep_fourth_type_branch != 0) hyp_quadlep_fourth_type();
	if (hyp_quadlep_second_type_branch != 0) hyp_quadlep_second_type();
	if (hyp_quadlep_third_type_branch != 0) hyp_quadlep_third_type();
	if (hyp_trilep_first_type_branch != 0) hyp_trilep_first_type();
	if (hyp_trilep_second_type_branch != 0) hyp_trilep_second_type();
	if (hyp_trilep_third_type_branch != 0) hyp_trilep_third_type();
	if (jets_closestElectron_branch != 0) jets_closestElectron();
	if (jets_closestMuon_branch != 0) jets_closestMuon();
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
	if (mus_met_flag_branch != 0) mus_met_flag();
	if (mus_closestEle_branch != 0) mus_closestEle();
	if (mus_closestJet_branch != 0) mus_closestJet();
	if (mus_pfmusidx_branch != 0) mus_pfmusidx();
	if (mus_charge_branch != 0) mus_charge();
	if (mus_chi2LocalMomentum_branch != 0) mus_chi2LocalMomentum();
	if (mus_chi2LocalPosition_branch != 0) mus_chi2LocalPosition();
	if (mus_gfit_validHits_branch != 0) mus_gfit_validHits();
	if (mus_gfit_validSTAHits_branch != 0) mus_gfit_validSTAHits();
	if (mus_gfit_validSiHits_branch != 0) mus_gfit_validSiHits();
	if (mus_glbKink_branch != 0) mus_glbKink();
	if (mus_glbTrackProbability_branch != 0) mus_glbTrackProbability();
	if (mus_globalDeltaEtaPhi_branch != 0) mus_globalDeltaEtaPhi();
	if (mus_goodmask_branch != 0) mus_goodmask();
	if (mus_iso03_ntrk_branch != 0) mus_iso03_ntrk();
	if (mus_iso05_ntrk_branch != 0) mus_iso05_ntrk();
	if (mus_localDistance_branch != 0) mus_localDistance();
	if (mus_lostHits_branch != 0) mus_lostHits();
	if (mus_nOverlaps_branch != 0) mus_nOverlaps();
	if (mus_nmatches_branch != 0) mus_nmatches();
	if (mus_overlap0_branch != 0) mus_overlap0();
	if (mus_overlap1_branch != 0) mus_overlap1();
	if (mus_pid_TM2DCompatibilityLoose_branch != 0) mus_pid_TM2DCompatibilityLoose();
	if (mus_pid_TM2DCompatibilityTight_branch != 0) mus_pid_TM2DCompatibilityTight();
	if (mus_pid_TMLastStationLoose_branch != 0) mus_pid_TMLastStationLoose();
	if (mus_pid_TMLastStationTight_branch != 0) mus_pid_TMLastStationTight();
	if (mus_staRelChi2_branch != 0) mus_staRelChi2();
	if (mus_sta_validHits_branch != 0) mus_sta_validHits();
	if (mus_timeDirection_branch != 0) mus_timeDirection();
	if (mus_timeNumStationsUsed_branch != 0) mus_timeNumStationsUsed();
	if (mus_trkKink_branch != 0) mus_trkKink();
	if (mus_trkRelChi2_branch != 0) mus_trkRelChi2();
	if (mus_trk_charge_branch != 0) mus_trk_charge();
	if (mus_trkidx_branch != 0) mus_trkidx();
	if (mus_type_branch != 0) mus_type();
	if (mus_validHits_branch != 0) mus_validHits();
	if (els_pat_genID_branch != 0) els_pat_genID();
	if (els_pat_genMotherID_branch != 0) els_pat_genMotherID();
	if (jets_pat_genPartonMother_id_branch != 0) jets_pat_genPartonMother_id();
	if (jets_pat_genParton_id_branch != 0) jets_pat_genParton_id();
	if (jets_pat_jetIDLoose_branch != 0) jets_pat_jetIDLoose();
	if (jets_pat_jetIDLooseAOD_branch != 0) jets_pat_jetIDLooseAOD();
	if (jets_pat_jetIDMinimal_branch != 0) jets_pat_jetIDMinimal();
	if (jets_pat_jetIDTight_branch != 0) jets_pat_jetIDTight();
	if (jets_pat_partonFlavour_branch != 0) jets_pat_partonFlavour();
	if (mus_pat_genID_branch != 0) mus_pat_genID();
	if (mus_pat_genMotherID_branch != 0) mus_pat_genMotherID();
	if (pfels_elsidx_branch != 0) pfels_elsidx();
	if (pfels_charge_branch != 0) pfels_charge();
	if (pfels_flag_branch != 0) pfels_flag();
	if (pfels_particleId_branch != 0) pfels_particleId();
	if (pfjets_chargedMultiplicity_branch != 0) pfjets_chargedMultiplicity();
	if (pfjets_muonMultiplicity_branch != 0) pfjets_muonMultiplicity();
	if (pfjets_neutralMultiplicity_branch != 0) pfjets_neutralMultiplicity();
	if (pfmus_musidx_branch != 0) pfmus_musidx();
	if (pfmus_charge_branch != 0) pfmus_charge();
	if (pfmus_flag_branch != 0) pfmus_flag();
	if (pfmus_particleId_branch != 0) pfmus_particleId();
	if (photons_fiduciality_branch != 0) photons_fiduciality();
	if (photons_scindex_branch != 0) photons_scindex();
	if (scs_detIdSeed_branch != 0) scs_detIdSeed();
	if (scs_elsidx_branch != 0) scs_elsidx();
	if (scs_severitySeed_branch != 0) scs_severitySeed();
	if (svs_isKs_branch != 0) svs_isKs();
	if (svs_isLambda_branch != 0) svs_isLambda();
	if (svs_mc3_id_branch != 0) svs_mc3_id();
	if (svs_nTrks_branch != 0) svs_nTrks();
	if (mus_tcmet_flag_branch != 0) mus_tcmet_flag();
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
	if (trks_elsidx_branch != 0) trks_elsidx();
	if (trk_musidx_branch != 0) trk_musidx();
	if (trkjets_ntrks_branch != 0) trkjets_ntrks();
	if (trkjets_vtxs_idx_branch != 0) trkjets_vtxs_idx();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (vtxs_tracksSize_branch != 0) vtxs_tracksSize();
	if (genps_lepdaughter_id_branch != 0) genps_lepdaughter_id();
	if (genps_lepdaughter_idx_branch != 0) genps_lepdaughter_idx();
	if (hlt_trigObjs_id_branch != 0) hlt_trigObjs_id();
	if (hyp_jets_idx_branch != 0) hyp_jets_idx();
	if (hyp_other_jets_idx_branch != 0) hyp_other_jets_idx();
	if (evt_nels_branch != 0) evt_nels();
	if (evt_detectorStatus_branch != 0) evt_detectorStatus();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (genps_flavorHistoryFilterResult_branch != 0) genps_flavorHistoryFilterResult();
	if (evt_ngenjets_branch != 0) evt_ngenjets();
	if (genps_signalProcessID_branch != 0) genps_signalProcessID();
	if (hlt_bits1_branch != 0) hlt_bits1();
	if (hlt_bits2_branch != 0) hlt_bits2();
	if (hlt_bits3_branch != 0) hlt_bits3();
	if (hlt_bits4_branch != 0) hlt_bits4();
	if (hlt_bits5_branch != 0) hlt_bits5();
	if (hlt_bits6_branch != 0) hlt_bits6();
	if (hlt_bits7_branch != 0) hlt_bits7();
	if (hlt_bits8_branch != 0) hlt_bits8();
	if (evt_njets_branch != 0) evt_njets();
	if (evt_njpts_branch != 0) evt_njpts();
	if (l1_bits1_branch != 0) l1_bits1();
	if (l1_bits2_branch != 0) l1_bits2();
	if (l1_bits3_branch != 0) l1_bits3();
	if (l1_bits4_branch != 0) l1_bits4();
	if (l1_techbits1_branch != 0) l1_techbits1();
	if (l1_techbits2_branch != 0) l1_techbits2();
	if (evt_nphotons_branch != 0) evt_nphotons();
	if (evt_ecalRecoStatus_branch != 0) evt_ecalRecoStatus();
	if (evt_nscs_branch != 0) evt_nscs();
	if (evt_ntrkjets_branch != 0) evt_ntrkjets();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (hlt_prescales_branch != 0) hlt_prescales();
	if (hyp_quadlep_bucket_branch != 0) hyp_quadlep_bucket();
	if (hyp_quadlep_first_index_branch != 0) hyp_quadlep_first_index();
	if (hyp_quadlep_fourth_index_branch != 0) hyp_quadlep_fourth_index();
	if (hyp_quadlep_second_index_branch != 0) hyp_quadlep_second_index();
	if (hyp_quadlep_third_index_branch != 0) hyp_quadlep_third_index();
	if (hyp_trilep_bucket_branch != 0) hyp_trilep_bucket();
	if (hyp_trilep_first_index_branch != 0) hyp_trilep_first_index();
	if (hyp_trilep_second_index_branch != 0) hyp_trilep_second_index();
	if (hyp_trilep_third_index_branch != 0) hyp_trilep_third_index();
	if (l1_prescales_branch != 0) l1_prescales();
	if (l1_techtrigprescales_branch != 0) l1_techtrigprescales();
	if (els_pat_flag_branch != 0) els_pat_flag();
	if (jets_pat_flag_branch != 0) jets_pat_flag();
	if (mus_pat_flag_branch != 0) mus_pat_flag();
	if (evt_nEvts_branch != 0) evt_nEvts();
	if (evt_filt_eff_branch != 0) evt_filt_eff();
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
	vector<TString> &hlt_trigNames()
	{
		if (not hlt_trigNames_isLoaded) {
			if (hlt_trigNames_branch != 0) {
				hlt_trigNames_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigNames_branch does not exist!\n");
				exit(1);
			}
			hlt_trigNames_isLoaded = true;
		}
		return hlt_trigNames_;
	}
	vector<TString> &l1_techtrigNames()
	{
		if (not l1_techtrigNames_isLoaded) {
			if (l1_techtrigNames_branch != 0) {
				l1_techtrigNames_branch->GetEntry(index);
			} else { 
				printf("branch l1_techtrigNames_branch does not exist!\n");
				exit(1);
			}
			l1_techtrigNames_isLoaded = true;
		}
		return l1_techtrigNames_;
	}
	vector<TString> &l1_trigNames()
	{
		if (not l1_trigNames_isLoaded) {
			if (l1_trigNames_branch != 0) {
				l1_trigNames_branch->GetEntry(index);
			} else { 
				printf("branch l1_trigNames_branch does not exist!\n");
				exit(1);
			}
			l1_trigNames_isLoaded = true;
		}
		return l1_trigNames_;
	}
	vector<TString> &evt_errCategory()
	{
		if (not evt_errCategory_isLoaded) {
			if (evt_errCategory_branch != 0) {
				evt_errCategory_branch->GetEntry(index);
			} else { 
				printf("branch evt_errCategory_branch does not exist!\n");
				exit(1);
			}
			evt_errCategory_isLoaded = true;
		}
		return evt_errCategory_;
	}
	vector<TString> &evt_errModule()
	{
		if (not evt_errModule_isLoaded) {
			if (evt_errModule_branch != 0) {
				evt_errModule_branch->GetEntry(index);
			} else { 
				printf("branch evt_errModule_branch does not exist!\n");
				exit(1);
			}
			evt_errModule_isLoaded = true;
		}
		return evt_errModule_;
	}
	vector<TString> &evt_errSeverity()
	{
		if (not evt_errSeverity_isLoaded) {
			if (evt_errSeverity_branch != 0) {
				evt_errSeverity_branch->GetEntry(index);
			} else { 
				printf("branch evt_errSeverity_branch does not exist!\n");
				exit(1);
			}
			evt_errSeverity_isLoaded = true;
		}
		return evt_errSeverity_;
	}
	bool &evt_eventHasHalo()
	{
		if (not evt_eventHasHalo_isLoaded) {
			if (evt_eventHasHalo_branch != 0) {
				evt_eventHasHalo_branch->GetEntry(index);
			} else { 
				printf("branch evt_eventHasHalo_branch does not exist!\n");
				exit(1);
			}
			evt_eventHasHalo_isLoaded = true;
		}
		return evt_eventHasHalo_;
	}
	bool &evt_hbheFilter()
	{
		if (not evt_hbheFilter_isLoaded) {
			if (evt_hbheFilter_branch != 0) {
				evt_hbheFilter_branch->GetEntry(index);
			} else { 
				printf("branch evt_hbheFilter_branch does not exist!\n");
				exit(1);
			}
			evt_hbheFilter_isLoaded = true;
		}
		return evt_hbheFilter_;
	}
	vector<bool> &mus_tightMatch()
	{
		if (not mus_tightMatch_isLoaded) {
			if (mus_tightMatch_branch != 0) {
				mus_tightMatch_branch->GetEntry(index);
			} else { 
				printf("branch mus_tightMatch_branch does not exist!\n");
				exit(1);
			}
			mus_tightMatch_isLoaded = true;
		}
		return mus_tightMatch_;
	}
	vector<bool> &mus_updatedSta()
	{
		if (not mus_updatedSta_isLoaded) {
			if (mus_updatedSta_branch != 0) {
				mus_updatedSta_branch->GetEntry(index);
			} else { 
				printf("branch mus_updatedSta_branch does not exist!\n");
				exit(1);
			}
			mus_updatedSta_isLoaded = true;
		}
		return mus_updatedSta_;
	}
	vector<bool> &photons_haspixelSeed()
	{
		if (not photons_haspixelSeed_isLoaded) {
			if (photons_haspixelSeed_branch != 0) {
				photons_haspixelSeed_branch->GetEntry(index);
			} else { 
				printf("branch photons_haspixelSeed_branch does not exist!\n");
				exit(1);
			}
			photons_haspixelSeed_isLoaded = true;
		}
		return photons_haspixelSeed_;
	}
	vector<double> &jets_closestElectron_DR()
	{
		if (not jets_closestElectron_DR_isLoaded) {
			if (jets_closestElectron_DR_branch != 0) {
				jets_closestElectron_DR_branch->GetEntry(index);
			} else { 
				printf("branch jets_closestElectron_DR_branch does not exist!\n");
				exit(1);
			}
			jets_closestElectron_DR_isLoaded = true;
		}
		return jets_closestElectron_DR_;
	}
	vector<double> &jets_closestMuon_DR()
	{
		if (not jets_closestMuon_DR_isLoaded) {
			if (jets_closestMuon_DR_branch != 0) {
				jets_closestMuon_DR_branch->GetEntry(index);
			} else { 
				printf("branch jets_closestMuon_DR_branch does not exist!\n");
				exit(1);
			}
			jets_closestMuon_DR_isLoaded = true;
		}
		return jets_closestMuon_DR_;
	}
	float &evt_bs_Xwidth()
	{
		if (not evt_bs_Xwidth_isLoaded) {
			if (evt_bs_Xwidth_branch != 0) {
				evt_bs_Xwidth_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_Xwidth_branch does not exist!\n");
				exit(1);
			}
			evt_bs_Xwidth_isLoaded = true;
		}
		return evt_bs_Xwidth_;
	}
	float &evt_bs_XwidthErr()
	{
		if (not evt_bs_XwidthErr_isLoaded) {
			if (evt_bs_XwidthErr_branch != 0) {
				evt_bs_XwidthErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_XwidthErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_XwidthErr_isLoaded = true;
		}
		return evt_bs_XwidthErr_;
	}
	float &evt_bs_Ywidth()
	{
		if (not evt_bs_Ywidth_isLoaded) {
			if (evt_bs_Ywidth_branch != 0) {
				evt_bs_Ywidth_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_Ywidth_branch does not exist!\n");
				exit(1);
			}
			evt_bs_Ywidth_isLoaded = true;
		}
		return evt_bs_Ywidth_;
	}
	float &evt_bs_YwidthErr()
	{
		if (not evt_bs_YwidthErr_isLoaded) {
			if (evt_bs_YwidthErr_branch != 0) {
				evt_bs_YwidthErr_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_YwidthErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_YwidthErr_isLoaded = true;
		}
		return evt_bs_YwidthErr_;
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
	float &evthcal_dmetx()
	{
		if (not evthcal_dmetx_isLoaded) {
			if (evthcal_dmetx_branch != 0) {
				evthcal_dmetx_branch->GetEntry(index);
			} else { 
				printf("branch evthcal_dmetx_branch does not exist!\n");
				exit(1);
			}
			evthcal_dmetx_isLoaded = true;
		}
		return evthcal_dmetx_;
	}
	float &evthcal_dmety()
	{
		if (not evthcal_dmety_isLoaded) {
			if (evthcal_dmety_branch != 0) {
				evthcal_dmety_branch->GetEntry(index);
			} else { 
				printf("branch evthcal_dmety_branch does not exist!\n");
				exit(1);
			}
			evthcal_dmety_isLoaded = true;
		}
		return evthcal_dmety_;
	}
	float &evthcal_dsumet()
	{
		if (not evthcal_dsumet_isLoaded) {
			if (evthcal_dsumet_branch != 0) {
				evthcal_dsumet_branch->GetEntry(index);
			} else { 
				printf("branch evthcal_dsumet_branch does not exist!\n");
				exit(1);
			}
			evthcal_dsumet_isLoaded = true;
		}
		return evthcal_dsumet_;
	}
	float &evthf_dmetx()
	{
		if (not evthf_dmetx_isLoaded) {
			if (evthf_dmetx_branch != 0) {
				evthf_dmetx_branch->GetEntry(index);
			} else { 
				printf("branch evthf_dmetx_branch does not exist!\n");
				exit(1);
			}
			evthf_dmetx_isLoaded = true;
		}
		return evthf_dmetx_;
	}
	float &evthf_dmety()
	{
		if (not evthf_dmety_isLoaded) {
			if (evthf_dmety_branch != 0) {
				evthf_dmety_branch->GetEntry(index);
			} else { 
				printf("branch evthf_dmety_branch does not exist!\n");
				exit(1);
			}
			evthf_dmety_isLoaded = true;
		}
		return evthf_dmety_;
	}
	float &evthf_dsumet()
	{
		if (not evthf_dsumet_isLoaded) {
			if (evthf_dsumet_branch != 0) {
				evthf_dsumet_branch->GetEntry(index);
			} else { 
				printf("branch evthf_dsumet_branch does not exist!\n");
				exit(1);
			}
			evthf_dsumet_isLoaded = true;
		}
		return evthf_dsumet_;
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
	float &evt_scale1fb()
	{
		if (not evt_scale1fb_isLoaded) {
			if (evt_scale1fb_branch != 0) {
				evt_scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch evt_scale1fb_branch does not exist!\n");
				exit(1);
			}
			evt_scale1fb_isLoaded = true;
		}
		return evt_scale1fb_;
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
	float &genps_alphaQCD()
	{
		if (not genps_alphaQCD_isLoaded) {
			if (genps_alphaQCD_branch != 0) {
				genps_alphaQCD_branch->GetEntry(index);
			} else { 
				printf("branch genps_alphaQCD_branch does not exist!\n");
				exit(1);
			}
			genps_alphaQCD_isLoaded = true;
		}
		return genps_alphaQCD_;
	}
	float &genps_pthat()
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
	float &genps_qScale()
	{
		if (not genps_qScale_isLoaded) {
			if (genps_qScale_branch != 0) {
				genps_qScale_branch->GetEntry(index);
			} else { 
				printf("branch genps_qScale_branch does not exist!\n");
				exit(1);
			}
			genps_qScale_isLoaded = true;
		}
		return genps_qScale_;
	}
	float &genps_weight()
	{
		if (not genps_weight_isLoaded) {
			if (genps_weight_branch != 0) {
				genps_weight_branch->GetEntry(index);
			} else { 
				printf("branch genps_weight_branch does not exist!\n");
				exit(1);
			}
			genps_weight_isLoaded = true;
		}
		return genps_weight_;
	}
	float &gen_sumEt()
	{
		if (not gen_sumEt_isLoaded) {
			if (gen_sumEt_branch != 0) {
				gen_sumEt_branch->GetEntry(index);
			} else { 
				printf("branch gen_sumEt_branch does not exist!\n");
				exit(1);
			}
			gen_sumEt_isLoaded = true;
		}
		return gen_sumEt_;
	}
	float &hcalnoise_eventChargeFraction()
	{
		if (not hcalnoise_eventChargeFraction_isLoaded) {
			if (hcalnoise_eventChargeFraction_branch != 0) {
				hcalnoise_eventChargeFraction_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_mht_mht_branch does not exist!\n");
				exit(1);
			}
			l1_mht_mht_isLoaded = true;
		}
		return l1_mht_mht_;
	}
	float &evt_ecalendcapm_met()
	{
		if (not evt_ecalendcapm_met_isLoaded) {
			if (evt_ecalendcapm_met_branch != 0) {
				evt_ecalendcapm_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalendcapm_met_branch does not exist!\n");
				exit(1);
			}
			evt_ecalendcapm_met_isLoaded = true;
		}
		return evt_ecalendcapm_met_;
	}
	float &evt_ecalendcapm_metPhi()
	{
		if (not evt_ecalendcapm_metPhi_isLoaded) {
			if (evt_ecalendcapm_metPhi_branch != 0) {
				evt_ecalendcapm_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalendcapm_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_ecalendcapm_metPhi_isLoaded = true;
		}
		return evt_ecalendcapm_metPhi_;
	}
	float &evt_ecalendcapp_met()
	{
		if (not evt_ecalendcapp_met_isLoaded) {
			if (evt_ecalendcapp_met_branch != 0) {
				evt_ecalendcapp_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalendcapp_met_branch does not exist!\n");
				exit(1);
			}
			evt_ecalendcapp_met_isLoaded = true;
		}
		return evt_ecalendcapp_met_;
	}
	float &evt_ecalendcapp_metPhi()
	{
		if (not evt_ecalendcapp_metPhi_isLoaded) {
			if (evt_ecalendcapp_metPhi_branch != 0) {
				evt_ecalendcapp_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalendcapp_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_ecalendcapp_metPhi_isLoaded = true;
		}
		return evt_ecalendcapp_metPhi_;
	}
	float &evt_ecalmet()
	{
		if (not evt_ecalmet_isLoaded) {
			if (evt_ecalmet_branch != 0) {
				evt_ecalmet_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalmet_branch does not exist!\n");
				exit(1);
			}
			evt_ecalmet_isLoaded = true;
		}
		return evt_ecalmet_;
	}
	float &evt_ecalmetPhi()
	{
		if (not evt_ecalmetPhi_isLoaded) {
			if (evt_ecalmetPhi_branch != 0) {
				evt_ecalmetPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_ecalmetPhi_isLoaded = true;
		}
		return evt_ecalmetPhi_;
	}
	float &evt_endcapm_met()
	{
		if (not evt_endcapm_met_isLoaded) {
			if (evt_endcapm_met_branch != 0) {
				evt_endcapm_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_endcapm_met_branch does not exist!\n");
				exit(1);
			}
			evt_endcapm_met_isLoaded = true;
		}
		return evt_endcapm_met_;
	}
	float &evt_endcapm_metPhi()
	{
		if (not evt_endcapm_metPhi_isLoaded) {
			if (evt_endcapm_metPhi_branch != 0) {
				evt_endcapm_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_endcapm_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_endcapm_metPhi_isLoaded = true;
		}
		return evt_endcapm_metPhi_;
	}
	float &evt_endcapp_met()
	{
		if (not evt_endcapp_met_isLoaded) {
			if (evt_endcapp_met_branch != 0) {
				evt_endcapp_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_endcapp_met_branch does not exist!\n");
				exit(1);
			}
			evt_endcapp_met_isLoaded = true;
		}
		return evt_endcapp_met_;
	}
	float &evt_endcapp_metPhi()
	{
		if (not evt_endcapp_metPhi_isLoaded) {
			if (evt_endcapp_metPhi_branch != 0) {
				evt_endcapp_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_endcapp_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_endcapp_metPhi_isLoaded = true;
		}
		return evt_endcapp_metPhi_;
	}
	float &evt_hcalendcapm_met()
	{
		if (not evt_hcalendcapm_met_isLoaded) {
			if (evt_hcalendcapm_met_branch != 0) {
				evt_hcalendcapm_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalendcapm_met_branch does not exist!\n");
				exit(1);
			}
			evt_hcalendcapm_met_isLoaded = true;
		}
		return evt_hcalendcapm_met_;
	}
	float &evt_hcalendcapm_metPhi()
	{
		if (not evt_hcalendcapm_metPhi_isLoaded) {
			if (evt_hcalendcapm_metPhi_branch != 0) {
				evt_hcalendcapm_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalendcapm_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_hcalendcapm_metPhi_isLoaded = true;
		}
		return evt_hcalendcapm_metPhi_;
	}
	float &evt_hcalendcapp_met()
	{
		if (not evt_hcalendcapp_met_isLoaded) {
			if (evt_hcalendcapp_met_branch != 0) {
				evt_hcalendcapp_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalendcapp_met_branch does not exist!\n");
				exit(1);
			}
			evt_hcalendcapp_met_isLoaded = true;
		}
		return evt_hcalendcapp_met_;
	}
	float &evt_hcalendcapp_metPhi()
	{
		if (not evt_hcalendcapp_metPhi_isLoaded) {
			if (evt_hcalendcapp_metPhi_branch != 0) {
				evt_hcalendcapp_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalendcapp_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_hcalendcapp_metPhi_isLoaded = true;
		}
		return evt_hcalendcapp_metPhi_;
	}
	float &evt_hcalmet()
	{
		if (not evt_hcalmet_isLoaded) {
			if (evt_hcalmet_branch != 0) {
				evt_hcalmet_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalmet_branch does not exist!\n");
				exit(1);
			}
			evt_hcalmet_isLoaded = true;
		}
		return evt_hcalmet_;
	}
	float &evt_hcalmetPhi()
	{
		if (not evt_hcalmetPhi_isLoaded) {
			if (evt_hcalmetPhi_branch != 0) {
				evt_hcalmetPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_hcalmetPhi_isLoaded = true;
		}
		return evt_hcalmetPhi_;
	}
	float &evt_met()
	{
		if (not evt_met_isLoaded) {
			if (evt_met_branch != 0) {
				evt_met_branch->GetEntry(index);
			} else { 
				printf("branch evt_met_branch does not exist!\n");
				exit(1);
			}
			evt_met_isLoaded = true;
		}
		return evt_met_;
	}
	float &evt_metHO()
	{
		if (not evt_metHO_isLoaded) {
			if (evt_metHO_branch != 0) {
				evt_metHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_metHO_branch does not exist!\n");
				exit(1);
			}
			evt_metHO_isLoaded = true;
		}
		return evt_metHO_;
	}
	float &evt_metHOPhi()
	{
		if (not evt_metHOPhi_isLoaded) {
			if (evt_metHOPhi_branch != 0) {
				evt_metHOPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metHOPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metHOPhi_isLoaded = true;
		}
		return evt_metHOPhi_;
	}
	float &evt_metHOSig()
	{
		if (not evt_metHOSig_isLoaded) {
			if (evt_metHOSig_branch != 0) {
				evt_metHOSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metHOSig_branch does not exist!\n");
				exit(1);
			}
			evt_metHOSig_isLoaded = true;
		}
		return evt_metHOSig_;
	}
	float &evt_metMuonCorr()
	{
		if (not evt_metMuonCorr_isLoaded) {
			if (evt_metMuonCorr_branch != 0) {
				evt_metMuonCorr_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonCorr_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonCorr_isLoaded = true;
		}
		return evt_metMuonCorr_;
	}
	float &evt_metMuonCorrPhi()
	{
		if (not evt_metMuonCorrPhi_isLoaded) {
			if (evt_metMuonCorrPhi_branch != 0) {
				evt_metMuonCorrPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonCorrPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonCorrPhi_isLoaded = true;
		}
		return evt_metMuonCorrPhi_;
	}
	float &evt_metMuonCorrSig()
	{
		if (not evt_metMuonCorrSig_isLoaded) {
			if (evt_metMuonCorrSig_branch != 0) {
				evt_metMuonCorrSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonCorrSig_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonCorrSig_isLoaded = true;
		}
		return evt_metMuonCorrSig_;
	}
	float &evt_metMuonJESCorr()
	{
		if (not evt_metMuonJESCorr_isLoaded) {
			if (evt_metMuonJESCorr_branch != 0) {
				evt_metMuonJESCorr_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonJESCorr_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonJESCorr_isLoaded = true;
		}
		return evt_metMuonJESCorr_;
	}
	float &evt_metMuonJESCorrPhi()
	{
		if (not evt_metMuonJESCorrPhi_isLoaded) {
			if (evt_metMuonJESCorrPhi_branch != 0) {
				evt_metMuonJESCorrPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonJESCorrPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonJESCorrPhi_isLoaded = true;
		}
		return evt_metMuonJESCorrPhi_;
	}
	float &evt_metMuonJESCorrSig()
	{
		if (not evt_metMuonJESCorrSig_isLoaded) {
			if (evt_metMuonJESCorrSig_branch != 0) {
				evt_metMuonJESCorrSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metMuonJESCorrSig_branch does not exist!\n");
				exit(1);
			}
			evt_metMuonJESCorrSig_isLoaded = true;
		}
		return evt_metMuonJESCorrSig_;
	}
	float &evt_metNoHF()
	{
		if (not evt_metNoHF_isLoaded) {
			if (evt_metNoHF_branch != 0) {
				evt_metNoHF_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHF_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHF_isLoaded = true;
		}
		return evt_metNoHF_;
	}
	float &evt_metNoHFHO()
	{
		if (not evt_metNoHFHO_isLoaded) {
			if (evt_metNoHFHO_branch != 0) {
				evt_metNoHFHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHFHO_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHFHO_isLoaded = true;
		}
		return evt_metNoHFHO_;
	}
	float &evt_metNoHFHOPhi()
	{
		if (not evt_metNoHFHOPhi_isLoaded) {
			if (evt_metNoHFHOPhi_branch != 0) {
				evt_metNoHFHOPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHFHOPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHFHOPhi_isLoaded = true;
		}
		return evt_metNoHFHOPhi_;
	}
	float &evt_metNoHFHOSig()
	{
		if (not evt_metNoHFHOSig_isLoaded) {
			if (evt_metNoHFHOSig_branch != 0) {
				evt_metNoHFHOSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHFHOSig_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHFHOSig_isLoaded = true;
		}
		return evt_metNoHFHOSig_;
	}
	float &evt_metNoHFPhi()
	{
		if (not evt_metNoHFPhi_isLoaded) {
			if (evt_metNoHFPhi_branch != 0) {
				evt_metNoHFPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHFPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHFPhi_isLoaded = true;
		}
		return evt_metNoHFPhi_;
	}
	float &evt_metNoHFSig()
	{
		if (not evt_metNoHFSig_isLoaded) {
			if (evt_metNoHFSig_branch != 0) {
				evt_metNoHFSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metNoHFSig_branch does not exist!\n");
				exit(1);
			}
			evt_metNoHFSig_isLoaded = true;
		}
		return evt_metNoHFSig_;
	}
	float &evt_metOpt()
	{
		if (not evt_metOpt_isLoaded) {
			if (evt_metOpt_branch != 0) {
				evt_metOpt_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOpt_branch does not exist!\n");
				exit(1);
			}
			evt_metOpt_isLoaded = true;
		}
		return evt_metOpt_;
	}
	float &evt_metOptHO()
	{
		if (not evt_metOptHO_isLoaded) {
			if (evt_metOptHO_branch != 0) {
				evt_metOptHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptHO_branch does not exist!\n");
				exit(1);
			}
			evt_metOptHO_isLoaded = true;
		}
		return evt_metOptHO_;
	}
	float &evt_metOptHOPhi()
	{
		if (not evt_metOptHOPhi_isLoaded) {
			if (evt_metOptHOPhi_branch != 0) {
				evt_metOptHOPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptHOPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metOptHOPhi_isLoaded = true;
		}
		return evt_metOptHOPhi_;
	}
	float &evt_metOptHOSig()
	{
		if (not evt_metOptHOSig_isLoaded) {
			if (evt_metOptHOSig_branch != 0) {
				evt_metOptHOSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptHOSig_branch does not exist!\n");
				exit(1);
			}
			evt_metOptHOSig_isLoaded = true;
		}
		return evt_metOptHOSig_;
	}
	float &evt_metOptNoHF()
	{
		if (not evt_metOptNoHF_isLoaded) {
			if (evt_metOptNoHF_branch != 0) {
				evt_metOptNoHF_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHF_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHF_isLoaded = true;
		}
		return evt_metOptNoHF_;
	}
	float &evt_metOptNoHFHO()
	{
		if (not evt_metOptNoHFHO_isLoaded) {
			if (evt_metOptNoHFHO_branch != 0) {
				evt_metOptNoHFHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHFHO_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHFHO_isLoaded = true;
		}
		return evt_metOptNoHFHO_;
	}
	float &evt_metOptNoHFHOPhi()
	{
		if (not evt_metOptNoHFHOPhi_isLoaded) {
			if (evt_metOptNoHFHOPhi_branch != 0) {
				evt_metOptNoHFHOPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHFHOPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHFHOPhi_isLoaded = true;
		}
		return evt_metOptNoHFHOPhi_;
	}
	float &evt_metOptNoHFHOSig()
	{
		if (not evt_metOptNoHFHOSig_isLoaded) {
			if (evt_metOptNoHFHOSig_branch != 0) {
				evt_metOptNoHFHOSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHFHOSig_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHFHOSig_isLoaded = true;
		}
		return evt_metOptNoHFHOSig_;
	}
	float &evt_metOptNoHFPhi()
	{
		if (not evt_metOptNoHFPhi_isLoaded) {
			if (evt_metOptNoHFPhi_branch != 0) {
				evt_metOptNoHFPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHFPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHFPhi_isLoaded = true;
		}
		return evt_metOptNoHFPhi_;
	}
	float &evt_metOptNoHFSig()
	{
		if (not evt_metOptNoHFSig_isLoaded) {
			if (evt_metOptNoHFSig_branch != 0) {
				evt_metOptNoHFSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptNoHFSig_branch does not exist!\n");
				exit(1);
			}
			evt_metOptNoHFSig_isLoaded = true;
		}
		return evt_metOptNoHFSig_;
	}
	float &evt_metOptPhi()
	{
		if (not evt_metOptPhi_isLoaded) {
			if (evt_metOptPhi_branch != 0) {
				evt_metOptPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metOptPhi_isLoaded = true;
		}
		return evt_metOptPhi_;
	}
	float &evt_metOptSig()
	{
		if (not evt_metOptSig_isLoaded) {
			if (evt_metOptSig_branch != 0) {
				evt_metOptSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metOptSig_branch does not exist!\n");
				exit(1);
			}
			evt_metOptSig_isLoaded = true;
		}
		return evt_metOptSig_;
	}
	float &evt_metPhi()
	{
		if (not evt_metPhi_isLoaded) {
			if (evt_metPhi_branch != 0) {
				evt_metPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_metPhi_branch does not exist!\n");
				exit(1);
			}
			evt_metPhi_isLoaded = true;
		}
		return evt_metPhi_;
	}
	float &evt_metSig()
	{
		if (not evt_metSig_isLoaded) {
			if (evt_metSig_branch != 0) {
				evt_metSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_metSig_branch does not exist!\n");
				exit(1);
			}
			evt_metSig_isLoaded = true;
		}
		return evt_metSig_;
	}
	float &evt_sumet()
	{
		if (not evt_sumet_isLoaded) {
			if (evt_sumet_branch != 0) {
				evt_sumet_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumet_branch does not exist!\n");
				exit(1);
			}
			evt_sumet_isLoaded = true;
		}
		return evt_sumet_;
	}
	float &evt_sumetHO()
	{
		if (not evt_sumetHO_isLoaded) {
			if (evt_sumetHO_branch != 0) {
				evt_sumetHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetHO_branch does not exist!\n");
				exit(1);
			}
			evt_sumetHO_isLoaded = true;
		}
		return evt_sumetHO_;
	}
	float &evt_sumetMuonCorr()
	{
		if (not evt_sumetMuonCorr_isLoaded) {
			if (evt_sumetMuonCorr_branch != 0) {
				evt_sumetMuonCorr_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetMuonCorr_branch does not exist!\n");
				exit(1);
			}
			evt_sumetMuonCorr_isLoaded = true;
		}
		return evt_sumetMuonCorr_;
	}
	float &evt_sumetNoHF()
	{
		if (not evt_sumetNoHF_isLoaded) {
			if (evt_sumetNoHF_branch != 0) {
				evt_sumetNoHF_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetNoHF_branch does not exist!\n");
				exit(1);
			}
			evt_sumetNoHF_isLoaded = true;
		}
		return evt_sumetNoHF_;
	}
	float &evt_sumetNoHFHO()
	{
		if (not evt_sumetNoHFHO_isLoaded) {
			if (evt_sumetNoHFHO_branch != 0) {
				evt_sumetNoHFHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetNoHFHO_branch does not exist!\n");
				exit(1);
			}
			evt_sumetNoHFHO_isLoaded = true;
		}
		return evt_sumetNoHFHO_;
	}
	float &evt_sumetOpt()
	{
		if (not evt_sumetOpt_isLoaded) {
			if (evt_sumetOpt_branch != 0) {
				evt_sumetOpt_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetOpt_branch does not exist!\n");
				exit(1);
			}
			evt_sumetOpt_isLoaded = true;
		}
		return evt_sumetOpt_;
	}
	float &evt_sumetOptHO()
	{
		if (not evt_sumetOptHO_isLoaded) {
			if (evt_sumetOptHO_branch != 0) {
				evt_sumetOptHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetOptHO_branch does not exist!\n");
				exit(1);
			}
			evt_sumetOptHO_isLoaded = true;
		}
		return evt_sumetOptHO_;
	}
	float &evt_sumetOptNoHF()
	{
		if (not evt_sumetOptNoHF_isLoaded) {
			if (evt_sumetOptNoHF_branch != 0) {
				evt_sumetOptNoHF_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetOptNoHF_branch does not exist!\n");
				exit(1);
			}
			evt_sumetOptNoHF_isLoaded = true;
		}
		return evt_sumetOptNoHF_;
	}
	float &evt_sumetOptNoHFHO()
	{
		if (not evt_sumetOptNoHFHO_isLoaded) {
			if (evt_sumetOptNoHFHO_branch != 0) {
				evt_sumetOptNoHFHO_branch->GetEntry(index);
			} else { 
				printf("branch evt_sumetOptNoHFHO_branch does not exist!\n");
				exit(1);
			}
			evt_sumetOptNoHFHO_isLoaded = true;
		}
		return evt_sumetOptNoHFHO_;
	}
	float &met_pat_metCor()
	{
		if (not met_pat_metCor_isLoaded) {
			if (met_pat_metCor_branch != 0) {
				met_pat_metCor_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metCor_branch does not exist!\n");
				exit(1);
			}
			met_pat_metCor_isLoaded = true;
		}
		return met_pat_metCor_;
	}
	float &met_pat_metPhiCor()
	{
		if (not met_pat_metPhiCor_isLoaded) {
			if (met_pat_metPhiCor_branch != 0) {
				met_pat_metPhiCor_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metPhiCor_branch does not exist!\n");
				exit(1);
			}
			met_pat_metPhiCor_isLoaded = true;
		}
		return met_pat_metPhiCor_;
	}
	float &met_pat_metPhiUncor()
	{
		if (not met_pat_metPhiUncor_isLoaded) {
			if (met_pat_metPhiUncor_branch != 0) {
				met_pat_metPhiUncor_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metPhiUncor_branch does not exist!\n");
				exit(1);
			}
			met_pat_metPhiUncor_isLoaded = true;
		}
		return met_pat_metPhiUncor_;
	}
	float &met_pat_metPhiUncorJES()
	{
		if (not met_pat_metPhiUncorJES_isLoaded) {
			if (met_pat_metPhiUncorJES_branch != 0) {
				met_pat_metPhiUncorJES_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metPhiUncorJES_branch does not exist!\n");
				exit(1);
			}
			met_pat_metPhiUncorJES_isLoaded = true;
		}
		return met_pat_metPhiUncorJES_;
	}
	float &met_pat_metPhiUncorMuon()
	{
		if (not met_pat_metPhiUncorMuon_isLoaded) {
			if (met_pat_metPhiUncorMuon_branch != 0) {
				met_pat_metPhiUncorMuon_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metPhiUncorMuon_branch does not exist!\n");
				exit(1);
			}
			met_pat_metPhiUncorMuon_isLoaded = true;
		}
		return met_pat_metPhiUncorMuon_;
	}
	float &met_pat_metUncor()
	{
		if (not met_pat_metUncor_isLoaded) {
			if (met_pat_metUncor_branch != 0) {
				met_pat_metUncor_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metUncor_branch does not exist!\n");
				exit(1);
			}
			met_pat_metUncor_isLoaded = true;
		}
		return met_pat_metUncor_;
	}
	float &met_pat_metUncorJES()
	{
		if (not met_pat_metUncorJES_isLoaded) {
			if (met_pat_metUncorJES_branch != 0) {
				met_pat_metUncorJES_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metUncorJES_branch does not exist!\n");
				exit(1);
			}
			met_pat_metUncorJES_isLoaded = true;
		}
		return met_pat_metUncorJES_;
	}
	float &met_pat_metUncorMuon()
	{
		if (not met_pat_metUncorMuon_isLoaded) {
			if (met_pat_metUncorMuon_branch != 0) {
				met_pat_metUncorMuon_branch->GetEntry(index);
			} else { 
				printf("branch met_pat_metUncorMuon_branch does not exist!\n");
				exit(1);
			}
			met_pat_metUncorMuon_isLoaded = true;
		}
		return met_pat_metUncorMuon_;
	}
	float &pdfinfo_scale()
	{
		if (not pdfinfo_scale_isLoaded) {
			if (pdfinfo_scale_branch != 0) {
				pdfinfo_scale_branch->GetEntry(index);
			} else { 
				printf("branch pdfinfo_scale_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_scale_isLoaded = true;
		}
		return pdfinfo_scale_;
	}
	float &pdfinfo_x1()
	{
		if (not pdfinfo_x1_isLoaded) {
			if (pdfinfo_x1_branch != 0) {
				pdfinfo_x1_branch->GetEntry(index);
			} else { 
				printf("branch pdfinfo_x1_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_x1_isLoaded = true;
		}
		return pdfinfo_x1_;
	}
	float &pdfinfo_x2()
	{
		if (not pdfinfo_x2_isLoaded) {
			if (pdfinfo_x2_branch != 0) {
				pdfinfo_x2_branch->GetEntry(index);
			} else { 
				printf("branch pdfinfo_x2_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_x2_isLoaded = true;
		}
		return pdfinfo_x2_;
	}
	float &evt_pfmet()
	{
		if (not evt_pfmet_isLoaded) {
			if (evt_pfmet_branch != 0) {
				evt_pfmet_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfmet_branch does not exist!\n");
				exit(1);
			}
			evt_pfmet_isLoaded = true;
		}
		return evt_pfmet_;
	}
	float &evt_pfmetPhi()
	{
		if (not evt_pfmetPhi_isLoaded) {
			if (evt_pfmetPhi_branch != 0) {
				evt_pfmetPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_pfmetPhi_isLoaded = true;
		}
		return evt_pfmetPhi_;
	}
	float &evt_pfmetSig()
	{
		if (not evt_pfmetSig_isLoaded) {
			if (evt_pfmetSig_branch != 0) {
				evt_pfmetSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfmetSig_branch does not exist!\n");
				exit(1);
			}
			evt_pfmetSig_isLoaded = true;
		}
		return evt_pfmetSig_;
	}
	float &evt_pfsumet()
	{
		if (not evt_pfsumet_isLoaded) {
			if (evt_pfsumet_branch != 0) {
				evt_pfsumet_branch->GetEntry(index);
			} else { 
				printf("branch evt_pfsumet_branch does not exist!\n");
				exit(1);
			}
			evt_pfsumet_isLoaded = true;
		}
		return evt_pfsumet_;
	}
	float &evt_tcmet()
	{
		if (not evt_tcmet_isLoaded) {
			if (evt_tcmet_branch != 0) {
				evt_tcmet_branch->GetEntry(index);
			} else { 
				printf("branch evt_tcmet_branch does not exist!\n");
				exit(1);
			}
			evt_tcmet_isLoaded = true;
		}
		return evt_tcmet_;
	}
	float &evt_tcmetPhi()
	{
		if (not evt_tcmetPhi_isLoaded) {
			if (evt_tcmetPhi_branch != 0) {
				evt_tcmetPhi_branch->GetEntry(index);
			} else { 
				printf("branch evt_tcmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_tcmetPhi_isLoaded = true;
		}
		return evt_tcmetPhi_;
	}
	float &evt_tcmetSig()
	{
		if (not evt_tcmetSig_isLoaded) {
			if (evt_tcmetSig_branch != 0) {
				evt_tcmetSig_branch->GetEntry(index);
			} else { 
				printf("branch evt_tcmetSig_branch does not exist!\n");
				exit(1);
			}
			evt_tcmetSig_isLoaded = true;
		}
		return evt_tcmetSig_;
	}
	float &evt_tcsumet()
	{
		if (not evt_tcsumet_isLoaded) {
			if (evt_tcsumet_branch != 0) {
				evt_tcsumet_branch->GetEntry(index);
			} else { 
				printf("branch evt_tcsumet_branch does not exist!\n");
				exit(1);
			}
			evt_tcsumet_isLoaded = true;
		}
		return evt_tcsumet_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4()
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
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_met_p4()
	{
		if (not l1_met_p4_isLoaded) {
			if (l1_met_p4_branch != 0) {
				l1_met_p4_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_mht_p4_branch does not exist!\n");
				exit(1);
			}
			l1_mht_p4_isLoaded = true;
		}
		return l1_mht_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_motherp4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_gp_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_motherp4()
	{
		if (not jets_mc_motherp4_isLoaded) {
			if (jets_mc_motherp4_branch != 0) {
				jets_mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			jets_mc_motherp4_isLoaded = true;
		}
		return jets_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_motherp4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_gp_p4()
	{
		if (not pfjets_mc_gp_p4_isLoaded) {
			if (pfjets_mc_gp_p4_branch != 0) {
				pfjets_mc_gp_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_gp_p4_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_gp_p4_isLoaded = true;
		}
		return pfjets_mc_gp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_motherp4()
	{
		if (not pfjets_mc_motherp4_isLoaded) {
			if (pfjets_mc_motherp4_branch != 0) {
				pfjets_mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_motherp4_isLoaded = true;
		}
		return pfjets_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_p4()
	{
		if (not pfjets_mc_p4_isLoaded) {
			if (pfjets_mc_p4_branch != 0) {
				pfjets_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_p4_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_p4_isLoaded = true;
		}
		return pfjets_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_motherp4()
	{
		if (not photons_mc_motherp4_isLoaded) {
			if (photons_mc_motherp4_branch != 0) {
				photons_mc_motherp4_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc_motherp4_branch does not exist!\n");
				exit(1);
			}
			photons_mc_motherp4_isLoaded = true;
		}
		return photons_mc_motherp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_p4()
	{
		if (not photons_mc_p4_isLoaded) {
			if (photons_mc_p4_branch != 0) {
				photons_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc_p4_branch does not exist!\n");
				exit(1);
			}
			photons_mc_p4_isLoaded = true;
		}
		return photons_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trk_mcp4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_conv_pos_p4()
	{
		if (not els_conv_pos_p4_isLoaded) {
			if (els_conv_pos_p4_branch != 0) {
				els_conv_pos_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_pos_p4_branch does not exist!\n");
				exit(1);
			}
			els_conv_pos_p4_isLoaded = true;
		}
		return els_conv_pos_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_inner_position()
	{
		if (not els_inner_position_isLoaded) {
			if (els_inner_position_branch != 0) {
				els_inner_position_branch->GetEntry(index);
			} else { 
				printf("branch els_inner_position_branch does not exist!\n");
				exit(1);
			}
			els_inner_position_isLoaded = true;
		}
		return els_inner_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_outer_position()
	{
		if (not els_outer_position_isLoaded) {
			if (els_outer_position_branch != 0) {
				els_outer_position_branch->GetEntry(index);
			} else { 
				printf("branch els_outer_position_branch does not exist!\n");
				exit(1);
			}
			els_outer_position_isLoaded = true;
		}
		return els_outer_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4()
	{
		if (not genjets_p4_isLoaded) {
			if (genjets_p4_branch != 0) {
				genjets_p4_branch->GetEntry(index);
			} else { 
				printf("branch genjets_p4_branch does not exist!\n");
				exit(1);
			}
			genjets_p4_isLoaded = true;
		}
		return genjets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_prod_vtx()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_inner_position()
	{
		if (not gsftrks_inner_position_isLoaded) {
			if (gsftrks_inner_position_branch != 0) {
				gsftrks_inner_position_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_inner_position_branch does not exist!\n");
				exit(1);
			}
			gsftrks_inner_position_isLoaded = true;
		}
		return gsftrks_inner_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_p4()
	{
		if (not gsftrks_outer_p4_isLoaded) {
			if (gsftrks_outer_p4_branch != 0) {
				gsftrks_outer_p4_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_outer_p4_branch does not exist!\n");
				exit(1);
			}
			gsftrks_outer_p4_isLoaded = true;
		}
		return gsftrks_outer_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_position()
	{
		if (not gsftrks_outer_position_isLoaded) {
			if (gsftrks_outer_position_branch != 0) {
				gsftrks_outer_position_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_outer_position_branch does not exist!\n");
				exit(1);
			}
			gsftrks_outer_position_isLoaded = true;
		}
		return gsftrks_outer_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_p4()
	{
		if (not gsftrks_p4_isLoaded) {
			if (gsftrks_p4_branch != 0) {
				gsftrks_p4_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_p4_branch does not exist!\n");
				exit(1);
			}
			gsftrks_p4_isLoaded = true;
		}
		return gsftrks_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_vertex_p4()
	{
		if (not gsftrks_vertex_p4_isLoaded) {
			if (gsftrks_vertex_p4_branch != 0) {
				gsftrks_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			gsftrks_vertex_p4_isLoaded = true;
		}
		return gsftrks_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_p4()
	{
		if (not hyp_ll_p4_isLoaded) {
			if (hyp_ll_p4_branch != 0) {
				hyp_ll_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_p4_isLoaded = true;
		}
		return hyp_ll_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_trk_p4()
	{
		if (not hyp_ll_trk_p4_isLoaded) {
			if (hyp_ll_trk_p4_branch != 0) {
				hyp_ll_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_trk_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_trk_p4_isLoaded = true;
		}
		return hyp_ll_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_p4()
	{
		if (not hyp_lt_p4_isLoaded) {
			if (hyp_lt_p4_branch != 0) {
				hyp_lt_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_p4_isLoaded = true;
		}
		return hyp_lt_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_trk_p4()
	{
		if (not hyp_lt_trk_p4_isLoaded) {
			if (hyp_lt_trk_p4_branch != 0) {
				hyp_lt_trk_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_trk_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_trk_p4_isLoaded = true;
		}
		return hyp_lt_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_p4()
	{
		if (not hyp_p4_isLoaded) {
			if (hyp_p4_branch != 0) {
				hyp_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_p4_isLoaded = true;
		}
		return hyp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_p4()
	{
		if (not hyp_FVFit_p4_isLoaded) {
			if (hyp_FVFit_p4_branch != 0) {
				hyp_FVFit_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_p4_isLoaded = true;
		}
		return hyp_FVFit_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_v4()
	{
		if (not hyp_FVFit_v4_isLoaded) {
			if (hyp_FVFit_v4_branch != 0) {
				hyp_FVFit_v4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4_isLoaded = true;
		}
		return hyp_FVFit_v4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_mc_p4()
	{
		if (not hyp_ll_mc_p4_isLoaded) {
			if (hyp_ll_mc_p4_branch != 0) {
				hyp_ll_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_mc_p4_isLoaded = true;
		}
		return hyp_ll_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_mc_p4()
	{
		if (not hyp_lt_mc_p4_isLoaded) {
			if (hyp_lt_mc_p4_branch != 0) {
				hyp_lt_mc_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_mc_p4_isLoaded = true;
		}
		return hyp_lt_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_p4()
	{
		if (not jets_p4_isLoaded) {
			if (jets_p4_branch != 0) {
				jets_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_p4_branch does not exist!\n");
				exit(1);
			}
			jets_p4_isLoaded = true;
		}
		return jets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_vertex_p4()
	{
		if (not jets_vertex_p4_isLoaded) {
			if (jets_vertex_p4_branch != 0) {
				jets_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			jets_vertex_p4_isLoaded = true;
		}
		return jets_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jpts_p4()
	{
		if (not jpts_p4_isLoaded) {
			if (jpts_p4_branch != 0) {
				jpts_p4_branch->GetEntry(index);
			} else { 
				printf("branch jpts_p4_branch does not exist!\n");
				exit(1);
			}
			jpts_p4_isLoaded = true;
		}
		return jpts_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emiso_p4()
	{
		if (not l1_emiso_p4_isLoaded) {
			if (l1_emiso_p4_branch != 0) {
				l1_emiso_p4_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_mus_p4_branch does not exist!\n");
				exit(1);
			}
			l1_mus_p4_isLoaded = true;
		}
		return l1_mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_ecalpos_p4()
	{
		if (not mus_ecalpos_p4_isLoaded) {
			if (mus_ecalpos_p4_branch != 0) {
				mus_ecalpos_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_ecalpos_p4_branch does not exist!\n");
				exit(1);
			}
			mus_ecalpos_p4_isLoaded = true;
		}
		return mus_ecalpos_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitdefault_p4()
	{
		if (not mus_fitdefault_p4_isLoaded) {
			if (mus_fitdefault_p4_branch != 0) {
				mus_fitdefault_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_fitdefault_p4_branch does not exist!\n");
				exit(1);
			}
			mus_fitdefault_p4_isLoaded = true;
		}
		return mus_fitdefault_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitfirsthit_p4()
	{
		if (not mus_fitfirsthit_p4_isLoaded) {
			if (mus_fitfirsthit_p4_branch != 0) {
				mus_fitfirsthit_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_fitfirsthit_p4_branch does not exist!\n");
				exit(1);
			}
			mus_fitfirsthit_p4_isLoaded = true;
		}
		return mus_fitfirsthit_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitpicky_p4()
	{
		if (not mus_fitpicky_p4_isLoaded) {
			if (mus_fitpicky_p4_branch != 0) {
				mus_fitpicky_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_fitpicky_p4_branch does not exist!\n");
				exit(1);
			}
			mus_fitpicky_p4_isLoaded = true;
		}
		return mus_fitpicky_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fittev_p4()
	{
		if (not mus_fittev_p4_isLoaded) {
			if (mus_fittev_p4_branch != 0) {
				mus_fittev_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_fittev_p4_branch does not exist!\n");
				exit(1);
			}
			mus_fittev_p4_isLoaded = true;
		}
		return mus_fittev_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_outerPos_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_p4()
	{
		if (not mus_gfit_p4_isLoaded) {
			if (mus_gfit_p4_branch != 0) {
				mus_gfit_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_p4_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_p4_isLoaded = true;
		}
		return mus_gfit_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_vertex_p4()
	{
		if (not mus_gfit_vertex_p4_isLoaded) {
			if (mus_gfit_vertex_p4_branch != 0) {
				mus_gfit_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_vertex_p4_isLoaded = true;
		}
		return mus_gfit_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_p4()
	{
		if (not mus_sta_p4_isLoaded) {
			if (mus_sta_p4_branch != 0) {
				mus_sta_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_p4_branch does not exist!\n");
				exit(1);
			}
			mus_sta_p4_isLoaded = true;
		}
		return mus_sta_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_vertex_p4()
	{
		if (not mus_sta_vertex_p4_isLoaded) {
			if (mus_sta_vertex_p4_branch != 0) {
				mus_sta_vertex_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			mus_sta_vertex_p4_isLoaded = true;
		}
		return mus_sta_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_trk_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_vertex_p4()
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genMotherP4()
	{
		if (not els_pat_genMotherP4_isLoaded) {
			if (els_pat_genMotherP4_branch != 0) {
				els_pat_genMotherP4_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_genMotherP4_branch does not exist!\n");
				exit(1);
			}
			els_pat_genMotherP4_isLoaded = true;
		}
		return els_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genP4()
	{
		if (not els_pat_genP4_isLoaded) {
			if (els_pat_genP4_branch != 0) {
				els_pat_genP4_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_genP4_branch does not exist!\n");
				exit(1);
			}
			els_pat_genP4_isLoaded = true;
		}
		return els_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_p4()
	{
		if (not els_pat_p4_isLoaded) {
			if (els_pat_p4_branch != 0) {
				els_pat_p4_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_p4_branch does not exist!\n");
				exit(1);
			}
			els_pat_p4_isLoaded = true;
		}
		return els_pat_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genJet_p4()
	{
		if (not jets_pat_genJet_p4_isLoaded) {
			if (jets_pat_genJet_p4_branch != 0) {
				jets_pat_genJet_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_genJet_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genJet_p4_isLoaded = true;
		}
		return jets_pat_genJet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genPartonMother_p4()
	{
		if (not jets_pat_genPartonMother_p4_isLoaded) {
			if (jets_pat_genPartonMother_p4_branch != 0) {
				jets_pat_genPartonMother_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_genPartonMother_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return jets_pat_genPartonMother_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genParton_p4()
	{
		if (not jets_pat_genParton_p4_isLoaded) {
			if (jets_pat_genParton_p4_branch != 0) {
				jets_pat_genParton_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_genParton_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genParton_p4_isLoaded = true;
		}
		return jets_pat_genParton_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_p4()
	{
		if (not jets_pat_jet_p4_isLoaded) {
			if (jets_pat_jet_p4_branch != 0) {
				jets_pat_jet_p4_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jet_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jet_p4_isLoaded = true;
		}
		return jets_pat_jet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_uncorp4()
	{
		if (not jets_pat_jet_uncorp4_isLoaded) {
			if (jets_pat_jet_uncorp4_branch != 0) {
				jets_pat_jet_uncorp4_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jet_uncorp4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jet_uncorp4_isLoaded = true;
		}
		return jets_pat_jet_uncorp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genMotherP4()
	{
		if (not mus_pat_genMotherP4_isLoaded) {
			if (mus_pat_genMotherP4_branch != 0) {
				mus_pat_genMotherP4_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_genMotherP4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genMotherP4_isLoaded = true;
		}
		return mus_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genP4()
	{
		if (not mus_pat_genP4_isLoaded) {
			if (mus_pat_genP4_branch != 0) {
				mus_pat_genP4_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_genP4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genP4_isLoaded = true;
		}
		return mus_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_p4()
	{
		if (not mus_pat_p4_isLoaded) {
			if (mus_pat_p4_branch != 0) {
				mus_pat_p4_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_p4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_p4_isLoaded = true;
		}
		return mus_pat_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_p4()
	{
		if (not pfels_p4_isLoaded) {
			if (pfels_p4_branch != 0) {
				pfels_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfels_p4_branch does not exist!\n");
				exit(1);
			}
			pfels_p4_isLoaded = true;
		}
		return pfels_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_posAtEcal_p4()
	{
		if (not pfels_posAtEcal_p4_isLoaded) {
			if (pfels_posAtEcal_p4_branch != 0) {
				pfels_posAtEcal_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfels_posAtEcal_p4_branch does not exist!\n");
				exit(1);
			}
			pfels_posAtEcal_p4_isLoaded = true;
		}
		return pfels_posAtEcal_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4()
	{
		if (not pfjets_p4_isLoaded) {
			if (pfjets_p4_branch != 0) {
				pfjets_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_p4_branch does not exist!\n");
				exit(1);
			}
			pfjets_p4_isLoaded = true;
		}
		return pfjets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_p4()
	{
		if (not pfmus_p4_isLoaded) {
			if (pfmus_p4_branch != 0) {
				pfmus_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_p4_branch does not exist!\n");
				exit(1);
			}
			pfmus_p4_isLoaded = true;
		}
		return pfmus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_posAtEcal_p4()
	{
		if (not pfmus_posAtEcal_p4_isLoaded) {
			if (pfmus_posAtEcal_p4_branch != 0) {
				pfmus_posAtEcal_p4_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_posAtEcal_p4_branch does not exist!\n");
				exit(1);
			}
			pfmus_posAtEcal_p4_isLoaded = true;
		}
		return pfmus_posAtEcal_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4()
	{
		if (not photons_p4_isLoaded) {
			if (photons_p4_branch != 0) {
				photons_p4_branch->GetEntry(index);
			} else { 
				printf("branch photons_p4_branch does not exist!\n");
				exit(1);
			}
			photons_p4_isLoaded = true;
		}
		return photons_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_p4()
	{
		if (not scs_p4_isLoaded) {
			if (scs_p4_branch != 0) {
				scs_p4_branch->GetEntry(index);
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
			} else { 
				printf("branch scs_vtx_p4_branch does not exist!\n");
				exit(1);
			}
			scs_vtx_p4_isLoaded = true;
		}
		return scs_vtx_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_flight()
	{
		if (not svs_flight_isLoaded) {
			if (svs_flight_branch != 0) {
				svs_flight_branch->GetEntry(index);
			} else { 
				printf("branch svs_flight_branch does not exist!\n");
				exit(1);
			}
			svs_flight_isLoaded = true;
		}
		return svs_flight_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_mc3_p4()
	{
		if (not svs_mc3_p4_isLoaded) {
			if (svs_mc3_p4_branch != 0) {
				svs_mc3_p4_branch->GetEntry(index);
			} else { 
				printf("branch svs_mc3_p4_branch does not exist!\n");
				exit(1);
			}
			svs_mc3_p4_isLoaded = true;
		}
		return svs_mc3_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_p4()
	{
		if (not svs_p4_isLoaded) {
			if (svs_p4_branch != 0) {
				svs_p4_branch->GetEntry(index);
			} else { 
				printf("branch svs_p4_branch does not exist!\n");
				exit(1);
			}
			svs_p4_isLoaded = true;
		}
		return svs_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_position()
	{
		if (not svs_position_isLoaded) {
			if (svs_position_branch != 0) {
				svs_position_branch->GetEntry(index);
			} else { 
				printf("branch svs_position_branch does not exist!\n");
				exit(1);
			}
			svs_position_isLoaded = true;
		}
		return svs_position_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_refitp4()
	{
		if (not svs_refitp4_isLoaded) {
			if (svs_refitp4_branch != 0) {
				svs_refitp4_branch->GetEntry(index);
			} else { 
				printf("branch svs_refitp4_branch does not exist!\n");
				exit(1);
			}
			svs_refitp4_isLoaded = true;
		}
		return svs_refitp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_inner_position()
	{
		if (not trks_inner_position_isLoaded) {
			if (trks_inner_position_branch != 0) {
				trks_inner_position_branch->GetEntry(index);
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
			} else { 
				printf("branch trks_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			trks_vertex_p4_isLoaded = true;
		}
		return trks_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trkjets_p4()
	{
		if (not trkjets_p4_isLoaded) {
			if (trkjets_p4_branch != 0) {
				trkjets_p4_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_p4_branch does not exist!\n");
				exit(1);
			}
			trkjets_p4_isLoaded = true;
		}
		return trkjets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position()
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &genps_lepdaughter_p4()
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4()
	{
		if (not hlt_trigObjs_p4_isLoaded) {
			if (hlt_trigObjs_p4_branch != 0) {
				hlt_trigObjs_p4_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigObjs_p4_branch does not exist!\n");
				exit(1);
			}
			hlt_trigObjs_p4_isLoaded = true;
		}
		return hlt_trigObjs_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_jets_p4()
	{
		if (not hyp_jets_p4_isLoaded) {
			if (hyp_jets_p4_branch != 0) {
				hyp_jets_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_jets_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_p4_isLoaded = true;
		}
		return hyp_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_other_jets_p4()
	{
		if (not hyp_other_jets_p4_isLoaded) {
			if (hyp_other_jets_p4_branch != 0) {
				hyp_other_jets_p4_branch->GetEntry(index);
			} else { 
				printf("branch hyp_other_jets_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_p4_isLoaded = true;
		}
		return hyp_other_jets_p4_;
	}
	vector<float> &jpts_combinedSecondaryVertexBJetTag()
	{
		if (not jpts_combinedSecondaryVertexBJetTag_isLoaded) {
			if (jpts_combinedSecondaryVertexBJetTag_branch != 0) {
				jpts_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return jpts_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &jpts_combinedSecondaryVertexMVABJetTag()
	{
		if (not jpts_combinedSecondaryVertexMVABJetTag_isLoaded) {
			if (jpts_combinedSecondaryVertexMVABJetTag_branch != 0) {
				jpts_combinedSecondaryVertexMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return jpts_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &jpts_jetBProbabilityBJetTag()
	{
		if (not jpts_jetBProbabilityBJetTag_isLoaded) {
			if (jpts_jetBProbabilityBJetTag_branch != 0) {
				jpts_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_jetBProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_jetBProbabilityBJetTag_isLoaded = true;
		}
		return jpts_jetBProbabilityBJetTag_;
	}
	vector<float> &jpts_jetProbabilityBJetTag()
	{
		if (not jpts_jetProbabilityBJetTag_isLoaded) {
			if (jpts_jetProbabilityBJetTag_branch != 0) {
				jpts_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_jetProbabilityBJetTag_isLoaded = true;
		}
		return jpts_jetProbabilityBJetTag_;
	}
	vector<float> &jpts_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not jpts_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (jpts_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				jpts_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_simpleSecondaryVertexHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return jpts_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &jpts_simpleSecondaryVertexHighPurBJetTags()
	{
		if (not jpts_simpleSecondaryVertexHighPurBJetTags_isLoaded) {
			if (jpts_simpleSecondaryVertexHighPurBJetTags_branch != 0) {
				jpts_simpleSecondaryVertexHighPurBJetTags_branch->GetEntry(index);
			} else { 
				printf("branch jpts_simpleSecondaryVertexHighPurBJetTags_branch does not exist!\n");
				exit(1);
			}
			jpts_simpleSecondaryVertexHighPurBJetTags_isLoaded = true;
		}
		return jpts_simpleSecondaryVertexHighPurBJetTags_;
	}
	vector<float> &jpts_softElectronByIP3dBJetTag()
	{
		if (not jpts_softElectronByIP3dBJetTag_isLoaded) {
			if (jpts_softElectronByIP3dBJetTag_branch != 0) {
				jpts_softElectronByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_softElectronByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_softElectronByIP3dBJetTag_isLoaded = true;
		}
		return jpts_softElectronByIP3dBJetTag_;
	}
	vector<float> &jpts_softElectronByPtBJetTag()
	{
		if (not jpts_softElectronByPtBJetTag_isLoaded) {
			if (jpts_softElectronByPtBJetTag_branch != 0) {
				jpts_softElectronByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_softElectronByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_softElectronByPtBJetTag_isLoaded = true;
		}
		return jpts_softElectronByPtBJetTag_;
	}
	vector<float> &jpts_softMuonBJetTag()
	{
		if (not jpts_softMuonBJetTag_isLoaded) {
			if (jpts_softMuonBJetTag_branch != 0) {
				jpts_softMuonBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_softMuonBJetTag_isLoaded = true;
		}
		return jpts_softMuonBJetTag_;
	}
	vector<float> &jpts_softMuonByIP3dBJetTag()
	{
		if (not jpts_softMuonByIP3dBJetTag_isLoaded) {
			if (jpts_softMuonByIP3dBJetTag_branch != 0) {
				jpts_softMuonByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_softMuonByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_softMuonByIP3dBJetTag_isLoaded = true;
		}
		return jpts_softMuonByIP3dBJetTag_;
	}
	vector<float> &jpts_softMuonByPtBJetTag()
	{
		if (not jpts_softMuonByPtBJetTag_isLoaded) {
			if (jpts_softMuonByPtBJetTag_branch != 0) {
				jpts_softMuonByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_softMuonByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_softMuonByPtBJetTag_isLoaded = true;
		}
		return jpts_softMuonByPtBJetTag_;
	}
	vector<float> &jpts_trackCountingHighEffBJetTag()
	{
		if (not jpts_trackCountingHighEffBJetTag_isLoaded) {
			if (jpts_trackCountingHighEffBJetTag_branch != 0) {
				jpts_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_trackCountingHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return jpts_trackCountingHighEffBJetTag_;
	}
	vector<float> &jpts_trackCountingHighPurBJetTag()
	{
		if (not jpts_trackCountingHighPurBJetTag_isLoaded) {
			if (jpts_trackCountingHighPurBJetTag_branch != 0) {
				jpts_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jpts_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jpts_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return jpts_trackCountingHighPurBJetTag_;
	}
	vector<float> &jets_combinedSecondaryVertexBJetTag()
	{
		if (not jets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (jets_combinedSecondaryVertexBJetTag_branch != 0) {
				jets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return jets_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &jets_combinedSecondaryVertexMVABJetTag()
	{
		if (not jets_combinedSecondaryVertexMVABJetTag_isLoaded) {
			if (jets_combinedSecondaryVertexMVABJetTag_branch != 0) {
				jets_combinedSecondaryVertexMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return jets_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &jets_jetBProbabilityBJetTag()
	{
		if (not jets_jetBProbabilityBJetTag_isLoaded) {
			if (jets_jetBProbabilityBJetTag_branch != 0) {
				jets_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_jetBProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_jetBProbabilityBJetTag_isLoaded = true;
		}
		return jets_jetBProbabilityBJetTag_;
	}
	vector<float> &jets_jetProbabilityBJetTag()
	{
		if (not jets_jetProbabilityBJetTag_isLoaded) {
			if (jets_jetProbabilityBJetTag_branch != 0) {
				jets_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_jetProbabilityBJetTag_isLoaded = true;
		}
		return jets_jetProbabilityBJetTag_;
	}
	vector<float> &jets_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not jets_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (jets_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				jets_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_simpleSecondaryVertexHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return jets_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &jets_simpleSecondaryVertexHighPurBJetTags()
	{
		if (not jets_simpleSecondaryVertexHighPurBJetTags_isLoaded) {
			if (jets_simpleSecondaryVertexHighPurBJetTags_branch != 0) {
				jets_simpleSecondaryVertexHighPurBJetTags_branch->GetEntry(index);
			} else { 
				printf("branch jets_simpleSecondaryVertexHighPurBJetTags_branch does not exist!\n");
				exit(1);
			}
			jets_simpleSecondaryVertexHighPurBJetTags_isLoaded = true;
		}
		return jets_simpleSecondaryVertexHighPurBJetTags_;
	}
	vector<float> &jets_softElectronByIP3dBJetTag()
	{
		if (not jets_softElectronByIP3dBJetTag_isLoaded) {
			if (jets_softElectronByIP3dBJetTag_branch != 0) {
				jets_softElectronByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_softElectronByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softElectronByIP3dBJetTag_isLoaded = true;
		}
		return jets_softElectronByIP3dBJetTag_;
	}
	vector<float> &jets_softElectronByPtBJetTag()
	{
		if (not jets_softElectronByPtBJetTag_isLoaded) {
			if (jets_softElectronByPtBJetTag_branch != 0) {
				jets_softElectronByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_softElectronByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softElectronByPtBJetTag_isLoaded = true;
		}
		return jets_softElectronByPtBJetTag_;
	}
	vector<float> &jets_softMuonBJetTag()
	{
		if (not jets_softMuonBJetTag_isLoaded) {
			if (jets_softMuonBJetTag_branch != 0) {
				jets_softMuonBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softMuonBJetTag_isLoaded = true;
		}
		return jets_softMuonBJetTag_;
	}
	vector<float> &jets_softMuonByIP3dBJetTag()
	{
		if (not jets_softMuonByIP3dBJetTag_isLoaded) {
			if (jets_softMuonByIP3dBJetTag_branch != 0) {
				jets_softMuonByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_softMuonByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softMuonByIP3dBJetTag_isLoaded = true;
		}
		return jets_softMuonByIP3dBJetTag_;
	}
	vector<float> &jets_softMuonByPtBJetTag()
	{
		if (not jets_softMuonByPtBJetTag_isLoaded) {
			if (jets_softMuonByPtBJetTag_branch != 0) {
				jets_softMuonByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_softMuonByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softMuonByPtBJetTag_isLoaded = true;
		}
		return jets_softMuonByPtBJetTag_;
	}
	vector<float> &jets_trackCountingHighEffBJetTag()
	{
		if (not jets_trackCountingHighEffBJetTag_isLoaded) {
			if (jets_trackCountingHighEffBJetTag_branch != 0) {
				jets_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_trackCountingHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return jets_trackCountingHighEffBJetTag_;
	}
	vector<float> &jets_trackCountingHighPurBJetTag()
	{
		if (not jets_trackCountingHighPurBJetTag_isLoaded) {
			if (jets_trackCountingHighPurBJetTag_branch != 0) {
				jets_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return jets_trackCountingHighPurBJetTag_;
	}
	vector<float> &pfjets_combinedSecondaryVertexBJetTag()
	{
		if (not pfjets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (pfjets_combinedSecondaryVertexBJetTag_branch != 0) {
				pfjets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return pfjets_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &pfjets_combinedSecondaryVertexMVABJetTag()
	{
		if (not pfjets_combinedSecondaryVertexMVABJetTag_isLoaded) {
			if (pfjets_combinedSecondaryVertexMVABJetTag_branch != 0) {
				pfjets_combinedSecondaryVertexMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return pfjets_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &pfjets_jetBProbabilityBJetTag()
	{
		if (not pfjets_jetBProbabilityBJetTag_isLoaded) {
			if (pfjets_jetBProbabilityBJetTag_branch != 0) {
				pfjets_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_jetBProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_jetBProbabilityBJetTag_isLoaded = true;
		}
		return pfjets_jetBProbabilityBJetTag_;
	}
	vector<float> &pfjets_jetProbabilityBJetTag()
	{
		if (not pfjets_jetProbabilityBJetTag_isLoaded) {
			if (pfjets_jetProbabilityBJetTag_branch != 0) {
				pfjets_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_jetProbabilityBJetTag_isLoaded = true;
		}
		return pfjets_jetProbabilityBJetTag_;
	}
	vector<float> &pfjets_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (pfjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				pfjets_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_simpleSecondaryVertexHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return pfjets_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &pfjets_simpleSecondaryVertexHighPurBJetTags()
	{
		if (not pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded) {
			if (pfjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) {
				pfjets_simpleSecondaryVertexHighPurBJetTags_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_simpleSecondaryVertexHighPurBJetTags_branch does not exist!\n");
				exit(1);
			}
			pfjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = true;
		}
		return pfjets_simpleSecondaryVertexHighPurBJetTags_;
	}
	vector<float> &pfjets_softElectronByIP3dBJetTag()
	{
		if (not pfjets_softElectronByIP3dBJetTag_isLoaded) {
			if (pfjets_softElectronByIP3dBJetTag_branch != 0) {
				pfjets_softElectronByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_softElectronByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_softElectronByIP3dBJetTag_isLoaded = true;
		}
		return pfjets_softElectronByIP3dBJetTag_;
	}
	vector<float> &pfjets_softElectronByPtBJetTag()
	{
		if (not pfjets_softElectronByPtBJetTag_isLoaded) {
			if (pfjets_softElectronByPtBJetTag_branch != 0) {
				pfjets_softElectronByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_softElectronByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_softElectronByPtBJetTag_isLoaded = true;
		}
		return pfjets_softElectronByPtBJetTag_;
	}
	vector<float> &pfjets_softMuonBJetTag()
	{
		if (not pfjets_softMuonBJetTag_isLoaded) {
			if (pfjets_softMuonBJetTag_branch != 0) {
				pfjets_softMuonBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_softMuonBJetTag_isLoaded = true;
		}
		return pfjets_softMuonBJetTag_;
	}
	vector<float> &pfjets_softMuonByIP3dBJetTag()
	{
		if (not pfjets_softMuonByIP3dBJetTag_isLoaded) {
			if (pfjets_softMuonByIP3dBJetTag_branch != 0) {
				pfjets_softMuonByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_softMuonByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_softMuonByIP3dBJetTag_isLoaded = true;
		}
		return pfjets_softMuonByIP3dBJetTag_;
	}
	vector<float> &pfjets_softMuonByPtBJetTag()
	{
		if (not pfjets_softMuonByPtBJetTag_isLoaded) {
			if (pfjets_softMuonByPtBJetTag_branch != 0) {
				pfjets_softMuonByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_softMuonByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_softMuonByPtBJetTag_isLoaded = true;
		}
		return pfjets_softMuonByPtBJetTag_;
	}
	vector<float> &pfjets_trackCountingHighEffBJetTag()
	{
		if (not pfjets_trackCountingHighEffBJetTag_isLoaded) {
			if (pfjets_trackCountingHighEffBJetTag_branch != 0) {
				pfjets_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_trackCountingHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return pfjets_trackCountingHighEffBJetTag_;
	}
	vector<float> &pfjets_trackCountingHighPurBJetTag()
	{
		if (not pfjets_trackCountingHighPurBJetTag_isLoaded) {
			if (pfjets_trackCountingHighPurBJetTag_branch != 0) {
				pfjets_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			pfjets_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return pfjets_trackCountingHighPurBJetTag_;
	}
	vector<float> &trkjets_combinedSecondaryVertexBJetTag()
	{
		if (not trkjets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (trkjets_combinedSecondaryVertexBJetTag_branch != 0) {
				trkjets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return trkjets_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &trkjets_combinedSecondaryVertexMVABJetTag()
	{
		if (not trkjets_combinedSecondaryVertexMVABJetTag_isLoaded) {
			if (trkjets_combinedSecondaryVertexMVABJetTag_branch != 0) {
				trkjets_combinedSecondaryVertexMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return trkjets_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &trkjets_jetBProbabilityBJetTag()
	{
		if (not trkjets_jetBProbabilityBJetTag_isLoaded) {
			if (trkjets_jetBProbabilityBJetTag_branch != 0) {
				trkjets_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_jetBProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_jetBProbabilityBJetTag_isLoaded = true;
		}
		return trkjets_jetBProbabilityBJetTag_;
	}
	vector<float> &trkjets_jetProbabilityBJetTag()
	{
		if (not trkjets_jetProbabilityBJetTag_isLoaded) {
			if (trkjets_jetProbabilityBJetTag_branch != 0) {
				trkjets_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_jetProbabilityBJetTag_isLoaded = true;
		}
		return trkjets_jetProbabilityBJetTag_;
	}
	vector<float> &trkjets_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not trkjets_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (trkjets_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				trkjets_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_simpleSecondaryVertexHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return trkjets_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &trkjets_simpleSecondaryVertexHighPurBJetTags()
	{
		if (not trkjets_simpleSecondaryVertexHighPurBJetTags_isLoaded) {
			if (trkjets_simpleSecondaryVertexHighPurBJetTags_branch != 0) {
				trkjets_simpleSecondaryVertexHighPurBJetTags_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_simpleSecondaryVertexHighPurBJetTags_branch does not exist!\n");
				exit(1);
			}
			trkjets_simpleSecondaryVertexHighPurBJetTags_isLoaded = true;
		}
		return trkjets_simpleSecondaryVertexHighPurBJetTags_;
	}
	vector<float> &trkjets_softElectronByIP3dBJetTag()
	{
		if (not trkjets_softElectronByIP3dBJetTag_isLoaded) {
			if (trkjets_softElectronByIP3dBJetTag_branch != 0) {
				trkjets_softElectronByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_softElectronByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_softElectronByIP3dBJetTag_isLoaded = true;
		}
		return trkjets_softElectronByIP3dBJetTag_;
	}
	vector<float> &trkjets_softElectronByPtBJetTag()
	{
		if (not trkjets_softElectronByPtBJetTag_isLoaded) {
			if (trkjets_softElectronByPtBJetTag_branch != 0) {
				trkjets_softElectronByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_softElectronByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_softElectronByPtBJetTag_isLoaded = true;
		}
		return trkjets_softElectronByPtBJetTag_;
	}
	vector<float> &trkjets_softMuonBJetTag()
	{
		if (not trkjets_softMuonBJetTag_isLoaded) {
			if (trkjets_softMuonBJetTag_branch != 0) {
				trkjets_softMuonBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_softMuonBJetTag_isLoaded = true;
		}
		return trkjets_softMuonBJetTag_;
	}
	vector<float> &trkjets_softMuonByIP3dBJetTag()
	{
		if (not trkjets_softMuonByIP3dBJetTag_isLoaded) {
			if (trkjets_softMuonByIP3dBJetTag_branch != 0) {
				trkjets_softMuonByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_softMuonByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_softMuonByIP3dBJetTag_isLoaded = true;
		}
		return trkjets_softMuonByIP3dBJetTag_;
	}
	vector<float> &trkjets_softMuonByPtBJetTag()
	{
		if (not trkjets_softMuonByPtBJetTag_isLoaded) {
			if (trkjets_softMuonByPtBJetTag_branch != 0) {
				trkjets_softMuonByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_softMuonByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_softMuonByPtBJetTag_isLoaded = true;
		}
		return trkjets_softMuonByPtBJetTag_;
	}
	vector<float> &trkjets_trackCountingHighEffBJetTag()
	{
		if (not trkjets_trackCountingHighEffBJetTag_isLoaded) {
			if (trkjets_trackCountingHighEffBJetTag_branch != 0) {
				trkjets_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_trackCountingHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return trkjets_trackCountingHighEffBJetTag_;
	}
	vector<float> &trkjets_trackCountingHighPurBJetTag()
	{
		if (not trkjets_trackCountingHighPurBJetTag_isLoaded) {
			if (trkjets_trackCountingHighPurBJetTag_branch != 0) {
				trkjets_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			trkjets_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return trkjets_trackCountingHighPurBJetTag_;
	}
	vector<float> &evt_bs_covMatrix()
	{
		if (not evt_bs_covMatrix_isLoaded) {
			if (evt_bs_covMatrix_branch != 0) {
				evt_bs_covMatrix_branch->GetEntry(index);
			} else { 
				printf("branch evt_bs_covMatrix_branch does not exist!\n");
				exit(1);
			}
			evt_bs_covMatrix_isLoaded = true;
		}
		return evt_bs_covMatrix_;
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
	vector<float> &jets_mc3dr()
	{
		if (not jets_mc3dr_isLoaded) {
			if (jets_mc3dr_branch != 0) {
				jets_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc3dr_branch does not exist!\n");
				exit(1);
			}
			jets_mc3dr_isLoaded = true;
		}
		return jets_mc3dr_;
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
	vector<float> &pfjets_mc3dr()
	{
		if (not pfjets_mc3dr_isLoaded) {
			if (pfjets_mc3dr_branch != 0) {
				pfjets_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc3dr_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc3dr_isLoaded = true;
		}
		return pfjets_mc3dr_;
	}
	vector<float> &pfjets_mcdr()
	{
		if (not pfjets_mcdr_isLoaded) {
			if (pfjets_mcdr_branch != 0) {
				pfjets_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mcdr_branch does not exist!\n");
				exit(1);
			}
			pfjets_mcdr_isLoaded = true;
		}
		return pfjets_mcdr_;
	}
	vector<float> &pfjets_mc_emEnergy()
	{
		if (not pfjets_mc_emEnergy_isLoaded) {
			if (pfjets_mc_emEnergy_branch != 0) {
				pfjets_mc_emEnergy_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_emEnergy_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_emEnergy_isLoaded = true;
		}
		return pfjets_mc_emEnergy_;
	}
	vector<float> &pfjets_mc_gpdr()
	{
		if (not pfjets_mc_gpdr_isLoaded) {
			if (pfjets_mc_gpdr_branch != 0) {
				pfjets_mc_gpdr_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_gpdr_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_gpdr_isLoaded = true;
		}
		return pfjets_mc_gpdr_;
	}
	vector<float> &pfjets_mc_hadEnergy()
	{
		if (not pfjets_mc_hadEnergy_isLoaded) {
			if (pfjets_mc_hadEnergy_branch != 0) {
				pfjets_mc_hadEnergy_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_hadEnergy_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_hadEnergy_isLoaded = true;
		}
		return pfjets_mc_hadEnergy_;
	}
	vector<float> &pfjets_mc_invEnergy()
	{
		if (not pfjets_mc_invEnergy_isLoaded) {
			if (pfjets_mc_invEnergy_branch != 0) {
				pfjets_mc_invEnergy_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_invEnergy_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_invEnergy_isLoaded = true;
		}
		return pfjets_mc_invEnergy_;
	}
	vector<float> &pfjets_mc_otherEnergy()
	{
		if (not pfjets_mc_otherEnergy_isLoaded) {
			if (pfjets_mc_otherEnergy_branch != 0) {
				pfjets_mc_otherEnergy_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_otherEnergy_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_otherEnergy_isLoaded = true;
		}
		return pfjets_mc_otherEnergy_;
	}
	vector<float> &photons_mc3dr()
	{
		if (not photons_mc3dr_isLoaded) {
			if (photons_mc3dr_branch != 0) {
				photons_mc3dr_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc3dr_branch does not exist!\n");
				exit(1);
			}
			photons_mc3dr_isLoaded = true;
		}
		return photons_mc3dr_;
	}
	vector<float> &photons_mcdr()
	{
		if (not photons_mcdr_isLoaded) {
			if (photons_mcdr_branch != 0) {
				photons_mcdr_branch->GetEntry(index);
			} else { 
				printf("branch photons_mcdr_branch does not exist!\n");
				exit(1);
			}
			photons_mcdr_isLoaded = true;
		}
		return photons_mcdr_;
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
	vector<float> &els_ecalJuraIso()
	{
		if (not els_ecalJuraIso_isLoaded) {
			if (els_ecalJuraIso_branch != 0) {
				els_ecalJuraIso_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalJuraIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalJuraIso_isLoaded = true;
		}
		return els_ecalJuraIso_;
	}
	vector<float> &els_ecalJuraTowerIso()
	{
		if (not els_ecalJuraTowerIso_isLoaded) {
			if (els_ecalJuraTowerIso_branch != 0) {
				els_ecalJuraTowerIso_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalJuraTowerIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalJuraTowerIso_isLoaded = true;
		}
		return els_ecalJuraTowerIso_;
	}
	vector<float> &els_hcalConeIso()
	{
		if (not els_hcalConeIso_isLoaded) {
			if (els_hcalConeIso_branch != 0) {
				els_hcalConeIso_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalConeIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalConeIso_isLoaded = true;
		}
		return els_hcalConeIso_;
	}
	vector<float> &els_tkJuraIso()
	{
		if (not els_tkJuraIso_isLoaded) {
			if (els_tkJuraIso_branch != 0) {
				els_tkJuraIso_branch->GetEntry(index);
			} else { 
				printf("branch els_tkJuraIso_branch does not exist!\n");
				exit(1);
			}
			els_tkJuraIso_isLoaded = true;
		}
		return els_tkJuraIso_;
	}
	vector<float> &els_jetdr()
	{
		if (not els_jetdr_isLoaded) {
			if (els_jetdr_branch != 0) {
				els_jetdr_branch->GetEntry(index);
			} else { 
				printf("branch els_jetdr_branch does not exist!\n");
				exit(1);
			}
			els_jetdr_isLoaded = true;
		}
		return els_jetdr_;
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
	vector<float> &els_conv_radius()
	{
		if (not els_conv_radius_isLoaded) {
			if (els_conv_radius_branch != 0) {
				els_conv_radius_branch->GetEntry(index);
			} else { 
				printf("branch els_conv_radius_branch does not exist!\n");
				exit(1);
			}
			els_conv_radius_isLoaded = true;
		}
		return els_conv_radius_;
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
	vector<float> &els_deltaEtaEleClusterTrackAtCalo()
	{
		if (not els_deltaEtaEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaEtaEleClusterTrackAtCalo_branch != 0) {
				els_deltaEtaEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { 
				printf("branch els_deltaEtaEleClusterTrackAtCalo_branch does not exist!\n");
				exit(1);
			}
			els_deltaEtaEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaEtaEleClusterTrackAtCalo_;
	}
	vector<float> &els_deltaPhiEleClusterTrackAtCalo()
	{
		if (not els_deltaPhiEleClusterTrackAtCalo_isLoaded) {
			if (els_deltaPhiEleClusterTrackAtCalo_branch != 0) {
				els_deltaPhiEleClusterTrackAtCalo_branch->GetEntry(index);
			} else { 
				printf("branch els_deltaPhiEleClusterTrackAtCalo_branch does not exist!\n");
				exit(1);
			}
			els_deltaPhiEleClusterTrackAtCalo_isLoaded = true;
		}
		return els_deltaPhiEleClusterTrackAtCalo_;
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
	vector<float> &els_eOverPOut()
	{
		if (not els_eOverPOut_isLoaded) {
			if (els_eOverPOut_branch != 0) {
				els_eOverPOut_branch->GetEntry(index);
			} else { 
				printf("branch els_eOverPOut_branch does not exist!\n");
				exit(1);
			}
			els_eOverPOut_isLoaded = true;
		}
		return els_eOverPOut_;
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
	vector<float> &els_ecalEnergy()
	{
		if (not els_ecalEnergy_isLoaded) {
			if (els_ecalEnergy_branch != 0) {
				els_ecalEnergy_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalEnergy_branch does not exist!\n");
				exit(1);
			}
			els_ecalEnergy_isLoaded = true;
		}
		return els_ecalEnergy_;
	}
	vector<float> &els_ecalEnergyError()
	{
		if (not els_ecalEnergyError_isLoaded) {
			if (els_ecalEnergyError_branch != 0) {
				els_ecalEnergyError_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalEnergyError_branch does not exist!\n");
				exit(1);
			}
			els_ecalEnergyError_isLoaded = true;
		}
		return els_ecalEnergyError_;
	}
	vector<float> &els_ecalIso()
	{
		if (not els_ecalIso_isLoaded) {
			if (els_ecalIso_branch != 0) {
				els_ecalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_ecalIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso_isLoaded = true;
		}
		return els_ecalIso_;
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
	vector<float> &els_electronMomentumError()
	{
		if (not els_electronMomentumError_isLoaded) {
			if (els_electronMomentumError_branch != 0) {
				els_electronMomentumError_branch->GetEntry(index);
			} else { 
				printf("branch els_electronMomentumError_branch does not exist!\n");
				exit(1);
			}
			els_electronMomentumError_isLoaded = true;
		}
		return els_electronMomentumError_;
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
	vector<float> &els_fbrem()
	{
		if (not els_fbrem_isLoaded) {
			if (els_fbrem_branch != 0) {
				els_fbrem_branch->GetEntry(index);
			} else { 
				printf("branch els_fbrem_branch does not exist!\n");
				exit(1);
			}
			els_fbrem_isLoaded = true;
		}
		return els_fbrem_;
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
	vector<float> &els_hcalDepth1OverEcal()
	{
		if (not els_hcalDepth1OverEcal_isLoaded) {
			if (els_hcalDepth1OverEcal_branch != 0) {
				els_hcalDepth1OverEcal_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1OverEcal_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1OverEcal_isLoaded = true;
		}
		return els_hcalDepth1OverEcal_;
	}
	vector<float> &els_hcalDepth1TowerSumEt()
	{
		if (not els_hcalDepth1TowerSumEt_isLoaded) {
			if (els_hcalDepth1TowerSumEt_branch != 0) {
				els_hcalDepth1TowerSumEt_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1TowerSumEt_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt_;
	}
	vector<float> &els_hcalDepth1TowerSumEt04()
	{
		if (not els_hcalDepth1TowerSumEt04_isLoaded) {
			if (els_hcalDepth1TowerSumEt04_branch != 0) {
				els_hcalDepth1TowerSumEt04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth1TowerSumEt04_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth1TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth1TowerSumEt04_;
	}
	vector<float> &els_hcalDepth2OverEcal()
	{
		if (not els_hcalDepth2OverEcal_isLoaded) {
			if (els_hcalDepth2OverEcal_branch != 0) {
				els_hcalDepth2OverEcal_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2OverEcal_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2OverEcal_isLoaded = true;
		}
		return els_hcalDepth2OverEcal_;
	}
	vector<float> &els_hcalDepth2TowerSumEt()
	{
		if (not els_hcalDepth2TowerSumEt_isLoaded) {
			if (els_hcalDepth2TowerSumEt_branch != 0) {
				els_hcalDepth2TowerSumEt_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2TowerSumEt_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2TowerSumEt_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt_;
	}
	vector<float> &els_hcalDepth2TowerSumEt04()
	{
		if (not els_hcalDepth2TowerSumEt04_isLoaded) {
			if (els_hcalDepth2TowerSumEt04_branch != 0) {
				els_hcalDepth2TowerSumEt04_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalDepth2TowerSumEt04_branch does not exist!\n");
				exit(1);
			}
			els_hcalDepth2TowerSumEt04_isLoaded = true;
		}
		return els_hcalDepth2TowerSumEt04_;
	}
	vector<float> &els_hcalIso()
	{
		if (not els_hcalIso_isLoaded) {
			if (els_hcalIso_branch != 0) {
				els_hcalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_hcalIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso_isLoaded = true;
		}
		return els_hcalIso_;
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
	vector<float> &els_mva()
	{
		if (not els_mva_isLoaded) {
			if (els_mva_branch != 0) {
				els_mva_branch->GetEntry(index);
			} else { 
				printf("branch els_mva_branch does not exist!\n");
				exit(1);
			}
			els_mva_isLoaded = true;
		}
		return els_mva_;
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
	vector<float> &els_sigmaIEtaIEtaSC()
	{
		if (not els_sigmaIEtaIEtaSC_isLoaded) {
			if (els_sigmaIEtaIEtaSC_branch != 0) {
				els_sigmaIEtaIEtaSC_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIEtaIEtaSC_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIEtaIEtaSC_isLoaded = true;
		}
		return els_sigmaIEtaIEtaSC_;
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
	vector<float> &els_sigmaIPhiIPhiSC()
	{
		if (not els_sigmaIPhiIPhiSC_isLoaded) {
			if (els_sigmaIPhiIPhiSC_branch != 0) {
				els_sigmaIPhiIPhiSC_branch->GetEntry(index);
			} else { 
				printf("branch els_sigmaIPhiIPhiSC_branch does not exist!\n");
				exit(1);
			}
			els_sigmaIPhiIPhiSC_isLoaded = true;
		}
		return els_sigmaIPhiIPhiSC_;
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
	vector<float> &els_tkIso()
	{
		if (not els_tkIso_isLoaded) {
			if (els_tkIso_branch != 0) {
				els_tkIso_branch->GetEntry(index);
			} else { 
				printf("branch els_tkIso_branch does not exist!\n");
				exit(1);
			}
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
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
	vector<float> &els_trackMomentumError()
	{
		if (not els_trackMomentumError_isLoaded) {
			if (els_trackMomentumError_branch != 0) {
				els_trackMomentumError_branch->GetEntry(index);
			} else { 
				printf("branch els_trackMomentumError_branch does not exist!\n");
				exit(1);
			}
			els_trackMomentumError_isLoaded = true;
		}
		return els_trackMomentumError_;
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
	vector<float> &gsftrks_chi2()
	{
		if (not gsftrks_chi2_isLoaded) {
			if (gsftrks_chi2_branch != 0) {
				gsftrks_chi2_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_chi2_branch does not exist!\n");
				exit(1);
			}
			gsftrks_chi2_isLoaded = true;
		}
		return gsftrks_chi2_;
	}
	vector<float> &gsftrks_d0()
	{
		if (not gsftrks_d0_isLoaded) {
			if (gsftrks_d0_branch != 0) {
				gsftrks_d0_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_d0_branch does not exist!\n");
				exit(1);
			}
			gsftrks_d0_isLoaded = true;
		}
		return gsftrks_d0_;
	}
	vector<float> &gsftrks_d0Err()
	{
		if (not gsftrks_d0Err_isLoaded) {
			if (gsftrks_d0Err_branch != 0) {
				gsftrks_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_d0Err_branch does not exist!\n");
				exit(1);
			}
			gsftrks_d0Err_isLoaded = true;
		}
		return gsftrks_d0Err_;
	}
	vector<float> &gsftrks_d0corr()
	{
		if (not gsftrks_d0corr_isLoaded) {
			if (gsftrks_d0corr_branch != 0) {
				gsftrks_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_d0corr_branch does not exist!\n");
				exit(1);
			}
			gsftrks_d0corr_isLoaded = true;
		}
		return gsftrks_d0corr_;
	}
	vector<float> &gsftrks_d0corrPhi()
	{
		if (not gsftrks_d0corrPhi_isLoaded) {
			if (gsftrks_d0corrPhi_branch != 0) {
				gsftrks_d0corrPhi_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_d0corrPhi_branch does not exist!\n");
				exit(1);
			}
			gsftrks_d0corrPhi_isLoaded = true;
		}
		return gsftrks_d0corrPhi_;
	}
	vector<float> &gsftrks_d0phiCov()
	{
		if (not gsftrks_d0phiCov_isLoaded) {
			if (gsftrks_d0phiCov_branch != 0) {
				gsftrks_d0phiCov_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_d0phiCov_branch does not exist!\n");
				exit(1);
			}
			gsftrks_d0phiCov_isLoaded = true;
		}
		return gsftrks_d0phiCov_;
	}
	vector<float> &gsftrks_etaErr()
	{
		if (not gsftrks_etaErr_isLoaded) {
			if (gsftrks_etaErr_branch != 0) {
				gsftrks_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_etaErr_branch does not exist!\n");
				exit(1);
			}
			gsftrks_etaErr_isLoaded = true;
		}
		return gsftrks_etaErr_;
	}
	vector<float> &gsftrks_layer1_charge()
	{
		if (not gsftrks_layer1_charge_isLoaded) {
			if (gsftrks_layer1_charge_branch != 0) {
				gsftrks_layer1_charge_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_layer1_charge_branch does not exist!\n");
				exit(1);
			}
			gsftrks_layer1_charge_isLoaded = true;
		}
		return gsftrks_layer1_charge_;
	}
	vector<float> &gsftrks_ndof()
	{
		if (not gsftrks_ndof_isLoaded) {
			if (gsftrks_ndof_branch != 0) {
				gsftrks_ndof_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_ndof_branch does not exist!\n");
				exit(1);
			}
			gsftrks_ndof_isLoaded = true;
		}
		return gsftrks_ndof_;
	}
	vector<float> &gsftrks_phiErr()
	{
		if (not gsftrks_phiErr_isLoaded) {
			if (gsftrks_phiErr_branch != 0) {
				gsftrks_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_phiErr_branch does not exist!\n");
				exit(1);
			}
			gsftrks_phiErr_isLoaded = true;
		}
		return gsftrks_phiErr_;
	}
	vector<float> &gsftrks_ptErr()
	{
		if (not gsftrks_ptErr_isLoaded) {
			if (gsftrks_ptErr_branch != 0) {
				gsftrks_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_ptErr_branch does not exist!\n");
				exit(1);
			}
			gsftrks_ptErr_isLoaded = true;
		}
		return gsftrks_ptErr_;
	}
	vector<float> &gsftrks_z0()
	{
		if (not gsftrks_z0_isLoaded) {
			if (gsftrks_z0_branch != 0) {
				gsftrks_z0_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_z0_branch does not exist!\n");
				exit(1);
			}
			gsftrks_z0_isLoaded = true;
		}
		return gsftrks_z0_;
	}
	vector<float> &gsftrks_z0Err()
	{
		if (not gsftrks_z0Err_isLoaded) {
			if (gsftrks_z0Err_branch != 0) {
				gsftrks_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_z0Err_branch does not exist!\n");
				exit(1);
			}
			gsftrks_z0Err_isLoaded = true;
		}
		return gsftrks_z0Err_;
	}
	vector<float> &gsftrks_z0corr()
	{
		if (not gsftrks_z0corr_isLoaded) {
			if (gsftrks_z0corr_branch != 0) {
				gsftrks_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_z0corr_branch does not exist!\n");
				exit(1);
			}
			gsftrks_z0corr_isLoaded = true;
		}
		return gsftrks_z0corr_;
	}
	vector<float> &hyp_Ht()
	{
		if (not hyp_Ht_isLoaded) {
			if (hyp_Ht_branch != 0) {
				hyp_Ht_branch->GetEntry(index);
			} else { 
				printf("branch hyp_Ht_branch does not exist!\n");
				exit(1);
			}
			hyp_Ht_isLoaded = true;
		}
		return hyp_Ht_;
	}
	vector<float> &hyp_dPhi_nJet_metMuonJESCorr()
	{
		if (not hyp_dPhi_nJet_metMuonJESCorr_isLoaded) {
			if (hyp_dPhi_nJet_metMuonJESCorr_branch != 0) {
				hyp_dPhi_nJet_metMuonJESCorr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_dPhi_nJet_metMuonJESCorr_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_metMuonJESCorr_isLoaded = true;
		}
		return hyp_dPhi_nJet_metMuonJESCorr_;
	}
	vector<float> &hyp_dPhi_nJet_muCorrMet()
	{
		if (not hyp_dPhi_nJet_muCorrMet_isLoaded) {
			if (hyp_dPhi_nJet_muCorrMet_branch != 0) {
				hyp_dPhi_nJet_muCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_dPhi_nJet_muCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_muCorrMet_isLoaded = true;
		}
		return hyp_dPhi_nJet_muCorrMet_;
	}
	vector<float> &hyp_dPhi_nJet_tcMet()
	{
		if (not hyp_dPhi_nJet_tcMet_isLoaded) {
			if (hyp_dPhi_nJet_tcMet_branch != 0) {
				hyp_dPhi_nJet_tcMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_dPhi_nJet_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_tcMet_isLoaded = true;
		}
		return hyp_dPhi_nJet_tcMet_;
	}
	vector<float> &hyp_dPhi_nJet_unCorrMet()
	{
		if (not hyp_dPhi_nJet_unCorrMet_isLoaded) {
			if (hyp_dPhi_nJet_unCorrMet_branch != 0) {
				hyp_dPhi_nJet_unCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_dPhi_nJet_unCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_unCorrMet_isLoaded = true;
		}
		return hyp_dPhi_nJet_unCorrMet_;
	}
	vector<float> &hyp_ll_chi2()
	{
		if (not hyp_ll_chi2_isLoaded) {
			if (hyp_ll_chi2_branch != 0) {
				hyp_ll_chi2_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_chi2_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_chi2_isLoaded = true;
		}
		return hyp_ll_chi2_;
	}
	vector<float> &hyp_ll_d0()
	{
		if (not hyp_ll_d0_isLoaded) {
			if (hyp_ll_d0_branch != 0) {
				hyp_ll_d0_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_d0_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_d0_isLoaded = true;
		}
		return hyp_ll_d0_;
	}
	vector<float> &hyp_ll_d0Err()
	{
		if (not hyp_ll_d0Err_isLoaded) {
			if (hyp_ll_d0Err_branch != 0) {
				hyp_ll_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_d0Err_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_d0Err_isLoaded = true;
		}
		return hyp_ll_d0Err_;
	}
	vector<float> &hyp_ll_d0corr()
	{
		if (not hyp_ll_d0corr_isLoaded) {
			if (hyp_ll_d0corr_branch != 0) {
				hyp_ll_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_d0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_d0corr_isLoaded = true;
		}
		return hyp_ll_d0corr_;
	}
	vector<float> &hyp_ll_dPhi_metMuonJESCorr()
	{
		if (not hyp_ll_dPhi_metMuonJESCorr_isLoaded) {
			if (hyp_ll_dPhi_metMuonJESCorr_branch != 0) {
				hyp_ll_dPhi_metMuonJESCorr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_dPhi_metMuonJESCorr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_metMuonJESCorr_isLoaded = true;
		}
		return hyp_ll_dPhi_metMuonJESCorr_;
	}
	vector<float> &hyp_ll_dPhi_muCorrMet()
	{
		if (not hyp_ll_dPhi_muCorrMet_isLoaded) {
			if (hyp_ll_dPhi_muCorrMet_branch != 0) {
				hyp_ll_dPhi_muCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_dPhi_muCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_muCorrMet_isLoaded = true;
		}
		return hyp_ll_dPhi_muCorrMet_;
	}
	vector<float> &hyp_ll_dPhi_tcMet()
	{
		if (not hyp_ll_dPhi_tcMet_isLoaded) {
			if (hyp_ll_dPhi_tcMet_branch != 0) {
				hyp_ll_dPhi_tcMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_dPhi_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_tcMet_isLoaded = true;
		}
		return hyp_ll_dPhi_tcMet_;
	}
	vector<float> &hyp_ll_dPhi_unCorrMet()
	{
		if (not hyp_ll_dPhi_unCorrMet_isLoaded) {
			if (hyp_ll_dPhi_unCorrMet_branch != 0) {
				hyp_ll_dPhi_unCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_dPhi_unCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_unCorrMet_isLoaded = true;
		}
		return hyp_ll_dPhi_unCorrMet_;
	}
	vector<float> &hyp_ll_etaErr()
	{
		if (not hyp_ll_etaErr_isLoaded) {
			if (hyp_ll_etaErr_branch != 0) {
				hyp_ll_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_etaErr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_etaErr_isLoaded = true;
		}
		return hyp_ll_etaErr_;
	}
	vector<float> &hyp_ll_ndof()
	{
		if (not hyp_ll_ndof_isLoaded) {
			if (hyp_ll_ndof_branch != 0) {
				hyp_ll_ndof_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_ndof_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_ndof_isLoaded = true;
		}
		return hyp_ll_ndof_;
	}
	vector<float> &hyp_ll_phiErr()
	{
		if (not hyp_ll_phiErr_isLoaded) {
			if (hyp_ll_phiErr_branch != 0) {
				hyp_ll_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_phiErr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_phiErr_isLoaded = true;
		}
		return hyp_ll_phiErr_;
	}
	vector<float> &hyp_ll_ptErr()
	{
		if (not hyp_ll_ptErr_isLoaded) {
			if (hyp_ll_ptErr_branch != 0) {
				hyp_ll_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_ptErr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_ptErr_isLoaded = true;
		}
		return hyp_ll_ptErr_;
	}
	vector<float> &hyp_ll_z0()
	{
		if (not hyp_ll_z0_isLoaded) {
			if (hyp_ll_z0_branch != 0) {
				hyp_ll_z0_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_z0_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_z0_isLoaded = true;
		}
		return hyp_ll_z0_;
	}
	vector<float> &hyp_ll_z0Err()
	{
		if (not hyp_ll_z0Err_isLoaded) {
			if (hyp_ll_z0Err_branch != 0) {
				hyp_ll_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_z0Err_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_z0Err_isLoaded = true;
		}
		return hyp_ll_z0Err_;
	}
	vector<float> &hyp_ll_z0corr()
	{
		if (not hyp_ll_z0corr_isLoaded) {
			if (hyp_ll_z0corr_branch != 0) {
				hyp_ll_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_z0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_z0corr_isLoaded = true;
		}
		return hyp_ll_z0corr_;
	}
	vector<float> &hyp_lt_chi2()
	{
		if (not hyp_lt_chi2_isLoaded) {
			if (hyp_lt_chi2_branch != 0) {
				hyp_lt_chi2_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_chi2_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_chi2_isLoaded = true;
		}
		return hyp_lt_chi2_;
	}
	vector<float> &hyp_lt_d0()
	{
		if (not hyp_lt_d0_isLoaded) {
			if (hyp_lt_d0_branch != 0) {
				hyp_lt_d0_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_d0_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_d0_isLoaded = true;
		}
		return hyp_lt_d0_;
	}
	vector<float> &hyp_lt_d0Err()
	{
		if (not hyp_lt_d0Err_isLoaded) {
			if (hyp_lt_d0Err_branch != 0) {
				hyp_lt_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_d0Err_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_d0Err_isLoaded = true;
		}
		return hyp_lt_d0Err_;
	}
	vector<float> &hyp_lt_d0corr()
	{
		if (not hyp_lt_d0corr_isLoaded) {
			if (hyp_lt_d0corr_branch != 0) {
				hyp_lt_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_d0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_d0corr_isLoaded = true;
		}
		return hyp_lt_d0corr_;
	}
	vector<float> &hyp_lt_dPhi_metMuonJESCorr()
	{
		if (not hyp_lt_dPhi_metMuonJESCorr_isLoaded) {
			if (hyp_lt_dPhi_metMuonJESCorr_branch != 0) {
				hyp_lt_dPhi_metMuonJESCorr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_dPhi_metMuonJESCorr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_metMuonJESCorr_isLoaded = true;
		}
		return hyp_lt_dPhi_metMuonJESCorr_;
	}
	vector<float> &hyp_lt_dPhi_muCorrMet()
	{
		if (not hyp_lt_dPhi_muCorrMet_isLoaded) {
			if (hyp_lt_dPhi_muCorrMet_branch != 0) {
				hyp_lt_dPhi_muCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_dPhi_muCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_muCorrMet_isLoaded = true;
		}
		return hyp_lt_dPhi_muCorrMet_;
	}
	vector<float> &hyp_lt_dPhi_tcMet()
	{
		if (not hyp_lt_dPhi_tcMet_isLoaded) {
			if (hyp_lt_dPhi_tcMet_branch != 0) {
				hyp_lt_dPhi_tcMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_dPhi_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_tcMet_isLoaded = true;
		}
		return hyp_lt_dPhi_tcMet_;
	}
	vector<float> &hyp_lt_dPhi_unCorrMet()
	{
		if (not hyp_lt_dPhi_unCorrMet_isLoaded) {
			if (hyp_lt_dPhi_unCorrMet_branch != 0) {
				hyp_lt_dPhi_unCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_dPhi_unCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_unCorrMet_isLoaded = true;
		}
		return hyp_lt_dPhi_unCorrMet_;
	}
	vector<float> &hyp_lt_etaErr()
	{
		if (not hyp_lt_etaErr_isLoaded) {
			if (hyp_lt_etaErr_branch != 0) {
				hyp_lt_etaErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_etaErr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_etaErr_isLoaded = true;
		}
		return hyp_lt_etaErr_;
	}
	vector<float> &hyp_lt_ndof()
	{
		if (not hyp_lt_ndof_isLoaded) {
			if (hyp_lt_ndof_branch != 0) {
				hyp_lt_ndof_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_ndof_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_ndof_isLoaded = true;
		}
		return hyp_lt_ndof_;
	}
	vector<float> &hyp_lt_phiErr()
	{
		if (not hyp_lt_phiErr_isLoaded) {
			if (hyp_lt_phiErr_branch != 0) {
				hyp_lt_phiErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_phiErr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_phiErr_isLoaded = true;
		}
		return hyp_lt_phiErr_;
	}
	vector<float> &hyp_lt_ptErr()
	{
		if (not hyp_lt_ptErr_isLoaded) {
			if (hyp_lt_ptErr_branch != 0) {
				hyp_lt_ptErr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_ptErr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_ptErr_isLoaded = true;
		}
		return hyp_lt_ptErr_;
	}
	vector<float> &hyp_lt_z0()
	{
		if (not hyp_lt_z0_isLoaded) {
			if (hyp_lt_z0_branch != 0) {
				hyp_lt_z0_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_z0_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_z0_isLoaded = true;
		}
		return hyp_lt_z0_;
	}
	vector<float> &hyp_lt_z0Err()
	{
		if (not hyp_lt_z0Err_isLoaded) {
			if (hyp_lt_z0Err_branch != 0) {
				hyp_lt_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_z0Err_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_z0Err_isLoaded = true;
		}
		return hyp_lt_z0Err_;
	}
	vector<float> &hyp_lt_z0corr()
	{
		if (not hyp_lt_z0corr_isLoaded) {
			if (hyp_lt_z0corr_branch != 0) {
				hyp_lt_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_z0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_z0corr_isLoaded = true;
		}
		return hyp_lt_z0corr_;
	}
	vector<float> &hyp_mt2_metMuonJESCorr()
	{
		if (not hyp_mt2_metMuonJESCorr_isLoaded) {
			if (hyp_mt2_metMuonJESCorr_branch != 0) {
				hyp_mt2_metMuonJESCorr_branch->GetEntry(index);
			} else { 
				printf("branch hyp_mt2_metMuonJESCorr_branch does not exist!\n");
				exit(1);
			}
			hyp_mt2_metMuonJESCorr_isLoaded = true;
		}
		return hyp_mt2_metMuonJESCorr_;
	}
	vector<float> &hyp_mt2_muCorrMet()
	{
		if (not hyp_mt2_muCorrMet_isLoaded) {
			if (hyp_mt2_muCorrMet_branch != 0) {
				hyp_mt2_muCorrMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_mt2_muCorrMet_branch does not exist!\n");
				exit(1);
			}
			hyp_mt2_muCorrMet_isLoaded = true;
		}
		return hyp_mt2_muCorrMet_;
	}
	vector<float> &hyp_mt2_tcMet()
	{
		if (not hyp_mt2_tcMet_isLoaded) {
			if (hyp_mt2_tcMet_branch != 0) {
				hyp_mt2_tcMet_branch->GetEntry(index);
			} else { 
				printf("branch hyp_mt2_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_mt2_tcMet_isLoaded = true;
		}
		return hyp_mt2_tcMet_;
	}
	vector<float> &hyp_sumJetPt()
	{
		if (not hyp_sumJetPt_isLoaded) {
			if (hyp_sumJetPt_branch != 0) {
				hyp_sumJetPt_branch->GetEntry(index);
			} else { 
				printf("branch hyp_sumJetPt_branch does not exist!\n");
				exit(1);
			}
			hyp_sumJetPt_isLoaded = true;
		}
		return hyp_sumJetPt_;
	}
	vector<float> &hyp_FVFit_chi2ndf()
	{
		if (not hyp_FVFit_chi2ndf_isLoaded) {
			if (hyp_FVFit_chi2ndf_branch != 0) {
				hyp_FVFit_chi2ndf_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_chi2ndf_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_chi2ndf_isLoaded = true;
		}
		return hyp_FVFit_chi2ndf_;
	}
	vector<float> &hyp_FVFit_prob()
	{
		if (not hyp_FVFit_prob_isLoaded) {
			if (hyp_FVFit_prob_branch != 0) {
				hyp_FVFit_prob_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_prob_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_prob_isLoaded = true;
		}
		return hyp_FVFit_prob_;
	}
	vector<float> &hyp_FVFit_v4cxx()
	{
		if (not hyp_FVFit_v4cxx_isLoaded) {
			if (hyp_FVFit_v4cxx_branch != 0) {
				hyp_FVFit_v4cxx_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4cxx_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4cxx_isLoaded = true;
		}
		return hyp_FVFit_v4cxx_;
	}
	vector<float> &hyp_FVFit_v4cxy()
	{
		if (not hyp_FVFit_v4cxy_isLoaded) {
			if (hyp_FVFit_v4cxy_branch != 0) {
				hyp_FVFit_v4cxy_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4cxy_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4cxy_isLoaded = true;
		}
		return hyp_FVFit_v4cxy_;
	}
	vector<float> &hyp_FVFit_v4cyy()
	{
		if (not hyp_FVFit_v4cyy_isLoaded) {
			if (hyp_FVFit_v4cyy_branch != 0) {
				hyp_FVFit_v4cyy_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4cyy_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4cyy_isLoaded = true;
		}
		return hyp_FVFit_v4cyy_;
	}
	vector<float> &hyp_FVFit_v4czx()
	{
		if (not hyp_FVFit_v4czx_isLoaded) {
			if (hyp_FVFit_v4czx_branch != 0) {
				hyp_FVFit_v4czx_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4czx_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4czx_isLoaded = true;
		}
		return hyp_FVFit_v4czx_;
	}
	vector<float> &hyp_FVFit_v4czy()
	{
		if (not hyp_FVFit_v4czy_isLoaded) {
			if (hyp_FVFit_v4czy_branch != 0) {
				hyp_FVFit_v4czy_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4czy_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4czy_isLoaded = true;
		}
		return hyp_FVFit_v4czy_;
	}
	vector<float> &hyp_FVFit_v4czz()
	{
		if (not hyp_FVFit_v4czz_isLoaded) {
			if (hyp_FVFit_v4czz_branch != 0) {
				hyp_FVFit_v4czz_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_v4czz_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_v4czz_isLoaded = true;
		}
		return hyp_FVFit_v4czz_;
	}
	vector<float> &hyp_ll_ecaliso()
	{
		if (not hyp_ll_ecaliso_isLoaded) {
			if (hyp_ll_ecaliso_branch != 0) {
				hyp_ll_ecaliso_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_ecaliso_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_ecaliso_isLoaded = true;
		}
		return hyp_ll_ecaliso_;
	}
	vector<float> &hyp_ll_trkiso()
	{
		if (not hyp_ll_trkiso_isLoaded) {
			if (hyp_ll_trkiso_branch != 0) {
				hyp_ll_trkiso_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_trkiso_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_trkiso_isLoaded = true;
		}
		return hyp_ll_trkiso_;
	}
	vector<float> &hyp_lt_ecaliso()
	{
		if (not hyp_lt_ecaliso_isLoaded) {
			if (hyp_lt_ecaliso_branch != 0) {
				hyp_lt_ecaliso_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_ecaliso_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_ecaliso_isLoaded = true;
		}
		return hyp_lt_ecaliso_;
	}
	vector<float> &hyp_lt_trkiso()
	{
		if (not hyp_lt_trkiso_isLoaded) {
			if (hyp_lt_trkiso_branch != 0) {
				hyp_lt_trkiso_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_trkiso_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_trkiso_isLoaded = true;
		}
		return hyp_lt_trkiso_;
	}
	vector<float> &jets_approximatefHPD()
	{
		if (not jets_approximatefHPD_isLoaded) {
			if (jets_approximatefHPD_branch != 0) {
				jets_approximatefHPD_branch->GetEntry(index);
			} else { 
				printf("branch jets_approximatefHPD_branch does not exist!\n");
				exit(1);
			}
			jets_approximatefHPD_isLoaded = true;
		}
		return jets_approximatefHPD_;
	}
	vector<float> &jets_approximatefRBX()
	{
		if (not jets_approximatefRBX_isLoaded) {
			if (jets_approximatefRBX_branch != 0) {
				jets_approximatefRBX_branch->GetEntry(index);
			} else { 
				printf("branch jets_approximatefRBX_branch does not exist!\n");
				exit(1);
			}
			jets_approximatefRBX_isLoaded = true;
		}
		return jets_approximatefRBX_;
	}
	vector<float> &jets_cor()
	{
		if (not jets_cor_isLoaded) {
			if (jets_cor_branch != 0) {
				jets_cor_branch->GetEntry(index);
			} else { 
				printf("branch jets_cor_branch does not exist!\n");
				exit(1);
			}
			jets_cor_isLoaded = true;
		}
		return jets_cor_;
	}
	vector<float> &jets_emFrac()
	{
		if (not jets_emFrac_isLoaded) {
			if (jets_emFrac_branch != 0) {
				jets_emFrac_branch->GetEntry(index);
			} else { 
				printf("branch jets_emFrac_branch does not exist!\n");
				exit(1);
			}
			jets_emFrac_isLoaded = true;
		}
		return jets_emFrac_;
	}
	vector<float> &jets_fHPD()
	{
		if (not jets_fHPD_isLoaded) {
			if (jets_fHPD_branch != 0) {
				jets_fHPD_branch->GetEntry(index);
			} else { 
				printf("branch jets_fHPD_branch does not exist!\n");
				exit(1);
			}
			jets_fHPD_isLoaded = true;
		}
		return jets_fHPD_;
	}
	vector<float> &jets_fRBX()
	{
		if (not jets_fRBX_isLoaded) {
			if (jets_fRBX_branch != 0) {
				jets_fRBX_branch->GetEntry(index);
			} else { 
				printf("branch jets_fRBX_branch does not exist!\n");
				exit(1);
			}
			jets_fRBX_isLoaded = true;
		}
		return jets_fRBX_;
	}
	vector<float> &jets_fSubDetector1()
	{
		if (not jets_fSubDetector1_isLoaded) {
			if (jets_fSubDetector1_branch != 0) {
				jets_fSubDetector1_branch->GetEntry(index);
			} else { 
				printf("branch jets_fSubDetector1_branch does not exist!\n");
				exit(1);
			}
			jets_fSubDetector1_isLoaded = true;
		}
		return jets_fSubDetector1_;
	}
	vector<float> &jets_fSubDetector2()
	{
		if (not jets_fSubDetector2_isLoaded) {
			if (jets_fSubDetector2_branch != 0) {
				jets_fSubDetector2_branch->GetEntry(index);
			} else { 
				printf("branch jets_fSubDetector2_branch does not exist!\n");
				exit(1);
			}
			jets_fSubDetector2_isLoaded = true;
		}
		return jets_fSubDetector2_;
	}
	vector<float> &jets_fSubDetector3()
	{
		if (not jets_fSubDetector3_isLoaded) {
			if (jets_fSubDetector3_branch != 0) {
				jets_fSubDetector3_branch->GetEntry(index);
			} else { 
				printf("branch jets_fSubDetector3_branch does not exist!\n");
				exit(1);
			}
			jets_fSubDetector3_isLoaded = true;
		}
		return jets_fSubDetector3_;
	}
	vector<float> &jets_fSubDetector4()
	{
		if (not jets_fSubDetector4_isLoaded) {
			if (jets_fSubDetector4_branch != 0) {
				jets_fSubDetector4_branch->GetEntry(index);
			} else { 
				printf("branch jets_fSubDetector4_branch does not exist!\n");
				exit(1);
			}
			jets_fSubDetector4_isLoaded = true;
		}
		return jets_fSubDetector4_;
	}
	vector<float> &jets_hitsInN90()
	{
		if (not jets_hitsInN90_isLoaded) {
			if (jets_hitsInN90_branch != 0) {
				jets_hitsInN90_branch->GetEntry(index);
			} else { 
				printf("branch jets_hitsInN90_branch does not exist!\n");
				exit(1);
			}
			jets_hitsInN90_isLoaded = true;
		}
		return jets_hitsInN90_;
	}
	vector<float> &jets_n90Hits()
	{
		if (not jets_n90Hits_isLoaded) {
			if (jets_n90Hits_branch != 0) {
				jets_n90Hits_branch->GetEntry(index);
			} else { 
				printf("branch jets_n90Hits_branch does not exist!\n");
				exit(1);
			}
			jets_n90Hits_isLoaded = true;
		}
		return jets_n90Hits_;
	}
	vector<float> &jets_nECALTowers()
	{
		if (not jets_nECALTowers_isLoaded) {
			if (jets_nECALTowers_branch != 0) {
				jets_nECALTowers_branch->GetEntry(index);
			} else { 
				printf("branch jets_nECALTowers_branch does not exist!\n");
				exit(1);
			}
			jets_nECALTowers_isLoaded = true;
		}
		return jets_nECALTowers_;
	}
	vector<float> &jets_nHCALTowers()
	{
		if (not jets_nHCALTowers_isLoaded) {
			if (jets_nHCALTowers_branch != 0) {
				jets_nHCALTowers_branch->GetEntry(index);
			} else { 
				printf("branch jets_nHCALTowers_branch does not exist!\n");
				exit(1);
			}
			jets_nHCALTowers_isLoaded = true;
		}
		return jets_nHCALTowers_;
	}
	vector<float> &jets_restrictedEMF()
	{
		if (not jets_restrictedEMF_isLoaded) {
			if (jets_restrictedEMF_branch != 0) {
				jets_restrictedEMF_branch->GetEntry(index);
			} else { 
				printf("branch jets_restrictedEMF_branch does not exist!\n");
				exit(1);
			}
			jets_restrictedEMF_isLoaded = true;
		}
		return jets_restrictedEMF_;
	}
	vector<float> &jpts_cor()
	{
		if (not jpts_cor_isLoaded) {
			if (jpts_cor_branch != 0) {
				jpts_cor_branch->GetEntry(index);
			} else { 
				printf("branch jpts_cor_branch does not exist!\n");
				exit(1);
			}
			jpts_cor_isLoaded = true;
		}
		return jpts_cor_;
	}
	vector<float> &jpts_emFrac()
	{
		if (not jpts_emFrac_isLoaded) {
			if (jpts_emFrac_branch != 0) {
				jpts_emFrac_branch->GetEntry(index);
			} else { 
				printf("branch jpts_emFrac_branch does not exist!\n");
				exit(1);
			}
			jpts_emFrac_isLoaded = true;
		}
		return jpts_emFrac_;
	}
	vector<float> &mus_met_deltax()
	{
		if (not mus_met_deltax_isLoaded) {
			if (mus_met_deltax_branch != 0) {
				mus_met_deltax_branch->GetEntry(index);
			} else { 
				printf("branch mus_met_deltax_branch does not exist!\n");
				exit(1);
			}
			mus_met_deltax_isLoaded = true;
		}
		return mus_met_deltax_;
	}
	vector<float> &mus_met_deltay()
	{
		if (not mus_met_deltay_isLoaded) {
			if (mus_met_deltay_branch != 0) {
				mus_met_deltay_branch->GetEntry(index);
			} else { 
				printf("branch mus_met_deltay_branch does not exist!\n");
				exit(1);
			}
			mus_met_deltay_isLoaded = true;
		}
		return mus_met_deltay_;
	}
	vector<float> &mus_eledr()
	{
		if (not mus_eledr_isLoaded) {
			if (mus_eledr_branch != 0) {
				mus_eledr_branch->GetEntry(index);
			} else { 
				printf("branch mus_eledr_branch does not exist!\n");
				exit(1);
			}
			mus_eledr_isLoaded = true;
		}
		return mus_eledr_;
	}
	vector<float> &mus_jetdr()
	{
		if (not mus_jetdr_isLoaded) {
			if (mus_jetdr_branch != 0) {
				mus_jetdr_branch->GetEntry(index);
			} else { 
				printf("branch mus_jetdr_branch does not exist!\n");
				exit(1);
			}
			mus_jetdr_isLoaded = true;
		}
		return mus_jetdr_;
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
	vector<float> &mus_gfit_d0()
	{
		if (not mus_gfit_d0_isLoaded) {
			if (mus_gfit_d0_branch != 0) {
				mus_gfit_d0_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_d0_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_d0_isLoaded = true;
		}
		return mus_gfit_d0_;
	}
	vector<float> &mus_gfit_d0Err()
	{
		if (not mus_gfit_d0Err_isLoaded) {
			if (mus_gfit_d0Err_branch != 0) {
				mus_gfit_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_d0Err_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_d0Err_isLoaded = true;
		}
		return mus_gfit_d0Err_;
	}
	vector<float> &mus_gfit_d0corr()
	{
		if (not mus_gfit_d0corr_isLoaded) {
			if (mus_gfit_d0corr_branch != 0) {
				mus_gfit_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_d0corr_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_d0corr_isLoaded = true;
		}
		return mus_gfit_d0corr_;
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
	vector<float> &mus_gfit_qoverp()
	{
		if (not mus_gfit_qoverp_isLoaded) {
			if (mus_gfit_qoverp_branch != 0) {
				mus_gfit_qoverp_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_qoverp_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_qoverp_isLoaded = true;
		}
		return mus_gfit_qoverp_;
	}
	vector<float> &mus_gfit_qoverpError()
	{
		if (not mus_gfit_qoverpError_isLoaded) {
			if (mus_gfit_qoverpError_branch != 0) {
				mus_gfit_qoverpError_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_qoverpError_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_qoverpError_isLoaded = true;
		}
		return mus_gfit_qoverpError_;
	}
	vector<float> &mus_gfit_z0()
	{
		if (not mus_gfit_z0_isLoaded) {
			if (mus_gfit_z0_branch != 0) {
				mus_gfit_z0_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_z0_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_z0_isLoaded = true;
		}
		return mus_gfit_z0_;
	}
	vector<float> &mus_gfit_z0Err()
	{
		if (not mus_gfit_z0Err_isLoaded) {
			if (mus_gfit_z0Err_branch != 0) {
				mus_gfit_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_z0Err_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_z0Err_isLoaded = true;
		}
		return mus_gfit_z0Err_;
	}
	vector<float> &mus_gfit_z0corr()
	{
		if (not mus_gfit_z0corr_isLoaded) {
			if (mus_gfit_z0corr_branch != 0) {
				mus_gfit_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_z0corr_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_z0corr_isLoaded = true;
		}
		return mus_gfit_z0corr_;
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
	vector<float> &mus_iso_ecalvetoDep()
	{
		if (not mus_iso_ecalvetoDep_isLoaded) {
			if (mus_iso_ecalvetoDep_branch != 0) {
				mus_iso_ecalvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso_ecalvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_iso_ecalvetoDep_isLoaded = true;
		}
		return mus_iso_ecalvetoDep_;
	}
	vector<float> &mus_iso_hcalvetoDep()
	{
		if (not mus_iso_hcalvetoDep_isLoaded) {
			if (mus_iso_hcalvetoDep_branch != 0) {
				mus_iso_hcalvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso_hcalvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_iso_hcalvetoDep_isLoaded = true;
		}
		return mus_iso_hcalvetoDep_;
	}
	vector<float> &mus_iso_hovetoDep()
	{
		if (not mus_iso_hovetoDep_isLoaded) {
			if (mus_iso_hovetoDep_branch != 0) {
				mus_iso_hovetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso_hovetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_iso_hovetoDep_isLoaded = true;
		}
		return mus_iso_hovetoDep_;
	}
	vector<float> &mus_iso_trckvetoDep()
	{
		if (not mus_iso_trckvetoDep_isLoaded) {
			if (mus_iso_trckvetoDep_branch != 0) {
				mus_iso_trckvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_iso_trckvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_iso_trckvetoDep_isLoaded = true;
		}
		return mus_iso_trckvetoDep_;
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
	vector<float> &mus_sta_chi2()
	{
		if (not mus_sta_chi2_isLoaded) {
			if (mus_sta_chi2_branch != 0) {
				mus_sta_chi2_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_chi2_branch does not exist!\n");
				exit(1);
			}
			mus_sta_chi2_isLoaded = true;
		}
		return mus_sta_chi2_;
	}
	vector<float> &mus_sta_d0()
	{
		if (not mus_sta_d0_isLoaded) {
			if (mus_sta_d0_branch != 0) {
				mus_sta_d0_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_d0_branch does not exist!\n");
				exit(1);
			}
			mus_sta_d0_isLoaded = true;
		}
		return mus_sta_d0_;
	}
	vector<float> &mus_sta_d0Err()
	{
		if (not mus_sta_d0Err_isLoaded) {
			if (mus_sta_d0Err_branch != 0) {
				mus_sta_d0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_d0Err_branch does not exist!\n");
				exit(1);
			}
			mus_sta_d0Err_isLoaded = true;
		}
		return mus_sta_d0Err_;
	}
	vector<float> &mus_sta_d0corr()
	{
		if (not mus_sta_d0corr_isLoaded) {
			if (mus_sta_d0corr_branch != 0) {
				mus_sta_d0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_d0corr_branch does not exist!\n");
				exit(1);
			}
			mus_sta_d0corr_isLoaded = true;
		}
		return mus_sta_d0corr_;
	}
	vector<float> &mus_sta_ndof()
	{
		if (not mus_sta_ndof_isLoaded) {
			if (mus_sta_ndof_branch != 0) {
				mus_sta_ndof_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_ndof_branch does not exist!\n");
				exit(1);
			}
			mus_sta_ndof_isLoaded = true;
		}
		return mus_sta_ndof_;
	}
	vector<float> &mus_sta_qoverp()
	{
		if (not mus_sta_qoverp_isLoaded) {
			if (mus_sta_qoverp_branch != 0) {
				mus_sta_qoverp_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_qoverp_branch does not exist!\n");
				exit(1);
			}
			mus_sta_qoverp_isLoaded = true;
		}
		return mus_sta_qoverp_;
	}
	vector<float> &mus_sta_qoverpError()
	{
		if (not mus_sta_qoverpError_isLoaded) {
			if (mus_sta_qoverpError_branch != 0) {
				mus_sta_qoverpError_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_qoverpError_branch does not exist!\n");
				exit(1);
			}
			mus_sta_qoverpError_isLoaded = true;
		}
		return mus_sta_qoverpError_;
	}
	vector<float> &mus_sta_z0()
	{
		if (not mus_sta_z0_isLoaded) {
			if (mus_sta_z0_branch != 0) {
				mus_sta_z0_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_z0_branch does not exist!\n");
				exit(1);
			}
			mus_sta_z0_isLoaded = true;
		}
		return mus_sta_z0_;
	}
	vector<float> &mus_sta_z0Err()
	{
		if (not mus_sta_z0Err_isLoaded) {
			if (mus_sta_z0Err_branch != 0) {
				mus_sta_z0Err_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_z0Err_branch does not exist!\n");
				exit(1);
			}
			mus_sta_z0Err_isLoaded = true;
		}
		return mus_sta_z0Err_;
	}
	vector<float> &mus_sta_z0corr()
	{
		if (not mus_sta_z0corr_isLoaded) {
			if (mus_sta_z0corr_branch != 0) {
				mus_sta_z0corr_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_z0corr_branch does not exist!\n");
				exit(1);
			}
			mus_sta_z0corr_isLoaded = true;
		}
		return mus_sta_z0corr_;
	}
	vector<float> &mus_timeAtIpInOut()
	{
		if (not mus_timeAtIpInOut_isLoaded) {
			if (mus_timeAtIpInOut_branch != 0) {
				mus_timeAtIpInOut_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeAtIpInOut_branch does not exist!\n");
				exit(1);
			}
			mus_timeAtIpInOut_isLoaded = true;
		}
		return mus_timeAtIpInOut_;
	}
	vector<float> &mus_timeAtIpInOutErr()
	{
		if (not mus_timeAtIpInOutErr_isLoaded) {
			if (mus_timeAtIpInOutErr_branch != 0) {
				mus_timeAtIpInOutErr_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeAtIpInOutErr_branch does not exist!\n");
				exit(1);
			}
			mus_timeAtIpInOutErr_isLoaded = true;
		}
		return mus_timeAtIpInOutErr_;
	}
	vector<float> &mus_timeAtIpOutIn()
	{
		if (not mus_timeAtIpOutIn_isLoaded) {
			if (mus_timeAtIpOutIn_branch != 0) {
				mus_timeAtIpOutIn_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeAtIpOutIn_branch does not exist!\n");
				exit(1);
			}
			mus_timeAtIpOutIn_isLoaded = true;
		}
		return mus_timeAtIpOutIn_;
	}
	vector<float> &mus_timeAtIpOutInErr()
	{
		if (not mus_timeAtIpOutInErr_isLoaded) {
			if (mus_timeAtIpOutInErr_branch != 0) {
				mus_timeAtIpOutInErr_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeAtIpOutInErr_branch does not exist!\n");
				exit(1);
			}
			mus_timeAtIpOutInErr_isLoaded = true;
		}
		return mus_timeAtIpOutInErr_;
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
	vector<float> &els_pat_caloIso()
	{
		if (not els_pat_caloIso_isLoaded) {
			if (els_pat_caloIso_branch != 0) {
				els_pat_caloIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_caloIso_branch does not exist!\n");
				exit(1);
			}
			els_pat_caloIso_isLoaded = true;
		}
		return els_pat_caloIso_;
	}
	vector<float> &els_pat_ecalIso()
	{
		if (not els_pat_ecalIso_isLoaded) {
			if (els_pat_ecalIso_branch != 0) {
				els_pat_ecalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_ecalIso_branch does not exist!\n");
				exit(1);
			}
			els_pat_ecalIso_isLoaded = true;
		}
		return els_pat_ecalIso_;
	}
	vector<float> &els_pat_hcalIso()
	{
		if (not els_pat_hcalIso_isLoaded) {
			if (els_pat_hcalIso_branch != 0) {
				els_pat_hcalIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_hcalIso_branch does not exist!\n");
				exit(1);
			}
			els_pat_hcalIso_isLoaded = true;
		}
		return els_pat_hcalIso_;
	}
	vector<float> &els_pat_looseId()
	{
		if (not els_pat_looseId_isLoaded) {
			if (els_pat_looseId_branch != 0) {
				els_pat_looseId_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_looseId_branch does not exist!\n");
				exit(1);
			}
			els_pat_looseId_isLoaded = true;
		}
		return els_pat_looseId_;
	}
	vector<float> &els_pat_robustHighEnergy()
	{
		if (not els_pat_robustHighEnergy_isLoaded) {
			if (els_pat_robustHighEnergy_branch != 0) {
				els_pat_robustHighEnergy_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_robustHighEnergy_branch does not exist!\n");
				exit(1);
			}
			els_pat_robustHighEnergy_isLoaded = true;
		}
		return els_pat_robustHighEnergy_;
	}
	vector<float> &els_pat_robustLooseId()
	{
		if (not els_pat_robustLooseId_isLoaded) {
			if (els_pat_robustLooseId_branch != 0) {
				els_pat_robustLooseId_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_robustLooseId_branch does not exist!\n");
				exit(1);
			}
			els_pat_robustLooseId_isLoaded = true;
		}
		return els_pat_robustLooseId_;
	}
	vector<float> &els_pat_robustTightId()
	{
		if (not els_pat_robustTightId_isLoaded) {
			if (els_pat_robustTightId_branch != 0) {
				els_pat_robustTightId_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_robustTightId_branch does not exist!\n");
				exit(1);
			}
			els_pat_robustTightId_isLoaded = true;
		}
		return els_pat_robustTightId_;
	}
	vector<float> &els_pat_scE1x5()
	{
		if (not els_pat_scE1x5_isLoaded) {
			if (els_pat_scE1x5_branch != 0) {
				els_pat_scE1x5_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_scE1x5_branch does not exist!\n");
				exit(1);
			}
			els_pat_scE1x5_isLoaded = true;
		}
		return els_pat_scE1x5_;
	}
	vector<float> &els_pat_scE2x5Max()
	{
		if (not els_pat_scE2x5Max_isLoaded) {
			if (els_pat_scE2x5Max_branch != 0) {
				els_pat_scE2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_scE2x5Max_branch does not exist!\n");
				exit(1);
			}
			els_pat_scE2x5Max_isLoaded = true;
		}
		return els_pat_scE2x5Max_;
	}
	vector<float> &els_pat_scE5x5()
	{
		if (not els_pat_scE5x5_isLoaded) {
			if (els_pat_scE5x5_branch != 0) {
				els_pat_scE5x5_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_scE5x5_branch does not exist!\n");
				exit(1);
			}
			els_pat_scE5x5_isLoaded = true;
		}
		return els_pat_scE5x5_;
	}
	vector<float> &els_pat_sigmaEtaEta()
	{
		if (not els_pat_sigmaEtaEta_isLoaded) {
			if (els_pat_sigmaEtaEta_branch != 0) {
				els_pat_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			els_pat_sigmaEtaEta_isLoaded = true;
		}
		return els_pat_sigmaEtaEta_;
	}
	vector<float> &els_pat_sigmaIEtaIEta()
	{
		if (not els_pat_sigmaIEtaIEta_isLoaded) {
			if (els_pat_sigmaIEtaIEta_branch != 0) {
				els_pat_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			els_pat_sigmaIEtaIEta_isLoaded = true;
		}
		return els_pat_sigmaIEtaIEta_;
	}
	vector<float> &els_pat_tightId()
	{
		if (not els_pat_tightId_isLoaded) {
			if (els_pat_tightId_branch != 0) {
				els_pat_tightId_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_tightId_branch does not exist!\n");
				exit(1);
			}
			els_pat_tightId_isLoaded = true;
		}
		return els_pat_tightId_;
	}
	vector<float> &els_pat_trackIso()
	{
		if (not els_pat_trackIso_isLoaded) {
			if (els_pat_trackIso_branch != 0) {
				els_pat_trackIso_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_trackIso_branch does not exist!\n");
				exit(1);
			}
			els_pat_trackIso_isLoaded = true;
		}
		return els_pat_trackIso_;
	}
	vector<float> &jets_pat_combinedSecondaryVertexBJetTag()
	{
		if (not jets_pat_combinedSecondaryVertexBJetTag_isLoaded) {
			if (jets_pat_combinedSecondaryVertexBJetTag_branch != 0) {
				jets_pat_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_combinedSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_combinedSecondaryVertexBJetTag_isLoaded = true;
		}
		return jets_pat_combinedSecondaryVertexBJetTag_;
	}
	vector<float> &jets_pat_combinedSecondaryVertexMVABJetTag()
	{
		if (not jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded) {
			if (jets_pat_combinedSecondaryVertexMVABJetTag_branch != 0) {
				jets_pat_combinedSecondaryVertexMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return jets_pat_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &jets_pat_coneIsolationTauJetTag()
	{
		if (not jets_pat_coneIsolationTauJetTag_isLoaded) {
			if (jets_pat_coneIsolationTauJetTag_branch != 0) {
				jets_pat_coneIsolationTauJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_coneIsolationTauJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_coneIsolationTauJetTag_isLoaded = true;
		}
		return jets_pat_coneIsolationTauJetTag_;
	}
	vector<float> &jets_pat_impactParameterMVABJetTag()
	{
		if (not jets_pat_impactParameterMVABJetTag_isLoaded) {
			if (jets_pat_impactParameterMVABJetTag_branch != 0) {
				jets_pat_impactParameterMVABJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_impactParameterMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_impactParameterMVABJetTag_isLoaded = true;
		}
		return jets_pat_impactParameterMVABJetTag_;
	}
	vector<float> &jets_pat_jetBProbabilityBJetTag()
	{
		if (not jets_pat_jetBProbabilityBJetTag_isLoaded) {
			if (jets_pat_jetBProbabilityBJetTag_branch != 0) {
				jets_pat_jetBProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetBProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetBProbabilityBJetTag_isLoaded = true;
		}
		return jets_pat_jetBProbabilityBJetTag_;
	}
	vector<float> &jets_pat_jetCharge()
	{
		if (not jets_pat_jetCharge_isLoaded) {
			if (jets_pat_jetCharge_branch != 0) {
				jets_pat_jetCharge_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetCharge_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetCharge_isLoaded = true;
		}
		return jets_pat_jetCharge_;
	}
	vector<float> &jets_pat_jetProbabilityBJetTag()
	{
		if (not jets_pat_jetProbabilityBJetTag_isLoaded) {
			if (jets_pat_jetProbabilityBJetTag_branch != 0) {
				jets_pat_jetProbabilityBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetProbabilityBJetTag_isLoaded = true;
		}
		return jets_pat_jetProbabilityBJetTag_;
	}
	vector<float> &jets_pat_noCorrF()
	{
		if (not jets_pat_noCorrF_isLoaded) {
			if (jets_pat_noCorrF_branch != 0) {
				jets_pat_noCorrF_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_noCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_noCorrF_isLoaded = true;
		}
		return jets_pat_noCorrF_;
	}
	vector<float> &jets_pat_simpleSecondaryVertexHighEffBJetTag()
	{
		if (not jets_pat_simpleSecondaryVertexHighEffBJetTag_isLoaded) {
			if (jets_pat_simpleSecondaryVertexHighEffBJetTag_branch != 0) {
				jets_pat_simpleSecondaryVertexHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_simpleSecondaryVertexHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_simpleSecondaryVertexHighEffBJetTag_isLoaded = true;
		}
		return jets_pat_simpleSecondaryVertexHighEffBJetTag_;
	}
	vector<float> &jets_pat_simpleSecondaryVertexHighPurBJetTag()
	{
		if (not jets_pat_simpleSecondaryVertexHighPurBJetTag_isLoaded) {
			if (jets_pat_simpleSecondaryVertexHighPurBJetTag_branch != 0) {
				jets_pat_simpleSecondaryVertexHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_simpleSecondaryVertexHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_simpleSecondaryVertexHighPurBJetTag_isLoaded = true;
		}
		return jets_pat_simpleSecondaryVertexHighPurBJetTag_;
	}
	vector<float> &jets_pat_softElectronByIP3dBJetTag()
	{
		if (not jets_pat_softElectronByIP3dBJetTag_isLoaded) {
			if (jets_pat_softElectronByIP3dBJetTag_branch != 0) {
				jets_pat_softElectronByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_softElectronByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softElectronByIP3dBJetTag_isLoaded = true;
		}
		return jets_pat_softElectronByIP3dBJetTag_;
	}
	vector<float> &jets_pat_softElectronByPtBJetTag()
	{
		if (not jets_pat_softElectronByPtBJetTag_isLoaded) {
			if (jets_pat_softElectronByPtBJetTag_branch != 0) {
				jets_pat_softElectronByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_softElectronByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softElectronByPtBJetTag_isLoaded = true;
		}
		return jets_pat_softElectronByPtBJetTag_;
	}
	vector<float> &jets_pat_softMuonBJetTag()
	{
		if (not jets_pat_softMuonBJetTag_isLoaded) {
			if (jets_pat_softMuonBJetTag_branch != 0) {
				jets_pat_softMuonBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softMuonBJetTag_isLoaded = true;
		}
		return jets_pat_softMuonBJetTag_;
	}
	vector<float> &jets_pat_softMuonByIP3dBJetTag()
	{
		if (not jets_pat_softMuonByIP3dBJetTag_isLoaded) {
			if (jets_pat_softMuonByIP3dBJetTag_branch != 0) {
				jets_pat_softMuonByIP3dBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_softMuonByIP3dBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softMuonByIP3dBJetTag_isLoaded = true;
		}
		return jets_pat_softMuonByIP3dBJetTag_;
	}
	vector<float> &jets_pat_softMuonByPtBJetTag()
	{
		if (not jets_pat_softMuonByPtBJetTag_isLoaded) {
			if (jets_pat_softMuonByPtBJetTag_branch != 0) {
				jets_pat_softMuonByPtBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_softMuonByPtBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softMuonByPtBJetTag_isLoaded = true;
		}
		return jets_pat_softMuonByPtBJetTag_;
	}
	vector<float> &jets_pat_trackCountingHighEffBJetTag()
	{
		if (not jets_pat_trackCountingHighEffBJetTag_isLoaded) {
			if (jets_pat_trackCountingHighEffBJetTag_branch != 0) {
				jets_pat_trackCountingHighEffBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_trackCountingHighEffBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_trackCountingHighEffBJetTag_isLoaded = true;
		}
		return jets_pat_trackCountingHighEffBJetTag_;
	}
	vector<float> &jets_pat_trackCountingHighPurBJetTag()
	{
		if (not jets_pat_trackCountingHighPurBJetTag_isLoaded) {
			if (jets_pat_trackCountingHighPurBJetTag_branch != 0) {
				jets_pat_trackCountingHighPurBJetTag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return jets_pat_trackCountingHighPurBJetTag_;
	}
	vector<float> &mus_pat_caloIso()
	{
		if (not mus_pat_caloIso_isLoaded) {
			if (mus_pat_caloIso_branch != 0) {
				mus_pat_caloIso_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_caloIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_caloIso_isLoaded = true;
		}
		return mus_pat_caloIso_;
	}
	vector<float> &mus_pat_calovetoDep()
	{
		if (not mus_pat_calovetoDep_isLoaded) {
			if (mus_pat_calovetoDep_branch != 0) {
				mus_pat_calovetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_calovetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_pat_calovetoDep_isLoaded = true;
		}
		return mus_pat_calovetoDep_;
	}
	vector<float> &mus_pat_ecalIso()
	{
		if (not mus_pat_ecalIso_isLoaded) {
			if (mus_pat_ecalIso_branch != 0) {
				mus_pat_ecalIso_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_ecalIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_ecalIso_isLoaded = true;
		}
		return mus_pat_ecalIso_;
	}
	vector<float> &mus_pat_ecalvetoDep()
	{
		if (not mus_pat_ecalvetoDep_isLoaded) {
			if (mus_pat_ecalvetoDep_branch != 0) {
				mus_pat_ecalvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_ecalvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_pat_ecalvetoDep_isLoaded = true;
		}
		return mus_pat_ecalvetoDep_;
	}
	vector<float> &mus_pat_hcalIso()
	{
		if (not mus_pat_hcalIso_isLoaded) {
			if (mus_pat_hcalIso_branch != 0) {
				mus_pat_hcalIso_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_hcalIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_hcalIso_isLoaded = true;
		}
		return mus_pat_hcalIso_;
	}
	vector<float> &mus_pat_hcalvetoDep()
	{
		if (not mus_pat_hcalvetoDep_isLoaded) {
			if (mus_pat_hcalvetoDep_branch != 0) {
				mus_pat_hcalvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_hcalvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_pat_hcalvetoDep_isLoaded = true;
		}
		return mus_pat_hcalvetoDep_;
	}
	vector<float> &mus_pat_trackIso()
	{
		if (not mus_pat_trackIso_isLoaded) {
			if (mus_pat_trackIso_branch != 0) {
				mus_pat_trackIso_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_trackIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_trackIso_isLoaded = true;
		}
		return mus_pat_trackIso_;
	}
	vector<float> &mus_pat_trckvetoDep()
	{
		if (not mus_pat_trckvetoDep_isLoaded) {
			if (mus_pat_trckvetoDep_branch != 0) {
				mus_pat_trckvetoDep_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_trckvetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_pat_trckvetoDep_isLoaded = true;
		}
		return mus_pat_trckvetoDep_;
	}
	vector<float> &pfels_deltaP()
	{
		if (not pfels_deltaP_isLoaded) {
			if (pfels_deltaP_branch != 0) {
				pfels_deltaP_branch->GetEntry(index);
			} else { 
				printf("branch pfels_deltaP_branch does not exist!\n");
				exit(1);
			}
			pfels_deltaP_isLoaded = true;
		}
		return pfels_deltaP_;
	}
	vector<float> &pfels_ecalE()
	{
		if (not pfels_ecalE_isLoaded) {
			if (pfels_ecalE_branch != 0) {
				pfels_ecalE_branch->GetEntry(index);
			} else { 
				printf("branch pfels_ecalE_branch does not exist!\n");
				exit(1);
			}
			pfels_ecalE_isLoaded = true;
		}
		return pfels_ecalE_;
	}
	vector<float> &pfels_hcalE()
	{
		if (not pfels_hcalE_isLoaded) {
			if (pfels_hcalE_branch != 0) {
				pfels_hcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfels_hcalE_branch does not exist!\n");
				exit(1);
			}
			pfels_hcalE_isLoaded = true;
		}
		return pfels_hcalE_;
	}
	vector<float> &pfels_isoChargedHadrons()
	{
		if (not pfels_isoChargedHadrons_isLoaded) {
			if (pfels_isoChargedHadrons_branch != 0) {
				pfels_isoChargedHadrons_branch->GetEntry(index);
			} else { 
				printf("branch pfels_isoChargedHadrons_branch does not exist!\n");
				exit(1);
			}
			pfels_isoChargedHadrons_isLoaded = true;
		}
		return pfels_isoChargedHadrons_;
	}
	vector<float> &pfels_isoNeutralHadrons()
	{
		if (not pfels_isoNeutralHadrons_isLoaded) {
			if (pfels_isoNeutralHadrons_branch != 0) {
				pfels_isoNeutralHadrons_branch->GetEntry(index);
			} else { 
				printf("branch pfels_isoNeutralHadrons_branch does not exist!\n");
				exit(1);
			}
			pfels_isoNeutralHadrons_isLoaded = true;
		}
		return pfels_isoNeutralHadrons_;
	}
	vector<float> &pfels_isoPhotons()
	{
		if (not pfels_isoPhotons_isLoaded) {
			if (pfels_isoPhotons_branch != 0) {
				pfels_isoPhotons_branch->GetEntry(index);
			} else { 
				printf("branch pfels_isoPhotons_branch does not exist!\n");
				exit(1);
			}
			pfels_isoPhotons_isLoaded = true;
		}
		return pfels_isoPhotons_;
	}
	vector<float> &pfels_mva_emu()
	{
		if (not pfels_mva_emu_isLoaded) {
			if (pfels_mva_emu_branch != 0) {
				pfels_mva_emu_branch->GetEntry(index);
			} else { 
				printf("branch pfels_mva_emu_branch does not exist!\n");
				exit(1);
			}
			pfels_mva_emu_isLoaded = true;
		}
		return pfels_mva_emu_;
	}
	vector<float> &pfels_mva_epi()
	{
		if (not pfels_mva_epi_isLoaded) {
			if (pfels_mva_epi_branch != 0) {
				pfels_mva_epi_branch->GetEntry(index);
			} else { 
				printf("branch pfels_mva_epi_branch does not exist!\n");
				exit(1);
			}
			pfels_mva_epi_isLoaded = true;
		}
		return pfels_mva_epi_;
	}
	vector<float> &pfels_mva_nothing_gamma()
	{
		if (not pfels_mva_nothing_gamma_isLoaded) {
			if (pfels_mva_nothing_gamma_branch != 0) {
				pfels_mva_nothing_gamma_branch->GetEntry(index);
			} else { 
				printf("branch pfels_mva_nothing_gamma_branch does not exist!\n");
				exit(1);
			}
			pfels_mva_nothing_gamma_isLoaded = true;
		}
		return pfels_mva_nothing_gamma_;
	}
	vector<float> &pfels_mva_nothing_nh()
	{
		if (not pfels_mva_nothing_nh_isLoaded) {
			if (pfels_mva_nothing_nh_branch != 0) {
				pfels_mva_nothing_nh_branch->GetEntry(index);
			} else { 
				printf("branch pfels_mva_nothing_nh_branch does not exist!\n");
				exit(1);
			}
			pfels_mva_nothing_nh_isLoaded = true;
		}
		return pfels_mva_nothing_nh_;
	}
	vector<float> &pfels_mva_pimu()
	{
		if (not pfels_mva_pimu_isLoaded) {
			if (pfels_mva_pimu_branch != 0) {
				pfels_mva_pimu_branch->GetEntry(index);
			} else { 
				printf("branch pfels_mva_pimu_branch does not exist!\n");
				exit(1);
			}
			pfels_mva_pimu_isLoaded = true;
		}
		return pfels_mva_pimu_;
	}
	vector<float> &pfels_pS1E()
	{
		if (not pfels_pS1E_isLoaded) {
			if (pfels_pS1E_branch != 0) {
				pfels_pS1E_branch->GetEntry(index);
			} else { 
				printf("branch pfels_pS1E_branch does not exist!\n");
				exit(1);
			}
			pfels_pS1E_isLoaded = true;
		}
		return pfels_pS1E_;
	}
	vector<float> &pfels_pS2E()
	{
		if (not pfels_pS2E_isLoaded) {
			if (pfels_pS2E_branch != 0) {
				pfels_pS2E_branch->GetEntry(index);
			} else { 
				printf("branch pfels_pS2E_branch does not exist!\n");
				exit(1);
			}
			pfels_pS2E_isLoaded = true;
		}
		return pfels_pS2E_;
	}
	vector<float> &pfels_rawEcalE()
	{
		if (not pfels_rawEcalE_isLoaded) {
			if (pfels_rawEcalE_branch != 0) {
				pfels_rawEcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfels_rawEcalE_branch does not exist!\n");
				exit(1);
			}
			pfels_rawEcalE_isLoaded = true;
		}
		return pfels_rawEcalE_;
	}
	vector<float> &pfels_rawHcalE()
	{
		if (not pfels_rawHcalE_isLoaded) {
			if (pfels_rawHcalE_branch != 0) {
				pfels_rawHcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfels_rawHcalE_branch does not exist!\n");
				exit(1);
			}
			pfels_rawHcalE_isLoaded = true;
		}
		return pfels_rawHcalE_;
	}
	vector<float> &pfjets_chargedEmE()
	{
		if (not pfjets_chargedEmE_isLoaded) {
			if (pfjets_chargedEmE_branch != 0) {
				pfjets_chargedEmE_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_chargedEmE_branch does not exist!\n");
				exit(1);
			}
			pfjets_chargedEmE_isLoaded = true;
		}
		return pfjets_chargedEmE_;
	}
	vector<float> &pfjets_chargedHadronE()
	{
		if (not pfjets_chargedHadronE_isLoaded) {
			if (pfjets_chargedHadronE_branch != 0) {
				pfjets_chargedHadronE_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_chargedHadronE_branch does not exist!\n");
				exit(1);
			}
			pfjets_chargedHadronE_isLoaded = true;
		}
		return pfjets_chargedHadronE_;
	}
	vector<float> &pfjets_cor()
	{
		if (not pfjets_cor_isLoaded) {
			if (pfjets_cor_branch != 0) {
				pfjets_cor_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_cor_branch does not exist!\n");
				exit(1);
			}
			pfjets_cor_isLoaded = true;
		}
		return pfjets_cor_;
	}
	vector<float> &pfjets_neutralEmE()
	{
		if (not pfjets_neutralEmE_isLoaded) {
			if (pfjets_neutralEmE_branch != 0) {
				pfjets_neutralEmE_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_neutralEmE_branch does not exist!\n");
				exit(1);
			}
			pfjets_neutralEmE_isLoaded = true;
		}
		return pfjets_neutralEmE_;
	}
	vector<float> &pfjets_neutralHadronE()
	{
		if (not pfjets_neutralHadronE_isLoaded) {
			if (pfjets_neutralHadronE_branch != 0) {
				pfjets_neutralHadronE_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_neutralHadronE_branch does not exist!\n");
				exit(1);
			}
			pfjets_neutralHadronE_isLoaded = true;
		}
		return pfjets_neutralHadronE_;
	}
	vector<float> &pfmus_deltaP()
	{
		if (not pfmus_deltaP_isLoaded) {
			if (pfmus_deltaP_branch != 0) {
				pfmus_deltaP_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_deltaP_branch does not exist!\n");
				exit(1);
			}
			pfmus_deltaP_isLoaded = true;
		}
		return pfmus_deltaP_;
	}
	vector<float> &pfmus_ecalE()
	{
		if (not pfmus_ecalE_isLoaded) {
			if (pfmus_ecalE_branch != 0) {
				pfmus_ecalE_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_ecalE_branch does not exist!\n");
				exit(1);
			}
			pfmus_ecalE_isLoaded = true;
		}
		return pfmus_ecalE_;
	}
	vector<float> &pfmus_hcalE()
	{
		if (not pfmus_hcalE_isLoaded) {
			if (pfmus_hcalE_branch != 0) {
				pfmus_hcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_hcalE_branch does not exist!\n");
				exit(1);
			}
			pfmus_hcalE_isLoaded = true;
		}
		return pfmus_hcalE_;
	}
	vector<float> &pfmus_isoChargedHadrons()
	{
		if (not pfmus_isoChargedHadrons_isLoaded) {
			if (pfmus_isoChargedHadrons_branch != 0) {
				pfmus_isoChargedHadrons_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_isoChargedHadrons_branch does not exist!\n");
				exit(1);
			}
			pfmus_isoChargedHadrons_isLoaded = true;
		}
		return pfmus_isoChargedHadrons_;
	}
	vector<float> &pfmus_isoNeutralHadrons()
	{
		if (not pfmus_isoNeutralHadrons_isLoaded) {
			if (pfmus_isoNeutralHadrons_branch != 0) {
				pfmus_isoNeutralHadrons_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_isoNeutralHadrons_branch does not exist!\n");
				exit(1);
			}
			pfmus_isoNeutralHadrons_isLoaded = true;
		}
		return pfmus_isoNeutralHadrons_;
	}
	vector<float> &pfmus_isoPhotons()
	{
		if (not pfmus_isoPhotons_isLoaded) {
			if (pfmus_isoPhotons_branch != 0) {
				pfmus_isoPhotons_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_isoPhotons_branch does not exist!\n");
				exit(1);
			}
			pfmus_isoPhotons_isLoaded = true;
		}
		return pfmus_isoPhotons_;
	}
	vector<float> &pfmus_mva_emu()
	{
		if (not pfmus_mva_emu_isLoaded) {
			if (pfmus_mva_emu_branch != 0) {
				pfmus_mva_emu_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_mva_emu_branch does not exist!\n");
				exit(1);
			}
			pfmus_mva_emu_isLoaded = true;
		}
		return pfmus_mva_emu_;
	}
	vector<float> &pfmus_mva_epi()
	{
		if (not pfmus_mva_epi_isLoaded) {
			if (pfmus_mva_epi_branch != 0) {
				pfmus_mva_epi_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_mva_epi_branch does not exist!\n");
				exit(1);
			}
			pfmus_mva_epi_isLoaded = true;
		}
		return pfmus_mva_epi_;
	}
	vector<float> &pfmus_mva_nothing_gamma()
	{
		if (not pfmus_mva_nothing_gamma_isLoaded) {
			if (pfmus_mva_nothing_gamma_branch != 0) {
				pfmus_mva_nothing_gamma_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_mva_nothing_gamma_branch does not exist!\n");
				exit(1);
			}
			pfmus_mva_nothing_gamma_isLoaded = true;
		}
		return pfmus_mva_nothing_gamma_;
	}
	vector<float> &pfmus_mva_nothing_nh()
	{
		if (not pfmus_mva_nothing_nh_isLoaded) {
			if (pfmus_mva_nothing_nh_branch != 0) {
				pfmus_mva_nothing_nh_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_mva_nothing_nh_branch does not exist!\n");
				exit(1);
			}
			pfmus_mva_nothing_nh_isLoaded = true;
		}
		return pfmus_mva_nothing_nh_;
	}
	vector<float> &pfmus_mva_pimu()
	{
		if (not pfmus_mva_pimu_isLoaded) {
			if (pfmus_mva_pimu_branch != 0) {
				pfmus_mva_pimu_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_mva_pimu_branch does not exist!\n");
				exit(1);
			}
			pfmus_mva_pimu_isLoaded = true;
		}
		return pfmus_mva_pimu_;
	}
	vector<float> &pfmus_pS1E()
	{
		if (not pfmus_pS1E_isLoaded) {
			if (pfmus_pS1E_branch != 0) {
				pfmus_pS1E_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_pS1E_branch does not exist!\n");
				exit(1);
			}
			pfmus_pS1E_isLoaded = true;
		}
		return pfmus_pS1E_;
	}
	vector<float> &pfmus_pS2E()
	{
		if (not pfmus_pS2E_isLoaded) {
			if (pfmus_pS2E_branch != 0) {
				pfmus_pS2E_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_pS2E_branch does not exist!\n");
				exit(1);
			}
			pfmus_pS2E_isLoaded = true;
		}
		return pfmus_pS2E_;
	}
	vector<float> &pfmus_rawEcalE()
	{
		if (not pfmus_rawEcalE_isLoaded) {
			if (pfmus_rawEcalE_branch != 0) {
				pfmus_rawEcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_rawEcalE_branch does not exist!\n");
				exit(1);
			}
			pfmus_rawEcalE_isLoaded = true;
		}
		return pfmus_rawEcalE_;
	}
	vector<float> &pfmus_rawHcalE()
	{
		if (not pfmus_rawHcalE_isLoaded) {
			if (pfmus_rawHcalE_branch != 0) {
				pfmus_rawHcalE_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_rawHcalE_branch does not exist!\n");
				exit(1);
			}
			pfmus_rawHcalE_isLoaded = true;
		}
		return pfmus_rawHcalE_;
	}
	vector<float> &photons_e1x5()
	{
		if (not photons_e1x5_isLoaded) {
			if (photons_e1x5_branch != 0) {
				photons_e1x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_e1x5_branch does not exist!\n");
				exit(1);
			}
			photons_e1x5_isLoaded = true;
		}
		return photons_e1x5_;
	}
	vector<float> &photons_e2x5Max()
	{
		if (not photons_e2x5Max_isLoaded) {
			if (photons_e2x5Max_branch != 0) {
				photons_e2x5Max_branch->GetEntry(index);
			} else { 
				printf("branch photons_e2x5Max_branch does not exist!\n");
				exit(1);
			}
			photons_e2x5Max_isLoaded = true;
		}
		return photons_e2x5Max_;
	}
	vector<float> &photons_e3x3()
	{
		if (not photons_e3x3_isLoaded) {
			if (photons_e3x3_branch != 0) {
				photons_e3x3_branch->GetEntry(index);
			} else { 
				printf("branch photons_e3x3_branch does not exist!\n");
				exit(1);
			}
			photons_e3x3_isLoaded = true;
		}
		return photons_e3x3_;
	}
	vector<float> &photons_e5x5()
	{
		if (not photons_e5x5_isLoaded) {
			if (photons_e5x5_branch != 0) {
				photons_e5x5_branch->GetEntry(index);
			} else { 
				printf("branch photons_e5x5_branch does not exist!\n");
				exit(1);
			}
			photons_e5x5_isLoaded = true;
		}
		return photons_e5x5_;
	}
	vector<float> &photons_ecalIso03()
	{
		if (not photons_ecalIso03_isLoaded) {
			if (photons_ecalIso03_branch != 0) {
				photons_ecalIso03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ecalIso03_branch does not exist!\n");
				exit(1);
			}
			photons_ecalIso03_isLoaded = true;
		}
		return photons_ecalIso03_;
	}
	vector<float> &photons_ecalIso04()
	{
		if (not photons_ecalIso04_isLoaded) {
			if (photons_ecalIso04_branch != 0) {
				photons_ecalIso04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ecalIso04_branch does not exist!\n");
				exit(1);
			}
			photons_ecalIso04_isLoaded = true;
		}
		return photons_ecalIso04_;
	}
	vector<float> &photons_hOverE()
	{
		if (not photons_hOverE_isLoaded) {
			if (photons_hOverE_branch != 0) {
				photons_hOverE_branch->GetEntry(index);
			} else { 
				printf("branch photons_hOverE_branch does not exist!\n");
				exit(1);
			}
			photons_hOverE_isLoaded = true;
		}
		return photons_hOverE_;
	}
	vector<float> &photons_hcalIso03()
	{
		if (not photons_hcalIso03_isLoaded) {
			if (photons_hcalIso03_branch != 0) {
				photons_hcalIso03_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalIso03_branch does not exist!\n");
				exit(1);
			}
			photons_hcalIso03_isLoaded = true;
		}
		return photons_hcalIso03_;
	}
	vector<float> &photons_hcalIso04()
	{
		if (not photons_hcalIso04_isLoaded) {
			if (photons_hcalIso04_branch != 0) {
				photons_hcalIso04_branch->GetEntry(index);
			} else { 
				printf("branch photons_hcalIso04_branch does not exist!\n");
				exit(1);
			}
			photons_hcalIso04_isLoaded = true;
		}
		return photons_hcalIso04_;
	}
	vector<float> &photons_ntkIsoHollow03()
	{
		if (not photons_ntkIsoHollow03_isLoaded) {
			if (photons_ntkIsoHollow03_branch != 0) {
				photons_ntkIsoHollow03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoHollow03_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoHollow03_isLoaded = true;
		}
		return photons_ntkIsoHollow03_;
	}
	vector<float> &photons_ntkIsoHollow04()
	{
		if (not photons_ntkIsoHollow04_isLoaded) {
			if (photons_ntkIsoHollow04_branch != 0) {
				photons_ntkIsoHollow04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoHollow04_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoHollow04_isLoaded = true;
		}
		return photons_ntkIsoHollow04_;
	}
	vector<float> &photons_ntkIsoSolid03()
	{
		if (not photons_ntkIsoSolid03_isLoaded) {
			if (photons_ntkIsoSolid03_branch != 0) {
				photons_ntkIsoSolid03_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoSolid03_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoSolid03_isLoaded = true;
		}
		return photons_ntkIsoSolid03_;
	}
	vector<float> &photons_ntkIsoSolid04()
	{
		if (not photons_ntkIsoSolid04_isLoaded) {
			if (photons_ntkIsoSolid04_branch != 0) {
				photons_ntkIsoSolid04_branch->GetEntry(index);
			} else { 
				printf("branch photons_ntkIsoSolid04_branch does not exist!\n");
				exit(1);
			}
			photons_ntkIsoSolid04_isLoaded = true;
		}
		return photons_ntkIsoSolid04_;
	}
	vector<float> &photons_sigmaEtaEta()
	{
		if (not photons_sigmaEtaEta_isLoaded) {
			if (photons_sigmaEtaEta_branch != 0) {
				photons_sigmaEtaEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_sigmaEtaEta_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaEtaEta_isLoaded = true;
		}
		return photons_sigmaEtaEta_;
	}
	vector<float> &photons_sigmaIEtaIEta()
	{
		if (not photons_sigmaIEtaIEta_isLoaded) {
			if (photons_sigmaIEtaIEta_branch != 0) {
				photons_sigmaIEtaIEta_branch->GetEntry(index);
			} else { 
				printf("branch photons_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_sigmaIEtaIEta_;
	}
	vector<float> &photons_swissSeed()
	{
		if (not photons_swissSeed_isLoaded) {
			if (photons_swissSeed_branch != 0) {
				photons_swissSeed_branch->GetEntry(index);
			} else { 
				printf("branch photons_swissSeed_branch does not exist!\n");
				exit(1);
			}
			photons_swissSeed_isLoaded = true;
		}
		return photons_swissSeed_;
	}
	vector<float> &photons_tkIsoHollow03()
	{
		if (not photons_tkIsoHollow03_isLoaded) {
			if (photons_tkIsoHollow03_branch != 0) {
				photons_tkIsoHollow03_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoHollow03_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoHollow03_isLoaded = true;
		}
		return photons_tkIsoHollow03_;
	}
	vector<float> &photons_tkIsoHollow04()
	{
		if (not photons_tkIsoHollow04_isLoaded) {
			if (photons_tkIsoHollow04_branch != 0) {
				photons_tkIsoHollow04_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoHollow04_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoHollow04_isLoaded = true;
		}
		return photons_tkIsoHollow04_;
	}
	vector<float> &photons_tkIsoSolid03()
	{
		if (not photons_tkIsoSolid03_isLoaded) {
			if (photons_tkIsoSolid03_branch != 0) {
				photons_tkIsoSolid03_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoSolid03_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoSolid03_isLoaded = true;
		}
		return photons_tkIsoSolid03_;
	}
	vector<float> &photons_tkIsoSolid04()
	{
		if (not photons_tkIsoSolid04_isLoaded) {
			if (photons_tkIsoSolid04_branch != 0) {
				photons_tkIsoSolid04_branch->GetEntry(index);
			} else { 
				printf("branch photons_tkIsoSolid04_branch does not exist!\n");
				exit(1);
			}
			photons_tkIsoSolid04_isLoaded = true;
		}
		return photons_tkIsoSolid04_;
	}
	vector<float> &scs_clustersSize()
	{
		if (not scs_clustersSize_isLoaded) {
			if (scs_clustersSize_branch != 0) {
				scs_clustersSize_branch->GetEntry(index);
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
			} else { 
				printf("branch scs_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIEta_isLoaded = true;
		}
		return scs_sigmaIEtaIEta_;
	}
	vector<float> &scs_sigmaIEtaIEtaSC()
	{
		if (not scs_sigmaIEtaIEtaSC_isLoaded) {
			if (scs_sigmaIEtaIEtaSC_branch != 0) {
				scs_sigmaIEtaIEtaSC_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaIEtaIEtaSC_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIEtaSC_isLoaded = true;
		}
		return scs_sigmaIEtaIEtaSC_;
	}
	vector<float> &scs_sigmaIEtaIPhi()
	{
		if (not scs_sigmaIEtaIPhi_isLoaded) {
			if (scs_sigmaIEtaIPhi_branch != 0) {
				scs_sigmaIEtaIPhi_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaIEtaIPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIPhi_isLoaded = true;
		}
		return scs_sigmaIEtaIPhi_;
	}
	vector<float> &scs_sigmaIEtaIPhiSC()
	{
		if (not scs_sigmaIEtaIPhiSC_isLoaded) {
			if (scs_sigmaIEtaIPhiSC_branch != 0) {
				scs_sigmaIEtaIPhiSC_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaIEtaIPhiSC_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIEtaIPhiSC_isLoaded = true;
		}
		return scs_sigmaIEtaIPhiSC_;
	}
	vector<float> &scs_sigmaIPhiIPhi()
	{
		if (not scs_sigmaIPhiIPhi_isLoaded) {
			if (scs_sigmaIPhiIPhi_branch != 0) {
				scs_sigmaIPhiIPhi_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaIPhiIPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIPhiIPhi_isLoaded = true;
		}
		return scs_sigmaIPhiIPhi_;
	}
	vector<float> &scs_sigmaIPhiIPhiSC()
	{
		if (not scs_sigmaIPhiIPhiSC_isLoaded) {
			if (scs_sigmaIPhiIPhiSC_branch != 0) {
				scs_sigmaIPhiIPhiSC_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaIPhiIPhiSC_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaIPhiIPhiSC_isLoaded = true;
		}
		return scs_sigmaIPhiIPhiSC_;
	}
	vector<float> &scs_sigmaPhiPhi()
	{
		if (not scs_sigmaPhiPhi_isLoaded) {
			if (scs_sigmaPhiPhi_branch != 0) {
				scs_sigmaPhiPhi_branch->GetEntry(index);
			} else { 
				printf("branch scs_sigmaPhiPhi_branch does not exist!\n");
				exit(1);
			}
			scs_sigmaPhiPhi_isLoaded = true;
		}
		return scs_sigmaPhiPhi_;
	}
	vector<float> &scs_timeSeed()
	{
		if (not scs_timeSeed_isLoaded) {
			if (scs_timeSeed_branch != 0) {
				scs_timeSeed_branch->GetEntry(index);
			} else { 
				printf("branch scs_timeSeed_branch does not exist!\n");
				exit(1);
			}
			scs_timeSeed_isLoaded = true;
		}
		return scs_timeSeed_;
	}
	vector<float> &svs_anglePV()
	{
		if (not svs_anglePV_isLoaded) {
			if (svs_anglePV_branch != 0) {
				svs_anglePV_branch->GetEntry(index);
			} else { 
				printf("branch svs_anglePV_branch does not exist!\n");
				exit(1);
			}
			svs_anglePV_isLoaded = true;
		}
		return svs_anglePV_;
	}
	vector<float> &svs_chi2()
	{
		if (not svs_chi2_isLoaded) {
			if (svs_chi2_branch != 0) {
				svs_chi2_branch->GetEntry(index);
			} else { 
				printf("branch svs_chi2_branch does not exist!\n");
				exit(1);
			}
			svs_chi2_isLoaded = true;
		}
		return svs_chi2_;
	}
	vector<float> &svs_dist3Dsig()
	{
		if (not svs_dist3Dsig_isLoaded) {
			if (svs_dist3Dsig_branch != 0) {
				svs_dist3Dsig_branch->GetEntry(index);
			} else { 
				printf("branch svs_dist3Dsig_branch does not exist!\n");
				exit(1);
			}
			svs_dist3Dsig_isLoaded = true;
		}
		return svs_dist3Dsig_;
	}
	vector<float> &svs_dist3Dval()
	{
		if (not svs_dist3Dval_isLoaded) {
			if (svs_dist3Dval_branch != 0) {
				svs_dist3Dval_branch->GetEntry(index);
			} else { 
				printf("branch svs_dist3Dval_branch does not exist!\n");
				exit(1);
			}
			svs_dist3Dval_isLoaded = true;
		}
		return svs_dist3Dval_;
	}
	vector<float> &svs_distXYsig()
	{
		if (not svs_distXYsig_isLoaded) {
			if (svs_distXYsig_branch != 0) {
				svs_distXYsig_branch->GetEntry(index);
			} else { 
				printf("branch svs_distXYsig_branch does not exist!\n");
				exit(1);
			}
			svs_distXYsig_isLoaded = true;
		}
		return svs_distXYsig_;
	}
	vector<float> &svs_distXYval()
	{
		if (not svs_distXYval_isLoaded) {
			if (svs_distXYval_branch != 0) {
				svs_distXYval_branch->GetEntry(index);
			} else { 
				printf("branch svs_distXYval_branch does not exist!\n");
				exit(1);
			}
			svs_distXYval_isLoaded = true;
		}
		return svs_distXYval_;
	}
	vector<float> &svs_ndof()
	{
		if (not svs_ndof_isLoaded) {
			if (svs_ndof_branch != 0) {
				svs_ndof_branch->GetEntry(index);
			} else { 
				printf("branch svs_ndof_branch does not exist!\n");
				exit(1);
			}
			svs_ndof_isLoaded = true;
		}
		return svs_ndof_;
	}
	vector<float> &svs_prob()
	{
		if (not svs_prob_isLoaded) {
			if (svs_prob_branch != 0) {
				svs_prob_branch->GetEntry(index);
			} else { 
				printf("branch svs_prob_branch does not exist!\n");
				exit(1);
			}
			svs_prob_isLoaded = true;
		}
		return svs_prob_;
	}
	vector<float> &svs_xError()
	{
		if (not svs_xError_isLoaded) {
			if (svs_xError_branch != 0) {
				svs_xError_branch->GetEntry(index);
			} else { 
				printf("branch svs_xError_branch does not exist!\n");
				exit(1);
			}
			svs_xError_isLoaded = true;
		}
		return svs_xError_;
	}
	vector<float> &svs_yError()
	{
		if (not svs_yError_isLoaded) {
			if (svs_yError_branch != 0) {
				svs_yError_branch->GetEntry(index);
			} else { 
				printf("branch svs_yError_branch does not exist!\n");
				exit(1);
			}
			svs_yError_isLoaded = true;
		}
		return svs_yError_;
	}
	vector<float> &svs_zError()
	{
		if (not svs_zError_isLoaded) {
			if (svs_zError_branch != 0) {
				svs_zError_branch->GetEntry(index);
			} else { 
				printf("branch svs_zError_branch does not exist!\n");
				exit(1);
			}
			svs_zError_isLoaded = true;
		}
		return svs_zError_;
	}
	vector<float> &mus_tcmet_deltax()
	{
		if (not mus_tcmet_deltax_isLoaded) {
			if (mus_tcmet_deltax_branch != 0) {
				mus_tcmet_deltax_branch->GetEntry(index);
			} else { 
				printf("branch mus_tcmet_deltax_branch does not exist!\n");
				exit(1);
			}
			mus_tcmet_deltax_isLoaded = true;
		}
		return mus_tcmet_deltax_;
	}
	vector<float> &mus_tcmet_deltay()
	{
		if (not mus_tcmet_deltay_isLoaded) {
			if (mus_tcmet_deltay_branch != 0) {
				mus_tcmet_deltay_branch->GetEntry(index);
			} else { 
				printf("branch mus_tcmet_deltay_branch does not exist!\n");
				exit(1);
			}
			mus_tcmet_deltay_isLoaded = true;
		}
		return mus_tcmet_deltay_;
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
	vector<float> &trks_d0phiCov()
	{
		if (not trks_d0phiCov_isLoaded) {
			if (trks_d0phiCov_branch != 0) {
				trks_d0phiCov_branch->GetEntry(index);
			} else { 
				printf("branch trks_d0phiCov_branch does not exist!\n");
				exit(1);
			}
			trks_d0phiCov_isLoaded = true;
		}
		return trks_d0phiCov_;
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
	vector<float> &trkjets_cor()
	{
		if (not trkjets_cor_isLoaded) {
			if (trkjets_cor_branch != 0) {
				trkjets_cor_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_cor_branch does not exist!\n");
				exit(1);
			}
			trkjets_cor_isLoaded = true;
		}
		return trkjets_cor_;
	}
	vector<float> &trks_d0Errvtx()
	{
		if (not trks_d0Errvtx_isLoaded) {
			if (trks_d0Errvtx_branch != 0) {
				trks_d0Errvtx_branch->GetEntry(index);
			} else { 
				printf("branch trks_d0Errvtx_branch does not exist!\n");
				exit(1);
			}
			trks_d0Errvtx_isLoaded = true;
		}
		return trks_d0Errvtx_;
	}
	vector<float> &trks_d0vtx()
	{
		if (not trks_d0vtx_isLoaded) {
			if (trks_d0vtx_branch != 0) {
				trks_d0vtx_branch->GetEntry(index);
			} else { 
				printf("branch trks_d0vtx_branch does not exist!\n");
				exit(1);
			}
			trks_d0vtx_isLoaded = true;
		}
		return trks_d0vtx_;
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
	vector<float> &vtxs_sumpt()
	{
		if (not vtxs_sumpt_isLoaded) {
			if (vtxs_sumpt_branch != 0) {
				vtxs_sumpt_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_sumpt_branch does not exist!\n");
				exit(1);
			}
			vtxs_sumpt_isLoaded = true;
		}
		return vtxs_sumpt_;
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
	vector<vector<float> > &vtxs_covMatrix()
	{
		if (not vtxs_covMatrix_isLoaded) {
			if (vtxs_covMatrix_branch != 0) {
				vtxs_covMatrix_branch->GetEntry(index);
			} else { 
				printf("branch vtxs_covMatrix_branch does not exist!\n");
				exit(1);
			}
			vtxs_covMatrix_isLoaded = true;
		}
		return vtxs_covMatrix_;
	}
	int &evt_cscLooseHaloId()
	{
		if (not evt_cscLooseHaloId_isLoaded) {
			if (evt_cscLooseHaloId_branch != 0) {
				evt_cscLooseHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_cscLooseHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_cscLooseHaloId_isLoaded = true;
		}
		return evt_cscLooseHaloId_;
	}
	int &evt_cscTightHaloId()
	{
		if (not evt_cscTightHaloId_isLoaded) {
			if (evt_cscTightHaloId_branch != 0) {
				evt_cscTightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_cscTightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_cscTightHaloId_isLoaded = true;
		}
		return evt_cscTightHaloId_;
	}
	int &evt_ecalLooseHaloId()
	{
		if (not evt_ecalLooseHaloId_isLoaded) {
			if (evt_ecalLooseHaloId_branch != 0) {
				evt_ecalLooseHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalLooseHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_ecalLooseHaloId_isLoaded = true;
		}
		return evt_ecalLooseHaloId_;
	}
	int &evt_ecalTightHaloId()
	{
		if (not evt_ecalTightHaloId_isLoaded) {
			if (evt_ecalTightHaloId_branch != 0) {
				evt_ecalTightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalTightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_ecalTightHaloId_isLoaded = true;
		}
		return evt_ecalTightHaloId_;
	}
	int &evt_extremeTightHaloId()
	{
		if (not evt_extremeTightHaloId_isLoaded) {
			if (evt_extremeTightHaloId_branch != 0) {
				evt_extremeTightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_extremeTightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_extremeTightHaloId_isLoaded = true;
		}
		return evt_extremeTightHaloId_;
	}
	int &evt_globalLooseHaloId()
	{
		if (not evt_globalLooseHaloId_isLoaded) {
			if (evt_globalLooseHaloId_branch != 0) {
				evt_globalLooseHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_globalLooseHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_globalLooseHaloId_isLoaded = true;
		}
		return evt_globalLooseHaloId_;
	}
	int &evt_globalTightHaloId()
	{
		if (not evt_globalTightHaloId_isLoaded) {
			if (evt_globalTightHaloId_branch != 0) {
				evt_globalTightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_globalTightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_globalTightHaloId_isLoaded = true;
		}
		return evt_globalTightHaloId_;
	}
	int &evt_hcalLooseHaloId()
	{
		if (not evt_hcalLooseHaloId_isLoaded) {
			if (evt_hcalLooseHaloId_branch != 0) {
				evt_hcalLooseHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalLooseHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_hcalLooseHaloId_isLoaded = true;
		}
		return evt_hcalLooseHaloId_;
	}
	int &evt_hcalTightHaloId()
	{
		if (not evt_hcalTightHaloId_isLoaded) {
			if (evt_hcalTightHaloId_branch != 0) {
				evt_hcalTightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcalTightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_hcalTightHaloId_isLoaded = true;
		}
		return evt_hcalTightHaloId_;
	}
	int &evt_looseHaloId()
	{
		if (not evt_looseHaloId_isLoaded) {
			if (evt_looseHaloId_branch != 0) {
				evt_looseHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_looseHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_looseHaloId_isLoaded = true;
		}
		return evt_looseHaloId_;
	}
	int &evt_nHaloLikeTracks()
	{
		if (not evt_nHaloLikeTracks_isLoaded) {
			if (evt_nHaloLikeTracks_branch != 0) {
				evt_nHaloLikeTracks_branch->GetEntry(index);
			} else { 
				printf("branch evt_nHaloLikeTracks_branch does not exist!\n");
				exit(1);
			}
			evt_nHaloLikeTracks_isLoaded = true;
		}
		return evt_nHaloLikeTracks_;
	}
	int &evt_nHaloTriggerCandidates()
	{
		if (not evt_nHaloTriggerCandidates_isLoaded) {
			if (evt_nHaloTriggerCandidates_branch != 0) {
				evt_nHaloTriggerCandidates_branch->GetEntry(index);
			} else { 
				printf("branch evt_nHaloTriggerCandidates_branch does not exist!\n");
				exit(1);
			}
			evt_nHaloTriggerCandidates_isLoaded = true;
		}
		return evt_nHaloTriggerCandidates_;
	}
	int &evt_tightHaloId()
	{
		if (not evt_tightHaloId_isLoaded) {
			if (evt_tightHaloId_branch != 0) {
				evt_tightHaloId_branch->GetEntry(index);
			} else { 
				printf("branch evt_tightHaloId_branch does not exist!\n");
				exit(1);
			}
			evt_tightHaloId_isLoaded = true;
		}
		return evt_tightHaloId_;
	}
	int &evt_bsType()
	{
		if (not evt_bsType_isLoaded) {
			if (evt_bsType_branch != 0) {
				evt_bsType_branch->GetEntry(index);
			} else { 
				printf("branch evt_bsType_branch does not exist!\n");
				exit(1);
			}
			evt_bsType_isLoaded = true;
		}
		return evt_bsType_;
	}
	int &evt_bunchCrossing()
	{
		if (not evt_bunchCrossing_isLoaded) {
			if (evt_bunchCrossing_branch != 0) {
				evt_bunchCrossing_branch->GetEntry(index);
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
			} else { 
				printf("branch evt_experimentType_branch does not exist!\n");
				exit(1);
			}
			evt_experimentType_isLoaded = true;
		}
		return evt_experimentType_;
	}
	int &evt_isRealData()
	{
		if (not evt_isRealData_isLoaded) {
			if (evt_isRealData_branch != 0) {
				evt_isRealData_branch->GetEntry(index);
			} else { 
				printf("branch evt_isRealData_branch does not exist!\n");
				exit(1);
			}
			evt_isRealData_isLoaded = true;
		}
		return evt_isRealData_;
	}
	int &evt_orbitNumber()
	{
		if (not evt_orbitNumber_isLoaded) {
			if (evt_orbitNumber_branch != 0) {
				evt_orbitNumber_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_nmus_branch does not exist!\n");
				exit(1);
			}
			l1_nmus_isLoaded = true;
		}
		return l1_nmus_;
	}
	int &pdfinfo_id1()
	{
		if (not pdfinfo_id1_isLoaded) {
			if (pdfinfo_id1_branch != 0) {
				pdfinfo_id1_branch->GetEntry(index);
			} else { 
				printf("branch pdfinfo_id1_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_id1_isLoaded = true;
		}
		return pdfinfo_id1_;
	}
	int &pdfinfo_id2()
	{
		if (not pdfinfo_id2_isLoaded) {
			if (pdfinfo_id2_branch != 0) {
				pdfinfo_id2_branch->GetEntry(index);
			} else { 
				printf("branch pdfinfo_id2_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_id2_isLoaded = true;
		}
		return pdfinfo_id2_;
	}
	vector<int> &evt_ecaliPhiSuspects()
	{
		if (not evt_ecaliPhiSuspects_isLoaded) {
			if (evt_ecaliPhiSuspects_branch != 0) {
				evt_ecaliPhiSuspects_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecaliPhiSuspects_branch does not exist!\n");
				exit(1);
			}
			evt_ecaliPhiSuspects_isLoaded = true;
		}
		return evt_ecaliPhiSuspects_;
	}
	vector<int> &evt_globaliPhiSuspects()
	{
		if (not evt_globaliPhiSuspects_isLoaded) {
			if (evt_globaliPhiSuspects_branch != 0) {
				evt_globaliPhiSuspects_branch->GetEntry(index);
			} else { 
				printf("branch evt_globaliPhiSuspects_branch does not exist!\n");
				exit(1);
			}
			evt_globaliPhiSuspects_isLoaded = true;
		}
		return evt_globaliPhiSuspects_;
	}
	vector<int> &evt_hcaliPhiSuspects()
	{
		if (not evt_hcaliPhiSuspects_isLoaded) {
			if (evt_hcaliPhiSuspects_branch != 0) {
				evt_hcaliPhiSuspects_branch->GetEntry(index);
			} else { 
				printf("branch evt_hcaliPhiSuspects_branch does not exist!\n");
				exit(1);
			}
			evt_hcaliPhiSuspects_isLoaded = true;
		}
		return evt_hcaliPhiSuspects_;
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
	vector<int> &els_mc3_motheridx()
	{
		if (not els_mc3_motheridx_isLoaded) {
			if (els_mc3_motheridx_branch != 0) {
				els_mc3_motheridx_branch->GetEntry(index);
			} else { 
				printf("branch els_mc3_motheridx_branch does not exist!\n");
				exit(1);
			}
			els_mc3_motheridx_isLoaded = true;
		}
		return els_mc3_motheridx_;
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
	vector<int> &jets_mc3_id()
	{
		if (not jets_mc3_id_isLoaded) {
			if (jets_mc3_id_branch != 0) {
				jets_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc3_id_branch does not exist!\n");
				exit(1);
			}
			jets_mc3_id_isLoaded = true;
		}
		return jets_mc3_id_;
	}
	vector<int> &jets_mc3idx()
	{
		if (not jets_mc3idx_isLoaded) {
			if (jets_mc3idx_branch != 0) {
				jets_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc3idx_branch does not exist!\n");
				exit(1);
			}
			jets_mc3idx_isLoaded = true;
		}
		return jets_mc3idx_;
	}
	vector<int> &jets_mc_gpidx()
	{
		if (not jets_mc_gpidx_isLoaded) {
			if (jets_mc_gpidx_branch != 0) {
				jets_mc_gpidx_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_gpidx_branch does not exist!\n");
				exit(1);
			}
			jets_mc_gpidx_isLoaded = true;
		}
		return jets_mc_gpidx_;
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
	vector<int> &jets_mcidx()
	{
		if (not jets_mcidx_isLoaded) {
			if (jets_mcidx_branch != 0) {
				jets_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch jets_mcidx_branch does not exist!\n");
				exit(1);
			}
			jets_mcidx_isLoaded = true;
		}
		return jets_mcidx_;
	}
	vector<int> &jets_mc_motherid()
	{
		if (not jets_mc_motherid_isLoaded) {
			if (jets_mc_motherid_branch != 0) {
				jets_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch jets_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			jets_mc_motherid_isLoaded = true;
		}
		return jets_mc_motherid_;
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
	vector<int> &mus_mc3_motheridx()
	{
		if (not mus_mc3_motheridx_isLoaded) {
			if (mus_mc3_motheridx_branch != 0) {
				mus_mc3_motheridx_branch->GetEntry(index);
			} else { 
				printf("branch mus_mc3_motheridx_branch does not exist!\n");
				exit(1);
			}
			mus_mc3_motheridx_isLoaded = true;
		}
		return mus_mc3_motheridx_;
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
	vector<int> &pfjets_mc3_id()
	{
		if (not pfjets_mc3_id_isLoaded) {
			if (pfjets_mc3_id_branch != 0) {
				pfjets_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc3_id_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc3_id_isLoaded = true;
		}
		return pfjets_mc3_id_;
	}
	vector<int> &pfjets_mc3idx()
	{
		if (not pfjets_mc3idx_isLoaded) {
			if (pfjets_mc3idx_branch != 0) {
				pfjets_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc3idx_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc3idx_isLoaded = true;
		}
		return pfjets_mc3idx_;
	}
	vector<int> &pfjets_mc_gpidx()
	{
		if (not pfjets_mc_gpidx_isLoaded) {
			if (pfjets_mc_gpidx_branch != 0) {
				pfjets_mc_gpidx_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_gpidx_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_gpidx_isLoaded = true;
		}
		return pfjets_mc_gpidx_;
	}
	vector<int> &pfjets_mc_id()
	{
		if (not pfjets_mc_id_isLoaded) {
			if (pfjets_mc_id_branch != 0) {
				pfjets_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_id_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_id_isLoaded = true;
		}
		return pfjets_mc_id_;
	}
	vector<int> &pfjets_mcidx()
	{
		if (not pfjets_mcidx_isLoaded) {
			if (pfjets_mcidx_branch != 0) {
				pfjets_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mcidx_branch does not exist!\n");
				exit(1);
			}
			pfjets_mcidx_isLoaded = true;
		}
		return pfjets_mcidx_;
	}
	vector<int> &pfjets_mc_motherid()
	{
		if (not pfjets_mc_motherid_isLoaded) {
			if (pfjets_mc_motherid_branch != 0) {
				pfjets_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			pfjets_mc_motherid_isLoaded = true;
		}
		return pfjets_mc_motherid_;
	}
	vector<int> &photons_mc3_id()
	{
		if (not photons_mc3_id_isLoaded) {
			if (photons_mc3_id_branch != 0) {
				photons_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc3_id_branch does not exist!\n");
				exit(1);
			}
			photons_mc3_id_isLoaded = true;
		}
		return photons_mc3_id_;
	}
	vector<int> &photons_mc3idx()
	{
		if (not photons_mc3idx_isLoaded) {
			if (photons_mc3idx_branch != 0) {
				photons_mc3idx_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc3idx_branch does not exist!\n");
				exit(1);
			}
			photons_mc3idx_isLoaded = true;
		}
		return photons_mc3idx_;
	}
	vector<int> &photons_mc3_motherid()
	{
		if (not photons_mc3_motherid_isLoaded) {
			if (photons_mc3_motherid_branch != 0) {
				photons_mc3_motherid_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc3_motherid_branch does not exist!\n");
				exit(1);
			}
			photons_mc3_motherid_isLoaded = true;
		}
		return photons_mc3_motherid_;
	}
	vector<int> &photons_mc3_motheridx()
	{
		if (not photons_mc3_motheridx_isLoaded) {
			if (photons_mc3_motheridx_branch != 0) {
				photons_mc3_motheridx_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc3_motheridx_branch does not exist!\n");
				exit(1);
			}
			photons_mc3_motheridx_isLoaded = true;
		}
		return photons_mc3_motheridx_;
	}
	vector<int> &photons_mc_id()
	{
		if (not photons_mc_id_isLoaded) {
			if (photons_mc_id_branch != 0) {
				photons_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc_id_branch does not exist!\n");
				exit(1);
			}
			photons_mc_id_isLoaded = true;
		}
		return photons_mc_id_;
	}
	vector<int> &photons_mcidx()
	{
		if (not photons_mcidx_isLoaded) {
			if (photons_mcidx_branch != 0) {
				photons_mcidx_branch->GetEntry(index);
			} else { 
				printf("branch photons_mcidx_branch does not exist!\n");
				exit(1);
			}
			photons_mcidx_isLoaded = true;
		}
		return photons_mcidx_;
	}
	vector<int> &photons_mc_motherid()
	{
		if (not photons_mc_motherid_isLoaded) {
			if (photons_mc_motherid_branch != 0) {
				photons_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch photons_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			photons_mc_motherid_isLoaded = true;
		}
		return photons_mc_motherid_;
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
	vector<int> &trk_mc3_motheridx()
	{
		if (not trk_mc3_motheridx_isLoaded) {
			if (trk_mc3_motheridx_branch != 0) {
				trk_mc3_motheridx_branch->GetEntry(index);
			} else { 
				printf("branch trk_mc3_motheridx_branch does not exist!\n");
				exit(1);
			}
			trk_mc3_motheridx_isLoaded = true;
		}
		return trk_mc3_motheridx_;
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
	vector<int> &els_exp_innerlayers39X()
	{
		if (not els_exp_innerlayers39X_isLoaded) {
			if (els_exp_innerlayers39X_branch != 0) {
				els_exp_innerlayers39X_branch->GetEntry(index);
			} else { 
				printf("branch els_exp_innerlayers39X_branch does not exist!\n");
				exit(1);
			}
			els_exp_innerlayers39X_isLoaded = true;
		}
		return els_exp_innerlayers39X_;
	}
	vector<int> &els_closestJet()
	{
		if (not els_closestJet_isLoaded) {
			if (els_closestJet_branch != 0) {
				els_closestJet_branch->GetEntry(index);
			} else { 
				printf("branch els_closestJet_branch does not exist!\n");
				exit(1);
			}
			els_closestJet_isLoaded = true;
		}
		return els_closestJet_;
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
	vector<int> &els_pfelsidx()
	{
		if (not els_pfelsidx_isLoaded) {
			if (els_pfelsidx_branch != 0) {
				els_pfelsidx_branch->GetEntry(index);
			} else { 
				printf("branch els_pfelsidx_branch does not exist!\n");
				exit(1);
			}
			els_pfelsidx_isLoaded = true;
		}
		return els_pfelsidx_;
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
	vector<int> &els_exp_innerlayers()
	{
		if (not els_exp_innerlayers_isLoaded) {
			if (els_exp_innerlayers_branch != 0) {
				els_exp_innerlayers_branch->GetEntry(index);
			} else { 
				printf("branch els_exp_innerlayers_branch does not exist!\n");
				exit(1);
			}
			els_exp_innerlayers_isLoaded = true;
		}
		return els_exp_innerlayers_;
	}
	vector<int> &els_exp_outerlayers()
	{
		if (not els_exp_outerlayers_isLoaded) {
			if (els_exp_outerlayers_branch != 0) {
				els_exp_outerlayers_branch->GetEntry(index);
			} else { 
				printf("branch els_exp_outerlayers_branch does not exist!\n");
				exit(1);
			}
			els_exp_outerlayers_isLoaded = true;
		}
		return els_exp_outerlayers_;
	}
	vector<int> &els_fiduciality()
	{
		if (not els_fiduciality_isLoaded) {
			if (els_fiduciality_branch != 0) {
				els_fiduciality_branch->GetEntry(index);
			} else { 
				printf("branch els_fiduciality_branch does not exist!\n");
				exit(1);
			}
			els_fiduciality_isLoaded = true;
		}
		return els_fiduciality_;
	}
	vector<int> &els_gsftrkidx()
	{
		if (not els_gsftrkidx_isLoaded) {
			if (els_gsftrkidx_branch != 0) {
				els_gsftrkidx_branch->GetEntry(index);
			} else { 
				printf("branch els_gsftrkidx_branch does not exist!\n");
				exit(1);
			}
			els_gsftrkidx_isLoaded = true;
		}
		return els_gsftrkidx_;
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
	vector<int> &els_sccharge()
	{
		if (not els_sccharge_isLoaded) {
			if (els_sccharge_branch != 0) {
				els_sccharge_branch->GetEntry(index);
			} else { 
				printf("branch els_sccharge_branch does not exist!\n");
				exit(1);
			}
			els_sccharge_isLoaded = true;
		}
		return els_sccharge_;
	}
	vector<int> &els_scindex()
	{
		if (not els_scindex_isLoaded) {
			if (els_scindex_branch != 0) {
				els_scindex_branch->GetEntry(index);
			} else { 
				printf("branch els_scindex_branch does not exist!\n");
				exit(1);
			}
			els_scindex_isLoaded = true;
		}
		return els_scindex_;
	}
	vector<int> &els_trk_charge()
	{
		if (not els_trk_charge_isLoaded) {
			if (els_trk_charge_branch != 0) {
				els_trk_charge_branch->GetEntry(index);
			} else { 
				printf("branch els_trk_charge_branch does not exist!\n");
				exit(1);
			}
			els_trk_charge_isLoaded = true;
		}
		return els_trk_charge_;
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
	vector<int> &els_type()
	{
		if (not els_type_isLoaded) {
			if (els_type_branch != 0) {
				els_type_branch->GetEntry(index);
			} else { 
				printf("branch els_type_branch does not exist!\n");
				exit(1);
			}
			els_type_isLoaded = true;
		}
		return els_type_;
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
	vector<int> &gsftrks_charge()
	{
		if (not gsftrks_charge_isLoaded) {
			if (gsftrks_charge_branch != 0) {
				gsftrks_charge_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_charge_branch does not exist!\n");
				exit(1);
			}
			gsftrks_charge_isLoaded = true;
		}
		return gsftrks_charge_;
	}
	vector<int> &gsftrks_exp_innerlayers()
	{
		if (not gsftrks_exp_innerlayers_isLoaded) {
			if (gsftrks_exp_innerlayers_branch != 0) {
				gsftrks_exp_innerlayers_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_exp_innerlayers_branch does not exist!\n");
				exit(1);
			}
			gsftrks_exp_innerlayers_isLoaded = true;
		}
		return gsftrks_exp_innerlayers_;
	}
	vector<int> &gsftrks_exp_outerlayers()
	{
		if (not gsftrks_exp_outerlayers_isLoaded) {
			if (gsftrks_exp_outerlayers_branch != 0) {
				gsftrks_exp_outerlayers_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_exp_outerlayers_branch does not exist!\n");
				exit(1);
			}
			gsftrks_exp_outerlayers_isLoaded = true;
		}
		return gsftrks_exp_outerlayers_;
	}
	vector<int> &gsftrks_layer1_det()
	{
		if (not gsftrks_layer1_det_isLoaded) {
			if (gsftrks_layer1_det_branch != 0) {
				gsftrks_layer1_det_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_layer1_det_branch does not exist!\n");
				exit(1);
			}
			gsftrks_layer1_det_isLoaded = true;
		}
		return gsftrks_layer1_det_;
	}
	vector<int> &gsftrks_layer1_layer()
	{
		if (not gsftrks_layer1_layer_isLoaded) {
			if (gsftrks_layer1_layer_branch != 0) {
				gsftrks_layer1_layer_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_layer1_layer_branch does not exist!\n");
				exit(1);
			}
			gsftrks_layer1_layer_isLoaded = true;
		}
		return gsftrks_layer1_layer_;
	}
	vector<int> &gsftrks_layer1_sizerphi()
	{
		if (not gsftrks_layer1_sizerphi_isLoaded) {
			if (gsftrks_layer1_sizerphi_branch != 0) {
				gsftrks_layer1_sizerphi_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_layer1_sizerphi_branch does not exist!\n");
				exit(1);
			}
			gsftrks_layer1_sizerphi_isLoaded = true;
		}
		return gsftrks_layer1_sizerphi_;
	}
	vector<int> &gsftrks_layer1_sizerz()
	{
		if (not gsftrks_layer1_sizerz_isLoaded) {
			if (gsftrks_layer1_sizerz_branch != 0) {
				gsftrks_layer1_sizerz_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_layer1_sizerz_branch does not exist!\n");
				exit(1);
			}
			gsftrks_layer1_sizerz_isLoaded = true;
		}
		return gsftrks_layer1_sizerz_;
	}
	vector<int> &gsftrks_lostHits()
	{
		if (not gsftrks_lostHits_isLoaded) {
			if (gsftrks_lostHits_branch != 0) {
				gsftrks_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_lostHits_branch does not exist!\n");
				exit(1);
			}
			gsftrks_lostHits_isLoaded = true;
		}
		return gsftrks_lostHits_;
	}
	vector<int> &gsftrks_lost_pixelhits()
	{
		if (not gsftrks_lost_pixelhits_isLoaded) {
			if (gsftrks_lost_pixelhits_branch != 0) {
				gsftrks_lost_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_lost_pixelhits_branch does not exist!\n");
				exit(1);
			}
			gsftrks_lost_pixelhits_isLoaded = true;
		}
		return gsftrks_lost_pixelhits_;
	}
	vector<int> &gsftrks_nlayers()
	{
		if (not gsftrks_nlayers_isLoaded) {
			if (gsftrks_nlayers_branch != 0) {
				gsftrks_nlayers_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_nlayers_branch does not exist!\n");
				exit(1);
			}
			gsftrks_nlayers_isLoaded = true;
		}
		return gsftrks_nlayers_;
	}
	vector<int> &gsftrks_nlayers3D()
	{
		if (not gsftrks_nlayers3D_isLoaded) {
			if (gsftrks_nlayers3D_branch != 0) {
				gsftrks_nlayers3D_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_nlayers3D_branch does not exist!\n");
				exit(1);
			}
			gsftrks_nlayers3D_isLoaded = true;
		}
		return gsftrks_nlayers3D_;
	}
	vector<int> &gsftrks_nlayersLost()
	{
		if (not gsftrks_nlayersLost_isLoaded) {
			if (gsftrks_nlayersLost_branch != 0) {
				gsftrks_nlayersLost_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_nlayersLost_branch does not exist!\n");
				exit(1);
			}
			gsftrks_nlayersLost_isLoaded = true;
		}
		return gsftrks_nlayersLost_;
	}
	vector<int> &gsftrks_validHits()
	{
		if (not gsftrks_validHits_isLoaded) {
			if (gsftrks_validHits_branch != 0) {
				gsftrks_validHits_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_validHits_branch does not exist!\n");
				exit(1);
			}
			gsftrks_validHits_isLoaded = true;
		}
		return gsftrks_validHits_;
	}
	vector<int> &gsftrks_valid_pixelhits()
	{
		if (not gsftrks_valid_pixelhits_isLoaded) {
			if (gsftrks_valid_pixelhits_branch != 0) {
				gsftrks_valid_pixelhits_branch->GetEntry(index);
			} else { 
				printf("branch gsftrks_valid_pixelhits_branch does not exist!\n");
				exit(1);
			}
			gsftrks_valid_pixelhits_isLoaded = true;
		}
		return gsftrks_valid_pixelhits_;
	}
	vector<int> &hyp_ll_charge()
	{
		if (not hyp_ll_charge_isLoaded) {
			if (hyp_ll_charge_branch != 0) {
				hyp_ll_charge_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_charge_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_charge_isLoaded = true;
		}
		return hyp_ll_charge_;
	}
	vector<int> &hyp_ll_id()
	{
		if (not hyp_ll_id_isLoaded) {
			if (hyp_ll_id_branch != 0) {
				hyp_ll_id_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_id_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_id_isLoaded = true;
		}
		return hyp_ll_id_;
	}
	vector<int> &hyp_ll_index()
	{
		if (not hyp_ll_index_isLoaded) {
			if (hyp_ll_index_branch != 0) {
				hyp_ll_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_index_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_index_isLoaded = true;
		}
		return hyp_ll_index_;
	}
	vector<int> &hyp_ll_lostHits()
	{
		if (not hyp_ll_lostHits_isLoaded) {
			if (hyp_ll_lostHits_branch != 0) {
				hyp_ll_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_lostHits_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_lostHits_isLoaded = true;
		}
		return hyp_ll_lostHits_;
	}
	vector<int> &hyp_ll_validHits()
	{
		if (not hyp_ll_validHits_isLoaded) {
			if (hyp_ll_validHits_branch != 0) {
				hyp_ll_validHits_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_validHits_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_validHits_isLoaded = true;
		}
		return hyp_ll_validHits_;
	}
	vector<int> &hyp_lt_charge()
	{
		if (not hyp_lt_charge_isLoaded) {
			if (hyp_lt_charge_branch != 0) {
				hyp_lt_charge_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_charge_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_charge_isLoaded = true;
		}
		return hyp_lt_charge_;
	}
	vector<int> &hyp_lt_id()
	{
		if (not hyp_lt_id_isLoaded) {
			if (hyp_lt_id_branch != 0) {
				hyp_lt_id_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_id_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_id_isLoaded = true;
		}
		return hyp_lt_id_;
	}
	vector<int> &hyp_lt_index()
	{
		if (not hyp_lt_index_isLoaded) {
			if (hyp_lt_index_branch != 0) {
				hyp_lt_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_index_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_index_isLoaded = true;
		}
		return hyp_lt_index_;
	}
	vector<int> &hyp_lt_lostHits()
	{
		if (not hyp_lt_lostHits_isLoaded) {
			if (hyp_lt_lostHits_branch != 0) {
				hyp_lt_lostHits_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_lostHits_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_lostHits_isLoaded = true;
		}
		return hyp_lt_lostHits_;
	}
	vector<int> &hyp_lt_validHits()
	{
		if (not hyp_lt_validHits_isLoaded) {
			if (hyp_lt_validHits_branch != 0) {
				hyp_lt_validHits_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_validHits_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_validHits_isLoaded = true;
		}
		return hyp_lt_validHits_;
	}
	vector<int> &hyp_njets()
	{
		if (not hyp_njets_isLoaded) {
			if (hyp_njets_branch != 0) {
				hyp_njets_branch->GetEntry(index);
			} else { 
				printf("branch hyp_njets_branch does not exist!\n");
				exit(1);
			}
			hyp_njets_isLoaded = true;
		}
		return hyp_njets_;
	}
	vector<int> &hyp_nojets()
	{
		if (not hyp_nojets_isLoaded) {
			if (hyp_nojets_branch != 0) {
				hyp_nojets_branch->GetEntry(index);
			} else { 
				printf("branch hyp_nojets_branch does not exist!\n");
				exit(1);
			}
			hyp_nojets_isLoaded = true;
		}
		return hyp_nojets_;
	}
	vector<int> &hyp_type()
	{
		if (not hyp_type_isLoaded) {
			if (hyp_type_branch != 0) {
				hyp_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_type_branch does not exist!\n");
				exit(1);
			}
			hyp_type_isLoaded = true;
		}
		return hyp_type_;
	}
	vector<int> &hyp_FVFit_ndf()
	{
		if (not hyp_FVFit_ndf_isLoaded) {
			if (hyp_FVFit_ndf_branch != 0) {
				hyp_FVFit_ndf_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_ndf_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_ndf_isLoaded = true;
		}
		return hyp_FVFit_ndf_;
	}
	vector<int> &hyp_FVFit_status()
	{
		if (not hyp_FVFit_status_isLoaded) {
			if (hyp_FVFit_status_branch != 0) {
				hyp_FVFit_status_branch->GetEntry(index);
			} else { 
				printf("branch hyp_FVFit_status_branch does not exist!\n");
				exit(1);
			}
			hyp_FVFit_status_isLoaded = true;
		}
		return hyp_FVFit_status_;
	}
	vector<int> &hyp_ll_mc_id()
	{
		if (not hyp_ll_mc_id_isLoaded) {
			if (hyp_ll_mc_id_branch != 0) {
				hyp_ll_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_mc_id_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_mc_id_isLoaded = true;
		}
		return hyp_ll_mc_id_;
	}
	vector<int> &hyp_ll_mc_motherid()
	{
		if (not hyp_ll_mc_motherid_isLoaded) {
			if (hyp_ll_mc_motherid_branch != 0) {
				hyp_ll_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch hyp_ll_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_mc_motherid_isLoaded = true;
		}
		return hyp_ll_mc_motherid_;
	}
	vector<int> &hyp_lt_mc_id()
	{
		if (not hyp_lt_mc_id_isLoaded) {
			if (hyp_lt_mc_id_branch != 0) {
				hyp_lt_mc_id_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_mc_id_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_mc_id_isLoaded = true;
		}
		return hyp_lt_mc_id_;
	}
	vector<int> &hyp_lt_mc_motherid()
	{
		if (not hyp_lt_mc_motherid_isLoaded) {
			if (hyp_lt_mc_motherid_branch != 0) {
				hyp_lt_mc_motherid_branch->GetEntry(index);
			} else { 
				printf("branch hyp_lt_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_mc_motherid_isLoaded = true;
		}
		return hyp_lt_mc_motherid_;
	}
	vector<int> &hyp_quadlep_first_type()
	{
		if (not hyp_quadlep_first_type_isLoaded) {
			if (hyp_quadlep_first_type_branch != 0) {
				hyp_quadlep_first_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_first_type_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_first_type_isLoaded = true;
		}
		return hyp_quadlep_first_type_;
	}
	vector<int> &hyp_quadlep_fourth_type()
	{
		if (not hyp_quadlep_fourth_type_isLoaded) {
			if (hyp_quadlep_fourth_type_branch != 0) {
				hyp_quadlep_fourth_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_fourth_type_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_fourth_type_isLoaded = true;
		}
		return hyp_quadlep_fourth_type_;
	}
	vector<int> &hyp_quadlep_second_type()
	{
		if (not hyp_quadlep_second_type_isLoaded) {
			if (hyp_quadlep_second_type_branch != 0) {
				hyp_quadlep_second_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_second_type_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_second_type_isLoaded = true;
		}
		return hyp_quadlep_second_type_;
	}
	vector<int> &hyp_quadlep_third_type()
	{
		if (not hyp_quadlep_third_type_isLoaded) {
			if (hyp_quadlep_third_type_branch != 0) {
				hyp_quadlep_third_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_third_type_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_third_type_isLoaded = true;
		}
		return hyp_quadlep_third_type_;
	}
	vector<int> &hyp_trilep_first_type()
	{
		if (not hyp_trilep_first_type_isLoaded) {
			if (hyp_trilep_first_type_branch != 0) {
				hyp_trilep_first_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_first_type_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_first_type_isLoaded = true;
		}
		return hyp_trilep_first_type_;
	}
	vector<int> &hyp_trilep_second_type()
	{
		if (not hyp_trilep_second_type_isLoaded) {
			if (hyp_trilep_second_type_branch != 0) {
				hyp_trilep_second_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_second_type_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_second_type_isLoaded = true;
		}
		return hyp_trilep_second_type_;
	}
	vector<int> &hyp_trilep_third_type()
	{
		if (not hyp_trilep_third_type_isLoaded) {
			if (hyp_trilep_third_type_branch != 0) {
				hyp_trilep_third_type_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_third_type_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_third_type_isLoaded = true;
		}
		return hyp_trilep_third_type_;
	}
	vector<int> &jets_closestElectron()
	{
		if (not jets_closestElectron_isLoaded) {
			if (jets_closestElectron_branch != 0) {
				jets_closestElectron_branch->GetEntry(index);
			} else { 
				printf("branch jets_closestElectron_branch does not exist!\n");
				exit(1);
			}
			jets_closestElectron_isLoaded = true;
		}
		return jets_closestElectron_;
	}
	vector<int> &jets_closestMuon()
	{
		if (not jets_closestMuon_isLoaded) {
			if (jets_closestMuon_branch != 0) {
				jets_closestMuon_branch->GetEntry(index);
			} else { 
				printf("branch jets_closestMuon_branch does not exist!\n");
				exit(1);
			}
			jets_closestMuon_isLoaded = true;
		}
		return jets_closestMuon_;
	}
	vector<int> &l1_emiso_ieta()
	{
		if (not l1_emiso_ieta_isLoaded) {
			if (l1_emiso_ieta_branch != 0) {
				l1_emiso_ieta_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_mus_qualFlags_branch does not exist!\n");
				exit(1);
			}
			l1_mus_qualFlags_isLoaded = true;
		}
		return l1_mus_qualFlags_;
	}
	vector<int> &mus_met_flag()
	{
		if (not mus_met_flag_isLoaded) {
			if (mus_met_flag_branch != 0) {
				mus_met_flag_branch->GetEntry(index);
			} else { 
				printf("branch mus_met_flag_branch does not exist!\n");
				exit(1);
			}
			mus_met_flag_isLoaded = true;
		}
		return mus_met_flag_;
	}
	vector<int> &mus_closestEle()
	{
		if (not mus_closestEle_isLoaded) {
			if (mus_closestEle_branch != 0) {
				mus_closestEle_branch->GetEntry(index);
			} else { 
				printf("branch mus_closestEle_branch does not exist!\n");
				exit(1);
			}
			mus_closestEle_isLoaded = true;
		}
		return mus_closestEle_;
	}
	vector<int> &mus_closestJet()
	{
		if (not mus_closestJet_isLoaded) {
			if (mus_closestJet_branch != 0) {
				mus_closestJet_branch->GetEntry(index);
			} else { 
				printf("branch mus_closestJet_branch does not exist!\n");
				exit(1);
			}
			mus_closestJet_isLoaded = true;
		}
		return mus_closestJet_;
	}
	vector<int> &mus_pfmusidx()
	{
		if (not mus_pfmusidx_isLoaded) {
			if (mus_pfmusidx_branch != 0) {
				mus_pfmusidx_branch->GetEntry(index);
			} else { 
				printf("branch mus_pfmusidx_branch does not exist!\n");
				exit(1);
			}
			mus_pfmusidx_isLoaded = true;
		}
		return mus_pfmusidx_;
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
	vector<int> &mus_chi2LocalMomentum()
	{
		if (not mus_chi2LocalMomentum_isLoaded) {
			if (mus_chi2LocalMomentum_branch != 0) {
				mus_chi2LocalMomentum_branch->GetEntry(index);
			} else { 
				printf("branch mus_chi2LocalMomentum_branch does not exist!\n");
				exit(1);
			}
			mus_chi2LocalMomentum_isLoaded = true;
		}
		return mus_chi2LocalMomentum_;
	}
	vector<int> &mus_chi2LocalPosition()
	{
		if (not mus_chi2LocalPosition_isLoaded) {
			if (mus_chi2LocalPosition_branch != 0) {
				mus_chi2LocalPosition_branch->GetEntry(index);
			} else { 
				printf("branch mus_chi2LocalPosition_branch does not exist!\n");
				exit(1);
			}
			mus_chi2LocalPosition_isLoaded = true;
		}
		return mus_chi2LocalPosition_;
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
	vector<int> &mus_gfit_validSTAHits()
	{
		if (not mus_gfit_validSTAHits_isLoaded) {
			if (mus_gfit_validSTAHits_branch != 0) {
				mus_gfit_validSTAHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_validSTAHits_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_validSTAHits_isLoaded = true;
		}
		return mus_gfit_validSTAHits_;
	}
	vector<int> &mus_gfit_validSiHits()
	{
		if (not mus_gfit_validSiHits_isLoaded) {
			if (mus_gfit_validSiHits_branch != 0) {
				mus_gfit_validSiHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_gfit_validSiHits_branch does not exist!\n");
				exit(1);
			}
			mus_gfit_validSiHits_isLoaded = true;
		}
		return mus_gfit_validSiHits_;
	}
	vector<int> &mus_glbKink()
	{
		if (not mus_glbKink_isLoaded) {
			if (mus_glbKink_branch != 0) {
				mus_glbKink_branch->GetEntry(index);
			} else { 
				printf("branch mus_glbKink_branch does not exist!\n");
				exit(1);
			}
			mus_glbKink_isLoaded = true;
		}
		return mus_glbKink_;
	}
	vector<int> &mus_glbTrackProbability()
	{
		if (not mus_glbTrackProbability_isLoaded) {
			if (mus_glbTrackProbability_branch != 0) {
				mus_glbTrackProbability_branch->GetEntry(index);
			} else { 
				printf("branch mus_glbTrackProbability_branch does not exist!\n");
				exit(1);
			}
			mus_glbTrackProbability_isLoaded = true;
		}
		return mus_glbTrackProbability_;
	}
	vector<int> &mus_globalDeltaEtaPhi()
	{
		if (not mus_globalDeltaEtaPhi_isLoaded) {
			if (mus_globalDeltaEtaPhi_branch != 0) {
				mus_globalDeltaEtaPhi_branch->GetEntry(index);
			} else { 
				printf("branch mus_globalDeltaEtaPhi_branch does not exist!\n");
				exit(1);
			}
			mus_globalDeltaEtaPhi_isLoaded = true;
		}
		return mus_globalDeltaEtaPhi_;
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
	vector<int> &mus_localDistance()
	{
		if (not mus_localDistance_isLoaded) {
			if (mus_localDistance_branch != 0) {
				mus_localDistance_branch->GetEntry(index);
			} else { 
				printf("branch mus_localDistance_branch does not exist!\n");
				exit(1);
			}
			mus_localDistance_isLoaded = true;
		}
		return mus_localDistance_;
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
	vector<int> &mus_nOverlaps()
	{
		if (not mus_nOverlaps_isLoaded) {
			if (mus_nOverlaps_branch != 0) {
				mus_nOverlaps_branch->GetEntry(index);
			} else { 
				printf("branch mus_nOverlaps_branch does not exist!\n");
				exit(1);
			}
			mus_nOverlaps_isLoaded = true;
		}
		return mus_nOverlaps_;
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
	vector<int> &mus_overlap0()
	{
		if (not mus_overlap0_isLoaded) {
			if (mus_overlap0_branch != 0) {
				mus_overlap0_branch->GetEntry(index);
			} else { 
				printf("branch mus_overlap0_branch does not exist!\n");
				exit(1);
			}
			mus_overlap0_isLoaded = true;
		}
		return mus_overlap0_;
	}
	vector<int> &mus_overlap1()
	{
		if (not mus_overlap1_isLoaded) {
			if (mus_overlap1_branch != 0) {
				mus_overlap1_branch->GetEntry(index);
			} else { 
				printf("branch mus_overlap1_branch does not exist!\n");
				exit(1);
			}
			mus_overlap1_isLoaded = true;
		}
		return mus_overlap1_;
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
	vector<int> &mus_staRelChi2()
	{
		if (not mus_staRelChi2_isLoaded) {
			if (mus_staRelChi2_branch != 0) {
				mus_staRelChi2_branch->GetEntry(index);
			} else { 
				printf("branch mus_staRelChi2_branch does not exist!\n");
				exit(1);
			}
			mus_staRelChi2_isLoaded = true;
		}
		return mus_staRelChi2_;
	}
	vector<int> &mus_sta_validHits()
	{
		if (not mus_sta_validHits_isLoaded) {
			if (mus_sta_validHits_branch != 0) {
				mus_sta_validHits_branch->GetEntry(index);
			} else { 
				printf("branch mus_sta_validHits_branch does not exist!\n");
				exit(1);
			}
			mus_sta_validHits_isLoaded = true;
		}
		return mus_sta_validHits_;
	}
	vector<int> &mus_timeDirection()
	{
		if (not mus_timeDirection_isLoaded) {
			if (mus_timeDirection_branch != 0) {
				mus_timeDirection_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeDirection_branch does not exist!\n");
				exit(1);
			}
			mus_timeDirection_isLoaded = true;
		}
		return mus_timeDirection_;
	}
	vector<int> &mus_timeNumStationsUsed()
	{
		if (not mus_timeNumStationsUsed_isLoaded) {
			if (mus_timeNumStationsUsed_branch != 0) {
				mus_timeNumStationsUsed_branch->GetEntry(index);
			} else { 
				printf("branch mus_timeNumStationsUsed_branch does not exist!\n");
				exit(1);
			}
			mus_timeNumStationsUsed_isLoaded = true;
		}
		return mus_timeNumStationsUsed_;
	}
	vector<int> &mus_trkKink()
	{
		if (not mus_trkKink_isLoaded) {
			if (mus_trkKink_branch != 0) {
				mus_trkKink_branch->GetEntry(index);
			} else { 
				printf("branch mus_trkKink_branch does not exist!\n");
				exit(1);
			}
			mus_trkKink_isLoaded = true;
		}
		return mus_trkKink_;
	}
	vector<int> &mus_trkRelChi2()
	{
		if (not mus_trkRelChi2_isLoaded) {
			if (mus_trkRelChi2_branch != 0) {
				mus_trkRelChi2_branch->GetEntry(index);
			} else { 
				printf("branch mus_trkRelChi2_branch does not exist!\n");
				exit(1);
			}
			mus_trkRelChi2_isLoaded = true;
		}
		return mus_trkRelChi2_;
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
	vector<int> &mus_trkidx()
	{
		if (not mus_trkidx_isLoaded) {
			if (mus_trkidx_branch != 0) {
				mus_trkidx_branch->GetEntry(index);
			} else { 
				printf("branch mus_trkidx_branch does not exist!\n");
				exit(1);
			}
			mus_trkidx_isLoaded = true;
		}
		return mus_trkidx_;
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
	vector<int> &els_pat_genID()
	{
		if (not els_pat_genID_isLoaded) {
			if (els_pat_genID_branch != 0) {
				els_pat_genID_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_genID_branch does not exist!\n");
				exit(1);
			}
			els_pat_genID_isLoaded = true;
		}
		return els_pat_genID_;
	}
	vector<int> &els_pat_genMotherID()
	{
		if (not els_pat_genMotherID_isLoaded) {
			if (els_pat_genMotherID_branch != 0) {
				els_pat_genMotherID_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_genMotherID_branch does not exist!\n");
				exit(1);
			}
			els_pat_genMotherID_isLoaded = true;
		}
		return els_pat_genMotherID_;
	}
	vector<int> &jets_pat_genPartonMother_id()
	{
		if (not jets_pat_genPartonMother_id_isLoaded) {
			if (jets_pat_genPartonMother_id_branch != 0) {
				jets_pat_genPartonMother_id_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_genPartonMother_id_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genPartonMother_id_isLoaded = true;
		}
		return jets_pat_genPartonMother_id_;
	}
	vector<int> &jets_pat_genParton_id()
	{
		if (not jets_pat_genParton_id_isLoaded) {
			if (jets_pat_genParton_id_branch != 0) {
				jets_pat_genParton_id_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_genParton_id_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genParton_id_isLoaded = true;
		}
		return jets_pat_genParton_id_;
	}
	vector<int> &jets_pat_jetIDLoose()
	{
		if (not jets_pat_jetIDLoose_isLoaded) {
			if (jets_pat_jetIDLoose_branch != 0) {
				jets_pat_jetIDLoose_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetIDLoose_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetIDLoose_isLoaded = true;
		}
		return jets_pat_jetIDLoose_;
	}
	vector<int> &jets_pat_jetIDLooseAOD()
	{
		if (not jets_pat_jetIDLooseAOD_isLoaded) {
			if (jets_pat_jetIDLooseAOD_branch != 0) {
				jets_pat_jetIDLooseAOD_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetIDLooseAOD_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetIDLooseAOD_isLoaded = true;
		}
		return jets_pat_jetIDLooseAOD_;
	}
	vector<int> &jets_pat_jetIDMinimal()
	{
		if (not jets_pat_jetIDMinimal_isLoaded) {
			if (jets_pat_jetIDMinimal_branch != 0) {
				jets_pat_jetIDMinimal_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetIDMinimal_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetIDMinimal_isLoaded = true;
		}
		return jets_pat_jetIDMinimal_;
	}
	vector<int> &jets_pat_jetIDTight()
	{
		if (not jets_pat_jetIDTight_isLoaded) {
			if (jets_pat_jetIDTight_branch != 0) {
				jets_pat_jetIDTight_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_jetIDTight_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jetIDTight_isLoaded = true;
		}
		return jets_pat_jetIDTight_;
	}
	vector<int> &jets_pat_partonFlavour()
	{
		if (not jets_pat_partonFlavour_isLoaded) {
			if (jets_pat_partonFlavour_branch != 0) {
				jets_pat_partonFlavour_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_partonFlavour_branch does not exist!\n");
				exit(1);
			}
			jets_pat_partonFlavour_isLoaded = true;
		}
		return jets_pat_partonFlavour_;
	}
	vector<int> &mus_pat_genID()
	{
		if (not mus_pat_genID_isLoaded) {
			if (mus_pat_genID_branch != 0) {
				mus_pat_genID_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_genID_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genID_isLoaded = true;
		}
		return mus_pat_genID_;
	}
	vector<int> &mus_pat_genMotherID()
	{
		if (not mus_pat_genMotherID_isLoaded) {
			if (mus_pat_genMotherID_branch != 0) {
				mus_pat_genMotherID_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_genMotherID_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genMotherID_isLoaded = true;
		}
		return mus_pat_genMotherID_;
	}
	vector<int> &pfels_elsidx()
	{
		if (not pfels_elsidx_isLoaded) {
			if (pfels_elsidx_branch != 0) {
				pfels_elsidx_branch->GetEntry(index);
			} else { 
				printf("branch pfels_elsidx_branch does not exist!\n");
				exit(1);
			}
			pfels_elsidx_isLoaded = true;
		}
		return pfels_elsidx_;
	}
	vector<int> &pfels_charge()
	{
		if (not pfels_charge_isLoaded) {
			if (pfels_charge_branch != 0) {
				pfels_charge_branch->GetEntry(index);
			} else { 
				printf("branch pfels_charge_branch does not exist!\n");
				exit(1);
			}
			pfels_charge_isLoaded = true;
		}
		return pfels_charge_;
	}
	vector<int> &pfels_flag()
	{
		if (not pfels_flag_isLoaded) {
			if (pfels_flag_branch != 0) {
				pfels_flag_branch->GetEntry(index);
			} else { 
				printf("branch pfels_flag_branch does not exist!\n");
				exit(1);
			}
			pfels_flag_isLoaded = true;
		}
		return pfels_flag_;
	}
	vector<int> &pfels_particleId()
	{
		if (not pfels_particleId_isLoaded) {
			if (pfels_particleId_branch != 0) {
				pfels_particleId_branch->GetEntry(index);
			} else { 
				printf("branch pfels_particleId_branch does not exist!\n");
				exit(1);
			}
			pfels_particleId_isLoaded = true;
		}
		return pfels_particleId_;
	}
	vector<int> &pfjets_chargedMultiplicity()
	{
		if (not pfjets_chargedMultiplicity_isLoaded) {
			if (pfjets_chargedMultiplicity_branch != 0) {
				pfjets_chargedMultiplicity_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_chargedMultiplicity_branch does not exist!\n");
				exit(1);
			}
			pfjets_chargedMultiplicity_isLoaded = true;
		}
		return pfjets_chargedMultiplicity_;
	}
	vector<int> &pfjets_muonMultiplicity()
	{
		if (not pfjets_muonMultiplicity_isLoaded) {
			if (pfjets_muonMultiplicity_branch != 0) {
				pfjets_muonMultiplicity_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_muonMultiplicity_branch does not exist!\n");
				exit(1);
			}
			pfjets_muonMultiplicity_isLoaded = true;
		}
		return pfjets_muonMultiplicity_;
	}
	vector<int> &pfjets_neutralMultiplicity()
	{
		if (not pfjets_neutralMultiplicity_isLoaded) {
			if (pfjets_neutralMultiplicity_branch != 0) {
				pfjets_neutralMultiplicity_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_neutralMultiplicity_branch does not exist!\n");
				exit(1);
			}
			pfjets_neutralMultiplicity_isLoaded = true;
		}
		return pfjets_neutralMultiplicity_;
	}
	vector<int> &pfmus_musidx()
	{
		if (not pfmus_musidx_isLoaded) {
			if (pfmus_musidx_branch != 0) {
				pfmus_musidx_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_musidx_branch does not exist!\n");
				exit(1);
			}
			pfmus_musidx_isLoaded = true;
		}
		return pfmus_musidx_;
	}
	vector<int> &pfmus_charge()
	{
		if (not pfmus_charge_isLoaded) {
			if (pfmus_charge_branch != 0) {
				pfmus_charge_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_charge_branch does not exist!\n");
				exit(1);
			}
			pfmus_charge_isLoaded = true;
		}
		return pfmus_charge_;
	}
	vector<int> &pfmus_flag()
	{
		if (not pfmus_flag_isLoaded) {
			if (pfmus_flag_branch != 0) {
				pfmus_flag_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_flag_branch does not exist!\n");
				exit(1);
			}
			pfmus_flag_isLoaded = true;
		}
		return pfmus_flag_;
	}
	vector<int> &pfmus_particleId()
	{
		if (not pfmus_particleId_isLoaded) {
			if (pfmus_particleId_branch != 0) {
				pfmus_particleId_branch->GetEntry(index);
			} else { 
				printf("branch pfmus_particleId_branch does not exist!\n");
				exit(1);
			}
			pfmus_particleId_isLoaded = true;
		}
		return pfmus_particleId_;
	}
	vector<int> &photons_fiduciality()
	{
		if (not photons_fiduciality_isLoaded) {
			if (photons_fiduciality_branch != 0) {
				photons_fiduciality_branch->GetEntry(index);
			} else { 
				printf("branch photons_fiduciality_branch does not exist!\n");
				exit(1);
			}
			photons_fiduciality_isLoaded = true;
		}
		return photons_fiduciality_;
	}
	vector<int> &photons_scindex()
	{
		if (not photons_scindex_isLoaded) {
			if (photons_scindex_branch != 0) {
				photons_scindex_branch->GetEntry(index);
			} else { 
				printf("branch photons_scindex_branch does not exist!\n");
				exit(1);
			}
			photons_scindex_isLoaded = true;
		}
		return photons_scindex_;
	}
	vector<int> &scs_detIdSeed()
	{
		if (not scs_detIdSeed_isLoaded) {
			if (scs_detIdSeed_branch != 0) {
				scs_detIdSeed_branch->GetEntry(index);
			} else { 
				printf("branch scs_detIdSeed_branch does not exist!\n");
				exit(1);
			}
			scs_detIdSeed_isLoaded = true;
		}
		return scs_detIdSeed_;
	}
	vector<int> &scs_elsidx()
	{
		if (not scs_elsidx_isLoaded) {
			if (scs_elsidx_branch != 0) {
				scs_elsidx_branch->GetEntry(index);
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
			} else { 
				printf("branch scs_severitySeed_branch does not exist!\n");
				exit(1);
			}
			scs_severitySeed_isLoaded = true;
		}
		return scs_severitySeed_;
	}
	vector<int> &svs_isKs()
	{
		if (not svs_isKs_isLoaded) {
			if (svs_isKs_branch != 0) {
				svs_isKs_branch->GetEntry(index);
			} else { 
				printf("branch svs_isKs_branch does not exist!\n");
				exit(1);
			}
			svs_isKs_isLoaded = true;
		}
		return svs_isKs_;
	}
	vector<int> &svs_isLambda()
	{
		if (not svs_isLambda_isLoaded) {
			if (svs_isLambda_branch != 0) {
				svs_isLambda_branch->GetEntry(index);
			} else { 
				printf("branch svs_isLambda_branch does not exist!\n");
				exit(1);
			}
			svs_isLambda_isLoaded = true;
		}
		return svs_isLambda_;
	}
	vector<int> &svs_mc3_id()
	{
		if (not svs_mc3_id_isLoaded) {
			if (svs_mc3_id_branch != 0) {
				svs_mc3_id_branch->GetEntry(index);
			} else { 
				printf("branch svs_mc3_id_branch does not exist!\n");
				exit(1);
			}
			svs_mc3_id_isLoaded = true;
		}
		return svs_mc3_id_;
	}
	vector<int> &svs_nTrks()
	{
		if (not svs_nTrks_isLoaded) {
			if (svs_nTrks_branch != 0) {
				svs_nTrks_branch->GetEntry(index);
			} else { 
				printf("branch svs_nTrks_branch does not exist!\n");
				exit(1);
			}
			svs_nTrks_isLoaded = true;
		}
		return svs_nTrks_;
	}
	vector<int> &mus_tcmet_flag()
	{
		if (not mus_tcmet_flag_isLoaded) {
			if (mus_tcmet_flag_branch != 0) {
				mus_tcmet_flag_branch->GetEntry(index);
			} else { 
				printf("branch mus_tcmet_flag_branch does not exist!\n");
				exit(1);
			}
			mus_tcmet_flag_isLoaded = true;
		}
		return mus_tcmet_flag_;
	}
	vector<int> &trks_algo()
	{
		if (not trks_algo_isLoaded) {
			if (trks_algo_branch != 0) {
				trks_algo_branch->GetEntry(index);
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
	vector<int> &trks_nlayers()
	{
		if (not trks_nlayers_isLoaded) {
			if (trks_nlayers_branch != 0) {
				trks_nlayers_branch->GetEntry(index);
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
	vector<int> &trkjets_ntrks()
	{
		if (not trkjets_ntrks_isLoaded) {
			if (trkjets_ntrks_branch != 0) {
				trkjets_ntrks_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_ntrks_branch does not exist!\n");
				exit(1);
			}
			trkjets_ntrks_isLoaded = true;
		}
		return trkjets_ntrks_;
	}
	vector<int> &trkjets_vtxs_idx()
	{
		if (not trkjets_vtxs_idx_isLoaded) {
			if (trkjets_vtxs_idx_branch != 0) {
				trkjets_vtxs_idx_branch->GetEntry(index);
			} else { 
				printf("branch trkjets_vtxs_idx_branch does not exist!\n");
				exit(1);
			}
			trkjets_vtxs_idx_isLoaded = true;
		}
		return trkjets_vtxs_idx_;
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
	vector<vector<int> > &genps_lepdaughter_id()
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
	vector<vector<int> > &genps_lepdaughter_idx()
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
	vector<vector<int> > &hlt_trigObjs_id()
	{
		if (not hlt_trigObjs_id_isLoaded) {
			if (hlt_trigObjs_id_branch != 0) {
				hlt_trigObjs_id_branch->GetEntry(index);
			} else { 
				printf("branch hlt_trigObjs_id_branch does not exist!\n");
				exit(1);
			}
			hlt_trigObjs_id_isLoaded = true;
		}
		return hlt_trigObjs_id_;
	}
	vector<vector<int> > &hyp_jets_idx()
	{
		if (not hyp_jets_idx_isLoaded) {
			if (hyp_jets_idx_branch != 0) {
				hyp_jets_idx_branch->GetEntry(index);
			} else { 
				printf("branch hyp_jets_idx_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_idx_isLoaded = true;
		}
		return hyp_jets_idx_;
	}
	vector<vector<int> > &hyp_other_jets_idx()
	{
		if (not hyp_other_jets_idx_isLoaded) {
			if (hyp_other_jets_idx_branch != 0) {
				hyp_other_jets_idx_branch->GetEntry(index);
			} else { 
				printf("branch hyp_other_jets_idx_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_idx_isLoaded = true;
		}
		return hyp_other_jets_idx_;
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
	unsigned int &evt_detectorStatus()
	{
		if (not evt_detectorStatus_isLoaded) {
			if (evt_detectorStatus_branch != 0) {
				evt_detectorStatus_branch->GetEntry(index);
			} else { 
				printf("branch evt_detectorStatus_branch does not exist!\n");
				exit(1);
			}
			evt_detectorStatus_isLoaded = true;
		}
		return evt_detectorStatus_;
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
	unsigned int &genps_flavorHistoryFilterResult()
	{
		if (not genps_flavorHistoryFilterResult_isLoaded) {
			if (genps_flavorHistoryFilterResult_branch != 0) {
				genps_flavorHistoryFilterResult_branch->GetEntry(index);
			} else { 
				printf("branch genps_flavorHistoryFilterResult_branch does not exist!\n");
				exit(1);
			}
			genps_flavorHistoryFilterResult_isLoaded = true;
		}
		return genps_flavorHistoryFilterResult_;
	}
	unsigned int &evt_ngenjets()
	{
		if (not evt_ngenjets_isLoaded) {
			if (evt_ngenjets_branch != 0) {
				evt_ngenjets_branch->GetEntry(index);
			} else { 
				printf("branch evt_ngenjets_branch does not exist!\n");
				exit(1);
			}
			evt_ngenjets_isLoaded = true;
		}
		return evt_ngenjets_;
	}
	unsigned int &genps_signalProcessID()
	{
		if (not genps_signalProcessID_isLoaded) {
			if (genps_signalProcessID_branch != 0) {
				genps_signalProcessID_branch->GetEntry(index);
			} else { 
				printf("branch genps_signalProcessID_branch does not exist!\n");
				exit(1);
			}
			genps_signalProcessID_isLoaded = true;
		}
		return genps_signalProcessID_;
	}
	unsigned int &hlt_bits1()
	{
		if (not hlt_bits1_isLoaded) {
			if (hlt_bits1_branch != 0) {
				hlt_bits1_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits1_branch does not exist!\n");
				exit(1);
			}
			hlt_bits1_isLoaded = true;
		}
		return hlt_bits1_;
	}
	unsigned int &hlt_bits2()
	{
		if (not hlt_bits2_isLoaded) {
			if (hlt_bits2_branch != 0) {
				hlt_bits2_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits2_branch does not exist!\n");
				exit(1);
			}
			hlt_bits2_isLoaded = true;
		}
		return hlt_bits2_;
	}
	unsigned int &hlt_bits3()
	{
		if (not hlt_bits3_isLoaded) {
			if (hlt_bits3_branch != 0) {
				hlt_bits3_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits3_branch does not exist!\n");
				exit(1);
			}
			hlt_bits3_isLoaded = true;
		}
		return hlt_bits3_;
	}
	unsigned int &hlt_bits4()
	{
		if (not hlt_bits4_isLoaded) {
			if (hlt_bits4_branch != 0) {
				hlt_bits4_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits4_branch does not exist!\n");
				exit(1);
			}
			hlt_bits4_isLoaded = true;
		}
		return hlt_bits4_;
	}
	unsigned int &hlt_bits5()
	{
		if (not hlt_bits5_isLoaded) {
			if (hlt_bits5_branch != 0) {
				hlt_bits5_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits5_branch does not exist!\n");
				exit(1);
			}
			hlt_bits5_isLoaded = true;
		}
		return hlt_bits5_;
	}
	unsigned int &hlt_bits6()
	{
		if (not hlt_bits6_isLoaded) {
			if (hlt_bits6_branch != 0) {
				hlt_bits6_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits6_branch does not exist!\n");
				exit(1);
			}
			hlt_bits6_isLoaded = true;
		}
		return hlt_bits6_;
	}
	unsigned int &hlt_bits7()
	{
		if (not hlt_bits7_isLoaded) {
			if (hlt_bits7_branch != 0) {
				hlt_bits7_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits7_branch does not exist!\n");
				exit(1);
			}
			hlt_bits7_isLoaded = true;
		}
		return hlt_bits7_;
	}
	unsigned int &hlt_bits8()
	{
		if (not hlt_bits8_isLoaded) {
			if (hlt_bits8_branch != 0) {
				hlt_bits8_branch->GetEntry(index);
			} else { 
				printf("branch hlt_bits8_branch does not exist!\n");
				exit(1);
			}
			hlt_bits8_isLoaded = true;
		}
		return hlt_bits8_;
	}
	unsigned int &evt_njets()
	{
		if (not evt_njets_isLoaded) {
			if (evt_njets_branch != 0) {
				evt_njets_branch->GetEntry(index);
			} else { 
				printf("branch evt_njets_branch does not exist!\n");
				exit(1);
			}
			evt_njets_isLoaded = true;
		}
		return evt_njets_;
	}
	unsigned int &evt_njpts()
	{
		if (not evt_njpts_isLoaded) {
			if (evt_njpts_branch != 0) {
				evt_njpts_branch->GetEntry(index);
			} else { 
				printf("branch evt_njpts_branch does not exist!\n");
				exit(1);
			}
			evt_njpts_isLoaded = true;
		}
		return evt_njpts_;
	}
	unsigned int &l1_bits1()
	{
		if (not l1_bits1_isLoaded) {
			if (l1_bits1_branch != 0) {
				l1_bits1_branch->GetEntry(index);
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
			} else { 
				printf("branch l1_techbits2_branch does not exist!\n");
				exit(1);
			}
			l1_techbits2_isLoaded = true;
		}
		return l1_techbits2_;
	}
	unsigned int &evt_nphotons()
	{
		if (not evt_nphotons_isLoaded) {
			if (evt_nphotons_branch != 0) {
				evt_nphotons_branch->GetEntry(index);
			} else { 
				printf("branch evt_nphotons_branch does not exist!\n");
				exit(1);
			}
			evt_nphotons_isLoaded = true;
		}
		return evt_nphotons_;
	}
	unsigned int &evt_ecalRecoStatus()
	{
		if (not evt_ecalRecoStatus_isLoaded) {
			if (evt_ecalRecoStatus_branch != 0) {
				evt_ecalRecoStatus_branch->GetEntry(index);
			} else { 
				printf("branch evt_ecalRecoStatus_branch does not exist!\n");
				exit(1);
			}
			evt_ecalRecoStatus_isLoaded = true;
		}
		return evt_ecalRecoStatus_;
	}
	unsigned int &evt_nscs()
	{
		if (not evt_nscs_isLoaded) {
			if (evt_nscs_branch != 0) {
				evt_nscs_branch->GetEntry(index);
			} else { 
				printf("branch evt_nscs_branch does not exist!\n");
				exit(1);
			}
			evt_nscs_isLoaded = true;
		}
		return evt_nscs_;
	}
	unsigned int &evt_ntrkjets()
	{
		if (not evt_ntrkjets_isLoaded) {
			if (evt_ntrkjets_branch != 0) {
				evt_ntrkjets_branch->GetEntry(index);
			} else { 
				printf("branch evt_ntrkjets_branch does not exist!\n");
				exit(1);
			}
			evt_ntrkjets_isLoaded = true;
		}
		return evt_ntrkjets_;
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
	vector<unsigned int> &hlt_prescales()
	{
		if (not hlt_prescales_isLoaded) {
			if (hlt_prescales_branch != 0) {
				hlt_prescales_branch->GetEntry(index);
			} else { 
				printf("branch hlt_prescales_branch does not exist!\n");
				exit(1);
			}
			hlt_prescales_isLoaded = true;
		}
		return hlt_prescales_;
	}
	vector<unsigned int> &hyp_quadlep_bucket()
	{
		if (not hyp_quadlep_bucket_isLoaded) {
			if (hyp_quadlep_bucket_branch != 0) {
				hyp_quadlep_bucket_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_bucket_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_bucket_isLoaded = true;
		}
		return hyp_quadlep_bucket_;
	}
	vector<unsigned int> &hyp_quadlep_first_index()
	{
		if (not hyp_quadlep_first_index_isLoaded) {
			if (hyp_quadlep_first_index_branch != 0) {
				hyp_quadlep_first_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_first_index_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_first_index_isLoaded = true;
		}
		return hyp_quadlep_first_index_;
	}
	vector<unsigned int> &hyp_quadlep_fourth_index()
	{
		if (not hyp_quadlep_fourth_index_isLoaded) {
			if (hyp_quadlep_fourth_index_branch != 0) {
				hyp_quadlep_fourth_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_fourth_index_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_fourth_index_isLoaded = true;
		}
		return hyp_quadlep_fourth_index_;
	}
	vector<unsigned int> &hyp_quadlep_second_index()
	{
		if (not hyp_quadlep_second_index_isLoaded) {
			if (hyp_quadlep_second_index_branch != 0) {
				hyp_quadlep_second_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_second_index_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_second_index_isLoaded = true;
		}
		return hyp_quadlep_second_index_;
	}
	vector<unsigned int> &hyp_quadlep_third_index()
	{
		if (not hyp_quadlep_third_index_isLoaded) {
			if (hyp_quadlep_third_index_branch != 0) {
				hyp_quadlep_third_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_quadlep_third_index_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_third_index_isLoaded = true;
		}
		return hyp_quadlep_third_index_;
	}
	vector<unsigned int> &hyp_trilep_bucket()
	{
		if (not hyp_trilep_bucket_isLoaded) {
			if (hyp_trilep_bucket_branch != 0) {
				hyp_trilep_bucket_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_bucket_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_bucket_isLoaded = true;
		}
		return hyp_trilep_bucket_;
	}
	vector<unsigned int> &hyp_trilep_first_index()
	{
		if (not hyp_trilep_first_index_isLoaded) {
			if (hyp_trilep_first_index_branch != 0) {
				hyp_trilep_first_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_first_index_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_first_index_isLoaded = true;
		}
		return hyp_trilep_first_index_;
	}
	vector<unsigned int> &hyp_trilep_second_index()
	{
		if (not hyp_trilep_second_index_isLoaded) {
			if (hyp_trilep_second_index_branch != 0) {
				hyp_trilep_second_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_second_index_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_second_index_isLoaded = true;
		}
		return hyp_trilep_second_index_;
	}
	vector<unsigned int> &hyp_trilep_third_index()
	{
		if (not hyp_trilep_third_index_isLoaded) {
			if (hyp_trilep_third_index_branch != 0) {
				hyp_trilep_third_index_branch->GetEntry(index);
			} else { 
				printf("branch hyp_trilep_third_index_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_third_index_isLoaded = true;
		}
		return hyp_trilep_third_index_;
	}
	vector<unsigned int> &l1_prescales()
	{
		if (not l1_prescales_isLoaded) {
			if (l1_prescales_branch != 0) {
				l1_prescales_branch->GetEntry(index);
			} else { 
				printf("branch l1_prescales_branch does not exist!\n");
				exit(1);
			}
			l1_prescales_isLoaded = true;
		}
		return l1_prescales_;
	}
	vector<unsigned int> &l1_techtrigprescales()
	{
		if (not l1_techtrigprescales_isLoaded) {
			if (l1_techtrigprescales_branch != 0) {
				l1_techtrigprescales_branch->GetEntry(index);
			} else { 
				printf("branch l1_techtrigprescales_branch does not exist!\n");
				exit(1);
			}
			l1_techtrigprescales_isLoaded = true;
		}
		return l1_techtrigprescales_;
	}
	vector<unsigned int> &els_pat_flag()
	{
		if (not els_pat_flag_isLoaded) {
			if (els_pat_flag_branch != 0) {
				els_pat_flag_branch->GetEntry(index);
			} else { 
				printf("branch els_pat_flag_branch does not exist!\n");
				exit(1);
			}
			els_pat_flag_isLoaded = true;
		}
		return els_pat_flag_;
	}
	vector<unsigned int> &jets_pat_flag()
	{
		if (not jets_pat_flag_isLoaded) {
			if (jets_pat_flag_branch != 0) {
				jets_pat_flag_branch->GetEntry(index);
			} else { 
				printf("branch jets_pat_flag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_flag_isLoaded = true;
		}
		return jets_pat_flag_;
	}
	vector<unsigned int> &mus_pat_flag()
	{
		if (not mus_pat_flag_isLoaded) {
			if (mus_pat_flag_branch != 0) {
				mus_pat_flag_branch->GetEntry(index);
			} else { 
				printf("branch mus_pat_flag_branch does not exist!\n");
				exit(1);
			}
			mus_pat_flag_isLoaded = true;
		}
		return mus_pat_flag_;
	}
	int &evt_nEvts()
	{
		if (not evt_nEvts_isLoaded) {
			if (evt_nEvts_branch != 0) {
				evt_nEvts_branch->GetEntry(index);
			} else { 
				printf("branch evt_nEvts_branch does not exist!\n");
				exit(1);
			}
			evt_nEvts_isLoaded = true;
		}
		return evt_nEvts_;
	}
	float &evt_filt_eff()
	{
		if (not evt_filt_eff_isLoaded) {
			if (evt_filt_eff_branch != 0) {
				evt_filt_eff_branch->GetEntry(index);
			} else { 
				printf("branch evt_filt_eff_branch does not exist!\n");
				exit(1);
			}
			evt_filt_eff_isLoaded = true;
		}
		return evt_filt_eff_;
	}
	bool passHLTTrigger(TString trigName) {
		int trigIndx;
		vector<TString>::const_iterator begin_it = hlt_trigNames().begin();
		vector<TString>::const_iterator end_it = hlt_trigNames().end();
		vector<TString>::const_iterator found_it = find(begin_it, end_it, trigName);
		if(found_it != end_it)
			trigIndx = found_it - begin_it;
		else {
			//cout << "Cannot find Trigger " << trigName << endl; 
			return 0;
		}

		if(trigIndx <= 31) {
			unsigned int bitmask = 1;
			bitmask <<= trigIndx;
			return hlt_bits1() & bitmask;
		}
		if(trigIndx >= 32 && trigIndx <= 63) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 32); 
			return hlt_bits2() & bitmask;
		}
		if(trigIndx >= 64 && trigIndx <= 95) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 64); 
			return hlt_bits3() & bitmask;
		}
		if(trigIndx >= 96 && trigIndx <= 127) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 96); 
			return hlt_bits4() & bitmask;
		}
		if(trigIndx >= 128 && trigIndx <= 159) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 128); 
			return hlt_bits5() & bitmask;
		}
		if(trigIndx >= 160 && trigIndx <= 191) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 160); 
			return hlt_bits6() & bitmask;
		}
		if(trigIndx >= 192 && trigIndx <= 223) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 192); 
			return hlt_bits7() & bitmask;
		}
		if(trigIndx >= 224 && trigIndx <= 255) {
			unsigned int bitmask = 1;
			bitmask <<= (trigIndx - 224); 
			return hlt_bits8() & bitmask;
		}
	return 0;
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
	TString &evt_CMS2tag();
	TString &evt_dataset();
	vector<TString> &hlt_trigNames();
	vector<TString> &l1_techtrigNames();
	vector<TString> &l1_trigNames();
	vector<TString> &evt_errCategory();
	vector<TString> &evt_errModule();
	vector<TString> &evt_errSeverity();
	bool &evt_eventHasHalo();
	bool &evt_hbheFilter();
	vector<bool> &mus_tightMatch();
	vector<bool> &mus_updatedSta();
	vector<bool> &photons_haspixelSeed();
	vector<double> &jets_closestElectron_DR();
	vector<double> &jets_closestMuon_DR();
	float &evt_bs_Xwidth();
	float &evt_bs_XwidthErr();
	float &evt_bs_Ywidth();
	float &evt_bs_YwidthErr();
	float &evt_bs_dxdz();
	float &evt_bs_dxdzErr();
	float &evt_bs_dydz();
	float &evt_bs_dydzErr();
	float &evt_bs_sigmaZ();
	float &evt_bs_sigmaZErr();
	float &evt_bs_xErr();
	float &evt_bs_yErr();
	float &evt_bs_zErr();
	float &evthcal_dmetx();
	float &evthcal_dmety();
	float &evthcal_dsumet();
	float &evthf_dmetx();
	float &evthf_dmety();
	float &evthf_dsumet();
	float &evt_bField();
	float &evt_kfactor();
	float &evt_scale1fb();
	float &evt_xsec_excl();
	float &evt_xsec_incl();
	float &gen_met();
	float &gen_metPhi();
	float &genps_alphaQCD();
	float &genps_pthat();
	float &genps_qScale();
	float &genps_weight();
	float &gen_sumEt();
	float &hcalnoise_eventChargeFraction();
	float &hcalnoise_eventEMEnergy();
	float &hcalnoise_eventEMFraction();
	float &hcalnoise_eventHadEnergy();
	float &hcalnoise_eventTrackEnergy();
	float &hcalnoise_max10GeVHitTime();
	float &hcalnoise_max25GeVHitTime();
	float &hcalnoise_min10GeVHitTime();
	float &hcalnoise_min25GeVHitTime();
	float &hcalnoise_minE10TS();
	float &hcalnoise_minE2Over10TS();
	float &hcalnoise_minE2TS();
	float &hcalnoise_minHPDEMF();
	float &hcalnoise_minRBXEMF();
	float &hcalnoise_rms10GeVHitTime();
	float &hcalnoise_rms25GeVHitTime();
	float &l1_met_etTot();
	float &l1_met_met();
	float &l1_mht_htTot();
	float &l1_mht_mht();
	float &evt_ecalendcapm_met();
	float &evt_ecalendcapm_metPhi();
	float &evt_ecalendcapp_met();
	float &evt_ecalendcapp_metPhi();
	float &evt_ecalmet();
	float &evt_ecalmetPhi();
	float &evt_endcapm_met();
	float &evt_endcapm_metPhi();
	float &evt_endcapp_met();
	float &evt_endcapp_metPhi();
	float &evt_hcalendcapm_met();
	float &evt_hcalendcapm_metPhi();
	float &evt_hcalendcapp_met();
	float &evt_hcalendcapp_metPhi();
	float &evt_hcalmet();
	float &evt_hcalmetPhi();
	float &evt_met();
	float &evt_metHO();
	float &evt_metHOPhi();
	float &evt_metHOSig();
	float &evt_metMuonCorr();
	float &evt_metMuonCorrPhi();
	float &evt_metMuonCorrSig();
	float &evt_metMuonJESCorr();
	float &evt_metMuonJESCorrPhi();
	float &evt_metMuonJESCorrSig();
	float &evt_metNoHF();
	float &evt_metNoHFHO();
	float &evt_metNoHFHOPhi();
	float &evt_metNoHFHOSig();
	float &evt_metNoHFPhi();
	float &evt_metNoHFSig();
	float &evt_metOpt();
	float &evt_metOptHO();
	float &evt_metOptHOPhi();
	float &evt_metOptHOSig();
	float &evt_metOptNoHF();
	float &evt_metOptNoHFHO();
	float &evt_metOptNoHFHOPhi();
	float &evt_metOptNoHFHOSig();
	float &evt_metOptNoHFPhi();
	float &evt_metOptNoHFSig();
	float &evt_metOptPhi();
	float &evt_metOptSig();
	float &evt_metPhi();
	float &evt_metSig();
	float &evt_sumet();
	float &evt_sumetHO();
	float &evt_sumetMuonCorr();
	float &evt_sumetNoHF();
	float &evt_sumetNoHFHO();
	float &evt_sumetOpt();
	float &evt_sumetOptHO();
	float &evt_sumetOptNoHF();
	float &evt_sumetOptNoHFHO();
	float &met_pat_metCor();
	float &met_pat_metPhiCor();
	float &met_pat_metPhiUncor();
	float &met_pat_metPhiUncorJES();
	float &met_pat_metPhiUncorMuon();
	float &met_pat_metUncor();
	float &met_pat_metUncorJES();
	float &met_pat_metUncorMuon();
	float &pdfinfo_scale();
	float &pdfinfo_x1();
	float &pdfinfo_x2();
	float &evt_pfmet();
	float &evt_pfmetPhi();
	float &evt_pfmetSig();
	float &evt_pfsumet();
	float &evt_tcmet();
	float &evt_tcmetPhi();
	float &evt_tcmetSig();
	float &evt_tcsumet();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_met_p4();
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_mht_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_gp_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_gp_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_motherp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trk_mcp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_conv_pos_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_inner_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_outer_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_prod_vtx();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_inner_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_v4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_mc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jpts_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emiso_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emnoiso_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsc_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsf_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetst_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_mus_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_ecalpos_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitdefault_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitfirsthit_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitpicky_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fittev_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_outerPos_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genMotherP4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genP4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genJet_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genPartonMother_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genParton_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_uncorp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genMotherP4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genP4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_posAtEcal_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_posAtEcal_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_pos_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_vtx_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_flight();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_mc3_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &svs_refitp4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_inner_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_position();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_trk_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_vertex_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trkjets_p4();
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &genps_lepdaughter_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_jets_p4();
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_other_jets_p4();
	vector<float> &jpts_combinedSecondaryVertexBJetTag();
	vector<float> &jpts_combinedSecondaryVertexMVABJetTag();
	vector<float> &jpts_jetBProbabilityBJetTag();
	vector<float> &jpts_jetProbabilityBJetTag();
	vector<float> &jpts_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &jpts_simpleSecondaryVertexHighPurBJetTags();
	vector<float> &jpts_softElectronByIP3dBJetTag();
	vector<float> &jpts_softElectronByPtBJetTag();
	vector<float> &jpts_softMuonBJetTag();
	vector<float> &jpts_softMuonByIP3dBJetTag();
	vector<float> &jpts_softMuonByPtBJetTag();
	vector<float> &jpts_trackCountingHighEffBJetTag();
	vector<float> &jpts_trackCountingHighPurBJetTag();
	vector<float> &jets_combinedSecondaryVertexBJetTag();
	vector<float> &jets_combinedSecondaryVertexMVABJetTag();
	vector<float> &jets_jetBProbabilityBJetTag();
	vector<float> &jets_jetProbabilityBJetTag();
	vector<float> &jets_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &jets_simpleSecondaryVertexHighPurBJetTags();
	vector<float> &jets_softElectronByIP3dBJetTag();
	vector<float> &jets_softElectronByPtBJetTag();
	vector<float> &jets_softMuonBJetTag();
	vector<float> &jets_softMuonByIP3dBJetTag();
	vector<float> &jets_softMuonByPtBJetTag();
	vector<float> &jets_trackCountingHighEffBJetTag();
	vector<float> &jets_trackCountingHighPurBJetTag();
	vector<float> &pfjets_combinedSecondaryVertexBJetTag();
	vector<float> &pfjets_combinedSecondaryVertexMVABJetTag();
	vector<float> &pfjets_jetBProbabilityBJetTag();
	vector<float> &pfjets_jetProbabilityBJetTag();
	vector<float> &pfjets_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &pfjets_simpleSecondaryVertexHighPurBJetTags();
	vector<float> &pfjets_softElectronByIP3dBJetTag();
	vector<float> &pfjets_softElectronByPtBJetTag();
	vector<float> &pfjets_softMuonBJetTag();
	vector<float> &pfjets_softMuonByIP3dBJetTag();
	vector<float> &pfjets_softMuonByPtBJetTag();
	vector<float> &pfjets_trackCountingHighEffBJetTag();
	vector<float> &pfjets_trackCountingHighPurBJetTag();
	vector<float> &trkjets_combinedSecondaryVertexBJetTag();
	vector<float> &trkjets_combinedSecondaryVertexMVABJetTag();
	vector<float> &trkjets_jetBProbabilityBJetTag();
	vector<float> &trkjets_jetProbabilityBJetTag();
	vector<float> &trkjets_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &trkjets_simpleSecondaryVertexHighPurBJetTags();
	vector<float> &trkjets_softElectronByIP3dBJetTag();
	vector<float> &trkjets_softElectronByPtBJetTag();
	vector<float> &trkjets_softMuonBJetTag();
	vector<float> &trkjets_softMuonByIP3dBJetTag();
	vector<float> &trkjets_softMuonByPtBJetTag();
	vector<float> &trkjets_trackCountingHighEffBJetTag();
	vector<float> &trkjets_trackCountingHighPurBJetTag();
	vector<float> &evt_bs_covMatrix();
	vector<float> &els_mc3dr();
	vector<float> &els_mcdr();
	vector<float> &jets_mc3dr();
	vector<float> &jets_mcdr();
	vector<float> &jets_mc_emEnergy();
	vector<float> &jets_mc_gpdr();
	vector<float> &jets_mc_hadEnergy();
	vector<float> &jets_mc_invEnergy();
	vector<float> &jets_mc_otherEnergy();
	vector<float> &mus_mc3dr();
	vector<float> &mus_mcdr();
	vector<float> &pfjets_mc3dr();
	vector<float> &pfjets_mcdr();
	vector<float> &pfjets_mc_emEnergy();
	vector<float> &pfjets_mc_gpdr();
	vector<float> &pfjets_mc_hadEnergy();
	vector<float> &pfjets_mc_invEnergy();
	vector<float> &pfjets_mc_otherEnergy();
	vector<float> &photons_mc3dr();
	vector<float> &photons_mcdr();
	vector<float> &trk_mc3dr();
	vector<float> &trk_mcdr();
	vector<float> &trks_conv_dcot();
	vector<float> &trks_conv_dist();
	vector<float> &els_ecalJuraIso();
	vector<float> &els_ecalJuraTowerIso();
	vector<float> &els_hcalConeIso();
	vector<float> &els_tkJuraIso();
	vector<float> &els_jetdr();
	vector<float> &els_musdr();
	vector<float> &els_chi2();
	vector<float> &els_conv_dcot();
	vector<float> &els_conv_dist();
	vector<float> &els_conv_radius();
	vector<float> &els_d0();
	vector<float> &els_d0Err();
	vector<float> &els_d0corr();
	vector<float> &els_dEtaIn();
	vector<float> &els_dEtaOut();
	vector<float> &els_dPhiIn();
	vector<float> &els_dPhiInPhiOut();
	vector<float> &els_dPhiOut();
	vector<float> &els_deltaEtaEleClusterTrackAtCalo();
	vector<float> &els_deltaPhiEleClusterTrackAtCalo();
	vector<float> &els_e1x5();
	vector<float> &els_e2x5Max();
	vector<float> &els_e3x3();
	vector<float> &els_e5x5();
	vector<float> &els_eMax();
	vector<float> &els_eOverPIn();
	vector<float> &els_eOverPOut();
	vector<float> &els_eSC();
	vector<float> &els_eSCPresh();
	vector<float> &els_eSCRaw();
	vector<float> &els_eSeed();
	vector<float> &els_eSeedOverPIn();
	vector<float> &els_eSeedOverPOut();
	vector<float> &els_ecalEnergy();
	vector<float> &els_ecalEnergyError();
	vector<float> &els_ecalIso();
	vector<float> &els_ecalIso04();
	vector<float> &els_egamma_looseId();
	vector<float> &els_egamma_robustHighEnergy();
	vector<float> &els_egamma_robustLooseId();
	vector<float> &els_egamma_robustTightId();
	vector<float> &els_egamma_tightId();
	vector<float> &els_electronMomentumError();
	vector<float> &els_etaErr();
	vector<float> &els_etaSC();
	vector<float> &els_fbrem();
	vector<float> &els_hOverE();
	vector<float> &els_hcalDepth1OverEcal();
	vector<float> &els_hcalDepth1TowerSumEt();
	vector<float> &els_hcalDepth1TowerSumEt04();
	vector<float> &els_hcalDepth2OverEcal();
	vector<float> &els_hcalDepth2TowerSumEt();
	vector<float> &els_hcalDepth2TowerSumEt04();
	vector<float> &els_hcalIso();
	vector<float> &els_hcalIso04();
	vector<float> &els_layer1_charge();
	vector<float> &els_mva();
	vector<float> &els_ndof();
	vector<float> &els_phiErr();
	vector<float> &els_phiSC();
	vector<float> &els_ptErr();
	vector<float> &els_sigmaEtaEta();
	vector<float> &els_sigmaIEtaIEta();
	vector<float> &els_sigmaIEtaIEtaSC();
	vector<float> &els_sigmaIPhiIPhi();
	vector<float> &els_sigmaIPhiIPhiSC();
	vector<float> &els_sigmaPhiPhi();
	vector<float> &els_tkIso();
	vector<float> &els_tkIso04();
	vector<float> &els_trackMomentumError();
	vector<float> &els_trkdr();
	vector<float> &els_trkshFrac();
	vector<float> &els_z0();
	vector<float> &els_z0Err();
	vector<float> &els_z0corr();
	vector<float> &gsftrks_chi2();
	vector<float> &gsftrks_d0();
	vector<float> &gsftrks_d0Err();
	vector<float> &gsftrks_d0corr();
	vector<float> &gsftrks_d0corrPhi();
	vector<float> &gsftrks_d0phiCov();
	vector<float> &gsftrks_etaErr();
	vector<float> &gsftrks_layer1_charge();
	vector<float> &gsftrks_ndof();
	vector<float> &gsftrks_phiErr();
	vector<float> &gsftrks_ptErr();
	vector<float> &gsftrks_z0();
	vector<float> &gsftrks_z0Err();
	vector<float> &gsftrks_z0corr();
	vector<float> &hyp_Ht();
	vector<float> &hyp_dPhi_nJet_metMuonJESCorr();
	vector<float> &hyp_dPhi_nJet_muCorrMet();
	vector<float> &hyp_dPhi_nJet_tcMet();
	vector<float> &hyp_dPhi_nJet_unCorrMet();
	vector<float> &hyp_ll_chi2();
	vector<float> &hyp_ll_d0();
	vector<float> &hyp_ll_d0Err();
	vector<float> &hyp_ll_d0corr();
	vector<float> &hyp_ll_dPhi_metMuonJESCorr();
	vector<float> &hyp_ll_dPhi_muCorrMet();
	vector<float> &hyp_ll_dPhi_tcMet();
	vector<float> &hyp_ll_dPhi_unCorrMet();
	vector<float> &hyp_ll_etaErr();
	vector<float> &hyp_ll_ndof();
	vector<float> &hyp_ll_phiErr();
	vector<float> &hyp_ll_ptErr();
	vector<float> &hyp_ll_z0();
	vector<float> &hyp_ll_z0Err();
	vector<float> &hyp_ll_z0corr();
	vector<float> &hyp_lt_chi2();
	vector<float> &hyp_lt_d0();
	vector<float> &hyp_lt_d0Err();
	vector<float> &hyp_lt_d0corr();
	vector<float> &hyp_lt_dPhi_metMuonJESCorr();
	vector<float> &hyp_lt_dPhi_muCorrMet();
	vector<float> &hyp_lt_dPhi_tcMet();
	vector<float> &hyp_lt_dPhi_unCorrMet();
	vector<float> &hyp_lt_etaErr();
	vector<float> &hyp_lt_ndof();
	vector<float> &hyp_lt_phiErr();
	vector<float> &hyp_lt_ptErr();
	vector<float> &hyp_lt_z0();
	vector<float> &hyp_lt_z0Err();
	vector<float> &hyp_lt_z0corr();
	vector<float> &hyp_mt2_metMuonJESCorr();
	vector<float> &hyp_mt2_muCorrMet();
	vector<float> &hyp_mt2_tcMet();
	vector<float> &hyp_sumJetPt();
	vector<float> &hyp_FVFit_chi2ndf();
	vector<float> &hyp_FVFit_prob();
	vector<float> &hyp_FVFit_v4cxx();
	vector<float> &hyp_FVFit_v4cxy();
	vector<float> &hyp_FVFit_v4cyy();
	vector<float> &hyp_FVFit_v4czx();
	vector<float> &hyp_FVFit_v4czy();
	vector<float> &hyp_FVFit_v4czz();
	vector<float> &hyp_ll_ecaliso();
	vector<float> &hyp_ll_trkiso();
	vector<float> &hyp_lt_ecaliso();
	vector<float> &hyp_lt_trkiso();
	vector<float> &jets_approximatefHPD();
	vector<float> &jets_approximatefRBX();
	vector<float> &jets_cor();
	vector<float> &jets_emFrac();
	vector<float> &jets_fHPD();
	vector<float> &jets_fRBX();
	vector<float> &jets_fSubDetector1();
	vector<float> &jets_fSubDetector2();
	vector<float> &jets_fSubDetector3();
	vector<float> &jets_fSubDetector4();
	vector<float> &jets_hitsInN90();
	vector<float> &jets_n90Hits();
	vector<float> &jets_nECALTowers();
	vector<float> &jets_nHCALTowers();
	vector<float> &jets_restrictedEMF();
	vector<float> &jpts_cor();
	vector<float> &jpts_emFrac();
	vector<float> &mus_met_deltax();
	vector<float> &mus_met_deltay();
	vector<float> &mus_eledr();
	vector<float> &mus_jetdr();
	vector<float> &mus_caloCompatibility();
	vector<float> &mus_chi2();
	vector<float> &mus_d0();
	vector<float> &mus_d0Err();
	vector<float> &mus_d0corr();
	vector<float> &mus_e_em();
	vector<float> &mus_e_emS9();
	vector<float> &mus_e_had();
	vector<float> &mus_e_hadS9();
	vector<float> &mus_e_ho();
	vector<float> &mus_e_hoS9();
	vector<float> &mus_etaErr();
	vector<float> &mus_gfit_chi2();
	vector<float> &mus_gfit_d0();
	vector<float> &mus_gfit_d0Err();
	vector<float> &mus_gfit_d0corr();
	vector<float> &mus_gfit_ndof();
	vector<float> &mus_gfit_qoverp();
	vector<float> &mus_gfit_qoverpError();
	vector<float> &mus_gfit_z0();
	vector<float> &mus_gfit_z0Err();
	vector<float> &mus_gfit_z0corr();
	vector<float> &mus_iso03_emEt();
	vector<float> &mus_iso03_hadEt();
	vector<float> &mus_iso03_hoEt();
	vector<float> &mus_iso03_sumPt();
	vector<float> &mus_iso05_emEt();
	vector<float> &mus_iso05_hadEt();
	vector<float> &mus_iso05_hoEt();
	vector<float> &mus_iso05_sumPt();
	vector<float> &mus_iso_ecalvetoDep();
	vector<float> &mus_iso_hcalvetoDep();
	vector<float> &mus_iso_hovetoDep();
	vector<float> &mus_iso_trckvetoDep();
	vector<float> &mus_ndof();
	vector<float> &mus_phiErr();
	vector<float> &mus_ptErr();
	vector<float> &mus_qoverp();
	vector<float> &mus_qoverpError();
	vector<float> &mus_sta_chi2();
	vector<float> &mus_sta_d0();
	vector<float> &mus_sta_d0Err();
	vector<float> &mus_sta_d0corr();
	vector<float> &mus_sta_ndof();
	vector<float> &mus_sta_qoverp();
	vector<float> &mus_sta_qoverpError();
	vector<float> &mus_sta_z0();
	vector<float> &mus_sta_z0Err();
	vector<float> &mus_sta_z0corr();
	vector<float> &mus_timeAtIpInOut();
	vector<float> &mus_timeAtIpInOutErr();
	vector<float> &mus_timeAtIpOutIn();
	vector<float> &mus_timeAtIpOutInErr();
	vector<float> &mus_vertexphi();
	vector<float> &mus_z0();
	vector<float> &mus_z0Err();
	vector<float> &mus_z0corr();
	vector<float> &els_pat_caloIso();
	vector<float> &els_pat_ecalIso();
	vector<float> &els_pat_hcalIso();
	vector<float> &els_pat_looseId();
	vector<float> &els_pat_robustHighEnergy();
	vector<float> &els_pat_robustLooseId();
	vector<float> &els_pat_robustTightId();
	vector<float> &els_pat_scE1x5();
	vector<float> &els_pat_scE2x5Max();
	vector<float> &els_pat_scE5x5();
	vector<float> &els_pat_sigmaEtaEta();
	vector<float> &els_pat_sigmaIEtaIEta();
	vector<float> &els_pat_tightId();
	vector<float> &els_pat_trackIso();
	vector<float> &jets_pat_combinedSecondaryVertexBJetTag();
	vector<float> &jets_pat_combinedSecondaryVertexMVABJetTag();
	vector<float> &jets_pat_coneIsolationTauJetTag();
	vector<float> &jets_pat_impactParameterMVABJetTag();
	vector<float> &jets_pat_jetBProbabilityBJetTag();
	vector<float> &jets_pat_jetCharge();
	vector<float> &jets_pat_jetProbabilityBJetTag();
	vector<float> &jets_pat_noCorrF();
	vector<float> &jets_pat_simpleSecondaryVertexHighEffBJetTag();
	vector<float> &jets_pat_simpleSecondaryVertexHighPurBJetTag();
	vector<float> &jets_pat_softElectronByIP3dBJetTag();
	vector<float> &jets_pat_softElectronByPtBJetTag();
	vector<float> &jets_pat_softMuonBJetTag();
	vector<float> &jets_pat_softMuonByIP3dBJetTag();
	vector<float> &jets_pat_softMuonByPtBJetTag();
	vector<float> &jets_pat_trackCountingHighEffBJetTag();
	vector<float> &jets_pat_trackCountingHighPurBJetTag();
	vector<float> &mus_pat_caloIso();
	vector<float> &mus_pat_calovetoDep();
	vector<float> &mus_pat_ecalIso();
	vector<float> &mus_pat_ecalvetoDep();
	vector<float> &mus_pat_hcalIso();
	vector<float> &mus_pat_hcalvetoDep();
	vector<float> &mus_pat_trackIso();
	vector<float> &mus_pat_trckvetoDep();
	vector<float> &pfels_deltaP();
	vector<float> &pfels_ecalE();
	vector<float> &pfels_hcalE();
	vector<float> &pfels_isoChargedHadrons();
	vector<float> &pfels_isoNeutralHadrons();
	vector<float> &pfels_isoPhotons();
	vector<float> &pfels_mva_emu();
	vector<float> &pfels_mva_epi();
	vector<float> &pfels_mva_nothing_gamma();
	vector<float> &pfels_mva_nothing_nh();
	vector<float> &pfels_mva_pimu();
	vector<float> &pfels_pS1E();
	vector<float> &pfels_pS2E();
	vector<float> &pfels_rawEcalE();
	vector<float> &pfels_rawHcalE();
	vector<float> &pfjets_chargedEmE();
	vector<float> &pfjets_chargedHadronE();
	vector<float> &pfjets_cor();
	vector<float> &pfjets_neutralEmE();
	vector<float> &pfjets_neutralHadronE();
	vector<float> &pfmus_deltaP();
	vector<float> &pfmus_ecalE();
	vector<float> &pfmus_hcalE();
	vector<float> &pfmus_isoChargedHadrons();
	vector<float> &pfmus_isoNeutralHadrons();
	vector<float> &pfmus_isoPhotons();
	vector<float> &pfmus_mva_emu();
	vector<float> &pfmus_mva_epi();
	vector<float> &pfmus_mva_nothing_gamma();
	vector<float> &pfmus_mva_nothing_nh();
	vector<float> &pfmus_mva_pimu();
	vector<float> &pfmus_pS1E();
	vector<float> &pfmus_pS2E();
	vector<float> &pfmus_rawEcalE();
	vector<float> &pfmus_rawHcalE();
	vector<float> &photons_e1x5();
	vector<float> &photons_e2x5Max();
	vector<float> &photons_e3x3();
	vector<float> &photons_e5x5();
	vector<float> &photons_ecalIso03();
	vector<float> &photons_ecalIso04();
	vector<float> &photons_hOverE();
	vector<float> &photons_hcalIso03();
	vector<float> &photons_hcalIso04();
	vector<float> &photons_ntkIsoHollow03();
	vector<float> &photons_ntkIsoHollow04();
	vector<float> &photons_ntkIsoSolid03();
	vector<float> &photons_ntkIsoSolid04();
	vector<float> &photons_sigmaEtaEta();
	vector<float> &photons_sigmaIEtaIEta();
	vector<float> &photons_swissSeed();
	vector<float> &photons_tkIsoHollow03();
	vector<float> &photons_tkIsoHollow04();
	vector<float> &photons_tkIsoSolid03();
	vector<float> &photons_tkIsoSolid04();
	vector<float> &scs_clustersSize();
	vector<float> &scs_crystalsSize();
	vector<float> &scs_e1x3();
	vector<float> &scs_e1x5();
	vector<float> &scs_e2nd();
	vector<float> &scs_e2x2();
	vector<float> &scs_e2x5Max();
	vector<float> &scs_e3x1();
	vector<float> &scs_e3x2();
	vector<float> &scs_e3x3();
	vector<float> &scs_e4x4();
	vector<float> &scs_e5x5();
	vector<float> &scs_eMax();
	vector<float> &scs_eSeed();
	vector<float> &scs_energy();
	vector<float> &scs_eta();
	vector<float> &scs_hoe();
	vector<float> &scs_phi();
	vector<float> &scs_preshowerEnergy();
	vector<float> &scs_rawEnergy();
	vector<float> &scs_sigmaEtaEta();
	vector<float> &scs_sigmaEtaPhi();
	vector<float> &scs_sigmaIEtaIEta();
	vector<float> &scs_sigmaIEtaIEtaSC();
	vector<float> &scs_sigmaIEtaIPhi();
	vector<float> &scs_sigmaIEtaIPhiSC();
	vector<float> &scs_sigmaIPhiIPhi();
	vector<float> &scs_sigmaIPhiIPhiSC();
	vector<float> &scs_sigmaPhiPhi();
	vector<float> &scs_timeSeed();
	vector<float> &svs_anglePV();
	vector<float> &svs_chi2();
	vector<float> &svs_dist3Dsig();
	vector<float> &svs_dist3Dval();
	vector<float> &svs_distXYsig();
	vector<float> &svs_distXYval();
	vector<float> &svs_ndof();
	vector<float> &svs_prob();
	vector<float> &svs_xError();
	vector<float> &svs_yError();
	vector<float> &svs_zError();
	vector<float> &mus_tcmet_deltax();
	vector<float> &mus_tcmet_deltay();
	vector<float> &trks_chi2();
	vector<float> &trks_d0();
	vector<float> &trks_d0Err();
	vector<float> &trks_d0corr();
	vector<float> &trks_d0corrPhi();
	vector<float> &trks_d0phiCov();
	vector<float> &trks_etaErr();
	vector<float> &trks_layer1_charge();
	vector<float> &trks_ndof();
	vector<float> &trks_phiErr();
	vector<float> &trks_ptErr();
	vector<float> &trks_z0();
	vector<float> &trks_z0Err();
	vector<float> &trks_z0corr();
	vector<float> &trkjets_cor();
	vector<float> &trks_d0Errvtx();
	vector<float> &trks_d0vtx();
	vector<float> &vtxs_chi2();
	vector<float> &vtxs_ndof();
	vector<float> &vtxs_sumpt();
	vector<float> &vtxs_xError();
	vector<float> &vtxs_yError();
	vector<float> &vtxs_zError();
	vector<vector<float> > &vtxs_covMatrix();
	int &evt_cscLooseHaloId();
	int &evt_cscTightHaloId();
	int &evt_ecalLooseHaloId();
	int &evt_ecalTightHaloId();
	int &evt_extremeTightHaloId();
	int &evt_globalLooseHaloId();
	int &evt_globalTightHaloId();
	int &evt_hcalLooseHaloId();
	int &evt_hcalTightHaloId();
	int &evt_looseHaloId();
	int &evt_nHaloLikeTracks();
	int &evt_nHaloTriggerCandidates();
	int &evt_tightHaloId();
	int &evt_bsType();
	int &evt_bunchCrossing();
	int &evt_experimentType();
	int &evt_isRealData();
	int &evt_orbitNumber();
	int &evt_storeNumber();
	int &hcalnoise_maxHPDHits();
	int &hcalnoise_maxRBXHits();
	int &hcalnoise_maxZeros();
	int &hcalnoise_noiseFilterStatus();
	int &hcalnoise_noiseType();
	int &hcalnoise_num10GeVHits();
	int &hcalnoise_num25GeVHits();
	int &hcalnoise_numProblematicRBXs();
	int &hcalnoise_passHighLevelNoiseFilter();
	int &hcalnoise_passLooseNoiseFilter();
	int &hcalnoise_passTightNoiseFilter();
	int &l1_nemiso();
	int &l1_nemnoiso();
	int &l1_njetsc();
	int &l1_njetsf();
	int &l1_njetst();
	int &l1_nmus();
	int &pdfinfo_id1();
	int &pdfinfo_id2();
	vector<int> &evt_ecaliPhiSuspects();
	vector<int> &evt_globaliPhiSuspects();
	vector<int> &evt_hcaliPhiSuspects();
	vector<int> &els_mc3_id();
	vector<int> &els_mc3idx();
	vector<int> &els_mc3_motherid();
	vector<int> &els_mc3_motheridx();
	vector<int> &els_mc_id();
	vector<int> &els_mcidx();
	vector<int> &els_mc_motherid();
	vector<int> &jets_mc3_id();
	vector<int> &jets_mc3idx();
	vector<int> &jets_mc_gpidx();
	vector<int> &jets_mc_id();
	vector<int> &jets_mcidx();
	vector<int> &jets_mc_motherid();
	vector<int> &mus_mc3_id();
	vector<int> &mus_mc3idx();
	vector<int> &mus_mc3_motherid();
	vector<int> &mus_mc3_motheridx();
	vector<int> &mus_mc_id();
	vector<int> &mus_mcidx();
	vector<int> &mus_mc_motherid();
	vector<int> &pfjets_mc3_id();
	vector<int> &pfjets_mc3idx();
	vector<int> &pfjets_mc_gpidx();
	vector<int> &pfjets_mc_id();
	vector<int> &pfjets_mcidx();
	vector<int> &pfjets_mc_motherid();
	vector<int> &photons_mc3_id();
	vector<int> &photons_mc3idx();
	vector<int> &photons_mc3_motherid();
	vector<int> &photons_mc3_motheridx();
	vector<int> &photons_mc_id();
	vector<int> &photons_mcidx();
	vector<int> &photons_mc_motherid();
	vector<int> &trk_mc3_id();
	vector<int> &trk_mc3idx();
	vector<int> &trk_mc3_motherid();
	vector<int> &trk_mc3_motheridx();
	vector<int> &trk_mc_id();
	vector<int> &trk_mcidx();
	vector<int> &trk_mc_motherid();
	vector<int> &trks_conv_tkidx();
	vector<int> &els_exp_innerlayers39X();
	vector<int> &els_closestJet();
	vector<int> &els_closestMuon();
	vector<int> &els_pfelsidx();
	vector<int> &els_category();
	vector<int> &els_charge();
	vector<int> &els_class();
	vector<int> &els_conv_tkidx();
	vector<int> &els_exp_innerlayers();
	vector<int> &els_exp_outerlayers();
	vector<int> &els_fiduciality();
	vector<int> &els_gsftrkidx();
	vector<int> &els_layer1_det();
	vector<int> &els_layer1_layer();
	vector<int> &els_layer1_sizerphi();
	vector<int> &els_layer1_sizerz();
	vector<int> &els_lostHits();
	vector<int> &els_lost_pixelhits();
	vector<int> &els_nSeed();
	vector<int> &els_sccharge();
	vector<int> &els_scindex();
	vector<int> &els_trk_charge();
	vector<int> &els_trkidx();
	vector<int> &els_type();
	vector<int> &els_validHits();
	vector<int> &els_valid_pixelhits();
	vector<int> &genps_id();
	vector<int> &genps_id_mother();
	vector<int> &genps_status();
	vector<int> &gsftrks_charge();
	vector<int> &gsftrks_exp_innerlayers();
	vector<int> &gsftrks_exp_outerlayers();
	vector<int> &gsftrks_layer1_det();
	vector<int> &gsftrks_layer1_layer();
	vector<int> &gsftrks_layer1_sizerphi();
	vector<int> &gsftrks_layer1_sizerz();
	vector<int> &gsftrks_lostHits();
	vector<int> &gsftrks_lost_pixelhits();
	vector<int> &gsftrks_nlayers();
	vector<int> &gsftrks_nlayers3D();
	vector<int> &gsftrks_nlayersLost();
	vector<int> &gsftrks_validHits();
	vector<int> &gsftrks_valid_pixelhits();
	vector<int> &hyp_ll_charge();
	vector<int> &hyp_ll_id();
	vector<int> &hyp_ll_index();
	vector<int> &hyp_ll_lostHits();
	vector<int> &hyp_ll_validHits();
	vector<int> &hyp_lt_charge();
	vector<int> &hyp_lt_id();
	vector<int> &hyp_lt_index();
	vector<int> &hyp_lt_lostHits();
	vector<int> &hyp_lt_validHits();
	vector<int> &hyp_njets();
	vector<int> &hyp_nojets();
	vector<int> &hyp_type();
	vector<int> &hyp_FVFit_ndf();
	vector<int> &hyp_FVFit_status();
	vector<int> &hyp_ll_mc_id();
	vector<int> &hyp_ll_mc_motherid();
	vector<int> &hyp_lt_mc_id();
	vector<int> &hyp_lt_mc_motherid();
	vector<int> &hyp_quadlep_first_type();
	vector<int> &hyp_quadlep_fourth_type();
	vector<int> &hyp_quadlep_second_type();
	vector<int> &hyp_quadlep_third_type();
	vector<int> &hyp_trilep_first_type();
	vector<int> &hyp_trilep_second_type();
	vector<int> &hyp_trilep_third_type();
	vector<int> &jets_closestElectron();
	vector<int> &jets_closestMuon();
	vector<int> &l1_emiso_ieta();
	vector<int> &l1_emiso_iphi();
	vector<int> &l1_emiso_rawId();
	vector<int> &l1_emiso_type();
	vector<int> &l1_emnoiso_ieta();
	vector<int> &l1_emnoiso_iphi();
	vector<int> &l1_emnoiso_rawId();
	vector<int> &l1_emnoiso_type();
	vector<int> &l1_jetsc_ieta();
	vector<int> &l1_jetsc_iphi();
	vector<int> &l1_jetsc_rawId();
	vector<int> &l1_jetsc_type();
	vector<int> &l1_jetsf_ieta();
	vector<int> &l1_jetsf_iphi();
	vector<int> &l1_jetsf_rawId();
	vector<int> &l1_jetsf_type();
	vector<int> &l1_jetst_ieta();
	vector<int> &l1_jetst_iphi();
	vector<int> &l1_jetst_rawId();
	vector<int> &l1_jetst_type();
	vector<int> &l1_mus_flags();
	vector<int> &l1_mus_q();
	vector<int> &l1_mus_qual();
	vector<int> &l1_mus_qualFlags();
	vector<int> &mus_met_flag();
	vector<int> &mus_closestEle();
	vector<int> &mus_closestJet();
	vector<int> &mus_pfmusidx();
	vector<int> &mus_charge();
	vector<int> &mus_chi2LocalMomentum();
	vector<int> &mus_chi2LocalPosition();
	vector<int> &mus_gfit_validHits();
	vector<int> &mus_gfit_validSTAHits();
	vector<int> &mus_gfit_validSiHits();
	vector<int> &mus_glbKink();
	vector<int> &mus_glbTrackProbability();
	vector<int> &mus_globalDeltaEtaPhi();
	vector<int> &mus_goodmask();
	vector<int> &mus_iso03_ntrk();
	vector<int> &mus_iso05_ntrk();
	vector<int> &mus_localDistance();
	vector<int> &mus_lostHits();
	vector<int> &mus_nOverlaps();
	vector<int> &mus_nmatches();
	vector<int> &mus_overlap0();
	vector<int> &mus_overlap1();
	vector<int> &mus_pid_TM2DCompatibilityLoose();
	vector<int> &mus_pid_TM2DCompatibilityTight();
	vector<int> &mus_pid_TMLastStationLoose();
	vector<int> &mus_pid_TMLastStationTight();
	vector<int> &mus_staRelChi2();
	vector<int> &mus_sta_validHits();
	vector<int> &mus_timeDirection();
	vector<int> &mus_timeNumStationsUsed();
	vector<int> &mus_trkKink();
	vector<int> &mus_trkRelChi2();
	vector<int> &mus_trk_charge();
	vector<int> &mus_trkidx();
	vector<int> &mus_type();
	vector<int> &mus_validHits();
	vector<int> &els_pat_genID();
	vector<int> &els_pat_genMotherID();
	vector<int> &jets_pat_genPartonMother_id();
	vector<int> &jets_pat_genParton_id();
	vector<int> &jets_pat_jetIDLoose();
	vector<int> &jets_pat_jetIDLooseAOD();
	vector<int> &jets_pat_jetIDMinimal();
	vector<int> &jets_pat_jetIDTight();
	vector<int> &jets_pat_partonFlavour();
	vector<int> &mus_pat_genID();
	vector<int> &mus_pat_genMotherID();
	vector<int> &pfels_elsidx();
	vector<int> &pfels_charge();
	vector<int> &pfels_flag();
	vector<int> &pfels_particleId();
	vector<int> &pfjets_chargedMultiplicity();
	vector<int> &pfjets_muonMultiplicity();
	vector<int> &pfjets_neutralMultiplicity();
	vector<int> &pfmus_musidx();
	vector<int> &pfmus_charge();
	vector<int> &pfmus_flag();
	vector<int> &pfmus_particleId();
	vector<int> &photons_fiduciality();
	vector<int> &photons_scindex();
	vector<int> &scs_detIdSeed();
	vector<int> &scs_elsidx();
	vector<int> &scs_severitySeed();
	vector<int> &svs_isKs();
	vector<int> &svs_isLambda();
	vector<int> &svs_mc3_id();
	vector<int> &svs_nTrks();
	vector<int> &mus_tcmet_flag();
	vector<int> &trks_algo();
	vector<int> &trks_charge();
	vector<int> &trks_exp_innerlayers();
	vector<int> &trks_exp_outerlayers();
	vector<int> &trks_layer1_det();
	vector<int> &trks_layer1_layer();
	vector<int> &trks_layer1_sizerphi();
	vector<int> &trks_layer1_sizerz();
	vector<int> &trks_lostHits();
	vector<int> &trks_lost_pixelhits();
	vector<int> &trks_nlayers();
	vector<int> &trks_nlayers3D();
	vector<int> &trks_nlayersLost();
	vector<int> &trks_qualityMask();
	vector<int> &trks_validHits();
	vector<int> &trks_valid_pixelhits();
	vector<int> &trks_elsidx();
	vector<int> &trk_musidx();
	vector<int> &trkjets_ntrks();
	vector<int> &trkjets_vtxs_idx();
	vector<int> &vtxs_isFake();
	vector<int> &vtxs_isValid();
	vector<int> &vtxs_tracksSize();
	vector<vector<int> > &genps_lepdaughter_id();
	vector<vector<int> > &genps_lepdaughter_idx();
	vector<vector<int> > &hlt_trigObjs_id();
	vector<vector<int> > &hyp_jets_idx();
	vector<vector<int> > &hyp_other_jets_idx();
	unsigned int &evt_nels();
	unsigned int &evt_detectorStatus();
	unsigned int &evt_event();
	unsigned int &evt_lumiBlock();
	unsigned int &evt_run();
	unsigned int &genps_flavorHistoryFilterResult();
	unsigned int &evt_ngenjets();
	unsigned int &genps_signalProcessID();
	unsigned int &hlt_bits1();
	unsigned int &hlt_bits2();
	unsigned int &hlt_bits3();
	unsigned int &hlt_bits4();
	unsigned int &hlt_bits5();
	unsigned int &hlt_bits6();
	unsigned int &hlt_bits7();
	unsigned int &hlt_bits8();
	unsigned int &evt_njets();
	unsigned int &evt_njpts();
	unsigned int &l1_bits1();
	unsigned int &l1_bits2();
	unsigned int &l1_bits3();
	unsigned int &l1_bits4();
	unsigned int &l1_techbits1();
	unsigned int &l1_techbits2();
	unsigned int &evt_nphotons();
	unsigned int &evt_ecalRecoStatus();
	unsigned int &evt_nscs();
	unsigned int &evt_ntrkjets();
	unsigned int &evt_nvtxs();
	vector<unsigned int> &hlt_prescales();
	vector<unsigned int> &hyp_quadlep_bucket();
	vector<unsigned int> &hyp_quadlep_first_index();
	vector<unsigned int> &hyp_quadlep_fourth_index();
	vector<unsigned int> &hyp_quadlep_second_index();
	vector<unsigned int> &hyp_quadlep_third_index();
	vector<unsigned int> &hyp_trilep_bucket();
	vector<unsigned int> &hyp_trilep_first_index();
	vector<unsigned int> &hyp_trilep_second_index();
	vector<unsigned int> &hyp_trilep_third_index();
	vector<unsigned int> &l1_prescales();
	vector<unsigned int> &l1_techtrigprescales();
	vector<unsigned int> &els_pat_flag();
	vector<unsigned int> &jets_pat_flag();
	vector<unsigned int> &mus_pat_flag();
	int &evt_nEvts();
	float &evt_filt_eff();
	bool passHLTTrigger(TString trigName);
	bool passL1Trigger(TString trigName);
}
#endif
