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
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  l1met_p4_;
	TBranch *l1met_p4_branch;
	bool l1met_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_calo_leadtrk_p4_;
	TBranch *taus_calo_leadtrk_p4_branch;
	bool taus_calo_leadtrk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_calo_p4_;
	TBranch *taus_calo_p4_branch;
	bool taus_calo_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genjets_p4_;
	TBranch *genjets_p4_branch;
	bool genjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_lepdaughter_p4_;
	TBranch *genps_lepdaughter_p4_branch;
	bool genps_lepdaughter_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > genps_prod_vtx_;
	TBranch *genps_prod_vtx_branch;
	bool genps_prod_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_ll_mc_p4_;
	TBranch *hyp_ll_mc_p4_branch;
	bool hyp_ll_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_ll_p4_;
	TBranch *hyp_ll_p4_branch;
	bool hyp_ll_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_ll_trk_p4_;
	TBranch *hyp_ll_trk_p4_branch;
	bool hyp_ll_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_lt_mc_p4_;
	TBranch *hyp_lt_mc_p4_branch;
	bool hyp_lt_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_lt_p4_;
	TBranch *hyp_lt_p4_branch;
	bool hyp_lt_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_lt_trk_p4_;
	TBranch *hyp_lt_trk_p4_branch;
	bool hyp_lt_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hyp_p4_;
	TBranch *hyp_p4_branch;
	bool hyp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_p4_;
	TBranch *jets_p4_branch;
	bool jets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jpts_p4_;
	TBranch *jpts_p4_branch;
	bool jpts_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1emiso_p4_;
	TBranch *l1emiso_p4_branch;
	bool l1emiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1emnoiso_p4_;
	TBranch *l1emnoiso_p4_branch;
	bool l1emnoiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1jetsc_p4_;
	TBranch *l1jetsc_p4_branch;
	bool l1jetsc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1jetsf_p4_;
	TBranch *l1jetsf_p4_branch;
	bool l1jetsf_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1jetst_p4_;
	TBranch *l1jetst_p4_branch;
	bool l1jetst_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > l1mus_p4_;
	TBranch *l1mus_p4_branch;
	bool l1mus_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_pat_genMotherP4_;
	TBranch *els_pat_genMotherP4_branch;
	bool els_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_pat_genP4_;
	TBranch *els_pat_genP4_branch;
	bool els_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > els_pat_p4_;
	TBranch *els_pat_p4_branch;
	bool els_pat_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_pat_genJet_p4_;
	TBranch *jets_pat_genJet_p4_branch;
	bool jets_pat_genJet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_pat_genPartonMother_p4_;
	TBranch *jets_pat_genPartonMother_p4_branch;
	bool jets_pat_genPartonMother_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_pat_genParton_p4_;
	TBranch *jets_pat_genParton_p4_branch;
	bool jets_pat_genParton_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_pat_jet_p4_;
	TBranch *jets_pat_jet_p4_branch;
	bool jets_pat_jet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > jets_pat_jet_uncorp4_;
	TBranch *jets_pat_jet_uncorp4_branch;
	bool jets_pat_jet_uncorp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_pat_genMotherP4_;
	TBranch *mus_pat_genMotherP4_branch;
	bool mus_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_pat_genP4_;
	TBranch *mus_pat_genP4_branch;
	bool mus_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > mus_pat_p4_;
	TBranch *mus_pat_p4_branch;
	bool mus_pat_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_pf_lead_chargecand_p4_;
	TBranch *taus_pf_lead_chargecand_p4_branch;
	bool taus_pf_lead_chargecand_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_pf_lead_neutrcand_p4_;
	TBranch *taus_pf_lead_neutrcand_p4_branch;
	bool taus_pf_lead_neutrcand_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_pf_leadtrk_p4_;
	TBranch *taus_pf_leadtrk_p4_branch;
	bool taus_pf_leadtrk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > taus_pf_p4_;
	TBranch *taus_pf_p4_branch;
	bool taus_pf_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > photons_p4_;
	TBranch *photons_p4_branch;
	bool photons_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > scs_p4_;
	TBranch *scs_p4_branch;
	bool scs_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > scs_pos_;
	TBranch *scs_pos_branch;
	bool scs_pos_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > scs_vtx_;
	TBranch *scs_vtx_branch;
	bool scs_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trks_trk_p4_;
	TBranch *trks_trk_p4_branch;
	bool trks_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trks_vertex_p4_;
	TBranch *trks_vertex_p4_branch;
	bool trks_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hlt2ele10LWR_p4_;
	TBranch *hlt2ele10LWR_p4_branch;
	bool hlt2ele10LWR_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hlt2mu3_p4_;
	TBranch *hlt2mu3_p4_branch;
	bool hlt2mu3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltLisoele18LWR_p4_;
	TBranch *hltLisoele18LWR_p4_branch;
	bool hltLisoele18LWR_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltisoele18R_p4_;
	TBranch *hltisoele18R_p4_branch;
	bool hltisoele18R_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltjet30_p4_;
	TBranch *hltjet30_p4_branch;
	bool hltjet30_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltl1jet15_p4_;
	TBranch *hltl1jet15_p4_branch;
	bool hltl1jet15_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltl1met20_p4_;
	TBranch *hltl1met20_p4_branch;
	bool hltl1met20_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltl1mu_p4_;
	TBranch *hltl1mu_p4_branch;
	bool hltl1mu_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltl2mu9_p4_;
	TBranch *hltl2mu9_p4_branch;
	bool hltl2mu9_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltmet25_p4_;
	TBranch *hltmet25_p4_branch;
	bool hltmet25_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltmu11_p4_;
	TBranch *hltmu11_p4_branch;
	bool hltmu11_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > hltmu9_p4_;
	TBranch *hltmu9_p4_branch;
	bool hltmu9_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > trkjets_p4_;
	TBranch *trkjets_p4_branch;
	bool trkjets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > hyp_jets_p4_;
	TBranch *hyp_jets_p4_branch;
	bool hyp_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > hyp_other_jets_p4_;
	TBranch *hyp_other_jets_p4_branch;
	bool hyp_other_jets_p4_isLoaded;
	vector<double> els_inner_positionx_;
	TBranch *els_inner_positionx_branch;
	bool els_inner_positionx_isLoaded;
	vector<double> els_inner_positiony_;
	TBranch *els_inner_positiony_branch;
	bool els_inner_positiony_isLoaded;
	vector<double> els_inner_positionz_;
	TBranch *els_inner_positionz_branch;
	bool els_inner_positionz_isLoaded;
	vector<double> jets_closestElectron_DR_;
	TBranch *jets_closestElectron_DR_branch;
	bool jets_closestElectron_DR_isLoaded;
	vector<double> jets_closestMuon_DR_;
	TBranch *jets_closestMuon_DR_branch;
	bool jets_closestMuon_DR_isLoaded;
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
	float evt_bs_width_;
	TBranch *evt_bs_width_branch;
	bool evt_bs_width_isLoaded;
	float evt_bs_widthErr_;
	TBranch *evt_bs_widthErr_branch;
	bool evt_bs_widthErr_isLoaded;
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
	float	evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float evt_weight_;
	TBranch *evt_weight_branch;
	bool evt_weight_isLoaded;
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
	float l1met_etHad_;
	TBranch *l1met_etHad_branch;
	bool l1met_etHad_isLoaded;
	float l1met_etTot_;
	TBranch *l1met_etTot_branch;
	bool l1met_etTot_isLoaded;
	float l1met_met_;
	TBranch *l1met_met_branch;
	bool l1met_met_isLoaded;
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
	float evt_pfsumet_;
	TBranch *evt_pfsumet_branch;
	bool evt_pfsumet_isLoaded;
	float evt_tcmet_;
	TBranch *evt_tcmet_branch;
	bool evt_tcmet_isLoaded;
	float evt_tcmetPhi_;
	TBranch *evt_tcmetPhi_branch;
	bool evt_tcmetPhi_isLoaded;
	float evt_tcsumet_;
	TBranch *evt_tcsumet_branch;
	bool evt_tcsumet_isLoaded;
	vector<float> jets_combinedSecondaryVertexBJetTag_;
	TBranch *jets_combinedSecondaryVertexBJetTag_branch;
	bool jets_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_combinedSecondaryVertexMVABJetTag_;
	TBranch *jets_combinedSecondaryVertexMVABJetTag_branch;
	bool jets_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> jets_impactParameterMVABJetTag_;
	TBranch *jets_impactParameterMVABJetTag_branch;
	bool jets_impactParameterMVABJetTag_isLoaded;
	vector<float> jets_jetBProbabilityBJetTag_;
	TBranch *jets_jetBProbabilityBJetTag_branch;
	bool jets_jetBProbabilityBJetTag_isLoaded;
	vector<float> jets_jetProbabilityBJetTag_;
	TBranch *jets_jetProbabilityBJetTag_branch;
	bool jets_jetProbabilityBJetTag_isLoaded;
	vector<float> jets_simpleSecondaryVertexBJetTag_;
	TBranch *jets_simpleSecondaryVertexBJetTag_branch;
	bool jets_simpleSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_softElectronBJetTag_;
	TBranch *jets_softElectronBJetTag_branch;
	bool jets_softElectronBJetTag_isLoaded;
	vector<float> jets_softMuonBJetTag_;
	TBranch *jets_softMuonBJetTag_branch;
	bool jets_softMuonBJetTag_isLoaded;
	vector<float> jets_softMuonNoIPBJetTag_;
	TBranch *jets_softMuonNoIPBJetTag_branch;
	bool jets_softMuonNoIPBJetTag_isLoaded;
	vector<float> jets_trackCountingHighEffBJetTag_;
	TBranch *jets_trackCountingHighEffBJetTag_branch;
	bool jets_trackCountingHighEffBJetTag_isLoaded;
	vector<float> jets_trackCountingHighPurBJetTag_;
	TBranch *jets_trackCountingHighPurBJetTag_branch;
	bool jets_trackCountingHighPurBJetTag_isLoaded;
	vector<float> taus_calo_isolationECALhitsEtSum_;
	TBranch *taus_calo_isolationECALhitsEtSum_branch;
	bool taus_calo_isolationECALhitsEtSum_isLoaded;
	vector<float> taus_calo_isolationtrksPtSum_;
	TBranch *taus_calo_isolationtrksPtSum_branch;
	bool taus_calo_isolationtrksPtSum_isLoaded;
	vector<float> taus_calo_leadtrk_HCAL3x3hitsEtSum_;
	TBranch *taus_calo_leadtrk_HCAL3x3hitsEtSum_branch;
	bool taus_calo_leadtrk_HCAL3x3hitsEtSum_isLoaded;
	vector<float> taus_calo_leadtrk_HCAL3x3hottesthitDEta_;
	TBranch *taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch;
	bool taus_calo_leadtrk_HCAL3x3hottesthitDEta_isLoaded;
	vector<float> taus_calo_leadtrk_Signed_Sipt_;
	TBranch *taus_calo_leadtrk_Signed_Sipt_branch;
	bool taus_calo_leadtrk_Signed_Sipt_isLoaded;
	vector<float> taus_calo_leadtrk_chi2_;
	TBranch *taus_calo_leadtrk_chi2_branch;
	bool taus_calo_leadtrk_chi2_isLoaded;
	vector<float> taus_calo_leadtrk_d0_;
	TBranch *taus_calo_leadtrk_d0_branch;
	bool taus_calo_leadtrk_d0_isLoaded;
	vector<float> taus_calo_leadtrk_lostHits_;
	TBranch *taus_calo_leadtrk_lostHits_branch;
	bool taus_calo_leadtrk_lostHits_isLoaded;
	vector<float> taus_calo_leadtrk_ndof_;
	TBranch *taus_calo_leadtrk_ndof_branch;
	bool taus_calo_leadtrk_ndof_isLoaded;
	vector<float> taus_calo_leadtrk_validHits_;
	TBranch *taus_calo_leadtrk_validHits_branch;
	bool taus_calo_leadtrk_validHits_isLoaded;
	vector<float> taus_calo_leadtrk_z0_;
	TBranch *taus_calo_leadtrk_z0_branch;
	bool taus_calo_leadtrk_z0_isLoaded;
	vector<float> taus_calo_maximumHCALhitEt_;
	TBranch *taus_calo_maximumHCALhitEt_branch;
	bool taus_calo_maximumHCALhitEt_isLoaded;
	vector<float> taus_calo_signaltrksInvariantMass_;
	TBranch *taus_calo_signaltrksInvariantMass_branch;
	bool taus_calo_signaltrksInvariantMass_isLoaded;
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
	vector<float> els_ecalJuraIso_;
	TBranch *els_ecalJuraIso_branch;
	bool els_ecalJuraIso_isLoaded;
	vector<float> els_ecalJuraTowerIso_;
	TBranch *els_ecalJuraTowerIso_branch;
	bool els_ecalJuraTowerIso_isLoaded;
	vector<float> els_hcalConeIso_;
	TBranch *els_hcalConeIso_branch;
	bool els_hcalConeIso_isLoaded;
	vector<float> els_jetdr_;
	TBranch *els_jetdr_branch;
	bool els_jetdr_isLoaded;
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
	vector<float> els_ecalIso_;
	TBranch *els_ecalIso_branch;
	bool els_ecalIso_isLoaded;
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
	vector<float> els_hcalIso_;
	TBranch *els_hcalIso_branch;
	bool els_hcalIso_isLoaded;
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
	vector<float> els_tkIso_;
	TBranch *els_tkIso_branch;
	bool els_tkIso_isLoaded;
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
	vector<float> hyp_Ht_;
	TBranch *hyp_Ht_branch;
	bool hyp_Ht_isLoaded;
	vector<float> hyp_dPhi_nJet_muCorrMet_;
	TBranch *hyp_dPhi_nJet_muCorrMet_branch;
	bool hyp_dPhi_nJet_muCorrMet_isLoaded;
	vector<float> hyp_dPhi_nJet_tcMet_;
	TBranch *hyp_dPhi_nJet_tcMet_branch;
	bool hyp_dPhi_nJet_tcMet_isLoaded;
	vector<float> hyp_dPhi_nJet_type1Met_;
	TBranch *hyp_dPhi_nJet_type1Met_branch;
	bool hyp_dPhi_nJet_type1Met_isLoaded;
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
	vector<float> hyp_ll_dPhi_muCorrMet_;
	TBranch *hyp_ll_dPhi_muCorrMet_branch;
	bool hyp_ll_dPhi_muCorrMet_isLoaded;
	vector<float> hyp_ll_dPhi_tcMet_;
	TBranch *hyp_ll_dPhi_tcMet_branch;
	bool hyp_ll_dPhi_tcMet_isLoaded;
	vector<float> hyp_ll_dPhi_type1Met_;
	TBranch *hyp_ll_dPhi_type1Met_branch;
	bool hyp_ll_dPhi_type1Met_isLoaded;
	vector<float> hyp_ll_dPhi_unCorrMet_;
	TBranch *hyp_ll_dPhi_unCorrMet_branch;
	bool hyp_ll_dPhi_unCorrMet_isLoaded;
	vector<float> hyp_ll_etaErr_;
	TBranch *hyp_ll_etaErr_branch;
	bool hyp_ll_etaErr_isLoaded;
	vector<float> hyp_ll_iso_;
	TBranch *hyp_ll_iso_branch;
	bool hyp_ll_iso_isLoaded;
	vector<float> hyp_ll_ndof_;
	TBranch *hyp_ll_ndof_branch;
	bool hyp_ll_ndof_isLoaded;
	vector<float> hyp_ll_outerEta_;
	TBranch *hyp_ll_outerEta_branch;
	bool hyp_ll_outerEta_isLoaded;
	vector<float> hyp_ll_outerPhi_;
	TBranch *hyp_ll_outerPhi_branch;
	bool hyp_ll_outerPhi_isLoaded;
	vector<float> hyp_ll_phiErr_;
	TBranch *hyp_ll_phiErr_branch;
	bool hyp_ll_phiErr_isLoaded;
	vector<float> hyp_ll_ptErr_;
	TBranch *hyp_ll_ptErr_branch;
	bool hyp_ll_ptErr_isLoaded;
	vector<float> hyp_ll_tkIso_;
	TBranch *hyp_ll_tkIso_branch;
	bool hyp_ll_tkIso_isLoaded;
	vector<float> hyp_ll_vertexphi_;
	TBranch *hyp_ll_vertexphi_branch;
	bool hyp_ll_vertexphi_isLoaded;
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
	vector<float> hyp_lt_dPhi_muCorrMet_;
	TBranch *hyp_lt_dPhi_muCorrMet_branch;
	bool hyp_lt_dPhi_muCorrMet_isLoaded;
	vector<float> hyp_lt_dPhi_tcMet_;
	TBranch *hyp_lt_dPhi_tcMet_branch;
	bool hyp_lt_dPhi_tcMet_isLoaded;
	vector<float> hyp_lt_dPhi_type1Met_;
	TBranch *hyp_lt_dPhi_type1Met_branch;
	bool hyp_lt_dPhi_type1Met_isLoaded;
	vector<float> hyp_lt_dPhi_unCorrMet_;
	TBranch *hyp_lt_dPhi_unCorrMet_branch;
	bool hyp_lt_dPhi_unCorrMet_isLoaded;
	vector<float> hyp_lt_etaErr_;
	TBranch *hyp_lt_etaErr_branch;
	bool hyp_lt_etaErr_isLoaded;
	vector<float> hyp_lt_iso_;
	TBranch *hyp_lt_iso_branch;
	bool hyp_lt_iso_isLoaded;
	vector<float> hyp_lt_ndof_;
	TBranch *hyp_lt_ndof_branch;
	bool hyp_lt_ndof_isLoaded;
	vector<float> hyp_lt_outerEta_;
	TBranch *hyp_lt_outerEta_branch;
	bool hyp_lt_outerEta_isLoaded;
	vector<float> hyp_lt_outerPhi_;
	TBranch *hyp_lt_outerPhi_branch;
	bool hyp_lt_outerPhi_isLoaded;
	vector<float> hyp_lt_phiErr_;
	TBranch *hyp_lt_phiErr_branch;
	bool hyp_lt_phiErr_isLoaded;
	vector<float> hyp_lt_ptErr_;
	TBranch *hyp_lt_ptErr_branch;
	bool hyp_lt_ptErr_isLoaded;
	vector<float> hyp_lt_tkIso_;
	TBranch *hyp_lt_tkIso_branch;
	bool hyp_lt_tkIso_isLoaded;
	vector<float> hyp_lt_vertexphi_;
	TBranch *hyp_lt_vertexphi_branch;
	bool hyp_lt_vertexphi_isLoaded;
	vector<float> hyp_lt_z0_;
	TBranch *hyp_lt_z0_branch;
	bool hyp_lt_z0_isLoaded;
	vector<float> hyp_lt_z0Err_;
	TBranch *hyp_lt_z0Err_branch;
	bool hyp_lt_z0Err_isLoaded;
	vector<float> hyp_lt_z0corr_;
	TBranch *hyp_lt_z0corr_branch;
	bool hyp_lt_z0corr_isLoaded;
	vector<float> hyp_sumJetPt_;
	TBranch *hyp_sumJetPt_branch;
	bool hyp_sumJetPt_isLoaded;
	vector<float> jets_EMFcor_;
	TBranch *jets_EMFcor_branch;
	bool jets_EMFcor_isLoaded;
	vector<float> jets_chFrac_;
	TBranch *jets_chFrac_branch;
	bool jets_chFrac_isLoaded;
	vector<float> jets_cor_;
	TBranch *jets_cor_branch;
	bool jets_cor_isLoaded;
	vector<float> jets_emFrac_;
	TBranch *jets_emFrac_branch;
	bool jets_emFrac_isLoaded;
	vector<float> jpts_cor_;
	TBranch *jpts_cor_branch;
	bool jpts_cor_isLoaded;
	vector<float> jpts_emFrac_;
	TBranch *jpts_emFrac_branch;
	bool jpts_emFrac_isLoaded;
	vector<float> jpts_jet_cor_;
	TBranch *jpts_jet_cor_branch;
	bool jpts_jet_cor_isLoaded;
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
	vector<float> mus_trkdr_;
	TBranch *mus_trkdr_branch;
	bool mus_trkdr_isLoaded;
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
	vector<float> jets_pat_bCorrF_;
	TBranch *jets_pat_bCorrF_branch;
	bool jets_pat_bCorrF_isLoaded;
	vector<float> jets_pat_cCorrF_;
	TBranch *jets_pat_cCorrF_branch;
	bool jets_pat_cCorrF_isLoaded;
	vector<float> jets_pat_combinedSecondaryVertexBJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexBJetTag_branch;
	bool jets_pat_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_pat_combinedSecondaryVertexMVABJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexMVABJetTag_branch;
	bool jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float> jets_pat_coneIsolationTauJetTag_;
	TBranch *jets_pat_coneIsolationTauJetTag_branch;
	bool jets_pat_coneIsolationTauJetTag_isLoaded;
	vector<float> jets_pat_gluCorrF_;
	TBranch *jets_pat_gluCorrF_branch;
	bool jets_pat_gluCorrF_isLoaded;
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
	vector<float> jets_pat_simpleSecondaryVertexBJetTag_;
	TBranch *jets_pat_simpleSecondaryVertexBJetTag_branch;
	bool jets_pat_simpleSecondaryVertexBJetTag_isLoaded;
	vector<float> jets_pat_softElectronBJetTag_;
	TBranch *jets_pat_softElectronBJetTag_branch;
	bool jets_pat_softElectronBJetTag_isLoaded;
	vector<float> jets_pat_softMuonBJetTag_;
	TBranch *jets_pat_softMuonBJetTag_branch;
	bool jets_pat_softMuonBJetTag_isLoaded;
	vector<float> jets_pat_softMuonNoIPBJetTag_;
	TBranch *jets_pat_softMuonNoIPBJetTag_branch;
	bool jets_pat_softMuonNoIPBJetTag_isLoaded;
	vector<float> jets_pat_trackCountingHighEffBJetTag_;
	TBranch *jets_pat_trackCountingHighEffBJetTag_branch;
	bool jets_pat_trackCountingHighEffBJetTag_isLoaded;
	vector<float> jets_pat_trackCountingHighPurBJetTag_;
	TBranch *jets_pat_trackCountingHighPurBJetTag_branch;
	bool jets_pat_trackCountingHighPurBJetTag_isLoaded;
	vector<float> jets_pat_udsCorrF_;
	TBranch *jets_pat_udsCorrF_branch;
	bool jets_pat_udsCorrF_isLoaded;
	vector<float> mus_pat_caloIso_;
	TBranch *mus_pat_caloIso_branch;
	bool mus_pat_caloIso_isLoaded;
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
	vector<float> mus_pat_vetoDep_;
	TBranch *mus_pat_vetoDep_branch;
	bool mus_pat_vetoDep_isLoaded;
	vector<float> taus_pf_bremsRecoveryEOverPLead_;
	TBranch *taus_pf_bremsRecoveryEOverPLead_branch;
	bool taus_pf_bremsRecoveryEOverPLead_isLoaded;
	vector<float> taus_pf_caloComp_;
	TBranch *taus_pf_caloComp_branch;
	bool taus_pf_caloComp_isLoaded;
	vector<float> taus_pf_ecalStripSumEOverPLead_;
	TBranch *taus_pf_ecalStripSumEOverPLead_branch;
	bool taus_pf_ecalStripSumEOverPLead_isLoaded;
	vector<float> taus_pf_emf_;
	TBranch *taus_pf_emf_branch;
	bool taus_pf_emf_isLoaded;
	vector<float> taus_pf_hcal3x3OverPLead_;
	TBranch *taus_pf_hcal3x3OverPLead_branch;
	bool taus_pf_hcal3x3OverPLead_isLoaded;
	vector<float> taus_pf_hcalMaxOverPLead_;
	TBranch *taus_pf_hcalMaxOverPLead_branch;
	bool taus_pf_hcalMaxOverPLead_isLoaded;
	vector<float> taus_pf_hcalTotOverPLead_;
	TBranch *taus_pf_hcalTotOverPLead_branch;
	bool taus_pf_hcalTotOverPLead_isLoaded;
	vector<float> taus_pf_isolationchargecandPtSum_;
	TBranch *taus_pf_isolationchargecandPtSum_branch;
	bool taus_pf_isolationchargecandPtSum_isLoaded;
	vector<float> taus_pf_isolationgammacandEtSum_;
	TBranch *taus_pf_isolationgammacandEtSum_branch;
	bool taus_pf_isolationgammacandEtSum_isLoaded;
	vector<float> taus_pf_lead_chargecand_Signed_Sipt_;
	TBranch *taus_pf_lead_chargecand_Signed_Sipt_branch;
	bool taus_pf_lead_chargecand_Signed_Sipt_isLoaded;
	vector<float> taus_pf_leadtrk_chi2_;
	TBranch *taus_pf_leadtrk_chi2_branch;
	bool taus_pf_leadtrk_chi2_isLoaded;
	vector<float> taus_pf_leadtrk_d0_;
	TBranch *taus_pf_leadtrk_d0_branch;
	bool taus_pf_leadtrk_d0_isLoaded;
	vector<float> taus_pf_leadtrk_lostHits_;
	TBranch *taus_pf_leadtrk_lostHits_branch;
	bool taus_pf_leadtrk_lostHits_isLoaded;
	vector<float> taus_pf_leadtrk_ndof_;
	TBranch *taus_pf_leadtrk_ndof_branch;
	bool taus_pf_leadtrk_ndof_isLoaded;
	vector<float> taus_pf_leadtrk_validHits_;
	TBranch *taus_pf_leadtrk_validHits_branch;
	bool taus_pf_leadtrk_validHits_isLoaded;
	vector<float> taus_pf_leadtrk_z0_;
	TBranch *taus_pf_leadtrk_z0_branch;
	bool taus_pf_leadtrk_z0_isLoaded;
	vector<float> taus_pf_maximumHCALPFClusterEt_;
	TBranch *taus_pf_maximumHCALPFClusterEt_branch;
	bool taus_pf_maximumHCALPFClusterEt_isLoaded;
	vector<float> taus_pf_segComp_;
	TBranch *taus_pf_segComp_branch;
	bool taus_pf_segComp_isLoaded;
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
	vector<float> photons_eMax_;
	TBranch *photons_eMax_branch;
	bool photons_eMax_isLoaded;
	vector<float> photons_eSC_;
	TBranch *photons_eSC_branch;
	bool photons_eSC_isLoaded;
	vector<float> photons_eSCPresh_;
	TBranch *photons_eSCPresh_branch;
	bool photons_eSCPresh_isLoaded;
	vector<float> photons_eSCRaw_;
	TBranch *photons_eSCRaw_branch;
	bool photons_eSCRaw_isLoaded;
	vector<float> photons_eSeed_;
	TBranch *photons_eSeed_branch;
	bool photons_eSeed_isLoaded;
	vector<float> photons_ecalIso_;
	TBranch *photons_ecalIso_branch;
	bool photons_ecalIso_isLoaded;
	vector<float> photons_hOverE_;
	TBranch *photons_hOverE_branch;
	bool photons_hOverE_isLoaded;
	vector<float> photons_hcalIso_;
	TBranch *photons_hcalIso_branch;
	bool photons_hcalIso_isLoaded;
	vector<float> photons_sigmaEtaEta_;
	TBranch *photons_sigmaEtaEta_branch;
	bool photons_sigmaEtaEta_isLoaded;
	vector<float> photons_sigmaIEtaIEta_;
	TBranch *photons_sigmaIEtaIEta_branch;
	bool photons_sigmaIEtaIEta_isLoaded;
	vector<float> photons_sigmaIPhiIPhi_;
	TBranch *photons_sigmaIPhiIPhi_branch;
	bool photons_sigmaIPhiIPhi_isLoaded;
	vector<float> photons_sigmaPhiPhi_;
	TBranch *photons_sigmaPhiPhi_branch;
	bool photons_sigmaPhiPhi_isLoaded;
	vector<float> photons_tkIso_;
	TBranch *photons_tkIso_branch;
	bool photons_tkIso_isLoaded;
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
	vector<float> trks_etaErr_;
	TBranch *trks_etaErr_branch;
	bool trks_etaErr_isLoaded;
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
	int evt_nl1emiso_;
	TBranch *evt_nl1emiso_branch;
	bool evt_nl1emiso_isLoaded;
	int evt_nl1emnoiso_;
	TBranch *evt_nl1emnoiso_branch;
	bool evt_nl1emnoiso_isLoaded;
	int evt_nl1jetsc_;
	TBranch *evt_nl1jetsc_branch;
	bool evt_nl1jetsc_isLoaded;
	int evt_nl1jetsf_;
	TBranch *evt_nl1jetsf_branch;
	bool evt_nl1jetsf_isLoaded;
	int evt_nl1jetst_;
	TBranch *evt_nl1jetst_branch;
	bool evt_nl1jetst_isLoaded;
	int evt_nl1mus_;
	TBranch *evt_nl1mus_branch;
	bool evt_nl1mus_isLoaded;
	int pdfinfo_id1_;
	TBranch *pdfinfo_id1_branch;
	bool pdfinfo_id1_isLoaded;
	int pdfinfo_id2_;
	TBranch *pdfinfo_id2_branch;
	bool pdfinfo_id2_isLoaded;
	vector<int> taus_calo_charge_;
	TBranch *taus_calo_charge_branch;
	bool taus_calo_charge_isLoaded;
	vector<int> taus_calo_iso_ntrks_;
	TBranch *taus_calo_iso_ntrks_branch;
	bool taus_calo_iso_ntrks_isLoaded;
	vector<int> taus_calo_sig_ntrks_;
	TBranch *taus_calo_sig_ntrks_branch;
	bool taus_calo_sig_ntrks_isLoaded;
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
	vector<int> els_closestJet_;
	TBranch *els_closestJet_branch;
	bool els_closestJet_isLoaded;
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
	vector<int> hyp_ll_mc_id_;
	TBranch *hyp_ll_mc_id_branch;
	bool hyp_ll_mc_id_isLoaded;
	vector<int> hyp_ll_mc_motherid_;
	TBranch *hyp_ll_mc_motherid_branch;
	bool hyp_ll_mc_motherid_isLoaded;
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
	vector<int> hyp_lt_mc_id_;
	TBranch *hyp_lt_mc_id_branch;
	bool hyp_lt_mc_id_isLoaded;
	vector<int> hyp_lt_mc_motherid_;
	TBranch *hyp_lt_mc_motherid_branch;
	bool hyp_lt_mc_motherid_isLoaded;
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
	vector<int> l1emiso_ieta_;
	TBranch *l1emiso_ieta_branch;
	bool l1emiso_ieta_isLoaded;
	vector<int> l1emiso_iphi_;
	TBranch *l1emiso_iphi_branch;
	bool l1emiso_iphi_isLoaded;
	vector<int> l1emiso_rawId_;
	TBranch *l1emiso_rawId_branch;
	bool l1emiso_rawId_isLoaded;
	vector<int> l1emiso_type_;
	TBranch *l1emiso_type_branch;
	bool l1emiso_type_isLoaded;
	vector<int> l1emnoiso_ieta_;
	TBranch *l1emnoiso_ieta_branch;
	bool l1emnoiso_ieta_isLoaded;
	vector<int> l1emnoiso_iphi_;
	TBranch *l1emnoiso_iphi_branch;
	bool l1emnoiso_iphi_isLoaded;
	vector<int> l1emnoiso_rawId_;
	TBranch *l1emnoiso_rawId_branch;
	bool l1emnoiso_rawId_isLoaded;
	vector<int> l1emnoiso_type_;
	TBranch *l1emnoiso_type_branch;
	bool l1emnoiso_type_isLoaded;
	vector<int> l1jetsc_ieta_;
	TBranch *l1jetsc_ieta_branch;
	bool l1jetsc_ieta_isLoaded;
	vector<int> l1jetsc_iphi_;
	TBranch *l1jetsc_iphi_branch;
	bool l1jetsc_iphi_isLoaded;
	vector<int> l1jetsc_rawId_;
	TBranch *l1jetsc_rawId_branch;
	bool l1jetsc_rawId_isLoaded;
	vector<int> l1jetsc_type_;
	TBranch *l1jetsc_type_branch;
	bool l1jetsc_type_isLoaded;
	vector<int> l1jetsf_ieta_;
	TBranch *l1jetsf_ieta_branch;
	bool l1jetsf_ieta_isLoaded;
	vector<int> l1jetsf_iphi_;
	TBranch *l1jetsf_iphi_branch;
	bool l1jetsf_iphi_isLoaded;
	vector<int> l1jetsf_rawId_;
	TBranch *l1jetsf_rawId_branch;
	bool l1jetsf_rawId_isLoaded;
	vector<int> l1jetsf_type_;
	TBranch *l1jetsf_type_branch;
	bool l1jetsf_type_isLoaded;
	vector<int> l1jetst_ieta_;
	TBranch *l1jetst_ieta_branch;
	bool l1jetst_ieta_isLoaded;
	vector<int> l1jetst_iphi_;
	TBranch *l1jetst_iphi_branch;
	bool l1jetst_iphi_isLoaded;
	vector<int> l1jetst_rawId_;
	TBranch *l1jetst_rawId_branch;
	bool l1jetst_rawId_isLoaded;
	vector<int> l1jetst_type_;
	TBranch *l1jetst_type_branch;
	bool l1jetst_type_isLoaded;
	vector<int> l1mus_flags_;
	TBranch *l1mus_flags_branch;
	bool l1mus_flags_isLoaded;
	vector<int> l1mus_q_;
	TBranch *l1mus_q_branch;
	bool l1mus_q_isLoaded;
	vector<int> l1mus_qual_;
	TBranch *l1mus_qual_branch;
	bool l1mus_qual_isLoaded;
	vector<int> l1mus_qualFlags_;
	TBranch *l1mus_qualFlags_branch;
	bool l1mus_qualFlags_isLoaded;
	vector<int> mus_met_flag_;
	TBranch *mus_met_flag_branch;
	bool mus_met_flag_isLoaded;
	vector<int> mus_closestEle_;
	TBranch *mus_closestEle_branch;
	bool mus_closestEle_isLoaded;
	vector<int> mus_closestJet_;
	TBranch *mus_closestJet_branch;
	bool mus_closestJet_isLoaded;
	vector<int> mus_trkidx_;
	TBranch *mus_trkidx_branch;
	bool mus_trkidx_isLoaded;
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
	vector<int> jets_pat_partonFlavour_;
	TBranch *jets_pat_partonFlavour_branch;
	bool jets_pat_partonFlavour_isLoaded;
	vector<int> mus_pat_genID_;
	TBranch *mus_pat_genID_branch;
	bool mus_pat_genID_isLoaded;
	vector<int> mus_pat_genMotherID_;
	TBranch *mus_pat_genMotherID_branch;
	bool mus_pat_genMotherID_isLoaded;
	vector<int> taus_pf_charge_;
	TBranch *taus_pf_charge_branch;
	bool taus_pf_charge_isLoaded;
	vector<int> taus_pf_electronPreID_;
	TBranch *taus_pf_electronPreID_branch;
	bool taus_pf_electronPreID_isLoaded;
	vector<int> taus_pf_hasMuonReference_;
	TBranch *taus_pf_hasMuonReference_branch;
	bool taus_pf_hasMuonReference_isLoaded;
	vector<int> taus_pf_iso_ncharge_cand_;
	TBranch *taus_pf_iso_ncharge_cand_branch;
	bool taus_pf_iso_ncharge_cand_isLoaded;
	vector<int> taus_pf_iso_ngamma_cand_;
	TBranch *taus_pf_iso_ngamma_cand_branch;
	bool taus_pf_iso_ngamma_cand_isLoaded;
	vector<int> taus_pf_iso_nneutr_cand_;
	TBranch *taus_pf_iso_nneutr_cand_branch;
	bool taus_pf_iso_nneutr_cand_isLoaded;
	vector<int> taus_pf_muonPreID_;
	TBranch *taus_pf_muonPreID_branch;
	bool taus_pf_muonPreID_isLoaded;
	vector<int> taus_pf_sig_ncharge_cand_;
	TBranch *taus_pf_sig_ncharge_cand_branch;
	bool taus_pf_sig_ncharge_cand_isLoaded;
	vector<int> taus_pf_sig_ngamma_cand_;
	TBranch *taus_pf_sig_ngamma_cand_branch;
	bool taus_pf_sig_ngamma_cand_isLoaded;
	vector<int> taus_pf_sig_nneutr_cand_;
	TBranch *taus_pf_sig_nneutr_cand_branch;
	bool taus_pf_sig_nneutr_cand_isLoaded;
	vector<int> mus_tcmet_flag_;
	TBranch *mus_tcmet_flag_branch;
	bool mus_tcmet_flag_isLoaded;
	vector<int> trks_charge_;
	TBranch *trks_charge_branch;
	bool trks_charge_isLoaded;
	vector<int> trks_lostHits_;
	TBranch *trks_lostHits_branch;
	bool trks_lostHits_isLoaded;
	vector<int> trks_validHits_;
	TBranch *trks_validHits_branch;
	bool trks_validHits_isLoaded;
	vector<int> trks_elsidx_;
	TBranch *trks_elsidx_branch;
	bool trks_elsidx_isLoaded;
	vector<int> trk_musidx_;
	TBranch *trk_musidx_branch;
	bool trk_musidx_isLoaded;
	vector<int> hlt2ele10LWR_id_;
	TBranch *hlt2ele10LWR_id_branch;
	bool hlt2ele10LWR_id_isLoaded;
	vector<int> hlt2ele10LWR_tid_;
	TBranch *hlt2ele10LWR_tid_branch;
	bool hlt2ele10LWR_tid_isLoaded;
	vector<int> hlt2mu3_id_;
	TBranch *hlt2mu3_id_branch;
	bool hlt2mu3_id_isLoaded;
	vector<int> hlt2mu3_tid_;
	TBranch *hlt2mu3_tid_branch;
	bool hlt2mu3_tid_isLoaded;
	vector<int> hltLisoele18LWR_id_;
	TBranch *hltLisoele18LWR_id_branch;
	bool hltLisoele18LWR_id_isLoaded;
	vector<int> hltLisoele18LWR_tid_;
	TBranch *hltLisoele18LWR_tid_branch;
	bool hltLisoele18LWR_tid_isLoaded;
	vector<int> hltisoele18R_id_;
	TBranch *hltisoele18R_id_branch;
	bool hltisoele18R_id_isLoaded;
	vector<int> hltisoele18R_tid_;
	TBranch *hltisoele18R_tid_branch;
	bool hltisoele18R_tid_isLoaded;
	vector<int> hltjet30_id_;
	TBranch *hltjet30_id_branch;
	bool hltjet30_id_isLoaded;
	vector<int> hltjet30_tid_;
	TBranch *hltjet30_tid_branch;
	bool hltjet30_tid_isLoaded;
	vector<int> hltl1jet15_id_;
	TBranch *hltl1jet15_id_branch;
	bool hltl1jet15_id_isLoaded;
	vector<int> hltl1jet15_tid_;
	TBranch *hltl1jet15_tid_branch;
	bool hltl1jet15_tid_isLoaded;
	vector<int> hltl1met20_id_;
	TBranch *hltl1met20_id_branch;
	bool hltl1met20_id_isLoaded;
	vector<int> hltl1met20_tid_;
	TBranch *hltl1met20_tid_branch;
	bool hltl1met20_tid_isLoaded;
	vector<int> hltl1mu_id_;
	TBranch *hltl1mu_id_branch;
	bool hltl1mu_id_isLoaded;
	vector<int> hltl1mu_tid_;
	TBranch *hltl1mu_tid_branch;
	bool hltl1mu_tid_isLoaded;
	vector<int> hltl2mu9_id_;
	TBranch *hltl2mu9_id_branch;
	bool hltl2mu9_id_isLoaded;
	vector<int> hltl2mu9_tid_;
	TBranch *hltl2mu9_tid_branch;
	bool hltl2mu9_tid_isLoaded;
	vector<int> hltmet25_id_;
	TBranch *hltmet25_id_branch;
	bool hltmet25_id_isLoaded;
	vector<int> hltmet25_tid_;
	TBranch *hltmet25_tid_branch;
	bool hltmet25_tid_isLoaded;
	vector<int> hltmu11_id_;
	TBranch *hltmu11_id_branch;
	bool hltmu11_id_isLoaded;
	vector<int> hltmu11_tid_;
	TBranch *hltmu11_tid_branch;
	bool hltmu11_tid_isLoaded;
	vector<int> hltmu9_id_;
	TBranch *hltmu9_id_branch;
	bool hltmu9_id_isLoaded;
	vector<int> hltmu9_tid_;
	TBranch *hltmu9_tid_branch;
	bool hltmu9_tid_isLoaded;
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
	vector<vector<int> > hyp_jets_idx_;
	TBranch *hyp_jets_idx_branch;
	bool hyp_jets_idx_isLoaded;
	vector<vector<int> > hyp_other_jets_idx_;
	TBranch *hyp_other_jets_idx_branch;
	bool hyp_other_jets_idx_isLoaded;
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
	unsigned int evt_ngenjets_;
	TBranch *evt_ngenjets_branch;
	bool evt_ngenjets_isLoaded;
	unsigned int evt_njets_;
	TBranch *evt_njets_branch;
	bool evt_njets_isLoaded;
	unsigned int evt_njpts_;
	TBranch *evt_njpts_branch;
	bool evt_njpts_isLoaded;
	unsigned int evt_nphotons_;
	TBranch *evt_nphotons_branch;
	bool evt_nphotons_isLoaded;
	unsigned int evt_nscs_;
	TBranch *evt_nscs_branch;
	bool evt_nscs_isLoaded;
	unsigned int evt_ntrkjets_;
	TBranch *evt_ntrkjets_branch;
	bool evt_ntrkjets_isLoaded;
	unsigned int evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
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
	vector<unsigned int> els_pat_flag_;
	TBranch *els_pat_flag_branch;
	bool els_pat_flag_isLoaded;
	vector<unsigned int> jets_pat_flag_;
	TBranch *jets_pat_flag_branch;
	bool jets_pat_flag_isLoaded;
	vector<unsigned int> mus_pat_flag_;
	TBranch *mus_pat_flag_branch;
	bool mus_pat_flag_isLoaded;
	float	evt_scale1fb_;
	TBranch *evt_scale1fb_branch;
	bool evt_scale1fb_isLoaded;
	int	evt_nEvts_;
	TBranch *evt_nEvts_branch;
	bool evt_nEvts_isLoaded;
	float	evt_filt_eff_;
	TBranch *evt_filt_eff_branch;
	bool evt_filt_eff_isLoaded;
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
	l1met_p4_branch = 0;
	if (tree->GetAlias("l1met_p4") != 0) {
		l1met_p4_branch = tree->GetBranch(tree->GetAlias("l1met_p4"));
		l1met_p4_branch->SetAddress(&l1met_p4_);
	}
	if(l1met_p4_branch == 0 ) {
	cout << "Branch l1met_p4 does not exist." << endl;
	}
	taus_calo_leadtrk_p4_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_p4") != 0) {
		taus_calo_leadtrk_p4_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_p4"));
		taus_calo_leadtrk_p4_branch->SetAddress(&taus_calo_leadtrk_p4_);
	}
	if(taus_calo_leadtrk_p4_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_p4 does not exist." << endl;
	}
	taus_calo_p4_branch = 0;
	if (tree->GetAlias("taus_calo_p4") != 0) {
		taus_calo_p4_branch = tree->GetBranch(tree->GetAlias("taus_calo_p4"));
		taus_calo_p4_branch->SetAddress(&taus_calo_p4_);
	}
	if(taus_calo_p4_branch == 0 ) {
	cout << "Branch taus_calo_p4 does not exist." << endl;
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
	genjets_p4_branch = 0;
	if (tree->GetAlias("genjets_p4") != 0) {
		genjets_p4_branch = tree->GetBranch(tree->GetAlias("genjets_p4"));
		genjets_p4_branch->SetAddress(&genjets_p4_);
	}
	if(genjets_p4_branch == 0 ) {
	cout << "Branch genjets_p4 does not exist." << endl;
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
	hyp_ll_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_p4") != 0) {
		hyp_ll_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_p4"));
		hyp_ll_mc_p4_branch->SetAddress(&hyp_ll_mc_p4_);
	}
	if(hyp_ll_mc_p4_branch == 0 ) {
	cout << "Branch hyp_ll_mc_p4 does not exist." << endl;
	}
	hyp_ll_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_p4") != 0) {
		hyp_ll_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_p4"));
		hyp_ll_p4_branch->SetAddress(&hyp_ll_p4_);
	}
	if(hyp_ll_p4_branch == 0 ) {
	cout << "Branch hyp_ll_p4 does not exist." << endl;
	}
	hyp_ll_trk_p4_branch = 0;
	if (tree->GetAlias("hyp_ll_trk_p4") != 0) {
		hyp_ll_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_trk_p4"));
		hyp_ll_trk_p4_branch->SetAddress(&hyp_ll_trk_p4_);
	}
	if(hyp_ll_trk_p4_branch == 0 ) {
	cout << "Branch hyp_ll_trk_p4 does not exist." << endl;
	}
	hyp_lt_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_p4") != 0) {
		hyp_lt_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_p4"));
		hyp_lt_mc_p4_branch->SetAddress(&hyp_lt_mc_p4_);
	}
	if(hyp_lt_mc_p4_branch == 0 ) {
	cout << "Branch hyp_lt_mc_p4 does not exist." << endl;
	}
	hyp_lt_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_p4") != 0) {
		hyp_lt_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_p4"));
		hyp_lt_p4_branch->SetAddress(&hyp_lt_p4_);
	}
	if(hyp_lt_p4_branch == 0 ) {
	cout << "Branch hyp_lt_p4 does not exist." << endl;
	}
	hyp_lt_trk_p4_branch = 0;
	if (tree->GetAlias("hyp_lt_trk_p4") != 0) {
		hyp_lt_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_trk_p4"));
		hyp_lt_trk_p4_branch->SetAddress(&hyp_lt_trk_p4_);
	}
	if(hyp_lt_trk_p4_branch == 0 ) {
	cout << "Branch hyp_lt_trk_p4 does not exist." << endl;
	}
	hyp_p4_branch = 0;
	if (tree->GetAlias("hyp_p4") != 0) {
		hyp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_p4"));
		hyp_p4_branch->SetAddress(&hyp_p4_);
	}
	if(hyp_p4_branch == 0 ) {
	cout << "Branch hyp_p4 does not exist." << endl;
	}
	jets_p4_branch = 0;
	if (tree->GetAlias("jets_p4") != 0) {
		jets_p4_branch = tree->GetBranch(tree->GetAlias("jets_p4"));
		jets_p4_branch->SetAddress(&jets_p4_);
	}
	if(jets_p4_branch == 0 ) {
	cout << "Branch jets_p4 does not exist." << endl;
	}
	jpts_p4_branch = 0;
	if (tree->GetAlias("jpts_p4") != 0) {
		jpts_p4_branch = tree->GetBranch(tree->GetAlias("jpts_p4"));
		jpts_p4_branch->SetAddress(&jpts_p4_);
	}
	if(jpts_p4_branch == 0 ) {
	cout << "Branch jpts_p4 does not exist." << endl;
	}
	l1emiso_p4_branch = 0;
	if (tree->GetAlias("l1emiso_p4") != 0) {
		l1emiso_p4_branch = tree->GetBranch(tree->GetAlias("l1emiso_p4"));
		l1emiso_p4_branch->SetAddress(&l1emiso_p4_);
	}
	if(l1emiso_p4_branch == 0 ) {
	cout << "Branch l1emiso_p4 does not exist." << endl;
	}
	l1emnoiso_p4_branch = 0;
	if (tree->GetAlias("l1emnoiso_p4") != 0) {
		l1emnoiso_p4_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_p4"));
		l1emnoiso_p4_branch->SetAddress(&l1emnoiso_p4_);
	}
	if(l1emnoiso_p4_branch == 0 ) {
	cout << "Branch l1emnoiso_p4 does not exist." << endl;
	}
	l1jetsc_p4_branch = 0;
	if (tree->GetAlias("l1jetsc_p4") != 0) {
		l1jetsc_p4_branch = tree->GetBranch(tree->GetAlias("l1jetsc_p4"));
		l1jetsc_p4_branch->SetAddress(&l1jetsc_p4_);
	}
	if(l1jetsc_p4_branch == 0 ) {
	cout << "Branch l1jetsc_p4 does not exist." << endl;
	}
	l1jetsf_p4_branch = 0;
	if (tree->GetAlias("l1jetsf_p4") != 0) {
		l1jetsf_p4_branch = tree->GetBranch(tree->GetAlias("l1jetsf_p4"));
		l1jetsf_p4_branch->SetAddress(&l1jetsf_p4_);
	}
	if(l1jetsf_p4_branch == 0 ) {
	cout << "Branch l1jetsf_p4 does not exist." << endl;
	}
	l1jetst_p4_branch = 0;
	if (tree->GetAlias("l1jetst_p4") != 0) {
		l1jetst_p4_branch = tree->GetBranch(tree->GetAlias("l1jetst_p4"));
		l1jetst_p4_branch->SetAddress(&l1jetst_p4_);
	}
	if(l1jetst_p4_branch == 0 ) {
	cout << "Branch l1jetst_p4 does not exist." << endl;
	}
	l1mus_p4_branch = 0;
	if (tree->GetAlias("l1mus_p4") != 0) {
		l1mus_p4_branch = tree->GetBranch(tree->GetAlias("l1mus_p4"));
		l1mus_p4_branch->SetAddress(&l1mus_p4_);
	}
	if(l1mus_p4_branch == 0 ) {
	cout << "Branch l1mus_p4 does not exist." << endl;
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
	els_pat_genMotherP4_branch = 0;
	if (tree->GetAlias("els_pat_genMotherP4") != 0) {
		els_pat_genMotherP4_branch = tree->GetBranch(tree->GetAlias("els_pat_genMotherP4"));
		els_pat_genMotherP4_branch->SetAddress(&els_pat_genMotherP4_);
	}
	if(els_pat_genMotherP4_branch == 0 ) {
	cout << "Branch els_pat_genMotherP4 does not exist." << endl;
	}
	els_pat_genP4_branch = 0;
	if (tree->GetAlias("els_pat_genP4") != 0) {
		els_pat_genP4_branch = tree->GetBranch(tree->GetAlias("els_pat_genP4"));
		els_pat_genP4_branch->SetAddress(&els_pat_genP4_);
	}
	if(els_pat_genP4_branch == 0 ) {
	cout << "Branch els_pat_genP4 does not exist." << endl;
	}
	els_pat_p4_branch = 0;
	if (tree->GetAlias("els_pat_p4") != 0) {
		els_pat_p4_branch = tree->GetBranch(tree->GetAlias("els_pat_p4"));
		els_pat_p4_branch->SetAddress(&els_pat_p4_);
	}
	if(els_pat_p4_branch == 0 ) {
	cout << "Branch els_pat_p4 does not exist." << endl;
	}
	jets_pat_genJet_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genJet_p4") != 0) {
		jets_pat_genJet_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genJet_p4"));
		jets_pat_genJet_p4_branch->SetAddress(&jets_pat_genJet_p4_);
	}
	if(jets_pat_genJet_p4_branch == 0 ) {
	cout << "Branch jets_pat_genJet_p4 does not exist." << endl;
	}
	jets_pat_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genPartonMother_p4") != 0) {
		jets_pat_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genPartonMother_p4"));
		jets_pat_genPartonMother_p4_branch->SetAddress(&jets_pat_genPartonMother_p4_);
	}
	if(jets_pat_genPartonMother_p4_branch == 0 ) {
	cout << "Branch jets_pat_genPartonMother_p4 does not exist." << endl;
	}
	jets_pat_genParton_p4_branch = 0;
	if (tree->GetAlias("jets_pat_genParton_p4") != 0) {
		jets_pat_genParton_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_genParton_p4"));
		jets_pat_genParton_p4_branch->SetAddress(&jets_pat_genParton_p4_);
	}
	if(jets_pat_genParton_p4_branch == 0 ) {
	cout << "Branch jets_pat_genParton_p4 does not exist." << endl;
	}
	jets_pat_jet_p4_branch = 0;
	if (tree->GetAlias("jets_pat_jet_p4") != 0) {
		jets_pat_jet_p4_branch = tree->GetBranch(tree->GetAlias("jets_pat_jet_p4"));
		jets_pat_jet_p4_branch->SetAddress(&jets_pat_jet_p4_);
	}
	if(jets_pat_jet_p4_branch == 0 ) {
	cout << "Branch jets_pat_jet_p4 does not exist." << endl;
	}
	jets_pat_jet_uncorp4_branch = 0;
	if (tree->GetAlias("jets_pat_jet_uncorp4") != 0) {
		jets_pat_jet_uncorp4_branch = tree->GetBranch(tree->GetAlias("jets_pat_jet_uncorp4"));
		jets_pat_jet_uncorp4_branch->SetAddress(&jets_pat_jet_uncorp4_);
	}
	if(jets_pat_jet_uncorp4_branch == 0 ) {
	cout << "Branch jets_pat_jet_uncorp4 does not exist." << endl;
	}
	mus_pat_genMotherP4_branch = 0;
	if (tree->GetAlias("mus_pat_genMotherP4") != 0) {
		mus_pat_genMotherP4_branch = tree->GetBranch(tree->GetAlias("mus_pat_genMotherP4"));
		mus_pat_genMotherP4_branch->SetAddress(&mus_pat_genMotherP4_);
	}
	if(mus_pat_genMotherP4_branch == 0 ) {
	cout << "Branch mus_pat_genMotherP4 does not exist." << endl;
	}
	mus_pat_genP4_branch = 0;
	if (tree->GetAlias("mus_pat_genP4") != 0) {
		mus_pat_genP4_branch = tree->GetBranch(tree->GetAlias("mus_pat_genP4"));
		mus_pat_genP4_branch->SetAddress(&mus_pat_genP4_);
	}
	if(mus_pat_genP4_branch == 0 ) {
	cout << "Branch mus_pat_genP4 does not exist." << endl;
	}
	mus_pat_p4_branch = 0;
	if (tree->GetAlias("mus_pat_p4") != 0) {
		mus_pat_p4_branch = tree->GetBranch(tree->GetAlias("mus_pat_p4"));
		mus_pat_p4_branch->SetAddress(&mus_pat_p4_);
	}
	if(mus_pat_p4_branch == 0 ) {
	cout << "Branch mus_pat_p4 does not exist." << endl;
	}
	taus_pf_lead_chargecand_p4_branch = 0;
	if (tree->GetAlias("taus_pf_lead_chargecand_p4") != 0) {
		taus_pf_lead_chargecand_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_lead_chargecand_p4"));
		taus_pf_lead_chargecand_p4_branch->SetAddress(&taus_pf_lead_chargecand_p4_);
	}
	if(taus_pf_lead_chargecand_p4_branch == 0 ) {
	cout << "Branch taus_pf_lead_chargecand_p4 does not exist." << endl;
	}
	taus_pf_lead_neutrcand_p4_branch = 0;
	if (tree->GetAlias("taus_pf_lead_neutrcand_p4") != 0) {
		taus_pf_lead_neutrcand_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_lead_neutrcand_p4"));
		taus_pf_lead_neutrcand_p4_branch->SetAddress(&taus_pf_lead_neutrcand_p4_);
	}
	if(taus_pf_lead_neutrcand_p4_branch == 0 ) {
	cout << "Branch taus_pf_lead_neutrcand_p4 does not exist." << endl;
	}
	taus_pf_leadtrk_p4_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_p4") != 0) {
		taus_pf_leadtrk_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_p4"));
		taus_pf_leadtrk_p4_branch->SetAddress(&taus_pf_leadtrk_p4_);
	}
	if(taus_pf_leadtrk_p4_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_p4 does not exist." << endl;
	}
	taus_pf_p4_branch = 0;
	if (tree->GetAlias("taus_pf_p4") != 0) {
		taus_pf_p4_branch = tree->GetBranch(tree->GetAlias("taus_pf_p4"));
		taus_pf_p4_branch->SetAddress(&taus_pf_p4_);
	}
	if(taus_pf_p4_branch == 0 ) {
	cout << "Branch taus_pf_p4 does not exist." << endl;
	}
	photons_p4_branch = 0;
	if (tree->GetAlias("photons_p4") != 0) {
		photons_p4_branch = tree->GetBranch(tree->GetAlias("photons_p4"));
		photons_p4_branch->SetAddress(&photons_p4_);
	}
	if(photons_p4_branch == 0 ) {
	cout << "Branch photons_p4 does not exist." << endl;
	}
	scs_p4_branch = 0;
	if (tree->GetAlias("scs_p4") != 0) {
		scs_p4_branch = tree->GetBranch(tree->GetAlias("scs_p4"));
		scs_p4_branch->SetAddress(&scs_p4_);
	}
	if(scs_p4_branch == 0 ) {
	cout << "Branch scs_p4 does not exist." << endl;
	}
	scs_pos_branch = 0;
	if (tree->GetAlias("scs_pos") != 0) {
		scs_pos_branch = tree->GetBranch(tree->GetAlias("scs_pos"));
		scs_pos_branch->SetAddress(&scs_pos_);
	}
	if(scs_pos_branch == 0 ) {
	cout << "Branch scs_pos does not exist." << endl;
	}
	scs_vtx_branch = 0;
	if (tree->GetAlias("scs_vtx") != 0) {
		scs_vtx_branch = tree->GetBranch(tree->GetAlias("scs_vtx"));
		scs_vtx_branch->SetAddress(&scs_vtx_);
	}
	if(scs_vtx_branch == 0 ) {
	cout << "Branch scs_vtx does not exist." << endl;
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
	hlt2ele10LWR_p4_branch = 0;
	if (tree->GetAlias("hlt2ele10LWR_p4") != 0) {
		hlt2ele10LWR_p4_branch = tree->GetBranch(tree->GetAlias("hlt2ele10LWR_p4"));
		hlt2ele10LWR_p4_branch->SetAddress(&hlt2ele10LWR_p4_);
	}
	if(hlt2ele10LWR_p4_branch == 0 ) {
	cout << "Branch hlt2ele10LWR_p4 does not exist." << endl;
	}
	hlt2mu3_p4_branch = 0;
	if (tree->GetAlias("hlt2mu3_p4") != 0) {
		hlt2mu3_p4_branch = tree->GetBranch(tree->GetAlias("hlt2mu3_p4"));
		hlt2mu3_p4_branch->SetAddress(&hlt2mu3_p4_);
	}
	if(hlt2mu3_p4_branch == 0 ) {
	cout << "Branch hlt2mu3_p4 does not exist." << endl;
	}
	hltLisoele18LWR_p4_branch = 0;
	if (tree->GetAlias("hltLisoele18LWR_p4") != 0) {
		hltLisoele18LWR_p4_branch = tree->GetBranch(tree->GetAlias("hltLisoele18LWR_p4"));
		hltLisoele18LWR_p4_branch->SetAddress(&hltLisoele18LWR_p4_);
	}
	if(hltLisoele18LWR_p4_branch == 0 ) {
	cout << "Branch hltLisoele18LWR_p4 does not exist." << endl;
	}
	hltisoele18R_p4_branch = 0;
	if (tree->GetAlias("hltisoele18R_p4") != 0) {
		hltisoele18R_p4_branch = tree->GetBranch(tree->GetAlias("hltisoele18R_p4"));
		hltisoele18R_p4_branch->SetAddress(&hltisoele18R_p4_);
	}
	if(hltisoele18R_p4_branch == 0 ) {
	cout << "Branch hltisoele18R_p4 does not exist." << endl;
	}
	hltjet30_p4_branch = 0;
	if (tree->GetAlias("hltjet30_p4") != 0) {
		hltjet30_p4_branch = tree->GetBranch(tree->GetAlias("hltjet30_p4"));
		hltjet30_p4_branch->SetAddress(&hltjet30_p4_);
	}
	if(hltjet30_p4_branch == 0 ) {
	cout << "Branch hltjet30_p4 does not exist." << endl;
	}
	hltl1jet15_p4_branch = 0;
	if (tree->GetAlias("hltl1jet15_p4") != 0) {
		hltl1jet15_p4_branch = tree->GetBranch(tree->GetAlias("hltl1jet15_p4"));
		hltl1jet15_p4_branch->SetAddress(&hltl1jet15_p4_);
	}
	if(hltl1jet15_p4_branch == 0 ) {
	cout << "Branch hltl1jet15_p4 does not exist." << endl;
	}
	hltl1met20_p4_branch = 0;
	if (tree->GetAlias("hltl1met20_p4") != 0) {
		hltl1met20_p4_branch = tree->GetBranch(tree->GetAlias("hltl1met20_p4"));
		hltl1met20_p4_branch->SetAddress(&hltl1met20_p4_);
	}
	if(hltl1met20_p4_branch == 0 ) {
	cout << "Branch hltl1met20_p4 does not exist." << endl;
	}
	hltl1mu_p4_branch = 0;
	if (tree->GetAlias("hltl1mu_p4") != 0) {
		hltl1mu_p4_branch = tree->GetBranch(tree->GetAlias("hltl1mu_p4"));
		hltl1mu_p4_branch->SetAddress(&hltl1mu_p4_);
	}
	if(hltl1mu_p4_branch == 0 ) {
	cout << "Branch hltl1mu_p4 does not exist." << endl;
	}
	hltl2mu9_p4_branch = 0;
	if (tree->GetAlias("hltl2mu9_p4") != 0) {
		hltl2mu9_p4_branch = tree->GetBranch(tree->GetAlias("hltl2mu9_p4"));
		hltl2mu9_p4_branch->SetAddress(&hltl2mu9_p4_);
	}
	if(hltl2mu9_p4_branch == 0 ) {
	cout << "Branch hltl2mu9_p4 does not exist." << endl;
	}
	hltmet25_p4_branch = 0;
	if (tree->GetAlias("hltmet25_p4") != 0) {
		hltmet25_p4_branch = tree->GetBranch(tree->GetAlias("hltmet25_p4"));
		hltmet25_p4_branch->SetAddress(&hltmet25_p4_);
	}
	if(hltmet25_p4_branch == 0 ) {
	cout << "Branch hltmet25_p4 does not exist." << endl;
	}
	hltmu11_p4_branch = 0;
	if (tree->GetAlias("hltmu11_p4") != 0) {
		hltmu11_p4_branch = tree->GetBranch(tree->GetAlias("hltmu11_p4"));
		hltmu11_p4_branch->SetAddress(&hltmu11_p4_);
	}
	if(hltmu11_p4_branch == 0 ) {
	cout << "Branch hltmu11_p4 does not exist." << endl;
	}
	hltmu9_p4_branch = 0;
	if (tree->GetAlias("hltmu9_p4") != 0) {
		hltmu9_p4_branch = tree->GetBranch(tree->GetAlias("hltmu9_p4"));
		hltmu9_p4_branch->SetAddress(&hltmu9_p4_);
	}
	if(hltmu9_p4_branch == 0 ) {
	cout << "Branch hltmu9_p4 does not exist." << endl;
	}
	trkjets_p4_branch = 0;
	if (tree->GetAlias("trkjets_p4") != 0) {
		trkjets_p4_branch = tree->GetBranch(tree->GetAlias("trkjets_p4"));
		trkjets_p4_branch->SetAddress(&trkjets_p4_);
	}
	if(trkjets_p4_branch == 0 ) {
	cout << "Branch trkjets_p4 does not exist." << endl;
	}
  tree->SetMakeClass(1);
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
	hyp_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_p4") != 0) {
		hyp_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_p4"));
		hyp_jets_p4_branch->SetAddress(&hyp_jets_p4_);
	}
	if(hyp_jets_p4_branch == 0 ) {
	cout << "Branch hyp_jets_p4 does not exist." << endl;
	}
	hyp_other_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_p4") != 0) {
		hyp_other_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_p4"));
		hyp_other_jets_p4_branch->SetAddress(&hyp_other_jets_p4_);
	}
	if(hyp_other_jets_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_p4 does not exist." << endl;
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
	jets_closestElectron_DR_branch = 0;
	if (tree->GetAlias("jets_closestElectron_DR") != 0) {
		jets_closestElectron_DR_branch = tree->GetBranch(tree->GetAlias("jets_closestElectron_DR"));
		jets_closestElectron_DR_branch->SetAddress(&jets_closestElectron_DR_);
	}
	if(jets_closestElectron_DR_branch == 0 ) {
	cout << "Branch jets_closestElectron_DR does not exist." << endl;
	}
	jets_closestMuon_DR_branch = 0;
	if (tree->GetAlias("jets_closestMuon_DR") != 0) {
		jets_closestMuon_DR_branch = tree->GetBranch(tree->GetAlias("jets_closestMuon_DR"));
		jets_closestMuon_DR_branch->SetAddress(&jets_closestMuon_DR_);
	}
	if(jets_closestMuon_DR_branch == 0 ) {
	cout << "Branch jets_closestMuon_DR does not exist." << endl;
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
	evt_bs_width_branch = 0;
	if (tree->GetAlias("evt_bs_width") != 0) {
		evt_bs_width_branch = tree->GetBranch(tree->GetAlias("evt_bs_width"));
		evt_bs_width_branch->SetAddress(&evt_bs_width_);
	}
	if(evt_bs_width_branch == 0 ) {
	cout << "Branch evt_bs_width does not exist." << endl;
	}
	evt_bs_widthErr_branch = 0;
	if (tree->GetAlias("evt_bs_widthErr") != 0) {
		evt_bs_widthErr_branch = tree->GetBranch(tree->GetAlias("evt_bs_widthErr"));
		evt_bs_widthErr_branch->SetAddress(&evt_bs_widthErr_);
	}
	if(evt_bs_widthErr_branch == 0 ) {
	cout << "Branch evt_bs_widthErr does not exist." << endl;
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
	l1met_etHad_branch = 0;
	if (tree->GetAlias("l1met_etHad") != 0) {
		l1met_etHad_branch = tree->GetBranch(tree->GetAlias("l1met_etHad"));
		l1met_etHad_branch->SetAddress(&l1met_etHad_);
	}
	if(l1met_etHad_branch == 0 ) {
	cout << "Branch l1met_etHad does not exist." << endl;
	}
	l1met_etTot_branch = 0;
	if (tree->GetAlias("l1met_etTot") != 0) {
		l1met_etTot_branch = tree->GetBranch(tree->GetAlias("l1met_etTot"));
		l1met_etTot_branch->SetAddress(&l1met_etTot_);
	}
	if(l1met_etTot_branch == 0 ) {
	cout << "Branch l1met_etTot does not exist." << endl;
	}
	l1met_met_branch = 0;
	if (tree->GetAlias("l1met_met") != 0) {
		l1met_met_branch = tree->GetBranch(tree->GetAlias("l1met_met"));
		l1met_met_branch->SetAddress(&l1met_met_);
	}
	if(l1met_met_branch == 0 ) {
	cout << "Branch l1met_met does not exist." << endl;
	}
	evt_met_branch = 0;
	if (tree->GetAlias("evt_met") != 0) {
		evt_met_branch = tree->GetBranch(tree->GetAlias("evt_met"));
		evt_met_branch->SetAddress(&evt_met_);
	}
	if(evt_met_branch == 0 ) {
	cout << "Branch evt_met does not exist." << endl;
	}
	evt_metHO_branch = 0;
	if (tree->GetAlias("evt_metHO") != 0) {
		evt_metHO_branch = tree->GetBranch(tree->GetAlias("evt_metHO"));
		evt_metHO_branch->SetAddress(&evt_metHO_);
	}
	if(evt_metHO_branch == 0 ) {
	cout << "Branch evt_metHO does not exist." << endl;
	}
	evt_metHOPhi_branch = 0;
	if (tree->GetAlias("evt_metHOPhi") != 0) {
		evt_metHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metHOPhi"));
		evt_metHOPhi_branch->SetAddress(&evt_metHOPhi_);
	}
	if(evt_metHOPhi_branch == 0 ) {
	cout << "Branch evt_metHOPhi does not exist." << endl;
	}
	evt_metHOSig_branch = 0;
	if (tree->GetAlias("evt_metHOSig") != 0) {
		evt_metHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metHOSig"));
		evt_metHOSig_branch->SetAddress(&evt_metHOSig_);
	}
	if(evt_metHOSig_branch == 0 ) {
	cout << "Branch evt_metHOSig does not exist." << endl;
	}
	evt_metMuonCorr_branch = 0;
	if (tree->GetAlias("evt_metMuonCorr") != 0) {
		evt_metMuonCorr_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorr"));
		evt_metMuonCorr_branch->SetAddress(&evt_metMuonCorr_);
	}
	if(evt_metMuonCorr_branch == 0 ) {
	cout << "Branch evt_metMuonCorr does not exist." << endl;
	}
	evt_metMuonCorrPhi_branch = 0;
	if (tree->GetAlias("evt_metMuonCorrPhi") != 0) {
		evt_metMuonCorrPhi_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorrPhi"));
		evt_metMuonCorrPhi_branch->SetAddress(&evt_metMuonCorrPhi_);
	}
	if(evt_metMuonCorrPhi_branch == 0 ) {
	cout << "Branch evt_metMuonCorrPhi does not exist." << endl;
	}
	evt_metMuonCorrSig_branch = 0;
	if (tree->GetAlias("evt_metMuonCorrSig") != 0) {
		evt_metMuonCorrSig_branch = tree->GetBranch(tree->GetAlias("evt_metMuonCorrSig"));
		evt_metMuonCorrSig_branch->SetAddress(&evt_metMuonCorrSig_);
	}
	if(evt_metMuonCorrSig_branch == 0 ) {
	cout << "Branch evt_metMuonCorrSig does not exist." << endl;
	}
	evt_metMuonJESCorr_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorr") != 0) {
		evt_metMuonJESCorr_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorr"));
		evt_metMuonJESCorr_branch->SetAddress(&evt_metMuonJESCorr_);
	}
	if(evt_metMuonJESCorr_branch == 0 ) {
	cout << "Branch evt_metMuonJESCorr does not exist." << endl;
	}
	evt_metMuonJESCorrPhi_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorrPhi") != 0) {
		evt_metMuonJESCorrPhi_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorrPhi"));
		evt_metMuonJESCorrPhi_branch->SetAddress(&evt_metMuonJESCorrPhi_);
	}
	if(evt_metMuonJESCorrPhi_branch == 0 ) {
	cout << "Branch evt_metMuonJESCorrPhi does not exist." << endl;
	}
	evt_metMuonJESCorrSig_branch = 0;
	if (tree->GetAlias("evt_metMuonJESCorrSig") != 0) {
		evt_metMuonJESCorrSig_branch = tree->GetBranch(tree->GetAlias("evt_metMuonJESCorrSig"));
		evt_metMuonJESCorrSig_branch->SetAddress(&evt_metMuonJESCorrSig_);
	}
	if(evt_metMuonJESCorrSig_branch == 0 ) {
	cout << "Branch evt_metMuonJESCorrSig does not exist." << endl;
	}
	evt_metNoHF_branch = 0;
	if (tree->GetAlias("evt_metNoHF") != 0) {
		evt_metNoHF_branch = tree->GetBranch(tree->GetAlias("evt_metNoHF"));
		evt_metNoHF_branch->SetAddress(&evt_metNoHF_);
	}
	if(evt_metNoHF_branch == 0 ) {
	cout << "Branch evt_metNoHF does not exist." << endl;
	}
	evt_metNoHFHO_branch = 0;
	if (tree->GetAlias("evt_metNoHFHO") != 0) {
		evt_metNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHO"));
		evt_metNoHFHO_branch->SetAddress(&evt_metNoHFHO_);
	}
	if(evt_metNoHFHO_branch == 0 ) {
	cout << "Branch evt_metNoHFHO does not exist." << endl;
	}
	evt_metNoHFHOPhi_branch = 0;
	if (tree->GetAlias("evt_metNoHFHOPhi") != 0) {
		evt_metNoHFHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHOPhi"));
		evt_metNoHFHOPhi_branch->SetAddress(&evt_metNoHFHOPhi_);
	}
	if(evt_metNoHFHOPhi_branch == 0 ) {
	cout << "Branch evt_metNoHFHOPhi does not exist." << endl;
	}
	evt_metNoHFHOSig_branch = 0;
	if (tree->GetAlias("evt_metNoHFHOSig") != 0) {
		evt_metNoHFHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFHOSig"));
		evt_metNoHFHOSig_branch->SetAddress(&evt_metNoHFHOSig_);
	}
	if(evt_metNoHFHOSig_branch == 0 ) {
	cout << "Branch evt_metNoHFHOSig does not exist." << endl;
	}
	evt_metNoHFPhi_branch = 0;
	if (tree->GetAlias("evt_metNoHFPhi") != 0) {
		evt_metNoHFPhi_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFPhi"));
		evt_metNoHFPhi_branch->SetAddress(&evt_metNoHFPhi_);
	}
	if(evt_metNoHFPhi_branch == 0 ) {
	cout << "Branch evt_metNoHFPhi does not exist." << endl;
	}
	evt_metNoHFSig_branch = 0;
	if (tree->GetAlias("evt_metNoHFSig") != 0) {
		evt_metNoHFSig_branch = tree->GetBranch(tree->GetAlias("evt_metNoHFSig"));
		evt_metNoHFSig_branch->SetAddress(&evt_metNoHFSig_);
	}
	if(evt_metNoHFSig_branch == 0 ) {
	cout << "Branch evt_metNoHFSig does not exist." << endl;
	}
	evt_metOpt_branch = 0;
	if (tree->GetAlias("evt_metOpt") != 0) {
		evt_metOpt_branch = tree->GetBranch(tree->GetAlias("evt_metOpt"));
		evt_metOpt_branch->SetAddress(&evt_metOpt_);
	}
	if(evt_metOpt_branch == 0 ) {
	cout << "Branch evt_metOpt does not exist." << endl;
	}
	evt_metOptHO_branch = 0;
	if (tree->GetAlias("evt_metOptHO") != 0) {
		evt_metOptHO_branch = tree->GetBranch(tree->GetAlias("evt_metOptHO"));
		evt_metOptHO_branch->SetAddress(&evt_metOptHO_);
	}
	if(evt_metOptHO_branch == 0 ) {
	cout << "Branch evt_metOptHO does not exist." << endl;
	}
	evt_metOptHOPhi_branch = 0;
	if (tree->GetAlias("evt_metOptHOPhi") != 0) {
		evt_metOptHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptHOPhi"));
		evt_metOptHOPhi_branch->SetAddress(&evt_metOptHOPhi_);
	}
	if(evt_metOptHOPhi_branch == 0 ) {
	cout << "Branch evt_metOptHOPhi does not exist." << endl;
	}
	evt_metOptHOSig_branch = 0;
	if (tree->GetAlias("evt_metOptHOSig") != 0) {
		evt_metOptHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptHOSig"));
		evt_metOptHOSig_branch->SetAddress(&evt_metOptHOSig_);
	}
	if(evt_metOptHOSig_branch == 0 ) {
	cout << "Branch evt_metOptHOSig does not exist." << endl;
	}
	evt_metOptNoHF_branch = 0;
	if (tree->GetAlias("evt_metOptNoHF") != 0) {
		evt_metOptNoHF_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHF"));
		evt_metOptNoHF_branch->SetAddress(&evt_metOptNoHF_);
	}
	if(evt_metOptNoHF_branch == 0 ) {
	cout << "Branch evt_metOptNoHF does not exist." << endl;
	}
	evt_metOptNoHFHO_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHO") != 0) {
		evt_metOptNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHO"));
		evt_metOptNoHFHO_branch->SetAddress(&evt_metOptNoHFHO_);
	}
	if(evt_metOptNoHFHO_branch == 0 ) {
	cout << "Branch evt_metOptNoHFHO does not exist." << endl;
	}
	evt_metOptNoHFHOPhi_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHOPhi") != 0) {
		evt_metOptNoHFHOPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHOPhi"));
		evt_metOptNoHFHOPhi_branch->SetAddress(&evt_metOptNoHFHOPhi_);
	}
	if(evt_metOptNoHFHOPhi_branch == 0 ) {
	cout << "Branch evt_metOptNoHFHOPhi does not exist." << endl;
	}
	evt_metOptNoHFHOSig_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFHOSig") != 0) {
		evt_metOptNoHFHOSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFHOSig"));
		evt_metOptNoHFHOSig_branch->SetAddress(&evt_metOptNoHFHOSig_);
	}
	if(evt_metOptNoHFHOSig_branch == 0 ) {
	cout << "Branch evt_metOptNoHFHOSig does not exist." << endl;
	}
	evt_metOptNoHFPhi_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFPhi") != 0) {
		evt_metOptNoHFPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFPhi"));
		evt_metOptNoHFPhi_branch->SetAddress(&evt_metOptNoHFPhi_);
	}
	if(evt_metOptNoHFPhi_branch == 0 ) {
	cout << "Branch evt_metOptNoHFPhi does not exist." << endl;
	}
	evt_metOptNoHFSig_branch = 0;
	if (tree->GetAlias("evt_metOptNoHFSig") != 0) {
		evt_metOptNoHFSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptNoHFSig"));
		evt_metOptNoHFSig_branch->SetAddress(&evt_metOptNoHFSig_);
	}
	if(evt_metOptNoHFSig_branch == 0 ) {
	cout << "Branch evt_metOptNoHFSig does not exist." << endl;
	}
	evt_metOptPhi_branch = 0;
	if (tree->GetAlias("evt_metOptPhi") != 0) {
		evt_metOptPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptPhi"));
		evt_metOptPhi_branch->SetAddress(&evt_metOptPhi_);
	}
	if(evt_metOptPhi_branch == 0 ) {
	cout << "Branch evt_metOptPhi does not exist." << endl;
	}
	evt_metOptSig_branch = 0;
	if (tree->GetAlias("evt_metOptSig") != 0) {
		evt_metOptSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptSig"));
		evt_metOptSig_branch->SetAddress(&evt_metOptSig_);
	}
	if(evt_metOptSig_branch == 0 ) {
	cout << "Branch evt_metOptSig does not exist." << endl;
	}
	evt_metPhi_branch = 0;
	if (tree->GetAlias("evt_metPhi") != 0) {
		evt_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_metPhi"));
		evt_metPhi_branch->SetAddress(&evt_metPhi_);
	}
	if(evt_metPhi_branch == 0 ) {
	cout << "Branch evt_metPhi does not exist." << endl;
	}
	evt_metSig_branch = 0;
	if (tree->GetAlias("evt_metSig") != 0) {
		evt_metSig_branch = tree->GetBranch(tree->GetAlias("evt_metSig"));
		evt_metSig_branch->SetAddress(&evt_metSig_);
	}
	if(evt_metSig_branch == 0 ) {
	cout << "Branch evt_metSig does not exist." << endl;
	}
	evt_sumet_branch = 0;
	if (tree->GetAlias("evt_sumet") != 0) {
		evt_sumet_branch = tree->GetBranch(tree->GetAlias("evt_sumet"));
		evt_sumet_branch->SetAddress(&evt_sumet_);
	}
	if(evt_sumet_branch == 0 ) {
	cout << "Branch evt_sumet does not exist." << endl;
	}
	evt_sumetHO_branch = 0;
	if (tree->GetAlias("evt_sumetHO") != 0) {
		evt_sumetHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetHO"));
		evt_sumetHO_branch->SetAddress(&evt_sumetHO_);
	}
	if(evt_sumetHO_branch == 0 ) {
	cout << "Branch evt_sumetHO does not exist." << endl;
	}
	evt_sumetMuonCorr_branch = 0;
	if (tree->GetAlias("evt_sumetMuonCorr") != 0) {
		evt_sumetMuonCorr_branch = tree->GetBranch(tree->GetAlias("evt_sumetMuonCorr"));
		evt_sumetMuonCorr_branch->SetAddress(&evt_sumetMuonCorr_);
	}
	if(evt_sumetMuonCorr_branch == 0 ) {
	cout << "Branch evt_sumetMuonCorr does not exist." << endl;
	}
	evt_sumetNoHF_branch = 0;
	if (tree->GetAlias("evt_sumetNoHF") != 0) {
		evt_sumetNoHF_branch = tree->GetBranch(tree->GetAlias("evt_sumetNoHF"));
		evt_sumetNoHF_branch->SetAddress(&evt_sumetNoHF_);
	}
	if(evt_sumetNoHF_branch == 0 ) {
	cout << "Branch evt_sumetNoHF does not exist." << endl;
	}
	evt_sumetNoHFHO_branch = 0;
	if (tree->GetAlias("evt_sumetNoHFHO") != 0) {
		evt_sumetNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetNoHFHO"));
		evt_sumetNoHFHO_branch->SetAddress(&evt_sumetNoHFHO_);
	}
	if(evt_sumetNoHFHO_branch == 0 ) {
	cout << "Branch evt_sumetNoHFHO does not exist." << endl;
	}
	evt_sumetOpt_branch = 0;
	if (tree->GetAlias("evt_sumetOpt") != 0) {
		evt_sumetOpt_branch = tree->GetBranch(tree->GetAlias("evt_sumetOpt"));
		evt_sumetOpt_branch->SetAddress(&evt_sumetOpt_);
	}
	if(evt_sumetOpt_branch == 0 ) {
	cout << "Branch evt_sumetOpt does not exist." << endl;
	}
	evt_sumetOptHO_branch = 0;
	if (tree->GetAlias("evt_sumetOptHO") != 0) {
		evt_sumetOptHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptHO"));
		evt_sumetOptHO_branch->SetAddress(&evt_sumetOptHO_);
	}
	if(evt_sumetOptHO_branch == 0 ) {
	cout << "Branch evt_sumetOptHO does not exist." << endl;
	}
	evt_sumetOptNoHF_branch = 0;
	if (tree->GetAlias("evt_sumetOptNoHF") != 0) {
		evt_sumetOptNoHF_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptNoHF"));
		evt_sumetOptNoHF_branch->SetAddress(&evt_sumetOptNoHF_);
	}
	if(evt_sumetOptNoHF_branch == 0 ) {
	cout << "Branch evt_sumetOptNoHF does not exist." << endl;
	}
	evt_sumetOptNoHFHO_branch = 0;
	if (tree->GetAlias("evt_sumetOptNoHFHO") != 0) {
		evt_sumetOptNoHFHO_branch = tree->GetBranch(tree->GetAlias("evt_sumetOptNoHFHO"));
		evt_sumetOptNoHFHO_branch->SetAddress(&evt_sumetOptNoHFHO_);
	}
	if(evt_sumetOptNoHFHO_branch == 0 ) {
	cout << "Branch evt_sumetOptNoHFHO does not exist." << endl;
	}
	met_pat_metCor_branch = 0;
	if (tree->GetAlias("met_pat_metCor") != 0) {
		met_pat_metCor_branch = tree->GetBranch(tree->GetAlias("met_pat_metCor"));
		met_pat_metCor_branch->SetAddress(&met_pat_metCor_);
	}
	if(met_pat_metCor_branch == 0 ) {
	cout << "Branch met_pat_metCor does not exist." << endl;
	}
	met_pat_metPhiCor_branch = 0;
	if (tree->GetAlias("met_pat_metPhiCor") != 0) {
		met_pat_metPhiCor_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiCor"));
		met_pat_metPhiCor_branch->SetAddress(&met_pat_metPhiCor_);
	}
	if(met_pat_metPhiCor_branch == 0 ) {
	cout << "Branch met_pat_metPhiCor does not exist." << endl;
	}
	met_pat_metPhiUncor_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncor") != 0) {
		met_pat_metPhiUncor_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncor"));
		met_pat_metPhiUncor_branch->SetAddress(&met_pat_metPhiUncor_);
	}
	if(met_pat_metPhiUncor_branch == 0 ) {
	cout << "Branch met_pat_metPhiUncor does not exist." << endl;
	}
	met_pat_metPhiUncorJES_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncorJES") != 0) {
		met_pat_metPhiUncorJES_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncorJES"));
		met_pat_metPhiUncorJES_branch->SetAddress(&met_pat_metPhiUncorJES_);
	}
	if(met_pat_metPhiUncorJES_branch == 0 ) {
	cout << "Branch met_pat_metPhiUncorJES does not exist." << endl;
	}
	met_pat_metPhiUncorMuon_branch = 0;
	if (tree->GetAlias("met_pat_metPhiUncorMuon") != 0) {
		met_pat_metPhiUncorMuon_branch = tree->GetBranch(tree->GetAlias("met_pat_metPhiUncorMuon"));
		met_pat_metPhiUncorMuon_branch->SetAddress(&met_pat_metPhiUncorMuon_);
	}
	if(met_pat_metPhiUncorMuon_branch == 0 ) {
	cout << "Branch met_pat_metPhiUncorMuon does not exist." << endl;
	}
	met_pat_metUncor_branch = 0;
	if (tree->GetAlias("met_pat_metUncor") != 0) {
		met_pat_metUncor_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncor"));
		met_pat_metUncor_branch->SetAddress(&met_pat_metUncor_);
	}
	if(met_pat_metUncor_branch == 0 ) {
	cout << "Branch met_pat_metUncor does not exist." << endl;
	}
	met_pat_metUncorJES_branch = 0;
	if (tree->GetAlias("met_pat_metUncorJES") != 0) {
		met_pat_metUncorJES_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncorJES"));
		met_pat_metUncorJES_branch->SetAddress(&met_pat_metUncorJES_);
	}
	if(met_pat_metUncorJES_branch == 0 ) {
	cout << "Branch met_pat_metUncorJES does not exist." << endl;
	}
	met_pat_metUncorMuon_branch = 0;
	if (tree->GetAlias("met_pat_metUncorMuon") != 0) {
		met_pat_metUncorMuon_branch = tree->GetBranch(tree->GetAlias("met_pat_metUncorMuon"));
		met_pat_metUncorMuon_branch->SetAddress(&met_pat_metUncorMuon_);
	}
	if(met_pat_metUncorMuon_branch == 0 ) {
	cout << "Branch met_pat_metUncorMuon does not exist." << endl;
	}
	pdfinfo_scale_branch = 0;
	if (tree->GetAlias("pdfinfo_scale") != 0) {
		pdfinfo_scale_branch = tree->GetBranch(tree->GetAlias("pdfinfo_scale"));
		pdfinfo_scale_branch->SetAddress(&pdfinfo_scale_);
	}
	if(pdfinfo_scale_branch == 0 ) {
	cout << "Branch pdfinfo_scale does not exist." << endl;
	}
	pdfinfo_x1_branch = 0;
	if (tree->GetAlias("pdfinfo_x1") != 0) {
		pdfinfo_x1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x1"));
		pdfinfo_x1_branch->SetAddress(&pdfinfo_x1_);
	}
	if(pdfinfo_x1_branch == 0 ) {
	cout << "Branch pdfinfo_x1 does not exist." << endl;
	}
	pdfinfo_x2_branch = 0;
	if (tree->GetAlias("pdfinfo_x2") != 0) {
		pdfinfo_x2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_x2"));
		pdfinfo_x2_branch->SetAddress(&pdfinfo_x2_);
	}
	if(pdfinfo_x2_branch == 0 ) {
	cout << "Branch pdfinfo_x2 does not exist." << endl;
	}
	evt_pfmet_branch = 0;
	if (tree->GetAlias("evt_pfmet") != 0) {
		evt_pfmet_branch = tree->GetBranch(tree->GetAlias("evt_pfmet"));
		evt_pfmet_branch->SetAddress(&evt_pfmet_);
	}
	if(evt_pfmet_branch == 0 ) {
	cout << "Branch evt_pfmet does not exist." << endl;
	}
	evt_pfmetPhi_branch = 0;
	if (tree->GetAlias("evt_pfmetPhi") != 0) {
		evt_pfmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_pfmetPhi"));
		evt_pfmetPhi_branch->SetAddress(&evt_pfmetPhi_);
	}
	if(evt_pfmetPhi_branch == 0 ) {
	cout << "Branch evt_pfmetPhi does not exist." << endl;
	}
	evt_pfsumet_branch = 0;
	if (tree->GetAlias("evt_pfsumet") != 0) {
		evt_pfsumet_branch = tree->GetBranch(tree->GetAlias("evt_pfsumet"));
		evt_pfsumet_branch->SetAddress(&evt_pfsumet_);
	}
	if(evt_pfsumet_branch == 0 ) {
	cout << "Branch evt_pfsumet does not exist." << endl;
	}
	evt_tcmet_branch = 0;
	if (tree->GetAlias("evt_tcmet") != 0) {
		evt_tcmet_branch = tree->GetBranch(tree->GetAlias("evt_tcmet"));
		evt_tcmet_branch->SetAddress(&evt_tcmet_);
	}
	if(evt_tcmet_branch == 0 ) {
	cout << "Branch evt_tcmet does not exist." << endl;
	}
	evt_tcmetPhi_branch = 0;
	if (tree->GetAlias("evt_tcmetPhi") != 0) {
		evt_tcmetPhi_branch = tree->GetBranch(tree->GetAlias("evt_tcmetPhi"));
		evt_tcmetPhi_branch->SetAddress(&evt_tcmetPhi_);
	}
	if(evt_tcmetPhi_branch == 0 ) {
	cout << "Branch evt_tcmetPhi does not exist." << endl;
	}
	evt_tcsumet_branch = 0;
	if (tree->GetAlias("evt_tcsumet") != 0) {
		evt_tcsumet_branch = tree->GetBranch(tree->GetAlias("evt_tcsumet"));
		evt_tcsumet_branch->SetAddress(&evt_tcsumet_);
	}
	if(evt_tcsumet_branch == 0 ) {
	cout << "Branch evt_tcsumet does not exist." << endl;
	}
	jets_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_combinedSecondaryVertexBJetTag") != 0) {
		jets_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_combinedSecondaryVertexBJetTag"));
		jets_combinedSecondaryVertexBJetTag_branch->SetAddress(&jets_combinedSecondaryVertexBJetTag_);
	}
	if(jets_combinedSecondaryVertexBJetTag_branch == 0 ) {
	cout << "Branch jets_combinedSecondaryVertexBJetTag does not exist." << endl;
	}
	jets_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_combinedSecondaryVertexMVABJetTag") != 0) {
		jets_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_combinedSecondaryVertexMVABJetTag"));
		jets_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&jets_combinedSecondaryVertexMVABJetTag_);
	}
	if(jets_combinedSecondaryVertexMVABJetTag_branch == 0 ) {
	cout << "Branch jets_combinedSecondaryVertexMVABJetTag does not exist." << endl;
	}
	jets_impactParameterMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_impactParameterMVABJetTag") != 0) {
		jets_impactParameterMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_impactParameterMVABJetTag"));
		jets_impactParameterMVABJetTag_branch->SetAddress(&jets_impactParameterMVABJetTag_);
	}
	if(jets_impactParameterMVABJetTag_branch == 0 ) {
	cout << "Branch jets_impactParameterMVABJetTag does not exist." << endl;
	}
	jets_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_jetBProbabilityBJetTag") != 0) {
		jets_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_jetBProbabilityBJetTag"));
		jets_jetBProbabilityBJetTag_branch->SetAddress(&jets_jetBProbabilityBJetTag_);
	}
	if(jets_jetBProbabilityBJetTag_branch == 0 ) {
	cout << "Branch jets_jetBProbabilityBJetTag does not exist." << endl;
	}
	jets_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_jetProbabilityBJetTag") != 0) {
		jets_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_jetProbabilityBJetTag"));
		jets_jetProbabilityBJetTag_branch->SetAddress(&jets_jetProbabilityBJetTag_);
	}
	if(jets_jetProbabilityBJetTag_branch == 0 ) {
	cout << "Branch jets_jetProbabilityBJetTag does not exist." << endl;
	}
	jets_simpleSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_simpleSecondaryVertexBJetTag") != 0) {
		jets_simpleSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_simpleSecondaryVertexBJetTag"));
		jets_simpleSecondaryVertexBJetTag_branch->SetAddress(&jets_simpleSecondaryVertexBJetTag_);
	}
	if(jets_simpleSecondaryVertexBJetTag_branch == 0 ) {
	cout << "Branch jets_simpleSecondaryVertexBJetTag does not exist." << endl;
	}
	jets_softElectronBJetTag_branch = 0;
	if (tree->GetAlias("jets_softElectronBJetTag") != 0) {
		jets_softElectronBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softElectronBJetTag"));
		jets_softElectronBJetTag_branch->SetAddress(&jets_softElectronBJetTag_);
	}
	if(jets_softElectronBJetTag_branch == 0 ) {
	cout << "Branch jets_softElectronBJetTag does not exist." << endl;
	}
	jets_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("jets_softMuonBJetTag") != 0) {
		jets_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softMuonBJetTag"));
		jets_softMuonBJetTag_branch->SetAddress(&jets_softMuonBJetTag_);
	}
	if(jets_softMuonBJetTag_branch == 0 ) {
	cout << "Branch jets_softMuonBJetTag does not exist." << endl;
	}
	jets_softMuonNoIPBJetTag_branch = 0;
	if (tree->GetAlias("jets_softMuonNoIPBJetTag") != 0) {
		jets_softMuonNoIPBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_softMuonNoIPBJetTag"));
		jets_softMuonNoIPBJetTag_branch->SetAddress(&jets_softMuonNoIPBJetTag_);
	}
	if(jets_softMuonNoIPBJetTag_branch == 0 ) {
	cout << "Branch jets_softMuonNoIPBJetTag does not exist." << endl;
	}
	jets_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_trackCountingHighEffBJetTag") != 0) {
		jets_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_trackCountingHighEffBJetTag"));
		jets_trackCountingHighEffBJetTag_branch->SetAddress(&jets_trackCountingHighEffBJetTag_);
	}
	if(jets_trackCountingHighEffBJetTag_branch == 0 ) {
	cout << "Branch jets_trackCountingHighEffBJetTag does not exist." << endl;
	}
	jets_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jets_trackCountingHighPurBJetTag") != 0) {
		jets_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_trackCountingHighPurBJetTag"));
		jets_trackCountingHighPurBJetTag_branch->SetAddress(&jets_trackCountingHighPurBJetTag_);
	}
	if(jets_trackCountingHighPurBJetTag_branch == 0 ) {
	cout << "Branch jets_trackCountingHighPurBJetTag does not exist." << endl;
	}
	taus_calo_isolationECALhitsEtSum_branch = 0;
	if (tree->GetAlias("taus_calo_isolationECALhitsEtSum") != 0) {
		taus_calo_isolationECALhitsEtSum_branch = tree->GetBranch(tree->GetAlias("taus_calo_isolationECALhitsEtSum"));
		taus_calo_isolationECALhitsEtSum_branch->SetAddress(&taus_calo_isolationECALhitsEtSum_);
	}
	if(taus_calo_isolationECALhitsEtSum_branch == 0 ) {
	cout << "Branch taus_calo_isolationECALhitsEtSum does not exist." << endl;
	}
	taus_calo_isolationtrksPtSum_branch = 0;
	if (tree->GetAlias("taus_calo_isolationtrksPtSum") != 0) {
		taus_calo_isolationtrksPtSum_branch = tree->GetBranch(tree->GetAlias("taus_calo_isolationtrksPtSum"));
		taus_calo_isolationtrksPtSum_branch->SetAddress(&taus_calo_isolationtrksPtSum_);
	}
	if(taus_calo_isolationtrksPtSum_branch == 0 ) {
	cout << "Branch taus_calo_isolationtrksPtSum does not exist." << endl;
	}
	taus_calo_leadtrk_HCAL3x3hitsEtSum_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_HCAL3x3hitsEtSum") != 0) {
		taus_calo_leadtrk_HCAL3x3hitsEtSum_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_HCAL3x3hitsEtSum"));
		taus_calo_leadtrk_HCAL3x3hitsEtSum_branch->SetAddress(&taus_calo_leadtrk_HCAL3x3hitsEtSum_);
	}
	if(taus_calo_leadtrk_HCAL3x3hitsEtSum_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_HCAL3x3hitsEtSum does not exist." << endl;
	}
	taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_HCAL3x3hottesthitDEta") != 0) {
		taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_HCAL3x3hottesthitDEta"));
		taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch->SetAddress(&taus_calo_leadtrk_HCAL3x3hottesthitDEta_);
	}
	if(taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_HCAL3x3hottesthitDEta does not exist." << endl;
	}
	taus_calo_leadtrk_Signed_Sipt_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_Signed_Sipt") != 0) {
		taus_calo_leadtrk_Signed_Sipt_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_Signed_Sipt"));
		taus_calo_leadtrk_Signed_Sipt_branch->SetAddress(&taus_calo_leadtrk_Signed_Sipt_);
	}
	if(taus_calo_leadtrk_Signed_Sipt_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_Signed_Sipt does not exist." << endl;
	}
	taus_calo_leadtrk_chi2_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_chi2") != 0) {
		taus_calo_leadtrk_chi2_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_chi2"));
		taus_calo_leadtrk_chi2_branch->SetAddress(&taus_calo_leadtrk_chi2_);
	}
	if(taus_calo_leadtrk_chi2_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_chi2 does not exist." << endl;
	}
	taus_calo_leadtrk_d0_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_d0") != 0) {
		taus_calo_leadtrk_d0_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_d0"));
		taus_calo_leadtrk_d0_branch->SetAddress(&taus_calo_leadtrk_d0_);
	}
	if(taus_calo_leadtrk_d0_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_d0 does not exist." << endl;
	}
	taus_calo_leadtrk_lostHits_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_lostHits") != 0) {
		taus_calo_leadtrk_lostHits_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_lostHits"));
		taus_calo_leadtrk_lostHits_branch->SetAddress(&taus_calo_leadtrk_lostHits_);
	}
	if(taus_calo_leadtrk_lostHits_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_lostHits does not exist." << endl;
	}
	taus_calo_leadtrk_ndof_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_ndof") != 0) {
		taus_calo_leadtrk_ndof_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_ndof"));
		taus_calo_leadtrk_ndof_branch->SetAddress(&taus_calo_leadtrk_ndof_);
	}
	if(taus_calo_leadtrk_ndof_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_ndof does not exist." << endl;
	}
	taus_calo_leadtrk_validHits_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_validHits") != 0) {
		taus_calo_leadtrk_validHits_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_validHits"));
		taus_calo_leadtrk_validHits_branch->SetAddress(&taus_calo_leadtrk_validHits_);
	}
	if(taus_calo_leadtrk_validHits_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_validHits does not exist." << endl;
	}
	taus_calo_leadtrk_z0_branch = 0;
	if (tree->GetAlias("taus_calo_leadtrk_z0") != 0) {
		taus_calo_leadtrk_z0_branch = tree->GetBranch(tree->GetAlias("taus_calo_leadtrk_z0"));
		taus_calo_leadtrk_z0_branch->SetAddress(&taus_calo_leadtrk_z0_);
	}
	if(taus_calo_leadtrk_z0_branch == 0 ) {
	cout << "Branch taus_calo_leadtrk_z0 does not exist." << endl;
	}
	taus_calo_maximumHCALhitEt_branch = 0;
	if (tree->GetAlias("taus_calo_maximumHCALhitEt") != 0) {
		taus_calo_maximumHCALhitEt_branch = tree->GetBranch(tree->GetAlias("taus_calo_maximumHCALhitEt"));
		taus_calo_maximumHCALhitEt_branch->SetAddress(&taus_calo_maximumHCALhitEt_);
	}
	if(taus_calo_maximumHCALhitEt_branch == 0 ) {
	cout << "Branch taus_calo_maximumHCALhitEt does not exist." << endl;
	}
	taus_calo_signaltrksInvariantMass_branch = 0;
	if (tree->GetAlias("taus_calo_signaltrksInvariantMass") != 0) {
		taus_calo_signaltrksInvariantMass_branch = tree->GetBranch(tree->GetAlias("taus_calo_signaltrksInvariantMass"));
		taus_calo_signaltrksInvariantMass_branch->SetAddress(&taus_calo_signaltrksInvariantMass_);
	}
	if(taus_calo_signaltrksInvariantMass_branch == 0 ) {
	cout << "Branch taus_calo_signaltrksInvariantMass does not exist." << endl;
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
	els_ecalJuraIso_branch = 0;
	if (tree->GetAlias("els_ecalJuraIso") != 0) {
		els_ecalJuraIso_branch = tree->GetBranch(tree->GetAlias("els_ecalJuraIso"));
		els_ecalJuraIso_branch->SetAddress(&els_ecalJuraIso_);
	}
	if(els_ecalJuraIso_branch == 0 ) {
	cout << "Branch els_ecalJuraIso does not exist." << endl;
	}
	els_ecalJuraTowerIso_branch = 0;
	if (tree->GetAlias("els_ecalJuraTowerIso") != 0) {
		els_ecalJuraTowerIso_branch = tree->GetBranch(tree->GetAlias("els_ecalJuraTowerIso"));
		els_ecalJuraTowerIso_branch->SetAddress(&els_ecalJuraTowerIso_);
	}
	if(els_ecalJuraTowerIso_branch == 0 ) {
	cout << "Branch els_ecalJuraTowerIso does not exist." << endl;
	}
	els_hcalConeIso_branch = 0;
	if (tree->GetAlias("els_hcalConeIso") != 0) {
		els_hcalConeIso_branch = tree->GetBranch(tree->GetAlias("els_hcalConeIso"));
		els_hcalConeIso_branch->SetAddress(&els_hcalConeIso_);
	}
	if(els_hcalConeIso_branch == 0 ) {
	cout << "Branch els_hcalConeIso does not exist." << endl;
	}
	els_jetdr_branch = 0;
	if (tree->GetAlias("els_jetdr") != 0) {
		els_jetdr_branch = tree->GetBranch(tree->GetAlias("els_jetdr"));
		els_jetdr_branch->SetAddress(&els_jetdr_);
	}
	if(els_jetdr_branch == 0 ) {
	cout << "Branch els_jetdr does not exist." << endl;
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
	els_ecalIso_branch = 0;
	if (tree->GetAlias("els_ecalIso") != 0) {
		els_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_ecalIso"));
		els_ecalIso_branch->SetAddress(&els_ecalIso_);
	}
	if(els_ecalIso_branch == 0 ) {
	cout << "Branch els_ecalIso does not exist." << endl;
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
	els_hcalIso_branch = 0;
	if (tree->GetAlias("els_hcalIso") != 0) {
		els_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_hcalIso"));
		els_hcalIso_branch->SetAddress(&els_hcalIso_);
	}
	if(els_hcalIso_branch == 0 ) {
	cout << "Branch els_hcalIso does not exist." << endl;
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
	els_tkIso_branch = 0;
	if (tree->GetAlias("els_tkIso") != 0) {
		els_tkIso_branch = tree->GetBranch(tree->GetAlias("els_tkIso"));
		els_tkIso_branch->SetAddress(&els_tkIso_);
	}
	if(els_tkIso_branch == 0 ) {
	cout << "Branch els_tkIso does not exist." << endl;
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
	hyp_Ht_branch = 0;
	if (tree->GetAlias("hyp_Ht") != 0) {
		hyp_Ht_branch = tree->GetBranch(tree->GetAlias("hyp_Ht"));
		hyp_Ht_branch->SetAddress(&hyp_Ht_);
	}
	if(hyp_Ht_branch == 0 ) {
	cout << "Branch hyp_Ht does not exist." << endl;
	}
	hyp_dPhi_nJet_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_muCorrMet") != 0) {
		hyp_dPhi_nJet_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_muCorrMet"));
		hyp_dPhi_nJet_muCorrMet_branch->SetAddress(&hyp_dPhi_nJet_muCorrMet_);
	}
	if(hyp_dPhi_nJet_muCorrMet_branch == 0 ) {
	cout << "Branch hyp_dPhi_nJet_muCorrMet does not exist." << endl;
	}
	hyp_dPhi_nJet_tcMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_tcMet") != 0) {
		hyp_dPhi_nJet_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_tcMet"));
		hyp_dPhi_nJet_tcMet_branch->SetAddress(&hyp_dPhi_nJet_tcMet_);
	}
	if(hyp_dPhi_nJet_tcMet_branch == 0 ) {
	cout << "Branch hyp_dPhi_nJet_tcMet does not exist." << endl;
	}
	hyp_dPhi_nJet_type1Met_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_type1Met") != 0) {
		hyp_dPhi_nJet_type1Met_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_type1Met"));
		hyp_dPhi_nJet_type1Met_branch->SetAddress(&hyp_dPhi_nJet_type1Met_);
	}
	if(hyp_dPhi_nJet_type1Met_branch == 0 ) {
	cout << "Branch hyp_dPhi_nJet_type1Met does not exist." << endl;
	}
	hyp_dPhi_nJet_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_dPhi_nJet_unCorrMet") != 0) {
		hyp_dPhi_nJet_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_dPhi_nJet_unCorrMet"));
		hyp_dPhi_nJet_unCorrMet_branch->SetAddress(&hyp_dPhi_nJet_unCorrMet_);
	}
	if(hyp_dPhi_nJet_unCorrMet_branch == 0 ) {
	cout << "Branch hyp_dPhi_nJet_unCorrMet does not exist." << endl;
	}
	hyp_ll_chi2_branch = 0;
	if (tree->GetAlias("hyp_ll_chi2") != 0) {
		hyp_ll_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_ll_chi2"));
		hyp_ll_chi2_branch->SetAddress(&hyp_ll_chi2_);
	}
	if(hyp_ll_chi2_branch == 0 ) {
	cout << "Branch hyp_ll_chi2 does not exist." << endl;
	}
	hyp_ll_d0_branch = 0;
	if (tree->GetAlias("hyp_ll_d0") != 0) {
		hyp_ll_d0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0"));
		hyp_ll_d0_branch->SetAddress(&hyp_ll_d0_);
	}
	if(hyp_ll_d0_branch == 0 ) {
	cout << "Branch hyp_ll_d0 does not exist." << endl;
	}
	hyp_ll_d0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_d0Err") != 0) {
		hyp_ll_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0Err"));
		hyp_ll_d0Err_branch->SetAddress(&hyp_ll_d0Err_);
	}
	if(hyp_ll_d0Err_branch == 0 ) {
	cout << "Branch hyp_ll_d0Err does not exist." << endl;
	}
	hyp_ll_d0corr_branch = 0;
	if (tree->GetAlias("hyp_ll_d0corr") != 0) {
		hyp_ll_d0corr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0corr"));
		hyp_ll_d0corr_branch->SetAddress(&hyp_ll_d0corr_);
	}
	if(hyp_ll_d0corr_branch == 0 ) {
	cout << "Branch hyp_ll_d0corr does not exist." << endl;
	}
	hyp_ll_dPhi_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_muCorrMet") != 0) {
		hyp_ll_dPhi_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_muCorrMet"));
		hyp_ll_dPhi_muCorrMet_branch->SetAddress(&hyp_ll_dPhi_muCorrMet_);
	}
	if(hyp_ll_dPhi_muCorrMet_branch == 0 ) {
	cout << "Branch hyp_ll_dPhi_muCorrMet does not exist." << endl;
	}
	hyp_ll_dPhi_tcMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_tcMet") != 0) {
		hyp_ll_dPhi_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_tcMet"));
		hyp_ll_dPhi_tcMet_branch->SetAddress(&hyp_ll_dPhi_tcMet_);
	}
	if(hyp_ll_dPhi_tcMet_branch == 0 ) {
	cout << "Branch hyp_ll_dPhi_tcMet does not exist." << endl;
	}
	hyp_ll_dPhi_type1Met_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_type1Met") != 0) {
		hyp_ll_dPhi_type1Met_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_type1Met"));
		hyp_ll_dPhi_type1Met_branch->SetAddress(&hyp_ll_dPhi_type1Met_);
	}
	if(hyp_ll_dPhi_type1Met_branch == 0 ) {
	cout << "Branch hyp_ll_dPhi_type1Met does not exist." << endl;
	}
	hyp_ll_dPhi_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_ll_dPhi_unCorrMet") != 0) {
		hyp_ll_dPhi_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_ll_dPhi_unCorrMet"));
		hyp_ll_dPhi_unCorrMet_branch->SetAddress(&hyp_ll_dPhi_unCorrMet_);
	}
	if(hyp_ll_dPhi_unCorrMet_branch == 0 ) {
	cout << "Branch hyp_ll_dPhi_unCorrMet does not exist." << endl;
	}
	hyp_ll_etaErr_branch = 0;
	if (tree->GetAlias("hyp_ll_etaErr") != 0) {
		hyp_ll_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_etaErr"));
		hyp_ll_etaErr_branch->SetAddress(&hyp_ll_etaErr_);
	}
	if(hyp_ll_etaErr_branch == 0 ) {
	cout << "Branch hyp_ll_etaErr does not exist." << endl;
	}
	hyp_ll_iso_branch = 0;
	if (tree->GetAlias("hyp_ll_iso") != 0) {
		hyp_ll_iso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_iso"));
		hyp_ll_iso_branch->SetAddress(&hyp_ll_iso_);
	}
	if(hyp_ll_iso_branch == 0 ) {
	cout << "Branch hyp_ll_iso does not exist." << endl;
	}
	hyp_ll_ndof_branch = 0;
	if (tree->GetAlias("hyp_ll_ndof") != 0) {
		hyp_ll_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ndof"));
		hyp_ll_ndof_branch->SetAddress(&hyp_ll_ndof_);
	}
	if(hyp_ll_ndof_branch == 0 ) {
	cout << "Branch hyp_ll_ndof does not exist." << endl;
	}
	hyp_ll_outerEta_branch = 0;
	if (tree->GetAlias("hyp_ll_outerEta") != 0) {
		hyp_ll_outerEta_branch = tree->GetBranch(tree->GetAlias("hyp_ll_outerEta"));
		hyp_ll_outerEta_branch->SetAddress(&hyp_ll_outerEta_);
	}
	if(hyp_ll_outerEta_branch == 0 ) {
	cout << "Branch hyp_ll_outerEta does not exist." << endl;
	}
	hyp_ll_outerPhi_branch = 0;
	if (tree->GetAlias("hyp_ll_outerPhi") != 0) {
		hyp_ll_outerPhi_branch = tree->GetBranch(tree->GetAlias("hyp_ll_outerPhi"));
		hyp_ll_outerPhi_branch->SetAddress(&hyp_ll_outerPhi_);
	}
	if(hyp_ll_outerPhi_branch == 0 ) {
	cout << "Branch hyp_ll_outerPhi does not exist." << endl;
	}
	hyp_ll_phiErr_branch = 0;
	if (tree->GetAlias("hyp_ll_phiErr") != 0) {
		hyp_ll_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_phiErr"));
		hyp_ll_phiErr_branch->SetAddress(&hyp_ll_phiErr_);
	}
	if(hyp_ll_phiErr_branch == 0 ) {
	cout << "Branch hyp_ll_phiErr does not exist." << endl;
	}
	hyp_ll_ptErr_branch = 0;
	if (tree->GetAlias("hyp_ll_ptErr") != 0) {
		hyp_ll_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ptErr"));
		hyp_ll_ptErr_branch->SetAddress(&hyp_ll_ptErr_);
	}
	if(hyp_ll_ptErr_branch == 0 ) {
	cout << "Branch hyp_ll_ptErr does not exist." << endl;
	}
	hyp_ll_tkIso_branch = 0;
	if (tree->GetAlias("hyp_ll_tkIso") != 0) {
		hyp_ll_tkIso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_tkIso"));
		hyp_ll_tkIso_branch->SetAddress(&hyp_ll_tkIso_);
	}
	if(hyp_ll_tkIso_branch == 0 ) {
	cout << "Branch hyp_ll_tkIso does not exist." << endl;
	}
	hyp_ll_vertexphi_branch = 0;
	if (tree->GetAlias("hyp_ll_vertexphi") != 0) {
		hyp_ll_vertexphi_branch = tree->GetBranch(tree->GetAlias("hyp_ll_vertexphi"));
		hyp_ll_vertexphi_branch->SetAddress(&hyp_ll_vertexphi_);
	}
	if(hyp_ll_vertexphi_branch == 0 ) {
	cout << "Branch hyp_ll_vertexphi does not exist." << endl;
	}
	hyp_ll_z0_branch = 0;
	if (tree->GetAlias("hyp_ll_z0") != 0) {
		hyp_ll_z0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0"));
		hyp_ll_z0_branch->SetAddress(&hyp_ll_z0_);
	}
	if(hyp_ll_z0_branch == 0 ) {
	cout << "Branch hyp_ll_z0 does not exist." << endl;
	}
	hyp_ll_z0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_z0Err") != 0) {
		hyp_ll_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0Err"));
		hyp_ll_z0Err_branch->SetAddress(&hyp_ll_z0Err_);
	}
	if(hyp_ll_z0Err_branch == 0 ) {
	cout << "Branch hyp_ll_z0Err does not exist." << endl;
	}
	hyp_ll_z0corr_branch = 0;
	if (tree->GetAlias("hyp_ll_z0corr") != 0) {
		hyp_ll_z0corr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0corr"));
		hyp_ll_z0corr_branch->SetAddress(&hyp_ll_z0corr_);
	}
	if(hyp_ll_z0corr_branch == 0 ) {
	cout << "Branch hyp_ll_z0corr does not exist." << endl;
	}
	hyp_lt_chi2_branch = 0;
	if (tree->GetAlias("hyp_lt_chi2") != 0) {
		hyp_lt_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_lt_chi2"));
		hyp_lt_chi2_branch->SetAddress(&hyp_lt_chi2_);
	}
	if(hyp_lt_chi2_branch == 0 ) {
	cout << "Branch hyp_lt_chi2 does not exist." << endl;
	}
	hyp_lt_d0_branch = 0;
	if (tree->GetAlias("hyp_lt_d0") != 0) {
		hyp_lt_d0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0"));
		hyp_lt_d0_branch->SetAddress(&hyp_lt_d0_);
	}
	if(hyp_lt_d0_branch == 0 ) {
	cout << "Branch hyp_lt_d0 does not exist." << endl;
	}
	hyp_lt_d0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_d0Err") != 0) {
		hyp_lt_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0Err"));
		hyp_lt_d0Err_branch->SetAddress(&hyp_lt_d0Err_);
	}
	if(hyp_lt_d0Err_branch == 0 ) {
	cout << "Branch hyp_lt_d0Err does not exist." << endl;
	}
	hyp_lt_d0corr_branch = 0;
	if (tree->GetAlias("hyp_lt_d0corr") != 0) {
		hyp_lt_d0corr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0corr"));
		hyp_lt_d0corr_branch->SetAddress(&hyp_lt_d0corr_);
	}
	if(hyp_lt_d0corr_branch == 0 ) {
	cout << "Branch hyp_lt_d0corr does not exist." << endl;
	}
	hyp_lt_dPhi_muCorrMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_muCorrMet") != 0) {
		hyp_lt_dPhi_muCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_muCorrMet"));
		hyp_lt_dPhi_muCorrMet_branch->SetAddress(&hyp_lt_dPhi_muCorrMet_);
	}
	if(hyp_lt_dPhi_muCorrMet_branch == 0 ) {
	cout << "Branch hyp_lt_dPhi_muCorrMet does not exist." << endl;
	}
	hyp_lt_dPhi_tcMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_tcMet") != 0) {
		hyp_lt_dPhi_tcMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_tcMet"));
		hyp_lt_dPhi_tcMet_branch->SetAddress(&hyp_lt_dPhi_tcMet_);
	}
	if(hyp_lt_dPhi_tcMet_branch == 0 ) {
	cout << "Branch hyp_lt_dPhi_tcMet does not exist." << endl;
	}
	hyp_lt_dPhi_type1Met_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_type1Met") != 0) {
		hyp_lt_dPhi_type1Met_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_type1Met"));
		hyp_lt_dPhi_type1Met_branch->SetAddress(&hyp_lt_dPhi_type1Met_);
	}
	if(hyp_lt_dPhi_type1Met_branch == 0 ) {
	cout << "Branch hyp_lt_dPhi_type1Met does not exist." << endl;
	}
	hyp_lt_dPhi_unCorrMet_branch = 0;
	if (tree->GetAlias("hyp_lt_dPhi_unCorrMet") != 0) {
		hyp_lt_dPhi_unCorrMet_branch = tree->GetBranch(tree->GetAlias("hyp_lt_dPhi_unCorrMet"));
		hyp_lt_dPhi_unCorrMet_branch->SetAddress(&hyp_lt_dPhi_unCorrMet_);
	}
	if(hyp_lt_dPhi_unCorrMet_branch == 0 ) {
	cout << "Branch hyp_lt_dPhi_unCorrMet does not exist." << endl;
	}
	hyp_lt_etaErr_branch = 0;
	if (tree->GetAlias("hyp_lt_etaErr") != 0) {
		hyp_lt_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_etaErr"));
		hyp_lt_etaErr_branch->SetAddress(&hyp_lt_etaErr_);
	}
	if(hyp_lt_etaErr_branch == 0 ) {
	cout << "Branch hyp_lt_etaErr does not exist." << endl;
	}
	hyp_lt_iso_branch = 0;
	if (tree->GetAlias("hyp_lt_iso") != 0) {
		hyp_lt_iso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_iso"));
		hyp_lt_iso_branch->SetAddress(&hyp_lt_iso_);
	}
	if(hyp_lt_iso_branch == 0 ) {
	cout << "Branch hyp_lt_iso does not exist." << endl;
	}
	hyp_lt_ndof_branch = 0;
	if (tree->GetAlias("hyp_lt_ndof") != 0) {
		hyp_lt_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ndof"));
		hyp_lt_ndof_branch->SetAddress(&hyp_lt_ndof_);
	}
	if(hyp_lt_ndof_branch == 0 ) {
	cout << "Branch hyp_lt_ndof does not exist." << endl;
	}
	hyp_lt_outerEta_branch = 0;
	if (tree->GetAlias("hyp_lt_outerEta") != 0) {
		hyp_lt_outerEta_branch = tree->GetBranch(tree->GetAlias("hyp_lt_outerEta"));
		hyp_lt_outerEta_branch->SetAddress(&hyp_lt_outerEta_);
	}
	if(hyp_lt_outerEta_branch == 0 ) {
	cout << "Branch hyp_lt_outerEta does not exist." << endl;
	}
	hyp_lt_outerPhi_branch = 0;
	if (tree->GetAlias("hyp_lt_outerPhi") != 0) {
		hyp_lt_outerPhi_branch = tree->GetBranch(tree->GetAlias("hyp_lt_outerPhi"));
		hyp_lt_outerPhi_branch->SetAddress(&hyp_lt_outerPhi_);
	}
	if(hyp_lt_outerPhi_branch == 0 ) {
	cout << "Branch hyp_lt_outerPhi does not exist." << endl;
	}
	hyp_lt_phiErr_branch = 0;
	if (tree->GetAlias("hyp_lt_phiErr") != 0) {
		hyp_lt_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_phiErr"));
		hyp_lt_phiErr_branch->SetAddress(&hyp_lt_phiErr_);
	}
	if(hyp_lt_phiErr_branch == 0 ) {
	cout << "Branch hyp_lt_phiErr does not exist." << endl;
	}
	hyp_lt_ptErr_branch = 0;
	if (tree->GetAlias("hyp_lt_ptErr") != 0) {
		hyp_lt_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ptErr"));
		hyp_lt_ptErr_branch->SetAddress(&hyp_lt_ptErr_);
	}
	if(hyp_lt_ptErr_branch == 0 ) {
	cout << "Branch hyp_lt_ptErr does not exist." << endl;
	}
	hyp_lt_tkIso_branch = 0;
	if (tree->GetAlias("hyp_lt_tkIso") != 0) {
		hyp_lt_tkIso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_tkIso"));
		hyp_lt_tkIso_branch->SetAddress(&hyp_lt_tkIso_);
	}
	if(hyp_lt_tkIso_branch == 0 ) {
	cout << "Branch hyp_lt_tkIso does not exist." << endl;
	}
	hyp_lt_vertexphi_branch = 0;
	if (tree->GetAlias("hyp_lt_vertexphi") != 0) {
		hyp_lt_vertexphi_branch = tree->GetBranch(tree->GetAlias("hyp_lt_vertexphi"));
		hyp_lt_vertexphi_branch->SetAddress(&hyp_lt_vertexphi_);
	}
	if(hyp_lt_vertexphi_branch == 0 ) {
	cout << "Branch hyp_lt_vertexphi does not exist." << endl;
	}
	hyp_lt_z0_branch = 0;
	if (tree->GetAlias("hyp_lt_z0") != 0) {
		hyp_lt_z0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0"));
		hyp_lt_z0_branch->SetAddress(&hyp_lt_z0_);
	}
	if(hyp_lt_z0_branch == 0 ) {
	cout << "Branch hyp_lt_z0 does not exist." << endl;
	}
	hyp_lt_z0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_z0Err") != 0) {
		hyp_lt_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0Err"));
		hyp_lt_z0Err_branch->SetAddress(&hyp_lt_z0Err_);
	}
	if(hyp_lt_z0Err_branch == 0 ) {
	cout << "Branch hyp_lt_z0Err does not exist." << endl;
	}
	hyp_lt_z0corr_branch = 0;
	if (tree->GetAlias("hyp_lt_z0corr") != 0) {
		hyp_lt_z0corr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0corr"));
		hyp_lt_z0corr_branch->SetAddress(&hyp_lt_z0corr_);
	}
	if(hyp_lt_z0corr_branch == 0 ) {
	cout << "Branch hyp_lt_z0corr does not exist." << endl;
	}
	hyp_sumJetPt_branch = 0;
	if (tree->GetAlias("hyp_sumJetPt") != 0) {
		hyp_sumJetPt_branch = tree->GetBranch(tree->GetAlias("hyp_sumJetPt"));
		hyp_sumJetPt_branch->SetAddress(&hyp_sumJetPt_);
	}
	if(hyp_sumJetPt_branch == 0 ) {
	cout << "Branch hyp_sumJetPt does not exist." << endl;
	}
	jets_EMFcor_branch = 0;
	if (tree->GetAlias("jets_EMFcor") != 0) {
		jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("jets_EMFcor"));
		jets_EMFcor_branch->SetAddress(&jets_EMFcor_);
	}
	if(jets_EMFcor_branch == 0 ) {
	cout << "Branch jets_EMFcor does not exist." << endl;
	}
	jets_chFrac_branch = 0;
	if (tree->GetAlias("jets_chFrac") != 0) {
		jets_chFrac_branch = tree->GetBranch(tree->GetAlias("jets_chFrac"));
		jets_chFrac_branch->SetAddress(&jets_chFrac_);
	}
	if(jets_chFrac_branch == 0 ) {
	cout << "Branch jets_chFrac does not exist." << endl;
	}
	jets_cor_branch = 0;
	if (tree->GetAlias("jets_cor") != 0) {
		jets_cor_branch = tree->GetBranch(tree->GetAlias("jets_cor"));
		jets_cor_branch->SetAddress(&jets_cor_);
	}
	if(jets_cor_branch == 0 ) {
	cout << "Branch jets_cor does not exist." << endl;
	}
	jets_emFrac_branch = 0;
	if (tree->GetAlias("jets_emFrac") != 0) {
		jets_emFrac_branch = tree->GetBranch(tree->GetAlias("jets_emFrac"));
		jets_emFrac_branch->SetAddress(&jets_emFrac_);
	}
	if(jets_emFrac_branch == 0 ) {
	cout << "Branch jets_emFrac does not exist." << endl;
	}
	jpts_cor_branch = 0;
	if (tree->GetAlias("jpts_cor") != 0) {
		jpts_cor_branch = tree->GetBranch(tree->GetAlias("jpts_cor"));
		jpts_cor_branch->SetAddress(&jpts_cor_);
	}
	if(jpts_cor_branch == 0 ) {
	cout << "Branch jpts_cor does not exist." << endl;
	}
	jpts_emFrac_branch = 0;
	if (tree->GetAlias("jpts_emFrac") != 0) {
		jpts_emFrac_branch = tree->GetBranch(tree->GetAlias("jpts_emFrac"));
		jpts_emFrac_branch->SetAddress(&jpts_emFrac_);
	}
	if(jpts_emFrac_branch == 0 ) {
	cout << "Branch jpts_emFrac does not exist." << endl;
	}
	jpts_jet_cor_branch = 0;
	if (tree->GetAlias("jpts_jet_cor") != 0) {
		jpts_jet_cor_branch = tree->GetBranch(tree->GetAlias("jpts_jet_cor"));
		jpts_jet_cor_branch->SetAddress(&jpts_jet_cor_);
	}
	if(jpts_jet_cor_branch == 0 ) {
	cout << "Branch jpts_jet_cor does not exist." << endl;
	}
	mus_met_deltax_branch = 0;
	if (tree->GetAlias("mus_met_deltax") != 0) {
		mus_met_deltax_branch = tree->GetBranch(tree->GetAlias("mus_met_deltax"));
		mus_met_deltax_branch->SetAddress(&mus_met_deltax_);
	}
	if(mus_met_deltax_branch == 0 ) {
	cout << "Branch mus_met_deltax does not exist." << endl;
	}
	mus_met_deltay_branch = 0;
	if (tree->GetAlias("mus_met_deltay") != 0) {
		mus_met_deltay_branch = tree->GetBranch(tree->GetAlias("mus_met_deltay"));
		mus_met_deltay_branch->SetAddress(&mus_met_deltay_);
	}
	if(mus_met_deltay_branch == 0 ) {
	cout << "Branch mus_met_deltay does not exist." << endl;
	}
	mus_eledr_branch = 0;
	if (tree->GetAlias("mus_eledr") != 0) {
		mus_eledr_branch = tree->GetBranch(tree->GetAlias("mus_eledr"));
		mus_eledr_branch->SetAddress(&mus_eledr_);
	}
	if(mus_eledr_branch == 0 ) {
	cout << "Branch mus_eledr does not exist." << endl;
	}
	mus_jetdr_branch = 0;
	if (tree->GetAlias("mus_jetdr") != 0) {
		mus_jetdr_branch = tree->GetBranch(tree->GetAlias("mus_jetdr"));
		mus_jetdr_branch->SetAddress(&mus_jetdr_);
	}
	if(mus_jetdr_branch == 0 ) {
	cout << "Branch mus_jetdr does not exist." << endl;
	}
	mus_trkdr_branch = 0;
	if (tree->GetAlias("mus_trkdr") != 0) {
		mus_trkdr_branch = tree->GetBranch(tree->GetAlias("mus_trkdr"));
		mus_trkdr_branch->SetAddress(&mus_trkdr_);
	}
	if(mus_trkdr_branch == 0 ) {
	cout << "Branch mus_trkdr does not exist." << endl;
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
	els_pat_caloIso_branch = 0;
	if (tree->GetAlias("els_pat_caloIso") != 0) {
		els_pat_caloIso_branch = tree->GetBranch(tree->GetAlias("els_pat_caloIso"));
		els_pat_caloIso_branch->SetAddress(&els_pat_caloIso_);
	}
	if(els_pat_caloIso_branch == 0 ) {
	cout << "Branch els_pat_caloIso does not exist." << endl;
	}
	els_pat_ecalIso_branch = 0;
	if (tree->GetAlias("els_pat_ecalIso") != 0) {
		els_pat_ecalIso_branch = tree->GetBranch(tree->GetAlias("els_pat_ecalIso"));
		els_pat_ecalIso_branch->SetAddress(&els_pat_ecalIso_);
	}
	if(els_pat_ecalIso_branch == 0 ) {
	cout << "Branch els_pat_ecalIso does not exist." << endl;
	}
	els_pat_hcalIso_branch = 0;
	if (tree->GetAlias("els_pat_hcalIso") != 0) {
		els_pat_hcalIso_branch = tree->GetBranch(tree->GetAlias("els_pat_hcalIso"));
		els_pat_hcalIso_branch->SetAddress(&els_pat_hcalIso_);
	}
	if(els_pat_hcalIso_branch == 0 ) {
	cout << "Branch els_pat_hcalIso does not exist." << endl;
	}
	els_pat_looseId_branch = 0;
	if (tree->GetAlias("els_pat_looseId") != 0) {
		els_pat_looseId_branch = tree->GetBranch(tree->GetAlias("els_pat_looseId"));
		els_pat_looseId_branch->SetAddress(&els_pat_looseId_);
	}
	if(els_pat_looseId_branch == 0 ) {
	cout << "Branch els_pat_looseId does not exist." << endl;
	}
	els_pat_robustHighEnergy_branch = 0;
	if (tree->GetAlias("els_pat_robustHighEnergy") != 0) {
		els_pat_robustHighEnergy_branch = tree->GetBranch(tree->GetAlias("els_pat_robustHighEnergy"));
		els_pat_robustHighEnergy_branch->SetAddress(&els_pat_robustHighEnergy_);
	}
	if(els_pat_robustHighEnergy_branch == 0 ) {
	cout << "Branch els_pat_robustHighEnergy does not exist." << endl;
	}
	els_pat_robustLooseId_branch = 0;
	if (tree->GetAlias("els_pat_robustLooseId") != 0) {
		els_pat_robustLooseId_branch = tree->GetBranch(tree->GetAlias("els_pat_robustLooseId"));
		els_pat_robustLooseId_branch->SetAddress(&els_pat_robustLooseId_);
	}
	if(els_pat_robustLooseId_branch == 0 ) {
	cout << "Branch els_pat_robustLooseId does not exist." << endl;
	}
	els_pat_robustTightId_branch = 0;
	if (tree->GetAlias("els_pat_robustTightId") != 0) {
		els_pat_robustTightId_branch = tree->GetBranch(tree->GetAlias("els_pat_robustTightId"));
		els_pat_robustTightId_branch->SetAddress(&els_pat_robustTightId_);
	}
	if(els_pat_robustTightId_branch == 0 ) {
	cout << "Branch els_pat_robustTightId does not exist." << endl;
	}
	els_pat_scE1x5_branch = 0;
	if (tree->GetAlias("els_pat_scE1x5") != 0) {
		els_pat_scE1x5_branch = tree->GetBranch(tree->GetAlias("els_pat_scE1x5"));
		els_pat_scE1x5_branch->SetAddress(&els_pat_scE1x5_);
	}
	if(els_pat_scE1x5_branch == 0 ) {
	cout << "Branch els_pat_scE1x5 does not exist." << endl;
	}
	els_pat_scE2x5Max_branch = 0;
	if (tree->GetAlias("els_pat_scE2x5Max") != 0) {
		els_pat_scE2x5Max_branch = tree->GetBranch(tree->GetAlias("els_pat_scE2x5Max"));
		els_pat_scE2x5Max_branch->SetAddress(&els_pat_scE2x5Max_);
	}
	if(els_pat_scE2x5Max_branch == 0 ) {
	cout << "Branch els_pat_scE2x5Max does not exist." << endl;
	}
	els_pat_scE5x5_branch = 0;
	if (tree->GetAlias("els_pat_scE5x5") != 0) {
		els_pat_scE5x5_branch = tree->GetBranch(tree->GetAlias("els_pat_scE5x5"));
		els_pat_scE5x5_branch->SetAddress(&els_pat_scE5x5_);
	}
	if(els_pat_scE5x5_branch == 0 ) {
	cout << "Branch els_pat_scE5x5 does not exist." << endl;
	}
	els_pat_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_pat_sigmaEtaEta") != 0) {
		els_pat_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_pat_sigmaEtaEta"));
		els_pat_sigmaEtaEta_branch->SetAddress(&els_pat_sigmaEtaEta_);
	}
	if(els_pat_sigmaEtaEta_branch == 0 ) {
	cout << "Branch els_pat_sigmaEtaEta does not exist." << endl;
	}
	els_pat_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("els_pat_sigmaIEtaIEta") != 0) {
		els_pat_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("els_pat_sigmaIEtaIEta"));
		els_pat_sigmaIEtaIEta_branch->SetAddress(&els_pat_sigmaIEtaIEta_);
	}
	if(els_pat_sigmaIEtaIEta_branch == 0 ) {
	cout << "Branch els_pat_sigmaIEtaIEta does not exist." << endl;
	}
	els_pat_tightId_branch = 0;
	if (tree->GetAlias("els_pat_tightId") != 0) {
		els_pat_tightId_branch = tree->GetBranch(tree->GetAlias("els_pat_tightId"));
		els_pat_tightId_branch->SetAddress(&els_pat_tightId_);
	}
	if(els_pat_tightId_branch == 0 ) {
	cout << "Branch els_pat_tightId does not exist." << endl;
	}
	els_pat_trackIso_branch = 0;
	if (tree->GetAlias("els_pat_trackIso") != 0) {
		els_pat_trackIso_branch = tree->GetBranch(tree->GetAlias("els_pat_trackIso"));
		els_pat_trackIso_branch->SetAddress(&els_pat_trackIso_);
	}
	if(els_pat_trackIso_branch == 0 ) {
	cout << "Branch els_pat_trackIso does not exist." << endl;
	}
	jets_pat_bCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_bCorrF") != 0) {
		jets_pat_bCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_bCorrF"));
		jets_pat_bCorrF_branch->SetAddress(&jets_pat_bCorrF_);
	}
	if(jets_pat_bCorrF_branch == 0 ) {
	cout << "Branch jets_pat_bCorrF does not exist." << endl;
	}
	jets_pat_cCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_cCorrF") != 0) {
		jets_pat_cCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_cCorrF"));
		jets_pat_cCorrF_branch->SetAddress(&jets_pat_cCorrF_);
	}
	if(jets_pat_cCorrF_branch == 0 ) {
	cout << "Branch jets_pat_cCorrF does not exist." << endl;
	}
	jets_pat_combinedSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_combinedSecondaryVertexBJetTag") != 0) {
		jets_pat_combinedSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_combinedSecondaryVertexBJetTag"));
		jets_pat_combinedSecondaryVertexBJetTag_branch->SetAddress(&jets_pat_combinedSecondaryVertexBJetTag_);
	}
	if(jets_pat_combinedSecondaryVertexBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_combinedSecondaryVertexBJetTag does not exist." << endl;
	}
	jets_pat_combinedSecondaryVertexMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_combinedSecondaryVertexMVABJetTag") != 0) {
		jets_pat_combinedSecondaryVertexMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_combinedSecondaryVertexMVABJetTag"));
		jets_pat_combinedSecondaryVertexMVABJetTag_branch->SetAddress(&jets_pat_combinedSecondaryVertexMVABJetTag_);
	}
	if(jets_pat_combinedSecondaryVertexMVABJetTag_branch == 0 ) {
	cout << "Branch jets_pat_combinedSecondaryVertexMVABJetTag does not exist." << endl;
	}
	jets_pat_coneIsolationTauJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_coneIsolationTauJetTag") != 0) {
		jets_pat_coneIsolationTauJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_coneIsolationTauJetTag"));
		jets_pat_coneIsolationTauJetTag_branch->SetAddress(&jets_pat_coneIsolationTauJetTag_);
	}
	if(jets_pat_coneIsolationTauJetTag_branch == 0 ) {
	cout << "Branch jets_pat_coneIsolationTauJetTag does not exist." << endl;
	}
	jets_pat_gluCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_gluCorrF") != 0) {
		jets_pat_gluCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_gluCorrF"));
		jets_pat_gluCorrF_branch->SetAddress(&jets_pat_gluCorrF_);
	}
	if(jets_pat_gluCorrF_branch == 0 ) {
	cout << "Branch jets_pat_gluCorrF does not exist." << endl;
	}
	jets_pat_impactParameterMVABJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_impactParameterMVABJetTag") != 0) {
		jets_pat_impactParameterMVABJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_impactParameterMVABJetTag"));
		jets_pat_impactParameterMVABJetTag_branch->SetAddress(&jets_pat_impactParameterMVABJetTag_);
	}
	if(jets_pat_impactParameterMVABJetTag_branch == 0 ) {
	cout << "Branch jets_pat_impactParameterMVABJetTag does not exist." << endl;
	}
	jets_pat_jetBProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_jetBProbabilityBJetTag") != 0) {
		jets_pat_jetBProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetBProbabilityBJetTag"));
		jets_pat_jetBProbabilityBJetTag_branch->SetAddress(&jets_pat_jetBProbabilityBJetTag_);
	}
	if(jets_pat_jetBProbabilityBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_jetBProbabilityBJetTag does not exist." << endl;
	}
	jets_pat_jetCharge_branch = 0;
	if (tree->GetAlias("jets_pat_jetCharge") != 0) {
		jets_pat_jetCharge_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetCharge"));
		jets_pat_jetCharge_branch->SetAddress(&jets_pat_jetCharge_);
	}
	if(jets_pat_jetCharge_branch == 0 ) {
	cout << "Branch jets_pat_jetCharge does not exist." << endl;
	}
	jets_pat_jetProbabilityBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_jetProbabilityBJetTag") != 0) {
		jets_pat_jetProbabilityBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetProbabilityBJetTag"));
		jets_pat_jetProbabilityBJetTag_branch->SetAddress(&jets_pat_jetProbabilityBJetTag_);
	}
	if(jets_pat_jetProbabilityBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_jetProbabilityBJetTag does not exist." << endl;
	}
	jets_pat_noCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_noCorrF") != 0) {
		jets_pat_noCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_noCorrF"));
		jets_pat_noCorrF_branch->SetAddress(&jets_pat_noCorrF_);
	}
	if(jets_pat_noCorrF_branch == 0 ) {
	cout << "Branch jets_pat_noCorrF does not exist." << endl;
	}
	jets_pat_simpleSecondaryVertexBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_simpleSecondaryVertexBJetTag") != 0) {
		jets_pat_simpleSecondaryVertexBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_simpleSecondaryVertexBJetTag"));
		jets_pat_simpleSecondaryVertexBJetTag_branch->SetAddress(&jets_pat_simpleSecondaryVertexBJetTag_);
	}
	if(jets_pat_simpleSecondaryVertexBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_simpleSecondaryVertexBJetTag does not exist." << endl;
	}
	jets_pat_softElectronBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softElectronBJetTag") != 0) {
		jets_pat_softElectronBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softElectronBJetTag"));
		jets_pat_softElectronBJetTag_branch->SetAddress(&jets_pat_softElectronBJetTag_);
	}
	if(jets_pat_softElectronBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_softElectronBJetTag does not exist." << endl;
	}
	jets_pat_softMuonBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softMuonBJetTag") != 0) {
		jets_pat_softMuonBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softMuonBJetTag"));
		jets_pat_softMuonBJetTag_branch->SetAddress(&jets_pat_softMuonBJetTag_);
	}
	if(jets_pat_softMuonBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_softMuonBJetTag does not exist." << endl;
	}
	jets_pat_softMuonNoIPBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_softMuonNoIPBJetTag") != 0) {
		jets_pat_softMuonNoIPBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_softMuonNoIPBJetTag"));
		jets_pat_softMuonNoIPBJetTag_branch->SetAddress(&jets_pat_softMuonNoIPBJetTag_);
	}
	if(jets_pat_softMuonNoIPBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_softMuonNoIPBJetTag does not exist." << endl;
	}
	jets_pat_trackCountingHighEffBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_trackCountingHighEffBJetTag") != 0) {
		jets_pat_trackCountingHighEffBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_trackCountingHighEffBJetTag"));
		jets_pat_trackCountingHighEffBJetTag_branch->SetAddress(&jets_pat_trackCountingHighEffBJetTag_);
	}
	if(jets_pat_trackCountingHighEffBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_trackCountingHighEffBJetTag does not exist." << endl;
	}
	jets_pat_trackCountingHighPurBJetTag_branch = 0;
	if (tree->GetAlias("jets_pat_trackCountingHighPurBJetTag") != 0) {
		jets_pat_trackCountingHighPurBJetTag_branch = tree->GetBranch(tree->GetAlias("jets_pat_trackCountingHighPurBJetTag"));
		jets_pat_trackCountingHighPurBJetTag_branch->SetAddress(&jets_pat_trackCountingHighPurBJetTag_);
	}
	if(jets_pat_trackCountingHighPurBJetTag_branch == 0 ) {
	cout << "Branch jets_pat_trackCountingHighPurBJetTag does not exist." << endl;
	}
	jets_pat_udsCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_udsCorrF") != 0) {
		jets_pat_udsCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_udsCorrF"));
		jets_pat_udsCorrF_branch->SetAddress(&jets_pat_udsCorrF_);
	}
	if(jets_pat_udsCorrF_branch == 0 ) {
	cout << "Branch jets_pat_udsCorrF does not exist." << endl;
	}
	mus_pat_caloIso_branch = 0;
	if (tree->GetAlias("mus_pat_caloIso") != 0) {
		mus_pat_caloIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_caloIso"));
		mus_pat_caloIso_branch->SetAddress(&mus_pat_caloIso_);
	}
	if(mus_pat_caloIso_branch == 0 ) {
	cout << "Branch mus_pat_caloIso does not exist." << endl;
	}
	mus_pat_ecalIso_branch = 0;
	if (tree->GetAlias("mus_pat_ecalIso") != 0) {
		mus_pat_ecalIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_ecalIso"));
		mus_pat_ecalIso_branch->SetAddress(&mus_pat_ecalIso_);
	}
	if(mus_pat_ecalIso_branch == 0 ) {
	cout << "Branch mus_pat_ecalIso does not exist." << endl;
	}
	mus_pat_ecalvetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_ecalvetoDep") != 0) {
		mus_pat_ecalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_ecalvetoDep"));
		mus_pat_ecalvetoDep_branch->SetAddress(&mus_pat_ecalvetoDep_);
	}
	if(mus_pat_ecalvetoDep_branch == 0 ) {
	cout << "Branch mus_pat_ecalvetoDep does not exist." << endl;
	}
	mus_pat_hcalIso_branch = 0;
	if (tree->GetAlias("mus_pat_hcalIso") != 0) {
		mus_pat_hcalIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_hcalIso"));
		mus_pat_hcalIso_branch->SetAddress(&mus_pat_hcalIso_);
	}
	if(mus_pat_hcalIso_branch == 0 ) {
	cout << "Branch mus_pat_hcalIso does not exist." << endl;
	}
	mus_pat_hcalvetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_hcalvetoDep") != 0) {
		mus_pat_hcalvetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_hcalvetoDep"));
		mus_pat_hcalvetoDep_branch->SetAddress(&mus_pat_hcalvetoDep_);
	}
	if(mus_pat_hcalvetoDep_branch == 0 ) {
	cout << "Branch mus_pat_hcalvetoDep does not exist." << endl;
	}
	mus_pat_trackIso_branch = 0;
	if (tree->GetAlias("mus_pat_trackIso") != 0) {
		mus_pat_trackIso_branch = tree->GetBranch(tree->GetAlias("mus_pat_trackIso"));
		mus_pat_trackIso_branch->SetAddress(&mus_pat_trackIso_);
	}
	if(mus_pat_trackIso_branch == 0 ) {
	cout << "Branch mus_pat_trackIso does not exist." << endl;
	}
	mus_pat_vetoDep_branch = 0;
	if (tree->GetAlias("mus_pat_vetoDep") != 0) {
		mus_pat_vetoDep_branch = tree->GetBranch(tree->GetAlias("mus_pat_vetoDep"));
		mus_pat_vetoDep_branch->SetAddress(&mus_pat_vetoDep_);
	}
	if(mus_pat_vetoDep_branch == 0 ) {
	cout << "Branch mus_pat_vetoDep does not exist." << endl;
	}
	taus_pf_bremsRecoveryEOverPLead_branch = 0;
	if (tree->GetAlias("taus_pf_bremsRecoveryEOverPLead") != 0) {
		taus_pf_bremsRecoveryEOverPLead_branch = tree->GetBranch(tree->GetAlias("taus_pf_bremsRecoveryEOverPLead"));
		taus_pf_bremsRecoveryEOverPLead_branch->SetAddress(&taus_pf_bremsRecoveryEOverPLead_);
	}
	if(taus_pf_bremsRecoveryEOverPLead_branch == 0 ) {
	cout << "Branch taus_pf_bremsRecoveryEOverPLead does not exist." << endl;
	}
	taus_pf_caloComp_branch = 0;
	if (tree->GetAlias("taus_pf_caloComp") != 0) {
		taus_pf_caloComp_branch = tree->GetBranch(tree->GetAlias("taus_pf_caloComp"));
		taus_pf_caloComp_branch->SetAddress(&taus_pf_caloComp_);
	}
	if(taus_pf_caloComp_branch == 0 ) {
	cout << "Branch taus_pf_caloComp does not exist." << endl;
	}
	taus_pf_ecalStripSumEOverPLead_branch = 0;
	if (tree->GetAlias("taus_pf_ecalStripSumEOverPLead") != 0) {
		taus_pf_ecalStripSumEOverPLead_branch = tree->GetBranch(tree->GetAlias("taus_pf_ecalStripSumEOverPLead"));
		taus_pf_ecalStripSumEOverPLead_branch->SetAddress(&taus_pf_ecalStripSumEOverPLead_);
	}
	if(taus_pf_ecalStripSumEOverPLead_branch == 0 ) {
	cout << "Branch taus_pf_ecalStripSumEOverPLead does not exist." << endl;
	}
	taus_pf_emf_branch = 0;
	if (tree->GetAlias("taus_pf_emf") != 0) {
		taus_pf_emf_branch = tree->GetBranch(tree->GetAlias("taus_pf_emf"));
		taus_pf_emf_branch->SetAddress(&taus_pf_emf_);
	}
	if(taus_pf_emf_branch == 0 ) {
	cout << "Branch taus_pf_emf does not exist." << endl;
	}
	taus_pf_hcal3x3OverPLead_branch = 0;
	if (tree->GetAlias("taus_pf_hcal3x3OverPLead") != 0) {
		taus_pf_hcal3x3OverPLead_branch = tree->GetBranch(tree->GetAlias("taus_pf_hcal3x3OverPLead"));
		taus_pf_hcal3x3OverPLead_branch->SetAddress(&taus_pf_hcal3x3OverPLead_);
	}
	if(taus_pf_hcal3x3OverPLead_branch == 0 ) {
	cout << "Branch taus_pf_hcal3x3OverPLead does not exist." << endl;
	}
	taus_pf_hcalMaxOverPLead_branch = 0;
	if (tree->GetAlias("taus_pf_hcalMaxOverPLead") != 0) {
		taus_pf_hcalMaxOverPLead_branch = tree->GetBranch(tree->GetAlias("taus_pf_hcalMaxOverPLead"));
		taus_pf_hcalMaxOverPLead_branch->SetAddress(&taus_pf_hcalMaxOverPLead_);
	}
	if(taus_pf_hcalMaxOverPLead_branch == 0 ) {
	cout << "Branch taus_pf_hcalMaxOverPLead does not exist." << endl;
	}
	taus_pf_hcalTotOverPLead_branch = 0;
	if (tree->GetAlias("taus_pf_hcalTotOverPLead") != 0) {
		taus_pf_hcalTotOverPLead_branch = tree->GetBranch(tree->GetAlias("taus_pf_hcalTotOverPLead"));
		taus_pf_hcalTotOverPLead_branch->SetAddress(&taus_pf_hcalTotOverPLead_);
	}
	if(taus_pf_hcalTotOverPLead_branch == 0 ) {
	cout << "Branch taus_pf_hcalTotOverPLead does not exist." << endl;
	}
	taus_pf_isolationchargecandPtSum_branch = 0;
	if (tree->GetAlias("taus_pf_isolationchargecandPtSum") != 0) {
		taus_pf_isolationchargecandPtSum_branch = tree->GetBranch(tree->GetAlias("taus_pf_isolationchargecandPtSum"));
		taus_pf_isolationchargecandPtSum_branch->SetAddress(&taus_pf_isolationchargecandPtSum_);
	}
	if(taus_pf_isolationchargecandPtSum_branch == 0 ) {
	cout << "Branch taus_pf_isolationchargecandPtSum does not exist." << endl;
	}
	taus_pf_isolationgammacandEtSum_branch = 0;
	if (tree->GetAlias("taus_pf_isolationgammacandEtSum") != 0) {
		taus_pf_isolationgammacandEtSum_branch = tree->GetBranch(tree->GetAlias("taus_pf_isolationgammacandEtSum"));
		taus_pf_isolationgammacandEtSum_branch->SetAddress(&taus_pf_isolationgammacandEtSum_);
	}
	if(taus_pf_isolationgammacandEtSum_branch == 0 ) {
	cout << "Branch taus_pf_isolationgammacandEtSum does not exist." << endl;
	}
	taus_pf_lead_chargecand_Signed_Sipt_branch = 0;
	if (tree->GetAlias("taus_pf_lead_chargecand_Signed_Sipt") != 0) {
		taus_pf_lead_chargecand_Signed_Sipt_branch = tree->GetBranch(tree->GetAlias("taus_pf_lead_chargecand_Signed_Sipt"));
		taus_pf_lead_chargecand_Signed_Sipt_branch->SetAddress(&taus_pf_lead_chargecand_Signed_Sipt_);
	}
	if(taus_pf_lead_chargecand_Signed_Sipt_branch == 0 ) {
	cout << "Branch taus_pf_lead_chargecand_Signed_Sipt does not exist." << endl;
	}
	taus_pf_leadtrk_chi2_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_chi2") != 0) {
		taus_pf_leadtrk_chi2_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_chi2"));
		taus_pf_leadtrk_chi2_branch->SetAddress(&taus_pf_leadtrk_chi2_);
	}
	if(taus_pf_leadtrk_chi2_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_chi2 does not exist." << endl;
	}
	taus_pf_leadtrk_d0_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_d0") != 0) {
		taus_pf_leadtrk_d0_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_d0"));
		taus_pf_leadtrk_d0_branch->SetAddress(&taus_pf_leadtrk_d0_);
	}
	if(taus_pf_leadtrk_d0_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_d0 does not exist." << endl;
	}
	taus_pf_leadtrk_lostHits_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_lostHits") != 0) {
		taus_pf_leadtrk_lostHits_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_lostHits"));
		taus_pf_leadtrk_lostHits_branch->SetAddress(&taus_pf_leadtrk_lostHits_);
	}
	if(taus_pf_leadtrk_lostHits_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_lostHits does not exist." << endl;
	}
	taus_pf_leadtrk_ndof_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_ndof") != 0) {
		taus_pf_leadtrk_ndof_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_ndof"));
		taus_pf_leadtrk_ndof_branch->SetAddress(&taus_pf_leadtrk_ndof_);
	}
	if(taus_pf_leadtrk_ndof_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_ndof does not exist." << endl;
	}
	taus_pf_leadtrk_validHits_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_validHits") != 0) {
		taus_pf_leadtrk_validHits_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_validHits"));
		taus_pf_leadtrk_validHits_branch->SetAddress(&taus_pf_leadtrk_validHits_);
	}
	if(taus_pf_leadtrk_validHits_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_validHits does not exist." << endl;
	}
	taus_pf_leadtrk_z0_branch = 0;
	if (tree->GetAlias("taus_pf_leadtrk_z0") != 0) {
		taus_pf_leadtrk_z0_branch = tree->GetBranch(tree->GetAlias("taus_pf_leadtrk_z0"));
		taus_pf_leadtrk_z0_branch->SetAddress(&taus_pf_leadtrk_z0_);
	}
	if(taus_pf_leadtrk_z0_branch == 0 ) {
	cout << "Branch taus_pf_leadtrk_z0 does not exist." << endl;
	}
	taus_pf_maximumHCALPFClusterEt_branch = 0;
	if (tree->GetAlias("taus_pf_maximumHCALPFClusterEt") != 0) {
		taus_pf_maximumHCALPFClusterEt_branch = tree->GetBranch(tree->GetAlias("taus_pf_maximumHCALPFClusterEt"));
		taus_pf_maximumHCALPFClusterEt_branch->SetAddress(&taus_pf_maximumHCALPFClusterEt_);
	}
	if(taus_pf_maximumHCALPFClusterEt_branch == 0 ) {
	cout << "Branch taus_pf_maximumHCALPFClusterEt does not exist." << endl;
	}
	taus_pf_segComp_branch = 0;
	if (tree->GetAlias("taus_pf_segComp") != 0) {
		taus_pf_segComp_branch = tree->GetBranch(tree->GetAlias("taus_pf_segComp"));
		taus_pf_segComp_branch->SetAddress(&taus_pf_segComp_);
	}
	if(taus_pf_segComp_branch == 0 ) {
	cout << "Branch taus_pf_segComp does not exist." << endl;
	}
	photons_e1x5_branch = 0;
	if (tree->GetAlias("photons_e1x5") != 0) {
		photons_e1x5_branch = tree->GetBranch(tree->GetAlias("photons_e1x5"));
		photons_e1x5_branch->SetAddress(&photons_e1x5_);
	}
	if(photons_e1x5_branch == 0 ) {
	cout << "Branch photons_e1x5 does not exist." << endl;
	}
	photons_e2x5Max_branch = 0;
	if (tree->GetAlias("photons_e2x5Max") != 0) {
		photons_e2x5Max_branch = tree->GetBranch(tree->GetAlias("photons_e2x5Max"));
		photons_e2x5Max_branch->SetAddress(&photons_e2x5Max_);
	}
	if(photons_e2x5Max_branch == 0 ) {
	cout << "Branch photons_e2x5Max does not exist." << endl;
	}
	photons_e3x3_branch = 0;
	if (tree->GetAlias("photons_e3x3") != 0) {
		photons_e3x3_branch = tree->GetBranch(tree->GetAlias("photons_e3x3"));
		photons_e3x3_branch->SetAddress(&photons_e3x3_);
	}
	if(photons_e3x3_branch == 0 ) {
	cout << "Branch photons_e3x3 does not exist." << endl;
	}
	photons_e5x5_branch = 0;
	if (tree->GetAlias("photons_e5x5") != 0) {
		photons_e5x5_branch = tree->GetBranch(tree->GetAlias("photons_e5x5"));
		photons_e5x5_branch->SetAddress(&photons_e5x5_);
	}
	if(photons_e5x5_branch == 0 ) {
	cout << "Branch photons_e5x5 does not exist." << endl;
	}
	photons_eMax_branch = 0;
	if (tree->GetAlias("photons_eMax") != 0) {
		photons_eMax_branch = tree->GetBranch(tree->GetAlias("photons_eMax"));
		photons_eMax_branch->SetAddress(&photons_eMax_);
	}
	if(photons_eMax_branch == 0 ) {
	cout << "Branch photons_eMax does not exist." << endl;
	}
	photons_eSC_branch = 0;
	if (tree->GetAlias("photons_eSC") != 0) {
		photons_eSC_branch = tree->GetBranch(tree->GetAlias("photons_eSC"));
		photons_eSC_branch->SetAddress(&photons_eSC_);
	}
	if(photons_eSC_branch == 0 ) {
	cout << "Branch photons_eSC does not exist." << endl;
	}
	photons_eSCPresh_branch = 0;
	if (tree->GetAlias("photons_eSCPresh") != 0) {
		photons_eSCPresh_branch = tree->GetBranch(tree->GetAlias("photons_eSCPresh"));
		photons_eSCPresh_branch->SetAddress(&photons_eSCPresh_);
	}
	if(photons_eSCPresh_branch == 0 ) {
	cout << "Branch photons_eSCPresh does not exist." << endl;
	}
	photons_eSCRaw_branch = 0;
	if (tree->GetAlias("photons_eSCRaw") != 0) {
		photons_eSCRaw_branch = tree->GetBranch(tree->GetAlias("photons_eSCRaw"));
		photons_eSCRaw_branch->SetAddress(&photons_eSCRaw_);
	}
	if(photons_eSCRaw_branch == 0 ) {
	cout << "Branch photons_eSCRaw does not exist." << endl;
	}
	photons_eSeed_branch = 0;
	if (tree->GetAlias("photons_eSeed") != 0) {
		photons_eSeed_branch = tree->GetBranch(tree->GetAlias("photons_eSeed"));
		photons_eSeed_branch->SetAddress(&photons_eSeed_);
	}
	if(photons_eSeed_branch == 0 ) {
	cout << "Branch photons_eSeed does not exist." << endl;
	}
	photons_ecalIso_branch = 0;
	if (tree->GetAlias("photons_ecalIso") != 0) {
		photons_ecalIso_branch = tree->GetBranch(tree->GetAlias("photons_ecalIso"));
		photons_ecalIso_branch->SetAddress(&photons_ecalIso_);
	}
	if(photons_ecalIso_branch == 0 ) {
	cout << "Branch photons_ecalIso does not exist." << endl;
	}
	photons_hOverE_branch = 0;
	if (tree->GetAlias("photons_hOverE") != 0) {
		photons_hOverE_branch = tree->GetBranch(tree->GetAlias("photons_hOverE"));
		photons_hOverE_branch->SetAddress(&photons_hOverE_);
	}
	if(photons_hOverE_branch == 0 ) {
	cout << "Branch photons_hOverE does not exist." << endl;
	}
	photons_hcalIso_branch = 0;
	if (tree->GetAlias("photons_hcalIso") != 0) {
		photons_hcalIso_branch = tree->GetBranch(tree->GetAlias("photons_hcalIso"));
		photons_hcalIso_branch->SetAddress(&photons_hcalIso_);
	}
	if(photons_hcalIso_branch == 0 ) {
	cout << "Branch photons_hcalIso does not exist." << endl;
	}
	photons_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("photons_sigmaEtaEta") != 0) {
		photons_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaEtaEta"));
		photons_sigmaEtaEta_branch->SetAddress(&photons_sigmaEtaEta_);
	}
	if(photons_sigmaEtaEta_branch == 0 ) {
	cout << "Branch photons_sigmaEtaEta does not exist." << endl;
	}
	photons_sigmaIEtaIEta_branch = 0;
	if (tree->GetAlias("photons_sigmaIEtaIEta") != 0) {
		photons_sigmaIEtaIEta_branch = tree->GetBranch(tree->GetAlias("photons_sigmaIEtaIEta"));
		photons_sigmaIEtaIEta_branch->SetAddress(&photons_sigmaIEtaIEta_);
	}
	if(photons_sigmaIEtaIEta_branch == 0 ) {
	cout << "Branch photons_sigmaIEtaIEta does not exist." << endl;
	}
	photons_sigmaIPhiIPhi_branch = 0;
	if (tree->GetAlias("photons_sigmaIPhiIPhi") != 0) {
		photons_sigmaIPhiIPhi_branch = tree->GetBranch(tree->GetAlias("photons_sigmaIPhiIPhi"));
		photons_sigmaIPhiIPhi_branch->SetAddress(&photons_sigmaIPhiIPhi_);
	}
	if(photons_sigmaIPhiIPhi_branch == 0 ) {
	cout << "Branch photons_sigmaIPhiIPhi does not exist." << endl;
	}
	photons_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("photons_sigmaPhiPhi") != 0) {
		photons_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("photons_sigmaPhiPhi"));
		photons_sigmaPhiPhi_branch->SetAddress(&photons_sigmaPhiPhi_);
	}
	if(photons_sigmaPhiPhi_branch == 0 ) {
	cout << "Branch photons_sigmaPhiPhi does not exist." << endl;
	}
	photons_tkIso_branch = 0;
	if (tree->GetAlias("photons_tkIso") != 0) {
		photons_tkIso_branch = tree->GetBranch(tree->GetAlias("photons_tkIso"));
		photons_tkIso_branch->SetAddress(&photons_tkIso_);
	}
	if(photons_tkIso_branch == 0 ) {
	cout << "Branch photons_tkIso does not exist." << endl;
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
	mus_tcmet_deltax_branch = 0;
	if (tree->GetAlias("mus_tcmet_deltax") != 0) {
		mus_tcmet_deltax_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_deltax"));
		mus_tcmet_deltax_branch->SetAddress(&mus_tcmet_deltax_);
	}
	if(mus_tcmet_deltax_branch == 0 ) {
	cout << "Branch mus_tcmet_deltax does not exist." << endl;
	}
	mus_tcmet_deltay_branch = 0;
	if (tree->GetAlias("mus_tcmet_deltay") != 0) {
		mus_tcmet_deltay_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_deltay"));
		mus_tcmet_deltay_branch->SetAddress(&mus_tcmet_deltay_);
	}
	if(mus_tcmet_deltay_branch == 0 ) {
	cout << "Branch mus_tcmet_deltay does not exist." << endl;
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
	evt_nl1emiso_branch = 0;
	if (tree->GetAlias("evt_nl1emiso") != 0) {
		evt_nl1emiso_branch = tree->GetBranch(tree->GetAlias("evt_nl1emiso"));
		evt_nl1emiso_branch->SetAddress(&evt_nl1emiso_);
	}
	if(evt_nl1emiso_branch == 0 ) {
	cout << "Branch evt_nl1emiso does not exist." << endl;
	}
	evt_nl1emnoiso_branch = 0;
	if (tree->GetAlias("evt_nl1emnoiso") != 0) {
		evt_nl1emnoiso_branch = tree->GetBranch(tree->GetAlias("evt_nl1emnoiso"));
		evt_nl1emnoiso_branch->SetAddress(&evt_nl1emnoiso_);
	}
	if(evt_nl1emnoiso_branch == 0 ) {
	cout << "Branch evt_nl1emnoiso does not exist." << endl;
	}
	evt_nl1jetsc_branch = 0;
	if (tree->GetAlias("evt_nl1jetsc") != 0) {
		evt_nl1jetsc_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetsc"));
		evt_nl1jetsc_branch->SetAddress(&evt_nl1jetsc_);
	}
	if(evt_nl1jetsc_branch == 0 ) {
	cout << "Branch evt_nl1jetsc does not exist." << endl;
	}
	evt_nl1jetsf_branch = 0;
	if (tree->GetAlias("evt_nl1jetsf") != 0) {
		evt_nl1jetsf_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetsf"));
		evt_nl1jetsf_branch->SetAddress(&evt_nl1jetsf_);
	}
	if(evt_nl1jetsf_branch == 0 ) {
	cout << "Branch evt_nl1jetsf does not exist." << endl;
	}
	evt_nl1jetst_branch = 0;
	if (tree->GetAlias("evt_nl1jetst") != 0) {
		evt_nl1jetst_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetst"));
		evt_nl1jetst_branch->SetAddress(&evt_nl1jetst_);
	}
	if(evt_nl1jetst_branch == 0 ) {
	cout << "Branch evt_nl1jetst does not exist." << endl;
	}
	evt_nl1mus_branch = 0;
	if (tree->GetAlias("evt_nl1mus") != 0) {
		evt_nl1mus_branch = tree->GetBranch(tree->GetAlias("evt_nl1mus"));
		evt_nl1mus_branch->SetAddress(&evt_nl1mus_);
	}
	if(evt_nl1mus_branch == 0 ) {
	cout << "Branch evt_nl1mus does not exist." << endl;
	}
	pdfinfo_id1_branch = 0;
	if (tree->GetAlias("pdfinfo_id1") != 0) {
		pdfinfo_id1_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id1"));
		pdfinfo_id1_branch->SetAddress(&pdfinfo_id1_);
	}
	if(pdfinfo_id1_branch == 0 ) {
	cout << "Branch pdfinfo_id1 does not exist." << endl;
	}
	pdfinfo_id2_branch = 0;
	if (tree->GetAlias("pdfinfo_id2") != 0) {
		pdfinfo_id2_branch = tree->GetBranch(tree->GetAlias("pdfinfo_id2"));
		pdfinfo_id2_branch->SetAddress(&pdfinfo_id2_);
	}
	if(pdfinfo_id2_branch == 0 ) {
	cout << "Branch pdfinfo_id2 does not exist." << endl;
	}
	taus_calo_charge_branch = 0;
	if (tree->GetAlias("taus_calo_charge") != 0) {
		taus_calo_charge_branch = tree->GetBranch(tree->GetAlias("taus_calo_charge"));
		taus_calo_charge_branch->SetAddress(&taus_calo_charge_);
	}
	if(taus_calo_charge_branch == 0 ) {
	cout << "Branch taus_calo_charge does not exist." << endl;
	}
	taus_calo_iso_ntrks_branch = 0;
	if (tree->GetAlias("taus_calo_iso_ntrks") != 0) {
		taus_calo_iso_ntrks_branch = tree->GetBranch(tree->GetAlias("taus_calo_iso_ntrks"));
		taus_calo_iso_ntrks_branch->SetAddress(&taus_calo_iso_ntrks_);
	}
	if(taus_calo_iso_ntrks_branch == 0 ) {
	cout << "Branch taus_calo_iso_ntrks does not exist." << endl;
	}
	taus_calo_sig_ntrks_branch = 0;
	if (tree->GetAlias("taus_calo_sig_ntrks") != 0) {
		taus_calo_sig_ntrks_branch = tree->GetBranch(tree->GetAlias("taus_calo_sig_ntrks"));
		taus_calo_sig_ntrks_branch->SetAddress(&taus_calo_sig_ntrks_);
	}
	if(taus_calo_sig_ntrks_branch == 0 ) {
	cout << "Branch taus_calo_sig_ntrks does not exist." << endl;
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
	els_closestJet_branch = 0;
	if (tree->GetAlias("els_closestJet") != 0) {
		els_closestJet_branch = tree->GetBranch(tree->GetAlias("els_closestJet"));
		els_closestJet_branch->SetAddress(&els_closestJet_);
	}
	if(els_closestJet_branch == 0 ) {
	cout << "Branch els_closestJet does not exist." << endl;
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
	hyp_ll_charge_branch = 0;
	if (tree->GetAlias("hyp_ll_charge") != 0) {
		hyp_ll_charge_branch = tree->GetBranch(tree->GetAlias("hyp_ll_charge"));
		hyp_ll_charge_branch->SetAddress(&hyp_ll_charge_);
	}
	if(hyp_ll_charge_branch == 0 ) {
	cout << "Branch hyp_ll_charge does not exist." << endl;
	}
	hyp_ll_id_branch = 0;
	if (tree->GetAlias("hyp_ll_id") != 0) {
		hyp_ll_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_id"));
		hyp_ll_id_branch->SetAddress(&hyp_ll_id_);
	}
	if(hyp_ll_id_branch == 0 ) {
	cout << "Branch hyp_ll_id does not exist." << endl;
	}
	hyp_ll_index_branch = 0;
	if (tree->GetAlias("hyp_ll_index") != 0) {
		hyp_ll_index_branch = tree->GetBranch(tree->GetAlias("hyp_ll_index"));
		hyp_ll_index_branch->SetAddress(&hyp_ll_index_);
	}
	if(hyp_ll_index_branch == 0 ) {
	cout << "Branch hyp_ll_index does not exist." << endl;
	}
	hyp_ll_lostHits_branch = 0;
	if (tree->GetAlias("hyp_ll_lostHits") != 0) {
		hyp_ll_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_lostHits"));
		hyp_ll_lostHits_branch->SetAddress(&hyp_ll_lostHits_);
	}
	if(hyp_ll_lostHits_branch == 0 ) {
	cout << "Branch hyp_ll_lostHits does not exist." << endl;
	}
	hyp_ll_mc_id_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_id") != 0) {
		hyp_ll_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_id"));
		hyp_ll_mc_id_branch->SetAddress(&hyp_ll_mc_id_);
	}
	if(hyp_ll_mc_id_branch == 0 ) {
	cout << "Branch hyp_ll_mc_id does not exist." << endl;
	}
	hyp_ll_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_motherid") != 0) {
		hyp_ll_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_motherid"));
		hyp_ll_mc_motherid_branch->SetAddress(&hyp_ll_mc_motherid_);
	}
	if(hyp_ll_mc_motherid_branch == 0 ) {
	cout << "Branch hyp_ll_mc_motherid does not exist." << endl;
	}
	hyp_ll_validHits_branch = 0;
	if (tree->GetAlias("hyp_ll_validHits") != 0) {
		hyp_ll_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_validHits"));
		hyp_ll_validHits_branch->SetAddress(&hyp_ll_validHits_);
	}
	if(hyp_ll_validHits_branch == 0 ) {
	cout << "Branch hyp_ll_validHits does not exist." << endl;
	}
	hyp_lt_charge_branch = 0;
	if (tree->GetAlias("hyp_lt_charge") != 0) {
		hyp_lt_charge_branch = tree->GetBranch(tree->GetAlias("hyp_lt_charge"));
		hyp_lt_charge_branch->SetAddress(&hyp_lt_charge_);
	}
	if(hyp_lt_charge_branch == 0 ) {
	cout << "Branch hyp_lt_charge does not exist." << endl;
	}
	hyp_lt_id_branch = 0;
	if (tree->GetAlias("hyp_lt_id") != 0) {
		hyp_lt_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_id"));
		hyp_lt_id_branch->SetAddress(&hyp_lt_id_);
	}
	if(hyp_lt_id_branch == 0 ) {
	cout << "Branch hyp_lt_id does not exist." << endl;
	}
	hyp_lt_index_branch = 0;
	if (tree->GetAlias("hyp_lt_index") != 0) {
		hyp_lt_index_branch = tree->GetBranch(tree->GetAlias("hyp_lt_index"));
		hyp_lt_index_branch->SetAddress(&hyp_lt_index_);
	}
	if(hyp_lt_index_branch == 0 ) {
	cout << "Branch hyp_lt_index does not exist." << endl;
	}
	hyp_lt_lostHits_branch = 0;
	if (tree->GetAlias("hyp_lt_lostHits") != 0) {
		hyp_lt_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_lostHits"));
		hyp_lt_lostHits_branch->SetAddress(&hyp_lt_lostHits_);
	}
	if(hyp_lt_lostHits_branch == 0 ) {
	cout << "Branch hyp_lt_lostHits does not exist." << endl;
	}
	hyp_lt_mc_id_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_id") != 0) {
		hyp_lt_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_id"));
		hyp_lt_mc_id_branch->SetAddress(&hyp_lt_mc_id_);
	}
	if(hyp_lt_mc_id_branch == 0 ) {
	cout << "Branch hyp_lt_mc_id does not exist." << endl;
	}
	hyp_lt_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_motherid") != 0) {
		hyp_lt_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_motherid"));
		hyp_lt_mc_motherid_branch->SetAddress(&hyp_lt_mc_motherid_);
	}
	if(hyp_lt_mc_motherid_branch == 0 ) {
	cout << "Branch hyp_lt_mc_motherid does not exist." << endl;
	}
	hyp_lt_validHits_branch = 0;
	if (tree->GetAlias("hyp_lt_validHits") != 0) {
		hyp_lt_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_validHits"));
		hyp_lt_validHits_branch->SetAddress(&hyp_lt_validHits_);
	}
	if(hyp_lt_validHits_branch == 0 ) {
	cout << "Branch hyp_lt_validHits does not exist." << endl;
	}
	hyp_njets_branch = 0;
	if (tree->GetAlias("hyp_njets") != 0) {
		hyp_njets_branch = tree->GetBranch(tree->GetAlias("hyp_njets"));
		hyp_njets_branch->SetAddress(&hyp_njets_);
	}
	if(hyp_njets_branch == 0 ) {
	cout << "Branch hyp_njets does not exist." << endl;
	}
	hyp_nojets_branch = 0;
	if (tree->GetAlias("hyp_nojets") != 0) {
		hyp_nojets_branch = tree->GetBranch(tree->GetAlias("hyp_nojets"));
		hyp_nojets_branch->SetAddress(&hyp_nojets_);
	}
	if(hyp_nojets_branch == 0 ) {
	cout << "Branch hyp_nojets does not exist." << endl;
	}
	hyp_type_branch = 0;
	if (tree->GetAlias("hyp_type") != 0) {
		hyp_type_branch = tree->GetBranch(tree->GetAlias("hyp_type"));
		hyp_type_branch->SetAddress(&hyp_type_);
	}
	if(hyp_type_branch == 0 ) {
	cout << "Branch hyp_type does not exist." << endl;
	}
	hyp_quadlep_first_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_type") != 0) {
		hyp_quadlep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_type"));
		hyp_quadlep_first_type_branch->SetAddress(&hyp_quadlep_first_type_);
	}
	if(hyp_quadlep_first_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_first_type does not exist." << endl;
	}
	hyp_quadlep_fourth_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_type") != 0) {
		hyp_quadlep_fourth_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_type"));
		hyp_quadlep_fourth_type_branch->SetAddress(&hyp_quadlep_fourth_type_);
	}
	if(hyp_quadlep_fourth_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_fourth_type does not exist." << endl;
	}
	hyp_quadlep_second_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_type") != 0) {
		hyp_quadlep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_type"));
		hyp_quadlep_second_type_branch->SetAddress(&hyp_quadlep_second_type_);
	}
	if(hyp_quadlep_second_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_second_type does not exist." << endl;
	}
	hyp_quadlep_third_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_type") != 0) {
		hyp_quadlep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_type"));
		hyp_quadlep_third_type_branch->SetAddress(&hyp_quadlep_third_type_);
	}
	if(hyp_quadlep_third_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_third_type does not exist." << endl;
	}
	hyp_trilep_first_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_type") != 0) {
		hyp_trilep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_type"));
		hyp_trilep_first_type_branch->SetAddress(&hyp_trilep_first_type_);
	}
	if(hyp_trilep_first_type_branch == 0 ) {
	cout << "Branch hyp_trilep_first_type does not exist." << endl;
	}
	hyp_trilep_second_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_type") != 0) {
		hyp_trilep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_type"));
		hyp_trilep_second_type_branch->SetAddress(&hyp_trilep_second_type_);
	}
	if(hyp_trilep_second_type_branch == 0 ) {
	cout << "Branch hyp_trilep_second_type does not exist." << endl;
	}
	hyp_trilep_third_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_type") != 0) {
		hyp_trilep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_type"));
		hyp_trilep_third_type_branch->SetAddress(&hyp_trilep_third_type_);
	}
	if(hyp_trilep_third_type_branch == 0 ) {
	cout << "Branch hyp_trilep_third_type does not exist." << endl;
	}
	jets_closestElectron_branch = 0;
	if (tree->GetAlias("jets_closestElectron") != 0) {
		jets_closestElectron_branch = tree->GetBranch(tree->GetAlias("jets_closestElectron"));
		jets_closestElectron_branch->SetAddress(&jets_closestElectron_);
	}
	if(jets_closestElectron_branch == 0 ) {
	cout << "Branch jets_closestElectron does not exist." << endl;
	}
	jets_closestMuon_branch = 0;
	if (tree->GetAlias("jets_closestMuon") != 0) {
		jets_closestMuon_branch = tree->GetBranch(tree->GetAlias("jets_closestMuon"));
		jets_closestMuon_branch->SetAddress(&jets_closestMuon_);
	}
	if(jets_closestMuon_branch == 0 ) {
	cout << "Branch jets_closestMuon does not exist." << endl;
	}
	l1emiso_ieta_branch = 0;
	if (tree->GetAlias("l1emiso_ieta") != 0) {
		l1emiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1emiso_ieta"));
		l1emiso_ieta_branch->SetAddress(&l1emiso_ieta_);
	}
	if(l1emiso_ieta_branch == 0 ) {
	cout << "Branch l1emiso_ieta does not exist." << endl;
	}
	l1emiso_iphi_branch = 0;
	if (tree->GetAlias("l1emiso_iphi") != 0) {
		l1emiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1emiso_iphi"));
		l1emiso_iphi_branch->SetAddress(&l1emiso_iphi_);
	}
	if(l1emiso_iphi_branch == 0 ) {
	cout << "Branch l1emiso_iphi does not exist." << endl;
	}
	l1emiso_rawId_branch = 0;
	if (tree->GetAlias("l1emiso_rawId") != 0) {
		l1emiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1emiso_rawId"));
		l1emiso_rawId_branch->SetAddress(&l1emiso_rawId_);
	}
	if(l1emiso_rawId_branch == 0 ) {
	cout << "Branch l1emiso_rawId does not exist." << endl;
	}
	l1emiso_type_branch = 0;
	if (tree->GetAlias("l1emiso_type") != 0) {
		l1emiso_type_branch = tree->GetBranch(tree->GetAlias("l1emiso_type"));
		l1emiso_type_branch->SetAddress(&l1emiso_type_);
	}
	if(l1emiso_type_branch == 0 ) {
	cout << "Branch l1emiso_type does not exist." << endl;
	}
	l1emnoiso_ieta_branch = 0;
	if (tree->GetAlias("l1emnoiso_ieta") != 0) {
		l1emnoiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_ieta"));
		l1emnoiso_ieta_branch->SetAddress(&l1emnoiso_ieta_);
	}
	if(l1emnoiso_ieta_branch == 0 ) {
	cout << "Branch l1emnoiso_ieta does not exist." << endl;
	}
	l1emnoiso_iphi_branch = 0;
	if (tree->GetAlias("l1emnoiso_iphi") != 0) {
		l1emnoiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_iphi"));
		l1emnoiso_iphi_branch->SetAddress(&l1emnoiso_iphi_);
	}
	if(l1emnoiso_iphi_branch == 0 ) {
	cout << "Branch l1emnoiso_iphi does not exist." << endl;
	}
	l1emnoiso_rawId_branch = 0;
	if (tree->GetAlias("l1emnoiso_rawId") != 0) {
		l1emnoiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_rawId"));
		l1emnoiso_rawId_branch->SetAddress(&l1emnoiso_rawId_);
	}
	if(l1emnoiso_rawId_branch == 0 ) {
	cout << "Branch l1emnoiso_rawId does not exist." << endl;
	}
	l1emnoiso_type_branch = 0;
	if (tree->GetAlias("l1emnoiso_type") != 0) {
		l1emnoiso_type_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_type"));
		l1emnoiso_type_branch->SetAddress(&l1emnoiso_type_);
	}
	if(l1emnoiso_type_branch == 0 ) {
	cout << "Branch l1emnoiso_type does not exist." << endl;
	}
	l1jetsc_ieta_branch = 0;
	if (tree->GetAlias("l1jetsc_ieta") != 0) {
		l1jetsc_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetsc_ieta"));
		l1jetsc_ieta_branch->SetAddress(&l1jetsc_ieta_);
	}
	if(l1jetsc_ieta_branch == 0 ) {
	cout << "Branch l1jetsc_ieta does not exist." << endl;
	}
	l1jetsc_iphi_branch = 0;
	if (tree->GetAlias("l1jetsc_iphi") != 0) {
		l1jetsc_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetsc_iphi"));
		l1jetsc_iphi_branch->SetAddress(&l1jetsc_iphi_);
	}
	if(l1jetsc_iphi_branch == 0 ) {
	cout << "Branch l1jetsc_iphi does not exist." << endl;
	}
	l1jetsc_rawId_branch = 0;
	if (tree->GetAlias("l1jetsc_rawId") != 0) {
		l1jetsc_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetsc_rawId"));
		l1jetsc_rawId_branch->SetAddress(&l1jetsc_rawId_);
	}
	if(l1jetsc_rawId_branch == 0 ) {
	cout << "Branch l1jetsc_rawId does not exist." << endl;
	}
	l1jetsc_type_branch = 0;
	if (tree->GetAlias("l1jetsc_type") != 0) {
		l1jetsc_type_branch = tree->GetBranch(tree->GetAlias("l1jetsc_type"));
		l1jetsc_type_branch->SetAddress(&l1jetsc_type_);
	}
	if(l1jetsc_type_branch == 0 ) {
	cout << "Branch l1jetsc_type does not exist." << endl;
	}
	l1jetsf_ieta_branch = 0;
	if (tree->GetAlias("l1jetsf_ieta") != 0) {
		l1jetsf_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetsf_ieta"));
		l1jetsf_ieta_branch->SetAddress(&l1jetsf_ieta_);
	}
	if(l1jetsf_ieta_branch == 0 ) {
	cout << "Branch l1jetsf_ieta does not exist." << endl;
	}
	l1jetsf_iphi_branch = 0;
	if (tree->GetAlias("l1jetsf_iphi") != 0) {
		l1jetsf_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetsf_iphi"));
		l1jetsf_iphi_branch->SetAddress(&l1jetsf_iphi_);
	}
	if(l1jetsf_iphi_branch == 0 ) {
	cout << "Branch l1jetsf_iphi does not exist." << endl;
	}
	l1jetsf_rawId_branch = 0;
	if (tree->GetAlias("l1jetsf_rawId") != 0) {
		l1jetsf_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetsf_rawId"));
		l1jetsf_rawId_branch->SetAddress(&l1jetsf_rawId_);
	}
	if(l1jetsf_rawId_branch == 0 ) {
	cout << "Branch l1jetsf_rawId does not exist." << endl;
	}
	l1jetsf_type_branch = 0;
	if (tree->GetAlias("l1jetsf_type") != 0) {
		l1jetsf_type_branch = tree->GetBranch(tree->GetAlias("l1jetsf_type"));
		l1jetsf_type_branch->SetAddress(&l1jetsf_type_);
	}
	if(l1jetsf_type_branch == 0 ) {
	cout << "Branch l1jetsf_type does not exist." << endl;
	}
	l1jetst_ieta_branch = 0;
	if (tree->GetAlias("l1jetst_ieta") != 0) {
		l1jetst_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetst_ieta"));
		l1jetst_ieta_branch->SetAddress(&l1jetst_ieta_);
	}
	if(l1jetst_ieta_branch == 0 ) {
	cout << "Branch l1jetst_ieta does not exist." << endl;
	}
	l1jetst_iphi_branch = 0;
	if (tree->GetAlias("l1jetst_iphi") != 0) {
		l1jetst_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetst_iphi"));
		l1jetst_iphi_branch->SetAddress(&l1jetst_iphi_);
	}
	if(l1jetst_iphi_branch == 0 ) {
	cout << "Branch l1jetst_iphi does not exist." << endl;
	}
	l1jetst_rawId_branch = 0;
	if (tree->GetAlias("l1jetst_rawId") != 0) {
		l1jetst_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetst_rawId"));
		l1jetst_rawId_branch->SetAddress(&l1jetst_rawId_);
	}
	if(l1jetst_rawId_branch == 0 ) {
	cout << "Branch l1jetst_rawId does not exist." << endl;
	}
	l1jetst_type_branch = 0;
	if (tree->GetAlias("l1jetst_type") != 0) {
		l1jetst_type_branch = tree->GetBranch(tree->GetAlias("l1jetst_type"));
		l1jetst_type_branch->SetAddress(&l1jetst_type_);
	}
	if(l1jetst_type_branch == 0 ) {
	cout << "Branch l1jetst_type does not exist." << endl;
	}
	l1mus_flags_branch = 0;
	if (tree->GetAlias("l1mus_flags") != 0) {
		l1mus_flags_branch = tree->GetBranch(tree->GetAlias("l1mus_flags"));
		l1mus_flags_branch->SetAddress(&l1mus_flags_);
	}
	if(l1mus_flags_branch == 0 ) {
	cout << "Branch l1mus_flags does not exist." << endl;
	}
	l1mus_q_branch = 0;
	if (tree->GetAlias("l1mus_q") != 0) {
		l1mus_q_branch = tree->GetBranch(tree->GetAlias("l1mus_q"));
		l1mus_q_branch->SetAddress(&l1mus_q_);
	}
	if(l1mus_q_branch == 0 ) {
	cout << "Branch l1mus_q does not exist." << endl;
	}
	l1mus_qual_branch = 0;
	if (tree->GetAlias("l1mus_qual") != 0) {
		l1mus_qual_branch = tree->GetBranch(tree->GetAlias("l1mus_qual"));
		l1mus_qual_branch->SetAddress(&l1mus_qual_);
	}
	if(l1mus_qual_branch == 0 ) {
	cout << "Branch l1mus_qual does not exist." << endl;
	}
	l1mus_qualFlags_branch = 0;
	if (tree->GetAlias("l1mus_qualFlags") != 0) {
		l1mus_qualFlags_branch = tree->GetBranch(tree->GetAlias("l1mus_qualFlags"));
		l1mus_qualFlags_branch->SetAddress(&l1mus_qualFlags_);
	}
	if(l1mus_qualFlags_branch == 0 ) {
	cout << "Branch l1mus_qualFlags does not exist." << endl;
	}
	mus_met_flag_branch = 0;
	if (tree->GetAlias("mus_met_flag") != 0) {
		mus_met_flag_branch = tree->GetBranch(tree->GetAlias("mus_met_flag"));
		mus_met_flag_branch->SetAddress(&mus_met_flag_);
	}
	if(mus_met_flag_branch == 0 ) {
	cout << "Branch mus_met_flag does not exist." << endl;
	}
	mus_closestEle_branch = 0;
	if (tree->GetAlias("mus_closestEle") != 0) {
		mus_closestEle_branch = tree->GetBranch(tree->GetAlias("mus_closestEle"));
		mus_closestEle_branch->SetAddress(&mus_closestEle_);
	}
	if(mus_closestEle_branch == 0 ) {
	cout << "Branch mus_closestEle does not exist." << endl;
	}
	mus_closestJet_branch = 0;
	if (tree->GetAlias("mus_closestJet") != 0) {
		mus_closestJet_branch = tree->GetBranch(tree->GetAlias("mus_closestJet"));
		mus_closestJet_branch->SetAddress(&mus_closestJet_);
	}
	if(mus_closestJet_branch == 0 ) {
	cout << "Branch mus_closestJet does not exist." << endl;
	}
	mus_trkidx_branch = 0;
	if (tree->GetAlias("mus_trkidx") != 0) {
		mus_trkidx_branch = tree->GetBranch(tree->GetAlias("mus_trkidx"));
		mus_trkidx_branch->SetAddress(&mus_trkidx_);
	}
	if(mus_trkidx_branch == 0 ) {
	cout << "Branch mus_trkidx does not exist." << endl;
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
	els_pat_genID_branch = 0;
	if (tree->GetAlias("els_pat_genID") != 0) {
		els_pat_genID_branch = tree->GetBranch(tree->GetAlias("els_pat_genID"));
		els_pat_genID_branch->SetAddress(&els_pat_genID_);
	}
	if(els_pat_genID_branch == 0 ) {
	cout << "Branch els_pat_genID does not exist." << endl;
	}
	els_pat_genMotherID_branch = 0;
	if (tree->GetAlias("els_pat_genMotherID") != 0) {
		els_pat_genMotherID_branch = tree->GetBranch(tree->GetAlias("els_pat_genMotherID"));
		els_pat_genMotherID_branch->SetAddress(&els_pat_genMotherID_);
	}
	if(els_pat_genMotherID_branch == 0 ) {
	cout << "Branch els_pat_genMotherID does not exist." << endl;
	}
	jets_pat_genPartonMother_id_branch = 0;
	if (tree->GetAlias("jets_pat_genPartonMother_id") != 0) {
		jets_pat_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("jets_pat_genPartonMother_id"));
		jets_pat_genPartonMother_id_branch->SetAddress(&jets_pat_genPartonMother_id_);
	}
	if(jets_pat_genPartonMother_id_branch == 0 ) {
	cout << "Branch jets_pat_genPartonMother_id does not exist." << endl;
	}
	jets_pat_genParton_id_branch = 0;
	if (tree->GetAlias("jets_pat_genParton_id") != 0) {
		jets_pat_genParton_id_branch = tree->GetBranch(tree->GetAlias("jets_pat_genParton_id"));
		jets_pat_genParton_id_branch->SetAddress(&jets_pat_genParton_id_);
	}
	if(jets_pat_genParton_id_branch == 0 ) {
	cout << "Branch jets_pat_genParton_id does not exist." << endl;
	}
	jets_pat_partonFlavour_branch = 0;
	if (tree->GetAlias("jets_pat_partonFlavour") != 0) {
		jets_pat_partonFlavour_branch = tree->GetBranch(tree->GetAlias("jets_pat_partonFlavour"));
		jets_pat_partonFlavour_branch->SetAddress(&jets_pat_partonFlavour_);
	}
	if(jets_pat_partonFlavour_branch == 0 ) {
	cout << "Branch jets_pat_partonFlavour does not exist." << endl;
	}
	mus_pat_genID_branch = 0;
	if (tree->GetAlias("mus_pat_genID") != 0) {
		mus_pat_genID_branch = tree->GetBranch(tree->GetAlias("mus_pat_genID"));
		mus_pat_genID_branch->SetAddress(&mus_pat_genID_);
	}
	if(mus_pat_genID_branch == 0 ) {
	cout << "Branch mus_pat_genID does not exist." << endl;
	}
	mus_pat_genMotherID_branch = 0;
	if (tree->GetAlias("mus_pat_genMotherID") != 0) {
		mus_pat_genMotherID_branch = tree->GetBranch(tree->GetAlias("mus_pat_genMotherID"));
		mus_pat_genMotherID_branch->SetAddress(&mus_pat_genMotherID_);
	}
	if(mus_pat_genMotherID_branch == 0 ) {
	cout << "Branch mus_pat_genMotherID does not exist." << endl;
	}
	taus_pf_charge_branch = 0;
	if (tree->GetAlias("taus_pf_charge") != 0) {
		taus_pf_charge_branch = tree->GetBranch(tree->GetAlias("taus_pf_charge"));
		taus_pf_charge_branch->SetAddress(&taus_pf_charge_);
	}
	if(taus_pf_charge_branch == 0 ) {
	cout << "Branch taus_pf_charge does not exist." << endl;
	}
	taus_pf_electronPreID_branch = 0;
	if (tree->GetAlias("taus_pf_electronPreID") != 0) {
		taus_pf_electronPreID_branch = tree->GetBranch(tree->GetAlias("taus_pf_electronPreID"));
		taus_pf_electronPreID_branch->SetAddress(&taus_pf_electronPreID_);
	}
	if(taus_pf_electronPreID_branch == 0 ) {
	cout << "Branch taus_pf_electronPreID does not exist." << endl;
	}
	taus_pf_hasMuonReference_branch = 0;
	if (tree->GetAlias("taus_pf_hasMuonReference") != 0) {
		taus_pf_hasMuonReference_branch = tree->GetBranch(tree->GetAlias("taus_pf_hasMuonReference"));
		taus_pf_hasMuonReference_branch->SetAddress(&taus_pf_hasMuonReference_);
	}
	if(taus_pf_hasMuonReference_branch == 0 ) {
	cout << "Branch taus_pf_hasMuonReference does not exist." << endl;
	}
	taus_pf_iso_ncharge_cand_branch = 0;
	if (tree->GetAlias("taus_pf_iso_ncharge_cand") != 0) {
		taus_pf_iso_ncharge_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_iso_ncharge_cand"));
		taus_pf_iso_ncharge_cand_branch->SetAddress(&taus_pf_iso_ncharge_cand_);
	}
	if(taus_pf_iso_ncharge_cand_branch == 0 ) {
	cout << "Branch taus_pf_iso_ncharge_cand does not exist." << endl;
	}
	taus_pf_iso_ngamma_cand_branch = 0;
	if (tree->GetAlias("taus_pf_iso_ngamma_cand") != 0) {
		taus_pf_iso_ngamma_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_iso_ngamma_cand"));
		taus_pf_iso_ngamma_cand_branch->SetAddress(&taus_pf_iso_ngamma_cand_);
	}
	if(taus_pf_iso_ngamma_cand_branch == 0 ) {
	cout << "Branch taus_pf_iso_ngamma_cand does not exist." << endl;
	}
	taus_pf_iso_nneutr_cand_branch = 0;
	if (tree->GetAlias("taus_pf_iso_nneutr_cand") != 0) {
		taus_pf_iso_nneutr_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_iso_nneutr_cand"));
		taus_pf_iso_nneutr_cand_branch->SetAddress(&taus_pf_iso_nneutr_cand_);
	}
	if(taus_pf_iso_nneutr_cand_branch == 0 ) {
	cout << "Branch taus_pf_iso_nneutr_cand does not exist." << endl;
	}
	taus_pf_muonPreID_branch = 0;
	if (tree->GetAlias("taus_pf_muonPreID") != 0) {
		taus_pf_muonPreID_branch = tree->GetBranch(tree->GetAlias("taus_pf_muonPreID"));
		taus_pf_muonPreID_branch->SetAddress(&taus_pf_muonPreID_);
	}
	if(taus_pf_muonPreID_branch == 0 ) {
	cout << "Branch taus_pf_muonPreID does not exist." << endl;
	}
	taus_pf_sig_ncharge_cand_branch = 0;
	if (tree->GetAlias("taus_pf_sig_ncharge_cand") != 0) {
		taus_pf_sig_ncharge_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_sig_ncharge_cand"));
		taus_pf_sig_ncharge_cand_branch->SetAddress(&taus_pf_sig_ncharge_cand_);
	}
	if(taus_pf_sig_ncharge_cand_branch == 0 ) {
	cout << "Branch taus_pf_sig_ncharge_cand does not exist." << endl;
	}
	taus_pf_sig_ngamma_cand_branch = 0;
	if (tree->GetAlias("taus_pf_sig_ngamma_cand") != 0) {
		taus_pf_sig_ngamma_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_sig_ngamma_cand"));
		taus_pf_sig_ngamma_cand_branch->SetAddress(&taus_pf_sig_ngamma_cand_);
	}
	if(taus_pf_sig_ngamma_cand_branch == 0 ) {
	cout << "Branch taus_pf_sig_ngamma_cand does not exist." << endl;
	}
	taus_pf_sig_nneutr_cand_branch = 0;
	if (tree->GetAlias("taus_pf_sig_nneutr_cand") != 0) {
		taus_pf_sig_nneutr_cand_branch = tree->GetBranch(tree->GetAlias("taus_pf_sig_nneutr_cand"));
		taus_pf_sig_nneutr_cand_branch->SetAddress(&taus_pf_sig_nneutr_cand_);
	}
	if(taus_pf_sig_nneutr_cand_branch == 0 ) {
	cout << "Branch taus_pf_sig_nneutr_cand does not exist." << endl;
	}
	mus_tcmet_flag_branch = 0;
	if (tree->GetAlias("mus_tcmet_flag") != 0) {
		mus_tcmet_flag_branch = tree->GetBranch(tree->GetAlias("mus_tcmet_flag"));
		mus_tcmet_flag_branch->SetAddress(&mus_tcmet_flag_);
	}
	if(mus_tcmet_flag_branch == 0 ) {
	cout << "Branch mus_tcmet_flag does not exist." << endl;
	}
	trks_charge_branch = 0;
	if (tree->GetAlias("trks_charge") != 0) {
		trks_charge_branch = tree->GetBranch(tree->GetAlias("trks_charge"));
		trks_charge_branch->SetAddress(&trks_charge_);
	}
	if(trks_charge_branch == 0 ) {
	cout << "Branch trks_charge does not exist." << endl;
	}
	trks_lostHits_branch = 0;
	if (tree->GetAlias("trks_lostHits") != 0) {
		trks_lostHits_branch = tree->GetBranch(tree->GetAlias("trks_lostHits"));
		trks_lostHits_branch->SetAddress(&trks_lostHits_);
	}
	if(trks_lostHits_branch == 0 ) {
	cout << "Branch trks_lostHits does not exist." << endl;
	}
	trks_validHits_branch = 0;
	if (tree->GetAlias("trks_validHits") != 0) {
		trks_validHits_branch = tree->GetBranch(tree->GetAlias("trks_validHits"));
		trks_validHits_branch->SetAddress(&trks_validHits_);
	}
	if(trks_validHits_branch == 0 ) {
	cout << "Branch trks_validHits does not exist." << endl;
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
	hlt2ele10LWR_id_branch = 0;
	if (tree->GetAlias("hlt2ele10LWR_id") != 0) {
		hlt2ele10LWR_id_branch = tree->GetBranch(tree->GetAlias("hlt2ele10LWR_id"));
		hlt2ele10LWR_id_branch->SetAddress(&hlt2ele10LWR_id_);
	}
	if(hlt2ele10LWR_id_branch == 0 ) {
	cout << "Branch hlt2ele10LWR_id does not exist." << endl;
	}
	hlt2ele10LWR_tid_branch = 0;
	if (tree->GetAlias("hlt2ele10LWR_tid") != 0) {
		hlt2ele10LWR_tid_branch = tree->GetBranch(tree->GetAlias("hlt2ele10LWR_tid"));
		hlt2ele10LWR_tid_branch->SetAddress(&hlt2ele10LWR_tid_);
	}
	if(hlt2ele10LWR_tid_branch == 0 ) {
	cout << "Branch hlt2ele10LWR_tid does not exist." << endl;
	}
	hlt2mu3_id_branch = 0;
	if (tree->GetAlias("hlt2mu3_id") != 0) {
		hlt2mu3_id_branch = tree->GetBranch(tree->GetAlias("hlt2mu3_id"));
		hlt2mu3_id_branch->SetAddress(&hlt2mu3_id_);
	}
	if(hlt2mu3_id_branch == 0 ) {
	cout << "Branch hlt2mu3_id does not exist." << endl;
	}
	hlt2mu3_tid_branch = 0;
	if (tree->GetAlias("hlt2mu3_tid") != 0) {
		hlt2mu3_tid_branch = tree->GetBranch(tree->GetAlias("hlt2mu3_tid"));
		hlt2mu3_tid_branch->SetAddress(&hlt2mu3_tid_);
	}
	if(hlt2mu3_tid_branch == 0 ) {
	cout << "Branch hlt2mu3_tid does not exist." << endl;
	}
	hltLisoele18LWR_id_branch = 0;
	if (tree->GetAlias("hltLisoele18LWR_id") != 0) {
		hltLisoele18LWR_id_branch = tree->GetBranch(tree->GetAlias("hltLisoele18LWR_id"));
		hltLisoele18LWR_id_branch->SetAddress(&hltLisoele18LWR_id_);
	}
	if(hltLisoele18LWR_id_branch == 0 ) {
	cout << "Branch hltLisoele18LWR_id does not exist." << endl;
	}
	hltLisoele18LWR_tid_branch = 0;
	if (tree->GetAlias("hltLisoele18LWR_tid") != 0) {
		hltLisoele18LWR_tid_branch = tree->GetBranch(tree->GetAlias("hltLisoele18LWR_tid"));
		hltLisoele18LWR_tid_branch->SetAddress(&hltLisoele18LWR_tid_);
	}
	if(hltLisoele18LWR_tid_branch == 0 ) {
	cout << "Branch hltLisoele18LWR_tid does not exist." << endl;
	}
	hltisoele18R_id_branch = 0;
	if (tree->GetAlias("hltisoele18R_id") != 0) {
		hltisoele18R_id_branch = tree->GetBranch(tree->GetAlias("hltisoele18R_id"));
		hltisoele18R_id_branch->SetAddress(&hltisoele18R_id_);
	}
	if(hltisoele18R_id_branch == 0 ) {
	cout << "Branch hltisoele18R_id does not exist." << endl;
	}
	hltisoele18R_tid_branch = 0;
	if (tree->GetAlias("hltisoele18R_tid") != 0) {
		hltisoele18R_tid_branch = tree->GetBranch(tree->GetAlias("hltisoele18R_tid"));
		hltisoele18R_tid_branch->SetAddress(&hltisoele18R_tid_);
	}
	if(hltisoele18R_tid_branch == 0 ) {
	cout << "Branch hltisoele18R_tid does not exist." << endl;
	}
	hltjet30_id_branch = 0;
	if (tree->GetAlias("hltjet30_id") != 0) {
		hltjet30_id_branch = tree->GetBranch(tree->GetAlias("hltjet30_id"));
		hltjet30_id_branch->SetAddress(&hltjet30_id_);
	}
	if(hltjet30_id_branch == 0 ) {
	cout << "Branch hltjet30_id does not exist." << endl;
	}
	hltjet30_tid_branch = 0;
	if (tree->GetAlias("hltjet30_tid") != 0) {
		hltjet30_tid_branch = tree->GetBranch(tree->GetAlias("hltjet30_tid"));
		hltjet30_tid_branch->SetAddress(&hltjet30_tid_);
	}
	if(hltjet30_tid_branch == 0 ) {
	cout << "Branch hltjet30_tid does not exist." << endl;
	}
	hltl1jet15_id_branch = 0;
	if (tree->GetAlias("hltl1jet15_id") != 0) {
		hltl1jet15_id_branch = tree->GetBranch(tree->GetAlias("hltl1jet15_id"));
		hltl1jet15_id_branch->SetAddress(&hltl1jet15_id_);
	}
	if(hltl1jet15_id_branch == 0 ) {
	cout << "Branch hltl1jet15_id does not exist." << endl;
	}
	hltl1jet15_tid_branch = 0;
	if (tree->GetAlias("hltl1jet15_tid") != 0) {
		hltl1jet15_tid_branch = tree->GetBranch(tree->GetAlias("hltl1jet15_tid"));
		hltl1jet15_tid_branch->SetAddress(&hltl1jet15_tid_);
	}
	if(hltl1jet15_tid_branch == 0 ) {
	cout << "Branch hltl1jet15_tid does not exist." << endl;
	}
	hltl1met20_id_branch = 0;
	if (tree->GetAlias("hltl1met20_id") != 0) {
		hltl1met20_id_branch = tree->GetBranch(tree->GetAlias("hltl1met20_id"));
		hltl1met20_id_branch->SetAddress(&hltl1met20_id_);
	}
	if(hltl1met20_id_branch == 0 ) {
	cout << "Branch hltl1met20_id does not exist." << endl;
	}
	hltl1met20_tid_branch = 0;
	if (tree->GetAlias("hltl1met20_tid") != 0) {
		hltl1met20_tid_branch = tree->GetBranch(tree->GetAlias("hltl1met20_tid"));
		hltl1met20_tid_branch->SetAddress(&hltl1met20_tid_);
	}
	if(hltl1met20_tid_branch == 0 ) {
	cout << "Branch hltl1met20_tid does not exist." << endl;
	}
	hltl1mu_id_branch = 0;
	if (tree->GetAlias("hltl1mu_id") != 0) {
		hltl1mu_id_branch = tree->GetBranch(tree->GetAlias("hltl1mu_id"));
		hltl1mu_id_branch->SetAddress(&hltl1mu_id_);
	}
	if(hltl1mu_id_branch == 0 ) {
	cout << "Branch hltl1mu_id does not exist." << endl;
	}
	hltl1mu_tid_branch = 0;
	if (tree->GetAlias("hltl1mu_tid") != 0) {
		hltl1mu_tid_branch = tree->GetBranch(tree->GetAlias("hltl1mu_tid"));
		hltl1mu_tid_branch->SetAddress(&hltl1mu_tid_);
	}
	if(hltl1mu_tid_branch == 0 ) {
	cout << "Branch hltl1mu_tid does not exist." << endl;
	}
	hltl2mu9_id_branch = 0;
	if (tree->GetAlias("hltl2mu9_id") != 0) {
		hltl2mu9_id_branch = tree->GetBranch(tree->GetAlias("hltl2mu9_id"));
		hltl2mu9_id_branch->SetAddress(&hltl2mu9_id_);
	}
	if(hltl2mu9_id_branch == 0 ) {
	cout << "Branch hltl2mu9_id does not exist." << endl;
	}
	hltl2mu9_tid_branch = 0;
	if (tree->GetAlias("hltl2mu9_tid") != 0) {
		hltl2mu9_tid_branch = tree->GetBranch(tree->GetAlias("hltl2mu9_tid"));
		hltl2mu9_tid_branch->SetAddress(&hltl2mu9_tid_);
	}
	if(hltl2mu9_tid_branch == 0 ) {
	cout << "Branch hltl2mu9_tid does not exist." << endl;
	}
	hltmet25_id_branch = 0;
	if (tree->GetAlias("hltmet25_id") != 0) {
		hltmet25_id_branch = tree->GetBranch(tree->GetAlias("hltmet25_id"));
		hltmet25_id_branch->SetAddress(&hltmet25_id_);
	}
	if(hltmet25_id_branch == 0 ) {
	cout << "Branch hltmet25_id does not exist." << endl;
	}
	hltmet25_tid_branch = 0;
	if (tree->GetAlias("hltmet25_tid") != 0) {
		hltmet25_tid_branch = tree->GetBranch(tree->GetAlias("hltmet25_tid"));
		hltmet25_tid_branch->SetAddress(&hltmet25_tid_);
	}
	if(hltmet25_tid_branch == 0 ) {
	cout << "Branch hltmet25_tid does not exist." << endl;
	}
	hltmu11_id_branch = 0;
	if (tree->GetAlias("hltmu11_id") != 0) {
		hltmu11_id_branch = tree->GetBranch(tree->GetAlias("hltmu11_id"));
		hltmu11_id_branch->SetAddress(&hltmu11_id_);
	}
	if(hltmu11_id_branch == 0 ) {
	cout << "Branch hltmu11_id does not exist." << endl;
	}
	hltmu11_tid_branch = 0;
	if (tree->GetAlias("hltmu11_tid") != 0) {
		hltmu11_tid_branch = tree->GetBranch(tree->GetAlias("hltmu11_tid"));
		hltmu11_tid_branch->SetAddress(&hltmu11_tid_);
	}
	if(hltmu11_tid_branch == 0 ) {
	cout << "Branch hltmu11_tid does not exist." << endl;
	}
	hltmu9_id_branch = 0;
	if (tree->GetAlias("hltmu9_id") != 0) {
		hltmu9_id_branch = tree->GetBranch(tree->GetAlias("hltmu9_id"));
		hltmu9_id_branch->SetAddress(&hltmu9_id_);
	}
	if(hltmu9_id_branch == 0 ) {
	cout << "Branch hltmu9_id does not exist." << endl;
	}
	hltmu9_tid_branch = 0;
	if (tree->GetAlias("hltmu9_tid") != 0) {
		hltmu9_tid_branch = tree->GetBranch(tree->GetAlias("hltmu9_tid"));
		hltmu9_tid_branch->SetAddress(&hltmu9_tid_);
	}
	if(hltmu9_tid_branch == 0 ) {
	cout << "Branch hltmu9_tid does not exist." << endl;
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
	hyp_jets_idx_branch = 0;
	if (tree->GetAlias("hyp_jets_idx") != 0) {
		hyp_jets_idx_branch = tree->GetBranch(tree->GetAlias("hyp_jets_idx"));
		hyp_jets_idx_branch->SetAddress(&hyp_jets_idx_);
	}
	if(hyp_jets_idx_branch == 0 ) {
	cout << "Branch hyp_jets_idx does not exist." << endl;
	}
	hyp_other_jets_idx_branch = 0;
	if (tree->GetAlias("hyp_other_jets_idx") != 0) {
		hyp_other_jets_idx_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_idx"));
		hyp_other_jets_idx_branch->SetAddress(&hyp_other_jets_idx_);
	}
	if(hyp_other_jets_idx_branch == 0 ) {
	cout << "Branch hyp_other_jets_idx does not exist." << endl;
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
	evt_ngenjets_branch = 0;
	if (tree->GetAlias("evt_ngenjets") != 0) {
		evt_ngenjets_branch = tree->GetBranch(tree->GetAlias("evt_ngenjets"));
		evt_ngenjets_branch->SetAddress(&evt_ngenjets_);
	}
	if(evt_ngenjets_branch == 0 ) {
	cout << "Branch evt_ngenjets does not exist." << endl;
	}
	evt_njets_branch = 0;
	if (tree->GetAlias("evt_njets") != 0) {
		evt_njets_branch = tree->GetBranch(tree->GetAlias("evt_njets"));
		evt_njets_branch->SetAddress(&evt_njets_);
	}
	if(evt_njets_branch == 0 ) {
	cout << "Branch evt_njets does not exist." << endl;
	}
	evt_njpts_branch = 0;
	if (tree->GetAlias("evt_njpts") != 0) {
		evt_njpts_branch = tree->GetBranch(tree->GetAlias("evt_njpts"));
		evt_njpts_branch->SetAddress(&evt_njpts_);
	}
	if(evt_njpts_branch == 0 ) {
	cout << "Branch evt_njpts does not exist." << endl;
	}
	evt_nphotons_branch = 0;
	if (tree->GetAlias("evt_nphotons") != 0) {
		evt_nphotons_branch = tree->GetBranch(tree->GetAlias("evt_nphotons"));
		evt_nphotons_branch->SetAddress(&evt_nphotons_);
	}
	if(evt_nphotons_branch == 0 ) {
	cout << "Branch evt_nphotons does not exist." << endl;
	}
	evt_nscs_branch = 0;
	if (tree->GetAlias("evt_nscs") != 0) {
		evt_nscs_branch = tree->GetBranch(tree->GetAlias("evt_nscs"));
		evt_nscs_branch->SetAddress(&evt_nscs_);
	}
	if(evt_nscs_branch == 0 ) {
	cout << "Branch evt_nscs does not exist." << endl;
	}
	evt_ntrkjets_branch = 0;
	if (tree->GetAlias("evt_ntrkjets") != 0) {
		evt_ntrkjets_branch = tree->GetBranch(tree->GetAlias("evt_ntrkjets"));
		evt_ntrkjets_branch->SetAddress(&evt_ntrkjets_);
	}
	if(evt_ntrkjets_branch == 0 ) {
	cout << "Branch evt_ntrkjets does not exist." << endl;
	}
	evt_nvtxs_branch = 0;
	if (tree->GetAlias("evt_nvtxs") != 0) {
		evt_nvtxs_branch = tree->GetBranch(tree->GetAlias("evt_nvtxs"));
		evt_nvtxs_branch->SetAddress(&evt_nvtxs_);
	}
	if(evt_nvtxs_branch == 0 ) {
	cout << "Branch evt_nvtxs does not exist." << endl;
	}
	hyp_quadlep_bucket_branch = 0;
	if (tree->GetAlias("hyp_quadlep_bucket") != 0) {
		hyp_quadlep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_bucket"));
		hyp_quadlep_bucket_branch->SetAddress(&hyp_quadlep_bucket_);
	}
	if(hyp_quadlep_bucket_branch == 0 ) {
	cout << "Branch hyp_quadlep_bucket does not exist." << endl;
	}
	hyp_quadlep_first_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_index") != 0) {
		hyp_quadlep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_index"));
		hyp_quadlep_first_index_branch->SetAddress(&hyp_quadlep_first_index_);
	}
	if(hyp_quadlep_first_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_first_index does not exist." << endl;
	}
	hyp_quadlep_fourth_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_index") != 0) {
		hyp_quadlep_fourth_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_index"));
		hyp_quadlep_fourth_index_branch->SetAddress(&hyp_quadlep_fourth_index_);
	}
	if(hyp_quadlep_fourth_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_fourth_index does not exist." << endl;
	}
	hyp_quadlep_second_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_index") != 0) {
		hyp_quadlep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_index"));
		hyp_quadlep_second_index_branch->SetAddress(&hyp_quadlep_second_index_);
	}
	if(hyp_quadlep_second_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_second_index does not exist." << endl;
	}
	hyp_quadlep_third_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_index") != 0) {
		hyp_quadlep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_index"));
		hyp_quadlep_third_index_branch->SetAddress(&hyp_quadlep_third_index_);
	}
	if(hyp_quadlep_third_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_third_index does not exist." << endl;
	}
	hyp_trilep_bucket_branch = 0;
	if (tree->GetAlias("hyp_trilep_bucket") != 0) {
		hyp_trilep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_bucket"));
		hyp_trilep_bucket_branch->SetAddress(&hyp_trilep_bucket_);
	}
	if(hyp_trilep_bucket_branch == 0 ) {
	cout << "Branch hyp_trilep_bucket does not exist." << endl;
	}
	hyp_trilep_first_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_index") != 0) {
		hyp_trilep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_index"));
		hyp_trilep_first_index_branch->SetAddress(&hyp_trilep_first_index_);
	}
	if(hyp_trilep_first_index_branch == 0 ) {
	cout << "Branch hyp_trilep_first_index does not exist." << endl;
	}
	hyp_trilep_second_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_index") != 0) {
		hyp_trilep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_index"));
		hyp_trilep_second_index_branch->SetAddress(&hyp_trilep_second_index_);
	}
	if(hyp_trilep_second_index_branch == 0 ) {
	cout << "Branch hyp_trilep_second_index does not exist." << endl;
	}
	hyp_trilep_third_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_index") != 0) {
		hyp_trilep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_index"));
		hyp_trilep_third_index_branch->SetAddress(&hyp_trilep_third_index_);
	}
	if(hyp_trilep_third_index_branch == 0 ) {
	cout << "Branch hyp_trilep_third_index does not exist." << endl;
	}
	els_pat_flag_branch = 0;
	if (tree->GetAlias("els_pat_flag") != 0) {
		els_pat_flag_branch = tree->GetBranch(tree->GetAlias("els_pat_flag"));
		els_pat_flag_branch->SetAddress(&els_pat_flag_);
	}
	if(els_pat_flag_branch == 0 ) {
	cout << "Branch els_pat_flag does not exist." << endl;
	}
	jets_pat_flag_branch = 0;
	if (tree->GetAlias("jets_pat_flag") != 0) {
		jets_pat_flag_branch = tree->GetBranch(tree->GetAlias("jets_pat_flag"));
		jets_pat_flag_branch->SetAddress(&jets_pat_flag_);
	}
	if(jets_pat_flag_branch == 0 ) {
	cout << "Branch jets_pat_flag does not exist." << endl;
	}
	mus_pat_flag_branch = 0;
	if (tree->GetAlias("mus_pat_flag") != 0) {
		mus_pat_flag_branch = tree->GetBranch(tree->GetAlias("mus_pat_flag"));
		mus_pat_flag_branch->SetAddress(&mus_pat_flag_);
	}
	if(mus_pat_flag_branch == 0 ) {
	cout << "Branch mus_pat_flag does not exist." << endl;
	}
	evt_scale1fb_branch = 0;
	if (tree->GetAlias("evt_scale1fb") != 0) {
		evt_scale1fb_branch = tree->GetBranch(tree->GetAlias("evt_scale1fb"));
		evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
	}
	if(evt_scale1fb_branch == 0 ) {
	cout << "Branch evt_scale1fb does not exist." << endl;
	}
	evt_nEvts_branch = 0;
	if (tree->GetAlias("evt_nEvts") != 0) {
		evt_nEvts_branch = tree->GetBranch(tree->GetAlias("evt_nEvts"));
		evt_nEvts_branch->SetAddress(&evt_nEvts_);
	}
	if(evt_nEvts_branch == 0 ) {
	cout << "Branch evt_nEvts does not exist." << endl;
	}
	evt_filt_eff_branch = 0;
	if (tree->GetAlias("evt_filt_eff") != 0) {
		evt_filt_eff_branch = tree->GetBranch(tree->GetAlias("evt_filt_eff"));
		evt_filt_eff_branch->SetAddress(&evt_filt_eff_);
	}
	if(evt_filt_eff_branch == 0 ) {
	cout << "Branch evt_filt_eff does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		evt_dataset_isLoaded = false;
		evt_HLT_trigNames_isLoaded = false;
		evt_L1_trigNames_isLoaded = false;
		genps_pthat_isLoaded = false;
		vtxs_position_isLoaded = false;
		evt_bsp4_isLoaded = false;
		l1met_p4_isLoaded = false;
		taus_calo_leadtrk_p4_isLoaded = false;
		taus_calo_p4_isLoaded = false;
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
		genjets_p4_isLoaded = false;
		genps_lepdaughter_p4_isLoaded = false;
		genps_p4_isLoaded = false;
		genps_prod_vtx_isLoaded = false;
		hyp_ll_mc_p4_isLoaded = false;
		hyp_ll_p4_isLoaded = false;
		hyp_ll_trk_p4_isLoaded = false;
		hyp_lt_mc_p4_isLoaded = false;
		hyp_lt_p4_isLoaded = false;
		hyp_lt_trk_p4_isLoaded = false;
		hyp_p4_isLoaded = false;
		jets_p4_isLoaded = false;
		jpts_p4_isLoaded = false;
		l1emiso_p4_isLoaded = false;
		l1emnoiso_p4_isLoaded = false;
		l1jetsc_p4_isLoaded = false;
		l1jetsf_p4_isLoaded = false;
		l1jetst_p4_isLoaded = false;
		l1mus_p4_isLoaded = false;
		mus_gfit_outerPos_p4_isLoaded = false;
		mus_p4_isLoaded = false;
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
		taus_pf_lead_chargecand_p4_isLoaded = false;
		taus_pf_lead_neutrcand_p4_isLoaded = false;
		taus_pf_leadtrk_p4_isLoaded = false;
		taus_pf_p4_isLoaded = false;
		photons_p4_isLoaded = false;
		scs_p4_isLoaded = false;
		scs_pos_isLoaded = false;
		scs_vtx_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		trks_vertex_p4_isLoaded = false;
		hlt2ele10LWR_p4_isLoaded = false;
		hlt2mu3_p4_isLoaded = false;
		hltLisoele18LWR_p4_isLoaded = false;
		hltisoele18R_p4_isLoaded = false;
		hltjet30_p4_isLoaded = false;
		hltl1jet15_p4_isLoaded = false;
		hltl1met20_p4_isLoaded = false;
		hltl1mu_p4_isLoaded = false;
		hltl2mu9_p4_isLoaded = false;
		hltmet25_p4_isLoaded = false;
		hltmu11_p4_isLoaded = false;
		hltmu9_p4_isLoaded = false;
		trkjets_p4_isLoaded = false;
		hyp_jets_p4_isLoaded = false;
		hyp_other_jets_p4_isLoaded = false;
		els_inner_positionx_isLoaded = false;
		els_inner_positiony_isLoaded = false;
		els_inner_positionz_isLoaded = false;
		jets_closestElectron_DR_isLoaded = false;
		jets_closestMuon_DR_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_width_isLoaded = false;
		evt_bs_widthErr_isLoaded = false;
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
		l1met_etHad_isLoaded = false;
		l1met_etTot_isLoaded = false;
		l1met_met_isLoaded = false;
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
		evt_pfsumet_isLoaded = false;
		evt_tcmet_isLoaded = false;
		evt_tcmetPhi_isLoaded = false;
		evt_tcsumet_isLoaded = false;
		jets_combinedSecondaryVertexBJetTag_isLoaded = false;
		jets_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		jets_impactParameterMVABJetTag_isLoaded = false;
		jets_jetBProbabilityBJetTag_isLoaded = false;
		jets_jetProbabilityBJetTag_isLoaded = false;
		jets_simpleSecondaryVertexBJetTag_isLoaded = false;
		jets_softElectronBJetTag_isLoaded = false;
		jets_softMuonBJetTag_isLoaded = false;
		jets_softMuonNoIPBJetTag_isLoaded = false;
		jets_trackCountingHighEffBJetTag_isLoaded = false;
		jets_trackCountingHighPurBJetTag_isLoaded = false;
		taus_calo_isolationECALhitsEtSum_isLoaded = false;
		taus_calo_isolationtrksPtSum_isLoaded = false;
		taus_calo_leadtrk_HCAL3x3hitsEtSum_isLoaded = false;
		taus_calo_leadtrk_HCAL3x3hottesthitDEta_isLoaded = false;
		taus_calo_leadtrk_Signed_Sipt_isLoaded = false;
		taus_calo_leadtrk_chi2_isLoaded = false;
		taus_calo_leadtrk_d0_isLoaded = false;
		taus_calo_leadtrk_lostHits_isLoaded = false;
		taus_calo_leadtrk_ndof_isLoaded = false;
		taus_calo_leadtrk_validHits_isLoaded = false;
		taus_calo_leadtrk_z0_isLoaded = false;
		taus_calo_maximumHCALhitEt_isLoaded = false;
		taus_calo_signaltrksInvariantMass_isLoaded = false;
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
		els_ecalJuraIso_isLoaded = false;
		els_ecalJuraTowerIso_isLoaded = false;
		els_hcalConeIso_isLoaded = false;
		els_jetdr_isLoaded = false;
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
		els_ecalIso_isLoaded = false;
		els_egamma_looseId_isLoaded = false;
		els_egamma_robustHighEnergy_isLoaded = false;
		els_egamma_robustLooseId_isLoaded = false;
		els_egamma_robustTightId_isLoaded = false;
		els_egamma_tightId_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_etaSC_isLoaded = false;
		els_fBrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_hcalIso_isLoaded = false;
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
		els_tkIso_isLoaded = false;
		els_vertexphi_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		els_z0corr_isLoaded = false;
		hyp_Ht_isLoaded = false;
		hyp_dPhi_nJet_muCorrMet_isLoaded = false;
		hyp_dPhi_nJet_tcMet_isLoaded = false;
		hyp_dPhi_nJet_type1Met_isLoaded = false;
		hyp_dPhi_nJet_unCorrMet_isLoaded = false;
		hyp_ll_chi2_isLoaded = false;
		hyp_ll_d0_isLoaded = false;
		hyp_ll_d0Err_isLoaded = false;
		hyp_ll_d0corr_isLoaded = false;
		hyp_ll_dPhi_muCorrMet_isLoaded = false;
		hyp_ll_dPhi_tcMet_isLoaded = false;
		hyp_ll_dPhi_type1Met_isLoaded = false;
		hyp_ll_dPhi_unCorrMet_isLoaded = false;
		hyp_ll_etaErr_isLoaded = false;
		hyp_ll_iso_isLoaded = false;
		hyp_ll_ndof_isLoaded = false;
		hyp_ll_outerEta_isLoaded = false;
		hyp_ll_outerPhi_isLoaded = false;
		hyp_ll_phiErr_isLoaded = false;
		hyp_ll_ptErr_isLoaded = false;
		hyp_ll_tkIso_isLoaded = false;
		hyp_ll_vertexphi_isLoaded = false;
		hyp_ll_z0_isLoaded = false;
		hyp_ll_z0Err_isLoaded = false;
		hyp_ll_z0corr_isLoaded = false;
		hyp_lt_chi2_isLoaded = false;
		hyp_lt_d0_isLoaded = false;
		hyp_lt_d0Err_isLoaded = false;
		hyp_lt_d0corr_isLoaded = false;
		hyp_lt_dPhi_muCorrMet_isLoaded = false;
		hyp_lt_dPhi_tcMet_isLoaded = false;
		hyp_lt_dPhi_type1Met_isLoaded = false;
		hyp_lt_dPhi_unCorrMet_isLoaded = false;
		hyp_lt_etaErr_isLoaded = false;
		hyp_lt_iso_isLoaded = false;
		hyp_lt_ndof_isLoaded = false;
		hyp_lt_outerEta_isLoaded = false;
		hyp_lt_outerPhi_isLoaded = false;
		hyp_lt_phiErr_isLoaded = false;
		hyp_lt_ptErr_isLoaded = false;
		hyp_lt_tkIso_isLoaded = false;
		hyp_lt_vertexphi_isLoaded = false;
		hyp_lt_z0_isLoaded = false;
		hyp_lt_z0Err_isLoaded = false;
		hyp_lt_z0corr_isLoaded = false;
		hyp_sumJetPt_isLoaded = false;
		jets_EMFcor_isLoaded = false;
		jets_chFrac_isLoaded = false;
		jets_cor_isLoaded = false;
		jets_emFrac_isLoaded = false;
		jpts_cor_isLoaded = false;
		jpts_emFrac_isLoaded = false;
		jpts_jet_cor_isLoaded = false;
		mus_met_deltax_isLoaded = false;
		mus_met_deltay_isLoaded = false;
		mus_eledr_isLoaded = false;
		mus_jetdr_isLoaded = false;
		mus_trkdr_isLoaded = false;
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
		jets_pat_bCorrF_isLoaded = false;
		jets_pat_cCorrF_isLoaded = false;
		jets_pat_combinedSecondaryVertexBJetTag_isLoaded = false;
		jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded = false;
		jets_pat_coneIsolationTauJetTag_isLoaded = false;
		jets_pat_gluCorrF_isLoaded = false;
		jets_pat_impactParameterMVABJetTag_isLoaded = false;
		jets_pat_jetBProbabilityBJetTag_isLoaded = false;
		jets_pat_jetCharge_isLoaded = false;
		jets_pat_jetProbabilityBJetTag_isLoaded = false;
		jets_pat_noCorrF_isLoaded = false;
		jets_pat_simpleSecondaryVertexBJetTag_isLoaded = false;
		jets_pat_softElectronBJetTag_isLoaded = false;
		jets_pat_softMuonBJetTag_isLoaded = false;
		jets_pat_softMuonNoIPBJetTag_isLoaded = false;
		jets_pat_trackCountingHighEffBJetTag_isLoaded = false;
		jets_pat_trackCountingHighPurBJetTag_isLoaded = false;
		jets_pat_udsCorrF_isLoaded = false;
		mus_pat_caloIso_isLoaded = false;
		mus_pat_ecalIso_isLoaded = false;
		mus_pat_ecalvetoDep_isLoaded = false;
		mus_pat_hcalIso_isLoaded = false;
		mus_pat_hcalvetoDep_isLoaded = false;
		mus_pat_trackIso_isLoaded = false;
		mus_pat_vetoDep_isLoaded = false;
		taus_pf_bremsRecoveryEOverPLead_isLoaded = false;
		taus_pf_caloComp_isLoaded = false;
		taus_pf_ecalStripSumEOverPLead_isLoaded = false;
		taus_pf_emf_isLoaded = false;
		taus_pf_hcal3x3OverPLead_isLoaded = false;
		taus_pf_hcalMaxOverPLead_isLoaded = false;
		taus_pf_hcalTotOverPLead_isLoaded = false;
		taus_pf_isolationchargecandPtSum_isLoaded = false;
		taus_pf_isolationgammacandEtSum_isLoaded = false;
		taus_pf_lead_chargecand_Signed_Sipt_isLoaded = false;
		taus_pf_leadtrk_chi2_isLoaded = false;
		taus_pf_leadtrk_d0_isLoaded = false;
		taus_pf_leadtrk_lostHits_isLoaded = false;
		taus_pf_leadtrk_ndof_isLoaded = false;
		taus_pf_leadtrk_validHits_isLoaded = false;
		taus_pf_leadtrk_z0_isLoaded = false;
		taus_pf_maximumHCALPFClusterEt_isLoaded = false;
		taus_pf_segComp_isLoaded = false;
		photons_e1x5_isLoaded = false;
		photons_e2x5Max_isLoaded = false;
		photons_e3x3_isLoaded = false;
		photons_e5x5_isLoaded = false;
		photons_eMax_isLoaded = false;
		photons_eSC_isLoaded = false;
		photons_eSCPresh_isLoaded = false;
		photons_eSCRaw_isLoaded = false;
		photons_eSeed_isLoaded = false;
		photons_ecalIso_isLoaded = false;
		photons_hOverE_isLoaded = false;
		photons_hcalIso_isLoaded = false;
		photons_sigmaEtaEta_isLoaded = false;
		photons_sigmaIEtaIEta_isLoaded = false;
		photons_sigmaIPhiIPhi_isLoaded = false;
		photons_sigmaPhiPhi_isLoaded = false;
		photons_tkIso_isLoaded = false;
		scs_clustersSize_isLoaded = false;
		scs_crystalsSize_isLoaded = false;
		scs_e1x3_isLoaded = false;
		scs_e1x5_isLoaded = false;
		scs_e2x2_isLoaded = false;
		scs_e2x5Max_isLoaded = false;
		scs_e3x1_isLoaded = false;
		scs_e3x2_isLoaded = false;
		scs_e3x3_isLoaded = false;
		scs_e4x4_isLoaded = false;
		scs_e5x5_isLoaded = false;
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
		mus_tcmet_deltax_isLoaded = false;
		mus_tcmet_deltay_isLoaded = false;
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_d0corr_isLoaded = false;
		trks_d0corrPhi_isLoaded = false;
		trks_etaErr_isLoaded = false;
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
		evt_nl1emiso_isLoaded = false;
		evt_nl1emnoiso_isLoaded = false;
		evt_nl1jetsc_isLoaded = false;
		evt_nl1jetsf_isLoaded = false;
		evt_nl1jetst_isLoaded = false;
		evt_nl1mus_isLoaded = false;
		pdfinfo_id1_isLoaded = false;
		pdfinfo_id2_isLoaded = false;
		taus_calo_charge_isLoaded = false;
		taus_calo_iso_ntrks_isLoaded = false;
		taus_calo_sig_ntrks_isLoaded = false;
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
		els_closestJet_isLoaded = false;
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
		hyp_ll_charge_isLoaded = false;
		hyp_ll_id_isLoaded = false;
		hyp_ll_index_isLoaded = false;
		hyp_ll_lostHits_isLoaded = false;
		hyp_ll_mc_id_isLoaded = false;
		hyp_ll_mc_motherid_isLoaded = false;
		hyp_ll_validHits_isLoaded = false;
		hyp_lt_charge_isLoaded = false;
		hyp_lt_id_isLoaded = false;
		hyp_lt_index_isLoaded = false;
		hyp_lt_lostHits_isLoaded = false;
		hyp_lt_mc_id_isLoaded = false;
		hyp_lt_mc_motherid_isLoaded = false;
		hyp_lt_validHits_isLoaded = false;
		hyp_njets_isLoaded = false;
		hyp_nojets_isLoaded = false;
		hyp_type_isLoaded = false;
		hyp_quadlep_first_type_isLoaded = false;
		hyp_quadlep_fourth_type_isLoaded = false;
		hyp_quadlep_second_type_isLoaded = false;
		hyp_quadlep_third_type_isLoaded = false;
		hyp_trilep_first_type_isLoaded = false;
		hyp_trilep_second_type_isLoaded = false;
		hyp_trilep_third_type_isLoaded = false;
		jets_closestElectron_isLoaded = false;
		jets_closestMuon_isLoaded = false;
		l1emiso_ieta_isLoaded = false;
		l1emiso_iphi_isLoaded = false;
		l1emiso_rawId_isLoaded = false;
		l1emiso_type_isLoaded = false;
		l1emnoiso_ieta_isLoaded = false;
		l1emnoiso_iphi_isLoaded = false;
		l1emnoiso_rawId_isLoaded = false;
		l1emnoiso_type_isLoaded = false;
		l1jetsc_ieta_isLoaded = false;
		l1jetsc_iphi_isLoaded = false;
		l1jetsc_rawId_isLoaded = false;
		l1jetsc_type_isLoaded = false;
		l1jetsf_ieta_isLoaded = false;
		l1jetsf_iphi_isLoaded = false;
		l1jetsf_rawId_isLoaded = false;
		l1jetsf_type_isLoaded = false;
		l1jetst_ieta_isLoaded = false;
		l1jetst_iphi_isLoaded = false;
		l1jetst_rawId_isLoaded = false;
		l1jetst_type_isLoaded = false;
		l1mus_flags_isLoaded = false;
		l1mus_q_isLoaded = false;
		l1mus_qual_isLoaded = false;
		l1mus_qualFlags_isLoaded = false;
		mus_met_flag_isLoaded = false;
		mus_closestEle_isLoaded = false;
		mus_closestJet_isLoaded = false;
		mus_trkidx_isLoaded = false;
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
		els_pat_genID_isLoaded = false;
		els_pat_genMotherID_isLoaded = false;
		jets_pat_genPartonMother_id_isLoaded = false;
		jets_pat_genParton_id_isLoaded = false;
		jets_pat_partonFlavour_isLoaded = false;
		mus_pat_genID_isLoaded = false;
		mus_pat_genMotherID_isLoaded = false;
		taus_pf_charge_isLoaded = false;
		taus_pf_electronPreID_isLoaded = false;
		taus_pf_hasMuonReference_isLoaded = false;
		taus_pf_iso_ncharge_cand_isLoaded = false;
		taus_pf_iso_ngamma_cand_isLoaded = false;
		taus_pf_iso_nneutr_cand_isLoaded = false;
		taus_pf_muonPreID_isLoaded = false;
		taus_pf_sig_ncharge_cand_isLoaded = false;
		taus_pf_sig_ngamma_cand_isLoaded = false;
		taus_pf_sig_nneutr_cand_isLoaded = false;
		mus_tcmet_flag_isLoaded = false;
		trks_charge_isLoaded = false;
		trks_lostHits_isLoaded = false;
		trks_validHits_isLoaded = false;
		trks_elsidx_isLoaded = false;
		trk_musidx_isLoaded = false;
		hlt2ele10LWR_id_isLoaded = false;
		hlt2ele10LWR_tid_isLoaded = false;
		hlt2mu3_id_isLoaded = false;
		hlt2mu3_tid_isLoaded = false;
		hltLisoele18LWR_id_isLoaded = false;
		hltLisoele18LWR_tid_isLoaded = false;
		hltisoele18R_id_isLoaded = false;
		hltisoele18R_tid_isLoaded = false;
		hltjet30_id_isLoaded = false;
		hltjet30_tid_isLoaded = false;
		hltl1jet15_id_isLoaded = false;
		hltl1jet15_tid_isLoaded = false;
		hltl1met20_id_isLoaded = false;
		hltl1met20_tid_isLoaded = false;
		hltl1mu_id_isLoaded = false;
		hltl1mu_tid_isLoaded = false;
		hltl2mu9_id_isLoaded = false;
		hltl2mu9_tid_isLoaded = false;
		hltmet25_id_isLoaded = false;
		hltmet25_tid_isLoaded = false;
		hltmu11_id_isLoaded = false;
		hltmu11_tid_isLoaded = false;
		hltmu9_id_isLoaded = false;
		hltmu9_tid_isLoaded = false;
		vtxs_hasRefittedTracks_isLoaded = false;
		vtxs_isFake_isLoaded = false;
		vtxs_isValid_isLoaded = false;
		vtxs_tracksSize_isLoaded = false;
		hyp_jets_idx_isLoaded = false;
		hyp_other_jets_idx_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		evt_ngenjets_isLoaded = false;
		evt_njets_isLoaded = false;
		evt_njpts_isLoaded = false;
		evt_nphotons_isLoaded = false;
		evt_nscs_isLoaded = false;
		evt_ntrkjets_isLoaded = false;
		evt_nvtxs_isLoaded = false;
		hyp_quadlep_bucket_isLoaded = false;
		hyp_quadlep_first_index_isLoaded = false;
		hyp_quadlep_fourth_index_isLoaded = false;
		hyp_quadlep_second_index_isLoaded = false;
		hyp_quadlep_third_index_isLoaded = false;
		hyp_trilep_bucket_isLoaded = false;
		hyp_trilep_first_index_isLoaded = false;
		hyp_trilep_second_index_isLoaded = false;
		hyp_trilep_third_index_isLoaded = false;
		els_pat_flag_isLoaded = false;
		jets_pat_flag_isLoaded = false;
		mus_pat_flag_isLoaded = false;
		evt_scale1fb_isLoaded = false;
		evt_nEvts_isLoaded = false;
		evt_filt_eff_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (evt_dataset_branch != 0) evt_dataset();
	if (evt_HLT_trigNames_branch != 0) evt_HLT_trigNames();
	if (evt_L1_trigNames_branch != 0) evt_L1_trigNames();
	if (genps_pthat_branch != 0) genps_pthat();
	if (vtxs_position_branch != 0) vtxs_position();
	if (evt_bsp4_branch != 0) evt_bsp4();
	if (l1met_p4_branch != 0) l1met_p4();
	if (taus_calo_leadtrk_p4_branch != 0) taus_calo_leadtrk_p4();
	if (taus_calo_p4_branch != 0) taus_calo_p4();
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
	if (genjets_p4_branch != 0) genjets_p4();
	if (genps_lepdaughter_p4_branch != 0) genps_lepdaughter_p4();
	if (genps_p4_branch != 0) genps_p4();
	if (genps_prod_vtx_branch != 0) genps_prod_vtx();
	if (hyp_ll_mc_p4_branch != 0) hyp_ll_mc_p4();
	if (hyp_ll_p4_branch != 0) hyp_ll_p4();
	if (hyp_ll_trk_p4_branch != 0) hyp_ll_trk_p4();
	if (hyp_lt_mc_p4_branch != 0) hyp_lt_mc_p4();
	if (hyp_lt_p4_branch != 0) hyp_lt_p4();
	if (hyp_lt_trk_p4_branch != 0) hyp_lt_trk_p4();
	if (hyp_p4_branch != 0) hyp_p4();
	if (jets_p4_branch != 0) jets_p4();
	if (jpts_p4_branch != 0) jpts_p4();
	if (l1emiso_p4_branch != 0) l1emiso_p4();
	if (l1emnoiso_p4_branch != 0) l1emnoiso_p4();
	if (l1jetsc_p4_branch != 0) l1jetsc_p4();
	if (l1jetsf_p4_branch != 0) l1jetsf_p4();
	if (l1jetst_p4_branch != 0) l1jetst_p4();
	if (l1mus_p4_branch != 0) l1mus_p4();
	if (mus_gfit_outerPos_p4_branch != 0) mus_gfit_outerPos_p4();
	if (mus_p4_branch != 0) mus_p4();
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
	if (taus_pf_lead_chargecand_p4_branch != 0) taus_pf_lead_chargecand_p4();
	if (taus_pf_lead_neutrcand_p4_branch != 0) taus_pf_lead_neutrcand_p4();
	if (taus_pf_leadtrk_p4_branch != 0) taus_pf_leadtrk_p4();
	if (taus_pf_p4_branch != 0) taus_pf_p4();
	if (photons_p4_branch != 0) photons_p4();
	if (scs_p4_branch != 0) scs_p4();
	if (scs_pos_branch != 0) scs_pos();
	if (scs_vtx_branch != 0) scs_vtx();
	if (trks_trk_p4_branch != 0) trks_trk_p4();
	if (trks_vertex_p4_branch != 0) trks_vertex_p4();
	if (hlt2ele10LWR_p4_branch != 0) hlt2ele10LWR_p4();
	if (hlt2mu3_p4_branch != 0) hlt2mu3_p4();
	if (hltLisoele18LWR_p4_branch != 0) hltLisoele18LWR_p4();
	if (hltisoele18R_p4_branch != 0) hltisoele18R_p4();
	if (hltjet30_p4_branch != 0) hltjet30_p4();
	if (hltl1jet15_p4_branch != 0) hltl1jet15_p4();
	if (hltl1met20_p4_branch != 0) hltl1met20_p4();
	if (hltl1mu_p4_branch != 0) hltl1mu_p4();
	if (hltl2mu9_p4_branch != 0) hltl2mu9_p4();
	if (hltmet25_p4_branch != 0) hltmet25_p4();
	if (hltmu11_p4_branch != 0) hltmu11_p4();
	if (hltmu9_p4_branch != 0) hltmu9_p4();
	if (trkjets_p4_branch != 0) trkjets_p4();
	if (hyp_jets_p4_branch != 0) hyp_jets_p4();
	if (hyp_other_jets_p4_branch != 0) hyp_other_jets_p4();
	if (els_inner_positionx_branch != 0) els_inner_positionx();
	if (els_inner_positiony_branch != 0) els_inner_positiony();
	if (els_inner_positionz_branch != 0) els_inner_positionz();
	if (jets_closestElectron_DR_branch != 0) jets_closestElectron_DR();
	if (jets_closestMuon_DR_branch != 0) jets_closestMuon_DR();
	if (evt_bs_dxdz_branch != 0) evt_bs_dxdz();
	if (evt_bs_dxdzErr_branch != 0) evt_bs_dxdzErr();
	if (evt_bs_dydz_branch != 0) evt_bs_dydz();
	if (evt_bs_dydzErr_branch != 0) evt_bs_dydzErr();
	if (evt_bs_sigmaZ_branch != 0) evt_bs_sigmaZ();
	if (evt_bs_sigmaZErr_branch != 0) evt_bs_sigmaZErr();
	if (evt_bs_width_branch != 0) evt_bs_width();
	if (evt_bs_widthErr_branch != 0) evt_bs_widthErr();
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
	if (l1met_etHad_branch != 0) l1met_etHad();
	if (l1met_etTot_branch != 0) l1met_etTot();
	if (l1met_met_branch != 0) l1met_met();
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
	if (evt_pfsumet_branch != 0) evt_pfsumet();
	if (evt_tcmet_branch != 0) evt_tcmet();
	if (evt_tcmetPhi_branch != 0) evt_tcmetPhi();
	if (evt_tcsumet_branch != 0) evt_tcsumet();
	if (jets_combinedSecondaryVertexBJetTag_branch != 0) jets_combinedSecondaryVertexBJetTag();
	if (jets_combinedSecondaryVertexMVABJetTag_branch != 0) jets_combinedSecondaryVertexMVABJetTag();
	if (jets_impactParameterMVABJetTag_branch != 0) jets_impactParameterMVABJetTag();
	if (jets_jetBProbabilityBJetTag_branch != 0) jets_jetBProbabilityBJetTag();
	if (jets_jetProbabilityBJetTag_branch != 0) jets_jetProbabilityBJetTag();
	if (jets_simpleSecondaryVertexBJetTag_branch != 0) jets_simpleSecondaryVertexBJetTag();
	if (jets_softElectronBJetTag_branch != 0) jets_softElectronBJetTag();
	if (jets_softMuonBJetTag_branch != 0) jets_softMuonBJetTag();
	if (jets_softMuonNoIPBJetTag_branch != 0) jets_softMuonNoIPBJetTag();
	if (jets_trackCountingHighEffBJetTag_branch != 0) jets_trackCountingHighEffBJetTag();
	if (jets_trackCountingHighPurBJetTag_branch != 0) jets_trackCountingHighPurBJetTag();
	if (taus_calo_isolationECALhitsEtSum_branch != 0) taus_calo_isolationECALhitsEtSum();
	if (taus_calo_isolationtrksPtSum_branch != 0) taus_calo_isolationtrksPtSum();
	if (taus_calo_leadtrk_HCAL3x3hitsEtSum_branch != 0) taus_calo_leadtrk_HCAL3x3hitsEtSum();
	if (taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch != 0) taus_calo_leadtrk_HCAL3x3hottesthitDEta();
	if (taus_calo_leadtrk_Signed_Sipt_branch != 0) taus_calo_leadtrk_Signed_Sipt();
	if (taus_calo_leadtrk_chi2_branch != 0) taus_calo_leadtrk_chi2();
	if (taus_calo_leadtrk_d0_branch != 0) taus_calo_leadtrk_d0();
	if (taus_calo_leadtrk_lostHits_branch != 0) taus_calo_leadtrk_lostHits();
	if (taus_calo_leadtrk_ndof_branch != 0) taus_calo_leadtrk_ndof();
	if (taus_calo_leadtrk_validHits_branch != 0) taus_calo_leadtrk_validHits();
	if (taus_calo_leadtrk_z0_branch != 0) taus_calo_leadtrk_z0();
	if (taus_calo_maximumHCALhitEt_branch != 0) taus_calo_maximumHCALhitEt();
	if (taus_calo_signaltrksInvariantMass_branch != 0) taus_calo_signaltrksInvariantMass();
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
	if (els_ecalJuraIso_branch != 0) els_ecalJuraIso();
	if (els_ecalJuraTowerIso_branch != 0) els_ecalJuraTowerIso();
	if (els_hcalConeIso_branch != 0) els_hcalConeIso();
	if (els_jetdr_branch != 0) els_jetdr();
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
	if (els_ecalIso_branch != 0) els_ecalIso();
	if (els_egamma_looseId_branch != 0) els_egamma_looseId();
	if (els_egamma_robustHighEnergy_branch != 0) els_egamma_robustHighEnergy();
	if (els_egamma_robustLooseId_branch != 0) els_egamma_robustLooseId();
	if (els_egamma_robustTightId_branch != 0) els_egamma_robustTightId();
	if (els_egamma_tightId_branch != 0) els_egamma_tightId();
	if (els_etaErr_branch != 0) els_etaErr();
	if (els_etaSC_branch != 0) els_etaSC();
	if (els_fBrem_branch != 0) els_fBrem();
	if (els_hOverE_branch != 0) els_hOverE();
	if (els_hcalIso_branch != 0) els_hcalIso();
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
	if (els_tkIso_branch != 0) els_tkIso();
	if (els_vertexphi_branch != 0) els_vertexphi();
	if (els_z0_branch != 0) els_z0();
	if (els_z0Err_branch != 0) els_z0Err();
	if (els_z0corr_branch != 0) els_z0corr();
	if (hyp_Ht_branch != 0) hyp_Ht();
	if (hyp_dPhi_nJet_muCorrMet_branch != 0) hyp_dPhi_nJet_muCorrMet();
	if (hyp_dPhi_nJet_tcMet_branch != 0) hyp_dPhi_nJet_tcMet();
	if (hyp_dPhi_nJet_type1Met_branch != 0) hyp_dPhi_nJet_type1Met();
	if (hyp_dPhi_nJet_unCorrMet_branch != 0) hyp_dPhi_nJet_unCorrMet();
	if (hyp_ll_chi2_branch != 0) hyp_ll_chi2();
	if (hyp_ll_d0_branch != 0) hyp_ll_d0();
	if (hyp_ll_d0Err_branch != 0) hyp_ll_d0Err();
	if (hyp_ll_d0corr_branch != 0) hyp_ll_d0corr();
	if (hyp_ll_dPhi_muCorrMet_branch != 0) hyp_ll_dPhi_muCorrMet();
	if (hyp_ll_dPhi_tcMet_branch != 0) hyp_ll_dPhi_tcMet();
	if (hyp_ll_dPhi_type1Met_branch != 0) hyp_ll_dPhi_type1Met();
	if (hyp_ll_dPhi_unCorrMet_branch != 0) hyp_ll_dPhi_unCorrMet();
	if (hyp_ll_etaErr_branch != 0) hyp_ll_etaErr();
	if (hyp_ll_iso_branch != 0) hyp_ll_iso();
	if (hyp_ll_ndof_branch != 0) hyp_ll_ndof();
	if (hyp_ll_outerEta_branch != 0) hyp_ll_outerEta();
	if (hyp_ll_outerPhi_branch != 0) hyp_ll_outerPhi();
	if (hyp_ll_phiErr_branch != 0) hyp_ll_phiErr();
	if (hyp_ll_ptErr_branch != 0) hyp_ll_ptErr();
	if (hyp_ll_tkIso_branch != 0) hyp_ll_tkIso();
	if (hyp_ll_vertexphi_branch != 0) hyp_ll_vertexphi();
	if (hyp_ll_z0_branch != 0) hyp_ll_z0();
	if (hyp_ll_z0Err_branch != 0) hyp_ll_z0Err();
	if (hyp_ll_z0corr_branch != 0) hyp_ll_z0corr();
	if (hyp_lt_chi2_branch != 0) hyp_lt_chi2();
	if (hyp_lt_d0_branch != 0) hyp_lt_d0();
	if (hyp_lt_d0Err_branch != 0) hyp_lt_d0Err();
	if (hyp_lt_d0corr_branch != 0) hyp_lt_d0corr();
	if (hyp_lt_dPhi_muCorrMet_branch != 0) hyp_lt_dPhi_muCorrMet();
	if (hyp_lt_dPhi_tcMet_branch != 0) hyp_lt_dPhi_tcMet();
	if (hyp_lt_dPhi_type1Met_branch != 0) hyp_lt_dPhi_type1Met();
	if (hyp_lt_dPhi_unCorrMet_branch != 0) hyp_lt_dPhi_unCorrMet();
	if (hyp_lt_etaErr_branch != 0) hyp_lt_etaErr();
	if (hyp_lt_iso_branch != 0) hyp_lt_iso();
	if (hyp_lt_ndof_branch != 0) hyp_lt_ndof();
	if (hyp_lt_outerEta_branch != 0) hyp_lt_outerEta();
	if (hyp_lt_outerPhi_branch != 0) hyp_lt_outerPhi();
	if (hyp_lt_phiErr_branch != 0) hyp_lt_phiErr();
	if (hyp_lt_ptErr_branch != 0) hyp_lt_ptErr();
	if (hyp_lt_tkIso_branch != 0) hyp_lt_tkIso();
	if (hyp_lt_vertexphi_branch != 0) hyp_lt_vertexphi();
	if (hyp_lt_z0_branch != 0) hyp_lt_z0();
	if (hyp_lt_z0Err_branch != 0) hyp_lt_z0Err();
	if (hyp_lt_z0corr_branch != 0) hyp_lt_z0corr();
	if (hyp_sumJetPt_branch != 0) hyp_sumJetPt();
	if (jets_EMFcor_branch != 0) jets_EMFcor();
	if (jets_chFrac_branch != 0) jets_chFrac();
	if (jets_cor_branch != 0) jets_cor();
	if (jets_emFrac_branch != 0) jets_emFrac();
	if (jpts_cor_branch != 0) jpts_cor();
	if (jpts_emFrac_branch != 0) jpts_emFrac();
	if (jpts_jet_cor_branch != 0) jpts_jet_cor();
	if (mus_met_deltax_branch != 0) mus_met_deltax();
	if (mus_met_deltay_branch != 0) mus_met_deltay();
	if (mus_eledr_branch != 0) mus_eledr();
	if (mus_jetdr_branch != 0) mus_jetdr();
	if (mus_trkdr_branch != 0) mus_trkdr();
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
	if (jets_pat_bCorrF_branch != 0) jets_pat_bCorrF();
	if (jets_pat_cCorrF_branch != 0) jets_pat_cCorrF();
	if (jets_pat_combinedSecondaryVertexBJetTag_branch != 0) jets_pat_combinedSecondaryVertexBJetTag();
	if (jets_pat_combinedSecondaryVertexMVABJetTag_branch != 0) jets_pat_combinedSecondaryVertexMVABJetTag();
	if (jets_pat_coneIsolationTauJetTag_branch != 0) jets_pat_coneIsolationTauJetTag();
	if (jets_pat_gluCorrF_branch != 0) jets_pat_gluCorrF();
	if (jets_pat_impactParameterMVABJetTag_branch != 0) jets_pat_impactParameterMVABJetTag();
	if (jets_pat_jetBProbabilityBJetTag_branch != 0) jets_pat_jetBProbabilityBJetTag();
	if (jets_pat_jetCharge_branch != 0) jets_pat_jetCharge();
	if (jets_pat_jetProbabilityBJetTag_branch != 0) jets_pat_jetProbabilityBJetTag();
	if (jets_pat_noCorrF_branch != 0) jets_pat_noCorrF();
	if (jets_pat_simpleSecondaryVertexBJetTag_branch != 0) jets_pat_simpleSecondaryVertexBJetTag();
	if (jets_pat_softElectronBJetTag_branch != 0) jets_pat_softElectronBJetTag();
	if (jets_pat_softMuonBJetTag_branch != 0) jets_pat_softMuonBJetTag();
	if (jets_pat_softMuonNoIPBJetTag_branch != 0) jets_pat_softMuonNoIPBJetTag();
	if (jets_pat_trackCountingHighEffBJetTag_branch != 0) jets_pat_trackCountingHighEffBJetTag();
	if (jets_pat_trackCountingHighPurBJetTag_branch != 0) jets_pat_trackCountingHighPurBJetTag();
	if (jets_pat_udsCorrF_branch != 0) jets_pat_udsCorrF();
	if (mus_pat_caloIso_branch != 0) mus_pat_caloIso();
	if (mus_pat_ecalIso_branch != 0) mus_pat_ecalIso();
	if (mus_pat_ecalvetoDep_branch != 0) mus_pat_ecalvetoDep();
	if (mus_pat_hcalIso_branch != 0) mus_pat_hcalIso();
	if (mus_pat_hcalvetoDep_branch != 0) mus_pat_hcalvetoDep();
	if (mus_pat_trackIso_branch != 0) mus_pat_trackIso();
	if (mus_pat_vetoDep_branch != 0) mus_pat_vetoDep();
	if (taus_pf_bremsRecoveryEOverPLead_branch != 0) taus_pf_bremsRecoveryEOverPLead();
	if (taus_pf_caloComp_branch != 0) taus_pf_caloComp();
	if (taus_pf_ecalStripSumEOverPLead_branch != 0) taus_pf_ecalStripSumEOverPLead();
	if (taus_pf_emf_branch != 0) taus_pf_emf();
	if (taus_pf_hcal3x3OverPLead_branch != 0) taus_pf_hcal3x3OverPLead();
	if (taus_pf_hcalMaxOverPLead_branch != 0) taus_pf_hcalMaxOverPLead();
	if (taus_pf_hcalTotOverPLead_branch != 0) taus_pf_hcalTotOverPLead();
	if (taus_pf_isolationchargecandPtSum_branch != 0) taus_pf_isolationchargecandPtSum();
	if (taus_pf_isolationgammacandEtSum_branch != 0) taus_pf_isolationgammacandEtSum();
	if (taus_pf_lead_chargecand_Signed_Sipt_branch != 0) taus_pf_lead_chargecand_Signed_Sipt();
	if (taus_pf_leadtrk_chi2_branch != 0) taus_pf_leadtrk_chi2();
	if (taus_pf_leadtrk_d0_branch != 0) taus_pf_leadtrk_d0();
	if (taus_pf_leadtrk_lostHits_branch != 0) taus_pf_leadtrk_lostHits();
	if (taus_pf_leadtrk_ndof_branch != 0) taus_pf_leadtrk_ndof();
	if (taus_pf_leadtrk_validHits_branch != 0) taus_pf_leadtrk_validHits();
	if (taus_pf_leadtrk_z0_branch != 0) taus_pf_leadtrk_z0();
	if (taus_pf_maximumHCALPFClusterEt_branch != 0) taus_pf_maximumHCALPFClusterEt();
	if (taus_pf_segComp_branch != 0) taus_pf_segComp();
	if (photons_e1x5_branch != 0) photons_e1x5();
	if (photons_e2x5Max_branch != 0) photons_e2x5Max();
	if (photons_e3x3_branch != 0) photons_e3x3();
	if (photons_e5x5_branch != 0) photons_e5x5();
	if (photons_eMax_branch != 0) photons_eMax();
	if (photons_eSC_branch != 0) photons_eSC();
	if (photons_eSCPresh_branch != 0) photons_eSCPresh();
	if (photons_eSCRaw_branch != 0) photons_eSCRaw();
	if (photons_eSeed_branch != 0) photons_eSeed();
	if (photons_ecalIso_branch != 0) photons_ecalIso();
	if (photons_hOverE_branch != 0) photons_hOverE();
	if (photons_hcalIso_branch != 0) photons_hcalIso();
	if (photons_sigmaEtaEta_branch != 0) photons_sigmaEtaEta();
	if (photons_sigmaIEtaIEta_branch != 0) photons_sigmaIEtaIEta();
	if (photons_sigmaIPhiIPhi_branch != 0) photons_sigmaIPhiIPhi();
	if (photons_sigmaPhiPhi_branch != 0) photons_sigmaPhiPhi();
	if (photons_tkIso_branch != 0) photons_tkIso();
	if (scs_clustersSize_branch != 0) scs_clustersSize();
	if (scs_crystalsSize_branch != 0) scs_crystalsSize();
	if (scs_e1x3_branch != 0) scs_e1x3();
	if (scs_e1x5_branch != 0) scs_e1x5();
	if (scs_e2x2_branch != 0) scs_e2x2();
	if (scs_e2x5Max_branch != 0) scs_e2x5Max();
	if (scs_e3x1_branch != 0) scs_e3x1();
	if (scs_e3x2_branch != 0) scs_e3x2();
	if (scs_e3x3_branch != 0) scs_e3x3();
	if (scs_e4x4_branch != 0) scs_e4x4();
	if (scs_e5x5_branch != 0) scs_e5x5();
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
	if (mus_tcmet_deltax_branch != 0) mus_tcmet_deltax();
	if (mus_tcmet_deltay_branch != 0) mus_tcmet_deltay();
	if (trks_chi2_branch != 0) trks_chi2();
	if (trks_d0_branch != 0) trks_d0();
	if (trks_d0Err_branch != 0) trks_d0Err();
	if (trks_d0corr_branch != 0) trks_d0corr();
	if (trks_d0corrPhi_branch != 0) trks_d0corrPhi();
	if (trks_etaErr_branch != 0) trks_etaErr();
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
	if (evt_nl1emiso_branch != 0) evt_nl1emiso();
	if (evt_nl1emnoiso_branch != 0) evt_nl1emnoiso();
	if (evt_nl1jetsc_branch != 0) evt_nl1jetsc();
	if (evt_nl1jetsf_branch != 0) evt_nl1jetsf();
	if (evt_nl1jetst_branch != 0) evt_nl1jetst();
	if (evt_nl1mus_branch != 0) evt_nl1mus();
	if (pdfinfo_id1_branch != 0) pdfinfo_id1();
	if (pdfinfo_id2_branch != 0) pdfinfo_id2();
	if (taus_calo_charge_branch != 0) taus_calo_charge();
	if (taus_calo_iso_ntrks_branch != 0) taus_calo_iso_ntrks();
	if (taus_calo_sig_ntrks_branch != 0) taus_calo_sig_ntrks();
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
	if (els_closestJet_branch != 0) els_closestJet();
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
	if (hyp_ll_charge_branch != 0) hyp_ll_charge();
	if (hyp_ll_id_branch != 0) hyp_ll_id();
	if (hyp_ll_index_branch != 0) hyp_ll_index();
	if (hyp_ll_lostHits_branch != 0) hyp_ll_lostHits();
	if (hyp_ll_mc_id_branch != 0) hyp_ll_mc_id();
	if (hyp_ll_mc_motherid_branch != 0) hyp_ll_mc_motherid();
	if (hyp_ll_validHits_branch != 0) hyp_ll_validHits();
	if (hyp_lt_charge_branch != 0) hyp_lt_charge();
	if (hyp_lt_id_branch != 0) hyp_lt_id();
	if (hyp_lt_index_branch != 0) hyp_lt_index();
	if (hyp_lt_lostHits_branch != 0) hyp_lt_lostHits();
	if (hyp_lt_mc_id_branch != 0) hyp_lt_mc_id();
	if (hyp_lt_mc_motherid_branch != 0) hyp_lt_mc_motherid();
	if (hyp_lt_validHits_branch != 0) hyp_lt_validHits();
	if (hyp_njets_branch != 0) hyp_njets();
	if (hyp_nojets_branch != 0) hyp_nojets();
	if (hyp_type_branch != 0) hyp_type();
	if (hyp_quadlep_first_type_branch != 0) hyp_quadlep_first_type();
	if (hyp_quadlep_fourth_type_branch != 0) hyp_quadlep_fourth_type();
	if (hyp_quadlep_second_type_branch != 0) hyp_quadlep_second_type();
	if (hyp_quadlep_third_type_branch != 0) hyp_quadlep_third_type();
	if (hyp_trilep_first_type_branch != 0) hyp_trilep_first_type();
	if (hyp_trilep_second_type_branch != 0) hyp_trilep_second_type();
	if (hyp_trilep_third_type_branch != 0) hyp_trilep_third_type();
	if (jets_closestElectron_branch != 0) jets_closestElectron();
	if (jets_closestMuon_branch != 0) jets_closestMuon();
	if (l1emiso_ieta_branch != 0) l1emiso_ieta();
	if (l1emiso_iphi_branch != 0) l1emiso_iphi();
	if (l1emiso_rawId_branch != 0) l1emiso_rawId();
	if (l1emiso_type_branch != 0) l1emiso_type();
	if (l1emnoiso_ieta_branch != 0) l1emnoiso_ieta();
	if (l1emnoiso_iphi_branch != 0) l1emnoiso_iphi();
	if (l1emnoiso_rawId_branch != 0) l1emnoiso_rawId();
	if (l1emnoiso_type_branch != 0) l1emnoiso_type();
	if (l1jetsc_ieta_branch != 0) l1jetsc_ieta();
	if (l1jetsc_iphi_branch != 0) l1jetsc_iphi();
	if (l1jetsc_rawId_branch != 0) l1jetsc_rawId();
	if (l1jetsc_type_branch != 0) l1jetsc_type();
	if (l1jetsf_ieta_branch != 0) l1jetsf_ieta();
	if (l1jetsf_iphi_branch != 0) l1jetsf_iphi();
	if (l1jetsf_rawId_branch != 0) l1jetsf_rawId();
	if (l1jetsf_type_branch != 0) l1jetsf_type();
	if (l1jetst_ieta_branch != 0) l1jetst_ieta();
	if (l1jetst_iphi_branch != 0) l1jetst_iphi();
	if (l1jetst_rawId_branch != 0) l1jetst_rawId();
	if (l1jetst_type_branch != 0) l1jetst_type();
	if (l1mus_flags_branch != 0) l1mus_flags();
	if (l1mus_q_branch != 0) l1mus_q();
	if (l1mus_qual_branch != 0) l1mus_qual();
	if (l1mus_qualFlags_branch != 0) l1mus_qualFlags();
	if (mus_met_flag_branch != 0) mus_met_flag();
	if (mus_closestEle_branch != 0) mus_closestEle();
	if (mus_closestJet_branch != 0) mus_closestJet();
	if (mus_trkidx_branch != 0) mus_trkidx();
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
	if (els_pat_genID_branch != 0) els_pat_genID();
	if (els_pat_genMotherID_branch != 0) els_pat_genMotherID();
	if (jets_pat_genPartonMother_id_branch != 0) jets_pat_genPartonMother_id();
	if (jets_pat_genParton_id_branch != 0) jets_pat_genParton_id();
	if (jets_pat_partonFlavour_branch != 0) jets_pat_partonFlavour();
	if (mus_pat_genID_branch != 0) mus_pat_genID();
	if (mus_pat_genMotherID_branch != 0) mus_pat_genMotherID();
	if (taus_pf_charge_branch != 0) taus_pf_charge();
	if (taus_pf_electronPreID_branch != 0) taus_pf_electronPreID();
	if (taus_pf_hasMuonReference_branch != 0) taus_pf_hasMuonReference();
	if (taus_pf_iso_ncharge_cand_branch != 0) taus_pf_iso_ncharge_cand();
	if (taus_pf_iso_ngamma_cand_branch != 0) taus_pf_iso_ngamma_cand();
	if (taus_pf_iso_nneutr_cand_branch != 0) taus_pf_iso_nneutr_cand();
	if (taus_pf_muonPreID_branch != 0) taus_pf_muonPreID();
	if (taus_pf_sig_ncharge_cand_branch != 0) taus_pf_sig_ncharge_cand();
	if (taus_pf_sig_ngamma_cand_branch != 0) taus_pf_sig_ngamma_cand();
	if (taus_pf_sig_nneutr_cand_branch != 0) taus_pf_sig_nneutr_cand();
	if (mus_tcmet_flag_branch != 0) mus_tcmet_flag();
	if (trks_charge_branch != 0) trks_charge();
	if (trks_lostHits_branch != 0) trks_lostHits();
	if (trks_validHits_branch != 0) trks_validHits();
	if (trks_elsidx_branch != 0) trks_elsidx();
	if (trk_musidx_branch != 0) trk_musidx();
	if (hlt2ele10LWR_id_branch != 0) hlt2ele10LWR_id();
	if (hlt2ele10LWR_tid_branch != 0) hlt2ele10LWR_tid();
	if (hlt2mu3_id_branch != 0) hlt2mu3_id();
	if (hlt2mu3_tid_branch != 0) hlt2mu3_tid();
	if (hltLisoele18LWR_id_branch != 0) hltLisoele18LWR_id();
	if (hltLisoele18LWR_tid_branch != 0) hltLisoele18LWR_tid();
	if (hltisoele18R_id_branch != 0) hltisoele18R_id();
	if (hltisoele18R_tid_branch != 0) hltisoele18R_tid();
	if (hltjet30_id_branch != 0) hltjet30_id();
	if (hltjet30_tid_branch != 0) hltjet30_tid();
	if (hltl1jet15_id_branch != 0) hltl1jet15_id();
	if (hltl1jet15_tid_branch != 0) hltl1jet15_tid();
	if (hltl1met20_id_branch != 0) hltl1met20_id();
	if (hltl1met20_tid_branch != 0) hltl1met20_tid();
	if (hltl1mu_id_branch != 0) hltl1mu_id();
	if (hltl1mu_tid_branch != 0) hltl1mu_tid();
	if (hltl2mu9_id_branch != 0) hltl2mu9_id();
	if (hltl2mu9_tid_branch != 0) hltl2mu9_tid();
	if (hltmet25_id_branch != 0) hltmet25_id();
	if (hltmet25_tid_branch != 0) hltmet25_tid();
	if (hltmu11_id_branch != 0) hltmu11_id();
	if (hltmu11_tid_branch != 0) hltmu11_tid();
	if (hltmu9_id_branch != 0) hltmu9_id();
	if (hltmu9_tid_branch != 0) hltmu9_tid();
	if (vtxs_hasRefittedTracks_branch != 0) vtxs_hasRefittedTracks();
	if (vtxs_isFake_branch != 0) vtxs_isFake();
	if (vtxs_isValid_branch != 0) vtxs_isValid();
	if (vtxs_tracksSize_branch != 0) vtxs_tracksSize();
	if (hyp_jets_idx_branch != 0) hyp_jets_idx();
	if (hyp_other_jets_idx_branch != 0) hyp_other_jets_idx();
	if (evt_nels_branch != 0) evt_nels();
	if (evt_event_branch != 0) evt_event();
	if (evt_lumiBlock_branch != 0) evt_lumiBlock();
	if (evt_run_branch != 0) evt_run();
	if (evt_ngenjets_branch != 0) evt_ngenjets();
	if (evt_njets_branch != 0) evt_njets();
	if (evt_njpts_branch != 0) evt_njpts();
	if (evt_nphotons_branch != 0) evt_nphotons();
	if (evt_nscs_branch != 0) evt_nscs();
	if (evt_ntrkjets_branch != 0) evt_ntrkjets();
	if (evt_nvtxs_branch != 0) evt_nvtxs();
	if (hyp_quadlep_bucket_branch != 0) hyp_quadlep_bucket();
	if (hyp_quadlep_first_index_branch != 0) hyp_quadlep_first_index();
	if (hyp_quadlep_fourth_index_branch != 0) hyp_quadlep_fourth_index();
	if (hyp_quadlep_second_index_branch != 0) hyp_quadlep_second_index();
	if (hyp_quadlep_third_index_branch != 0) hyp_quadlep_third_index();
	if (hyp_trilep_bucket_branch != 0) hyp_trilep_bucket();
	if (hyp_trilep_first_index_branch != 0) hyp_trilep_first_index();
	if (hyp_trilep_second_index_branch != 0) hyp_trilep_second_index();
	if (hyp_trilep_third_index_branch != 0) hyp_trilep_third_index();
	if (els_pat_flag_branch != 0) els_pat_flag();
	if (jets_pat_flag_branch != 0) jets_pat_flag();
	if (mus_pat_flag_branch != 0) mus_pat_flag();
	if (evt_scale1fb_branch != 0) evt_scale1fb();
	if (evt_nEvts_branch != 0) evt_nEvts();
	if (evt_filt_eff_branch != 0) evt_filt_eff();
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
	vector<TString> &evt_HLT_trigNames()
	{
		if (not evt_HLT_trigNames_isLoaded) {
			if (evt_HLT_trigNames_branch != 0) {
				evt_HLT_trigNames_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				int e;
				frexp(evt_bsp4_.pt(), &e);
				if (not isfinite(evt_bsp4_.pt()) || e > 30) {
					printf("branch evt_bsp4_branch contains a bad float: %f\n", evt_bsp4_.pt());
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bsp4_branch does not exist!\n");
				exit(1);
			}
			evt_bsp4_isLoaded = true;
		}
		return evt_bsp4_;
	}
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &l1met_p4()
	{
		if (not l1met_p4_isLoaded) {
			if (l1met_p4_branch != 0) {
				l1met_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				int e;
				frexp(l1met_p4_.pt(), &e);
				if (not isfinite(l1met_p4_.pt()) || e > 30) {
					printf("branch l1met_p4_branch contains a bad float: %f\n", l1met_p4_.pt());
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1met_p4_branch does not exist!\n");
				exit(1);
			}
			l1met_p4_isLoaded = true;
		}
		return l1met_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_calo_leadtrk_p4()
	{
		if (not taus_calo_leadtrk_p4_isLoaded) {
			if (taus_calo_leadtrk_p4_branch != 0) {
				taus_calo_leadtrk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_calo_leadtrk_p4_.begin(); i != taus_calo_leadtrk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_calo_leadtrk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_p4_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_p4_isLoaded = true;
		}
		return taus_calo_leadtrk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_calo_p4()
	{
		if (not taus_calo_p4_isLoaded) {
			if (taus_calo_p4_branch != 0) {
				taus_calo_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_calo_p4_.begin(); i != taus_calo_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_calo_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_p4_branch does not exist!\n");
				exit(1);
			}
			taus_calo_p4_isLoaded = true;
		}
		return taus_calo_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_motherp4()
	{
		if (not els_mc3_motherp4_isLoaded) {
			if (els_mc3_motherp4_branch != 0) {
				els_mc3_motherp4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_mc3_motherp4_.begin(); i != els_mc3_motherp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_mc3_motherp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_mc3_p4_.begin(); i != els_mc3_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_mc3_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_mc_motherp4_.begin(); i != els_mc_motherp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_mc_motherp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_mc_p4_.begin(); i != els_mc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_mc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_mc_gp_p4_.begin(); i != jets_mc_gp_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_mc_gp_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_mc_p4_.begin(); i != jets_mc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_mc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_mc3_motherp4_.begin(); i != mus_mc3_motherp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_mc3_motherp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_mc3_p4_.begin(); i != mus_mc3_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_mc3_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_mc_motherp4_.begin(); i != mus_mc_motherp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_mc_motherp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_mc_p4_.begin(); i != mus_mc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_mc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = trk_mc3p4_.begin(); i != trk_mc3p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trk_mc3p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = trk_mcp4_.begin(); i != trk_mcp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trk_mcp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_p4_.begin(); i != els_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_p4In_.begin(); i != els_p4In_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_p4In_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_p4Out_.begin(); i != els_p4Out_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_p4Out_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_trk_p4_.begin(); i != els_trk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_trk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_vertex_p4_.begin(); i != els_vertex_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_vertex_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			els_vertex_p4_isLoaded = true;
		}
		return els_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genjets_p4()
	{
		if (not genjets_p4_isLoaded) {
			if (genjets_p4_branch != 0) {
				genjets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = genjets_p4_.begin(); i != genjets_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch genjets_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genjets_p4_branch does not exist!\n");
				exit(1);
			}
			genjets_p4_isLoaded = true;
		}
		return genjets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_lepdaughter_p4()
	{
		if (not genps_lepdaughter_p4_isLoaded) {
			if (genps_lepdaughter_p4_branch != 0) {
				genps_lepdaughter_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = genps_lepdaughter_p4_.begin(); i != genps_lepdaughter_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch genps_lepdaughter_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = genps_p4_.begin(); i != genps_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch genps_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = genps_prod_vtx_.begin(); i != genps_prod_vtx_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch genps_prod_vtx_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genps_prod_vtx_branch does not exist!\n");
				exit(1);
			}
			genps_prod_vtx_isLoaded = true;
		}
		return genps_prod_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_mc_p4()
	{
		if (not hyp_ll_mc_p4_isLoaded) {
			if (hyp_ll_mc_p4_branch != 0) {
				hyp_ll_mc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_ll_mc_p4_.begin(); i != hyp_ll_mc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_ll_mc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_mc_p4_isLoaded = true;
		}
		return hyp_ll_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_p4()
	{
		if (not hyp_ll_p4_isLoaded) {
			if (hyp_ll_p4_branch != 0) {
				hyp_ll_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_ll_p4_.begin(); i != hyp_ll_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_ll_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_p4_isLoaded = true;
		}
		return hyp_ll_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_trk_p4()
	{
		if (not hyp_ll_trk_p4_isLoaded) {
			if (hyp_ll_trk_p4_branch != 0) {
				hyp_ll_trk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_ll_trk_p4_.begin(); i != hyp_ll_trk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_ll_trk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_trk_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_trk_p4_isLoaded = true;
		}
		return hyp_ll_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_mc_p4()
	{
		if (not hyp_lt_mc_p4_isLoaded) {
			if (hyp_lt_mc_p4_branch != 0) {
				hyp_lt_mc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_lt_mc_p4_.begin(); i != hyp_lt_mc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_lt_mc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_mc_p4_isLoaded = true;
		}
		return hyp_lt_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_p4()
	{
		if (not hyp_lt_p4_isLoaded) {
			if (hyp_lt_p4_branch != 0) {
				hyp_lt_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_lt_p4_.begin(); i != hyp_lt_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_lt_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_p4_isLoaded = true;
		}
		return hyp_lt_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_trk_p4()
	{
		if (not hyp_lt_trk_p4_isLoaded) {
			if (hyp_lt_trk_p4_branch != 0) {
				hyp_lt_trk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_lt_trk_p4_.begin(); i != hyp_lt_trk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_lt_trk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_trk_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_trk_p4_isLoaded = true;
		}
		return hyp_lt_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_p4()
	{
		if (not hyp_p4_isLoaded) {
			if (hyp_p4_branch != 0) {
				hyp_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hyp_p4_.begin(); i != hyp_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hyp_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_p4_isLoaded = true;
		}
		return hyp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_p4()
	{
		if (not jets_p4_isLoaded) {
			if (jets_p4_branch != 0) {
				jets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_p4_.begin(); i != jets_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_p4_branch does not exist!\n");
				exit(1);
			}
			jets_p4_isLoaded = true;
		}
		return jets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jpts_p4()
	{
		if (not jpts_p4_isLoaded) {
			if (jpts_p4_branch != 0) {
				jpts_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jpts_p4_.begin(); i != jpts_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jpts_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jpts_p4_branch does not exist!\n");
				exit(1);
			}
			jpts_p4_isLoaded = true;
		}
		return jpts_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emiso_p4()
	{
		if (not l1emiso_p4_isLoaded) {
			if (l1emiso_p4_branch != 0) {
				l1emiso_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1emiso_p4_.begin(); i != l1emiso_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1emiso_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emiso_p4_branch does not exist!\n");
				exit(1);
			}
			l1emiso_p4_isLoaded = true;
		}
		return l1emiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emnoiso_p4()
	{
		if (not l1emnoiso_p4_isLoaded) {
			if (l1emnoiso_p4_branch != 0) {
				l1emnoiso_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1emnoiso_p4_.begin(); i != l1emnoiso_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1emnoiso_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emnoiso_p4_branch does not exist!\n");
				exit(1);
			}
			l1emnoiso_p4_isLoaded = true;
		}
		return l1emnoiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsc_p4()
	{
		if (not l1jetsc_p4_isLoaded) {
			if (l1jetsc_p4_branch != 0) {
				l1jetsc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1jetsc_p4_.begin(); i != l1jetsc_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1jetsc_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsc_p4_branch does not exist!\n");
				exit(1);
			}
			l1jetsc_p4_isLoaded = true;
		}
		return l1jetsc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsf_p4()
	{
		if (not l1jetsf_p4_isLoaded) {
			if (l1jetsf_p4_branch != 0) {
				l1jetsf_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1jetsf_p4_.begin(); i != l1jetsf_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1jetsf_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsf_p4_branch does not exist!\n");
				exit(1);
			}
			l1jetsf_p4_isLoaded = true;
		}
		return l1jetsf_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetst_p4()
	{
		if (not l1jetst_p4_isLoaded) {
			if (l1jetst_p4_branch != 0) {
				l1jetst_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1jetst_p4_.begin(); i != l1jetst_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1jetst_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetst_p4_branch does not exist!\n");
				exit(1);
			}
			l1jetst_p4_isLoaded = true;
		}
		return l1jetst_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1mus_p4()
	{
		if (not l1mus_p4_isLoaded) {
			if (l1mus_p4_branch != 0) {
				l1mus_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = l1mus_p4_.begin(); i != l1mus_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch l1mus_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1mus_p4_branch does not exist!\n");
				exit(1);
			}
			l1mus_p4_isLoaded = true;
		}
		return l1mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_gfit_outerPos_p4()
	{
		if (not mus_gfit_outerPos_p4_isLoaded) {
			if (mus_gfit_outerPos_p4_branch != 0) {
				mus_gfit_outerPos_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_gfit_outerPos_p4_.begin(); i != mus_gfit_outerPos_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_gfit_outerPos_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_p4_.begin(); i != mus_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_trk_p4_.begin(); i != mus_trk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_trk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_vertex_p4_.begin(); i != mus_vertex_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_vertex_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_vertex_p4_branch does not exist!\n");
				exit(1);
			}
			mus_vertex_p4_isLoaded = true;
		}
		return mus_vertex_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genMotherP4()
	{
		if (not els_pat_genMotherP4_isLoaded) {
			if (els_pat_genMotherP4_branch != 0) {
				els_pat_genMotherP4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_pat_genMotherP4_.begin(); i != els_pat_genMotherP4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_pat_genMotherP4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_pat_genMotherP4_branch does not exist!\n");
				exit(1);
			}
			els_pat_genMotherP4_isLoaded = true;
		}
		return els_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genP4()
	{
		if (not els_pat_genP4_isLoaded) {
			if (els_pat_genP4_branch != 0) {
				els_pat_genP4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_pat_genP4_.begin(); i != els_pat_genP4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_pat_genP4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_pat_genP4_branch does not exist!\n");
				exit(1);
			}
			els_pat_genP4_isLoaded = true;
		}
		return els_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_p4()
	{
		if (not els_pat_p4_isLoaded) {
			if (els_pat_p4_branch != 0) {
				els_pat_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = els_pat_p4_.begin(); i != els_pat_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch els_pat_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_pat_p4_branch does not exist!\n");
				exit(1);
			}
			els_pat_p4_isLoaded = true;
		}
		return els_pat_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genJet_p4()
	{
		if (not jets_pat_genJet_p4_isLoaded) {
			if (jets_pat_genJet_p4_branch != 0) {
				jets_pat_genJet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_pat_genJet_p4_.begin(); i != jets_pat_genJet_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_pat_genJet_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_genJet_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genJet_p4_isLoaded = true;
		}
		return jets_pat_genJet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genPartonMother_p4()
	{
		if (not jets_pat_genPartonMother_p4_isLoaded) {
			if (jets_pat_genPartonMother_p4_branch != 0) {
				jets_pat_genPartonMother_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_pat_genPartonMother_p4_.begin(); i != jets_pat_genPartonMother_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_pat_genPartonMother_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_genPartonMother_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return jets_pat_genPartonMother_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genParton_p4()
	{
		if (not jets_pat_genParton_p4_isLoaded) {
			if (jets_pat_genParton_p4_branch != 0) {
				jets_pat_genParton_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_pat_genParton_p4_.begin(); i != jets_pat_genParton_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_pat_genParton_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_genParton_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genParton_p4_isLoaded = true;
		}
		return jets_pat_genParton_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_p4()
	{
		if (not jets_pat_jet_p4_isLoaded) {
			if (jets_pat_jet_p4_branch != 0) {
				jets_pat_jet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_pat_jet_p4_.begin(); i != jets_pat_jet_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_pat_jet_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_jet_p4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jet_p4_isLoaded = true;
		}
		return jets_pat_jet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_uncorp4()
	{
		if (not jets_pat_jet_uncorp4_isLoaded) {
			if (jets_pat_jet_uncorp4_branch != 0) {
				jets_pat_jet_uncorp4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = jets_pat_jet_uncorp4_.begin(); i != jets_pat_jet_uncorp4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch jets_pat_jet_uncorp4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_jet_uncorp4_branch does not exist!\n");
				exit(1);
			}
			jets_pat_jet_uncorp4_isLoaded = true;
		}
		return jets_pat_jet_uncorp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genMotherP4()
	{
		if (not mus_pat_genMotherP4_isLoaded) {
			if (mus_pat_genMotherP4_branch != 0) {
				mus_pat_genMotherP4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_pat_genMotherP4_.begin(); i != mus_pat_genMotherP4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_pat_genMotherP4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_genMotherP4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genMotherP4_isLoaded = true;
		}
		return mus_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genP4()
	{
		if (not mus_pat_genP4_isLoaded) {
			if (mus_pat_genP4_branch != 0) {
				mus_pat_genP4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_pat_genP4_.begin(); i != mus_pat_genP4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_pat_genP4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_genP4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genP4_isLoaded = true;
		}
		return mus_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_p4()
	{
		if (not mus_pat_p4_isLoaded) {
			if (mus_pat_p4_branch != 0) {
				mus_pat_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = mus_pat_p4_.begin(); i != mus_pat_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch mus_pat_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_p4_branch does not exist!\n");
				exit(1);
			}
			mus_pat_p4_isLoaded = true;
		}
		return mus_pat_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_lead_chargecand_p4()
	{
		if (not taus_pf_lead_chargecand_p4_isLoaded) {
			if (taus_pf_lead_chargecand_p4_branch != 0) {
				taus_pf_lead_chargecand_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_pf_lead_chargecand_p4_.begin(); i != taus_pf_lead_chargecand_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_pf_lead_chargecand_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_lead_chargecand_p4_branch does not exist!\n");
				exit(1);
			}
			taus_pf_lead_chargecand_p4_isLoaded = true;
		}
		return taus_pf_lead_chargecand_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_lead_neutrcand_p4()
	{
		if (not taus_pf_lead_neutrcand_p4_isLoaded) {
			if (taus_pf_lead_neutrcand_p4_branch != 0) {
				taus_pf_lead_neutrcand_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_pf_lead_neutrcand_p4_.begin(); i != taus_pf_lead_neutrcand_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_pf_lead_neutrcand_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_lead_neutrcand_p4_branch does not exist!\n");
				exit(1);
			}
			taus_pf_lead_neutrcand_p4_isLoaded = true;
		}
		return taus_pf_lead_neutrcand_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_leadtrk_p4()
	{
		if (not taus_pf_leadtrk_p4_isLoaded) {
			if (taus_pf_leadtrk_p4_branch != 0) {
				taus_pf_leadtrk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_pf_leadtrk_p4_.begin(); i != taus_pf_leadtrk_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_pf_leadtrk_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_p4_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_p4_isLoaded = true;
		}
		return taus_pf_leadtrk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_p4()
	{
		if (not taus_pf_p4_isLoaded) {
			if (taus_pf_p4_branch != 0) {
				taus_pf_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = taus_pf_p4_.begin(); i != taus_pf_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch taus_pf_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_p4_branch does not exist!\n");
				exit(1);
			}
			taus_pf_p4_isLoaded = true;
		}
		return taus_pf_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &photons_p4()
	{
		if (not photons_p4_isLoaded) {
			if (photons_p4_branch != 0) {
				photons_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = photons_p4_.begin(); i != photons_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch photons_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_p4_branch does not exist!\n");
				exit(1);
			}
			photons_p4_isLoaded = true;
		}
		return photons_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_p4()
	{
		if (not scs_p4_isLoaded) {
			if (scs_p4_branch != 0) {
				scs_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = scs_p4_.begin(); i != scs_p4_.end(); ++i) {
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_pos()
	{
		if (not scs_pos_isLoaded) {
			if (scs_pos_branch != 0) {
				scs_pos_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = scs_pos_.begin(); i != scs_pos_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch scs_pos_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_pos_branch does not exist!\n");
				exit(1);
			}
			scs_pos_isLoaded = true;
		}
		return scs_pos_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_vtx()
	{
		if (not scs_vtx_isLoaded) {
			if (scs_vtx_branch != 0) {
				scs_vtx_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = scs_vtx_.begin(); i != scs_vtx_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch scs_vtx_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_vtx_branch does not exist!\n");
				exit(1);
			}
			scs_vtx_isLoaded = true;
		}
		return scs_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4()
	{
		if (not trks_trk_p4_isLoaded) {
			if (trks_trk_p4_branch != 0) {
				trks_trk_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = trks_trk_p4_.begin(); i != trks_trk_p4_.end(); ++i) {
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_vertex_p4()
	{
		if (not trks_vertex_p4_isLoaded) {
			if (trks_vertex_p4_branch != 0) {
				trks_vertex_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = trks_vertex_p4_.begin(); i != trks_vertex_p4_.end(); ++i) {
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hlt2ele10LWR_p4()
	{
		if (not hlt2ele10LWR_p4_isLoaded) {
			if (hlt2ele10LWR_p4_branch != 0) {
				hlt2ele10LWR_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hlt2ele10LWR_p4_.begin(); i != hlt2ele10LWR_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hlt2ele10LWR_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2ele10LWR_p4_branch does not exist!\n");
				exit(1);
			}
			hlt2ele10LWR_p4_isLoaded = true;
		}
		return hlt2ele10LWR_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hlt2mu3_p4()
	{
		if (not hlt2mu3_p4_isLoaded) {
			if (hlt2mu3_p4_branch != 0) {
				hlt2mu3_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hlt2mu3_p4_.begin(); i != hlt2mu3_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hlt2mu3_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2mu3_p4_branch does not exist!\n");
				exit(1);
			}
			hlt2mu3_p4_isLoaded = true;
		}
		return hlt2mu3_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltLisoele18LWR_p4()
	{
		if (not hltLisoele18LWR_p4_isLoaded) {
			if (hltLisoele18LWR_p4_branch != 0) {
				hltLisoele18LWR_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltLisoele18LWR_p4_.begin(); i != hltLisoele18LWR_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltLisoele18LWR_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltLisoele18LWR_p4_branch does not exist!\n");
				exit(1);
			}
			hltLisoele18LWR_p4_isLoaded = true;
		}
		return hltLisoele18LWR_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltisoele18R_p4()
	{
		if (not hltisoele18R_p4_isLoaded) {
			if (hltisoele18R_p4_branch != 0) {
				hltisoele18R_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltisoele18R_p4_.begin(); i != hltisoele18R_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltisoele18R_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltisoele18R_p4_branch does not exist!\n");
				exit(1);
			}
			hltisoele18R_p4_isLoaded = true;
		}
		return hltisoele18R_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltjet30_p4()
	{
		if (not hltjet30_p4_isLoaded) {
			if (hltjet30_p4_branch != 0) {
				hltjet30_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltjet30_p4_.begin(); i != hltjet30_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltjet30_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltjet30_p4_branch does not exist!\n");
				exit(1);
			}
			hltjet30_p4_isLoaded = true;
		}
		return hltjet30_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1jet15_p4()
	{
		if (not hltl1jet15_p4_isLoaded) {
			if (hltl1jet15_p4_branch != 0) {
				hltl1jet15_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltl1jet15_p4_.begin(); i != hltl1jet15_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltl1jet15_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1jet15_p4_branch does not exist!\n");
				exit(1);
			}
			hltl1jet15_p4_isLoaded = true;
		}
		return hltl1jet15_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1met20_p4()
	{
		if (not hltl1met20_p4_isLoaded) {
			if (hltl1met20_p4_branch != 0) {
				hltl1met20_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltl1met20_p4_.begin(); i != hltl1met20_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltl1met20_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1met20_p4_branch does not exist!\n");
				exit(1);
			}
			hltl1met20_p4_isLoaded = true;
		}
		return hltl1met20_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1mu_p4()
	{
		if (not hltl1mu_p4_isLoaded) {
			if (hltl1mu_p4_branch != 0) {
				hltl1mu_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltl1mu_p4_.begin(); i != hltl1mu_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltl1mu_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1mu_p4_branch does not exist!\n");
				exit(1);
			}
			hltl1mu_p4_isLoaded = true;
		}
		return hltl1mu_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl2mu9_p4()
	{
		if (not hltl2mu9_p4_isLoaded) {
			if (hltl2mu9_p4_branch != 0) {
				hltl2mu9_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltl2mu9_p4_.begin(); i != hltl2mu9_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltl2mu9_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl2mu9_p4_branch does not exist!\n");
				exit(1);
			}
			hltl2mu9_p4_isLoaded = true;
		}
		return hltl2mu9_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmet25_p4()
	{
		if (not hltmet25_p4_isLoaded) {
			if (hltmet25_p4_branch != 0) {
				hltmet25_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltmet25_p4_.begin(); i != hltmet25_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltmet25_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmet25_p4_branch does not exist!\n");
				exit(1);
			}
			hltmet25_p4_isLoaded = true;
		}
		return hltmet25_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmu11_p4()
	{
		if (not hltmu11_p4_isLoaded) {
			if (hltmu11_p4_branch != 0) {
				hltmu11_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltmu11_p4_.begin(); i != hltmu11_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltmu11_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu11_p4_branch does not exist!\n");
				exit(1);
			}
			hltmu11_p4_isLoaded = true;
		}
		return hltmu11_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmu9_p4()
	{
		if (not hltmu9_p4_isLoaded) {
			if (hltmu9_p4_branch != 0) {
				hltmu9_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = hltmu9_p4_.begin(); i != hltmu9_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch hltmu9_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu9_p4_branch does not exist!\n");
				exit(1);
			}
			hltmu9_p4_isLoaded = true;
		}
		return hltmu9_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trkjets_p4()
	{
		if (not trkjets_p4_isLoaded) {
			if (trkjets_p4_branch != 0) {
				trkjets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = trkjets_p4_.begin(); i != trkjets_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch trkjets_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trkjets_p4_branch does not exist!\n");
				exit(1);
			}
			trkjets_p4_isLoaded = true;
		}
		return trkjets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_p4()
	{
		if (not hyp_jets_p4_isLoaded) {
			if (hyp_jets_p4_branch != 0) {
				hyp_jets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_p4_.begin(); i != hyp_jets_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_p4_isLoaded = true;
		}
		return hyp_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_p4()
	{
		if (not hyp_other_jets_p4_isLoaded) {
			if (hyp_other_jets_p4_branch != 0) {
				hyp_other_jets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_p4_.begin(); i != hyp_other_jets_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_p4_isLoaded = true;
		}
		return hyp_other_jets_p4_;
	}
	vector<double> &els_inner_positionx()
	{
		if (not els_inner_positionx_isLoaded) {
			if (els_inner_positionx_branch != 0) {
				els_inner_positionx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_inner_positionz_branch does not exist!\n");
				exit(1);
			}
			els_inner_positionz_isLoaded = true;
		}
		return els_inner_positionz_;
	}
	vector<double> &jets_closestElectron_DR()
	{
		if (not jets_closestElectron_DR_isLoaded) {
			if (jets_closestElectron_DR_branch != 0) {
				jets_closestElectron_DR_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_closestMuon_DR_branch does not exist!\n");
				exit(1);
			}
			jets_closestMuon_DR_isLoaded = true;
		}
		return jets_closestMuon_DR_;
	}
	float &evt_bs_dxdz()
	{
		if (not evt_bs_dxdz_isLoaded) {
			if (evt_bs_dxdz_branch != 0) {
				evt_bs_dxdz_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_bs_dxdz_)) {
					printf("branch evt_bs_dxdz_branch contains a bad float: %f\n", evt_bs_dxdz_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_dxdzErr_)) {
					printf("branch evt_bs_dxdzErr_branch contains a bad float: %f\n", evt_bs_dxdzErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_dydz_)) {
					printf("branch evt_bs_dydz_branch contains a bad float: %f\n", evt_bs_dydz_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_dydzErr_)) {
					printf("branch evt_bs_dydzErr_branch contains a bad float: %f\n", evt_bs_dydzErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_sigmaZ_)) {
					printf("branch evt_bs_sigmaZ_branch contains a bad float: %f\n", evt_bs_sigmaZ_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_sigmaZErr_)) {
					printf("branch evt_bs_sigmaZErr_branch contains a bad float: %f\n", evt_bs_sigmaZErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bs_sigmaZErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_sigmaZErr_isLoaded = true;
		}
		return evt_bs_sigmaZErr_;
	}
	float &evt_bs_width()
	{
		if (not evt_bs_width_isLoaded) {
			if (evt_bs_width_branch != 0) {
				evt_bs_width_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_bs_width_)) {
					printf("branch evt_bs_width_branch contains a bad float: %f\n", evt_bs_width_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bs_width_branch does not exist!\n");
				exit(1);
			}
			evt_bs_width_isLoaded = true;
		}
		return evt_bs_width_;
	}
	float &evt_bs_widthErr()
	{
		if (not evt_bs_widthErr_isLoaded) {
			if (evt_bs_widthErr_branch != 0) {
				evt_bs_widthErr_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_bs_widthErr_)) {
					printf("branch evt_bs_widthErr_branch contains a bad float: %f\n", evt_bs_widthErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_bs_widthErr_branch does not exist!\n");
				exit(1);
			}
			evt_bs_widthErr_isLoaded = true;
		}
		return evt_bs_widthErr_;
	}
	float &evt_bs_xErr()
	{
		if (not evt_bs_xErr_isLoaded) {
			if (evt_bs_xErr_branch != 0) {
				evt_bs_xErr_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_bs_xErr_)) {
					printf("branch evt_bs_xErr_branch contains a bad float: %f\n", evt_bs_xErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_yErr_)) {
					printf("branch evt_bs_yErr_branch contains a bad float: %f\n", evt_bs_yErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_bs_zErr_)) {
					printf("branch evt_bs_zErr_branch contains a bad float: %f\n", evt_bs_zErr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) {
				evt_kfactor_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_weight_)) {
					printf("branch evt_weight_branch contains a bad float: %f\n", evt_weight_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(gen_met_)) {
					printf("branch gen_met_branch contains a bad float: %f\n", gen_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(gen_metPhi_)) {
					printf("branch gen_metPhi_branch contains a bad float: %f\n", gen_metPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch gen_metPhi_branch does not exist!\n");
				exit(1);
			}
			gen_metPhi_isLoaded = true;
		}
		return gen_metPhi_;
	}
	float &l1met_etHad()
	{
		if (not l1met_etHad_isLoaded) {
			if (l1met_etHad_branch != 0) {
				l1met_etHad_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1met_etHad_)) {
					printf("branch l1met_etHad_branch contains a bad float: %f\n", l1met_etHad_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1met_etHad_branch does not exist!\n");
				exit(1);
			}
			l1met_etHad_isLoaded = true;
		}
		return l1met_etHad_;
	}
	float &l1met_etTot()
	{
		if (not l1met_etTot_isLoaded) {
			if (l1met_etTot_branch != 0) {
				l1met_etTot_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1met_etTot_)) {
					printf("branch l1met_etTot_branch contains a bad float: %f\n", l1met_etTot_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1met_etTot_branch does not exist!\n");
				exit(1);
			}
			l1met_etTot_isLoaded = true;
		}
		return l1met_etTot_;
	}
	float &l1met_met()
	{
		if (not l1met_met_isLoaded) {
			if (l1met_met_branch != 0) {
				l1met_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(l1met_met_)) {
					printf("branch l1met_met_branch contains a bad float: %f\n", l1met_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1met_met_branch does not exist!\n");
				exit(1);
			}
			l1met_met_isLoaded = true;
		}
		return l1met_met_;
	}
	float &evt_met()
	{
		if (not evt_met_isLoaded) {
			if (evt_met_branch != 0) {
				evt_met_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_met_)) {
					printf("branch evt_met_branch contains a bad float: %f\n", evt_met_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metHO_)) {
					printf("branch evt_metHO_branch contains a bad float: %f\n", evt_metHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metHOPhi_)) {
					printf("branch evt_metHOPhi_branch contains a bad float: %f\n", evt_metHOPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metHOSig_)) {
					printf("branch evt_metHOSig_branch contains a bad float: %f\n", evt_metHOSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonCorr_)) {
					printf("branch evt_metMuonCorr_branch contains a bad float: %f\n", evt_metMuonCorr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonCorrPhi_)) {
					printf("branch evt_metMuonCorrPhi_branch contains a bad float: %f\n", evt_metMuonCorrPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonCorrSig_)) {
					printf("branch evt_metMuonCorrSig_branch contains a bad float: %f\n", evt_metMuonCorrSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonJESCorr_)) {
					printf("branch evt_metMuonJESCorr_branch contains a bad float: %f\n", evt_metMuonJESCorr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonJESCorrPhi_)) {
					printf("branch evt_metMuonJESCorrPhi_branch contains a bad float: %f\n", evt_metMuonJESCorrPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metMuonJESCorrSig_)) {
					printf("branch evt_metMuonJESCorrSig_branch contains a bad float: %f\n", evt_metMuonJESCorrSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHF_)) {
					printf("branch evt_metNoHF_branch contains a bad float: %f\n", evt_metNoHF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHFHO_)) {
					printf("branch evt_metNoHFHO_branch contains a bad float: %f\n", evt_metNoHFHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHFHOPhi_)) {
					printf("branch evt_metNoHFHOPhi_branch contains a bad float: %f\n", evt_metNoHFHOPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHFHOSig_)) {
					printf("branch evt_metNoHFHOSig_branch contains a bad float: %f\n", evt_metNoHFHOSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHFPhi_)) {
					printf("branch evt_metNoHFPhi_branch contains a bad float: %f\n", evt_metNoHFPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metNoHFSig_)) {
					printf("branch evt_metNoHFSig_branch contains a bad float: %f\n", evt_metNoHFSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOpt_)) {
					printf("branch evt_metOpt_branch contains a bad float: %f\n", evt_metOpt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptHO_)) {
					printf("branch evt_metOptHO_branch contains a bad float: %f\n", evt_metOptHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptHOPhi_)) {
					printf("branch evt_metOptHOPhi_branch contains a bad float: %f\n", evt_metOptHOPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptHOSig_)) {
					printf("branch evt_metOptHOSig_branch contains a bad float: %f\n", evt_metOptHOSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHF_)) {
					printf("branch evt_metOptNoHF_branch contains a bad float: %f\n", evt_metOptNoHF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHFHO_)) {
					printf("branch evt_metOptNoHFHO_branch contains a bad float: %f\n", evt_metOptNoHFHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHFHOPhi_)) {
					printf("branch evt_metOptNoHFHOPhi_branch contains a bad float: %f\n", evt_metOptNoHFHOPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHFHOSig_)) {
					printf("branch evt_metOptNoHFHOSig_branch contains a bad float: %f\n", evt_metOptNoHFHOSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHFPhi_)) {
					printf("branch evt_metOptNoHFPhi_branch contains a bad float: %f\n", evt_metOptNoHFPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptNoHFSig_)) {
					printf("branch evt_metOptNoHFSig_branch contains a bad float: %f\n", evt_metOptNoHFSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptPhi_)) {
					printf("branch evt_metOptPhi_branch contains a bad float: %f\n", evt_metOptPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metOptSig_)) {
					printf("branch evt_metOptSig_branch contains a bad float: %f\n", evt_metOptSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metPhi_)) {
					printf("branch evt_metPhi_branch contains a bad float: %f\n", evt_metPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_metSig_)) {
					printf("branch evt_metSig_branch contains a bad float: %f\n", evt_metSig_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumet_)) {
					printf("branch evt_sumet_branch contains a bad float: %f\n", evt_sumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetHO_)) {
					printf("branch evt_sumetHO_branch contains a bad float: %f\n", evt_sumetHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetMuonCorr_)) {
					printf("branch evt_sumetMuonCorr_branch contains a bad float: %f\n", evt_sumetMuonCorr_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetNoHF_)) {
					printf("branch evt_sumetNoHF_branch contains a bad float: %f\n", evt_sumetNoHF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetNoHFHO_)) {
					printf("branch evt_sumetNoHFHO_branch contains a bad float: %f\n", evt_sumetNoHFHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetOpt_)) {
					printf("branch evt_sumetOpt_branch contains a bad float: %f\n", evt_sumetOpt_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetOptHO_)) {
					printf("branch evt_sumetOptHO_branch contains a bad float: %f\n", evt_sumetOptHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetOptNoHF_)) {
					printf("branch evt_sumetOptNoHF_branch contains a bad float: %f\n", evt_sumetOptNoHF_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_sumetOptNoHFHO_)) {
					printf("branch evt_sumetOptNoHFHO_branch contains a bad float: %f\n", evt_sumetOptNoHFHO_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metCor_)) {
					printf("branch met_pat_metCor_branch contains a bad float: %f\n", met_pat_metCor_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metPhiCor_)) {
					printf("branch met_pat_metPhiCor_branch contains a bad float: %f\n", met_pat_metPhiCor_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metPhiUncor_)) {
					printf("branch met_pat_metPhiUncor_branch contains a bad float: %f\n", met_pat_metPhiUncor_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metPhiUncorJES_)) {
					printf("branch met_pat_metPhiUncorJES_branch contains a bad float: %f\n", met_pat_metPhiUncorJES_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metPhiUncorMuon_)) {
					printf("branch met_pat_metPhiUncorMuon_branch contains a bad float: %f\n", met_pat_metPhiUncorMuon_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metUncor_)) {
					printf("branch met_pat_metUncor_branch contains a bad float: %f\n", met_pat_metUncor_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metUncorJES_)) {
					printf("branch met_pat_metUncorJES_branch contains a bad float: %f\n", met_pat_metUncorJES_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(met_pat_metUncorMuon_)) {
					printf("branch met_pat_metUncorMuon_branch contains a bad float: %f\n", met_pat_metUncorMuon_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(pdfinfo_scale_)) {
					printf("branch pdfinfo_scale_branch contains a bad float: %f\n", pdfinfo_scale_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(pdfinfo_x1_)) {
					printf("branch pdfinfo_x1_branch contains a bad float: %f\n", pdfinfo_x1_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(pdfinfo_x2_)) {
					printf("branch pdfinfo_x2_branch contains a bad float: %f\n", pdfinfo_x2_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_pfmet_)) {
					printf("branch evt_pfmet_branch contains a bad float: %f\n", evt_pfmet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_pfmetPhi_)) {
					printf("branch evt_pfmetPhi_branch contains a bad float: %f\n", evt_pfmetPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_pfmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_pfmetPhi_isLoaded = true;
		}
		return evt_pfmetPhi_;
	}
	float &evt_pfsumet()
	{
		if (not evt_pfsumet_isLoaded) {
			if (evt_pfsumet_branch != 0) {
				evt_pfsumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_pfsumet_)) {
					printf("branch evt_pfsumet_branch contains a bad float: %f\n", evt_pfsumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_tcmet_)) {
					printf("branch evt_tcmet_branch contains a bad float: %f\n", evt_tcmet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				if (not isfinite(evt_tcmetPhi_)) {
					printf("branch evt_tcmetPhi_branch contains a bad float: %f\n", evt_tcmetPhi_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_tcmetPhi_branch does not exist!\n");
				exit(1);
			}
			evt_tcmetPhi_isLoaded = true;
		}
		return evt_tcmetPhi_;
	}
	float &evt_tcsumet()
	{
		if (not evt_tcsumet_isLoaded) {
			if (evt_tcsumet_branch != 0) {
				evt_tcsumet_branch->GetEntry(index);
				#ifdef PARANOIA
				if (not isfinite(evt_tcsumet_)) {
					printf("branch evt_tcsumet_branch contains a bad float: %f\n", evt_tcsumet_);
					exit(1);
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_tcsumet_branch does not exist!\n");
				exit(1);
			}
			evt_tcsumet_isLoaded = true;
		}
		return evt_tcsumet_;
	}
	vector<float> &jets_combinedSecondaryVertexBJetTag()
	{
		if (not jets_combinedSecondaryVertexBJetTag_isLoaded) {
			if (jets_combinedSecondaryVertexBJetTag_branch != 0) {
				jets_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_combinedSecondaryVertexBJetTag_.begin(); i != jets_combinedSecondaryVertexBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_combinedSecondaryVertexBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_combinedSecondaryVertexMVABJetTag_.begin(); i != jets_combinedSecondaryVertexMVABJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_combinedSecondaryVertexMVABJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_combinedSecondaryVertexMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_combinedSecondaryVertexMVABJetTag_isLoaded = true;
		}
		return jets_combinedSecondaryVertexMVABJetTag_;
	}
	vector<float> &jets_impactParameterMVABJetTag()
	{
		if (not jets_impactParameterMVABJetTag_isLoaded) {
			if (jets_impactParameterMVABJetTag_branch != 0) {
				jets_impactParameterMVABJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_impactParameterMVABJetTag_.begin(); i != jets_impactParameterMVABJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_impactParameterMVABJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_impactParameterMVABJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_impactParameterMVABJetTag_isLoaded = true;
		}
		return jets_impactParameterMVABJetTag_;
	}
	vector<float> &jets_jetBProbabilityBJetTag()
	{
		if (not jets_jetBProbabilityBJetTag_isLoaded) {
			if (jets_jetBProbabilityBJetTag_branch != 0) {
				jets_jetBProbabilityBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_jetBProbabilityBJetTag_.begin(); i != jets_jetBProbabilityBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_jetBProbabilityBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_jetProbabilityBJetTag_.begin(); i != jets_jetProbabilityBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_jetProbabilityBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_jetProbabilityBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_jetProbabilityBJetTag_isLoaded = true;
		}
		return jets_jetProbabilityBJetTag_;
	}
	vector<float> &jets_simpleSecondaryVertexBJetTag()
	{
		if (not jets_simpleSecondaryVertexBJetTag_isLoaded) {
			if (jets_simpleSecondaryVertexBJetTag_branch != 0) {
				jets_simpleSecondaryVertexBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_simpleSecondaryVertexBJetTag_.begin(); i != jets_simpleSecondaryVertexBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_simpleSecondaryVertexBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_simpleSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_simpleSecondaryVertexBJetTag_isLoaded = true;
		}
		return jets_simpleSecondaryVertexBJetTag_;
	}
	vector<float> &jets_softElectronBJetTag()
	{
		if (not jets_softElectronBJetTag_isLoaded) {
			if (jets_softElectronBJetTag_branch != 0) {
				jets_softElectronBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_softElectronBJetTag_.begin(); i != jets_softElectronBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_softElectronBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_softElectronBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softElectronBJetTag_isLoaded = true;
		}
		return jets_softElectronBJetTag_;
	}
	vector<float> &jets_softMuonBJetTag()
	{
		if (not jets_softMuonBJetTag_isLoaded) {
			if (jets_softMuonBJetTag_branch != 0) {
				jets_softMuonBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_softMuonBJetTag_.begin(); i != jets_softMuonBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_softMuonBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softMuonBJetTag_isLoaded = true;
		}
		return jets_softMuonBJetTag_;
	}
	vector<float> &jets_softMuonNoIPBJetTag()
	{
		if (not jets_softMuonNoIPBJetTag_isLoaded) {
			if (jets_softMuonNoIPBJetTag_branch != 0) {
				jets_softMuonNoIPBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_softMuonNoIPBJetTag_.begin(); i != jets_softMuonNoIPBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_softMuonNoIPBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_softMuonNoIPBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_softMuonNoIPBJetTag_isLoaded = true;
		}
		return jets_softMuonNoIPBJetTag_;
	}
	vector<float> &jets_trackCountingHighEffBJetTag()
	{
		if (not jets_trackCountingHighEffBJetTag_isLoaded) {
			if (jets_trackCountingHighEffBJetTag_branch != 0) {
				jets_trackCountingHighEffBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_trackCountingHighEffBJetTag_.begin(); i != jets_trackCountingHighEffBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_trackCountingHighEffBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_trackCountingHighPurBJetTag_.begin(); i != jets_trackCountingHighPurBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_trackCountingHighPurBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return jets_trackCountingHighPurBJetTag_;
	}
	vector<float> &taus_calo_isolationECALhitsEtSum()
	{
		if (not taus_calo_isolationECALhitsEtSum_isLoaded) {
			if (taus_calo_isolationECALhitsEtSum_branch != 0) {
				taus_calo_isolationECALhitsEtSum_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_isolationECALhitsEtSum_.begin(); i != taus_calo_isolationECALhitsEtSum_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_isolationECALhitsEtSum_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_isolationECALhitsEtSum_branch does not exist!\n");
				exit(1);
			}
			taus_calo_isolationECALhitsEtSum_isLoaded = true;
		}
		return taus_calo_isolationECALhitsEtSum_;
	}
	vector<float> &taus_calo_isolationtrksPtSum()
	{
		if (not taus_calo_isolationtrksPtSum_isLoaded) {
			if (taus_calo_isolationtrksPtSum_branch != 0) {
				taus_calo_isolationtrksPtSum_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_isolationtrksPtSum_.begin(); i != taus_calo_isolationtrksPtSum_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_isolationtrksPtSum_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_isolationtrksPtSum_branch does not exist!\n");
				exit(1);
			}
			taus_calo_isolationtrksPtSum_isLoaded = true;
		}
		return taus_calo_isolationtrksPtSum_;
	}
	vector<float> &taus_calo_leadtrk_HCAL3x3hitsEtSum()
	{
		if (not taus_calo_leadtrk_HCAL3x3hitsEtSum_isLoaded) {
			if (taus_calo_leadtrk_HCAL3x3hitsEtSum_branch != 0) {
				taus_calo_leadtrk_HCAL3x3hitsEtSum_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_HCAL3x3hitsEtSum_.begin(); i != taus_calo_leadtrk_HCAL3x3hitsEtSum_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_HCAL3x3hitsEtSum_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_HCAL3x3hitsEtSum_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_HCAL3x3hitsEtSum_isLoaded = true;
		}
		return taus_calo_leadtrk_HCAL3x3hitsEtSum_;
	}
	vector<float> &taus_calo_leadtrk_HCAL3x3hottesthitDEta()
	{
		if (not taus_calo_leadtrk_HCAL3x3hottesthitDEta_isLoaded) {
			if (taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch != 0) {
				taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_HCAL3x3hottesthitDEta_.begin(); i != taus_calo_leadtrk_HCAL3x3hottesthitDEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_HCAL3x3hottesthitDEta_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_HCAL3x3hottesthitDEta_isLoaded = true;
		}
		return taus_calo_leadtrk_HCAL3x3hottesthitDEta_;
	}
	vector<float> &taus_calo_leadtrk_Signed_Sipt()
	{
		if (not taus_calo_leadtrk_Signed_Sipt_isLoaded) {
			if (taus_calo_leadtrk_Signed_Sipt_branch != 0) {
				taus_calo_leadtrk_Signed_Sipt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_Signed_Sipt_.begin(); i != taus_calo_leadtrk_Signed_Sipt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_Signed_Sipt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_Signed_Sipt_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_Signed_Sipt_isLoaded = true;
		}
		return taus_calo_leadtrk_Signed_Sipt_;
	}
	vector<float> &taus_calo_leadtrk_chi2()
	{
		if (not taus_calo_leadtrk_chi2_isLoaded) {
			if (taus_calo_leadtrk_chi2_branch != 0) {
				taus_calo_leadtrk_chi2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_chi2_.begin(); i != taus_calo_leadtrk_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_chi2_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_chi2_isLoaded = true;
		}
		return taus_calo_leadtrk_chi2_;
	}
	vector<float> &taus_calo_leadtrk_d0()
	{
		if (not taus_calo_leadtrk_d0_isLoaded) {
			if (taus_calo_leadtrk_d0_branch != 0) {
				taus_calo_leadtrk_d0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_d0_.begin(); i != taus_calo_leadtrk_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_d0_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_d0_isLoaded = true;
		}
		return taus_calo_leadtrk_d0_;
	}
	vector<float> &taus_calo_leadtrk_lostHits()
	{
		if (not taus_calo_leadtrk_lostHits_isLoaded) {
			if (taus_calo_leadtrk_lostHits_branch != 0) {
				taus_calo_leadtrk_lostHits_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_lostHits_.begin(); i != taus_calo_leadtrk_lostHits_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_lostHits_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_lostHits_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_lostHits_isLoaded = true;
		}
		return taus_calo_leadtrk_lostHits_;
	}
	vector<float> &taus_calo_leadtrk_ndof()
	{
		if (not taus_calo_leadtrk_ndof_isLoaded) {
			if (taus_calo_leadtrk_ndof_branch != 0) {
				taus_calo_leadtrk_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_ndof_.begin(); i != taus_calo_leadtrk_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_ndof_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_ndof_isLoaded = true;
		}
		return taus_calo_leadtrk_ndof_;
	}
	vector<float> &taus_calo_leadtrk_validHits()
	{
		if (not taus_calo_leadtrk_validHits_isLoaded) {
			if (taus_calo_leadtrk_validHits_branch != 0) {
				taus_calo_leadtrk_validHits_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_validHits_.begin(); i != taus_calo_leadtrk_validHits_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_validHits_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_validHits_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_validHits_isLoaded = true;
		}
		return taus_calo_leadtrk_validHits_;
	}
	vector<float> &taus_calo_leadtrk_z0()
	{
		if (not taus_calo_leadtrk_z0_isLoaded) {
			if (taus_calo_leadtrk_z0_branch != 0) {
				taus_calo_leadtrk_z0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_leadtrk_z0_.begin(); i != taus_calo_leadtrk_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_leadtrk_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_leadtrk_z0_branch does not exist!\n");
				exit(1);
			}
			taus_calo_leadtrk_z0_isLoaded = true;
		}
		return taus_calo_leadtrk_z0_;
	}
	vector<float> &taus_calo_maximumHCALhitEt()
	{
		if (not taus_calo_maximumHCALhitEt_isLoaded) {
			if (taus_calo_maximumHCALhitEt_branch != 0) {
				taus_calo_maximumHCALhitEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_maximumHCALhitEt_.begin(); i != taus_calo_maximumHCALhitEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_maximumHCALhitEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_maximumHCALhitEt_branch does not exist!\n");
				exit(1);
			}
			taus_calo_maximumHCALhitEt_isLoaded = true;
		}
		return taus_calo_maximumHCALhitEt_;
	}
	vector<float> &taus_calo_signaltrksInvariantMass()
	{
		if (not taus_calo_signaltrksInvariantMass_isLoaded) {
			if (taus_calo_signaltrksInvariantMass_branch != 0) {
				taus_calo_signaltrksInvariantMass_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_calo_signaltrksInvariantMass_.begin(); i != taus_calo_signaltrksInvariantMass_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_calo_signaltrksInvariantMass_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_signaltrksInvariantMass_branch does not exist!\n");
				exit(1);
			}
			taus_calo_signaltrksInvariantMass_isLoaded = true;
		}
		return taus_calo_signaltrksInvariantMass_;
	}
	vector<float> &els_mc3dr()
	{
		if (not els_mc3dr_isLoaded) {
			if (els_mc3dr_branch != 0) {
				els_mc3dr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_mc3dr_.begin(); i != els_mc3dr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_mc3dr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_mcdr_.begin(); i != els_mcdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_mcdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mcdr_.begin(); i != jets_mcdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mcdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mc_emEnergy_.begin(); i != jets_mc_emEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mc_emEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mc_gpdr_.begin(); i != jets_mc_gpdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mc_gpdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mc_hadEnergy_.begin(); i != jets_mc_hadEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mc_hadEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mc_invEnergy_.begin(); i != jets_mc_invEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mc_invEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_mc_otherEnergy_.begin(); i != jets_mc_otherEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_mc_otherEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_mc3dr_.begin(); i != mus_mc3dr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_mc3dr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_mcdr_.begin(); i != mus_mcdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_mcdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trk_mc3dr_.begin(); i != trk_mc3dr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trk_mc3dr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trk_mcdr_.begin(); i != trk_mcdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trk_mcdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_conv_dcot_.begin(); i != els_conv_dcot_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_conv_dcot_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_conv_dist_.begin(); i != els_conv_dist_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_conv_dist_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_conv_dcot_.begin(); i != trks_conv_dcot_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_conv_dcot_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_conv_dist_.begin(); i != trks_conv_dist_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_conv_dist_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ecalJuraIso_.begin(); i != els_ecalJuraIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ecalJuraIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ecalJuraTowerIso_.begin(); i != els_ecalJuraTowerIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ecalJuraTowerIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_hcalConeIso_.begin(); i != els_hcalConeIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_hcalConeIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_hcalConeIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalConeIso_isLoaded = true;
		}
		return els_hcalConeIso_;
	}
	vector<float> &els_jetdr()
	{
		if (not els_jetdr_isLoaded) {
			if (els_jetdr_branch != 0) {
				els_jetdr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_jetdr_.begin(); i != els_jetdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_jetdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_musdr_.begin(); i != els_musdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_musdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_trkdr_.begin(); i != els_trkdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_trkdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_trkshFrac_.begin(); i != els_trkshFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_trkshFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_chi2_.begin(); i != els_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_d0_.begin(); i != els_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_d0Err_.begin(); i != els_d0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_d0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_d0corr_.begin(); i != els_d0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_d0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dEtaIn_.begin(); i != els_dEtaIn_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dEtaIn_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dEtaOut_.begin(); i != els_dEtaOut_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dEtaOut_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dPhiIn_.begin(); i != els_dPhiIn_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dPhiIn_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dPhiInPhiOut_.begin(); i != els_dPhiInPhiOut_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dPhiInPhiOut_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_dPhiOut_.begin(); i != els_dPhiOut_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_dPhiOut_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_e1x5_.begin(); i != els_e1x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_e1x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_e2x5Max_.begin(); i != els_e2x5Max_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_e2x5Max_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_e3x3_.begin(); i != els_e3x3_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_e3x3_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_e5x5_.begin(); i != els_e5x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_e5x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eMax_.begin(); i != els_eMax_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eMax_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eOverPIn_.begin(); i != els_eOverPIn_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eOverPIn_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSC_.begin(); i != els_eSC_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSC_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSCPresh_.begin(); i != els_eSCPresh_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSCPresh_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSCRaw_.begin(); i != els_eSCRaw_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSCRaw_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSeed_.begin(); i != els_eSeed_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSeed_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSeedOverPIn_.begin(); i != els_eSeedOverPIn_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSeedOverPIn_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_eSeedOverPOut_.begin(); i != els_eSeedOverPOut_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_eSeedOverPOut_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_eSeedOverPOut_branch does not exist!\n");
				exit(1);
			}
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	vector<float> &els_ecalIso()
	{
		if (not els_ecalIso_isLoaded) {
			if (els_ecalIso_branch != 0) {
				els_ecalIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ecalIso_.begin(); i != els_ecalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ecalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_ecalIso_branch does not exist!\n");
				exit(1);
			}
			els_ecalIso_isLoaded = true;
		}
		return els_ecalIso_;
	}
	vector<float> &els_egamma_looseId()
	{
		if (not els_egamma_looseId_isLoaded) {
			if (els_egamma_looseId_branch != 0) {
				els_egamma_looseId_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_egamma_looseId_.begin(); i != els_egamma_looseId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_egamma_looseId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_egamma_robustHighEnergy_.begin(); i != els_egamma_robustHighEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_egamma_robustHighEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_egamma_robustLooseId_.begin(); i != els_egamma_robustLooseId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_egamma_robustLooseId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_egamma_robustTightId_.begin(); i != els_egamma_robustTightId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_egamma_robustTightId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_egamma_tightId_.begin(); i != els_egamma_tightId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_egamma_tightId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_etaErr_.begin(); i != els_etaErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_etaErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_etaSC_.begin(); i != els_etaSC_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_etaSC_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_fBrem_.begin(); i != els_fBrem_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_fBrem_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_hOverE_.begin(); i != els_hOverE_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_hOverE_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_hOverE_branch does not exist!\n");
				exit(1);
			}
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	vector<float> &els_hcalIso()
	{
		if (not els_hcalIso_isLoaded) {
			if (els_hcalIso_branch != 0) {
				els_hcalIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_hcalIso_.begin(); i != els_hcalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_hcalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_hcalIso_branch does not exist!\n");
				exit(1);
			}
			els_hcalIso_isLoaded = true;
		}
		return els_hcalIso_;
	}
	vector<float> &els_layer1_charge()
	{
		if (not els_layer1_charge_isLoaded) {
			if (els_layer1_charge_branch != 0) {
				els_layer1_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_layer1_charge_.begin(); i != els_layer1_charge_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_layer1_charge_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ndof_.begin(); i != els_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_outerEta_.begin(); i != els_outerEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_outerEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_outerPhi_.begin(); i != els_outerPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_outerPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_phiErr_.begin(); i != els_phiErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_phiErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_phiSC_.begin(); i != els_phiSC_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_phiSC_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pin_.begin(); i != els_pin_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pin_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pout_.begin(); i != els_pout_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pout_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_ptErr_.begin(); i != els_ptErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_ptErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_sigmaEtaEta_.begin(); i != els_sigmaEtaEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_sigmaEtaEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_sigmaIEtaIEta_.begin(); i != els_sigmaIEtaIEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_sigmaIEtaIEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_sigmaIPhiIPhi_.begin(); i != els_sigmaIPhiIPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_sigmaIPhiIPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_sigmaPhiPhi_.begin(); i != els_sigmaPhiPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_sigmaPhiPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_tkIso_.begin(); i != els_tkIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_tkIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_tkIso_branch does not exist!\n");
				exit(1);
			}
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
	}
	vector<float> &els_vertexphi()
	{
		if (not els_vertexphi_isLoaded) {
			if (els_vertexphi_branch != 0) {
				els_vertexphi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_vertexphi_.begin(); i != els_vertexphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_vertexphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_z0_.begin(); i != els_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_z0Err_.begin(); i != els_z0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_z0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_z0corr_.begin(); i != els_z0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_z0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_z0corr_branch does not exist!\n");
				exit(1);
			}
			els_z0corr_isLoaded = true;
		}
		return els_z0corr_;
	}
	vector<float> &hyp_Ht()
	{
		if (not hyp_Ht_isLoaded) {
			if (hyp_Ht_branch != 0) {
				hyp_Ht_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_Ht_.begin(); i != hyp_Ht_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_Ht_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_Ht_branch does not exist!\n");
				exit(1);
			}
			hyp_Ht_isLoaded = true;
		}
		return hyp_Ht_;
	}
	vector<float> &hyp_dPhi_nJet_muCorrMet()
	{
		if (not hyp_dPhi_nJet_muCorrMet_isLoaded) {
			if (hyp_dPhi_nJet_muCorrMet_branch != 0) {
				hyp_dPhi_nJet_muCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_dPhi_nJet_muCorrMet_.begin(); i != hyp_dPhi_nJet_muCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_dPhi_nJet_muCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_dPhi_nJet_tcMet_.begin(); i != hyp_dPhi_nJet_tcMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_dPhi_nJet_tcMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_dPhi_nJet_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_tcMet_isLoaded = true;
		}
		return hyp_dPhi_nJet_tcMet_;
	}
	vector<float> &hyp_dPhi_nJet_type1Met()
	{
		if (not hyp_dPhi_nJet_type1Met_isLoaded) {
			if (hyp_dPhi_nJet_type1Met_branch != 0) {
				hyp_dPhi_nJet_type1Met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_dPhi_nJet_type1Met_.begin(); i != hyp_dPhi_nJet_type1Met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_dPhi_nJet_type1Met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_dPhi_nJet_type1Met_branch does not exist!\n");
				exit(1);
			}
			hyp_dPhi_nJet_type1Met_isLoaded = true;
		}
		return hyp_dPhi_nJet_type1Met_;
	}
	vector<float> &hyp_dPhi_nJet_unCorrMet()
	{
		if (not hyp_dPhi_nJet_unCorrMet_isLoaded) {
			if (hyp_dPhi_nJet_unCorrMet_branch != 0) {
				hyp_dPhi_nJet_unCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_dPhi_nJet_unCorrMet_.begin(); i != hyp_dPhi_nJet_unCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_dPhi_nJet_unCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_chi2_.begin(); i != hyp_ll_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_d0_.begin(); i != hyp_ll_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_d0Err_.begin(); i != hyp_ll_d0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_d0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_d0corr_.begin(); i != hyp_ll_d0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_d0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_d0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_d0corr_isLoaded = true;
		}
		return hyp_ll_d0corr_;
	}
	vector<float> &hyp_ll_dPhi_muCorrMet()
	{
		if (not hyp_ll_dPhi_muCorrMet_isLoaded) {
			if (hyp_ll_dPhi_muCorrMet_branch != 0) {
				hyp_ll_dPhi_muCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_dPhi_muCorrMet_.begin(); i != hyp_ll_dPhi_muCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_dPhi_muCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_dPhi_tcMet_.begin(); i != hyp_ll_dPhi_tcMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_dPhi_tcMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_dPhi_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_tcMet_isLoaded = true;
		}
		return hyp_ll_dPhi_tcMet_;
	}
	vector<float> &hyp_ll_dPhi_type1Met()
	{
		if (not hyp_ll_dPhi_type1Met_isLoaded) {
			if (hyp_ll_dPhi_type1Met_branch != 0) {
				hyp_ll_dPhi_type1Met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_dPhi_type1Met_.begin(); i != hyp_ll_dPhi_type1Met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_dPhi_type1Met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_dPhi_type1Met_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_dPhi_type1Met_isLoaded = true;
		}
		return hyp_ll_dPhi_type1Met_;
	}
	vector<float> &hyp_ll_dPhi_unCorrMet()
	{
		if (not hyp_ll_dPhi_unCorrMet_isLoaded) {
			if (hyp_ll_dPhi_unCorrMet_branch != 0) {
				hyp_ll_dPhi_unCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_dPhi_unCorrMet_.begin(); i != hyp_ll_dPhi_unCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_dPhi_unCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_etaErr_.begin(); i != hyp_ll_etaErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_etaErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_etaErr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_etaErr_isLoaded = true;
		}
		return hyp_ll_etaErr_;
	}
	vector<float> &hyp_ll_iso()
	{
		if (not hyp_ll_iso_isLoaded) {
			if (hyp_ll_iso_branch != 0) {
				hyp_ll_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_iso_.begin(); i != hyp_ll_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_iso_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_iso_isLoaded = true;
		}
		return hyp_ll_iso_;
	}
	vector<float> &hyp_ll_ndof()
	{
		if (not hyp_ll_ndof_isLoaded) {
			if (hyp_ll_ndof_branch != 0) {
				hyp_ll_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_ndof_.begin(); i != hyp_ll_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_ndof_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_ndof_isLoaded = true;
		}
		return hyp_ll_ndof_;
	}
	vector<float> &hyp_ll_outerEta()
	{
		if (not hyp_ll_outerEta_isLoaded) {
			if (hyp_ll_outerEta_branch != 0) {
				hyp_ll_outerEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_outerEta_.begin(); i != hyp_ll_outerEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_outerEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_outerEta_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_outerEta_isLoaded = true;
		}
		return hyp_ll_outerEta_;
	}
	vector<float> &hyp_ll_outerPhi()
	{
		if (not hyp_ll_outerPhi_isLoaded) {
			if (hyp_ll_outerPhi_branch != 0) {
				hyp_ll_outerPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_outerPhi_.begin(); i != hyp_ll_outerPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_outerPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_outerPhi_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_outerPhi_isLoaded = true;
		}
		return hyp_ll_outerPhi_;
	}
	vector<float> &hyp_ll_phiErr()
	{
		if (not hyp_ll_phiErr_isLoaded) {
			if (hyp_ll_phiErr_branch != 0) {
				hyp_ll_phiErr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_phiErr_.begin(); i != hyp_ll_phiErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_phiErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_ptErr_.begin(); i != hyp_ll_ptErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_ptErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_ptErr_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_ptErr_isLoaded = true;
		}
		return hyp_ll_ptErr_;
	}
	vector<float> &hyp_ll_tkIso()
	{
		if (not hyp_ll_tkIso_isLoaded) {
			if (hyp_ll_tkIso_branch != 0) {
				hyp_ll_tkIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_tkIso_.begin(); i != hyp_ll_tkIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_tkIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_tkIso_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_tkIso_isLoaded = true;
		}
		return hyp_ll_tkIso_;
	}
	vector<float> &hyp_ll_vertexphi()
	{
		if (not hyp_ll_vertexphi_isLoaded) {
			if (hyp_ll_vertexphi_branch != 0) {
				hyp_ll_vertexphi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_vertexphi_.begin(); i != hyp_ll_vertexphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_vertexphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_vertexphi_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_vertexphi_isLoaded = true;
		}
		return hyp_ll_vertexphi_;
	}
	vector<float> &hyp_ll_z0()
	{
		if (not hyp_ll_z0_isLoaded) {
			if (hyp_ll_z0_branch != 0) {
				hyp_ll_z0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_z0_.begin(); i != hyp_ll_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_z0Err_.begin(); i != hyp_ll_z0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_z0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_ll_z0corr_.begin(); i != hyp_ll_z0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_ll_z0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_chi2_.begin(); i != hyp_lt_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_d0_.begin(); i != hyp_lt_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_d0Err_.begin(); i != hyp_lt_d0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_d0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_d0corr_.begin(); i != hyp_lt_d0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_d0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_d0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_d0corr_isLoaded = true;
		}
		return hyp_lt_d0corr_;
	}
	vector<float> &hyp_lt_dPhi_muCorrMet()
	{
		if (not hyp_lt_dPhi_muCorrMet_isLoaded) {
			if (hyp_lt_dPhi_muCorrMet_branch != 0) {
				hyp_lt_dPhi_muCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_dPhi_muCorrMet_.begin(); i != hyp_lt_dPhi_muCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_dPhi_muCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_dPhi_tcMet_.begin(); i != hyp_lt_dPhi_tcMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_dPhi_tcMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_dPhi_tcMet_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_tcMet_isLoaded = true;
		}
		return hyp_lt_dPhi_tcMet_;
	}
	vector<float> &hyp_lt_dPhi_type1Met()
	{
		if (not hyp_lt_dPhi_type1Met_isLoaded) {
			if (hyp_lt_dPhi_type1Met_branch != 0) {
				hyp_lt_dPhi_type1Met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_dPhi_type1Met_.begin(); i != hyp_lt_dPhi_type1Met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_dPhi_type1Met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_dPhi_type1Met_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_dPhi_type1Met_isLoaded = true;
		}
		return hyp_lt_dPhi_type1Met_;
	}
	vector<float> &hyp_lt_dPhi_unCorrMet()
	{
		if (not hyp_lt_dPhi_unCorrMet_isLoaded) {
			if (hyp_lt_dPhi_unCorrMet_branch != 0) {
				hyp_lt_dPhi_unCorrMet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_dPhi_unCorrMet_.begin(); i != hyp_lt_dPhi_unCorrMet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_dPhi_unCorrMet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_etaErr_.begin(); i != hyp_lt_etaErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_etaErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_etaErr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_etaErr_isLoaded = true;
		}
		return hyp_lt_etaErr_;
	}
	vector<float> &hyp_lt_iso()
	{
		if (not hyp_lt_iso_isLoaded) {
			if (hyp_lt_iso_branch != 0) {
				hyp_lt_iso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_iso_.begin(); i != hyp_lt_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_iso_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_iso_isLoaded = true;
		}
		return hyp_lt_iso_;
	}
	vector<float> &hyp_lt_ndof()
	{
		if (not hyp_lt_ndof_isLoaded) {
			if (hyp_lt_ndof_branch != 0) {
				hyp_lt_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_ndof_.begin(); i != hyp_lt_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_ndof_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_ndof_isLoaded = true;
		}
		return hyp_lt_ndof_;
	}
	vector<float> &hyp_lt_outerEta()
	{
		if (not hyp_lt_outerEta_isLoaded) {
			if (hyp_lt_outerEta_branch != 0) {
				hyp_lt_outerEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_outerEta_.begin(); i != hyp_lt_outerEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_outerEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_outerEta_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_outerEta_isLoaded = true;
		}
		return hyp_lt_outerEta_;
	}
	vector<float> &hyp_lt_outerPhi()
	{
		if (not hyp_lt_outerPhi_isLoaded) {
			if (hyp_lt_outerPhi_branch != 0) {
				hyp_lt_outerPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_outerPhi_.begin(); i != hyp_lt_outerPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_outerPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_outerPhi_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_outerPhi_isLoaded = true;
		}
		return hyp_lt_outerPhi_;
	}
	vector<float> &hyp_lt_phiErr()
	{
		if (not hyp_lt_phiErr_isLoaded) {
			if (hyp_lt_phiErr_branch != 0) {
				hyp_lt_phiErr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_phiErr_.begin(); i != hyp_lt_phiErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_phiErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_ptErr_.begin(); i != hyp_lt_ptErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_ptErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_ptErr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_ptErr_isLoaded = true;
		}
		return hyp_lt_ptErr_;
	}
	vector<float> &hyp_lt_tkIso()
	{
		if (not hyp_lt_tkIso_isLoaded) {
			if (hyp_lt_tkIso_branch != 0) {
				hyp_lt_tkIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_tkIso_.begin(); i != hyp_lt_tkIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_tkIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_tkIso_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_tkIso_isLoaded = true;
		}
		return hyp_lt_tkIso_;
	}
	vector<float> &hyp_lt_vertexphi()
	{
		if (not hyp_lt_vertexphi_isLoaded) {
			if (hyp_lt_vertexphi_branch != 0) {
				hyp_lt_vertexphi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_vertexphi_.begin(); i != hyp_lt_vertexphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_vertexphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_vertexphi_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_vertexphi_isLoaded = true;
		}
		return hyp_lt_vertexphi_;
	}
	vector<float> &hyp_lt_z0()
	{
		if (not hyp_lt_z0_isLoaded) {
			if (hyp_lt_z0_branch != 0) {
				hyp_lt_z0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_z0_.begin(); i != hyp_lt_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_z0Err_.begin(); i != hyp_lt_z0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_z0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_lt_z0corr_.begin(); i != hyp_lt_z0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_lt_z0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_z0corr_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_z0corr_isLoaded = true;
		}
		return hyp_lt_z0corr_;
	}
	vector<float> &hyp_sumJetPt()
	{
		if (not hyp_sumJetPt_isLoaded) {
			if (hyp_sumJetPt_branch != 0) {
				hyp_sumJetPt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_sumJetPt_.begin(); i != hyp_sumJetPt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_sumJetPt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_sumJetPt_branch does not exist!\n");
				exit(1);
			}
			hyp_sumJetPt_isLoaded = true;
		}
		return hyp_sumJetPt_;
	}
	vector<float> &jets_EMFcor()
	{
		if (not jets_EMFcor_isLoaded) {
			if (jets_EMFcor_branch != 0) {
				jets_EMFcor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_EMFcor_.begin(); i != jets_EMFcor_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_EMFcor_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_EMFcor_branch does not exist!\n");
				exit(1);
			}
			jets_EMFcor_isLoaded = true;
		}
		return jets_EMFcor_;
	}
	vector<float> &jets_chFrac()
	{
		if (not jets_chFrac_isLoaded) {
			if (jets_chFrac_branch != 0) {
				jets_chFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_chFrac_.begin(); i != jets_chFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_chFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_chFrac_branch does not exist!\n");
				exit(1);
			}
			jets_chFrac_isLoaded = true;
		}
		return jets_chFrac_;
	}
	vector<float> &jets_cor()
	{
		if (not jets_cor_isLoaded) {
			if (jets_cor_branch != 0) {
				jets_cor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_cor_.begin(); i != jets_cor_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_cor_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_emFrac_.begin(); i != jets_emFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_emFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_emFrac_branch does not exist!\n");
				exit(1);
			}
			jets_emFrac_isLoaded = true;
		}
		return jets_emFrac_;
	}
	vector<float> &jpts_cor()
	{
		if (not jpts_cor_isLoaded) {
			if (jpts_cor_branch != 0) {
				jpts_cor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jpts_cor_.begin(); i != jpts_cor_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jpts_cor_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jpts_emFrac_.begin(); i != jpts_emFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jpts_emFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jpts_emFrac_branch does not exist!\n");
				exit(1);
			}
			jpts_emFrac_isLoaded = true;
		}
		return jpts_emFrac_;
	}
	vector<float> &jpts_jet_cor()
	{
		if (not jpts_jet_cor_isLoaded) {
			if (jpts_jet_cor_branch != 0) {
				jpts_jet_cor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jpts_jet_cor_.begin(); i != jpts_jet_cor_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jpts_jet_cor_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jpts_jet_cor_branch does not exist!\n");
				exit(1);
			}
			jpts_jet_cor_isLoaded = true;
		}
		return jpts_jet_cor_;
	}
	vector<float> &mus_met_deltax()
	{
		if (not mus_met_deltax_isLoaded) {
			if (mus_met_deltax_branch != 0) {
				mus_met_deltax_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_met_deltax_.begin(); i != mus_met_deltax_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_met_deltax_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_met_deltay_.begin(); i != mus_met_deltay_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_met_deltay_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_eledr_.begin(); i != mus_eledr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_eledr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_jetdr_.begin(); i != mus_jetdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_jetdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_jetdr_branch does not exist!\n");
				exit(1);
			}
			mus_jetdr_isLoaded = true;
		}
		return mus_jetdr_;
	}
	vector<float> &mus_trkdr()
	{
		if (not mus_trkdr_isLoaded) {
			if (mus_trkdr_branch != 0) {
				mus_trkdr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_trkdr_.begin(); i != mus_trkdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_trkdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_trkdr_branch does not exist!\n");
				exit(1);
			}
			mus_trkdr_isLoaded = true;
		}
		return mus_trkdr_;
	}
	vector<float> &mus_caloCompatibility()
	{
		if (not mus_caloCompatibility_isLoaded) {
			if (mus_caloCompatibility_branch != 0) {
				mus_caloCompatibility_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_caloCompatibility_.begin(); i != mus_caloCompatibility_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_caloCompatibility_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_chi2_.begin(); i != mus_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_d0_.begin(); i != mus_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_d0Err_.begin(); i != mus_d0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_d0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_d0corr_.begin(); i != mus_d0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_d0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_em_.begin(); i != mus_e_em_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_em_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_emS9_.begin(); i != mus_e_emS9_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_emS9_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_had_.begin(); i != mus_e_had_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_had_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_hadS9_.begin(); i != mus_e_hadS9_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_hadS9_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_ho_.begin(); i != mus_e_ho_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_ho_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_e_hoS9_.begin(); i != mus_e_hoS9_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_e_hoS9_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_etaErr_.begin(); i != mus_etaErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_etaErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_gfit_chi2_.begin(); i != mus_gfit_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_gfit_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_gfit_ndof_.begin(); i != mus_gfit_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_gfit_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso_.begin(); i != mus_iso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso03_emEt_.begin(); i != mus_iso03_emEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso03_emEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso03_hadEt_.begin(); i != mus_iso03_hadEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso03_hadEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso03_hoEt_.begin(); i != mus_iso03_hoEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso03_hoEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso03_sumPt_.begin(); i != mus_iso03_sumPt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso03_sumPt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso05_emEt_.begin(); i != mus_iso05_emEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso05_emEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso05_hadEt_.begin(); i != mus_iso05_hadEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso05_hadEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso05_hoEt_.begin(); i != mus_iso05_hoEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso05_hoEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_iso05_sumPt_.begin(); i != mus_iso05_sumPt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_iso05_sumPt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_ndof_.begin(); i != mus_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_outerEta_.begin(); i != mus_outerEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_outerEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_outerPhi_.begin(); i != mus_outerPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_outerPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_phiErr_.begin(); i != mus_phiErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_phiErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_ptErr_.begin(); i != mus_ptErr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_ptErr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_qoverp_.begin(); i != mus_qoverp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_qoverp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_qoverpError_.begin(); i != mus_qoverpError_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_qoverpError_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_vertexphi_.begin(); i != mus_vertexphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_vertexphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_z0_.begin(); i != mus_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_z0Err_.begin(); i != mus_z0Err_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_z0Err_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_z0corr_.begin(); i != mus_z0corr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_z0corr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_caloIso_.begin(); i != els_pat_caloIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_caloIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_ecalIso_.begin(); i != els_pat_ecalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_ecalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_hcalIso_.begin(); i != els_pat_hcalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_hcalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_looseId_.begin(); i != els_pat_looseId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_looseId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_robustHighEnergy_.begin(); i != els_pat_robustHighEnergy_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_robustHighEnergy_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_robustLooseId_.begin(); i != els_pat_robustLooseId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_robustLooseId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_robustTightId_.begin(); i != els_pat_robustTightId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_robustTightId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_scE1x5_.begin(); i != els_pat_scE1x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_scE1x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_scE2x5Max_.begin(); i != els_pat_scE2x5Max_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_scE2x5Max_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_scE5x5_.begin(); i != els_pat_scE5x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_scE5x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_sigmaEtaEta_.begin(); i != els_pat_sigmaEtaEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_sigmaEtaEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_sigmaIEtaIEta_.begin(); i != els_pat_sigmaIEtaIEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_sigmaIEtaIEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_tightId_.begin(); i != els_pat_tightId_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_tightId_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = els_pat_trackIso_.begin(); i != els_pat_trackIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch els_pat_trackIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_pat_trackIso_branch does not exist!\n");
				exit(1);
			}
			els_pat_trackIso_isLoaded = true;
		}
		return els_pat_trackIso_;
	}
	vector<float> &jets_pat_bCorrF()
	{
		if (not jets_pat_bCorrF_isLoaded) {
			if (jets_pat_bCorrF_branch != 0) {
				jets_pat_bCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_bCorrF_.begin(); i != jets_pat_bCorrF_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_bCorrF_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_bCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_bCorrF_isLoaded = true;
		}
		return jets_pat_bCorrF_;
	}
	vector<float> &jets_pat_cCorrF()
	{
		if (not jets_pat_cCorrF_isLoaded) {
			if (jets_pat_cCorrF_branch != 0) {
				jets_pat_cCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_cCorrF_.begin(); i != jets_pat_cCorrF_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_cCorrF_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_cCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_cCorrF_isLoaded = true;
		}
		return jets_pat_cCorrF_;
	}
	vector<float> &jets_pat_combinedSecondaryVertexBJetTag()
	{
		if (not jets_pat_combinedSecondaryVertexBJetTag_isLoaded) {
			if (jets_pat_combinedSecondaryVertexBJetTag_branch != 0) {
				jets_pat_combinedSecondaryVertexBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_combinedSecondaryVertexBJetTag_.begin(); i != jets_pat_combinedSecondaryVertexBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_combinedSecondaryVertexBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_combinedSecondaryVertexMVABJetTag_.begin(); i != jets_pat_combinedSecondaryVertexMVABJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_combinedSecondaryVertexMVABJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_coneIsolationTauJetTag_.begin(); i != jets_pat_coneIsolationTauJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_coneIsolationTauJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_coneIsolationTauJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_coneIsolationTauJetTag_isLoaded = true;
		}
		return jets_pat_coneIsolationTauJetTag_;
	}
	vector<float> &jets_pat_gluCorrF()
	{
		if (not jets_pat_gluCorrF_isLoaded) {
			if (jets_pat_gluCorrF_branch != 0) {
				jets_pat_gluCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_gluCorrF_.begin(); i != jets_pat_gluCorrF_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_gluCorrF_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_gluCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_gluCorrF_isLoaded = true;
		}
		return jets_pat_gluCorrF_;
	}
	vector<float> &jets_pat_impactParameterMVABJetTag()
	{
		if (not jets_pat_impactParameterMVABJetTag_isLoaded) {
			if (jets_pat_impactParameterMVABJetTag_branch != 0) {
				jets_pat_impactParameterMVABJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_impactParameterMVABJetTag_.begin(); i != jets_pat_impactParameterMVABJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_impactParameterMVABJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_jetBProbabilityBJetTag_.begin(); i != jets_pat_jetBProbabilityBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_jetBProbabilityBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_jetCharge_.begin(); i != jets_pat_jetCharge_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_jetCharge_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_jetProbabilityBJetTag_.begin(); i != jets_pat_jetProbabilityBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_jetProbabilityBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_noCorrF_.begin(); i != jets_pat_noCorrF_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_noCorrF_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_noCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_noCorrF_isLoaded = true;
		}
		return jets_pat_noCorrF_;
	}
	vector<float> &jets_pat_simpleSecondaryVertexBJetTag()
	{
		if (not jets_pat_simpleSecondaryVertexBJetTag_isLoaded) {
			if (jets_pat_simpleSecondaryVertexBJetTag_branch != 0) {
				jets_pat_simpleSecondaryVertexBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_simpleSecondaryVertexBJetTag_.begin(); i != jets_pat_simpleSecondaryVertexBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_simpleSecondaryVertexBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_simpleSecondaryVertexBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_simpleSecondaryVertexBJetTag_isLoaded = true;
		}
		return jets_pat_simpleSecondaryVertexBJetTag_;
	}
	vector<float> &jets_pat_softElectronBJetTag()
	{
		if (not jets_pat_softElectronBJetTag_isLoaded) {
			if (jets_pat_softElectronBJetTag_branch != 0) {
				jets_pat_softElectronBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_softElectronBJetTag_.begin(); i != jets_pat_softElectronBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_softElectronBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_softElectronBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softElectronBJetTag_isLoaded = true;
		}
		return jets_pat_softElectronBJetTag_;
	}
	vector<float> &jets_pat_softMuonBJetTag()
	{
		if (not jets_pat_softMuonBJetTag_isLoaded) {
			if (jets_pat_softMuonBJetTag_branch != 0) {
				jets_pat_softMuonBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_softMuonBJetTag_.begin(); i != jets_pat_softMuonBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_softMuonBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_softMuonBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softMuonBJetTag_isLoaded = true;
		}
		return jets_pat_softMuonBJetTag_;
	}
	vector<float> &jets_pat_softMuonNoIPBJetTag()
	{
		if (not jets_pat_softMuonNoIPBJetTag_isLoaded) {
			if (jets_pat_softMuonNoIPBJetTag_branch != 0) {
				jets_pat_softMuonNoIPBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_softMuonNoIPBJetTag_.begin(); i != jets_pat_softMuonNoIPBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_softMuonNoIPBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_softMuonNoIPBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_softMuonNoIPBJetTag_isLoaded = true;
		}
		return jets_pat_softMuonNoIPBJetTag_;
	}
	vector<float> &jets_pat_trackCountingHighEffBJetTag()
	{
		if (not jets_pat_trackCountingHighEffBJetTag_isLoaded) {
			if (jets_pat_trackCountingHighEffBJetTag_branch != 0) {
				jets_pat_trackCountingHighEffBJetTag_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_trackCountingHighEffBJetTag_.begin(); i != jets_pat_trackCountingHighEffBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_trackCountingHighEffBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_trackCountingHighPurBJetTag_.begin(); i != jets_pat_trackCountingHighPurBJetTag_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_trackCountingHighPurBJetTag_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_trackCountingHighPurBJetTag_branch does not exist!\n");
				exit(1);
			}
			jets_pat_trackCountingHighPurBJetTag_isLoaded = true;
		}
		return jets_pat_trackCountingHighPurBJetTag_;
	}
	vector<float> &jets_pat_udsCorrF()
	{
		if (not jets_pat_udsCorrF_isLoaded) {
			if (jets_pat_udsCorrF_branch != 0) {
				jets_pat_udsCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jets_pat_udsCorrF_.begin(); i != jets_pat_udsCorrF_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jets_pat_udsCorrF_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_udsCorrF_branch does not exist!\n");
				exit(1);
			}
			jets_pat_udsCorrF_isLoaded = true;
		}
		return jets_pat_udsCorrF_;
	}
	vector<float> &mus_pat_caloIso()
	{
		if (not mus_pat_caloIso_isLoaded) {
			if (mus_pat_caloIso_branch != 0) {
				mus_pat_caloIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_caloIso_.begin(); i != mus_pat_caloIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_caloIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_caloIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_caloIso_isLoaded = true;
		}
		return mus_pat_caloIso_;
	}
	vector<float> &mus_pat_ecalIso()
	{
		if (not mus_pat_ecalIso_isLoaded) {
			if (mus_pat_ecalIso_branch != 0) {
				mus_pat_ecalIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_ecalIso_.begin(); i != mus_pat_ecalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_ecalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_ecalvetoDep_.begin(); i != mus_pat_ecalvetoDep_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_ecalvetoDep_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_hcalIso_.begin(); i != mus_pat_hcalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_hcalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_hcalvetoDep_.begin(); i != mus_pat_hcalvetoDep_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_hcalvetoDep_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_trackIso_.begin(); i != mus_pat_trackIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_trackIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_trackIso_branch does not exist!\n");
				exit(1);
			}
			mus_pat_trackIso_isLoaded = true;
		}
		return mus_pat_trackIso_;
	}
	vector<float> &mus_pat_vetoDep()
	{
		if (not mus_pat_vetoDep_isLoaded) {
			if (mus_pat_vetoDep_branch != 0) {
				mus_pat_vetoDep_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_pat_vetoDep_.begin(); i != mus_pat_vetoDep_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_pat_vetoDep_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_vetoDep_branch does not exist!\n");
				exit(1);
			}
			mus_pat_vetoDep_isLoaded = true;
		}
		return mus_pat_vetoDep_;
	}
	vector<float> &taus_pf_bremsRecoveryEOverPLead()
	{
		if (not taus_pf_bremsRecoveryEOverPLead_isLoaded) {
			if (taus_pf_bremsRecoveryEOverPLead_branch != 0) {
				taus_pf_bremsRecoveryEOverPLead_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_bremsRecoveryEOverPLead_.begin(); i != taus_pf_bremsRecoveryEOverPLead_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_bremsRecoveryEOverPLead_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_bremsRecoveryEOverPLead_branch does not exist!\n");
				exit(1);
			}
			taus_pf_bremsRecoveryEOverPLead_isLoaded = true;
		}
		return taus_pf_bremsRecoveryEOverPLead_;
	}
	vector<float> &taus_pf_caloComp()
	{
		if (not taus_pf_caloComp_isLoaded) {
			if (taus_pf_caloComp_branch != 0) {
				taus_pf_caloComp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_caloComp_.begin(); i != taus_pf_caloComp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_caloComp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_caloComp_branch does not exist!\n");
				exit(1);
			}
			taus_pf_caloComp_isLoaded = true;
		}
		return taus_pf_caloComp_;
	}
	vector<float> &taus_pf_ecalStripSumEOverPLead()
	{
		if (not taus_pf_ecalStripSumEOverPLead_isLoaded) {
			if (taus_pf_ecalStripSumEOverPLead_branch != 0) {
				taus_pf_ecalStripSumEOverPLead_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_ecalStripSumEOverPLead_.begin(); i != taus_pf_ecalStripSumEOverPLead_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_ecalStripSumEOverPLead_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_ecalStripSumEOverPLead_branch does not exist!\n");
				exit(1);
			}
			taus_pf_ecalStripSumEOverPLead_isLoaded = true;
		}
		return taus_pf_ecalStripSumEOverPLead_;
	}
	vector<float> &taus_pf_emf()
	{
		if (not taus_pf_emf_isLoaded) {
			if (taus_pf_emf_branch != 0) {
				taus_pf_emf_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_emf_.begin(); i != taus_pf_emf_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_emf_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_emf_branch does not exist!\n");
				exit(1);
			}
			taus_pf_emf_isLoaded = true;
		}
		return taus_pf_emf_;
	}
	vector<float> &taus_pf_hcal3x3OverPLead()
	{
		if (not taus_pf_hcal3x3OverPLead_isLoaded) {
			if (taus_pf_hcal3x3OverPLead_branch != 0) {
				taus_pf_hcal3x3OverPLead_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_hcal3x3OverPLead_.begin(); i != taus_pf_hcal3x3OverPLead_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_hcal3x3OverPLead_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_hcal3x3OverPLead_branch does not exist!\n");
				exit(1);
			}
			taus_pf_hcal3x3OverPLead_isLoaded = true;
		}
		return taus_pf_hcal3x3OverPLead_;
	}
	vector<float> &taus_pf_hcalMaxOverPLead()
	{
		if (not taus_pf_hcalMaxOverPLead_isLoaded) {
			if (taus_pf_hcalMaxOverPLead_branch != 0) {
				taus_pf_hcalMaxOverPLead_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_hcalMaxOverPLead_.begin(); i != taus_pf_hcalMaxOverPLead_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_hcalMaxOverPLead_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_hcalMaxOverPLead_branch does not exist!\n");
				exit(1);
			}
			taus_pf_hcalMaxOverPLead_isLoaded = true;
		}
		return taus_pf_hcalMaxOverPLead_;
	}
	vector<float> &taus_pf_hcalTotOverPLead()
	{
		if (not taus_pf_hcalTotOverPLead_isLoaded) {
			if (taus_pf_hcalTotOverPLead_branch != 0) {
				taus_pf_hcalTotOverPLead_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_hcalTotOverPLead_.begin(); i != taus_pf_hcalTotOverPLead_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_hcalTotOverPLead_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_hcalTotOverPLead_branch does not exist!\n");
				exit(1);
			}
			taus_pf_hcalTotOverPLead_isLoaded = true;
		}
		return taus_pf_hcalTotOverPLead_;
	}
	vector<float> &taus_pf_isolationchargecandPtSum()
	{
		if (not taus_pf_isolationchargecandPtSum_isLoaded) {
			if (taus_pf_isolationchargecandPtSum_branch != 0) {
				taus_pf_isolationchargecandPtSum_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_isolationchargecandPtSum_.begin(); i != taus_pf_isolationchargecandPtSum_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_isolationchargecandPtSum_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_isolationchargecandPtSum_branch does not exist!\n");
				exit(1);
			}
			taus_pf_isolationchargecandPtSum_isLoaded = true;
		}
		return taus_pf_isolationchargecandPtSum_;
	}
	vector<float> &taus_pf_isolationgammacandEtSum()
	{
		if (not taus_pf_isolationgammacandEtSum_isLoaded) {
			if (taus_pf_isolationgammacandEtSum_branch != 0) {
				taus_pf_isolationgammacandEtSum_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_isolationgammacandEtSum_.begin(); i != taus_pf_isolationgammacandEtSum_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_isolationgammacandEtSum_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_isolationgammacandEtSum_branch does not exist!\n");
				exit(1);
			}
			taus_pf_isolationgammacandEtSum_isLoaded = true;
		}
		return taus_pf_isolationgammacandEtSum_;
	}
	vector<float> &taus_pf_lead_chargecand_Signed_Sipt()
	{
		if (not taus_pf_lead_chargecand_Signed_Sipt_isLoaded) {
			if (taus_pf_lead_chargecand_Signed_Sipt_branch != 0) {
				taus_pf_lead_chargecand_Signed_Sipt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_lead_chargecand_Signed_Sipt_.begin(); i != taus_pf_lead_chargecand_Signed_Sipt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_lead_chargecand_Signed_Sipt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_lead_chargecand_Signed_Sipt_branch does not exist!\n");
				exit(1);
			}
			taus_pf_lead_chargecand_Signed_Sipt_isLoaded = true;
		}
		return taus_pf_lead_chargecand_Signed_Sipt_;
	}
	vector<float> &taus_pf_leadtrk_chi2()
	{
		if (not taus_pf_leadtrk_chi2_isLoaded) {
			if (taus_pf_leadtrk_chi2_branch != 0) {
				taus_pf_leadtrk_chi2_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_chi2_.begin(); i != taus_pf_leadtrk_chi2_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_chi2_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_chi2_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_chi2_isLoaded = true;
		}
		return taus_pf_leadtrk_chi2_;
	}
	vector<float> &taus_pf_leadtrk_d0()
	{
		if (not taus_pf_leadtrk_d0_isLoaded) {
			if (taus_pf_leadtrk_d0_branch != 0) {
				taus_pf_leadtrk_d0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_d0_.begin(); i != taus_pf_leadtrk_d0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_d0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_d0_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_d0_isLoaded = true;
		}
		return taus_pf_leadtrk_d0_;
	}
	vector<float> &taus_pf_leadtrk_lostHits()
	{
		if (not taus_pf_leadtrk_lostHits_isLoaded) {
			if (taus_pf_leadtrk_lostHits_branch != 0) {
				taus_pf_leadtrk_lostHits_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_lostHits_.begin(); i != taus_pf_leadtrk_lostHits_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_lostHits_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_lostHits_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_lostHits_isLoaded = true;
		}
		return taus_pf_leadtrk_lostHits_;
	}
	vector<float> &taus_pf_leadtrk_ndof()
	{
		if (not taus_pf_leadtrk_ndof_isLoaded) {
			if (taus_pf_leadtrk_ndof_branch != 0) {
				taus_pf_leadtrk_ndof_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_ndof_.begin(); i != taus_pf_leadtrk_ndof_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_ndof_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_ndof_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_ndof_isLoaded = true;
		}
		return taus_pf_leadtrk_ndof_;
	}
	vector<float> &taus_pf_leadtrk_validHits()
	{
		if (not taus_pf_leadtrk_validHits_isLoaded) {
			if (taus_pf_leadtrk_validHits_branch != 0) {
				taus_pf_leadtrk_validHits_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_validHits_.begin(); i != taus_pf_leadtrk_validHits_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_validHits_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_validHits_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_validHits_isLoaded = true;
		}
		return taus_pf_leadtrk_validHits_;
	}
	vector<float> &taus_pf_leadtrk_z0()
	{
		if (not taus_pf_leadtrk_z0_isLoaded) {
			if (taus_pf_leadtrk_z0_branch != 0) {
				taus_pf_leadtrk_z0_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_leadtrk_z0_.begin(); i != taus_pf_leadtrk_z0_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_leadtrk_z0_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_leadtrk_z0_branch does not exist!\n");
				exit(1);
			}
			taus_pf_leadtrk_z0_isLoaded = true;
		}
		return taus_pf_leadtrk_z0_;
	}
	vector<float> &taus_pf_maximumHCALPFClusterEt()
	{
		if (not taus_pf_maximumHCALPFClusterEt_isLoaded) {
			if (taus_pf_maximumHCALPFClusterEt_branch != 0) {
				taus_pf_maximumHCALPFClusterEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_maximumHCALPFClusterEt_.begin(); i != taus_pf_maximumHCALPFClusterEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_maximumHCALPFClusterEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_maximumHCALPFClusterEt_branch does not exist!\n");
				exit(1);
			}
			taus_pf_maximumHCALPFClusterEt_isLoaded = true;
		}
		return taus_pf_maximumHCALPFClusterEt_;
	}
	vector<float> &taus_pf_segComp()
	{
		if (not taus_pf_segComp_isLoaded) {
			if (taus_pf_segComp_branch != 0) {
				taus_pf_segComp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = taus_pf_segComp_.begin(); i != taus_pf_segComp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch taus_pf_segComp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_segComp_branch does not exist!\n");
				exit(1);
			}
			taus_pf_segComp_isLoaded = true;
		}
		return taus_pf_segComp_;
	}
	vector<float> &photons_e1x5()
	{
		if (not photons_e1x5_isLoaded) {
			if (photons_e1x5_branch != 0) {
				photons_e1x5_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_e1x5_.begin(); i != photons_e1x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_e1x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_e2x5Max_.begin(); i != photons_e2x5Max_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_e2x5Max_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_e3x3_.begin(); i != photons_e3x3_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_e3x3_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_e5x5_.begin(); i != photons_e5x5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_e5x5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_e5x5_branch does not exist!\n");
				exit(1);
			}
			photons_e5x5_isLoaded = true;
		}
		return photons_e5x5_;
	}
	vector<float> &photons_eMax()
	{
		if (not photons_eMax_isLoaded) {
			if (photons_eMax_branch != 0) {
				photons_eMax_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_eMax_.begin(); i != photons_eMax_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_eMax_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_eMax_branch does not exist!\n");
				exit(1);
			}
			photons_eMax_isLoaded = true;
		}
		return photons_eMax_;
	}
	vector<float> &photons_eSC()
	{
		if (not photons_eSC_isLoaded) {
			if (photons_eSC_branch != 0) {
				photons_eSC_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_eSC_.begin(); i != photons_eSC_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_eSC_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_eSC_branch does not exist!\n");
				exit(1);
			}
			photons_eSC_isLoaded = true;
		}
		return photons_eSC_;
	}
	vector<float> &photons_eSCPresh()
	{
		if (not photons_eSCPresh_isLoaded) {
			if (photons_eSCPresh_branch != 0) {
				photons_eSCPresh_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_eSCPresh_.begin(); i != photons_eSCPresh_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_eSCPresh_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_eSCPresh_branch does not exist!\n");
				exit(1);
			}
			photons_eSCPresh_isLoaded = true;
		}
		return photons_eSCPresh_;
	}
	vector<float> &photons_eSCRaw()
	{
		if (not photons_eSCRaw_isLoaded) {
			if (photons_eSCRaw_branch != 0) {
				photons_eSCRaw_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_eSCRaw_.begin(); i != photons_eSCRaw_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_eSCRaw_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_eSCRaw_branch does not exist!\n");
				exit(1);
			}
			photons_eSCRaw_isLoaded = true;
		}
		return photons_eSCRaw_;
	}
	vector<float> &photons_eSeed()
	{
		if (not photons_eSeed_isLoaded) {
			if (photons_eSeed_branch != 0) {
				photons_eSeed_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_eSeed_.begin(); i != photons_eSeed_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_eSeed_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_eSeed_branch does not exist!\n");
				exit(1);
			}
			photons_eSeed_isLoaded = true;
		}
		return photons_eSeed_;
	}
	vector<float> &photons_ecalIso()
	{
		if (not photons_ecalIso_isLoaded) {
			if (photons_ecalIso_branch != 0) {
				photons_ecalIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_ecalIso_.begin(); i != photons_ecalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_ecalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_ecalIso_branch does not exist!\n");
				exit(1);
			}
			photons_ecalIso_isLoaded = true;
		}
		return photons_ecalIso_;
	}
	vector<float> &photons_hOverE()
	{
		if (not photons_hOverE_isLoaded) {
			if (photons_hOverE_branch != 0) {
				photons_hOverE_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_hOverE_.begin(); i != photons_hOverE_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_hOverE_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_hOverE_branch does not exist!\n");
				exit(1);
			}
			photons_hOverE_isLoaded = true;
		}
		return photons_hOverE_;
	}
	vector<float> &photons_hcalIso()
	{
		if (not photons_hcalIso_isLoaded) {
			if (photons_hcalIso_branch != 0) {
				photons_hcalIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_hcalIso_.begin(); i != photons_hcalIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_hcalIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_hcalIso_branch does not exist!\n");
				exit(1);
			}
			photons_hcalIso_isLoaded = true;
		}
		return photons_hcalIso_;
	}
	vector<float> &photons_sigmaEtaEta()
	{
		if (not photons_sigmaEtaEta_isLoaded) {
			if (photons_sigmaEtaEta_branch != 0) {
				photons_sigmaEtaEta_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_sigmaEtaEta_.begin(); i != photons_sigmaEtaEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_sigmaEtaEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_sigmaIEtaIEta_.begin(); i != photons_sigmaIEtaIEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_sigmaIEtaIEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_sigmaIEtaIEta_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaIEtaIEta_isLoaded = true;
		}
		return photons_sigmaIEtaIEta_;
	}
	vector<float> &photons_sigmaIPhiIPhi()
	{
		if (not photons_sigmaIPhiIPhi_isLoaded) {
			if (photons_sigmaIPhiIPhi_branch != 0) {
				photons_sigmaIPhiIPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_sigmaIPhiIPhi_.begin(); i != photons_sigmaIPhiIPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_sigmaIPhiIPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_sigmaIPhiIPhi_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaIPhiIPhi_isLoaded = true;
		}
		return photons_sigmaIPhiIPhi_;
	}
	vector<float> &photons_sigmaPhiPhi()
	{
		if (not photons_sigmaPhiPhi_isLoaded) {
			if (photons_sigmaPhiPhi_branch != 0) {
				photons_sigmaPhiPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_sigmaPhiPhi_.begin(); i != photons_sigmaPhiPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_sigmaPhiPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_sigmaPhiPhi_branch does not exist!\n");
				exit(1);
			}
			photons_sigmaPhiPhi_isLoaded = true;
		}
		return photons_sigmaPhiPhi_;
	}
	vector<float> &photons_tkIso()
	{
		if (not photons_tkIso_isLoaded) {
			if (photons_tkIso_branch != 0) {
				photons_tkIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = photons_tkIso_.begin(); i != photons_tkIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch photons_tkIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch photons_tkIso_branch does not exist!\n");
				exit(1);
			}
			photons_tkIso_isLoaded = true;
		}
		return photons_tkIso_;
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
	vector<float> &mus_tcmet_deltax()
	{
		if (not mus_tcmet_deltax_isLoaded) {
			if (mus_tcmet_deltax_branch != 0) {
				mus_tcmet_deltax_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_tcmet_deltax_.begin(); i != mus_tcmet_deltax_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_tcmet_deltax_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = mus_tcmet_deltay_.begin(); i != mus_tcmet_deltay_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch mus_tcmet_deltay_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
	vector<float> &trks_outerEt()
	{
		if (not trks_outerEt_isLoaded) {
			if (trks_outerEt_branch != 0) {
				trks_outerEt_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_outerEt_.begin(); i != trks_outerEt_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_outerEt_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_outerEta_.begin(); i != trks_outerEta_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_outerEta_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_outerPhi_.begin(); i != trks_outerPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_outerPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
	vector<float> &trks_tkIso()
	{
		if (not trks_tkIso_isLoaded) {
			if (trks_tkIso_branch != 0) {
				trks_tkIso_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_tkIso_.begin(); i != trks_tkIso_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_tkIso_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_vertexphi_.begin(); i != trks_vertexphi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_vertexphi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
	vector<float> &trks_elsdr()
	{
		if (not trks_elsdr_isLoaded) {
			if (trks_elsdr_branch != 0) {
				trks_elsdr_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_elsdr_.begin(); i != trks_elsdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_elsdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trks_elsshFrac_.begin(); i != trks_elsshFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trks_elsshFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = trk_musdr_.begin(); i != trk_musdr_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch trk_musdr_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
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
	int &evt_HLT1()
	{
		if (not evt_HLT1_isLoaded) {
			if (evt_HLT1_branch != 0) {
				evt_HLT1_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_L1_4_branch does not exist!\n");
				exit(1);
			}
			evt_L1_4_isLoaded = true;
		}
		return evt_L1_4_;
	}
	int &evt_nl1emiso()
	{
		if (not evt_nl1emiso_isLoaded) {
			if (evt_nl1emiso_branch != 0) {
				evt_nl1emiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1emiso_branch does not exist!\n");
				exit(1);
			}
			evt_nl1emiso_isLoaded = true;
		}
		return evt_nl1emiso_;
	}
	int &evt_nl1emnoiso()
	{
		if (not evt_nl1emnoiso_isLoaded) {
			if (evt_nl1emnoiso_branch != 0) {
				evt_nl1emnoiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1emnoiso_branch does not exist!\n");
				exit(1);
			}
			evt_nl1emnoiso_isLoaded = true;
		}
		return evt_nl1emnoiso_;
	}
	int &evt_nl1jetsc()
	{
		if (not evt_nl1jetsc_isLoaded) {
			if (evt_nl1jetsc_branch != 0) {
				evt_nl1jetsc_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1jetsc_branch does not exist!\n");
				exit(1);
			}
			evt_nl1jetsc_isLoaded = true;
		}
		return evt_nl1jetsc_;
	}
	int &evt_nl1jetsf()
	{
		if (not evt_nl1jetsf_isLoaded) {
			if (evt_nl1jetsf_branch != 0) {
				evt_nl1jetsf_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1jetsf_branch does not exist!\n");
				exit(1);
			}
			evt_nl1jetsf_isLoaded = true;
		}
		return evt_nl1jetsf_;
	}
	int &evt_nl1jetst()
	{
		if (not evt_nl1jetst_isLoaded) {
			if (evt_nl1jetst_branch != 0) {
				evt_nl1jetst_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1jetst_branch does not exist!\n");
				exit(1);
			}
			evt_nl1jetst_isLoaded = true;
		}
		return evt_nl1jetst_;
	}
	int &evt_nl1mus()
	{
		if (not evt_nl1mus_isLoaded) {
			if (evt_nl1mus_branch != 0) {
				evt_nl1mus_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nl1mus_branch does not exist!\n");
				exit(1);
			}
			evt_nl1mus_isLoaded = true;
		}
		return evt_nl1mus_;
	}
	int &pdfinfo_id1()
	{
		if (not pdfinfo_id1_isLoaded) {
			if (pdfinfo_id1_branch != 0) {
				pdfinfo_id1_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch pdfinfo_id2_branch does not exist!\n");
				exit(1);
			}
			pdfinfo_id2_isLoaded = true;
		}
		return pdfinfo_id2_;
	}
	vector<int> &taus_calo_charge()
	{
		if (not taus_calo_charge_isLoaded) {
			if (taus_calo_charge_branch != 0) {
				taus_calo_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_charge_branch does not exist!\n");
				exit(1);
			}
			taus_calo_charge_isLoaded = true;
		}
		return taus_calo_charge_;
	}
	vector<int> &taus_calo_iso_ntrks()
	{
		if (not taus_calo_iso_ntrks_isLoaded) {
			if (taus_calo_iso_ntrks_branch != 0) {
				taus_calo_iso_ntrks_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_iso_ntrks_branch does not exist!\n");
				exit(1);
			}
			taus_calo_iso_ntrks_isLoaded = true;
		}
		return taus_calo_iso_ntrks_;
	}
	vector<int> &taus_calo_sig_ntrks()
	{
		if (not taus_calo_sig_ntrks_isLoaded) {
			if (taus_calo_sig_ntrks_branch != 0) {
				taus_calo_sig_ntrks_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_calo_sig_ntrks_branch does not exist!\n");
				exit(1);
			}
			taus_calo_sig_ntrks_isLoaded = true;
		}
		return taus_calo_sig_ntrks_;
	}
	vector<int> &els_mc3_id()
	{
		if (not els_mc3_id_isLoaded) {
			if (els_mc3_id_branch != 0) {
				els_mc3_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trks_conv_tkidx_branch does not exist!\n");
				exit(1);
			}
			trks_conv_tkidx_isLoaded = true;
		}
		return trks_conv_tkidx_;
	}
	vector<int> &els_closestJet()
	{
		if (not els_closestJet_isLoaded) {
			if (els_closestJet_branch != 0) {
				els_closestJet_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch els_nSeed_branch does not exist!\n");
				exit(1);
			}
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	vector<int> &els_pass3looseId()
	{
		if (not els_pass3looseId_isLoaded) {
			if (els_pass3looseId_branch != 0) {
				els_pass3looseId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch genps_status_branch does not exist!\n");
				exit(1);
			}
			genps_status_isLoaded = true;
		}
		return genps_status_;
	}
	vector<int> &hyp_ll_charge()
	{
		if (not hyp_ll_charge_isLoaded) {
			if (hyp_ll_charge_branch != 0) {
				hyp_ll_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_lostHits_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_lostHits_isLoaded = true;
		}
		return hyp_ll_lostHits_;
	}
	vector<int> &hyp_ll_mc_id()
	{
		if (not hyp_ll_mc_id_isLoaded) {
			if (hyp_ll_mc_id_branch != 0) {
				hyp_ll_mc_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_ll_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			hyp_ll_mc_motherid_isLoaded = true;
		}
		return hyp_ll_mc_motherid_;
	}
	vector<int> &hyp_ll_validHits()
	{
		if (not hyp_ll_validHits_isLoaded) {
			if (hyp_ll_validHits_branch != 0) {
				hyp_ll_validHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_lostHits_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_lostHits_isLoaded = true;
		}
		return hyp_lt_lostHits_;
	}
	vector<int> &hyp_lt_mc_id()
	{
		if (not hyp_lt_mc_id_isLoaded) {
			if (hyp_lt_mc_id_branch != 0) {
				hyp_lt_mc_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_lt_mc_motherid_branch does not exist!\n");
				exit(1);
			}
			hyp_lt_mc_motherid_isLoaded = true;
		}
		return hyp_lt_mc_motherid_;
	}
	vector<int> &hyp_lt_validHits()
	{
		if (not hyp_lt_validHits_isLoaded) {
			if (hyp_lt_validHits_branch != 0) {
				hyp_lt_validHits_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_type_branch does not exist!\n");
				exit(1);
			}
			hyp_type_isLoaded = true;
		}
		return hyp_type_;
	}
	vector<int> &hyp_quadlep_first_type()
	{
		if (not hyp_quadlep_first_type_isLoaded) {
			if (hyp_quadlep_first_type_branch != 0) {
				hyp_quadlep_first_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_closestMuon_branch does not exist!\n");
				exit(1);
			}
			jets_closestMuon_isLoaded = true;
		}
		return jets_closestMuon_;
	}
	vector<int> &l1emiso_ieta()
	{
		if (not l1emiso_ieta_isLoaded) {
			if (l1emiso_ieta_branch != 0) {
				l1emiso_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emiso_ieta_branch does not exist!\n");
				exit(1);
			}
			l1emiso_ieta_isLoaded = true;
		}
		return l1emiso_ieta_;
	}
	vector<int> &l1emiso_iphi()
	{
		if (not l1emiso_iphi_isLoaded) {
			if (l1emiso_iphi_branch != 0) {
				l1emiso_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emiso_iphi_branch does not exist!\n");
				exit(1);
			}
			l1emiso_iphi_isLoaded = true;
		}
		return l1emiso_iphi_;
	}
	vector<int> &l1emiso_rawId()
	{
		if (not l1emiso_rawId_isLoaded) {
			if (l1emiso_rawId_branch != 0) {
				l1emiso_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emiso_rawId_branch does not exist!\n");
				exit(1);
			}
			l1emiso_rawId_isLoaded = true;
		}
		return l1emiso_rawId_;
	}
	vector<int> &l1emiso_type()
	{
		if (not l1emiso_type_isLoaded) {
			if (l1emiso_type_branch != 0) {
				l1emiso_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emiso_type_branch does not exist!\n");
				exit(1);
			}
			l1emiso_type_isLoaded = true;
		}
		return l1emiso_type_;
	}
	vector<int> &l1emnoiso_ieta()
	{
		if (not l1emnoiso_ieta_isLoaded) {
			if (l1emnoiso_ieta_branch != 0) {
				l1emnoiso_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emnoiso_ieta_branch does not exist!\n");
				exit(1);
			}
			l1emnoiso_ieta_isLoaded = true;
		}
		return l1emnoiso_ieta_;
	}
	vector<int> &l1emnoiso_iphi()
	{
		if (not l1emnoiso_iphi_isLoaded) {
			if (l1emnoiso_iphi_branch != 0) {
				l1emnoiso_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emnoiso_iphi_branch does not exist!\n");
				exit(1);
			}
			l1emnoiso_iphi_isLoaded = true;
		}
		return l1emnoiso_iphi_;
	}
	vector<int> &l1emnoiso_rawId()
	{
		if (not l1emnoiso_rawId_isLoaded) {
			if (l1emnoiso_rawId_branch != 0) {
				l1emnoiso_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emnoiso_rawId_branch does not exist!\n");
				exit(1);
			}
			l1emnoiso_rawId_isLoaded = true;
		}
		return l1emnoiso_rawId_;
	}
	vector<int> &l1emnoiso_type()
	{
		if (not l1emnoiso_type_isLoaded) {
			if (l1emnoiso_type_branch != 0) {
				l1emnoiso_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1emnoiso_type_branch does not exist!\n");
				exit(1);
			}
			l1emnoiso_type_isLoaded = true;
		}
		return l1emnoiso_type_;
	}
	vector<int> &l1jetsc_ieta()
	{
		if (not l1jetsc_ieta_isLoaded) {
			if (l1jetsc_ieta_branch != 0) {
				l1jetsc_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsc_ieta_branch does not exist!\n");
				exit(1);
			}
			l1jetsc_ieta_isLoaded = true;
		}
		return l1jetsc_ieta_;
	}
	vector<int> &l1jetsc_iphi()
	{
		if (not l1jetsc_iphi_isLoaded) {
			if (l1jetsc_iphi_branch != 0) {
				l1jetsc_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsc_iphi_branch does not exist!\n");
				exit(1);
			}
			l1jetsc_iphi_isLoaded = true;
		}
		return l1jetsc_iphi_;
	}
	vector<int> &l1jetsc_rawId()
	{
		if (not l1jetsc_rawId_isLoaded) {
			if (l1jetsc_rawId_branch != 0) {
				l1jetsc_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsc_rawId_branch does not exist!\n");
				exit(1);
			}
			l1jetsc_rawId_isLoaded = true;
		}
		return l1jetsc_rawId_;
	}
	vector<int> &l1jetsc_type()
	{
		if (not l1jetsc_type_isLoaded) {
			if (l1jetsc_type_branch != 0) {
				l1jetsc_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsc_type_branch does not exist!\n");
				exit(1);
			}
			l1jetsc_type_isLoaded = true;
		}
		return l1jetsc_type_;
	}
	vector<int> &l1jetsf_ieta()
	{
		if (not l1jetsf_ieta_isLoaded) {
			if (l1jetsf_ieta_branch != 0) {
				l1jetsf_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsf_ieta_branch does not exist!\n");
				exit(1);
			}
			l1jetsf_ieta_isLoaded = true;
		}
		return l1jetsf_ieta_;
	}
	vector<int> &l1jetsf_iphi()
	{
		if (not l1jetsf_iphi_isLoaded) {
			if (l1jetsf_iphi_branch != 0) {
				l1jetsf_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsf_iphi_branch does not exist!\n");
				exit(1);
			}
			l1jetsf_iphi_isLoaded = true;
		}
		return l1jetsf_iphi_;
	}
	vector<int> &l1jetsf_rawId()
	{
		if (not l1jetsf_rawId_isLoaded) {
			if (l1jetsf_rawId_branch != 0) {
				l1jetsf_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsf_rawId_branch does not exist!\n");
				exit(1);
			}
			l1jetsf_rawId_isLoaded = true;
		}
		return l1jetsf_rawId_;
	}
	vector<int> &l1jetsf_type()
	{
		if (not l1jetsf_type_isLoaded) {
			if (l1jetsf_type_branch != 0) {
				l1jetsf_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetsf_type_branch does not exist!\n");
				exit(1);
			}
			l1jetsf_type_isLoaded = true;
		}
		return l1jetsf_type_;
	}
	vector<int> &l1jetst_ieta()
	{
		if (not l1jetst_ieta_isLoaded) {
			if (l1jetst_ieta_branch != 0) {
				l1jetst_ieta_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetst_ieta_branch does not exist!\n");
				exit(1);
			}
			l1jetst_ieta_isLoaded = true;
		}
		return l1jetst_ieta_;
	}
	vector<int> &l1jetst_iphi()
	{
		if (not l1jetst_iphi_isLoaded) {
			if (l1jetst_iphi_branch != 0) {
				l1jetst_iphi_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetst_iphi_branch does not exist!\n");
				exit(1);
			}
			l1jetst_iphi_isLoaded = true;
		}
		return l1jetst_iphi_;
	}
	vector<int> &l1jetst_rawId()
	{
		if (not l1jetst_rawId_isLoaded) {
			if (l1jetst_rawId_branch != 0) {
				l1jetst_rawId_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetst_rawId_branch does not exist!\n");
				exit(1);
			}
			l1jetst_rawId_isLoaded = true;
		}
		return l1jetst_rawId_;
	}
	vector<int> &l1jetst_type()
	{
		if (not l1jetst_type_isLoaded) {
			if (l1jetst_type_branch != 0) {
				l1jetst_type_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1jetst_type_branch does not exist!\n");
				exit(1);
			}
			l1jetst_type_isLoaded = true;
		}
		return l1jetst_type_;
	}
	vector<int> &l1mus_flags()
	{
		if (not l1mus_flags_isLoaded) {
			if (l1mus_flags_branch != 0) {
				l1mus_flags_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1mus_flags_branch does not exist!\n");
				exit(1);
			}
			l1mus_flags_isLoaded = true;
		}
		return l1mus_flags_;
	}
	vector<int> &l1mus_q()
	{
		if (not l1mus_q_isLoaded) {
			if (l1mus_q_branch != 0) {
				l1mus_q_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1mus_q_branch does not exist!\n");
				exit(1);
			}
			l1mus_q_isLoaded = true;
		}
		return l1mus_q_;
	}
	vector<int> &l1mus_qual()
	{
		if (not l1mus_qual_isLoaded) {
			if (l1mus_qual_branch != 0) {
				l1mus_qual_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1mus_qual_branch does not exist!\n");
				exit(1);
			}
			l1mus_qual_isLoaded = true;
		}
		return l1mus_qual_;
	}
	vector<int> &l1mus_qualFlags()
	{
		if (not l1mus_qualFlags_isLoaded) {
			if (l1mus_qualFlags_branch != 0) {
				l1mus_qualFlags_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch l1mus_qualFlags_branch does not exist!\n");
				exit(1);
			}
			l1mus_qualFlags_isLoaded = true;
		}
		return l1mus_qualFlags_;
	}
	vector<int> &mus_met_flag()
	{
		if (not mus_met_flag_isLoaded) {
			if (mus_met_flag_branch != 0) {
				mus_met_flag_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_closestJet_branch does not exist!\n");
				exit(1);
			}
			mus_closestJet_isLoaded = true;
		}
		return mus_closestJet_;
	}
	vector<int> &mus_trkidx()
	{
		if (not mus_trkidx_isLoaded) {
			if (mus_trkidx_branch != 0) {
				mus_trkidx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_trkidx_branch does not exist!\n");
				exit(1);
			}
			mus_trkidx_isLoaded = true;
		}
		return mus_trkidx_;
	}
	vector<int> &mus_charge()
	{
		if (not mus_charge_isLoaded) {
			if (mus_charge_branch != 0) {
				mus_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jets_pat_genParton_id_branch does not exist!\n");
				exit(1);
			}
			jets_pat_genParton_id_isLoaded = true;
		}
		return jets_pat_genParton_id_;
	}
	vector<int> &jets_pat_partonFlavour()
	{
		if (not jets_pat_partonFlavour_isLoaded) {
			if (jets_pat_partonFlavour_branch != 0) {
				jets_pat_partonFlavour_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_genMotherID_branch does not exist!\n");
				exit(1);
			}
			mus_pat_genMotherID_isLoaded = true;
		}
		return mus_pat_genMotherID_;
	}
	vector<int> &taus_pf_charge()
	{
		if (not taus_pf_charge_isLoaded) {
			if (taus_pf_charge_branch != 0) {
				taus_pf_charge_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_charge_branch does not exist!\n");
				exit(1);
			}
			taus_pf_charge_isLoaded = true;
		}
		return taus_pf_charge_;
	}
	vector<int> &taus_pf_electronPreID()
	{
		if (not taus_pf_electronPreID_isLoaded) {
			if (taus_pf_electronPreID_branch != 0) {
				taus_pf_electronPreID_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_electronPreID_branch does not exist!\n");
				exit(1);
			}
			taus_pf_electronPreID_isLoaded = true;
		}
		return taus_pf_electronPreID_;
	}
	vector<int> &taus_pf_hasMuonReference()
	{
		if (not taus_pf_hasMuonReference_isLoaded) {
			if (taus_pf_hasMuonReference_branch != 0) {
				taus_pf_hasMuonReference_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_hasMuonReference_branch does not exist!\n");
				exit(1);
			}
			taus_pf_hasMuonReference_isLoaded = true;
		}
		return taus_pf_hasMuonReference_;
	}
	vector<int> &taus_pf_iso_ncharge_cand()
	{
		if (not taus_pf_iso_ncharge_cand_isLoaded) {
			if (taus_pf_iso_ncharge_cand_branch != 0) {
				taus_pf_iso_ncharge_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_iso_ncharge_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_iso_ncharge_cand_isLoaded = true;
		}
		return taus_pf_iso_ncharge_cand_;
	}
	vector<int> &taus_pf_iso_ngamma_cand()
	{
		if (not taus_pf_iso_ngamma_cand_isLoaded) {
			if (taus_pf_iso_ngamma_cand_branch != 0) {
				taus_pf_iso_ngamma_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_iso_ngamma_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_iso_ngamma_cand_isLoaded = true;
		}
		return taus_pf_iso_ngamma_cand_;
	}
	vector<int> &taus_pf_iso_nneutr_cand()
	{
		if (not taus_pf_iso_nneutr_cand_isLoaded) {
			if (taus_pf_iso_nneutr_cand_branch != 0) {
				taus_pf_iso_nneutr_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_iso_nneutr_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_iso_nneutr_cand_isLoaded = true;
		}
		return taus_pf_iso_nneutr_cand_;
	}
	vector<int> &taus_pf_muonPreID()
	{
		if (not taus_pf_muonPreID_isLoaded) {
			if (taus_pf_muonPreID_branch != 0) {
				taus_pf_muonPreID_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_muonPreID_branch does not exist!\n");
				exit(1);
			}
			taus_pf_muonPreID_isLoaded = true;
		}
		return taus_pf_muonPreID_;
	}
	vector<int> &taus_pf_sig_ncharge_cand()
	{
		if (not taus_pf_sig_ncharge_cand_isLoaded) {
			if (taus_pf_sig_ncharge_cand_branch != 0) {
				taus_pf_sig_ncharge_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_sig_ncharge_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_sig_ncharge_cand_isLoaded = true;
		}
		return taus_pf_sig_ncharge_cand_;
	}
	vector<int> &taus_pf_sig_ngamma_cand()
	{
		if (not taus_pf_sig_ngamma_cand_isLoaded) {
			if (taus_pf_sig_ngamma_cand_branch != 0) {
				taus_pf_sig_ngamma_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_sig_ngamma_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_sig_ngamma_cand_isLoaded = true;
		}
		return taus_pf_sig_ngamma_cand_;
	}
	vector<int> &taus_pf_sig_nneutr_cand()
	{
		if (not taus_pf_sig_nneutr_cand_isLoaded) {
			if (taus_pf_sig_nneutr_cand_branch != 0) {
				taus_pf_sig_nneutr_cand_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch taus_pf_sig_nneutr_cand_branch does not exist!\n");
				exit(1);
			}
			taus_pf_sig_nneutr_cand_isLoaded = true;
		}
		return taus_pf_sig_nneutr_cand_;
	}
	vector<int> &mus_tcmet_flag()
	{
		if (not mus_tcmet_flag_isLoaded) {
			if (mus_tcmet_flag_branch != 0) {
				mus_tcmet_flag_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_tcmet_flag_branch does not exist!\n");
				exit(1);
			}
			mus_tcmet_flag_isLoaded = true;
		}
		return mus_tcmet_flag_;
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
	vector<int> &trks_elsidx()
	{
		if (not trks_elsidx_isLoaded) {
			if (trks_elsidx_branch != 0) {
				trks_elsidx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch trk_musidx_branch does not exist!\n");
				exit(1);
			}
			trk_musidx_isLoaded = true;
		}
		return trk_musidx_;
	}
	vector<int> &hlt2ele10LWR_id()
	{
		if (not hlt2ele10LWR_id_isLoaded) {
			if (hlt2ele10LWR_id_branch != 0) {
				hlt2ele10LWR_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2ele10LWR_id_branch does not exist!\n");
				exit(1);
			}
			hlt2ele10LWR_id_isLoaded = true;
		}
		return hlt2ele10LWR_id_;
	}
	vector<int> &hlt2ele10LWR_tid()
	{
		if (not hlt2ele10LWR_tid_isLoaded) {
			if (hlt2ele10LWR_tid_branch != 0) {
				hlt2ele10LWR_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2ele10LWR_tid_branch does not exist!\n");
				exit(1);
			}
			hlt2ele10LWR_tid_isLoaded = true;
		}
		return hlt2ele10LWR_tid_;
	}
	vector<int> &hlt2mu3_id()
	{
		if (not hlt2mu3_id_isLoaded) {
			if (hlt2mu3_id_branch != 0) {
				hlt2mu3_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2mu3_id_branch does not exist!\n");
				exit(1);
			}
			hlt2mu3_id_isLoaded = true;
		}
		return hlt2mu3_id_;
	}
	vector<int> &hlt2mu3_tid()
	{
		if (not hlt2mu3_tid_isLoaded) {
			if (hlt2mu3_tid_branch != 0) {
				hlt2mu3_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hlt2mu3_tid_branch does not exist!\n");
				exit(1);
			}
			hlt2mu3_tid_isLoaded = true;
		}
		return hlt2mu3_tid_;
	}
	vector<int> &hltLisoele18LWR_id()
	{
		if (not hltLisoele18LWR_id_isLoaded) {
			if (hltLisoele18LWR_id_branch != 0) {
				hltLisoele18LWR_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltLisoele18LWR_id_branch does not exist!\n");
				exit(1);
			}
			hltLisoele18LWR_id_isLoaded = true;
		}
		return hltLisoele18LWR_id_;
	}
	vector<int> &hltLisoele18LWR_tid()
	{
		if (not hltLisoele18LWR_tid_isLoaded) {
			if (hltLisoele18LWR_tid_branch != 0) {
				hltLisoele18LWR_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltLisoele18LWR_tid_branch does not exist!\n");
				exit(1);
			}
			hltLisoele18LWR_tid_isLoaded = true;
		}
		return hltLisoele18LWR_tid_;
	}
	vector<int> &hltisoele18R_id()
	{
		if (not hltisoele18R_id_isLoaded) {
			if (hltisoele18R_id_branch != 0) {
				hltisoele18R_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltisoele18R_id_branch does not exist!\n");
				exit(1);
			}
			hltisoele18R_id_isLoaded = true;
		}
		return hltisoele18R_id_;
	}
	vector<int> &hltisoele18R_tid()
	{
		if (not hltisoele18R_tid_isLoaded) {
			if (hltisoele18R_tid_branch != 0) {
				hltisoele18R_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltisoele18R_tid_branch does not exist!\n");
				exit(1);
			}
			hltisoele18R_tid_isLoaded = true;
		}
		return hltisoele18R_tid_;
	}
	vector<int> &hltjet30_id()
	{
		if (not hltjet30_id_isLoaded) {
			if (hltjet30_id_branch != 0) {
				hltjet30_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltjet30_id_branch does not exist!\n");
				exit(1);
			}
			hltjet30_id_isLoaded = true;
		}
		return hltjet30_id_;
	}
	vector<int> &hltjet30_tid()
	{
		if (not hltjet30_tid_isLoaded) {
			if (hltjet30_tid_branch != 0) {
				hltjet30_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltjet30_tid_branch does not exist!\n");
				exit(1);
			}
			hltjet30_tid_isLoaded = true;
		}
		return hltjet30_tid_;
	}
	vector<int> &hltl1jet15_id()
	{
		if (not hltl1jet15_id_isLoaded) {
			if (hltl1jet15_id_branch != 0) {
				hltl1jet15_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1jet15_id_branch does not exist!\n");
				exit(1);
			}
			hltl1jet15_id_isLoaded = true;
		}
		return hltl1jet15_id_;
	}
	vector<int> &hltl1jet15_tid()
	{
		if (not hltl1jet15_tid_isLoaded) {
			if (hltl1jet15_tid_branch != 0) {
				hltl1jet15_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1jet15_tid_branch does not exist!\n");
				exit(1);
			}
			hltl1jet15_tid_isLoaded = true;
		}
		return hltl1jet15_tid_;
	}
	vector<int> &hltl1met20_id()
	{
		if (not hltl1met20_id_isLoaded) {
			if (hltl1met20_id_branch != 0) {
				hltl1met20_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1met20_id_branch does not exist!\n");
				exit(1);
			}
			hltl1met20_id_isLoaded = true;
		}
		return hltl1met20_id_;
	}
	vector<int> &hltl1met20_tid()
	{
		if (not hltl1met20_tid_isLoaded) {
			if (hltl1met20_tid_branch != 0) {
				hltl1met20_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1met20_tid_branch does not exist!\n");
				exit(1);
			}
			hltl1met20_tid_isLoaded = true;
		}
		return hltl1met20_tid_;
	}
	vector<int> &hltl1mu_id()
	{
		if (not hltl1mu_id_isLoaded) {
			if (hltl1mu_id_branch != 0) {
				hltl1mu_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1mu_id_branch does not exist!\n");
				exit(1);
			}
			hltl1mu_id_isLoaded = true;
		}
		return hltl1mu_id_;
	}
	vector<int> &hltl1mu_tid()
	{
		if (not hltl1mu_tid_isLoaded) {
			if (hltl1mu_tid_branch != 0) {
				hltl1mu_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl1mu_tid_branch does not exist!\n");
				exit(1);
			}
			hltl1mu_tid_isLoaded = true;
		}
		return hltl1mu_tid_;
	}
	vector<int> &hltl2mu9_id()
	{
		if (not hltl2mu9_id_isLoaded) {
			if (hltl2mu9_id_branch != 0) {
				hltl2mu9_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl2mu9_id_branch does not exist!\n");
				exit(1);
			}
			hltl2mu9_id_isLoaded = true;
		}
		return hltl2mu9_id_;
	}
	vector<int> &hltl2mu9_tid()
	{
		if (not hltl2mu9_tid_isLoaded) {
			if (hltl2mu9_tid_branch != 0) {
				hltl2mu9_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltl2mu9_tid_branch does not exist!\n");
				exit(1);
			}
			hltl2mu9_tid_isLoaded = true;
		}
		return hltl2mu9_tid_;
	}
	vector<int> &hltmet25_id()
	{
		if (not hltmet25_id_isLoaded) {
			if (hltmet25_id_branch != 0) {
				hltmet25_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmet25_id_branch does not exist!\n");
				exit(1);
			}
			hltmet25_id_isLoaded = true;
		}
		return hltmet25_id_;
	}
	vector<int> &hltmet25_tid()
	{
		if (not hltmet25_tid_isLoaded) {
			if (hltmet25_tid_branch != 0) {
				hltmet25_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmet25_tid_branch does not exist!\n");
				exit(1);
			}
			hltmet25_tid_isLoaded = true;
		}
		return hltmet25_tid_;
	}
	vector<int> &hltmu11_id()
	{
		if (not hltmu11_id_isLoaded) {
			if (hltmu11_id_branch != 0) {
				hltmu11_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu11_id_branch does not exist!\n");
				exit(1);
			}
			hltmu11_id_isLoaded = true;
		}
		return hltmu11_id_;
	}
	vector<int> &hltmu11_tid()
	{
		if (not hltmu11_tid_isLoaded) {
			if (hltmu11_tid_branch != 0) {
				hltmu11_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu11_tid_branch does not exist!\n");
				exit(1);
			}
			hltmu11_tid_isLoaded = true;
		}
		return hltmu11_tid_;
	}
	vector<int> &hltmu9_id()
	{
		if (not hltmu9_id_isLoaded) {
			if (hltmu9_id_branch != 0) {
				hltmu9_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu9_id_branch does not exist!\n");
				exit(1);
			}
			hltmu9_id_isLoaded = true;
		}
		return hltmu9_id_;
	}
	vector<int> &hltmu9_tid()
	{
		if (not hltmu9_tid_isLoaded) {
			if (hltmu9_tid_branch != 0) {
				hltmu9_tid_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hltmu9_tid_branch does not exist!\n");
				exit(1);
			}
			hltmu9_tid_isLoaded = true;
		}
		return hltmu9_tid_;
	}
	vector<int> &vtxs_hasRefittedTracks()
	{
		if (not vtxs_hasRefittedTracks_isLoaded) {
			if (vtxs_hasRefittedTracks_branch != 0) {
				vtxs_hasRefittedTracks_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
	vector<vector<int> > &hyp_jets_idx()
	{
		if (not hyp_jets_idx_isLoaded) {
			if (hyp_jets_idx_branch != 0) {
				hyp_jets_idx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
	unsigned int &evt_ngenjets()
	{
		if (not evt_ngenjets_isLoaded) {
			if (evt_ngenjets_branch != 0) {
				evt_ngenjets_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_ngenjets_branch does not exist!\n");
				exit(1);
			}
			evt_ngenjets_isLoaded = true;
		}
		return evt_ngenjets_;
	}
	unsigned int &evt_njets()
	{
		if (not evt_njets_isLoaded) {
			if (evt_njets_branch != 0) {
				evt_njets_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_njpts_branch does not exist!\n");
				exit(1);
			}
			evt_njpts_isLoaded = true;
		}
		return evt_njpts_;
	}
	unsigned int &evt_nphotons()
	{
		if (not evt_nphotons_isLoaded) {
			if (evt_nphotons_branch != 0) {
				evt_nphotons_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nphotons_branch does not exist!\n");
				exit(1);
			}
			evt_nphotons_isLoaded = true;
		}
		return evt_nphotons_;
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
	unsigned int &evt_ntrkjets()
	{
		if (not evt_ntrkjets_isLoaded) {
			if (evt_ntrkjets_branch != 0) {
				evt_ntrkjets_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
	vector<unsigned int> &hyp_quadlep_bucket()
	{
		if (not hyp_quadlep_bucket_isLoaded) {
			if (hyp_quadlep_bucket_branch != 0) {
				hyp_quadlep_bucket_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_trilep_third_index_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_third_index_isLoaded = true;
		}
		return hyp_trilep_third_index_;
	}
	vector<unsigned int> &els_pat_flag()
	{
		if (not els_pat_flag_isLoaded) {
			if (els_pat_flag_branch != 0) {
				els_pat_flag_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch mus_pat_flag_branch does not exist!\n");
				exit(1);
			}
			mus_pat_flag_isLoaded = true;
		}
		return mus_pat_flag_;
	}
	float &evt_scale1fb()
	{
		if (not evt_scale1fb_isLoaded) {
			if (evt_scale1fb_branch != 0) {
				evt_scale1fb_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_scale1fb_branch does not exist!\n");
				exit(1);
			}
			evt_scale1fb_isLoaded = true;
		}
		return evt_scale1fb_;
	}
	int &evt_nEvts()
	{
		if (not evt_nEvts_isLoaded) {
			if (evt_nEvts_branch != 0) {
				evt_nEvts_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
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
	TString &evt_dataset() { return cms2.evt_dataset(); }
	vector<TString> &evt_HLT_trigNames() { return cms2.evt_HLT_trigNames(); }
	vector<TString> &evt_L1_trigNames() { return cms2.evt_L1_trigNames(); }
	double &genps_pthat() { return cms2.genps_pthat(); }
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &vtxs_position() { return cms2.vtxs_position(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &evt_bsp4() { return cms2.evt_bsp4(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &l1met_p4() { return cms2.l1met_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_calo_leadtrk_p4() { return cms2.taus_calo_leadtrk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_calo_p4() { return cms2.taus_calo_p4(); }
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genjets_p4() { return cms2.genjets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_lepdaughter_p4() { return cms2.genps_lepdaughter_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_p4() { return cms2.genps_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_prod_vtx() { return cms2.genps_prod_vtx(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_mc_p4() { return cms2.hyp_ll_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_p4() { return cms2.hyp_ll_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_trk_p4() { return cms2.hyp_ll_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_mc_p4() { return cms2.hyp_lt_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_p4() { return cms2.hyp_lt_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_trk_p4() { return cms2.hyp_lt_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_p4() { return cms2.hyp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_p4() { return cms2.jets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jpts_p4() { return cms2.jpts_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emiso_p4() { return cms2.l1emiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emnoiso_p4() { return cms2.l1emnoiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsc_p4() { return cms2.l1jetsc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsf_p4() { return cms2.l1jetsf_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetst_p4() { return cms2.l1jetst_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1mus_p4() { return cms2.l1mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_gfit_outerPos_p4() { return cms2.mus_gfit_outerPos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4() { return cms2.mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4() { return cms2.mus_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_vertex_p4() { return cms2.mus_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genMotherP4() { return cms2.els_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genP4() { return cms2.els_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_p4() { return cms2.els_pat_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genJet_p4() { return cms2.jets_pat_genJet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genPartonMother_p4() { return cms2.jets_pat_genPartonMother_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genParton_p4() { return cms2.jets_pat_genParton_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_p4() { return cms2.jets_pat_jet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_uncorp4() { return cms2.jets_pat_jet_uncorp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genMotherP4() { return cms2.mus_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genP4() { return cms2.mus_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_p4() { return cms2.mus_pat_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_lead_chargecand_p4() { return cms2.taus_pf_lead_chargecand_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_lead_neutrcand_p4() { return cms2.taus_pf_lead_neutrcand_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_leadtrk_p4() { return cms2.taus_pf_leadtrk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &taus_pf_p4() { return cms2.taus_pf_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &photons_p4() { return cms2.photons_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_p4() { return cms2.scs_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_pos() { return cms2.scs_pos(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_vtx() { return cms2.scs_vtx(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4() { return cms2.trks_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_vertex_p4() { return cms2.trks_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hlt2ele10LWR_p4() { return cms2.hlt2ele10LWR_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hlt2mu3_p4() { return cms2.hlt2mu3_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltLisoele18LWR_p4() { return cms2.hltLisoele18LWR_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltisoele18R_p4() { return cms2.hltisoele18R_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltjet30_p4() { return cms2.hltjet30_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1jet15_p4() { return cms2.hltl1jet15_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1met20_p4() { return cms2.hltl1met20_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl1mu_p4() { return cms2.hltl1mu_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltl2mu9_p4() { return cms2.hltl2mu9_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmet25_p4() { return cms2.hltmet25_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmu11_p4() { return cms2.hltmu11_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hltmu9_p4() { return cms2.hltmu9_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trkjets_p4() { return cms2.trkjets_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_p4() { return cms2.hyp_jets_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_p4() { return cms2.hyp_other_jets_p4(); }
	vector<double> &els_inner_positionx() { return cms2.els_inner_positionx(); }
	vector<double> &els_inner_positiony() { return cms2.els_inner_positiony(); }
	vector<double> &els_inner_positionz() { return cms2.els_inner_positionz(); }
	vector<double> &jets_closestElectron_DR() { return cms2.jets_closestElectron_DR(); }
	vector<double> &jets_closestMuon_DR() { return cms2.jets_closestMuon_DR(); }
	float &evt_bs_dxdz() { return cms2.evt_bs_dxdz(); }
	float &evt_bs_dxdzErr() { return cms2.evt_bs_dxdzErr(); }
	float &evt_bs_dydz() { return cms2.evt_bs_dydz(); }
	float &evt_bs_dydzErr() { return cms2.evt_bs_dydzErr(); }
	float &evt_bs_sigmaZ() { return cms2.evt_bs_sigmaZ(); }
	float &evt_bs_sigmaZErr() { return cms2.evt_bs_sigmaZErr(); }
	float &evt_bs_width() { return cms2.evt_bs_width(); }
	float &evt_bs_widthErr() { return cms2.evt_bs_widthErr(); }
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
	float &l1met_etHad() { return cms2.l1met_etHad(); }
	float &l1met_etTot() { return cms2.l1met_etTot(); }
	float &l1met_met() { return cms2.l1met_met(); }
	float &evt_met() { return cms2.evt_met(); }
	float &evt_metHO() { return cms2.evt_metHO(); }
	float &evt_metHOPhi() { return cms2.evt_metHOPhi(); }
	float &evt_metHOSig() { return cms2.evt_metHOSig(); }
	float &evt_metMuonCorr() { return cms2.evt_metMuonCorr(); }
	float &evt_metMuonCorrPhi() { return cms2.evt_metMuonCorrPhi(); }
	float &evt_metMuonCorrSig() { return cms2.evt_metMuonCorrSig(); }
	float &evt_metMuonJESCorr() { return cms2.evt_metMuonJESCorr(); }
	float &evt_metMuonJESCorrPhi() { return cms2.evt_metMuonJESCorrPhi(); }
	float &evt_metMuonJESCorrSig() { return cms2.evt_metMuonJESCorrSig(); }
	float &evt_metNoHF() { return cms2.evt_metNoHF(); }
	float &evt_metNoHFHO() { return cms2.evt_metNoHFHO(); }
	float &evt_metNoHFHOPhi() { return cms2.evt_metNoHFHOPhi(); }
	float &evt_metNoHFHOSig() { return cms2.evt_metNoHFHOSig(); }
	float &evt_metNoHFPhi() { return cms2.evt_metNoHFPhi(); }
	float &evt_metNoHFSig() { return cms2.evt_metNoHFSig(); }
	float &evt_metOpt() { return cms2.evt_metOpt(); }
	float &evt_metOptHO() { return cms2.evt_metOptHO(); }
	float &evt_metOptHOPhi() { return cms2.evt_metOptHOPhi(); }
	float &evt_metOptHOSig() { return cms2.evt_metOptHOSig(); }
	float &evt_metOptNoHF() { return cms2.evt_metOptNoHF(); }
	float &evt_metOptNoHFHO() { return cms2.evt_metOptNoHFHO(); }
	float &evt_metOptNoHFHOPhi() { return cms2.evt_metOptNoHFHOPhi(); }
	float &evt_metOptNoHFHOSig() { return cms2.evt_metOptNoHFHOSig(); }
	float &evt_metOptNoHFPhi() { return cms2.evt_metOptNoHFPhi(); }
	float &evt_metOptNoHFSig() { return cms2.evt_metOptNoHFSig(); }
	float &evt_metOptPhi() { return cms2.evt_metOptPhi(); }
	float &evt_metOptSig() { return cms2.evt_metOptSig(); }
	float &evt_metPhi() { return cms2.evt_metPhi(); }
	float &evt_metSig() { return cms2.evt_metSig(); }
	float &evt_sumet() { return cms2.evt_sumet(); }
	float &evt_sumetHO() { return cms2.evt_sumetHO(); }
	float &evt_sumetMuonCorr() { return cms2.evt_sumetMuonCorr(); }
	float &evt_sumetNoHF() { return cms2.evt_sumetNoHF(); }
	float &evt_sumetNoHFHO() { return cms2.evt_sumetNoHFHO(); }
	float &evt_sumetOpt() { return cms2.evt_sumetOpt(); }
	float &evt_sumetOptHO() { return cms2.evt_sumetOptHO(); }
	float &evt_sumetOptNoHF() { return cms2.evt_sumetOptNoHF(); }
	float &evt_sumetOptNoHFHO() { return cms2.evt_sumetOptNoHFHO(); }
	float &met_pat_metCor() { return cms2.met_pat_metCor(); }
	float &met_pat_metPhiCor() { return cms2.met_pat_metPhiCor(); }
	float &met_pat_metPhiUncor() { return cms2.met_pat_metPhiUncor(); }
	float &met_pat_metPhiUncorJES() { return cms2.met_pat_metPhiUncorJES(); }
	float &met_pat_metPhiUncorMuon() { return cms2.met_pat_metPhiUncorMuon(); }
	float &met_pat_metUncor() { return cms2.met_pat_metUncor(); }
	float &met_pat_metUncorJES() { return cms2.met_pat_metUncorJES(); }
	float &met_pat_metUncorMuon() { return cms2.met_pat_metUncorMuon(); }
	float &pdfinfo_scale() { return cms2.pdfinfo_scale(); }
	float &pdfinfo_x1() { return cms2.pdfinfo_x1(); }
	float &pdfinfo_x2() { return cms2.pdfinfo_x2(); }
	float &evt_pfmet() { return cms2.evt_pfmet(); }
	float &evt_pfmetPhi() { return cms2.evt_pfmetPhi(); }
	float &evt_pfsumet() { return cms2.evt_pfsumet(); }
	float &evt_tcmet() { return cms2.evt_tcmet(); }
	float &evt_tcmetPhi() { return cms2.evt_tcmetPhi(); }
	float &evt_tcsumet() { return cms2.evt_tcsumet(); }
	vector<float> &jets_combinedSecondaryVertexBJetTag() { return cms2.jets_combinedSecondaryVertexBJetTag(); }
	vector<float> &jets_combinedSecondaryVertexMVABJetTag() { return cms2.jets_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &jets_impactParameterMVABJetTag() { return cms2.jets_impactParameterMVABJetTag(); }
	vector<float> &jets_jetBProbabilityBJetTag() { return cms2.jets_jetBProbabilityBJetTag(); }
	vector<float> &jets_jetProbabilityBJetTag() { return cms2.jets_jetProbabilityBJetTag(); }
	vector<float> &jets_simpleSecondaryVertexBJetTag() { return cms2.jets_simpleSecondaryVertexBJetTag(); }
	vector<float> &jets_softElectronBJetTag() { return cms2.jets_softElectronBJetTag(); }
	vector<float> &jets_softMuonBJetTag() { return cms2.jets_softMuonBJetTag(); }
	vector<float> &jets_softMuonNoIPBJetTag() { return cms2.jets_softMuonNoIPBJetTag(); }
	vector<float> &jets_trackCountingHighEffBJetTag() { return cms2.jets_trackCountingHighEffBJetTag(); }
	vector<float> &jets_trackCountingHighPurBJetTag() { return cms2.jets_trackCountingHighPurBJetTag(); }
	vector<float> &taus_calo_isolationECALhitsEtSum() { return cms2.taus_calo_isolationECALhitsEtSum(); }
	vector<float> &taus_calo_isolationtrksPtSum() { return cms2.taus_calo_isolationtrksPtSum(); }
	vector<float> &taus_calo_leadtrk_HCAL3x3hitsEtSum() { return cms2.taus_calo_leadtrk_HCAL3x3hitsEtSum(); }
	vector<float> &taus_calo_leadtrk_HCAL3x3hottesthitDEta() { return cms2.taus_calo_leadtrk_HCAL3x3hottesthitDEta(); }
	vector<float> &taus_calo_leadtrk_Signed_Sipt() { return cms2.taus_calo_leadtrk_Signed_Sipt(); }
	vector<float> &taus_calo_leadtrk_chi2() { return cms2.taus_calo_leadtrk_chi2(); }
	vector<float> &taus_calo_leadtrk_d0() { return cms2.taus_calo_leadtrk_d0(); }
	vector<float> &taus_calo_leadtrk_lostHits() { return cms2.taus_calo_leadtrk_lostHits(); }
	vector<float> &taus_calo_leadtrk_ndof() { return cms2.taus_calo_leadtrk_ndof(); }
	vector<float> &taus_calo_leadtrk_validHits() { return cms2.taus_calo_leadtrk_validHits(); }
	vector<float> &taus_calo_leadtrk_z0() { return cms2.taus_calo_leadtrk_z0(); }
	vector<float> &taus_calo_maximumHCALhitEt() { return cms2.taus_calo_maximumHCALhitEt(); }
	vector<float> &taus_calo_signaltrksInvariantMass() { return cms2.taus_calo_signaltrksInvariantMass(); }
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
	vector<float> &els_ecalJuraIso() { return cms2.els_ecalJuraIso(); }
	vector<float> &els_ecalJuraTowerIso() { return cms2.els_ecalJuraTowerIso(); }
	vector<float> &els_hcalConeIso() { return cms2.els_hcalConeIso(); }
	vector<float> &els_jetdr() { return cms2.els_jetdr(); }
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
	vector<float> &els_ecalIso() { return cms2.els_ecalIso(); }
	vector<float> &els_egamma_looseId() { return cms2.els_egamma_looseId(); }
	vector<float> &els_egamma_robustHighEnergy() { return cms2.els_egamma_robustHighEnergy(); }
	vector<float> &els_egamma_robustLooseId() { return cms2.els_egamma_robustLooseId(); }
	vector<float> &els_egamma_robustTightId() { return cms2.els_egamma_robustTightId(); }
	vector<float> &els_egamma_tightId() { return cms2.els_egamma_tightId(); }
	vector<float> &els_etaErr() { return cms2.els_etaErr(); }
	vector<float> &els_etaSC() { return cms2.els_etaSC(); }
	vector<float> &els_fBrem() { return cms2.els_fBrem(); }
	vector<float> &els_hOverE() { return cms2.els_hOverE(); }
	vector<float> &els_hcalIso() { return cms2.els_hcalIso(); }
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
	vector<float> &els_tkIso() { return cms2.els_tkIso(); }
	vector<float> &els_vertexphi() { return cms2.els_vertexphi(); }
	vector<float> &els_z0() { return cms2.els_z0(); }
	vector<float> &els_z0Err() { return cms2.els_z0Err(); }
	vector<float> &els_z0corr() { return cms2.els_z0corr(); }
	vector<float> &hyp_Ht() { return cms2.hyp_Ht(); }
	vector<float> &hyp_dPhi_nJet_muCorrMet() { return cms2.hyp_dPhi_nJet_muCorrMet(); }
	vector<float> &hyp_dPhi_nJet_tcMet() { return cms2.hyp_dPhi_nJet_tcMet(); }
	vector<float> &hyp_dPhi_nJet_type1Met() { return cms2.hyp_dPhi_nJet_type1Met(); }
	vector<float> &hyp_dPhi_nJet_unCorrMet() { return cms2.hyp_dPhi_nJet_unCorrMet(); }
	vector<float> &hyp_ll_chi2() { return cms2.hyp_ll_chi2(); }
	vector<float> &hyp_ll_d0() { return cms2.hyp_ll_d0(); }
	vector<float> &hyp_ll_d0Err() { return cms2.hyp_ll_d0Err(); }
	vector<float> &hyp_ll_d0corr() { return cms2.hyp_ll_d0corr(); }
	vector<float> &hyp_ll_dPhi_muCorrMet() { return cms2.hyp_ll_dPhi_muCorrMet(); }
	vector<float> &hyp_ll_dPhi_tcMet() { return cms2.hyp_ll_dPhi_tcMet(); }
	vector<float> &hyp_ll_dPhi_type1Met() { return cms2.hyp_ll_dPhi_type1Met(); }
	vector<float> &hyp_ll_dPhi_unCorrMet() { return cms2.hyp_ll_dPhi_unCorrMet(); }
	vector<float> &hyp_ll_etaErr() { return cms2.hyp_ll_etaErr(); }
	vector<float> &hyp_ll_iso() { return cms2.hyp_ll_iso(); }
	vector<float> &hyp_ll_ndof() { return cms2.hyp_ll_ndof(); }
	vector<float> &hyp_ll_outerEta() { return cms2.hyp_ll_outerEta(); }
	vector<float> &hyp_ll_outerPhi() { return cms2.hyp_ll_outerPhi(); }
	vector<float> &hyp_ll_phiErr() { return cms2.hyp_ll_phiErr(); }
	vector<float> &hyp_ll_ptErr() { return cms2.hyp_ll_ptErr(); }
	vector<float> &hyp_ll_tkIso() { return cms2.hyp_ll_tkIso(); }
	vector<float> &hyp_ll_vertexphi() { return cms2.hyp_ll_vertexphi(); }
	vector<float> &hyp_ll_z0() { return cms2.hyp_ll_z0(); }
	vector<float> &hyp_ll_z0Err() { return cms2.hyp_ll_z0Err(); }
	vector<float> &hyp_ll_z0corr() { return cms2.hyp_ll_z0corr(); }
	vector<float> &hyp_lt_chi2() { return cms2.hyp_lt_chi2(); }
	vector<float> &hyp_lt_d0() { return cms2.hyp_lt_d0(); }
	vector<float> &hyp_lt_d0Err() { return cms2.hyp_lt_d0Err(); }
	vector<float> &hyp_lt_d0corr() { return cms2.hyp_lt_d0corr(); }
	vector<float> &hyp_lt_dPhi_muCorrMet() { return cms2.hyp_lt_dPhi_muCorrMet(); }
	vector<float> &hyp_lt_dPhi_tcMet() { return cms2.hyp_lt_dPhi_tcMet(); }
	vector<float> &hyp_lt_dPhi_type1Met() { return cms2.hyp_lt_dPhi_type1Met(); }
	vector<float> &hyp_lt_dPhi_unCorrMet() { return cms2.hyp_lt_dPhi_unCorrMet(); }
	vector<float> &hyp_lt_etaErr() { return cms2.hyp_lt_etaErr(); }
	vector<float> &hyp_lt_iso() { return cms2.hyp_lt_iso(); }
	vector<float> &hyp_lt_ndof() { return cms2.hyp_lt_ndof(); }
	vector<float> &hyp_lt_outerEta() { return cms2.hyp_lt_outerEta(); }
	vector<float> &hyp_lt_outerPhi() { return cms2.hyp_lt_outerPhi(); }
	vector<float> &hyp_lt_phiErr() { return cms2.hyp_lt_phiErr(); }
	vector<float> &hyp_lt_ptErr() { return cms2.hyp_lt_ptErr(); }
	vector<float> &hyp_lt_tkIso() { return cms2.hyp_lt_tkIso(); }
	vector<float> &hyp_lt_vertexphi() { return cms2.hyp_lt_vertexphi(); }
	vector<float> &hyp_lt_z0() { return cms2.hyp_lt_z0(); }
	vector<float> &hyp_lt_z0Err() { return cms2.hyp_lt_z0Err(); }
	vector<float> &hyp_lt_z0corr() { return cms2.hyp_lt_z0corr(); }
	vector<float> &hyp_sumJetPt() { return cms2.hyp_sumJetPt(); }
	vector<float> &jets_EMFcor() { return cms2.jets_EMFcor(); }
	vector<float> &jets_chFrac() { return cms2.jets_chFrac(); }
	vector<float> &jets_cor() { return cms2.jets_cor(); }
	vector<float> &jets_emFrac() { return cms2.jets_emFrac(); }
	vector<float> &jpts_cor() { return cms2.jpts_cor(); }
	vector<float> &jpts_emFrac() { return cms2.jpts_emFrac(); }
	vector<float> &jpts_jet_cor() { return cms2.jpts_jet_cor(); }
	vector<float> &mus_met_deltax() { return cms2.mus_met_deltax(); }
	vector<float> &mus_met_deltay() { return cms2.mus_met_deltay(); }
	vector<float> &mus_eledr() { return cms2.mus_eledr(); }
	vector<float> &mus_jetdr() { return cms2.mus_jetdr(); }
	vector<float> &mus_trkdr() { return cms2.mus_trkdr(); }
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
	vector<float> &els_pat_caloIso() { return cms2.els_pat_caloIso(); }
	vector<float> &els_pat_ecalIso() { return cms2.els_pat_ecalIso(); }
	vector<float> &els_pat_hcalIso() { return cms2.els_pat_hcalIso(); }
	vector<float> &els_pat_looseId() { return cms2.els_pat_looseId(); }
	vector<float> &els_pat_robustHighEnergy() { return cms2.els_pat_robustHighEnergy(); }
	vector<float> &els_pat_robustLooseId() { return cms2.els_pat_robustLooseId(); }
	vector<float> &els_pat_robustTightId() { return cms2.els_pat_robustTightId(); }
	vector<float> &els_pat_scE1x5() { return cms2.els_pat_scE1x5(); }
	vector<float> &els_pat_scE2x5Max() { return cms2.els_pat_scE2x5Max(); }
	vector<float> &els_pat_scE5x5() { return cms2.els_pat_scE5x5(); }
	vector<float> &els_pat_sigmaEtaEta() { return cms2.els_pat_sigmaEtaEta(); }
	vector<float> &els_pat_sigmaIEtaIEta() { return cms2.els_pat_sigmaIEtaIEta(); }
	vector<float> &els_pat_tightId() { return cms2.els_pat_tightId(); }
	vector<float> &els_pat_trackIso() { return cms2.els_pat_trackIso(); }
	vector<float> &jets_pat_bCorrF() { return cms2.jets_pat_bCorrF(); }
	vector<float> &jets_pat_cCorrF() { return cms2.jets_pat_cCorrF(); }
	vector<float> &jets_pat_combinedSecondaryVertexBJetTag() { return cms2.jets_pat_combinedSecondaryVertexBJetTag(); }
	vector<float> &jets_pat_combinedSecondaryVertexMVABJetTag() { return cms2.jets_pat_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &jets_pat_coneIsolationTauJetTag() { return cms2.jets_pat_coneIsolationTauJetTag(); }
	vector<float> &jets_pat_gluCorrF() { return cms2.jets_pat_gluCorrF(); }
	vector<float> &jets_pat_impactParameterMVABJetTag() { return cms2.jets_pat_impactParameterMVABJetTag(); }
	vector<float> &jets_pat_jetBProbabilityBJetTag() { return cms2.jets_pat_jetBProbabilityBJetTag(); }
	vector<float> &jets_pat_jetCharge() { return cms2.jets_pat_jetCharge(); }
	vector<float> &jets_pat_jetProbabilityBJetTag() { return cms2.jets_pat_jetProbabilityBJetTag(); }
	vector<float> &jets_pat_noCorrF() { return cms2.jets_pat_noCorrF(); }
	vector<float> &jets_pat_simpleSecondaryVertexBJetTag() { return cms2.jets_pat_simpleSecondaryVertexBJetTag(); }
	vector<float> &jets_pat_softElectronBJetTag() { return cms2.jets_pat_softElectronBJetTag(); }
	vector<float> &jets_pat_softMuonBJetTag() { return cms2.jets_pat_softMuonBJetTag(); }
	vector<float> &jets_pat_softMuonNoIPBJetTag() { return cms2.jets_pat_softMuonNoIPBJetTag(); }
	vector<float> &jets_pat_trackCountingHighEffBJetTag() { return cms2.jets_pat_trackCountingHighEffBJetTag(); }
	vector<float> &jets_pat_trackCountingHighPurBJetTag() { return cms2.jets_pat_trackCountingHighPurBJetTag(); }
	vector<float> &jets_pat_udsCorrF() { return cms2.jets_pat_udsCorrF(); }
	vector<float> &mus_pat_caloIso() { return cms2.mus_pat_caloIso(); }
	vector<float> &mus_pat_ecalIso() { return cms2.mus_pat_ecalIso(); }
	vector<float> &mus_pat_ecalvetoDep() { return cms2.mus_pat_ecalvetoDep(); }
	vector<float> &mus_pat_hcalIso() { return cms2.mus_pat_hcalIso(); }
	vector<float> &mus_pat_hcalvetoDep() { return cms2.mus_pat_hcalvetoDep(); }
	vector<float> &mus_pat_trackIso() { return cms2.mus_pat_trackIso(); }
	vector<float> &mus_pat_vetoDep() { return cms2.mus_pat_vetoDep(); }
	vector<float> &taus_pf_bremsRecoveryEOverPLead() { return cms2.taus_pf_bremsRecoveryEOverPLead(); }
	vector<float> &taus_pf_caloComp() { return cms2.taus_pf_caloComp(); }
	vector<float> &taus_pf_ecalStripSumEOverPLead() { return cms2.taus_pf_ecalStripSumEOverPLead(); }
	vector<float> &taus_pf_emf() { return cms2.taus_pf_emf(); }
	vector<float> &taus_pf_hcal3x3OverPLead() { return cms2.taus_pf_hcal3x3OverPLead(); }
	vector<float> &taus_pf_hcalMaxOverPLead() { return cms2.taus_pf_hcalMaxOverPLead(); }
	vector<float> &taus_pf_hcalTotOverPLead() { return cms2.taus_pf_hcalTotOverPLead(); }
	vector<float> &taus_pf_isolationchargecandPtSum() { return cms2.taus_pf_isolationchargecandPtSum(); }
	vector<float> &taus_pf_isolationgammacandEtSum() { return cms2.taus_pf_isolationgammacandEtSum(); }
	vector<float> &taus_pf_lead_chargecand_Signed_Sipt() { return cms2.taus_pf_lead_chargecand_Signed_Sipt(); }
	vector<float> &taus_pf_leadtrk_chi2() { return cms2.taus_pf_leadtrk_chi2(); }
	vector<float> &taus_pf_leadtrk_d0() { return cms2.taus_pf_leadtrk_d0(); }
	vector<float> &taus_pf_leadtrk_lostHits() { return cms2.taus_pf_leadtrk_lostHits(); }
	vector<float> &taus_pf_leadtrk_ndof() { return cms2.taus_pf_leadtrk_ndof(); }
	vector<float> &taus_pf_leadtrk_validHits() { return cms2.taus_pf_leadtrk_validHits(); }
	vector<float> &taus_pf_leadtrk_z0() { return cms2.taus_pf_leadtrk_z0(); }
	vector<float> &taus_pf_maximumHCALPFClusterEt() { return cms2.taus_pf_maximumHCALPFClusterEt(); }
	vector<float> &taus_pf_segComp() { return cms2.taus_pf_segComp(); }
	vector<float> &photons_e1x5() { return cms2.photons_e1x5(); }
	vector<float> &photons_e2x5Max() { return cms2.photons_e2x5Max(); }
	vector<float> &photons_e3x3() { return cms2.photons_e3x3(); }
	vector<float> &photons_e5x5() { return cms2.photons_e5x5(); }
	vector<float> &photons_eMax() { return cms2.photons_eMax(); }
	vector<float> &photons_eSC() { return cms2.photons_eSC(); }
	vector<float> &photons_eSCPresh() { return cms2.photons_eSCPresh(); }
	vector<float> &photons_eSCRaw() { return cms2.photons_eSCRaw(); }
	vector<float> &photons_eSeed() { return cms2.photons_eSeed(); }
	vector<float> &photons_ecalIso() { return cms2.photons_ecalIso(); }
	vector<float> &photons_hOverE() { return cms2.photons_hOverE(); }
	vector<float> &photons_hcalIso() { return cms2.photons_hcalIso(); }
	vector<float> &photons_sigmaEtaEta() { return cms2.photons_sigmaEtaEta(); }
	vector<float> &photons_sigmaIEtaIEta() { return cms2.photons_sigmaIEtaIEta(); }
	vector<float> &photons_sigmaIPhiIPhi() { return cms2.photons_sigmaIPhiIPhi(); }
	vector<float> &photons_sigmaPhiPhi() { return cms2.photons_sigmaPhiPhi(); }
	vector<float> &photons_tkIso() { return cms2.photons_tkIso(); }
	vector<float> &scs_clustersSize() { return cms2.scs_clustersSize(); }
	vector<float> &scs_crystalsSize() { return cms2.scs_crystalsSize(); }
	vector<float> &scs_e1x3() { return cms2.scs_e1x3(); }
	vector<float> &scs_e1x5() { return cms2.scs_e1x5(); }
	vector<float> &scs_e2x2() { return cms2.scs_e2x2(); }
	vector<float> &scs_e2x5Max() { return cms2.scs_e2x5Max(); }
	vector<float> &scs_e3x1() { return cms2.scs_e3x1(); }
	vector<float> &scs_e3x2() { return cms2.scs_e3x2(); }
	vector<float> &scs_e3x3() { return cms2.scs_e3x3(); }
	vector<float> &scs_e4x4() { return cms2.scs_e4x4(); }
	vector<float> &scs_e5x5() { return cms2.scs_e5x5(); }
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
	vector<float> &mus_tcmet_deltax() { return cms2.mus_tcmet_deltax(); }
	vector<float> &mus_tcmet_deltay() { return cms2.mus_tcmet_deltay(); }
	vector<float> &trks_chi2() { return cms2.trks_chi2(); }
	vector<float> &trks_d0() { return cms2.trks_d0(); }
	vector<float> &trks_d0Err() { return cms2.trks_d0Err(); }
	vector<float> &trks_d0corr() { return cms2.trks_d0corr(); }
	vector<float> &trks_d0corrPhi() { return cms2.trks_d0corrPhi(); }
	vector<float> &trks_etaErr() { return cms2.trks_etaErr(); }
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
	int &evt_nl1emiso() { return cms2.evt_nl1emiso(); }
	int &evt_nl1emnoiso() { return cms2.evt_nl1emnoiso(); }
	int &evt_nl1jetsc() { return cms2.evt_nl1jetsc(); }
	int &evt_nl1jetsf() { return cms2.evt_nl1jetsf(); }
	int &evt_nl1jetst() { return cms2.evt_nl1jetst(); }
	int &evt_nl1mus() { return cms2.evt_nl1mus(); }
	int &pdfinfo_id1() { return cms2.pdfinfo_id1(); }
	int &pdfinfo_id2() { return cms2.pdfinfo_id2(); }
	vector<int> &taus_calo_charge() { return cms2.taus_calo_charge(); }
	vector<int> &taus_calo_iso_ntrks() { return cms2.taus_calo_iso_ntrks(); }
	vector<int> &taus_calo_sig_ntrks() { return cms2.taus_calo_sig_ntrks(); }
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
	vector<int> &els_closestJet() { return cms2.els_closestJet(); }
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
	vector<int> &hyp_ll_charge() { return cms2.hyp_ll_charge(); }
	vector<int> &hyp_ll_id() { return cms2.hyp_ll_id(); }
	vector<int> &hyp_ll_index() { return cms2.hyp_ll_index(); }
	vector<int> &hyp_ll_lostHits() { return cms2.hyp_ll_lostHits(); }
	vector<int> &hyp_ll_mc_id() { return cms2.hyp_ll_mc_id(); }
	vector<int> &hyp_ll_mc_motherid() { return cms2.hyp_ll_mc_motherid(); }
	vector<int> &hyp_ll_validHits() { return cms2.hyp_ll_validHits(); }
	vector<int> &hyp_lt_charge() { return cms2.hyp_lt_charge(); }
	vector<int> &hyp_lt_id() { return cms2.hyp_lt_id(); }
	vector<int> &hyp_lt_index() { return cms2.hyp_lt_index(); }
	vector<int> &hyp_lt_lostHits() { return cms2.hyp_lt_lostHits(); }
	vector<int> &hyp_lt_mc_id() { return cms2.hyp_lt_mc_id(); }
	vector<int> &hyp_lt_mc_motherid() { return cms2.hyp_lt_mc_motherid(); }
	vector<int> &hyp_lt_validHits() { return cms2.hyp_lt_validHits(); }
	vector<int> &hyp_njets() { return cms2.hyp_njets(); }
	vector<int> &hyp_nojets() { return cms2.hyp_nojets(); }
	vector<int> &hyp_type() { return cms2.hyp_type(); }
	vector<int> &hyp_quadlep_first_type() { return cms2.hyp_quadlep_first_type(); }
	vector<int> &hyp_quadlep_fourth_type() { return cms2.hyp_quadlep_fourth_type(); }
	vector<int> &hyp_quadlep_second_type() { return cms2.hyp_quadlep_second_type(); }
	vector<int> &hyp_quadlep_third_type() { return cms2.hyp_quadlep_third_type(); }
	vector<int> &hyp_trilep_first_type() { return cms2.hyp_trilep_first_type(); }
	vector<int> &hyp_trilep_second_type() { return cms2.hyp_trilep_second_type(); }
	vector<int> &hyp_trilep_third_type() { return cms2.hyp_trilep_third_type(); }
	vector<int> &jets_closestElectron() { return cms2.jets_closestElectron(); }
	vector<int> &jets_closestMuon() { return cms2.jets_closestMuon(); }
	vector<int> &l1emiso_ieta() { return cms2.l1emiso_ieta(); }
	vector<int> &l1emiso_iphi() { return cms2.l1emiso_iphi(); }
	vector<int> &l1emiso_rawId() { return cms2.l1emiso_rawId(); }
	vector<int> &l1emiso_type() { return cms2.l1emiso_type(); }
	vector<int> &l1emnoiso_ieta() { return cms2.l1emnoiso_ieta(); }
	vector<int> &l1emnoiso_iphi() { return cms2.l1emnoiso_iphi(); }
	vector<int> &l1emnoiso_rawId() { return cms2.l1emnoiso_rawId(); }
	vector<int> &l1emnoiso_type() { return cms2.l1emnoiso_type(); }
	vector<int> &l1jetsc_ieta() { return cms2.l1jetsc_ieta(); }
	vector<int> &l1jetsc_iphi() { return cms2.l1jetsc_iphi(); }
	vector<int> &l1jetsc_rawId() { return cms2.l1jetsc_rawId(); }
	vector<int> &l1jetsc_type() { return cms2.l1jetsc_type(); }
	vector<int> &l1jetsf_ieta() { return cms2.l1jetsf_ieta(); }
	vector<int> &l1jetsf_iphi() { return cms2.l1jetsf_iphi(); }
	vector<int> &l1jetsf_rawId() { return cms2.l1jetsf_rawId(); }
	vector<int> &l1jetsf_type() { return cms2.l1jetsf_type(); }
	vector<int> &l1jetst_ieta() { return cms2.l1jetst_ieta(); }
	vector<int> &l1jetst_iphi() { return cms2.l1jetst_iphi(); }
	vector<int> &l1jetst_rawId() { return cms2.l1jetst_rawId(); }
	vector<int> &l1jetst_type() { return cms2.l1jetst_type(); }
	vector<int> &l1mus_flags() { return cms2.l1mus_flags(); }
	vector<int> &l1mus_q() { return cms2.l1mus_q(); }
	vector<int> &l1mus_qual() { return cms2.l1mus_qual(); }
	vector<int> &l1mus_qualFlags() { return cms2.l1mus_qualFlags(); }
	vector<int> &mus_met_flag() { return cms2.mus_met_flag(); }
	vector<int> &mus_closestEle() { return cms2.mus_closestEle(); }
	vector<int> &mus_closestJet() { return cms2.mus_closestJet(); }
	vector<int> &mus_trkidx() { return cms2.mus_trkidx(); }
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
	vector<int> &els_pat_genID() { return cms2.els_pat_genID(); }
	vector<int> &els_pat_genMotherID() { return cms2.els_pat_genMotherID(); }
	vector<int> &jets_pat_genPartonMother_id() { return cms2.jets_pat_genPartonMother_id(); }
	vector<int> &jets_pat_genParton_id() { return cms2.jets_pat_genParton_id(); }
	vector<int> &jets_pat_partonFlavour() { return cms2.jets_pat_partonFlavour(); }
	vector<int> &mus_pat_genID() { return cms2.mus_pat_genID(); }
	vector<int> &mus_pat_genMotherID() { return cms2.mus_pat_genMotherID(); }
	vector<int> &taus_pf_charge() { return cms2.taus_pf_charge(); }
	vector<int> &taus_pf_electronPreID() { return cms2.taus_pf_electronPreID(); }
	vector<int> &taus_pf_hasMuonReference() { return cms2.taus_pf_hasMuonReference(); }
	vector<int> &taus_pf_iso_ncharge_cand() { return cms2.taus_pf_iso_ncharge_cand(); }
	vector<int> &taus_pf_iso_ngamma_cand() { return cms2.taus_pf_iso_ngamma_cand(); }
	vector<int> &taus_pf_iso_nneutr_cand() { return cms2.taus_pf_iso_nneutr_cand(); }
	vector<int> &taus_pf_muonPreID() { return cms2.taus_pf_muonPreID(); }
	vector<int> &taus_pf_sig_ncharge_cand() { return cms2.taus_pf_sig_ncharge_cand(); }
	vector<int> &taus_pf_sig_ngamma_cand() { return cms2.taus_pf_sig_ngamma_cand(); }
	vector<int> &taus_pf_sig_nneutr_cand() { return cms2.taus_pf_sig_nneutr_cand(); }
	vector<int> &mus_tcmet_flag() { return cms2.mus_tcmet_flag(); }
	vector<int> &trks_charge() { return cms2.trks_charge(); }
	vector<int> &trks_lostHits() { return cms2.trks_lostHits(); }
	vector<int> &trks_validHits() { return cms2.trks_validHits(); }
	vector<int> &trks_elsidx() { return cms2.trks_elsidx(); }
	vector<int> &trk_musidx() { return cms2.trk_musidx(); }
	vector<int> &hlt2ele10LWR_id() { return cms2.hlt2ele10LWR_id(); }
	vector<int> &hlt2ele10LWR_tid() { return cms2.hlt2ele10LWR_tid(); }
	vector<int> &hlt2mu3_id() { return cms2.hlt2mu3_id(); }
	vector<int> &hlt2mu3_tid() { return cms2.hlt2mu3_tid(); }
	vector<int> &hltLisoele18LWR_id() { return cms2.hltLisoele18LWR_id(); }
	vector<int> &hltLisoele18LWR_tid() { return cms2.hltLisoele18LWR_tid(); }
	vector<int> &hltisoele18R_id() { return cms2.hltisoele18R_id(); }
	vector<int> &hltisoele18R_tid() { return cms2.hltisoele18R_tid(); }
	vector<int> &hltjet30_id() { return cms2.hltjet30_id(); }
	vector<int> &hltjet30_tid() { return cms2.hltjet30_tid(); }
	vector<int> &hltl1jet15_id() { return cms2.hltl1jet15_id(); }
	vector<int> &hltl1jet15_tid() { return cms2.hltl1jet15_tid(); }
	vector<int> &hltl1met20_id() { return cms2.hltl1met20_id(); }
	vector<int> &hltl1met20_tid() { return cms2.hltl1met20_tid(); }
	vector<int> &hltl1mu_id() { return cms2.hltl1mu_id(); }
	vector<int> &hltl1mu_tid() { return cms2.hltl1mu_tid(); }
	vector<int> &hltl2mu9_id() { return cms2.hltl2mu9_id(); }
	vector<int> &hltl2mu9_tid() { return cms2.hltl2mu9_tid(); }
	vector<int> &hltmet25_id() { return cms2.hltmet25_id(); }
	vector<int> &hltmet25_tid() { return cms2.hltmet25_tid(); }
	vector<int> &hltmu11_id() { return cms2.hltmu11_id(); }
	vector<int> &hltmu11_tid() { return cms2.hltmu11_tid(); }
	vector<int> &hltmu9_id() { return cms2.hltmu9_id(); }
	vector<int> &hltmu9_tid() { return cms2.hltmu9_tid(); }
	vector<int> &vtxs_hasRefittedTracks() { return cms2.vtxs_hasRefittedTracks(); }
	vector<int> &vtxs_isFake() { return cms2.vtxs_isFake(); }
	vector<int> &vtxs_isValid() { return cms2.vtxs_isValid(); }
	vector<int> &vtxs_tracksSize() { return cms2.vtxs_tracksSize(); }
	vector<vector<int> > &hyp_jets_idx() { return cms2.hyp_jets_idx(); }
	vector<vector<int> > &hyp_other_jets_idx() { return cms2.hyp_other_jets_idx(); }
	unsigned int &evt_nels() { return cms2.evt_nels(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	unsigned int &evt_ngenjets() { return cms2.evt_ngenjets(); }
	unsigned int &evt_njets() { return cms2.evt_njets(); }
	unsigned int &evt_njpts() { return cms2.evt_njpts(); }
	unsigned int &evt_nphotons() { return cms2.evt_nphotons(); }
	unsigned int &evt_nscs() { return cms2.evt_nscs(); }
	unsigned int &evt_ntrkjets() { return cms2.evt_ntrkjets(); }
	unsigned int &evt_nvtxs() { return cms2.evt_nvtxs(); }
	vector<unsigned int> &hyp_quadlep_bucket() { return cms2.hyp_quadlep_bucket(); }
	vector<unsigned int> &hyp_quadlep_first_index() { return cms2.hyp_quadlep_first_index(); }
	vector<unsigned int> &hyp_quadlep_fourth_index() { return cms2.hyp_quadlep_fourth_index(); }
	vector<unsigned int> &hyp_quadlep_second_index() { return cms2.hyp_quadlep_second_index(); }
	vector<unsigned int> &hyp_quadlep_third_index() { return cms2.hyp_quadlep_third_index(); }
	vector<unsigned int> &hyp_trilep_bucket() { return cms2.hyp_trilep_bucket(); }
	vector<unsigned int> &hyp_trilep_first_index() { return cms2.hyp_trilep_first_index(); }
	vector<unsigned int> &hyp_trilep_second_index() { return cms2.hyp_trilep_second_index(); }
	vector<unsigned int> &hyp_trilep_third_index() { return cms2.hyp_trilep_third_index(); }
	vector<unsigned int> &els_pat_flag() { return cms2.els_pat_flag(); }
	vector<unsigned int> &jets_pat_flag() { return cms2.jets_pat_flag(); }
	vector<unsigned int> &mus_pat_flag() { return cms2.mus_pat_flag(); }
	float &evt_scale1fb() { return cms2.evt_scale1fb(); }
	int &evt_nEvts() { return cms2.evt_nEvts(); }
	float &evt_filt_eff() { return cms2.evt_filt_eff(); }
	static bool passHLTTrigger(TString trigName) { return cms2.passHLTTrigger(trigName); }
	static bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }
}
#endif
