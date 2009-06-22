// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include <vector> 

#define PARANOIA

using namespace std; 
class CMS2 { 
private: 
	 TH1F *samplehisto;
protected: 
	unsigned int index;
	TString	evt_dataset_;
	TBranch *evt_dataset_branch;
	bool evt_dataset_isLoaded;
	vector<TString>	evt_HLT_trigNames_;
	TBranch *evt_HLT_trigNames_branch;
	bool evt_HLT_trigNames_isLoaded;
	vector<TString>	evt_L1_trigNames_;
	TBranch *evt_L1_trigNames_branch;
	bool evt_L1_trigNames_isLoaded;
	double	genps_pthat_;
	TBranch *genps_pthat_branch;
	bool genps_pthat_isLoaded;
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >	scs_pos_;
	TBranch *scs_pos_branch;
	bool scs_pos_isLoaded;
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >	scs_vtx_;
	TBranch *scs_vtx_branch;
	bool scs_vtx_isLoaded;
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >	vtxs_position_;
	TBranch *vtxs_position_branch;
	bool vtxs_position_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 	evt_bsp4_;
	TBranch *evt_bsp4_branch;
	bool evt_bsp4_isLoaded;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 	l1met_p4_;
	TBranch *l1met_p4_branch;
	bool l1met_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	alltrkjets_p4_;
	TBranch *alltrkjets_p4_branch;
	bool alltrkjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_mc3_p4_;
	TBranch *els_mc3_p4_branch;
	bool els_mc3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_mc_p4_;
	TBranch *els_mc_p4_branch;
	bool els_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_mc_gp_p4_;
	TBranch *jets_mc_gp_p4_branch;
	bool jets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_mc_p4_;
	TBranch *jets_mc_p4_branch;
	bool jets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_mc3_p4_;
	TBranch *mus_mc3_p4_branch;
	bool mus_mc3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_mc_p4_;
	TBranch *mus_mc_p4_branch;
	bool mus_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trk_mc3p4_;
	TBranch *trk_mc3p4_branch;
	bool trk_mc3p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trk_mcp4_;
	TBranch *trk_mcp4_branch;
	bool trk_mcp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4_;
	TBranch *els_p4_branch;
	bool els_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4In_;
	TBranch *els_p4In_branch;
	bool els_p4In_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4Out_;
	TBranch *els_p4Out_branch;
	bool els_p4Out_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_trk_p4_;
	TBranch *els_trk_p4_branch;
	bool els_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_vertex_p4_;
	TBranch *els_vertex_p4_branch;
	bool els_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	genjets_p4_;
	TBranch *genjets_p4_branch;
	bool genjets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	genps_p4_;
	TBranch *genps_p4_branch;
	bool genps_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	genps_prod_vtx_;
	TBranch *genps_prod_vtx_branch;
	bool genps_prod_vtx_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_mc_p4_;
	TBranch *hyp_ll_mc_p4_branch;
	bool hyp_ll_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_p4_;
	TBranch *hyp_ll_p4_branch;
	bool hyp_ll_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_trk_p4_;
	TBranch *hyp_ll_trk_p4_branch;
	bool hyp_ll_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_mc_p4_;
	TBranch *hyp_lt_mc_p4_branch;
	bool hyp_lt_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_p4_;
	TBranch *hyp_lt_p4_branch;
	bool hyp_lt_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_trk_p4_;
	TBranch *hyp_lt_trk_p4_branch;
	bool hyp_lt_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_p4_;
	TBranch *hyp_p4_branch;
	bool hyp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_p4_;
	TBranch *jets_p4_branch;
	bool jets_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jpts_p4_;
	TBranch *jpts_p4_branch;
	bool jpts_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1emiso_p4_;
	TBranch *l1emiso_p4_branch;
	bool l1emiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1emnoiso_p4_;
	TBranch *l1emnoiso_p4_branch;
	bool l1emnoiso_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetsc_p4_;
	TBranch *l1jetsc_p4_branch;
	bool l1jetsc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetsf_p4_;
	TBranch *l1jetsf_p4_branch;
	bool l1jetsf_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetst_p4_;
	TBranch *l1jetst_p4_branch;
	bool l1jetst_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1mus_p4_;
	TBranch *l1mus_p4_branch;
	bool l1mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_trk_p4_;
	TBranch *mus_trk_p4_branch;
	bool mus_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_vertex_p4_;
	TBranch *mus_vertex_p4_branch;
	bool mus_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_pat_genMotherP4_;
	TBranch *els_pat_genMotherP4_branch;
	bool els_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_pat_genP4_;
	TBranch *els_pat_genP4_branch;
	bool els_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_pat_genJet_p4_;
	TBranch *jets_pat_genJet_p4_branch;
	bool jets_pat_genJet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_pat_genPartonMother_p4_;
	TBranch *jets_pat_genPartonMother_p4_branch;
	bool jets_pat_genPartonMother_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_pat_genParton_p4_;
	TBranch *jets_pat_genParton_p4_branch;
	bool jets_pat_genParton_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_pat_jet_p4_;
	TBranch *jets_pat_jet_p4_branch;
	bool jets_pat_jet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_pat_jet_uncorp4_;
	TBranch *jets_pat_jet_uncorp4_branch;
	bool jets_pat_jet_uncorp4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_pat_genMotherP4_;
	TBranch *mus_pat_genMotherP4_branch;
	bool mus_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_pat_genP4_;
	TBranch *mus_pat_genP4_branch;
	bool mus_pat_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	scs_p4_;
	TBranch *scs_p4_branch;
	bool scs_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trks_trk_p4_;
	TBranch *trks_trk_p4_branch;
	bool trks_trk_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trks_vertex_p4_;
	TBranch *trks_vertex_p4_branch;
	bool trks_vertex_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hlt2ele10LWR_p4_;
	TBranch *hlt2ele10LWR_p4_branch;
	bool hlt2ele10LWR_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hlt2mu3_p4_;
	TBranch *hlt2mu3_p4_branch;
	bool hlt2mu3_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltLisoele18LWR_p4_;
	TBranch *hltLisoele18LWR_p4_branch;
	bool hltLisoele18LWR_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltisoele18R_p4_;
	TBranch *hltisoele18R_p4_branch;
	bool hltisoele18R_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltjet30_p4_;
	TBranch *hltjet30_p4_branch;
	bool hltjet30_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltl1jet15_p4_;
	TBranch *hltl1jet15_p4_branch;
	bool hltl1jet15_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltl1met20_p4_;
	TBranch *hltl1met20_p4_branch;
	bool hltl1met20_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltl1mu_p4_;
	TBranch *hltl1mu_p4_branch;
	bool hltl1mu_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltl2mu9_p4_;
	TBranch *hltl2mu9_p4_branch;
	bool hltl2mu9_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltmet25_p4_;
	TBranch *hltmet25_p4_branch;
	bool hltmet25_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltmu11_p4_;
	TBranch *hltmu11_p4_branch;
	bool hltmu11_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hltmu9_p4_;
	TBranch *hltmu9_p4_branch;
	bool hltmu9_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trkjets_p4_;
	TBranch *trkjets_p4_branch;
	bool trkjets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_gp_p4_;
	TBranch *hyp_jets_mc_gp_p4_branch;
	bool hyp_jets_mc_gp_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_p4_;
	TBranch *hyp_jets_mc_p4_branch;
	bool hyp_jets_mc_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_p4_;
	TBranch *hyp_jets_p4_branch;
	bool hyp_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_pat_genJet_p4_;
	TBranch *hyp_jets_pat_genJet_p4_branch;
	bool hyp_jets_pat_genJet_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_pat_genPartonMother_p4_;
	TBranch *hyp_jets_pat_genPartonMother_p4_branch;
	bool hyp_jets_pat_genPartonMother_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_pat_genParton_p4_;
	TBranch *hyp_jets_pat_genParton_p4_branch;
	bool hyp_jets_pat_genParton_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_pat_jet_p4_;
	TBranch *hyp_jets_pat_jet_p4_branch;
	bool hyp_jets_pat_jet_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_gp_p4_;
	TBranch *hyp_other_jets_mc_gp_p4_branch;
	bool hyp_other_jets_mc_gp_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_p4_;
	TBranch *hyp_other_jets_mc_p4_branch;
	bool hyp_other_jets_mc_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_p4_;
	TBranch *hyp_other_jets_p4_branch;
	bool hyp_other_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_pat_genJet_p4_;
	TBranch *hyp_other_jets_pat_genJet_p4_branch;
	bool hyp_other_jets_pat_genJet_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_pat_genPartonMother_p4_;
	TBranch *hyp_other_jets_pat_genPartonMother_p4_branch;
	bool hyp_other_jets_pat_genPartonMother_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_pat_genParton_p4_;
	TBranch *hyp_other_jets_pat_genParton_p4_branch;
	bool hyp_other_jets_pat_genParton_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_pat_jet_p4_;
	TBranch *hyp_other_jets_pat_jet_p4_branch;
	bool hyp_other_jets_pat_jet_p4_isLoaded;
	vector<double>	jets_closestElectron_DR_;
	TBranch *jets_closestElectron_DR_branch;
	bool jets_closestElectron_DR_isLoaded;
	vector<double>	jets_closestMuon_DR_;
	TBranch *jets_closestMuon_DR_branch;
	bool jets_closestMuon_DR_isLoaded;
	float	evt_bs_dxdz_;
	TBranch *evt_bs_dxdz_branch;
	bool evt_bs_dxdz_isLoaded;
	float	evt_bs_dxdzErr_;
	TBranch *evt_bs_dxdzErr_branch;
	bool evt_bs_dxdzErr_isLoaded;
	float	evt_bs_dydz_;
	TBranch *evt_bs_dydz_branch;
	bool evt_bs_dydz_isLoaded;
	float	evt_bs_dydzErr_;
	TBranch *evt_bs_dydzErr_branch;
	bool evt_bs_dydzErr_isLoaded;
	float	evt_bs_sigmaZ_;
	TBranch *evt_bs_sigmaZ_branch;
	bool evt_bs_sigmaZ_isLoaded;
	float	evt_bs_sigmaZErr_;
	TBranch *evt_bs_sigmaZErr_branch;
	bool evt_bs_sigmaZErr_isLoaded;
	float	evt_bs_widthErr_;
	TBranch *evt_bs_widthErr_branch;
	bool evt_bs_widthErr_isLoaded;
	float	evt_bs_xErr_;
	TBranch *evt_bs_xErr_branch;
	bool evt_bs_xErr_isLoaded;
	float	evt_bs_yErr_;
	TBranch *evt_bs_yErr_branch;
	bool evt_bs_yErr_isLoaded;
	float	evt_bs_zErr_;
	TBranch *evt_bs_zErr_branch;
	bool evt_bs_zErr_isLoaded;
	float	gen_met_;
	TBranch *gen_met_branch;
	bool gen_met_isLoaded;
	float	gen_metPhi_;
	TBranch *gen_metPhi_branch;
	bool gen_metPhi_isLoaded;
	float	evt_bField_;
	TBranch *evt_bField_branch;
	bool evt_bField_isLoaded;
	float	evt_kfactor_;
	TBranch *evt_kfactor_branch;
	bool evt_kfactor_isLoaded;
	float	evt_weight_;
	TBranch *evt_weight_branch;
	bool evt_weight_isLoaded;
	float	evt_xsec_excl_;
	TBranch *evt_xsec_excl_branch;
	bool evt_xsec_excl_isLoaded;
	float	evt_xsec_incl_;
	TBranch *evt_xsec_incl_branch;
	bool evt_xsec_incl_isLoaded;
	float	l1met_etHad_;
	TBranch *l1met_etHad_branch;
	bool l1met_etHad_isLoaded;
	float	l1met_etTot_;
	TBranch *l1met_etTot_branch;
	bool l1met_etTot_isLoaded;
	float	l1met_met_;
	TBranch *l1met_met_branch;
	bool l1met_met_isLoaded;
	float	evt_met_;
	TBranch *evt_met_branch;
	bool evt_met_isLoaded;
	float	evt_metHO_;
	TBranch *evt_metHO_branch;
	bool evt_metHO_isLoaded;
	float	evt_metHOPhi_;
	TBranch *evt_metHOPhi_branch;
	bool evt_metHOPhi_isLoaded;
	float	evt_metHOSig_;
	TBranch *evt_metHOSig_branch;
	bool evt_metHOSig_isLoaded;
	float	evt_metMuonCorr_;
	TBranch *evt_metMuonCorr_branch;
	bool evt_metMuonCorr_isLoaded;
	float	evt_metMuonCorrPhi_;
	TBranch *evt_metMuonCorrPhi_branch;
	bool evt_metMuonCorrPhi_isLoaded;
	float	evt_metMuonCorrSig_;
	TBranch *evt_metMuonCorrSig_branch;
	bool evt_metMuonCorrSig_isLoaded;
	float	evt_metNoHF_;
	TBranch *evt_metNoHF_branch;
	bool evt_metNoHF_isLoaded;
	float	evt_metNoHFHO_;
	TBranch *evt_metNoHFHO_branch;
	bool evt_metNoHFHO_isLoaded;
	float	evt_metNoHFHOPhi_;
	TBranch *evt_metNoHFHOPhi_branch;
	bool evt_metNoHFHOPhi_isLoaded;
	float	evt_metNoHFHOSig_;
	TBranch *evt_metNoHFHOSig_branch;
	bool evt_metNoHFHOSig_isLoaded;
	float	evt_metNoHFPhi_;
	TBranch *evt_metNoHFPhi_branch;
	bool evt_metNoHFPhi_isLoaded;
	float	evt_metSig_;
	TBranch *evt_metSig_branch;
	bool evt_metSig_isLoaded;
	float	evt_metOpt_;
	TBranch *evt_metOpt_branch;
	bool evt_metOpt_isLoaded;
	float	evt_metOptHO_;
	TBranch *evt_metOptHO_branch;
	bool evt_metOptHO_isLoaded;
	float	evt_metOptHOPhi_;
	TBranch *evt_metOptHOPhi_branch;
	bool evt_metOptHOPhi_isLoaded;
	float	evt_metOptHOSig_;
	TBranch *evt_metOptHOSig_branch;
	bool evt_metOptHOSig_isLoaded;
	float	evt_metOptNoHF_;
	TBranch *evt_metOptNoHF_branch;
	bool evt_metOptNoHF_isLoaded;
	float	evt_metOptNoHFHO_;
	TBranch *evt_metOptNoHFHO_branch;
	bool evt_metOptNoHFHO_isLoaded;
	float	evt_metOptNoHFHOPhi_;
	TBranch *evt_metOptNoHFHOPhi_branch;
	bool evt_metOptNoHFHOPhi_isLoaded;
	float	evt_metOptNoHFHOSig_;
	TBranch *evt_metOptNoHFHOSig_branch;
	bool evt_metOptNoHFHOSig_isLoaded;
	float	evt_metOptNoHFPhi_;
	TBranch *evt_metOptNoHFPhi_branch;
	bool evt_metOptNoHFPhi_isLoaded;
	float	evt_metOptSig_;
	TBranch *evt_metOptSig_branch;
	bool evt_metOptSig_isLoaded;
	float	evt_metOptPhi_;
	TBranch *evt_metOptPhi_branch;
	bool evt_metOptPhi_isLoaded;
	float	evt_metPhi_;
	TBranch *evt_metPhi_branch;
	bool evt_metPhi_isLoaded;
	float	evt_sumet_;
	TBranch *evt_sumet_branch;
	bool evt_sumet_isLoaded;
	float	evt_sumetHO_;
	TBranch *evt_sumetHO_branch;
	bool evt_sumetHO_isLoaded;
	float	evt_sumetMuonCorr_;
	TBranch *evt_sumetMuonCorr_branch;
	bool evt_sumetMuonCorr_isLoaded;
	float	evt_sumetNoHF_;
	TBranch *evt_sumetNoHF_branch;
	bool evt_sumetNoHF_isLoaded;
	float	evt_sumetNoHFHO_;
	TBranch *evt_sumetNoHFHO_branch;
	bool evt_sumetNoHFHO_isLoaded;
	float	evt_sumetOpt_;
	TBranch *evt_sumetOpt_branch;
	bool evt_sumetOpt_isLoaded;
	float	evt_sumetOptHO_;
	TBranch *evt_sumetOptHO_branch;
	bool evt_sumetOptHO_isLoaded;
	float	evt_sumetOptNoHF_;
	TBranch *evt_sumetOptNoHF_branch;
	bool evt_sumetOptNoHF_isLoaded;
	float	evt_sumetOptNoHFHO_;
	TBranch *evt_sumetOptNoHFHO_branch;
	bool evt_sumetOptNoHFHO_isLoaded;
	float	met_pat_metCor_;
	TBranch *met_pat_metCor_branch;
	bool met_pat_metCor_isLoaded;
	float	met_pat_metPhiCor_;
	TBranch *met_pat_metPhiCor_branch;
	bool met_pat_metPhiCor_isLoaded;
	float	met_pat_metPhiUncor_;
	TBranch *met_pat_metPhiUncor_branch;
	bool met_pat_metPhiUncor_isLoaded;
	float	met_pat_metPhiUncorJES_;
	TBranch *met_pat_metPhiUncorJES_branch;
	bool met_pat_metPhiUncorJES_isLoaded;
	float	met_pat_metPhiUncorMuon_;
	TBranch *met_pat_metPhiUncorMuon_branch;
	bool met_pat_metPhiUncorMuon_isLoaded;
	float	met_pat_metUncor_;
	TBranch *met_pat_metUncor_branch;
	bool met_pat_metUncor_isLoaded;
	float	met_pat_metUncorJES_;
	TBranch *met_pat_metUncorJES_branch;
	bool met_pat_metUncorJES_isLoaded;
	float	met_pat_metUncorMuon_;
	TBranch *met_pat_metUncorMuon_branch;
	bool met_pat_metUncorMuon_isLoaded;
	float	evt_tcmet_;
	TBranch *evt_tcmet_branch;
	bool evt_tcmet_isLoaded;
	float	evt_tcmetPhi_;
	TBranch *evt_tcmetPhi_branch;
	bool evt_tcmetPhi_isLoaded;
	float	evt_tcsumet_;
	TBranch *evt_tcsumet_branch;
	bool evt_tcsumet_isLoaded;
	vector<float>	els_mc3dr_;
	TBranch *els_mc3dr_branch;
	bool els_mc3dr_isLoaded;
	vector<float>	els_mcdr_;
	TBranch *els_mcdr_branch;
	bool els_mcdr_isLoaded;
	vector<float>	jets_mcdr_;
	TBranch *jets_mcdr_branch;
	bool jets_mcdr_isLoaded;
	vector<float>	jets_mc_emEnergy_;
	TBranch *jets_mc_emEnergy_branch;
	bool jets_mc_emEnergy_isLoaded;
	vector<float>	jets_mc_gpdr_;
	TBranch *jets_mc_gpdr_branch;
	bool jets_mc_gpdr_isLoaded;
	vector<float>	jets_mc_hadEnergy_;
	TBranch *jets_mc_hadEnergy_branch;
	bool jets_mc_hadEnergy_isLoaded;
	vector<float>	jets_mc_invEnergy_;
	TBranch *jets_mc_invEnergy_branch;
	bool jets_mc_invEnergy_isLoaded;
	vector<float>	jets_mc_otherEnergy_;
	TBranch *jets_mc_otherEnergy_branch;
	bool jets_mc_otherEnergy_isLoaded;
	vector<float>	mus_mc3dr_;
	TBranch *mus_mc3dr_branch;
	bool mus_mc3dr_isLoaded;
	vector<float>	mus_mcdr_;
	TBranch *mus_mcdr_branch;
	bool mus_mcdr_isLoaded;
	vector<float>	trk_mc3dr_;
	TBranch *trk_mc3dr_branch;
	bool trk_mc3dr_isLoaded;
	vector<float>	trk_mcdr_;
	TBranch *trk_mcdr_branch;
	bool trk_mcdr_isLoaded;
	vector<float>	els_conv_dcot_;
	TBranch *els_conv_dcot_branch;
	bool els_conv_dcot_isLoaded;
	vector<float>	els_conv_dist_;
	TBranch *els_conv_dist_branch;
	bool els_conv_dist_isLoaded;
	vector<float>	trks_conv_dcot_;
	TBranch *trks_conv_dcot_branch;
	bool trks_conv_dcot_isLoaded;
	vector<float>	trks_conv_dist_;
	TBranch *trks_conv_dist_branch;
	bool trks_conv_dist_isLoaded;
	vector<float>	els_ecalJuraIso_;
	TBranch *els_ecalJuraIso_branch;
	bool els_ecalJuraIso_isLoaded;
	vector<float>	els_ecalJuraTowerIso_;
	TBranch *els_ecalJuraTowerIso_branch;
	bool els_ecalJuraTowerIso_isLoaded;
	vector<float>	els_hcalConeIso_;
	TBranch *els_hcalConeIso_branch;
	bool els_hcalConeIso_isLoaded;
	vector<float>	els_jetdr_;
	TBranch *els_jetdr_branch;
	bool els_jetdr_isLoaded;
	vector<float>	els_musdr_;
	TBranch *els_musdr_branch;
	bool els_musdr_isLoaded;
	vector<float>	els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float>	els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
	vector<float>	els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float>	els_d0_;
	TBranch *els_d0_branch;
	bool els_d0_isLoaded;
	vector<float>	els_d0Err_;
	TBranch *els_d0Err_branch;
	bool els_d0Err_isLoaded;
	vector<float>	els_d0corr_;
	TBranch *els_d0corr_branch;
	bool els_d0corr_isLoaded;
	vector<float>	els_dEtaIn_;
	TBranch *els_dEtaIn_branch;
	bool els_dEtaIn_isLoaded;
	vector<float>	els_dEtaOut_;
	TBranch *els_dEtaOut_branch;
	bool els_dEtaOut_isLoaded;
	vector<float>	els_dPhiIn_;
	TBranch *els_dPhiIn_branch;
	bool els_dPhiIn_isLoaded;
	vector<float>	els_dPhiInPhiOut_;
	TBranch *els_dPhiInPhiOut_branch;
	bool els_dPhiInPhiOut_isLoaded;
	vector<float>	els_dPhiOut_;
	TBranch *els_dPhiOut_branch;
	bool els_dPhiOut_isLoaded;
	vector<float>	els_e2x5Max_;
	TBranch *els_e2x5Max_branch;
	bool els_e2x5Max_isLoaded;
	vector<float>	els_e3x3_;
	TBranch *els_e3x3_branch;
	bool els_e3x3_isLoaded;
	vector<float>	els_e5x5_;
	TBranch *els_e5x5_branch;
	bool els_e5x5_isLoaded;
	vector<float>	els_eMax_;
	TBranch *els_eMax_branch;
	bool els_eMax_isLoaded;
	vector<float>	els_eOverPIn_;
	TBranch *els_eOverPIn_branch;
	bool els_eOverPIn_isLoaded;
	vector<float>	els_eSC_;
	TBranch *els_eSC_branch;
	bool els_eSC_isLoaded;
	vector<float>	els_eSCPresh_;
	TBranch *els_eSCPresh_branch;
	bool els_eSCPresh_isLoaded;
	vector<float>	els_eSCRaw_;
	TBranch *els_eSCRaw_branch;
	bool els_eSCRaw_isLoaded;
	vector<float>	els_eSeed_;
	TBranch *els_eSeed_branch;
	bool els_eSeed_isLoaded;
	vector<float>	els_eSeedOverPOut_;
	TBranch *els_eSeedOverPOut_branch;
	bool els_eSeedOverPOut_isLoaded;
	vector<float>	els_etaErr_;
	TBranch *els_etaErr_branch;
	bool els_etaErr_isLoaded;
	vector<float>	els_fBrem_;
	TBranch *els_fBrem_branch;
	bool els_fBrem_isLoaded;
	vector<float>	els_hOverE_;
	TBranch *els_hOverE_branch;
	bool els_hOverE_isLoaded;
	vector<float>	els_ndof_;
	TBranch *els_ndof_branch;
	bool els_ndof_isLoaded;
	vector<float>	els_outerEta_;
	TBranch *els_outerEta_branch;
	bool els_outerEta_isLoaded;
	vector<float>	els_outerPhi_;
	TBranch *els_outerPhi_branch;
	bool els_outerPhi_isLoaded;
	vector<float>	els_phiErr_;
	TBranch *els_phiErr_branch;
	bool els_phiErr_isLoaded;
	vector<float>	els_ptErr_;
	TBranch *els_ptErr_branch;
	bool els_ptErr_isLoaded;
	vector<float>	els_sigmaEtaEta_;
	TBranch *els_sigmaEtaEta_branch;
	bool els_sigmaEtaEta_isLoaded;
	vector<float>	els_sigmaIEtaIEta_;
	TBranch *els_sigmaIEtaIEta_branch;
	bool els_sigmaIEtaIEta_isLoaded;
	vector<float>	els_sigmaIPhiIPhi_;
	TBranch *els_sigmaIPhiIPhi_branch;
	bool els_sigmaIPhiIPhi_isLoaded;
	vector<float>	els_sigmaPhiPhi_;
	TBranch *els_sigmaPhiPhi_branch;
	bool els_sigmaPhiPhi_isLoaded;
	vector<float>	els_tkIso_;
	TBranch *els_tkIso_branch;
	bool els_tkIso_isLoaded;
	vector<float>	els_vertexphi_;
	TBranch *els_vertexphi_branch;
	bool els_vertexphi_isLoaded;
	vector<float>	els_z0_;
	TBranch *els_z0_branch;
	bool els_z0_isLoaded;
	vector<float>	els_z0Err_;
	TBranch *els_z0Err_branch;
	bool els_z0Err_isLoaded;
	vector<float>	els_z0corr_;
	TBranch *els_z0corr_branch;
	bool els_z0corr_isLoaded;
	vector<float>	hyp_ll_chi2_;
	TBranch *hyp_ll_chi2_branch;
	bool hyp_ll_chi2_isLoaded;
	vector<float>	hyp_ll_d0_;
	TBranch *hyp_ll_d0_branch;
	bool hyp_ll_d0_isLoaded;
	vector<float>	hyp_ll_d0Err_;
	TBranch *hyp_ll_d0Err_branch;
	bool hyp_ll_d0Err_isLoaded;
	vector<float>	hyp_ll_d0corr_;
	TBranch *hyp_ll_d0corr_branch;
	bool hyp_ll_d0corr_isLoaded;
	vector<float>	hyp_ll_etaErr_;
	TBranch *hyp_ll_etaErr_branch;
	bool hyp_ll_etaErr_isLoaded;
	vector<float>	hyp_ll_iso_;
	TBranch *hyp_ll_iso_branch;
	bool hyp_ll_iso_isLoaded;
	vector<float>	hyp_ll_ndof_;
	TBranch *hyp_ll_ndof_branch;
	bool hyp_ll_ndof_isLoaded;
	vector<float>	hyp_ll_outerEta_;
	TBranch *hyp_ll_outerEta_branch;
	bool hyp_ll_outerEta_isLoaded;
	vector<float>	hyp_ll_outerPhi_;
	TBranch *hyp_ll_outerPhi_branch;
	bool hyp_ll_outerPhi_isLoaded;
	vector<float>	hyp_ll_phiErr_;
	TBranch *hyp_ll_phiErr_branch;
	bool hyp_ll_phiErr_isLoaded;
	vector<float>	hyp_ll_ptErr_;
	TBranch *hyp_ll_ptErr_branch;
	bool hyp_ll_ptErr_isLoaded;
	vector<float>	hyp_ll_tkIso_;
	TBranch *hyp_ll_tkIso_branch;
	bool hyp_ll_tkIso_isLoaded;
	vector<float>	hyp_ll_vertexphi_;
	TBranch *hyp_ll_vertexphi_branch;
	bool hyp_ll_vertexphi_isLoaded;
	vector<float>	hyp_ll_z0_;
	TBranch *hyp_ll_z0_branch;
	bool hyp_ll_z0_isLoaded;
	vector<float>	hyp_ll_z0Err_;
	TBranch *hyp_ll_z0Err_branch;
	bool hyp_ll_z0Err_isLoaded;
	vector<float>	hyp_ll_z0corr_;
	TBranch *hyp_ll_z0corr_branch;
	bool hyp_ll_z0corr_isLoaded;
	vector<float>	hyp_lt_chi2_;
	TBranch *hyp_lt_chi2_branch;
	bool hyp_lt_chi2_isLoaded;
	vector<float>	hyp_lt_d0_;
	TBranch *hyp_lt_d0_branch;
	bool hyp_lt_d0_isLoaded;
	vector<float>	hyp_lt_d0Err_;
	TBranch *hyp_lt_d0Err_branch;
	bool hyp_lt_d0Err_isLoaded;
	vector<float>	hyp_lt_d0corr_;
	TBranch *hyp_lt_d0corr_branch;
	bool hyp_lt_d0corr_isLoaded;
	vector<float>	hyp_lt_etaErr_;
	TBranch *hyp_lt_etaErr_branch;
	bool hyp_lt_etaErr_isLoaded;
	vector<float>	hyp_lt_iso_;
	TBranch *hyp_lt_iso_branch;
	bool hyp_lt_iso_isLoaded;
	vector<float>	hyp_lt_ndof_;
	TBranch *hyp_lt_ndof_branch;
	bool hyp_lt_ndof_isLoaded;
	vector<float>	hyp_lt_outerEta_;
	TBranch *hyp_lt_outerEta_branch;
	bool hyp_lt_outerEta_isLoaded;
	vector<float>	hyp_lt_outerPhi_;
	TBranch *hyp_lt_outerPhi_branch;
	bool hyp_lt_outerPhi_isLoaded;
	vector<float>	hyp_lt_phiErr_;
	TBranch *hyp_lt_phiErr_branch;
	bool hyp_lt_phiErr_isLoaded;
	vector<float>	hyp_lt_ptErr_;
	TBranch *hyp_lt_ptErr_branch;
	bool hyp_lt_ptErr_isLoaded;
	vector<float>	hyp_lt_tkIso_;
	TBranch *hyp_lt_tkIso_branch;
	bool hyp_lt_tkIso_isLoaded;
	vector<float>	hyp_lt_vertexphi_;
	TBranch *hyp_lt_vertexphi_branch;
	bool hyp_lt_vertexphi_isLoaded;
	vector<float>	hyp_lt_z0_;
	TBranch *hyp_lt_z0_branch;
	bool hyp_lt_z0_isLoaded;
	vector<float>	hyp_lt_z0Err_;
	TBranch *hyp_lt_z0Err_branch;
	bool hyp_lt_z0Err_isLoaded;
	vector<float>	hyp_lt_z0corr_;
	TBranch *hyp_lt_z0corr_branch;
	bool hyp_lt_z0corr_isLoaded;
	vector<float>	hyp_met_;
	TBranch *hyp_met_branch;
	bool hyp_met_isLoaded;
	vector<float>	hyp_metAll_;
	TBranch *hyp_metAll_branch;
	bool hyp_metAll_isLoaded;
	vector<float>	hyp_metAllCaloExp_;
	TBranch *hyp_metAllCaloExp_branch;
	bool hyp_metAllCaloExp_isLoaded;
	vector<float>	hyp_metCaloExp_;
	TBranch *hyp_metCaloExp_branch;
	bool hyp_metCaloExp_isLoaded;
	vector<float>	hyp_metCone_;
	TBranch *hyp_metCone_branch;
	bool hyp_metCone_isLoaded;
	vector<float>	hyp_metDPhiJet10_;
	TBranch *hyp_metDPhiJet10_branch;
	bool hyp_metDPhiJet10_isLoaded;
	vector<float>	hyp_metDPhiJet15_;
	TBranch *hyp_metDPhiJet15_branch;
	bool hyp_metDPhiJet15_isLoaded;
	vector<float>	hyp_metDPhiJet20_;
	TBranch *hyp_metDPhiJet20_branch;
	bool hyp_metDPhiJet20_isLoaded;
	vector<float>	hyp_metDPhiTrk10_;
	TBranch *hyp_metDPhiTrk10_branch;
	bool hyp_metDPhiTrk10_isLoaded;
	vector<float>	hyp_metDPhiTrk25_;
	TBranch *hyp_metDPhiTrk25_branch;
	bool hyp_metDPhiTrk25_isLoaded;
	vector<float>	hyp_metDPhiTrk50_;
	TBranch *hyp_metDPhiTrk50_branch;
	bool hyp_metDPhiTrk50_isLoaded;
	vector<float>	hyp_metJes10_;
	TBranch *hyp_metJes10_branch;
	bool hyp_metJes10_isLoaded;
	vector<float>	hyp_metJes15_;
	TBranch *hyp_metJes15_branch;
	bool hyp_metJes15_isLoaded;
	vector<float>	hyp_metJes30_;
	TBranch *hyp_metJes30_branch;
	bool hyp_metJes30_isLoaded;
	vector<float>	hyp_metJes5_;
	TBranch *hyp_metJes5_branch;
	bool hyp_metJes5_isLoaded;
	vector<float>	hyp_metJes50_;
	TBranch *hyp_metJes50_branch;
	bool hyp_metJes50_isLoaded;
	vector<float>	hyp_metNoCalo_;
	TBranch *hyp_metNoCalo_branch;
	bool hyp_metNoCalo_isLoaded;
	vector<float>	hyp_metPhi_;
	TBranch *hyp_metPhi_branch;
	bool hyp_metPhi_isLoaded;
	vector<float>	hyp_metPhiAll_;
	TBranch *hyp_metPhiAll_branch;
	bool hyp_metPhiAll_isLoaded;
	vector<float>	hyp_metPhiAllCaloExp_;
	TBranch *hyp_metPhiAllCaloExp_branch;
	bool hyp_metPhiAllCaloExp_isLoaded;
	vector<float>	hyp_metPhiCaloExp_;
	TBranch *hyp_metPhiCaloExp_branch;
	bool hyp_metPhiCaloExp_isLoaded;
	vector<float>	hyp_metPhiCone_;
	TBranch *hyp_metPhiCone_branch;
	bool hyp_metPhiCone_isLoaded;
	vector<float>	hyp_metPhiJes10_;
	TBranch *hyp_metPhiJes10_branch;
	bool hyp_metPhiJes10_isLoaded;
	vector<float>	hyp_metPhiJes15_;
	TBranch *hyp_metPhiJes15_branch;
	bool hyp_metPhiJes15_isLoaded;
	vector<float>	hyp_metPhiJes30_;
	TBranch *hyp_metPhiJes30_branch;
	bool hyp_metPhiJes30_isLoaded;
	vector<float>	hyp_metPhiJes5_;
	TBranch *hyp_metPhiJes5_branch;
	bool hyp_metPhiJes5_isLoaded;
	vector<float>	hyp_metPhiJes50_;
	TBranch *hyp_metPhiJes50_branch;
	bool hyp_metPhiJes50_isLoaded;
	vector<float>	hyp_metPhiNoCalo_;
	TBranch *hyp_metPhiNoCalo_branch;
	bool hyp_metPhiNoCalo_isLoaded;
	vector<float>	hyp_quadlep_met_;
	TBranch *hyp_quadlep_met_branch;
	bool hyp_quadlep_met_isLoaded;
	vector<float>	hyp_quadlep_metAll_;
	TBranch *hyp_quadlep_metAll_branch;
	bool hyp_quadlep_metAll_isLoaded;
	vector<float>	hyp_trilep_met_;
	TBranch *hyp_trilep_met_branch;
	bool hyp_trilep_met_isLoaded;
	vector<float>	hyp_trilep_metAll_;
	TBranch *hyp_trilep_metAll_branch;
	bool hyp_trilep_metAll_isLoaded;
	vector<float>	jets_EMFcor_;
	TBranch *jets_EMFcor_branch;
	bool jets_EMFcor_isLoaded;
	vector<float>	jets_chFrac_;
	TBranch *jets_chFrac_branch;
	bool jets_chFrac_isLoaded;
	vector<float>	jets_cor_;
	TBranch *jets_cor_branch;
	bool jets_cor_isLoaded;
	vector<float>	jets_emFrac_;
	TBranch *jets_emFrac_branch;
	bool jets_emFrac_isLoaded;
	vector<float>	jpts_chFrac_;
	TBranch *jpts_chFrac_branch;
	bool jpts_chFrac_isLoaded;
	vector<float>	jpts_cor_;
	TBranch *jpts_cor_branch;
	bool jpts_cor_isLoaded;
	vector<float>	jpts_emFrac_;
	TBranch *jpts_emFrac_branch;
	bool jpts_emFrac_isLoaded;
	vector<float>	mus_eledr_;
	TBranch *mus_eledr_branch;
	bool mus_eledr_isLoaded;
	vector<float>	mus_jetdr_;
	TBranch *mus_jetdr_branch;
	bool mus_jetdr_isLoaded;
	vector<float>	mus_trkdr_;
	TBranch *mus_trkdr_branch;
	bool mus_trkdr_isLoaded;
	vector<float>	mus_caloCompatibility_;
	TBranch *mus_caloCompatibility_branch;
	bool mus_caloCompatibility_isLoaded;
	vector<float>	mus_chi2_;
	TBranch *mus_chi2_branch;
	bool mus_chi2_isLoaded;
	vector<float>	mus_d0_;
	TBranch *mus_d0_branch;
	bool mus_d0_isLoaded;
	vector<float>	mus_d0Err_;
	TBranch *mus_d0Err_branch;
	bool mus_d0Err_isLoaded;
	vector<float>	mus_d0corr_;
	TBranch *mus_d0corr_branch;
	bool mus_d0corr_isLoaded;
	vector<float>	mus_e_em_;
	TBranch *mus_e_em_branch;
	bool mus_e_em_isLoaded;
	vector<float>	mus_e_emS9_;
	TBranch *mus_e_emS9_branch;
	bool mus_e_emS9_isLoaded;
	vector<float>	mus_e_had_;
	TBranch *mus_e_had_branch;
	bool mus_e_had_isLoaded;
	vector<float>	mus_e_hadS9_;
	TBranch *mus_e_hadS9_branch;
	bool mus_e_hadS9_isLoaded;
	vector<float>	mus_e_ho_;
	TBranch *mus_e_ho_branch;
	bool mus_e_ho_isLoaded;
	vector<float>	mus_e_hoS9_;
	TBranch *mus_e_hoS9_branch;
	bool mus_e_hoS9_isLoaded;
	vector<float>	mus_etaErr_;
	TBranch *mus_etaErr_branch;
	bool mus_etaErr_isLoaded;
	vector<float>	mus_gfit_chi2_;
	TBranch *mus_gfit_chi2_branch;
	bool mus_gfit_chi2_isLoaded;
	vector<float>	mus_gfit_ndof_;
	TBranch *mus_gfit_ndof_branch;
	bool mus_gfit_ndof_isLoaded;
	vector<float>	mus_iso_;
	TBranch *mus_iso_branch;
	bool mus_iso_isLoaded;
	vector<float>	mus_iso03_emEt_;
	TBranch *mus_iso03_emEt_branch;
	bool mus_iso03_emEt_isLoaded;
	vector<float>	mus_iso03_hadEt_;
	TBranch *mus_iso03_hadEt_branch;
	bool mus_iso03_hadEt_isLoaded;
	vector<float>	mus_iso03_hoEt_;
	TBranch *mus_iso03_hoEt_branch;
	bool mus_iso03_hoEt_isLoaded;
	vector<float>	mus_iso03_sumPt_;
	TBranch *mus_iso03_sumPt_branch;
	bool mus_iso03_sumPt_isLoaded;
	vector<float>	mus_iso05_emEt_;
	TBranch *mus_iso05_emEt_branch;
	bool mus_iso05_emEt_isLoaded;
	vector<float>	mus_iso05_hadEt_;
	TBranch *mus_iso05_hadEt_branch;
	bool mus_iso05_hadEt_isLoaded;
	vector<float>	mus_iso05_hoEt_;
	TBranch *mus_iso05_hoEt_branch;
	bool mus_iso05_hoEt_isLoaded;
	vector<float>	mus_iso05_sumPt_;
	TBranch *mus_iso05_sumPt_branch;
	bool mus_iso05_sumPt_isLoaded;
	vector<float>	mus_ndof_;
	TBranch *mus_ndof_branch;
	bool mus_ndof_isLoaded;
	vector<float>	mus_outerEta_;
	TBranch *mus_outerEta_branch;
	bool mus_outerEta_isLoaded;
	vector<float>	mus_outerPhi_;
	TBranch *mus_outerPhi_branch;
	bool mus_outerPhi_isLoaded;
	vector<float>	mus_phiErr_;
	TBranch *mus_phiErr_branch;
	bool mus_phiErr_isLoaded;
	vector<float>	mus_ptErr_;
	TBranch *mus_ptErr_branch;
	bool mus_ptErr_isLoaded;
	vector<float>	mus_vertexphi_;
	TBranch *mus_vertexphi_branch;
	bool mus_vertexphi_isLoaded;
	vector<float>	mus_z0_;
	TBranch *mus_z0_branch;
	bool mus_z0_isLoaded;
	vector<float>	mus_z0Err_;
	TBranch *mus_z0Err_branch;
	bool mus_z0Err_isLoaded;
	vector<float>	mus_z0corr_;
	TBranch *mus_z0corr_branch;
	bool mus_z0corr_isLoaded;
	vector<float>	els_pat_caloIso_;
	TBranch *els_pat_caloIso_branch;
	bool els_pat_caloIso_isLoaded;
	vector<float>	els_pat_ecalIso_;
	TBranch *els_pat_ecalIso_branch;
	bool els_pat_ecalIso_isLoaded;
	vector<float>	els_pat_hcalIso_;
	TBranch *els_pat_hcalIso_branch;
	bool els_pat_hcalIso_isLoaded;
	vector<float>	els_pat_looseId_;
	TBranch *els_pat_looseId_branch;
	bool els_pat_looseId_isLoaded;
	vector<float>	els_pat_robustHighEnergy_;
	TBranch *els_pat_robustHighEnergy_branch;
	bool els_pat_robustHighEnergy_isLoaded;
	vector<float>	els_pat_robustLooseId_;
	TBranch *els_pat_robustLooseId_branch;
	bool els_pat_robustLooseId_isLoaded;
	vector<float>	els_pat_robustTightId_;
	TBranch *els_pat_robustTightId_branch;
	bool els_pat_robustTightId_isLoaded;
	vector<float>	els_pat_scE1x5_;
	TBranch *els_pat_scE1x5_branch;
	bool els_pat_scE1x5_isLoaded;
	vector<float>	els_pat_scE2x5Max_;
	TBranch *els_pat_scE2x5Max_branch;
	bool els_pat_scE2x5Max_isLoaded;
	vector<float>	els_pat_scE5x5_;
	TBranch *els_pat_scE5x5_branch;
	bool els_pat_scE5x5_isLoaded;
	vector<float>	els_pat_sigmaEtaEta_;
	TBranch *els_pat_sigmaEtaEta_branch;
	bool els_pat_sigmaEtaEta_isLoaded;
	vector<float>	els_pat_sigmaIEtaIEta_;
	TBranch *els_pat_sigmaIEtaIEta_branch;
	bool els_pat_sigmaIEtaIEta_isLoaded;
	vector<float>	els_pat_tightId_;
	TBranch *els_pat_tightId_branch;
	bool els_pat_tightId_isLoaded;
	vector<float>	els_pat_trackIso_;
	TBranch *els_pat_trackIso_branch;
	bool els_pat_trackIso_isLoaded;
	vector<float>	jets_pat_bCorrF_;
	TBranch *jets_pat_bCorrF_branch;
	bool jets_pat_bCorrF_isLoaded;
	vector<float>	jets_pat_cCorrF_;
	TBranch *jets_pat_cCorrF_branch;
	bool jets_pat_cCorrF_isLoaded;
	vector<float>	jets_pat_combinedSecondaryVertexBJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexBJetTag_branch;
	bool jets_pat_combinedSecondaryVertexBJetTag_isLoaded;
	vector<float>	jets_pat_combinedSecondaryVertexMVABJetTag_;
	TBranch *jets_pat_combinedSecondaryVertexMVABJetTag_branch;
	bool jets_pat_combinedSecondaryVertexMVABJetTag_isLoaded;
	vector<float>	jets_pat_coneIsolationTauJetTag_;
	TBranch *jets_pat_coneIsolationTauJetTag_branch;
	bool jets_pat_coneIsolationTauJetTag_isLoaded;
	vector<float>	jets_pat_gluCorrF_;
	TBranch *jets_pat_gluCorrF_branch;
	bool jets_pat_gluCorrF_isLoaded;
	vector<float>	jets_pat_impactParameterMVABJetTag_;
	TBranch *jets_pat_impactParameterMVABJetTag_branch;
	bool jets_pat_impactParameterMVABJetTag_isLoaded;
	vector<float>	jets_pat_jetBProbabilityBJetTag_;
	TBranch *jets_pat_jetBProbabilityBJetTag_branch;
	bool jets_pat_jetBProbabilityBJetTag_isLoaded;
	vector<float>	jets_pat_jetCharge_;
	TBranch *jets_pat_jetCharge_branch;
	bool jets_pat_jetCharge_isLoaded;
	vector<float>	jets_pat_jetProbabilityBJetTag_;
	TBranch *jets_pat_jetProbabilityBJetTag_branch;
	bool jets_pat_jetProbabilityBJetTag_isLoaded;
	vector<float>	jets_pat_noCorrF_;
	TBranch *jets_pat_noCorrF_branch;
	bool jets_pat_noCorrF_isLoaded;
	vector<float>	jets_pat_simpleSecondaryVertexBJetTag_;
	TBranch *jets_pat_simpleSecondaryVertexBJetTag_branch;
	bool jets_pat_simpleSecondaryVertexBJetTag_isLoaded;
	vector<float>	jets_pat_softElectronBJetTag_;
	TBranch *jets_pat_softElectronBJetTag_branch;
	bool jets_pat_softElectronBJetTag_isLoaded;
	vector<float>	jets_pat_softMuonBJetTag_;
	TBranch *jets_pat_softMuonBJetTag_branch;
	bool jets_pat_softMuonBJetTag_isLoaded;
	vector<float>	jets_pat_softMuonNoIPBJetTag_;
	TBranch *jets_pat_softMuonNoIPBJetTag_branch;
	bool jets_pat_softMuonNoIPBJetTag_isLoaded;
	vector<float>	jets_pat_trackCountingHighEffBJetTag_;
	TBranch *jets_pat_trackCountingHighEffBJetTag_branch;
	bool jets_pat_trackCountingHighEffBJetTag_isLoaded;
	vector<float>	jets_pat_trackCountingHighPurBJetTag_;
	TBranch *jets_pat_trackCountingHighPurBJetTag_branch;
	bool jets_pat_trackCountingHighPurBJetTag_isLoaded;
	vector<float>	jets_pat_udsCorrF_;
	TBranch *jets_pat_udsCorrF_branch;
	bool jets_pat_udsCorrF_isLoaded;
	vector<float>	mus_pat_caloIso_;
	TBranch *mus_pat_caloIso_branch;
	bool mus_pat_caloIso_isLoaded;
	vector<float>	mus_pat_ecalIso_;
	TBranch *mus_pat_ecalIso_branch;
	bool mus_pat_ecalIso_isLoaded;
	vector<float>	mus_pat_ecalvetoDep_;
	TBranch *mus_pat_ecalvetoDep_branch;
	bool mus_pat_ecalvetoDep_isLoaded;
	vector<float>	mus_pat_hcalIso_;
	TBranch *mus_pat_hcalIso_branch;
	bool mus_pat_hcalIso_isLoaded;
	vector<float>	mus_pat_hcalvetoDep_;
	TBranch *mus_pat_hcalvetoDep_branch;
	bool mus_pat_hcalvetoDep_isLoaded;
	vector<float>	mus_pat_trackIso_;
	TBranch *mus_pat_trackIso_branch;
	bool mus_pat_trackIso_isLoaded;
	vector<float>	mus_pat_vetoDep_;
	TBranch *mus_pat_vetoDep_branch;
	bool mus_pat_vetoDep_isLoaded;
	vector<float>	scs_clustersSize_;
	TBranch *scs_clustersSize_branch;
	bool scs_clustersSize_isLoaded;
	vector<float>	scs_crystalsSize_;
	TBranch *scs_crystalsSize_branch;
	bool scs_crystalsSize_isLoaded;
	vector<float>	scs_e1x3_;
	TBranch *scs_e1x3_branch;
	bool scs_e1x3_isLoaded;
	vector<float>	scs_e1x5_;
	TBranch *scs_e1x5_branch;
	bool scs_e1x5_isLoaded;
	vector<float>	scs_e2x2_;
	TBranch *scs_e2x2_branch;
	bool scs_e2x2_isLoaded;
	vector<float>	scs_e2x5Max_;
	TBranch *scs_e2x5Max_branch;
	bool scs_e2x5Max_isLoaded;
	vector<float>	scs_e3x1_;
	TBranch *scs_e3x1_branch;
	bool scs_e3x1_isLoaded;
	vector<float>	scs_e3x2_;
	TBranch *scs_e3x2_branch;
	bool scs_e3x2_isLoaded;
	vector<float>	scs_e3x3_;
	TBranch *scs_e3x3_branch;
	bool scs_e3x3_isLoaded;
	vector<float>	scs_e4x4_;
	TBranch *scs_e4x4_branch;
	bool scs_e4x4_isLoaded;
	vector<float>	scs_e5x5_;
	TBranch *scs_e5x5_branch;
	bool scs_e5x5_isLoaded;
	vector<float>	scs_energy_;
	TBranch *scs_energy_branch;
	bool scs_energy_isLoaded;
	vector<float>	scs_eta_;
	TBranch *scs_eta_branch;
	bool scs_eta_isLoaded;
	vector<float>	scs_hoe_;
	TBranch *scs_hoe_branch;
	bool scs_hoe_isLoaded;
	vector<float>	scs_phi_;
	TBranch *scs_phi_branch;
	bool scs_phi_isLoaded;
	vector<float>	scs_preshowerEnergy_;
	TBranch *scs_preshowerEnergy_branch;
	bool scs_preshowerEnergy_isLoaded;
	vector<float>	scs_rawEnergy_;
	TBranch *scs_rawEnergy_branch;
	bool scs_rawEnergy_isLoaded;
	vector<float>	scs_sigmaEtaEta_;
	TBranch *scs_sigmaEtaEta_branch;
	bool scs_sigmaEtaEta_isLoaded;
	vector<float>	scs_sigmaEtaPhi_;
	TBranch *scs_sigmaEtaPhi_branch;
	bool scs_sigmaEtaPhi_isLoaded;
	vector<float>	scs_sigmaIEtaIEta_;
	TBranch *scs_sigmaIEtaIEta_branch;
	bool scs_sigmaIEtaIEta_isLoaded;
	vector<float>	scs_sigmaIEtaIPhi_;
	TBranch *scs_sigmaIEtaIPhi_branch;
	bool scs_sigmaIEtaIPhi_isLoaded;
	vector<float>	scs_sigmaIPhiIPhi_;
	TBranch *scs_sigmaIPhiIPhi_branch;
	bool scs_sigmaIPhiIPhi_isLoaded;
	vector<float>	scs_sigmaPhiPhi_;
	TBranch *scs_sigmaPhiPhi_branch;
	bool scs_sigmaPhiPhi_isLoaded;
	vector<float>	trks_chi2_;
	TBranch *trks_chi2_branch;
	bool trks_chi2_isLoaded;
	vector<float>	trks_d0_;
	TBranch *trks_d0_branch;
	bool trks_d0_isLoaded;
	vector<float>	trks_d0Err_;
	TBranch *trks_d0Err_branch;
	bool trks_d0Err_isLoaded;
	vector<float>	trks_d0corr_;
	TBranch *trks_d0corr_branch;
	bool trks_d0corr_isLoaded;
	vector<float>	trks_etaErr_;
	TBranch *trks_etaErr_branch;
	bool trks_etaErr_isLoaded;
	vector<float>	trks_ndof_;
	TBranch *trks_ndof_branch;
	bool trks_ndof_isLoaded;
	vector<float>	trks_outerEta_;
	TBranch *trks_outerEta_branch;
	bool trks_outerEta_isLoaded;
	vector<float>	trks_outerPhi_;
	TBranch *trks_outerPhi_branch;
	bool trks_outerPhi_isLoaded;
	vector<float>	trks_phiErr_;
	TBranch *trks_phiErr_branch;
	bool trks_phiErr_isLoaded;
	vector<float>	trks_ptErr_;
	TBranch *trks_ptErr_branch;
	bool trks_ptErr_isLoaded;
	vector<float>	trks_tkIso_;
	TBranch *trks_tkIso_branch;
	bool trks_tkIso_isLoaded;
	vector<float>	trks_vertexphi_;
	TBranch *trks_vertexphi_branch;
	bool trks_vertexphi_isLoaded;
	vector<float>	trks_z0_;
	TBranch *trks_z0_branch;
	bool trks_z0_isLoaded;
	vector<float>	trks_z0Err_;
	TBranch *trks_z0Err_branch;
	bool trks_z0Err_isLoaded;
	vector<float>	trks_z0corr_;
	TBranch *trks_z0corr_branch;
	bool trks_z0corr_isLoaded;
	vector<float>	trks_elsdr_;
	TBranch *trks_elsdr_branch;
	bool trks_elsdr_isLoaded;
	vector<float>	trks_elsshFrac_;
	TBranch *trks_elsshFrac_branch;
	bool trks_elsshFrac_isLoaded;
	vector<float>	trk_musdr_;
	TBranch *trk_musdr_branch;
	bool trk_musdr_isLoaded;
	vector<float>	vtxs_chi2_;
	TBranch *vtxs_chi2_branch;
	bool vtxs_chi2_isLoaded;
	vector<float>	vtxs_ndof_;
	TBranch *vtxs_ndof_branch;
	bool vtxs_ndof_isLoaded;
	vector<float>	vtxs_xError_;
	TBranch *vtxs_xError_branch;
	bool vtxs_xError_isLoaded;
	vector<float>	vtxs_yError_;
	TBranch *vtxs_yError_branch;
	bool vtxs_yError_isLoaded;
	vector<float>	vtxs_zError_;
	TBranch *vtxs_zError_branch;
	bool vtxs_zError_isLoaded;
	vector<float>	ww_pmet_;
	TBranch *ww_pmet_branch;
	bool ww_pmet_isLoaded;
	vector<vector<float> >	hyp_jets_EMFcor_;
	TBranch *hyp_jets_EMFcor_branch;
	bool hyp_jets_EMFcor_isLoaded;
	vector<vector<float> >	hyp_jets_chFrac_;
	TBranch *hyp_jets_chFrac_branch;
	bool hyp_jets_chFrac_isLoaded;
	vector<vector<float> >	hyp_jets_cor_;
	TBranch *hyp_jets_cor_branch;
	bool hyp_jets_cor_isLoaded;
	vector<vector<float> >	hyp_jets_emFrac_;
	TBranch *hyp_jets_emFrac_branch;
	bool hyp_jets_emFrac_isLoaded;
	vector<vector<float> >	hyp_jets_mc_emEnergy_;
	TBranch *hyp_jets_mc_emEnergy_branch;
	bool hyp_jets_mc_emEnergy_isLoaded;
	vector<vector<float> >	hyp_jets_mc_hadEnergy_;
	TBranch *hyp_jets_mc_hadEnergy_branch;
	bool hyp_jets_mc_hadEnergy_isLoaded;
	vector<vector<float> >	hyp_jets_mc_invEnergy_;
	TBranch *hyp_jets_mc_invEnergy_branch;
	bool hyp_jets_mc_invEnergy_isLoaded;
	vector<vector<float> >	hyp_jets_mc_otherEnergy_;
	TBranch *hyp_jets_mc_otherEnergy_branch;
	bool hyp_jets_mc_otherEnergy_isLoaded;
	vector<vector<float> >	hyp_jets_pat_bCorrF_;
	TBranch *hyp_jets_pat_bCorrF_branch;
	bool hyp_jets_pat_bCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_pat_cCorrF_;
	TBranch *hyp_jets_pat_cCorrF_branch;
	bool hyp_jets_pat_cCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_pat_gluCorrF_;
	TBranch *hyp_jets_pat_gluCorrF_branch;
	bool hyp_jets_pat_gluCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_pat_jetCharge_;
	TBranch *hyp_jets_pat_jetCharge_branch;
	bool hyp_jets_pat_jetCharge_isLoaded;
	vector<vector<float> >	hyp_jets_pat_noCorrF_;
	TBranch *hyp_jets_pat_noCorrF_branch;
	bool hyp_jets_pat_noCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_pat_udsCorrF_;
	TBranch *hyp_jets_pat_udsCorrF_branch;
	bool hyp_jets_pat_udsCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_EMFcor_;
	TBranch *hyp_other_jets_EMFcor_branch;
	bool hyp_other_jets_EMFcor_isLoaded;
	vector<vector<float> >	hyp_other_jets_chFrac_;
	TBranch *hyp_other_jets_chFrac_branch;
	bool hyp_other_jets_chFrac_isLoaded;
	vector<vector<float> >	hyp_other_jets_cor_;
	TBranch *hyp_other_jets_cor_branch;
	bool hyp_other_jets_cor_isLoaded;
	vector<vector<float> >	hyp_other_jets_emFrac_;
	TBranch *hyp_other_jets_emFrac_branch;
	bool hyp_other_jets_emFrac_isLoaded;
	vector<vector<float> >	hyp_other_jets_mc_emEnergy_;
	TBranch *hyp_other_jets_mc_emEnergy_branch;
	bool hyp_other_jets_mc_emEnergy_isLoaded;
	vector<vector<float> >	hyp_other_jets_mc_hadEnergy_;
	TBranch *hyp_other_jets_mc_hadEnergy_branch;
	bool hyp_other_jets_mc_hadEnergy_isLoaded;
	vector<vector<float> >	hyp_other_jets_mc_invEnergy_;
	TBranch *hyp_other_jets_mc_invEnergy_branch;
	bool hyp_other_jets_mc_invEnergy_isLoaded;
	vector<vector<float> >	hyp_other_jets_mc_otherEnergy_;
	TBranch *hyp_other_jets_mc_otherEnergy_branch;
	bool hyp_other_jets_mc_otherEnergy_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_bCorrF_;
	TBranch *hyp_other_jets_pat_bCorrF_branch;
	bool hyp_other_jets_pat_bCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_cCorrF_;
	TBranch *hyp_other_jets_pat_cCorrF_branch;
	bool hyp_other_jets_pat_cCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_gluCorrF_;
	TBranch *hyp_other_jets_pat_gluCorrF_branch;
	bool hyp_other_jets_pat_gluCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_jetCharge_;
	TBranch *hyp_other_jets_pat_jetCharge_branch;
	bool hyp_other_jets_pat_jetCharge_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_noCorrF_;
	TBranch *hyp_other_jets_pat_noCorrF_branch;
	bool hyp_other_jets_pat_noCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_pat_udsCorrF_;
	TBranch *hyp_other_jets_pat_udsCorrF_branch;
	bool hyp_other_jets_pat_udsCorrF_isLoaded;
	int	evt_HLT1_;
	TBranch *evt_HLT1_branch;
	bool evt_HLT1_isLoaded;
	int	evt_HLT2_;
	TBranch *evt_HLT2_branch;
	bool evt_HLT2_isLoaded;
	int	evt_HLT3_;
	TBranch *evt_HLT3_branch;
	bool evt_HLT3_isLoaded;
	int	evt_HLT4_;
	TBranch *evt_HLT4_branch;
	bool evt_HLT4_isLoaded;
	int	evt_HLT5_;
	TBranch *evt_HLT5_branch;
	bool evt_HLT5_isLoaded;
	int	evt_HLT6_;
	TBranch *evt_HLT6_branch;
	bool evt_HLT6_isLoaded;
	int	evt_HLT7_;
	TBranch *evt_HLT7_branch;
	bool evt_HLT7_isLoaded;
	int	evt_HLT8_;
	TBranch *evt_HLT8_branch;
	bool evt_HLT8_isLoaded;
	int	evt_L1_1_;
	TBranch *evt_L1_1_branch;
	bool evt_L1_1_isLoaded;
	int	evt_L1_2_;
	TBranch *evt_L1_2_branch;
	bool evt_L1_2_isLoaded;
	int	evt_L1_3_;
	TBranch *evt_L1_3_branch;
	bool evt_L1_3_isLoaded;
	int	evt_L1_4_;
	TBranch *evt_L1_4_branch;
	bool evt_L1_4_isLoaded;
	int	evt_nl1emiso_;
	TBranch *evt_nl1emiso_branch;
	bool evt_nl1emiso_isLoaded;
	int	evt_nl1emnoiso_;
	TBranch *evt_nl1emnoiso_branch;
	bool evt_nl1emnoiso_isLoaded;
	int	evt_nl1jetsc_;
	TBranch *evt_nl1jetsc_branch;
	bool evt_nl1jetsc_isLoaded;
	int	evt_nl1jetsf_;
	TBranch *evt_nl1jetsf_branch;
	bool evt_nl1jetsf_isLoaded;
	int	evt_nl1jetst_;
	TBranch *evt_nl1jetst_branch;
	bool evt_nl1jetst_isLoaded;
	int	evt_nl1mus_;
	TBranch *evt_nl1mus_branch;
	bool evt_nl1mus_isLoaded;
	vector<int>	els_mc3_id_;
	TBranch *els_mc3_id_branch;
	bool els_mc3_id_isLoaded;
	vector<int>	els_mc3idx_;
	TBranch *els_mc3idx_branch;
	bool els_mc3idx_isLoaded;
	vector<int>	els_mc3_motherid_;
	TBranch *els_mc3_motherid_branch;
	bool els_mc3_motherid_isLoaded;
	vector<int>	els_mc_id_;
	TBranch *els_mc_id_branch;
	bool els_mc_id_isLoaded;
	vector<int>	els_mcidx_;
	TBranch *els_mcidx_branch;
	bool els_mcidx_isLoaded;
	vector<int>	els_mc_motherid_;
	TBranch *els_mc_motherid_branch;
	bool els_mc_motherid_isLoaded;
	vector<int>	jets_mc_id_;
	TBranch *jets_mc_id_branch;
	bool jets_mc_id_isLoaded;
	vector<int>	mus_mc3_id_;
	TBranch *mus_mc3_id_branch;
	bool mus_mc3_id_isLoaded;
	vector<int>	mus_mc3idx_;
	TBranch *mus_mc3idx_branch;
	bool mus_mc3idx_isLoaded;
	vector<int>	mus_mc3_motherid_;
	TBranch *mus_mc3_motherid_branch;
	bool mus_mc3_motherid_isLoaded;
	vector<int>	mus_mc_id_;
	TBranch *mus_mc_id_branch;
	bool mus_mc_id_isLoaded;
	vector<int>	mus_mcidx_;
	TBranch *mus_mcidx_branch;
	bool mus_mcidx_isLoaded;
	vector<int>	mus_mc_motherid_;
	TBranch *mus_mc_motherid_branch;
	bool mus_mc_motherid_isLoaded;
	vector<int>	trk_mc3_id_;
	TBranch *trk_mc3_id_branch;
	bool trk_mc3_id_isLoaded;
	vector<int>	trk_mc3idx_;
	TBranch *trk_mc3idx_branch;
	bool trk_mc3idx_isLoaded;
	vector<int>	trk_mc3_motherid_;
	TBranch *trk_mc3_motherid_branch;
	bool trk_mc3_motherid_isLoaded;
	vector<int>	trk_mc_id_;
	TBranch *trk_mc_id_branch;
	bool trk_mc_id_isLoaded;
	vector<int>	trk_mcidx_;
	TBranch *trk_mcidx_branch;
	bool trk_mcidx_isLoaded;
	vector<int>	trk_mc_motherid_;
	TBranch *trk_mc_motherid_branch;
	bool trk_mc_motherid_isLoaded;
	vector<int>	els_conv_tkidx_;
	TBranch *els_conv_tkidx_branch;
	bool els_conv_tkidx_isLoaded;
	vector<int>	trks_conv_tkidx_;
	TBranch *trks_conv_tkidx_branch;
	bool trks_conv_tkidx_isLoaded;
	vector<int>	els_closestJet_;
	TBranch *els_closestJet_branch;
	bool els_closestJet_isLoaded;
	vector<int>	els_closestMuon_;
	TBranch *els_closestMuon_branch;
	bool els_closestMuon_isLoaded;
	vector<int>	els_trkidx_;
	TBranch *els_trkidx_branch;
	bool els_trkidx_isLoaded;
	vector<int>	els_category_;
	TBranch *els_category_branch;
	bool els_category_isLoaded;
	vector<int>	els_categoryold_;
	TBranch *els_categoryold_branch;
	bool els_categoryold_isLoaded;
	vector<int>	els_charge_;
	TBranch *els_charge_branch;
	bool els_charge_isLoaded;
	vector<int>	els_class_;
	TBranch *els_class_branch;
	bool els_class_isLoaded;
	vector<int>	els_looseId_;
	TBranch *els_looseId_branch;
	bool els_looseId_isLoaded;
	vector<int>	els_lostHits_;
	TBranch *els_lostHits_branch;
	bool els_lostHits_isLoaded;
	vector<int>	els_nSeed_;
	TBranch *els_nSeed_branch;
	bool els_nSeed_isLoaded;
	vector<int>	els_pass3looseId_;
	TBranch *els_pass3looseId_branch;
	bool els_pass3looseId_isLoaded;
	vector<int>	els_pass3simpleId_;
	TBranch *els_pass3simpleId_branch;
	bool els_pass3simpleId_isLoaded;
	vector<int>	els_pass3tightId_;
	TBranch *els_pass3tightId_branch;
	bool els_pass3tightId_isLoaded;
	vector<int>	els_robustId_;
	TBranch *els_robustId_branch;
	bool els_robustId_isLoaded;
	vector<int>	els_simpleIdPlus_;
	TBranch *els_simpleIdPlus_branch;
	bool els_simpleIdPlus_isLoaded;
	vector<int>	els_tightId_;
	TBranch *els_tightId_branch;
	bool els_tightId_isLoaded;
	vector<int>	els_tightId22XMaxMatteo_;
	TBranch *els_tightId22XMaxMatteo_branch;
	bool els_tightId22XMaxMatteo_isLoaded;
	vector<int>	els_tightId22XMinMatteo_;
	TBranch *els_tightId22XMinMatteo_branch;
	bool els_tightId22XMinMatteo_isLoaded;
	vector<int>	els_validHits_;
	TBranch *els_validHits_branch;
	bool els_validHits_isLoaded;
	vector<int>	genps_id_;
	TBranch *genps_id_branch;
	bool genps_id_isLoaded;
	vector<int>	genps_id_mother_;
	TBranch *genps_id_mother_branch;
	bool genps_id_mother_isLoaded;
	vector<int>	genps_status_;
	TBranch *genps_status_branch;
	bool genps_status_isLoaded;
	vector<int>	hyp_ll_charge_;
	TBranch *hyp_ll_charge_branch;
	bool hyp_ll_charge_isLoaded;
	vector<int>	hyp_ll_id_;
	TBranch *hyp_ll_id_branch;
	bool hyp_ll_id_isLoaded;
	vector<int>	hyp_ll_index_;
	TBranch *hyp_ll_index_branch;
	bool hyp_ll_index_isLoaded;
	vector<int>	hyp_ll_lostHits_;
	TBranch *hyp_ll_lostHits_branch;
	bool hyp_ll_lostHits_isLoaded;
	vector<int>	hyp_ll_mc_id_;
	TBranch *hyp_ll_mc_id_branch;
	bool hyp_ll_mc_id_isLoaded;
	vector<int>	hyp_ll_mc_motherid_;
	TBranch *hyp_ll_mc_motherid_branch;
	bool hyp_ll_mc_motherid_isLoaded;
	vector<int>	hyp_ll_validHits_;
	TBranch *hyp_ll_validHits_branch;
	bool hyp_ll_validHits_isLoaded;
	vector<int>	hyp_lt_charge_;
	TBranch *hyp_lt_charge_branch;
	bool hyp_lt_charge_isLoaded;
	vector<int>	hyp_lt_id_;
	TBranch *hyp_lt_id_branch;
	bool hyp_lt_id_isLoaded;
	vector<int>	hyp_lt_index_;
	TBranch *hyp_lt_index_branch;
	bool hyp_lt_index_isLoaded;
	vector<int>	hyp_lt_lostHits_;
	TBranch *hyp_lt_lostHits_branch;
	bool hyp_lt_lostHits_isLoaded;
	vector<int>	hyp_lt_mc_id_;
	TBranch *hyp_lt_mc_id_branch;
	bool hyp_lt_mc_id_isLoaded;
	vector<int>	hyp_lt_mc_motherid_;
	TBranch *hyp_lt_mc_motherid_branch;
	bool hyp_lt_mc_motherid_isLoaded;
	vector<int>	hyp_lt_validHits_;
	TBranch *hyp_lt_validHits_branch;
	bool hyp_lt_validHits_isLoaded;
	vector<int>	hyp_njets_;
	TBranch *hyp_njets_branch;
	bool hyp_njets_isLoaded;
	vector<int>	hyp_nojets_;
	TBranch *hyp_nojets_branch;
	bool hyp_nojets_isLoaded;
	vector<int>	hyp_type_;
	TBranch *hyp_type_branch;
	bool hyp_type_isLoaded;
	vector<int>	hyp_quadlep_first_type_;
	TBranch *hyp_quadlep_first_type_branch;
	bool hyp_quadlep_first_type_isLoaded;
	vector<int>	hyp_quadlep_fourth_type_;
	TBranch *hyp_quadlep_fourth_type_branch;
	bool hyp_quadlep_fourth_type_isLoaded;
	vector<int>	hyp_quadlep_second_type_;
	TBranch *hyp_quadlep_second_type_branch;
	bool hyp_quadlep_second_type_isLoaded;
	vector<int>	hyp_quadlep_third_type_;
	TBranch *hyp_quadlep_third_type_branch;
	bool hyp_quadlep_third_type_isLoaded;
	vector<int>	hyp_trilep_first_type_;
	TBranch *hyp_trilep_first_type_branch;
	bool hyp_trilep_first_type_isLoaded;
	vector<int>	hyp_trilep_second_type_;
	TBranch *hyp_trilep_second_type_branch;
	bool hyp_trilep_second_type_isLoaded;
	vector<int>	hyp_trilep_third_type_;
	TBranch *hyp_trilep_third_type_branch;
	bool hyp_trilep_third_type_isLoaded;
	vector<int>	jets_closestElectron_;
	TBranch *jets_closestElectron_branch;
	bool jets_closestElectron_isLoaded;
	vector<int>	jets_closestMuon_;
	TBranch *jets_closestMuon_branch;
	bool jets_closestMuon_isLoaded;
	vector<int>	l1emiso_ieta_;
	TBranch *l1emiso_ieta_branch;
	bool l1emiso_ieta_isLoaded;
	vector<int>	l1emiso_iphi_;
	TBranch *l1emiso_iphi_branch;
	bool l1emiso_iphi_isLoaded;
	vector<int>	l1emiso_rawId_;
	TBranch *l1emiso_rawId_branch;
	bool l1emiso_rawId_isLoaded;
	vector<int>	l1emiso_type_;
	TBranch *l1emiso_type_branch;
	bool l1emiso_type_isLoaded;
	vector<int>	l1emnoiso_ieta_;
	TBranch *l1emnoiso_ieta_branch;
	bool l1emnoiso_ieta_isLoaded;
	vector<int>	l1emnoiso_iphi_;
	TBranch *l1emnoiso_iphi_branch;
	bool l1emnoiso_iphi_isLoaded;
	vector<int>	l1emnoiso_rawId_;
	TBranch *l1emnoiso_rawId_branch;
	bool l1emnoiso_rawId_isLoaded;
	vector<int>	l1emnoiso_type_;
	TBranch *l1emnoiso_type_branch;
	bool l1emnoiso_type_isLoaded;
	vector<int>	l1jetsc_ieta_;
	TBranch *l1jetsc_ieta_branch;
	bool l1jetsc_ieta_isLoaded;
	vector<int>	l1jetsc_iphi_;
	TBranch *l1jetsc_iphi_branch;
	bool l1jetsc_iphi_isLoaded;
	vector<int>	l1jetsc_rawId_;
	TBranch *l1jetsc_rawId_branch;
	bool l1jetsc_rawId_isLoaded;
	vector<int>	l1jetsc_type_;
	TBranch *l1jetsc_type_branch;
	bool l1jetsc_type_isLoaded;
	vector<int>	l1jetsf_ieta_;
	TBranch *l1jetsf_ieta_branch;
	bool l1jetsf_ieta_isLoaded;
	vector<int>	l1jetsf_iphi_;
	TBranch *l1jetsf_iphi_branch;
	bool l1jetsf_iphi_isLoaded;
	vector<int>	l1jetsf_rawId_;
	TBranch *l1jetsf_rawId_branch;
	bool l1jetsf_rawId_isLoaded;
	vector<int>	l1jetsf_type_;
	TBranch *l1jetsf_type_branch;
	bool l1jetsf_type_isLoaded;
	vector<int>	l1jetst_ieta_;
	TBranch *l1jetst_ieta_branch;
	bool l1jetst_ieta_isLoaded;
	vector<int>	l1jetst_iphi_;
	TBranch *l1jetst_iphi_branch;
	bool l1jetst_iphi_isLoaded;
	vector<int>	l1jetst_rawId_;
	TBranch *l1jetst_rawId_branch;
	bool l1jetst_rawId_isLoaded;
	vector<int>	l1jetst_type_;
	TBranch *l1jetst_type_branch;
	bool l1jetst_type_isLoaded;
	vector<int>	l1mus_flags_;
	TBranch *l1mus_flags_branch;
	bool l1mus_flags_isLoaded;
	vector<int>	l1mus_q_;
	TBranch *l1mus_q_branch;
	bool l1mus_q_isLoaded;
	vector<int>	l1mus_qual_;
	TBranch *l1mus_qual_branch;
	bool l1mus_qual_isLoaded;
	vector<int>	l1mus_qualFlags_;
	TBranch *l1mus_qualFlags_branch;
	bool l1mus_qualFlags_isLoaded;
	vector<int>	mus_closestEle_;
	TBranch *mus_closestEle_branch;
	bool mus_closestEle_isLoaded;
	vector<int>	mus_closestJet_;
	TBranch *mus_closestJet_branch;
	bool mus_closestJet_isLoaded;
	vector<int>	mus_trkidx_;
	TBranch *mus_trkidx_branch;
	bool mus_trkidx_isLoaded;
	vector<int>	mus_charge_;
	TBranch *mus_charge_branch;
	bool mus_charge_isLoaded;
	vector<int>	mus_gfit_validHits_;
	TBranch *mus_gfit_validHits_branch;
	bool mus_gfit_validHits_isLoaded;
	vector<int>	mus_goodmask_;
	TBranch *mus_goodmask_branch;
	bool mus_goodmask_isLoaded;
	vector<int>	mus_iso03_ntrk_;
	TBranch *mus_iso03_ntrk_branch;
	bool mus_iso03_ntrk_isLoaded;
	vector<int>	mus_iso05_ntrk_;
	TBranch *mus_iso05_ntrk_branch;
	bool mus_iso05_ntrk_isLoaded;
	vector<int>	mus_lostHits_;
	TBranch *mus_lostHits_branch;
	bool mus_lostHits_isLoaded;
	vector<int>	mus_nmatches_;
	TBranch *mus_nmatches_branch;
	bool mus_nmatches_isLoaded;
	vector<int>	mus_pid_TM2DCompatibilityLoose_;
	TBranch *mus_pid_TM2DCompatibilityLoose_branch;
	bool mus_pid_TM2DCompatibilityLoose_isLoaded;
	vector<int>	mus_pid_TM2DCompatibilityTight_;
	TBranch *mus_pid_TM2DCompatibilityTight_branch;
	bool mus_pid_TM2DCompatibilityTight_isLoaded;
	vector<int>	mus_pid_TMLastStationLoose_;
	TBranch *mus_pid_TMLastStationLoose_branch;
	bool mus_pid_TMLastStationLoose_isLoaded;
	vector<int>	mus_pid_TMLastStationTight_;
	TBranch *mus_pid_TMLastStationTight_branch;
	bool mus_pid_TMLastStationTight_isLoaded;
	vector<int>	mus_trk_charge_;
	TBranch *mus_trk_charge_branch;
	bool mus_trk_charge_isLoaded;
	vector<int>	mus_trkrefkey_;
	TBranch *mus_trkrefkey_branch;
	bool mus_trkrefkey_isLoaded;
	vector<int>	mus_type_;
	TBranch *mus_type_branch;
	bool mus_type_isLoaded;
	vector<int>	mus_validHits_;
	TBranch *mus_validHits_branch;
	bool mus_validHits_isLoaded;
	vector<int>	els_pat_genID_;
	TBranch *els_pat_genID_branch;
	bool els_pat_genID_isLoaded;
	vector<int>	els_pat_genMotherID_;
	TBranch *els_pat_genMotherID_branch;
	bool els_pat_genMotherID_isLoaded;
	vector<int>	jets_pat_genPartonMother_id_;
	TBranch *jets_pat_genPartonMother_id_branch;
	bool jets_pat_genPartonMother_id_isLoaded;
	vector<int>	jets_pat_genParton_id_;
	TBranch *jets_pat_genParton_id_branch;
	bool jets_pat_genParton_id_isLoaded;
	vector<int>	jets_pat_partonFlavour_;
	TBranch *jets_pat_partonFlavour_branch;
	bool jets_pat_partonFlavour_isLoaded;
	vector<int>	mus_pat_genID_;
	TBranch *mus_pat_genID_branch;
	bool mus_pat_genID_isLoaded;
	vector<int>	mus_pat_genMotherID_;
	TBranch *mus_pat_genMotherID_branch;
	bool mus_pat_genMotherID_isLoaded;
	vector<int>	trks_charge_;
	TBranch *trks_charge_branch;
	bool trks_charge_isLoaded;
	vector<int>	trks_lostHits_;
	TBranch *trks_lostHits_branch;
	bool trks_lostHits_isLoaded;
	vector<int>	trks_validHits_;
	TBranch *trks_validHits_branch;
	bool trks_validHits_isLoaded;
	vector<int>	trks_elsidx_;
	TBranch *trks_elsidx_branch;
	bool trks_elsidx_isLoaded;
	vector<int>	trk_musidx_;
	TBranch *trk_musidx_branch;
	bool trk_musidx_isLoaded;
	vector<int>	hlt2ele10LWR_id_;
	TBranch *hlt2ele10LWR_id_branch;
	bool hlt2ele10LWR_id_isLoaded;
	vector<int>	hlt2ele10LWR_tid_;
	TBranch *hlt2ele10LWR_tid_branch;
	bool hlt2ele10LWR_tid_isLoaded;
	vector<int>	hlt2mu3_id_;
	TBranch *hlt2mu3_id_branch;
	bool hlt2mu3_id_isLoaded;
	vector<int>	hlt2mu3_tid_;
	TBranch *hlt2mu3_tid_branch;
	bool hlt2mu3_tid_isLoaded;
	vector<int>	hltLisoele18LWR_id_;
	TBranch *hltLisoele18LWR_id_branch;
	bool hltLisoele18LWR_id_isLoaded;
	vector<int>	hltLisoele18LWR_tid_;
	TBranch *hltLisoele18LWR_tid_branch;
	bool hltLisoele18LWR_tid_isLoaded;
	vector<int>	hltisoele18R_id_;
	TBranch *hltisoele18R_id_branch;
	bool hltisoele18R_id_isLoaded;
	vector<int>	hltisoele18R_tid_;
	TBranch *hltisoele18R_tid_branch;
	bool hltisoele18R_tid_isLoaded;
	vector<int>	hltjet30_id_;
	TBranch *hltjet30_id_branch;
	bool hltjet30_id_isLoaded;
	vector<int>	hltjet30_tid_;
	TBranch *hltjet30_tid_branch;
	bool hltjet30_tid_isLoaded;
	vector<int>	hltl1jet15_id_;
	TBranch *hltl1jet15_id_branch;
	bool hltl1jet15_id_isLoaded;
	vector<int>	hltl1jet15_tid_;
	TBranch *hltl1jet15_tid_branch;
	bool hltl1jet15_tid_isLoaded;
	vector<int>	hltl1met20_id_;
	TBranch *hltl1met20_id_branch;
	bool hltl1met20_id_isLoaded;
	vector<int>	hltl1met20_tid_;
	TBranch *hltl1met20_tid_branch;
	bool hltl1met20_tid_isLoaded;
	vector<int>	hltl1mu_id_;
	TBranch *hltl1mu_id_branch;
	bool hltl1mu_id_isLoaded;
	vector<int>	hltl1mu_tid_;
	TBranch *hltl1mu_tid_branch;
	bool hltl1mu_tid_isLoaded;
	vector<int>	hltl2mu9_id_;
	TBranch *hltl2mu9_id_branch;
	bool hltl2mu9_id_isLoaded;
	vector<int>	hltl2mu9_tid_;
	TBranch *hltl2mu9_tid_branch;
	bool hltl2mu9_tid_isLoaded;
	vector<int>	hltmet25_id_;
	TBranch *hltmet25_id_branch;
	bool hltmet25_id_isLoaded;
	vector<int>	hltmet25_tid_;
	TBranch *hltmet25_tid_branch;
	bool hltmet25_tid_isLoaded;
	vector<int>	hltmu11_id_;
	TBranch *hltmu11_id_branch;
	bool hltmu11_id_isLoaded;
	vector<int>	hltmu11_tid_;
	TBranch *hltmu11_tid_branch;
	bool hltmu11_tid_isLoaded;
	vector<int>	hltmu9_id_;
	TBranch *hltmu9_id_branch;
	bool hltmu9_id_isLoaded;
	vector<int>	hltmu9_tid_;
	TBranch *hltmu9_tid_branch;
	bool hltmu9_tid_isLoaded;
	vector<int>	vtxs_hasRefittedTracks_;
	TBranch *vtxs_hasRefittedTracks_branch;
	bool vtxs_hasRefittedTracks_isLoaded;
	vector<int>	vtxs_isFake_;
	TBranch *vtxs_isFake_branch;
	bool vtxs_isFake_isLoaded;
	vector<int>	vtxs_isValid_;
	TBranch *vtxs_isValid_branch;
	bool vtxs_isValid_isLoaded;
	vector<int>	vtxs_tracksSize_;
	TBranch *vtxs_tracksSize_branch;
	bool vtxs_tracksSize_isLoaded;
	vector<int>	ww_isdyee_;
	TBranch *ww_isdyee_branch;
	bool ww_isdyee_isLoaded;
	vector<int>	ww_isdymm_;
	TBranch *ww_isdymm_branch;
	bool ww_isdymm_isLoaded;
	vector<int>	ww_isdytt_;
	TBranch *ww_isdytt_branch;
	bool ww_isdytt_isLoaded;
	vector<int>	ww_llgoodel_;
	TBranch *ww_llgoodel_branch;
	bool ww_llgoodel_isLoaded;
	vector<int>	ww_llgoodeliso_;
	TBranch *ww_llgoodeliso_branch;
	bool ww_llgoodeliso_isLoaded;
	vector<int>	ww_llgoodmu_;
	TBranch *ww_llgoodmu_branch;
	bool ww_llgoodmu_isLoaded;
	vector<int>	ww_llgoodmuiso_;
	TBranch *ww_llgoodmuiso_branch;
	bool ww_llgoodmuiso_isLoaded;
	vector<int>	ww_ltgoodel_;
	TBranch *ww_ltgoodel_branch;
	bool ww_ltgoodel_isLoaded;
	vector<int>	ww_ltgoodeliso_;
	TBranch *ww_ltgoodeliso_branch;
	bool ww_ltgoodeliso_isLoaded;
	vector<int>	ww_ltgoodmu_;
	TBranch *ww_ltgoodmu_branch;
	bool ww_ltgoodmu_isLoaded;
	vector<int>	ww_ltgoodmuiso_;
	TBranch *ww_ltgoodmuiso_branch;
	bool ww_ltgoodmuiso_isLoaded;
	vector<int>	ww_oppsign_;
	TBranch *ww_oppsign_branch;
	bool ww_oppsign_isLoaded;
	vector<int>	ww_pass2met_;
	TBranch *ww_pass2met_branch;
	bool ww_pass2met_isLoaded;
	vector<int>	ww_pass4met_;
	TBranch *ww_pass4met_branch;
	bool ww_pass4met_isLoaded;
	vector<int>	ww_passFebsel_no_jetveto_;
	TBranch *ww_passFebsel_no_jetveto_branch;
	bool ww_passFebsel_no_jetveto_isLoaded;
	vector<int>	ww_passFebsel_with_jetveto_;
	TBranch *ww_passFebsel_with_jetveto_branch;
	bool ww_passFebsel_with_jetveto_isLoaded;
	vector<int>	ww_passaddzveto_;
	TBranch *ww_passaddzveto_branch;
	bool ww_passaddzveto_isLoaded;
	vector<int>	ww_passzveto_;
	TBranch *ww_passzveto_branch;
	bool ww_passzveto_isLoaded;
	vector<vector<int> >	hyp_jets_mc_id_;
	TBranch *hyp_jets_mc_id_branch;
	bool hyp_jets_mc_id_isLoaded;
	vector<vector<int> >	hyp_jets_pat_genPartonMother_id_;
	TBranch *hyp_jets_pat_genPartonMother_id_branch;
	bool hyp_jets_pat_genPartonMother_id_isLoaded;
	vector<vector<int> >	hyp_jets_pat_genParton_id_;
	TBranch *hyp_jets_pat_genParton_id_branch;
	bool hyp_jets_pat_genParton_id_isLoaded;
	vector<vector<int> >	hyp_jets_pat_partonFlavour_;
	TBranch *hyp_jets_pat_partonFlavour_branch;
	bool hyp_jets_pat_partonFlavour_isLoaded;
	vector<vector<int> >	hyp_other_jets_mc_id_;
	TBranch *hyp_other_jets_mc_id_branch;
	bool hyp_other_jets_mc_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_pat_genPartonMother_id_;
	TBranch *hyp_other_jets_pat_genPartonMother_id_branch;
	bool hyp_other_jets_pat_genPartonMother_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_pat_genParton_id_;
	TBranch *hyp_other_jets_pat_genParton_id_branch;
	bool hyp_other_jets_pat_genParton_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_pat_partonFlavour_;
	TBranch *hyp_other_jets_pat_partonFlavour_branch;
	bool hyp_other_jets_pat_partonFlavour_isLoaded;
	vector<vector<int> >	hyp_quadlep_jets_index_;
	TBranch *hyp_quadlep_jets_index_branch;
	bool hyp_quadlep_jets_index_isLoaded;
	vector<vector<int> >	hyp_trilep_jets_index_;
	TBranch *hyp_trilep_jets_index_branch;
	bool hyp_trilep_jets_index_isLoaded;
	unsigned int	evt_nalltrkjets_;
	TBranch *evt_nalltrkjets_branch;
	bool evt_nalltrkjets_isLoaded;
	unsigned int	evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int	evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int	evt_lumiBlock_;
	TBranch *evt_lumiBlock_branch;
	bool evt_lumiBlock_isLoaded;
	unsigned int	evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int	evt_ngenjets_;
	TBranch *evt_ngenjets_branch;
	bool evt_ngenjets_isLoaded;
	unsigned int	evt_njets_;
	TBranch *evt_njets_branch;
	bool evt_njets_isLoaded;
	unsigned int	evt_njpts_;
	TBranch *evt_njpts_branch;
	bool evt_njpts_isLoaded;
	unsigned int	evt_nscs_;
	TBranch *evt_nscs_branch;
	bool evt_nscs_isLoaded;
	unsigned int	evt_ntrkjets_;
	TBranch *evt_ntrkjets_branch;
	bool evt_ntrkjets_isLoaded;
	unsigned int	evt_nvtxs_;
	TBranch *evt_nvtxs_branch;
	bool evt_nvtxs_isLoaded;
	vector<unsigned int>	hyp_quadlep_bucket_;
	TBranch *hyp_quadlep_bucket_branch;
	bool hyp_quadlep_bucket_isLoaded;
	vector<unsigned int>	hyp_quadlep_first_index_;
	TBranch *hyp_quadlep_first_index_branch;
	bool hyp_quadlep_first_index_isLoaded;
	vector<unsigned int>	hyp_quadlep_fourth_index_;
	TBranch *hyp_quadlep_fourth_index_branch;
	bool hyp_quadlep_fourth_index_isLoaded;
	vector<unsigned int>	hyp_quadlep_second_index_;
	TBranch *hyp_quadlep_second_index_branch;
	bool hyp_quadlep_second_index_isLoaded;
	vector<unsigned int>	hyp_quadlep_third_index_;
	TBranch *hyp_quadlep_third_index_branch;
	bool hyp_quadlep_third_index_isLoaded;
	vector<unsigned int>	hyp_trilep_bucket_;
	TBranch *hyp_trilep_bucket_branch;
	bool hyp_trilep_bucket_isLoaded;
	vector<unsigned int>	hyp_trilep_first_index_;
	TBranch *hyp_trilep_first_index_branch;
	bool hyp_trilep_first_index_isLoaded;
	vector<unsigned int>	hyp_trilep_second_index_;
	TBranch *hyp_trilep_second_index_branch;
	bool hyp_trilep_second_index_isLoaded;
	vector<unsigned int>	hyp_trilep_third_index_;
	TBranch *hyp_trilep_third_index_branch;
	bool hyp_trilep_third_index_isLoaded;
	vector<unsigned int>	els_pat_flag_;
	TBranch *els_pat_flag_branch;
	bool els_pat_flag_isLoaded;
	vector<unsigned int>	jets_pat_flag_;
	TBranch *jets_pat_flag_branch;
	bool jets_pat_flag_isLoaded;
	vector<unsigned int>	mus_pat_flag_;
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
int ScanChain(class TChain* chain, int nEvents=-1);
void Init(TTree *tree) {
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
	alltrkjets_p4_branch = 0;
	if (tree->GetAlias("alltrkjets_p4") != 0) {
		alltrkjets_p4_branch = tree->GetBranch(tree->GetAlias("alltrkjets_p4"));
		alltrkjets_p4_branch->SetAddress(&alltrkjets_p4_);
	}
	if(alltrkjets_p4_branch == 0 ) {
	cout << "Branch alltrkjets_p4 does not exist." << endl;
	}
	els_mc3_p4_branch = 0;
	if (tree->GetAlias("els_mc3_p4") != 0) {
		els_mc3_p4_branch = tree->GetBranch(tree->GetAlias("els_mc3_p4"));
		els_mc3_p4_branch->SetAddress(&els_mc3_p4_);
	}
	if(els_mc3_p4_branch == 0 ) {
	cout << "Branch els_mc3_p4 does not exist." << endl;
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
	mus_mc3_p4_branch = 0;
	if (tree->GetAlias("mus_mc3_p4") != 0) {
		mus_mc3_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc3_p4"));
		mus_mc3_p4_branch->SetAddress(&mus_mc3_p4_);
	}
	if(mus_mc3_p4_branch == 0 ) {
	cout << "Branch mus_mc3_p4 does not exist." << endl;
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
	scs_p4_branch = 0;
	if (tree->GetAlias("scs_p4") != 0) {
		scs_p4_branch = tree->GetBranch(tree->GetAlias("scs_p4"));
		scs_p4_branch->SetAddress(&scs_p4_);
	}
	if(scs_p4_branch == 0 ) {
	cout << "Branch scs_p4 does not exist." << endl;
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
	vtxs_position_branch = 0;
	if (tree->GetAlias("vtxs_position") != 0) {
		vtxs_position_branch = tree->GetBranch(tree->GetAlias("vtxs_position"));
		vtxs_position_branch->SetAddress(&vtxs_position_);
	}
	if(vtxs_position_branch == 0 ) {
	cout << "Branch vtxs_position does not exist." << endl;
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
	hyp_jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_gp_p4") != 0) {
		hyp_jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_gp_p4"));
		hyp_jets_mc_gp_p4_branch->SetAddress(&hyp_jets_mc_gp_p4_);
	}
	if(hyp_jets_mc_gp_p4_branch == 0 ) {
	cout << "Branch hyp_jets_mc_gp_p4 does not exist." << endl;
	}
	hyp_jets_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_p4") != 0) {
		hyp_jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_p4"));
		hyp_jets_mc_p4_branch->SetAddress(&hyp_jets_mc_p4_);
	}
	if(hyp_jets_mc_p4_branch == 0 ) {
	cout << "Branch hyp_jets_mc_p4 does not exist." << endl;
	}
	hyp_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_p4") != 0) {
		hyp_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_p4"));
		hyp_jets_p4_branch->SetAddress(&hyp_jets_p4_);
	}
	if(hyp_jets_p4_branch == 0 ) {
	cout << "Branch hyp_jets_p4 does not exist." << endl;
	}
	hyp_jets_pat_genJet_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_genJet_p4") != 0) {
		hyp_jets_pat_genJet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_genJet_p4"));
		hyp_jets_pat_genJet_p4_branch->SetAddress(&hyp_jets_pat_genJet_p4_);
	}
	if(hyp_jets_pat_genJet_p4_branch == 0 ) {
	cout << "Branch hyp_jets_pat_genJet_p4 does not exist." << endl;
	}
	hyp_jets_pat_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_genPartonMother_p4") != 0) {
		hyp_jets_pat_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_genPartonMother_p4"));
		hyp_jets_pat_genPartonMother_p4_branch->SetAddress(&hyp_jets_pat_genPartonMother_p4_);
	}
	if(hyp_jets_pat_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_jets_pat_genPartonMother_p4 does not exist." << endl;
	}
	hyp_jets_pat_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_genParton_p4") != 0) {
		hyp_jets_pat_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_genParton_p4"));
		hyp_jets_pat_genParton_p4_branch->SetAddress(&hyp_jets_pat_genParton_p4_);
	}
	if(hyp_jets_pat_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_jets_pat_genParton_p4 does not exist." << endl;
	}
	hyp_jets_pat_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_jet_p4") != 0) {
		hyp_jets_pat_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_jet_p4"));
		hyp_jets_pat_jet_p4_branch->SetAddress(&hyp_jets_pat_jet_p4_);
	}
	if(hyp_jets_pat_jet_p4_branch == 0 ) {
	cout << "Branch hyp_jets_pat_jet_p4 does not exist." << endl;
	}
	hyp_other_jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_gp_p4") != 0) {
		hyp_other_jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_gp_p4"));
		hyp_other_jets_mc_gp_p4_branch->SetAddress(&hyp_other_jets_mc_gp_p4_);
	}
	if(hyp_other_jets_mc_gp_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_gp_p4 does not exist." << endl;
	}
	hyp_other_jets_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_p4") != 0) {
		hyp_other_jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_p4"));
		hyp_other_jets_mc_p4_branch->SetAddress(&hyp_other_jets_mc_p4_);
	}
	if(hyp_other_jets_mc_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_p4 does not exist." << endl;
	}
	hyp_other_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_p4") != 0) {
		hyp_other_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_p4"));
		hyp_other_jets_p4_branch->SetAddress(&hyp_other_jets_p4_);
	}
	if(hyp_other_jets_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_p4 does not exist." << endl;
	}
	hyp_other_jets_pat_genJet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_genJet_p4") != 0) {
		hyp_other_jets_pat_genJet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_genJet_p4"));
		hyp_other_jets_pat_genJet_p4_branch->SetAddress(&hyp_other_jets_pat_genJet_p4_);
	}
	if(hyp_other_jets_pat_genJet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_genJet_p4 does not exist." << endl;
	}
	hyp_other_jets_pat_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_genPartonMother_p4") != 0) {
		hyp_other_jets_pat_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_genPartonMother_p4"));
		hyp_other_jets_pat_genPartonMother_p4_branch->SetAddress(&hyp_other_jets_pat_genPartonMother_p4_);
	}
	if(hyp_other_jets_pat_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_genPartonMother_p4 does not exist." << endl;
	}
	hyp_other_jets_pat_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_genParton_p4") != 0) {
		hyp_other_jets_pat_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_genParton_p4"));
		hyp_other_jets_pat_genParton_p4_branch->SetAddress(&hyp_other_jets_pat_genParton_p4_);
	}
	if(hyp_other_jets_pat_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_genParton_p4 does not exist." << endl;
	}
	hyp_other_jets_pat_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_jet_p4") != 0) {
		hyp_other_jets_pat_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_jet_p4"));
		hyp_other_jets_pat_jet_p4_branch->SetAddress(&hyp_other_jets_pat_jet_p4_);
	}
	if(hyp_other_jets_pat_jet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_jet_p4 does not exist." << endl;
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
	evt_metSig_branch = 0;
	if (tree->GetAlias("evt_metSig") != 0) {
		evt_metSig_branch = tree->GetBranch(tree->GetAlias("evt_metSig"));
		evt_metSig_branch->SetAddress(&evt_metSig_);
	}
	if(evt_metSig_branch == 0 ) {
	cout << "Branch evt_metSig does not exist." << endl;
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
	evt_metOptSig_branch = 0;
	if (tree->GetAlias("evt_metOptSig") != 0) {
		evt_metOptSig_branch = tree->GetBranch(tree->GetAlias("evt_metOptSig"));
		evt_metOptSig_branch->SetAddress(&evt_metOptSig_);
	}
	if(evt_metOptSig_branch == 0 ) {
	cout << "Branch evt_metOptSig does not exist." << endl;
	}
	evt_metOptPhi_branch = 0;
	if (tree->GetAlias("evt_metOptPhi") != 0) {
		evt_metOptPhi_branch = tree->GetBranch(tree->GetAlias("evt_metOptPhi"));
		evt_metOptPhi_branch->SetAddress(&evt_metOptPhi_);
	}
	if(evt_metOptPhi_branch == 0 ) {
	cout << "Branch evt_metOptPhi does not exist." << endl;
	}
	evt_metPhi_branch = 0;
	if (tree->GetAlias("evt_metPhi") != 0) {
		evt_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_metPhi"));
		evt_metPhi_branch->SetAddress(&evt_metPhi_);
	}
	if(evt_metPhi_branch == 0 ) {
	cout << "Branch evt_metPhi does not exist." << endl;
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
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut_);
	}
	if(els_eSeedOverPOut_branch == 0 ) {
	cout << "Branch els_eSeedOverPOut does not exist." << endl;
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		els_etaErr_branch->SetAddress(&els_etaErr_);
	}
	if(els_etaErr_branch == 0 ) {
	cout << "Branch els_etaErr does not exist." << endl;
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
	hyp_met_branch = 0;
	if (tree->GetAlias("hyp_met") != 0) {
		hyp_met_branch = tree->GetBranch(tree->GetAlias("hyp_met"));
		hyp_met_branch->SetAddress(&hyp_met_);
	}
	if(hyp_met_branch == 0 ) {
	cout << "Branch hyp_met does not exist." << endl;
	}
	hyp_metAll_branch = 0;
	if (tree->GetAlias("hyp_metAll") != 0) {
		hyp_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_metAll"));
		hyp_metAll_branch->SetAddress(&hyp_metAll_);
	}
	if(hyp_metAll_branch == 0 ) {
	cout << "Branch hyp_metAll does not exist." << endl;
	}
	hyp_metAllCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metAllCaloExp") != 0) {
		hyp_metAllCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metAllCaloExp"));
		hyp_metAllCaloExp_branch->SetAddress(&hyp_metAllCaloExp_);
	}
	if(hyp_metAllCaloExp_branch == 0 ) {
	cout << "Branch hyp_metAllCaloExp does not exist." << endl;
	}
	hyp_metCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metCaloExp") != 0) {
		hyp_metCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metCaloExp"));
		hyp_metCaloExp_branch->SetAddress(&hyp_metCaloExp_);
	}
	if(hyp_metCaloExp_branch == 0 ) {
	cout << "Branch hyp_metCaloExp does not exist." << endl;
	}
	hyp_metCone_branch = 0;
	if (tree->GetAlias("hyp_metCone") != 0) {
		hyp_metCone_branch = tree->GetBranch(tree->GetAlias("hyp_metCone"));
		hyp_metCone_branch->SetAddress(&hyp_metCone_);
	}
	if(hyp_metCone_branch == 0 ) {
	cout << "Branch hyp_metCone does not exist." << endl;
	}
	hyp_metDPhiJet10_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet10") != 0) {
		hyp_metDPhiJet10_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet10"));
		hyp_metDPhiJet10_branch->SetAddress(&hyp_metDPhiJet10_);
	}
	if(hyp_metDPhiJet10_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet10 does not exist." << endl;
	}
	hyp_metDPhiJet15_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet15") != 0) {
		hyp_metDPhiJet15_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet15"));
		hyp_metDPhiJet15_branch->SetAddress(&hyp_metDPhiJet15_);
	}
	if(hyp_metDPhiJet15_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet15 does not exist." << endl;
	}
	hyp_metDPhiJet20_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet20") != 0) {
		hyp_metDPhiJet20_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet20"));
		hyp_metDPhiJet20_branch->SetAddress(&hyp_metDPhiJet20_);
	}
	if(hyp_metDPhiJet20_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet20 does not exist." << endl;
	}
	hyp_metDPhiTrk10_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk10") != 0) {
		hyp_metDPhiTrk10_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk10"));
		hyp_metDPhiTrk10_branch->SetAddress(&hyp_metDPhiTrk10_);
	}
	if(hyp_metDPhiTrk10_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk10 does not exist." << endl;
	}
	hyp_metDPhiTrk25_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk25") != 0) {
		hyp_metDPhiTrk25_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk25"));
		hyp_metDPhiTrk25_branch->SetAddress(&hyp_metDPhiTrk25_);
	}
	if(hyp_metDPhiTrk25_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk25 does not exist." << endl;
	}
	hyp_metDPhiTrk50_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk50") != 0) {
		hyp_metDPhiTrk50_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk50"));
		hyp_metDPhiTrk50_branch->SetAddress(&hyp_metDPhiTrk50_);
	}
	if(hyp_metDPhiTrk50_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk50 does not exist." << endl;
	}
	hyp_metJes10_branch = 0;
	if (tree->GetAlias("hyp_metJes10") != 0) {
		hyp_metJes10_branch = tree->GetBranch(tree->GetAlias("hyp_metJes10"));
		hyp_metJes10_branch->SetAddress(&hyp_metJes10_);
	}
	if(hyp_metJes10_branch == 0 ) {
	cout << "Branch hyp_metJes10 does not exist." << endl;
	}
	hyp_metJes15_branch = 0;
	if (tree->GetAlias("hyp_metJes15") != 0) {
		hyp_metJes15_branch = tree->GetBranch(tree->GetAlias("hyp_metJes15"));
		hyp_metJes15_branch->SetAddress(&hyp_metJes15_);
	}
	if(hyp_metJes15_branch == 0 ) {
	cout << "Branch hyp_metJes15 does not exist." << endl;
	}
	hyp_metJes30_branch = 0;
	if (tree->GetAlias("hyp_metJes30") != 0) {
		hyp_metJes30_branch = tree->GetBranch(tree->GetAlias("hyp_metJes30"));
		hyp_metJes30_branch->SetAddress(&hyp_metJes30_);
	}
	if(hyp_metJes30_branch == 0 ) {
	cout << "Branch hyp_metJes30 does not exist." << endl;
	}
	hyp_metJes5_branch = 0;
	if (tree->GetAlias("hyp_metJes5") != 0) {
		hyp_metJes5_branch = tree->GetBranch(tree->GetAlias("hyp_metJes5"));
		hyp_metJes5_branch->SetAddress(&hyp_metJes5_);
	}
	if(hyp_metJes5_branch == 0 ) {
	cout << "Branch hyp_metJes5 does not exist." << endl;
	}
	hyp_metJes50_branch = 0;
	if (tree->GetAlias("hyp_metJes50") != 0) {
		hyp_metJes50_branch = tree->GetBranch(tree->GetAlias("hyp_metJes50"));
		hyp_metJes50_branch->SetAddress(&hyp_metJes50_);
	}
	if(hyp_metJes50_branch == 0 ) {
	cout << "Branch hyp_metJes50 does not exist." << endl;
	}
	hyp_metNoCalo_branch = 0;
	if (tree->GetAlias("hyp_metNoCalo") != 0) {
		hyp_metNoCalo_branch = tree->GetBranch(tree->GetAlias("hyp_metNoCalo"));
		hyp_metNoCalo_branch->SetAddress(&hyp_metNoCalo_);
	}
	if(hyp_metNoCalo_branch == 0 ) {
	cout << "Branch hyp_metNoCalo does not exist." << endl;
	}
	hyp_metPhi_branch = 0;
	if (tree->GetAlias("hyp_metPhi") != 0) {
		hyp_metPhi_branch = tree->GetBranch(tree->GetAlias("hyp_metPhi"));
		hyp_metPhi_branch->SetAddress(&hyp_metPhi_);
	}
	if(hyp_metPhi_branch == 0 ) {
	cout << "Branch hyp_metPhi does not exist." << endl;
	}
	hyp_metPhiAll_branch = 0;
	if (tree->GetAlias("hyp_metPhiAll") != 0) {
		hyp_metPhiAll_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiAll"));
		hyp_metPhiAll_branch->SetAddress(&hyp_metPhiAll_);
	}
	if(hyp_metPhiAll_branch == 0 ) {
	cout << "Branch hyp_metPhiAll does not exist." << endl;
	}
	hyp_metPhiAllCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metPhiAllCaloExp") != 0) {
		hyp_metPhiAllCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiAllCaloExp"));
		hyp_metPhiAllCaloExp_branch->SetAddress(&hyp_metPhiAllCaloExp_);
	}
	if(hyp_metPhiAllCaloExp_branch == 0 ) {
	cout << "Branch hyp_metPhiAllCaloExp does not exist." << endl;
	}
	hyp_metPhiCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metPhiCaloExp") != 0) {
		hyp_metPhiCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiCaloExp"));
		hyp_metPhiCaloExp_branch->SetAddress(&hyp_metPhiCaloExp_);
	}
	if(hyp_metPhiCaloExp_branch == 0 ) {
	cout << "Branch hyp_metPhiCaloExp does not exist." << endl;
	}
	hyp_metPhiCone_branch = 0;
	if (tree->GetAlias("hyp_metPhiCone") != 0) {
		hyp_metPhiCone_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiCone"));
		hyp_metPhiCone_branch->SetAddress(&hyp_metPhiCone_);
	}
	if(hyp_metPhiCone_branch == 0 ) {
	cout << "Branch hyp_metPhiCone does not exist." << endl;
	}
	hyp_metPhiJes10_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes10") != 0) {
		hyp_metPhiJes10_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes10"));
		hyp_metPhiJes10_branch->SetAddress(&hyp_metPhiJes10_);
	}
	if(hyp_metPhiJes10_branch == 0 ) {
	cout << "Branch hyp_metPhiJes10 does not exist." << endl;
	}
	hyp_metPhiJes15_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes15") != 0) {
		hyp_metPhiJes15_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes15"));
		hyp_metPhiJes15_branch->SetAddress(&hyp_metPhiJes15_);
	}
	if(hyp_metPhiJes15_branch == 0 ) {
	cout << "Branch hyp_metPhiJes15 does not exist." << endl;
	}
	hyp_metPhiJes30_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes30") != 0) {
		hyp_metPhiJes30_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes30"));
		hyp_metPhiJes30_branch->SetAddress(&hyp_metPhiJes30_);
	}
	if(hyp_metPhiJes30_branch == 0 ) {
	cout << "Branch hyp_metPhiJes30 does not exist." << endl;
	}
	hyp_metPhiJes5_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes5") != 0) {
		hyp_metPhiJes5_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes5"));
		hyp_metPhiJes5_branch->SetAddress(&hyp_metPhiJes5_);
	}
	if(hyp_metPhiJes5_branch == 0 ) {
	cout << "Branch hyp_metPhiJes5 does not exist." << endl;
	}
	hyp_metPhiJes50_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes50") != 0) {
		hyp_metPhiJes50_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes50"));
		hyp_metPhiJes50_branch->SetAddress(&hyp_metPhiJes50_);
	}
	if(hyp_metPhiJes50_branch == 0 ) {
	cout << "Branch hyp_metPhiJes50 does not exist." << endl;
	}
	hyp_metPhiNoCalo_branch = 0;
	if (tree->GetAlias("hyp_metPhiNoCalo") != 0) {
		hyp_metPhiNoCalo_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiNoCalo"));
		hyp_metPhiNoCalo_branch->SetAddress(&hyp_metPhiNoCalo_);
	}
	if(hyp_metPhiNoCalo_branch == 0 ) {
	cout << "Branch hyp_metPhiNoCalo does not exist." << endl;
	}
	hyp_quadlep_met_branch = 0;
	if (tree->GetAlias("hyp_quadlep_met") != 0) {
		hyp_quadlep_met_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_met"));
		hyp_quadlep_met_branch->SetAddress(&hyp_quadlep_met_);
	}
	if(hyp_quadlep_met_branch == 0 ) {
	cout << "Branch hyp_quadlep_met does not exist." << endl;
	}
	hyp_quadlep_metAll_branch = 0;
	if (tree->GetAlias("hyp_quadlep_metAll") != 0) {
		hyp_quadlep_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_metAll"));
		hyp_quadlep_metAll_branch->SetAddress(&hyp_quadlep_metAll_);
	}
	if(hyp_quadlep_metAll_branch == 0 ) {
	cout << "Branch hyp_quadlep_metAll does not exist." << endl;
	}
	hyp_trilep_met_branch = 0;
	if (tree->GetAlias("hyp_trilep_met") != 0) {
		hyp_trilep_met_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_met"));
		hyp_trilep_met_branch->SetAddress(&hyp_trilep_met_);
	}
	if(hyp_trilep_met_branch == 0 ) {
	cout << "Branch hyp_trilep_met does not exist." << endl;
	}
	hyp_trilep_metAll_branch = 0;
	if (tree->GetAlias("hyp_trilep_metAll") != 0) {
		hyp_trilep_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_metAll"));
		hyp_trilep_metAll_branch->SetAddress(&hyp_trilep_metAll_);
	}
	if(hyp_trilep_metAll_branch == 0 ) {
	cout << "Branch hyp_trilep_metAll does not exist." << endl;
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
	jpts_chFrac_branch = 0;
	if (tree->GetAlias("jpts_chFrac") != 0) {
		jpts_chFrac_branch = tree->GetBranch(tree->GetAlias("jpts_chFrac"));
		jpts_chFrac_branch->SetAddress(&jpts_chFrac_);
	}
	if(jpts_chFrac_branch == 0 ) {
	cout << "Branch jpts_chFrac does not exist." << endl;
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
	ww_pmet_branch = 0;
	if (tree->GetAlias("ww_pmet") != 0) {
		ww_pmet_branch = tree->GetBranch(tree->GetAlias("ww_pmet"));
		ww_pmet_branch->SetAddress(&ww_pmet_);
	}
	if(ww_pmet_branch == 0 ) {
	cout << "Branch ww_pmet does not exist." << endl;
	}
	hyp_jets_EMFcor_branch = 0;
	if (tree->GetAlias("hyp_jets_EMFcor") != 0) {
		hyp_jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("hyp_jets_EMFcor"));
		hyp_jets_EMFcor_branch->SetAddress(&hyp_jets_EMFcor_);
	}
	if(hyp_jets_EMFcor_branch == 0 ) {
	cout << "Branch hyp_jets_EMFcor does not exist." << endl;
	}
	hyp_jets_chFrac_branch = 0;
	if (tree->GetAlias("hyp_jets_chFrac") != 0) {
		hyp_jets_chFrac_branch = tree->GetBranch(tree->GetAlias("hyp_jets_chFrac"));
		hyp_jets_chFrac_branch->SetAddress(&hyp_jets_chFrac_);
	}
	if(hyp_jets_chFrac_branch == 0 ) {
	cout << "Branch hyp_jets_chFrac does not exist." << endl;
	}
	hyp_jets_cor_branch = 0;
	if (tree->GetAlias("hyp_jets_cor") != 0) {
		hyp_jets_cor_branch = tree->GetBranch(tree->GetAlias("hyp_jets_cor"));
		hyp_jets_cor_branch->SetAddress(&hyp_jets_cor_);
	}
	if(hyp_jets_cor_branch == 0 ) {
	cout << "Branch hyp_jets_cor does not exist." << endl;
	}
	hyp_jets_emFrac_branch = 0;
	if (tree->GetAlias("hyp_jets_emFrac") != 0) {
		hyp_jets_emFrac_branch = tree->GetBranch(tree->GetAlias("hyp_jets_emFrac"));
		hyp_jets_emFrac_branch->SetAddress(&hyp_jets_emFrac_);
	}
	if(hyp_jets_emFrac_branch == 0 ) {
	cout << "Branch hyp_jets_emFrac does not exist." << endl;
	}
	hyp_jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_emEnergy") != 0) {
		hyp_jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_emEnergy"));
		hyp_jets_mc_emEnergy_branch->SetAddress(&hyp_jets_mc_emEnergy_);
	}
	if(hyp_jets_mc_emEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_emEnergy does not exist." << endl;
	}
	hyp_jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_hadEnergy") != 0) {
		hyp_jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_hadEnergy"));
		hyp_jets_mc_hadEnergy_branch->SetAddress(&hyp_jets_mc_hadEnergy_);
	}
	if(hyp_jets_mc_hadEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_hadEnergy does not exist." << endl;
	}
	hyp_jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_invEnergy") != 0) {
		hyp_jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_invEnergy"));
		hyp_jets_mc_invEnergy_branch->SetAddress(&hyp_jets_mc_invEnergy_);
	}
	if(hyp_jets_mc_invEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_invEnergy does not exist." << endl;
	}
	hyp_jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_otherEnergy") != 0) {
		hyp_jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_otherEnergy"));
		hyp_jets_mc_otherEnergy_branch->SetAddress(&hyp_jets_mc_otherEnergy_);
	}
	if(hyp_jets_mc_otherEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_otherEnergy does not exist." << endl;
	}
	hyp_jets_pat_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_bCorrF") != 0) {
		hyp_jets_pat_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_bCorrF"));
		hyp_jets_pat_bCorrF_branch->SetAddress(&hyp_jets_pat_bCorrF_);
	}
	if(hyp_jets_pat_bCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_pat_bCorrF does not exist." << endl;
	}
	hyp_jets_pat_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_cCorrF") != 0) {
		hyp_jets_pat_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_cCorrF"));
		hyp_jets_pat_cCorrF_branch->SetAddress(&hyp_jets_pat_cCorrF_);
	}
	if(hyp_jets_pat_cCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_pat_cCorrF does not exist." << endl;
	}
	hyp_jets_pat_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_gluCorrF") != 0) {
		hyp_jets_pat_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_gluCorrF"));
		hyp_jets_pat_gluCorrF_branch->SetAddress(&hyp_jets_pat_gluCorrF_);
	}
	if(hyp_jets_pat_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_pat_gluCorrF does not exist." << endl;
	}
	hyp_jets_pat_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_jetCharge") != 0) {
		hyp_jets_pat_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_jetCharge"));
		hyp_jets_pat_jetCharge_branch->SetAddress(&hyp_jets_pat_jetCharge_);
	}
	if(hyp_jets_pat_jetCharge_branch == 0 ) {
	cout << "Branch hyp_jets_pat_jetCharge does not exist." << endl;
	}
	hyp_jets_pat_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_noCorrF") != 0) {
		hyp_jets_pat_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_noCorrF"));
		hyp_jets_pat_noCorrF_branch->SetAddress(&hyp_jets_pat_noCorrF_);
	}
	if(hyp_jets_pat_noCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_pat_noCorrF does not exist." << endl;
	}
	hyp_jets_pat_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_udsCorrF") != 0) {
		hyp_jets_pat_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_udsCorrF"));
		hyp_jets_pat_udsCorrF_branch->SetAddress(&hyp_jets_pat_udsCorrF_);
	}
	if(hyp_jets_pat_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_pat_udsCorrF does not exist." << endl;
	}
	hyp_other_jets_EMFcor_branch = 0;
	if (tree->GetAlias("hyp_other_jets_EMFcor") != 0) {
		hyp_other_jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_EMFcor"));
		hyp_other_jets_EMFcor_branch->SetAddress(&hyp_other_jets_EMFcor_);
	}
	if(hyp_other_jets_EMFcor_branch == 0 ) {
	cout << "Branch hyp_other_jets_EMFcor does not exist." << endl;
	}
	hyp_other_jets_chFrac_branch = 0;
	if (tree->GetAlias("hyp_other_jets_chFrac") != 0) {
		hyp_other_jets_chFrac_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_chFrac"));
		hyp_other_jets_chFrac_branch->SetAddress(&hyp_other_jets_chFrac_);
	}
	if(hyp_other_jets_chFrac_branch == 0 ) {
	cout << "Branch hyp_other_jets_chFrac does not exist." << endl;
	}
	hyp_other_jets_cor_branch = 0;
	if (tree->GetAlias("hyp_other_jets_cor") != 0) {
		hyp_other_jets_cor_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_cor"));
		hyp_other_jets_cor_branch->SetAddress(&hyp_other_jets_cor_);
	}
	if(hyp_other_jets_cor_branch == 0 ) {
	cout << "Branch hyp_other_jets_cor does not exist." << endl;
	}
	hyp_other_jets_emFrac_branch = 0;
	if (tree->GetAlias("hyp_other_jets_emFrac") != 0) {
		hyp_other_jets_emFrac_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_emFrac"));
		hyp_other_jets_emFrac_branch->SetAddress(&hyp_other_jets_emFrac_);
	}
	if(hyp_other_jets_emFrac_branch == 0 ) {
	cout << "Branch hyp_other_jets_emFrac does not exist." << endl;
	}
	hyp_other_jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_emEnergy") != 0) {
		hyp_other_jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_emEnergy"));
		hyp_other_jets_mc_emEnergy_branch->SetAddress(&hyp_other_jets_mc_emEnergy_);
	}
	if(hyp_other_jets_mc_emEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_emEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_hadEnergy") != 0) {
		hyp_other_jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_hadEnergy"));
		hyp_other_jets_mc_hadEnergy_branch->SetAddress(&hyp_other_jets_mc_hadEnergy_);
	}
	if(hyp_other_jets_mc_hadEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_hadEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_invEnergy") != 0) {
		hyp_other_jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_invEnergy"));
		hyp_other_jets_mc_invEnergy_branch->SetAddress(&hyp_other_jets_mc_invEnergy_);
	}
	if(hyp_other_jets_mc_invEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_invEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_otherEnergy") != 0) {
		hyp_other_jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_otherEnergy"));
		hyp_other_jets_mc_otherEnergy_branch->SetAddress(&hyp_other_jets_mc_otherEnergy_);
	}
	if(hyp_other_jets_mc_otherEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_otherEnergy does not exist." << endl;
	}
	hyp_other_jets_pat_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_bCorrF") != 0) {
		hyp_other_jets_pat_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_bCorrF"));
		hyp_other_jets_pat_bCorrF_branch->SetAddress(&hyp_other_jets_pat_bCorrF_);
	}
	if(hyp_other_jets_pat_bCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_bCorrF does not exist." << endl;
	}
	hyp_other_jets_pat_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_cCorrF") != 0) {
		hyp_other_jets_pat_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_cCorrF"));
		hyp_other_jets_pat_cCorrF_branch->SetAddress(&hyp_other_jets_pat_cCorrF_);
	}
	if(hyp_other_jets_pat_cCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_cCorrF does not exist." << endl;
	}
	hyp_other_jets_pat_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_gluCorrF") != 0) {
		hyp_other_jets_pat_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_gluCorrF"));
		hyp_other_jets_pat_gluCorrF_branch->SetAddress(&hyp_other_jets_pat_gluCorrF_);
	}
	if(hyp_other_jets_pat_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_gluCorrF does not exist." << endl;
	}
	hyp_other_jets_pat_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_jetCharge") != 0) {
		hyp_other_jets_pat_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_jetCharge"));
		hyp_other_jets_pat_jetCharge_branch->SetAddress(&hyp_other_jets_pat_jetCharge_);
	}
	if(hyp_other_jets_pat_jetCharge_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_jetCharge does not exist." << endl;
	}
	hyp_other_jets_pat_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_noCorrF") != 0) {
		hyp_other_jets_pat_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_noCorrF"));
		hyp_other_jets_pat_noCorrF_branch->SetAddress(&hyp_other_jets_pat_noCorrF_);
	}
	if(hyp_other_jets_pat_noCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_noCorrF does not exist." << endl;
	}
	hyp_other_jets_pat_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_udsCorrF") != 0) {
		hyp_other_jets_pat_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_udsCorrF"));
		hyp_other_jets_pat_udsCorrF_branch->SetAddress(&hyp_other_jets_pat_udsCorrF_);
	}
	if(hyp_other_jets_pat_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_udsCorrF does not exist." << endl;
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
	ww_isdyee_branch = 0;
	if (tree->GetAlias("ww_isdyee") != 0) {
		ww_isdyee_branch = tree->GetBranch(tree->GetAlias("ww_isdyee"));
		ww_isdyee_branch->SetAddress(&ww_isdyee_);
	}
	if(ww_isdyee_branch == 0 ) {
	cout << "Branch ww_isdyee does not exist." << endl;
	}
	ww_isdymm_branch = 0;
	if (tree->GetAlias("ww_isdymm") != 0) {
		ww_isdymm_branch = tree->GetBranch(tree->GetAlias("ww_isdymm"));
		ww_isdymm_branch->SetAddress(&ww_isdymm_);
	}
	if(ww_isdymm_branch == 0 ) {
	cout << "Branch ww_isdymm does not exist." << endl;
	}
	ww_isdytt_branch = 0;
	if (tree->GetAlias("ww_isdytt") != 0) {
		ww_isdytt_branch = tree->GetBranch(tree->GetAlias("ww_isdytt"));
		ww_isdytt_branch->SetAddress(&ww_isdytt_);
	}
	if(ww_isdytt_branch == 0 ) {
	cout << "Branch ww_isdytt does not exist." << endl;
	}
	ww_llgoodel_branch = 0;
	if (tree->GetAlias("ww_llgoodel") != 0) {
		ww_llgoodel_branch = tree->GetBranch(tree->GetAlias("ww_llgoodel"));
		ww_llgoodel_branch->SetAddress(&ww_llgoodel_);
	}
	if(ww_llgoodel_branch == 0 ) {
	cout << "Branch ww_llgoodel does not exist." << endl;
	}
	ww_llgoodeliso_branch = 0;
	if (tree->GetAlias("ww_llgoodeliso") != 0) {
		ww_llgoodeliso_branch = tree->GetBranch(tree->GetAlias("ww_llgoodeliso"));
		ww_llgoodeliso_branch->SetAddress(&ww_llgoodeliso_);
	}
	if(ww_llgoodeliso_branch == 0 ) {
	cout << "Branch ww_llgoodeliso does not exist." << endl;
	}
	ww_llgoodmu_branch = 0;
	if (tree->GetAlias("ww_llgoodmu") != 0) {
		ww_llgoodmu_branch = tree->GetBranch(tree->GetAlias("ww_llgoodmu"));
		ww_llgoodmu_branch->SetAddress(&ww_llgoodmu_);
	}
	if(ww_llgoodmu_branch == 0 ) {
	cout << "Branch ww_llgoodmu does not exist." << endl;
	}
	ww_llgoodmuiso_branch = 0;
	if (tree->GetAlias("ww_llgoodmuiso") != 0) {
		ww_llgoodmuiso_branch = tree->GetBranch(tree->GetAlias("ww_llgoodmuiso"));
		ww_llgoodmuiso_branch->SetAddress(&ww_llgoodmuiso_);
	}
	if(ww_llgoodmuiso_branch == 0 ) {
	cout << "Branch ww_llgoodmuiso does not exist." << endl;
	}
	ww_ltgoodel_branch = 0;
	if (tree->GetAlias("ww_ltgoodel") != 0) {
		ww_ltgoodel_branch = tree->GetBranch(tree->GetAlias("ww_ltgoodel"));
		ww_ltgoodel_branch->SetAddress(&ww_ltgoodel_);
	}
	if(ww_ltgoodel_branch == 0 ) {
	cout << "Branch ww_ltgoodel does not exist." << endl;
	}
	ww_ltgoodeliso_branch = 0;
	if (tree->GetAlias("ww_ltgoodeliso") != 0) {
		ww_ltgoodeliso_branch = tree->GetBranch(tree->GetAlias("ww_ltgoodeliso"));
		ww_ltgoodeliso_branch->SetAddress(&ww_ltgoodeliso_);
	}
	if(ww_ltgoodeliso_branch == 0 ) {
	cout << "Branch ww_ltgoodeliso does not exist." << endl;
	}
	ww_ltgoodmu_branch = 0;
	if (tree->GetAlias("ww_ltgoodmu") != 0) {
		ww_ltgoodmu_branch = tree->GetBranch(tree->GetAlias("ww_ltgoodmu"));
		ww_ltgoodmu_branch->SetAddress(&ww_ltgoodmu_);
	}
	if(ww_ltgoodmu_branch == 0 ) {
	cout << "Branch ww_ltgoodmu does not exist." << endl;
	}
	ww_ltgoodmuiso_branch = 0;
	if (tree->GetAlias("ww_ltgoodmuiso") != 0) {
		ww_ltgoodmuiso_branch = tree->GetBranch(tree->GetAlias("ww_ltgoodmuiso"));
		ww_ltgoodmuiso_branch->SetAddress(&ww_ltgoodmuiso_);
	}
	if(ww_ltgoodmuiso_branch == 0 ) {
	cout << "Branch ww_ltgoodmuiso does not exist." << endl;
	}
	ww_oppsign_branch = 0;
	if (tree->GetAlias("ww_oppsign") != 0) {
		ww_oppsign_branch = tree->GetBranch(tree->GetAlias("ww_oppsign"));
		ww_oppsign_branch->SetAddress(&ww_oppsign_);
	}
	if(ww_oppsign_branch == 0 ) {
	cout << "Branch ww_oppsign does not exist." << endl;
	}
	ww_pass2met_branch = 0;
	if (tree->GetAlias("ww_pass2met") != 0) {
		ww_pass2met_branch = tree->GetBranch(tree->GetAlias("ww_pass2met"));
		ww_pass2met_branch->SetAddress(&ww_pass2met_);
	}
	if(ww_pass2met_branch == 0 ) {
	cout << "Branch ww_pass2met does not exist." << endl;
	}
	ww_pass4met_branch = 0;
	if (tree->GetAlias("ww_pass4met") != 0) {
		ww_pass4met_branch = tree->GetBranch(tree->GetAlias("ww_pass4met"));
		ww_pass4met_branch->SetAddress(&ww_pass4met_);
	}
	if(ww_pass4met_branch == 0 ) {
	cout << "Branch ww_pass4met does not exist." << endl;
	}
	ww_passFebsel_no_jetveto_branch = 0;
	if (tree->GetAlias("ww_passFebsel_no_jetveto") != 0) {
		ww_passFebsel_no_jetveto_branch = tree->GetBranch(tree->GetAlias("ww_passFebsel_no_jetveto"));
		ww_passFebsel_no_jetveto_branch->SetAddress(&ww_passFebsel_no_jetveto_);
	}
	if(ww_passFebsel_no_jetveto_branch == 0 ) {
	cout << "Branch ww_passFebsel_no_jetveto does not exist." << endl;
	}
	ww_passFebsel_with_jetveto_branch = 0;
	if (tree->GetAlias("ww_passFebsel_with_jetveto") != 0) {
		ww_passFebsel_with_jetveto_branch = tree->GetBranch(tree->GetAlias("ww_passFebsel_with_jetveto"));
		ww_passFebsel_with_jetveto_branch->SetAddress(&ww_passFebsel_with_jetveto_);
	}
	if(ww_passFebsel_with_jetveto_branch == 0 ) {
	cout << "Branch ww_passFebsel_with_jetveto does not exist." << endl;
	}
	ww_passaddzveto_branch = 0;
	if (tree->GetAlias("ww_passaddzveto") != 0) {
		ww_passaddzveto_branch = tree->GetBranch(tree->GetAlias("ww_passaddzveto"));
		ww_passaddzveto_branch->SetAddress(&ww_passaddzveto_);
	}
	if(ww_passaddzveto_branch == 0 ) {
	cout << "Branch ww_passaddzveto does not exist." << endl;
	}
	ww_passzveto_branch = 0;
	if (tree->GetAlias("ww_passzveto") != 0) {
		ww_passzveto_branch = tree->GetBranch(tree->GetAlias("ww_passzveto"));
		ww_passzveto_branch->SetAddress(&ww_passzveto_);
	}
	if(ww_passzveto_branch == 0 ) {
	cout << "Branch ww_passzveto does not exist." << endl;
	}
	hyp_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_id") != 0) {
		hyp_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_id"));
		hyp_jets_mc_id_branch->SetAddress(&hyp_jets_mc_id_);
	}
	if(hyp_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_jets_mc_id does not exist." << endl;
	}
	hyp_jets_pat_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_genPartonMother_id") != 0) {
		hyp_jets_pat_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_genPartonMother_id"));
		hyp_jets_pat_genPartonMother_id_branch->SetAddress(&hyp_jets_pat_genPartonMother_id_);
	}
	if(hyp_jets_pat_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_jets_pat_genPartonMother_id does not exist." << endl;
	}
	hyp_jets_pat_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_genParton_id") != 0) {
		hyp_jets_pat_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_genParton_id"));
		hyp_jets_pat_genParton_id_branch->SetAddress(&hyp_jets_pat_genParton_id_);
	}
	if(hyp_jets_pat_genParton_id_branch == 0 ) {
	cout << "Branch hyp_jets_pat_genParton_id does not exist." << endl;
	}
	hyp_jets_pat_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_partonFlavour") != 0) {
		hyp_jets_pat_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_partonFlavour"));
		hyp_jets_pat_partonFlavour_branch->SetAddress(&hyp_jets_pat_partonFlavour_);
	}
	if(hyp_jets_pat_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_jets_pat_partonFlavour does not exist." << endl;
	}
	hyp_other_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_id") != 0) {
		hyp_other_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_id"));
		hyp_other_jets_mc_id_branch->SetAddress(&hyp_other_jets_mc_id_);
	}
	if(hyp_other_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_id does not exist." << endl;
	}
	hyp_other_jets_pat_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_genPartonMother_id") != 0) {
		hyp_other_jets_pat_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_genPartonMother_id"));
		hyp_other_jets_pat_genPartonMother_id_branch->SetAddress(&hyp_other_jets_pat_genPartonMother_id_);
	}
	if(hyp_other_jets_pat_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_genPartonMother_id does not exist." << endl;
	}
	hyp_other_jets_pat_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_genParton_id") != 0) {
		hyp_other_jets_pat_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_genParton_id"));
		hyp_other_jets_pat_genParton_id_branch->SetAddress(&hyp_other_jets_pat_genParton_id_);
	}
	if(hyp_other_jets_pat_genParton_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_genParton_id does not exist." << endl;
	}
	hyp_other_jets_pat_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_other_jets_pat_partonFlavour") != 0) {
		hyp_other_jets_pat_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_pat_partonFlavour"));
		hyp_other_jets_pat_partonFlavour_branch->SetAddress(&hyp_other_jets_pat_partonFlavour_);
	}
	if(hyp_other_jets_pat_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_other_jets_pat_partonFlavour does not exist." << endl;
	}
	hyp_quadlep_jets_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_jets_index") != 0) {
		hyp_quadlep_jets_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_jets_index"));
		hyp_quadlep_jets_index_branch->SetAddress(&hyp_quadlep_jets_index_);
	}
	if(hyp_quadlep_jets_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_jets_index does not exist." << endl;
	}
	hyp_trilep_jets_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_jets_index") != 0) {
		hyp_trilep_jets_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_jets_index"));
		hyp_trilep_jets_index_branch->SetAddress(&hyp_trilep_jets_index_);
	}
	if(hyp_trilep_jets_index_branch == 0 ) {
	cout << "Branch hyp_trilep_jets_index does not exist." << endl;
	}
	evt_nalltrkjets_branch = 0;
	if (tree->GetAlias("evt_nalltrkjets") != 0) {
		evt_nalltrkjets_branch = tree->GetBranch(tree->GetAlias("evt_nalltrkjets"));
		evt_nalltrkjets_branch->SetAddress(&evt_nalltrkjets_);
	}
	if(evt_nalltrkjets_branch == 0 ) {
	cout << "Branch evt_nalltrkjets does not exist." << endl;
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
		scs_pos_isLoaded = false;
		scs_vtx_isLoaded = false;
		vtxs_position_isLoaded = false;
		evt_bsp4_isLoaded = false;
		l1met_p4_isLoaded = false;
		alltrkjets_p4_isLoaded = false;
		els_mc3_p4_isLoaded = false;
		els_mc_p4_isLoaded = false;
		jets_mc_gp_p4_isLoaded = false;
		jets_mc_p4_isLoaded = false;
		mus_mc3_p4_isLoaded = false;
		mus_mc_p4_isLoaded = false;
		trk_mc3p4_isLoaded = false;
		trk_mcp4_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
		els_vertex_p4_isLoaded = false;
		genjets_p4_isLoaded = false;
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
		mus_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		mus_vertex_p4_isLoaded = false;
		els_pat_genMotherP4_isLoaded = false;
		els_pat_genP4_isLoaded = false;
		jets_pat_genJet_p4_isLoaded = false;
		jets_pat_genPartonMother_p4_isLoaded = false;
		jets_pat_genParton_p4_isLoaded = false;
		jets_pat_jet_p4_isLoaded = false;
		jets_pat_jet_uncorp4_isLoaded = false;
		mus_pat_genMotherP4_isLoaded = false;
		mus_pat_genP4_isLoaded = false;
		scs_p4_isLoaded = false;
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
		hyp_jets_mc_gp_p4_isLoaded = false;
		hyp_jets_mc_p4_isLoaded = false;
		hyp_jets_p4_isLoaded = false;
		hyp_jets_pat_genJet_p4_isLoaded = false;
		hyp_jets_pat_genPartonMother_p4_isLoaded = false;
		hyp_jets_pat_genParton_p4_isLoaded = false;
		hyp_jets_pat_jet_p4_isLoaded = false;
		hyp_other_jets_mc_gp_p4_isLoaded = false;
		hyp_other_jets_mc_p4_isLoaded = false;
		hyp_other_jets_p4_isLoaded = false;
		hyp_other_jets_pat_genJet_p4_isLoaded = false;
		hyp_other_jets_pat_genPartonMother_p4_isLoaded = false;
		hyp_other_jets_pat_genParton_p4_isLoaded = false;
		hyp_other_jets_pat_jet_p4_isLoaded = false;
		jets_closestElectron_DR_isLoaded = false;
		jets_closestMuon_DR_isLoaded = false;
		evt_bs_dxdz_isLoaded = false;
		evt_bs_dxdzErr_isLoaded = false;
		evt_bs_dydz_isLoaded = false;
		evt_bs_dydzErr_isLoaded = false;
		evt_bs_sigmaZ_isLoaded = false;
		evt_bs_sigmaZErr_isLoaded = false;
		evt_bs_widthErr_isLoaded = false;
		evt_bs_xErr_isLoaded = false;
		evt_bs_yErr_isLoaded = false;
		evt_bs_zErr_isLoaded = false;
		gen_met_isLoaded = false;
		gen_metPhi_isLoaded = false;
		evt_bField_isLoaded = false;
		evt_kfactor_isLoaded = false;
		evt_weight_isLoaded = false;
		evt_xsec_excl_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
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
		evt_metNoHF_isLoaded = false;
		evt_metNoHFHO_isLoaded = false;
		evt_metNoHFHOPhi_isLoaded = false;
		evt_metNoHFHOSig_isLoaded = false;
		evt_metNoHFPhi_isLoaded = false;
		evt_metSig_isLoaded = false;
		evt_metOpt_isLoaded = false;
		evt_metOptHO_isLoaded = false;
		evt_metOptHOPhi_isLoaded = false;
		evt_metOptHOSig_isLoaded = false;
		evt_metOptNoHF_isLoaded = false;
		evt_metOptNoHFHO_isLoaded = false;
		evt_metOptNoHFHOPhi_isLoaded = false;
		evt_metOptNoHFHOSig_isLoaded = false;
		evt_metOptNoHFPhi_isLoaded = false;
		evt_metOptSig_isLoaded = false;
		evt_metOptPhi_isLoaded = false;
		evt_metPhi_isLoaded = false;
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
		evt_tcmet_isLoaded = false;
		evt_tcmetPhi_isLoaded = false;
		evt_tcsumet_isLoaded = false;
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
		els_e2x5Max_isLoaded = false;
		els_e3x3_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_eMax_isLoaded = false;
		els_eOverPIn_isLoaded = false;
		els_eSC_isLoaded = false;
		els_eSCPresh_isLoaded = false;
		els_eSCRaw_isLoaded = false;
		els_eSeed_isLoaded = false;
		els_eSeedOverPOut_isLoaded = false;
		els_etaErr_isLoaded = false;
		els_fBrem_isLoaded = false;
		els_hOverE_isLoaded = false;
		els_ndof_isLoaded = false;
		els_outerEta_isLoaded = false;
		els_outerPhi_isLoaded = false;
		els_phiErr_isLoaded = false;
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
		hyp_ll_chi2_isLoaded = false;
		hyp_ll_d0_isLoaded = false;
		hyp_ll_d0Err_isLoaded = false;
		hyp_ll_d0corr_isLoaded = false;
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
		hyp_met_isLoaded = false;
		hyp_metAll_isLoaded = false;
		hyp_metAllCaloExp_isLoaded = false;
		hyp_metCaloExp_isLoaded = false;
		hyp_metCone_isLoaded = false;
		hyp_metDPhiJet10_isLoaded = false;
		hyp_metDPhiJet15_isLoaded = false;
		hyp_metDPhiJet20_isLoaded = false;
		hyp_metDPhiTrk10_isLoaded = false;
		hyp_metDPhiTrk25_isLoaded = false;
		hyp_metDPhiTrk50_isLoaded = false;
		hyp_metJes10_isLoaded = false;
		hyp_metJes15_isLoaded = false;
		hyp_metJes30_isLoaded = false;
		hyp_metJes5_isLoaded = false;
		hyp_metJes50_isLoaded = false;
		hyp_metNoCalo_isLoaded = false;
		hyp_metPhi_isLoaded = false;
		hyp_metPhiAll_isLoaded = false;
		hyp_metPhiAllCaloExp_isLoaded = false;
		hyp_metPhiCaloExp_isLoaded = false;
		hyp_metPhiCone_isLoaded = false;
		hyp_metPhiJes10_isLoaded = false;
		hyp_metPhiJes15_isLoaded = false;
		hyp_metPhiJes30_isLoaded = false;
		hyp_metPhiJes5_isLoaded = false;
		hyp_metPhiJes50_isLoaded = false;
		hyp_metPhiNoCalo_isLoaded = false;
		hyp_quadlep_met_isLoaded = false;
		hyp_quadlep_metAll_isLoaded = false;
		hyp_trilep_met_isLoaded = false;
		hyp_trilep_metAll_isLoaded = false;
		jets_EMFcor_isLoaded = false;
		jets_chFrac_isLoaded = false;
		jets_cor_isLoaded = false;
		jets_emFrac_isLoaded = false;
		jpts_chFrac_isLoaded = false;
		jpts_cor_isLoaded = false;
		jpts_emFrac_isLoaded = false;
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
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_d0corr_isLoaded = false;
		trks_etaErr_isLoaded = false;
		trks_ndof_isLoaded = false;
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
		ww_pmet_isLoaded = false;
		hyp_jets_EMFcor_isLoaded = false;
		hyp_jets_chFrac_isLoaded = false;
		hyp_jets_cor_isLoaded = false;
		hyp_jets_emFrac_isLoaded = false;
		hyp_jets_mc_emEnergy_isLoaded = false;
		hyp_jets_mc_hadEnergy_isLoaded = false;
		hyp_jets_mc_invEnergy_isLoaded = false;
		hyp_jets_mc_otherEnergy_isLoaded = false;
		hyp_jets_pat_bCorrF_isLoaded = false;
		hyp_jets_pat_cCorrF_isLoaded = false;
		hyp_jets_pat_gluCorrF_isLoaded = false;
		hyp_jets_pat_jetCharge_isLoaded = false;
		hyp_jets_pat_noCorrF_isLoaded = false;
		hyp_jets_pat_udsCorrF_isLoaded = false;
		hyp_other_jets_EMFcor_isLoaded = false;
		hyp_other_jets_chFrac_isLoaded = false;
		hyp_other_jets_cor_isLoaded = false;
		hyp_other_jets_emFrac_isLoaded = false;
		hyp_other_jets_mc_emEnergy_isLoaded = false;
		hyp_other_jets_mc_hadEnergy_isLoaded = false;
		hyp_other_jets_mc_invEnergy_isLoaded = false;
		hyp_other_jets_mc_otherEnergy_isLoaded = false;
		hyp_other_jets_pat_bCorrF_isLoaded = false;
		hyp_other_jets_pat_cCorrF_isLoaded = false;
		hyp_other_jets_pat_gluCorrF_isLoaded = false;
		hyp_other_jets_pat_jetCharge_isLoaded = false;
		hyp_other_jets_pat_noCorrF_isLoaded = false;
		hyp_other_jets_pat_udsCorrF_isLoaded = false;
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
		els_looseId_isLoaded = false;
		els_lostHits_isLoaded = false;
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
		genps_id_isLoaded = false;
		genps_id_mother_isLoaded = false;
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
		ww_isdyee_isLoaded = false;
		ww_isdymm_isLoaded = false;
		ww_isdytt_isLoaded = false;
		ww_llgoodel_isLoaded = false;
		ww_llgoodeliso_isLoaded = false;
		ww_llgoodmu_isLoaded = false;
		ww_llgoodmuiso_isLoaded = false;
		ww_ltgoodel_isLoaded = false;
		ww_ltgoodeliso_isLoaded = false;
		ww_ltgoodmu_isLoaded = false;
		ww_ltgoodmuiso_isLoaded = false;
		ww_oppsign_isLoaded = false;
		ww_pass2met_isLoaded = false;
		ww_pass4met_isLoaded = false;
		ww_passFebsel_no_jetveto_isLoaded = false;
		ww_passFebsel_with_jetveto_isLoaded = false;
		ww_passaddzveto_isLoaded = false;
		ww_passzveto_isLoaded = false;
		hyp_jets_mc_id_isLoaded = false;
		hyp_jets_pat_genPartonMother_id_isLoaded = false;
		hyp_jets_pat_genParton_id_isLoaded = false;
		hyp_jets_pat_partonFlavour_isLoaded = false;
		hyp_other_jets_mc_id_isLoaded = false;
		hyp_other_jets_pat_genPartonMother_id_isLoaded = false;
		hyp_other_jets_pat_genParton_id_isLoaded = false;
		hyp_other_jets_pat_partonFlavour_isLoaded = false;
		hyp_quadlep_jets_index_isLoaded = false;
		hyp_trilep_jets_index_isLoaded = false;
		evt_nalltrkjets_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_event_isLoaded = false;
		evt_lumiBlock_isLoaded = false;
		evt_run_isLoaded = false;
		evt_ngenjets_isLoaded = false;
		evt_njets_isLoaded = false;
		evt_njpts_isLoaded = false;
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
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &scs_pos()
	{
		if (not scs_pos_isLoaded) {
			if (scs_pos_branch != 0) {
				scs_pos_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_pos_branch does not exist!\n");
				exit(1);
			}
			scs_pos_isLoaded = true;
		}
		return scs_pos_;
	}
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &scs_vtx()
	{
		if (not scs_vtx_isLoaded) {
			if (scs_vtx_branch != 0) {
				scs_vtx_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch scs_vtx_branch does not exist!\n");
				exit(1);
			}
			scs_vtx_isLoaded = true;
		}
		return scs_vtx_;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &alltrkjets_p4()
	{
		if (not alltrkjets_p4_isLoaded) {
			if (alltrkjets_p4_branch != 0) {
				alltrkjets_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator i = alltrkjets_p4_.begin(); i != alltrkjets_p4_.end(); ++i) {
					int e;
					frexp(i->pt(), &e);
					if (not isfinite(i->pt()) || e > 30) {
						printf("branch alltrkjets_p4_branch contains a bad float: %f\n", i->pt());
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch alltrkjets_p4_branch does not exist!\n");
				exit(1);
			}
			alltrkjets_p4_isLoaded = true;
		}
		return alltrkjets_p4_;
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_gp_p4()
	{
		if (not hyp_jets_mc_gp_p4_isLoaded) {
			if (hyp_jets_mc_gp_p4_branch != 0) {
				hyp_jets_mc_gp_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_mc_gp_p4_.begin(); i != hyp_jets_mc_gp_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_mc_gp_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_gp_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_p4()
	{
		if (not hyp_jets_mc_p4_isLoaded) {
			if (hyp_jets_mc_p4_branch != 0) {
				hyp_jets_mc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_mc_p4_.begin(); i != hyp_jets_mc_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_mc_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_p4_isLoaded = true;
		}
		return hyp_jets_mc_p4_;
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genJet_p4()
	{
		if (not hyp_jets_pat_genJet_p4_isLoaded) {
			if (hyp_jets_pat_genJet_p4_branch != 0) {
				hyp_jets_pat_genJet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_pat_genJet_p4_.begin(); i != hyp_jets_pat_genJet_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_pat_genJet_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_genJet_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_genJet_p4_isLoaded = true;
		}
		return hyp_jets_pat_genJet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genPartonMother_p4()
	{
		if (not hyp_jets_pat_genPartonMother_p4_isLoaded) {
			if (hyp_jets_pat_genPartonMother_p4_branch != 0) {
				hyp_jets_pat_genPartonMother_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_pat_genPartonMother_p4_.begin(); i != hyp_jets_pat_genPartonMother_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_pat_genPartonMother_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_genPartonMother_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return hyp_jets_pat_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genParton_p4()
	{
		if (not hyp_jets_pat_genParton_p4_isLoaded) {
			if (hyp_jets_pat_genParton_p4_branch != 0) {
				hyp_jets_pat_genParton_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_pat_genParton_p4_.begin(); i != hyp_jets_pat_genParton_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_pat_genParton_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_genParton_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_genParton_p4_isLoaded = true;
		}
		return hyp_jets_pat_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_jet_p4()
	{
		if (not hyp_jets_pat_jet_p4_isLoaded) {
			if (hyp_jets_pat_jet_p4_branch != 0) {
				hyp_jets_pat_jet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_jets_pat_jet_p4_.begin(); i != hyp_jets_pat_jet_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_jets_pat_jet_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_jet_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_jet_p4_isLoaded = true;
		}
		return hyp_jets_pat_jet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_gp_p4()
	{
		if (not hyp_other_jets_mc_gp_p4_isLoaded) {
			if (hyp_other_jets_mc_gp_p4_branch != 0) {
				hyp_other_jets_mc_gp_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_mc_gp_p4_.begin(); i != hyp_other_jets_mc_gp_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_mc_gp_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_gp_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_p4()
	{
		if (not hyp_other_jets_mc_p4_isLoaded) {
			if (hyp_other_jets_mc_p4_branch != 0) {
				hyp_other_jets_mc_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_mc_p4_.begin(); i != hyp_other_jets_mc_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_mc_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_p4_;
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genJet_p4()
	{
		if (not hyp_other_jets_pat_genJet_p4_isLoaded) {
			if (hyp_other_jets_pat_genJet_p4_branch != 0) {
				hyp_other_jets_pat_genJet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_pat_genJet_p4_.begin(); i != hyp_other_jets_pat_genJet_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_pat_genJet_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_genJet_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_genJet_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genJet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genPartonMother_p4()
	{
		if (not hyp_other_jets_pat_genPartonMother_p4_isLoaded) {
			if (hyp_other_jets_pat_genPartonMother_p4_branch != 0) {
				hyp_other_jets_pat_genPartonMother_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_pat_genPartonMother_p4_.begin(); i != hyp_other_jets_pat_genPartonMother_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_pat_genPartonMother_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_genPartonMother_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genParton_p4()
	{
		if (not hyp_other_jets_pat_genParton_p4_isLoaded) {
			if (hyp_other_jets_pat_genParton_p4_branch != 0) {
				hyp_other_jets_pat_genParton_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_pat_genParton_p4_.begin(); i != hyp_other_jets_pat_genParton_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_pat_genParton_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_genParton_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_genParton_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_jet_p4()
	{
		if (not hyp_other_jets_pat_jet_p4_isLoaded) {
			if (hyp_other_jets_pat_jet_p4_branch != 0) {
				hyp_other_jets_pat_jet_p4_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >::const_iterator i = hyp_other_jets_pat_jet_p4_.begin(); i != hyp_other_jets_pat_jet_p4_.end(); ++i) {
					for (vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >::const_iterator j = i->begin(); j != i->end(); ++j) {
						int e;
						frexp(j->pt(), &e);
						if (not isfinite(j->pt()) || e > 30) {
							printf("branch hyp_other_jets_pat_jet_p4_branch contains a bad float: %f\n", j->pt());
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_jet_p4_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_jet_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_jet_p4_;
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
	vector<float> &hyp_met()
	{
		if (not hyp_met_isLoaded) {
			if (hyp_met_branch != 0) {
				hyp_met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_met_.begin(); i != hyp_met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_met_branch does not exist!\n");
				exit(1);
			}
			hyp_met_isLoaded = true;
		}
		return hyp_met_;
	}
	vector<float> &hyp_metAll()
	{
		if (not hyp_metAll_isLoaded) {
			if (hyp_metAll_branch != 0) {
				hyp_metAll_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metAll_.begin(); i != hyp_metAll_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metAll_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metAll_branch does not exist!\n");
				exit(1);
			}
			hyp_metAll_isLoaded = true;
		}
		return hyp_metAll_;
	}
	vector<float> &hyp_metAllCaloExp()
	{
		if (not hyp_metAllCaloExp_isLoaded) {
			if (hyp_metAllCaloExp_branch != 0) {
				hyp_metAllCaloExp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metAllCaloExp_.begin(); i != hyp_metAllCaloExp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metAllCaloExp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metAllCaloExp_branch does not exist!\n");
				exit(1);
			}
			hyp_metAllCaloExp_isLoaded = true;
		}
		return hyp_metAllCaloExp_;
	}
	vector<float> &hyp_metCaloExp()
	{
		if (not hyp_metCaloExp_isLoaded) {
			if (hyp_metCaloExp_branch != 0) {
				hyp_metCaloExp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metCaloExp_.begin(); i != hyp_metCaloExp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metCaloExp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metCaloExp_branch does not exist!\n");
				exit(1);
			}
			hyp_metCaloExp_isLoaded = true;
		}
		return hyp_metCaloExp_;
	}
	vector<float> &hyp_metCone()
	{
		if (not hyp_metCone_isLoaded) {
			if (hyp_metCone_branch != 0) {
				hyp_metCone_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metCone_.begin(); i != hyp_metCone_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metCone_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metCone_branch does not exist!\n");
				exit(1);
			}
			hyp_metCone_isLoaded = true;
		}
		return hyp_metCone_;
	}
	vector<float> &hyp_metDPhiJet10()
	{
		if (not hyp_metDPhiJet10_isLoaded) {
			if (hyp_metDPhiJet10_branch != 0) {
				hyp_metDPhiJet10_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiJet10_.begin(); i != hyp_metDPhiJet10_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiJet10_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiJet10_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiJet10_isLoaded = true;
		}
		return hyp_metDPhiJet10_;
	}
	vector<float> &hyp_metDPhiJet15()
	{
		if (not hyp_metDPhiJet15_isLoaded) {
			if (hyp_metDPhiJet15_branch != 0) {
				hyp_metDPhiJet15_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiJet15_.begin(); i != hyp_metDPhiJet15_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiJet15_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiJet15_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiJet15_isLoaded = true;
		}
		return hyp_metDPhiJet15_;
	}
	vector<float> &hyp_metDPhiJet20()
	{
		if (not hyp_metDPhiJet20_isLoaded) {
			if (hyp_metDPhiJet20_branch != 0) {
				hyp_metDPhiJet20_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiJet20_.begin(); i != hyp_metDPhiJet20_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiJet20_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiJet20_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiJet20_isLoaded = true;
		}
		return hyp_metDPhiJet20_;
	}
	vector<float> &hyp_metDPhiTrk10()
	{
		if (not hyp_metDPhiTrk10_isLoaded) {
			if (hyp_metDPhiTrk10_branch != 0) {
				hyp_metDPhiTrk10_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiTrk10_.begin(); i != hyp_metDPhiTrk10_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiTrk10_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiTrk10_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiTrk10_isLoaded = true;
		}
		return hyp_metDPhiTrk10_;
	}
	vector<float> &hyp_metDPhiTrk25()
	{
		if (not hyp_metDPhiTrk25_isLoaded) {
			if (hyp_metDPhiTrk25_branch != 0) {
				hyp_metDPhiTrk25_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiTrk25_.begin(); i != hyp_metDPhiTrk25_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiTrk25_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiTrk25_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiTrk25_isLoaded = true;
		}
		return hyp_metDPhiTrk25_;
	}
	vector<float> &hyp_metDPhiTrk50()
	{
		if (not hyp_metDPhiTrk50_isLoaded) {
			if (hyp_metDPhiTrk50_branch != 0) {
				hyp_metDPhiTrk50_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metDPhiTrk50_.begin(); i != hyp_metDPhiTrk50_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metDPhiTrk50_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metDPhiTrk50_branch does not exist!\n");
				exit(1);
			}
			hyp_metDPhiTrk50_isLoaded = true;
		}
		return hyp_metDPhiTrk50_;
	}
	vector<float> &hyp_metJes10()
	{
		if (not hyp_metJes10_isLoaded) {
			if (hyp_metJes10_branch != 0) {
				hyp_metJes10_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metJes10_.begin(); i != hyp_metJes10_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metJes10_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metJes10_branch does not exist!\n");
				exit(1);
			}
			hyp_metJes10_isLoaded = true;
		}
		return hyp_metJes10_;
	}
	vector<float> &hyp_metJes15()
	{
		if (not hyp_metJes15_isLoaded) {
			if (hyp_metJes15_branch != 0) {
				hyp_metJes15_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metJes15_.begin(); i != hyp_metJes15_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metJes15_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metJes15_branch does not exist!\n");
				exit(1);
			}
			hyp_metJes15_isLoaded = true;
		}
		return hyp_metJes15_;
	}
	vector<float> &hyp_metJes30()
	{
		if (not hyp_metJes30_isLoaded) {
			if (hyp_metJes30_branch != 0) {
				hyp_metJes30_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metJes30_.begin(); i != hyp_metJes30_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metJes30_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metJes30_branch does not exist!\n");
				exit(1);
			}
			hyp_metJes30_isLoaded = true;
		}
		return hyp_metJes30_;
	}
	vector<float> &hyp_metJes5()
	{
		if (not hyp_metJes5_isLoaded) {
			if (hyp_metJes5_branch != 0) {
				hyp_metJes5_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metJes5_.begin(); i != hyp_metJes5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metJes5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metJes5_branch does not exist!\n");
				exit(1);
			}
			hyp_metJes5_isLoaded = true;
		}
		return hyp_metJes5_;
	}
	vector<float> &hyp_metJes50()
	{
		if (not hyp_metJes50_isLoaded) {
			if (hyp_metJes50_branch != 0) {
				hyp_metJes50_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metJes50_.begin(); i != hyp_metJes50_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metJes50_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metJes50_branch does not exist!\n");
				exit(1);
			}
			hyp_metJes50_isLoaded = true;
		}
		return hyp_metJes50_;
	}
	vector<float> &hyp_metNoCalo()
	{
		if (not hyp_metNoCalo_isLoaded) {
			if (hyp_metNoCalo_branch != 0) {
				hyp_metNoCalo_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metNoCalo_.begin(); i != hyp_metNoCalo_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metNoCalo_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metNoCalo_branch does not exist!\n");
				exit(1);
			}
			hyp_metNoCalo_isLoaded = true;
		}
		return hyp_metNoCalo_;
	}
	vector<float> &hyp_metPhi()
	{
		if (not hyp_metPhi_isLoaded) {
			if (hyp_metPhi_branch != 0) {
				hyp_metPhi_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhi_.begin(); i != hyp_metPhi_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhi_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhi_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhi_isLoaded = true;
		}
		return hyp_metPhi_;
	}
	vector<float> &hyp_metPhiAll()
	{
		if (not hyp_metPhiAll_isLoaded) {
			if (hyp_metPhiAll_branch != 0) {
				hyp_metPhiAll_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiAll_.begin(); i != hyp_metPhiAll_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiAll_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiAll_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiAll_isLoaded = true;
		}
		return hyp_metPhiAll_;
	}
	vector<float> &hyp_metPhiAllCaloExp()
	{
		if (not hyp_metPhiAllCaloExp_isLoaded) {
			if (hyp_metPhiAllCaloExp_branch != 0) {
				hyp_metPhiAllCaloExp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiAllCaloExp_.begin(); i != hyp_metPhiAllCaloExp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiAllCaloExp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiAllCaloExp_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiAllCaloExp_isLoaded = true;
		}
		return hyp_metPhiAllCaloExp_;
	}
	vector<float> &hyp_metPhiCaloExp()
	{
		if (not hyp_metPhiCaloExp_isLoaded) {
			if (hyp_metPhiCaloExp_branch != 0) {
				hyp_metPhiCaloExp_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiCaloExp_.begin(); i != hyp_metPhiCaloExp_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiCaloExp_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiCaloExp_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiCaloExp_isLoaded = true;
		}
		return hyp_metPhiCaloExp_;
	}
	vector<float> &hyp_metPhiCone()
	{
		if (not hyp_metPhiCone_isLoaded) {
			if (hyp_metPhiCone_branch != 0) {
				hyp_metPhiCone_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiCone_.begin(); i != hyp_metPhiCone_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiCone_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiCone_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiCone_isLoaded = true;
		}
		return hyp_metPhiCone_;
	}
	vector<float> &hyp_metPhiJes10()
	{
		if (not hyp_metPhiJes10_isLoaded) {
			if (hyp_metPhiJes10_branch != 0) {
				hyp_metPhiJes10_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiJes10_.begin(); i != hyp_metPhiJes10_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiJes10_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiJes10_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiJes10_isLoaded = true;
		}
		return hyp_metPhiJes10_;
	}
	vector<float> &hyp_metPhiJes15()
	{
		if (not hyp_metPhiJes15_isLoaded) {
			if (hyp_metPhiJes15_branch != 0) {
				hyp_metPhiJes15_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiJes15_.begin(); i != hyp_metPhiJes15_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiJes15_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiJes15_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiJes15_isLoaded = true;
		}
		return hyp_metPhiJes15_;
	}
	vector<float> &hyp_metPhiJes30()
	{
		if (not hyp_metPhiJes30_isLoaded) {
			if (hyp_metPhiJes30_branch != 0) {
				hyp_metPhiJes30_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiJes30_.begin(); i != hyp_metPhiJes30_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiJes30_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiJes30_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiJes30_isLoaded = true;
		}
		return hyp_metPhiJes30_;
	}
	vector<float> &hyp_metPhiJes5()
	{
		if (not hyp_metPhiJes5_isLoaded) {
			if (hyp_metPhiJes5_branch != 0) {
				hyp_metPhiJes5_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiJes5_.begin(); i != hyp_metPhiJes5_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiJes5_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiJes5_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiJes5_isLoaded = true;
		}
		return hyp_metPhiJes5_;
	}
	vector<float> &hyp_metPhiJes50()
	{
		if (not hyp_metPhiJes50_isLoaded) {
			if (hyp_metPhiJes50_branch != 0) {
				hyp_metPhiJes50_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiJes50_.begin(); i != hyp_metPhiJes50_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiJes50_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiJes50_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiJes50_isLoaded = true;
		}
		return hyp_metPhiJes50_;
	}
	vector<float> &hyp_metPhiNoCalo()
	{
		if (not hyp_metPhiNoCalo_isLoaded) {
			if (hyp_metPhiNoCalo_branch != 0) {
				hyp_metPhiNoCalo_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_metPhiNoCalo_.begin(); i != hyp_metPhiNoCalo_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_metPhiNoCalo_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_metPhiNoCalo_branch does not exist!\n");
				exit(1);
			}
			hyp_metPhiNoCalo_isLoaded = true;
		}
		return hyp_metPhiNoCalo_;
	}
	vector<float> &hyp_quadlep_met()
	{
		if (not hyp_quadlep_met_isLoaded) {
			if (hyp_quadlep_met_branch != 0) {
				hyp_quadlep_met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_quadlep_met_.begin(); i != hyp_quadlep_met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_quadlep_met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_quadlep_met_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_met_isLoaded = true;
		}
		return hyp_quadlep_met_;
	}
	vector<float> &hyp_quadlep_metAll()
	{
		if (not hyp_quadlep_metAll_isLoaded) {
			if (hyp_quadlep_metAll_branch != 0) {
				hyp_quadlep_metAll_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_quadlep_metAll_.begin(); i != hyp_quadlep_metAll_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_quadlep_metAll_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_quadlep_metAll_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_metAll_isLoaded = true;
		}
		return hyp_quadlep_metAll_;
	}
	vector<float> &hyp_trilep_met()
	{
		if (not hyp_trilep_met_isLoaded) {
			if (hyp_trilep_met_branch != 0) {
				hyp_trilep_met_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_trilep_met_.begin(); i != hyp_trilep_met_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_trilep_met_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_trilep_met_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_met_isLoaded = true;
		}
		return hyp_trilep_met_;
	}
	vector<float> &hyp_trilep_metAll()
	{
		if (not hyp_trilep_metAll_isLoaded) {
			if (hyp_trilep_metAll_branch != 0) {
				hyp_trilep_metAll_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = hyp_trilep_metAll_.begin(); i != hyp_trilep_metAll_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch hyp_trilep_metAll_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_trilep_metAll_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_metAll_isLoaded = true;
		}
		return hyp_trilep_metAll_;
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
	vector<float> &jpts_chFrac()
	{
		if (not jpts_chFrac_isLoaded) {
			if (jpts_chFrac_branch != 0) {
				jpts_chFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = jpts_chFrac_.begin(); i != jpts_chFrac_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch jpts_chFrac_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch jpts_chFrac_branch does not exist!\n");
				exit(1);
			}
			jpts_chFrac_isLoaded = true;
		}
		return jpts_chFrac_;
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
	vector<float> &ww_pmet()
	{
		if (not ww_pmet_isLoaded) {
			if (ww_pmet_branch != 0) {
				ww_pmet_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<float>::const_iterator i = ww_pmet_.begin(); i != ww_pmet_.end(); ++i) {
					if (not isfinite(*i)) {
						printf("branch ww_pmet_branch contains a bad float: %f\n", *i);
						exit(1);
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_pmet_branch does not exist!\n");
				exit(1);
			}
			ww_pmet_isLoaded = true;
		}
		return ww_pmet_;
	}
	vector<vector<float> > &hyp_jets_EMFcor()
	{
		if (not hyp_jets_EMFcor_isLoaded) {
			if (hyp_jets_EMFcor_branch != 0) {
				hyp_jets_EMFcor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_EMFcor_.begin(); i != hyp_jets_EMFcor_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_EMFcor_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_EMFcor_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_EMFcor_isLoaded = true;
		}
		return hyp_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_jets_chFrac()
	{
		if (not hyp_jets_chFrac_isLoaded) {
			if (hyp_jets_chFrac_branch != 0) {
				hyp_jets_chFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_chFrac_.begin(); i != hyp_jets_chFrac_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_chFrac_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_chFrac_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_chFrac_isLoaded = true;
		}
		return hyp_jets_chFrac_;
	}
	vector<vector<float> > &hyp_jets_cor()
	{
		if (not hyp_jets_cor_isLoaded) {
			if (hyp_jets_cor_branch != 0) {
				hyp_jets_cor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_cor_.begin(); i != hyp_jets_cor_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_cor_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_cor_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_cor_isLoaded = true;
		}
		return hyp_jets_cor_;
	}
	vector<vector<float> > &hyp_jets_emFrac()
	{
		if (not hyp_jets_emFrac_isLoaded) {
			if (hyp_jets_emFrac_branch != 0) {
				hyp_jets_emFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_emFrac_.begin(); i != hyp_jets_emFrac_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_emFrac_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_emFrac_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_emFrac_isLoaded = true;
		}
		return hyp_jets_emFrac_;
	}
	vector<vector<float> > &hyp_jets_mc_emEnergy()
	{
		if (not hyp_jets_mc_emEnergy_isLoaded) {
			if (hyp_jets_mc_emEnergy_branch != 0) {
				hyp_jets_mc_emEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_mc_emEnergy_.begin(); i != hyp_jets_mc_emEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_mc_emEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_emEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_hadEnergy()
	{
		if (not hyp_jets_mc_hadEnergy_isLoaded) {
			if (hyp_jets_mc_hadEnergy_branch != 0) {
				hyp_jets_mc_hadEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_mc_hadEnergy_.begin(); i != hyp_jets_mc_hadEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_mc_hadEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_hadEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_invEnergy()
	{
		if (not hyp_jets_mc_invEnergy_isLoaded) {
			if (hyp_jets_mc_invEnergy_branch != 0) {
				hyp_jets_mc_invEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_mc_invEnergy_.begin(); i != hyp_jets_mc_invEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_mc_invEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_invEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_otherEnergy()
	{
		if (not hyp_jets_mc_otherEnergy_isLoaded) {
			if (hyp_jets_mc_otherEnergy_branch != 0) {
				hyp_jets_mc_otherEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_mc_otherEnergy_.begin(); i != hyp_jets_mc_otherEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_mc_otherEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_otherEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_jets_pat_bCorrF()
	{
		if (not hyp_jets_pat_bCorrF_isLoaded) {
			if (hyp_jets_pat_bCorrF_branch != 0) {
				hyp_jets_pat_bCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_bCorrF_.begin(); i != hyp_jets_pat_bCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_bCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_bCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_bCorrF_isLoaded = true;
		}
		return hyp_jets_pat_bCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_cCorrF()
	{
		if (not hyp_jets_pat_cCorrF_isLoaded) {
			if (hyp_jets_pat_cCorrF_branch != 0) {
				hyp_jets_pat_cCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_cCorrF_.begin(); i != hyp_jets_pat_cCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_cCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_cCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_cCorrF_isLoaded = true;
		}
		return hyp_jets_pat_cCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_gluCorrF()
	{
		if (not hyp_jets_pat_gluCorrF_isLoaded) {
			if (hyp_jets_pat_gluCorrF_branch != 0) {
				hyp_jets_pat_gluCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_gluCorrF_.begin(); i != hyp_jets_pat_gluCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_gluCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_gluCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_gluCorrF_isLoaded = true;
		}
		return hyp_jets_pat_gluCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_jetCharge()
	{
		if (not hyp_jets_pat_jetCharge_isLoaded) {
			if (hyp_jets_pat_jetCharge_branch != 0) {
				hyp_jets_pat_jetCharge_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_jetCharge_.begin(); i != hyp_jets_pat_jetCharge_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_jetCharge_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_jetCharge_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_jetCharge_isLoaded = true;
		}
		return hyp_jets_pat_jetCharge_;
	}
	vector<vector<float> > &hyp_jets_pat_noCorrF()
	{
		if (not hyp_jets_pat_noCorrF_isLoaded) {
			if (hyp_jets_pat_noCorrF_branch != 0) {
				hyp_jets_pat_noCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_noCorrF_.begin(); i != hyp_jets_pat_noCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_noCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_noCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_noCorrF_isLoaded = true;
		}
		return hyp_jets_pat_noCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_udsCorrF()
	{
		if (not hyp_jets_pat_udsCorrF_isLoaded) {
			if (hyp_jets_pat_udsCorrF_branch != 0) {
				hyp_jets_pat_udsCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_jets_pat_udsCorrF_.begin(); i != hyp_jets_pat_udsCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_jets_pat_udsCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_udsCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_udsCorrF_isLoaded = true;
		}
		return hyp_jets_pat_udsCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_EMFcor()
	{
		if (not hyp_other_jets_EMFcor_isLoaded) {
			if (hyp_other_jets_EMFcor_branch != 0) {
				hyp_other_jets_EMFcor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_EMFcor_.begin(); i != hyp_other_jets_EMFcor_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_EMFcor_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_EMFcor_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_EMFcor_isLoaded = true;
		}
		return hyp_other_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_other_jets_chFrac()
	{
		if (not hyp_other_jets_chFrac_isLoaded) {
			if (hyp_other_jets_chFrac_branch != 0) {
				hyp_other_jets_chFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_chFrac_.begin(); i != hyp_other_jets_chFrac_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_chFrac_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_chFrac_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_chFrac_isLoaded = true;
		}
		return hyp_other_jets_chFrac_;
	}
	vector<vector<float> > &hyp_other_jets_cor()
	{
		if (not hyp_other_jets_cor_isLoaded) {
			if (hyp_other_jets_cor_branch != 0) {
				hyp_other_jets_cor_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_cor_.begin(); i != hyp_other_jets_cor_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_cor_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_cor_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_cor_isLoaded = true;
		}
		return hyp_other_jets_cor_;
	}
	vector<vector<float> > &hyp_other_jets_emFrac()
	{
		if (not hyp_other_jets_emFrac_isLoaded) {
			if (hyp_other_jets_emFrac_branch != 0) {
				hyp_other_jets_emFrac_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_emFrac_.begin(); i != hyp_other_jets_emFrac_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_emFrac_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_emFrac_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_emFrac_isLoaded = true;
		}
		return hyp_other_jets_emFrac_;
	}
	vector<vector<float> > &hyp_other_jets_mc_emEnergy()
	{
		if (not hyp_other_jets_mc_emEnergy_isLoaded) {
			if (hyp_other_jets_mc_emEnergy_branch != 0) {
				hyp_other_jets_mc_emEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_mc_emEnergy_.begin(); i != hyp_other_jets_mc_emEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_mc_emEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_emEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_hadEnergy()
	{
		if (not hyp_other_jets_mc_hadEnergy_isLoaded) {
			if (hyp_other_jets_mc_hadEnergy_branch != 0) {
				hyp_other_jets_mc_hadEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_mc_hadEnergy_.begin(); i != hyp_other_jets_mc_hadEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_mc_hadEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_hadEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_invEnergy()
	{
		if (not hyp_other_jets_mc_invEnergy_isLoaded) {
			if (hyp_other_jets_mc_invEnergy_branch != 0) {
				hyp_other_jets_mc_invEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_mc_invEnergy_.begin(); i != hyp_other_jets_mc_invEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_mc_invEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_invEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_otherEnergy()
	{
		if (not hyp_other_jets_mc_otherEnergy_isLoaded) {
			if (hyp_other_jets_mc_otherEnergy_branch != 0) {
				hyp_other_jets_mc_otherEnergy_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_mc_otherEnergy_.begin(); i != hyp_other_jets_mc_otherEnergy_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_mc_otherEnergy_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_otherEnergy_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_pat_bCorrF()
	{
		if (not hyp_other_jets_pat_bCorrF_isLoaded) {
			if (hyp_other_jets_pat_bCorrF_branch != 0) {
				hyp_other_jets_pat_bCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_bCorrF_.begin(); i != hyp_other_jets_pat_bCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_bCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_bCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_bCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_bCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_cCorrF()
	{
		if (not hyp_other_jets_pat_cCorrF_isLoaded) {
			if (hyp_other_jets_pat_cCorrF_branch != 0) {
				hyp_other_jets_pat_cCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_cCorrF_.begin(); i != hyp_other_jets_pat_cCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_cCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_cCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_cCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_cCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_gluCorrF()
	{
		if (not hyp_other_jets_pat_gluCorrF_isLoaded) {
			if (hyp_other_jets_pat_gluCorrF_branch != 0) {
				hyp_other_jets_pat_gluCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_gluCorrF_.begin(); i != hyp_other_jets_pat_gluCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_gluCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_gluCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_gluCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_gluCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_jetCharge()
	{
		if (not hyp_other_jets_pat_jetCharge_isLoaded) {
			if (hyp_other_jets_pat_jetCharge_branch != 0) {
				hyp_other_jets_pat_jetCharge_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_jetCharge_.begin(); i != hyp_other_jets_pat_jetCharge_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_jetCharge_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_jetCharge_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_jetCharge_isLoaded = true;
		}
		return hyp_other_jets_pat_jetCharge_;
	}
	vector<vector<float> > &hyp_other_jets_pat_noCorrF()
	{
		if (not hyp_other_jets_pat_noCorrF_isLoaded) {
			if (hyp_other_jets_pat_noCorrF_branch != 0) {
				hyp_other_jets_pat_noCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_noCorrF_.begin(); i != hyp_other_jets_pat_noCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_noCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_noCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_noCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_noCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_udsCorrF()
	{
		if (not hyp_other_jets_pat_udsCorrF_isLoaded) {
			if (hyp_other_jets_pat_udsCorrF_branch != 0) {
				hyp_other_jets_pat_udsCorrF_branch->GetEntry(index);
				#ifdef PARANOIA
				for (vector<vector<float> >::const_iterator i = hyp_other_jets_pat_udsCorrF_.begin(); i != hyp_other_jets_pat_udsCorrF_.end(); ++i) {
					for (vector<float>::const_iterator j = i->begin(); j != i->end(); ++j) {
						if (not isfinite(*j)) {
							printf("branch hyp_other_jets_pat_udsCorrF_branch contains a bad float: %f\n", *j);
							exit(1);
						}
					}
				}
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_udsCorrF_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_udsCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_udsCorrF_;
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
	vector<int> &ww_isdyee()
	{
		if (not ww_isdyee_isLoaded) {
			if (ww_isdyee_branch != 0) {
				ww_isdyee_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_isdyee_branch does not exist!\n");
				exit(1);
			}
			ww_isdyee_isLoaded = true;
		}
		return ww_isdyee_;
	}
	vector<int> &ww_isdymm()
	{
		if (not ww_isdymm_isLoaded) {
			if (ww_isdymm_branch != 0) {
				ww_isdymm_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_isdymm_branch does not exist!\n");
				exit(1);
			}
			ww_isdymm_isLoaded = true;
		}
		return ww_isdymm_;
	}
	vector<int> &ww_isdytt()
	{
		if (not ww_isdytt_isLoaded) {
			if (ww_isdytt_branch != 0) {
				ww_isdytt_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_isdytt_branch does not exist!\n");
				exit(1);
			}
			ww_isdytt_isLoaded = true;
		}
		return ww_isdytt_;
	}
	vector<int> &ww_llgoodel()
	{
		if (not ww_llgoodel_isLoaded) {
			if (ww_llgoodel_branch != 0) {
				ww_llgoodel_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_llgoodel_branch does not exist!\n");
				exit(1);
			}
			ww_llgoodel_isLoaded = true;
		}
		return ww_llgoodel_;
	}
	vector<int> &ww_llgoodeliso()
	{
		if (not ww_llgoodeliso_isLoaded) {
			if (ww_llgoodeliso_branch != 0) {
				ww_llgoodeliso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_llgoodeliso_branch does not exist!\n");
				exit(1);
			}
			ww_llgoodeliso_isLoaded = true;
		}
		return ww_llgoodeliso_;
	}
	vector<int> &ww_llgoodmu()
	{
		if (not ww_llgoodmu_isLoaded) {
			if (ww_llgoodmu_branch != 0) {
				ww_llgoodmu_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_llgoodmu_branch does not exist!\n");
				exit(1);
			}
			ww_llgoodmu_isLoaded = true;
		}
		return ww_llgoodmu_;
	}
	vector<int> &ww_llgoodmuiso()
	{
		if (not ww_llgoodmuiso_isLoaded) {
			if (ww_llgoodmuiso_branch != 0) {
				ww_llgoodmuiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_llgoodmuiso_branch does not exist!\n");
				exit(1);
			}
			ww_llgoodmuiso_isLoaded = true;
		}
		return ww_llgoodmuiso_;
	}
	vector<int> &ww_ltgoodel()
	{
		if (not ww_ltgoodel_isLoaded) {
			if (ww_ltgoodel_branch != 0) {
				ww_ltgoodel_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_ltgoodel_branch does not exist!\n");
				exit(1);
			}
			ww_ltgoodel_isLoaded = true;
		}
		return ww_ltgoodel_;
	}
	vector<int> &ww_ltgoodeliso()
	{
		if (not ww_ltgoodeliso_isLoaded) {
			if (ww_ltgoodeliso_branch != 0) {
				ww_ltgoodeliso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_ltgoodeliso_branch does not exist!\n");
				exit(1);
			}
			ww_ltgoodeliso_isLoaded = true;
		}
		return ww_ltgoodeliso_;
	}
	vector<int> &ww_ltgoodmu()
	{
		if (not ww_ltgoodmu_isLoaded) {
			if (ww_ltgoodmu_branch != 0) {
				ww_ltgoodmu_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_ltgoodmu_branch does not exist!\n");
				exit(1);
			}
			ww_ltgoodmu_isLoaded = true;
		}
		return ww_ltgoodmu_;
	}
	vector<int> &ww_ltgoodmuiso()
	{
		if (not ww_ltgoodmuiso_isLoaded) {
			if (ww_ltgoodmuiso_branch != 0) {
				ww_ltgoodmuiso_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_ltgoodmuiso_branch does not exist!\n");
				exit(1);
			}
			ww_ltgoodmuiso_isLoaded = true;
		}
		return ww_ltgoodmuiso_;
	}
	vector<int> &ww_oppsign()
	{
		if (not ww_oppsign_isLoaded) {
			if (ww_oppsign_branch != 0) {
				ww_oppsign_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_oppsign_branch does not exist!\n");
				exit(1);
			}
			ww_oppsign_isLoaded = true;
		}
		return ww_oppsign_;
	}
	vector<int> &ww_pass2met()
	{
		if (not ww_pass2met_isLoaded) {
			if (ww_pass2met_branch != 0) {
				ww_pass2met_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_pass2met_branch does not exist!\n");
				exit(1);
			}
			ww_pass2met_isLoaded = true;
		}
		return ww_pass2met_;
	}
	vector<int> &ww_pass4met()
	{
		if (not ww_pass4met_isLoaded) {
			if (ww_pass4met_branch != 0) {
				ww_pass4met_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_pass4met_branch does not exist!\n");
				exit(1);
			}
			ww_pass4met_isLoaded = true;
		}
		return ww_pass4met_;
	}
	vector<int> &ww_passFebsel_no_jetveto()
	{
		if (not ww_passFebsel_no_jetveto_isLoaded) {
			if (ww_passFebsel_no_jetveto_branch != 0) {
				ww_passFebsel_no_jetveto_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_passFebsel_no_jetveto_branch does not exist!\n");
				exit(1);
			}
			ww_passFebsel_no_jetveto_isLoaded = true;
		}
		return ww_passFebsel_no_jetveto_;
	}
	vector<int> &ww_passFebsel_with_jetveto()
	{
		if (not ww_passFebsel_with_jetveto_isLoaded) {
			if (ww_passFebsel_with_jetveto_branch != 0) {
				ww_passFebsel_with_jetveto_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_passFebsel_with_jetveto_branch does not exist!\n");
				exit(1);
			}
			ww_passFebsel_with_jetveto_isLoaded = true;
		}
		return ww_passFebsel_with_jetveto_;
	}
	vector<int> &ww_passaddzveto()
	{
		if (not ww_passaddzveto_isLoaded) {
			if (ww_passaddzveto_branch != 0) {
				ww_passaddzveto_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_passaddzveto_branch does not exist!\n");
				exit(1);
			}
			ww_passaddzveto_isLoaded = true;
		}
		return ww_passaddzveto_;
	}
	vector<int> &ww_passzveto()
	{
		if (not ww_passzveto_isLoaded) {
			if (ww_passzveto_branch != 0) {
				ww_passzveto_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch ww_passzveto_branch does not exist!\n");
				exit(1);
			}
			ww_passzveto_isLoaded = true;
		}
		return ww_passzveto_;
	}
	vector<vector<int> > &hyp_jets_mc_id()
	{
		if (not hyp_jets_mc_id_isLoaded) {
			if (hyp_jets_mc_id_branch != 0) {
				hyp_jets_mc_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_mc_id_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_mc_id_isLoaded = true;
		}
		return hyp_jets_mc_id_;
	}
	vector<vector<int> > &hyp_jets_pat_genPartonMother_id()
	{
		if (not hyp_jets_pat_genPartonMother_id_isLoaded) {
			if (hyp_jets_pat_genPartonMother_id_branch != 0) {
				hyp_jets_pat_genPartonMother_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_genPartonMother_id_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_genPartonMother_id_isLoaded = true;
		}
		return hyp_jets_pat_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_jets_pat_genParton_id()
	{
		if (not hyp_jets_pat_genParton_id_isLoaded) {
			if (hyp_jets_pat_genParton_id_branch != 0) {
				hyp_jets_pat_genParton_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_genParton_id_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_genParton_id_isLoaded = true;
		}
		return hyp_jets_pat_genParton_id_;
	}
	vector<vector<int> > &hyp_jets_pat_partonFlavour()
	{
		if (not hyp_jets_pat_partonFlavour_isLoaded) {
			if (hyp_jets_pat_partonFlavour_branch != 0) {
				hyp_jets_pat_partonFlavour_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_jets_pat_partonFlavour_branch does not exist!\n");
				exit(1);
			}
			hyp_jets_pat_partonFlavour_isLoaded = true;
		}
		return hyp_jets_pat_partonFlavour_;
	}
	vector<vector<int> > &hyp_other_jets_mc_id()
	{
		if (not hyp_other_jets_mc_id_isLoaded) {
			if (hyp_other_jets_mc_id_branch != 0) {
				hyp_other_jets_mc_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_mc_id_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_mc_id_isLoaded = true;
		}
		return hyp_other_jets_mc_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_genPartonMother_id()
	{
		if (not hyp_other_jets_pat_genPartonMother_id_isLoaded) {
			if (hyp_other_jets_pat_genPartonMother_id_branch != 0) {
				hyp_other_jets_pat_genPartonMother_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_genPartonMother_id_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_genPartonMother_id_isLoaded = true;
		}
		return hyp_other_jets_pat_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_genParton_id()
	{
		if (not hyp_other_jets_pat_genParton_id_isLoaded) {
			if (hyp_other_jets_pat_genParton_id_branch != 0) {
				hyp_other_jets_pat_genParton_id_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_genParton_id_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_genParton_id_isLoaded = true;
		}
		return hyp_other_jets_pat_genParton_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_partonFlavour()
	{
		if (not hyp_other_jets_pat_partonFlavour_isLoaded) {
			if (hyp_other_jets_pat_partonFlavour_branch != 0) {
				hyp_other_jets_pat_partonFlavour_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_other_jets_pat_partonFlavour_branch does not exist!\n");
				exit(1);
			}
			hyp_other_jets_pat_partonFlavour_isLoaded = true;
		}
		return hyp_other_jets_pat_partonFlavour_;
	}
	vector<vector<int> > &hyp_quadlep_jets_index()
	{
		if (not hyp_quadlep_jets_index_isLoaded) {
			if (hyp_quadlep_jets_index_branch != 0) {
				hyp_quadlep_jets_index_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_quadlep_jets_index_branch does not exist!\n");
				exit(1);
			}
			hyp_quadlep_jets_index_isLoaded = true;
		}
		return hyp_quadlep_jets_index_;
	}
	vector<vector<int> > &hyp_trilep_jets_index()
	{
		if (not hyp_trilep_jets_index_isLoaded) {
			if (hyp_trilep_jets_index_branch != 0) {
				hyp_trilep_jets_index_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch hyp_trilep_jets_index_branch does not exist!\n");
				exit(1);
			}
			hyp_trilep_jets_index_isLoaded = true;
		}
		return hyp_trilep_jets_index_;
	}
	unsigned int &evt_nalltrkjets()
	{
		if (not evt_nalltrkjets_isLoaded) {
			if (evt_nalltrkjets_branch != 0) {
				evt_nalltrkjets_branch->GetEntry(index);
				#ifdef PARANOIA
				#endif // #ifdef PARANOIA
			} else { 
				printf("branch evt_nalltrkjets_branch does not exist!\n");
				exit(1);
			}
			evt_nalltrkjets_isLoaded = true;
		}
		return evt_nalltrkjets_;
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
#endif
