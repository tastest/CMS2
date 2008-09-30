// -*- C++ -*-
#ifndef CMS2_H
#define CMS2_H
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"

#include "TH1F.h"

#include <vector> 
using namespace std; 
class CMS2 { 
private: 
	 TH1F *samplehisto;
protected: 
	unsigned int index;
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 	l1met_p4_;
	TBranch *l1met_p4_branch;
	bool l1met_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_tq_genMotherP4_;
	TBranch *els_tq_genMotherP4_branch;
	bool els_tq_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_tq_genP4_;
	TBranch *els_tq_genP4_branch;
	bool els_tq_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genJet_p4_;
	TBranch *jets_tq_genJet_p4_branch;
	bool jets_tq_genJet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genPartonMother_p4_;
	TBranch *jets_tq_genPartonMother_p4_branch;
	bool jets_tq_genPartonMother_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genParton_p4_;
	TBranch *jets_tq_genParton_p4_branch;
	bool jets_tq_genParton_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_jet_p4_;
	TBranch *jets_tq_jet_p4_branch;
	bool jets_tq_jet_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_tq_genMotherP4_;
	TBranch *mus_tq_genMotherP4_branch;
	bool mus_tq_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_tq_genP4_;
	TBranch *mus_tq_genP4_branch;
	bool mus_tq_genP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trks_trk_p4_;
	TBranch *trks_trk_p4_branch;
	bool trks_trk_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_gp_p4_;
	TBranch *hyp_jets_mc_gp_p4_branch;
	bool hyp_jets_mc_gp_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_p4_;
	TBranch *hyp_jets_mc_p4_branch;
	bool hyp_jets_mc_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_p4_;
	TBranch *hyp_jets_p4_branch;
	bool hyp_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_p4_;
	TBranch *hyp_jets_tq_p4_branch;
	bool hyp_jets_tq_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_genPartonMother_p4_;
	TBranch *hyp_jets_tq_genPartonMother_p4_branch;
	bool hyp_jets_tq_genPartonMother_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_genParton_p4_;
	TBranch *hyp_jets_tq_genParton_p4_branch;
	bool hyp_jets_tq_genParton_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_jet_p4_;
	TBranch *hyp_jets_tq_jet_p4_branch;
	bool hyp_jets_tq_jet_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_gp_p4_;
	TBranch *hyp_other_jets_mc_gp_p4_branch;
	bool hyp_other_jets_mc_gp_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_p4_;
	TBranch *hyp_other_jets_mc_p4_branch;
	bool hyp_other_jets_mc_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_p4_;
	TBranch *hyp_other_jets_p4_branch;
	bool hyp_other_jets_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genJet_p4_;
	TBranch *hyp_other_jets_tq_genJet_p4_branch;
	bool hyp_other_jets_tq_genJet_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genPartonMother_p4_;
	TBranch *hyp_other_jets_tq_genPartonMother_p4_branch;
	bool hyp_other_jets_tq_genPartonMother_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genParton_p4_;
	TBranch *hyp_other_jets_tq_genParton_p4_branch;
	bool hyp_other_jets_tq_genParton_p4_isLoaded;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_jet_p4_;
	TBranch *hyp_other_jets_tq_jet_p4_branch;
	bool hyp_other_jets_tq_jet_p4_isLoaded;
	float	evt_CSA07FilterEff_;
	TBranch *evt_CSA07FilterEff_branch;
	bool evt_CSA07FilterEff_isLoaded;
	float	evt_CSA07Pthat_;
	TBranch *evt_CSA07Pthat_branch;
	bool evt_CSA07Pthat_isLoaded;
	float	evt_CSA07Weight_;
	TBranch *evt_CSA07Weight_branch;
	bool evt_CSA07Weight_isLoaded;
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
	float	evt_scale1fb_;
	TBranch *evt_scale1fb_branch;
	bool evt_scale1fb_isLoaded;
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
	float	evt_metPhi_;
	TBranch *evt_metPhi_branch;
	bool evt_metPhi_isLoaded;
	float	evt_metSig_;
	TBranch *evt_metSig_branch;
	bool evt_metSig_isLoaded;
	float	evt_met_jetcorr_;
	TBranch *evt_met_jetcorr_branch;
	bool evt_met_jetcorr_isLoaded;
	float	metphi_jetcorr_;
	TBranch *metphi_jetcorr_branch;
	bool metphi_jetcorr_isLoaded;
	vector<float>	els_musdr_;
	TBranch *els_musdr_branch;
	bool els_musdr_isLoaded;
	vector<float>	els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float>	els_ESc_;
	TBranch *els_ESc_branch;
	bool els_ESc_isLoaded;
	vector<float>	els_ESc_raw_;
	TBranch *els_ESc_raw_branch;
	bool els_ESc_raw_isLoaded;
	vector<float>	els_ESeed_;
	TBranch *els_ESeed_branch;
	bool els_ESeed_isLoaded;
	vector<float>	els_chi2_;
	TBranch *els_chi2_branch;
	bool els_chi2_isLoaded;
	vector<float>	els_d0_;
	TBranch *els_d0_branch;
	bool els_d0_isLoaded;
	vector<float>	els_d0Err_;
	TBranch *els_d0Err_branch;
	bool els_d0Err_isLoaded;
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
	vector<float>	els_e3x3_;
	TBranch *els_e3x3_branch;
	bool els_e3x3_isLoaded;
	vector<float>	els_e5x5_;
	TBranch *els_e5x5_branch;
	bool els_e5x5_isLoaded;
	vector<float>	els_eOverPIn_;
	TBranch *els_eOverPIn_branch;
	bool els_eOverPIn_isLoaded;
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
	vector<float>	hyp_ll_chi2_;
	TBranch *hyp_ll_chi2_branch;
	bool hyp_ll_chi2_isLoaded;
	vector<float>	hyp_ll_d0_;
	TBranch *hyp_ll_d0_branch;
	bool hyp_ll_d0_isLoaded;
	vector<float>	hyp_ll_d0Err_;
	TBranch *hyp_ll_d0Err_branch;
	bool hyp_ll_d0Err_isLoaded;
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
	vector<float>	hyp_lt_chi2_;
	TBranch *hyp_lt_chi2_branch;
	bool hyp_lt_chi2_isLoaded;
	vector<float>	hyp_lt_d0_;
	TBranch *hyp_lt_d0_branch;
	bool hyp_lt_d0_isLoaded;
	vector<float>	hyp_lt_d0Err_;
	TBranch *hyp_lt_d0Err_branch;
	bool hyp_lt_d0Err_isLoaded;
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
	vector<float>	mus_eledr_;
	TBranch *mus_eledr_branch;
	bool mus_eledr_isLoaded;
	vector<float>	mus_jetdr_;
	TBranch *mus_jetdr_branch;
	bool mus_jetdr_isLoaded;
	vector<float>	mus_trkdr_;
	TBranch *mus_trkdr_branch;
	bool mus_trkdr_isLoaded;
	vector<float>	mus_chi2_;
	TBranch *mus_chi2_branch;
	bool mus_chi2_isLoaded;
	vector<float>	mus_d0_;
	TBranch *mus_d0_branch;
	bool mus_d0_isLoaded;
	vector<float>	mus_d0Err_;
	TBranch *mus_d0Err_branch;
	bool mus_d0Err_isLoaded;
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
	vector<float>	els_tq_LRComb_;
	TBranch *els_tq_LRComb_branch;
	bool els_tq_LRComb_isLoaded;
	vector<float>	els_tq_caloIso_;
	TBranch *els_tq_caloIso_branch;
	bool els_tq_caloIso_isLoaded;
	vector<float>	els_tq_egammaEcalIso_;
	TBranch *els_tq_egammaEcalIso_branch;
	bool els_tq_egammaEcalIso_isLoaded;
	vector<float>	els_tq_egammaHcalIso_;
	TBranch *els_tq_egammaHcalIso_branch;
	bool els_tq_egammaHcalIso_isLoaded;
	vector<float>	els_tq_egammaTkIso_;
	TBranch *els_tq_egammaTkIso_branch;
	bool els_tq_egammaTkIso_isLoaded;
	vector<float>	els_tq_electronIDRobust_;
	TBranch *els_tq_electronIDRobust_branch;
	bool els_tq_electronIDRobust_isLoaded;
	vector<float>	els_tq_leptonID_;
	TBranch *els_tq_leptonID_branch;
	bool els_tq_leptonID_isLoaded;
	vector<float>	els_tq_trackIso_;
	TBranch *els_tq_trackIso_branch;
	bool els_tq_trackIso_isLoaded;
	vector<float>	jets_tq_bCorrF_;
	TBranch *jets_tq_bCorrF_branch;
	bool jets_tq_bCorrF_isLoaded;
	vector<float>	jets_tq_cCorrF_;
	TBranch *jets_tq_cCorrF_branch;
	bool jets_tq_cCorrF_isLoaded;
	vector<float>	jets_tq_gluCorrF_;
	TBranch *jets_tq_gluCorrF_branch;
	bool jets_tq_gluCorrF_isLoaded;
	vector<float>	jets_tq_jetCharge_;
	TBranch *jets_tq_jetCharge_branch;
	bool jets_tq_jetCharge_isLoaded;
	vector<float>	jets_tq_noCorrF_;
	TBranch *jets_tq_noCorrF_branch;
	bool jets_tq_noCorrF_isLoaded;
	vector<float>	jets_tq_udsCorrF_;
	TBranch *jets_tq_udsCorrF_branch;
	bool jets_tq_udsCorrF_isLoaded;
	vector<float>	mus_tq_caloIso_;
	TBranch *mus_tq_caloIso_branch;
	bool mus_tq_caloIso_isLoaded;
	vector<float>	mus_tq_leptonID_;
	TBranch *mus_tq_leptonID_branch;
	bool mus_tq_leptonID_isLoaded;
	vector<float>	mus_tq_lrComb_;
	TBranch *mus_tq_lrComb_branch;
	bool mus_tq_lrComb_isLoaded;
	vector<float>	mus_tq_trackIso_;
	TBranch *mus_tq_trackIso_branch;
	bool mus_tq_trackIso_isLoaded;
	vector<float>	trks_chi2_;
	TBranch *trks_chi2_branch;
	bool trks_chi2_isLoaded;
	vector<float>	trks_d0_;
	TBranch *trks_d0_branch;
	bool trks_d0_isLoaded;
	vector<float>	trks_d0Err_;
	TBranch *trks_d0Err_branch;
	bool trks_d0Err_isLoaded;
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
	vector<float>	trks_vertexphi_;
	TBranch *trks_vertexphi_branch;
	bool trks_vertexphi_isLoaded;
	vector<float>	trks_z0_;
	TBranch *trks_z0_branch;
	bool trks_z0_isLoaded;
	vector<float>	trks_z0Err_;
	TBranch *trks_z0Err_branch;
	bool trks_z0Err_isLoaded;
	vector<float>	trk_elsdr_;
	TBranch *trk_elsdr_branch;
	bool trk_elsdr_isLoaded;
	vector<float>	trk_musdr_;
	TBranch *trk_musdr_branch;
	bool trk_musdr_isLoaded;
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
	vector<vector<float> >	hyp_jets_tq_bCorrF_;
	TBranch *hyp_jets_tq_bCorrF_branch;
	bool hyp_jets_tq_bCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_tq_cCorrF_;
	TBranch *hyp_jets_tq_cCorrF_branch;
	bool hyp_jets_tq_cCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_tq_gluCorrF_;
	TBranch *hyp_jets_tq_gluCorrF_branch;
	bool hyp_jets_tq_gluCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_tq_jetCharge_;
	TBranch *hyp_jets_tq_jetCharge_branch;
	bool hyp_jets_tq_jetCharge_isLoaded;
	vector<vector<float> >	hyp_jets_tq_noCorrF_;
	TBranch *hyp_jets_tq_noCorrF_branch;
	bool hyp_jets_tq_noCorrF_isLoaded;
	vector<vector<float> >	hyp_jets_tq_udsCorrF_;
	TBranch *hyp_jets_tq_udsCorrF_branch;
	bool hyp_jets_tq_udsCorrF_isLoaded;
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
	vector<vector<float> >	hyp_other_jets_tq_bCorrF_;
	TBranch *hyp_other_jets_tq_bCorrF_branch;
	bool hyp_other_jets_tq_bCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_tq_cCorrF_;
	TBranch *hyp_other_jets_tq_cCorrF_branch;
	bool hyp_other_jets_tq_cCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_tq_gluCorrF_;
	TBranch *hyp_other_jets_tq_gluCorrF_branch;
	bool hyp_other_jets_tq_gluCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_tq_jetCharge_;
	TBranch *hyp_other_jets_tq_jetCharge_branch;
	bool hyp_other_jets_tq_jetCharge_isLoaded;
	vector<vector<float> >	hyp_other_jets_tq_noCorrF_;
	TBranch *hyp_other_jets_tq_noCorrF_branch;
	bool hyp_other_jets_tq_noCorrF_isLoaded;
	vector<vector<float> >	hyp_other_jets_tq_udsCorrF_;
	TBranch *hyp_other_jets_tq_udsCorrF_branch;
	bool hyp_other_jets_tq_udsCorrF_isLoaded;
	int	evt_CSA07Process_;
	TBranch *evt_CSA07Process_branch;
	bool evt_CSA07Process_isLoaded;
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
	vector<int>	els_closestMuon_;
	TBranch *els_closestMuon_branch;
	bool els_closestMuon_isLoaded;
	vector<int>	els_trkidx_;
	TBranch *els_trkidx_branch;
	bool els_trkidx_isLoaded;
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
	vector<int>	mus_trkrefkey_;
	TBranch *mus_trkrefkey_branch;
	bool mus_trkrefkey_isLoaded;
	vector<int>	mus_validHits_;
	TBranch *mus_validHits_branch;
	bool mus_validHits_isLoaded;
	vector<int>	els_tq_egammaTkNumIso_;
	TBranch *els_tq_egammaTkNumIso_branch;
	bool els_tq_egammaTkNumIso_isLoaded;
	vector<int>	els_tq_genID_;
	TBranch *els_tq_genID_branch;
	bool els_tq_genID_isLoaded;
	vector<int>	els_tq_genMotherID_;
	TBranch *els_tq_genMotherID_branch;
	bool els_tq_genMotherID_isLoaded;
	vector<int>	jets_tq_genPartonMother_id_;
	TBranch *jets_tq_genPartonMother_id_branch;
	bool jets_tq_genPartonMother_id_isLoaded;
	vector<int>	jets_tq_genParton_id_;
	TBranch *jets_tq_genParton_id_branch;
	bool jets_tq_genParton_id_isLoaded;
	vector<int>	jets_tq_partonFlavour_;
	TBranch *jets_tq_partonFlavour_branch;
	bool jets_tq_partonFlavour_isLoaded;
	vector<int>	mus_tq_genID_;
	TBranch *mus_tq_genID_branch;
	bool mus_tq_genID_isLoaded;
	vector<int>	mus_tq_genMotherID_;
	TBranch *mus_tq_genMotherID_branch;
	bool mus_tq_genMotherID_isLoaded;
	vector<int>	trks_charge_;
	TBranch *trks_charge_branch;
	bool trks_charge_isLoaded;
	vector<int>	trks_lostHits_;
	TBranch *trks_lostHits_branch;
	bool trks_lostHits_isLoaded;
	vector<int>	trks_validHits_;
	TBranch *trks_validHits_branch;
	bool trks_validHits_isLoaded;
	vector<int>	trk_elsidx_;
	TBranch *trk_elsidx_branch;
	bool trk_elsidx_isLoaded;
	vector<int>	trk_musidx_;
	TBranch *trk_musidx_branch;
	bool trk_musidx_isLoaded;
	vector<vector<int> >	hyp_jets_mc_id_;
	TBranch *hyp_jets_mc_id_branch;
	bool hyp_jets_mc_id_isLoaded;
	vector<vector<int> >	hyp_jets_tq_genPartonMother_id_;
	TBranch *hyp_jets_tq_genPartonMother_id_branch;
	bool hyp_jets_tq_genPartonMother_id_isLoaded;
	vector<vector<int> >	hyp_jets_tq_genParton_id_;
	TBranch *hyp_jets_tq_genParton_id_branch;
	bool hyp_jets_tq_genParton_id_isLoaded;
	vector<vector<int> >	hyp_jets_tq_partonFlavour_;
	TBranch *hyp_jets_tq_partonFlavour_branch;
	bool hyp_jets_tq_partonFlavour_isLoaded;
	vector<vector<int> >	hyp_other_jets_mc_id_;
	TBranch *hyp_other_jets_mc_id_branch;
	bool hyp_other_jets_mc_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_tq_genPartonMother_id_;
	TBranch *hyp_other_jets_tq_genPartonMother_id_branch;
	bool hyp_other_jets_tq_genPartonMother_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_tq_genParton_id_;
	TBranch *hyp_other_jets_tq_genParton_id_branch;
	bool hyp_other_jets_tq_genParton_id_isLoaded;
	vector<vector<int> >	hyp_other_jets_tq_partonFlavour_;
	TBranch *hyp_other_jets_tq_partonFlavour_branch;
	bool hyp_other_jets_tq_partonFlavour_isLoaded;
	vector<vector<int> >	hyp_quadlep_jets_index_;
	TBranch *hyp_quadlep_jets_index_branch;
	bool hyp_quadlep_jets_index_isLoaded;
	vector<vector<int> >	hyp_trilep_jets_index_;
	TBranch *hyp_trilep_jets_index_branch;
	bool hyp_trilep_jets_index_isLoaded;
	unsigned int	evt_nels_;
	TBranch *evt_nels_branch;
	bool evt_nels_isLoaded;
	unsigned int	evt_event_;
	TBranch *evt_event_branch;
	bool evt_event_isLoaded;
	unsigned int	evt_run_;
	TBranch *evt_run_branch;
	bool evt_run_isLoaded;
	unsigned int	evt_njets_;
	TBranch *evt_njets_branch;
	bool evt_njets_isLoaded;
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
public: 
int ScanChain( TChain* chain, int nEvents=-1);
void Init(TTree *tree) {
	l1met_p4_branch = 0;
	if (tree->GetAlias("l1met_p4") != 0) {
		l1met_p4_branch = tree->GetBranch(tree->GetAlias("l1met_p4"));
		l1met_p4_branch->SetAddress(&l1met_p4_);
	}
	if(l1met_p4_branch == 0 ) {
	cout << "Branch l1met_p4 does not exist." << endl;
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
	els_tq_genMotherP4_branch = 0;
	if (tree->GetAlias("els_tq_genMotherP4") != 0) {
		els_tq_genMotherP4_branch = tree->GetBranch(tree->GetAlias("els_tq_genMotherP4"));
		els_tq_genMotherP4_branch->SetAddress(&els_tq_genMotherP4_);
	}
	if(els_tq_genMotherP4_branch == 0 ) {
	cout << "Branch els_tq_genMotherP4 does not exist." << endl;
	}
	els_tq_genP4_branch = 0;
	if (tree->GetAlias("els_tq_genP4") != 0) {
		els_tq_genP4_branch = tree->GetBranch(tree->GetAlias("els_tq_genP4"));
		els_tq_genP4_branch->SetAddress(&els_tq_genP4_);
	}
	if(els_tq_genP4_branch == 0 ) {
	cout << "Branch els_tq_genP4 does not exist." << endl;
	}
	jets_tq_genJet_p4_branch = 0;
	if (tree->GetAlias("jets_tq_genJet_p4") != 0) {
		jets_tq_genJet_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genJet_p4"));
		jets_tq_genJet_p4_branch->SetAddress(&jets_tq_genJet_p4_);
	}
	if(jets_tq_genJet_p4_branch == 0 ) {
	cout << "Branch jets_tq_genJet_p4 does not exist." << endl;
	}
	jets_tq_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("jets_tq_genPartonMother_p4") != 0) {
		jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genPartonMother_p4"));
		jets_tq_genPartonMother_p4_branch->SetAddress(&jets_tq_genPartonMother_p4_);
	}
	if(jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	jets_tq_genParton_p4_branch = 0;
	if (tree->GetAlias("jets_tq_genParton_p4") != 0) {
		jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genParton_p4"));
		jets_tq_genParton_p4_branch->SetAddress(&jets_tq_genParton_p4_);
	}
	if(jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch jets_tq_genParton_p4 does not exist." << endl;
	}
	jets_tq_jet_p4_branch = 0;
	if (tree->GetAlias("jets_tq_jet_p4") != 0) {
		jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_jet_p4"));
		jets_tq_jet_p4_branch->SetAddress(&jets_tq_jet_p4_);
	}
	if(jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch jets_tq_jet_p4 does not exist." << endl;
	}
	mus_tq_genMotherP4_branch = 0;
	if (tree->GetAlias("mus_tq_genMotherP4") != 0) {
		mus_tq_genMotherP4_branch = tree->GetBranch(tree->GetAlias("mus_tq_genMotherP4"));
		mus_tq_genMotherP4_branch->SetAddress(&mus_tq_genMotherP4_);
	}
	if(mus_tq_genMotherP4_branch == 0 ) {
	cout << "Branch mus_tq_genMotherP4 does not exist." << endl;
	}
	mus_tq_genP4_branch = 0;
	if (tree->GetAlias("mus_tq_genP4") != 0) {
		mus_tq_genP4_branch = tree->GetBranch(tree->GetAlias("mus_tq_genP4"));
		mus_tq_genP4_branch->SetAddress(&mus_tq_genP4_);
	}
	if(mus_tq_genP4_branch == 0 ) {
	cout << "Branch mus_tq_genP4 does not exist." << endl;
	}
	trks_trk_p4_branch = 0;
	if (tree->GetAlias("trks_trk_p4") != 0) {
		trks_trk_p4_branch = tree->GetBranch(tree->GetAlias("trks_trk_p4"));
		trks_trk_p4_branch->SetAddress(&trks_trk_p4_);
	}
	if(trks_trk_p4_branch == 0 ) {
	cout << "Branch trks_trk_p4 does not exist." << endl;
	}
  tree->SetMakeClass(1);
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
	hyp_jets_tq_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_p4") != 0) {
		hyp_jets_tq_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_p4"));
		hyp_jets_tq_p4_branch->SetAddress(&hyp_jets_tq_p4_);
	}
	if(hyp_jets_tq_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_p4 does not exist." << endl;
	}
	hyp_jets_tq_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genPartonMother_p4") != 0) {
		hyp_jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genPartonMother_p4"));
		hyp_jets_tq_genPartonMother_p4_branch->SetAddress(&hyp_jets_tq_genPartonMother_p4_);
	}
	if(hyp_jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	hyp_jets_tq_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genParton_p4") != 0) {
		hyp_jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genParton_p4"));
		hyp_jets_tq_genParton_p4_branch->SetAddress(&hyp_jets_tq_genParton_p4_);
	}
	if(hyp_jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genParton_p4 does not exist." << endl;
	}
	hyp_jets_tq_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_jet_p4") != 0) {
		hyp_jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_jet_p4"));
		hyp_jets_tq_jet_p4_branch->SetAddress(&hyp_jets_tq_jet_p4_);
	}
	if(hyp_jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_jet_p4 does not exist." << endl;
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
	hyp_other_jets_tq_genJet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genJet_p4") != 0) {
		hyp_other_jets_tq_genJet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genJet_p4"));
		hyp_other_jets_tq_genJet_p4_branch->SetAddress(&hyp_other_jets_tq_genJet_p4_);
	}
	if(hyp_other_jets_tq_genJet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genJet_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genPartonMother_p4") != 0) {
		hyp_other_jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genPartonMother_p4"));
		hyp_other_jets_tq_genPartonMother_p4_branch->SetAddress(&hyp_other_jets_tq_genPartonMother_p4_);
	}
	if(hyp_other_jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genParton_p4") != 0) {
		hyp_other_jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genParton_p4"));
		hyp_other_jets_tq_genParton_p4_branch->SetAddress(&hyp_other_jets_tq_genParton_p4_);
	}
	if(hyp_other_jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genParton_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_jet_p4") != 0) {
		hyp_other_jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_jet_p4"));
		hyp_other_jets_tq_jet_p4_branch->SetAddress(&hyp_other_jets_tq_jet_p4_);
	}
	if(hyp_other_jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_jet_p4 does not exist." << endl;
	}
	evt_CSA07FilterEff_branch = 0;
	if (tree->GetAlias("evt_CSA07FilterEff") != 0) {
		evt_CSA07FilterEff_branch = tree->GetBranch(tree->GetAlias("evt_CSA07FilterEff"));
		evt_CSA07FilterEff_branch->SetAddress(&evt_CSA07FilterEff_);
	}
	if(evt_CSA07FilterEff_branch == 0 ) {
	cout << "Branch evt_CSA07FilterEff does not exist." << endl;
	}
	evt_CSA07Pthat_branch = 0;
	if (tree->GetAlias("evt_CSA07Pthat") != 0) {
		evt_CSA07Pthat_branch = tree->GetBranch(tree->GetAlias("evt_CSA07Pthat"));
		evt_CSA07Pthat_branch->SetAddress(&evt_CSA07Pthat_);
	}
	if(evt_CSA07Pthat_branch == 0 ) {
	cout << "Branch evt_CSA07Pthat does not exist." << endl;
	}
	evt_CSA07Weight_branch = 0;
	if (tree->GetAlias("evt_CSA07Weight") != 0) {
		evt_CSA07Weight_branch = tree->GetBranch(tree->GetAlias("evt_CSA07Weight"));
		evt_CSA07Weight_branch->SetAddress(&evt_CSA07Weight_);
	}
	if(evt_CSA07Weight_branch == 0 ) {
	cout << "Branch evt_CSA07Weight does not exist." << endl;
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
	evt_scale1fb_branch = 0;
	if (tree->GetAlias("evt_scale1fb") != 0) {
		evt_scale1fb_branch = tree->GetBranch(tree->GetAlias("evt_scale1fb"));
		evt_scale1fb_branch->SetAddress(&evt_scale1fb_);
	}
	if(evt_scale1fb_branch == 0 ) {
	cout << "Branch evt_scale1fb does not exist." << endl;
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
	evt_met_jetcorr_branch = 0;
	if (tree->GetAlias("evt_met_jetcorr") != 0) {
		evt_met_jetcorr_branch = tree->GetBranch(tree->GetAlias("evt_met_jetcorr"));
		evt_met_jetcorr_branch->SetAddress(&evt_met_jetcorr_);
	}
	if(evt_met_jetcorr_branch == 0 ) {
	cout << "Branch evt_met_jetcorr does not exist." << endl;
	}
	metphi_jetcorr_branch = 0;
	if (tree->GetAlias("metphi_jetcorr") != 0) {
		metphi_jetcorr_branch = tree->GetBranch(tree->GetAlias("metphi_jetcorr"));
		metphi_jetcorr_branch->SetAddress(&metphi_jetcorr_);
	}
	if(metphi_jetcorr_branch == 0 ) {
	cout << "Branch metphi_jetcorr does not exist." << endl;
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
	els_ESc_branch = 0;
	if (tree->GetAlias("els_ESc") != 0) {
		els_ESc_branch = tree->GetBranch(tree->GetAlias("els_ESc"));
		els_ESc_branch->SetAddress(&els_ESc_);
	}
	if(els_ESc_branch == 0 ) {
	cout << "Branch els_ESc does not exist." << endl;
	}
	els_ESc_raw_branch = 0;
	if (tree->GetAlias("els_ESc_raw") != 0) {
		els_ESc_raw_branch = tree->GetBranch(tree->GetAlias("els_ESc_raw"));
		els_ESc_raw_branch->SetAddress(&els_ESc_raw_);
	}
	if(els_ESc_raw_branch == 0 ) {
	cout << "Branch els_ESc_raw does not exist." << endl;
	}
	els_ESeed_branch = 0;
	if (tree->GetAlias("els_ESeed") != 0) {
		els_ESeed_branch = tree->GetBranch(tree->GetAlias("els_ESeed"));
		els_ESeed_branch->SetAddress(&els_ESeed_);
	}
	if(els_ESeed_branch == 0 ) {
	cout << "Branch els_ESeed does not exist." << endl;
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
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		els_eOverPIn_branch->SetAddress(&els_eOverPIn_);
	}
	if(els_eOverPIn_branch == 0 ) {
	cout << "Branch els_eOverPIn does not exist." << endl;
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
	els_tq_LRComb_branch = 0;
	if (tree->GetAlias("els_tq_LRComb") != 0) {
		els_tq_LRComb_branch = tree->GetBranch(tree->GetAlias("els_tq_LRComb"));
		els_tq_LRComb_branch->SetAddress(&els_tq_LRComb_);
	}
	if(els_tq_LRComb_branch == 0 ) {
	cout << "Branch els_tq_LRComb does not exist." << endl;
	}
	els_tq_caloIso_branch = 0;
	if (tree->GetAlias("els_tq_caloIso") != 0) {
		els_tq_caloIso_branch = tree->GetBranch(tree->GetAlias("els_tq_caloIso"));
		els_tq_caloIso_branch->SetAddress(&els_tq_caloIso_);
	}
	if(els_tq_caloIso_branch == 0 ) {
	cout << "Branch els_tq_caloIso does not exist." << endl;
	}
	els_tq_egammaEcalIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaEcalIso") != 0) {
		els_tq_egammaEcalIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaEcalIso"));
		els_tq_egammaEcalIso_branch->SetAddress(&els_tq_egammaEcalIso_);
	}
	if(els_tq_egammaEcalIso_branch == 0 ) {
	cout << "Branch els_tq_egammaEcalIso does not exist." << endl;
	}
	els_tq_egammaHcalIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaHcalIso") != 0) {
		els_tq_egammaHcalIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaHcalIso"));
		els_tq_egammaHcalIso_branch->SetAddress(&els_tq_egammaHcalIso_);
	}
	if(els_tq_egammaHcalIso_branch == 0 ) {
	cout << "Branch els_tq_egammaHcalIso does not exist." << endl;
	}
	els_tq_egammaTkIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaTkIso") != 0) {
		els_tq_egammaTkIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaTkIso"));
		els_tq_egammaTkIso_branch->SetAddress(&els_tq_egammaTkIso_);
	}
	if(els_tq_egammaTkIso_branch == 0 ) {
	cout << "Branch els_tq_egammaTkIso does not exist." << endl;
	}
	els_tq_electronIDRobust_branch = 0;
	if (tree->GetAlias("els_tq_electronIDRobust") != 0) {
		els_tq_electronIDRobust_branch = tree->GetBranch(tree->GetAlias("els_tq_electronIDRobust"));
		els_tq_electronIDRobust_branch->SetAddress(&els_tq_electronIDRobust_);
	}
	if(els_tq_electronIDRobust_branch == 0 ) {
	cout << "Branch els_tq_electronIDRobust does not exist." << endl;
	}
	els_tq_leptonID_branch = 0;
	if (tree->GetAlias("els_tq_leptonID") != 0) {
		els_tq_leptonID_branch = tree->GetBranch(tree->GetAlias("els_tq_leptonID"));
		els_tq_leptonID_branch->SetAddress(&els_tq_leptonID_);
	}
	if(els_tq_leptonID_branch == 0 ) {
	cout << "Branch els_tq_leptonID does not exist." << endl;
	}
	els_tq_trackIso_branch = 0;
	if (tree->GetAlias("els_tq_trackIso") != 0) {
		els_tq_trackIso_branch = tree->GetBranch(tree->GetAlias("els_tq_trackIso"));
		els_tq_trackIso_branch->SetAddress(&els_tq_trackIso_);
	}
	if(els_tq_trackIso_branch == 0 ) {
	cout << "Branch els_tq_trackIso does not exist." << endl;
	}
	jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_bCorrF") != 0) {
		jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_bCorrF"));
		jets_tq_bCorrF_branch->SetAddress(&jets_tq_bCorrF_);
	}
	if(jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch jets_tq_bCorrF does not exist." << endl;
	}
	jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_cCorrF") != 0) {
		jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_cCorrF"));
		jets_tq_cCorrF_branch->SetAddress(&jets_tq_cCorrF_);
	}
	if(jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch jets_tq_cCorrF does not exist." << endl;
	}
	jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_gluCorrF") != 0) {
		jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_gluCorrF"));
		jets_tq_gluCorrF_branch->SetAddress(&jets_tq_gluCorrF_);
	}
	if(jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch jets_tq_gluCorrF does not exist." << endl;
	}
	jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("jets_tq_jetCharge") != 0) {
		jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("jets_tq_jetCharge"));
		jets_tq_jetCharge_branch->SetAddress(&jets_tq_jetCharge_);
	}
	if(jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch jets_tq_jetCharge does not exist." << endl;
	}
	jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_noCorrF") != 0) {
		jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_noCorrF"));
		jets_tq_noCorrF_branch->SetAddress(&jets_tq_noCorrF_);
	}
	if(jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch jets_tq_noCorrF does not exist." << endl;
	}
	jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_udsCorrF") != 0) {
		jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_udsCorrF"));
		jets_tq_udsCorrF_branch->SetAddress(&jets_tq_udsCorrF_);
	}
	if(jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch jets_tq_udsCorrF does not exist." << endl;
	}
	mus_tq_caloIso_branch = 0;
	if (tree->GetAlias("mus_tq_caloIso") != 0) {
		mus_tq_caloIso_branch = tree->GetBranch(tree->GetAlias("mus_tq_caloIso"));
		mus_tq_caloIso_branch->SetAddress(&mus_tq_caloIso_);
	}
	if(mus_tq_caloIso_branch == 0 ) {
	cout << "Branch mus_tq_caloIso does not exist." << endl;
	}
	mus_tq_leptonID_branch = 0;
	if (tree->GetAlias("mus_tq_leptonID") != 0) {
		mus_tq_leptonID_branch = tree->GetBranch(tree->GetAlias("mus_tq_leptonID"));
		mus_tq_leptonID_branch->SetAddress(&mus_tq_leptonID_);
	}
	if(mus_tq_leptonID_branch == 0 ) {
	cout << "Branch mus_tq_leptonID does not exist." << endl;
	}
	mus_tq_lrComb_branch = 0;
	if (tree->GetAlias("mus_tq_lrComb") != 0) {
		mus_tq_lrComb_branch = tree->GetBranch(tree->GetAlias("mus_tq_lrComb"));
		mus_tq_lrComb_branch->SetAddress(&mus_tq_lrComb_);
	}
	if(mus_tq_lrComb_branch == 0 ) {
	cout << "Branch mus_tq_lrComb does not exist." << endl;
	}
	mus_tq_trackIso_branch = 0;
	if (tree->GetAlias("mus_tq_trackIso") != 0) {
		mus_tq_trackIso_branch = tree->GetBranch(tree->GetAlias("mus_tq_trackIso"));
		mus_tq_trackIso_branch->SetAddress(&mus_tq_trackIso_);
	}
	if(mus_tq_trackIso_branch == 0 ) {
	cout << "Branch mus_tq_trackIso does not exist." << endl;
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
	trk_elsdr_branch = 0;
	if (tree->GetAlias("trk_elsdr") != 0) {
		trk_elsdr_branch = tree->GetBranch(tree->GetAlias("trk_elsdr"));
		trk_elsdr_branch->SetAddress(&trk_elsdr_);
	}
	if(trk_elsdr_branch == 0 ) {
	cout << "Branch trk_elsdr does not exist." << endl;
	}
	trk_musdr_branch = 0;
	if (tree->GetAlias("trk_musdr") != 0) {
		trk_musdr_branch = tree->GetBranch(tree->GetAlias("trk_musdr"));
		trk_musdr_branch->SetAddress(&trk_musdr_);
	}
	if(trk_musdr_branch == 0 ) {
	cout << "Branch trk_musdr does not exist." << endl;
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
	hyp_jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_bCorrF") != 0) {
		hyp_jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_bCorrF"));
		hyp_jets_tq_bCorrF_branch->SetAddress(&hyp_jets_tq_bCorrF_);
	}
	if(hyp_jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_bCorrF does not exist." << endl;
	}
	hyp_jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_cCorrF") != 0) {
		hyp_jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_cCorrF"));
		hyp_jets_tq_cCorrF_branch->SetAddress(&hyp_jets_tq_cCorrF_);
	}
	if(hyp_jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_cCorrF does not exist." << endl;
	}
	hyp_jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_gluCorrF") != 0) {
		hyp_jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_gluCorrF"));
		hyp_jets_tq_gluCorrF_branch->SetAddress(&hyp_jets_tq_gluCorrF_);
	}
	if(hyp_jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_gluCorrF does not exist." << endl;
	}
	hyp_jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_jetCharge") != 0) {
		hyp_jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_jetCharge"));
		hyp_jets_tq_jetCharge_branch->SetAddress(&hyp_jets_tq_jetCharge_);
	}
	if(hyp_jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch hyp_jets_tq_jetCharge does not exist." << endl;
	}
	hyp_jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_noCorrF") != 0) {
		hyp_jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_noCorrF"));
		hyp_jets_tq_noCorrF_branch->SetAddress(&hyp_jets_tq_noCorrF_);
	}
	if(hyp_jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_noCorrF does not exist." << endl;
	}
	hyp_jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_udsCorrF") != 0) {
		hyp_jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_udsCorrF"));
		hyp_jets_tq_udsCorrF_branch->SetAddress(&hyp_jets_tq_udsCorrF_);
	}
	if(hyp_jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_udsCorrF does not exist." << endl;
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
	hyp_other_jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_bCorrF") != 0) {
		hyp_other_jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_bCorrF"));
		hyp_other_jets_tq_bCorrF_branch->SetAddress(&hyp_other_jets_tq_bCorrF_);
	}
	if(hyp_other_jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_bCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_cCorrF") != 0) {
		hyp_other_jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_cCorrF"));
		hyp_other_jets_tq_cCorrF_branch->SetAddress(&hyp_other_jets_tq_cCorrF_);
	}
	if(hyp_other_jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_cCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_gluCorrF") != 0) {
		hyp_other_jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_gluCorrF"));
		hyp_other_jets_tq_gluCorrF_branch->SetAddress(&hyp_other_jets_tq_gluCorrF_);
	}
	if(hyp_other_jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_gluCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_jetCharge") != 0) {
		hyp_other_jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_jetCharge"));
		hyp_other_jets_tq_jetCharge_branch->SetAddress(&hyp_other_jets_tq_jetCharge_);
	}
	if(hyp_other_jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_jetCharge does not exist." << endl;
	}
	hyp_other_jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_noCorrF") != 0) {
		hyp_other_jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_noCorrF"));
		hyp_other_jets_tq_noCorrF_branch->SetAddress(&hyp_other_jets_tq_noCorrF_);
	}
	if(hyp_other_jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_noCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_udsCorrF") != 0) {
		hyp_other_jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_udsCorrF"));
		hyp_other_jets_tq_udsCorrF_branch->SetAddress(&hyp_other_jets_tq_udsCorrF_);
	}
	if(hyp_other_jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_udsCorrF does not exist." << endl;
	}
	evt_CSA07Process_branch = 0;
	if (tree->GetAlias("evt_CSA07Process") != 0) {
		evt_CSA07Process_branch = tree->GetBranch(tree->GetAlias("evt_CSA07Process"));
		evt_CSA07Process_branch->SetAddress(&evt_CSA07Process_);
	}
	if(evt_CSA07Process_branch == 0 ) {
	cout << "Branch evt_CSA07Process does not exist." << endl;
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
	mus_trkrefkey_branch = 0;
	if (tree->GetAlias("mus_trkrefkey") != 0) {
		mus_trkrefkey_branch = tree->GetBranch(tree->GetAlias("mus_trkrefkey"));
		mus_trkrefkey_branch->SetAddress(&mus_trkrefkey_);
	}
	if(mus_trkrefkey_branch == 0 ) {
	cout << "Branch mus_trkrefkey does not exist." << endl;
	}
	mus_validHits_branch = 0;
	if (tree->GetAlias("mus_validHits") != 0) {
		mus_validHits_branch = tree->GetBranch(tree->GetAlias("mus_validHits"));
		mus_validHits_branch->SetAddress(&mus_validHits_);
	}
	if(mus_validHits_branch == 0 ) {
	cout << "Branch mus_validHits does not exist." << endl;
	}
	els_tq_egammaTkNumIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaTkNumIso") != 0) {
		els_tq_egammaTkNumIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaTkNumIso"));
		els_tq_egammaTkNumIso_branch->SetAddress(&els_tq_egammaTkNumIso_);
	}
	if(els_tq_egammaTkNumIso_branch == 0 ) {
	cout << "Branch els_tq_egammaTkNumIso does not exist." << endl;
	}
	els_tq_genID_branch = 0;
	if (tree->GetAlias("els_tq_genID") != 0) {
		els_tq_genID_branch = tree->GetBranch(tree->GetAlias("els_tq_genID"));
		els_tq_genID_branch->SetAddress(&els_tq_genID_);
	}
	if(els_tq_genID_branch == 0 ) {
	cout << "Branch els_tq_genID does not exist." << endl;
	}
	els_tq_genMotherID_branch = 0;
	if (tree->GetAlias("els_tq_genMotherID") != 0) {
		els_tq_genMotherID_branch = tree->GetBranch(tree->GetAlias("els_tq_genMotherID"));
		els_tq_genMotherID_branch->SetAddress(&els_tq_genMotherID_);
	}
	if(els_tq_genMotherID_branch == 0 ) {
	cout << "Branch els_tq_genMotherID does not exist." << endl;
	}
	jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("jets_tq_genPartonMother_id") != 0) {
		jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("jets_tq_genPartonMother_id"));
		jets_tq_genPartonMother_id_branch->SetAddress(&jets_tq_genPartonMother_id_);
	}
	if(jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch jets_tq_genPartonMother_id does not exist." << endl;
	}
	jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("jets_tq_genParton_id") != 0) {
		jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("jets_tq_genParton_id"));
		jets_tq_genParton_id_branch->SetAddress(&jets_tq_genParton_id_);
	}
	if(jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch jets_tq_genParton_id does not exist." << endl;
	}
	jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("jets_tq_partonFlavour") != 0) {
		jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("jets_tq_partonFlavour"));
		jets_tq_partonFlavour_branch->SetAddress(&jets_tq_partonFlavour_);
	}
	if(jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch jets_tq_partonFlavour does not exist." << endl;
	}
	mus_tq_genID_branch = 0;
	if (tree->GetAlias("mus_tq_genID") != 0) {
		mus_tq_genID_branch = tree->GetBranch(tree->GetAlias("mus_tq_genID"));
		mus_tq_genID_branch->SetAddress(&mus_tq_genID_);
	}
	if(mus_tq_genID_branch == 0 ) {
	cout << "Branch mus_tq_genID does not exist." << endl;
	}
	mus_tq_genMotherID_branch = 0;
	if (tree->GetAlias("mus_tq_genMotherID") != 0) {
		mus_tq_genMotherID_branch = tree->GetBranch(tree->GetAlias("mus_tq_genMotherID"));
		mus_tq_genMotherID_branch->SetAddress(&mus_tq_genMotherID_);
	}
	if(mus_tq_genMotherID_branch == 0 ) {
	cout << "Branch mus_tq_genMotherID does not exist." << endl;
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
	trk_elsidx_branch = 0;
	if (tree->GetAlias("trk_elsidx") != 0) {
		trk_elsidx_branch = tree->GetBranch(tree->GetAlias("trk_elsidx"));
		trk_elsidx_branch->SetAddress(&trk_elsidx_);
	}
	if(trk_elsidx_branch == 0 ) {
	cout << "Branch trk_elsidx does not exist." << endl;
	}
	trk_musidx_branch = 0;
	if (tree->GetAlias("trk_musidx") != 0) {
		trk_musidx_branch = tree->GetBranch(tree->GetAlias("trk_musidx"));
		trk_musidx_branch->SetAddress(&trk_musidx_);
	}
	if(trk_musidx_branch == 0 ) {
	cout << "Branch trk_musidx does not exist." << endl;
	}
	hyp_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_id") != 0) {
		hyp_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_id"));
		hyp_jets_mc_id_branch->SetAddress(&hyp_jets_mc_id_);
	}
	if(hyp_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_jets_mc_id does not exist." << endl;
	}
	hyp_jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genPartonMother_id") != 0) {
		hyp_jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genPartonMother_id"));
		hyp_jets_tq_genPartonMother_id_branch->SetAddress(&hyp_jets_tq_genPartonMother_id_);
	}
	if(hyp_jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genPartonMother_id does not exist." << endl;
	}
	hyp_jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genParton_id") != 0) {
		hyp_jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genParton_id"));
		hyp_jets_tq_genParton_id_branch->SetAddress(&hyp_jets_tq_genParton_id_);
	}
	if(hyp_jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genParton_id does not exist." << endl;
	}
	hyp_jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_partonFlavour") != 0) {
		hyp_jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_partonFlavour"));
		hyp_jets_tq_partonFlavour_branch->SetAddress(&hyp_jets_tq_partonFlavour_);
	}
	if(hyp_jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_jets_tq_partonFlavour does not exist." << endl;
	}
	hyp_other_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_id") != 0) {
		hyp_other_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_id"));
		hyp_other_jets_mc_id_branch->SetAddress(&hyp_other_jets_mc_id_);
	}
	if(hyp_other_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_id does not exist." << endl;
	}
	hyp_other_jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genPartonMother_id") != 0) {
		hyp_other_jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genPartonMother_id"));
		hyp_other_jets_tq_genPartonMother_id_branch->SetAddress(&hyp_other_jets_tq_genPartonMother_id_);
	}
	if(hyp_other_jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genPartonMother_id does not exist." << endl;
	}
	hyp_other_jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genParton_id") != 0) {
		hyp_other_jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genParton_id"));
		hyp_other_jets_tq_genParton_id_branch->SetAddress(&hyp_other_jets_tq_genParton_id_);
	}
	if(hyp_other_jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genParton_id does not exist." << endl;
	}
	hyp_other_jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_partonFlavour") != 0) {
		hyp_other_jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_partonFlavour"));
		hyp_other_jets_tq_partonFlavour_branch->SetAddress(&hyp_other_jets_tq_partonFlavour_);
	}
	if(hyp_other_jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_partonFlavour does not exist." << endl;
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
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		evt_run_branch->SetAddress(&evt_run_);
	}
	if(evt_run_branch == 0 ) {
	cout << "Branch evt_run does not exist." << endl;
	}
	evt_njets_branch = 0;
	if (tree->GetAlias("evt_njets") != 0) {
		evt_njets_branch = tree->GetBranch(tree->GetAlias("evt_njets"));
		evt_njets_branch->SetAddress(&evt_njets_);
	}
	if(evt_njets_branch == 0 ) {
	cout << "Branch evt_njets does not exist." << endl;
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
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		l1met_p4_isLoaded = false;
		els_p4_isLoaded = false;
		els_p4In_isLoaded = false;
		els_p4Out_isLoaded = false;
		els_trk_p4_isLoaded = false;
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
		l1emiso_p4_isLoaded = false;
		l1emnoiso_p4_isLoaded = false;
		l1jetsc_p4_isLoaded = false;
		l1jetsf_p4_isLoaded = false;
		l1jetst_p4_isLoaded = false;
		l1mus_p4_isLoaded = false;
		mus_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		els_tq_genMotherP4_isLoaded = false;
		els_tq_genP4_isLoaded = false;
		jets_tq_genJet_p4_isLoaded = false;
		jets_tq_genPartonMother_p4_isLoaded = false;
		jets_tq_genParton_p4_isLoaded = false;
		jets_tq_jet_p4_isLoaded = false;
		mus_tq_genMotherP4_isLoaded = false;
		mus_tq_genP4_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		hyp_jets_mc_gp_p4_isLoaded = false;
		hyp_jets_mc_p4_isLoaded = false;
		hyp_jets_p4_isLoaded = false;
		hyp_jets_tq_p4_isLoaded = false;
		hyp_jets_tq_genPartonMother_p4_isLoaded = false;
		hyp_jets_tq_genParton_p4_isLoaded = false;
		hyp_jets_tq_jet_p4_isLoaded = false;
		hyp_other_jets_mc_gp_p4_isLoaded = false;
		hyp_other_jets_mc_p4_isLoaded = false;
		hyp_other_jets_p4_isLoaded = false;
		hyp_other_jets_tq_genJet_p4_isLoaded = false;
		hyp_other_jets_tq_genPartonMother_p4_isLoaded = false;
		hyp_other_jets_tq_genParton_p4_isLoaded = false;
		hyp_other_jets_tq_jet_p4_isLoaded = false;
		evt_CSA07FilterEff_isLoaded = false;
		evt_CSA07Pthat_isLoaded = false;
		evt_CSA07Weight_isLoaded = false;
		evt_kfactor_isLoaded = false;
		evt_weight_isLoaded = false;
		evt_xsec_excl_isLoaded = false;
		evt_xsec_incl_isLoaded = false;
		evt_scale1fb_isLoaded = false;
		l1met_etHad_isLoaded = false;
		l1met_etTot_isLoaded = false;
		l1met_met_isLoaded = false;
		evt_met_isLoaded = false;
		evt_metPhi_isLoaded = false;
		evt_metSig_isLoaded = false;
		evt_met_jetcorr_isLoaded = false;
		metphi_jetcorr_isLoaded = false;
		els_musdr_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_ESc_isLoaded = false;
		els_ESc_raw_isLoaded = false;
		els_ESeed_isLoaded = false;
		els_chi2_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_dEtaIn_isLoaded = false;
		els_dEtaOut_isLoaded = false;
		els_dPhiIn_isLoaded = false;
		els_dPhiInPhiOut_isLoaded = false;
		els_dPhiOut_isLoaded = false;
		els_e3x3_isLoaded = false;
		els_e5x5_isLoaded = false;
		els_eOverPIn_isLoaded = false;
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
		els_sigmaPhiPhi_isLoaded = false;
		els_tkIso_isLoaded = false;
		els_vertexphi_isLoaded = false;
		els_z0_isLoaded = false;
		els_z0Err_isLoaded = false;
		hyp_ll_chi2_isLoaded = false;
		hyp_ll_d0_isLoaded = false;
		hyp_ll_d0Err_isLoaded = false;
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
		hyp_lt_chi2_isLoaded = false;
		hyp_lt_d0_isLoaded = false;
		hyp_lt_d0Err_isLoaded = false;
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
		mus_eledr_isLoaded = false;
		mus_jetdr_isLoaded = false;
		mus_trkdr_isLoaded = false;
		mus_chi2_isLoaded = false;
		mus_d0_isLoaded = false;
		mus_d0Err_isLoaded = false;
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
		els_tq_LRComb_isLoaded = false;
		els_tq_caloIso_isLoaded = false;
		els_tq_egammaEcalIso_isLoaded = false;
		els_tq_egammaHcalIso_isLoaded = false;
		els_tq_egammaTkIso_isLoaded = false;
		els_tq_electronIDRobust_isLoaded = false;
		els_tq_leptonID_isLoaded = false;
		els_tq_trackIso_isLoaded = false;
		jets_tq_bCorrF_isLoaded = false;
		jets_tq_cCorrF_isLoaded = false;
		jets_tq_gluCorrF_isLoaded = false;
		jets_tq_jetCharge_isLoaded = false;
		jets_tq_noCorrF_isLoaded = false;
		jets_tq_udsCorrF_isLoaded = false;
		mus_tq_caloIso_isLoaded = false;
		mus_tq_leptonID_isLoaded = false;
		mus_tq_lrComb_isLoaded = false;
		mus_tq_trackIso_isLoaded = false;
		trks_chi2_isLoaded = false;
		trks_d0_isLoaded = false;
		trks_d0Err_isLoaded = false;
		trks_etaErr_isLoaded = false;
		trks_ndof_isLoaded = false;
		trks_outerEta_isLoaded = false;
		trks_outerPhi_isLoaded = false;
		trks_phiErr_isLoaded = false;
		trks_ptErr_isLoaded = false;
		trks_vertexphi_isLoaded = false;
		trks_z0_isLoaded = false;
		trks_z0Err_isLoaded = false;
		trk_elsdr_isLoaded = false;
		trk_musdr_isLoaded = false;
		hyp_jets_EMFcor_isLoaded = false;
		hyp_jets_chFrac_isLoaded = false;
		hyp_jets_cor_isLoaded = false;
		hyp_jets_emFrac_isLoaded = false;
		hyp_jets_mc_emEnergy_isLoaded = false;
		hyp_jets_mc_hadEnergy_isLoaded = false;
		hyp_jets_mc_invEnergy_isLoaded = false;
		hyp_jets_mc_otherEnergy_isLoaded = false;
		hyp_jets_tq_bCorrF_isLoaded = false;
		hyp_jets_tq_cCorrF_isLoaded = false;
		hyp_jets_tq_gluCorrF_isLoaded = false;
		hyp_jets_tq_jetCharge_isLoaded = false;
		hyp_jets_tq_noCorrF_isLoaded = false;
		hyp_jets_tq_udsCorrF_isLoaded = false;
		hyp_other_jets_EMFcor_isLoaded = false;
		hyp_other_jets_chFrac_isLoaded = false;
		hyp_other_jets_cor_isLoaded = false;
		hyp_other_jets_emFrac_isLoaded = false;
		hyp_other_jets_mc_emEnergy_isLoaded = false;
		hyp_other_jets_mc_hadEnergy_isLoaded = false;
		hyp_other_jets_mc_invEnergy_isLoaded = false;
		hyp_other_jets_mc_otherEnergy_isLoaded = false;
		hyp_other_jets_tq_bCorrF_isLoaded = false;
		hyp_other_jets_tq_cCorrF_isLoaded = false;
		hyp_other_jets_tq_gluCorrF_isLoaded = false;
		hyp_other_jets_tq_jetCharge_isLoaded = false;
		hyp_other_jets_tq_noCorrF_isLoaded = false;
		hyp_other_jets_tq_udsCorrF_isLoaded = false;
		evt_CSA07Process_isLoaded = false;
		evt_HLT1_isLoaded = false;
		evt_HLT2_isLoaded = false;
		evt_HLT3_isLoaded = false;
		evt_HLT4_isLoaded = false;
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
		els_closestMuon_isLoaded = false;
		els_trkidx_isLoaded = false;
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
		mus_iso03_ntrk_isLoaded = false;
		mus_iso05_ntrk_isLoaded = false;
		mus_lostHits_isLoaded = false;
		mus_nmatches_isLoaded = false;
		mus_pid_TM2DCompatibilityLoose_isLoaded = false;
		mus_pid_TM2DCompatibilityTight_isLoaded = false;
		mus_pid_TMLastStationLoose_isLoaded = false;
		mus_pid_TMLastStationTight_isLoaded = false;
		mus_trkrefkey_isLoaded = false;
		mus_validHits_isLoaded = false;
		els_tq_egammaTkNumIso_isLoaded = false;
		els_tq_genID_isLoaded = false;
		els_tq_genMotherID_isLoaded = false;
		jets_tq_genPartonMother_id_isLoaded = false;
		jets_tq_genParton_id_isLoaded = false;
		jets_tq_partonFlavour_isLoaded = false;
		mus_tq_genID_isLoaded = false;
		mus_tq_genMotherID_isLoaded = false;
		trks_charge_isLoaded = false;
		trks_lostHits_isLoaded = false;
		trks_validHits_isLoaded = false;
		trk_elsidx_isLoaded = false;
		trk_musidx_isLoaded = false;
		hyp_jets_mc_id_isLoaded = false;
		hyp_jets_tq_genPartonMother_id_isLoaded = false;
		hyp_jets_tq_genParton_id_isLoaded = false;
		hyp_jets_tq_partonFlavour_isLoaded = false;
		hyp_other_jets_mc_id_isLoaded = false;
		hyp_other_jets_tq_genPartonMother_id_isLoaded = false;
		hyp_other_jets_tq_genParton_id_isLoaded = false;
		hyp_other_jets_tq_partonFlavour_isLoaded = false;
		hyp_quadlep_jets_index_isLoaded = false;
		hyp_trilep_jets_index_isLoaded = false;
		evt_nels_isLoaded = false;
		evt_event_isLoaded = false;
		evt_run_isLoaded = false;
		evt_njets_isLoaded = false;
		hyp_quadlep_bucket_isLoaded = false;
		hyp_quadlep_first_index_isLoaded = false;
		hyp_quadlep_fourth_index_isLoaded = false;
		hyp_quadlep_second_index_isLoaded = false;
		hyp_quadlep_third_index_isLoaded = false;
		hyp_trilep_bucket_isLoaded = false;
		hyp_trilep_first_index_isLoaded = false;
		hyp_trilep_second_index_isLoaded = false;
		hyp_trilep_third_index_isLoaded = false;
	}

	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &l1met_p4()
	{
		if (not l1met_p4_isLoaded) {
			if (l1met_p4_branch != 0) 
				l1met_p4_branch->GetEntry(index);
			l1met_p4_isLoaded = true;
		}
		return l1met_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) 
				els_p4_branch->GetEntry(index);
			els_p4_isLoaded = true;
		}
		return els_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4In()
	{
		if (not els_p4In_isLoaded) {
			if (els_p4In_branch != 0) 
				els_p4In_branch->GetEntry(index);
			els_p4In_isLoaded = true;
		}
		return els_p4In_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4Out()
	{
		if (not els_p4Out_isLoaded) {
			if (els_p4Out_branch != 0) 
				els_p4Out_branch->GetEntry(index);
			els_p4Out_isLoaded = true;
		}
		return els_p4Out_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_trk_p4()
	{
		if (not els_trk_p4_isLoaded) {
			if (els_trk_p4_branch != 0) 
				els_trk_p4_branch->GetEntry(index);
			els_trk_p4_isLoaded = true;
		}
		return els_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_p4()
	{
		if (not genps_p4_isLoaded) {
			if (genps_p4_branch != 0) 
				genps_p4_branch->GetEntry(index);
			genps_p4_isLoaded = true;
		}
		return genps_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_prod_vtx()
	{
		if (not genps_prod_vtx_isLoaded) {
			if (genps_prod_vtx_branch != 0) 
				genps_prod_vtx_branch->GetEntry(index);
			genps_prod_vtx_isLoaded = true;
		}
		return genps_prod_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_mc_p4()
	{
		if (not hyp_ll_mc_p4_isLoaded) {
			if (hyp_ll_mc_p4_branch != 0) 
				hyp_ll_mc_p4_branch->GetEntry(index);
			hyp_ll_mc_p4_isLoaded = true;
		}
		return hyp_ll_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_p4()
	{
		if (not hyp_ll_p4_isLoaded) {
			if (hyp_ll_p4_branch != 0) 
				hyp_ll_p4_branch->GetEntry(index);
			hyp_ll_p4_isLoaded = true;
		}
		return hyp_ll_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_trk_p4()
	{
		if (not hyp_ll_trk_p4_isLoaded) {
			if (hyp_ll_trk_p4_branch != 0) 
				hyp_ll_trk_p4_branch->GetEntry(index);
			hyp_ll_trk_p4_isLoaded = true;
		}
		return hyp_ll_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_mc_p4()
	{
		if (not hyp_lt_mc_p4_isLoaded) {
			if (hyp_lt_mc_p4_branch != 0) 
				hyp_lt_mc_p4_branch->GetEntry(index);
			hyp_lt_mc_p4_isLoaded = true;
		}
		return hyp_lt_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_p4()
	{
		if (not hyp_lt_p4_isLoaded) {
			if (hyp_lt_p4_branch != 0) 
				hyp_lt_p4_branch->GetEntry(index);
			hyp_lt_p4_isLoaded = true;
		}
		return hyp_lt_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_trk_p4()
	{
		if (not hyp_lt_trk_p4_isLoaded) {
			if (hyp_lt_trk_p4_branch != 0) 
				hyp_lt_trk_p4_branch->GetEntry(index);
			hyp_lt_trk_p4_isLoaded = true;
		}
		return hyp_lt_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_p4()
	{
		if (not hyp_p4_isLoaded) {
			if (hyp_p4_branch != 0) 
				hyp_p4_branch->GetEntry(index);
			hyp_p4_isLoaded = true;
		}
		return hyp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_p4()
	{
		if (not jets_p4_isLoaded) {
			if (jets_p4_branch != 0) 
				jets_p4_branch->GetEntry(index);
			jets_p4_isLoaded = true;
		}
		return jets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emiso_p4()
	{
		if (not l1emiso_p4_isLoaded) {
			if (l1emiso_p4_branch != 0) 
				l1emiso_p4_branch->GetEntry(index);
			l1emiso_p4_isLoaded = true;
		}
		return l1emiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1emnoiso_p4()
	{
		if (not l1emnoiso_p4_isLoaded) {
			if (l1emnoiso_p4_branch != 0) 
				l1emnoiso_p4_branch->GetEntry(index);
			l1emnoiso_p4_isLoaded = true;
		}
		return l1emnoiso_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsc_p4()
	{
		if (not l1jetsc_p4_isLoaded) {
			if (l1jetsc_p4_branch != 0) 
				l1jetsc_p4_branch->GetEntry(index);
			l1jetsc_p4_isLoaded = true;
		}
		return l1jetsc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetsf_p4()
	{
		if (not l1jetsf_p4_isLoaded) {
			if (l1jetsf_p4_branch != 0) 
				l1jetsf_p4_branch->GetEntry(index);
			l1jetsf_p4_isLoaded = true;
		}
		return l1jetsf_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1jetst_p4()
	{
		if (not l1jetst_p4_isLoaded) {
			if (l1jetst_p4_branch != 0) 
				l1jetst_p4_branch->GetEntry(index);
			l1jetst_p4_isLoaded = true;
		}
		return l1jetst_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &l1mus_p4()
	{
		if (not l1mus_p4_isLoaded) {
			if (l1mus_p4_branch != 0) 
				l1mus_p4_branch->GetEntry(index);
			l1mus_p4_isLoaded = true;
		}
		return l1mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4()
	{
		if (not mus_p4_isLoaded) {
			if (mus_p4_branch != 0) 
				mus_p4_branch->GetEntry(index);
			mus_p4_isLoaded = true;
		}
		return mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4()
	{
		if (not mus_trk_p4_isLoaded) {
			if (mus_trk_p4_branch != 0) 
				mus_trk_p4_branch->GetEntry(index);
			mus_trk_p4_isLoaded = true;
		}
		return mus_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_tq_genMotherP4()
	{
		if (not els_tq_genMotherP4_isLoaded) {
			if (els_tq_genMotherP4_branch != 0) 
				els_tq_genMotherP4_branch->GetEntry(index);
			els_tq_genMotherP4_isLoaded = true;
		}
		return els_tq_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_tq_genP4()
	{
		if (not els_tq_genP4_isLoaded) {
			if (els_tq_genP4_branch != 0) 
				els_tq_genP4_branch->GetEntry(index);
			els_tq_genP4_isLoaded = true;
		}
		return els_tq_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_tq_genJet_p4()
	{
		if (not jets_tq_genJet_p4_isLoaded) {
			if (jets_tq_genJet_p4_branch != 0) 
				jets_tq_genJet_p4_branch->GetEntry(index);
			jets_tq_genJet_p4_isLoaded = true;
		}
		return jets_tq_genJet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_tq_genPartonMother_p4()
	{
		if (not jets_tq_genPartonMother_p4_isLoaded) {
			if (jets_tq_genPartonMother_p4_branch != 0) 
				jets_tq_genPartonMother_p4_branch->GetEntry(index);
			jets_tq_genPartonMother_p4_isLoaded = true;
		}
		return jets_tq_genPartonMother_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_tq_genParton_p4()
	{
		if (not jets_tq_genParton_p4_isLoaded) {
			if (jets_tq_genParton_p4_branch != 0) 
				jets_tq_genParton_p4_branch->GetEntry(index);
			jets_tq_genParton_p4_isLoaded = true;
		}
		return jets_tq_genParton_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_tq_jet_p4()
	{
		if (not jets_tq_jet_p4_isLoaded) {
			if (jets_tq_jet_p4_branch != 0) 
				jets_tq_jet_p4_branch->GetEntry(index);
			jets_tq_jet_p4_isLoaded = true;
		}
		return jets_tq_jet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_tq_genMotherP4()
	{
		if (not mus_tq_genMotherP4_isLoaded) {
			if (mus_tq_genMotherP4_branch != 0) 
				mus_tq_genMotherP4_branch->GetEntry(index);
			mus_tq_genMotherP4_isLoaded = true;
		}
		return mus_tq_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_tq_genP4()
	{
		if (not mus_tq_genP4_isLoaded) {
			if (mus_tq_genP4_branch != 0) 
				mus_tq_genP4_branch->GetEntry(index);
			mus_tq_genP4_isLoaded = true;
		}
		return mus_tq_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4()
	{
		if (not trks_trk_p4_isLoaded) {
			if (trks_trk_p4_branch != 0) 
				trks_trk_p4_branch->GetEntry(index);
			trks_trk_p4_isLoaded = true;
		}
		return trks_trk_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_gp_p4()
	{
		if (not hyp_jets_mc_gp_p4_isLoaded) {
			if (hyp_jets_mc_gp_p4_branch != 0) 
				hyp_jets_mc_gp_p4_branch->GetEntry(index);
			hyp_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_p4()
	{
		if (not hyp_jets_mc_p4_isLoaded) {
			if (hyp_jets_mc_p4_branch != 0) 
				hyp_jets_mc_p4_branch->GetEntry(index);
			hyp_jets_mc_p4_isLoaded = true;
		}
		return hyp_jets_mc_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_p4()
	{
		if (not hyp_jets_p4_isLoaded) {
			if (hyp_jets_p4_branch != 0) 
				hyp_jets_p4_branch->GetEntry(index);
			hyp_jets_p4_isLoaded = true;
		}
		return hyp_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_tq_p4()
	{
		if (not hyp_jets_tq_p4_isLoaded) {
			if (hyp_jets_tq_p4_branch != 0) 
				hyp_jets_tq_p4_branch->GetEntry(index);
			hyp_jets_tq_p4_isLoaded = true;
		}
		return hyp_jets_tq_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_tq_genPartonMother_p4()
	{
		if (not hyp_jets_tq_genPartonMother_p4_isLoaded) {
			if (hyp_jets_tq_genPartonMother_p4_branch != 0) 
				hyp_jets_tq_genPartonMother_p4_branch->GetEntry(index);
			hyp_jets_tq_genPartonMother_p4_isLoaded = true;
		}
		return hyp_jets_tq_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_tq_genParton_p4()
	{
		if (not hyp_jets_tq_genParton_p4_isLoaded) {
			if (hyp_jets_tq_genParton_p4_branch != 0) 
				hyp_jets_tq_genParton_p4_branch->GetEntry(index);
			hyp_jets_tq_genParton_p4_isLoaded = true;
		}
		return hyp_jets_tq_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_tq_jet_p4()
	{
		if (not hyp_jets_tq_jet_p4_isLoaded) {
			if (hyp_jets_tq_jet_p4_branch != 0) 
				hyp_jets_tq_jet_p4_branch->GetEntry(index);
			hyp_jets_tq_jet_p4_isLoaded = true;
		}
		return hyp_jets_tq_jet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_gp_p4()
	{
		if (not hyp_other_jets_mc_gp_p4_isLoaded) {
			if (hyp_other_jets_mc_gp_p4_branch != 0) 
				hyp_other_jets_mc_gp_p4_branch->GetEntry(index);
			hyp_other_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_p4()
	{
		if (not hyp_other_jets_mc_p4_isLoaded) {
			if (hyp_other_jets_mc_p4_branch != 0) 
				hyp_other_jets_mc_p4_branch->GetEntry(index);
			hyp_other_jets_mc_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_p4()
	{
		if (not hyp_other_jets_p4_isLoaded) {
			if (hyp_other_jets_p4_branch != 0) 
				hyp_other_jets_p4_branch->GetEntry(index);
			hyp_other_jets_p4_isLoaded = true;
		}
		return hyp_other_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_tq_genJet_p4()
	{
		if (not hyp_other_jets_tq_genJet_p4_isLoaded) {
			if (hyp_other_jets_tq_genJet_p4_branch != 0) 
				hyp_other_jets_tq_genJet_p4_branch->GetEntry(index);
			hyp_other_jets_tq_genJet_p4_isLoaded = true;
		}
		return hyp_other_jets_tq_genJet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_tq_genPartonMother_p4()
	{
		if (not hyp_other_jets_tq_genPartonMother_p4_isLoaded) {
			if (hyp_other_jets_tq_genPartonMother_p4_branch != 0) 
				hyp_other_jets_tq_genPartonMother_p4_branch->GetEntry(index);
			hyp_other_jets_tq_genPartonMother_p4_isLoaded = true;
		}
		return hyp_other_jets_tq_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_tq_genParton_p4()
	{
		if (not hyp_other_jets_tq_genParton_p4_isLoaded) {
			if (hyp_other_jets_tq_genParton_p4_branch != 0) 
				hyp_other_jets_tq_genParton_p4_branch->GetEntry(index);
			hyp_other_jets_tq_genParton_p4_isLoaded = true;
		}
		return hyp_other_jets_tq_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_tq_jet_p4()
	{
		if (not hyp_other_jets_tq_jet_p4_isLoaded) {
			if (hyp_other_jets_tq_jet_p4_branch != 0) 
				hyp_other_jets_tq_jet_p4_branch->GetEntry(index);
			hyp_other_jets_tq_jet_p4_isLoaded = true;
		}
		return hyp_other_jets_tq_jet_p4_;
	}
	float &evt_CSA07FilterEff()
	{
		if (not evt_CSA07FilterEff_isLoaded) {
			if (evt_CSA07FilterEff_branch != 0) 
				evt_CSA07FilterEff_branch->GetEntry(index);
			evt_CSA07FilterEff_isLoaded = true;
		}
		return evt_CSA07FilterEff_;
	}
	float &evt_CSA07Pthat()
	{
		if (not evt_CSA07Pthat_isLoaded) {
			if (evt_CSA07Pthat_branch != 0) 
				evt_CSA07Pthat_branch->GetEntry(index);
			evt_CSA07Pthat_isLoaded = true;
		}
		return evt_CSA07Pthat_;
	}
	float &evt_CSA07Weight()
	{
		if (not evt_CSA07Weight_isLoaded) {
			if (evt_CSA07Weight_branch != 0) 
				evt_CSA07Weight_branch->GetEntry(index);
			evt_CSA07Weight_isLoaded = true;
		}
		return evt_CSA07Weight_;
	}
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) 
				evt_kfactor_branch->GetEntry(index);
			evt_kfactor_isLoaded = true;
		}
		return evt_kfactor_;
	}
	float &evt_weight()
	{
		if (not evt_weight_isLoaded) {
			if (evt_weight_branch != 0) 
				evt_weight_branch->GetEntry(index);
			evt_weight_isLoaded = true;
		}
		return evt_weight_;
	}
	float &evt_xsec_excl()
	{
		if (not evt_xsec_excl_isLoaded) {
			if (evt_xsec_excl_branch != 0) 
				evt_xsec_excl_branch->GetEntry(index);
			evt_xsec_excl_isLoaded = true;
		}
		return evt_xsec_excl_;
	}
	float &evt_xsec_incl()
	{
		if (not evt_xsec_incl_isLoaded) {
			if (evt_xsec_incl_branch != 0) 
				evt_xsec_incl_branch->GetEntry(index);
			evt_xsec_incl_isLoaded = true;
		}
		return evt_xsec_incl_;
	}
	float &evt_scale1fb()
	{
		if (not evt_scale1fb_isLoaded) {
			if (evt_scale1fb_branch != 0) 
				evt_scale1fb_branch->GetEntry(index);
			evt_scale1fb_isLoaded = true;
		}
		return evt_scale1fb_;
	}
	float &l1met_etHad()
	{
		if (not l1met_etHad_isLoaded) {
			if (l1met_etHad_branch != 0) 
				l1met_etHad_branch->GetEntry(index);
			l1met_etHad_isLoaded = true;
		}
		return l1met_etHad_;
	}
	float &l1met_etTot()
	{
		if (not l1met_etTot_isLoaded) {
			if (l1met_etTot_branch != 0) 
				l1met_etTot_branch->GetEntry(index);
			l1met_etTot_isLoaded = true;
		}
		return l1met_etTot_;
	}
	float &l1met_met()
	{
		if (not l1met_met_isLoaded) {
			if (l1met_met_branch != 0) 
				l1met_met_branch->GetEntry(index);
			l1met_met_isLoaded = true;
		}
		return l1met_met_;
	}
	float &evt_met()
	{
		if (not evt_met_isLoaded) {
			if (evt_met_branch != 0) 
				evt_met_branch->GetEntry(index);
			evt_met_isLoaded = true;
		}
		return evt_met_;
	}
	float &evt_metPhi()
	{
		if (not evt_metPhi_isLoaded) {
			if (evt_metPhi_branch != 0) 
				evt_metPhi_branch->GetEntry(index);
			evt_metPhi_isLoaded = true;
		}
		return evt_metPhi_;
	}
	float &evt_metSig()
	{
		if (not evt_metSig_isLoaded) {
			if (evt_metSig_branch != 0) 
				evt_metSig_branch->GetEntry(index);
			evt_metSig_isLoaded = true;
		}
		return evt_metSig_;
	}
	float &evt_met_jetcorr()
	{
		if (not evt_met_jetcorr_isLoaded) {
			if (evt_met_jetcorr_branch != 0) 
				evt_met_jetcorr_branch->GetEntry(index);
			evt_met_jetcorr_isLoaded = true;
		}
		return evt_met_jetcorr_;
	}
	float &metphi_jetcorr()
	{
		if (not metphi_jetcorr_isLoaded) {
			if (metphi_jetcorr_branch != 0) 
				metphi_jetcorr_branch->GetEntry(index);
			metphi_jetcorr_isLoaded = true;
		}
		return metphi_jetcorr_;
	}
	vector<float> &els_musdr()
	{
		if (not els_musdr_isLoaded) {
			if (els_musdr_branch != 0) 
				els_musdr_branch->GetEntry(index);
			els_musdr_isLoaded = true;
		}
		return els_musdr_;
	}
	vector<float> &els_trkdr()
	{
		if (not els_trkdr_isLoaded) {
			if (els_trkdr_branch != 0) 
				els_trkdr_branch->GetEntry(index);
			els_trkdr_isLoaded = true;
		}
		return els_trkdr_;
	}
	vector<float> &els_ESc()
	{
		if (not els_ESc_isLoaded) {
			if (els_ESc_branch != 0) 
				els_ESc_branch->GetEntry(index);
			els_ESc_isLoaded = true;
		}
		return els_ESc_;
	}
	vector<float> &els_ESc_raw()
	{
		if (not els_ESc_raw_isLoaded) {
			if (els_ESc_raw_branch != 0) 
				els_ESc_raw_branch->GetEntry(index);
			els_ESc_raw_isLoaded = true;
		}
		return els_ESc_raw_;
	}
	vector<float> &els_ESeed()
	{
		if (not els_ESeed_isLoaded) {
			if (els_ESeed_branch != 0) 
				els_ESeed_branch->GetEntry(index);
			els_ESeed_isLoaded = true;
		}
		return els_ESeed_;
	}
	vector<float> &els_chi2()
	{
		if (not els_chi2_isLoaded) {
			if (els_chi2_branch != 0) 
				els_chi2_branch->GetEntry(index);
			els_chi2_isLoaded = true;
		}
		return els_chi2_;
	}
	vector<float> &els_d0()
	{
		if (not els_d0_isLoaded) {
			if (els_d0_branch != 0) 
				els_d0_branch->GetEntry(index);
			els_d0_isLoaded = true;
		}
		return els_d0_;
	}
	vector<float> &els_d0Err()
	{
		if (not els_d0Err_isLoaded) {
			if (els_d0Err_branch != 0) 
				els_d0Err_branch->GetEntry(index);
			els_d0Err_isLoaded = true;
		}
		return els_d0Err_;
	}
	vector<float> &els_dEtaIn()
	{
		if (not els_dEtaIn_isLoaded) {
			if (els_dEtaIn_branch != 0) 
				els_dEtaIn_branch->GetEntry(index);
			els_dEtaIn_isLoaded = true;
		}
		return els_dEtaIn_;
	}
	vector<float> &els_dEtaOut()
	{
		if (not els_dEtaOut_isLoaded) {
			if (els_dEtaOut_branch != 0) 
				els_dEtaOut_branch->GetEntry(index);
			els_dEtaOut_isLoaded = true;
		}
		return els_dEtaOut_;
	}
	vector<float> &els_dPhiIn()
	{
		if (not els_dPhiIn_isLoaded) {
			if (els_dPhiIn_branch != 0) 
				els_dPhiIn_branch->GetEntry(index);
			els_dPhiIn_isLoaded = true;
		}
		return els_dPhiIn_;
	}
	vector<float> &els_dPhiInPhiOut()
	{
		if (not els_dPhiInPhiOut_isLoaded) {
			if (els_dPhiInPhiOut_branch != 0) 
				els_dPhiInPhiOut_branch->GetEntry(index);
			els_dPhiInPhiOut_isLoaded = true;
		}
		return els_dPhiInPhiOut_;
	}
	vector<float> &els_dPhiOut()
	{
		if (not els_dPhiOut_isLoaded) {
			if (els_dPhiOut_branch != 0) 
				els_dPhiOut_branch->GetEntry(index);
			els_dPhiOut_isLoaded = true;
		}
		return els_dPhiOut_;
	}
	vector<float> &els_e3x3()
	{
		if (not els_e3x3_isLoaded) {
			if (els_e3x3_branch != 0) 
				els_e3x3_branch->GetEntry(index);
			els_e3x3_isLoaded = true;
		}
		return els_e3x3_;
	}
	vector<float> &els_e5x5()
	{
		if (not els_e5x5_isLoaded) {
			if (els_e5x5_branch != 0) 
				els_e5x5_branch->GetEntry(index);
			els_e5x5_isLoaded = true;
		}
		return els_e5x5_;
	}
	vector<float> &els_eOverPIn()
	{
		if (not els_eOverPIn_isLoaded) {
			if (els_eOverPIn_branch != 0) 
				els_eOverPIn_branch->GetEntry(index);
			els_eOverPIn_isLoaded = true;
		}
		return els_eOverPIn_;
	}
	vector<float> &els_eSeedOverPOut()
	{
		if (not els_eSeedOverPOut_isLoaded) {
			if (els_eSeedOverPOut_branch != 0) 
				els_eSeedOverPOut_branch->GetEntry(index);
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	vector<float> &els_etaErr()
	{
		if (not els_etaErr_isLoaded) {
			if (els_etaErr_branch != 0) 
				els_etaErr_branch->GetEntry(index);
			els_etaErr_isLoaded = true;
		}
		return els_etaErr_;
	}
	vector<float> &els_fBrem()
	{
		if (not els_fBrem_isLoaded) {
			if (els_fBrem_branch != 0) 
				els_fBrem_branch->GetEntry(index);
			els_fBrem_isLoaded = true;
		}
		return els_fBrem_;
	}
	vector<float> &els_hOverE()
	{
		if (not els_hOverE_isLoaded) {
			if (els_hOverE_branch != 0) 
				els_hOverE_branch->GetEntry(index);
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	vector<float> &els_ndof()
	{
		if (not els_ndof_isLoaded) {
			if (els_ndof_branch != 0) 
				els_ndof_branch->GetEntry(index);
			els_ndof_isLoaded = true;
		}
		return els_ndof_;
	}
	vector<float> &els_outerEta()
	{
		if (not els_outerEta_isLoaded) {
			if (els_outerEta_branch != 0) 
				els_outerEta_branch->GetEntry(index);
			els_outerEta_isLoaded = true;
		}
		return els_outerEta_;
	}
	vector<float> &els_outerPhi()
	{
		if (not els_outerPhi_isLoaded) {
			if (els_outerPhi_branch != 0) 
				els_outerPhi_branch->GetEntry(index);
			els_outerPhi_isLoaded = true;
		}
		return els_outerPhi_;
	}
	vector<float> &els_phiErr()
	{
		if (not els_phiErr_isLoaded) {
			if (els_phiErr_branch != 0) 
				els_phiErr_branch->GetEntry(index);
			els_phiErr_isLoaded = true;
		}
		return els_phiErr_;
	}
	vector<float> &els_ptErr()
	{
		if (not els_ptErr_isLoaded) {
			if (els_ptErr_branch != 0) 
				els_ptErr_branch->GetEntry(index);
			els_ptErr_isLoaded = true;
		}
		return els_ptErr_;
	}
	vector<float> &els_sigmaEtaEta()
	{
		if (not els_sigmaEtaEta_isLoaded) {
			if (els_sigmaEtaEta_branch != 0) 
				els_sigmaEtaEta_branch->GetEntry(index);
			els_sigmaEtaEta_isLoaded = true;
		}
		return els_sigmaEtaEta_;
	}
	vector<float> &els_sigmaPhiPhi()
	{
		if (not els_sigmaPhiPhi_isLoaded) {
			if (els_sigmaPhiPhi_branch != 0) 
				els_sigmaPhiPhi_branch->GetEntry(index);
			els_sigmaPhiPhi_isLoaded = true;
		}
		return els_sigmaPhiPhi_;
	}
	vector<float> &els_tkIso()
	{
		if (not els_tkIso_isLoaded) {
			if (els_tkIso_branch != 0) 
				els_tkIso_branch->GetEntry(index);
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
	}
	vector<float> &els_vertexphi()
	{
		if (not els_vertexphi_isLoaded) {
			if (els_vertexphi_branch != 0) 
				els_vertexphi_branch->GetEntry(index);
			els_vertexphi_isLoaded = true;
		}
		return els_vertexphi_;
	}
	vector<float> &els_z0()
	{
		if (not els_z0_isLoaded) {
			if (els_z0_branch != 0) 
				els_z0_branch->GetEntry(index);
			els_z0_isLoaded = true;
		}
		return els_z0_;
	}
	vector<float> &els_z0Err()
	{
		if (not els_z0Err_isLoaded) {
			if (els_z0Err_branch != 0) 
				els_z0Err_branch->GetEntry(index);
			els_z0Err_isLoaded = true;
		}
		return els_z0Err_;
	}
	vector<float> &hyp_ll_chi2()
	{
		if (not hyp_ll_chi2_isLoaded) {
			if (hyp_ll_chi2_branch != 0) 
				hyp_ll_chi2_branch->GetEntry(index);
			hyp_ll_chi2_isLoaded = true;
		}
		return hyp_ll_chi2_;
	}
	vector<float> &hyp_ll_d0()
	{
		if (not hyp_ll_d0_isLoaded) {
			if (hyp_ll_d0_branch != 0) 
				hyp_ll_d0_branch->GetEntry(index);
			hyp_ll_d0_isLoaded = true;
		}
		return hyp_ll_d0_;
	}
	vector<float> &hyp_ll_d0Err()
	{
		if (not hyp_ll_d0Err_isLoaded) {
			if (hyp_ll_d0Err_branch != 0) 
				hyp_ll_d0Err_branch->GetEntry(index);
			hyp_ll_d0Err_isLoaded = true;
		}
		return hyp_ll_d0Err_;
	}
	vector<float> &hyp_ll_etaErr()
	{
		if (not hyp_ll_etaErr_isLoaded) {
			if (hyp_ll_etaErr_branch != 0) 
				hyp_ll_etaErr_branch->GetEntry(index);
			hyp_ll_etaErr_isLoaded = true;
		}
		return hyp_ll_etaErr_;
	}
	vector<float> &hyp_ll_iso()
	{
		if (not hyp_ll_iso_isLoaded) {
			if (hyp_ll_iso_branch != 0) 
				hyp_ll_iso_branch->GetEntry(index);
			hyp_ll_iso_isLoaded = true;
		}
		return hyp_ll_iso_;
	}
	vector<float> &hyp_ll_ndof()
	{
		if (not hyp_ll_ndof_isLoaded) {
			if (hyp_ll_ndof_branch != 0) 
				hyp_ll_ndof_branch->GetEntry(index);
			hyp_ll_ndof_isLoaded = true;
		}
		return hyp_ll_ndof_;
	}
	vector<float> &hyp_ll_outerEta()
	{
		if (not hyp_ll_outerEta_isLoaded) {
			if (hyp_ll_outerEta_branch != 0) 
				hyp_ll_outerEta_branch->GetEntry(index);
			hyp_ll_outerEta_isLoaded = true;
		}
		return hyp_ll_outerEta_;
	}
	vector<float> &hyp_ll_outerPhi()
	{
		if (not hyp_ll_outerPhi_isLoaded) {
			if (hyp_ll_outerPhi_branch != 0) 
				hyp_ll_outerPhi_branch->GetEntry(index);
			hyp_ll_outerPhi_isLoaded = true;
		}
		return hyp_ll_outerPhi_;
	}
	vector<float> &hyp_ll_phiErr()
	{
		if (not hyp_ll_phiErr_isLoaded) {
			if (hyp_ll_phiErr_branch != 0) 
				hyp_ll_phiErr_branch->GetEntry(index);
			hyp_ll_phiErr_isLoaded = true;
		}
		return hyp_ll_phiErr_;
	}
	vector<float> &hyp_ll_ptErr()
	{
		if (not hyp_ll_ptErr_isLoaded) {
			if (hyp_ll_ptErr_branch != 0) 
				hyp_ll_ptErr_branch->GetEntry(index);
			hyp_ll_ptErr_isLoaded = true;
		}
		return hyp_ll_ptErr_;
	}
	vector<float> &hyp_ll_tkIso()
	{
		if (not hyp_ll_tkIso_isLoaded) {
			if (hyp_ll_tkIso_branch != 0) 
				hyp_ll_tkIso_branch->GetEntry(index);
			hyp_ll_tkIso_isLoaded = true;
		}
		return hyp_ll_tkIso_;
	}
	vector<float> &hyp_ll_vertexphi()
	{
		if (not hyp_ll_vertexphi_isLoaded) {
			if (hyp_ll_vertexphi_branch != 0) 
				hyp_ll_vertexphi_branch->GetEntry(index);
			hyp_ll_vertexphi_isLoaded = true;
		}
		return hyp_ll_vertexphi_;
	}
	vector<float> &hyp_ll_z0()
	{
		if (not hyp_ll_z0_isLoaded) {
			if (hyp_ll_z0_branch != 0) 
				hyp_ll_z0_branch->GetEntry(index);
			hyp_ll_z0_isLoaded = true;
		}
		return hyp_ll_z0_;
	}
	vector<float> &hyp_ll_z0Err()
	{
		if (not hyp_ll_z0Err_isLoaded) {
			if (hyp_ll_z0Err_branch != 0) 
				hyp_ll_z0Err_branch->GetEntry(index);
			hyp_ll_z0Err_isLoaded = true;
		}
		return hyp_ll_z0Err_;
	}
	vector<float> &hyp_lt_chi2()
	{
		if (not hyp_lt_chi2_isLoaded) {
			if (hyp_lt_chi2_branch != 0) 
				hyp_lt_chi2_branch->GetEntry(index);
			hyp_lt_chi2_isLoaded = true;
		}
		return hyp_lt_chi2_;
	}
	vector<float> &hyp_lt_d0()
	{
		if (not hyp_lt_d0_isLoaded) {
			if (hyp_lt_d0_branch != 0) 
				hyp_lt_d0_branch->GetEntry(index);
			hyp_lt_d0_isLoaded = true;
		}
		return hyp_lt_d0_;
	}
	vector<float> &hyp_lt_d0Err()
	{
		if (not hyp_lt_d0Err_isLoaded) {
			if (hyp_lt_d0Err_branch != 0) 
				hyp_lt_d0Err_branch->GetEntry(index);
			hyp_lt_d0Err_isLoaded = true;
		}
		return hyp_lt_d0Err_;
	}
	vector<float> &hyp_lt_etaErr()
	{
		if (not hyp_lt_etaErr_isLoaded) {
			if (hyp_lt_etaErr_branch != 0) 
				hyp_lt_etaErr_branch->GetEntry(index);
			hyp_lt_etaErr_isLoaded = true;
		}
		return hyp_lt_etaErr_;
	}
	vector<float> &hyp_lt_iso()
	{
		if (not hyp_lt_iso_isLoaded) {
			if (hyp_lt_iso_branch != 0) 
				hyp_lt_iso_branch->GetEntry(index);
			hyp_lt_iso_isLoaded = true;
		}
		return hyp_lt_iso_;
	}
	vector<float> &hyp_lt_ndof()
	{
		if (not hyp_lt_ndof_isLoaded) {
			if (hyp_lt_ndof_branch != 0) 
				hyp_lt_ndof_branch->GetEntry(index);
			hyp_lt_ndof_isLoaded = true;
		}
		return hyp_lt_ndof_;
	}
	vector<float> &hyp_lt_outerEta()
	{
		if (not hyp_lt_outerEta_isLoaded) {
			if (hyp_lt_outerEta_branch != 0) 
				hyp_lt_outerEta_branch->GetEntry(index);
			hyp_lt_outerEta_isLoaded = true;
		}
		return hyp_lt_outerEta_;
	}
	vector<float> &hyp_lt_outerPhi()
	{
		if (not hyp_lt_outerPhi_isLoaded) {
			if (hyp_lt_outerPhi_branch != 0) 
				hyp_lt_outerPhi_branch->GetEntry(index);
			hyp_lt_outerPhi_isLoaded = true;
		}
		return hyp_lt_outerPhi_;
	}
	vector<float> &hyp_lt_phiErr()
	{
		if (not hyp_lt_phiErr_isLoaded) {
			if (hyp_lt_phiErr_branch != 0) 
				hyp_lt_phiErr_branch->GetEntry(index);
			hyp_lt_phiErr_isLoaded = true;
		}
		return hyp_lt_phiErr_;
	}
	vector<float> &hyp_lt_ptErr()
	{
		if (not hyp_lt_ptErr_isLoaded) {
			if (hyp_lt_ptErr_branch != 0) 
				hyp_lt_ptErr_branch->GetEntry(index);
			hyp_lt_ptErr_isLoaded = true;
		}
		return hyp_lt_ptErr_;
	}
	vector<float> &hyp_lt_tkIso()
	{
		if (not hyp_lt_tkIso_isLoaded) {
			if (hyp_lt_tkIso_branch != 0) 
				hyp_lt_tkIso_branch->GetEntry(index);
			hyp_lt_tkIso_isLoaded = true;
		}
		return hyp_lt_tkIso_;
	}
	vector<float> &hyp_lt_vertexphi()
	{
		if (not hyp_lt_vertexphi_isLoaded) {
			if (hyp_lt_vertexphi_branch != 0) 
				hyp_lt_vertexphi_branch->GetEntry(index);
			hyp_lt_vertexphi_isLoaded = true;
		}
		return hyp_lt_vertexphi_;
	}
	vector<float> &hyp_lt_z0()
	{
		if (not hyp_lt_z0_isLoaded) {
			if (hyp_lt_z0_branch != 0) 
				hyp_lt_z0_branch->GetEntry(index);
			hyp_lt_z0_isLoaded = true;
		}
		return hyp_lt_z0_;
	}
	vector<float> &hyp_lt_z0Err()
	{
		if (not hyp_lt_z0Err_isLoaded) {
			if (hyp_lt_z0Err_branch != 0) 
				hyp_lt_z0Err_branch->GetEntry(index);
			hyp_lt_z0Err_isLoaded = true;
		}
		return hyp_lt_z0Err_;
	}
	vector<float> &hyp_met()
	{
		if (not hyp_met_isLoaded) {
			if (hyp_met_branch != 0) 
				hyp_met_branch->GetEntry(index);
			hyp_met_isLoaded = true;
		}
		return hyp_met_;
	}
	vector<float> &hyp_metAll()
	{
		if (not hyp_metAll_isLoaded) {
			if (hyp_metAll_branch != 0) 
				hyp_metAll_branch->GetEntry(index);
			hyp_metAll_isLoaded = true;
		}
		return hyp_metAll_;
	}
	vector<float> &hyp_metAllCaloExp()
	{
		if (not hyp_metAllCaloExp_isLoaded) {
			if (hyp_metAllCaloExp_branch != 0) 
				hyp_metAllCaloExp_branch->GetEntry(index);
			hyp_metAllCaloExp_isLoaded = true;
		}
		return hyp_metAllCaloExp_;
	}
	vector<float> &hyp_metCaloExp()
	{
		if (not hyp_metCaloExp_isLoaded) {
			if (hyp_metCaloExp_branch != 0) 
				hyp_metCaloExp_branch->GetEntry(index);
			hyp_metCaloExp_isLoaded = true;
		}
		return hyp_metCaloExp_;
	}
	vector<float> &hyp_metCone()
	{
		if (not hyp_metCone_isLoaded) {
			if (hyp_metCone_branch != 0) 
				hyp_metCone_branch->GetEntry(index);
			hyp_metCone_isLoaded = true;
		}
		return hyp_metCone_;
	}
	vector<float> &hyp_metDPhiJet10()
	{
		if (not hyp_metDPhiJet10_isLoaded) {
			if (hyp_metDPhiJet10_branch != 0) 
				hyp_metDPhiJet10_branch->GetEntry(index);
			hyp_metDPhiJet10_isLoaded = true;
		}
		return hyp_metDPhiJet10_;
	}
	vector<float> &hyp_metDPhiJet15()
	{
		if (not hyp_metDPhiJet15_isLoaded) {
			if (hyp_metDPhiJet15_branch != 0) 
				hyp_metDPhiJet15_branch->GetEntry(index);
			hyp_metDPhiJet15_isLoaded = true;
		}
		return hyp_metDPhiJet15_;
	}
	vector<float> &hyp_metDPhiJet20()
	{
		if (not hyp_metDPhiJet20_isLoaded) {
			if (hyp_metDPhiJet20_branch != 0) 
				hyp_metDPhiJet20_branch->GetEntry(index);
			hyp_metDPhiJet20_isLoaded = true;
		}
		return hyp_metDPhiJet20_;
	}
	vector<float> &hyp_metDPhiTrk10()
	{
		if (not hyp_metDPhiTrk10_isLoaded) {
			if (hyp_metDPhiTrk10_branch != 0) 
				hyp_metDPhiTrk10_branch->GetEntry(index);
			hyp_metDPhiTrk10_isLoaded = true;
		}
		return hyp_metDPhiTrk10_;
	}
	vector<float> &hyp_metDPhiTrk25()
	{
		if (not hyp_metDPhiTrk25_isLoaded) {
			if (hyp_metDPhiTrk25_branch != 0) 
				hyp_metDPhiTrk25_branch->GetEntry(index);
			hyp_metDPhiTrk25_isLoaded = true;
		}
		return hyp_metDPhiTrk25_;
	}
	vector<float> &hyp_metDPhiTrk50()
	{
		if (not hyp_metDPhiTrk50_isLoaded) {
			if (hyp_metDPhiTrk50_branch != 0) 
				hyp_metDPhiTrk50_branch->GetEntry(index);
			hyp_metDPhiTrk50_isLoaded = true;
		}
		return hyp_metDPhiTrk50_;
	}
	vector<float> &hyp_metJes10()
	{
		if (not hyp_metJes10_isLoaded) {
			if (hyp_metJes10_branch != 0) 
				hyp_metJes10_branch->GetEntry(index);
			hyp_metJes10_isLoaded = true;
		}
		return hyp_metJes10_;
	}
	vector<float> &hyp_metJes15()
	{
		if (not hyp_metJes15_isLoaded) {
			if (hyp_metJes15_branch != 0) 
				hyp_metJes15_branch->GetEntry(index);
			hyp_metJes15_isLoaded = true;
		}
		return hyp_metJes15_;
	}
	vector<float> &hyp_metJes30()
	{
		if (not hyp_metJes30_isLoaded) {
			if (hyp_metJes30_branch != 0) 
				hyp_metJes30_branch->GetEntry(index);
			hyp_metJes30_isLoaded = true;
		}
		return hyp_metJes30_;
	}
	vector<float> &hyp_metJes5()
	{
		if (not hyp_metJes5_isLoaded) {
			if (hyp_metJes5_branch != 0) 
				hyp_metJes5_branch->GetEntry(index);
			hyp_metJes5_isLoaded = true;
		}
		return hyp_metJes5_;
	}
	vector<float> &hyp_metJes50()
	{
		if (not hyp_metJes50_isLoaded) {
			if (hyp_metJes50_branch != 0) 
				hyp_metJes50_branch->GetEntry(index);
			hyp_metJes50_isLoaded = true;
		}
		return hyp_metJes50_;
	}
	vector<float> &hyp_metNoCalo()
	{
		if (not hyp_metNoCalo_isLoaded) {
			if (hyp_metNoCalo_branch != 0) 
				hyp_metNoCalo_branch->GetEntry(index);
			hyp_metNoCalo_isLoaded = true;
		}
		return hyp_metNoCalo_;
	}
	vector<float> &hyp_metPhi()
	{
		if (not hyp_metPhi_isLoaded) {
			if (hyp_metPhi_branch != 0) 
				hyp_metPhi_branch->GetEntry(index);
			hyp_metPhi_isLoaded = true;
		}
		return hyp_metPhi_;
	}
	vector<float> &hyp_metPhiAll()
	{
		if (not hyp_metPhiAll_isLoaded) {
			if (hyp_metPhiAll_branch != 0) 
				hyp_metPhiAll_branch->GetEntry(index);
			hyp_metPhiAll_isLoaded = true;
		}
		return hyp_metPhiAll_;
	}
	vector<float> &hyp_metPhiAllCaloExp()
	{
		if (not hyp_metPhiAllCaloExp_isLoaded) {
			if (hyp_metPhiAllCaloExp_branch != 0) 
				hyp_metPhiAllCaloExp_branch->GetEntry(index);
			hyp_metPhiAllCaloExp_isLoaded = true;
		}
		return hyp_metPhiAllCaloExp_;
	}
	vector<float> &hyp_metPhiCaloExp()
	{
		if (not hyp_metPhiCaloExp_isLoaded) {
			if (hyp_metPhiCaloExp_branch != 0) 
				hyp_metPhiCaloExp_branch->GetEntry(index);
			hyp_metPhiCaloExp_isLoaded = true;
		}
		return hyp_metPhiCaloExp_;
	}
	vector<float> &hyp_metPhiCone()
	{
		if (not hyp_metPhiCone_isLoaded) {
			if (hyp_metPhiCone_branch != 0) 
				hyp_metPhiCone_branch->GetEntry(index);
			hyp_metPhiCone_isLoaded = true;
		}
		return hyp_metPhiCone_;
	}
	vector<float> &hyp_metPhiJes10()
	{
		if (not hyp_metPhiJes10_isLoaded) {
			if (hyp_metPhiJes10_branch != 0) 
				hyp_metPhiJes10_branch->GetEntry(index);
			hyp_metPhiJes10_isLoaded = true;
		}
		return hyp_metPhiJes10_;
	}
	vector<float> &hyp_metPhiJes15()
	{
		if (not hyp_metPhiJes15_isLoaded) {
			if (hyp_metPhiJes15_branch != 0) 
				hyp_metPhiJes15_branch->GetEntry(index);
			hyp_metPhiJes15_isLoaded = true;
		}
		return hyp_metPhiJes15_;
	}
	vector<float> &hyp_metPhiJes30()
	{
		if (not hyp_metPhiJes30_isLoaded) {
			if (hyp_metPhiJes30_branch != 0) 
				hyp_metPhiJes30_branch->GetEntry(index);
			hyp_metPhiJes30_isLoaded = true;
		}
		return hyp_metPhiJes30_;
	}
	vector<float> &hyp_metPhiJes5()
	{
		if (not hyp_metPhiJes5_isLoaded) {
			if (hyp_metPhiJes5_branch != 0) 
				hyp_metPhiJes5_branch->GetEntry(index);
			hyp_metPhiJes5_isLoaded = true;
		}
		return hyp_metPhiJes5_;
	}
	vector<float> &hyp_metPhiJes50()
	{
		if (not hyp_metPhiJes50_isLoaded) {
			if (hyp_metPhiJes50_branch != 0) 
				hyp_metPhiJes50_branch->GetEntry(index);
			hyp_metPhiJes50_isLoaded = true;
		}
		return hyp_metPhiJes50_;
	}
	vector<float> &hyp_metPhiNoCalo()
	{
		if (not hyp_metPhiNoCalo_isLoaded) {
			if (hyp_metPhiNoCalo_branch != 0) 
				hyp_metPhiNoCalo_branch->GetEntry(index);
			hyp_metPhiNoCalo_isLoaded = true;
		}
		return hyp_metPhiNoCalo_;
	}
	vector<float> &hyp_quadlep_met()
	{
		if (not hyp_quadlep_met_isLoaded) {
			if (hyp_quadlep_met_branch != 0) 
				hyp_quadlep_met_branch->GetEntry(index);
			hyp_quadlep_met_isLoaded = true;
		}
		return hyp_quadlep_met_;
	}
	vector<float> &hyp_quadlep_metAll()
	{
		if (not hyp_quadlep_metAll_isLoaded) {
			if (hyp_quadlep_metAll_branch != 0) 
				hyp_quadlep_metAll_branch->GetEntry(index);
			hyp_quadlep_metAll_isLoaded = true;
		}
		return hyp_quadlep_metAll_;
	}
	vector<float> &hyp_trilep_met()
	{
		if (not hyp_trilep_met_isLoaded) {
			if (hyp_trilep_met_branch != 0) 
				hyp_trilep_met_branch->GetEntry(index);
			hyp_trilep_met_isLoaded = true;
		}
		return hyp_trilep_met_;
	}
	vector<float> &hyp_trilep_metAll()
	{
		if (not hyp_trilep_metAll_isLoaded) {
			if (hyp_trilep_metAll_branch != 0) 
				hyp_trilep_metAll_branch->GetEntry(index);
			hyp_trilep_metAll_isLoaded = true;
		}
		return hyp_trilep_metAll_;
	}
	vector<float> &jets_EMFcor()
	{
		if (not jets_EMFcor_isLoaded) {
			if (jets_EMFcor_branch != 0) 
				jets_EMFcor_branch->GetEntry(index);
			jets_EMFcor_isLoaded = true;
		}
		return jets_EMFcor_;
	}
	vector<float> &jets_chFrac()
	{
		if (not jets_chFrac_isLoaded) {
			if (jets_chFrac_branch != 0) 
				jets_chFrac_branch->GetEntry(index);
			jets_chFrac_isLoaded = true;
		}
		return jets_chFrac_;
	}
	vector<float> &jets_cor()
	{
		if (not jets_cor_isLoaded) {
			if (jets_cor_branch != 0) 
				jets_cor_branch->GetEntry(index);
			jets_cor_isLoaded = true;
		}
		return jets_cor_;
	}
	vector<float> &jets_emFrac()
	{
		if (not jets_emFrac_isLoaded) {
			if (jets_emFrac_branch != 0) 
				jets_emFrac_branch->GetEntry(index);
			jets_emFrac_isLoaded = true;
		}
		return jets_emFrac_;
	}
	vector<float> &mus_eledr()
	{
		if (not mus_eledr_isLoaded) {
			if (mus_eledr_branch != 0) 
				mus_eledr_branch->GetEntry(index);
			mus_eledr_isLoaded = true;
		}
		return mus_eledr_;
	}
	vector<float> &mus_jetdr()
	{
		if (not mus_jetdr_isLoaded) {
			if (mus_jetdr_branch != 0) 
				mus_jetdr_branch->GetEntry(index);
			mus_jetdr_isLoaded = true;
		}
		return mus_jetdr_;
	}
	vector<float> &mus_trkdr()
	{
		if (not mus_trkdr_isLoaded) {
			if (mus_trkdr_branch != 0) 
				mus_trkdr_branch->GetEntry(index);
			mus_trkdr_isLoaded = true;
		}
		return mus_trkdr_;
	}
	vector<float> &mus_chi2()
	{
		if (not mus_chi2_isLoaded) {
			if (mus_chi2_branch != 0) 
				mus_chi2_branch->GetEntry(index);
			mus_chi2_isLoaded = true;
		}
		return mus_chi2_;
	}
	vector<float> &mus_d0()
	{
		if (not mus_d0_isLoaded) {
			if (mus_d0_branch != 0) 
				mus_d0_branch->GetEntry(index);
			mus_d0_isLoaded = true;
		}
		return mus_d0_;
	}
	vector<float> &mus_d0Err()
	{
		if (not mus_d0Err_isLoaded) {
			if (mus_d0Err_branch != 0) 
				mus_d0Err_branch->GetEntry(index);
			mus_d0Err_isLoaded = true;
		}
		return mus_d0Err_;
	}
	vector<float> &mus_e_em()
	{
		if (not mus_e_em_isLoaded) {
			if (mus_e_em_branch != 0) 
				mus_e_em_branch->GetEntry(index);
			mus_e_em_isLoaded = true;
		}
		return mus_e_em_;
	}
	vector<float> &mus_e_emS9()
	{
		if (not mus_e_emS9_isLoaded) {
			if (mus_e_emS9_branch != 0) 
				mus_e_emS9_branch->GetEntry(index);
			mus_e_emS9_isLoaded = true;
		}
		return mus_e_emS9_;
	}
	vector<float> &mus_e_had()
	{
		if (not mus_e_had_isLoaded) {
			if (mus_e_had_branch != 0) 
				mus_e_had_branch->GetEntry(index);
			mus_e_had_isLoaded = true;
		}
		return mus_e_had_;
	}
	vector<float> &mus_e_hadS9()
	{
		if (not mus_e_hadS9_isLoaded) {
			if (mus_e_hadS9_branch != 0) 
				mus_e_hadS9_branch->GetEntry(index);
			mus_e_hadS9_isLoaded = true;
		}
		return mus_e_hadS9_;
	}
	vector<float> &mus_e_ho()
	{
		if (not mus_e_ho_isLoaded) {
			if (mus_e_ho_branch != 0) 
				mus_e_ho_branch->GetEntry(index);
			mus_e_ho_isLoaded = true;
		}
		return mus_e_ho_;
	}
	vector<float> &mus_e_hoS9()
	{
		if (not mus_e_hoS9_isLoaded) {
			if (mus_e_hoS9_branch != 0) 
				mus_e_hoS9_branch->GetEntry(index);
			mus_e_hoS9_isLoaded = true;
		}
		return mus_e_hoS9_;
	}
	vector<float> &mus_etaErr()
	{
		if (not mus_etaErr_isLoaded) {
			if (mus_etaErr_branch != 0) 
				mus_etaErr_branch->GetEntry(index);
			mus_etaErr_isLoaded = true;
		}
		return mus_etaErr_;
	}
	vector<float> &mus_gfit_chi2()
	{
		if (not mus_gfit_chi2_isLoaded) {
			if (mus_gfit_chi2_branch != 0) 
				mus_gfit_chi2_branch->GetEntry(index);
			mus_gfit_chi2_isLoaded = true;
		}
		return mus_gfit_chi2_;
	}
	vector<float> &mus_gfit_ndof()
	{
		if (not mus_gfit_ndof_isLoaded) {
			if (mus_gfit_ndof_branch != 0) 
				mus_gfit_ndof_branch->GetEntry(index);
			mus_gfit_ndof_isLoaded = true;
		}
		return mus_gfit_ndof_;
	}
	vector<float> &mus_iso()
	{
		if (not mus_iso_isLoaded) {
			if (mus_iso_branch != 0) 
				mus_iso_branch->GetEntry(index);
			mus_iso_isLoaded = true;
		}
		return mus_iso_;
	}
	vector<float> &mus_iso03_emEt()
	{
		if (not mus_iso03_emEt_isLoaded) {
			if (mus_iso03_emEt_branch != 0) 
				mus_iso03_emEt_branch->GetEntry(index);
			mus_iso03_emEt_isLoaded = true;
		}
		return mus_iso03_emEt_;
	}
	vector<float> &mus_iso03_hadEt()
	{
		if (not mus_iso03_hadEt_isLoaded) {
			if (mus_iso03_hadEt_branch != 0) 
				mus_iso03_hadEt_branch->GetEntry(index);
			mus_iso03_hadEt_isLoaded = true;
		}
		return mus_iso03_hadEt_;
	}
	vector<float> &mus_iso03_hoEt()
	{
		if (not mus_iso03_hoEt_isLoaded) {
			if (mus_iso03_hoEt_branch != 0) 
				mus_iso03_hoEt_branch->GetEntry(index);
			mus_iso03_hoEt_isLoaded = true;
		}
		return mus_iso03_hoEt_;
	}
	vector<float> &mus_iso03_sumPt()
	{
		if (not mus_iso03_sumPt_isLoaded) {
			if (mus_iso03_sumPt_branch != 0) 
				mus_iso03_sumPt_branch->GetEntry(index);
			mus_iso03_sumPt_isLoaded = true;
		}
		return mus_iso03_sumPt_;
	}
	vector<float> &mus_iso05_emEt()
	{
		if (not mus_iso05_emEt_isLoaded) {
			if (mus_iso05_emEt_branch != 0) 
				mus_iso05_emEt_branch->GetEntry(index);
			mus_iso05_emEt_isLoaded = true;
		}
		return mus_iso05_emEt_;
	}
	vector<float> &mus_iso05_hadEt()
	{
		if (not mus_iso05_hadEt_isLoaded) {
			if (mus_iso05_hadEt_branch != 0) 
				mus_iso05_hadEt_branch->GetEntry(index);
			mus_iso05_hadEt_isLoaded = true;
		}
		return mus_iso05_hadEt_;
	}
	vector<float> &mus_iso05_hoEt()
	{
		if (not mus_iso05_hoEt_isLoaded) {
			if (mus_iso05_hoEt_branch != 0) 
				mus_iso05_hoEt_branch->GetEntry(index);
			mus_iso05_hoEt_isLoaded = true;
		}
		return mus_iso05_hoEt_;
	}
	vector<float> &mus_iso05_sumPt()
	{
		if (not mus_iso05_sumPt_isLoaded) {
			if (mus_iso05_sumPt_branch != 0) 
				mus_iso05_sumPt_branch->GetEntry(index);
			mus_iso05_sumPt_isLoaded = true;
		}
		return mus_iso05_sumPt_;
	}
	vector<float> &mus_ndof()
	{
		if (not mus_ndof_isLoaded) {
			if (mus_ndof_branch != 0) 
				mus_ndof_branch->GetEntry(index);
			mus_ndof_isLoaded = true;
		}
		return mus_ndof_;
	}
	vector<float> &mus_outerEta()
	{
		if (not mus_outerEta_isLoaded) {
			if (mus_outerEta_branch != 0) 
				mus_outerEta_branch->GetEntry(index);
			mus_outerEta_isLoaded = true;
		}
		return mus_outerEta_;
	}
	vector<float> &mus_outerPhi()
	{
		if (not mus_outerPhi_isLoaded) {
			if (mus_outerPhi_branch != 0) 
				mus_outerPhi_branch->GetEntry(index);
			mus_outerPhi_isLoaded = true;
		}
		return mus_outerPhi_;
	}
	vector<float> &mus_phiErr()
	{
		if (not mus_phiErr_isLoaded) {
			if (mus_phiErr_branch != 0) 
				mus_phiErr_branch->GetEntry(index);
			mus_phiErr_isLoaded = true;
		}
		return mus_phiErr_;
	}
	vector<float> &mus_ptErr()
	{
		if (not mus_ptErr_isLoaded) {
			if (mus_ptErr_branch != 0) 
				mus_ptErr_branch->GetEntry(index);
			mus_ptErr_isLoaded = true;
		}
		return mus_ptErr_;
	}
	vector<float> &mus_vertexphi()
	{
		if (not mus_vertexphi_isLoaded) {
			if (mus_vertexphi_branch != 0) 
				mus_vertexphi_branch->GetEntry(index);
			mus_vertexphi_isLoaded = true;
		}
		return mus_vertexphi_;
	}
	vector<float> &mus_z0()
	{
		if (not mus_z0_isLoaded) {
			if (mus_z0_branch != 0) 
				mus_z0_branch->GetEntry(index);
			mus_z0_isLoaded = true;
		}
		return mus_z0_;
	}
	vector<float> &mus_z0Err()
	{
		if (not mus_z0Err_isLoaded) {
			if (mus_z0Err_branch != 0) 
				mus_z0Err_branch->GetEntry(index);
			mus_z0Err_isLoaded = true;
		}
		return mus_z0Err_;
	}
	vector<float> &els_tq_LRComb()
	{
		if (not els_tq_LRComb_isLoaded) {
			if (els_tq_LRComb_branch != 0) 
				els_tq_LRComb_branch->GetEntry(index);
			els_tq_LRComb_isLoaded = true;
		}
		return els_tq_LRComb_;
	}
	vector<float> &els_tq_caloIso()
	{
		if (not els_tq_caloIso_isLoaded) {
			if (els_tq_caloIso_branch != 0) 
				els_tq_caloIso_branch->GetEntry(index);
			els_tq_caloIso_isLoaded = true;
		}
		return els_tq_caloIso_;
	}
	vector<float> &els_tq_egammaEcalIso()
	{
		if (not els_tq_egammaEcalIso_isLoaded) {
			if (els_tq_egammaEcalIso_branch != 0) 
				els_tq_egammaEcalIso_branch->GetEntry(index);
			els_tq_egammaEcalIso_isLoaded = true;
		}
		return els_tq_egammaEcalIso_;
	}
	vector<float> &els_tq_egammaHcalIso()
	{
		if (not els_tq_egammaHcalIso_isLoaded) {
			if (els_tq_egammaHcalIso_branch != 0) 
				els_tq_egammaHcalIso_branch->GetEntry(index);
			els_tq_egammaHcalIso_isLoaded = true;
		}
		return els_tq_egammaHcalIso_;
	}
	vector<float> &els_tq_egammaTkIso()
	{
		if (not els_tq_egammaTkIso_isLoaded) {
			if (els_tq_egammaTkIso_branch != 0) 
				els_tq_egammaTkIso_branch->GetEntry(index);
			els_tq_egammaTkIso_isLoaded = true;
		}
		return els_tq_egammaTkIso_;
	}
	vector<float> &els_tq_electronIDRobust()
	{
		if (not els_tq_electronIDRobust_isLoaded) {
			if (els_tq_electronIDRobust_branch != 0) 
				els_tq_electronIDRobust_branch->GetEntry(index);
			els_tq_electronIDRobust_isLoaded = true;
		}
		return els_tq_electronIDRobust_;
	}
	vector<float> &els_tq_leptonID()
	{
		if (not els_tq_leptonID_isLoaded) {
			if (els_tq_leptonID_branch != 0) 
				els_tq_leptonID_branch->GetEntry(index);
			els_tq_leptonID_isLoaded = true;
		}
		return els_tq_leptonID_;
	}
	vector<float> &els_tq_trackIso()
	{
		if (not els_tq_trackIso_isLoaded) {
			if (els_tq_trackIso_branch != 0) 
				els_tq_trackIso_branch->GetEntry(index);
			els_tq_trackIso_isLoaded = true;
		}
		return els_tq_trackIso_;
	}
	vector<float> &jets_tq_bCorrF()
	{
		if (not jets_tq_bCorrF_isLoaded) {
			if (jets_tq_bCorrF_branch != 0) 
				jets_tq_bCorrF_branch->GetEntry(index);
			jets_tq_bCorrF_isLoaded = true;
		}
		return jets_tq_bCorrF_;
	}
	vector<float> &jets_tq_cCorrF()
	{
		if (not jets_tq_cCorrF_isLoaded) {
			if (jets_tq_cCorrF_branch != 0) 
				jets_tq_cCorrF_branch->GetEntry(index);
			jets_tq_cCorrF_isLoaded = true;
		}
		return jets_tq_cCorrF_;
	}
	vector<float> &jets_tq_gluCorrF()
	{
		if (not jets_tq_gluCorrF_isLoaded) {
			if (jets_tq_gluCorrF_branch != 0) 
				jets_tq_gluCorrF_branch->GetEntry(index);
			jets_tq_gluCorrF_isLoaded = true;
		}
		return jets_tq_gluCorrF_;
	}
	vector<float> &jets_tq_jetCharge()
	{
		if (not jets_tq_jetCharge_isLoaded) {
			if (jets_tq_jetCharge_branch != 0) 
				jets_tq_jetCharge_branch->GetEntry(index);
			jets_tq_jetCharge_isLoaded = true;
		}
		return jets_tq_jetCharge_;
	}
	vector<float> &jets_tq_noCorrF()
	{
		if (not jets_tq_noCorrF_isLoaded) {
			if (jets_tq_noCorrF_branch != 0) 
				jets_tq_noCorrF_branch->GetEntry(index);
			jets_tq_noCorrF_isLoaded = true;
		}
		return jets_tq_noCorrF_;
	}
	vector<float> &jets_tq_udsCorrF()
	{
		if (not jets_tq_udsCorrF_isLoaded) {
			if (jets_tq_udsCorrF_branch != 0) 
				jets_tq_udsCorrF_branch->GetEntry(index);
			jets_tq_udsCorrF_isLoaded = true;
		}
		return jets_tq_udsCorrF_;
	}
	vector<float> &mus_tq_caloIso()
	{
		if (not mus_tq_caloIso_isLoaded) {
			if (mus_tq_caloIso_branch != 0) 
				mus_tq_caloIso_branch->GetEntry(index);
			mus_tq_caloIso_isLoaded = true;
		}
		return mus_tq_caloIso_;
	}
	vector<float> &mus_tq_leptonID()
	{
		if (not mus_tq_leptonID_isLoaded) {
			if (mus_tq_leptonID_branch != 0) 
				mus_tq_leptonID_branch->GetEntry(index);
			mus_tq_leptonID_isLoaded = true;
		}
		return mus_tq_leptonID_;
	}
	vector<float> &mus_tq_lrComb()
	{
		if (not mus_tq_lrComb_isLoaded) {
			if (mus_tq_lrComb_branch != 0) 
				mus_tq_lrComb_branch->GetEntry(index);
			mus_tq_lrComb_isLoaded = true;
		}
		return mus_tq_lrComb_;
	}
	vector<float> &mus_tq_trackIso()
	{
		if (not mus_tq_trackIso_isLoaded) {
			if (mus_tq_trackIso_branch != 0) 
				mus_tq_trackIso_branch->GetEntry(index);
			mus_tq_trackIso_isLoaded = true;
		}
		return mus_tq_trackIso_;
	}
	vector<float> &trks_chi2()
	{
		if (not trks_chi2_isLoaded) {
			if (trks_chi2_branch != 0) 
				trks_chi2_branch->GetEntry(index);
			trks_chi2_isLoaded = true;
		}
		return trks_chi2_;
	}
	vector<float> &trks_d0()
	{
		if (not trks_d0_isLoaded) {
			if (trks_d0_branch != 0) 
				trks_d0_branch->GetEntry(index);
			trks_d0_isLoaded = true;
		}
		return trks_d0_;
	}
	vector<float> &trks_d0Err()
	{
		if (not trks_d0Err_isLoaded) {
			if (trks_d0Err_branch != 0) 
				trks_d0Err_branch->GetEntry(index);
			trks_d0Err_isLoaded = true;
		}
		return trks_d0Err_;
	}
	vector<float> &trks_etaErr()
	{
		if (not trks_etaErr_isLoaded) {
			if (trks_etaErr_branch != 0) 
				trks_etaErr_branch->GetEntry(index);
			trks_etaErr_isLoaded = true;
		}
		return trks_etaErr_;
	}
	vector<float> &trks_ndof()
	{
		if (not trks_ndof_isLoaded) {
			if (trks_ndof_branch != 0) 
				trks_ndof_branch->GetEntry(index);
			trks_ndof_isLoaded = true;
		}
		return trks_ndof_;
	}
	vector<float> &trks_outerEta()
	{
		if (not trks_outerEta_isLoaded) {
			if (trks_outerEta_branch != 0) 
				trks_outerEta_branch->GetEntry(index);
			trks_outerEta_isLoaded = true;
		}
		return trks_outerEta_;
	}
	vector<float> &trks_outerPhi()
	{
		if (not trks_outerPhi_isLoaded) {
			if (trks_outerPhi_branch != 0) 
				trks_outerPhi_branch->GetEntry(index);
			trks_outerPhi_isLoaded = true;
		}
		return trks_outerPhi_;
	}
	vector<float> &trks_phiErr()
	{
		if (not trks_phiErr_isLoaded) {
			if (trks_phiErr_branch != 0) 
				trks_phiErr_branch->GetEntry(index);
			trks_phiErr_isLoaded = true;
		}
		return trks_phiErr_;
	}
	vector<float> &trks_ptErr()
	{
		if (not trks_ptErr_isLoaded) {
			if (trks_ptErr_branch != 0) 
				trks_ptErr_branch->GetEntry(index);
			trks_ptErr_isLoaded = true;
		}
		return trks_ptErr_;
	}
	vector<float> &trks_vertexphi()
	{
		if (not trks_vertexphi_isLoaded) {
			if (trks_vertexphi_branch != 0) 
				trks_vertexphi_branch->GetEntry(index);
			trks_vertexphi_isLoaded = true;
		}
		return trks_vertexphi_;
	}
	vector<float> &trks_z0()
	{
		if (not trks_z0_isLoaded) {
			if (trks_z0_branch != 0) 
				trks_z0_branch->GetEntry(index);
			trks_z0_isLoaded = true;
		}
		return trks_z0_;
	}
	vector<float> &trks_z0Err()
	{
		if (not trks_z0Err_isLoaded) {
			if (trks_z0Err_branch != 0) 
				trks_z0Err_branch->GetEntry(index);
			trks_z0Err_isLoaded = true;
		}
		return trks_z0Err_;
	}
	vector<float> &trk_elsdr()
	{
		if (not trk_elsdr_isLoaded) {
			if (trk_elsdr_branch != 0) 
				trk_elsdr_branch->GetEntry(index);
			trk_elsdr_isLoaded = true;
		}
		return trk_elsdr_;
	}
	vector<float> &trk_musdr()
	{
		if (not trk_musdr_isLoaded) {
			if (trk_musdr_branch != 0) 
				trk_musdr_branch->GetEntry(index);
			trk_musdr_isLoaded = true;
		}
		return trk_musdr_;
	}
	vector<vector<float> > &hyp_jets_EMFcor()
	{
		if (not hyp_jets_EMFcor_isLoaded) {
			if (hyp_jets_EMFcor_branch != 0) 
				hyp_jets_EMFcor_branch->GetEntry(index);
			hyp_jets_EMFcor_isLoaded = true;
		}
		return hyp_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_jets_chFrac()
	{
		if (not hyp_jets_chFrac_isLoaded) {
			if (hyp_jets_chFrac_branch != 0) 
				hyp_jets_chFrac_branch->GetEntry(index);
			hyp_jets_chFrac_isLoaded = true;
		}
		return hyp_jets_chFrac_;
	}
	vector<vector<float> > &hyp_jets_cor()
	{
		if (not hyp_jets_cor_isLoaded) {
			if (hyp_jets_cor_branch != 0) 
				hyp_jets_cor_branch->GetEntry(index);
			hyp_jets_cor_isLoaded = true;
		}
		return hyp_jets_cor_;
	}
	vector<vector<float> > &hyp_jets_emFrac()
	{
		if (not hyp_jets_emFrac_isLoaded) {
			if (hyp_jets_emFrac_branch != 0) 
				hyp_jets_emFrac_branch->GetEntry(index);
			hyp_jets_emFrac_isLoaded = true;
		}
		return hyp_jets_emFrac_;
	}
	vector<vector<float> > &hyp_jets_mc_emEnergy()
	{
		if (not hyp_jets_mc_emEnergy_isLoaded) {
			if (hyp_jets_mc_emEnergy_branch != 0) 
				hyp_jets_mc_emEnergy_branch->GetEntry(index);
			hyp_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_hadEnergy()
	{
		if (not hyp_jets_mc_hadEnergy_isLoaded) {
			if (hyp_jets_mc_hadEnergy_branch != 0) 
				hyp_jets_mc_hadEnergy_branch->GetEntry(index);
			hyp_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_invEnergy()
	{
		if (not hyp_jets_mc_invEnergy_isLoaded) {
			if (hyp_jets_mc_invEnergy_branch != 0) 
				hyp_jets_mc_invEnergy_branch->GetEntry(index);
			hyp_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_otherEnergy()
	{
		if (not hyp_jets_mc_otherEnergy_isLoaded) {
			if (hyp_jets_mc_otherEnergy_branch != 0) 
				hyp_jets_mc_otherEnergy_branch->GetEntry(index);
			hyp_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_jets_tq_bCorrF()
	{
		if (not hyp_jets_tq_bCorrF_isLoaded) {
			if (hyp_jets_tq_bCorrF_branch != 0) 
				hyp_jets_tq_bCorrF_branch->GetEntry(index);
			hyp_jets_tq_bCorrF_isLoaded = true;
		}
		return hyp_jets_tq_bCorrF_;
	}
	vector<vector<float> > &hyp_jets_tq_cCorrF()
	{
		if (not hyp_jets_tq_cCorrF_isLoaded) {
			if (hyp_jets_tq_cCorrF_branch != 0) 
				hyp_jets_tq_cCorrF_branch->GetEntry(index);
			hyp_jets_tq_cCorrF_isLoaded = true;
		}
		return hyp_jets_tq_cCorrF_;
	}
	vector<vector<float> > &hyp_jets_tq_gluCorrF()
	{
		if (not hyp_jets_tq_gluCorrF_isLoaded) {
			if (hyp_jets_tq_gluCorrF_branch != 0) 
				hyp_jets_tq_gluCorrF_branch->GetEntry(index);
			hyp_jets_tq_gluCorrF_isLoaded = true;
		}
		return hyp_jets_tq_gluCorrF_;
	}
	vector<vector<float> > &hyp_jets_tq_jetCharge()
	{
		if (not hyp_jets_tq_jetCharge_isLoaded) {
			if (hyp_jets_tq_jetCharge_branch != 0) 
				hyp_jets_tq_jetCharge_branch->GetEntry(index);
			hyp_jets_tq_jetCharge_isLoaded = true;
		}
		return hyp_jets_tq_jetCharge_;
	}
	vector<vector<float> > &hyp_jets_tq_noCorrF()
	{
		if (not hyp_jets_tq_noCorrF_isLoaded) {
			if (hyp_jets_tq_noCorrF_branch != 0) 
				hyp_jets_tq_noCorrF_branch->GetEntry(index);
			hyp_jets_tq_noCorrF_isLoaded = true;
		}
		return hyp_jets_tq_noCorrF_;
	}
	vector<vector<float> > &hyp_jets_tq_udsCorrF()
	{
		if (not hyp_jets_tq_udsCorrF_isLoaded) {
			if (hyp_jets_tq_udsCorrF_branch != 0) 
				hyp_jets_tq_udsCorrF_branch->GetEntry(index);
			hyp_jets_tq_udsCorrF_isLoaded = true;
		}
		return hyp_jets_tq_udsCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_EMFcor()
	{
		if (not hyp_other_jets_EMFcor_isLoaded) {
			if (hyp_other_jets_EMFcor_branch != 0) 
				hyp_other_jets_EMFcor_branch->GetEntry(index);
			hyp_other_jets_EMFcor_isLoaded = true;
		}
		return hyp_other_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_other_jets_chFrac()
	{
		if (not hyp_other_jets_chFrac_isLoaded) {
			if (hyp_other_jets_chFrac_branch != 0) 
				hyp_other_jets_chFrac_branch->GetEntry(index);
			hyp_other_jets_chFrac_isLoaded = true;
		}
		return hyp_other_jets_chFrac_;
	}
	vector<vector<float> > &hyp_other_jets_cor()
	{
		if (not hyp_other_jets_cor_isLoaded) {
			if (hyp_other_jets_cor_branch != 0) 
				hyp_other_jets_cor_branch->GetEntry(index);
			hyp_other_jets_cor_isLoaded = true;
		}
		return hyp_other_jets_cor_;
	}
	vector<vector<float> > &hyp_other_jets_emFrac()
	{
		if (not hyp_other_jets_emFrac_isLoaded) {
			if (hyp_other_jets_emFrac_branch != 0) 
				hyp_other_jets_emFrac_branch->GetEntry(index);
			hyp_other_jets_emFrac_isLoaded = true;
		}
		return hyp_other_jets_emFrac_;
	}
	vector<vector<float> > &hyp_other_jets_mc_emEnergy()
	{
		if (not hyp_other_jets_mc_emEnergy_isLoaded) {
			if (hyp_other_jets_mc_emEnergy_branch != 0) 
				hyp_other_jets_mc_emEnergy_branch->GetEntry(index);
			hyp_other_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_hadEnergy()
	{
		if (not hyp_other_jets_mc_hadEnergy_isLoaded) {
			if (hyp_other_jets_mc_hadEnergy_branch != 0) 
				hyp_other_jets_mc_hadEnergy_branch->GetEntry(index);
			hyp_other_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_invEnergy()
	{
		if (not hyp_other_jets_mc_invEnergy_isLoaded) {
			if (hyp_other_jets_mc_invEnergy_branch != 0) 
				hyp_other_jets_mc_invEnergy_branch->GetEntry(index);
			hyp_other_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_otherEnergy()
	{
		if (not hyp_other_jets_mc_otherEnergy_isLoaded) {
			if (hyp_other_jets_mc_otherEnergy_branch != 0) 
				hyp_other_jets_mc_otherEnergy_branch->GetEntry(index);
			hyp_other_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_tq_bCorrF()
	{
		if (not hyp_other_jets_tq_bCorrF_isLoaded) {
			if (hyp_other_jets_tq_bCorrF_branch != 0) 
				hyp_other_jets_tq_bCorrF_branch->GetEntry(index);
			hyp_other_jets_tq_bCorrF_isLoaded = true;
		}
		return hyp_other_jets_tq_bCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_tq_cCorrF()
	{
		if (not hyp_other_jets_tq_cCorrF_isLoaded) {
			if (hyp_other_jets_tq_cCorrF_branch != 0) 
				hyp_other_jets_tq_cCorrF_branch->GetEntry(index);
			hyp_other_jets_tq_cCorrF_isLoaded = true;
		}
		return hyp_other_jets_tq_cCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_tq_gluCorrF()
	{
		if (not hyp_other_jets_tq_gluCorrF_isLoaded) {
			if (hyp_other_jets_tq_gluCorrF_branch != 0) 
				hyp_other_jets_tq_gluCorrF_branch->GetEntry(index);
			hyp_other_jets_tq_gluCorrF_isLoaded = true;
		}
		return hyp_other_jets_tq_gluCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_tq_jetCharge()
	{
		if (not hyp_other_jets_tq_jetCharge_isLoaded) {
			if (hyp_other_jets_tq_jetCharge_branch != 0) 
				hyp_other_jets_tq_jetCharge_branch->GetEntry(index);
			hyp_other_jets_tq_jetCharge_isLoaded = true;
		}
		return hyp_other_jets_tq_jetCharge_;
	}
	vector<vector<float> > &hyp_other_jets_tq_noCorrF()
	{
		if (not hyp_other_jets_tq_noCorrF_isLoaded) {
			if (hyp_other_jets_tq_noCorrF_branch != 0) 
				hyp_other_jets_tq_noCorrF_branch->GetEntry(index);
			hyp_other_jets_tq_noCorrF_isLoaded = true;
		}
		return hyp_other_jets_tq_noCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_tq_udsCorrF()
	{
		if (not hyp_other_jets_tq_udsCorrF_isLoaded) {
			if (hyp_other_jets_tq_udsCorrF_branch != 0) 
				hyp_other_jets_tq_udsCorrF_branch->GetEntry(index);
			hyp_other_jets_tq_udsCorrF_isLoaded = true;
		}
		return hyp_other_jets_tq_udsCorrF_;
	}
	int &evt_CSA07Process()
	{
		if (not evt_CSA07Process_isLoaded) {
			if (evt_CSA07Process_branch != 0) 
				evt_CSA07Process_branch->GetEntry(index);
			evt_CSA07Process_isLoaded = true;
		}
		return evt_CSA07Process_;
	}
	int &evt_HLT1()
	{
		if (not evt_HLT1_isLoaded) {
			if (evt_HLT1_branch != 0) 
				evt_HLT1_branch->GetEntry(index);
			evt_HLT1_isLoaded = true;
		}
		return evt_HLT1_;
	}
	int &evt_HLT2()
	{
		if (not evt_HLT2_isLoaded) {
			if (evt_HLT2_branch != 0) 
				evt_HLT2_branch->GetEntry(index);
			evt_HLT2_isLoaded = true;
		}
		return evt_HLT2_;
	}
	int &evt_HLT3()
	{
		if (not evt_HLT3_isLoaded) {
			if (evt_HLT3_branch != 0) 
				evt_HLT3_branch->GetEntry(index);
			evt_HLT3_isLoaded = true;
		}
		return evt_HLT3_;
	}
	int &evt_HLT4()
	{
		if (not evt_HLT4_isLoaded) {
			if (evt_HLT4_branch != 0) 
				evt_HLT4_branch->GetEntry(index);
			evt_HLT4_isLoaded = true;
		}
		return evt_HLT4_;
	}
	int &evt_L1_1()
	{
		if (not evt_L1_1_isLoaded) {
			if (evt_L1_1_branch != 0) 
				evt_L1_1_branch->GetEntry(index);
			evt_L1_1_isLoaded = true;
		}
		return evt_L1_1_;
	}
	int &evt_L1_2()
	{
		if (not evt_L1_2_isLoaded) {
			if (evt_L1_2_branch != 0) 
				evt_L1_2_branch->GetEntry(index);
			evt_L1_2_isLoaded = true;
		}
		return evt_L1_2_;
	}
	int &evt_L1_3()
	{
		if (not evt_L1_3_isLoaded) {
			if (evt_L1_3_branch != 0) 
				evt_L1_3_branch->GetEntry(index);
			evt_L1_3_isLoaded = true;
		}
		return evt_L1_3_;
	}
	int &evt_L1_4()
	{
		if (not evt_L1_4_isLoaded) {
			if (evt_L1_4_branch != 0) 
				evt_L1_4_branch->GetEntry(index);
			evt_L1_4_isLoaded = true;
		}
		return evt_L1_4_;
	}
	int &evt_nl1emiso()
	{
		if (not evt_nl1emiso_isLoaded) {
			if (evt_nl1emiso_branch != 0) 
				evt_nl1emiso_branch->GetEntry(index);
			evt_nl1emiso_isLoaded = true;
		}
		return evt_nl1emiso_;
	}
	int &evt_nl1emnoiso()
	{
		if (not evt_nl1emnoiso_isLoaded) {
			if (evt_nl1emnoiso_branch != 0) 
				evt_nl1emnoiso_branch->GetEntry(index);
			evt_nl1emnoiso_isLoaded = true;
		}
		return evt_nl1emnoiso_;
	}
	int &evt_nl1jetsc()
	{
		if (not evt_nl1jetsc_isLoaded) {
			if (evt_nl1jetsc_branch != 0) 
				evt_nl1jetsc_branch->GetEntry(index);
			evt_nl1jetsc_isLoaded = true;
		}
		return evt_nl1jetsc_;
	}
	int &evt_nl1jetsf()
	{
		if (not evt_nl1jetsf_isLoaded) {
			if (evt_nl1jetsf_branch != 0) 
				evt_nl1jetsf_branch->GetEntry(index);
			evt_nl1jetsf_isLoaded = true;
		}
		return evt_nl1jetsf_;
	}
	int &evt_nl1jetst()
	{
		if (not evt_nl1jetst_isLoaded) {
			if (evt_nl1jetst_branch != 0) 
				evt_nl1jetst_branch->GetEntry(index);
			evt_nl1jetst_isLoaded = true;
		}
		return evt_nl1jetst_;
	}
	int &evt_nl1mus()
	{
		if (not evt_nl1mus_isLoaded) {
			if (evt_nl1mus_branch != 0) 
				evt_nl1mus_branch->GetEntry(index);
			evt_nl1mus_isLoaded = true;
		}
		return evt_nl1mus_;
	}
	vector<int> &els_closestMuon()
	{
		if (not els_closestMuon_isLoaded) {
			if (els_closestMuon_branch != 0) 
				els_closestMuon_branch->GetEntry(index);
			els_closestMuon_isLoaded = true;
		}
		return els_closestMuon_;
	}
	vector<int> &els_trkidx()
	{
		if (not els_trkidx_isLoaded) {
			if (els_trkidx_branch != 0) 
				els_trkidx_branch->GetEntry(index);
			els_trkidx_isLoaded = true;
		}
		return els_trkidx_;
	}
	vector<int> &els_charge()
	{
		if (not els_charge_isLoaded) {
			if (els_charge_branch != 0) 
				els_charge_branch->GetEntry(index);
			els_charge_isLoaded = true;
		}
		return els_charge_;
	}
	vector<int> &els_class()
	{
		if (not els_class_isLoaded) {
			if (els_class_branch != 0) 
				els_class_branch->GetEntry(index);
			els_class_isLoaded = true;
		}
		return els_class_;
	}
	vector<int> &els_looseId()
	{
		if (not els_looseId_isLoaded) {
			if (els_looseId_branch != 0) 
				els_looseId_branch->GetEntry(index);
			els_looseId_isLoaded = true;
		}
		return els_looseId_;
	}
	vector<int> &els_lostHits()
	{
		if (not els_lostHits_isLoaded) {
			if (els_lostHits_branch != 0) 
				els_lostHits_branch->GetEntry(index);
			els_lostHits_isLoaded = true;
		}
		return els_lostHits_;
	}
	vector<int> &els_nSeed()
	{
		if (not els_nSeed_isLoaded) {
			if (els_nSeed_branch != 0) 
				els_nSeed_branch->GetEntry(index);
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	vector<int> &els_pass3looseId()
	{
		if (not els_pass3looseId_isLoaded) {
			if (els_pass3looseId_branch != 0) 
				els_pass3looseId_branch->GetEntry(index);
			els_pass3looseId_isLoaded = true;
		}
		return els_pass3looseId_;
	}
	vector<int> &els_pass3simpleId()
	{
		if (not els_pass3simpleId_isLoaded) {
			if (els_pass3simpleId_branch != 0) 
				els_pass3simpleId_branch->GetEntry(index);
			els_pass3simpleId_isLoaded = true;
		}
		return els_pass3simpleId_;
	}
	vector<int> &els_pass3tightId()
	{
		if (not els_pass3tightId_isLoaded) {
			if (els_pass3tightId_branch != 0) 
				els_pass3tightId_branch->GetEntry(index);
			els_pass3tightId_isLoaded = true;
		}
		return els_pass3tightId_;
	}
	vector<int> &els_robustId()
	{
		if (not els_robustId_isLoaded) {
			if (els_robustId_branch != 0) 
				els_robustId_branch->GetEntry(index);
			els_robustId_isLoaded = true;
		}
		return els_robustId_;
	}
	vector<int> &els_simpleIdPlus()
	{
		if (not els_simpleIdPlus_isLoaded) {
			if (els_simpleIdPlus_branch != 0) 
				els_simpleIdPlus_branch->GetEntry(index);
			els_simpleIdPlus_isLoaded = true;
		}
		return els_simpleIdPlus_;
	}
	vector<int> &els_tightId()
	{
		if (not els_tightId_isLoaded) {
			if (els_tightId_branch != 0) 
				els_tightId_branch->GetEntry(index);
			els_tightId_isLoaded = true;
		}
		return els_tightId_;
	}
	vector<int> &els_validHits()
	{
		if (not els_validHits_isLoaded) {
			if (els_validHits_branch != 0) 
				els_validHits_branch->GetEntry(index);
			els_validHits_isLoaded = true;
		}
		return els_validHits_;
	}
	vector<int> &genps_id()
	{
		if (not genps_id_isLoaded) {
			if (genps_id_branch != 0) 
				genps_id_branch->GetEntry(index);
			genps_id_isLoaded = true;
		}
		return genps_id_;
	}
	vector<int> &genps_id_mother()
	{
		if (not genps_id_mother_isLoaded) {
			if (genps_id_mother_branch != 0) 
				genps_id_mother_branch->GetEntry(index);
			genps_id_mother_isLoaded = true;
		}
		return genps_id_mother_;
	}
	vector<int> &genps_status()
	{
		if (not genps_status_isLoaded) {
			if (genps_status_branch != 0) 
				genps_status_branch->GetEntry(index);
			genps_status_isLoaded = true;
		}
		return genps_status_;
	}
	vector<int> &hyp_ll_charge()
	{
		if (not hyp_ll_charge_isLoaded) {
			if (hyp_ll_charge_branch != 0) 
				hyp_ll_charge_branch->GetEntry(index);
			hyp_ll_charge_isLoaded = true;
		}
		return hyp_ll_charge_;
	}
	vector<int> &hyp_ll_id()
	{
		if (not hyp_ll_id_isLoaded) {
			if (hyp_ll_id_branch != 0) 
				hyp_ll_id_branch->GetEntry(index);
			hyp_ll_id_isLoaded = true;
		}
		return hyp_ll_id_;
	}
	vector<int> &hyp_ll_index()
	{
		if (not hyp_ll_index_isLoaded) {
			if (hyp_ll_index_branch != 0) 
				hyp_ll_index_branch->GetEntry(index);
			hyp_ll_index_isLoaded = true;
		}
		return hyp_ll_index_;
	}
	vector<int> &hyp_ll_lostHits()
	{
		if (not hyp_ll_lostHits_isLoaded) {
			if (hyp_ll_lostHits_branch != 0) 
				hyp_ll_lostHits_branch->GetEntry(index);
			hyp_ll_lostHits_isLoaded = true;
		}
		return hyp_ll_lostHits_;
	}
	vector<int> &hyp_ll_mc_id()
	{
		if (not hyp_ll_mc_id_isLoaded) {
			if (hyp_ll_mc_id_branch != 0) 
				hyp_ll_mc_id_branch->GetEntry(index);
			hyp_ll_mc_id_isLoaded = true;
		}
		return hyp_ll_mc_id_;
	}
	vector<int> &hyp_ll_mc_motherid()
	{
		if (not hyp_ll_mc_motherid_isLoaded) {
			if (hyp_ll_mc_motherid_branch != 0) 
				hyp_ll_mc_motherid_branch->GetEntry(index);
			hyp_ll_mc_motherid_isLoaded = true;
		}
		return hyp_ll_mc_motherid_;
	}
	vector<int> &hyp_ll_validHits()
	{
		if (not hyp_ll_validHits_isLoaded) {
			if (hyp_ll_validHits_branch != 0) 
				hyp_ll_validHits_branch->GetEntry(index);
			hyp_ll_validHits_isLoaded = true;
		}
		return hyp_ll_validHits_;
	}
	vector<int> &hyp_lt_charge()
	{
		if (not hyp_lt_charge_isLoaded) {
			if (hyp_lt_charge_branch != 0) 
				hyp_lt_charge_branch->GetEntry(index);
			hyp_lt_charge_isLoaded = true;
		}
		return hyp_lt_charge_;
	}
	vector<int> &hyp_lt_id()
	{
		if (not hyp_lt_id_isLoaded) {
			if (hyp_lt_id_branch != 0) 
				hyp_lt_id_branch->GetEntry(index);
			hyp_lt_id_isLoaded = true;
		}
		return hyp_lt_id_;
	}
	vector<int> &hyp_lt_index()
	{
		if (not hyp_lt_index_isLoaded) {
			if (hyp_lt_index_branch != 0) 
				hyp_lt_index_branch->GetEntry(index);
			hyp_lt_index_isLoaded = true;
		}
		return hyp_lt_index_;
	}
	vector<int> &hyp_lt_lostHits()
	{
		if (not hyp_lt_lostHits_isLoaded) {
			if (hyp_lt_lostHits_branch != 0) 
				hyp_lt_lostHits_branch->GetEntry(index);
			hyp_lt_lostHits_isLoaded = true;
		}
		return hyp_lt_lostHits_;
	}
	vector<int> &hyp_lt_mc_id()
	{
		if (not hyp_lt_mc_id_isLoaded) {
			if (hyp_lt_mc_id_branch != 0) 
				hyp_lt_mc_id_branch->GetEntry(index);
			hyp_lt_mc_id_isLoaded = true;
		}
		return hyp_lt_mc_id_;
	}
	vector<int> &hyp_lt_mc_motherid()
	{
		if (not hyp_lt_mc_motherid_isLoaded) {
			if (hyp_lt_mc_motherid_branch != 0) 
				hyp_lt_mc_motherid_branch->GetEntry(index);
			hyp_lt_mc_motherid_isLoaded = true;
		}
		return hyp_lt_mc_motherid_;
	}
	vector<int> &hyp_lt_validHits()
	{
		if (not hyp_lt_validHits_isLoaded) {
			if (hyp_lt_validHits_branch != 0) 
				hyp_lt_validHits_branch->GetEntry(index);
			hyp_lt_validHits_isLoaded = true;
		}
		return hyp_lt_validHits_;
	}
	vector<int> &hyp_njets()
	{
		if (not hyp_njets_isLoaded) {
			if (hyp_njets_branch != 0) 
				hyp_njets_branch->GetEntry(index);
			hyp_njets_isLoaded = true;
		}
		return hyp_njets_;
	}
	vector<int> &hyp_nojets()
	{
		if (not hyp_nojets_isLoaded) {
			if (hyp_nojets_branch != 0) 
				hyp_nojets_branch->GetEntry(index);
			hyp_nojets_isLoaded = true;
		}
		return hyp_nojets_;
	}
	vector<int> &hyp_type()
	{
		if (not hyp_type_isLoaded) {
			if (hyp_type_branch != 0) 
				hyp_type_branch->GetEntry(index);
			hyp_type_isLoaded = true;
		}
		return hyp_type_;
	}
	vector<int> &hyp_quadlep_first_type()
	{
		if (not hyp_quadlep_first_type_isLoaded) {
			if (hyp_quadlep_first_type_branch != 0) 
				hyp_quadlep_first_type_branch->GetEntry(index);
			hyp_quadlep_first_type_isLoaded = true;
		}
		return hyp_quadlep_first_type_;
	}
	vector<int> &hyp_quadlep_fourth_type()
	{
		if (not hyp_quadlep_fourth_type_isLoaded) {
			if (hyp_quadlep_fourth_type_branch != 0) 
				hyp_quadlep_fourth_type_branch->GetEntry(index);
			hyp_quadlep_fourth_type_isLoaded = true;
		}
		return hyp_quadlep_fourth_type_;
	}
	vector<int> &hyp_quadlep_second_type()
	{
		if (not hyp_quadlep_second_type_isLoaded) {
			if (hyp_quadlep_second_type_branch != 0) 
				hyp_quadlep_second_type_branch->GetEntry(index);
			hyp_quadlep_second_type_isLoaded = true;
		}
		return hyp_quadlep_second_type_;
	}
	vector<int> &hyp_quadlep_third_type()
	{
		if (not hyp_quadlep_third_type_isLoaded) {
			if (hyp_quadlep_third_type_branch != 0) 
				hyp_quadlep_third_type_branch->GetEntry(index);
			hyp_quadlep_third_type_isLoaded = true;
		}
		return hyp_quadlep_third_type_;
	}
	vector<int> &hyp_trilep_first_type()
	{
		if (not hyp_trilep_first_type_isLoaded) {
			if (hyp_trilep_first_type_branch != 0) 
				hyp_trilep_first_type_branch->GetEntry(index);
			hyp_trilep_first_type_isLoaded = true;
		}
		return hyp_trilep_first_type_;
	}
	vector<int> &hyp_trilep_second_type()
	{
		if (not hyp_trilep_second_type_isLoaded) {
			if (hyp_trilep_second_type_branch != 0) 
				hyp_trilep_second_type_branch->GetEntry(index);
			hyp_trilep_second_type_isLoaded = true;
		}
		return hyp_trilep_second_type_;
	}
	vector<int> &hyp_trilep_third_type()
	{
		if (not hyp_trilep_third_type_isLoaded) {
			if (hyp_trilep_third_type_branch != 0) 
				hyp_trilep_third_type_branch->GetEntry(index);
			hyp_trilep_third_type_isLoaded = true;
		}
		return hyp_trilep_third_type_;
	}
	vector<int> &jets_closestElectron()
	{
		if (not jets_closestElectron_isLoaded) {
			if (jets_closestElectron_branch != 0) 
				jets_closestElectron_branch->GetEntry(index);
			jets_closestElectron_isLoaded = true;
		}
		return jets_closestElectron_;
	}
	vector<int> &jets_closestMuon()
	{
		if (not jets_closestMuon_isLoaded) {
			if (jets_closestMuon_branch != 0) 
				jets_closestMuon_branch->GetEntry(index);
			jets_closestMuon_isLoaded = true;
		}
		return jets_closestMuon_;
	}
	vector<int> &l1emiso_ieta()
	{
		if (not l1emiso_ieta_isLoaded) {
			if (l1emiso_ieta_branch != 0) 
				l1emiso_ieta_branch->GetEntry(index);
			l1emiso_ieta_isLoaded = true;
		}
		return l1emiso_ieta_;
	}
	vector<int> &l1emiso_iphi()
	{
		if (not l1emiso_iphi_isLoaded) {
			if (l1emiso_iphi_branch != 0) 
				l1emiso_iphi_branch->GetEntry(index);
			l1emiso_iphi_isLoaded = true;
		}
		return l1emiso_iphi_;
	}
	vector<int> &l1emiso_rawId()
	{
		if (not l1emiso_rawId_isLoaded) {
			if (l1emiso_rawId_branch != 0) 
				l1emiso_rawId_branch->GetEntry(index);
			l1emiso_rawId_isLoaded = true;
		}
		return l1emiso_rawId_;
	}
	vector<int> &l1emiso_type()
	{
		if (not l1emiso_type_isLoaded) {
			if (l1emiso_type_branch != 0) 
				l1emiso_type_branch->GetEntry(index);
			l1emiso_type_isLoaded = true;
		}
		return l1emiso_type_;
	}
	vector<int> &l1emnoiso_ieta()
	{
		if (not l1emnoiso_ieta_isLoaded) {
			if (l1emnoiso_ieta_branch != 0) 
				l1emnoiso_ieta_branch->GetEntry(index);
			l1emnoiso_ieta_isLoaded = true;
		}
		return l1emnoiso_ieta_;
	}
	vector<int> &l1emnoiso_iphi()
	{
		if (not l1emnoiso_iphi_isLoaded) {
			if (l1emnoiso_iphi_branch != 0) 
				l1emnoiso_iphi_branch->GetEntry(index);
			l1emnoiso_iphi_isLoaded = true;
		}
		return l1emnoiso_iphi_;
	}
	vector<int> &l1emnoiso_rawId()
	{
		if (not l1emnoiso_rawId_isLoaded) {
			if (l1emnoiso_rawId_branch != 0) 
				l1emnoiso_rawId_branch->GetEntry(index);
			l1emnoiso_rawId_isLoaded = true;
		}
		return l1emnoiso_rawId_;
	}
	vector<int> &l1emnoiso_type()
	{
		if (not l1emnoiso_type_isLoaded) {
			if (l1emnoiso_type_branch != 0) 
				l1emnoiso_type_branch->GetEntry(index);
			l1emnoiso_type_isLoaded = true;
		}
		return l1emnoiso_type_;
	}
	vector<int> &l1jetsc_ieta()
	{
		if (not l1jetsc_ieta_isLoaded) {
			if (l1jetsc_ieta_branch != 0) 
				l1jetsc_ieta_branch->GetEntry(index);
			l1jetsc_ieta_isLoaded = true;
		}
		return l1jetsc_ieta_;
	}
	vector<int> &l1jetsc_iphi()
	{
		if (not l1jetsc_iphi_isLoaded) {
			if (l1jetsc_iphi_branch != 0) 
				l1jetsc_iphi_branch->GetEntry(index);
			l1jetsc_iphi_isLoaded = true;
		}
		return l1jetsc_iphi_;
	}
	vector<int> &l1jetsc_rawId()
	{
		if (not l1jetsc_rawId_isLoaded) {
			if (l1jetsc_rawId_branch != 0) 
				l1jetsc_rawId_branch->GetEntry(index);
			l1jetsc_rawId_isLoaded = true;
		}
		return l1jetsc_rawId_;
	}
	vector<int> &l1jetsc_type()
	{
		if (not l1jetsc_type_isLoaded) {
			if (l1jetsc_type_branch != 0) 
				l1jetsc_type_branch->GetEntry(index);
			l1jetsc_type_isLoaded = true;
		}
		return l1jetsc_type_;
	}
	vector<int> &l1jetsf_ieta()
	{
		if (not l1jetsf_ieta_isLoaded) {
			if (l1jetsf_ieta_branch != 0) 
				l1jetsf_ieta_branch->GetEntry(index);
			l1jetsf_ieta_isLoaded = true;
		}
		return l1jetsf_ieta_;
	}
	vector<int> &l1jetsf_iphi()
	{
		if (not l1jetsf_iphi_isLoaded) {
			if (l1jetsf_iphi_branch != 0) 
				l1jetsf_iphi_branch->GetEntry(index);
			l1jetsf_iphi_isLoaded = true;
		}
		return l1jetsf_iphi_;
	}
	vector<int> &l1jetsf_rawId()
	{
		if (not l1jetsf_rawId_isLoaded) {
			if (l1jetsf_rawId_branch != 0) 
				l1jetsf_rawId_branch->GetEntry(index);
			l1jetsf_rawId_isLoaded = true;
		}
		return l1jetsf_rawId_;
	}
	vector<int> &l1jetsf_type()
	{
		if (not l1jetsf_type_isLoaded) {
			if (l1jetsf_type_branch != 0) 
				l1jetsf_type_branch->GetEntry(index);
			l1jetsf_type_isLoaded = true;
		}
		return l1jetsf_type_;
	}
	vector<int> &l1jetst_ieta()
	{
		if (not l1jetst_ieta_isLoaded) {
			if (l1jetst_ieta_branch != 0) 
				l1jetst_ieta_branch->GetEntry(index);
			l1jetst_ieta_isLoaded = true;
		}
		return l1jetst_ieta_;
	}
	vector<int> &l1jetst_iphi()
	{
		if (not l1jetst_iphi_isLoaded) {
			if (l1jetst_iphi_branch != 0) 
				l1jetst_iphi_branch->GetEntry(index);
			l1jetst_iphi_isLoaded = true;
		}
		return l1jetst_iphi_;
	}
	vector<int> &l1jetst_rawId()
	{
		if (not l1jetst_rawId_isLoaded) {
			if (l1jetst_rawId_branch != 0) 
				l1jetst_rawId_branch->GetEntry(index);
			l1jetst_rawId_isLoaded = true;
		}
		return l1jetst_rawId_;
	}
	vector<int> &l1jetst_type()
	{
		if (not l1jetst_type_isLoaded) {
			if (l1jetst_type_branch != 0) 
				l1jetst_type_branch->GetEntry(index);
			l1jetst_type_isLoaded = true;
		}
		return l1jetst_type_;
	}
	vector<int> &l1mus_flags()
	{
		if (not l1mus_flags_isLoaded) {
			if (l1mus_flags_branch != 0) 
				l1mus_flags_branch->GetEntry(index);
			l1mus_flags_isLoaded = true;
		}
		return l1mus_flags_;
	}
	vector<int> &l1mus_q()
	{
		if (not l1mus_q_isLoaded) {
			if (l1mus_q_branch != 0) 
				l1mus_q_branch->GetEntry(index);
			l1mus_q_isLoaded = true;
		}
		return l1mus_q_;
	}
	vector<int> &l1mus_qual()
	{
		if (not l1mus_qual_isLoaded) {
			if (l1mus_qual_branch != 0) 
				l1mus_qual_branch->GetEntry(index);
			l1mus_qual_isLoaded = true;
		}
		return l1mus_qual_;
	}
	vector<int> &l1mus_qualFlags()
	{
		if (not l1mus_qualFlags_isLoaded) {
			if (l1mus_qualFlags_branch != 0) 
				l1mus_qualFlags_branch->GetEntry(index);
			l1mus_qualFlags_isLoaded = true;
		}
		return l1mus_qualFlags_;
	}
	vector<int> &mus_closestEle()
	{
		if (not mus_closestEle_isLoaded) {
			if (mus_closestEle_branch != 0) 
				mus_closestEle_branch->GetEntry(index);
			mus_closestEle_isLoaded = true;
		}
		return mus_closestEle_;
	}
	vector<int> &mus_closestJet()
	{
		if (not mus_closestJet_isLoaded) {
			if (mus_closestJet_branch != 0) 
				mus_closestJet_branch->GetEntry(index);
			mus_closestJet_isLoaded = true;
		}
		return mus_closestJet_;
	}
	vector<int> &mus_trkidx()
	{
		if (not mus_trkidx_isLoaded) {
			if (mus_trkidx_branch != 0) 
				mus_trkidx_branch->GetEntry(index);
			mus_trkidx_isLoaded = true;
		}
		return mus_trkidx_;
	}
	vector<int> &mus_charge()
	{
		if (not mus_charge_isLoaded) {
			if (mus_charge_branch != 0) 
				mus_charge_branch->GetEntry(index);
			mus_charge_isLoaded = true;
		}
		return mus_charge_;
	}
	vector<int> &mus_gfit_validHits()
	{
		if (not mus_gfit_validHits_isLoaded) {
			if (mus_gfit_validHits_branch != 0) 
				mus_gfit_validHits_branch->GetEntry(index);
			mus_gfit_validHits_isLoaded = true;
		}
		return mus_gfit_validHits_;
	}
	vector<int> &mus_iso03_ntrk()
	{
		if (not mus_iso03_ntrk_isLoaded) {
			if (mus_iso03_ntrk_branch != 0) 
				mus_iso03_ntrk_branch->GetEntry(index);
			mus_iso03_ntrk_isLoaded = true;
		}
		return mus_iso03_ntrk_;
	}
	vector<int> &mus_iso05_ntrk()
	{
		if (not mus_iso05_ntrk_isLoaded) {
			if (mus_iso05_ntrk_branch != 0) 
				mus_iso05_ntrk_branch->GetEntry(index);
			mus_iso05_ntrk_isLoaded = true;
		}
		return mus_iso05_ntrk_;
	}
	vector<int> &mus_lostHits()
	{
		if (not mus_lostHits_isLoaded) {
			if (mus_lostHits_branch != 0) 
				mus_lostHits_branch->GetEntry(index);
			mus_lostHits_isLoaded = true;
		}
		return mus_lostHits_;
	}
	vector<int> &mus_nmatches()
	{
		if (not mus_nmatches_isLoaded) {
			if (mus_nmatches_branch != 0) 
				mus_nmatches_branch->GetEntry(index);
			mus_nmatches_isLoaded = true;
		}
		return mus_nmatches_;
	}
	vector<int> &mus_pid_TM2DCompatibilityLoose()
	{
		if (not mus_pid_TM2DCompatibilityLoose_isLoaded) {
			if (mus_pid_TM2DCompatibilityLoose_branch != 0) 
				mus_pid_TM2DCompatibilityLoose_branch->GetEntry(index);
			mus_pid_TM2DCompatibilityLoose_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityLoose_;
	}
	vector<int> &mus_pid_TM2DCompatibilityTight()
	{
		if (not mus_pid_TM2DCompatibilityTight_isLoaded) {
			if (mus_pid_TM2DCompatibilityTight_branch != 0) 
				mus_pid_TM2DCompatibilityTight_branch->GetEntry(index);
			mus_pid_TM2DCompatibilityTight_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityTight_;
	}
	vector<int> &mus_pid_TMLastStationLoose()
	{
		if (not mus_pid_TMLastStationLoose_isLoaded) {
			if (mus_pid_TMLastStationLoose_branch != 0) 
				mus_pid_TMLastStationLoose_branch->GetEntry(index);
			mus_pid_TMLastStationLoose_isLoaded = true;
		}
		return mus_pid_TMLastStationLoose_;
	}
	vector<int> &mus_pid_TMLastStationTight()
	{
		if (not mus_pid_TMLastStationTight_isLoaded) {
			if (mus_pid_TMLastStationTight_branch != 0) 
				mus_pid_TMLastStationTight_branch->GetEntry(index);
			mus_pid_TMLastStationTight_isLoaded = true;
		}
		return mus_pid_TMLastStationTight_;
	}
	vector<int> &mus_trkrefkey()
	{
		if (not mus_trkrefkey_isLoaded) {
			if (mus_trkrefkey_branch != 0) 
				mus_trkrefkey_branch->GetEntry(index);
			mus_trkrefkey_isLoaded = true;
		}
		return mus_trkrefkey_;
	}
	vector<int> &mus_validHits()
	{
		if (not mus_validHits_isLoaded) {
			if (mus_validHits_branch != 0) 
				mus_validHits_branch->GetEntry(index);
			mus_validHits_isLoaded = true;
		}
		return mus_validHits_;
	}
	vector<int> &els_tq_egammaTkNumIso()
	{
		if (not els_tq_egammaTkNumIso_isLoaded) {
			if (els_tq_egammaTkNumIso_branch != 0) 
				els_tq_egammaTkNumIso_branch->GetEntry(index);
			els_tq_egammaTkNumIso_isLoaded = true;
		}
		return els_tq_egammaTkNumIso_;
	}
	vector<int> &els_tq_genID()
	{
		if (not els_tq_genID_isLoaded) {
			if (els_tq_genID_branch != 0) 
				els_tq_genID_branch->GetEntry(index);
			els_tq_genID_isLoaded = true;
		}
		return els_tq_genID_;
	}
	vector<int> &els_tq_genMotherID()
	{
		if (not els_tq_genMotherID_isLoaded) {
			if (els_tq_genMotherID_branch != 0) 
				els_tq_genMotherID_branch->GetEntry(index);
			els_tq_genMotherID_isLoaded = true;
		}
		return els_tq_genMotherID_;
	}
	vector<int> &jets_tq_genPartonMother_id()
	{
		if (not jets_tq_genPartonMother_id_isLoaded) {
			if (jets_tq_genPartonMother_id_branch != 0) 
				jets_tq_genPartonMother_id_branch->GetEntry(index);
			jets_tq_genPartonMother_id_isLoaded = true;
		}
		return jets_tq_genPartonMother_id_;
	}
	vector<int> &jets_tq_genParton_id()
	{
		if (not jets_tq_genParton_id_isLoaded) {
			if (jets_tq_genParton_id_branch != 0) 
				jets_tq_genParton_id_branch->GetEntry(index);
			jets_tq_genParton_id_isLoaded = true;
		}
		return jets_tq_genParton_id_;
	}
	vector<int> &jets_tq_partonFlavour()
	{
		if (not jets_tq_partonFlavour_isLoaded) {
			if (jets_tq_partonFlavour_branch != 0) 
				jets_tq_partonFlavour_branch->GetEntry(index);
			jets_tq_partonFlavour_isLoaded = true;
		}
		return jets_tq_partonFlavour_;
	}
	vector<int> &mus_tq_genID()
	{
		if (not mus_tq_genID_isLoaded) {
			if (mus_tq_genID_branch != 0) 
				mus_tq_genID_branch->GetEntry(index);
			mus_tq_genID_isLoaded = true;
		}
		return mus_tq_genID_;
	}
	vector<int> &mus_tq_genMotherID()
	{
		if (not mus_tq_genMotherID_isLoaded) {
			if (mus_tq_genMotherID_branch != 0) 
				mus_tq_genMotherID_branch->GetEntry(index);
			mus_tq_genMotherID_isLoaded = true;
		}
		return mus_tq_genMotherID_;
	}
	vector<int> &trks_charge()
	{
		if (not trks_charge_isLoaded) {
			if (trks_charge_branch != 0) 
				trks_charge_branch->GetEntry(index);
			trks_charge_isLoaded = true;
		}
		return trks_charge_;
	}
	vector<int> &trks_lostHits()
	{
		if (not trks_lostHits_isLoaded) {
			if (trks_lostHits_branch != 0) 
				trks_lostHits_branch->GetEntry(index);
			trks_lostHits_isLoaded = true;
		}
		return trks_lostHits_;
	}
	vector<int> &trks_validHits()
	{
		if (not trks_validHits_isLoaded) {
			if (trks_validHits_branch != 0) 
				trks_validHits_branch->GetEntry(index);
			trks_validHits_isLoaded = true;
		}
		return trks_validHits_;
	}
	vector<int> &trk_elsidx()
	{
		if (not trk_elsidx_isLoaded) {
			if (trk_elsidx_branch != 0) 
				trk_elsidx_branch->GetEntry(index);
			trk_elsidx_isLoaded = true;
		}
		return trk_elsidx_;
	}
	vector<int> &trk_musidx()
	{
		if (not trk_musidx_isLoaded) {
			if (trk_musidx_branch != 0) 
				trk_musidx_branch->GetEntry(index);
			trk_musidx_isLoaded = true;
		}
		return trk_musidx_;
	}
	vector<vector<int> > &hyp_jets_mc_id()
	{
		if (not hyp_jets_mc_id_isLoaded) {
			if (hyp_jets_mc_id_branch != 0) 
				hyp_jets_mc_id_branch->GetEntry(index);
			hyp_jets_mc_id_isLoaded = true;
		}
		return hyp_jets_mc_id_;
	}
	vector<vector<int> > &hyp_jets_tq_genPartonMother_id()
	{
		if (not hyp_jets_tq_genPartonMother_id_isLoaded) {
			if (hyp_jets_tq_genPartonMother_id_branch != 0) 
				hyp_jets_tq_genPartonMother_id_branch->GetEntry(index);
			hyp_jets_tq_genPartonMother_id_isLoaded = true;
		}
		return hyp_jets_tq_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_jets_tq_genParton_id()
	{
		if (not hyp_jets_tq_genParton_id_isLoaded) {
			if (hyp_jets_tq_genParton_id_branch != 0) 
				hyp_jets_tq_genParton_id_branch->GetEntry(index);
			hyp_jets_tq_genParton_id_isLoaded = true;
		}
		return hyp_jets_tq_genParton_id_;
	}
	vector<vector<int> > &hyp_jets_tq_partonFlavour()
	{
		if (not hyp_jets_tq_partonFlavour_isLoaded) {
			if (hyp_jets_tq_partonFlavour_branch != 0) 
				hyp_jets_tq_partonFlavour_branch->GetEntry(index);
			hyp_jets_tq_partonFlavour_isLoaded = true;
		}
		return hyp_jets_tq_partonFlavour_;
	}
	vector<vector<int> > &hyp_other_jets_mc_id()
	{
		if (not hyp_other_jets_mc_id_isLoaded) {
			if (hyp_other_jets_mc_id_branch != 0) 
				hyp_other_jets_mc_id_branch->GetEntry(index);
			hyp_other_jets_mc_id_isLoaded = true;
		}
		return hyp_other_jets_mc_id_;
	}
	vector<vector<int> > &hyp_other_jets_tq_genPartonMother_id()
	{
		if (not hyp_other_jets_tq_genPartonMother_id_isLoaded) {
			if (hyp_other_jets_tq_genPartonMother_id_branch != 0) 
				hyp_other_jets_tq_genPartonMother_id_branch->GetEntry(index);
			hyp_other_jets_tq_genPartonMother_id_isLoaded = true;
		}
		return hyp_other_jets_tq_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_other_jets_tq_genParton_id()
	{
		if (not hyp_other_jets_tq_genParton_id_isLoaded) {
			if (hyp_other_jets_tq_genParton_id_branch != 0) 
				hyp_other_jets_tq_genParton_id_branch->GetEntry(index);
			hyp_other_jets_tq_genParton_id_isLoaded = true;
		}
		return hyp_other_jets_tq_genParton_id_;
	}
	vector<vector<int> > &hyp_other_jets_tq_partonFlavour()
	{
		if (not hyp_other_jets_tq_partonFlavour_isLoaded) {
			if (hyp_other_jets_tq_partonFlavour_branch != 0) 
				hyp_other_jets_tq_partonFlavour_branch->GetEntry(index);
			hyp_other_jets_tq_partonFlavour_isLoaded = true;
		}
		return hyp_other_jets_tq_partonFlavour_;
	}
	vector<vector<int> > &hyp_quadlep_jets_index()
	{
		if (not hyp_quadlep_jets_index_isLoaded) {
			if (hyp_quadlep_jets_index_branch != 0) 
				hyp_quadlep_jets_index_branch->GetEntry(index);
			hyp_quadlep_jets_index_isLoaded = true;
		}
		return hyp_quadlep_jets_index_;
	}
	vector<vector<int> > &hyp_trilep_jets_index()
	{
		if (not hyp_trilep_jets_index_isLoaded) {
			if (hyp_trilep_jets_index_branch != 0) 
				hyp_trilep_jets_index_branch->GetEntry(index);
			hyp_trilep_jets_index_isLoaded = true;
		}
		return hyp_trilep_jets_index_;
	}
	unsigned int &evt_nels()
	{
		if (not evt_nels_isLoaded) {
			if (evt_nels_branch != 0) 
				evt_nels_branch->GetEntry(index);
			evt_nels_isLoaded = true;
		}
		return evt_nels_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) 
				evt_event_branch->GetEntry(index);
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) 
				evt_run_branch->GetEntry(index);
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_njets()
	{
		if (not evt_njets_isLoaded) {
			if (evt_njets_branch != 0) 
				evt_njets_branch->GetEntry(index);
			evt_njets_isLoaded = true;
		}
		return evt_njets_;
	}
	vector<unsigned int> &hyp_quadlep_bucket()
	{
		if (not hyp_quadlep_bucket_isLoaded) {
			if (hyp_quadlep_bucket_branch != 0) 
				hyp_quadlep_bucket_branch->GetEntry(index);
			hyp_quadlep_bucket_isLoaded = true;
		}
		return hyp_quadlep_bucket_;
	}
	vector<unsigned int> &hyp_quadlep_first_index()
	{
		if (not hyp_quadlep_first_index_isLoaded) {
			if (hyp_quadlep_first_index_branch != 0) 
				hyp_quadlep_first_index_branch->GetEntry(index);
			hyp_quadlep_first_index_isLoaded = true;
		}
		return hyp_quadlep_first_index_;
	}
	vector<unsigned int> &hyp_quadlep_fourth_index()
	{
		if (not hyp_quadlep_fourth_index_isLoaded) {
			if (hyp_quadlep_fourth_index_branch != 0) 
				hyp_quadlep_fourth_index_branch->GetEntry(index);
			hyp_quadlep_fourth_index_isLoaded = true;
		}
		return hyp_quadlep_fourth_index_;
	}
	vector<unsigned int> &hyp_quadlep_second_index()
	{
		if (not hyp_quadlep_second_index_isLoaded) {
			if (hyp_quadlep_second_index_branch != 0) 
				hyp_quadlep_second_index_branch->GetEntry(index);
			hyp_quadlep_second_index_isLoaded = true;
		}
		return hyp_quadlep_second_index_;
	}
	vector<unsigned int> &hyp_quadlep_third_index()
	{
		if (not hyp_quadlep_third_index_isLoaded) {
			if (hyp_quadlep_third_index_branch != 0) 
				hyp_quadlep_third_index_branch->GetEntry(index);
			hyp_quadlep_third_index_isLoaded = true;
		}
		return hyp_quadlep_third_index_;
	}
	vector<unsigned int> &hyp_trilep_bucket()
	{
		if (not hyp_trilep_bucket_isLoaded) {
			if (hyp_trilep_bucket_branch != 0) 
				hyp_trilep_bucket_branch->GetEntry(index);
			hyp_trilep_bucket_isLoaded = true;
		}
		return hyp_trilep_bucket_;
	}
	vector<unsigned int> &hyp_trilep_first_index()
	{
		if (not hyp_trilep_first_index_isLoaded) {
			if (hyp_trilep_first_index_branch != 0) 
				hyp_trilep_first_index_branch->GetEntry(index);
			hyp_trilep_first_index_isLoaded = true;
		}
		return hyp_trilep_first_index_;
	}
	vector<unsigned int> &hyp_trilep_second_index()
	{
		if (not hyp_trilep_second_index_isLoaded) {
			if (hyp_trilep_second_index_branch != 0) 
				hyp_trilep_second_index_branch->GetEntry(index);
			hyp_trilep_second_index_isLoaded = true;
		}
		return hyp_trilep_second_index_;
	}
	vector<unsigned int> &hyp_trilep_third_index()
	{
		if (not hyp_trilep_third_index_isLoaded) {
			if (hyp_trilep_third_index_branch != 0) 
				hyp_trilep_third_index_branch->GetEntry(index);
			hyp_trilep_third_index_isLoaded = true;
		}
		return hyp_trilep_third_index_;
	}
};

#endif
