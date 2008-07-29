// -*- C++ -*-
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > 	l1met_p4;
	TBranch *l1met_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4;
	TBranch *els_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4In;
	TBranch *els_p4In_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_p4Out;
	TBranch *els_p4Out_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_trk_p4;
	TBranch *els_trk_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	genps_p4;
	TBranch *genps_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	genps_prod_vtx;
	TBranch *genps_prod_vtx_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_mc_p4;
	TBranch *hyp_ll_mc_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_p4;
	TBranch *hyp_ll_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_ll_trk_p4;
	TBranch *hyp_ll_trk_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_mc_p4;
	TBranch *hyp_lt_mc_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_p4;
	TBranch *hyp_lt_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_lt_trk_p4;
	TBranch *hyp_lt_trk_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	hyp_p4;
	TBranch *hyp_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_p4;
	TBranch *jets_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1emiso_p4;
	TBranch *l1emiso_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1emnoiso_p4;
	TBranch *l1emnoiso_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetsc_p4;
	TBranch *l1jetsc_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetsf_p4;
	TBranch *l1jetsf_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1jetst_p4;
	TBranch *l1jetst_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	l1mus_p4;
	TBranch *l1mus_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_p4;
	TBranch *mus_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_trk_p4;
	TBranch *mus_trk_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_tq_genMotherP4;
	TBranch *els_tq_genMotherP4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_tq_genP4;
	TBranch *els_tq_genP4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genJet_p4;
	TBranch *jets_tq_genJet_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genPartonMother_p4;
	TBranch *jets_tq_genPartonMother_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_genParton_p4;
	TBranch *jets_tq_genParton_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_tq_jet_p4;
	TBranch *jets_tq_jet_p4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_tq_genMotherP4;
	TBranch *mus_tq_genMotherP4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_tq_genP4;
	TBranch *mus_tq_genP4_branch;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	trks_trk_p4;
	TBranch *trks_trk_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_gp_p4;
	TBranch *hyp_jets_mc_gp_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_mc_p4;
	TBranch *hyp_jets_mc_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_p4;
	TBranch *hyp_jets_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_p4;
	TBranch *hyp_jets_tq_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_genPartonMother_p4;
	TBranch *hyp_jets_tq_genPartonMother_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_genParton_p4;
	TBranch *hyp_jets_tq_genParton_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_tq_jet_p4;
	TBranch *hyp_jets_tq_jet_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_gp_p4;
	TBranch *hyp_other_jets_mc_gp_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_mc_p4;
	TBranch *hyp_other_jets_mc_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_p4;
	TBranch *hyp_other_jets_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genJet_p4;
	TBranch *hyp_other_jets_tq_genJet_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genPartonMother_p4;
	TBranch *hyp_other_jets_tq_genPartonMother_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_genParton_p4;
	TBranch *hyp_other_jets_tq_genParton_p4_branch;
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_other_jets_tq_jet_p4;
	TBranch *hyp_other_jets_tq_jet_p4_branch;
	float	evt_kfactor;
	TBranch *evt_kfactor_branch;
	float	evt_weight;
	TBranch *evt_weight_branch;
	float	evt_xsec_excl;
	TBranch *evt_xsec_excl_branch;
	float	evt_xsec_incl;
	TBranch *evt_xsec_incl_branch;
	float	l1met_etHad;
	TBranch *l1met_etHad_branch;
	float	l1met_etTot;
	TBranch *l1met_etTot_branch;
	float	l1met_met;
	TBranch *l1met_met_branch;
	float	evt_met;
	TBranch *evt_met_branch;
	float	evt_metPhi;
	TBranch *evt_metPhi_branch;
	float	evt_metSig;
	TBranch *evt_metSig_branch;
	float	evt_met_jetcorr;
	TBranch *evt_met_jetcorr_branch;
	float	metphi_jetcorr;
	TBranch *metphi_jetcorr_branch;
	vector<float>	els_musdr;
	TBranch *els_musdr_branch;
	vector<float>	els_trkdr;
	TBranch *els_trkdr_branch;
	vector<float>	els_ESc;
	TBranch *els_ESc_branch;
	vector<float>	els_ESc_raw;
	TBranch *els_ESc_raw_branch;
	vector<float>	els_ESeed;
	TBranch *els_ESeed_branch;
	vector<float>	els_chi2;
	TBranch *els_chi2_branch;
	vector<float>	els_d0;
	TBranch *els_d0_branch;
	vector<float>	els_d0Err;
	TBranch *els_d0Err_branch;
	vector<float>	els_dEtaIn;
	TBranch *els_dEtaIn_branch;
	vector<float>	els_dEtaOut;
	TBranch *els_dEtaOut_branch;
	vector<float>	els_dPhiIn;
	TBranch *els_dPhiIn_branch;
	vector<float>	els_dPhiInPhiOut;
	TBranch *els_dPhiInPhiOut_branch;
	vector<float>	els_dPhiOut;
	TBranch *els_dPhiOut_branch;
	vector<float>	els_e3x3;
	TBranch *els_e3x3_branch;
	vector<float>	els_e5x5;
	TBranch *els_e5x5_branch;
	vector<float>	els_eOverPIn;
	TBranch *els_eOverPIn_branch;
	vector<float>	els_eSeedOverPOut;
	TBranch *els_eSeedOverPOut_branch;
	vector<float>	els_etaErr;
	TBranch *els_etaErr_branch;
	vector<float>	els_fBrem;
	TBranch *els_fBrem_branch;
	vector<float>	els_hOverE;
	TBranch *els_hOverE_branch;
	vector<float>	els_ndof;
	TBranch *els_ndof_branch;
	vector<float>	els_outerEta;
	TBranch *els_outerEta_branch;
	vector<float>	els_outerPhi;
	TBranch *els_outerPhi_branch;
	vector<float>	els_phiErr;
	TBranch *els_phiErr_branch;
	vector<float>	els_ptErr;
	TBranch *els_ptErr_branch;
	vector<float>	els_sigmaEtaEta;
	TBranch *els_sigmaEtaEta_branch;
	vector<float>	els_sigmaPhiPhi;
	TBranch *els_sigmaPhiPhi_branch;
	vector<float>	els_tkIso;
	TBranch *els_tkIso_branch;
	vector<float>	els_vertexphi;
	TBranch *els_vertexphi_branch;
	vector<float>	els_z0;
	TBranch *els_z0_branch;
	vector<float>	els_z0Err;
	TBranch *els_z0Err_branch;
	vector<float>	hyp_ll_chi2;
	TBranch *hyp_ll_chi2_branch;
	vector<float>	hyp_ll_d0;
	TBranch *hyp_ll_d0_branch;
	vector<float>	hyp_ll_d0Err;
	TBranch *hyp_ll_d0Err_branch;
	vector<float>	hyp_ll_etaErr;
	TBranch *hyp_ll_etaErr_branch;
	vector<float>	hyp_ll_iso;
	TBranch *hyp_ll_iso_branch;
	vector<float>	hyp_ll_ndof;
	TBranch *hyp_ll_ndof_branch;
	vector<float>	hyp_ll_outerEta;
	TBranch *hyp_ll_outerEta_branch;
	vector<float>	hyp_ll_outerPhi;
	TBranch *hyp_ll_outerPhi_branch;
	vector<float>	hyp_ll_phiErr;
	TBranch *hyp_ll_phiErr_branch;
	vector<float>	hyp_ll_ptErr;
	TBranch *hyp_ll_ptErr_branch;
	vector<float>	hyp_ll_tkIso;
	TBranch *hyp_ll_tkIso_branch;
	vector<float>	hyp_ll_vertexphi;
	TBranch *hyp_ll_vertexphi_branch;
	vector<float>	hyp_ll_z0;
	TBranch *hyp_ll_z0_branch;
	vector<float>	hyp_ll_z0Err;
	TBranch *hyp_ll_z0Err_branch;
	vector<float>	hyp_lt_chi2;
	TBranch *hyp_lt_chi2_branch;
	vector<float>	hyp_lt_d0;
	TBranch *hyp_lt_d0_branch;
	vector<float>	hyp_lt_d0Err;
	TBranch *hyp_lt_d0Err_branch;
	vector<float>	hyp_lt_etaErr;
	TBranch *hyp_lt_etaErr_branch;
	vector<float>	hyp_lt_iso;
	TBranch *hyp_lt_iso_branch;
	vector<float>	hyp_lt_ndof;
	TBranch *hyp_lt_ndof_branch;
	vector<float>	hyp_lt_outerEta;
	TBranch *hyp_lt_outerEta_branch;
	vector<float>	hyp_lt_outerPhi;
	TBranch *hyp_lt_outerPhi_branch;
	vector<float>	hyp_lt_phiErr;
	TBranch *hyp_lt_phiErr_branch;
	vector<float>	hyp_lt_ptErr;
	TBranch *hyp_lt_ptErr_branch;
	vector<float>	hyp_lt_tkIso;
	TBranch *hyp_lt_tkIso_branch;
	vector<float>	hyp_lt_vertexphi;
	TBranch *hyp_lt_vertexphi_branch;
	vector<float>	hyp_lt_z0;
	TBranch *hyp_lt_z0_branch;
	vector<float>	hyp_lt_z0Err;
	TBranch *hyp_lt_z0Err_branch;
	vector<float>	hyp_met;
	TBranch *hyp_met_branch;
	vector<float>	hyp_metAll;
	TBranch *hyp_metAll_branch;
	vector<float>	hyp_metAllCaloExp;
	TBranch *hyp_metAllCaloExp_branch;
	vector<float>	hyp_metCaloExp;
	TBranch *hyp_metCaloExp_branch;
	vector<float>	hyp_metCone;
	TBranch *hyp_metCone_branch;
	vector<float>	hyp_metDPhiJet10;
	TBranch *hyp_metDPhiJet10_branch;
	vector<float>	hyp_metDPhiJet15;
	TBranch *hyp_metDPhiJet15_branch;
	vector<float>	hyp_metDPhiJet20;
	TBranch *hyp_metDPhiJet20_branch;
	vector<float>	hyp_metDPhiTrk10;
	TBranch *hyp_metDPhiTrk10_branch;
	vector<float>	hyp_metDPhiTrk25;
	TBranch *hyp_metDPhiTrk25_branch;
	vector<float>	hyp_metDPhiTrk50;
	TBranch *hyp_metDPhiTrk50_branch;
	vector<float>	hyp_metJes10;
	TBranch *hyp_metJes10_branch;
	vector<float>	hyp_metJes15;
	TBranch *hyp_metJes15_branch;
	vector<float>	hyp_metJes30;
	TBranch *hyp_metJes30_branch;
	vector<float>	hyp_metJes5;
	TBranch *hyp_metJes5_branch;
	vector<float>	hyp_metJes50;
	TBranch *hyp_metJes50_branch;
	vector<float>	hyp_metNoCalo;
	TBranch *hyp_metNoCalo_branch;
	vector<float>	hyp_metPhi;
	TBranch *hyp_metPhi_branch;
	vector<float>	hyp_metPhiAll;
	TBranch *hyp_metPhiAll_branch;
	vector<float>	hyp_metPhiAllCaloExp;
	TBranch *hyp_metPhiAllCaloExp_branch;
	vector<float>	hyp_metPhiCaloExp;
	TBranch *hyp_metPhiCaloExp_branch;
	vector<float>	hyp_metPhiCone;
	TBranch *hyp_metPhiCone_branch;
	vector<float>	hyp_metPhiJes10;
	TBranch *hyp_metPhiJes10_branch;
	vector<float>	hyp_metPhiJes15;
	TBranch *hyp_metPhiJes15_branch;
	vector<float>	hyp_metPhiJes30;
	TBranch *hyp_metPhiJes30_branch;
	vector<float>	hyp_metPhiJes5;
	TBranch *hyp_metPhiJes5_branch;
	vector<float>	hyp_metPhiJes50;
	TBranch *hyp_metPhiJes50_branch;
	vector<float>	hyp_metPhiNoCalo;
	TBranch *hyp_metPhiNoCalo_branch;
	vector<float>	hyp_quadlep_met;
	TBranch *hyp_quadlep_met_branch;
	vector<float>	hyp_quadlep_metAll;
	TBranch *hyp_quadlep_metAll_branch;
	vector<float>	hyp_trilep_met;
	TBranch *hyp_trilep_met_branch;
	vector<float>	hyp_trilep_metAll;
	TBranch *hyp_trilep_metAll_branch;
	vector<float>	jets_EMFcor;
	TBranch *jets_EMFcor_branch;
	vector<float>	jets_chFrac;
	TBranch *jets_chFrac_branch;
	vector<float>	jets_cor;
	TBranch *jets_cor_branch;
	vector<float>	jets_emFrac;
	TBranch *jets_emFrac_branch;
	vector<float>	mus_eledr;
	TBranch *mus_eledr_branch;
	vector<float>	mus_jetdr;
	TBranch *mus_jetdr_branch;
	vector<float>	mus_trkdr;
	TBranch *mus_trkdr_branch;
	vector<float>	mus_chi2;
	TBranch *mus_chi2_branch;
	vector<float>	mus_d0;
	TBranch *mus_d0_branch;
	vector<float>	mus_d0Err;
	TBranch *mus_d0Err_branch;
	vector<float>	mus_e_em;
	TBranch *mus_e_em_branch;
	vector<float>	mus_e_emS9;
	TBranch *mus_e_emS9_branch;
	vector<float>	mus_e_had;
	TBranch *mus_e_had_branch;
	vector<float>	mus_e_hadS9;
	TBranch *mus_e_hadS9_branch;
	vector<float>	mus_e_ho;
	TBranch *mus_e_ho_branch;
	vector<float>	mus_e_hoS9;
	TBranch *mus_e_hoS9_branch;
	vector<float>	mus_etaErr;
	TBranch *mus_etaErr_branch;
	vector<float>	mus_gfit_chi2;
	TBranch *mus_gfit_chi2_branch;
	vector<float>	mus_gfit_ndof;
	TBranch *mus_gfit_ndof_branch;
	vector<float>	mus_iso;
	TBranch *mus_iso_branch;
	vector<float>	mus_iso03_emEt;
	TBranch *mus_iso03_emEt_branch;
	vector<float>	mus_iso03_hadEt;
	TBranch *mus_iso03_hadEt_branch;
	vector<float>	mus_iso03_hoEt;
	TBranch *mus_iso03_hoEt_branch;
	vector<float>	mus_iso03_sumPt;
	TBranch *mus_iso03_sumPt_branch;
	vector<float>	mus_iso05_emEt;
	TBranch *mus_iso05_emEt_branch;
	vector<float>	mus_iso05_hadEt;
	TBranch *mus_iso05_hadEt_branch;
	vector<float>	mus_iso05_hoEt;
	TBranch *mus_iso05_hoEt_branch;
	vector<float>	mus_iso05_sumPt;
	TBranch *mus_iso05_sumPt_branch;
	vector<float>	mus_ndof;
	TBranch *mus_ndof_branch;
	vector<float>	mus_outerEta;
	TBranch *mus_outerEta_branch;
	vector<float>	mus_outerPhi;
	TBranch *mus_outerPhi_branch;
	vector<float>	mus_phiErr;
	TBranch *mus_phiErr_branch;
	vector<float>	mus_ptErr;
	TBranch *mus_ptErr_branch;
	vector<float>	mus_vertexphi;
	TBranch *mus_vertexphi_branch;
	vector<float>	mus_z0;
	TBranch *mus_z0_branch;
	vector<float>	mus_z0Err;
	TBranch *mus_z0Err_branch;
	vector<float>	els_tq_LRComb;
	TBranch *els_tq_LRComb_branch;
	vector<float>	els_tq_caloIso;
	TBranch *els_tq_caloIso_branch;
	vector<float>	els_tq_egammaEcalIso;
	TBranch *els_tq_egammaEcalIso_branch;
	vector<float>	els_tq_egammaHcalIso;
	TBranch *els_tq_egammaHcalIso_branch;
	vector<float>	els_tq_egammaTkIso;
	TBranch *els_tq_egammaTkIso_branch;
	vector<float>	els_tq_electronIDRobust;
	TBranch *els_tq_electronIDRobust_branch;
	vector<float>	els_tq_leptonID;
	TBranch *els_tq_leptonID_branch;
	vector<float>	els_tq_trackIso;
	TBranch *els_tq_trackIso_branch;
	vector<float>	jets_tq_bCorrF;
	TBranch *jets_tq_bCorrF_branch;
	vector<float>	jets_tq_cCorrF;
	TBranch *jets_tq_cCorrF_branch;
	vector<float>	jets_tq_gluCorrF;
	TBranch *jets_tq_gluCorrF_branch;
	vector<float>	jets_tq_jetCharge;
	TBranch *jets_tq_jetCharge_branch;
	vector<float>	jets_tq_noCorrF;
	TBranch *jets_tq_noCorrF_branch;
	vector<float>	jets_tq_udsCorrF;
	TBranch *jets_tq_udsCorrF_branch;
	vector<float>	mus_tq_caloIso;
	TBranch *mus_tq_caloIso_branch;
	vector<float>	mus_tq_leptonID;
	TBranch *mus_tq_leptonID_branch;
	vector<float>	mus_tq_lrComb;
	TBranch *mus_tq_lrComb_branch;
	vector<float>	mus_tq_trackIso;
	TBranch *mus_tq_trackIso_branch;
	vector<float>	trks_chi2;
	TBranch *trks_chi2_branch;
	vector<float>	trks_d0;
	TBranch *trks_d0_branch;
	vector<float>	trks_d0Err;
	TBranch *trks_d0Err_branch;
	vector<float>	trks_etaErr;
	TBranch *trks_etaErr_branch;
	vector<float>	trks_ndof;
	TBranch *trks_ndof_branch;
	vector<float>	trks_outerEta;
	TBranch *trks_outerEta_branch;
	vector<float>	trks_outerPhi;
	TBranch *trks_outerPhi_branch;
	vector<float>	trks_phiErr;
	TBranch *trks_phiErr_branch;
	vector<float>	trks_ptErr;
	TBranch *trks_ptErr_branch;
	vector<float>	trks_vertexphi;
	TBranch *trks_vertexphi_branch;
	vector<float>	trks_z0;
	TBranch *trks_z0_branch;
	vector<float>	trks_z0Err;
	TBranch *trks_z0Err_branch;
	vector<vector<float> >	hyp_jets_EMFcor;
	TBranch *hyp_jets_EMFcor_branch;
	vector<vector<float> >	hyp_jets_chFrac;
	TBranch *hyp_jets_chFrac_branch;
	vector<vector<float> >	hyp_jets_cor;
	TBranch *hyp_jets_cor_branch;
	vector<vector<float> >	hyp_jets_emFrac;
	TBranch *hyp_jets_emFrac_branch;
	vector<vector<float> >	hyp_jets_mc_emEnergy;
	TBranch *hyp_jets_mc_emEnergy_branch;
	vector<vector<float> >	hyp_jets_mc_hadEnergy;
	TBranch *hyp_jets_mc_hadEnergy_branch;
	vector<vector<float> >	hyp_jets_mc_invEnergy;
	TBranch *hyp_jets_mc_invEnergy_branch;
	vector<vector<float> >	hyp_jets_mc_otherEnergy;
	TBranch *hyp_jets_mc_otherEnergy_branch;
	vector<vector<float> >	hyp_jets_tq_bCorrF;
	TBranch *hyp_jets_tq_bCorrF_branch;
	vector<vector<float> >	hyp_jets_tq_cCorrF;
	TBranch *hyp_jets_tq_cCorrF_branch;
	vector<vector<float> >	hyp_jets_tq_gluCorrF;
	TBranch *hyp_jets_tq_gluCorrF_branch;
	vector<vector<float> >	hyp_jets_tq_jetCharge;
	TBranch *hyp_jets_tq_jetCharge_branch;
	vector<vector<float> >	hyp_jets_tq_noCorrF;
	TBranch *hyp_jets_tq_noCorrF_branch;
	vector<vector<float> >	hyp_jets_tq_udsCorrF;
	TBranch *hyp_jets_tq_udsCorrF_branch;
	vector<vector<float> >	hyp_other_jets_EMFcor;
	TBranch *hyp_other_jets_EMFcor_branch;
	vector<vector<float> >	hyp_other_jets_chFrac;
	TBranch *hyp_other_jets_chFrac_branch;
	vector<vector<float> >	hyp_other_jets_cor;
	TBranch *hyp_other_jets_cor_branch;
	vector<vector<float> >	hyp_other_jets_emFrac;
	TBranch *hyp_other_jets_emFrac_branch;
	vector<vector<float> >	hyp_other_jets_mc_emEnergy;
	TBranch *hyp_other_jets_mc_emEnergy_branch;
	vector<vector<float> >	hyp_other_jets_mc_hadEnergy;
	TBranch *hyp_other_jets_mc_hadEnergy_branch;
	vector<vector<float> >	hyp_other_jets_mc_invEnergy;
	TBranch *hyp_other_jets_mc_invEnergy_branch;
	vector<vector<float> >	hyp_other_jets_mc_otherEnergy;
	TBranch *hyp_other_jets_mc_otherEnergy_branch;
	vector<vector<float> >	hyp_other_jets_tq_bCorrF;
	TBranch *hyp_other_jets_tq_bCorrF_branch;
	vector<vector<float> >	hyp_other_jets_tq_cCorrF;
	TBranch *hyp_other_jets_tq_cCorrF_branch;
	vector<vector<float> >	hyp_other_jets_tq_gluCorrF;
	TBranch *hyp_other_jets_tq_gluCorrF_branch;
	vector<vector<float> >	hyp_other_jets_tq_jetCharge;
	TBranch *hyp_other_jets_tq_jetCharge_branch;
	vector<vector<float> >	hyp_other_jets_tq_noCorrF;
	TBranch *hyp_other_jets_tq_noCorrF_branch;
	vector<vector<float> >	hyp_other_jets_tq_udsCorrF;
	TBranch *hyp_other_jets_tq_udsCorrF_branch;
	int	evt_HLT1;
	TBranch *evt_HLT1_branch;
	int	evt_HLT2;
	TBranch *evt_HLT2_branch;
	int	evt_HLT3;
	TBranch *evt_HLT3_branch;
	int	evt_HLT4;
	TBranch *evt_HLT4_branch;
	int	evt_L1_1;
	TBranch *evt_L1_1_branch;
	int	evt_L1_2;
	TBranch *evt_L1_2_branch;
	int	evt_L1_3;
	TBranch *evt_L1_3_branch;
	int	evt_L1_4;
	TBranch *evt_L1_4_branch;
	int	evt_event;
	TBranch *evt_event_branch;
	int	evt_run;
	TBranch *evt_run_branch;
	int	evt_nl1emiso;
	TBranch *evt_nl1emiso_branch;
	int	evt_nl1emnoiso;
	TBranch *evt_nl1emnoiso_branch;
	int	evt_nl1jetsc;
	TBranch *evt_nl1jetsc_branch;
	int	evt_nl1jetsf;
	TBranch *evt_nl1jetsf_branch;
	int	evt_nl1jetst;
	TBranch *evt_nl1jetst_branch;
	int	evt_nl1mus;
	TBranch *evt_nl1mus_branch;
	vector<int>	els_closestMuon;
	TBranch *els_closestMuon_branch;
	vector<int>	els_trkidx;
	TBranch *els_trkidx_branch;
	vector<int>	els_charge;
	TBranch *els_charge_branch;
	vector<int>	els_class;
	TBranch *els_class_branch;
	vector<int>	els_looseId;
	TBranch *els_looseId_branch;
	vector<int>	els_lostHits;
	TBranch *els_lostHits_branch;
	vector<int>	els_nSeed;
	TBranch *els_nSeed_branch;
	vector<int>	els_pass3looseId;
	TBranch *els_pass3looseId_branch;
	vector<int>	els_pass3simpleId;
	TBranch *els_pass3simpleId_branch;
	vector<int>	els_pass3tightId;
	TBranch *els_pass3tightId_branch;
	vector<int>	els_robustId;
	TBranch *els_robustId_branch;
	vector<int>	els_simpleIdPlus;
	TBranch *els_simpleIdPlus_branch;
	vector<int>	els_tightId;
	TBranch *els_tightId_branch;
	vector<int>	els_validHits;
	TBranch *els_validHits_branch;
	vector<int>	genps_id;
	TBranch *genps_id_branch;
	vector<int>	genps_id_mother;
	TBranch *genps_id_mother_branch;
	vector<int>	genps_status;
	TBranch *genps_status_branch;
	vector<int>	hyp_ll_charge;
	TBranch *hyp_ll_charge_branch;
	vector<int>	hyp_ll_id;
	TBranch *hyp_ll_id_branch;
	vector<int>	hyp_ll_index;
	TBranch *hyp_ll_index_branch;
	vector<int>	hyp_ll_lostHits;
	TBranch *hyp_ll_lostHits_branch;
	vector<int>	hyp_ll_mc_id;
	TBranch *hyp_ll_mc_id_branch;
	vector<int>	hyp_ll_mc_motherid;
	TBranch *hyp_ll_mc_motherid_branch;
	vector<int>	hyp_ll_validHits;
	TBranch *hyp_ll_validHits_branch;
	vector<int>	hyp_lt_charge;
	TBranch *hyp_lt_charge_branch;
	vector<int>	hyp_lt_id;
	TBranch *hyp_lt_id_branch;
	vector<int>	hyp_lt_index;
	TBranch *hyp_lt_index_branch;
	vector<int>	hyp_lt_lostHits;
	TBranch *hyp_lt_lostHits_branch;
	vector<int>	hyp_lt_mc_id;
	TBranch *hyp_lt_mc_id_branch;
	vector<int>	hyp_lt_mc_motherid;
	TBranch *hyp_lt_mc_motherid_branch;
	vector<int>	hyp_lt_validHits;
	TBranch *hyp_lt_validHits_branch;
	vector<int>	hyp_njets;
	TBranch *hyp_njets_branch;
	vector<int>	hyp_nojets;
	TBranch *hyp_nojets_branch;
	vector<int>	hyp_type;
	TBranch *hyp_type_branch;
	vector<int>	hyp_quadlep_first_type;
	TBranch *hyp_quadlep_first_type_branch;
	vector<int>	hyp_quadlep_fourth_type;
	TBranch *hyp_quadlep_fourth_type_branch;
	vector<int>	hyp_quadlep_second_type;
	TBranch *hyp_quadlep_second_type_branch;
	vector<int>	hyp_quadlep_third_type;
	TBranch *hyp_quadlep_third_type_branch;
	vector<int>	hyp_trilep_first_type;
	TBranch *hyp_trilep_first_type_branch;
	vector<int>	hyp_trilep_second_type;
	TBranch *hyp_trilep_second_type_branch;
	vector<int>	hyp_trilep_third_type;
	TBranch *hyp_trilep_third_type_branch;
	vector<int>	jets_closestElectron;
	TBranch *jets_closestElectron_branch;
	vector<int>	jets_closestMuon;
	TBranch *jets_closestMuon_branch;
	vector<int>	l1emiso_ieta;
	TBranch *l1emiso_ieta_branch;
	vector<int>	l1emiso_iphi;
	TBranch *l1emiso_iphi_branch;
	vector<int>	l1emiso_rawId;
	TBranch *l1emiso_rawId_branch;
	vector<int>	l1emiso_type;
	TBranch *l1emiso_type_branch;
	vector<int>	l1emnoiso_ieta;
	TBranch *l1emnoiso_ieta_branch;
	vector<int>	l1emnoiso_iphi;
	TBranch *l1emnoiso_iphi_branch;
	vector<int>	l1emnoiso_rawId;
	TBranch *l1emnoiso_rawId_branch;
	vector<int>	l1emnoiso_type;
	TBranch *l1emnoiso_type_branch;
	vector<int>	l1jetsc_ieta;
	TBranch *l1jetsc_ieta_branch;
	vector<int>	l1jetsc_iphi;
	TBranch *l1jetsc_iphi_branch;
	vector<int>	l1jetsc_rawId;
	TBranch *l1jetsc_rawId_branch;
	vector<int>	l1jetsc_type;
	TBranch *l1jetsc_type_branch;
	vector<int>	l1jetsf_ieta;
	TBranch *l1jetsf_ieta_branch;
	vector<int>	l1jetsf_iphi;
	TBranch *l1jetsf_iphi_branch;
	vector<int>	l1jetsf_rawId;
	TBranch *l1jetsf_rawId_branch;
	vector<int>	l1jetsf_type;
	TBranch *l1jetsf_type_branch;
	vector<int>	l1jetst_ieta;
	TBranch *l1jetst_ieta_branch;
	vector<int>	l1jetst_iphi;
	TBranch *l1jetst_iphi_branch;
	vector<int>	l1jetst_rawId;
	TBranch *l1jetst_rawId_branch;
	vector<int>	l1jetst_type;
	TBranch *l1jetst_type_branch;
	vector<int>	l1mus_flags;
	TBranch *l1mus_flags_branch;
	vector<int>	l1mus_q;
	TBranch *l1mus_q_branch;
	vector<int>	l1mus_qual;
	TBranch *l1mus_qual_branch;
	vector<int>	l1mus_qualFlags;
	TBranch *l1mus_qualFlags_branch;
	vector<int>	mus_closestEle;
	TBranch *mus_closestEle_branch;
	vector<int>	mus_closestJet;
	TBranch *mus_closestJet_branch;
	vector<int>	mus_trkidx;
	TBranch *mus_trkidx_branch;
	vector<int>	mus_charge;
	TBranch *mus_charge_branch;
	vector<int>	mus_gfit_validHits;
	TBranch *mus_gfit_validHits_branch;
	vector<int>	mus_iso03_ntrk;
	TBranch *mus_iso03_ntrk_branch;
	vector<int>	mus_iso05_ntrk;
	TBranch *mus_iso05_ntrk_branch;
	vector<int>	mus_lostHits;
	TBranch *mus_lostHits_branch;
	vector<int>	mus_nmatches;
	TBranch *mus_nmatches_branch;
	vector<int>	mus_pid_TM2DCompatibilityLoose;
	TBranch *mus_pid_TM2DCompatibilityLoose_branch;
	vector<int>	mus_pid_TM2DCompatibilityTight;
	TBranch *mus_pid_TM2DCompatibilityTight_branch;
	vector<int>	mus_pid_TMLastStationLoose;
	TBranch *mus_pid_TMLastStationLoose_branch;
	vector<int>	mus_pid_TMLastStationTight;
	TBranch *mus_pid_TMLastStationTight_branch;
	vector<int>	mus_trkrefkey;
	TBranch *mus_trkrefkey_branch;
	vector<int>	mus_validHits;
	TBranch *mus_validHits_branch;
	vector<int>	els_tq_egammaTkNumIso;
	TBranch *els_tq_egammaTkNumIso_branch;
	vector<int>	els_tq_genID;
	TBranch *els_tq_genID_branch;
	vector<int>	els_tq_genMotherID;
	TBranch *els_tq_genMotherID_branch;
	vector<int>	jets_tq_genPartonMother_id;
	TBranch *jets_tq_genPartonMother_id_branch;
	vector<int>	jets_tq_genParton_id;
	TBranch *jets_tq_genParton_id_branch;
	vector<int>	jets_tq_partonFlavour;
	TBranch *jets_tq_partonFlavour_branch;
	vector<int>	mus_tq_genID;
	TBranch *mus_tq_genID_branch;
	vector<int>	mus_tq_genMotherID;
	TBranch *mus_tq_genMotherID_branch;
	vector<int>	trks_charge;
	TBranch *trks_charge_branch;
	vector<int>	trks_lostHits;
	TBranch *trks_lostHits_branch;
	vector<int>	trks_validHits;
	TBranch *trks_validHits_branch;
	vector<vector<int> >	hyp_jets_mc_id;
	TBranch *hyp_jets_mc_id_branch;
	vector<vector<int> >	hyp_jets_tq_genPartonMother_id;
	TBranch *hyp_jets_tq_genPartonMother_id_branch;
	vector<vector<int> >	hyp_jets_tq_genParton_id;
	TBranch *hyp_jets_tq_genParton_id_branch;
	vector<vector<int> >	hyp_jets_tq_partonFlavour;
	TBranch *hyp_jets_tq_partonFlavour_branch;
	vector<vector<int> >	hyp_other_jets_mc_id;
	TBranch *hyp_other_jets_mc_id_branch;
	vector<vector<int> >	hyp_other_jets_tq_genPartonMother_id;
	TBranch *hyp_other_jets_tq_genPartonMother_id_branch;
	vector<vector<int> >	hyp_other_jets_tq_genParton_id;
	TBranch *hyp_other_jets_tq_genParton_id_branch;
	vector<vector<int> >	hyp_other_jets_tq_partonFlavour;
	TBranch *hyp_other_jets_tq_partonFlavour_branch;
	vector<vector<int> >	hyp_quadlep_jets_index;
	TBranch *hyp_quadlep_jets_index_branch;
	vector<vector<int> >	hyp_trilep_jets_index;
	TBranch *hyp_trilep_jets_index_branch;
	unsigned int	evt_nels;
	TBranch *evt_nels_branch;
	unsigned int	evt_njets;
	TBranch *evt_njets_branch;
	vector<unsigned int>	hyp_quadlep_bucket;
	TBranch *hyp_quadlep_bucket_branch;
	vector<unsigned int>	hyp_quadlep_first_index;
	TBranch *hyp_quadlep_first_index_branch;
	vector<unsigned int>	hyp_quadlep_fourth_index;
	TBranch *hyp_quadlep_fourth_index_branch;
	vector<unsigned int>	hyp_quadlep_second_index;
	TBranch *hyp_quadlep_second_index_branch;
	vector<unsigned int>	hyp_quadlep_third_index;
	TBranch *hyp_quadlep_third_index_branch;
	vector<unsigned int>	hyp_trilep_bucket;
	TBranch *hyp_trilep_bucket_branch;
	vector<unsigned int>	hyp_trilep_first_index;
	TBranch *hyp_trilep_first_index_branch;
	vector<unsigned int>	hyp_trilep_second_index;
	TBranch *hyp_trilep_second_index_branch;
	vector<unsigned int>	hyp_trilep_third_index;
	TBranch *hyp_trilep_third_index_branch;
void Init(TTree *tree) {
	l1met_p4_branch = tree->GetBranch(tree->GetAlias("l1met_p4"));
	l1met_p4_branch->SetAddress(&l1met_p4);
	if(l1met_p4_branch == 0 ) {
	cout << "Branch l1met_p4 does not exist." << endl;
	}
	els_p4_branch = tree->GetBranch(tree->GetAlias("els_p4"));
	els_p4_branch->SetAddress(&els_p4);
	if(els_p4_branch == 0 ) {
	cout << "Branch els_p4 does not exist." << endl;
	}
	els_p4In_branch = tree->GetBranch(tree->GetAlias("els_p4In"));
	els_p4In_branch->SetAddress(&els_p4In);
	if(els_p4In_branch == 0 ) {
	cout << "Branch els_p4In does not exist." << endl;
	}
	els_p4Out_branch = tree->GetBranch(tree->GetAlias("els_p4Out"));
	els_p4Out_branch->SetAddress(&els_p4Out);
	if(els_p4Out_branch == 0 ) {
	cout << "Branch els_p4Out does not exist." << endl;
	}
	els_trk_p4_branch = tree->GetBranch(tree->GetAlias("els_trk_p4"));
	els_trk_p4_branch->SetAddress(&els_trk_p4);
	if(els_trk_p4_branch == 0 ) {
	cout << "Branch els_trk_p4 does not exist." << endl;
	}
	genps_p4_branch = tree->GetBranch(tree->GetAlias("genps_p4"));
	genps_p4_branch->SetAddress(&genps_p4);
	if(genps_p4_branch == 0 ) {
	cout << "Branch genps_p4 does not exist." << endl;
	}
	genps_prod_vtx_branch = tree->GetBranch(tree->GetAlias("genps_prod_vtx"));
	genps_prod_vtx_branch->SetAddress(&genps_prod_vtx);
	if(genps_prod_vtx_branch == 0 ) {
	cout << "Branch genps_prod_vtx does not exist." << endl;
	}
	hyp_ll_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_p4"));
	hyp_ll_mc_p4_branch->SetAddress(&hyp_ll_mc_p4);
	if(hyp_ll_mc_p4_branch == 0 ) {
	cout << "Branch hyp_ll_mc_p4 does not exist." << endl;
	}
	hyp_ll_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_p4"));
	hyp_ll_p4_branch->SetAddress(&hyp_ll_p4);
	if(hyp_ll_p4_branch == 0 ) {
	cout << "Branch hyp_ll_p4 does not exist." << endl;
	}
	hyp_ll_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_ll_trk_p4"));
	hyp_ll_trk_p4_branch->SetAddress(&hyp_ll_trk_p4);
	if(hyp_ll_trk_p4_branch == 0 ) {
	cout << "Branch hyp_ll_trk_p4 does not exist." << endl;
	}
	hyp_lt_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_p4"));
	hyp_lt_mc_p4_branch->SetAddress(&hyp_lt_mc_p4);
	if(hyp_lt_mc_p4_branch == 0 ) {
	cout << "Branch hyp_lt_mc_p4 does not exist." << endl;
	}
	hyp_lt_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_p4"));
	hyp_lt_p4_branch->SetAddress(&hyp_lt_p4);
	if(hyp_lt_p4_branch == 0 ) {
	cout << "Branch hyp_lt_p4 does not exist." << endl;
	}
	hyp_lt_trk_p4_branch = tree->GetBranch(tree->GetAlias("hyp_lt_trk_p4"));
	hyp_lt_trk_p4_branch->SetAddress(&hyp_lt_trk_p4);
	if(hyp_lt_trk_p4_branch == 0 ) {
	cout << "Branch hyp_lt_trk_p4 does not exist." << endl;
	}
	hyp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_p4"));
	hyp_p4_branch->SetAddress(&hyp_p4);
	if(hyp_p4_branch == 0 ) {
	cout << "Branch hyp_p4 does not exist." << endl;
	}
	jets_p4_branch = tree->GetBranch(tree->GetAlias("jets_p4"));
	jets_p4_branch->SetAddress(&jets_p4);
	if(jets_p4_branch == 0 ) {
	cout << "Branch jets_p4 does not exist." << endl;
	}
	l1emiso_p4_branch = tree->GetBranch(tree->GetAlias("l1emiso_p4"));
	l1emiso_p4_branch->SetAddress(&l1emiso_p4);
	if(l1emiso_p4_branch == 0 ) {
	cout << "Branch l1emiso_p4 does not exist." << endl;
	}
	l1emnoiso_p4_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_p4"));
	l1emnoiso_p4_branch->SetAddress(&l1emnoiso_p4);
	if(l1emnoiso_p4_branch == 0 ) {
	cout << "Branch l1emnoiso_p4 does not exist." << endl;
	}
	l1jetsc_p4_branch = tree->GetBranch(tree->GetAlias("l1jetsc_p4"));
	l1jetsc_p4_branch->SetAddress(&l1jetsc_p4);
	if(l1jetsc_p4_branch == 0 ) {
	cout << "Branch l1jetsc_p4 does not exist." << endl;
	}
	l1jetsf_p4_branch = tree->GetBranch(tree->GetAlias("l1jetsf_p4"));
	l1jetsf_p4_branch->SetAddress(&l1jetsf_p4);
	if(l1jetsf_p4_branch == 0 ) {
	cout << "Branch l1jetsf_p4 does not exist." << endl;
	}
	l1jetst_p4_branch = tree->GetBranch(tree->GetAlias("l1jetst_p4"));
	l1jetst_p4_branch->SetAddress(&l1jetst_p4);
	if(l1jetst_p4_branch == 0 ) {
	cout << "Branch l1jetst_p4 does not exist." << endl;
	}
	l1mus_p4_branch = tree->GetBranch(tree->GetAlias("l1mus_p4"));
	l1mus_p4_branch->SetAddress(&l1mus_p4);
	if(l1mus_p4_branch == 0 ) {
	cout << "Branch l1mus_p4 does not exist." << endl;
	}
	mus_p4_branch = tree->GetBranch(tree->GetAlias("mus_p4"));
	mus_p4_branch->SetAddress(&mus_p4);
	if(mus_p4_branch == 0 ) {
	cout << "Branch mus_p4 does not exist." << endl;
	}
	mus_trk_p4_branch = tree->GetBranch(tree->GetAlias("mus_trk_p4"));
	mus_trk_p4_branch->SetAddress(&mus_trk_p4);
	if(mus_trk_p4_branch == 0 ) {
	cout << "Branch mus_trk_p4 does not exist." << endl;
	}
	els_tq_genMotherP4_branch = tree->GetBranch(tree->GetAlias("els_tq_genMotherP4"));
	els_tq_genMotherP4_branch->SetAddress(&els_tq_genMotherP4);
	if(els_tq_genMotherP4_branch == 0 ) {
	cout << "Branch els_tq_genMotherP4 does not exist." << endl;
	}
	els_tq_genP4_branch = tree->GetBranch(tree->GetAlias("els_tq_genP4"));
	els_tq_genP4_branch->SetAddress(&els_tq_genP4);
	if(els_tq_genP4_branch == 0 ) {
	cout << "Branch els_tq_genP4 does not exist." << endl;
	}
	jets_tq_genJet_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genJet_p4"));
	jets_tq_genJet_p4_branch->SetAddress(&jets_tq_genJet_p4);
	if(jets_tq_genJet_p4_branch == 0 ) {
	cout << "Branch jets_tq_genJet_p4 does not exist." << endl;
	}
	jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genPartonMother_p4"));
	jets_tq_genPartonMother_p4_branch->SetAddress(&jets_tq_genPartonMother_p4);
	if(jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_genParton_p4"));
	jets_tq_genParton_p4_branch->SetAddress(&jets_tq_genParton_p4);
	if(jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch jets_tq_genParton_p4 does not exist." << endl;
	}
	jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("jets_tq_jet_p4"));
	jets_tq_jet_p4_branch->SetAddress(&jets_tq_jet_p4);
	if(jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch jets_tq_jet_p4 does not exist." << endl;
	}
	mus_tq_genMotherP4_branch = tree->GetBranch(tree->GetAlias("mus_tq_genMotherP4"));
	mus_tq_genMotherP4_branch->SetAddress(&mus_tq_genMotherP4);
	if(mus_tq_genMotherP4_branch == 0 ) {
	cout << "Branch mus_tq_genMotherP4 does not exist." << endl;
	}
	mus_tq_genP4_branch = tree->GetBranch(tree->GetAlias("mus_tq_genP4"));
	mus_tq_genP4_branch->SetAddress(&mus_tq_genP4);
	if(mus_tq_genP4_branch == 0 ) {
	cout << "Branch mus_tq_genP4 does not exist." << endl;
	}
	trks_trk_p4_branch = tree->GetBranch(tree->GetAlias("trks_trk_p4"));
	trks_trk_p4_branch->SetAddress(&trks_trk_p4);
	if(trks_trk_p4_branch == 0 ) {
	cout << "Branch trks_trk_p4 does not exist." << endl;
	}
  tree->SetMakeClass(1);
	hyp_jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_gp_p4") != 0) {
		hyp_jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_gp_p4"));
		hyp_jets_mc_gp_p4_branch->SetAddress(&hyp_jets_mc_gp_p4);
	}
	if(hyp_jets_mc_gp_p4_branch == 0 ) {
	cout << "Branch hyp_jets_mc_gp_p4 does not exist." << endl;
	}
	hyp_jets_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_p4") != 0) {
		hyp_jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_p4"));
		hyp_jets_mc_p4_branch->SetAddress(&hyp_jets_mc_p4);
	}
	if(hyp_jets_mc_p4_branch == 0 ) {
	cout << "Branch hyp_jets_mc_p4 does not exist." << endl;
	}
	hyp_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_p4") != 0) {
		hyp_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_p4"));
		hyp_jets_p4_branch->SetAddress(&hyp_jets_p4);
	}
	if(hyp_jets_p4_branch == 0 ) {
	cout << "Branch hyp_jets_p4 does not exist." << endl;
	}
	hyp_jets_tq_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_p4") != 0) {
		hyp_jets_tq_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_p4"));
		hyp_jets_tq_p4_branch->SetAddress(&hyp_jets_tq_p4);
	}
	if(hyp_jets_tq_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_p4 does not exist." << endl;
	}
	hyp_jets_tq_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genPartonMother_p4") != 0) {
		hyp_jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genPartonMother_p4"));
		hyp_jets_tq_genPartonMother_p4_branch->SetAddress(&hyp_jets_tq_genPartonMother_p4);
	}
	if(hyp_jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	hyp_jets_tq_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genParton_p4") != 0) {
		hyp_jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genParton_p4"));
		hyp_jets_tq_genParton_p4_branch->SetAddress(&hyp_jets_tq_genParton_p4);
	}
	if(hyp_jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genParton_p4 does not exist." << endl;
	}
	hyp_jets_tq_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_jet_p4") != 0) {
		hyp_jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_jet_p4"));
		hyp_jets_tq_jet_p4_branch->SetAddress(&hyp_jets_tq_jet_p4);
	}
	if(hyp_jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch hyp_jets_tq_jet_p4 does not exist." << endl;
	}
	hyp_other_jets_mc_gp_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_gp_p4") != 0) {
		hyp_other_jets_mc_gp_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_gp_p4"));
		hyp_other_jets_mc_gp_p4_branch->SetAddress(&hyp_other_jets_mc_gp_p4);
	}
	if(hyp_other_jets_mc_gp_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_gp_p4 does not exist." << endl;
	}
	hyp_other_jets_mc_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_p4") != 0) {
		hyp_other_jets_mc_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_p4"));
		hyp_other_jets_mc_p4_branch->SetAddress(&hyp_other_jets_mc_p4);
	}
	if(hyp_other_jets_mc_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_p4 does not exist." << endl;
	}
	hyp_other_jets_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_p4") != 0) {
		hyp_other_jets_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_p4"));
		hyp_other_jets_p4_branch->SetAddress(&hyp_other_jets_p4);
	}
	if(hyp_other_jets_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_genJet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genJet_p4") != 0) {
		hyp_other_jets_tq_genJet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genJet_p4"));
		hyp_other_jets_tq_genJet_p4_branch->SetAddress(&hyp_other_jets_tq_genJet_p4);
	}
	if(hyp_other_jets_tq_genJet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genJet_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_genPartonMother_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genPartonMother_p4") != 0) {
		hyp_other_jets_tq_genPartonMother_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genPartonMother_p4"));
		hyp_other_jets_tq_genPartonMother_p4_branch->SetAddress(&hyp_other_jets_tq_genPartonMother_p4);
	}
	if(hyp_other_jets_tq_genPartonMother_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genPartonMother_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_genParton_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genParton_p4") != 0) {
		hyp_other_jets_tq_genParton_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genParton_p4"));
		hyp_other_jets_tq_genParton_p4_branch->SetAddress(&hyp_other_jets_tq_genParton_p4);
	}
	if(hyp_other_jets_tq_genParton_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genParton_p4 does not exist." << endl;
	}
	hyp_other_jets_tq_jet_p4_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_jet_p4") != 0) {
		hyp_other_jets_tq_jet_p4_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_jet_p4"));
		hyp_other_jets_tq_jet_p4_branch->SetAddress(&hyp_other_jets_tq_jet_p4);
	}
	if(hyp_other_jets_tq_jet_p4_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_jet_p4 does not exist." << endl;
	}
	evt_kfactor_branch = 0;
	if (tree->GetAlias("evt_kfactor") != 0) {
		evt_kfactor_branch = tree->GetBranch(tree->GetAlias("evt_kfactor"));
		evt_kfactor_branch->SetAddress(&evt_kfactor);
	}
	if(evt_kfactor_branch == 0 ) {
	cout << "Branch evt_kfactor does not exist." << endl;
	}
	evt_weight_branch = 0;
	if (tree->GetAlias("evt_weight") != 0) {
		evt_weight_branch = tree->GetBranch(tree->GetAlias("evt_weight"));
		evt_weight_branch->SetAddress(&evt_weight);
	}
	if(evt_weight_branch == 0 ) {
	cout << "Branch evt_weight does not exist." << endl;
	}
	evt_xsec_excl_branch = 0;
	if (tree->GetAlias("evt_xsec_excl") != 0) {
		evt_xsec_excl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_excl"));
		evt_xsec_excl_branch->SetAddress(&evt_xsec_excl);
	}
	if(evt_xsec_excl_branch == 0 ) {
	cout << "Branch evt_xsec_excl does not exist." << endl;
	}
	evt_xsec_incl_branch = 0;
	if (tree->GetAlias("evt_xsec_incl") != 0) {
		evt_xsec_incl_branch = tree->GetBranch(tree->GetAlias("evt_xsec_incl"));
		evt_xsec_incl_branch->SetAddress(&evt_xsec_incl);
	}
	if(evt_xsec_incl_branch == 0 ) {
	cout << "Branch evt_xsec_incl does not exist." << endl;
	}
	l1met_etHad_branch = 0;
	if (tree->GetAlias("l1met_etHad") != 0) {
		l1met_etHad_branch = tree->GetBranch(tree->GetAlias("l1met_etHad"));
		l1met_etHad_branch->SetAddress(&l1met_etHad);
	}
	if(l1met_etHad_branch == 0 ) {
	cout << "Branch l1met_etHad does not exist." << endl;
	}
	l1met_etTot_branch = 0;
	if (tree->GetAlias("l1met_etTot") != 0) {
		l1met_etTot_branch = tree->GetBranch(tree->GetAlias("l1met_etTot"));
		l1met_etTot_branch->SetAddress(&l1met_etTot);
	}
	if(l1met_etTot_branch == 0 ) {
	cout << "Branch l1met_etTot does not exist." << endl;
	}
	l1met_met_branch = 0;
	if (tree->GetAlias("l1met_met") != 0) {
		l1met_met_branch = tree->GetBranch(tree->GetAlias("l1met_met"));
		l1met_met_branch->SetAddress(&l1met_met);
	}
	if(l1met_met_branch == 0 ) {
	cout << "Branch l1met_met does not exist." << endl;
	}
	evt_met_branch = 0;
	if (tree->GetAlias("evt_met") != 0) {
		evt_met_branch = tree->GetBranch(tree->GetAlias("evt_met"));
		evt_met_branch->SetAddress(&evt_met);
	}
	if(evt_met_branch == 0 ) {
	cout << "Branch evt_met does not exist." << endl;
	}
	evt_metPhi_branch = 0;
	if (tree->GetAlias("evt_metPhi") != 0) {
		evt_metPhi_branch = tree->GetBranch(tree->GetAlias("evt_metPhi"));
		evt_metPhi_branch->SetAddress(&evt_metPhi);
	}
	if(evt_metPhi_branch == 0 ) {
	cout << "Branch evt_metPhi does not exist." << endl;
	}
	evt_metSig_branch = 0;
	if (tree->GetAlias("evt_metSig") != 0) {
		evt_metSig_branch = tree->GetBranch(tree->GetAlias("evt_metSig"));
		evt_metSig_branch->SetAddress(&evt_metSig);
	}
	if(evt_metSig_branch == 0 ) {
	cout << "Branch evt_metSig does not exist." << endl;
	}
	evt_met_jetcorr_branch = 0;
	if (tree->GetAlias("evt_met_jetcorr") != 0) {
		evt_met_jetcorr_branch = tree->GetBranch(tree->GetAlias("evt_met_jetcorr"));
		evt_met_jetcorr_branch->SetAddress(&evt_met_jetcorr);
	}
	if(evt_met_jetcorr_branch == 0 ) {
	cout << "Branch evt_met_jetcorr does not exist." << endl;
	}
	metphi_jetcorr_branch = 0;
	if (tree->GetAlias("metphi_jetcorr") != 0) {
		metphi_jetcorr_branch = tree->GetBranch(tree->GetAlias("metphi_jetcorr"));
		metphi_jetcorr_branch->SetAddress(&metphi_jetcorr);
	}
	if(metphi_jetcorr_branch == 0 ) {
	cout << "Branch metphi_jetcorr does not exist." << endl;
	}
	els_musdr_branch = 0;
	if (tree->GetAlias("els_musdr") != 0) {
		els_musdr_branch = tree->GetBranch(tree->GetAlias("els_musdr"));
		els_musdr_branch->SetAddress(&els_musdr);
	}
	if(els_musdr_branch == 0 ) {
	cout << "Branch els_musdr does not exist." << endl;
	}
	els_trkdr_branch = 0;
	if (tree->GetAlias("els_trkdr") != 0) {
		els_trkdr_branch = tree->GetBranch(tree->GetAlias("els_trkdr"));
		els_trkdr_branch->SetAddress(&els_trkdr);
	}
	if(els_trkdr_branch == 0 ) {
	cout << "Branch els_trkdr does not exist." << endl;
	}
	els_ESc_branch = 0;
	if (tree->GetAlias("els_ESc") != 0) {
		els_ESc_branch = tree->GetBranch(tree->GetAlias("els_ESc"));
		els_ESc_branch->SetAddress(&els_ESc);
	}
	if(els_ESc_branch == 0 ) {
	cout << "Branch els_ESc does not exist." << endl;
	}
	els_ESc_raw_branch = 0;
	if (tree->GetAlias("els_ESc_raw") != 0) {
		els_ESc_raw_branch = tree->GetBranch(tree->GetAlias("els_ESc_raw"));
		els_ESc_raw_branch->SetAddress(&els_ESc_raw);
	}
	if(els_ESc_raw_branch == 0 ) {
	cout << "Branch els_ESc_raw does not exist." << endl;
	}
	els_ESeed_branch = 0;
	if (tree->GetAlias("els_ESeed") != 0) {
		els_ESeed_branch = tree->GetBranch(tree->GetAlias("els_ESeed"));
		els_ESeed_branch->SetAddress(&els_ESeed);
	}
	if(els_ESeed_branch == 0 ) {
	cout << "Branch els_ESeed does not exist." << endl;
	}
	els_chi2_branch = 0;
	if (tree->GetAlias("els_chi2") != 0) {
		els_chi2_branch = tree->GetBranch(tree->GetAlias("els_chi2"));
		els_chi2_branch->SetAddress(&els_chi2);
	}
	if(els_chi2_branch == 0 ) {
	cout << "Branch els_chi2 does not exist." << endl;
	}
	els_d0_branch = 0;
	if (tree->GetAlias("els_d0") != 0) {
		els_d0_branch = tree->GetBranch(tree->GetAlias("els_d0"));
		els_d0_branch->SetAddress(&els_d0);
	}
	if(els_d0_branch == 0 ) {
	cout << "Branch els_d0 does not exist." << endl;
	}
	els_d0Err_branch = 0;
	if (tree->GetAlias("els_d0Err") != 0) {
		els_d0Err_branch = tree->GetBranch(tree->GetAlias("els_d0Err"));
		els_d0Err_branch->SetAddress(&els_d0Err);
	}
	if(els_d0Err_branch == 0 ) {
	cout << "Branch els_d0Err does not exist." << endl;
	}
	els_dEtaIn_branch = 0;
	if (tree->GetAlias("els_dEtaIn") != 0) {
		els_dEtaIn_branch = tree->GetBranch(tree->GetAlias("els_dEtaIn"));
		els_dEtaIn_branch->SetAddress(&els_dEtaIn);
	}
	if(els_dEtaIn_branch == 0 ) {
	cout << "Branch els_dEtaIn does not exist." << endl;
	}
	els_dEtaOut_branch = 0;
	if (tree->GetAlias("els_dEtaOut") != 0) {
		els_dEtaOut_branch = tree->GetBranch(tree->GetAlias("els_dEtaOut"));
		els_dEtaOut_branch->SetAddress(&els_dEtaOut);
	}
	if(els_dEtaOut_branch == 0 ) {
	cout << "Branch els_dEtaOut does not exist." << endl;
	}
	els_dPhiIn_branch = 0;
	if (tree->GetAlias("els_dPhiIn") != 0) {
		els_dPhiIn_branch = tree->GetBranch(tree->GetAlias("els_dPhiIn"));
		els_dPhiIn_branch->SetAddress(&els_dPhiIn);
	}
	if(els_dPhiIn_branch == 0 ) {
	cout << "Branch els_dPhiIn does not exist." << endl;
	}
	els_dPhiInPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiInPhiOut") != 0) {
		els_dPhiInPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiInPhiOut"));
		els_dPhiInPhiOut_branch->SetAddress(&els_dPhiInPhiOut);
	}
	if(els_dPhiInPhiOut_branch == 0 ) {
	cout << "Branch els_dPhiInPhiOut does not exist." << endl;
	}
	els_dPhiOut_branch = 0;
	if (tree->GetAlias("els_dPhiOut") != 0) {
		els_dPhiOut_branch = tree->GetBranch(tree->GetAlias("els_dPhiOut"));
		els_dPhiOut_branch->SetAddress(&els_dPhiOut);
	}
	if(els_dPhiOut_branch == 0 ) {
	cout << "Branch els_dPhiOut does not exist." << endl;
	}
	els_e3x3_branch = 0;
	if (tree->GetAlias("els_e3x3") != 0) {
		els_e3x3_branch = tree->GetBranch(tree->GetAlias("els_e3x3"));
		els_e3x3_branch->SetAddress(&els_e3x3);
	}
	if(els_e3x3_branch == 0 ) {
	cout << "Branch els_e3x3 does not exist." << endl;
	}
	els_e5x5_branch = 0;
	if (tree->GetAlias("els_e5x5") != 0) {
		els_e5x5_branch = tree->GetBranch(tree->GetAlias("els_e5x5"));
		els_e5x5_branch->SetAddress(&els_e5x5);
	}
	if(els_e5x5_branch == 0 ) {
	cout << "Branch els_e5x5 does not exist." << endl;
	}
	els_eOverPIn_branch = 0;
	if (tree->GetAlias("els_eOverPIn") != 0) {
		els_eOverPIn_branch = tree->GetBranch(tree->GetAlias("els_eOverPIn"));
		els_eOverPIn_branch->SetAddress(&els_eOverPIn);
	}
	if(els_eOverPIn_branch == 0 ) {
	cout << "Branch els_eOverPIn does not exist." << endl;
	}
	els_eSeedOverPOut_branch = 0;
	if (tree->GetAlias("els_eSeedOverPOut") != 0) {
		els_eSeedOverPOut_branch = tree->GetBranch(tree->GetAlias("els_eSeedOverPOut"));
		els_eSeedOverPOut_branch->SetAddress(&els_eSeedOverPOut);
	}
	if(els_eSeedOverPOut_branch == 0 ) {
	cout << "Branch els_eSeedOverPOut does not exist." << endl;
	}
	els_etaErr_branch = 0;
	if (tree->GetAlias("els_etaErr") != 0) {
		els_etaErr_branch = tree->GetBranch(tree->GetAlias("els_etaErr"));
		els_etaErr_branch->SetAddress(&els_etaErr);
	}
	if(els_etaErr_branch == 0 ) {
	cout << "Branch els_etaErr does not exist." << endl;
	}
	els_fBrem_branch = 0;
	if (tree->GetAlias("els_fBrem") != 0) {
		els_fBrem_branch = tree->GetBranch(tree->GetAlias("els_fBrem"));
		els_fBrem_branch->SetAddress(&els_fBrem);
	}
	if(els_fBrem_branch == 0 ) {
	cout << "Branch els_fBrem does not exist." << endl;
	}
	els_hOverE_branch = 0;
	if (tree->GetAlias("els_hOverE") != 0) {
		els_hOverE_branch = tree->GetBranch(tree->GetAlias("els_hOverE"));
		els_hOverE_branch->SetAddress(&els_hOverE);
	}
	if(els_hOverE_branch == 0 ) {
	cout << "Branch els_hOverE does not exist." << endl;
	}
	els_ndof_branch = 0;
	if (tree->GetAlias("els_ndof") != 0) {
		els_ndof_branch = tree->GetBranch(tree->GetAlias("els_ndof"));
		els_ndof_branch->SetAddress(&els_ndof);
	}
	if(els_ndof_branch == 0 ) {
	cout << "Branch els_ndof does not exist." << endl;
	}
	els_outerEta_branch = 0;
	if (tree->GetAlias("els_outerEta") != 0) {
		els_outerEta_branch = tree->GetBranch(tree->GetAlias("els_outerEta"));
		els_outerEta_branch->SetAddress(&els_outerEta);
	}
	if(els_outerEta_branch == 0 ) {
	cout << "Branch els_outerEta does not exist." << endl;
	}
	els_outerPhi_branch = 0;
	if (tree->GetAlias("els_outerPhi") != 0) {
		els_outerPhi_branch = tree->GetBranch(tree->GetAlias("els_outerPhi"));
		els_outerPhi_branch->SetAddress(&els_outerPhi);
	}
	if(els_outerPhi_branch == 0 ) {
	cout << "Branch els_outerPhi does not exist." << endl;
	}
	els_phiErr_branch = 0;
	if (tree->GetAlias("els_phiErr") != 0) {
		els_phiErr_branch = tree->GetBranch(tree->GetAlias("els_phiErr"));
		els_phiErr_branch->SetAddress(&els_phiErr);
	}
	if(els_phiErr_branch == 0 ) {
	cout << "Branch els_phiErr does not exist." << endl;
	}
	els_ptErr_branch = 0;
	if (tree->GetAlias("els_ptErr") != 0) {
		els_ptErr_branch = tree->GetBranch(tree->GetAlias("els_ptErr"));
		els_ptErr_branch->SetAddress(&els_ptErr);
	}
	if(els_ptErr_branch == 0 ) {
	cout << "Branch els_ptErr does not exist." << endl;
	}
	els_sigmaEtaEta_branch = 0;
	if (tree->GetAlias("els_sigmaEtaEta") != 0) {
		els_sigmaEtaEta_branch = tree->GetBranch(tree->GetAlias("els_sigmaEtaEta"));
		els_sigmaEtaEta_branch->SetAddress(&els_sigmaEtaEta);
	}
	if(els_sigmaEtaEta_branch == 0 ) {
	cout << "Branch els_sigmaEtaEta does not exist." << endl;
	}
	els_sigmaPhiPhi_branch = 0;
	if (tree->GetAlias("els_sigmaPhiPhi") != 0) {
		els_sigmaPhiPhi_branch = tree->GetBranch(tree->GetAlias("els_sigmaPhiPhi"));
		els_sigmaPhiPhi_branch->SetAddress(&els_sigmaPhiPhi);
	}
	if(els_sigmaPhiPhi_branch == 0 ) {
	cout << "Branch els_sigmaPhiPhi does not exist." << endl;
	}
	els_tkIso_branch = 0;
	if (tree->GetAlias("els_tkIso") != 0) {
		els_tkIso_branch = tree->GetBranch(tree->GetAlias("els_tkIso"));
		els_tkIso_branch->SetAddress(&els_tkIso);
	}
	if(els_tkIso_branch == 0 ) {
	cout << "Branch els_tkIso does not exist." << endl;
	}
	els_vertexphi_branch = 0;
	if (tree->GetAlias("els_vertexphi") != 0) {
		els_vertexphi_branch = tree->GetBranch(tree->GetAlias("els_vertexphi"));
		els_vertexphi_branch->SetAddress(&els_vertexphi);
	}
	if(els_vertexphi_branch == 0 ) {
	cout << "Branch els_vertexphi does not exist." << endl;
	}
	els_z0_branch = 0;
	if (tree->GetAlias("els_z0") != 0) {
		els_z0_branch = tree->GetBranch(tree->GetAlias("els_z0"));
		els_z0_branch->SetAddress(&els_z0);
	}
	if(els_z0_branch == 0 ) {
	cout << "Branch els_z0 does not exist." << endl;
	}
	els_z0Err_branch = 0;
	if (tree->GetAlias("els_z0Err") != 0) {
		els_z0Err_branch = tree->GetBranch(tree->GetAlias("els_z0Err"));
		els_z0Err_branch->SetAddress(&els_z0Err);
	}
	if(els_z0Err_branch == 0 ) {
	cout << "Branch els_z0Err does not exist." << endl;
	}
	hyp_ll_chi2_branch = 0;
	if (tree->GetAlias("hyp_ll_chi2") != 0) {
		hyp_ll_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_ll_chi2"));
		hyp_ll_chi2_branch->SetAddress(&hyp_ll_chi2);
	}
	if(hyp_ll_chi2_branch == 0 ) {
	cout << "Branch hyp_ll_chi2 does not exist." << endl;
	}
	hyp_ll_d0_branch = 0;
	if (tree->GetAlias("hyp_ll_d0") != 0) {
		hyp_ll_d0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0"));
		hyp_ll_d0_branch->SetAddress(&hyp_ll_d0);
	}
	if(hyp_ll_d0_branch == 0 ) {
	cout << "Branch hyp_ll_d0 does not exist." << endl;
	}
	hyp_ll_d0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_d0Err") != 0) {
		hyp_ll_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_d0Err"));
		hyp_ll_d0Err_branch->SetAddress(&hyp_ll_d0Err);
	}
	if(hyp_ll_d0Err_branch == 0 ) {
	cout << "Branch hyp_ll_d0Err does not exist." << endl;
	}
	hyp_ll_etaErr_branch = 0;
	if (tree->GetAlias("hyp_ll_etaErr") != 0) {
		hyp_ll_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_etaErr"));
		hyp_ll_etaErr_branch->SetAddress(&hyp_ll_etaErr);
	}
	if(hyp_ll_etaErr_branch == 0 ) {
	cout << "Branch hyp_ll_etaErr does not exist." << endl;
	}
	hyp_ll_iso_branch = 0;
	if (tree->GetAlias("hyp_ll_iso") != 0) {
		hyp_ll_iso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_iso"));
		hyp_ll_iso_branch->SetAddress(&hyp_ll_iso);
	}
	if(hyp_ll_iso_branch == 0 ) {
	cout << "Branch hyp_ll_iso does not exist." << endl;
	}
	hyp_ll_ndof_branch = 0;
	if (tree->GetAlias("hyp_ll_ndof") != 0) {
		hyp_ll_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ndof"));
		hyp_ll_ndof_branch->SetAddress(&hyp_ll_ndof);
	}
	if(hyp_ll_ndof_branch == 0 ) {
	cout << "Branch hyp_ll_ndof does not exist." << endl;
	}
	hyp_ll_outerEta_branch = 0;
	if (tree->GetAlias("hyp_ll_outerEta") != 0) {
		hyp_ll_outerEta_branch = tree->GetBranch(tree->GetAlias("hyp_ll_outerEta"));
		hyp_ll_outerEta_branch->SetAddress(&hyp_ll_outerEta);
	}
	if(hyp_ll_outerEta_branch == 0 ) {
	cout << "Branch hyp_ll_outerEta does not exist." << endl;
	}
	hyp_ll_outerPhi_branch = 0;
	if (tree->GetAlias("hyp_ll_outerPhi") != 0) {
		hyp_ll_outerPhi_branch = tree->GetBranch(tree->GetAlias("hyp_ll_outerPhi"));
		hyp_ll_outerPhi_branch->SetAddress(&hyp_ll_outerPhi);
	}
	if(hyp_ll_outerPhi_branch == 0 ) {
	cout << "Branch hyp_ll_outerPhi does not exist." << endl;
	}
	hyp_ll_phiErr_branch = 0;
	if (tree->GetAlias("hyp_ll_phiErr") != 0) {
		hyp_ll_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_phiErr"));
		hyp_ll_phiErr_branch->SetAddress(&hyp_ll_phiErr);
	}
	if(hyp_ll_phiErr_branch == 0 ) {
	cout << "Branch hyp_ll_phiErr does not exist." << endl;
	}
	hyp_ll_ptErr_branch = 0;
	if (tree->GetAlias("hyp_ll_ptErr") != 0) {
		hyp_ll_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_ll_ptErr"));
		hyp_ll_ptErr_branch->SetAddress(&hyp_ll_ptErr);
	}
	if(hyp_ll_ptErr_branch == 0 ) {
	cout << "Branch hyp_ll_ptErr does not exist." << endl;
	}
	hyp_ll_tkIso_branch = 0;
	if (tree->GetAlias("hyp_ll_tkIso") != 0) {
		hyp_ll_tkIso_branch = tree->GetBranch(tree->GetAlias("hyp_ll_tkIso"));
		hyp_ll_tkIso_branch->SetAddress(&hyp_ll_tkIso);
	}
	if(hyp_ll_tkIso_branch == 0 ) {
	cout << "Branch hyp_ll_tkIso does not exist." << endl;
	}
	hyp_ll_vertexphi_branch = 0;
	if (tree->GetAlias("hyp_ll_vertexphi") != 0) {
		hyp_ll_vertexphi_branch = tree->GetBranch(tree->GetAlias("hyp_ll_vertexphi"));
		hyp_ll_vertexphi_branch->SetAddress(&hyp_ll_vertexphi);
	}
	if(hyp_ll_vertexphi_branch == 0 ) {
	cout << "Branch hyp_ll_vertexphi does not exist." << endl;
	}
	hyp_ll_z0_branch = 0;
	if (tree->GetAlias("hyp_ll_z0") != 0) {
		hyp_ll_z0_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0"));
		hyp_ll_z0_branch->SetAddress(&hyp_ll_z0);
	}
	if(hyp_ll_z0_branch == 0 ) {
	cout << "Branch hyp_ll_z0 does not exist." << endl;
	}
	hyp_ll_z0Err_branch = 0;
	if (tree->GetAlias("hyp_ll_z0Err") != 0) {
		hyp_ll_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_ll_z0Err"));
		hyp_ll_z0Err_branch->SetAddress(&hyp_ll_z0Err);
	}
	if(hyp_ll_z0Err_branch == 0 ) {
	cout << "Branch hyp_ll_z0Err does not exist." << endl;
	}
	hyp_lt_chi2_branch = 0;
	if (tree->GetAlias("hyp_lt_chi2") != 0) {
		hyp_lt_chi2_branch = tree->GetBranch(tree->GetAlias("hyp_lt_chi2"));
		hyp_lt_chi2_branch->SetAddress(&hyp_lt_chi2);
	}
	if(hyp_lt_chi2_branch == 0 ) {
	cout << "Branch hyp_lt_chi2 does not exist." << endl;
	}
	hyp_lt_d0_branch = 0;
	if (tree->GetAlias("hyp_lt_d0") != 0) {
		hyp_lt_d0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0"));
		hyp_lt_d0_branch->SetAddress(&hyp_lt_d0);
	}
	if(hyp_lt_d0_branch == 0 ) {
	cout << "Branch hyp_lt_d0 does not exist." << endl;
	}
	hyp_lt_d0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_d0Err") != 0) {
		hyp_lt_d0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_d0Err"));
		hyp_lt_d0Err_branch->SetAddress(&hyp_lt_d0Err);
	}
	if(hyp_lt_d0Err_branch == 0 ) {
	cout << "Branch hyp_lt_d0Err does not exist." << endl;
	}
	hyp_lt_etaErr_branch = 0;
	if (tree->GetAlias("hyp_lt_etaErr") != 0) {
		hyp_lt_etaErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_etaErr"));
		hyp_lt_etaErr_branch->SetAddress(&hyp_lt_etaErr);
	}
	if(hyp_lt_etaErr_branch == 0 ) {
	cout << "Branch hyp_lt_etaErr does not exist." << endl;
	}
	hyp_lt_iso_branch = 0;
	if (tree->GetAlias("hyp_lt_iso") != 0) {
		hyp_lt_iso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_iso"));
		hyp_lt_iso_branch->SetAddress(&hyp_lt_iso);
	}
	if(hyp_lt_iso_branch == 0 ) {
	cout << "Branch hyp_lt_iso does not exist." << endl;
	}
	hyp_lt_ndof_branch = 0;
	if (tree->GetAlias("hyp_lt_ndof") != 0) {
		hyp_lt_ndof_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ndof"));
		hyp_lt_ndof_branch->SetAddress(&hyp_lt_ndof);
	}
	if(hyp_lt_ndof_branch == 0 ) {
	cout << "Branch hyp_lt_ndof does not exist." << endl;
	}
	hyp_lt_outerEta_branch = 0;
	if (tree->GetAlias("hyp_lt_outerEta") != 0) {
		hyp_lt_outerEta_branch = tree->GetBranch(tree->GetAlias("hyp_lt_outerEta"));
		hyp_lt_outerEta_branch->SetAddress(&hyp_lt_outerEta);
	}
	if(hyp_lt_outerEta_branch == 0 ) {
	cout << "Branch hyp_lt_outerEta does not exist." << endl;
	}
	hyp_lt_outerPhi_branch = 0;
	if (tree->GetAlias("hyp_lt_outerPhi") != 0) {
		hyp_lt_outerPhi_branch = tree->GetBranch(tree->GetAlias("hyp_lt_outerPhi"));
		hyp_lt_outerPhi_branch->SetAddress(&hyp_lt_outerPhi);
	}
	if(hyp_lt_outerPhi_branch == 0 ) {
	cout << "Branch hyp_lt_outerPhi does not exist." << endl;
	}
	hyp_lt_phiErr_branch = 0;
	if (tree->GetAlias("hyp_lt_phiErr") != 0) {
		hyp_lt_phiErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_phiErr"));
		hyp_lt_phiErr_branch->SetAddress(&hyp_lt_phiErr);
	}
	if(hyp_lt_phiErr_branch == 0 ) {
	cout << "Branch hyp_lt_phiErr does not exist." << endl;
	}
	hyp_lt_ptErr_branch = 0;
	if (tree->GetAlias("hyp_lt_ptErr") != 0) {
		hyp_lt_ptErr_branch = tree->GetBranch(tree->GetAlias("hyp_lt_ptErr"));
		hyp_lt_ptErr_branch->SetAddress(&hyp_lt_ptErr);
	}
	if(hyp_lt_ptErr_branch == 0 ) {
	cout << "Branch hyp_lt_ptErr does not exist." << endl;
	}
	hyp_lt_tkIso_branch = 0;
	if (tree->GetAlias("hyp_lt_tkIso") != 0) {
		hyp_lt_tkIso_branch = tree->GetBranch(tree->GetAlias("hyp_lt_tkIso"));
		hyp_lt_tkIso_branch->SetAddress(&hyp_lt_tkIso);
	}
	if(hyp_lt_tkIso_branch == 0 ) {
	cout << "Branch hyp_lt_tkIso does not exist." << endl;
	}
	hyp_lt_vertexphi_branch = 0;
	if (tree->GetAlias("hyp_lt_vertexphi") != 0) {
		hyp_lt_vertexphi_branch = tree->GetBranch(tree->GetAlias("hyp_lt_vertexphi"));
		hyp_lt_vertexphi_branch->SetAddress(&hyp_lt_vertexphi);
	}
	if(hyp_lt_vertexphi_branch == 0 ) {
	cout << "Branch hyp_lt_vertexphi does not exist." << endl;
	}
	hyp_lt_z0_branch = 0;
	if (tree->GetAlias("hyp_lt_z0") != 0) {
		hyp_lt_z0_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0"));
		hyp_lt_z0_branch->SetAddress(&hyp_lt_z0);
	}
	if(hyp_lt_z0_branch == 0 ) {
	cout << "Branch hyp_lt_z0 does not exist." << endl;
	}
	hyp_lt_z0Err_branch = 0;
	if (tree->GetAlias("hyp_lt_z0Err") != 0) {
		hyp_lt_z0Err_branch = tree->GetBranch(tree->GetAlias("hyp_lt_z0Err"));
		hyp_lt_z0Err_branch->SetAddress(&hyp_lt_z0Err);
	}
	if(hyp_lt_z0Err_branch == 0 ) {
	cout << "Branch hyp_lt_z0Err does not exist." << endl;
	}
	hyp_met_branch = 0;
	if (tree->GetAlias("hyp_met") != 0) {
		hyp_met_branch = tree->GetBranch(tree->GetAlias("hyp_met"));
		hyp_met_branch->SetAddress(&hyp_met);
	}
	if(hyp_met_branch == 0 ) {
	cout << "Branch hyp_met does not exist." << endl;
	}
	hyp_metAll_branch = 0;
	if (tree->GetAlias("hyp_metAll") != 0) {
		hyp_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_metAll"));
		hyp_metAll_branch->SetAddress(&hyp_metAll);
	}
	if(hyp_metAll_branch == 0 ) {
	cout << "Branch hyp_metAll does not exist." << endl;
	}
	hyp_metAllCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metAllCaloExp") != 0) {
		hyp_metAllCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metAllCaloExp"));
		hyp_metAllCaloExp_branch->SetAddress(&hyp_metAllCaloExp);
	}
	if(hyp_metAllCaloExp_branch == 0 ) {
	cout << "Branch hyp_metAllCaloExp does not exist." << endl;
	}
	hyp_metCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metCaloExp") != 0) {
		hyp_metCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metCaloExp"));
		hyp_metCaloExp_branch->SetAddress(&hyp_metCaloExp);
	}
	if(hyp_metCaloExp_branch == 0 ) {
	cout << "Branch hyp_metCaloExp does not exist." << endl;
	}
	hyp_metCone_branch = 0;
	if (tree->GetAlias("hyp_metCone") != 0) {
		hyp_metCone_branch = tree->GetBranch(tree->GetAlias("hyp_metCone"));
		hyp_metCone_branch->SetAddress(&hyp_metCone);
	}
	if(hyp_metCone_branch == 0 ) {
	cout << "Branch hyp_metCone does not exist." << endl;
	}
	hyp_metDPhiJet10_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet10") != 0) {
		hyp_metDPhiJet10_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet10"));
		hyp_metDPhiJet10_branch->SetAddress(&hyp_metDPhiJet10);
	}
	if(hyp_metDPhiJet10_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet10 does not exist." << endl;
	}
	hyp_metDPhiJet15_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet15") != 0) {
		hyp_metDPhiJet15_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet15"));
		hyp_metDPhiJet15_branch->SetAddress(&hyp_metDPhiJet15);
	}
	if(hyp_metDPhiJet15_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet15 does not exist." << endl;
	}
	hyp_metDPhiJet20_branch = 0;
	if (tree->GetAlias("hyp_metDPhiJet20") != 0) {
		hyp_metDPhiJet20_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiJet20"));
		hyp_metDPhiJet20_branch->SetAddress(&hyp_metDPhiJet20);
	}
	if(hyp_metDPhiJet20_branch == 0 ) {
	cout << "Branch hyp_metDPhiJet20 does not exist." << endl;
	}
	hyp_metDPhiTrk10_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk10") != 0) {
		hyp_metDPhiTrk10_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk10"));
		hyp_metDPhiTrk10_branch->SetAddress(&hyp_metDPhiTrk10);
	}
	if(hyp_metDPhiTrk10_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk10 does not exist." << endl;
	}
	hyp_metDPhiTrk25_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk25") != 0) {
		hyp_metDPhiTrk25_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk25"));
		hyp_metDPhiTrk25_branch->SetAddress(&hyp_metDPhiTrk25);
	}
	if(hyp_metDPhiTrk25_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk25 does not exist." << endl;
	}
	hyp_metDPhiTrk50_branch = 0;
	if (tree->GetAlias("hyp_metDPhiTrk50") != 0) {
		hyp_metDPhiTrk50_branch = tree->GetBranch(tree->GetAlias("hyp_metDPhiTrk50"));
		hyp_metDPhiTrk50_branch->SetAddress(&hyp_metDPhiTrk50);
	}
	if(hyp_metDPhiTrk50_branch == 0 ) {
	cout << "Branch hyp_metDPhiTrk50 does not exist." << endl;
	}
	hyp_metJes10_branch = 0;
	if (tree->GetAlias("hyp_metJes10") != 0) {
		hyp_metJes10_branch = tree->GetBranch(tree->GetAlias("hyp_metJes10"));
		hyp_metJes10_branch->SetAddress(&hyp_metJes10);
	}
	if(hyp_metJes10_branch == 0 ) {
	cout << "Branch hyp_metJes10 does not exist." << endl;
	}
	hyp_metJes15_branch = 0;
	if (tree->GetAlias("hyp_metJes15") != 0) {
		hyp_metJes15_branch = tree->GetBranch(tree->GetAlias("hyp_metJes15"));
		hyp_metJes15_branch->SetAddress(&hyp_metJes15);
	}
	if(hyp_metJes15_branch == 0 ) {
	cout << "Branch hyp_metJes15 does not exist." << endl;
	}
	hyp_metJes30_branch = 0;
	if (tree->GetAlias("hyp_metJes30") != 0) {
		hyp_metJes30_branch = tree->GetBranch(tree->GetAlias("hyp_metJes30"));
		hyp_metJes30_branch->SetAddress(&hyp_metJes30);
	}
	if(hyp_metJes30_branch == 0 ) {
	cout << "Branch hyp_metJes30 does not exist." << endl;
	}
	hyp_metJes5_branch = 0;
	if (tree->GetAlias("hyp_metJes5") != 0) {
		hyp_metJes5_branch = tree->GetBranch(tree->GetAlias("hyp_metJes5"));
		hyp_metJes5_branch->SetAddress(&hyp_metJes5);
	}
	if(hyp_metJes5_branch == 0 ) {
	cout << "Branch hyp_metJes5 does not exist." << endl;
	}
	hyp_metJes50_branch = 0;
	if (tree->GetAlias("hyp_metJes50") != 0) {
		hyp_metJes50_branch = tree->GetBranch(tree->GetAlias("hyp_metJes50"));
		hyp_metJes50_branch->SetAddress(&hyp_metJes50);
	}
	if(hyp_metJes50_branch == 0 ) {
	cout << "Branch hyp_metJes50 does not exist." << endl;
	}
	hyp_metNoCalo_branch = 0;
	if (tree->GetAlias("hyp_metNoCalo") != 0) {
		hyp_metNoCalo_branch = tree->GetBranch(tree->GetAlias("hyp_metNoCalo"));
		hyp_metNoCalo_branch->SetAddress(&hyp_metNoCalo);
	}
	if(hyp_metNoCalo_branch == 0 ) {
	cout << "Branch hyp_metNoCalo does not exist." << endl;
	}
	hyp_metPhi_branch = 0;
	if (tree->GetAlias("hyp_metPhi") != 0) {
		hyp_metPhi_branch = tree->GetBranch(tree->GetAlias("hyp_metPhi"));
		hyp_metPhi_branch->SetAddress(&hyp_metPhi);
	}
	if(hyp_metPhi_branch == 0 ) {
	cout << "Branch hyp_metPhi does not exist." << endl;
	}
	hyp_metPhiAll_branch = 0;
	if (tree->GetAlias("hyp_metPhiAll") != 0) {
		hyp_metPhiAll_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiAll"));
		hyp_metPhiAll_branch->SetAddress(&hyp_metPhiAll);
	}
	if(hyp_metPhiAll_branch == 0 ) {
	cout << "Branch hyp_metPhiAll does not exist." << endl;
	}
	hyp_metPhiAllCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metPhiAllCaloExp") != 0) {
		hyp_metPhiAllCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiAllCaloExp"));
		hyp_metPhiAllCaloExp_branch->SetAddress(&hyp_metPhiAllCaloExp);
	}
	if(hyp_metPhiAllCaloExp_branch == 0 ) {
	cout << "Branch hyp_metPhiAllCaloExp does not exist." << endl;
	}
	hyp_metPhiCaloExp_branch = 0;
	if (tree->GetAlias("hyp_metPhiCaloExp") != 0) {
		hyp_metPhiCaloExp_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiCaloExp"));
		hyp_metPhiCaloExp_branch->SetAddress(&hyp_metPhiCaloExp);
	}
	if(hyp_metPhiCaloExp_branch == 0 ) {
	cout << "Branch hyp_metPhiCaloExp does not exist." << endl;
	}
	hyp_metPhiCone_branch = 0;
	if (tree->GetAlias("hyp_metPhiCone") != 0) {
		hyp_metPhiCone_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiCone"));
		hyp_metPhiCone_branch->SetAddress(&hyp_metPhiCone);
	}
	if(hyp_metPhiCone_branch == 0 ) {
	cout << "Branch hyp_metPhiCone does not exist." << endl;
	}
	hyp_metPhiJes10_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes10") != 0) {
		hyp_metPhiJes10_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes10"));
		hyp_metPhiJes10_branch->SetAddress(&hyp_metPhiJes10);
	}
	if(hyp_metPhiJes10_branch == 0 ) {
	cout << "Branch hyp_metPhiJes10 does not exist." << endl;
	}
	hyp_metPhiJes15_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes15") != 0) {
		hyp_metPhiJes15_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes15"));
		hyp_metPhiJes15_branch->SetAddress(&hyp_metPhiJes15);
	}
	if(hyp_metPhiJes15_branch == 0 ) {
	cout << "Branch hyp_metPhiJes15 does not exist." << endl;
	}
	hyp_metPhiJes30_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes30") != 0) {
		hyp_metPhiJes30_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes30"));
		hyp_metPhiJes30_branch->SetAddress(&hyp_metPhiJes30);
	}
	if(hyp_metPhiJes30_branch == 0 ) {
	cout << "Branch hyp_metPhiJes30 does not exist." << endl;
	}
	hyp_metPhiJes5_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes5") != 0) {
		hyp_metPhiJes5_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes5"));
		hyp_metPhiJes5_branch->SetAddress(&hyp_metPhiJes5);
	}
	if(hyp_metPhiJes5_branch == 0 ) {
	cout << "Branch hyp_metPhiJes5 does not exist." << endl;
	}
	hyp_metPhiJes50_branch = 0;
	if (tree->GetAlias("hyp_metPhiJes50") != 0) {
		hyp_metPhiJes50_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiJes50"));
		hyp_metPhiJes50_branch->SetAddress(&hyp_metPhiJes50);
	}
	if(hyp_metPhiJes50_branch == 0 ) {
	cout << "Branch hyp_metPhiJes50 does not exist." << endl;
	}
	hyp_metPhiNoCalo_branch = 0;
	if (tree->GetAlias("hyp_metPhiNoCalo") != 0) {
		hyp_metPhiNoCalo_branch = tree->GetBranch(tree->GetAlias("hyp_metPhiNoCalo"));
		hyp_metPhiNoCalo_branch->SetAddress(&hyp_metPhiNoCalo);
	}
	if(hyp_metPhiNoCalo_branch == 0 ) {
	cout << "Branch hyp_metPhiNoCalo does not exist." << endl;
	}
	hyp_quadlep_met_branch = 0;
	if (tree->GetAlias("hyp_quadlep_met") != 0) {
		hyp_quadlep_met_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_met"));
		hyp_quadlep_met_branch->SetAddress(&hyp_quadlep_met);
	}
	if(hyp_quadlep_met_branch == 0 ) {
	cout << "Branch hyp_quadlep_met does not exist." << endl;
	}
	hyp_quadlep_metAll_branch = 0;
	if (tree->GetAlias("hyp_quadlep_metAll") != 0) {
		hyp_quadlep_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_metAll"));
		hyp_quadlep_metAll_branch->SetAddress(&hyp_quadlep_metAll);
	}
	if(hyp_quadlep_metAll_branch == 0 ) {
	cout << "Branch hyp_quadlep_metAll does not exist." << endl;
	}
	hyp_trilep_met_branch = 0;
	if (tree->GetAlias("hyp_trilep_met") != 0) {
		hyp_trilep_met_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_met"));
		hyp_trilep_met_branch->SetAddress(&hyp_trilep_met);
	}
	if(hyp_trilep_met_branch == 0 ) {
	cout << "Branch hyp_trilep_met does not exist." << endl;
	}
	hyp_trilep_metAll_branch = 0;
	if (tree->GetAlias("hyp_trilep_metAll") != 0) {
		hyp_trilep_metAll_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_metAll"));
		hyp_trilep_metAll_branch->SetAddress(&hyp_trilep_metAll);
	}
	if(hyp_trilep_metAll_branch == 0 ) {
	cout << "Branch hyp_trilep_metAll does not exist." << endl;
	}
	jets_EMFcor_branch = 0;
	if (tree->GetAlias("jets_EMFcor") != 0) {
		jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("jets_EMFcor"));
		jets_EMFcor_branch->SetAddress(&jets_EMFcor);
	}
	if(jets_EMFcor_branch == 0 ) {
	cout << "Branch jets_EMFcor does not exist." << endl;
	}
	jets_chFrac_branch = 0;
	if (tree->GetAlias("jets_chFrac") != 0) {
		jets_chFrac_branch = tree->GetBranch(tree->GetAlias("jets_chFrac"));
		jets_chFrac_branch->SetAddress(&jets_chFrac);
	}
	if(jets_chFrac_branch == 0 ) {
	cout << "Branch jets_chFrac does not exist." << endl;
	}
	jets_cor_branch = 0;
	if (tree->GetAlias("jets_cor") != 0) {
		jets_cor_branch = tree->GetBranch(tree->GetAlias("jets_cor"));
		jets_cor_branch->SetAddress(&jets_cor);
	}
	if(jets_cor_branch == 0 ) {
	cout << "Branch jets_cor does not exist." << endl;
	}
	jets_emFrac_branch = 0;
	if (tree->GetAlias("jets_emFrac") != 0) {
		jets_emFrac_branch = tree->GetBranch(tree->GetAlias("jets_emFrac"));
		jets_emFrac_branch->SetAddress(&jets_emFrac);
	}
	if(jets_emFrac_branch == 0 ) {
	cout << "Branch jets_emFrac does not exist." << endl;
	}
	mus_eledr_branch = 0;
	if (tree->GetAlias("mus_eledr") != 0) {
		mus_eledr_branch = tree->GetBranch(tree->GetAlias("mus_eledr"));
		mus_eledr_branch->SetAddress(&mus_eledr);
	}
	if(mus_eledr_branch == 0 ) {
	cout << "Branch mus_eledr does not exist." << endl;
	}
	mus_jetdr_branch = 0;
	if (tree->GetAlias("mus_jetdr") != 0) {
		mus_jetdr_branch = tree->GetBranch(tree->GetAlias("mus_jetdr"));
		mus_jetdr_branch->SetAddress(&mus_jetdr);
	}
	if(mus_jetdr_branch == 0 ) {
	cout << "Branch mus_jetdr does not exist." << endl;
	}
	mus_trkdr_branch = 0;
	if (tree->GetAlias("mus_trkdr") != 0) {
		mus_trkdr_branch = tree->GetBranch(tree->GetAlias("mus_trkdr"));
		mus_trkdr_branch->SetAddress(&mus_trkdr);
	}
	if(mus_trkdr_branch == 0 ) {
	cout << "Branch mus_trkdr does not exist." << endl;
	}
	mus_chi2_branch = 0;
	if (tree->GetAlias("mus_chi2") != 0) {
		mus_chi2_branch = tree->GetBranch(tree->GetAlias("mus_chi2"));
		mus_chi2_branch->SetAddress(&mus_chi2);
	}
	if(mus_chi2_branch == 0 ) {
	cout << "Branch mus_chi2 does not exist." << endl;
	}
	mus_d0_branch = 0;
	if (tree->GetAlias("mus_d0") != 0) {
		mus_d0_branch = tree->GetBranch(tree->GetAlias("mus_d0"));
		mus_d0_branch->SetAddress(&mus_d0);
	}
	if(mus_d0_branch == 0 ) {
	cout << "Branch mus_d0 does not exist." << endl;
	}
	mus_d0Err_branch = 0;
	if (tree->GetAlias("mus_d0Err") != 0) {
		mus_d0Err_branch = tree->GetBranch(tree->GetAlias("mus_d0Err"));
		mus_d0Err_branch->SetAddress(&mus_d0Err);
	}
	if(mus_d0Err_branch == 0 ) {
	cout << "Branch mus_d0Err does not exist." << endl;
	}
	mus_e_em_branch = 0;
	if (tree->GetAlias("mus_e_em") != 0) {
		mus_e_em_branch = tree->GetBranch(tree->GetAlias("mus_e_em"));
		mus_e_em_branch->SetAddress(&mus_e_em);
	}
	if(mus_e_em_branch == 0 ) {
	cout << "Branch mus_e_em does not exist." << endl;
	}
	mus_e_emS9_branch = 0;
	if (tree->GetAlias("mus_e_emS9") != 0) {
		mus_e_emS9_branch = tree->GetBranch(tree->GetAlias("mus_e_emS9"));
		mus_e_emS9_branch->SetAddress(&mus_e_emS9);
	}
	if(mus_e_emS9_branch == 0 ) {
	cout << "Branch mus_e_emS9 does not exist." << endl;
	}
	mus_e_had_branch = 0;
	if (tree->GetAlias("mus_e_had") != 0) {
		mus_e_had_branch = tree->GetBranch(tree->GetAlias("mus_e_had"));
		mus_e_had_branch->SetAddress(&mus_e_had);
	}
	if(mus_e_had_branch == 0 ) {
	cout << "Branch mus_e_had does not exist." << endl;
	}
	mus_e_hadS9_branch = 0;
	if (tree->GetAlias("mus_e_hadS9") != 0) {
		mus_e_hadS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hadS9"));
		mus_e_hadS9_branch->SetAddress(&mus_e_hadS9);
	}
	if(mus_e_hadS9_branch == 0 ) {
	cout << "Branch mus_e_hadS9 does not exist." << endl;
	}
	mus_e_ho_branch = 0;
	if (tree->GetAlias("mus_e_ho") != 0) {
		mus_e_ho_branch = tree->GetBranch(tree->GetAlias("mus_e_ho"));
		mus_e_ho_branch->SetAddress(&mus_e_ho);
	}
	if(mus_e_ho_branch == 0 ) {
	cout << "Branch mus_e_ho does not exist." << endl;
	}
	mus_e_hoS9_branch = 0;
	if (tree->GetAlias("mus_e_hoS9") != 0) {
		mus_e_hoS9_branch = tree->GetBranch(tree->GetAlias("mus_e_hoS9"));
		mus_e_hoS9_branch->SetAddress(&mus_e_hoS9);
	}
	if(mus_e_hoS9_branch == 0 ) {
	cout << "Branch mus_e_hoS9 does not exist." << endl;
	}
	mus_etaErr_branch = 0;
	if (tree->GetAlias("mus_etaErr") != 0) {
		mus_etaErr_branch = tree->GetBranch(tree->GetAlias("mus_etaErr"));
		mus_etaErr_branch->SetAddress(&mus_etaErr);
	}
	if(mus_etaErr_branch == 0 ) {
	cout << "Branch mus_etaErr does not exist." << endl;
	}
	mus_gfit_chi2_branch = 0;
	if (tree->GetAlias("mus_gfit_chi2") != 0) {
		mus_gfit_chi2_branch = tree->GetBranch(tree->GetAlias("mus_gfit_chi2"));
		mus_gfit_chi2_branch->SetAddress(&mus_gfit_chi2);
	}
	if(mus_gfit_chi2_branch == 0 ) {
	cout << "Branch mus_gfit_chi2 does not exist." << endl;
	}
	mus_gfit_ndof_branch = 0;
	if (tree->GetAlias("mus_gfit_ndof") != 0) {
		mus_gfit_ndof_branch = tree->GetBranch(tree->GetAlias("mus_gfit_ndof"));
		mus_gfit_ndof_branch->SetAddress(&mus_gfit_ndof);
	}
	if(mus_gfit_ndof_branch == 0 ) {
	cout << "Branch mus_gfit_ndof does not exist." << endl;
	}
	mus_iso_branch = 0;
	if (tree->GetAlias("mus_iso") != 0) {
		mus_iso_branch = tree->GetBranch(tree->GetAlias("mus_iso"));
		mus_iso_branch->SetAddress(&mus_iso);
	}
	if(mus_iso_branch == 0 ) {
	cout << "Branch mus_iso does not exist." << endl;
	}
	mus_iso03_emEt_branch = 0;
	if (tree->GetAlias("mus_iso03_emEt") != 0) {
		mus_iso03_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_emEt"));
		mus_iso03_emEt_branch->SetAddress(&mus_iso03_emEt);
	}
	if(mus_iso03_emEt_branch == 0 ) {
	cout << "Branch mus_iso03_emEt does not exist." << endl;
	}
	mus_iso03_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hadEt") != 0) {
		mus_iso03_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hadEt"));
		mus_iso03_hadEt_branch->SetAddress(&mus_iso03_hadEt);
	}
	if(mus_iso03_hadEt_branch == 0 ) {
	cout << "Branch mus_iso03_hadEt does not exist." << endl;
	}
	mus_iso03_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso03_hoEt") != 0) {
		mus_iso03_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_hoEt"));
		mus_iso03_hoEt_branch->SetAddress(&mus_iso03_hoEt);
	}
	if(mus_iso03_hoEt_branch == 0 ) {
	cout << "Branch mus_iso03_hoEt does not exist." << endl;
	}
	mus_iso03_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso03_sumPt") != 0) {
		mus_iso03_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso03_sumPt"));
		mus_iso03_sumPt_branch->SetAddress(&mus_iso03_sumPt);
	}
	if(mus_iso03_sumPt_branch == 0 ) {
	cout << "Branch mus_iso03_sumPt does not exist." << endl;
	}
	mus_iso05_emEt_branch = 0;
	if (tree->GetAlias("mus_iso05_emEt") != 0) {
		mus_iso05_emEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_emEt"));
		mus_iso05_emEt_branch->SetAddress(&mus_iso05_emEt);
	}
	if(mus_iso05_emEt_branch == 0 ) {
	cout << "Branch mus_iso05_emEt does not exist." << endl;
	}
	mus_iso05_hadEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hadEt") != 0) {
		mus_iso05_hadEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hadEt"));
		mus_iso05_hadEt_branch->SetAddress(&mus_iso05_hadEt);
	}
	if(mus_iso05_hadEt_branch == 0 ) {
	cout << "Branch mus_iso05_hadEt does not exist." << endl;
	}
	mus_iso05_hoEt_branch = 0;
	if (tree->GetAlias("mus_iso05_hoEt") != 0) {
		mus_iso05_hoEt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_hoEt"));
		mus_iso05_hoEt_branch->SetAddress(&mus_iso05_hoEt);
	}
	if(mus_iso05_hoEt_branch == 0 ) {
	cout << "Branch mus_iso05_hoEt does not exist." << endl;
	}
	mus_iso05_sumPt_branch = 0;
	if (tree->GetAlias("mus_iso05_sumPt") != 0) {
		mus_iso05_sumPt_branch = tree->GetBranch(tree->GetAlias("mus_iso05_sumPt"));
		mus_iso05_sumPt_branch->SetAddress(&mus_iso05_sumPt);
	}
	if(mus_iso05_sumPt_branch == 0 ) {
	cout << "Branch mus_iso05_sumPt does not exist." << endl;
	}
	mus_ndof_branch = 0;
	if (tree->GetAlias("mus_ndof") != 0) {
		mus_ndof_branch = tree->GetBranch(tree->GetAlias("mus_ndof"));
		mus_ndof_branch->SetAddress(&mus_ndof);
	}
	if(mus_ndof_branch == 0 ) {
	cout << "Branch mus_ndof does not exist." << endl;
	}
	mus_outerEta_branch = 0;
	if (tree->GetAlias("mus_outerEta") != 0) {
		mus_outerEta_branch = tree->GetBranch(tree->GetAlias("mus_outerEta"));
		mus_outerEta_branch->SetAddress(&mus_outerEta);
	}
	if(mus_outerEta_branch == 0 ) {
	cout << "Branch mus_outerEta does not exist." << endl;
	}
	mus_outerPhi_branch = 0;
	if (tree->GetAlias("mus_outerPhi") != 0) {
		mus_outerPhi_branch = tree->GetBranch(tree->GetAlias("mus_outerPhi"));
		mus_outerPhi_branch->SetAddress(&mus_outerPhi);
	}
	if(mus_outerPhi_branch == 0 ) {
	cout << "Branch mus_outerPhi does not exist." << endl;
	}
	mus_phiErr_branch = 0;
	if (tree->GetAlias("mus_phiErr") != 0) {
		mus_phiErr_branch = tree->GetBranch(tree->GetAlias("mus_phiErr"));
		mus_phiErr_branch->SetAddress(&mus_phiErr);
	}
	if(mus_phiErr_branch == 0 ) {
	cout << "Branch mus_phiErr does not exist." << endl;
	}
	mus_ptErr_branch = 0;
	if (tree->GetAlias("mus_ptErr") != 0) {
		mus_ptErr_branch = tree->GetBranch(tree->GetAlias("mus_ptErr"));
		mus_ptErr_branch->SetAddress(&mus_ptErr);
	}
	if(mus_ptErr_branch == 0 ) {
	cout << "Branch mus_ptErr does not exist." << endl;
	}
	mus_vertexphi_branch = 0;
	if (tree->GetAlias("mus_vertexphi") != 0) {
		mus_vertexphi_branch = tree->GetBranch(tree->GetAlias("mus_vertexphi"));
		mus_vertexphi_branch->SetAddress(&mus_vertexphi);
	}
	if(mus_vertexphi_branch == 0 ) {
	cout << "Branch mus_vertexphi does not exist." << endl;
	}
	mus_z0_branch = 0;
	if (tree->GetAlias("mus_z0") != 0) {
		mus_z0_branch = tree->GetBranch(tree->GetAlias("mus_z0"));
		mus_z0_branch->SetAddress(&mus_z0);
	}
	if(mus_z0_branch == 0 ) {
	cout << "Branch mus_z0 does not exist." << endl;
	}
	mus_z0Err_branch = 0;
	if (tree->GetAlias("mus_z0Err") != 0) {
		mus_z0Err_branch = tree->GetBranch(tree->GetAlias("mus_z0Err"));
		mus_z0Err_branch->SetAddress(&mus_z0Err);
	}
	if(mus_z0Err_branch == 0 ) {
	cout << "Branch mus_z0Err does not exist." << endl;
	}
	els_tq_LRComb_branch = 0;
	if (tree->GetAlias("els_tq_LRComb") != 0) {
		els_tq_LRComb_branch = tree->GetBranch(tree->GetAlias("els_tq_LRComb"));
		els_tq_LRComb_branch->SetAddress(&els_tq_LRComb);
	}
	if(els_tq_LRComb_branch == 0 ) {
	cout << "Branch els_tq_LRComb does not exist." << endl;
	}
	els_tq_caloIso_branch = 0;
	if (tree->GetAlias("els_tq_caloIso") != 0) {
		els_tq_caloIso_branch = tree->GetBranch(tree->GetAlias("els_tq_caloIso"));
		els_tq_caloIso_branch->SetAddress(&els_tq_caloIso);
	}
	if(els_tq_caloIso_branch == 0 ) {
	cout << "Branch els_tq_caloIso does not exist." << endl;
	}
	els_tq_egammaEcalIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaEcalIso") != 0) {
		els_tq_egammaEcalIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaEcalIso"));
		els_tq_egammaEcalIso_branch->SetAddress(&els_tq_egammaEcalIso);
	}
	if(els_tq_egammaEcalIso_branch == 0 ) {
	cout << "Branch els_tq_egammaEcalIso does not exist." << endl;
	}
	els_tq_egammaHcalIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaHcalIso") != 0) {
		els_tq_egammaHcalIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaHcalIso"));
		els_tq_egammaHcalIso_branch->SetAddress(&els_tq_egammaHcalIso);
	}
	if(els_tq_egammaHcalIso_branch == 0 ) {
	cout << "Branch els_tq_egammaHcalIso does not exist." << endl;
	}
	els_tq_egammaTkIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaTkIso") != 0) {
		els_tq_egammaTkIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaTkIso"));
		els_tq_egammaTkIso_branch->SetAddress(&els_tq_egammaTkIso);
	}
	if(els_tq_egammaTkIso_branch == 0 ) {
	cout << "Branch els_tq_egammaTkIso does not exist." << endl;
	}
	els_tq_electronIDRobust_branch = 0;
	if (tree->GetAlias("els_tq_electronIDRobust") != 0) {
		els_tq_electronIDRobust_branch = tree->GetBranch(tree->GetAlias("els_tq_electronIDRobust"));
		els_tq_electronIDRobust_branch->SetAddress(&els_tq_electronIDRobust);
	}
	if(els_tq_electronIDRobust_branch == 0 ) {
	cout << "Branch els_tq_electronIDRobust does not exist." << endl;
	}
	els_tq_leptonID_branch = 0;
	if (tree->GetAlias("els_tq_leptonID") != 0) {
		els_tq_leptonID_branch = tree->GetBranch(tree->GetAlias("els_tq_leptonID"));
		els_tq_leptonID_branch->SetAddress(&els_tq_leptonID);
	}
	if(els_tq_leptonID_branch == 0 ) {
	cout << "Branch els_tq_leptonID does not exist." << endl;
	}
	els_tq_trackIso_branch = 0;
	if (tree->GetAlias("els_tq_trackIso") != 0) {
		els_tq_trackIso_branch = tree->GetBranch(tree->GetAlias("els_tq_trackIso"));
		els_tq_trackIso_branch->SetAddress(&els_tq_trackIso);
	}
	if(els_tq_trackIso_branch == 0 ) {
	cout << "Branch els_tq_trackIso does not exist." << endl;
	}
	jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_bCorrF") != 0) {
		jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_bCorrF"));
		jets_tq_bCorrF_branch->SetAddress(&jets_tq_bCorrF);
	}
	if(jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch jets_tq_bCorrF does not exist." << endl;
	}
	jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_cCorrF") != 0) {
		jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_cCorrF"));
		jets_tq_cCorrF_branch->SetAddress(&jets_tq_cCorrF);
	}
	if(jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch jets_tq_cCorrF does not exist." << endl;
	}
	jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_gluCorrF") != 0) {
		jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_gluCorrF"));
		jets_tq_gluCorrF_branch->SetAddress(&jets_tq_gluCorrF);
	}
	if(jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch jets_tq_gluCorrF does not exist." << endl;
	}
	jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("jets_tq_jetCharge") != 0) {
		jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("jets_tq_jetCharge"));
		jets_tq_jetCharge_branch->SetAddress(&jets_tq_jetCharge);
	}
	if(jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch jets_tq_jetCharge does not exist." << endl;
	}
	jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_noCorrF") != 0) {
		jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_noCorrF"));
		jets_tq_noCorrF_branch->SetAddress(&jets_tq_noCorrF);
	}
	if(jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch jets_tq_noCorrF does not exist." << endl;
	}
	jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("jets_tq_udsCorrF") != 0) {
		jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("jets_tq_udsCorrF"));
		jets_tq_udsCorrF_branch->SetAddress(&jets_tq_udsCorrF);
	}
	if(jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch jets_tq_udsCorrF does not exist." << endl;
	}
	mus_tq_caloIso_branch = 0;
	if (tree->GetAlias("mus_tq_caloIso") != 0) {
		mus_tq_caloIso_branch = tree->GetBranch(tree->GetAlias("mus_tq_caloIso"));
		mus_tq_caloIso_branch->SetAddress(&mus_tq_caloIso);
	}
	if(mus_tq_caloIso_branch == 0 ) {
	cout << "Branch mus_tq_caloIso does not exist." << endl;
	}
	mus_tq_leptonID_branch = 0;
	if (tree->GetAlias("mus_tq_leptonID") != 0) {
		mus_tq_leptonID_branch = tree->GetBranch(tree->GetAlias("mus_tq_leptonID"));
		mus_tq_leptonID_branch->SetAddress(&mus_tq_leptonID);
	}
	if(mus_tq_leptonID_branch == 0 ) {
	cout << "Branch mus_tq_leptonID does not exist." << endl;
	}
	mus_tq_lrComb_branch = 0;
	if (tree->GetAlias("mus_tq_lrComb") != 0) {
		mus_tq_lrComb_branch = tree->GetBranch(tree->GetAlias("mus_tq_lrComb"));
		mus_tq_lrComb_branch->SetAddress(&mus_tq_lrComb);
	}
	if(mus_tq_lrComb_branch == 0 ) {
	cout << "Branch mus_tq_lrComb does not exist." << endl;
	}
	mus_tq_trackIso_branch = 0;
	if (tree->GetAlias("mus_tq_trackIso") != 0) {
		mus_tq_trackIso_branch = tree->GetBranch(tree->GetAlias("mus_tq_trackIso"));
		mus_tq_trackIso_branch->SetAddress(&mus_tq_trackIso);
	}
	if(mus_tq_trackIso_branch == 0 ) {
	cout << "Branch mus_tq_trackIso does not exist." << endl;
	}
	trks_chi2_branch = 0;
	if (tree->GetAlias("trks_chi2") != 0) {
		trks_chi2_branch = tree->GetBranch(tree->GetAlias("trks_chi2"));
		trks_chi2_branch->SetAddress(&trks_chi2);
	}
	if(trks_chi2_branch == 0 ) {
	cout << "Branch trks_chi2 does not exist." << endl;
	}
	trks_d0_branch = 0;
	if (tree->GetAlias("trks_d0") != 0) {
		trks_d0_branch = tree->GetBranch(tree->GetAlias("trks_d0"));
		trks_d0_branch->SetAddress(&trks_d0);
	}
	if(trks_d0_branch == 0 ) {
	cout << "Branch trks_d0 does not exist." << endl;
	}
	trks_d0Err_branch = 0;
	if (tree->GetAlias("trks_d0Err") != 0) {
		trks_d0Err_branch = tree->GetBranch(tree->GetAlias("trks_d0Err"));
		trks_d0Err_branch->SetAddress(&trks_d0Err);
	}
	if(trks_d0Err_branch == 0 ) {
	cout << "Branch trks_d0Err does not exist." << endl;
	}
	trks_etaErr_branch = 0;
	if (tree->GetAlias("trks_etaErr") != 0) {
		trks_etaErr_branch = tree->GetBranch(tree->GetAlias("trks_etaErr"));
		trks_etaErr_branch->SetAddress(&trks_etaErr);
	}
	if(trks_etaErr_branch == 0 ) {
	cout << "Branch trks_etaErr does not exist." << endl;
	}
	trks_ndof_branch = 0;
	if (tree->GetAlias("trks_ndof") != 0) {
		trks_ndof_branch = tree->GetBranch(tree->GetAlias("trks_ndof"));
		trks_ndof_branch->SetAddress(&trks_ndof);
	}
	if(trks_ndof_branch == 0 ) {
	cout << "Branch trks_ndof does not exist." << endl;
	}
	trks_outerEta_branch = 0;
	if (tree->GetAlias("trks_outerEta") != 0) {
		trks_outerEta_branch = tree->GetBranch(tree->GetAlias("trks_outerEta"));
		trks_outerEta_branch->SetAddress(&trks_outerEta);
	}
	if(trks_outerEta_branch == 0 ) {
	cout << "Branch trks_outerEta does not exist." << endl;
	}
	trks_outerPhi_branch = 0;
	if (tree->GetAlias("trks_outerPhi") != 0) {
		trks_outerPhi_branch = tree->GetBranch(tree->GetAlias("trks_outerPhi"));
		trks_outerPhi_branch->SetAddress(&trks_outerPhi);
	}
	if(trks_outerPhi_branch == 0 ) {
	cout << "Branch trks_outerPhi does not exist." << endl;
	}
	trks_phiErr_branch = 0;
	if (tree->GetAlias("trks_phiErr") != 0) {
		trks_phiErr_branch = tree->GetBranch(tree->GetAlias("trks_phiErr"));
		trks_phiErr_branch->SetAddress(&trks_phiErr);
	}
	if(trks_phiErr_branch == 0 ) {
	cout << "Branch trks_phiErr does not exist." << endl;
	}
	trks_ptErr_branch = 0;
	if (tree->GetAlias("trks_ptErr") != 0) {
		trks_ptErr_branch = tree->GetBranch(tree->GetAlias("trks_ptErr"));
		trks_ptErr_branch->SetAddress(&trks_ptErr);
	}
	if(trks_ptErr_branch == 0 ) {
	cout << "Branch trks_ptErr does not exist." << endl;
	}
	trks_vertexphi_branch = 0;
	if (tree->GetAlias("trks_vertexphi") != 0) {
		trks_vertexphi_branch = tree->GetBranch(tree->GetAlias("trks_vertexphi"));
		trks_vertexphi_branch->SetAddress(&trks_vertexphi);
	}
	if(trks_vertexphi_branch == 0 ) {
	cout << "Branch trks_vertexphi does not exist." << endl;
	}
	trks_z0_branch = 0;
	if (tree->GetAlias("trks_z0") != 0) {
		trks_z0_branch = tree->GetBranch(tree->GetAlias("trks_z0"));
		trks_z0_branch->SetAddress(&trks_z0);
	}
	if(trks_z0_branch == 0 ) {
	cout << "Branch trks_z0 does not exist." << endl;
	}
	trks_z0Err_branch = 0;
	if (tree->GetAlias("trks_z0Err") != 0) {
		trks_z0Err_branch = tree->GetBranch(tree->GetAlias("trks_z0Err"));
		trks_z0Err_branch->SetAddress(&trks_z0Err);
	}
	if(trks_z0Err_branch == 0 ) {
	cout << "Branch trks_z0Err does not exist." << endl;
	}
	hyp_jets_EMFcor_branch = 0;
	if (tree->GetAlias("hyp_jets_EMFcor") != 0) {
		hyp_jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("hyp_jets_EMFcor"));
		hyp_jets_EMFcor_branch->SetAddress(&hyp_jets_EMFcor);
	}
	if(hyp_jets_EMFcor_branch == 0 ) {
	cout << "Branch hyp_jets_EMFcor does not exist." << endl;
	}
	hyp_jets_chFrac_branch = 0;
	if (tree->GetAlias("hyp_jets_chFrac") != 0) {
		hyp_jets_chFrac_branch = tree->GetBranch(tree->GetAlias("hyp_jets_chFrac"));
		hyp_jets_chFrac_branch->SetAddress(&hyp_jets_chFrac);
	}
	if(hyp_jets_chFrac_branch == 0 ) {
	cout << "Branch hyp_jets_chFrac does not exist." << endl;
	}
	hyp_jets_cor_branch = 0;
	if (tree->GetAlias("hyp_jets_cor") != 0) {
		hyp_jets_cor_branch = tree->GetBranch(tree->GetAlias("hyp_jets_cor"));
		hyp_jets_cor_branch->SetAddress(&hyp_jets_cor);
	}
	if(hyp_jets_cor_branch == 0 ) {
	cout << "Branch hyp_jets_cor does not exist." << endl;
	}
	hyp_jets_emFrac_branch = 0;
	if (tree->GetAlias("hyp_jets_emFrac") != 0) {
		hyp_jets_emFrac_branch = tree->GetBranch(tree->GetAlias("hyp_jets_emFrac"));
		hyp_jets_emFrac_branch->SetAddress(&hyp_jets_emFrac);
	}
	if(hyp_jets_emFrac_branch == 0 ) {
	cout << "Branch hyp_jets_emFrac does not exist." << endl;
	}
	hyp_jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_emEnergy") != 0) {
		hyp_jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_emEnergy"));
		hyp_jets_mc_emEnergy_branch->SetAddress(&hyp_jets_mc_emEnergy);
	}
	if(hyp_jets_mc_emEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_emEnergy does not exist." << endl;
	}
	hyp_jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_hadEnergy") != 0) {
		hyp_jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_hadEnergy"));
		hyp_jets_mc_hadEnergy_branch->SetAddress(&hyp_jets_mc_hadEnergy);
	}
	if(hyp_jets_mc_hadEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_hadEnergy does not exist." << endl;
	}
	hyp_jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_invEnergy") != 0) {
		hyp_jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_invEnergy"));
		hyp_jets_mc_invEnergy_branch->SetAddress(&hyp_jets_mc_invEnergy);
	}
	if(hyp_jets_mc_invEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_invEnergy does not exist." << endl;
	}
	hyp_jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_otherEnergy") != 0) {
		hyp_jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_otherEnergy"));
		hyp_jets_mc_otherEnergy_branch->SetAddress(&hyp_jets_mc_otherEnergy);
	}
	if(hyp_jets_mc_otherEnergy_branch == 0 ) {
	cout << "Branch hyp_jets_mc_otherEnergy does not exist." << endl;
	}
	hyp_jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_bCorrF") != 0) {
		hyp_jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_bCorrF"));
		hyp_jets_tq_bCorrF_branch->SetAddress(&hyp_jets_tq_bCorrF);
	}
	if(hyp_jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_bCorrF does not exist." << endl;
	}
	hyp_jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_cCorrF") != 0) {
		hyp_jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_cCorrF"));
		hyp_jets_tq_cCorrF_branch->SetAddress(&hyp_jets_tq_cCorrF);
	}
	if(hyp_jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_cCorrF does not exist." << endl;
	}
	hyp_jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_gluCorrF") != 0) {
		hyp_jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_gluCorrF"));
		hyp_jets_tq_gluCorrF_branch->SetAddress(&hyp_jets_tq_gluCorrF);
	}
	if(hyp_jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_gluCorrF does not exist." << endl;
	}
	hyp_jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_jetCharge") != 0) {
		hyp_jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_jetCharge"));
		hyp_jets_tq_jetCharge_branch->SetAddress(&hyp_jets_tq_jetCharge);
	}
	if(hyp_jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch hyp_jets_tq_jetCharge does not exist." << endl;
	}
	hyp_jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_noCorrF") != 0) {
		hyp_jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_noCorrF"));
		hyp_jets_tq_noCorrF_branch->SetAddress(&hyp_jets_tq_noCorrF);
	}
	if(hyp_jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_noCorrF does not exist." << endl;
	}
	hyp_jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_udsCorrF") != 0) {
		hyp_jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_udsCorrF"));
		hyp_jets_tq_udsCorrF_branch->SetAddress(&hyp_jets_tq_udsCorrF);
	}
	if(hyp_jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_jets_tq_udsCorrF does not exist." << endl;
	}
	hyp_other_jets_EMFcor_branch = 0;
	if (tree->GetAlias("hyp_other_jets_EMFcor") != 0) {
		hyp_other_jets_EMFcor_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_EMFcor"));
		hyp_other_jets_EMFcor_branch->SetAddress(&hyp_other_jets_EMFcor);
	}
	if(hyp_other_jets_EMFcor_branch == 0 ) {
	cout << "Branch hyp_other_jets_EMFcor does not exist." << endl;
	}
	hyp_other_jets_chFrac_branch = 0;
	if (tree->GetAlias("hyp_other_jets_chFrac") != 0) {
		hyp_other_jets_chFrac_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_chFrac"));
		hyp_other_jets_chFrac_branch->SetAddress(&hyp_other_jets_chFrac);
	}
	if(hyp_other_jets_chFrac_branch == 0 ) {
	cout << "Branch hyp_other_jets_chFrac does not exist." << endl;
	}
	hyp_other_jets_cor_branch = 0;
	if (tree->GetAlias("hyp_other_jets_cor") != 0) {
		hyp_other_jets_cor_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_cor"));
		hyp_other_jets_cor_branch->SetAddress(&hyp_other_jets_cor);
	}
	if(hyp_other_jets_cor_branch == 0 ) {
	cout << "Branch hyp_other_jets_cor does not exist." << endl;
	}
	hyp_other_jets_emFrac_branch = 0;
	if (tree->GetAlias("hyp_other_jets_emFrac") != 0) {
		hyp_other_jets_emFrac_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_emFrac"));
		hyp_other_jets_emFrac_branch->SetAddress(&hyp_other_jets_emFrac);
	}
	if(hyp_other_jets_emFrac_branch == 0 ) {
	cout << "Branch hyp_other_jets_emFrac does not exist." << endl;
	}
	hyp_other_jets_mc_emEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_emEnergy") != 0) {
		hyp_other_jets_mc_emEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_emEnergy"));
		hyp_other_jets_mc_emEnergy_branch->SetAddress(&hyp_other_jets_mc_emEnergy);
	}
	if(hyp_other_jets_mc_emEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_emEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_hadEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_hadEnergy") != 0) {
		hyp_other_jets_mc_hadEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_hadEnergy"));
		hyp_other_jets_mc_hadEnergy_branch->SetAddress(&hyp_other_jets_mc_hadEnergy);
	}
	if(hyp_other_jets_mc_hadEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_hadEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_invEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_invEnergy") != 0) {
		hyp_other_jets_mc_invEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_invEnergy"));
		hyp_other_jets_mc_invEnergy_branch->SetAddress(&hyp_other_jets_mc_invEnergy);
	}
	if(hyp_other_jets_mc_invEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_invEnergy does not exist." << endl;
	}
	hyp_other_jets_mc_otherEnergy_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_otherEnergy") != 0) {
		hyp_other_jets_mc_otherEnergy_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_otherEnergy"));
		hyp_other_jets_mc_otherEnergy_branch->SetAddress(&hyp_other_jets_mc_otherEnergy);
	}
	if(hyp_other_jets_mc_otherEnergy_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_otherEnergy does not exist." << endl;
	}
	hyp_other_jets_tq_bCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_bCorrF") != 0) {
		hyp_other_jets_tq_bCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_bCorrF"));
		hyp_other_jets_tq_bCorrF_branch->SetAddress(&hyp_other_jets_tq_bCorrF);
	}
	if(hyp_other_jets_tq_bCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_bCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_cCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_cCorrF") != 0) {
		hyp_other_jets_tq_cCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_cCorrF"));
		hyp_other_jets_tq_cCorrF_branch->SetAddress(&hyp_other_jets_tq_cCorrF);
	}
	if(hyp_other_jets_tq_cCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_cCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_gluCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_gluCorrF") != 0) {
		hyp_other_jets_tq_gluCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_gluCorrF"));
		hyp_other_jets_tq_gluCorrF_branch->SetAddress(&hyp_other_jets_tq_gluCorrF);
	}
	if(hyp_other_jets_tq_gluCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_gluCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_jetCharge_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_jetCharge") != 0) {
		hyp_other_jets_tq_jetCharge_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_jetCharge"));
		hyp_other_jets_tq_jetCharge_branch->SetAddress(&hyp_other_jets_tq_jetCharge);
	}
	if(hyp_other_jets_tq_jetCharge_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_jetCharge does not exist." << endl;
	}
	hyp_other_jets_tq_noCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_noCorrF") != 0) {
		hyp_other_jets_tq_noCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_noCorrF"));
		hyp_other_jets_tq_noCorrF_branch->SetAddress(&hyp_other_jets_tq_noCorrF);
	}
	if(hyp_other_jets_tq_noCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_noCorrF does not exist." << endl;
	}
	hyp_other_jets_tq_udsCorrF_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_udsCorrF") != 0) {
		hyp_other_jets_tq_udsCorrF_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_udsCorrF"));
		hyp_other_jets_tq_udsCorrF_branch->SetAddress(&hyp_other_jets_tq_udsCorrF);
	}
	if(hyp_other_jets_tq_udsCorrF_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_udsCorrF does not exist." << endl;
	}
	evt_HLT1_branch = 0;
	if (tree->GetAlias("evt_HLT1") != 0) {
		evt_HLT1_branch = tree->GetBranch(tree->GetAlias("evt_HLT1"));
		evt_HLT1_branch->SetAddress(&evt_HLT1);
	}
	if(evt_HLT1_branch == 0 ) {
	cout << "Branch evt_HLT1 does not exist." << endl;
	}
	evt_HLT2_branch = 0;
	if (tree->GetAlias("evt_HLT2") != 0) {
		evt_HLT2_branch = tree->GetBranch(tree->GetAlias("evt_HLT2"));
		evt_HLT2_branch->SetAddress(&evt_HLT2);
	}
	if(evt_HLT2_branch == 0 ) {
	cout << "Branch evt_HLT2 does not exist." << endl;
	}
	evt_HLT3_branch = 0;
	if (tree->GetAlias("evt_HLT3") != 0) {
		evt_HLT3_branch = tree->GetBranch(tree->GetAlias("evt_HLT3"));
		evt_HLT3_branch->SetAddress(&evt_HLT3);
	}
	if(evt_HLT3_branch == 0 ) {
	cout << "Branch evt_HLT3 does not exist." << endl;
	}
	evt_HLT4_branch = 0;
	if (tree->GetAlias("evt_HLT4") != 0) {
		evt_HLT4_branch = tree->GetBranch(tree->GetAlias("evt_HLT4"));
		evt_HLT4_branch->SetAddress(&evt_HLT4);
	}
	if(evt_HLT4_branch == 0 ) {
	cout << "Branch evt_HLT4 does not exist." << endl;
	}
	evt_L1_1_branch = 0;
	if (tree->GetAlias("evt_L1_1") != 0) {
		evt_L1_1_branch = tree->GetBranch(tree->GetAlias("evt_L1_1"));
		evt_L1_1_branch->SetAddress(&evt_L1_1);
	}
	if(evt_L1_1_branch == 0 ) {
	cout << "Branch evt_L1_1 does not exist." << endl;
	}
	evt_L1_2_branch = 0;
	if (tree->GetAlias("evt_L1_2") != 0) {
		evt_L1_2_branch = tree->GetBranch(tree->GetAlias("evt_L1_2"));
		evt_L1_2_branch->SetAddress(&evt_L1_2);
	}
	if(evt_L1_2_branch == 0 ) {
	cout << "Branch evt_L1_2 does not exist." << endl;
	}
	evt_L1_3_branch = 0;
	if (tree->GetAlias("evt_L1_3") != 0) {
		evt_L1_3_branch = tree->GetBranch(tree->GetAlias("evt_L1_3"));
		evt_L1_3_branch->SetAddress(&evt_L1_3);
	}
	if(evt_L1_3_branch == 0 ) {
	cout << "Branch evt_L1_3 does not exist." << endl;
	}
	evt_L1_4_branch = 0;
	if (tree->GetAlias("evt_L1_4") != 0) {
		evt_L1_4_branch = tree->GetBranch(tree->GetAlias("evt_L1_4"));
		evt_L1_4_branch->SetAddress(&evt_L1_4);
	}
	if(evt_L1_4_branch == 0 ) {
	cout << "Branch evt_L1_4 does not exist." << endl;
	}
	evt_event_branch = 0;
	if (tree->GetAlias("evt_event") != 0) {
		evt_event_branch = tree->GetBranch(tree->GetAlias("evt_event"));
		evt_event_branch->SetAddress(&evt_event);
	}
	if(evt_event_branch == 0 ) {
	cout << "Branch evt_event does not exist." << endl;
	}
	evt_run_branch = 0;
	if (tree->GetAlias("evt_run") != 0) {
		evt_run_branch = tree->GetBranch(tree->GetAlias("evt_run"));
		evt_run_branch->SetAddress(&evt_run);
	}
	if(evt_run_branch == 0 ) {
	cout << "Branch evt_run does not exist." << endl;
	}
	evt_nl1emiso_branch = 0;
	if (tree->GetAlias("evt_nl1emiso") != 0) {
		evt_nl1emiso_branch = tree->GetBranch(tree->GetAlias("evt_nl1emiso"));
		evt_nl1emiso_branch->SetAddress(&evt_nl1emiso);
	}
	if(evt_nl1emiso_branch == 0 ) {
	cout << "Branch evt_nl1emiso does not exist." << endl;
	}
	evt_nl1emnoiso_branch = 0;
	if (tree->GetAlias("evt_nl1emnoiso") != 0) {
		evt_nl1emnoiso_branch = tree->GetBranch(tree->GetAlias("evt_nl1emnoiso"));
		evt_nl1emnoiso_branch->SetAddress(&evt_nl1emnoiso);
	}
	if(evt_nl1emnoiso_branch == 0 ) {
	cout << "Branch evt_nl1emnoiso does not exist." << endl;
	}
	evt_nl1jetsc_branch = 0;
	if (tree->GetAlias("evt_nl1jetsc") != 0) {
		evt_nl1jetsc_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetsc"));
		evt_nl1jetsc_branch->SetAddress(&evt_nl1jetsc);
	}
	if(evt_nl1jetsc_branch == 0 ) {
	cout << "Branch evt_nl1jetsc does not exist." << endl;
	}
	evt_nl1jetsf_branch = 0;
	if (tree->GetAlias("evt_nl1jetsf") != 0) {
		evt_nl1jetsf_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetsf"));
		evt_nl1jetsf_branch->SetAddress(&evt_nl1jetsf);
	}
	if(evt_nl1jetsf_branch == 0 ) {
	cout << "Branch evt_nl1jetsf does not exist." << endl;
	}
	evt_nl1jetst_branch = 0;
	if (tree->GetAlias("evt_nl1jetst") != 0) {
		evt_nl1jetst_branch = tree->GetBranch(tree->GetAlias("evt_nl1jetst"));
		evt_nl1jetst_branch->SetAddress(&evt_nl1jetst);
	}
	if(evt_nl1jetst_branch == 0 ) {
	cout << "Branch evt_nl1jetst does not exist." << endl;
	}
	evt_nl1mus_branch = 0;
	if (tree->GetAlias("evt_nl1mus") != 0) {
		evt_nl1mus_branch = tree->GetBranch(tree->GetAlias("evt_nl1mus"));
		evt_nl1mus_branch->SetAddress(&evt_nl1mus);
	}
	if(evt_nl1mus_branch == 0 ) {
	cout << "Branch evt_nl1mus does not exist." << endl;
	}
	els_closestMuon_branch = 0;
	if (tree->GetAlias("els_closestMuon") != 0) {
		els_closestMuon_branch = tree->GetBranch(tree->GetAlias("els_closestMuon"));
		els_closestMuon_branch->SetAddress(&els_closestMuon);
	}
	if(els_closestMuon_branch == 0 ) {
	cout << "Branch els_closestMuon does not exist." << endl;
	}
	els_trkidx_branch = 0;
	if (tree->GetAlias("els_trkidx") != 0) {
		els_trkidx_branch = tree->GetBranch(tree->GetAlias("els_trkidx"));
		els_trkidx_branch->SetAddress(&els_trkidx);
	}
	if(els_trkidx_branch == 0 ) {
	cout << "Branch els_trkidx does not exist." << endl;
	}
	els_charge_branch = 0;
	if (tree->GetAlias("els_charge") != 0) {
		els_charge_branch = tree->GetBranch(tree->GetAlias("els_charge"));
		els_charge_branch->SetAddress(&els_charge);
	}
	if(els_charge_branch == 0 ) {
	cout << "Branch els_charge does not exist." << endl;
	}
	els_class_branch = 0;
	if (tree->GetAlias("els_class") != 0) {
		els_class_branch = tree->GetBranch(tree->GetAlias("els_class"));
		els_class_branch->SetAddress(&els_class);
	}
	if(els_class_branch == 0 ) {
	cout << "Branch els_class does not exist." << endl;
	}
	els_looseId_branch = 0;
	if (tree->GetAlias("els_looseId") != 0) {
		els_looseId_branch = tree->GetBranch(tree->GetAlias("els_looseId"));
		els_looseId_branch->SetAddress(&els_looseId);
	}
	if(els_looseId_branch == 0 ) {
	cout << "Branch els_looseId does not exist." << endl;
	}
	els_lostHits_branch = 0;
	if (tree->GetAlias("els_lostHits") != 0) {
		els_lostHits_branch = tree->GetBranch(tree->GetAlias("els_lostHits"));
		els_lostHits_branch->SetAddress(&els_lostHits);
	}
	if(els_lostHits_branch == 0 ) {
	cout << "Branch els_lostHits does not exist." << endl;
	}
	els_nSeed_branch = 0;
	if (tree->GetAlias("els_nSeed") != 0) {
		els_nSeed_branch = tree->GetBranch(tree->GetAlias("els_nSeed"));
		els_nSeed_branch->SetAddress(&els_nSeed);
	}
	if(els_nSeed_branch == 0 ) {
	cout << "Branch els_nSeed does not exist." << endl;
	}
	els_pass3looseId_branch = 0;
	if (tree->GetAlias("els_pass3looseId") != 0) {
		els_pass3looseId_branch = tree->GetBranch(tree->GetAlias("els_pass3looseId"));
		els_pass3looseId_branch->SetAddress(&els_pass3looseId);
	}
	if(els_pass3looseId_branch == 0 ) {
	cout << "Branch els_pass3looseId does not exist." << endl;
	}
	els_pass3simpleId_branch = 0;
	if (tree->GetAlias("els_pass3simpleId") != 0) {
		els_pass3simpleId_branch = tree->GetBranch(tree->GetAlias("els_pass3simpleId"));
		els_pass3simpleId_branch->SetAddress(&els_pass3simpleId);
	}
	if(els_pass3simpleId_branch == 0 ) {
	cout << "Branch els_pass3simpleId does not exist." << endl;
	}
	els_pass3tightId_branch = 0;
	if (tree->GetAlias("els_pass3tightId") != 0) {
		els_pass3tightId_branch = tree->GetBranch(tree->GetAlias("els_pass3tightId"));
		els_pass3tightId_branch->SetAddress(&els_pass3tightId);
	}
	if(els_pass3tightId_branch == 0 ) {
	cout << "Branch els_pass3tightId does not exist." << endl;
	}
	els_robustId_branch = 0;
	if (tree->GetAlias("els_robustId") != 0) {
		els_robustId_branch = tree->GetBranch(tree->GetAlias("els_robustId"));
		els_robustId_branch->SetAddress(&els_robustId);
	}
	if(els_robustId_branch == 0 ) {
	cout << "Branch els_robustId does not exist." << endl;
	}
	els_simpleIdPlus_branch = 0;
	if (tree->GetAlias("els_simpleIdPlus") != 0) {
		els_simpleIdPlus_branch = tree->GetBranch(tree->GetAlias("els_simpleIdPlus"));
		els_simpleIdPlus_branch->SetAddress(&els_simpleIdPlus);
	}
	if(els_simpleIdPlus_branch == 0 ) {
	cout << "Branch els_simpleIdPlus does not exist." << endl;
	}
	els_tightId_branch = 0;
	if (tree->GetAlias("els_tightId") != 0) {
		els_tightId_branch = tree->GetBranch(tree->GetAlias("els_tightId"));
		els_tightId_branch->SetAddress(&els_tightId);
	}
	if(els_tightId_branch == 0 ) {
	cout << "Branch els_tightId does not exist." << endl;
	}
	els_validHits_branch = 0;
	if (tree->GetAlias("els_validHits") != 0) {
		els_validHits_branch = tree->GetBranch(tree->GetAlias("els_validHits"));
		els_validHits_branch->SetAddress(&els_validHits);
	}
	if(els_validHits_branch == 0 ) {
	cout << "Branch els_validHits does not exist." << endl;
	}
	genps_id_branch = 0;
	if (tree->GetAlias("genps_id") != 0) {
		genps_id_branch = tree->GetBranch(tree->GetAlias("genps_id"));
		genps_id_branch->SetAddress(&genps_id);
	}
	if(genps_id_branch == 0 ) {
	cout << "Branch genps_id does not exist." << endl;
	}
	genps_id_mother_branch = 0;
	if (tree->GetAlias("genps_id_mother") != 0) {
		genps_id_mother_branch = tree->GetBranch(tree->GetAlias("genps_id_mother"));
		genps_id_mother_branch->SetAddress(&genps_id_mother);
	}
	if(genps_id_mother_branch == 0 ) {
	cout << "Branch genps_id_mother does not exist." << endl;
	}
	genps_status_branch = 0;
	if (tree->GetAlias("genps_status") != 0) {
		genps_status_branch = tree->GetBranch(tree->GetAlias("genps_status"));
		genps_status_branch->SetAddress(&genps_status);
	}
	if(genps_status_branch == 0 ) {
	cout << "Branch genps_status does not exist." << endl;
	}
	hyp_ll_charge_branch = 0;
	if (tree->GetAlias("hyp_ll_charge") != 0) {
		hyp_ll_charge_branch = tree->GetBranch(tree->GetAlias("hyp_ll_charge"));
		hyp_ll_charge_branch->SetAddress(&hyp_ll_charge);
	}
	if(hyp_ll_charge_branch == 0 ) {
	cout << "Branch hyp_ll_charge does not exist." << endl;
	}
	hyp_ll_id_branch = 0;
	if (tree->GetAlias("hyp_ll_id") != 0) {
		hyp_ll_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_id"));
		hyp_ll_id_branch->SetAddress(&hyp_ll_id);
	}
	if(hyp_ll_id_branch == 0 ) {
	cout << "Branch hyp_ll_id does not exist." << endl;
	}
	hyp_ll_index_branch = 0;
	if (tree->GetAlias("hyp_ll_index") != 0) {
		hyp_ll_index_branch = tree->GetBranch(tree->GetAlias("hyp_ll_index"));
		hyp_ll_index_branch->SetAddress(&hyp_ll_index);
	}
	if(hyp_ll_index_branch == 0 ) {
	cout << "Branch hyp_ll_index does not exist." << endl;
	}
	hyp_ll_lostHits_branch = 0;
	if (tree->GetAlias("hyp_ll_lostHits") != 0) {
		hyp_ll_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_lostHits"));
		hyp_ll_lostHits_branch->SetAddress(&hyp_ll_lostHits);
	}
	if(hyp_ll_lostHits_branch == 0 ) {
	cout << "Branch hyp_ll_lostHits does not exist." << endl;
	}
	hyp_ll_mc_id_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_id") != 0) {
		hyp_ll_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_id"));
		hyp_ll_mc_id_branch->SetAddress(&hyp_ll_mc_id);
	}
	if(hyp_ll_mc_id_branch == 0 ) {
	cout << "Branch hyp_ll_mc_id does not exist." << endl;
	}
	hyp_ll_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_ll_mc_motherid") != 0) {
		hyp_ll_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_ll_mc_motherid"));
		hyp_ll_mc_motherid_branch->SetAddress(&hyp_ll_mc_motherid);
	}
	if(hyp_ll_mc_motherid_branch == 0 ) {
	cout << "Branch hyp_ll_mc_motherid does not exist." << endl;
	}
	hyp_ll_validHits_branch = 0;
	if (tree->GetAlias("hyp_ll_validHits") != 0) {
		hyp_ll_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_ll_validHits"));
		hyp_ll_validHits_branch->SetAddress(&hyp_ll_validHits);
	}
	if(hyp_ll_validHits_branch == 0 ) {
	cout << "Branch hyp_ll_validHits does not exist." << endl;
	}
	hyp_lt_charge_branch = 0;
	if (tree->GetAlias("hyp_lt_charge") != 0) {
		hyp_lt_charge_branch = tree->GetBranch(tree->GetAlias("hyp_lt_charge"));
		hyp_lt_charge_branch->SetAddress(&hyp_lt_charge);
	}
	if(hyp_lt_charge_branch == 0 ) {
	cout << "Branch hyp_lt_charge does not exist." << endl;
	}
	hyp_lt_id_branch = 0;
	if (tree->GetAlias("hyp_lt_id") != 0) {
		hyp_lt_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_id"));
		hyp_lt_id_branch->SetAddress(&hyp_lt_id);
	}
	if(hyp_lt_id_branch == 0 ) {
	cout << "Branch hyp_lt_id does not exist." << endl;
	}
	hyp_lt_index_branch = 0;
	if (tree->GetAlias("hyp_lt_index") != 0) {
		hyp_lt_index_branch = tree->GetBranch(tree->GetAlias("hyp_lt_index"));
		hyp_lt_index_branch->SetAddress(&hyp_lt_index);
	}
	if(hyp_lt_index_branch == 0 ) {
	cout << "Branch hyp_lt_index does not exist." << endl;
	}
	hyp_lt_lostHits_branch = 0;
	if (tree->GetAlias("hyp_lt_lostHits") != 0) {
		hyp_lt_lostHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_lostHits"));
		hyp_lt_lostHits_branch->SetAddress(&hyp_lt_lostHits);
	}
	if(hyp_lt_lostHits_branch == 0 ) {
	cout << "Branch hyp_lt_lostHits does not exist." << endl;
	}
	hyp_lt_mc_id_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_id") != 0) {
		hyp_lt_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_id"));
		hyp_lt_mc_id_branch->SetAddress(&hyp_lt_mc_id);
	}
	if(hyp_lt_mc_id_branch == 0 ) {
	cout << "Branch hyp_lt_mc_id does not exist." << endl;
	}
	hyp_lt_mc_motherid_branch = 0;
	if (tree->GetAlias("hyp_lt_mc_motherid") != 0) {
		hyp_lt_mc_motherid_branch = tree->GetBranch(tree->GetAlias("hyp_lt_mc_motherid"));
		hyp_lt_mc_motherid_branch->SetAddress(&hyp_lt_mc_motherid);
	}
	if(hyp_lt_mc_motherid_branch == 0 ) {
	cout << "Branch hyp_lt_mc_motherid does not exist." << endl;
	}
	hyp_lt_validHits_branch = 0;
	if (tree->GetAlias("hyp_lt_validHits") != 0) {
		hyp_lt_validHits_branch = tree->GetBranch(tree->GetAlias("hyp_lt_validHits"));
		hyp_lt_validHits_branch->SetAddress(&hyp_lt_validHits);
	}
	if(hyp_lt_validHits_branch == 0 ) {
	cout << "Branch hyp_lt_validHits does not exist." << endl;
	}
	hyp_njets_branch = 0;
	if (tree->GetAlias("hyp_njets") != 0) {
		hyp_njets_branch = tree->GetBranch(tree->GetAlias("hyp_njets"));
		hyp_njets_branch->SetAddress(&hyp_njets);
	}
	if(hyp_njets_branch == 0 ) {
	cout << "Branch hyp_njets does not exist." << endl;
	}
	hyp_nojets_branch = 0;
	if (tree->GetAlias("hyp_nojets") != 0) {
		hyp_nojets_branch = tree->GetBranch(tree->GetAlias("hyp_nojets"));
		hyp_nojets_branch->SetAddress(&hyp_nojets);
	}
	if(hyp_nojets_branch == 0 ) {
	cout << "Branch hyp_nojets does not exist." << endl;
	}
	hyp_type_branch = 0;
	if (tree->GetAlias("hyp_type") != 0) {
		hyp_type_branch = tree->GetBranch(tree->GetAlias("hyp_type"));
		hyp_type_branch->SetAddress(&hyp_type);
	}
	if(hyp_type_branch == 0 ) {
	cout << "Branch hyp_type does not exist." << endl;
	}
	hyp_quadlep_first_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_type") != 0) {
		hyp_quadlep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_type"));
		hyp_quadlep_first_type_branch->SetAddress(&hyp_quadlep_first_type);
	}
	if(hyp_quadlep_first_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_first_type does not exist." << endl;
	}
	hyp_quadlep_fourth_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_type") != 0) {
		hyp_quadlep_fourth_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_type"));
		hyp_quadlep_fourth_type_branch->SetAddress(&hyp_quadlep_fourth_type);
	}
	if(hyp_quadlep_fourth_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_fourth_type does not exist." << endl;
	}
	hyp_quadlep_second_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_type") != 0) {
		hyp_quadlep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_type"));
		hyp_quadlep_second_type_branch->SetAddress(&hyp_quadlep_second_type);
	}
	if(hyp_quadlep_second_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_second_type does not exist." << endl;
	}
	hyp_quadlep_third_type_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_type") != 0) {
		hyp_quadlep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_type"));
		hyp_quadlep_third_type_branch->SetAddress(&hyp_quadlep_third_type);
	}
	if(hyp_quadlep_third_type_branch == 0 ) {
	cout << "Branch hyp_quadlep_third_type does not exist." << endl;
	}
	hyp_trilep_first_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_type") != 0) {
		hyp_trilep_first_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_type"));
		hyp_trilep_first_type_branch->SetAddress(&hyp_trilep_first_type);
	}
	if(hyp_trilep_first_type_branch == 0 ) {
	cout << "Branch hyp_trilep_first_type does not exist." << endl;
	}
	hyp_trilep_second_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_type") != 0) {
		hyp_trilep_second_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_type"));
		hyp_trilep_second_type_branch->SetAddress(&hyp_trilep_second_type);
	}
	if(hyp_trilep_second_type_branch == 0 ) {
	cout << "Branch hyp_trilep_second_type does not exist." << endl;
	}
	hyp_trilep_third_type_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_type") != 0) {
		hyp_trilep_third_type_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_type"));
		hyp_trilep_third_type_branch->SetAddress(&hyp_trilep_third_type);
	}
	if(hyp_trilep_third_type_branch == 0 ) {
	cout << "Branch hyp_trilep_third_type does not exist." << endl;
	}
	jets_closestElectron_branch = 0;
	if (tree->GetAlias("jets_closestElectron") != 0) {
		jets_closestElectron_branch = tree->GetBranch(tree->GetAlias("jets_closestElectron"));
		jets_closestElectron_branch->SetAddress(&jets_closestElectron);
	}
	if(jets_closestElectron_branch == 0 ) {
	cout << "Branch jets_closestElectron does not exist." << endl;
	}
	jets_closestMuon_branch = 0;
	if (tree->GetAlias("jets_closestMuon") != 0) {
		jets_closestMuon_branch = tree->GetBranch(tree->GetAlias("jets_closestMuon"));
		jets_closestMuon_branch->SetAddress(&jets_closestMuon);
	}
	if(jets_closestMuon_branch == 0 ) {
	cout << "Branch jets_closestMuon does not exist." << endl;
	}
	l1emiso_ieta_branch = 0;
	if (tree->GetAlias("l1emiso_ieta") != 0) {
		l1emiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1emiso_ieta"));
		l1emiso_ieta_branch->SetAddress(&l1emiso_ieta);
	}
	if(l1emiso_ieta_branch == 0 ) {
	cout << "Branch l1emiso_ieta does not exist." << endl;
	}
	l1emiso_iphi_branch = 0;
	if (tree->GetAlias("l1emiso_iphi") != 0) {
		l1emiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1emiso_iphi"));
		l1emiso_iphi_branch->SetAddress(&l1emiso_iphi);
	}
	if(l1emiso_iphi_branch == 0 ) {
	cout << "Branch l1emiso_iphi does not exist." << endl;
	}
	l1emiso_rawId_branch = 0;
	if (tree->GetAlias("l1emiso_rawId") != 0) {
		l1emiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1emiso_rawId"));
		l1emiso_rawId_branch->SetAddress(&l1emiso_rawId);
	}
	if(l1emiso_rawId_branch == 0 ) {
	cout << "Branch l1emiso_rawId does not exist." << endl;
	}
	l1emiso_type_branch = 0;
	if (tree->GetAlias("l1emiso_type") != 0) {
		l1emiso_type_branch = tree->GetBranch(tree->GetAlias("l1emiso_type"));
		l1emiso_type_branch->SetAddress(&l1emiso_type);
	}
	if(l1emiso_type_branch == 0 ) {
	cout << "Branch l1emiso_type does not exist." << endl;
	}
	l1emnoiso_ieta_branch = 0;
	if (tree->GetAlias("l1emnoiso_ieta") != 0) {
		l1emnoiso_ieta_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_ieta"));
		l1emnoiso_ieta_branch->SetAddress(&l1emnoiso_ieta);
	}
	if(l1emnoiso_ieta_branch == 0 ) {
	cout << "Branch l1emnoiso_ieta does not exist." << endl;
	}
	l1emnoiso_iphi_branch = 0;
	if (tree->GetAlias("l1emnoiso_iphi") != 0) {
		l1emnoiso_iphi_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_iphi"));
		l1emnoiso_iphi_branch->SetAddress(&l1emnoiso_iphi);
	}
	if(l1emnoiso_iphi_branch == 0 ) {
	cout << "Branch l1emnoiso_iphi does not exist." << endl;
	}
	l1emnoiso_rawId_branch = 0;
	if (tree->GetAlias("l1emnoiso_rawId") != 0) {
		l1emnoiso_rawId_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_rawId"));
		l1emnoiso_rawId_branch->SetAddress(&l1emnoiso_rawId);
	}
	if(l1emnoiso_rawId_branch == 0 ) {
	cout << "Branch l1emnoiso_rawId does not exist." << endl;
	}
	l1emnoiso_type_branch = 0;
	if (tree->GetAlias("l1emnoiso_type") != 0) {
		l1emnoiso_type_branch = tree->GetBranch(tree->GetAlias("l1emnoiso_type"));
		l1emnoiso_type_branch->SetAddress(&l1emnoiso_type);
	}
	if(l1emnoiso_type_branch == 0 ) {
	cout << "Branch l1emnoiso_type does not exist." << endl;
	}
	l1jetsc_ieta_branch = 0;
	if (tree->GetAlias("l1jetsc_ieta") != 0) {
		l1jetsc_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetsc_ieta"));
		l1jetsc_ieta_branch->SetAddress(&l1jetsc_ieta);
	}
	if(l1jetsc_ieta_branch == 0 ) {
	cout << "Branch l1jetsc_ieta does not exist." << endl;
	}
	l1jetsc_iphi_branch = 0;
	if (tree->GetAlias("l1jetsc_iphi") != 0) {
		l1jetsc_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetsc_iphi"));
		l1jetsc_iphi_branch->SetAddress(&l1jetsc_iphi);
	}
	if(l1jetsc_iphi_branch == 0 ) {
	cout << "Branch l1jetsc_iphi does not exist." << endl;
	}
	l1jetsc_rawId_branch = 0;
	if (tree->GetAlias("l1jetsc_rawId") != 0) {
		l1jetsc_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetsc_rawId"));
		l1jetsc_rawId_branch->SetAddress(&l1jetsc_rawId);
	}
	if(l1jetsc_rawId_branch == 0 ) {
	cout << "Branch l1jetsc_rawId does not exist." << endl;
	}
	l1jetsc_type_branch = 0;
	if (tree->GetAlias("l1jetsc_type") != 0) {
		l1jetsc_type_branch = tree->GetBranch(tree->GetAlias("l1jetsc_type"));
		l1jetsc_type_branch->SetAddress(&l1jetsc_type);
	}
	if(l1jetsc_type_branch == 0 ) {
	cout << "Branch l1jetsc_type does not exist." << endl;
	}
	l1jetsf_ieta_branch = 0;
	if (tree->GetAlias("l1jetsf_ieta") != 0) {
		l1jetsf_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetsf_ieta"));
		l1jetsf_ieta_branch->SetAddress(&l1jetsf_ieta);
	}
	if(l1jetsf_ieta_branch == 0 ) {
	cout << "Branch l1jetsf_ieta does not exist." << endl;
	}
	l1jetsf_iphi_branch = 0;
	if (tree->GetAlias("l1jetsf_iphi") != 0) {
		l1jetsf_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetsf_iphi"));
		l1jetsf_iphi_branch->SetAddress(&l1jetsf_iphi);
	}
	if(l1jetsf_iphi_branch == 0 ) {
	cout << "Branch l1jetsf_iphi does not exist." << endl;
	}
	l1jetsf_rawId_branch = 0;
	if (tree->GetAlias("l1jetsf_rawId") != 0) {
		l1jetsf_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetsf_rawId"));
		l1jetsf_rawId_branch->SetAddress(&l1jetsf_rawId);
	}
	if(l1jetsf_rawId_branch == 0 ) {
	cout << "Branch l1jetsf_rawId does not exist." << endl;
	}
	l1jetsf_type_branch = 0;
	if (tree->GetAlias("l1jetsf_type") != 0) {
		l1jetsf_type_branch = tree->GetBranch(tree->GetAlias("l1jetsf_type"));
		l1jetsf_type_branch->SetAddress(&l1jetsf_type);
	}
	if(l1jetsf_type_branch == 0 ) {
	cout << "Branch l1jetsf_type does not exist." << endl;
	}
	l1jetst_ieta_branch = 0;
	if (tree->GetAlias("l1jetst_ieta") != 0) {
		l1jetst_ieta_branch = tree->GetBranch(tree->GetAlias("l1jetst_ieta"));
		l1jetst_ieta_branch->SetAddress(&l1jetst_ieta);
	}
	if(l1jetst_ieta_branch == 0 ) {
	cout << "Branch l1jetst_ieta does not exist." << endl;
	}
	l1jetst_iphi_branch = 0;
	if (tree->GetAlias("l1jetst_iphi") != 0) {
		l1jetst_iphi_branch = tree->GetBranch(tree->GetAlias("l1jetst_iphi"));
		l1jetst_iphi_branch->SetAddress(&l1jetst_iphi);
	}
	if(l1jetst_iphi_branch == 0 ) {
	cout << "Branch l1jetst_iphi does not exist." << endl;
	}
	l1jetst_rawId_branch = 0;
	if (tree->GetAlias("l1jetst_rawId") != 0) {
		l1jetst_rawId_branch = tree->GetBranch(tree->GetAlias("l1jetst_rawId"));
		l1jetst_rawId_branch->SetAddress(&l1jetst_rawId);
	}
	if(l1jetst_rawId_branch == 0 ) {
	cout << "Branch l1jetst_rawId does not exist." << endl;
	}
	l1jetst_type_branch = 0;
	if (tree->GetAlias("l1jetst_type") != 0) {
		l1jetst_type_branch = tree->GetBranch(tree->GetAlias("l1jetst_type"));
		l1jetst_type_branch->SetAddress(&l1jetst_type);
	}
	if(l1jetst_type_branch == 0 ) {
	cout << "Branch l1jetst_type does not exist." << endl;
	}
	l1mus_flags_branch = 0;
	if (tree->GetAlias("l1mus_flags") != 0) {
		l1mus_flags_branch = tree->GetBranch(tree->GetAlias("l1mus_flags"));
		l1mus_flags_branch->SetAddress(&l1mus_flags);
	}
	if(l1mus_flags_branch == 0 ) {
	cout << "Branch l1mus_flags does not exist." << endl;
	}
	l1mus_q_branch = 0;
	if (tree->GetAlias("l1mus_q") != 0) {
		l1mus_q_branch = tree->GetBranch(tree->GetAlias("l1mus_q"));
		l1mus_q_branch->SetAddress(&l1mus_q);
	}
	if(l1mus_q_branch == 0 ) {
	cout << "Branch l1mus_q does not exist." << endl;
	}
	l1mus_qual_branch = 0;
	if (tree->GetAlias("l1mus_qual") != 0) {
		l1mus_qual_branch = tree->GetBranch(tree->GetAlias("l1mus_qual"));
		l1mus_qual_branch->SetAddress(&l1mus_qual);
	}
	if(l1mus_qual_branch == 0 ) {
	cout << "Branch l1mus_qual does not exist." << endl;
	}
	l1mus_qualFlags_branch = 0;
	if (tree->GetAlias("l1mus_qualFlags") != 0) {
		l1mus_qualFlags_branch = tree->GetBranch(tree->GetAlias("l1mus_qualFlags"));
		l1mus_qualFlags_branch->SetAddress(&l1mus_qualFlags);
	}
	if(l1mus_qualFlags_branch == 0 ) {
	cout << "Branch l1mus_qualFlags does not exist." << endl;
	}
	mus_closestEle_branch = 0;
	if (tree->GetAlias("mus_closestEle") != 0) {
		mus_closestEle_branch = tree->GetBranch(tree->GetAlias("mus_closestEle"));
		mus_closestEle_branch->SetAddress(&mus_closestEle);
	}
	if(mus_closestEle_branch == 0 ) {
	cout << "Branch mus_closestEle does not exist." << endl;
	}
	mus_closestJet_branch = 0;
	if (tree->GetAlias("mus_closestJet") != 0) {
		mus_closestJet_branch = tree->GetBranch(tree->GetAlias("mus_closestJet"));
		mus_closestJet_branch->SetAddress(&mus_closestJet);
	}
	if(mus_closestJet_branch == 0 ) {
	cout << "Branch mus_closestJet does not exist." << endl;
	}
	mus_trkidx_branch = 0;
	if (tree->GetAlias("mus_trkidx") != 0) {
		mus_trkidx_branch = tree->GetBranch(tree->GetAlias("mus_trkidx"));
		mus_trkidx_branch->SetAddress(&mus_trkidx);
	}
	if(mus_trkidx_branch == 0 ) {
	cout << "Branch mus_trkidx does not exist." << endl;
	}
	mus_charge_branch = 0;
	if (tree->GetAlias("mus_charge") != 0) {
		mus_charge_branch = tree->GetBranch(tree->GetAlias("mus_charge"));
		mus_charge_branch->SetAddress(&mus_charge);
	}
	if(mus_charge_branch == 0 ) {
	cout << "Branch mus_charge does not exist." << endl;
	}
	mus_gfit_validHits_branch = 0;
	if (tree->GetAlias("mus_gfit_validHits") != 0) {
		mus_gfit_validHits_branch = tree->GetBranch(tree->GetAlias("mus_gfit_validHits"));
		mus_gfit_validHits_branch->SetAddress(&mus_gfit_validHits);
	}
	if(mus_gfit_validHits_branch == 0 ) {
	cout << "Branch mus_gfit_validHits does not exist." << endl;
	}
	mus_iso03_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso03_ntrk") != 0) {
		mus_iso03_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso03_ntrk"));
		mus_iso03_ntrk_branch->SetAddress(&mus_iso03_ntrk);
	}
	if(mus_iso03_ntrk_branch == 0 ) {
	cout << "Branch mus_iso03_ntrk does not exist." << endl;
	}
	mus_iso05_ntrk_branch = 0;
	if (tree->GetAlias("mus_iso05_ntrk") != 0) {
		mus_iso05_ntrk_branch = tree->GetBranch(tree->GetAlias("mus_iso05_ntrk"));
		mus_iso05_ntrk_branch->SetAddress(&mus_iso05_ntrk);
	}
	if(mus_iso05_ntrk_branch == 0 ) {
	cout << "Branch mus_iso05_ntrk does not exist." << endl;
	}
	mus_lostHits_branch = 0;
	if (tree->GetAlias("mus_lostHits") != 0) {
		mus_lostHits_branch = tree->GetBranch(tree->GetAlias("mus_lostHits"));
		mus_lostHits_branch->SetAddress(&mus_lostHits);
	}
	if(mus_lostHits_branch == 0 ) {
	cout << "Branch mus_lostHits does not exist." << endl;
	}
	mus_nmatches_branch = 0;
	if (tree->GetAlias("mus_nmatches") != 0) {
		mus_nmatches_branch = tree->GetBranch(tree->GetAlias("mus_nmatches"));
		mus_nmatches_branch->SetAddress(&mus_nmatches);
	}
	if(mus_nmatches_branch == 0 ) {
	cout << "Branch mus_nmatches does not exist." << endl;
	}
	mus_pid_TM2DCompatibilityLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityLoose") != 0) {
		mus_pid_TM2DCompatibilityLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityLoose"));
		mus_pid_TM2DCompatibilityLoose_branch->SetAddress(&mus_pid_TM2DCompatibilityLoose);
	}
	if(mus_pid_TM2DCompatibilityLoose_branch == 0 ) {
	cout << "Branch mus_pid_TM2DCompatibilityLoose does not exist." << endl;
	}
	mus_pid_TM2DCompatibilityTight_branch = 0;
	if (tree->GetAlias("mus_pid_TM2DCompatibilityTight") != 0) {
		mus_pid_TM2DCompatibilityTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TM2DCompatibilityTight"));
		mus_pid_TM2DCompatibilityTight_branch->SetAddress(&mus_pid_TM2DCompatibilityTight);
	}
	if(mus_pid_TM2DCompatibilityTight_branch == 0 ) {
	cout << "Branch mus_pid_TM2DCompatibilityTight does not exist." << endl;
	}
	mus_pid_TMLastStationLoose_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationLoose") != 0) {
		mus_pid_TMLastStationLoose_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationLoose"));
		mus_pid_TMLastStationLoose_branch->SetAddress(&mus_pid_TMLastStationLoose);
	}
	if(mus_pid_TMLastStationLoose_branch == 0 ) {
	cout << "Branch mus_pid_TMLastStationLoose does not exist." << endl;
	}
	mus_pid_TMLastStationTight_branch = 0;
	if (tree->GetAlias("mus_pid_TMLastStationTight") != 0) {
		mus_pid_TMLastStationTight_branch = tree->GetBranch(tree->GetAlias("mus_pid_TMLastStationTight"));
		mus_pid_TMLastStationTight_branch->SetAddress(&mus_pid_TMLastStationTight);
	}
	if(mus_pid_TMLastStationTight_branch == 0 ) {
	cout << "Branch mus_pid_TMLastStationTight does not exist." << endl;
	}
	mus_trkrefkey_branch = 0;
	if (tree->GetAlias("mus_trkrefkey") != 0) {
		mus_trkrefkey_branch = tree->GetBranch(tree->GetAlias("mus_trkrefkey"));
		mus_trkrefkey_branch->SetAddress(&mus_trkrefkey);
	}
	if(mus_trkrefkey_branch == 0 ) {
	cout << "Branch mus_trkrefkey does not exist." << endl;
	}
	mus_validHits_branch = 0;
	if (tree->GetAlias("mus_validHits") != 0) {
		mus_validHits_branch = tree->GetBranch(tree->GetAlias("mus_validHits"));
		mus_validHits_branch->SetAddress(&mus_validHits);
	}
	if(mus_validHits_branch == 0 ) {
	cout << "Branch mus_validHits does not exist." << endl;
	}
	els_tq_egammaTkNumIso_branch = 0;
	if (tree->GetAlias("els_tq_egammaTkNumIso") != 0) {
		els_tq_egammaTkNumIso_branch = tree->GetBranch(tree->GetAlias("els_tq_egammaTkNumIso"));
		els_tq_egammaTkNumIso_branch->SetAddress(&els_tq_egammaTkNumIso);
	}
	if(els_tq_egammaTkNumIso_branch == 0 ) {
	cout << "Branch els_tq_egammaTkNumIso does not exist." << endl;
	}
	els_tq_genID_branch = 0;
	if (tree->GetAlias("els_tq_genID") != 0) {
		els_tq_genID_branch = tree->GetBranch(tree->GetAlias("els_tq_genID"));
		els_tq_genID_branch->SetAddress(&els_tq_genID);
	}
	if(els_tq_genID_branch == 0 ) {
	cout << "Branch els_tq_genID does not exist." << endl;
	}
	els_tq_genMotherID_branch = 0;
	if (tree->GetAlias("els_tq_genMotherID") != 0) {
		els_tq_genMotherID_branch = tree->GetBranch(tree->GetAlias("els_tq_genMotherID"));
		els_tq_genMotherID_branch->SetAddress(&els_tq_genMotherID);
	}
	if(els_tq_genMotherID_branch == 0 ) {
	cout << "Branch els_tq_genMotherID does not exist." << endl;
	}
	jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("jets_tq_genPartonMother_id") != 0) {
		jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("jets_tq_genPartonMother_id"));
		jets_tq_genPartonMother_id_branch->SetAddress(&jets_tq_genPartonMother_id);
	}
	if(jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch jets_tq_genPartonMother_id does not exist." << endl;
	}
	jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("jets_tq_genParton_id") != 0) {
		jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("jets_tq_genParton_id"));
		jets_tq_genParton_id_branch->SetAddress(&jets_tq_genParton_id);
	}
	if(jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch jets_tq_genParton_id does not exist." << endl;
	}
	jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("jets_tq_partonFlavour") != 0) {
		jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("jets_tq_partonFlavour"));
		jets_tq_partonFlavour_branch->SetAddress(&jets_tq_partonFlavour);
	}
	if(jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch jets_tq_partonFlavour does not exist." << endl;
	}
	mus_tq_genID_branch = 0;
	if (tree->GetAlias("mus_tq_genID") != 0) {
		mus_tq_genID_branch = tree->GetBranch(tree->GetAlias("mus_tq_genID"));
		mus_tq_genID_branch->SetAddress(&mus_tq_genID);
	}
	if(mus_tq_genID_branch == 0 ) {
	cout << "Branch mus_tq_genID does not exist." << endl;
	}
	mus_tq_genMotherID_branch = 0;
	if (tree->GetAlias("mus_tq_genMotherID") != 0) {
		mus_tq_genMotherID_branch = tree->GetBranch(tree->GetAlias("mus_tq_genMotherID"));
		mus_tq_genMotherID_branch->SetAddress(&mus_tq_genMotherID);
	}
	if(mus_tq_genMotherID_branch == 0 ) {
	cout << "Branch mus_tq_genMotherID does not exist." << endl;
	}
	trks_charge_branch = 0;
	if (tree->GetAlias("trks_charge") != 0) {
		trks_charge_branch = tree->GetBranch(tree->GetAlias("trks_charge"));
		trks_charge_branch->SetAddress(&trks_charge);
	}
	if(trks_charge_branch == 0 ) {
	cout << "Branch trks_charge does not exist." << endl;
	}
	trks_lostHits_branch = 0;
	if (tree->GetAlias("trks_lostHits") != 0) {
		trks_lostHits_branch = tree->GetBranch(tree->GetAlias("trks_lostHits"));
		trks_lostHits_branch->SetAddress(&trks_lostHits);
	}
	if(trks_lostHits_branch == 0 ) {
	cout << "Branch trks_lostHits does not exist." << endl;
	}
	trks_validHits_branch = 0;
	if (tree->GetAlias("trks_validHits") != 0) {
		trks_validHits_branch = tree->GetBranch(tree->GetAlias("trks_validHits"));
		trks_validHits_branch->SetAddress(&trks_validHits);
	}
	if(trks_validHits_branch == 0 ) {
	cout << "Branch trks_validHits does not exist." << endl;
	}
	hyp_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_jets_mc_id") != 0) {
		hyp_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_mc_id"));
		hyp_jets_mc_id_branch->SetAddress(&hyp_jets_mc_id);
	}
	if(hyp_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_jets_mc_id does not exist." << endl;
	}
	hyp_jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genPartonMother_id") != 0) {
		hyp_jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genPartonMother_id"));
		hyp_jets_tq_genPartonMother_id_branch->SetAddress(&hyp_jets_tq_genPartonMother_id);
	}
	if(hyp_jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genPartonMother_id does not exist." << endl;
	}
	hyp_jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_genParton_id") != 0) {
		hyp_jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_genParton_id"));
		hyp_jets_tq_genParton_id_branch->SetAddress(&hyp_jets_tq_genParton_id);
	}
	if(hyp_jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch hyp_jets_tq_genParton_id does not exist." << endl;
	}
	hyp_jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_jets_tq_partonFlavour") != 0) {
		hyp_jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_jets_tq_partonFlavour"));
		hyp_jets_tq_partonFlavour_branch->SetAddress(&hyp_jets_tq_partonFlavour);
	}
	if(hyp_jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_jets_tq_partonFlavour does not exist." << endl;
	}
	hyp_other_jets_mc_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_mc_id") != 0) {
		hyp_other_jets_mc_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_mc_id"));
		hyp_other_jets_mc_id_branch->SetAddress(&hyp_other_jets_mc_id);
	}
	if(hyp_other_jets_mc_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_mc_id does not exist." << endl;
	}
	hyp_other_jets_tq_genPartonMother_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genPartonMother_id") != 0) {
		hyp_other_jets_tq_genPartonMother_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genPartonMother_id"));
		hyp_other_jets_tq_genPartonMother_id_branch->SetAddress(&hyp_other_jets_tq_genPartonMother_id);
	}
	if(hyp_other_jets_tq_genPartonMother_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genPartonMother_id does not exist." << endl;
	}
	hyp_other_jets_tq_genParton_id_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_genParton_id") != 0) {
		hyp_other_jets_tq_genParton_id_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_genParton_id"));
		hyp_other_jets_tq_genParton_id_branch->SetAddress(&hyp_other_jets_tq_genParton_id);
	}
	if(hyp_other_jets_tq_genParton_id_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_genParton_id does not exist." << endl;
	}
	hyp_other_jets_tq_partonFlavour_branch = 0;
	if (tree->GetAlias("hyp_other_jets_tq_partonFlavour") != 0) {
		hyp_other_jets_tq_partonFlavour_branch = tree->GetBranch(tree->GetAlias("hyp_other_jets_tq_partonFlavour"));
		hyp_other_jets_tq_partonFlavour_branch->SetAddress(&hyp_other_jets_tq_partonFlavour);
	}
	if(hyp_other_jets_tq_partonFlavour_branch == 0 ) {
	cout << "Branch hyp_other_jets_tq_partonFlavour does not exist." << endl;
	}
	hyp_quadlep_jets_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_jets_index") != 0) {
		hyp_quadlep_jets_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_jets_index"));
		hyp_quadlep_jets_index_branch->SetAddress(&hyp_quadlep_jets_index);
	}
	if(hyp_quadlep_jets_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_jets_index does not exist." << endl;
	}
	hyp_trilep_jets_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_jets_index") != 0) {
		hyp_trilep_jets_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_jets_index"));
		hyp_trilep_jets_index_branch->SetAddress(&hyp_trilep_jets_index);
	}
	if(hyp_trilep_jets_index_branch == 0 ) {
	cout << "Branch hyp_trilep_jets_index does not exist." << endl;
	}
	evt_nels_branch = 0;
	if (tree->GetAlias("evt_nels") != 0) {
		evt_nels_branch = tree->GetBranch(tree->GetAlias("evt_nels"));
		evt_nels_branch->SetAddress(&evt_nels);
	}
	if(evt_nels_branch == 0 ) {
	cout << "Branch evt_nels does not exist." << endl;
	}
	evt_njets_branch = 0;
	if (tree->GetAlias("evt_njets") != 0) {
		evt_njets_branch = tree->GetBranch(tree->GetAlias("evt_njets"));
		evt_njets_branch->SetAddress(&evt_njets);
	}
	if(evt_njets_branch == 0 ) {
	cout << "Branch evt_njets does not exist." << endl;
	}
	hyp_quadlep_bucket_branch = 0;
	if (tree->GetAlias("hyp_quadlep_bucket") != 0) {
		hyp_quadlep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_bucket"));
		hyp_quadlep_bucket_branch->SetAddress(&hyp_quadlep_bucket);
	}
	if(hyp_quadlep_bucket_branch == 0 ) {
	cout << "Branch hyp_quadlep_bucket does not exist." << endl;
	}
	hyp_quadlep_first_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_first_index") != 0) {
		hyp_quadlep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_first_index"));
		hyp_quadlep_first_index_branch->SetAddress(&hyp_quadlep_first_index);
	}
	if(hyp_quadlep_first_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_first_index does not exist." << endl;
	}
	hyp_quadlep_fourth_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_fourth_index") != 0) {
		hyp_quadlep_fourth_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_fourth_index"));
		hyp_quadlep_fourth_index_branch->SetAddress(&hyp_quadlep_fourth_index);
	}
	if(hyp_quadlep_fourth_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_fourth_index does not exist." << endl;
	}
	hyp_quadlep_second_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_second_index") != 0) {
		hyp_quadlep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_second_index"));
		hyp_quadlep_second_index_branch->SetAddress(&hyp_quadlep_second_index);
	}
	if(hyp_quadlep_second_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_second_index does not exist." << endl;
	}
	hyp_quadlep_third_index_branch = 0;
	if (tree->GetAlias("hyp_quadlep_third_index") != 0) {
		hyp_quadlep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_quadlep_third_index"));
		hyp_quadlep_third_index_branch->SetAddress(&hyp_quadlep_third_index);
	}
	if(hyp_quadlep_third_index_branch == 0 ) {
	cout << "Branch hyp_quadlep_third_index does not exist." << endl;
	}
	hyp_trilep_bucket_branch = 0;
	if (tree->GetAlias("hyp_trilep_bucket") != 0) {
		hyp_trilep_bucket_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_bucket"));
		hyp_trilep_bucket_branch->SetAddress(&hyp_trilep_bucket);
	}
	if(hyp_trilep_bucket_branch == 0 ) {
	cout << "Branch hyp_trilep_bucket does not exist." << endl;
	}
	hyp_trilep_first_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_first_index") != 0) {
		hyp_trilep_first_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_first_index"));
		hyp_trilep_first_index_branch->SetAddress(&hyp_trilep_first_index);
	}
	if(hyp_trilep_first_index_branch == 0 ) {
	cout << "Branch hyp_trilep_first_index does not exist." << endl;
	}
	hyp_trilep_second_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_second_index") != 0) {
		hyp_trilep_second_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_second_index"));
		hyp_trilep_second_index_branch->SetAddress(&hyp_trilep_second_index);
	}
	if(hyp_trilep_second_index_branch == 0 ) {
	cout << "Branch hyp_trilep_second_index does not exist." << endl;
	}
	hyp_trilep_third_index_branch = 0;
	if (tree->GetAlias("hyp_trilep_third_index") != 0) {
		hyp_trilep_third_index_branch = tree->GetBranch(tree->GetAlias("hyp_trilep_third_index"));
		hyp_trilep_third_index_branch->SetAddress(&hyp_trilep_third_index);
	}
	if(hyp_trilep_third_index_branch == 0 ) {
	cout << "Branch hyp_trilep_third_index does not exist." << endl;
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int index) {
	if (l1met_p4_branch != 0) {
	l1met_p4_branch->GetEntry(index);
	}
	if (els_p4_branch != 0) {
	els_p4_branch->GetEntry(index);
	}
	if (els_p4In_branch != 0) {
	els_p4In_branch->GetEntry(index);
	}
	if (els_p4Out_branch != 0) {
	els_p4Out_branch->GetEntry(index);
	}
	if (els_trk_p4_branch != 0) {
	els_trk_p4_branch->GetEntry(index);
	}
	if (genps_p4_branch != 0) {
	genps_p4_branch->GetEntry(index);
	}
	if (genps_prod_vtx_branch != 0) {
	genps_prod_vtx_branch->GetEntry(index);
	}
	if (hyp_ll_mc_p4_branch != 0) {
	hyp_ll_mc_p4_branch->GetEntry(index);
	}
	if (hyp_ll_p4_branch != 0) {
	hyp_ll_p4_branch->GetEntry(index);
	}
	if (hyp_ll_trk_p4_branch != 0) {
	hyp_ll_trk_p4_branch->GetEntry(index);
	}
	if (hyp_lt_mc_p4_branch != 0) {
	hyp_lt_mc_p4_branch->GetEntry(index);
	}
	if (hyp_lt_p4_branch != 0) {
	hyp_lt_p4_branch->GetEntry(index);
	}
	if (hyp_lt_trk_p4_branch != 0) {
	hyp_lt_trk_p4_branch->GetEntry(index);
	}
	if (hyp_p4_branch != 0) {
	hyp_p4_branch->GetEntry(index);
	}
	if (jets_p4_branch != 0) {
	jets_p4_branch->GetEntry(index);
	}
	if (l1emiso_p4_branch != 0) {
	l1emiso_p4_branch->GetEntry(index);
	}
	if (l1emnoiso_p4_branch != 0) {
	l1emnoiso_p4_branch->GetEntry(index);
	}
	if (l1jetsc_p4_branch != 0) {
	l1jetsc_p4_branch->GetEntry(index);
	}
	if (l1jetsf_p4_branch != 0) {
	l1jetsf_p4_branch->GetEntry(index);
	}
	if (l1jetst_p4_branch != 0) {
	l1jetst_p4_branch->GetEntry(index);
	}
	if (l1mus_p4_branch != 0) {
	l1mus_p4_branch->GetEntry(index);
	}
	if (mus_p4_branch != 0) {
	mus_p4_branch->GetEntry(index);
	}
	if (mus_trk_p4_branch != 0) {
	mus_trk_p4_branch->GetEntry(index);
	}
	if (els_tq_genMotherP4_branch != 0) {
	els_tq_genMotherP4_branch->GetEntry(index);
	}
	if (els_tq_genP4_branch != 0) {
	els_tq_genP4_branch->GetEntry(index);
	}
	if (jets_tq_genJet_p4_branch != 0) {
	jets_tq_genJet_p4_branch->GetEntry(index);
	}
	if (jets_tq_genPartonMother_p4_branch != 0) {
	jets_tq_genPartonMother_p4_branch->GetEntry(index);
	}
	if (jets_tq_genParton_p4_branch != 0) {
	jets_tq_genParton_p4_branch->GetEntry(index);
	}
	if (jets_tq_jet_p4_branch != 0) {
	jets_tq_jet_p4_branch->GetEntry(index);
	}
	if (mus_tq_genMotherP4_branch != 0) {
	mus_tq_genMotherP4_branch->GetEntry(index);
	}
	if (mus_tq_genP4_branch != 0) {
	mus_tq_genP4_branch->GetEntry(index);
	}
	if (trks_trk_p4_branch != 0) {
	trks_trk_p4_branch->GetEntry(index);
	}
	if (hyp_jets_mc_gp_p4_branch != 0) {
	hyp_jets_mc_gp_p4_branch->GetEntry(index);
	}
	if (hyp_jets_mc_p4_branch != 0) {
	hyp_jets_mc_p4_branch->GetEntry(index);
	}
	if (hyp_jets_p4_branch != 0) {
	hyp_jets_p4_branch->GetEntry(index);
	}
	if (hyp_jets_tq_p4_branch != 0) {
	hyp_jets_tq_p4_branch->GetEntry(index);
	}
	if (hyp_jets_tq_genPartonMother_p4_branch != 0) {
	hyp_jets_tq_genPartonMother_p4_branch->GetEntry(index);
	}
	if (hyp_jets_tq_genParton_p4_branch != 0) {
	hyp_jets_tq_genParton_p4_branch->GetEntry(index);
	}
	if (hyp_jets_tq_jet_p4_branch != 0) {
	hyp_jets_tq_jet_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_gp_p4_branch != 0) {
	hyp_other_jets_mc_gp_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_p4_branch != 0) {
	hyp_other_jets_mc_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_p4_branch != 0) {
	hyp_other_jets_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_genJet_p4_branch != 0) {
	hyp_other_jets_tq_genJet_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_genPartonMother_p4_branch != 0) {
	hyp_other_jets_tq_genPartonMother_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_genParton_p4_branch != 0) {
	hyp_other_jets_tq_genParton_p4_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_jet_p4_branch != 0) {
	hyp_other_jets_tq_jet_p4_branch->GetEntry(index);
	}
	if (evt_kfactor_branch != 0) {
	evt_kfactor_branch->GetEntry(index);
	}
	if (evt_weight_branch != 0) {
	evt_weight_branch->GetEntry(index);
	}
	if (evt_xsec_excl_branch != 0) {
	evt_xsec_excl_branch->GetEntry(index);
	}
	if (evt_xsec_incl_branch != 0) {
	evt_xsec_incl_branch->GetEntry(index);
	}
	if (l1met_etHad_branch != 0) {
	l1met_etHad_branch->GetEntry(index);
	}
	if (l1met_etTot_branch != 0) {
	l1met_etTot_branch->GetEntry(index);
	}
	if (l1met_met_branch != 0) {
	l1met_met_branch->GetEntry(index);
	}
	if (evt_met_branch != 0) {
	evt_met_branch->GetEntry(index);
	}
	if (evt_metPhi_branch != 0) {
	evt_metPhi_branch->GetEntry(index);
	}
	if (evt_metSig_branch != 0) {
	evt_metSig_branch->GetEntry(index);
	}
	if (evt_met_jetcorr_branch != 0) {
	evt_met_jetcorr_branch->GetEntry(index);
	}
	if (metphi_jetcorr_branch != 0) {
	metphi_jetcorr_branch->GetEntry(index);
	}
	if (els_musdr_branch != 0) {
	els_musdr_branch->GetEntry(index);
	}
	if (els_trkdr_branch != 0) {
	els_trkdr_branch->GetEntry(index);
	}
	if (els_ESc_branch != 0) {
	els_ESc_branch->GetEntry(index);
	}
	if (els_ESc_raw_branch != 0) {
	els_ESc_raw_branch->GetEntry(index);
	}
	if (els_ESeed_branch != 0) {
	els_ESeed_branch->GetEntry(index);
	}
	if (els_chi2_branch != 0) {
	els_chi2_branch->GetEntry(index);
	}
	if (els_d0_branch != 0) {
	els_d0_branch->GetEntry(index);
	}
	if (els_d0Err_branch != 0) {
	els_d0Err_branch->GetEntry(index);
	}
	if (els_dEtaIn_branch != 0) {
	els_dEtaIn_branch->GetEntry(index);
	}
	if (els_dEtaOut_branch != 0) {
	els_dEtaOut_branch->GetEntry(index);
	}
	if (els_dPhiIn_branch != 0) {
	els_dPhiIn_branch->GetEntry(index);
	}
	if (els_dPhiInPhiOut_branch != 0) {
	els_dPhiInPhiOut_branch->GetEntry(index);
	}
	if (els_dPhiOut_branch != 0) {
	els_dPhiOut_branch->GetEntry(index);
	}
	if (els_e3x3_branch != 0) {
	els_e3x3_branch->GetEntry(index);
	}
	if (els_e5x5_branch != 0) {
	els_e5x5_branch->GetEntry(index);
	}
	if (els_eOverPIn_branch != 0) {
	els_eOverPIn_branch->GetEntry(index);
	}
	if (els_eSeedOverPOut_branch != 0) {
	els_eSeedOverPOut_branch->GetEntry(index);
	}
	if (els_etaErr_branch != 0) {
	els_etaErr_branch->GetEntry(index);
	}
	if (els_fBrem_branch != 0) {
	els_fBrem_branch->GetEntry(index);
	}
	if (els_hOverE_branch != 0) {
	els_hOverE_branch->GetEntry(index);
	}
	if (els_ndof_branch != 0) {
	els_ndof_branch->GetEntry(index);
	}
	if (els_outerEta_branch != 0) {
	els_outerEta_branch->GetEntry(index);
	}
	if (els_outerPhi_branch != 0) {
	els_outerPhi_branch->GetEntry(index);
	}
	if (els_phiErr_branch != 0) {
	els_phiErr_branch->GetEntry(index);
	}
	if (els_ptErr_branch != 0) {
	els_ptErr_branch->GetEntry(index);
	}
	if (els_sigmaEtaEta_branch != 0) {
	els_sigmaEtaEta_branch->GetEntry(index);
	}
	if (els_sigmaPhiPhi_branch != 0) {
	els_sigmaPhiPhi_branch->GetEntry(index);
	}
	if (els_tkIso_branch != 0) {
	els_tkIso_branch->GetEntry(index);
	}
	if (els_vertexphi_branch != 0) {
	els_vertexphi_branch->GetEntry(index);
	}
	if (els_z0_branch != 0) {
	els_z0_branch->GetEntry(index);
	}
	if (els_z0Err_branch != 0) {
	els_z0Err_branch->GetEntry(index);
	}
	if (hyp_ll_chi2_branch != 0) {
	hyp_ll_chi2_branch->GetEntry(index);
	}
	if (hyp_ll_d0_branch != 0) {
	hyp_ll_d0_branch->GetEntry(index);
	}
	if (hyp_ll_d0Err_branch != 0) {
	hyp_ll_d0Err_branch->GetEntry(index);
	}
	if (hyp_ll_etaErr_branch != 0) {
	hyp_ll_etaErr_branch->GetEntry(index);
	}
	if (hyp_ll_iso_branch != 0) {
	hyp_ll_iso_branch->GetEntry(index);
	}
	if (hyp_ll_ndof_branch != 0) {
	hyp_ll_ndof_branch->GetEntry(index);
	}
	if (hyp_ll_outerEta_branch != 0) {
	hyp_ll_outerEta_branch->GetEntry(index);
	}
	if (hyp_ll_outerPhi_branch != 0) {
	hyp_ll_outerPhi_branch->GetEntry(index);
	}
	if (hyp_ll_phiErr_branch != 0) {
	hyp_ll_phiErr_branch->GetEntry(index);
	}
	if (hyp_ll_ptErr_branch != 0) {
	hyp_ll_ptErr_branch->GetEntry(index);
	}
	if (hyp_ll_tkIso_branch != 0) {
	hyp_ll_tkIso_branch->GetEntry(index);
	}
	if (hyp_ll_vertexphi_branch != 0) {
	hyp_ll_vertexphi_branch->GetEntry(index);
	}
	if (hyp_ll_z0_branch != 0) {
	hyp_ll_z0_branch->GetEntry(index);
	}
	if (hyp_ll_z0Err_branch != 0) {
	hyp_ll_z0Err_branch->GetEntry(index);
	}
	if (hyp_lt_chi2_branch != 0) {
	hyp_lt_chi2_branch->GetEntry(index);
	}
	if (hyp_lt_d0_branch != 0) {
	hyp_lt_d0_branch->GetEntry(index);
	}
	if (hyp_lt_d0Err_branch != 0) {
	hyp_lt_d0Err_branch->GetEntry(index);
	}
	if (hyp_lt_etaErr_branch != 0) {
	hyp_lt_etaErr_branch->GetEntry(index);
	}
	if (hyp_lt_iso_branch != 0) {
	hyp_lt_iso_branch->GetEntry(index);
	}
	if (hyp_lt_ndof_branch != 0) {
	hyp_lt_ndof_branch->GetEntry(index);
	}
	if (hyp_lt_outerEta_branch != 0) {
	hyp_lt_outerEta_branch->GetEntry(index);
	}
	if (hyp_lt_outerPhi_branch != 0) {
	hyp_lt_outerPhi_branch->GetEntry(index);
	}
	if (hyp_lt_phiErr_branch != 0) {
	hyp_lt_phiErr_branch->GetEntry(index);
	}
	if (hyp_lt_ptErr_branch != 0) {
	hyp_lt_ptErr_branch->GetEntry(index);
	}
	if (hyp_lt_tkIso_branch != 0) {
	hyp_lt_tkIso_branch->GetEntry(index);
	}
	if (hyp_lt_vertexphi_branch != 0) {
	hyp_lt_vertexphi_branch->GetEntry(index);
	}
	if (hyp_lt_z0_branch != 0) {
	hyp_lt_z0_branch->GetEntry(index);
	}
	if (hyp_lt_z0Err_branch != 0) {
	hyp_lt_z0Err_branch->GetEntry(index);
	}
	if (hyp_met_branch != 0) {
	hyp_met_branch->GetEntry(index);
	}
	if (hyp_metAll_branch != 0) {
	hyp_metAll_branch->GetEntry(index);
	}
	if (hyp_metAllCaloExp_branch != 0) {
	hyp_metAllCaloExp_branch->GetEntry(index);
	}
	if (hyp_metCaloExp_branch != 0) {
	hyp_metCaloExp_branch->GetEntry(index);
	}
	if (hyp_metCone_branch != 0) {
	hyp_metCone_branch->GetEntry(index);
	}
	if (hyp_metDPhiJet10_branch != 0) {
	hyp_metDPhiJet10_branch->GetEntry(index);
	}
	if (hyp_metDPhiJet15_branch != 0) {
	hyp_metDPhiJet15_branch->GetEntry(index);
	}
	if (hyp_metDPhiJet20_branch != 0) {
	hyp_metDPhiJet20_branch->GetEntry(index);
	}
	if (hyp_metDPhiTrk10_branch != 0) {
	hyp_metDPhiTrk10_branch->GetEntry(index);
	}
	if (hyp_metDPhiTrk25_branch != 0) {
	hyp_metDPhiTrk25_branch->GetEntry(index);
	}
	if (hyp_metDPhiTrk50_branch != 0) {
	hyp_metDPhiTrk50_branch->GetEntry(index);
	}
	if (hyp_metJes10_branch != 0) {
	hyp_metJes10_branch->GetEntry(index);
	}
	if (hyp_metJes15_branch != 0) {
	hyp_metJes15_branch->GetEntry(index);
	}
	if (hyp_metJes30_branch != 0) {
	hyp_metJes30_branch->GetEntry(index);
	}
	if (hyp_metJes5_branch != 0) {
	hyp_metJes5_branch->GetEntry(index);
	}
	if (hyp_metJes50_branch != 0) {
	hyp_metJes50_branch->GetEntry(index);
	}
	if (hyp_metNoCalo_branch != 0) {
	hyp_metNoCalo_branch->GetEntry(index);
	}
	if (hyp_metPhi_branch != 0) {
	hyp_metPhi_branch->GetEntry(index);
	}
	if (hyp_metPhiAll_branch != 0) {
	hyp_metPhiAll_branch->GetEntry(index);
	}
	if (hyp_metPhiAllCaloExp_branch != 0) {
	hyp_metPhiAllCaloExp_branch->GetEntry(index);
	}
	if (hyp_metPhiCaloExp_branch != 0) {
	hyp_metPhiCaloExp_branch->GetEntry(index);
	}
	if (hyp_metPhiCone_branch != 0) {
	hyp_metPhiCone_branch->GetEntry(index);
	}
	if (hyp_metPhiJes10_branch != 0) {
	hyp_metPhiJes10_branch->GetEntry(index);
	}
	if (hyp_metPhiJes15_branch != 0) {
	hyp_metPhiJes15_branch->GetEntry(index);
	}
	if (hyp_metPhiJes30_branch != 0) {
	hyp_metPhiJes30_branch->GetEntry(index);
	}
	if (hyp_metPhiJes5_branch != 0) {
	hyp_metPhiJes5_branch->GetEntry(index);
	}
	if (hyp_metPhiJes50_branch != 0) {
	hyp_metPhiJes50_branch->GetEntry(index);
	}
	if (hyp_metPhiNoCalo_branch != 0) {
	hyp_metPhiNoCalo_branch->GetEntry(index);
	}
	if (hyp_quadlep_met_branch != 0) {
	hyp_quadlep_met_branch->GetEntry(index);
	}
	if (hyp_quadlep_metAll_branch != 0) {
	hyp_quadlep_metAll_branch->GetEntry(index);
	}
	if (hyp_trilep_met_branch != 0) {
	hyp_trilep_met_branch->GetEntry(index);
	}
	if (hyp_trilep_metAll_branch != 0) {
	hyp_trilep_metAll_branch->GetEntry(index);
	}
	if (jets_EMFcor_branch != 0) {
	jets_EMFcor_branch->GetEntry(index);
	}
	if (jets_chFrac_branch != 0) {
	jets_chFrac_branch->GetEntry(index);
	}
	if (jets_cor_branch != 0) {
	jets_cor_branch->GetEntry(index);
	}
	if (jets_emFrac_branch != 0) {
	jets_emFrac_branch->GetEntry(index);
	}
	if (mus_eledr_branch != 0) {
	mus_eledr_branch->GetEntry(index);
	}
	if (mus_jetdr_branch != 0) {
	mus_jetdr_branch->GetEntry(index);
	}
	if (mus_trkdr_branch != 0) {
	mus_trkdr_branch->GetEntry(index);
	}
	if (mus_chi2_branch != 0) {
	mus_chi2_branch->GetEntry(index);
	}
	if (mus_d0_branch != 0) {
	mus_d0_branch->GetEntry(index);
	}
	if (mus_d0Err_branch != 0) {
	mus_d0Err_branch->GetEntry(index);
	}
	if (mus_e_em_branch != 0) {
	mus_e_em_branch->GetEntry(index);
	}
	if (mus_e_emS9_branch != 0) {
	mus_e_emS9_branch->GetEntry(index);
	}
	if (mus_e_had_branch != 0) {
	mus_e_had_branch->GetEntry(index);
	}
	if (mus_e_hadS9_branch != 0) {
	mus_e_hadS9_branch->GetEntry(index);
	}
	if (mus_e_ho_branch != 0) {
	mus_e_ho_branch->GetEntry(index);
	}
	if (mus_e_hoS9_branch != 0) {
	mus_e_hoS9_branch->GetEntry(index);
	}
	if (mus_etaErr_branch != 0) {
	mus_etaErr_branch->GetEntry(index);
	}
	if (mus_gfit_chi2_branch != 0) {
	mus_gfit_chi2_branch->GetEntry(index);
	}
	if (mus_gfit_ndof_branch != 0) {
	mus_gfit_ndof_branch->GetEntry(index);
	}
	if (mus_iso_branch != 0) {
	mus_iso_branch->GetEntry(index);
	}
	if (mus_iso03_emEt_branch != 0) {
	mus_iso03_emEt_branch->GetEntry(index);
	}
	if (mus_iso03_hadEt_branch != 0) {
	mus_iso03_hadEt_branch->GetEntry(index);
	}
	if (mus_iso03_hoEt_branch != 0) {
	mus_iso03_hoEt_branch->GetEntry(index);
	}
	if (mus_iso03_sumPt_branch != 0) {
	mus_iso03_sumPt_branch->GetEntry(index);
	}
	if (mus_iso05_emEt_branch != 0) {
	mus_iso05_emEt_branch->GetEntry(index);
	}
	if (mus_iso05_hadEt_branch != 0) {
	mus_iso05_hadEt_branch->GetEntry(index);
	}
	if (mus_iso05_hoEt_branch != 0) {
	mus_iso05_hoEt_branch->GetEntry(index);
	}
	if (mus_iso05_sumPt_branch != 0) {
	mus_iso05_sumPt_branch->GetEntry(index);
	}
	if (mus_ndof_branch != 0) {
	mus_ndof_branch->GetEntry(index);
	}
	if (mus_outerEta_branch != 0) {
	mus_outerEta_branch->GetEntry(index);
	}
	if (mus_outerPhi_branch != 0) {
	mus_outerPhi_branch->GetEntry(index);
	}
	if (mus_phiErr_branch != 0) {
	mus_phiErr_branch->GetEntry(index);
	}
	if (mus_ptErr_branch != 0) {
	mus_ptErr_branch->GetEntry(index);
	}
	if (mus_vertexphi_branch != 0) {
	mus_vertexphi_branch->GetEntry(index);
	}
	if (mus_z0_branch != 0) {
	mus_z0_branch->GetEntry(index);
	}
	if (mus_z0Err_branch != 0) {
	mus_z0Err_branch->GetEntry(index);
	}
	if (els_tq_LRComb_branch != 0) {
	els_tq_LRComb_branch->GetEntry(index);
	}
	if (els_tq_caloIso_branch != 0) {
	els_tq_caloIso_branch->GetEntry(index);
	}
	if (els_tq_egammaEcalIso_branch != 0) {
	els_tq_egammaEcalIso_branch->GetEntry(index);
	}
	if (els_tq_egammaHcalIso_branch != 0) {
	els_tq_egammaHcalIso_branch->GetEntry(index);
	}
	if (els_tq_egammaTkIso_branch != 0) {
	els_tq_egammaTkIso_branch->GetEntry(index);
	}
	if (els_tq_electronIDRobust_branch != 0) {
	els_tq_electronIDRobust_branch->GetEntry(index);
	}
	if (els_tq_leptonID_branch != 0) {
	els_tq_leptonID_branch->GetEntry(index);
	}
	if (els_tq_trackIso_branch != 0) {
	els_tq_trackIso_branch->GetEntry(index);
	}
	if (jets_tq_bCorrF_branch != 0) {
	jets_tq_bCorrF_branch->GetEntry(index);
	}
	if (jets_tq_cCorrF_branch != 0) {
	jets_tq_cCorrF_branch->GetEntry(index);
	}
	if (jets_tq_gluCorrF_branch != 0) {
	jets_tq_gluCorrF_branch->GetEntry(index);
	}
	if (jets_tq_jetCharge_branch != 0) {
	jets_tq_jetCharge_branch->GetEntry(index);
	}
	if (jets_tq_noCorrF_branch != 0) {
	jets_tq_noCorrF_branch->GetEntry(index);
	}
	if (jets_tq_udsCorrF_branch != 0) {
	jets_tq_udsCorrF_branch->GetEntry(index);
	}
	if (mus_tq_caloIso_branch != 0) {
	mus_tq_caloIso_branch->GetEntry(index);
	}
	if (mus_tq_leptonID_branch != 0) {
	mus_tq_leptonID_branch->GetEntry(index);
	}
	if (mus_tq_lrComb_branch != 0) {
	mus_tq_lrComb_branch->GetEntry(index);
	}
	if (mus_tq_trackIso_branch != 0) {
	mus_tq_trackIso_branch->GetEntry(index);
	}
	if (trks_chi2_branch != 0) {
	trks_chi2_branch->GetEntry(index);
	}
	if (trks_d0_branch != 0) {
	trks_d0_branch->GetEntry(index);
	}
	if (trks_d0Err_branch != 0) {
	trks_d0Err_branch->GetEntry(index);
	}
	if (trks_etaErr_branch != 0) {
	trks_etaErr_branch->GetEntry(index);
	}
	if (trks_ndof_branch != 0) {
	trks_ndof_branch->GetEntry(index);
	}
	if (trks_outerEta_branch != 0) {
	trks_outerEta_branch->GetEntry(index);
	}
	if (trks_outerPhi_branch != 0) {
	trks_outerPhi_branch->GetEntry(index);
	}
	if (trks_phiErr_branch != 0) {
	trks_phiErr_branch->GetEntry(index);
	}
	if (trks_ptErr_branch != 0) {
	trks_ptErr_branch->GetEntry(index);
	}
	if (trks_vertexphi_branch != 0) {
	trks_vertexphi_branch->GetEntry(index);
	}
	if (trks_z0_branch != 0) {
	trks_z0_branch->GetEntry(index);
	}
	if (trks_z0Err_branch != 0) {
	trks_z0Err_branch->GetEntry(index);
	}
	if (hyp_jets_EMFcor_branch != 0) {
	hyp_jets_EMFcor_branch->GetEntry(index);
	}
	if (hyp_jets_chFrac_branch != 0) {
	hyp_jets_chFrac_branch->GetEntry(index);
	}
	if (hyp_jets_cor_branch != 0) {
	hyp_jets_cor_branch->GetEntry(index);
	}
	if (hyp_jets_emFrac_branch != 0) {
	hyp_jets_emFrac_branch->GetEntry(index);
	}
	if (hyp_jets_mc_emEnergy_branch != 0) {
	hyp_jets_mc_emEnergy_branch->GetEntry(index);
	}
	if (hyp_jets_mc_hadEnergy_branch != 0) {
	hyp_jets_mc_hadEnergy_branch->GetEntry(index);
	}
	if (hyp_jets_mc_invEnergy_branch != 0) {
	hyp_jets_mc_invEnergy_branch->GetEntry(index);
	}
	if (hyp_jets_mc_otherEnergy_branch != 0) {
	hyp_jets_mc_otherEnergy_branch->GetEntry(index);
	}
	if (hyp_jets_tq_bCorrF_branch != 0) {
	hyp_jets_tq_bCorrF_branch->GetEntry(index);
	}
	if (hyp_jets_tq_cCorrF_branch != 0) {
	hyp_jets_tq_cCorrF_branch->GetEntry(index);
	}
	if (hyp_jets_tq_gluCorrF_branch != 0) {
	hyp_jets_tq_gluCorrF_branch->GetEntry(index);
	}
	if (hyp_jets_tq_jetCharge_branch != 0) {
	hyp_jets_tq_jetCharge_branch->GetEntry(index);
	}
	if (hyp_jets_tq_noCorrF_branch != 0) {
	hyp_jets_tq_noCorrF_branch->GetEntry(index);
	}
	if (hyp_jets_tq_udsCorrF_branch != 0) {
	hyp_jets_tq_udsCorrF_branch->GetEntry(index);
	}
	if (hyp_other_jets_EMFcor_branch != 0) {
	hyp_other_jets_EMFcor_branch->GetEntry(index);
	}
	if (hyp_other_jets_chFrac_branch != 0) {
	hyp_other_jets_chFrac_branch->GetEntry(index);
	}
	if (hyp_other_jets_cor_branch != 0) {
	hyp_other_jets_cor_branch->GetEntry(index);
	}
	if (hyp_other_jets_emFrac_branch != 0) {
	hyp_other_jets_emFrac_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_emEnergy_branch != 0) {
	hyp_other_jets_mc_emEnergy_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_hadEnergy_branch != 0) {
	hyp_other_jets_mc_hadEnergy_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_invEnergy_branch != 0) {
	hyp_other_jets_mc_invEnergy_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_otherEnergy_branch != 0) {
	hyp_other_jets_mc_otherEnergy_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_bCorrF_branch != 0) {
	hyp_other_jets_tq_bCorrF_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_cCorrF_branch != 0) {
	hyp_other_jets_tq_cCorrF_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_gluCorrF_branch != 0) {
	hyp_other_jets_tq_gluCorrF_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_jetCharge_branch != 0) {
	hyp_other_jets_tq_jetCharge_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_noCorrF_branch != 0) {
	hyp_other_jets_tq_noCorrF_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_udsCorrF_branch != 0) {
	hyp_other_jets_tq_udsCorrF_branch->GetEntry(index);
	}
	if (evt_HLT1_branch != 0) {
	evt_HLT1_branch->GetEntry(index);
	}
	if (evt_HLT2_branch != 0) {
	evt_HLT2_branch->GetEntry(index);
	}
	if (evt_HLT3_branch != 0) {
	evt_HLT3_branch->GetEntry(index);
	}
	if (evt_HLT4_branch != 0) {
	evt_HLT4_branch->GetEntry(index);
	}
	if (evt_L1_1_branch != 0) {
	evt_L1_1_branch->GetEntry(index);
	}
	if (evt_L1_2_branch != 0) {
	evt_L1_2_branch->GetEntry(index);
	}
	if (evt_L1_3_branch != 0) {
	evt_L1_3_branch->GetEntry(index);
	}
	if (evt_L1_4_branch != 0) {
	evt_L1_4_branch->GetEntry(index);
	}
	if (evt_event_branch != 0) {
	evt_event_branch->GetEntry(index);
	}
	if (evt_run_branch != 0) {
	evt_run_branch->GetEntry(index);
	}
	if (evt_nl1emiso_branch != 0) {
	evt_nl1emiso_branch->GetEntry(index);
	}
	if (evt_nl1emnoiso_branch != 0) {
	evt_nl1emnoiso_branch->GetEntry(index);
	}
	if (evt_nl1jetsc_branch != 0) {
	evt_nl1jetsc_branch->GetEntry(index);
	}
	if (evt_nl1jetsf_branch != 0) {
	evt_nl1jetsf_branch->GetEntry(index);
	}
	if (evt_nl1jetst_branch != 0) {
	evt_nl1jetst_branch->GetEntry(index);
	}
	if (evt_nl1mus_branch != 0) {
	evt_nl1mus_branch->GetEntry(index);
	}
	if (els_closestMuon_branch != 0) {
	els_closestMuon_branch->GetEntry(index);
	}
	if (els_trkidx_branch != 0) {
	els_trkidx_branch->GetEntry(index);
	}
	if (els_charge_branch != 0) {
	els_charge_branch->GetEntry(index);
	}
	if (els_class_branch != 0) {
	els_class_branch->GetEntry(index);
	}
	if (els_looseId_branch != 0) {
	els_looseId_branch->GetEntry(index);
	}
	if (els_lostHits_branch != 0) {
	els_lostHits_branch->GetEntry(index);
	}
	if (els_nSeed_branch != 0) {
	els_nSeed_branch->GetEntry(index);
	}
	if (els_pass3looseId_branch != 0) {
	els_pass3looseId_branch->GetEntry(index);
	}
	if (els_pass3simpleId_branch != 0) {
	els_pass3simpleId_branch->GetEntry(index);
	}
	if (els_pass3tightId_branch != 0) {
	els_pass3tightId_branch->GetEntry(index);
	}
	if (els_robustId_branch != 0) {
	els_robustId_branch->GetEntry(index);
	}
	if (els_simpleIdPlus_branch != 0) {
	els_simpleIdPlus_branch->GetEntry(index);
	}
	if (els_tightId_branch != 0) {
	els_tightId_branch->GetEntry(index);
	}
	if (els_validHits_branch != 0) {
	els_validHits_branch->GetEntry(index);
	}
	if (genps_id_branch != 0) {
	genps_id_branch->GetEntry(index);
	}
	if (genps_id_mother_branch != 0) {
	genps_id_mother_branch->GetEntry(index);
	}
	if (genps_status_branch != 0) {
	genps_status_branch->GetEntry(index);
	}
	if (hyp_ll_charge_branch != 0) {
	hyp_ll_charge_branch->GetEntry(index);
	}
	if (hyp_ll_id_branch != 0) {
	hyp_ll_id_branch->GetEntry(index);
	}
	if (hyp_ll_index_branch != 0) {
	hyp_ll_index_branch->GetEntry(index);
	}
	if (hyp_ll_lostHits_branch != 0) {
	hyp_ll_lostHits_branch->GetEntry(index);
	}
	if (hyp_ll_mc_id_branch != 0) {
	hyp_ll_mc_id_branch->GetEntry(index);
	}
	if (hyp_ll_mc_motherid_branch != 0) {
	hyp_ll_mc_motherid_branch->GetEntry(index);
	}
	if (hyp_ll_validHits_branch != 0) {
	hyp_ll_validHits_branch->GetEntry(index);
	}
	if (hyp_lt_charge_branch != 0) {
	hyp_lt_charge_branch->GetEntry(index);
	}
	if (hyp_lt_id_branch != 0) {
	hyp_lt_id_branch->GetEntry(index);
	}
	if (hyp_lt_index_branch != 0) {
	hyp_lt_index_branch->GetEntry(index);
	}
	if (hyp_lt_lostHits_branch != 0) {
	hyp_lt_lostHits_branch->GetEntry(index);
	}
	if (hyp_lt_mc_id_branch != 0) {
	hyp_lt_mc_id_branch->GetEntry(index);
	}
	if (hyp_lt_mc_motherid_branch != 0) {
	hyp_lt_mc_motherid_branch->GetEntry(index);
	}
	if (hyp_lt_validHits_branch != 0) {
	hyp_lt_validHits_branch->GetEntry(index);
	}
	if (hyp_njets_branch != 0) {
	hyp_njets_branch->GetEntry(index);
	}
	if (hyp_nojets_branch != 0) {
	hyp_nojets_branch->GetEntry(index);
	}
	if (hyp_type_branch != 0) {
	hyp_type_branch->GetEntry(index);
	}
	if (hyp_quadlep_first_type_branch != 0) {
	hyp_quadlep_first_type_branch->GetEntry(index);
	}
	if (hyp_quadlep_fourth_type_branch != 0) {
	hyp_quadlep_fourth_type_branch->GetEntry(index);
	}
	if (hyp_quadlep_second_type_branch != 0) {
	hyp_quadlep_second_type_branch->GetEntry(index);
	}
	if (hyp_quadlep_third_type_branch != 0) {
	hyp_quadlep_third_type_branch->GetEntry(index);
	}
	if (hyp_trilep_first_type_branch != 0) {
	hyp_trilep_first_type_branch->GetEntry(index);
	}
	if (hyp_trilep_second_type_branch != 0) {
	hyp_trilep_second_type_branch->GetEntry(index);
	}
	if (hyp_trilep_third_type_branch != 0) {
	hyp_trilep_third_type_branch->GetEntry(index);
	}
	if (jets_closestElectron_branch != 0) {
	jets_closestElectron_branch->GetEntry(index);
	}
	if (jets_closestMuon_branch != 0) {
	jets_closestMuon_branch->GetEntry(index);
	}
	if (l1emiso_ieta_branch != 0) {
	l1emiso_ieta_branch->GetEntry(index);
	}
	if (l1emiso_iphi_branch != 0) {
	l1emiso_iphi_branch->GetEntry(index);
	}
	if (l1emiso_rawId_branch != 0) {
	l1emiso_rawId_branch->GetEntry(index);
	}
	if (l1emiso_type_branch != 0) {
	l1emiso_type_branch->GetEntry(index);
	}
	if (l1emnoiso_ieta_branch != 0) {
	l1emnoiso_ieta_branch->GetEntry(index);
	}
	if (l1emnoiso_iphi_branch != 0) {
	l1emnoiso_iphi_branch->GetEntry(index);
	}
	if (l1emnoiso_rawId_branch != 0) {
	l1emnoiso_rawId_branch->GetEntry(index);
	}
	if (l1emnoiso_type_branch != 0) {
	l1emnoiso_type_branch->GetEntry(index);
	}
	if (l1jetsc_ieta_branch != 0) {
	l1jetsc_ieta_branch->GetEntry(index);
	}
	if (l1jetsc_iphi_branch != 0) {
	l1jetsc_iphi_branch->GetEntry(index);
	}
	if (l1jetsc_rawId_branch != 0) {
	l1jetsc_rawId_branch->GetEntry(index);
	}
	if (l1jetsc_type_branch != 0) {
	l1jetsc_type_branch->GetEntry(index);
	}
	if (l1jetsf_ieta_branch != 0) {
	l1jetsf_ieta_branch->GetEntry(index);
	}
	if (l1jetsf_iphi_branch != 0) {
	l1jetsf_iphi_branch->GetEntry(index);
	}
	if (l1jetsf_rawId_branch != 0) {
	l1jetsf_rawId_branch->GetEntry(index);
	}
	if (l1jetsf_type_branch != 0) {
	l1jetsf_type_branch->GetEntry(index);
	}
	if (l1jetst_ieta_branch != 0) {
	l1jetst_ieta_branch->GetEntry(index);
	}
	if (l1jetst_iphi_branch != 0) {
	l1jetst_iphi_branch->GetEntry(index);
	}
	if (l1jetst_rawId_branch != 0) {
	l1jetst_rawId_branch->GetEntry(index);
	}
	if (l1jetst_type_branch != 0) {
	l1jetst_type_branch->GetEntry(index);
	}
	if (l1mus_flags_branch != 0) {
	l1mus_flags_branch->GetEntry(index);
	}
	if (l1mus_q_branch != 0) {
	l1mus_q_branch->GetEntry(index);
	}
	if (l1mus_qual_branch != 0) {
	l1mus_qual_branch->GetEntry(index);
	}
	if (l1mus_qualFlags_branch != 0) {
	l1mus_qualFlags_branch->GetEntry(index);
	}
	if (mus_closestEle_branch != 0) {
	mus_closestEle_branch->GetEntry(index);
	}
	if (mus_closestJet_branch != 0) {
	mus_closestJet_branch->GetEntry(index);
	}
	if (mus_trkidx_branch != 0) {
	mus_trkidx_branch->GetEntry(index);
	}
	if (mus_charge_branch != 0) {
	mus_charge_branch->GetEntry(index);
	}
	if (mus_gfit_validHits_branch != 0) {
	mus_gfit_validHits_branch->GetEntry(index);
	}
	if (mus_iso03_ntrk_branch != 0) {
	mus_iso03_ntrk_branch->GetEntry(index);
	}
	if (mus_iso05_ntrk_branch != 0) {
	mus_iso05_ntrk_branch->GetEntry(index);
	}
	if (mus_lostHits_branch != 0) {
	mus_lostHits_branch->GetEntry(index);
	}
	if (mus_nmatches_branch != 0) {
	mus_nmatches_branch->GetEntry(index);
	}
	if (mus_pid_TM2DCompatibilityLoose_branch != 0) {
	mus_pid_TM2DCompatibilityLoose_branch->GetEntry(index);
	}
	if (mus_pid_TM2DCompatibilityTight_branch != 0) {
	mus_pid_TM2DCompatibilityTight_branch->GetEntry(index);
	}
	if (mus_pid_TMLastStationLoose_branch != 0) {
	mus_pid_TMLastStationLoose_branch->GetEntry(index);
	}
	if (mus_pid_TMLastStationTight_branch != 0) {
	mus_pid_TMLastStationTight_branch->GetEntry(index);
	}
	if (mus_trkrefkey_branch != 0) {
	mus_trkrefkey_branch->GetEntry(index);
	}
	if (mus_validHits_branch != 0) {
	mus_validHits_branch->GetEntry(index);
	}
	if (els_tq_egammaTkNumIso_branch != 0) {
	els_tq_egammaTkNumIso_branch->GetEntry(index);
	}
	if (els_tq_genID_branch != 0) {
	els_tq_genID_branch->GetEntry(index);
	}
	if (els_tq_genMotherID_branch != 0) {
	els_tq_genMotherID_branch->GetEntry(index);
	}
	if (jets_tq_genPartonMother_id_branch != 0) {
	jets_tq_genPartonMother_id_branch->GetEntry(index);
	}
	if (jets_tq_genParton_id_branch != 0) {
	jets_tq_genParton_id_branch->GetEntry(index);
	}
	if (jets_tq_partonFlavour_branch != 0) {
	jets_tq_partonFlavour_branch->GetEntry(index);
	}
	if (mus_tq_genID_branch != 0) {
	mus_tq_genID_branch->GetEntry(index);
	}
	if (mus_tq_genMotherID_branch != 0) {
	mus_tq_genMotherID_branch->GetEntry(index);
	}
	if (trks_charge_branch != 0) {
	trks_charge_branch->GetEntry(index);
	}
	if (trks_lostHits_branch != 0) {
	trks_lostHits_branch->GetEntry(index);
	}
	if (trks_validHits_branch != 0) {
	trks_validHits_branch->GetEntry(index);
	}
	if (hyp_jets_mc_id_branch != 0) {
	hyp_jets_mc_id_branch->GetEntry(index);
	}
	if (hyp_jets_tq_genPartonMother_id_branch != 0) {
	hyp_jets_tq_genPartonMother_id_branch->GetEntry(index);
	}
	if (hyp_jets_tq_genParton_id_branch != 0) {
	hyp_jets_tq_genParton_id_branch->GetEntry(index);
	}
	if (hyp_jets_tq_partonFlavour_branch != 0) {
	hyp_jets_tq_partonFlavour_branch->GetEntry(index);
	}
	if (hyp_other_jets_mc_id_branch != 0) {
	hyp_other_jets_mc_id_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_genPartonMother_id_branch != 0) {
	hyp_other_jets_tq_genPartonMother_id_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_genParton_id_branch != 0) {
	hyp_other_jets_tq_genParton_id_branch->GetEntry(index);
	}
	if (hyp_other_jets_tq_partonFlavour_branch != 0) {
	hyp_other_jets_tq_partonFlavour_branch->GetEntry(index);
	}
	if (hyp_quadlep_jets_index_branch != 0) {
	hyp_quadlep_jets_index_branch->GetEntry(index);
	}
	if (hyp_trilep_jets_index_branch != 0) {
	hyp_trilep_jets_index_branch->GetEntry(index);
	}
	if (evt_nels_branch != 0) {
	evt_nels_branch->GetEntry(index);
	}
	if (evt_njets_branch != 0) {
	evt_njets_branch->GetEntry(index);
	}
	if (hyp_quadlep_bucket_branch != 0) {
	hyp_quadlep_bucket_branch->GetEntry(index);
	}
	if (hyp_quadlep_first_index_branch != 0) {
	hyp_quadlep_first_index_branch->GetEntry(index);
	}
	if (hyp_quadlep_fourth_index_branch != 0) {
	hyp_quadlep_fourth_index_branch->GetEntry(index);
	}
	if (hyp_quadlep_second_index_branch != 0) {
	hyp_quadlep_second_index_branch->GetEntry(index);
	}
	if (hyp_quadlep_third_index_branch != 0) {
	hyp_quadlep_third_index_branch->GetEntry(index);
	}
	if (hyp_trilep_bucket_branch != 0) {
	hyp_trilep_bucket_branch->GetEntry(index);
	}
	if (hyp_trilep_first_index_branch != 0) {
	hyp_trilep_first_index_branch->GetEntry(index);
	}
	if (hyp_trilep_second_index_branch != 0) {
	hyp_trilep_second_index_branch->GetEntry(index);
	}
	if (hyp_trilep_third_index_branch != 0) {
	hyp_trilep_third_index_branch->GetEntry(index);
	}
}
