	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> 	evt_bs_;
	TBranch *evt_bs_branch;
	bool evt_bs_isLoaded;
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> >	mus_gfit_outerPos_;
	TBranch *mus_gfit_outerPos_branch;
	bool mus_gfit_outerPos_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	els_mc_p4_;
	TBranch *els_mc_p4_branch;
	bool els_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_mc_gp_p4_;
	TBranch *jets_mc_gp_p4_branch;
	bool jets_mc_gp_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	jets_mc_p4_;
	TBranch *jets_mc_p4_branch;
	bool jets_mc_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_mc_p4_;
	TBranch *mus_mc_p4_branch;
	bool mus_mc_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_p4_;
	TBranch *mus_p4_branch;
	bool mus_p4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_trk_p4_;
	TBranch *mus_trk_p4_branch;
	bool mus_trk_p4_isLoaded;
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_pat_genMotherP4_;
	TBranch *mus_pat_genMotherP4_branch;
	bool mus_pat_genMotherP4_isLoaded;
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > >	mus_pat_genP4_;
	TBranch *mus_pat_genP4_branch;
	bool mus_pat_genP4_isLoaded;
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > >	hyp_jets_pat_p4_;
	TBranch *hyp_jets_pat_p4_branch;
	bool hyp_jets_pat_p4_isLoaded;
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
	vector<float>	mus_mcdr_;
	TBranch *mus_mcdr_branch;
	bool mus_mcdr_isLoaded;
	vector<float>	trk_mcdr_;
	TBranch *trk_mcdr_branch;
	bool trk_mcdr_isLoaded;
	vector<float>	els_musdr_;
	TBranch *els_musdr_branch;
	bool els_musdr_isLoaded;
	vector<float>	els_trkdr_;
	TBranch *els_trkdr_branch;
	bool els_trkdr_isLoaded;
	vector<float>	els_trkshFrac_;
	TBranch *els_trkshFrac_branch;
	bool els_trkshFrac_isLoaded;
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
	vector<float>	els_pat_robustLooseId_;
	TBranch *els_pat_robustLooseId_branch;
	bool els_pat_robustLooseId_isLoaded;
	vector<float>	els_pat_robustTightId_;
	TBranch *els_pat_robustTightId_branch;
	bool els_pat_robustTightId_isLoaded;
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
	vector<float>	jets_pat_gluCorrF_;
	TBranch *jets_pat_gluCorrF_branch;
	bool jets_pat_gluCorrF_isLoaded;
	vector<float>	jets_pat_jetCharge_;
	TBranch *jets_pat_jetCharge_branch;
	bool jets_pat_jetCharge_isLoaded;
	vector<float>	jets_pat_noCorrF_;
	TBranch *jets_pat_noCorrF_branch;
	bool jets_pat_noCorrF_isLoaded;
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
	vector<float>	mus_pat_leptonID_;
	TBranch *mus_pat_leptonID_branch;
	bool mus_pat_leptonID_isLoaded;
	vector<float>	mus_pat_trackIso_;
	TBranch *mus_pat_trackIso_branch;
	bool mus_pat_trackIso_isLoaded;
	vector<float>	mus_pat_vetoDep_;
	TBranch *mus_pat_vetoDep_branch;
	bool mus_pat_vetoDep_isLoaded;
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
	vector<int>	mus_mc_id_;
	TBranch *mus_mc_id_branch;
	bool mus_mc_id_isLoaded;
	vector<int>	mus_mcidx_;
	TBranch *mus_mcidx_branch;
	bool mus_mcidx_isLoaded;
	vector<int>	mus_mc_motherid_;
	TBranch *mus_mc_motherid_branch;
	bool mus_mc_motherid_isLoaded;
	vector<int>	trk_mc_id_;
	TBranch *trk_mc_id_branch;
	bool trk_mc_id_isLoaded;
	vector<int>	trk_mcidx_;
	TBranch *trk_mcidx_branch;
	bool trk_mcidx_isLoaded;
	vector<int>	trk_mc_motherid_;
	TBranch *trk_mc_motherid_branch;
	bool trk_mc_motherid_isLoaded;
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
	vector<unsigned int>	els_pat_flag_;
	TBranch *els_pat_flag_branch;
	bool els_pat_flag_isLoaded;
	vector<unsigned int>	jets_pat_flag_;
	TBranch *jets_pat_flag_branch;
	bool jets_pat_flag_isLoaded;
	vector<unsigned int>	mus_pat_flag_;
	TBranch *mus_pat_flag_branch;
	bool mus_pat_flag_isLoaded;
void Init(TTree *tree) {
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
	mus_mc_p4_branch = 0;
	if (tree->GetAlias("mus_mc_p4") != 0) {
		mus_mc_p4_branch = tree->GetBranch(tree->GetAlias("mus_mc_p4"));
		mus_mc_p4_branch->SetAddress(&mus_mc_p4_);
	}
	if(mus_mc_p4_branch == 0 ) {
	cout << "Branch mus_mc_p4 does not exist." << endl;
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
	trks_trk_p4_branch = 0;
	if (tree->GetAlias("trks_trk_p4") != 0) {
		trks_trk_p4_branch = tree->GetBranch(tree->GetAlias("trks_trk_p4"));
		trks_trk_p4_branch->SetAddress(&trks_trk_p4_);
	}
	if(trks_trk_p4_branch == 0 ) {
	cout << "Branch trks_trk_p4 does not exist." << endl;
	}
  tree->SetMakeClass(1);
	evt_bs_branch = 0;
	if (tree->GetAlias("evt_bs") != 0) {
		evt_bs_branch = tree->GetBranch(tree->GetAlias("evt_bs"));
		evt_bs_branch->SetAddress(&evt_bs_);
	}
	if(evt_bs_branch == 0 ) {
	cout << "Branch evt_bs does not exist." << endl;
	}
	mus_gfit_outerPos_branch = 0;
	if (tree->GetAlias("mus_gfit_outerPos") != 0) {
		mus_gfit_outerPos_branch = tree->GetBranch(tree->GetAlias("mus_gfit_outerPos"));
		mus_gfit_outerPos_branch->SetAddress(&mus_gfit_outerPos_);
	}
	if(mus_gfit_outerPos_branch == 0 ) {
	cout << "Branch mus_gfit_outerPos does not exist." << endl;
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
	hyp_jets_pat_p4_branch = 0;
	if (tree->GetAlias("hyp_jets_pat_p4") != 0) {
		hyp_jets_pat_p4_branch = tree->GetBranch(tree->GetAlias("hyp_jets_pat_p4"));
		hyp_jets_pat_p4_branch->SetAddress(&hyp_jets_pat_p4_);
	}
	if(hyp_jets_pat_p4_branch == 0 ) {
	cout << "Branch hyp_jets_pat_p4 does not exist." << endl;
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
	mus_mcdr_branch = 0;
	if (tree->GetAlias("mus_mcdr") != 0) {
		mus_mcdr_branch = tree->GetBranch(tree->GetAlias("mus_mcdr"));
		mus_mcdr_branch->SetAddress(&mus_mcdr_);
	}
	if(mus_mcdr_branch == 0 ) {
	cout << "Branch mus_mcdr does not exist." << endl;
	}
	trk_mcdr_branch = 0;
	if (tree->GetAlias("trk_mcdr") != 0) {
		trk_mcdr_branch = tree->GetBranch(tree->GetAlias("trk_mcdr"));
		trk_mcdr_branch->SetAddress(&trk_mcdr_);
	}
	if(trk_mcdr_branch == 0 ) {
	cout << "Branch trk_mcdr does not exist." << endl;
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
	jets_pat_gluCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_gluCorrF") != 0) {
		jets_pat_gluCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_gluCorrF"));
		jets_pat_gluCorrF_branch->SetAddress(&jets_pat_gluCorrF_);
	}
	if(jets_pat_gluCorrF_branch == 0 ) {
	cout << "Branch jets_pat_gluCorrF does not exist." << endl;
	}
	jets_pat_jetCharge_branch = 0;
	if (tree->GetAlias("jets_pat_jetCharge") != 0) {
		jets_pat_jetCharge_branch = tree->GetBranch(tree->GetAlias("jets_pat_jetCharge"));
		jets_pat_jetCharge_branch->SetAddress(&jets_pat_jetCharge_);
	}
	if(jets_pat_jetCharge_branch == 0 ) {
	cout << "Branch jets_pat_jetCharge does not exist." << endl;
	}
	jets_pat_noCorrF_branch = 0;
	if (tree->GetAlias("jets_pat_noCorrF") != 0) {
		jets_pat_noCorrF_branch = tree->GetBranch(tree->GetAlias("jets_pat_noCorrF"));
		jets_pat_noCorrF_branch->SetAddress(&jets_pat_noCorrF_);
	}
	if(jets_pat_noCorrF_branch == 0 ) {
	cout << "Branch jets_pat_noCorrF does not exist." << endl;
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
	mus_pat_leptonID_branch = 0;
	if (tree->GetAlias("mus_pat_leptonID") != 0) {
		mus_pat_leptonID_branch = tree->GetBranch(tree->GetAlias("mus_pat_leptonID"));
		mus_pat_leptonID_branch->SetAddress(&mus_pat_leptonID_);
	}
	if(mus_pat_leptonID_branch == 0 ) {
	cout << "Branch mus_pat_leptonID does not exist." << endl;
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
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index_ = idx;
		evt_bs_isLoaded = false;
		mus_gfit_outerPos_isLoaded = false;
		els_mc_p4_isLoaded = false;
		jets_mc_gp_p4_isLoaded = false;
		jets_mc_p4_isLoaded = false;
		mus_mc_p4_isLoaded = false;
		trk_mcp4_isLoaded = false;
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
		mus_p4_isLoaded = false;
		mus_trk_p4_isLoaded = false;
		els_pat_genMotherP4_isLoaded = false;
		els_pat_genP4_isLoaded = false;
		jets_pat_genJet_p4_isLoaded = false;
		jets_pat_genPartonMother_p4_isLoaded = false;
		jets_pat_genParton_p4_isLoaded = false;
		jets_pat_jet_p4_isLoaded = false;
		mus_pat_genMotherP4_isLoaded = false;
		mus_pat_genP4_isLoaded = false;
		trks_trk_p4_isLoaded = false;
		hyp_jets_mc_gp_p4_isLoaded = false;
		hyp_jets_mc_p4_isLoaded = false;
		hyp_jets_p4_isLoaded = false;
		hyp_jets_pat_p4_isLoaded = false;
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
		evt_met_isLoaded = false;
		evt_metHO_isLoaded = false;
		evt_metHOPhi_isLoaded = false;
		evt_metHOSig_isLoaded = false;
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
		met_pat_metCor_isLoaded = false;
		met_pat_metPhiCor_isLoaded = false;
		met_pat_metPhiUncor_isLoaded = false;
		met_pat_metPhiUncorJES_isLoaded = false;
		met_pat_metPhiUncorMuon_isLoaded = false;
		met_pat_metUncor_isLoaded = false;
		met_pat_metUncorJES_isLoaded = false;
		met_pat_metUncorMuon_isLoaded = false;
		els_mcdr_isLoaded = false;
		jets_mcdr_isLoaded = false;
		jets_mc_emEnergy_isLoaded = false;
		jets_mc_gpdr_isLoaded = false;
		jets_mc_hadEnergy_isLoaded = false;
		jets_mc_invEnergy_isLoaded = false;
		jets_mc_otherEnergy_isLoaded = false;
		mus_mcdr_isLoaded = false;
		trk_mcdr_isLoaded = false;
		els_musdr_isLoaded = false;
		els_trkdr_isLoaded = false;
		els_trkshFrac_isLoaded = false;
		els_ESc_isLoaded = false;
		els_ESc_raw_isLoaded = false;
		els_ESeed_isLoaded = false;
		els_chi2_isLoaded = false;
		els_d0_isLoaded = false;
		els_d0Err_isLoaded = false;
		els_d0corr_isLoaded = false;
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
		mus_eledr_isLoaded = false;
		mus_jetdr_isLoaded = false;
		mus_trkdr_isLoaded = false;
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
		els_pat_robustLooseId_isLoaded = false;
		els_pat_robustTightId_isLoaded = false;
		els_pat_tightId_isLoaded = false;
		els_pat_trackIso_isLoaded = false;
		jets_pat_bCorrF_isLoaded = false;
		jets_pat_cCorrF_isLoaded = false;
		jets_pat_gluCorrF_isLoaded = false;
		jets_pat_jetCharge_isLoaded = false;
		jets_pat_noCorrF_isLoaded = false;
		jets_pat_udsCorrF_isLoaded = false;
		mus_pat_caloIso_isLoaded = false;
		mus_pat_ecalIso_isLoaded = false;
		mus_pat_ecalvetoDep_isLoaded = false;
		mus_pat_hcalIso_isLoaded = false;
		mus_pat_hcalvetoDep_isLoaded = false;
		mus_pat_leptonID_isLoaded = false;
		mus_pat_trackIso_isLoaded = false;
		mus_pat_vetoDep_isLoaded = false;
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
		trks_vertexphi_isLoaded = false;
		trks_z0_isLoaded = false;
		trks_z0Err_isLoaded = false;
		trks_z0corr_isLoaded = false;
		trks_elsdr_isLoaded = false;
		trks_elsshFrac_isLoaded = false;
		trk_musdr_isLoaded = false;
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
		els_mc_id_isLoaded = false;
		els_mcidx_isLoaded = false;
		els_mc_motherid_isLoaded = false;
		jets_mc_id_isLoaded = false;
		mus_mc_id_isLoaded = false;
		mus_mcidx_isLoaded = false;
		mus_mc_motherid_isLoaded = false;
		trk_mc_id_isLoaded = false;
		trk_mcidx_isLoaded = false;
		trk_mc_motherid_isLoaded = false;
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
		els_pat_flag_isLoaded = false;
		jets_pat_flag_isLoaded = false;
		mus_pat_flag_isLoaded = false;
	}

	ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag>  &evt_bs()
	{
		if (not evt_bs_isLoaded) {
			if (evt_bs_branch != 0) 
				evt_bs_branch->GetEntry(index_);
			evt_bs_isLoaded = true;
		}
		return evt_bs_;
	}
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &mus_gfit_outerPos()
	{
		if (not mus_gfit_outerPos_isLoaded) {
			if (mus_gfit_outerPos_branch != 0) 
				mus_gfit_outerPos_branch->GetEntry(index_);
			mus_gfit_outerPos_isLoaded = true;
		}
		return mus_gfit_outerPos_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_p4()
	{
		if (not els_mc_p4_isLoaded) {
			if (els_mc_p4_branch != 0) 
				els_mc_p4_branch->GetEntry(index_);
			els_mc_p4_isLoaded = true;
		}
		return els_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_gp_p4()
	{
		if (not jets_mc_gp_p4_isLoaded) {
			if (jets_mc_gp_p4_branch != 0) 
				jets_mc_gp_p4_branch->GetEntry(index_);
			jets_mc_gp_p4_isLoaded = true;
		}
		return jets_mc_gp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_p4()
	{
		if (not jets_mc_p4_isLoaded) {
			if (jets_mc_p4_branch != 0) 
				jets_mc_p4_branch->GetEntry(index_);
			jets_mc_p4_isLoaded = true;
		}
		return jets_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc_p4()
	{
		if (not mus_mc_p4_isLoaded) {
			if (mus_mc_p4_branch != 0) 
				mus_mc_p4_branch->GetEntry(index_);
			mus_mc_p4_isLoaded = true;
		}
		return mus_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trk_mcp4()
	{
		if (not trk_mcp4_isLoaded) {
			if (trk_mcp4_branch != 0) 
				trk_mcp4_branch->GetEntry(index_);
			trk_mcp4_isLoaded = true;
		}
		return trk_mcp4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4()
	{
		if (not els_p4_isLoaded) {
			if (els_p4_branch != 0) 
				els_p4_branch->GetEntry(index_);
			els_p4_isLoaded = true;
		}
		return els_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4In()
	{
		if (not els_p4In_isLoaded) {
			if (els_p4In_branch != 0) 
				els_p4In_branch->GetEntry(index_);
			els_p4In_isLoaded = true;
		}
		return els_p4In_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_p4Out()
	{
		if (not els_p4Out_isLoaded) {
			if (els_p4Out_branch != 0) 
				els_p4Out_branch->GetEntry(index_);
			els_p4Out_isLoaded = true;
		}
		return els_p4Out_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_trk_p4()
	{
		if (not els_trk_p4_isLoaded) {
			if (els_trk_p4_branch != 0) 
				els_trk_p4_branch->GetEntry(index_);
			els_trk_p4_isLoaded = true;
		}
		return els_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_p4()
	{
		if (not genps_p4_isLoaded) {
			if (genps_p4_branch != 0) 
				genps_p4_branch->GetEntry(index_);
			genps_p4_isLoaded = true;
		}
		return genps_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &genps_prod_vtx()
	{
		if (not genps_prod_vtx_isLoaded) {
			if (genps_prod_vtx_branch != 0) 
				genps_prod_vtx_branch->GetEntry(index_);
			genps_prod_vtx_isLoaded = true;
		}
		return genps_prod_vtx_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_mc_p4()
	{
		if (not hyp_ll_mc_p4_isLoaded) {
			if (hyp_ll_mc_p4_branch != 0) 
				hyp_ll_mc_p4_branch->GetEntry(index_);
			hyp_ll_mc_p4_isLoaded = true;
		}
		return hyp_ll_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_p4()
	{
		if (not hyp_ll_p4_isLoaded) {
			if (hyp_ll_p4_branch != 0) 
				hyp_ll_p4_branch->GetEntry(index_);
			hyp_ll_p4_isLoaded = true;
		}
		return hyp_ll_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_ll_trk_p4()
	{
		if (not hyp_ll_trk_p4_isLoaded) {
			if (hyp_ll_trk_p4_branch != 0) 
				hyp_ll_trk_p4_branch->GetEntry(index_);
			hyp_ll_trk_p4_isLoaded = true;
		}
		return hyp_ll_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_mc_p4()
	{
		if (not hyp_lt_mc_p4_isLoaded) {
			if (hyp_lt_mc_p4_branch != 0) 
				hyp_lt_mc_p4_branch->GetEntry(index_);
			hyp_lt_mc_p4_isLoaded = true;
		}
		return hyp_lt_mc_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_p4()
	{
		if (not hyp_lt_p4_isLoaded) {
			if (hyp_lt_p4_branch != 0) 
				hyp_lt_p4_branch->GetEntry(index_);
			hyp_lt_p4_isLoaded = true;
		}
		return hyp_lt_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_lt_trk_p4()
	{
		if (not hyp_lt_trk_p4_isLoaded) {
			if (hyp_lt_trk_p4_branch != 0) 
				hyp_lt_trk_p4_branch->GetEntry(index_);
			hyp_lt_trk_p4_isLoaded = true;
		}
		return hyp_lt_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &hyp_p4()
	{
		if (not hyp_p4_isLoaded) {
			if (hyp_p4_branch != 0) 
				hyp_p4_branch->GetEntry(index_);
			hyp_p4_isLoaded = true;
		}
		return hyp_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_p4()
	{
		if (not jets_p4_isLoaded) {
			if (jets_p4_branch != 0) 
				jets_p4_branch->GetEntry(index_);
			jets_p4_isLoaded = true;
		}
		return jets_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4()
	{
		if (not mus_p4_isLoaded) {
			if (mus_p4_branch != 0) 
				mus_p4_branch->GetEntry(index_);
			mus_p4_isLoaded = true;
		}
		return mus_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4()
	{
		if (not mus_trk_p4_isLoaded) {
			if (mus_trk_p4_branch != 0) 
				mus_trk_p4_branch->GetEntry(index_);
			mus_trk_p4_isLoaded = true;
		}
		return mus_trk_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genMotherP4()
	{
		if (not els_pat_genMotherP4_isLoaded) {
			if (els_pat_genMotherP4_branch != 0) 
				els_pat_genMotherP4_branch->GetEntry(index_);
			els_pat_genMotherP4_isLoaded = true;
		}
		return els_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genP4()
	{
		if (not els_pat_genP4_isLoaded) {
			if (els_pat_genP4_branch != 0) 
				els_pat_genP4_branch->GetEntry(index_);
			els_pat_genP4_isLoaded = true;
		}
		return els_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genJet_p4()
	{
		if (not jets_pat_genJet_p4_isLoaded) {
			if (jets_pat_genJet_p4_branch != 0) 
				jets_pat_genJet_p4_branch->GetEntry(index_);
			jets_pat_genJet_p4_isLoaded = true;
		}
		return jets_pat_genJet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genPartonMother_p4()
	{
		if (not jets_pat_genPartonMother_p4_isLoaded) {
			if (jets_pat_genPartonMother_p4_branch != 0) 
				jets_pat_genPartonMother_p4_branch->GetEntry(index_);
			jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return jets_pat_genPartonMother_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genParton_p4()
	{
		if (not jets_pat_genParton_p4_isLoaded) {
			if (jets_pat_genParton_p4_branch != 0) 
				jets_pat_genParton_p4_branch->GetEntry(index_);
			jets_pat_genParton_p4_isLoaded = true;
		}
		return jets_pat_genParton_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_p4()
	{
		if (not jets_pat_jet_p4_isLoaded) {
			if (jets_pat_jet_p4_branch != 0) 
				jets_pat_jet_p4_branch->GetEntry(index_);
			jets_pat_jet_p4_isLoaded = true;
		}
		return jets_pat_jet_p4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genMotherP4()
	{
		if (not mus_pat_genMotherP4_isLoaded) {
			if (mus_pat_genMotherP4_branch != 0) 
				mus_pat_genMotherP4_branch->GetEntry(index_);
			mus_pat_genMotherP4_isLoaded = true;
		}
		return mus_pat_genMotherP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genP4()
	{
		if (not mus_pat_genP4_isLoaded) {
			if (mus_pat_genP4_branch != 0) 
				mus_pat_genP4_branch->GetEntry(index_);
			mus_pat_genP4_isLoaded = true;
		}
		return mus_pat_genP4_;
	}
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &trks_trk_p4()
	{
		if (not trks_trk_p4_isLoaded) {
			if (trks_trk_p4_branch != 0) 
				trks_trk_p4_branch->GetEntry(index_);
			trks_trk_p4_isLoaded = true;
		}
		return trks_trk_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_gp_p4()
	{
		if (not hyp_jets_mc_gp_p4_isLoaded) {
			if (hyp_jets_mc_gp_p4_branch != 0) 
				hyp_jets_mc_gp_p4_branch->GetEntry(index_);
			hyp_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_p4()
	{
		if (not hyp_jets_mc_p4_isLoaded) {
			if (hyp_jets_mc_p4_branch != 0) 
				hyp_jets_mc_p4_branch->GetEntry(index_);
			hyp_jets_mc_p4_isLoaded = true;
		}
		return hyp_jets_mc_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_p4()
	{
		if (not hyp_jets_p4_isLoaded) {
			if (hyp_jets_p4_branch != 0) 
				hyp_jets_p4_branch->GetEntry(index_);
			hyp_jets_p4_isLoaded = true;
		}
		return hyp_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_p4()
	{
		if (not hyp_jets_pat_p4_isLoaded) {
			if (hyp_jets_pat_p4_branch != 0) 
				hyp_jets_pat_p4_branch->GetEntry(index_);
			hyp_jets_pat_p4_isLoaded = true;
		}
		return hyp_jets_pat_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genPartonMother_p4()
	{
		if (not hyp_jets_pat_genPartonMother_p4_isLoaded) {
			if (hyp_jets_pat_genPartonMother_p4_branch != 0) 
				hyp_jets_pat_genPartonMother_p4_branch->GetEntry(index_);
			hyp_jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return hyp_jets_pat_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genParton_p4()
	{
		if (not hyp_jets_pat_genParton_p4_isLoaded) {
			if (hyp_jets_pat_genParton_p4_branch != 0) 
				hyp_jets_pat_genParton_p4_branch->GetEntry(index_);
			hyp_jets_pat_genParton_p4_isLoaded = true;
		}
		return hyp_jets_pat_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_jet_p4()
	{
		if (not hyp_jets_pat_jet_p4_isLoaded) {
			if (hyp_jets_pat_jet_p4_branch != 0) 
				hyp_jets_pat_jet_p4_branch->GetEntry(index_);
			hyp_jets_pat_jet_p4_isLoaded = true;
		}
		return hyp_jets_pat_jet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_gp_p4()
	{
		if (not hyp_other_jets_mc_gp_p4_isLoaded) {
			if (hyp_other_jets_mc_gp_p4_branch != 0) 
				hyp_other_jets_mc_gp_p4_branch->GetEntry(index_);
			hyp_other_jets_mc_gp_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_gp_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_p4()
	{
		if (not hyp_other_jets_mc_p4_isLoaded) {
			if (hyp_other_jets_mc_p4_branch != 0) 
				hyp_other_jets_mc_p4_branch->GetEntry(index_);
			hyp_other_jets_mc_p4_isLoaded = true;
		}
		return hyp_other_jets_mc_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_p4()
	{
		if (not hyp_other_jets_p4_isLoaded) {
			if (hyp_other_jets_p4_branch != 0) 
				hyp_other_jets_p4_branch->GetEntry(index_);
			hyp_other_jets_p4_isLoaded = true;
		}
		return hyp_other_jets_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genJet_p4()
	{
		if (not hyp_other_jets_pat_genJet_p4_isLoaded) {
			if (hyp_other_jets_pat_genJet_p4_branch != 0) 
				hyp_other_jets_pat_genJet_p4_branch->GetEntry(index_);
			hyp_other_jets_pat_genJet_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genJet_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genPartonMother_p4()
	{
		if (not hyp_other_jets_pat_genPartonMother_p4_isLoaded) {
			if (hyp_other_jets_pat_genPartonMother_p4_branch != 0) 
				hyp_other_jets_pat_genPartonMother_p4_branch->GetEntry(index_);
			hyp_other_jets_pat_genPartonMother_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genPartonMother_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genParton_p4()
	{
		if (not hyp_other_jets_pat_genParton_p4_isLoaded) {
			if (hyp_other_jets_pat_genParton_p4_branch != 0) 
				hyp_other_jets_pat_genParton_p4_branch->GetEntry(index_);
			hyp_other_jets_pat_genParton_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_genParton_p4_;
	}
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_jet_p4()
	{
		if (not hyp_other_jets_pat_jet_p4_isLoaded) {
			if (hyp_other_jets_pat_jet_p4_branch != 0) 
				hyp_other_jets_pat_jet_p4_branch->GetEntry(index_);
			hyp_other_jets_pat_jet_p4_isLoaded = true;
		}
		return hyp_other_jets_pat_jet_p4_;
	}
	vector<double> &jets_closestElectron_DR()
	{
		if (not jets_closestElectron_DR_isLoaded) {
			if (jets_closestElectron_DR_branch != 0) 
				jets_closestElectron_DR_branch->GetEntry(index_);
			jets_closestElectron_DR_isLoaded = true;
		}
		return jets_closestElectron_DR_;
	}
	vector<double> &jets_closestMuon_DR()
	{
		if (not jets_closestMuon_DR_isLoaded) {
			if (jets_closestMuon_DR_branch != 0) 
				jets_closestMuon_DR_branch->GetEntry(index_);
			jets_closestMuon_DR_isLoaded = true;
		}
		return jets_closestMuon_DR_;
	}
	float &evt_bs_dxdz()
	{
		if (not evt_bs_dxdz_isLoaded) {
			if (evt_bs_dxdz_branch != 0) 
				evt_bs_dxdz_branch->GetEntry(index_);
			evt_bs_dxdz_isLoaded = true;
		}
		return evt_bs_dxdz_;
	}
	float &evt_bs_dxdzErr()
	{
		if (not evt_bs_dxdzErr_isLoaded) {
			if (evt_bs_dxdzErr_branch != 0) 
				evt_bs_dxdzErr_branch->GetEntry(index_);
			evt_bs_dxdzErr_isLoaded = true;
		}
		return evt_bs_dxdzErr_;
	}
	float &evt_bs_dydz()
	{
		if (not evt_bs_dydz_isLoaded) {
			if (evt_bs_dydz_branch != 0) 
				evt_bs_dydz_branch->GetEntry(index_);
			evt_bs_dydz_isLoaded = true;
		}
		return evt_bs_dydz_;
	}
	float &evt_bs_dydzErr()
	{
		if (not evt_bs_dydzErr_isLoaded) {
			if (evt_bs_dydzErr_branch != 0) 
				evt_bs_dydzErr_branch->GetEntry(index_);
			evt_bs_dydzErr_isLoaded = true;
		}
		return evt_bs_dydzErr_;
	}
	float &evt_bs_sigmaZ()
	{
		if (not evt_bs_sigmaZ_isLoaded) {
			if (evt_bs_sigmaZ_branch != 0) 
				evt_bs_sigmaZ_branch->GetEntry(index_);
			evt_bs_sigmaZ_isLoaded = true;
		}
		return evt_bs_sigmaZ_;
	}
	float &evt_bs_sigmaZErr()
	{
		if (not evt_bs_sigmaZErr_isLoaded) {
			if (evt_bs_sigmaZErr_branch != 0) 
				evt_bs_sigmaZErr_branch->GetEntry(index_);
			evt_bs_sigmaZErr_isLoaded = true;
		}
		return evt_bs_sigmaZErr_;
	}
	float &evt_bs_widthErr()
	{
		if (not evt_bs_widthErr_isLoaded) {
			if (evt_bs_widthErr_branch != 0) 
				evt_bs_widthErr_branch->GetEntry(index_);
			evt_bs_widthErr_isLoaded = true;
		}
		return evt_bs_widthErr_;
	}
	float &evt_bs_xErr()
	{
		if (not evt_bs_xErr_isLoaded) {
			if (evt_bs_xErr_branch != 0) 
				evt_bs_xErr_branch->GetEntry(index_);
			evt_bs_xErr_isLoaded = true;
		}
		return evt_bs_xErr_;
	}
	float &evt_bs_yErr()
	{
		if (not evt_bs_yErr_isLoaded) {
			if (evt_bs_yErr_branch != 0) 
				evt_bs_yErr_branch->GetEntry(index_);
			evt_bs_yErr_isLoaded = true;
		}
		return evt_bs_yErr_;
	}
	float &evt_bs_zErr()
	{
		if (not evt_bs_zErr_isLoaded) {
			if (evt_bs_zErr_branch != 0) 
				evt_bs_zErr_branch->GetEntry(index_);
			evt_bs_zErr_isLoaded = true;
		}
		return evt_bs_zErr_;
	}
	float &gen_met()
	{
		if (not gen_met_isLoaded) {
			if (gen_met_branch != 0) 
				gen_met_branch->GetEntry(index_);
			gen_met_isLoaded = true;
		}
		return gen_met_;
	}
	float &gen_metPhi()
	{
		if (not gen_metPhi_isLoaded) {
			if (gen_metPhi_branch != 0) 
				gen_metPhi_branch->GetEntry(index_);
			gen_metPhi_isLoaded = true;
		}
		return gen_metPhi_;
	}
	float &evt_bField()
	{
		if (not evt_bField_isLoaded) {
			if (evt_bField_branch != 0) 
				evt_bField_branch->GetEntry(index_);
			evt_bField_isLoaded = true;
		}
		return evt_bField_;
	}
	float &evt_kfactor()
	{
		if (not evt_kfactor_isLoaded) {
			if (evt_kfactor_branch != 0) 
				evt_kfactor_branch->GetEntry(index_);
			evt_kfactor_isLoaded = true;
		}
		return evt_kfactor_;
	}
	float &evt_weight()
	{
		if (not evt_weight_isLoaded) {
			if (evt_weight_branch != 0) 
				evt_weight_branch->GetEntry(index_);
			evt_weight_isLoaded = true;
		}
		return evt_weight_;
	}
	float &evt_xsec_excl()
	{
		if (not evt_xsec_excl_isLoaded) {
			if (evt_xsec_excl_branch != 0) 
				evt_xsec_excl_branch->GetEntry(index_);
			evt_xsec_excl_isLoaded = true;
		}
		return evt_xsec_excl_;
	}
	float &evt_xsec_incl()
	{
		if (not evt_xsec_incl_isLoaded) {
			if (evt_xsec_incl_branch != 0) 
				evt_xsec_incl_branch->GetEntry(index_);
			evt_xsec_incl_isLoaded = true;
		}
		return evt_xsec_incl_;
	}
	float &evt_met()
	{
		if (not evt_met_isLoaded) {
			if (evt_met_branch != 0) 
				evt_met_branch->GetEntry(index_);
			evt_met_isLoaded = true;
		}
		return evt_met_;
	}
	float &evt_metHO()
	{
		if (not evt_metHO_isLoaded) {
			if (evt_metHO_branch != 0) 
				evt_metHO_branch->GetEntry(index_);
			evt_metHO_isLoaded = true;
		}
		return evt_metHO_;
	}
	float &evt_metHOPhi()
	{
		if (not evt_metHOPhi_isLoaded) {
			if (evt_metHOPhi_branch != 0) 
				evt_metHOPhi_branch->GetEntry(index_);
			evt_metHOPhi_isLoaded = true;
		}
		return evt_metHOPhi_;
	}
	float &evt_metHOSig()
	{
		if (not evt_metHOSig_isLoaded) {
			if (evt_metHOSig_branch != 0) 
				evt_metHOSig_branch->GetEntry(index_);
			evt_metHOSig_isLoaded = true;
		}
		return evt_metHOSig_;
	}
	float &evt_metNoHF()
	{
		if (not evt_metNoHF_isLoaded) {
			if (evt_metNoHF_branch != 0) 
				evt_metNoHF_branch->GetEntry(index_);
			evt_metNoHF_isLoaded = true;
		}
		return evt_metNoHF_;
	}
	float &evt_metNoHFHO()
	{
		if (not evt_metNoHFHO_isLoaded) {
			if (evt_metNoHFHO_branch != 0) 
				evt_metNoHFHO_branch->GetEntry(index_);
			evt_metNoHFHO_isLoaded = true;
		}
		return evt_metNoHFHO_;
	}
	float &evt_metNoHFHOPhi()
	{
		if (not evt_metNoHFHOPhi_isLoaded) {
			if (evt_metNoHFHOPhi_branch != 0) 
				evt_metNoHFHOPhi_branch->GetEntry(index_);
			evt_metNoHFHOPhi_isLoaded = true;
		}
		return evt_metNoHFHOPhi_;
	}
	float &evt_metNoHFHOSig()
	{
		if (not evt_metNoHFHOSig_isLoaded) {
			if (evt_metNoHFHOSig_branch != 0) 
				evt_metNoHFHOSig_branch->GetEntry(index_);
			evt_metNoHFHOSig_isLoaded = true;
		}
		return evt_metNoHFHOSig_;
	}
	float &evt_metNoHFPhi()
	{
		if (not evt_metNoHFPhi_isLoaded) {
			if (evt_metNoHFPhi_branch != 0) 
				evt_metNoHFPhi_branch->GetEntry(index_);
			evt_metNoHFPhi_isLoaded = true;
		}
		return evt_metNoHFPhi_;
	}
	float &evt_metSig()
	{
		if (not evt_metSig_isLoaded) {
			if (evt_metSig_branch != 0) 
				evt_metSig_branch->GetEntry(index_);
			evt_metSig_isLoaded = true;
		}
		return evt_metSig_;
	}
	float &evt_metOpt()
	{
		if (not evt_metOpt_isLoaded) {
			if (evt_metOpt_branch != 0) 
				evt_metOpt_branch->GetEntry(index_);
			evt_metOpt_isLoaded = true;
		}
		return evt_metOpt_;
	}
	float &evt_metOptHO()
	{
		if (not evt_metOptHO_isLoaded) {
			if (evt_metOptHO_branch != 0) 
				evt_metOptHO_branch->GetEntry(index_);
			evt_metOptHO_isLoaded = true;
		}
		return evt_metOptHO_;
	}
	float &evt_metOptHOPhi()
	{
		if (not evt_metOptHOPhi_isLoaded) {
			if (evt_metOptHOPhi_branch != 0) 
				evt_metOptHOPhi_branch->GetEntry(index_);
			evt_metOptHOPhi_isLoaded = true;
		}
		return evt_metOptHOPhi_;
	}
	float &evt_metOptHOSig()
	{
		if (not evt_metOptHOSig_isLoaded) {
			if (evt_metOptHOSig_branch != 0) 
				evt_metOptHOSig_branch->GetEntry(index_);
			evt_metOptHOSig_isLoaded = true;
		}
		return evt_metOptHOSig_;
	}
	float &evt_metOptNoHF()
	{
		if (not evt_metOptNoHF_isLoaded) {
			if (evt_metOptNoHF_branch != 0) 
				evt_metOptNoHF_branch->GetEntry(index_);
			evt_metOptNoHF_isLoaded = true;
		}
		return evt_metOptNoHF_;
	}
	float &evt_metOptNoHFHO()
	{
		if (not evt_metOptNoHFHO_isLoaded) {
			if (evt_metOptNoHFHO_branch != 0) 
				evt_metOptNoHFHO_branch->GetEntry(index_);
			evt_metOptNoHFHO_isLoaded = true;
		}
		return evt_metOptNoHFHO_;
	}
	float &evt_metOptNoHFHOPhi()
	{
		if (not evt_metOptNoHFHOPhi_isLoaded) {
			if (evt_metOptNoHFHOPhi_branch != 0) 
				evt_metOptNoHFHOPhi_branch->GetEntry(index_);
			evt_metOptNoHFHOPhi_isLoaded = true;
		}
		return evt_metOptNoHFHOPhi_;
	}
	float &evt_metOptNoHFHOSig()
	{
		if (not evt_metOptNoHFHOSig_isLoaded) {
			if (evt_metOptNoHFHOSig_branch != 0) 
				evt_metOptNoHFHOSig_branch->GetEntry(index_);
			evt_metOptNoHFHOSig_isLoaded = true;
		}
		return evt_metOptNoHFHOSig_;
	}
	float &evt_metOptNoHFPhi()
	{
		if (not evt_metOptNoHFPhi_isLoaded) {
			if (evt_metOptNoHFPhi_branch != 0) 
				evt_metOptNoHFPhi_branch->GetEntry(index_);
			evt_metOptNoHFPhi_isLoaded = true;
		}
		return evt_metOptNoHFPhi_;
	}
	float &evt_metOptSig()
	{
		if (not evt_metOptSig_isLoaded) {
			if (evt_metOptSig_branch != 0) 
				evt_metOptSig_branch->GetEntry(index_);
			evt_metOptSig_isLoaded = true;
		}
		return evt_metOptSig_;
	}
	float &evt_metOptPhi()
	{
		if (not evt_metOptPhi_isLoaded) {
			if (evt_metOptPhi_branch != 0) 
				evt_metOptPhi_branch->GetEntry(index_);
			evt_metOptPhi_isLoaded = true;
		}
		return evt_metOptPhi_;
	}
	float &evt_metPhi()
	{
		if (not evt_metPhi_isLoaded) {
			if (evt_metPhi_branch != 0) 
				evt_metPhi_branch->GetEntry(index_);
			evt_metPhi_isLoaded = true;
		}
		return evt_metPhi_;
	}
	float &met_pat_metCor()
	{
		if (not met_pat_metCor_isLoaded) {
			if (met_pat_metCor_branch != 0) 
				met_pat_metCor_branch->GetEntry(index_);
			met_pat_metCor_isLoaded = true;
		}
		return met_pat_metCor_;
	}
	float &met_pat_metPhiCor()
	{
		if (not met_pat_metPhiCor_isLoaded) {
			if (met_pat_metPhiCor_branch != 0) 
				met_pat_metPhiCor_branch->GetEntry(index_);
			met_pat_metPhiCor_isLoaded = true;
		}
		return met_pat_metPhiCor_;
	}
	float &met_pat_metPhiUncor()
	{
		if (not met_pat_metPhiUncor_isLoaded) {
			if (met_pat_metPhiUncor_branch != 0) 
				met_pat_metPhiUncor_branch->GetEntry(index_);
			met_pat_metPhiUncor_isLoaded = true;
		}
		return met_pat_metPhiUncor_;
	}
	float &met_pat_metPhiUncorJES()
	{
		if (not met_pat_metPhiUncorJES_isLoaded) {
			if (met_pat_metPhiUncorJES_branch != 0) 
				met_pat_metPhiUncorJES_branch->GetEntry(index_);
			met_pat_metPhiUncorJES_isLoaded = true;
		}
		return met_pat_metPhiUncorJES_;
	}
	float &met_pat_metPhiUncorMuon()
	{
		if (not met_pat_metPhiUncorMuon_isLoaded) {
			if (met_pat_metPhiUncorMuon_branch != 0) 
				met_pat_metPhiUncorMuon_branch->GetEntry(index_);
			met_pat_metPhiUncorMuon_isLoaded = true;
		}
		return met_pat_metPhiUncorMuon_;
	}
	float &met_pat_metUncor()
	{
		if (not met_pat_metUncor_isLoaded) {
			if (met_pat_metUncor_branch != 0) 
				met_pat_metUncor_branch->GetEntry(index_);
			met_pat_metUncor_isLoaded = true;
		}
		return met_pat_metUncor_;
	}
	float &met_pat_metUncorJES()
	{
		if (not met_pat_metUncorJES_isLoaded) {
			if (met_pat_metUncorJES_branch != 0) 
				met_pat_metUncorJES_branch->GetEntry(index_);
			met_pat_metUncorJES_isLoaded = true;
		}
		return met_pat_metUncorJES_;
	}
	float &met_pat_metUncorMuon()
	{
		if (not met_pat_metUncorMuon_isLoaded) {
			if (met_pat_metUncorMuon_branch != 0) 
				met_pat_metUncorMuon_branch->GetEntry(index_);
			met_pat_metUncorMuon_isLoaded = true;
		}
		return met_pat_metUncorMuon_;
	}
	vector<float> &els_mcdr()
	{
		if (not els_mcdr_isLoaded) {
			if (els_mcdr_branch != 0) 
				els_mcdr_branch->GetEntry(index_);
			els_mcdr_isLoaded = true;
		}
		return els_mcdr_;
	}
	vector<float> &jets_mcdr()
	{
		if (not jets_mcdr_isLoaded) {
			if (jets_mcdr_branch != 0) 
				jets_mcdr_branch->GetEntry(index_);
			jets_mcdr_isLoaded = true;
		}
		return jets_mcdr_;
	}
	vector<float> &jets_mc_emEnergy()
	{
		if (not jets_mc_emEnergy_isLoaded) {
			if (jets_mc_emEnergy_branch != 0) 
				jets_mc_emEnergy_branch->GetEntry(index_);
			jets_mc_emEnergy_isLoaded = true;
		}
		return jets_mc_emEnergy_;
	}
	vector<float> &jets_mc_gpdr()
	{
		if (not jets_mc_gpdr_isLoaded) {
			if (jets_mc_gpdr_branch != 0) 
				jets_mc_gpdr_branch->GetEntry(index_);
			jets_mc_gpdr_isLoaded = true;
		}
		return jets_mc_gpdr_;
	}
	vector<float> &jets_mc_hadEnergy()
	{
		if (not jets_mc_hadEnergy_isLoaded) {
			if (jets_mc_hadEnergy_branch != 0) 
				jets_mc_hadEnergy_branch->GetEntry(index_);
			jets_mc_hadEnergy_isLoaded = true;
		}
		return jets_mc_hadEnergy_;
	}
	vector<float> &jets_mc_invEnergy()
	{
		if (not jets_mc_invEnergy_isLoaded) {
			if (jets_mc_invEnergy_branch != 0) 
				jets_mc_invEnergy_branch->GetEntry(index_);
			jets_mc_invEnergy_isLoaded = true;
		}
		return jets_mc_invEnergy_;
	}
	vector<float> &jets_mc_otherEnergy()
	{
		if (not jets_mc_otherEnergy_isLoaded) {
			if (jets_mc_otherEnergy_branch != 0) 
				jets_mc_otherEnergy_branch->GetEntry(index_);
			jets_mc_otherEnergy_isLoaded = true;
		}
		return jets_mc_otherEnergy_;
	}
	vector<float> &mus_mcdr()
	{
		if (not mus_mcdr_isLoaded) {
			if (mus_mcdr_branch != 0) 
				mus_mcdr_branch->GetEntry(index_);
			mus_mcdr_isLoaded = true;
		}
		return mus_mcdr_;
	}
	vector<float> &trk_mcdr()
	{
		if (not trk_mcdr_isLoaded) {
			if (trk_mcdr_branch != 0) 
				trk_mcdr_branch->GetEntry(index_);
			trk_mcdr_isLoaded = true;
		}
		return trk_mcdr_;
	}
	vector<float> &els_musdr()
	{
		if (not els_musdr_isLoaded) {
			if (els_musdr_branch != 0) 
				els_musdr_branch->GetEntry(index_);
			els_musdr_isLoaded = true;
		}
		return els_musdr_;
	}
	vector<float> &els_trkdr()
	{
		if (not els_trkdr_isLoaded) {
			if (els_trkdr_branch != 0) 
				els_trkdr_branch->GetEntry(index_);
			els_trkdr_isLoaded = true;
		}
		return els_trkdr_;
	}
	vector<float> &els_trkshFrac()
	{
		if (not els_trkshFrac_isLoaded) {
			if (els_trkshFrac_branch != 0) 
				els_trkshFrac_branch->GetEntry(index_);
			els_trkshFrac_isLoaded = true;
		}
		return els_trkshFrac_;
	}
	vector<float> &els_ESc()
	{
		if (not els_ESc_isLoaded) {
			if (els_ESc_branch != 0) 
				els_ESc_branch->GetEntry(index_);
			els_ESc_isLoaded = true;
		}
		return els_ESc_;
	}
	vector<float> &els_ESc_raw()
	{
		if (not els_ESc_raw_isLoaded) {
			if (els_ESc_raw_branch != 0) 
				els_ESc_raw_branch->GetEntry(index_);
			els_ESc_raw_isLoaded = true;
		}
		return els_ESc_raw_;
	}
	vector<float> &els_ESeed()
	{
		if (not els_ESeed_isLoaded) {
			if (els_ESeed_branch != 0) 
				els_ESeed_branch->GetEntry(index_);
			els_ESeed_isLoaded = true;
		}
		return els_ESeed_;
	}
	vector<float> &els_chi2()
	{
		if (not els_chi2_isLoaded) {
			if (els_chi2_branch != 0) 
				els_chi2_branch->GetEntry(index_);
			els_chi2_isLoaded = true;
		}
		return els_chi2_;
	}
	vector<float> &els_d0()
	{
		if (not els_d0_isLoaded) {
			if (els_d0_branch != 0) 
				els_d0_branch->GetEntry(index_);
			els_d0_isLoaded = true;
		}
		return els_d0_;
	}
	vector<float> &els_d0Err()
	{
		if (not els_d0Err_isLoaded) {
			if (els_d0Err_branch != 0) 
				els_d0Err_branch->GetEntry(index_);
			els_d0Err_isLoaded = true;
		}
		return els_d0Err_;
	}
	vector<float> &els_d0corr()
	{
		if (not els_d0corr_isLoaded) {
			if (els_d0corr_branch != 0) 
				els_d0corr_branch->GetEntry(index_);
			els_d0corr_isLoaded = true;
		}
		return els_d0corr_;
	}
	vector<float> &els_dEtaIn()
	{
		if (not els_dEtaIn_isLoaded) {
			if (els_dEtaIn_branch != 0) 
				els_dEtaIn_branch->GetEntry(index_);
			els_dEtaIn_isLoaded = true;
		}
		return els_dEtaIn_;
	}
	vector<float> &els_dEtaOut()
	{
		if (not els_dEtaOut_isLoaded) {
			if (els_dEtaOut_branch != 0) 
				els_dEtaOut_branch->GetEntry(index_);
			els_dEtaOut_isLoaded = true;
		}
		return els_dEtaOut_;
	}
	vector<float> &els_dPhiIn()
	{
		if (not els_dPhiIn_isLoaded) {
			if (els_dPhiIn_branch != 0) 
				els_dPhiIn_branch->GetEntry(index_);
			els_dPhiIn_isLoaded = true;
		}
		return els_dPhiIn_;
	}
	vector<float> &els_dPhiInPhiOut()
	{
		if (not els_dPhiInPhiOut_isLoaded) {
			if (els_dPhiInPhiOut_branch != 0) 
				els_dPhiInPhiOut_branch->GetEntry(index_);
			els_dPhiInPhiOut_isLoaded = true;
		}
		return els_dPhiInPhiOut_;
	}
	vector<float> &els_dPhiOut()
	{
		if (not els_dPhiOut_isLoaded) {
			if (els_dPhiOut_branch != 0) 
				els_dPhiOut_branch->GetEntry(index_);
			els_dPhiOut_isLoaded = true;
		}
		return els_dPhiOut_;
	}
	vector<float> &els_e3x3()
	{
		if (not els_e3x3_isLoaded) {
			if (els_e3x3_branch != 0) 
				els_e3x3_branch->GetEntry(index_);
			els_e3x3_isLoaded = true;
		}
		return els_e3x3_;
	}
	vector<float> &els_e5x5()
	{
		if (not els_e5x5_isLoaded) {
			if (els_e5x5_branch != 0) 
				els_e5x5_branch->GetEntry(index_);
			els_e5x5_isLoaded = true;
		}
		return els_e5x5_;
	}
	vector<float> &els_eOverPIn()
	{
		if (not els_eOverPIn_isLoaded) {
			if (els_eOverPIn_branch != 0) 
				els_eOverPIn_branch->GetEntry(index_);
			els_eOverPIn_isLoaded = true;
		}
		return els_eOverPIn_;
	}
	vector<float> &els_eSeedOverPOut()
	{
		if (not els_eSeedOverPOut_isLoaded) {
			if (els_eSeedOverPOut_branch != 0) 
				els_eSeedOverPOut_branch->GetEntry(index_);
			els_eSeedOverPOut_isLoaded = true;
		}
		return els_eSeedOverPOut_;
	}
	vector<float> &els_etaErr()
	{
		if (not els_etaErr_isLoaded) {
			if (els_etaErr_branch != 0) 
				els_etaErr_branch->GetEntry(index_);
			els_etaErr_isLoaded = true;
		}
		return els_etaErr_;
	}
	vector<float> &els_fBrem()
	{
		if (not els_fBrem_isLoaded) {
			if (els_fBrem_branch != 0) 
				els_fBrem_branch->GetEntry(index_);
			els_fBrem_isLoaded = true;
		}
		return els_fBrem_;
	}
	vector<float> &els_hOverE()
	{
		if (not els_hOverE_isLoaded) {
			if (els_hOverE_branch != 0) 
				els_hOverE_branch->GetEntry(index_);
			els_hOverE_isLoaded = true;
		}
		return els_hOverE_;
	}
	vector<float> &els_ndof()
	{
		if (not els_ndof_isLoaded) {
			if (els_ndof_branch != 0) 
				els_ndof_branch->GetEntry(index_);
			els_ndof_isLoaded = true;
		}
		return els_ndof_;
	}
	vector<float> &els_outerEta()
	{
		if (not els_outerEta_isLoaded) {
			if (els_outerEta_branch != 0) 
				els_outerEta_branch->GetEntry(index_);
			els_outerEta_isLoaded = true;
		}
		return els_outerEta_;
	}
	vector<float> &els_outerPhi()
	{
		if (not els_outerPhi_isLoaded) {
			if (els_outerPhi_branch != 0) 
				els_outerPhi_branch->GetEntry(index_);
			els_outerPhi_isLoaded = true;
		}
		return els_outerPhi_;
	}
	vector<float> &els_phiErr()
	{
		if (not els_phiErr_isLoaded) {
			if (els_phiErr_branch != 0) 
				els_phiErr_branch->GetEntry(index_);
			els_phiErr_isLoaded = true;
		}
		return els_phiErr_;
	}
	vector<float> &els_ptErr()
	{
		if (not els_ptErr_isLoaded) {
			if (els_ptErr_branch != 0) 
				els_ptErr_branch->GetEntry(index_);
			els_ptErr_isLoaded = true;
		}
		return els_ptErr_;
	}
	vector<float> &els_sigmaEtaEta()
	{
		if (not els_sigmaEtaEta_isLoaded) {
			if (els_sigmaEtaEta_branch != 0) 
				els_sigmaEtaEta_branch->GetEntry(index_);
			els_sigmaEtaEta_isLoaded = true;
		}
		return els_sigmaEtaEta_;
	}
	vector<float> &els_sigmaPhiPhi()
	{
		if (not els_sigmaPhiPhi_isLoaded) {
			if (els_sigmaPhiPhi_branch != 0) 
				els_sigmaPhiPhi_branch->GetEntry(index_);
			els_sigmaPhiPhi_isLoaded = true;
		}
		return els_sigmaPhiPhi_;
	}
	vector<float> &els_tkIso()
	{
		if (not els_tkIso_isLoaded) {
			if (els_tkIso_branch != 0) 
				els_tkIso_branch->GetEntry(index_);
			els_tkIso_isLoaded = true;
		}
		return els_tkIso_;
	}
	vector<float> &els_vertexphi()
	{
		if (not els_vertexphi_isLoaded) {
			if (els_vertexphi_branch != 0) 
				els_vertexphi_branch->GetEntry(index_);
			els_vertexphi_isLoaded = true;
		}
		return els_vertexphi_;
	}
	vector<float> &els_z0()
	{
		if (not els_z0_isLoaded) {
			if (els_z0_branch != 0) 
				els_z0_branch->GetEntry(index_);
			els_z0_isLoaded = true;
		}
		return els_z0_;
	}
	vector<float> &els_z0Err()
	{
		if (not els_z0Err_isLoaded) {
			if (els_z0Err_branch != 0) 
				els_z0Err_branch->GetEntry(index_);
			els_z0Err_isLoaded = true;
		}
		return els_z0Err_;
	}
	vector<float> &els_z0corr()
	{
		if (not els_z0corr_isLoaded) {
			if (els_z0corr_branch != 0) 
				els_z0corr_branch->GetEntry(index_);
			els_z0corr_isLoaded = true;
		}
		return els_z0corr_;
	}
	vector<float> &hyp_ll_chi2()
	{
		if (not hyp_ll_chi2_isLoaded) {
			if (hyp_ll_chi2_branch != 0) 
				hyp_ll_chi2_branch->GetEntry(index_);
			hyp_ll_chi2_isLoaded = true;
		}
		return hyp_ll_chi2_;
	}
	vector<float> &hyp_ll_d0()
	{
		if (not hyp_ll_d0_isLoaded) {
			if (hyp_ll_d0_branch != 0) 
				hyp_ll_d0_branch->GetEntry(index_);
			hyp_ll_d0_isLoaded = true;
		}
		return hyp_ll_d0_;
	}
	vector<float> &hyp_ll_d0Err()
	{
		if (not hyp_ll_d0Err_isLoaded) {
			if (hyp_ll_d0Err_branch != 0) 
				hyp_ll_d0Err_branch->GetEntry(index_);
			hyp_ll_d0Err_isLoaded = true;
		}
		return hyp_ll_d0Err_;
	}
	vector<float> &hyp_ll_d0corr()
	{
		if (not hyp_ll_d0corr_isLoaded) {
			if (hyp_ll_d0corr_branch != 0) 
				hyp_ll_d0corr_branch->GetEntry(index_);
			hyp_ll_d0corr_isLoaded = true;
		}
		return hyp_ll_d0corr_;
	}
	vector<float> &hyp_ll_etaErr()
	{
		if (not hyp_ll_etaErr_isLoaded) {
			if (hyp_ll_etaErr_branch != 0) 
				hyp_ll_etaErr_branch->GetEntry(index_);
			hyp_ll_etaErr_isLoaded = true;
		}
		return hyp_ll_etaErr_;
	}
	vector<float> &hyp_ll_iso()
	{
		if (not hyp_ll_iso_isLoaded) {
			if (hyp_ll_iso_branch != 0) 
				hyp_ll_iso_branch->GetEntry(index_);
			hyp_ll_iso_isLoaded = true;
		}
		return hyp_ll_iso_;
	}
	vector<float> &hyp_ll_ndof()
	{
		if (not hyp_ll_ndof_isLoaded) {
			if (hyp_ll_ndof_branch != 0) 
				hyp_ll_ndof_branch->GetEntry(index_);
			hyp_ll_ndof_isLoaded = true;
		}
		return hyp_ll_ndof_;
	}
	vector<float> &hyp_ll_outerEta()
	{
		if (not hyp_ll_outerEta_isLoaded) {
			if (hyp_ll_outerEta_branch != 0) 
				hyp_ll_outerEta_branch->GetEntry(index_);
			hyp_ll_outerEta_isLoaded = true;
		}
		return hyp_ll_outerEta_;
	}
	vector<float> &hyp_ll_outerPhi()
	{
		if (not hyp_ll_outerPhi_isLoaded) {
			if (hyp_ll_outerPhi_branch != 0) 
				hyp_ll_outerPhi_branch->GetEntry(index_);
			hyp_ll_outerPhi_isLoaded = true;
		}
		return hyp_ll_outerPhi_;
	}
	vector<float> &hyp_ll_phiErr()
	{
		if (not hyp_ll_phiErr_isLoaded) {
			if (hyp_ll_phiErr_branch != 0) 
				hyp_ll_phiErr_branch->GetEntry(index_);
			hyp_ll_phiErr_isLoaded = true;
		}
		return hyp_ll_phiErr_;
	}
	vector<float> &hyp_ll_ptErr()
	{
		if (not hyp_ll_ptErr_isLoaded) {
			if (hyp_ll_ptErr_branch != 0) 
				hyp_ll_ptErr_branch->GetEntry(index_);
			hyp_ll_ptErr_isLoaded = true;
		}
		return hyp_ll_ptErr_;
	}
	vector<float> &hyp_ll_tkIso()
	{
		if (not hyp_ll_tkIso_isLoaded) {
			if (hyp_ll_tkIso_branch != 0) 
				hyp_ll_tkIso_branch->GetEntry(index_);
			hyp_ll_tkIso_isLoaded = true;
		}
		return hyp_ll_tkIso_;
	}
	vector<float> &hyp_ll_vertexphi()
	{
		if (not hyp_ll_vertexphi_isLoaded) {
			if (hyp_ll_vertexphi_branch != 0) 
				hyp_ll_vertexphi_branch->GetEntry(index_);
			hyp_ll_vertexphi_isLoaded = true;
		}
		return hyp_ll_vertexphi_;
	}
	vector<float> &hyp_ll_z0()
	{
		if (not hyp_ll_z0_isLoaded) {
			if (hyp_ll_z0_branch != 0) 
				hyp_ll_z0_branch->GetEntry(index_);
			hyp_ll_z0_isLoaded = true;
		}
		return hyp_ll_z0_;
	}
	vector<float> &hyp_ll_z0Err()
	{
		if (not hyp_ll_z0Err_isLoaded) {
			if (hyp_ll_z0Err_branch != 0) 
				hyp_ll_z0Err_branch->GetEntry(index_);
			hyp_ll_z0Err_isLoaded = true;
		}
		return hyp_ll_z0Err_;
	}
	vector<float> &hyp_ll_z0corr()
	{
		if (not hyp_ll_z0corr_isLoaded) {
			if (hyp_ll_z0corr_branch != 0) 
				hyp_ll_z0corr_branch->GetEntry(index_);
			hyp_ll_z0corr_isLoaded = true;
		}
		return hyp_ll_z0corr_;
	}
	vector<float> &hyp_lt_chi2()
	{
		if (not hyp_lt_chi2_isLoaded) {
			if (hyp_lt_chi2_branch != 0) 
				hyp_lt_chi2_branch->GetEntry(index_);
			hyp_lt_chi2_isLoaded = true;
		}
		return hyp_lt_chi2_;
	}
	vector<float> &hyp_lt_d0()
	{
		if (not hyp_lt_d0_isLoaded) {
			if (hyp_lt_d0_branch != 0) 
				hyp_lt_d0_branch->GetEntry(index_);
			hyp_lt_d0_isLoaded = true;
		}
		return hyp_lt_d0_;
	}
	vector<float> &hyp_lt_d0Err()
	{
		if (not hyp_lt_d0Err_isLoaded) {
			if (hyp_lt_d0Err_branch != 0) 
				hyp_lt_d0Err_branch->GetEntry(index_);
			hyp_lt_d0Err_isLoaded = true;
		}
		return hyp_lt_d0Err_;
	}
	vector<float> &hyp_lt_d0corr()
	{
		if (not hyp_lt_d0corr_isLoaded) {
			if (hyp_lt_d0corr_branch != 0) 
				hyp_lt_d0corr_branch->GetEntry(index_);
			hyp_lt_d0corr_isLoaded = true;
		}
		return hyp_lt_d0corr_;
	}
	vector<float> &hyp_lt_etaErr()
	{
		if (not hyp_lt_etaErr_isLoaded) {
			if (hyp_lt_etaErr_branch != 0) 
				hyp_lt_etaErr_branch->GetEntry(index_);
			hyp_lt_etaErr_isLoaded = true;
		}
		return hyp_lt_etaErr_;
	}
	vector<float> &hyp_lt_iso()
	{
		if (not hyp_lt_iso_isLoaded) {
			if (hyp_lt_iso_branch != 0) 
				hyp_lt_iso_branch->GetEntry(index_);
			hyp_lt_iso_isLoaded = true;
		}
		return hyp_lt_iso_;
	}
	vector<float> &hyp_lt_ndof()
	{
		if (not hyp_lt_ndof_isLoaded) {
			if (hyp_lt_ndof_branch != 0) 
				hyp_lt_ndof_branch->GetEntry(index_);
			hyp_lt_ndof_isLoaded = true;
		}
		return hyp_lt_ndof_;
	}
	vector<float> &hyp_lt_outerEta()
	{
		if (not hyp_lt_outerEta_isLoaded) {
			if (hyp_lt_outerEta_branch != 0) 
				hyp_lt_outerEta_branch->GetEntry(index_);
			hyp_lt_outerEta_isLoaded = true;
		}
		return hyp_lt_outerEta_;
	}
	vector<float> &hyp_lt_outerPhi()
	{
		if (not hyp_lt_outerPhi_isLoaded) {
			if (hyp_lt_outerPhi_branch != 0) 
				hyp_lt_outerPhi_branch->GetEntry(index_);
			hyp_lt_outerPhi_isLoaded = true;
		}
		return hyp_lt_outerPhi_;
	}
	vector<float> &hyp_lt_phiErr()
	{
		if (not hyp_lt_phiErr_isLoaded) {
			if (hyp_lt_phiErr_branch != 0) 
				hyp_lt_phiErr_branch->GetEntry(index_);
			hyp_lt_phiErr_isLoaded = true;
		}
		return hyp_lt_phiErr_;
	}
	vector<float> &hyp_lt_ptErr()
	{
		if (not hyp_lt_ptErr_isLoaded) {
			if (hyp_lt_ptErr_branch != 0) 
				hyp_lt_ptErr_branch->GetEntry(index_);
			hyp_lt_ptErr_isLoaded = true;
		}
		return hyp_lt_ptErr_;
	}
	vector<float> &hyp_lt_tkIso()
	{
		if (not hyp_lt_tkIso_isLoaded) {
			if (hyp_lt_tkIso_branch != 0) 
				hyp_lt_tkIso_branch->GetEntry(index_);
			hyp_lt_tkIso_isLoaded = true;
		}
		return hyp_lt_tkIso_;
	}
	vector<float> &hyp_lt_vertexphi()
	{
		if (not hyp_lt_vertexphi_isLoaded) {
			if (hyp_lt_vertexphi_branch != 0) 
				hyp_lt_vertexphi_branch->GetEntry(index_);
			hyp_lt_vertexphi_isLoaded = true;
		}
		return hyp_lt_vertexphi_;
	}
	vector<float> &hyp_lt_z0()
	{
		if (not hyp_lt_z0_isLoaded) {
			if (hyp_lt_z0_branch != 0) 
				hyp_lt_z0_branch->GetEntry(index_);
			hyp_lt_z0_isLoaded = true;
		}
		return hyp_lt_z0_;
	}
	vector<float> &hyp_lt_z0Err()
	{
		if (not hyp_lt_z0Err_isLoaded) {
			if (hyp_lt_z0Err_branch != 0) 
				hyp_lt_z0Err_branch->GetEntry(index_);
			hyp_lt_z0Err_isLoaded = true;
		}
		return hyp_lt_z0Err_;
	}
	vector<float> &hyp_lt_z0corr()
	{
		if (not hyp_lt_z0corr_isLoaded) {
			if (hyp_lt_z0corr_branch != 0) 
				hyp_lt_z0corr_branch->GetEntry(index_);
			hyp_lt_z0corr_isLoaded = true;
		}
		return hyp_lt_z0corr_;
	}
	vector<float> &hyp_met()
	{
		if (not hyp_met_isLoaded) {
			if (hyp_met_branch != 0) 
				hyp_met_branch->GetEntry(index_);
			hyp_met_isLoaded = true;
		}
		return hyp_met_;
	}
	vector<float> &hyp_metAll()
	{
		if (not hyp_metAll_isLoaded) {
			if (hyp_metAll_branch != 0) 
				hyp_metAll_branch->GetEntry(index_);
			hyp_metAll_isLoaded = true;
		}
		return hyp_metAll_;
	}
	vector<float> &hyp_metAllCaloExp()
	{
		if (not hyp_metAllCaloExp_isLoaded) {
			if (hyp_metAllCaloExp_branch != 0) 
				hyp_metAllCaloExp_branch->GetEntry(index_);
			hyp_metAllCaloExp_isLoaded = true;
		}
		return hyp_metAllCaloExp_;
	}
	vector<float> &hyp_metCaloExp()
	{
		if (not hyp_metCaloExp_isLoaded) {
			if (hyp_metCaloExp_branch != 0) 
				hyp_metCaloExp_branch->GetEntry(index_);
			hyp_metCaloExp_isLoaded = true;
		}
		return hyp_metCaloExp_;
	}
	vector<float> &hyp_metCone()
	{
		if (not hyp_metCone_isLoaded) {
			if (hyp_metCone_branch != 0) 
				hyp_metCone_branch->GetEntry(index_);
			hyp_metCone_isLoaded = true;
		}
		return hyp_metCone_;
	}
	vector<float> &hyp_metDPhiJet10()
	{
		if (not hyp_metDPhiJet10_isLoaded) {
			if (hyp_metDPhiJet10_branch != 0) 
				hyp_metDPhiJet10_branch->GetEntry(index_);
			hyp_metDPhiJet10_isLoaded = true;
		}
		return hyp_metDPhiJet10_;
	}
	vector<float> &hyp_metDPhiJet15()
	{
		if (not hyp_metDPhiJet15_isLoaded) {
			if (hyp_metDPhiJet15_branch != 0) 
				hyp_metDPhiJet15_branch->GetEntry(index_);
			hyp_metDPhiJet15_isLoaded = true;
		}
		return hyp_metDPhiJet15_;
	}
	vector<float> &hyp_metDPhiJet20()
	{
		if (not hyp_metDPhiJet20_isLoaded) {
			if (hyp_metDPhiJet20_branch != 0) 
				hyp_metDPhiJet20_branch->GetEntry(index_);
			hyp_metDPhiJet20_isLoaded = true;
		}
		return hyp_metDPhiJet20_;
	}
	vector<float> &hyp_metDPhiTrk10()
	{
		if (not hyp_metDPhiTrk10_isLoaded) {
			if (hyp_metDPhiTrk10_branch != 0) 
				hyp_metDPhiTrk10_branch->GetEntry(index_);
			hyp_metDPhiTrk10_isLoaded = true;
		}
		return hyp_metDPhiTrk10_;
	}
	vector<float> &hyp_metDPhiTrk25()
	{
		if (not hyp_metDPhiTrk25_isLoaded) {
			if (hyp_metDPhiTrk25_branch != 0) 
				hyp_metDPhiTrk25_branch->GetEntry(index_);
			hyp_metDPhiTrk25_isLoaded = true;
		}
		return hyp_metDPhiTrk25_;
	}
	vector<float> &hyp_metDPhiTrk50()
	{
		if (not hyp_metDPhiTrk50_isLoaded) {
			if (hyp_metDPhiTrk50_branch != 0) 
				hyp_metDPhiTrk50_branch->GetEntry(index_);
			hyp_metDPhiTrk50_isLoaded = true;
		}
		return hyp_metDPhiTrk50_;
	}
	vector<float> &hyp_metJes10()
	{
		if (not hyp_metJes10_isLoaded) {
			if (hyp_metJes10_branch != 0) 
				hyp_metJes10_branch->GetEntry(index_);
			hyp_metJes10_isLoaded = true;
		}
		return hyp_metJes10_;
	}
	vector<float> &hyp_metJes15()
	{
		if (not hyp_metJes15_isLoaded) {
			if (hyp_metJes15_branch != 0) 
				hyp_metJes15_branch->GetEntry(index_);
			hyp_metJes15_isLoaded = true;
		}
		return hyp_metJes15_;
	}
	vector<float> &hyp_metJes30()
	{
		if (not hyp_metJes30_isLoaded) {
			if (hyp_metJes30_branch != 0) 
				hyp_metJes30_branch->GetEntry(index_);
			hyp_metJes30_isLoaded = true;
		}
		return hyp_metJes30_;
	}
	vector<float> &hyp_metJes5()
	{
		if (not hyp_metJes5_isLoaded) {
			if (hyp_metJes5_branch != 0) 
				hyp_metJes5_branch->GetEntry(index_);
			hyp_metJes5_isLoaded = true;
		}
		return hyp_metJes5_;
	}
	vector<float> &hyp_metJes50()
	{
		if (not hyp_metJes50_isLoaded) {
			if (hyp_metJes50_branch != 0) 
				hyp_metJes50_branch->GetEntry(index_);
			hyp_metJes50_isLoaded = true;
		}
		return hyp_metJes50_;
	}
	vector<float> &hyp_metNoCalo()
	{
		if (not hyp_metNoCalo_isLoaded) {
			if (hyp_metNoCalo_branch != 0) 
				hyp_metNoCalo_branch->GetEntry(index_);
			hyp_metNoCalo_isLoaded = true;
		}
		return hyp_metNoCalo_;
	}
	vector<float> &hyp_metPhi()
	{
		if (not hyp_metPhi_isLoaded) {
			if (hyp_metPhi_branch != 0) 
				hyp_metPhi_branch->GetEntry(index_);
			hyp_metPhi_isLoaded = true;
		}
		return hyp_metPhi_;
	}
	vector<float> &hyp_metPhiAll()
	{
		if (not hyp_metPhiAll_isLoaded) {
			if (hyp_metPhiAll_branch != 0) 
				hyp_metPhiAll_branch->GetEntry(index_);
			hyp_metPhiAll_isLoaded = true;
		}
		return hyp_metPhiAll_;
	}
	vector<float> &hyp_metPhiAllCaloExp()
	{
		if (not hyp_metPhiAllCaloExp_isLoaded) {
			if (hyp_metPhiAllCaloExp_branch != 0) 
				hyp_metPhiAllCaloExp_branch->GetEntry(index_);
			hyp_metPhiAllCaloExp_isLoaded = true;
		}
		return hyp_metPhiAllCaloExp_;
	}
	vector<float> &hyp_metPhiCaloExp()
	{
		if (not hyp_metPhiCaloExp_isLoaded) {
			if (hyp_metPhiCaloExp_branch != 0) 
				hyp_metPhiCaloExp_branch->GetEntry(index_);
			hyp_metPhiCaloExp_isLoaded = true;
		}
		return hyp_metPhiCaloExp_;
	}
	vector<float> &hyp_metPhiCone()
	{
		if (not hyp_metPhiCone_isLoaded) {
			if (hyp_metPhiCone_branch != 0) 
				hyp_metPhiCone_branch->GetEntry(index_);
			hyp_metPhiCone_isLoaded = true;
		}
		return hyp_metPhiCone_;
	}
	vector<float> &hyp_metPhiJes10()
	{
		if (not hyp_metPhiJes10_isLoaded) {
			if (hyp_metPhiJes10_branch != 0) 
				hyp_metPhiJes10_branch->GetEntry(index_);
			hyp_metPhiJes10_isLoaded = true;
		}
		return hyp_metPhiJes10_;
	}
	vector<float> &hyp_metPhiJes15()
	{
		if (not hyp_metPhiJes15_isLoaded) {
			if (hyp_metPhiJes15_branch != 0) 
				hyp_metPhiJes15_branch->GetEntry(index_);
			hyp_metPhiJes15_isLoaded = true;
		}
		return hyp_metPhiJes15_;
	}
	vector<float> &hyp_metPhiJes30()
	{
		if (not hyp_metPhiJes30_isLoaded) {
			if (hyp_metPhiJes30_branch != 0) 
				hyp_metPhiJes30_branch->GetEntry(index_);
			hyp_metPhiJes30_isLoaded = true;
		}
		return hyp_metPhiJes30_;
	}
	vector<float> &hyp_metPhiJes5()
	{
		if (not hyp_metPhiJes5_isLoaded) {
			if (hyp_metPhiJes5_branch != 0) 
				hyp_metPhiJes5_branch->GetEntry(index_);
			hyp_metPhiJes5_isLoaded = true;
		}
		return hyp_metPhiJes5_;
	}
	vector<float> &hyp_metPhiJes50()
	{
		if (not hyp_metPhiJes50_isLoaded) {
			if (hyp_metPhiJes50_branch != 0) 
				hyp_metPhiJes50_branch->GetEntry(index_);
			hyp_metPhiJes50_isLoaded = true;
		}
		return hyp_metPhiJes50_;
	}
	vector<float> &hyp_metPhiNoCalo()
	{
		if (not hyp_metPhiNoCalo_isLoaded) {
			if (hyp_metPhiNoCalo_branch != 0) 
				hyp_metPhiNoCalo_branch->GetEntry(index_);
			hyp_metPhiNoCalo_isLoaded = true;
		}
		return hyp_metPhiNoCalo_;
	}
	vector<float> &hyp_quadlep_met()
	{
		if (not hyp_quadlep_met_isLoaded) {
			if (hyp_quadlep_met_branch != 0) 
				hyp_quadlep_met_branch->GetEntry(index_);
			hyp_quadlep_met_isLoaded = true;
		}
		return hyp_quadlep_met_;
	}
	vector<float> &hyp_quadlep_metAll()
	{
		if (not hyp_quadlep_metAll_isLoaded) {
			if (hyp_quadlep_metAll_branch != 0) 
				hyp_quadlep_metAll_branch->GetEntry(index_);
			hyp_quadlep_metAll_isLoaded = true;
		}
		return hyp_quadlep_metAll_;
	}
	vector<float> &hyp_trilep_met()
	{
		if (not hyp_trilep_met_isLoaded) {
			if (hyp_trilep_met_branch != 0) 
				hyp_trilep_met_branch->GetEntry(index_);
			hyp_trilep_met_isLoaded = true;
		}
		return hyp_trilep_met_;
	}
	vector<float> &hyp_trilep_metAll()
	{
		if (not hyp_trilep_metAll_isLoaded) {
			if (hyp_trilep_metAll_branch != 0) 
				hyp_trilep_metAll_branch->GetEntry(index_);
			hyp_trilep_metAll_isLoaded = true;
		}
		return hyp_trilep_metAll_;
	}
	vector<float> &jets_EMFcor()
	{
		if (not jets_EMFcor_isLoaded) {
			if (jets_EMFcor_branch != 0) 
				jets_EMFcor_branch->GetEntry(index_);
			jets_EMFcor_isLoaded = true;
		}
		return jets_EMFcor_;
	}
	vector<float> &jets_chFrac()
	{
		if (not jets_chFrac_isLoaded) {
			if (jets_chFrac_branch != 0) 
				jets_chFrac_branch->GetEntry(index_);
			jets_chFrac_isLoaded = true;
		}
		return jets_chFrac_;
	}
	vector<float> &jets_cor()
	{
		if (not jets_cor_isLoaded) {
			if (jets_cor_branch != 0) 
				jets_cor_branch->GetEntry(index_);
			jets_cor_isLoaded = true;
		}
		return jets_cor_;
	}
	vector<float> &jets_emFrac()
	{
		if (not jets_emFrac_isLoaded) {
			if (jets_emFrac_branch != 0) 
				jets_emFrac_branch->GetEntry(index_);
			jets_emFrac_isLoaded = true;
		}
		return jets_emFrac_;
	}
	vector<float> &mus_eledr()
	{
		if (not mus_eledr_isLoaded) {
			if (mus_eledr_branch != 0) 
				mus_eledr_branch->GetEntry(index_);
			mus_eledr_isLoaded = true;
		}
		return mus_eledr_;
	}
	vector<float> &mus_jetdr()
	{
		if (not mus_jetdr_isLoaded) {
			if (mus_jetdr_branch != 0) 
				mus_jetdr_branch->GetEntry(index_);
			mus_jetdr_isLoaded = true;
		}
		return mus_jetdr_;
	}
	vector<float> &mus_trkdr()
	{
		if (not mus_trkdr_isLoaded) {
			if (mus_trkdr_branch != 0) 
				mus_trkdr_branch->GetEntry(index_);
			mus_trkdr_isLoaded = true;
		}
		return mus_trkdr_;
	}
	vector<float> &mus_chi2()
	{
		if (not mus_chi2_isLoaded) {
			if (mus_chi2_branch != 0) 
				mus_chi2_branch->GetEntry(index_);
			mus_chi2_isLoaded = true;
		}
		return mus_chi2_;
	}
	vector<float> &mus_d0()
	{
		if (not mus_d0_isLoaded) {
			if (mus_d0_branch != 0) 
				mus_d0_branch->GetEntry(index_);
			mus_d0_isLoaded = true;
		}
		return mus_d0_;
	}
	vector<float> &mus_d0Err()
	{
		if (not mus_d0Err_isLoaded) {
			if (mus_d0Err_branch != 0) 
				mus_d0Err_branch->GetEntry(index_);
			mus_d0Err_isLoaded = true;
		}
		return mus_d0Err_;
	}
	vector<float> &mus_d0corr()
	{
		if (not mus_d0corr_isLoaded) {
			if (mus_d0corr_branch != 0) 
				mus_d0corr_branch->GetEntry(index_);
			mus_d0corr_isLoaded = true;
		}
		return mus_d0corr_;
	}
	vector<float> &mus_e_em()
	{
		if (not mus_e_em_isLoaded) {
			if (mus_e_em_branch != 0) 
				mus_e_em_branch->GetEntry(index_);
			mus_e_em_isLoaded = true;
		}
		return mus_e_em_;
	}
	vector<float> &mus_e_emS9()
	{
		if (not mus_e_emS9_isLoaded) {
			if (mus_e_emS9_branch != 0) 
				mus_e_emS9_branch->GetEntry(index_);
			mus_e_emS9_isLoaded = true;
		}
		return mus_e_emS9_;
	}
	vector<float> &mus_e_had()
	{
		if (not mus_e_had_isLoaded) {
			if (mus_e_had_branch != 0) 
				mus_e_had_branch->GetEntry(index_);
			mus_e_had_isLoaded = true;
		}
		return mus_e_had_;
	}
	vector<float> &mus_e_hadS9()
	{
		if (not mus_e_hadS9_isLoaded) {
			if (mus_e_hadS9_branch != 0) 
				mus_e_hadS9_branch->GetEntry(index_);
			mus_e_hadS9_isLoaded = true;
		}
		return mus_e_hadS9_;
	}
	vector<float> &mus_e_ho()
	{
		if (not mus_e_ho_isLoaded) {
			if (mus_e_ho_branch != 0) 
				mus_e_ho_branch->GetEntry(index_);
			mus_e_ho_isLoaded = true;
		}
		return mus_e_ho_;
	}
	vector<float> &mus_e_hoS9()
	{
		if (not mus_e_hoS9_isLoaded) {
			if (mus_e_hoS9_branch != 0) 
				mus_e_hoS9_branch->GetEntry(index_);
			mus_e_hoS9_isLoaded = true;
		}
		return mus_e_hoS9_;
	}
	vector<float> &mus_etaErr()
	{
		if (not mus_etaErr_isLoaded) {
			if (mus_etaErr_branch != 0) 
				mus_etaErr_branch->GetEntry(index_);
			mus_etaErr_isLoaded = true;
		}
		return mus_etaErr_;
	}
	vector<float> &mus_gfit_chi2()
	{
		if (not mus_gfit_chi2_isLoaded) {
			if (mus_gfit_chi2_branch != 0) 
				mus_gfit_chi2_branch->GetEntry(index_);
			mus_gfit_chi2_isLoaded = true;
		}
		return mus_gfit_chi2_;
	}
	vector<float> &mus_gfit_ndof()
	{
		if (not mus_gfit_ndof_isLoaded) {
			if (mus_gfit_ndof_branch != 0) 
				mus_gfit_ndof_branch->GetEntry(index_);
			mus_gfit_ndof_isLoaded = true;
		}
		return mus_gfit_ndof_;
	}
	vector<float> &mus_iso()
	{
		if (not mus_iso_isLoaded) {
			if (mus_iso_branch != 0) 
				mus_iso_branch->GetEntry(index_);
			mus_iso_isLoaded = true;
		}
		return mus_iso_;
	}
	vector<float> &mus_iso03_emEt()
	{
		if (not mus_iso03_emEt_isLoaded) {
			if (mus_iso03_emEt_branch != 0) 
				mus_iso03_emEt_branch->GetEntry(index_);
			mus_iso03_emEt_isLoaded = true;
		}
		return mus_iso03_emEt_;
	}
	vector<float> &mus_iso03_hadEt()
	{
		if (not mus_iso03_hadEt_isLoaded) {
			if (mus_iso03_hadEt_branch != 0) 
				mus_iso03_hadEt_branch->GetEntry(index_);
			mus_iso03_hadEt_isLoaded = true;
		}
		return mus_iso03_hadEt_;
	}
	vector<float> &mus_iso03_hoEt()
	{
		if (not mus_iso03_hoEt_isLoaded) {
			if (mus_iso03_hoEt_branch != 0) 
				mus_iso03_hoEt_branch->GetEntry(index_);
			mus_iso03_hoEt_isLoaded = true;
		}
		return mus_iso03_hoEt_;
	}
	vector<float> &mus_iso03_sumPt()
	{
		if (not mus_iso03_sumPt_isLoaded) {
			if (mus_iso03_sumPt_branch != 0) 
				mus_iso03_sumPt_branch->GetEntry(index_);
			mus_iso03_sumPt_isLoaded = true;
		}
		return mus_iso03_sumPt_;
	}
	vector<float> &mus_iso05_emEt()
	{
		if (not mus_iso05_emEt_isLoaded) {
			if (mus_iso05_emEt_branch != 0) 
				mus_iso05_emEt_branch->GetEntry(index_);
			mus_iso05_emEt_isLoaded = true;
		}
		return mus_iso05_emEt_;
	}
	vector<float> &mus_iso05_hadEt()
	{
		if (not mus_iso05_hadEt_isLoaded) {
			if (mus_iso05_hadEt_branch != 0) 
				mus_iso05_hadEt_branch->GetEntry(index_);
			mus_iso05_hadEt_isLoaded = true;
		}
		return mus_iso05_hadEt_;
	}
	vector<float> &mus_iso05_hoEt()
	{
		if (not mus_iso05_hoEt_isLoaded) {
			if (mus_iso05_hoEt_branch != 0) 
				mus_iso05_hoEt_branch->GetEntry(index_);
			mus_iso05_hoEt_isLoaded = true;
		}
		return mus_iso05_hoEt_;
	}
	vector<float> &mus_iso05_sumPt()
	{
		if (not mus_iso05_sumPt_isLoaded) {
			if (mus_iso05_sumPt_branch != 0) 
				mus_iso05_sumPt_branch->GetEntry(index_);
			mus_iso05_sumPt_isLoaded = true;
		}
		return mus_iso05_sumPt_;
	}
	vector<float> &mus_ndof()
	{
		if (not mus_ndof_isLoaded) {
			if (mus_ndof_branch != 0) 
				mus_ndof_branch->GetEntry(index_);
			mus_ndof_isLoaded = true;
		}
		return mus_ndof_;
	}
	vector<float> &mus_outerEta()
	{
		if (not mus_outerEta_isLoaded) {
			if (mus_outerEta_branch != 0) 
				mus_outerEta_branch->GetEntry(index_);
			mus_outerEta_isLoaded = true;
		}
		return mus_outerEta_;
	}
	vector<float> &mus_outerPhi()
	{
		if (not mus_outerPhi_isLoaded) {
			if (mus_outerPhi_branch != 0) 
				mus_outerPhi_branch->GetEntry(index_);
			mus_outerPhi_isLoaded = true;
		}
		return mus_outerPhi_;
	}
	vector<float> &mus_phiErr()
	{
		if (not mus_phiErr_isLoaded) {
			if (mus_phiErr_branch != 0) 
				mus_phiErr_branch->GetEntry(index_);
			mus_phiErr_isLoaded = true;
		}
		return mus_phiErr_;
	}
	vector<float> &mus_ptErr()
	{
		if (not mus_ptErr_isLoaded) {
			if (mus_ptErr_branch != 0) 
				mus_ptErr_branch->GetEntry(index_);
			mus_ptErr_isLoaded = true;
		}
		return mus_ptErr_;
	}
	vector<float> &mus_vertexphi()
	{
		if (not mus_vertexphi_isLoaded) {
			if (mus_vertexphi_branch != 0) 
				mus_vertexphi_branch->GetEntry(index_);
			mus_vertexphi_isLoaded = true;
		}
		return mus_vertexphi_;
	}
	vector<float> &mus_z0()
	{
		if (not mus_z0_isLoaded) {
			if (mus_z0_branch != 0) 
				mus_z0_branch->GetEntry(index_);
			mus_z0_isLoaded = true;
		}
		return mus_z0_;
	}
	vector<float> &mus_z0Err()
	{
		if (not mus_z0Err_isLoaded) {
			if (mus_z0Err_branch != 0) 
				mus_z0Err_branch->GetEntry(index_);
			mus_z0Err_isLoaded = true;
		}
		return mus_z0Err_;
	}
	vector<float> &mus_z0corr()
	{
		if (not mus_z0corr_isLoaded) {
			if (mus_z0corr_branch != 0) 
				mus_z0corr_branch->GetEntry(index_);
			mus_z0corr_isLoaded = true;
		}
		return mus_z0corr_;
	}
	vector<float> &els_pat_caloIso()
	{
		if (not els_pat_caloIso_isLoaded) {
			if (els_pat_caloIso_branch != 0) 
				els_pat_caloIso_branch->GetEntry(index_);
			els_pat_caloIso_isLoaded = true;
		}
		return els_pat_caloIso_;
	}
	vector<float> &els_pat_ecalIso()
	{
		if (not els_pat_ecalIso_isLoaded) {
			if (els_pat_ecalIso_branch != 0) 
				els_pat_ecalIso_branch->GetEntry(index_);
			els_pat_ecalIso_isLoaded = true;
		}
		return els_pat_ecalIso_;
	}
	vector<float> &els_pat_hcalIso()
	{
		if (not els_pat_hcalIso_isLoaded) {
			if (els_pat_hcalIso_branch != 0) 
				els_pat_hcalIso_branch->GetEntry(index_);
			els_pat_hcalIso_isLoaded = true;
		}
		return els_pat_hcalIso_;
	}
	vector<float> &els_pat_looseId()
	{
		if (not els_pat_looseId_isLoaded) {
			if (els_pat_looseId_branch != 0) 
				els_pat_looseId_branch->GetEntry(index_);
			els_pat_looseId_isLoaded = true;
		}
		return els_pat_looseId_;
	}
	vector<float> &els_pat_robustLooseId()
	{
		if (not els_pat_robustLooseId_isLoaded) {
			if (els_pat_robustLooseId_branch != 0) 
				els_pat_robustLooseId_branch->GetEntry(index_);
			els_pat_robustLooseId_isLoaded = true;
		}
		return els_pat_robustLooseId_;
	}
	vector<float> &els_pat_robustTightId()
	{
		if (not els_pat_robustTightId_isLoaded) {
			if (els_pat_robustTightId_branch != 0) 
				els_pat_robustTightId_branch->GetEntry(index_);
			els_pat_robustTightId_isLoaded = true;
		}
		return els_pat_robustTightId_;
	}
	vector<float> &els_pat_tightId()
	{
		if (not els_pat_tightId_isLoaded) {
			if (els_pat_tightId_branch != 0) 
				els_pat_tightId_branch->GetEntry(index_);
			els_pat_tightId_isLoaded = true;
		}
		return els_pat_tightId_;
	}
	vector<float> &els_pat_trackIso()
	{
		if (not els_pat_trackIso_isLoaded) {
			if (els_pat_trackIso_branch != 0) 
				els_pat_trackIso_branch->GetEntry(index_);
			els_pat_trackIso_isLoaded = true;
		}
		return els_pat_trackIso_;
	}
	vector<float> &jets_pat_bCorrF()
	{
		if (not jets_pat_bCorrF_isLoaded) {
			if (jets_pat_bCorrF_branch != 0) 
				jets_pat_bCorrF_branch->GetEntry(index_);
			jets_pat_bCorrF_isLoaded = true;
		}
		return jets_pat_bCorrF_;
	}
	vector<float> &jets_pat_cCorrF()
	{
		if (not jets_pat_cCorrF_isLoaded) {
			if (jets_pat_cCorrF_branch != 0) 
				jets_pat_cCorrF_branch->GetEntry(index_);
			jets_pat_cCorrF_isLoaded = true;
		}
		return jets_pat_cCorrF_;
	}
	vector<float> &jets_pat_gluCorrF()
	{
		if (not jets_pat_gluCorrF_isLoaded) {
			if (jets_pat_gluCorrF_branch != 0) 
				jets_pat_gluCorrF_branch->GetEntry(index_);
			jets_pat_gluCorrF_isLoaded = true;
		}
		return jets_pat_gluCorrF_;
	}
	vector<float> &jets_pat_jetCharge()
	{
		if (not jets_pat_jetCharge_isLoaded) {
			if (jets_pat_jetCharge_branch != 0) 
				jets_pat_jetCharge_branch->GetEntry(index_);
			jets_pat_jetCharge_isLoaded = true;
		}
		return jets_pat_jetCharge_;
	}
	vector<float> &jets_pat_noCorrF()
	{
		if (not jets_pat_noCorrF_isLoaded) {
			if (jets_pat_noCorrF_branch != 0) 
				jets_pat_noCorrF_branch->GetEntry(index_);
			jets_pat_noCorrF_isLoaded = true;
		}
		return jets_pat_noCorrF_;
	}
	vector<float> &jets_pat_udsCorrF()
	{
		if (not jets_pat_udsCorrF_isLoaded) {
			if (jets_pat_udsCorrF_branch != 0) 
				jets_pat_udsCorrF_branch->GetEntry(index_);
			jets_pat_udsCorrF_isLoaded = true;
		}
		return jets_pat_udsCorrF_;
	}
	vector<float> &mus_pat_caloIso()
	{
		if (not mus_pat_caloIso_isLoaded) {
			if (mus_pat_caloIso_branch != 0) 
				mus_pat_caloIso_branch->GetEntry(index_);
			mus_pat_caloIso_isLoaded = true;
		}
		return mus_pat_caloIso_;
	}
	vector<float> &mus_pat_ecalIso()
	{
		if (not mus_pat_ecalIso_isLoaded) {
			if (mus_pat_ecalIso_branch != 0) 
				mus_pat_ecalIso_branch->GetEntry(index_);
			mus_pat_ecalIso_isLoaded = true;
		}
		return mus_pat_ecalIso_;
	}
	vector<float> &mus_pat_ecalvetoDep()
	{
		if (not mus_pat_ecalvetoDep_isLoaded) {
			if (mus_pat_ecalvetoDep_branch != 0) 
				mus_pat_ecalvetoDep_branch->GetEntry(index_);
			mus_pat_ecalvetoDep_isLoaded = true;
		}
		return mus_pat_ecalvetoDep_;
	}
	vector<float> &mus_pat_hcalIso()
	{
		if (not mus_pat_hcalIso_isLoaded) {
			if (mus_pat_hcalIso_branch != 0) 
				mus_pat_hcalIso_branch->GetEntry(index_);
			mus_pat_hcalIso_isLoaded = true;
		}
		return mus_pat_hcalIso_;
	}
	vector<float> &mus_pat_hcalvetoDep()
	{
		if (not mus_pat_hcalvetoDep_isLoaded) {
			if (mus_pat_hcalvetoDep_branch != 0) 
				mus_pat_hcalvetoDep_branch->GetEntry(index_);
			mus_pat_hcalvetoDep_isLoaded = true;
		}
		return mus_pat_hcalvetoDep_;
	}
	vector<float> &mus_pat_leptonID()
	{
		if (not mus_pat_leptonID_isLoaded) {
			if (mus_pat_leptonID_branch != 0) 
				mus_pat_leptonID_branch->GetEntry(index_);
			mus_pat_leptonID_isLoaded = true;
		}
		return mus_pat_leptonID_;
	}
	vector<float> &mus_pat_trackIso()
	{
		if (not mus_pat_trackIso_isLoaded) {
			if (mus_pat_trackIso_branch != 0) 
				mus_pat_trackIso_branch->GetEntry(index_);
			mus_pat_trackIso_isLoaded = true;
		}
		return mus_pat_trackIso_;
	}
	vector<float> &mus_pat_vetoDep()
	{
		if (not mus_pat_vetoDep_isLoaded) {
			if (mus_pat_vetoDep_branch != 0) 
				mus_pat_vetoDep_branch->GetEntry(index_);
			mus_pat_vetoDep_isLoaded = true;
		}
		return mus_pat_vetoDep_;
	}
	vector<float> &trks_chi2()
	{
		if (not trks_chi2_isLoaded) {
			if (trks_chi2_branch != 0) 
				trks_chi2_branch->GetEntry(index_);
			trks_chi2_isLoaded = true;
		}
		return trks_chi2_;
	}
	vector<float> &trks_d0()
	{
		if (not trks_d0_isLoaded) {
			if (trks_d0_branch != 0) 
				trks_d0_branch->GetEntry(index_);
			trks_d0_isLoaded = true;
		}
		return trks_d0_;
	}
	vector<float> &trks_d0Err()
	{
		if (not trks_d0Err_isLoaded) {
			if (trks_d0Err_branch != 0) 
				trks_d0Err_branch->GetEntry(index_);
			trks_d0Err_isLoaded = true;
		}
		return trks_d0Err_;
	}
	vector<float> &trks_d0corr()
	{
		if (not trks_d0corr_isLoaded) {
			if (trks_d0corr_branch != 0) 
				trks_d0corr_branch->GetEntry(index_);
			trks_d0corr_isLoaded = true;
		}
		return trks_d0corr_;
	}
	vector<float> &trks_etaErr()
	{
		if (not trks_etaErr_isLoaded) {
			if (trks_etaErr_branch != 0) 
				trks_etaErr_branch->GetEntry(index_);
			trks_etaErr_isLoaded = true;
		}
		return trks_etaErr_;
	}
	vector<float> &trks_ndof()
	{
		if (not trks_ndof_isLoaded) {
			if (trks_ndof_branch != 0) 
				trks_ndof_branch->GetEntry(index_);
			trks_ndof_isLoaded = true;
		}
		return trks_ndof_;
	}
	vector<float> &trks_outerEta()
	{
		if (not trks_outerEta_isLoaded) {
			if (trks_outerEta_branch != 0) 
				trks_outerEta_branch->GetEntry(index_);
			trks_outerEta_isLoaded = true;
		}
		return trks_outerEta_;
	}
	vector<float> &trks_outerPhi()
	{
		if (not trks_outerPhi_isLoaded) {
			if (trks_outerPhi_branch != 0) 
				trks_outerPhi_branch->GetEntry(index_);
			trks_outerPhi_isLoaded = true;
		}
		return trks_outerPhi_;
	}
	vector<float> &trks_phiErr()
	{
		if (not trks_phiErr_isLoaded) {
			if (trks_phiErr_branch != 0) 
				trks_phiErr_branch->GetEntry(index_);
			trks_phiErr_isLoaded = true;
		}
		return trks_phiErr_;
	}
	vector<float> &trks_ptErr()
	{
		if (not trks_ptErr_isLoaded) {
			if (trks_ptErr_branch != 0) 
				trks_ptErr_branch->GetEntry(index_);
			trks_ptErr_isLoaded = true;
		}
		return trks_ptErr_;
	}
	vector<float> &trks_vertexphi()
	{
		if (not trks_vertexphi_isLoaded) {
			if (trks_vertexphi_branch != 0) 
				trks_vertexphi_branch->GetEntry(index_);
			trks_vertexphi_isLoaded = true;
		}
		return trks_vertexphi_;
	}
	vector<float> &trks_z0()
	{
		if (not trks_z0_isLoaded) {
			if (trks_z0_branch != 0) 
				trks_z0_branch->GetEntry(index_);
			trks_z0_isLoaded = true;
		}
		return trks_z0_;
	}
	vector<float> &trks_z0Err()
	{
		if (not trks_z0Err_isLoaded) {
			if (trks_z0Err_branch != 0) 
				trks_z0Err_branch->GetEntry(index_);
			trks_z0Err_isLoaded = true;
		}
		return trks_z0Err_;
	}
	vector<float> &trks_z0corr()
	{
		if (not trks_z0corr_isLoaded) {
			if (trks_z0corr_branch != 0) 
				trks_z0corr_branch->GetEntry(index_);
			trks_z0corr_isLoaded = true;
		}
		return trks_z0corr_;
	}
	vector<float> &trks_elsdr()
	{
		if (not trks_elsdr_isLoaded) {
			if (trks_elsdr_branch != 0) 
				trks_elsdr_branch->GetEntry(index_);
			trks_elsdr_isLoaded = true;
		}
		return trks_elsdr_;
	}
	vector<float> &trks_elsshFrac()
	{
		if (not trks_elsshFrac_isLoaded) {
			if (trks_elsshFrac_branch != 0) 
				trks_elsshFrac_branch->GetEntry(index_);
			trks_elsshFrac_isLoaded = true;
		}
		return trks_elsshFrac_;
	}
	vector<float> &trk_musdr()
	{
		if (not trk_musdr_isLoaded) {
			if (trk_musdr_branch != 0) 
				trk_musdr_branch->GetEntry(index_);
			trk_musdr_isLoaded = true;
		}
		return trk_musdr_;
	}
	vector<vector<float> > &hyp_jets_EMFcor()
	{
		if (not hyp_jets_EMFcor_isLoaded) {
			if (hyp_jets_EMFcor_branch != 0) 
				hyp_jets_EMFcor_branch->GetEntry(index_);
			hyp_jets_EMFcor_isLoaded = true;
		}
		return hyp_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_jets_chFrac()
	{
		if (not hyp_jets_chFrac_isLoaded) {
			if (hyp_jets_chFrac_branch != 0) 
				hyp_jets_chFrac_branch->GetEntry(index_);
			hyp_jets_chFrac_isLoaded = true;
		}
		return hyp_jets_chFrac_;
	}
	vector<vector<float> > &hyp_jets_cor()
	{
		if (not hyp_jets_cor_isLoaded) {
			if (hyp_jets_cor_branch != 0) 
				hyp_jets_cor_branch->GetEntry(index_);
			hyp_jets_cor_isLoaded = true;
		}
		return hyp_jets_cor_;
	}
	vector<vector<float> > &hyp_jets_emFrac()
	{
		if (not hyp_jets_emFrac_isLoaded) {
			if (hyp_jets_emFrac_branch != 0) 
				hyp_jets_emFrac_branch->GetEntry(index_);
			hyp_jets_emFrac_isLoaded = true;
		}
		return hyp_jets_emFrac_;
	}
	vector<vector<float> > &hyp_jets_mc_emEnergy()
	{
		if (not hyp_jets_mc_emEnergy_isLoaded) {
			if (hyp_jets_mc_emEnergy_branch != 0) 
				hyp_jets_mc_emEnergy_branch->GetEntry(index_);
			hyp_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_hadEnergy()
	{
		if (not hyp_jets_mc_hadEnergy_isLoaded) {
			if (hyp_jets_mc_hadEnergy_branch != 0) 
				hyp_jets_mc_hadEnergy_branch->GetEntry(index_);
			hyp_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_invEnergy()
	{
		if (not hyp_jets_mc_invEnergy_isLoaded) {
			if (hyp_jets_mc_invEnergy_branch != 0) 
				hyp_jets_mc_invEnergy_branch->GetEntry(index_);
			hyp_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_jets_mc_otherEnergy()
	{
		if (not hyp_jets_mc_otherEnergy_isLoaded) {
			if (hyp_jets_mc_otherEnergy_branch != 0) 
				hyp_jets_mc_otherEnergy_branch->GetEntry(index_);
			hyp_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_jets_pat_bCorrF()
	{
		if (not hyp_jets_pat_bCorrF_isLoaded) {
			if (hyp_jets_pat_bCorrF_branch != 0) 
				hyp_jets_pat_bCorrF_branch->GetEntry(index_);
			hyp_jets_pat_bCorrF_isLoaded = true;
		}
		return hyp_jets_pat_bCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_cCorrF()
	{
		if (not hyp_jets_pat_cCorrF_isLoaded) {
			if (hyp_jets_pat_cCorrF_branch != 0) 
				hyp_jets_pat_cCorrF_branch->GetEntry(index_);
			hyp_jets_pat_cCorrF_isLoaded = true;
		}
		return hyp_jets_pat_cCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_gluCorrF()
	{
		if (not hyp_jets_pat_gluCorrF_isLoaded) {
			if (hyp_jets_pat_gluCorrF_branch != 0) 
				hyp_jets_pat_gluCorrF_branch->GetEntry(index_);
			hyp_jets_pat_gluCorrF_isLoaded = true;
		}
		return hyp_jets_pat_gluCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_jetCharge()
	{
		if (not hyp_jets_pat_jetCharge_isLoaded) {
			if (hyp_jets_pat_jetCharge_branch != 0) 
				hyp_jets_pat_jetCharge_branch->GetEntry(index_);
			hyp_jets_pat_jetCharge_isLoaded = true;
		}
		return hyp_jets_pat_jetCharge_;
	}
	vector<vector<float> > &hyp_jets_pat_noCorrF()
	{
		if (not hyp_jets_pat_noCorrF_isLoaded) {
			if (hyp_jets_pat_noCorrF_branch != 0) 
				hyp_jets_pat_noCorrF_branch->GetEntry(index_);
			hyp_jets_pat_noCorrF_isLoaded = true;
		}
		return hyp_jets_pat_noCorrF_;
	}
	vector<vector<float> > &hyp_jets_pat_udsCorrF()
	{
		if (not hyp_jets_pat_udsCorrF_isLoaded) {
			if (hyp_jets_pat_udsCorrF_branch != 0) 
				hyp_jets_pat_udsCorrF_branch->GetEntry(index_);
			hyp_jets_pat_udsCorrF_isLoaded = true;
		}
		return hyp_jets_pat_udsCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_EMFcor()
	{
		if (not hyp_other_jets_EMFcor_isLoaded) {
			if (hyp_other_jets_EMFcor_branch != 0) 
				hyp_other_jets_EMFcor_branch->GetEntry(index_);
			hyp_other_jets_EMFcor_isLoaded = true;
		}
		return hyp_other_jets_EMFcor_;
	}
	vector<vector<float> > &hyp_other_jets_chFrac()
	{
		if (not hyp_other_jets_chFrac_isLoaded) {
			if (hyp_other_jets_chFrac_branch != 0) 
				hyp_other_jets_chFrac_branch->GetEntry(index_);
			hyp_other_jets_chFrac_isLoaded = true;
		}
		return hyp_other_jets_chFrac_;
	}
	vector<vector<float> > &hyp_other_jets_cor()
	{
		if (not hyp_other_jets_cor_isLoaded) {
			if (hyp_other_jets_cor_branch != 0) 
				hyp_other_jets_cor_branch->GetEntry(index_);
			hyp_other_jets_cor_isLoaded = true;
		}
		return hyp_other_jets_cor_;
	}
	vector<vector<float> > &hyp_other_jets_emFrac()
	{
		if (not hyp_other_jets_emFrac_isLoaded) {
			if (hyp_other_jets_emFrac_branch != 0) 
				hyp_other_jets_emFrac_branch->GetEntry(index_);
			hyp_other_jets_emFrac_isLoaded = true;
		}
		return hyp_other_jets_emFrac_;
	}
	vector<vector<float> > &hyp_other_jets_mc_emEnergy()
	{
		if (not hyp_other_jets_mc_emEnergy_isLoaded) {
			if (hyp_other_jets_mc_emEnergy_branch != 0) 
				hyp_other_jets_mc_emEnergy_branch->GetEntry(index_);
			hyp_other_jets_mc_emEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_emEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_hadEnergy()
	{
		if (not hyp_other_jets_mc_hadEnergy_isLoaded) {
			if (hyp_other_jets_mc_hadEnergy_branch != 0) 
				hyp_other_jets_mc_hadEnergy_branch->GetEntry(index_);
			hyp_other_jets_mc_hadEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_hadEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_invEnergy()
	{
		if (not hyp_other_jets_mc_invEnergy_isLoaded) {
			if (hyp_other_jets_mc_invEnergy_branch != 0) 
				hyp_other_jets_mc_invEnergy_branch->GetEntry(index_);
			hyp_other_jets_mc_invEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_invEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_mc_otherEnergy()
	{
		if (not hyp_other_jets_mc_otherEnergy_isLoaded) {
			if (hyp_other_jets_mc_otherEnergy_branch != 0) 
				hyp_other_jets_mc_otherEnergy_branch->GetEntry(index_);
			hyp_other_jets_mc_otherEnergy_isLoaded = true;
		}
		return hyp_other_jets_mc_otherEnergy_;
	}
	vector<vector<float> > &hyp_other_jets_pat_bCorrF()
	{
		if (not hyp_other_jets_pat_bCorrF_isLoaded) {
			if (hyp_other_jets_pat_bCorrF_branch != 0) 
				hyp_other_jets_pat_bCorrF_branch->GetEntry(index_);
			hyp_other_jets_pat_bCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_bCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_cCorrF()
	{
		if (not hyp_other_jets_pat_cCorrF_isLoaded) {
			if (hyp_other_jets_pat_cCorrF_branch != 0) 
				hyp_other_jets_pat_cCorrF_branch->GetEntry(index_);
			hyp_other_jets_pat_cCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_cCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_gluCorrF()
	{
		if (not hyp_other_jets_pat_gluCorrF_isLoaded) {
			if (hyp_other_jets_pat_gluCorrF_branch != 0) 
				hyp_other_jets_pat_gluCorrF_branch->GetEntry(index_);
			hyp_other_jets_pat_gluCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_gluCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_jetCharge()
	{
		if (not hyp_other_jets_pat_jetCharge_isLoaded) {
			if (hyp_other_jets_pat_jetCharge_branch != 0) 
				hyp_other_jets_pat_jetCharge_branch->GetEntry(index_);
			hyp_other_jets_pat_jetCharge_isLoaded = true;
		}
		return hyp_other_jets_pat_jetCharge_;
	}
	vector<vector<float> > &hyp_other_jets_pat_noCorrF()
	{
		if (not hyp_other_jets_pat_noCorrF_isLoaded) {
			if (hyp_other_jets_pat_noCorrF_branch != 0) 
				hyp_other_jets_pat_noCorrF_branch->GetEntry(index_);
			hyp_other_jets_pat_noCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_noCorrF_;
	}
	vector<vector<float> > &hyp_other_jets_pat_udsCorrF()
	{
		if (not hyp_other_jets_pat_udsCorrF_isLoaded) {
			if (hyp_other_jets_pat_udsCorrF_branch != 0) 
				hyp_other_jets_pat_udsCorrF_branch->GetEntry(index_);
			hyp_other_jets_pat_udsCorrF_isLoaded = true;
		}
		return hyp_other_jets_pat_udsCorrF_;
	}
	int &evt_HLT1()
	{
		if (not evt_HLT1_isLoaded) {
			if (evt_HLT1_branch != 0) 
				evt_HLT1_branch->GetEntry(index_);
			evt_HLT1_isLoaded = true;
		}
		return evt_HLT1_;
	}
	int &evt_HLT2()
	{
		if (not evt_HLT2_isLoaded) {
			if (evt_HLT2_branch != 0) 
				evt_HLT2_branch->GetEntry(index_);
			evt_HLT2_isLoaded = true;
		}
		return evt_HLT2_;
	}
	int &evt_HLT3()
	{
		if (not evt_HLT3_isLoaded) {
			if (evt_HLT3_branch != 0) 
				evt_HLT3_branch->GetEntry(index_);
			evt_HLT3_isLoaded = true;
		}
		return evt_HLT3_;
	}
	int &evt_HLT4()
	{
		if (not evt_HLT4_isLoaded) {
			if (evt_HLT4_branch != 0) 
				evt_HLT4_branch->GetEntry(index_);
			evt_HLT4_isLoaded = true;
		}
		return evt_HLT4_;
	}
	int &evt_HLT5()
	{
		if (not evt_HLT5_isLoaded) {
			if (evt_HLT5_branch != 0) 
				evt_HLT5_branch->GetEntry(index_);
			evt_HLT5_isLoaded = true;
		}
		return evt_HLT5_;
	}
	int &evt_HLT6()
	{
		if (not evt_HLT6_isLoaded) {
			if (evt_HLT6_branch != 0) 
				evt_HLT6_branch->GetEntry(index_);
			evt_HLT6_isLoaded = true;
		}
		return evt_HLT6_;
	}
	int &evt_HLT7()
	{
		if (not evt_HLT7_isLoaded) {
			if (evt_HLT7_branch != 0) 
				evt_HLT7_branch->GetEntry(index_);
			evt_HLT7_isLoaded = true;
		}
		return evt_HLT7_;
	}
	int &evt_HLT8()
	{
		if (not evt_HLT8_isLoaded) {
			if (evt_HLT8_branch != 0) 
				evt_HLT8_branch->GetEntry(index_);
			evt_HLT8_isLoaded = true;
		}
		return evt_HLT8_;
	}
	int &evt_L1_1()
	{
		if (not evt_L1_1_isLoaded) {
			if (evt_L1_1_branch != 0) 
				evt_L1_1_branch->GetEntry(index_);
			evt_L1_1_isLoaded = true;
		}
		return evt_L1_1_;
	}
	int &evt_L1_2()
	{
		if (not evt_L1_2_isLoaded) {
			if (evt_L1_2_branch != 0) 
				evt_L1_2_branch->GetEntry(index_);
			evt_L1_2_isLoaded = true;
		}
		return evt_L1_2_;
	}
	int &evt_L1_3()
	{
		if (not evt_L1_3_isLoaded) {
			if (evt_L1_3_branch != 0) 
				evt_L1_3_branch->GetEntry(index_);
			evt_L1_3_isLoaded = true;
		}
		return evt_L1_3_;
	}
	int &evt_L1_4()
	{
		if (not evt_L1_4_isLoaded) {
			if (evt_L1_4_branch != 0) 
				evt_L1_4_branch->GetEntry(index_);
			evt_L1_4_isLoaded = true;
		}
		return evt_L1_4_;
	}
	vector<int> &els_mc_id()
	{
		if (not els_mc_id_isLoaded) {
			if (els_mc_id_branch != 0) 
				els_mc_id_branch->GetEntry(index_);
			els_mc_id_isLoaded = true;
		}
		return els_mc_id_;
	}
	vector<int> &els_mcidx()
	{
		if (not els_mcidx_isLoaded) {
			if (els_mcidx_branch != 0) 
				els_mcidx_branch->GetEntry(index_);
			els_mcidx_isLoaded = true;
		}
		return els_mcidx_;
	}
	vector<int> &els_mc_motherid()
	{
		if (not els_mc_motherid_isLoaded) {
			if (els_mc_motherid_branch != 0) 
				els_mc_motherid_branch->GetEntry(index_);
			els_mc_motherid_isLoaded = true;
		}
		return els_mc_motherid_;
	}
	vector<int> &jets_mc_id()
	{
		if (not jets_mc_id_isLoaded) {
			if (jets_mc_id_branch != 0) 
				jets_mc_id_branch->GetEntry(index_);
			jets_mc_id_isLoaded = true;
		}
		return jets_mc_id_;
	}
	vector<int> &mus_mc_id()
	{
		if (not mus_mc_id_isLoaded) {
			if (mus_mc_id_branch != 0) 
				mus_mc_id_branch->GetEntry(index_);
			mus_mc_id_isLoaded = true;
		}
		return mus_mc_id_;
	}
	vector<int> &mus_mcidx()
	{
		if (not mus_mcidx_isLoaded) {
			if (mus_mcidx_branch != 0) 
				mus_mcidx_branch->GetEntry(index_);
			mus_mcidx_isLoaded = true;
		}
		return mus_mcidx_;
	}
	vector<int> &mus_mc_motherid()
	{
		if (not mus_mc_motherid_isLoaded) {
			if (mus_mc_motherid_branch != 0) 
				mus_mc_motherid_branch->GetEntry(index_);
			mus_mc_motherid_isLoaded = true;
		}
		return mus_mc_motherid_;
	}
	vector<int> &trk_mc_id()
	{
		if (not trk_mc_id_isLoaded) {
			if (trk_mc_id_branch != 0) 
				trk_mc_id_branch->GetEntry(index_);
			trk_mc_id_isLoaded = true;
		}
		return trk_mc_id_;
	}
	vector<int> &trk_mcidx()
	{
		if (not trk_mcidx_isLoaded) {
			if (trk_mcidx_branch != 0) 
				trk_mcidx_branch->GetEntry(index_);
			trk_mcidx_isLoaded = true;
		}
		return trk_mcidx_;
	}
	vector<int> &trk_mc_motherid()
	{
		if (not trk_mc_motherid_isLoaded) {
			if (trk_mc_motherid_branch != 0) 
				trk_mc_motherid_branch->GetEntry(index_);
			trk_mc_motherid_isLoaded = true;
		}
		return trk_mc_motherid_;
	}
	vector<int> &els_closestMuon()
	{
		if (not els_closestMuon_isLoaded) {
			if (els_closestMuon_branch != 0) 
				els_closestMuon_branch->GetEntry(index_);
			els_closestMuon_isLoaded = true;
		}
		return els_closestMuon_;
	}
	vector<int> &els_trkidx()
	{
		if (not els_trkidx_isLoaded) {
			if (els_trkidx_branch != 0) 
				els_trkidx_branch->GetEntry(index_);
			els_trkidx_isLoaded = true;
		}
		return els_trkidx_;
	}
	vector<int> &els_category()
	{
		if (not els_category_isLoaded) {
			if (els_category_branch != 0) 
				els_category_branch->GetEntry(index_);
			els_category_isLoaded = true;
		}
		return els_category_;
	}
	vector<int> &els_categoryold()
	{
		if (not els_categoryold_isLoaded) {
			if (els_categoryold_branch != 0) 
				els_categoryold_branch->GetEntry(index_);
			els_categoryold_isLoaded = true;
		}
		return els_categoryold_;
	}
	vector<int> &els_charge()
	{
		if (not els_charge_isLoaded) {
			if (els_charge_branch != 0) 
				els_charge_branch->GetEntry(index_);
			els_charge_isLoaded = true;
		}
		return els_charge_;
	}
	vector<int> &els_class()
	{
		if (not els_class_isLoaded) {
			if (els_class_branch != 0) 
				els_class_branch->GetEntry(index_);
			els_class_isLoaded = true;
		}
		return els_class_;
	}
	vector<int> &els_looseId()
	{
		if (not els_looseId_isLoaded) {
			if (els_looseId_branch != 0) 
				els_looseId_branch->GetEntry(index_);
			els_looseId_isLoaded = true;
		}
		return els_looseId_;
	}
	vector<int> &els_lostHits()
	{
		if (not els_lostHits_isLoaded) {
			if (els_lostHits_branch != 0) 
				els_lostHits_branch->GetEntry(index_);
			els_lostHits_isLoaded = true;
		}
		return els_lostHits_;
	}
	vector<int> &els_nSeed()
	{
		if (not els_nSeed_isLoaded) {
			if (els_nSeed_branch != 0) 
				els_nSeed_branch->GetEntry(index_);
			els_nSeed_isLoaded = true;
		}
		return els_nSeed_;
	}
	vector<int> &els_pass3looseId()
	{
		if (not els_pass3looseId_isLoaded) {
			if (els_pass3looseId_branch != 0) 
				els_pass3looseId_branch->GetEntry(index_);
			els_pass3looseId_isLoaded = true;
		}
		return els_pass3looseId_;
	}
	vector<int> &els_pass3simpleId()
	{
		if (not els_pass3simpleId_isLoaded) {
			if (els_pass3simpleId_branch != 0) 
				els_pass3simpleId_branch->GetEntry(index_);
			els_pass3simpleId_isLoaded = true;
		}
		return els_pass3simpleId_;
	}
	vector<int> &els_pass3tightId()
	{
		if (not els_pass3tightId_isLoaded) {
			if (els_pass3tightId_branch != 0) 
				els_pass3tightId_branch->GetEntry(index_);
			els_pass3tightId_isLoaded = true;
		}
		return els_pass3tightId_;
	}
	vector<int> &els_robustId()
	{
		if (not els_robustId_isLoaded) {
			if (els_robustId_branch != 0) 
				els_robustId_branch->GetEntry(index_);
			els_robustId_isLoaded = true;
		}
		return els_robustId_;
	}
	vector<int> &els_simpleIdPlus()
	{
		if (not els_simpleIdPlus_isLoaded) {
			if (els_simpleIdPlus_branch != 0) 
				els_simpleIdPlus_branch->GetEntry(index_);
			els_simpleIdPlus_isLoaded = true;
		}
		return els_simpleIdPlus_;
	}
	vector<int> &els_tightId()
	{
		if (not els_tightId_isLoaded) {
			if (els_tightId_branch != 0) 
				els_tightId_branch->GetEntry(index_);
			els_tightId_isLoaded = true;
		}
		return els_tightId_;
	}
	vector<int> &els_validHits()
	{
		if (not els_validHits_isLoaded) {
			if (els_validHits_branch != 0) 
				els_validHits_branch->GetEntry(index_);
			els_validHits_isLoaded = true;
		}
		return els_validHits_;
	}
	vector<int> &genps_id()
	{
		if (not genps_id_isLoaded) {
			if (genps_id_branch != 0) 
				genps_id_branch->GetEntry(index_);
			genps_id_isLoaded = true;
		}
		return genps_id_;
	}
	vector<int> &genps_id_mother()
	{
		if (not genps_id_mother_isLoaded) {
			if (genps_id_mother_branch != 0) 
				genps_id_mother_branch->GetEntry(index_);
			genps_id_mother_isLoaded = true;
		}
		return genps_id_mother_;
	}
	vector<int> &genps_status()
	{
		if (not genps_status_isLoaded) {
			if (genps_status_branch != 0) 
				genps_status_branch->GetEntry(index_);
			genps_status_isLoaded = true;
		}
		return genps_status_;
	}
	vector<int> &hyp_ll_charge()
	{
		if (not hyp_ll_charge_isLoaded) {
			if (hyp_ll_charge_branch != 0) 
				hyp_ll_charge_branch->GetEntry(index_);
			hyp_ll_charge_isLoaded = true;
		}
		return hyp_ll_charge_;
	}
	vector<int> &hyp_ll_id()
	{
		if (not hyp_ll_id_isLoaded) {
			if (hyp_ll_id_branch != 0) 
				hyp_ll_id_branch->GetEntry(index_);
			hyp_ll_id_isLoaded = true;
		}
		return hyp_ll_id_;
	}
	vector<int> &hyp_ll_index()
	{
		if (not hyp_ll_index_isLoaded) {
			if (hyp_ll_index_branch != 0) 
				hyp_ll_index_branch->GetEntry(index_);
			hyp_ll_index_isLoaded = true;
		}
		return hyp_ll_index_;
	}
	vector<int> &hyp_ll_lostHits()
	{
		if (not hyp_ll_lostHits_isLoaded) {
			if (hyp_ll_lostHits_branch != 0) 
				hyp_ll_lostHits_branch->GetEntry(index_);
			hyp_ll_lostHits_isLoaded = true;
		}
		return hyp_ll_lostHits_;
	}
	vector<int> &hyp_ll_mc_id()
	{
		if (not hyp_ll_mc_id_isLoaded) {
			if (hyp_ll_mc_id_branch != 0) 
				hyp_ll_mc_id_branch->GetEntry(index_);
			hyp_ll_mc_id_isLoaded = true;
		}
		return hyp_ll_mc_id_;
	}
	vector<int> &hyp_ll_mc_motherid()
	{
		if (not hyp_ll_mc_motherid_isLoaded) {
			if (hyp_ll_mc_motherid_branch != 0) 
				hyp_ll_mc_motherid_branch->GetEntry(index_);
			hyp_ll_mc_motherid_isLoaded = true;
		}
		return hyp_ll_mc_motherid_;
	}
	vector<int> &hyp_ll_validHits()
	{
		if (not hyp_ll_validHits_isLoaded) {
			if (hyp_ll_validHits_branch != 0) 
				hyp_ll_validHits_branch->GetEntry(index_);
			hyp_ll_validHits_isLoaded = true;
		}
		return hyp_ll_validHits_;
	}
	vector<int> &hyp_lt_charge()
	{
		if (not hyp_lt_charge_isLoaded) {
			if (hyp_lt_charge_branch != 0) 
				hyp_lt_charge_branch->GetEntry(index_);
			hyp_lt_charge_isLoaded = true;
		}
		return hyp_lt_charge_;
	}
	vector<int> &hyp_lt_id()
	{
		if (not hyp_lt_id_isLoaded) {
			if (hyp_lt_id_branch != 0) 
				hyp_lt_id_branch->GetEntry(index_);
			hyp_lt_id_isLoaded = true;
		}
		return hyp_lt_id_;
	}
	vector<int> &hyp_lt_index()
	{
		if (not hyp_lt_index_isLoaded) {
			if (hyp_lt_index_branch != 0) 
				hyp_lt_index_branch->GetEntry(index_);
			hyp_lt_index_isLoaded = true;
		}
		return hyp_lt_index_;
	}
	vector<int> &hyp_lt_lostHits()
	{
		if (not hyp_lt_lostHits_isLoaded) {
			if (hyp_lt_lostHits_branch != 0) 
				hyp_lt_lostHits_branch->GetEntry(index_);
			hyp_lt_lostHits_isLoaded = true;
		}
		return hyp_lt_lostHits_;
	}
	vector<int> &hyp_lt_mc_id()
	{
		if (not hyp_lt_mc_id_isLoaded) {
			if (hyp_lt_mc_id_branch != 0) 
				hyp_lt_mc_id_branch->GetEntry(index_);
			hyp_lt_mc_id_isLoaded = true;
		}
		return hyp_lt_mc_id_;
	}
	vector<int> &hyp_lt_mc_motherid()
	{
		if (not hyp_lt_mc_motherid_isLoaded) {
			if (hyp_lt_mc_motherid_branch != 0) 
				hyp_lt_mc_motherid_branch->GetEntry(index_);
			hyp_lt_mc_motherid_isLoaded = true;
		}
		return hyp_lt_mc_motherid_;
	}
	vector<int> &hyp_lt_validHits()
	{
		if (not hyp_lt_validHits_isLoaded) {
			if (hyp_lt_validHits_branch != 0) 
				hyp_lt_validHits_branch->GetEntry(index_);
			hyp_lt_validHits_isLoaded = true;
		}
		return hyp_lt_validHits_;
	}
	vector<int> &hyp_njets()
	{
		if (not hyp_njets_isLoaded) {
			if (hyp_njets_branch != 0) 
				hyp_njets_branch->GetEntry(index_);
			hyp_njets_isLoaded = true;
		}
		return hyp_njets_;
	}
	vector<int> &hyp_nojets()
	{
		if (not hyp_nojets_isLoaded) {
			if (hyp_nojets_branch != 0) 
				hyp_nojets_branch->GetEntry(index_);
			hyp_nojets_isLoaded = true;
		}
		return hyp_nojets_;
	}
	vector<int> &hyp_type()
	{
		if (not hyp_type_isLoaded) {
			if (hyp_type_branch != 0) 
				hyp_type_branch->GetEntry(index_);
			hyp_type_isLoaded = true;
		}
		return hyp_type_;
	}
	vector<int> &hyp_quadlep_first_type()
	{
		if (not hyp_quadlep_first_type_isLoaded) {
			if (hyp_quadlep_first_type_branch != 0) 
				hyp_quadlep_first_type_branch->GetEntry(index_);
			hyp_quadlep_first_type_isLoaded = true;
		}
		return hyp_quadlep_first_type_;
	}
	vector<int> &hyp_quadlep_fourth_type()
	{
		if (not hyp_quadlep_fourth_type_isLoaded) {
			if (hyp_quadlep_fourth_type_branch != 0) 
				hyp_quadlep_fourth_type_branch->GetEntry(index_);
			hyp_quadlep_fourth_type_isLoaded = true;
		}
		return hyp_quadlep_fourth_type_;
	}
	vector<int> &hyp_quadlep_second_type()
	{
		if (not hyp_quadlep_second_type_isLoaded) {
			if (hyp_quadlep_second_type_branch != 0) 
				hyp_quadlep_second_type_branch->GetEntry(index_);
			hyp_quadlep_second_type_isLoaded = true;
		}
		return hyp_quadlep_second_type_;
	}
	vector<int> &hyp_quadlep_third_type()
	{
		if (not hyp_quadlep_third_type_isLoaded) {
			if (hyp_quadlep_third_type_branch != 0) 
				hyp_quadlep_third_type_branch->GetEntry(index_);
			hyp_quadlep_third_type_isLoaded = true;
		}
		return hyp_quadlep_third_type_;
	}
	vector<int> &hyp_trilep_first_type()
	{
		if (not hyp_trilep_first_type_isLoaded) {
			if (hyp_trilep_first_type_branch != 0) 
				hyp_trilep_first_type_branch->GetEntry(index_);
			hyp_trilep_first_type_isLoaded = true;
		}
		return hyp_trilep_first_type_;
	}
	vector<int> &hyp_trilep_second_type()
	{
		if (not hyp_trilep_second_type_isLoaded) {
			if (hyp_trilep_second_type_branch != 0) 
				hyp_trilep_second_type_branch->GetEntry(index_);
			hyp_trilep_second_type_isLoaded = true;
		}
		return hyp_trilep_second_type_;
	}
	vector<int> &hyp_trilep_third_type()
	{
		if (not hyp_trilep_third_type_isLoaded) {
			if (hyp_trilep_third_type_branch != 0) 
				hyp_trilep_third_type_branch->GetEntry(index_);
			hyp_trilep_third_type_isLoaded = true;
		}
		return hyp_trilep_third_type_;
	}
	vector<int> &jets_closestElectron()
	{
		if (not jets_closestElectron_isLoaded) {
			if (jets_closestElectron_branch != 0) 
				jets_closestElectron_branch->GetEntry(index_);
			jets_closestElectron_isLoaded = true;
		}
		return jets_closestElectron_;
	}
	vector<int> &jets_closestMuon()
	{
		if (not jets_closestMuon_isLoaded) {
			if (jets_closestMuon_branch != 0) 
				jets_closestMuon_branch->GetEntry(index_);
			jets_closestMuon_isLoaded = true;
		}
		return jets_closestMuon_;
	}
	vector<int> &mus_closestEle()
	{
		if (not mus_closestEle_isLoaded) {
			if (mus_closestEle_branch != 0) 
				mus_closestEle_branch->GetEntry(index_);
			mus_closestEle_isLoaded = true;
		}
		return mus_closestEle_;
	}
	vector<int> &mus_closestJet()
	{
		if (not mus_closestJet_isLoaded) {
			if (mus_closestJet_branch != 0) 
				mus_closestJet_branch->GetEntry(index_);
			mus_closestJet_isLoaded = true;
		}
		return mus_closestJet_;
	}
	vector<int> &mus_trkidx()
	{
		if (not mus_trkidx_isLoaded) {
			if (mus_trkidx_branch != 0) 
				mus_trkidx_branch->GetEntry(index_);
			mus_trkidx_isLoaded = true;
		}
		return mus_trkidx_;
	}
	vector<int> &mus_charge()
	{
		if (not mus_charge_isLoaded) {
			if (mus_charge_branch != 0) 
				mus_charge_branch->GetEntry(index_);
			mus_charge_isLoaded = true;
		}
		return mus_charge_;
	}
	vector<int> &mus_gfit_validHits()
	{
		if (not mus_gfit_validHits_isLoaded) {
			if (mus_gfit_validHits_branch != 0) 
				mus_gfit_validHits_branch->GetEntry(index_);
			mus_gfit_validHits_isLoaded = true;
		}
		return mus_gfit_validHits_;
	}
	vector<int> &mus_goodmask()
	{
		if (not mus_goodmask_isLoaded) {
			if (mus_goodmask_branch != 0) 
				mus_goodmask_branch->GetEntry(index_);
			mus_goodmask_isLoaded = true;
		}
		return mus_goodmask_;
	}
	vector<int> &mus_iso03_ntrk()
	{
		if (not mus_iso03_ntrk_isLoaded) {
			if (mus_iso03_ntrk_branch != 0) 
				mus_iso03_ntrk_branch->GetEntry(index_);
			mus_iso03_ntrk_isLoaded = true;
		}
		return mus_iso03_ntrk_;
	}
	vector<int> &mus_iso05_ntrk()
	{
		if (not mus_iso05_ntrk_isLoaded) {
			if (mus_iso05_ntrk_branch != 0) 
				mus_iso05_ntrk_branch->GetEntry(index_);
			mus_iso05_ntrk_isLoaded = true;
		}
		return mus_iso05_ntrk_;
	}
	vector<int> &mus_lostHits()
	{
		if (not mus_lostHits_isLoaded) {
			if (mus_lostHits_branch != 0) 
				mus_lostHits_branch->GetEntry(index_);
			mus_lostHits_isLoaded = true;
		}
		return mus_lostHits_;
	}
	vector<int> &mus_nmatches()
	{
		if (not mus_nmatches_isLoaded) {
			if (mus_nmatches_branch != 0) 
				mus_nmatches_branch->GetEntry(index_);
			mus_nmatches_isLoaded = true;
		}
		return mus_nmatches_;
	}
	vector<int> &mus_pid_TM2DCompatibilityLoose()
	{
		if (not mus_pid_TM2DCompatibilityLoose_isLoaded) {
			if (mus_pid_TM2DCompatibilityLoose_branch != 0) 
				mus_pid_TM2DCompatibilityLoose_branch->GetEntry(index_);
			mus_pid_TM2DCompatibilityLoose_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityLoose_;
	}
	vector<int> &mus_pid_TM2DCompatibilityTight()
	{
		if (not mus_pid_TM2DCompatibilityTight_isLoaded) {
			if (mus_pid_TM2DCompatibilityTight_branch != 0) 
				mus_pid_TM2DCompatibilityTight_branch->GetEntry(index_);
			mus_pid_TM2DCompatibilityTight_isLoaded = true;
		}
		return mus_pid_TM2DCompatibilityTight_;
	}
	vector<int> &mus_pid_TMLastStationLoose()
	{
		if (not mus_pid_TMLastStationLoose_isLoaded) {
			if (mus_pid_TMLastStationLoose_branch != 0) 
				mus_pid_TMLastStationLoose_branch->GetEntry(index_);
			mus_pid_TMLastStationLoose_isLoaded = true;
		}
		return mus_pid_TMLastStationLoose_;
	}
	vector<int> &mus_pid_TMLastStationTight()
	{
		if (not mus_pid_TMLastStationTight_isLoaded) {
			if (mus_pid_TMLastStationTight_branch != 0) 
				mus_pid_TMLastStationTight_branch->GetEntry(index_);
			mus_pid_TMLastStationTight_isLoaded = true;
		}
		return mus_pid_TMLastStationTight_;
	}
	vector<int> &mus_trk_charge()
	{
		if (not mus_trk_charge_isLoaded) {
			if (mus_trk_charge_branch != 0) 
				mus_trk_charge_branch->GetEntry(index_);
			mus_trk_charge_isLoaded = true;
		}
		return mus_trk_charge_;
	}
	vector<int> &mus_trkrefkey()
	{
		if (not mus_trkrefkey_isLoaded) {
			if (mus_trkrefkey_branch != 0) 
				mus_trkrefkey_branch->GetEntry(index_);
			mus_trkrefkey_isLoaded = true;
		}
		return mus_trkrefkey_;
	}
	vector<int> &mus_type()
	{
		if (not mus_type_isLoaded) {
			if (mus_type_branch != 0) 
				mus_type_branch->GetEntry(index_);
			mus_type_isLoaded = true;
		}
		return mus_type_;
	}
	vector<int> &mus_validHits()
	{
		if (not mus_validHits_isLoaded) {
			if (mus_validHits_branch != 0) 
				mus_validHits_branch->GetEntry(index_);
			mus_validHits_isLoaded = true;
		}
		return mus_validHits_;
	}
	vector<int> &els_pat_genID()
	{
		if (not els_pat_genID_isLoaded) {
			if (els_pat_genID_branch != 0) 
				els_pat_genID_branch->GetEntry(index_);
			els_pat_genID_isLoaded = true;
		}
		return els_pat_genID_;
	}
	vector<int> &els_pat_genMotherID()
	{
		if (not els_pat_genMotherID_isLoaded) {
			if (els_pat_genMotherID_branch != 0) 
				els_pat_genMotherID_branch->GetEntry(index_);
			els_pat_genMotherID_isLoaded = true;
		}
		return els_pat_genMotherID_;
	}
	vector<int> &jets_pat_genPartonMother_id()
	{
		if (not jets_pat_genPartonMother_id_isLoaded) {
			if (jets_pat_genPartonMother_id_branch != 0) 
				jets_pat_genPartonMother_id_branch->GetEntry(index_);
			jets_pat_genPartonMother_id_isLoaded = true;
		}
		return jets_pat_genPartonMother_id_;
	}
	vector<int> &jets_pat_genParton_id()
	{
		if (not jets_pat_genParton_id_isLoaded) {
			if (jets_pat_genParton_id_branch != 0) 
				jets_pat_genParton_id_branch->GetEntry(index_);
			jets_pat_genParton_id_isLoaded = true;
		}
		return jets_pat_genParton_id_;
	}
	vector<int> &jets_pat_partonFlavour()
	{
		if (not jets_pat_partonFlavour_isLoaded) {
			if (jets_pat_partonFlavour_branch != 0) 
				jets_pat_partonFlavour_branch->GetEntry(index_);
			jets_pat_partonFlavour_isLoaded = true;
		}
		return jets_pat_partonFlavour_;
	}
	vector<int> &mus_pat_genID()
	{
		if (not mus_pat_genID_isLoaded) {
			if (mus_pat_genID_branch != 0) 
				mus_pat_genID_branch->GetEntry(index_);
			mus_pat_genID_isLoaded = true;
		}
		return mus_pat_genID_;
	}
	vector<int> &mus_pat_genMotherID()
	{
		if (not mus_pat_genMotherID_isLoaded) {
			if (mus_pat_genMotherID_branch != 0) 
				mus_pat_genMotherID_branch->GetEntry(index_);
			mus_pat_genMotherID_isLoaded = true;
		}
		return mus_pat_genMotherID_;
	}
	vector<int> &trks_charge()
	{
		if (not trks_charge_isLoaded) {
			if (trks_charge_branch != 0) 
				trks_charge_branch->GetEntry(index_);
			trks_charge_isLoaded = true;
		}
		return trks_charge_;
	}
	vector<int> &trks_lostHits()
	{
		if (not trks_lostHits_isLoaded) {
			if (trks_lostHits_branch != 0) 
				trks_lostHits_branch->GetEntry(index_);
			trks_lostHits_isLoaded = true;
		}
		return trks_lostHits_;
	}
	vector<int> &trks_validHits()
	{
		if (not trks_validHits_isLoaded) {
			if (trks_validHits_branch != 0) 
				trks_validHits_branch->GetEntry(index_);
			trks_validHits_isLoaded = true;
		}
		return trks_validHits_;
	}
	vector<int> &trks_elsidx()
	{
		if (not trks_elsidx_isLoaded) {
			if (trks_elsidx_branch != 0) 
				trks_elsidx_branch->GetEntry(index_);
			trks_elsidx_isLoaded = true;
		}
		return trks_elsidx_;
	}
	vector<int> &trk_musidx()
	{
		if (not trk_musidx_isLoaded) {
			if (trk_musidx_branch != 0) 
				trk_musidx_branch->GetEntry(index_);
			trk_musidx_isLoaded = true;
		}
		return trk_musidx_;
	}
	vector<vector<int> > &hyp_jets_mc_id()
	{
		if (not hyp_jets_mc_id_isLoaded) {
			if (hyp_jets_mc_id_branch != 0) 
				hyp_jets_mc_id_branch->GetEntry(index_);
			hyp_jets_mc_id_isLoaded = true;
		}
		return hyp_jets_mc_id_;
	}
	vector<vector<int> > &hyp_jets_pat_genPartonMother_id()
	{
		if (not hyp_jets_pat_genPartonMother_id_isLoaded) {
			if (hyp_jets_pat_genPartonMother_id_branch != 0) 
				hyp_jets_pat_genPartonMother_id_branch->GetEntry(index_);
			hyp_jets_pat_genPartonMother_id_isLoaded = true;
		}
		return hyp_jets_pat_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_jets_pat_genParton_id()
	{
		if (not hyp_jets_pat_genParton_id_isLoaded) {
			if (hyp_jets_pat_genParton_id_branch != 0) 
				hyp_jets_pat_genParton_id_branch->GetEntry(index_);
			hyp_jets_pat_genParton_id_isLoaded = true;
		}
		return hyp_jets_pat_genParton_id_;
	}
	vector<vector<int> > &hyp_jets_pat_partonFlavour()
	{
		if (not hyp_jets_pat_partonFlavour_isLoaded) {
			if (hyp_jets_pat_partonFlavour_branch != 0) 
				hyp_jets_pat_partonFlavour_branch->GetEntry(index_);
			hyp_jets_pat_partonFlavour_isLoaded = true;
		}
		return hyp_jets_pat_partonFlavour_;
	}
	vector<vector<int> > &hyp_other_jets_mc_id()
	{
		if (not hyp_other_jets_mc_id_isLoaded) {
			if (hyp_other_jets_mc_id_branch != 0) 
				hyp_other_jets_mc_id_branch->GetEntry(index_);
			hyp_other_jets_mc_id_isLoaded = true;
		}
		return hyp_other_jets_mc_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_genPartonMother_id()
	{
		if (not hyp_other_jets_pat_genPartonMother_id_isLoaded) {
			if (hyp_other_jets_pat_genPartonMother_id_branch != 0) 
				hyp_other_jets_pat_genPartonMother_id_branch->GetEntry(index_);
			hyp_other_jets_pat_genPartonMother_id_isLoaded = true;
		}
		return hyp_other_jets_pat_genPartonMother_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_genParton_id()
	{
		if (not hyp_other_jets_pat_genParton_id_isLoaded) {
			if (hyp_other_jets_pat_genParton_id_branch != 0) 
				hyp_other_jets_pat_genParton_id_branch->GetEntry(index_);
			hyp_other_jets_pat_genParton_id_isLoaded = true;
		}
		return hyp_other_jets_pat_genParton_id_;
	}
	vector<vector<int> > &hyp_other_jets_pat_partonFlavour()
	{
		if (not hyp_other_jets_pat_partonFlavour_isLoaded) {
			if (hyp_other_jets_pat_partonFlavour_branch != 0) 
				hyp_other_jets_pat_partonFlavour_branch->GetEntry(index_);
			hyp_other_jets_pat_partonFlavour_isLoaded = true;
		}
		return hyp_other_jets_pat_partonFlavour_;
	}
	vector<vector<int> > &hyp_quadlep_jets_index()
	{
		if (not hyp_quadlep_jets_index_isLoaded) {
			if (hyp_quadlep_jets_index_branch != 0) 
				hyp_quadlep_jets_index_branch->GetEntry(index_);
			hyp_quadlep_jets_index_isLoaded = true;
		}
		return hyp_quadlep_jets_index_;
	}
	vector<vector<int> > &hyp_trilep_jets_index()
	{
		if (not hyp_trilep_jets_index_isLoaded) {
			if (hyp_trilep_jets_index_branch != 0) 
				hyp_trilep_jets_index_branch->GetEntry(index_);
			hyp_trilep_jets_index_isLoaded = true;
		}
		return hyp_trilep_jets_index_;
	}
	unsigned int &evt_nels()
	{
		if (not evt_nels_isLoaded) {
			if (evt_nels_branch != 0) 
				evt_nels_branch->GetEntry(index_);
			evt_nels_isLoaded = true;
		}
		return evt_nels_;
	}
	unsigned int &evt_event()
	{
		if (not evt_event_isLoaded) {
			if (evt_event_branch != 0) 
				evt_event_branch->GetEntry(index_);
			evt_event_isLoaded = true;
		}
		return evt_event_;
	}
	unsigned int &evt_run()
	{
		if (not evt_run_isLoaded) {
			if (evt_run_branch != 0) 
				evt_run_branch->GetEntry(index_);
			evt_run_isLoaded = true;
		}
		return evt_run_;
	}
	unsigned int &evt_njets()
	{
		if (not evt_njets_isLoaded) {
			if (evt_njets_branch != 0) 
				evt_njets_branch->GetEntry(index_);
			evt_njets_isLoaded = true;
		}
		return evt_njets_;
	}
	vector<unsigned int> &hyp_quadlep_bucket()
	{
		if (not hyp_quadlep_bucket_isLoaded) {
			if (hyp_quadlep_bucket_branch != 0) 
				hyp_quadlep_bucket_branch->GetEntry(index_);
			hyp_quadlep_bucket_isLoaded = true;
		}
		return hyp_quadlep_bucket_;
	}
	vector<unsigned int> &hyp_quadlep_first_index()
	{
		if (not hyp_quadlep_first_index_isLoaded) {
			if (hyp_quadlep_first_index_branch != 0) 
				hyp_quadlep_first_index_branch->GetEntry(index_);
			hyp_quadlep_first_index_isLoaded = true;
		}
		return hyp_quadlep_first_index_;
	}
	vector<unsigned int> &hyp_quadlep_fourth_index()
	{
		if (not hyp_quadlep_fourth_index_isLoaded) {
			if (hyp_quadlep_fourth_index_branch != 0) 
				hyp_quadlep_fourth_index_branch->GetEntry(index_);
			hyp_quadlep_fourth_index_isLoaded = true;
		}
		return hyp_quadlep_fourth_index_;
	}
	vector<unsigned int> &hyp_quadlep_second_index()
	{
		if (not hyp_quadlep_second_index_isLoaded) {
			if (hyp_quadlep_second_index_branch != 0) 
				hyp_quadlep_second_index_branch->GetEntry(index_);
			hyp_quadlep_second_index_isLoaded = true;
		}
		return hyp_quadlep_second_index_;
	}
	vector<unsigned int> &hyp_quadlep_third_index()
	{
		if (not hyp_quadlep_third_index_isLoaded) {
			if (hyp_quadlep_third_index_branch != 0) 
				hyp_quadlep_third_index_branch->GetEntry(index_);
			hyp_quadlep_third_index_isLoaded = true;
		}
		return hyp_quadlep_third_index_;
	}
	vector<unsigned int> &hyp_trilep_bucket()
	{
		if (not hyp_trilep_bucket_isLoaded) {
			if (hyp_trilep_bucket_branch != 0) 
				hyp_trilep_bucket_branch->GetEntry(index_);
			hyp_trilep_bucket_isLoaded = true;
		}
		return hyp_trilep_bucket_;
	}
	vector<unsigned int> &hyp_trilep_first_index()
	{
		if (not hyp_trilep_first_index_isLoaded) {
			if (hyp_trilep_first_index_branch != 0) 
				hyp_trilep_first_index_branch->GetEntry(index_);
			hyp_trilep_first_index_isLoaded = true;
		}
		return hyp_trilep_first_index_;
	}
	vector<unsigned int> &hyp_trilep_second_index()
	{
		if (not hyp_trilep_second_index_isLoaded) {
			if (hyp_trilep_second_index_branch != 0) 
				hyp_trilep_second_index_branch->GetEntry(index_);
			hyp_trilep_second_index_isLoaded = true;
		}
		return hyp_trilep_second_index_;
	}
	vector<unsigned int> &hyp_trilep_third_index()
	{
		if (not hyp_trilep_third_index_isLoaded) {
			if (hyp_trilep_third_index_branch != 0) 
				hyp_trilep_third_index_branch->GetEntry(index_);
			hyp_trilep_third_index_isLoaded = true;
		}
		return hyp_trilep_third_index_;
	}
	vector<unsigned int> &els_pat_flag()
	{
		if (not els_pat_flag_isLoaded) {
			if (els_pat_flag_branch != 0) 
				els_pat_flag_branch->GetEntry(index_);
			els_pat_flag_isLoaded = true;
		}
		return els_pat_flag_;
	}
	vector<unsigned int> &jets_pat_flag()
	{
		if (not jets_pat_flag_isLoaded) {
			if (jets_pat_flag_branch != 0) 
				jets_pat_flag_branch->GetEntry(index_);
			jets_pat_flag_isLoaded = true;
		}
		return jets_pat_flag_;
	}
	vector<unsigned int> &mus_pat_flag()
	{
		if (not mus_pat_flag_isLoaded) {
			if (mus_pat_flag_branch != 0) 
				mus_pat_flag_branch->GetEntry(index_);
			mus_pat_flag_isLoaded = true;
		}
		return mus_pat_flag_;
	}
