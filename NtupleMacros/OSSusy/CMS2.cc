#include "CMS2.h"
CMS2 cms2;
namespace tas {
	TString &evt_CMS2tag() { return cms2.evt_CMS2tag(); }
	TString &evt_dataset() { return cms2.evt_dataset(); }
	vector<TString> &hlt1e31_trigNames() { return cms2.hlt1e31_trigNames(); }
	vector<TString> &hlt8e29_trigNames() { return cms2.hlt8e29_trigNames(); }
	vector<TString> &hlt_trigNames() { return cms2.hlt_trigNames(); }
	vector<TString> &l1_techtrigNames() { return cms2.l1_techtrigNames(); }
	vector<TString> &l1_trigNames() { return cms2.l1_trigNames(); }
	bool &evt_eventHasHalo() { return cms2.evt_eventHasHalo(); }
	vector<bool> &photons_haspixelSeed() { return cms2.photons_haspixelSeed(); }
	vector<double> &jets_closestElectron_DR() { return cms2.jets_closestElectron_DR(); }
	vector<double> &jets_closestMuon_DR() { return cms2.jets_closestMuon_DR(); }
	float &evt_bs_Xwidth() { return cms2.evt_bs_Xwidth(); }
	float &evt_bs_XwidthErr() { return cms2.evt_bs_XwidthErr(); }
	float &evt_bs_Ywidth() { return cms2.evt_bs_Ywidth(); }
	float &evt_bs_YwidthErr() { return cms2.evt_bs_YwidthErr(); }
	float &evt_bs_dxdz() { return cms2.evt_bs_dxdz(); }
	float &evt_bs_dxdzErr() { return cms2.evt_bs_dxdzErr(); }
	float &evt_bs_dydz() { return cms2.evt_bs_dydz(); }
	float &evt_bs_dydzErr() { return cms2.evt_bs_dydzErr(); }
	float &evt_bs_sigmaZ() { return cms2.evt_bs_sigmaZ(); }
	float &evt_bs_sigmaZErr() { return cms2.evt_bs_sigmaZErr(); }
	float &evt_bs_xErr() { return cms2.evt_bs_xErr(); }
	float &evt_bs_yErr() { return cms2.evt_bs_yErr(); }
	float &evt_bs_zErr() { return cms2.evt_bs_zErr(); }
	float &evt_bField() { return cms2.evt_bField(); }
	float &evt_kfactor() { return cms2.evt_kfactor(); }
	float &evt_scale1fb() { return cms2.evt_scale1fb(); }
	float &evt_xsec_excl() { return cms2.evt_xsec_excl(); }
	float &evt_xsec_incl() { return cms2.evt_xsec_incl(); }
	float &gen_met() { return cms2.gen_met(); }
	float &gen_metPhi() { return cms2.gen_metPhi(); }
	float &genps_alphaQCD() { return cms2.genps_alphaQCD(); }
	float &genps_pthat() { return cms2.genps_pthat(); }
	float &genps_qScale() { return cms2.genps_qScale(); }
	float &genps_weight() { return cms2.genps_weight(); }
	float &gen_sumEt() { return cms2.gen_sumEt(); }
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
	float &evt35x_met() { return cms2.evt35x_met(); }
	float &evt35x_metPhi() { return cms2.evt35x_metPhi(); }
	float &evt35x_metSig() { return cms2.evt35x_metSig(); }
	float &evt35x_sumet() { return cms2.evt35x_sumet(); }
	float &evt_ecalendcapm_met() { return cms2.evt_ecalendcapm_met(); }
	float &evt_ecalendcapm_metPhi() { return cms2.evt_ecalendcapm_metPhi(); }
	float &evt_ecalendcapp_met() { return cms2.evt_ecalendcapp_met(); }
	float &evt_ecalendcapp_metPhi() { return cms2.evt_ecalendcapp_metPhi(); }
	float &evt_ecalmet() { return cms2.evt_ecalmet(); }
	float &evt_ecalmetPhi() { return cms2.evt_ecalmetPhi(); }
	float &evt_endcapm_met() { return cms2.evt_endcapm_met(); }
	float &evt_endcapm_metPhi() { return cms2.evt_endcapm_metPhi(); }
	float &evt_endcapp_met() { return cms2.evt_endcapp_met(); }
	float &evt_endcapp_metPhi() { return cms2.evt_endcapp_metPhi(); }
	float &evt_hbheFilter() { return cms2.evt_hbheFilter(); }
	float &evt_hcalendcapm_met() { return cms2.evt_hcalendcapm_met(); }
	float &evt_hcalendcapm_metPhi() { return cms2.evt_hcalendcapm_metPhi(); }
	float &evt_hcalendcapp_met() { return cms2.evt_hcalendcapp_met(); }
	float &evt_hcalendcapp_metPhi() { return cms2.evt_hcalendcapp_metPhi(); }
	float &evt_hcalmet() { return cms2.evt_hcalmet(); }
	float &evt_hcalmetPhi() { return cms2.evt_hcalmetPhi(); }
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
	float &evt_pfmetSig() { return cms2.evt_pfmetSig(); }
	float &evt_pfsumet() { return cms2.evt_pfsumet(); }
	float &evt_tcmet() { return cms2.evt_tcmet(); }
	float &evt_tcmetPhi() { return cms2.evt_tcmetPhi(); }
	float &evt_tcmetSig() { return cms2.evt_tcmetSig(); }
	float &evt_tcsumet() { return cms2.evt_tcsumet(); }
	float &evt35X_tcmet() { return cms2.evt35X_tcmet(); }
	float &evt35X_tcmetPhi() { return cms2.evt35X_tcmetPhi(); }
	float &evt35X_tcmetSig() { return cms2.evt35X_tcmetSig(); }
	float &evt35X_tcsumet() { return cms2.evt35X_tcsumet(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &evt_bsp4() { return cms2.evt_bsp4(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_met_p4() { return cms2.l1_met_p4(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> >  &l1_mht_p4() { return cms2.l1_mht_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_motherp4() { return cms2.els_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_mc_p4() { return cms2.els_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_gp_p4() { return cms2.jets_mc_gp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_motherp4() { return cms2.jets_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_mc_p4() { return cms2.jets_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_motherp4() { return cms2.mus_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_mc_p4() { return cms2.mus_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_gp_p4() { return cms2.pfjets_mc_gp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_motherp4() { return cms2.pfjets_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_mc_p4() { return cms2.pfjets_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_motherp4() { return cms2.photons_mc_motherp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_mc_p4() { return cms2.photons_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trk_mcp4() { return cms2.trk_mcp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_conv_pos_p4() { return cms2.els_conv_pos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_inner_position() { return cms2.els_inner_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_outer_position() { return cms2.els_outer_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4() { return cms2.els_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4In() { return cms2.els_p4In(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_p4Out() { return cms2.els_p4Out(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_trk_p4() { return cms2.els_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_vertex_p4() { return cms2.els_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genjets_p4() { return cms2.genjets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_p4() { return cms2.genps_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &genps_prod_vtx() { return cms2.genps_prod_vtx(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_inner_position() { return cms2.gsftrks_inner_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_p4() { return cms2.gsftrks_outer_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_outer_position() { return cms2.gsftrks_outer_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_p4() { return cms2.gsftrks_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &gsftrks_vertex_p4() { return cms2.gsftrks_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_p4() { return cms2.hyp_ll_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_trk_p4() { return cms2.hyp_ll_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_p4() { return cms2.hyp_lt_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_trk_p4() { return cms2.hyp_lt_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_p4() { return cms2.hyp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_p4() { return cms2.hyp_FVFit_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_FVFit_v4() { return cms2.hyp_FVFit_v4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_ll_mc_p4() { return cms2.hyp_ll_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &hyp_lt_mc_p4() { return cms2.hyp_lt_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_p4() { return cms2.jets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_vertex_p4() { return cms2.jets_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jpts_p4() { return cms2.jpts_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emiso_p4() { return cms2.l1_emiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_emnoiso_p4() { return cms2.l1_emnoiso_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsc_p4() { return cms2.l1_jetsc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetsf_p4() { return cms2.l1_jetsf_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_jetst_p4() { return cms2.l1_jetst_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &l1_mus_p4() { return cms2.l1_mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_ecalpos_p4() { return cms2.mus_ecalpos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitdefault_p4() { return cms2.mus_fitdefault_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitfirsthit_p4() { return cms2.mus_fitfirsthit_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fitpicky_p4() { return cms2.mus_fitpicky_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_fittev_p4() { return cms2.mus_fittev_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_outerPos_p4() { return cms2.mus_gfit_outerPos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_p4() { return cms2.mus_gfit_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_gfit_vertex_p4() { return cms2.mus_gfit_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_p4() { return cms2.mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_p4() { return cms2.mus_sta_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_sta_vertex_p4() { return cms2.mus_sta_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_trk_p4() { return cms2.mus_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_vertex_p4() { return cms2.mus_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genMotherP4() { return cms2.els_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_genP4() { return cms2.els_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &els_pat_p4() { return cms2.els_pat_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genJet_p4() { return cms2.jets_pat_genJet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genPartonMother_p4() { return cms2.jets_pat_genPartonMother_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_genParton_p4() { return cms2.jets_pat_genParton_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_p4() { return cms2.jets_pat_jet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &jets_pat_jet_uncorp4() { return cms2.jets_pat_jet_uncorp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genMotherP4() { return cms2.mus_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_genP4() { return cms2.mus_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &mus_pat_p4() { return cms2.mus_pat_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_p4() { return cms2.pfels_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfels_posAtEcal_p4() { return cms2.pfels_posAtEcal_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfjets_p4() { return cms2.pfjets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_p4() { return cms2.pfmus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &pfmus_posAtEcal_p4() { return cms2.pfmus_posAtEcal_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &photons_p4() { return cms2.photons_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_p4() { return cms2.scs_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_pos_p4() { return cms2.scs_pos_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &scs_vtx_p4() { return cms2.scs_vtx_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_inner_position() { return cms2.trks_inner_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_p4() { return cms2.trks_outer_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_outer_position() { return cms2.trks_outer_position(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_trk_p4() { return cms2.trks_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trks_vertex_p4() { return cms2.trks_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &trkjets_p4() { return cms2.trkjets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > &vtxs_position() { return cms2.vtxs_position(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &genps_lepdaughter_p4() { return cms2.genps_lepdaughter_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt1e31_trigObjs_p4() { return cms2.hlt1e31_trigObjs_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt8e29_trigObjs_p4() { return cms2.hlt8e29_trigObjs_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hlt_trigObjs_p4() { return cms2.hlt_trigObjs_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_jets_p4() { return cms2.hyp_jets_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > > &hyp_other_jets_p4() { return cms2.hyp_other_jets_p4(); }
	vector<float> &jpts_combinedSecondaryVertexBJetTag() { return cms2.jpts_combinedSecondaryVertexBJetTag(); }
	vector<float> &jpts_combinedSecondaryVertexMVABJetTag() { return cms2.jpts_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &jpts_jetBProbabilityBJetTag() { return cms2.jpts_jetBProbabilityBJetTag(); }
	vector<float> &jpts_jetProbabilityBJetTag() { return cms2.jpts_jetProbabilityBJetTag(); }
	vector<float> &jpts_simpleSecondaryVertexHighEffBJetTag() { return cms2.jpts_simpleSecondaryVertexHighEffBJetTag(); }
	vector<float> &jpts_simpleSecondaryVertexHighPurBJetTags() { return cms2.jpts_simpleSecondaryVertexHighPurBJetTags(); }
	vector<float> &jpts_softElectronByIP3dBJetTag() { return cms2.jpts_softElectronByIP3dBJetTag(); }
	vector<float> &jpts_softElectronByPtBJetTag() { return cms2.jpts_softElectronByPtBJetTag(); }
	vector<float> &jpts_softMuonBJetTag() { return cms2.jpts_softMuonBJetTag(); }
	vector<float> &jpts_softMuonByIP3dBJetTag() { return cms2.jpts_softMuonByIP3dBJetTag(); }
	vector<float> &jpts_softMuonByPtBJetTag() { return cms2.jpts_softMuonByPtBJetTag(); }
	vector<float> &jpts_trackCountingHighEffBJetTag() { return cms2.jpts_trackCountingHighEffBJetTag(); }
	vector<float> &jpts_trackCountingHighPurBJetTag() { return cms2.jpts_trackCountingHighPurBJetTag(); }
	vector<float> &jets_combinedSecondaryVertexBJetTag() { return cms2.jets_combinedSecondaryVertexBJetTag(); }
	vector<float> &jets_combinedSecondaryVertexMVABJetTag() { return cms2.jets_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &jets_jetBProbabilityBJetTag() { return cms2.jets_jetBProbabilityBJetTag(); }
	vector<float> &jets_jetProbabilityBJetTag() { return cms2.jets_jetProbabilityBJetTag(); }
	vector<float> &jets_simpleSecondaryVertexHighEffBJetTag() { return cms2.jets_simpleSecondaryVertexHighEffBJetTag(); }
	vector<float> &jets_simpleSecondaryVertexHighPurBJetTags() { return cms2.jets_simpleSecondaryVertexHighPurBJetTags(); }
	vector<float> &jets_softElectronByIP3dBJetTag() { return cms2.jets_softElectronByIP3dBJetTag(); }
	vector<float> &jets_softElectronByPtBJetTag() { return cms2.jets_softElectronByPtBJetTag(); }
	vector<float> &jets_softMuonBJetTag() { return cms2.jets_softMuonBJetTag(); }
	vector<float> &jets_softMuonByIP3dBJetTag() { return cms2.jets_softMuonByIP3dBJetTag(); }
	vector<float> &jets_softMuonByPtBJetTag() { return cms2.jets_softMuonByPtBJetTag(); }
	vector<float> &jets_trackCountingHighEffBJetTag() { return cms2.jets_trackCountingHighEffBJetTag(); }
	vector<float> &jets_trackCountingHighPurBJetTag() { return cms2.jets_trackCountingHighPurBJetTag(); }
	vector<float> &pfjets_combinedSecondaryVertexBJetTag() { return cms2.pfjets_combinedSecondaryVertexBJetTag(); }
	vector<float> &pfjets_combinedSecondaryVertexMVABJetTag() { return cms2.pfjets_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &pfjets_jetBProbabilityBJetTag() { return cms2.pfjets_jetBProbabilityBJetTag(); }
	vector<float> &pfjets_jetProbabilityBJetTag() { return cms2.pfjets_jetProbabilityBJetTag(); }
	vector<float> &pfjets_simpleSecondaryVertexHighEffBJetTag() { return cms2.pfjets_simpleSecondaryVertexHighEffBJetTag(); }
	vector<float> &pfjets_simpleSecondaryVertexHighPurBJetTags() { return cms2.pfjets_simpleSecondaryVertexHighPurBJetTags(); }
	vector<float> &pfjets_softElectronByIP3dBJetTag() { return cms2.pfjets_softElectronByIP3dBJetTag(); }
	vector<float> &pfjets_softElectronByPtBJetTag() { return cms2.pfjets_softElectronByPtBJetTag(); }
	vector<float> &pfjets_softMuonBJetTag() { return cms2.pfjets_softMuonBJetTag(); }
	vector<float> &pfjets_softMuonByIP3dBJetTag() { return cms2.pfjets_softMuonByIP3dBJetTag(); }
	vector<float> &pfjets_softMuonByPtBJetTag() { return cms2.pfjets_softMuonByPtBJetTag(); }
	vector<float> &pfjets_trackCountingHighEffBJetTag() { return cms2.pfjets_trackCountingHighEffBJetTag(); }
	vector<float> &pfjets_trackCountingHighPurBJetTag() { return cms2.pfjets_trackCountingHighPurBJetTag(); }
	vector<float> &trkjets_combinedSecondaryVertexBJetTag() { return cms2.trkjets_combinedSecondaryVertexBJetTag(); }
	vector<float> &trkjets_combinedSecondaryVertexMVABJetTag() { return cms2.trkjets_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &trkjets_jetBProbabilityBJetTag() { return cms2.trkjets_jetBProbabilityBJetTag(); }
	vector<float> &trkjets_jetProbabilityBJetTag() { return cms2.trkjets_jetProbabilityBJetTag(); }
	vector<float> &trkjets_simpleSecondaryVertexHighEffBJetTag() { return cms2.trkjets_simpleSecondaryVertexHighEffBJetTag(); }
	vector<float> &trkjets_simpleSecondaryVertexHighPurBJetTags() { return cms2.trkjets_simpleSecondaryVertexHighPurBJetTags(); }
	vector<float> &trkjets_softElectronByIP3dBJetTag() { return cms2.trkjets_softElectronByIP3dBJetTag(); }
	vector<float> &trkjets_softElectronByPtBJetTag() { return cms2.trkjets_softElectronByPtBJetTag(); }
	vector<float> &trkjets_softMuonBJetTag() { return cms2.trkjets_softMuonBJetTag(); }
	vector<float> &trkjets_softMuonByIP3dBJetTag() { return cms2.trkjets_softMuonByIP3dBJetTag(); }
	vector<float> &trkjets_softMuonByPtBJetTag() { return cms2.trkjets_softMuonByPtBJetTag(); }
	vector<float> &trkjets_trackCountingHighEffBJetTag() { return cms2.trkjets_trackCountingHighEffBJetTag(); }
	vector<float> &trkjets_trackCountingHighPurBJetTag() { return cms2.trkjets_trackCountingHighPurBJetTag(); }
	vector<float> &evt_bs_covMatrix() { return cms2.evt_bs_covMatrix(); }
	vector<float> &els_mc3dr() { return cms2.els_mc3dr(); }
	vector<float> &els_mcdr() { return cms2.els_mcdr(); }
	vector<float> &jets_mc3dr() { return cms2.jets_mc3dr(); }
	vector<float> &jets_mcdr() { return cms2.jets_mcdr(); }
	vector<float> &jets_mc_emEnergy() { return cms2.jets_mc_emEnergy(); }
	vector<float> &jets_mc_gpdr() { return cms2.jets_mc_gpdr(); }
	vector<float> &jets_mc_hadEnergy() { return cms2.jets_mc_hadEnergy(); }
	vector<float> &jets_mc_invEnergy() { return cms2.jets_mc_invEnergy(); }
	vector<float> &jets_mc_otherEnergy() { return cms2.jets_mc_otherEnergy(); }
	vector<float> &mus_mc3dr() { return cms2.mus_mc3dr(); }
	vector<float> &mus_mcdr() { return cms2.mus_mcdr(); }
	vector<float> &pfjets_mc3dr() { return cms2.pfjets_mc3dr(); }
	vector<float> &pfjets_mcdr() { return cms2.pfjets_mcdr(); }
	vector<float> &pfjets_mc_emEnergy() { return cms2.pfjets_mc_emEnergy(); }
	vector<float> &pfjets_mc_gpdr() { return cms2.pfjets_mc_gpdr(); }
	vector<float> &pfjets_mc_hadEnergy() { return cms2.pfjets_mc_hadEnergy(); }
	vector<float> &pfjets_mc_invEnergy() { return cms2.pfjets_mc_invEnergy(); }
	vector<float> &pfjets_mc_otherEnergy() { return cms2.pfjets_mc_otherEnergy(); }
	vector<float> &photons_mc3dr() { return cms2.photons_mc3dr(); }
	vector<float> &photons_mcdr() { return cms2.photons_mcdr(); }
	vector<float> &trk_mc3dr() { return cms2.trk_mc3dr(); }
	vector<float> &trk_mcdr() { return cms2.trk_mcdr(); }
	vector<float> &trks_conv_dcot() { return cms2.trks_conv_dcot(); }
	vector<float> &trks_conv_dist() { return cms2.trks_conv_dist(); }
	vector<float> &els_ecalJuraIso() { return cms2.els_ecalJuraIso(); }
	vector<float> &els_ecalJuraTowerIso() { return cms2.els_ecalJuraTowerIso(); }
	vector<float> &els_hcalConeIso() { return cms2.els_hcalConeIso(); }
	vector<float> &els_tkJuraIso() { return cms2.els_tkJuraIso(); }
	vector<float> &els_jetdr() { return cms2.els_jetdr(); }
	vector<float> &els_musdr() { return cms2.els_musdr(); }
	vector<float> &els_chi2() { return cms2.els_chi2(); }
	vector<float> &els_conv_dcot() { return cms2.els_conv_dcot(); }
	vector<float> &els_conv_dist() { return cms2.els_conv_dist(); }
	vector<float> &els_conv_radius() { return cms2.els_conv_radius(); }
	vector<float> &els_d0() { return cms2.els_d0(); }
	vector<float> &els_d0Err() { return cms2.els_d0Err(); }
	vector<float> &els_d0corr() { return cms2.els_d0corr(); }
	vector<float> &els_dEtaIn() { return cms2.els_dEtaIn(); }
	vector<float> &els_dEtaOut() { return cms2.els_dEtaOut(); }
	vector<float> &els_dPhiIn() { return cms2.els_dPhiIn(); }
	vector<float> &els_dPhiInPhiOut() { return cms2.els_dPhiInPhiOut(); }
	vector<float> &els_dPhiOut() { return cms2.els_dPhiOut(); }
	vector<float> &els_deltaEtaEleClusterTrackAtCalo() { return cms2.els_deltaEtaEleClusterTrackAtCalo(); }
	vector<float> &els_deltaPhiEleClusterTrackAtCalo() { return cms2.els_deltaPhiEleClusterTrackAtCalo(); }
	vector<float> &els_e1x5() { return cms2.els_e1x5(); }
	vector<float> &els_e2x5Max() { return cms2.els_e2x5Max(); }
	vector<float> &els_e3x3() { return cms2.els_e3x3(); }
	vector<float> &els_e5x5() { return cms2.els_e5x5(); }
	vector<float> &els_eMax() { return cms2.els_eMax(); }
	vector<float> &els_eOverPIn() { return cms2.els_eOverPIn(); }
	vector<float> &els_eOverPOut() { return cms2.els_eOverPOut(); }
	vector<float> &els_eSC() { return cms2.els_eSC(); }
	vector<float> &els_eSCPresh() { return cms2.els_eSCPresh(); }
	vector<float> &els_eSCRaw() { return cms2.els_eSCRaw(); }
	vector<float> &els_eSeed() { return cms2.els_eSeed(); }
	vector<float> &els_eSeedOverPIn() { return cms2.els_eSeedOverPIn(); }
	vector<float> &els_eSeedOverPOut() { return cms2.els_eSeedOverPOut(); }
	vector<float> &els_ecalEnergy() { return cms2.els_ecalEnergy(); }
	vector<float> &els_ecalEnergyError() { return cms2.els_ecalEnergyError(); }
	vector<float> &els_ecalIso() { return cms2.els_ecalIso(); }
	vector<float> &els_ecalIso04() { return cms2.els_ecalIso04(); }
	vector<float> &els_egamma_looseId() { return cms2.els_egamma_looseId(); }
	vector<float> &els_egamma_robustHighEnergy() { return cms2.els_egamma_robustHighEnergy(); }
	vector<float> &els_egamma_robustLooseId() { return cms2.els_egamma_robustLooseId(); }
	vector<float> &els_egamma_robustTightId() { return cms2.els_egamma_robustTightId(); }
	vector<float> &els_egamma_tightId() { return cms2.els_egamma_tightId(); }
	vector<float> &els_electronMomentumError() { return cms2.els_electronMomentumError(); }
	vector<float> &els_etaErr() { return cms2.els_etaErr(); }
	vector<float> &els_etaSC() { return cms2.els_etaSC(); }
	vector<float> &els_fbrem() { return cms2.els_fbrem(); }
	vector<float> &els_hOverE() { return cms2.els_hOverE(); }
	vector<float> &els_hcalDepth1OverEcal() { return cms2.els_hcalDepth1OverEcal(); }
	vector<float> &els_hcalDepth1TowerSumEt() { return cms2.els_hcalDepth1TowerSumEt(); }
	vector<float> &els_hcalDepth1TowerSumEt04() { return cms2.els_hcalDepth1TowerSumEt04(); }
	vector<float> &els_hcalDepth2OverEcal() { return cms2.els_hcalDepth2OverEcal(); }
	vector<float> &els_hcalDepth2TowerSumEt() { return cms2.els_hcalDepth2TowerSumEt(); }
	vector<float> &els_hcalDepth2TowerSumEt04() { return cms2.els_hcalDepth2TowerSumEt04(); }
	vector<float> &els_hcalIso() { return cms2.els_hcalIso(); }
	vector<float> &els_hcalIso04() { return cms2.els_hcalIso04(); }
	vector<float> &els_layer1_charge() { return cms2.els_layer1_charge(); }
	vector<float> &els_mva() { return cms2.els_mva(); }
	vector<float> &els_ndof() { return cms2.els_ndof(); }
	vector<float> &els_phiErr() { return cms2.els_phiErr(); }
	vector<float> &els_phiSC() { return cms2.els_phiSC(); }
	vector<float> &els_ptErr() { return cms2.els_ptErr(); }
	vector<float> &els_sigmaEtaEta() { return cms2.els_sigmaEtaEta(); }
	vector<float> &els_sigmaIEtaIEta() { return cms2.els_sigmaIEtaIEta(); }
	vector<float> &els_sigmaIEtaIEtaSC() { return cms2.els_sigmaIEtaIEtaSC(); }
	vector<float> &els_sigmaIPhiIPhi() { return cms2.els_sigmaIPhiIPhi(); }
	vector<float> &els_sigmaIPhiIPhiSC() { return cms2.els_sigmaIPhiIPhiSC(); }
	vector<float> &els_sigmaPhiPhi() { return cms2.els_sigmaPhiPhi(); }
	vector<float> &els_tkIso() { return cms2.els_tkIso(); }
	vector<float> &els_tkIso04() { return cms2.els_tkIso04(); }
	vector<float> &els_trackMomentumError() { return cms2.els_trackMomentumError(); }
	vector<float> &els_trkdr() { return cms2.els_trkdr(); }
	vector<float> &els_trkshFrac() { return cms2.els_trkshFrac(); }
	vector<float> &els_z0() { return cms2.els_z0(); }
	vector<float> &els_z0Err() { return cms2.els_z0Err(); }
	vector<float> &els_z0corr() { return cms2.els_z0corr(); }
	vector<float> &gsftrks_chi2() { return cms2.gsftrks_chi2(); }
	vector<float> &gsftrks_d0() { return cms2.gsftrks_d0(); }
	vector<float> &gsftrks_d0Err() { return cms2.gsftrks_d0Err(); }
	vector<float> &gsftrks_d0corr() { return cms2.gsftrks_d0corr(); }
	vector<float> &gsftrks_d0corrPhi() { return cms2.gsftrks_d0corrPhi(); }
	vector<float> &gsftrks_d0phiCov() { return cms2.gsftrks_d0phiCov(); }
	vector<float> &gsftrks_etaErr() { return cms2.gsftrks_etaErr(); }
	vector<float> &gsftrks_layer1_charge() { return cms2.gsftrks_layer1_charge(); }
	vector<float> &gsftrks_ndof() { return cms2.gsftrks_ndof(); }
	vector<float> &gsftrks_phiErr() { return cms2.gsftrks_phiErr(); }
	vector<float> &gsftrks_ptErr() { return cms2.gsftrks_ptErr(); }
	vector<float> &gsftrks_z0() { return cms2.gsftrks_z0(); }
	vector<float> &gsftrks_z0Err() { return cms2.gsftrks_z0Err(); }
	vector<float> &gsftrks_z0corr() { return cms2.gsftrks_z0corr(); }
	vector<float> &hyp_Ht() { return cms2.hyp_Ht(); }
	vector<float> &hyp_dPhi_nJet_metMuonJESCorr() { return cms2.hyp_dPhi_nJet_metMuonJESCorr(); }
	vector<float> &hyp_dPhi_nJet_muCorrMet() { return cms2.hyp_dPhi_nJet_muCorrMet(); }
	vector<float> &hyp_dPhi_nJet_tcMet() { return cms2.hyp_dPhi_nJet_tcMet(); }
	vector<float> &hyp_dPhi_nJet_unCorrMet() { return cms2.hyp_dPhi_nJet_unCorrMet(); }
	vector<float> &hyp_ll_chi2() { return cms2.hyp_ll_chi2(); }
	vector<float> &hyp_ll_d0() { return cms2.hyp_ll_d0(); }
	vector<float> &hyp_ll_d0Err() { return cms2.hyp_ll_d0Err(); }
	vector<float> &hyp_ll_d0corr() { return cms2.hyp_ll_d0corr(); }
	vector<float> &hyp_ll_dPhi_metMuonJESCorr() { return cms2.hyp_ll_dPhi_metMuonJESCorr(); }
	vector<float> &hyp_ll_dPhi_muCorrMet() { return cms2.hyp_ll_dPhi_muCorrMet(); }
	vector<float> &hyp_ll_dPhi_tcMet() { return cms2.hyp_ll_dPhi_tcMet(); }
	vector<float> &hyp_ll_dPhi_unCorrMet() { return cms2.hyp_ll_dPhi_unCorrMet(); }
	vector<float> &hyp_ll_etaErr() { return cms2.hyp_ll_etaErr(); }
	vector<float> &hyp_ll_ndof() { return cms2.hyp_ll_ndof(); }
	vector<float> &hyp_ll_phiErr() { return cms2.hyp_ll_phiErr(); }
	vector<float> &hyp_ll_ptErr() { return cms2.hyp_ll_ptErr(); }
	vector<float> &hyp_ll_z0() { return cms2.hyp_ll_z0(); }
	vector<float> &hyp_ll_z0Err() { return cms2.hyp_ll_z0Err(); }
	vector<float> &hyp_ll_z0corr() { return cms2.hyp_ll_z0corr(); }
	vector<float> &hyp_lt_chi2() { return cms2.hyp_lt_chi2(); }
	vector<float> &hyp_lt_d0() { return cms2.hyp_lt_d0(); }
	vector<float> &hyp_lt_d0Err() { return cms2.hyp_lt_d0Err(); }
	vector<float> &hyp_lt_d0corr() { return cms2.hyp_lt_d0corr(); }
	vector<float> &hyp_lt_dPhi_metMuonJESCorr() { return cms2.hyp_lt_dPhi_metMuonJESCorr(); }
	vector<float> &hyp_lt_dPhi_muCorrMet() { return cms2.hyp_lt_dPhi_muCorrMet(); }
	vector<float> &hyp_lt_dPhi_tcMet() { return cms2.hyp_lt_dPhi_tcMet(); }
	vector<float> &hyp_lt_dPhi_unCorrMet() { return cms2.hyp_lt_dPhi_unCorrMet(); }
	vector<float> &hyp_lt_etaErr() { return cms2.hyp_lt_etaErr(); }
	vector<float> &hyp_lt_ndof() { return cms2.hyp_lt_ndof(); }
	vector<float> &hyp_lt_phiErr() { return cms2.hyp_lt_phiErr(); }
	vector<float> &hyp_lt_ptErr() { return cms2.hyp_lt_ptErr(); }
	vector<float> &hyp_lt_z0() { return cms2.hyp_lt_z0(); }
	vector<float> &hyp_lt_z0Err() { return cms2.hyp_lt_z0Err(); }
	vector<float> &hyp_lt_z0corr() { return cms2.hyp_lt_z0corr(); }
	vector<float> &hyp_mt2_metMuonJESCorr() { return cms2.hyp_mt2_metMuonJESCorr(); }
	vector<float> &hyp_mt2_muCorrMet() { return cms2.hyp_mt2_muCorrMet(); }
	vector<float> &hyp_mt2_tcMet() { return cms2.hyp_mt2_tcMet(); }
	vector<float> &hyp_sumJetPt() { return cms2.hyp_sumJetPt(); }
	vector<float> &hyp_FVFit_chi2ndf() { return cms2.hyp_FVFit_chi2ndf(); }
	vector<float> &hyp_FVFit_prob() { return cms2.hyp_FVFit_prob(); }
	vector<float> &hyp_FVFit_v4cxx() { return cms2.hyp_FVFit_v4cxx(); }
	vector<float> &hyp_FVFit_v4cxy() { return cms2.hyp_FVFit_v4cxy(); }
	vector<float> &hyp_FVFit_v4cyy() { return cms2.hyp_FVFit_v4cyy(); }
	vector<float> &hyp_FVFit_v4czx() { return cms2.hyp_FVFit_v4czx(); }
	vector<float> &hyp_FVFit_v4czy() { return cms2.hyp_FVFit_v4czy(); }
	vector<float> &hyp_FVFit_v4czz() { return cms2.hyp_FVFit_v4czz(); }
	vector<float> &hyp_ll_ecaliso() { return cms2.hyp_ll_ecaliso(); }
	vector<float> &hyp_ll_trkiso() { return cms2.hyp_ll_trkiso(); }
	vector<float> &hyp_lt_ecaliso() { return cms2.hyp_lt_ecaliso(); }
	vector<float> &hyp_lt_trkiso() { return cms2.hyp_lt_trkiso(); }
	vector<float> &jets_approximatefHPD() { return cms2.jets_approximatefHPD(); }
	vector<float> &jets_approximatefRBX() { return cms2.jets_approximatefRBX(); }
	vector<float> &jets_cor() { return cms2.jets_cor(); }
	vector<float> &jets_emFrac() { return cms2.jets_emFrac(); }
	vector<float> &jets_fHPD() { return cms2.jets_fHPD(); }
	vector<float> &jets_fRBX() { return cms2.jets_fRBX(); }
	vector<float> &jets_fSubDetector1() { return cms2.jets_fSubDetector1(); }
	vector<float> &jets_fSubDetector2() { return cms2.jets_fSubDetector2(); }
	vector<float> &jets_fSubDetector3() { return cms2.jets_fSubDetector3(); }
	vector<float> &jets_fSubDetector4() { return cms2.jets_fSubDetector4(); }
	vector<float> &jets_hitsInN90() { return cms2.jets_hitsInN90(); }
	vector<float> &jets_n90Hits() { return cms2.jets_n90Hits(); }
	vector<float> &jets_nECALTowers() { return cms2.jets_nECALTowers(); }
	vector<float> &jets_nHCALTowers() { return cms2.jets_nHCALTowers(); }
	vector<float> &jets_restrictedEMF() { return cms2.jets_restrictedEMF(); }
	vector<float> &jpts_cor() { return cms2.jpts_cor(); }
	vector<float> &jpts_emFrac() { return cms2.jpts_emFrac(); }
	vector<float> &mus_met_deltax() { return cms2.mus_met_deltax(); }
	vector<float> &mus_met_deltay() { return cms2.mus_met_deltay(); }
	vector<float> &mus_eledr() { return cms2.mus_eledr(); }
	vector<float> &mus_jetdr() { return cms2.mus_jetdr(); }
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
	vector<float> &mus_gfit_d0() { return cms2.mus_gfit_d0(); }
	vector<float> &mus_gfit_d0Err() { return cms2.mus_gfit_d0Err(); }
	vector<float> &mus_gfit_d0corr() { return cms2.mus_gfit_d0corr(); }
	vector<float> &mus_gfit_ndof() { return cms2.mus_gfit_ndof(); }
	vector<float> &mus_gfit_qoverp() { return cms2.mus_gfit_qoverp(); }
	vector<float> &mus_gfit_qoverpError() { return cms2.mus_gfit_qoverpError(); }
	vector<float> &mus_gfit_z0() { return cms2.mus_gfit_z0(); }
	vector<float> &mus_gfit_z0Err() { return cms2.mus_gfit_z0Err(); }
	vector<float> &mus_gfit_z0corr() { return cms2.mus_gfit_z0corr(); }
	vector<float> &mus_iso03_emEt() { return cms2.mus_iso03_emEt(); }
	vector<float> &mus_iso03_hadEt() { return cms2.mus_iso03_hadEt(); }
	vector<float> &mus_iso03_hoEt() { return cms2.mus_iso03_hoEt(); }
	vector<float> &mus_iso03_sumPt() { return cms2.mus_iso03_sumPt(); }
	vector<float> &mus_iso05_emEt() { return cms2.mus_iso05_emEt(); }
	vector<float> &mus_iso05_hadEt() { return cms2.mus_iso05_hadEt(); }
	vector<float> &mus_iso05_hoEt() { return cms2.mus_iso05_hoEt(); }
	vector<float> &mus_iso05_sumPt() { return cms2.mus_iso05_sumPt(); }
	vector<float> &mus_iso_ecalvetoDep() { return cms2.mus_iso_ecalvetoDep(); }
	vector<float> &mus_iso_hcalvetoDep() { return cms2.mus_iso_hcalvetoDep(); }
	vector<float> &mus_iso_hovetoDep() { return cms2.mus_iso_hovetoDep(); }
	vector<float> &mus_iso_trckvetoDep() { return cms2.mus_iso_trckvetoDep(); }
	vector<float> &mus_ndof() { return cms2.mus_ndof(); }
	vector<float> &mus_phiErr() { return cms2.mus_phiErr(); }
	vector<float> &mus_ptErr() { return cms2.mus_ptErr(); }
	vector<float> &mus_qoverp() { return cms2.mus_qoverp(); }
	vector<float> &mus_qoverpError() { return cms2.mus_qoverpError(); }
	vector<float> &mus_sta_chi2() { return cms2.mus_sta_chi2(); }
	vector<float> &mus_sta_d0() { return cms2.mus_sta_d0(); }
	vector<float> &mus_sta_d0Err() { return cms2.mus_sta_d0Err(); }
	vector<float> &mus_sta_d0corr() { return cms2.mus_sta_d0corr(); }
	vector<float> &mus_sta_ndof() { return cms2.mus_sta_ndof(); }
	vector<float> &mus_sta_qoverp() { return cms2.mus_sta_qoverp(); }
	vector<float> &mus_sta_qoverpError() { return cms2.mus_sta_qoverpError(); }
	vector<float> &mus_sta_z0() { return cms2.mus_sta_z0(); }
	vector<float> &mus_sta_z0Err() { return cms2.mus_sta_z0Err(); }
	vector<float> &mus_sta_z0corr() { return cms2.mus_sta_z0corr(); }
	vector<float> &mus_timeAtIpInOut() { return cms2.mus_timeAtIpInOut(); }
	vector<float> &mus_timeAtIpInOutErr() { return cms2.mus_timeAtIpInOutErr(); }
	vector<float> &mus_timeAtIpOutIn() { return cms2.mus_timeAtIpOutIn(); }
	vector<float> &mus_timeAtIpOutInErr() { return cms2.mus_timeAtIpOutInErr(); }
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
	vector<float> &jets_pat_combinedSecondaryVertexBJetTag() { return cms2.jets_pat_combinedSecondaryVertexBJetTag(); }
	vector<float> &jets_pat_combinedSecondaryVertexMVABJetTag() { return cms2.jets_pat_combinedSecondaryVertexMVABJetTag(); }
	vector<float> &jets_pat_coneIsolationTauJetTag() { return cms2.jets_pat_coneIsolationTauJetTag(); }
	vector<float> &jets_pat_impactParameterMVABJetTag() { return cms2.jets_pat_impactParameterMVABJetTag(); }
	vector<float> &jets_pat_jetBProbabilityBJetTag() { return cms2.jets_pat_jetBProbabilityBJetTag(); }
	vector<float> &jets_pat_jetCharge() { return cms2.jets_pat_jetCharge(); }
	vector<float> &jets_pat_jetProbabilityBJetTag() { return cms2.jets_pat_jetProbabilityBJetTag(); }
	vector<float> &jets_pat_noCorrF() { return cms2.jets_pat_noCorrF(); }
	vector<float> &jets_pat_simpleSecondaryVertexHighEffBJetTag() { return cms2.jets_pat_simpleSecondaryVertexHighEffBJetTag(); }
	vector<float> &jets_pat_simpleSecondaryVertexHighPurBJetTag() { return cms2.jets_pat_simpleSecondaryVertexHighPurBJetTag(); }
	vector<float> &jets_pat_softElectronByIP3dBJetTag() { return cms2.jets_pat_softElectronByIP3dBJetTag(); }
	vector<float> &jets_pat_softElectronByPtBJetTag() { return cms2.jets_pat_softElectronByPtBJetTag(); }
	vector<float> &jets_pat_softMuonBJetTag() { return cms2.jets_pat_softMuonBJetTag(); }
	vector<float> &jets_pat_softMuonByIP3dBJetTag() { return cms2.jets_pat_softMuonByIP3dBJetTag(); }
	vector<float> &jets_pat_softMuonByPtBJetTag() { return cms2.jets_pat_softMuonByPtBJetTag(); }
	vector<float> &jets_pat_trackCountingHighEffBJetTag() { return cms2.jets_pat_trackCountingHighEffBJetTag(); }
	vector<float> &jets_pat_trackCountingHighPurBJetTag() { return cms2.jets_pat_trackCountingHighPurBJetTag(); }
	vector<float> &mus_pat_caloIso() { return cms2.mus_pat_caloIso(); }
	vector<float> &mus_pat_calovetoDep() { return cms2.mus_pat_calovetoDep(); }
	vector<float> &mus_pat_ecalIso() { return cms2.mus_pat_ecalIso(); }
	vector<float> &mus_pat_ecalvetoDep() { return cms2.mus_pat_ecalvetoDep(); }
	vector<float> &mus_pat_hcalIso() { return cms2.mus_pat_hcalIso(); }
	vector<float> &mus_pat_hcalvetoDep() { return cms2.mus_pat_hcalvetoDep(); }
	vector<float> &mus_pat_trackIso() { return cms2.mus_pat_trackIso(); }
	vector<float> &mus_pat_trckvetoDep() { return cms2.mus_pat_trckvetoDep(); }
	vector<float> &pfels_deltaP() { return cms2.pfels_deltaP(); }
	vector<float> &pfels_ecalE() { return cms2.pfels_ecalE(); }
	vector<float> &pfels_hcalE() { return cms2.pfels_hcalE(); }
	vector<float> &pfels_isoChargedHadrons() { return cms2.pfels_isoChargedHadrons(); }
	vector<float> &pfels_isoNeutralHadrons() { return cms2.pfels_isoNeutralHadrons(); }
	vector<float> &pfels_isoPhotons() { return cms2.pfels_isoPhotons(); }
	vector<float> &pfels_mva_emu() { return cms2.pfels_mva_emu(); }
	vector<float> &pfels_mva_epi() { return cms2.pfels_mva_epi(); }
	vector<float> &pfels_mva_nothing_gamma() { return cms2.pfels_mva_nothing_gamma(); }
	vector<float> &pfels_mva_nothing_nh() { return cms2.pfels_mva_nothing_nh(); }
	vector<float> &pfels_mva_pimu() { return cms2.pfels_mva_pimu(); }
	vector<float> &pfels_pS1E() { return cms2.pfels_pS1E(); }
	vector<float> &pfels_pS2E() { return cms2.pfels_pS2E(); }
	vector<float> &pfels_rawEcalE() { return cms2.pfels_rawEcalE(); }
	vector<float> &pfels_rawHcalE() { return cms2.pfels_rawHcalE(); }
	vector<float> &pfjets_chargedEmE() { return cms2.pfjets_chargedEmE(); }
	vector<float> &pfjets_chargedHadronE() { return cms2.pfjets_chargedHadronE(); }
	vector<float> &pfjets_cor() { return cms2.pfjets_cor(); }
	vector<float> &pfjets_neutralEmE() { return cms2.pfjets_neutralEmE(); }
	vector<float> &pfjets_neutralHadronE() { return cms2.pfjets_neutralHadronE(); }
	vector<float> &pfmus_deltaP() { return cms2.pfmus_deltaP(); }
	vector<float> &pfmus_ecalE() { return cms2.pfmus_ecalE(); }
	vector<float> &pfmus_hcalE() { return cms2.pfmus_hcalE(); }
	vector<float> &pfmus_isoChargedHadrons() { return cms2.pfmus_isoChargedHadrons(); }
	vector<float> &pfmus_isoNeutralHadrons() { return cms2.pfmus_isoNeutralHadrons(); }
	vector<float> &pfmus_isoPhotons() { return cms2.pfmus_isoPhotons(); }
	vector<float> &pfmus_mva_emu() { return cms2.pfmus_mva_emu(); }
	vector<float> &pfmus_mva_epi() { return cms2.pfmus_mva_epi(); }
	vector<float> &pfmus_mva_nothing_gamma() { return cms2.pfmus_mva_nothing_gamma(); }
	vector<float> &pfmus_mva_nothing_nh() { return cms2.pfmus_mva_nothing_nh(); }
	vector<float> &pfmus_mva_pimu() { return cms2.pfmus_mva_pimu(); }
	vector<float> &pfmus_pS1E() { return cms2.pfmus_pS1E(); }
	vector<float> &pfmus_pS2E() { return cms2.pfmus_pS2E(); }
	vector<float> &pfmus_rawEcalE() { return cms2.pfmus_rawEcalE(); }
	vector<float> &pfmus_rawHcalE() { return cms2.pfmus_rawHcalE(); }
	vector<float> &photons_e1x5() { return cms2.photons_e1x5(); }
	vector<float> &photons_e2x5Max() { return cms2.photons_e2x5Max(); }
	vector<float> &photons_e3x3() { return cms2.photons_e3x3(); }
	vector<float> &photons_e5x5() { return cms2.photons_e5x5(); }
	vector<float> &photons_ecalIso03() { return cms2.photons_ecalIso03(); }
	vector<float> &photons_ecalIso04() { return cms2.photons_ecalIso04(); }
	vector<float> &photons_hOverE() { return cms2.photons_hOverE(); }
	vector<float> &photons_hcalIso03() { return cms2.photons_hcalIso03(); }
	vector<float> &photons_hcalIso04() { return cms2.photons_hcalIso04(); }
	vector<float> &photons_ntkIsoHollow03() { return cms2.photons_ntkIsoHollow03(); }
	vector<float> &photons_ntkIsoHollow04() { return cms2.photons_ntkIsoHollow04(); }
	vector<float> &photons_ntkIsoSolid03() { return cms2.photons_ntkIsoSolid03(); }
	vector<float> &photons_ntkIsoSolid04() { return cms2.photons_ntkIsoSolid04(); }
	vector<float> &photons_sigmaEtaEta() { return cms2.photons_sigmaEtaEta(); }
	vector<float> &photons_sigmaIEtaIEta() { return cms2.photons_sigmaIEtaIEta(); }
	vector<float> &photons_swissSeed() { return cms2.photons_swissSeed(); }
	vector<float> &photons_tkIsoHollow03() { return cms2.photons_tkIsoHollow03(); }
	vector<float> &photons_tkIsoHollow04() { return cms2.photons_tkIsoHollow04(); }
	vector<float> &photons_tkIsoSolid03() { return cms2.photons_tkIsoSolid03(); }
	vector<float> &photons_tkIsoSolid04() { return cms2.photons_tkIsoSolid04(); }
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
	vector<float> &scs_sigmaIEtaIEtaSC() { return cms2.scs_sigmaIEtaIEtaSC(); }
	vector<float> &scs_sigmaIEtaIPhi() { return cms2.scs_sigmaIEtaIPhi(); }
	vector<float> &scs_sigmaIEtaIPhiSC() { return cms2.scs_sigmaIEtaIPhiSC(); }
	vector<float> &scs_sigmaIPhiIPhi() { return cms2.scs_sigmaIPhiIPhi(); }
	vector<float> &scs_sigmaIPhiIPhiSC() { return cms2.scs_sigmaIPhiIPhiSC(); }
	vector<float> &scs_sigmaPhiPhi() { return cms2.scs_sigmaPhiPhi(); }
	vector<float> &scs_timeSeed() { return cms2.scs_timeSeed(); }
	vector<float> &mus_tcmet_deltax() { return cms2.mus_tcmet_deltax(); }
	vector<float> &mus_tcmet_deltay() { return cms2.mus_tcmet_deltay(); }
	vector<float> &evt35X_mus_tcmet_deltax() { return cms2.evt35X_mus_tcmet_deltax(); }
	vector<float> &evt35X_mus_tcmet_deltay() { return cms2.evt35X_mus_tcmet_deltay(); }
	vector<float> &trks_chi2() { return cms2.trks_chi2(); }
	vector<float> &trks_d0() { return cms2.trks_d0(); }
	vector<float> &trks_d0Err() { return cms2.trks_d0Err(); }
	vector<float> &trks_d0corr() { return cms2.trks_d0corr(); }
	vector<float> &trks_d0corrPhi() { return cms2.trks_d0corrPhi(); }
	vector<float> &trks_d0phiCov() { return cms2.trks_d0phiCov(); }
	vector<float> &trks_etaErr() { return cms2.trks_etaErr(); }
	vector<float> &trks_layer1_charge() { return cms2.trks_layer1_charge(); }
	vector<float> &trks_ndof() { return cms2.trks_ndof(); }
	vector<float> &trks_phiErr() { return cms2.trks_phiErr(); }
	vector<float> &trks_ptErr() { return cms2.trks_ptErr(); }
	vector<float> &trks_z0() { return cms2.trks_z0(); }
	vector<float> &trks_z0Err() { return cms2.trks_z0Err(); }
	vector<float> &trks_z0corr() { return cms2.trks_z0corr(); }
	vector<float> &trkjets_cor() { return cms2.trkjets_cor(); }
	vector<float> &trks_d0Errvtx() { return cms2.trks_d0Errvtx(); }
	vector<float> &trks_d0vtx() { return cms2.trks_d0vtx(); }
	vector<float> &vtxs_chi2() { return cms2.vtxs_chi2(); }
	vector<float> &vtxs_ndof() { return cms2.vtxs_ndof(); }
	vector<float> &vtxs_sumpt() { return cms2.vtxs_sumpt(); }
	vector<float> &vtxs_xError() { return cms2.vtxs_xError(); }
	vector<float> &vtxs_yError() { return cms2.vtxs_yError(); }
	vector<float> &vtxs_zError() { return cms2.vtxs_zError(); }
	vector<vector<float> > &vtxs_covMatrix() { return cms2.vtxs_covMatrix(); }
	int &evt_cscLooseHaloId() { return cms2.evt_cscLooseHaloId(); }
	int &evt_cscTightHaloId() { return cms2.evt_cscTightHaloId(); }
	int &evt_ecalLooseHaloId() { return cms2.evt_ecalLooseHaloId(); }
	int &evt_ecalTightHaloId() { return cms2.evt_ecalTightHaloId(); }
	int &evt_extremeTightHaloId() { return cms2.evt_extremeTightHaloId(); }
	int &evt_globalLooseHaloId() { return cms2.evt_globalLooseHaloId(); }
	int &evt_globalTightHaloId() { return cms2.evt_globalTightHaloId(); }
	int &evt_hcalLooseHaloId() { return cms2.evt_hcalLooseHaloId(); }
	int &evt_hcalTightHaloId() { return cms2.evt_hcalTightHaloId(); }
	int &evt_looseHaloId() { return cms2.evt_looseHaloId(); }
	int &evt_nHaloLikeTracks() { return cms2.evt_nHaloLikeTracks(); }
	int &evt_nHaloTriggerCandidates() { return cms2.evt_nHaloTriggerCandidates(); }
	int &evt_tightHaloId() { return cms2.evt_tightHaloId(); }
	int &evt_bsType() { return cms2.evt_bsType(); }
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
	int &pdfinfo_id1() { return cms2.pdfinfo_id1(); }
	int &pdfinfo_id2() { return cms2.pdfinfo_id2(); }
	vector<int> &evt_ecaliPhiSuspects() { return cms2.evt_ecaliPhiSuspects(); }
	vector<int> &evt_globaliPhiSuspects() { return cms2.evt_globaliPhiSuspects(); }
	vector<int> &evt_hcaliPhiSuspects() { return cms2.evt_hcaliPhiSuspects(); }
	vector<int> &els_mc3_id() { return cms2.els_mc3_id(); }
	vector<int> &els_mc3idx() { return cms2.els_mc3idx(); }
	vector<int> &els_mc3_motherid() { return cms2.els_mc3_motherid(); }
	vector<int> &els_mc3_motheridx() { return cms2.els_mc3_motheridx(); }
	vector<int> &els_mc_id() { return cms2.els_mc_id(); }
	vector<int> &els_mcidx() { return cms2.els_mcidx(); }
	vector<int> &els_mc_motherid() { return cms2.els_mc_motherid(); }
	vector<int> &jets_mc3_id() { return cms2.jets_mc3_id(); }
	vector<int> &jets_mc3idx() { return cms2.jets_mc3idx(); }
	vector<int> &jets_mc_gpidx() { return cms2.jets_mc_gpidx(); }
	vector<int> &jets_mc_id() { return cms2.jets_mc_id(); }
	vector<int> &jets_mcidx() { return cms2.jets_mcidx(); }
	vector<int> &jets_mc_motherid() { return cms2.jets_mc_motherid(); }
	vector<int> &mus_mc3_id() { return cms2.mus_mc3_id(); }
	vector<int> &mus_mc3idx() { return cms2.mus_mc3idx(); }
	vector<int> &mus_mc3_motherid() { return cms2.mus_mc3_motherid(); }
	vector<int> &mus_mc3_motheridx() { return cms2.mus_mc3_motheridx(); }
	vector<int> &mus_mc_id() { return cms2.mus_mc_id(); }
	vector<int> &mus_mcidx() { return cms2.mus_mcidx(); }
	vector<int> &mus_mc_motherid() { return cms2.mus_mc_motherid(); }
	vector<int> &pfjets_mc3_id() { return cms2.pfjets_mc3_id(); }
	vector<int> &pfjets_mc3idx() { return cms2.pfjets_mc3idx(); }
	vector<int> &pfjets_mc_gpidx() { return cms2.pfjets_mc_gpidx(); }
	vector<int> &pfjets_mc_id() { return cms2.pfjets_mc_id(); }
	vector<int> &pfjets_mcidx() { return cms2.pfjets_mcidx(); }
	vector<int> &pfjets_mc_motherid() { return cms2.pfjets_mc_motherid(); }
	vector<int> &photons_mc3_id() { return cms2.photons_mc3_id(); }
	vector<int> &photons_mc3idx() { return cms2.photons_mc3idx(); }
	vector<int> &photons_mc3_motherid() { return cms2.photons_mc3_motherid(); }
	vector<int> &photons_mc3_motheridx() { return cms2.photons_mc3_motheridx(); }
	vector<int> &photons_mc_id() { return cms2.photons_mc_id(); }
	vector<int> &photons_mcidx() { return cms2.photons_mcidx(); }
	vector<int> &photons_mc_motherid() { return cms2.photons_mc_motherid(); }
	vector<int> &trk_mc3_id() { return cms2.trk_mc3_id(); }
	vector<int> &trk_mc3idx() { return cms2.trk_mc3idx(); }
	vector<int> &trk_mc3_motherid() { return cms2.trk_mc3_motherid(); }
	vector<int> &trk_mc3_motheridx() { return cms2.trk_mc3_motheridx(); }
	vector<int> &trk_mc_id() { return cms2.trk_mc_id(); }
	vector<int> &trk_mcidx() { return cms2.trk_mcidx(); }
	vector<int> &trk_mc_motherid() { return cms2.trk_mc_motherid(); }
	vector<int> &trks_conv_tkidx() { return cms2.trks_conv_tkidx(); }
	vector<int> &els_closestJet() { return cms2.els_closestJet(); }
	vector<int> &els_closestMuon() { return cms2.els_closestMuon(); }
	vector<int> &els_pfelsidx() { return cms2.els_pfelsidx(); }
	vector<int> &els_category() { return cms2.els_category(); }
	vector<int> &els_charge() { return cms2.els_charge(); }
	vector<int> &els_class() { return cms2.els_class(); }
	vector<int> &els_conv_tkidx() { return cms2.els_conv_tkidx(); }
	vector<int> &els_exp_innerlayers() { return cms2.els_exp_innerlayers(); }
	vector<int> &els_exp_outerlayers() { return cms2.els_exp_outerlayers(); }
	vector<int> &els_fiduciality() { return cms2.els_fiduciality(); }
	vector<int> &els_layer1_det() { return cms2.els_layer1_det(); }
	vector<int> &els_layer1_layer() { return cms2.els_layer1_layer(); }
	vector<int> &els_layer1_sizerphi() { return cms2.els_layer1_sizerphi(); }
	vector<int> &els_layer1_sizerz() { return cms2.els_layer1_sizerz(); }
	vector<int> &els_lostHits() { return cms2.els_lostHits(); }
	vector<int> &els_lost_pixelhits() { return cms2.els_lost_pixelhits(); }
	vector<int> &els_nSeed() { return cms2.els_nSeed(); }
	vector<int> &els_sccharge() { return cms2.els_sccharge(); }
	vector<int> &els_scindex() { return cms2.els_scindex(); }
	vector<int> &els_trk_charge() { return cms2.els_trk_charge(); }
	vector<int> &els_trkidx() { return cms2.els_trkidx(); }
	vector<int> &els_type() { return cms2.els_type(); }
	vector<int> &els_validHits() { return cms2.els_validHits(); }
	vector<int> &els_valid_pixelhits() { return cms2.els_valid_pixelhits(); }
	vector<int> &genps_id() { return cms2.genps_id(); }
	vector<int> &genps_id_mother() { return cms2.genps_id_mother(); }
	vector<int> &genps_status() { return cms2.genps_status(); }
	vector<int> &gsftrks_charge() { return cms2.gsftrks_charge(); }
	vector<int> &gsftrks_exp_innerlayers() { return cms2.gsftrks_exp_innerlayers(); }
	vector<int> &gsftrks_exp_outerlayers() { return cms2.gsftrks_exp_outerlayers(); }
	vector<int> &gsftrks_layer1_det() { return cms2.gsftrks_layer1_det(); }
	vector<int> &gsftrks_layer1_layer() { return cms2.gsftrks_layer1_layer(); }
	vector<int> &gsftrks_layer1_sizerphi() { return cms2.gsftrks_layer1_sizerphi(); }
	vector<int> &gsftrks_layer1_sizerz() { return cms2.gsftrks_layer1_sizerz(); }
	vector<int> &gsftrks_lostHits() { return cms2.gsftrks_lostHits(); }
	vector<int> &gsftrks_lost_pixelhits() { return cms2.gsftrks_lost_pixelhits(); }
	vector<int> &gsftrks_nlayers() { return cms2.gsftrks_nlayers(); }
	vector<int> &gsftrks_nlayers3D() { return cms2.gsftrks_nlayers3D(); }
	vector<int> &gsftrks_nlayersLost() { return cms2.gsftrks_nlayersLost(); }
	vector<int> &gsftrks_validHits() { return cms2.gsftrks_validHits(); }
	vector<int> &gsftrks_valid_pixelhits() { return cms2.gsftrks_valid_pixelhits(); }
	vector<int> &hyp_ll_charge() { return cms2.hyp_ll_charge(); }
	vector<int> &hyp_ll_id() { return cms2.hyp_ll_id(); }
	vector<int> &hyp_ll_index() { return cms2.hyp_ll_index(); }
	vector<int> &hyp_ll_lostHits() { return cms2.hyp_ll_lostHits(); }
	vector<int> &hyp_ll_validHits() { return cms2.hyp_ll_validHits(); }
	vector<int> &hyp_lt_charge() { return cms2.hyp_lt_charge(); }
	vector<int> &hyp_lt_id() { return cms2.hyp_lt_id(); }
	vector<int> &hyp_lt_index() { return cms2.hyp_lt_index(); }
	vector<int> &hyp_lt_lostHits() { return cms2.hyp_lt_lostHits(); }
	vector<int> &hyp_lt_validHits() { return cms2.hyp_lt_validHits(); }
	vector<int> &hyp_njets() { return cms2.hyp_njets(); }
	vector<int> &hyp_nojets() { return cms2.hyp_nojets(); }
	vector<int> &hyp_type() { return cms2.hyp_type(); }
	vector<int> &hyp_FVFit_ndf() { return cms2.hyp_FVFit_ndf(); }
	vector<int> &hyp_FVFit_status() { return cms2.hyp_FVFit_status(); }
	vector<int> &hyp_ll_mc_id() { return cms2.hyp_ll_mc_id(); }
	vector<int> &hyp_ll_mc_motherid() { return cms2.hyp_ll_mc_motherid(); }
	vector<int> &hyp_lt_mc_id() { return cms2.hyp_lt_mc_id(); }
	vector<int> &hyp_lt_mc_motherid() { return cms2.hyp_lt_mc_motherid(); }
	vector<int> &hyp_quadlep_first_type() { return cms2.hyp_quadlep_first_type(); }
	vector<int> &hyp_quadlep_fourth_type() { return cms2.hyp_quadlep_fourth_type(); }
	vector<int> &hyp_quadlep_second_type() { return cms2.hyp_quadlep_second_type(); }
	vector<int> &hyp_quadlep_third_type() { return cms2.hyp_quadlep_third_type(); }
	vector<int> &hyp_trilep_first_type() { return cms2.hyp_trilep_first_type(); }
	vector<int> &hyp_trilep_second_type() { return cms2.hyp_trilep_second_type(); }
	vector<int> &hyp_trilep_third_type() { return cms2.hyp_trilep_third_type(); }
	vector<int> &jets_closestElectron() { return cms2.jets_closestElectron(); }
	vector<int> &jets_closestMuon() { return cms2.jets_closestMuon(); }
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
	vector<int> &mus_met_flag() { return cms2.mus_met_flag(); }
	vector<int> &mus_closestEle() { return cms2.mus_closestEle(); }
	vector<int> &mus_closestJet() { return cms2.mus_closestJet(); }
	vector<int> &mus_pfmusidx() { return cms2.mus_pfmusidx(); }
	vector<int> &mus_charge() { return cms2.mus_charge(); }
	vector<int> &mus_gfit_validHits() { return cms2.mus_gfit_validHits(); }
	vector<int> &mus_gfit_validSTAHits() { return cms2.mus_gfit_validSTAHits(); }
	vector<int> &mus_gfit_validSiHits() { return cms2.mus_gfit_validSiHits(); }
	vector<int> &mus_goodmask() { return cms2.mus_goodmask(); }
	vector<int> &mus_iso03_ntrk() { return cms2.mus_iso03_ntrk(); }
	vector<int> &mus_iso05_ntrk() { return cms2.mus_iso05_ntrk(); }
	vector<int> &mus_lostHits() { return cms2.mus_lostHits(); }
	vector<int> &mus_nOverlaps() { return cms2.mus_nOverlaps(); }
	vector<int> &mus_nmatches() { return cms2.mus_nmatches(); }
	vector<int> &mus_overlap0() { return cms2.mus_overlap0(); }
	vector<int> &mus_overlap1() { return cms2.mus_overlap1(); }
	vector<int> &mus_pid_TM2DCompatibilityLoose() { return cms2.mus_pid_TM2DCompatibilityLoose(); }
	vector<int> &mus_pid_TM2DCompatibilityTight() { return cms2.mus_pid_TM2DCompatibilityTight(); }
	vector<int> &mus_pid_TMLastStationLoose() { return cms2.mus_pid_TMLastStationLoose(); }
	vector<int> &mus_pid_TMLastStationTight() { return cms2.mus_pid_TMLastStationTight(); }
	vector<int> &mus_sta_validHits() { return cms2.mus_sta_validHits(); }
	vector<int> &mus_timeDirection() { return cms2.mus_timeDirection(); }
	vector<int> &mus_timeNumStationsUsed() { return cms2.mus_timeNumStationsUsed(); }
	vector<int> &mus_trk_charge() { return cms2.mus_trk_charge(); }
	vector<int> &mus_trkidx() { return cms2.mus_trkidx(); }
	vector<int> &mus_type() { return cms2.mus_type(); }
	vector<int> &mus_validHits() { return cms2.mus_validHits(); }
	vector<int> &els_pat_genID() { return cms2.els_pat_genID(); }
	vector<int> &els_pat_genMotherID() { return cms2.els_pat_genMotherID(); }
	vector<int> &jets_pat_genPartonMother_id() { return cms2.jets_pat_genPartonMother_id(); }
	vector<int> &jets_pat_genParton_id() { return cms2.jets_pat_genParton_id(); }
	vector<int> &jets_pat_jetIDLoose() { return cms2.jets_pat_jetIDLoose(); }
	vector<int> &jets_pat_jetIDLooseAOD() { return cms2.jets_pat_jetIDLooseAOD(); }
	vector<int> &jets_pat_jetIDMinimal() { return cms2.jets_pat_jetIDMinimal(); }
	vector<int> &jets_pat_jetIDTight() { return cms2.jets_pat_jetIDTight(); }
	vector<int> &jets_pat_partonFlavour() { return cms2.jets_pat_partonFlavour(); }
	vector<int> &mus_pat_genID() { return cms2.mus_pat_genID(); }
	vector<int> &mus_pat_genMotherID() { return cms2.mus_pat_genMotherID(); }
	vector<int> &pfels_elsidx() { return cms2.pfels_elsidx(); }
	vector<int> &pfels_charge() { return cms2.pfels_charge(); }
	vector<int> &pfels_flag() { return cms2.pfels_flag(); }
	vector<int> &pfels_particleId() { return cms2.pfels_particleId(); }
	vector<int> &pfjets_chargedMultiplicity() { return cms2.pfjets_chargedMultiplicity(); }
	vector<int> &pfjets_muonMultiplicity() { return cms2.pfjets_muonMultiplicity(); }
	vector<int> &pfjets_neutralMultiplicity() { return cms2.pfjets_neutralMultiplicity(); }
	vector<int> &pfmus_musidx() { return cms2.pfmus_musidx(); }
	vector<int> &pfmus_charge() { return cms2.pfmus_charge(); }
	vector<int> &pfmus_flag() { return cms2.pfmus_flag(); }
	vector<int> &pfmus_particleId() { return cms2.pfmus_particleId(); }
	vector<int> &photons_fiduciality() { return cms2.photons_fiduciality(); }
	vector<int> &photons_scindex() { return cms2.photons_scindex(); }
	vector<int> &scs_detIdSeed() { return cms2.scs_detIdSeed(); }
	vector<int> &scs_elsidx() { return cms2.scs_elsidx(); }
	vector<int> &scs_severitySeed() { return cms2.scs_severitySeed(); }
	vector<int> &mus_tcmet_flag() { return cms2.mus_tcmet_flag(); }
	vector<int> &evt35X_mus_tcmet_flag() { return cms2.evt35X_mus_tcmet_flag(); }
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
	vector<int> &trks_elsidx() { return cms2.trks_elsidx(); }
	vector<int> &trk_musidx() { return cms2.trk_musidx(); }
	vector<int> &vtxs_isFake() { return cms2.vtxs_isFake(); }
	vector<int> &vtxs_isValid() { return cms2.vtxs_isValid(); }
	vector<int> &vtxs_tracksSize() { return cms2.vtxs_tracksSize(); }
	vector<vector<int> > &genps_lepdaughter_id() { return cms2.genps_lepdaughter_id(); }
	vector<vector<int> > &genps_lepdaughter_idx() { return cms2.genps_lepdaughter_idx(); }
	vector<vector<int> > &hlt1e31_trigObjs_id() { return cms2.hlt1e31_trigObjs_id(); }
	vector<vector<int> > &hlt8e29_trigObjs_id() { return cms2.hlt8e29_trigObjs_id(); }
	vector<vector<int> > &hlt_trigObjs_id() { return cms2.hlt_trigObjs_id(); }
	vector<vector<int> > &hyp_jets_idx() { return cms2.hyp_jets_idx(); }
	vector<vector<int> > &hyp_other_jets_idx() { return cms2.hyp_other_jets_idx(); }
	unsigned int &evt_nels() { return cms2.evt_nels(); }
	unsigned int &evt_detectorStatus() { return cms2.evt_detectorStatus(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	unsigned int &genps_flavorHistoryFilterResult() { return cms2.genps_flavorHistoryFilterResult(); }
	unsigned int &evt_ngenjets() { return cms2.evt_ngenjets(); }
	unsigned int &genps_signalProcessID() { return cms2.genps_signalProcessID(); }
	unsigned int &hlt1e31_bits1() { return cms2.hlt1e31_bits1(); }
	unsigned int &hlt1e31_bits2() { return cms2.hlt1e31_bits2(); }
	unsigned int &hlt1e31_bits3() { return cms2.hlt1e31_bits3(); }
	unsigned int &hlt1e31_bits4() { return cms2.hlt1e31_bits4(); }
	unsigned int &hlt1e31_bits5() { return cms2.hlt1e31_bits5(); }
	unsigned int &hlt1e31_bits6() { return cms2.hlt1e31_bits6(); }
	unsigned int &hlt1e31_bits7() { return cms2.hlt1e31_bits7(); }
	unsigned int &hlt1e31_bits8() { return cms2.hlt1e31_bits8(); }
	unsigned int &hlt8e29_bits1() { return cms2.hlt8e29_bits1(); }
	unsigned int &hlt8e29_bits2() { return cms2.hlt8e29_bits2(); }
	unsigned int &hlt8e29_bits3() { return cms2.hlt8e29_bits3(); }
	unsigned int &hlt8e29_bits4() { return cms2.hlt8e29_bits4(); }
	unsigned int &hlt8e29_bits5() { return cms2.hlt8e29_bits5(); }
	unsigned int &hlt8e29_bits6() { return cms2.hlt8e29_bits6(); }
	unsigned int &hlt8e29_bits7() { return cms2.hlt8e29_bits7(); }
	unsigned int &hlt8e29_bits8() { return cms2.hlt8e29_bits8(); }
	unsigned int &hlt_bits1() { return cms2.hlt_bits1(); }
	unsigned int &hlt_bits2() { return cms2.hlt_bits2(); }
	unsigned int &hlt_bits3() { return cms2.hlt_bits3(); }
	unsigned int &hlt_bits4() { return cms2.hlt_bits4(); }
	unsigned int &hlt_bits5() { return cms2.hlt_bits5(); }
	unsigned int &hlt_bits6() { return cms2.hlt_bits6(); }
	unsigned int &hlt_bits7() { return cms2.hlt_bits7(); }
	unsigned int &hlt_bits8() { return cms2.hlt_bits8(); }
	unsigned int &evt_njets() { return cms2.evt_njets(); }
	unsigned int &evt_njpts() { return cms2.evt_njpts(); }
	unsigned int &l1_bits1() { return cms2.l1_bits1(); }
	unsigned int &l1_bits2() { return cms2.l1_bits2(); }
	unsigned int &l1_bits3() { return cms2.l1_bits3(); }
	unsigned int &l1_bits4() { return cms2.l1_bits4(); }
	unsigned int &l1_techbits1() { return cms2.l1_techbits1(); }
	unsigned int &l1_techbits2() { return cms2.l1_techbits2(); }
	unsigned int &evt_nphotons() { return cms2.evt_nphotons(); }
	unsigned int &evt_ecalRecoStatus() { return cms2.evt_ecalRecoStatus(); }
	unsigned int &evt_nscs() { return cms2.evt_nscs(); }
	unsigned int &evt_ntrkjets() { return cms2.evt_ntrkjets(); }
	unsigned int &evt_nvtxs() { return cms2.evt_nvtxs(); }
	vector<unsigned int> &hlt1e31_prescales() { return cms2.hlt1e31_prescales(); }
	vector<unsigned int> &hlt8e29_prescales() { return cms2.hlt8e29_prescales(); }
	vector<unsigned int> &hlt_prescales() { return cms2.hlt_prescales(); }
	vector<unsigned int> &hyp_quadlep_bucket() { return cms2.hyp_quadlep_bucket(); }
	vector<unsigned int> &hyp_quadlep_first_index() { return cms2.hyp_quadlep_first_index(); }
	vector<unsigned int> &hyp_quadlep_fourth_index() { return cms2.hyp_quadlep_fourth_index(); }
	vector<unsigned int> &hyp_quadlep_second_index() { return cms2.hyp_quadlep_second_index(); }
	vector<unsigned int> &hyp_quadlep_third_index() { return cms2.hyp_quadlep_third_index(); }
	vector<unsigned int> &hyp_trilep_bucket() { return cms2.hyp_trilep_bucket(); }
	vector<unsigned int> &hyp_trilep_first_index() { return cms2.hyp_trilep_first_index(); }
	vector<unsigned int> &hyp_trilep_second_index() { return cms2.hyp_trilep_second_index(); }
	vector<unsigned int> &hyp_trilep_third_index() { return cms2.hyp_trilep_third_index(); }
	vector<unsigned int> &l1_prescales() { return cms2.l1_prescales(); }
	vector<unsigned int> &l1_techtrigprescales() { return cms2.l1_techtrigprescales(); }
	vector<unsigned int> &els_pat_flag() { return cms2.els_pat_flag(); }
	vector<unsigned int> &jets_pat_flag() { return cms2.jets_pat_flag(); }
	vector<unsigned int> &mus_pat_flag() { return cms2.mus_pat_flag(); }
	int &evt_nEvts() { return cms2.evt_nEvts(); }
	float &evt_filt_eff() { return cms2.evt_filt_eff(); }
	bool passHLTTrigger(TString trigName) { return cms2.passHLTTrigger(trigName); }
	bool passHLT8E29Trigger(TString trigName) { return cms2.passHLT8E29Trigger(trigName); }
	bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }
}
