#include "CMS2.h"

// this is the global fellow who's used extern everywhere
CMS2 cms2;
namespace tas {
	vector<TString> &evt_HLT_trigNames() { return cms2.evt_HLT_trigNames(); }
	vector<TString> &evt_L1_trigNames() { return cms2.evt_L1_trigNames(); }
	double &genps_pthat() { return cms2.genps_pthat(); }
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &scs_pos() { return cms2.scs_pos(); }
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &scs_vtx() { return cms2.scs_vtx(); }
	vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::DefaultCoordinateSystemTag> > &vtxs_position() { return cms2.vtxs_position(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &evt_bsp4() { return cms2.evt_bsp4(); }
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  &l1met_p4() { return cms2.l1met_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &alltrkjets_p4() { return cms2.alltrkjets_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc3_p4() { return cms2.els_mc3_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_mc_p4() { return cms2.els_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_gp_p4() { return cms2.jets_mc_gp_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_mc_p4() { return cms2.jets_mc_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_mc3_p4() { return cms2.mus_mc3_p4(); }
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
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_p4() { return cms2.mus_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_trk_p4() { return cms2.mus_trk_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_vertex_p4() { return cms2.mus_vertex_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genMotherP4() { return cms2.els_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &els_pat_genP4() { return cms2.els_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genJet_p4() { return cms2.jets_pat_genJet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genPartonMother_p4() { return cms2.jets_pat_genPartonMother_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_genParton_p4() { return cms2.jets_pat_genParton_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_p4() { return cms2.jets_pat_jet_p4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &jets_pat_jet_uncorp4() { return cms2.jets_pat_jet_uncorp4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genMotherP4() { return cms2.mus_pat_genMotherP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &mus_pat_genP4() { return cms2.mus_pat_genP4(); }
	vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > &scs_p4() { return cms2.scs_p4(); }
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
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_gp_p4() { return cms2.hyp_jets_mc_gp_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_mc_p4() { return cms2.hyp_jets_mc_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_p4() { return cms2.hyp_jets_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genJet_p4() { return cms2.hyp_jets_pat_genJet_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genPartonMother_p4() { return cms2.hyp_jets_pat_genPartonMother_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_genParton_p4() { return cms2.hyp_jets_pat_genParton_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_jets_pat_jet_p4() { return cms2.hyp_jets_pat_jet_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_gp_p4() { return cms2.hyp_other_jets_mc_gp_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_mc_p4() { return cms2.hyp_other_jets_mc_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_p4() { return cms2.hyp_other_jets_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genJet_p4() { return cms2.hyp_other_jets_pat_genJet_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genPartonMother_p4() { return cms2.hyp_other_jets_pat_genPartonMother_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_genParton_p4() { return cms2.hyp_other_jets_pat_genParton_p4(); }
	vector<vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > &hyp_other_jets_pat_jet_p4() { return cms2.hyp_other_jets_pat_jet_p4(); }
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
	float &gen_met() { return cms2.gen_met(); }
	float &gen_metPhi() { return cms2.gen_metPhi(); }
	float &evt_bField() { return cms2.evt_bField(); }
	float &evt_kfactor() { return cms2.evt_kfactor(); }
	float &evt_weight() { return cms2.evt_weight(); }
	float &evt_xsec_excl() { return cms2.evt_xsec_excl(); }
	float &evt_xsec_incl() { return cms2.evt_xsec_incl(); }
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
	float &evt_metNoHF() { return cms2.evt_metNoHF(); }
	float &evt_metNoHFHO() { return cms2.evt_metNoHFHO(); }
	float &evt_metNoHFHOPhi() { return cms2.evt_metNoHFHOPhi(); }
	float &evt_metNoHFHOSig() { return cms2.evt_metNoHFHOSig(); }
	float &evt_metNoHFPhi() { return cms2.evt_metNoHFPhi(); }
	float &evt_metSig() { return cms2.evt_metSig(); }
	float &evt_metOpt() { return cms2.evt_metOpt(); }
	float &evt_metOptHO() { return cms2.evt_metOptHO(); }
	float &evt_metOptHOPhi() { return cms2.evt_metOptHOPhi(); }
	float &evt_metOptHOSig() { return cms2.evt_metOptHOSig(); }
	float &evt_metOptNoHF() { return cms2.evt_metOptNoHF(); }
	float &evt_metOptNoHFHO() { return cms2.evt_metOptNoHFHO(); }
	float &evt_metOptNoHFHOPhi() { return cms2.evt_metOptNoHFHOPhi(); }
	float &evt_metOptNoHFHOSig() { return cms2.evt_metOptNoHFHOSig(); }
	float &evt_metOptNoHFPhi() { return cms2.evt_metOptNoHFPhi(); }
	float &evt_metOptSig() { return cms2.evt_metOptSig(); }
	float &evt_metOptPhi() { return cms2.evt_metOptPhi(); }
	float &evt_metPhi() { return cms2.evt_metPhi(); }
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
	float &evt_tcmet() { return cms2.evt_tcmet(); }
	float &evt_tcmetPhi() { return cms2.evt_tcmetPhi(); }
	float &evt_tcsumet() { return cms2.evt_tcsumet(); }
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
	vector<float> &els_e2x5Max() { return cms2.els_e2x5Max(); }
	vector<float> &els_e3x3() { return cms2.els_e3x3(); }
	vector<float> &els_e5x5() { return cms2.els_e5x5(); }
	vector<float> &els_eMax() { return cms2.els_eMax(); }
	vector<float> &els_eOverPIn() { return cms2.els_eOverPIn(); }
	vector<float> &els_eSC() { return cms2.els_eSC(); }
	vector<float> &els_eSCPresh() { return cms2.els_eSCPresh(); }
	vector<float> &els_eSCRaw() { return cms2.els_eSCRaw(); }
	vector<float> &els_eSeed() { return cms2.els_eSeed(); }
	vector<float> &els_eSeedOverPOut() { return cms2.els_eSeedOverPOut(); }
	vector<float> &els_etaErr() { return cms2.els_etaErr(); }
	vector<float> &els_fBrem() { return cms2.els_fBrem(); }
	vector<float> &els_hOverE() { return cms2.els_hOverE(); }
	vector<float> &els_ndof() { return cms2.els_ndof(); }
	vector<float> &els_outerEta() { return cms2.els_outerEta(); }
	vector<float> &els_outerPhi() { return cms2.els_outerPhi(); }
	vector<float> &els_phiErr() { return cms2.els_phiErr(); }
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
	vector<float> &hyp_ll_chi2() { return cms2.hyp_ll_chi2(); }
	vector<float> &hyp_ll_d0() { return cms2.hyp_ll_d0(); }
	vector<float> &hyp_ll_d0Err() { return cms2.hyp_ll_d0Err(); }
	vector<float> &hyp_ll_d0corr() { return cms2.hyp_ll_d0corr(); }
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
	vector<float> &hyp_met() { return cms2.hyp_met(); }
	vector<float> &hyp_metAll() { return cms2.hyp_metAll(); }
	vector<float> &hyp_metAllCaloExp() { return cms2.hyp_metAllCaloExp(); }
	vector<float> &hyp_metCaloExp() { return cms2.hyp_metCaloExp(); }
	vector<float> &hyp_metCone() { return cms2.hyp_metCone(); }
	vector<float> &hyp_metDPhiJet10() { return cms2.hyp_metDPhiJet10(); }
	vector<float> &hyp_metDPhiJet15() { return cms2.hyp_metDPhiJet15(); }
	vector<float> &hyp_metDPhiJet20() { return cms2.hyp_metDPhiJet20(); }
	vector<float> &hyp_metDPhiTrk10() { return cms2.hyp_metDPhiTrk10(); }
	vector<float> &hyp_metDPhiTrk25() { return cms2.hyp_metDPhiTrk25(); }
	vector<float> &hyp_metDPhiTrk50() { return cms2.hyp_metDPhiTrk50(); }
	vector<float> &hyp_metJes10() { return cms2.hyp_metJes10(); }
	vector<float> &hyp_metJes15() { return cms2.hyp_metJes15(); }
	vector<float> &hyp_metJes30() { return cms2.hyp_metJes30(); }
	vector<float> &hyp_metJes5() { return cms2.hyp_metJes5(); }
	vector<float> &hyp_metJes50() { return cms2.hyp_metJes50(); }
	vector<float> &hyp_metNoCalo() { return cms2.hyp_metNoCalo(); }
	vector<float> &hyp_metPhi() { return cms2.hyp_metPhi(); }
	vector<float> &hyp_metPhiAll() { return cms2.hyp_metPhiAll(); }
	vector<float> &hyp_metPhiAllCaloExp() { return cms2.hyp_metPhiAllCaloExp(); }
	vector<float> &hyp_metPhiCaloExp() { return cms2.hyp_metPhiCaloExp(); }
	vector<float> &hyp_metPhiCone() { return cms2.hyp_metPhiCone(); }
	vector<float> &hyp_metPhiJes10() { return cms2.hyp_metPhiJes10(); }
	vector<float> &hyp_metPhiJes15() { return cms2.hyp_metPhiJes15(); }
	vector<float> &hyp_metPhiJes30() { return cms2.hyp_metPhiJes30(); }
	vector<float> &hyp_metPhiJes5() { return cms2.hyp_metPhiJes5(); }
	vector<float> &hyp_metPhiJes50() { return cms2.hyp_metPhiJes50(); }
	vector<float> &hyp_metPhiNoCalo() { return cms2.hyp_metPhiNoCalo(); }
	vector<float> &hyp_quadlep_met() { return cms2.hyp_quadlep_met(); }
	vector<float> &hyp_quadlep_metAll() { return cms2.hyp_quadlep_metAll(); }
	vector<float> &hyp_trilep_met() { return cms2.hyp_trilep_met(); }
	vector<float> &hyp_trilep_metAll() { return cms2.hyp_trilep_metAll(); }
	vector<float> &jets_EMFcor() { return cms2.jets_EMFcor(); }
	vector<float> &jets_chFrac() { return cms2.jets_chFrac(); }
	vector<float> &jets_cor() { return cms2.jets_cor(); }
	vector<float> &jets_emFrac() { return cms2.jets_emFrac(); }
	vector<float> &jpts_chFrac() { return cms2.jpts_chFrac(); }
	vector<float> &jpts_cor() { return cms2.jpts_cor(); }
	vector<float> &jpts_emFrac() { return cms2.jpts_emFrac(); }
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
	vector<float> &trks_chi2() { return cms2.trks_chi2(); }
	vector<float> &trks_d0() { return cms2.trks_d0(); }
	vector<float> &trks_d0Err() { return cms2.trks_d0Err(); }
	vector<float> &trks_d0corr() { return cms2.trks_d0corr(); }
	vector<float> &trks_d0corrPhi() { return cms2.trks_d0corrPhi(); }
	vector<float> &trks_etaErr() { return cms2.trks_etaErr(); }
	vector<float> &trks_ndof() { return cms2.trks_ndof(); }
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
	vector<float> &ww_pmet() { return cms2.ww_pmet(); }
	vector<vector<float> > &hyp_jets_EMFcor() { return cms2.hyp_jets_EMFcor(); }
	vector<vector<float> > &hyp_jets_chFrac() { return cms2.hyp_jets_chFrac(); }
	vector<vector<float> > &hyp_jets_cor() { return cms2.hyp_jets_cor(); }
	vector<vector<float> > &hyp_jets_emFrac() { return cms2.hyp_jets_emFrac(); }
	vector<vector<float> > &hyp_jets_mc_emEnergy() { return cms2.hyp_jets_mc_emEnergy(); }
	vector<vector<float> > &hyp_jets_mc_hadEnergy() { return cms2.hyp_jets_mc_hadEnergy(); }
	vector<vector<float> > &hyp_jets_mc_invEnergy() { return cms2.hyp_jets_mc_invEnergy(); }
	vector<vector<float> > &hyp_jets_mc_otherEnergy() { return cms2.hyp_jets_mc_otherEnergy(); }
	vector<vector<float> > &hyp_jets_pat_bCorrF() { return cms2.hyp_jets_pat_bCorrF(); }
	vector<vector<float> > &hyp_jets_pat_cCorrF() { return cms2.hyp_jets_pat_cCorrF(); }
	vector<vector<float> > &hyp_jets_pat_gluCorrF() { return cms2.hyp_jets_pat_gluCorrF(); }
	vector<vector<float> > &hyp_jets_pat_jetCharge() { return cms2.hyp_jets_pat_jetCharge(); }
	vector<vector<float> > &hyp_jets_pat_noCorrF() { return cms2.hyp_jets_pat_noCorrF(); }
	vector<vector<float> > &hyp_jets_pat_udsCorrF() { return cms2.hyp_jets_pat_udsCorrF(); }
	vector<vector<float> > &hyp_other_jets_EMFcor() { return cms2.hyp_other_jets_EMFcor(); }
	vector<vector<float> > &hyp_other_jets_chFrac() { return cms2.hyp_other_jets_chFrac(); }
	vector<vector<float> > &hyp_other_jets_cor() { return cms2.hyp_other_jets_cor(); }
	vector<vector<float> > &hyp_other_jets_emFrac() { return cms2.hyp_other_jets_emFrac(); }
	vector<vector<float> > &hyp_other_jets_mc_emEnergy() { return cms2.hyp_other_jets_mc_emEnergy(); }
	vector<vector<float> > &hyp_other_jets_mc_hadEnergy() { return cms2.hyp_other_jets_mc_hadEnergy(); }
	vector<vector<float> > &hyp_other_jets_mc_invEnergy() { return cms2.hyp_other_jets_mc_invEnergy(); }
	vector<vector<float> > &hyp_other_jets_mc_otherEnergy() { return cms2.hyp_other_jets_mc_otherEnergy(); }
	vector<vector<float> > &hyp_other_jets_pat_bCorrF() { return cms2.hyp_other_jets_pat_bCorrF(); }
	vector<vector<float> > &hyp_other_jets_pat_cCorrF() { return cms2.hyp_other_jets_pat_cCorrF(); }
	vector<vector<float> > &hyp_other_jets_pat_gluCorrF() { return cms2.hyp_other_jets_pat_gluCorrF(); }
	vector<vector<float> > &hyp_other_jets_pat_jetCharge() { return cms2.hyp_other_jets_pat_jetCharge(); }
	vector<vector<float> > &hyp_other_jets_pat_noCorrF() { return cms2.hyp_other_jets_pat_noCorrF(); }
	vector<vector<float> > &hyp_other_jets_pat_udsCorrF() { return cms2.hyp_other_jets_pat_udsCorrF(); }
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
	vector<int> &els_looseId() { return cms2.els_looseId(); }
	vector<int> &els_lostHits() { return cms2.els_lostHits(); }
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
	vector<int> &ww_isdyee() { return cms2.ww_isdyee(); }
	vector<int> &ww_isdymm() { return cms2.ww_isdymm(); }
	vector<int> &ww_isdytt() { return cms2.ww_isdytt(); }
	vector<int> &ww_llgoodel() { return cms2.ww_llgoodel(); }
	vector<int> &ww_llgoodeliso() { return cms2.ww_llgoodeliso(); }
	vector<int> &ww_llgoodmu() { return cms2.ww_llgoodmu(); }
	vector<int> &ww_llgoodmuiso() { return cms2.ww_llgoodmuiso(); }
	vector<int> &ww_ltgoodel() { return cms2.ww_ltgoodel(); }
	vector<int> &ww_ltgoodeliso() { return cms2.ww_ltgoodeliso(); }
	vector<int> &ww_ltgoodmu() { return cms2.ww_ltgoodmu(); }
	vector<int> &ww_ltgoodmuiso() { return cms2.ww_ltgoodmuiso(); }
	vector<int> &ww_oppsign() { return cms2.ww_oppsign(); }
	vector<int> &ww_pass2met() { return cms2.ww_pass2met(); }
	vector<int> &ww_pass4met() { return cms2.ww_pass4met(); }
	vector<int> &ww_passFebsel_no_jetveto() { return cms2.ww_passFebsel_no_jetveto(); }
	vector<int> &ww_passFebsel_with_jetveto() { return cms2.ww_passFebsel_with_jetveto(); }
	vector<int> &ww_passaddzveto() { return cms2.ww_passaddzveto(); }
	vector<int> &ww_passzveto() { return cms2.ww_passzveto(); }
	vector<vector<int> > &hyp_jets_mc_id() { return cms2.hyp_jets_mc_id(); }
	vector<vector<int> > &hyp_jets_pat_genPartonMother_id() { return cms2.hyp_jets_pat_genPartonMother_id(); }
	vector<vector<int> > &hyp_jets_pat_genParton_id() { return cms2.hyp_jets_pat_genParton_id(); }
	vector<vector<int> > &hyp_jets_pat_partonFlavour() { return cms2.hyp_jets_pat_partonFlavour(); }
	vector<vector<int> > &hyp_other_jets_mc_id() { return cms2.hyp_other_jets_mc_id(); }
	vector<vector<int> > &hyp_other_jets_pat_genPartonMother_id() { return cms2.hyp_other_jets_pat_genPartonMother_id(); }
	vector<vector<int> > &hyp_other_jets_pat_genParton_id() { return cms2.hyp_other_jets_pat_genParton_id(); }
	vector<vector<int> > &hyp_other_jets_pat_partonFlavour() { return cms2.hyp_other_jets_pat_partonFlavour(); }
	vector<vector<int> > &hyp_quadlep_jets_index() { return cms2.hyp_quadlep_jets_index(); }
	vector<vector<int> > &hyp_trilep_jets_index() { return cms2.hyp_trilep_jets_index(); }
	unsigned int &evt_nalltrkjets() { return cms2.evt_nalltrkjets(); }
	unsigned int &evt_nels() { return cms2.evt_nels(); }
	unsigned int &evt_event() { return cms2.evt_event(); }
	unsigned int &evt_lumiBlock() { return cms2.evt_lumiBlock(); }
	unsigned int &evt_run() { return cms2.evt_run(); }
	unsigned int &evt_ngenjets() { return cms2.evt_ngenjets(); }
	unsigned int &evt_njets() { return cms2.evt_njets(); }
	unsigned int &evt_njpts() { return cms2.evt_njpts(); }
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
	bool passHLTTrigger(TString trigName) { return cms2.passHLTTrigger(trigName); }
	bool passL1Trigger(TString trigName) { return cms2.passL1Trigger(trigName); }
}
