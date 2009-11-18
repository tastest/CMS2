
//#include "TH2F.h"
#include "plotResults.h"
#include "Tools/HistogramUtilities.h"
//#include "Tools/DataSource.h"
#include "Tools/Utilities.h"
#include "Tools/histtools.cc"

#include "TROOT.h"
//#include "TGraphAsymmErrors.h"

using namespace std;



// for 2_2_X
//const static sources_t &theSources = theSources_22X;

TString rootfilename = "wzratio.root";


void plotResults() {

  //HistogramUtilities* h1 = new HistogramUtilities("Results.root");  //normalization is in weight
  HistogramUtilities* h1 = new HistogramUtilities("ABCDResults.root");  //normalization is in weight

  vector<DataSource> vSources;
  vSources.push_back(   fH_WEJET_ALP()   );	
  vSources.push_back(   fH_WMJET_ALP()   );	
  vSources.push_back(   fH_WTJET_ALP()   );	
  vSources.push_back(	fH_ZEEJET_ALP()	);
  vSources.push_back(   fH_ZMMJET_ALP() );
  vSources.push_back(   fH_ZTTJET_ALP() );
  vSources.push_back(	fH_QCD30()	);
  vSources.push_back(   fH_QCD80()      );
  vSources.push_back(	fH_QCDEM()	);
  vSources.push_back(	fH_PHOTONJET()	);
  vSources.push_back(	fH_QCDBCTOE()	);
  vSources.push_back(   fH_MU15_SINGLE() );
  vSources.push_back(	fH_TTBAR()	);
  //sources_.push_back(   fH_TTBAR_SINGLE() ); //it's single, but not called that
  h1->setOrder(vSources);
  
  TLegend *lg_all = h1->getLegend(theSources, "lep_pt", "", "all");
  /*
  makeStack(h1, lg_all, theSources, "lep_pt",				"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_transmass",		"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_tcmet",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_calomet_muon",		"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_met_dphi",			"", "all", true);
  //isos
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "all", true);
  //els
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_nlep_nod0iso",		"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_nlep_nometiso",	"", "e", true);
  makeStack(h1, lg_all, theSources, "lep_conversions",		"", "e", true);
  //mus
  makeStack(h1, lg_all, theSources, "lep_trckIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_ecalIso",			"", "m", true);
  //makeStack(h1, lg_all, theSources, "lep_hcalIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_relIso",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "m", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "m", true);

  makeStack(h1, lg_all, theSources, "lep_nlep",				"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "lep_njet30",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_0_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_1_pt",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_pt",				"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_mass",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_tcmet",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_calomet_muon",	"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet20",			"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_njet30",			"", "all", true);
  */ /*
  makeStack(h1, lg_all, theSources, "dilep_mass_ll20",		"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_pt_mgen",		"", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_0_genpt", "", "all", true);
  makeStack(h1, lg_all, theSources, "dilep_1_genpt", "", "all", true);
  makeStack(h1, lg_all, theSources, "lep_genpt", "", "all", true);
	 */
  makeStack(h1, lg_all, theSources, "dilep_mass_ll20",		"", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_pt_mgen",		"", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_0_genpt", "", "ee", true);
  //makeStack(h1, lg_all, theSources, "dilep_1_genpt", "", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_0_genpt_mch", "", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_1_genpt_mch", "", "ee", true);
  makeStack(h1, lg_all, theSources, "dilep_genpt", "", "ee", true);
  makeStack(h1, lg_all, theSources, "lep_genpt", "", "e", true);
  makeStack(h1, lg_all, theSources, "lep_genpt_mch", "", "e", true);
  TH1F* dilep_genpt = h1->getHistogram( sigdSources, "dilep_genpt", "", "ee");
  //TH1F* dilep_1_genpt = h1->getHistogram( sigdSources, "dilep_1_genpt", "", "ee");
  TH1F* lep_genpt = h1->getHistogram( sigsSources, "lep_genpt", "", "e");
  TH1F* lep_genmet = h1->getHistogram( sigsSources, "lep_genmet", "", "e");
  TH1F* lep_accgenmet = h1->getHistogram( sigsSources, "lep_accgenmet", "", "e");

  TString opt = "colz";

  //ABCD
  //tcmet vs reliso
  //e+m
  /*
  TH2F* tcmet_reliso_sig = h1->get2dHistogram(sigsSources, "lep_tcMet_relIso", "", "all", 1, "_sig"); //all or e/mu?
  //tcmet_reliso_sig->SetName("tcMet_relIso_signal");
  tcmet_reliso_sig->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_sig->GetName()+".png");

  TH2F* tcmet_reliso_bkg = h1->get2dHistogram(bkgSources, "lep_tcMet_relIso", "", "all", 1, "_bkg"); //all or e/mu?
  //tcmet_reliso_bkg->SetName("tcMet_relIso_background");
  tcmet_reliso_bkg->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_bkg->GetName()+".png");
  */

  //e
  TH2F* tcmet_reliso_all_e = getTH2F(h1, theSources,  "lep_tcMet_relIso", "", "e", 1, "_alls", opt);
  TH2F* tcmet_reliso_sig_e = getTH2F(h1, sigsSources, "lep_tcMet_relIso", "", "e", 1, "_sigs", opt);
  //this one has neutrino acceptance: require nu.eta to be in 2.4 to see if changes shape more similiar to z lep pt
  //TH2F* tcmetacc_reliso_all_e = getTH2F(h1, theSources,  "lep_acc_tcMet_relIso", "", "e", 1, "_alls", opt);
  TH2F* tcmetacc_reliso_sig_e = getTH2F(h1, sigsSources, "lep_acc_tcMet_relIso", "", "e", 1, "_sigs", opt);
  TH2F* tcmet_reliso_bkg_e = getTH2F(h1,bkgSources, "lep_tcMet_relIso", "", "e", 1, "_bkg", opt);
  TH2F* tcmet_reliso_zfake_e = getTH2F(h1, sigdSources, "lep_tcMet_relIso", "", "e", 1, "_zbkg", opt);

  /*
  TH2F* tcmet_reliso_mu15_e = h1->get2dHistogram((1ll << H_MU15_SINGLE), "lep_tcMet_relIso", "", "e", 1, "_mu15");
  tcmet_reliso_mu15_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_mu15_e->GetName()+".png");

  TH2F* tcmet_reliso_qcdem_e = h1->get2dHistogram((1ll << H_QCDEM), "lep_tcMet_relIso", "", "e", 1, "_qcdem");
  tcmet_reliso_qcdem_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_qcdem_e->GetName()+".png");

  TH2F* tcmet_reliso_qcdbc_e = h1->get2dHistogram((1ll << H_QCDBCTOE), "lep_tcMet_relIso", "", "e", 1, "_qcdbc");
  tcmet_reliso_qcdbc_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_qcdbc_e->GetName()+".png");

  TH2F* tcmet_reliso_photn_e = h1->get2dHistogram((1ll << H_PHOTONJET), "lep_tcMet_relIso", "", "e", 1, "_photn");
  tcmet_reliso_photn_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_photn_e->GetName()+".png");

  TH2F* tcmet_reliso_ttbar_e = h1->get2dHistogram((1ll << H_TTBAR), "lep_tcMet_relIso", "", "e", 1, "_ttbar");
  tcmet_reliso_ttbar_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_ttbar_e->GetName()+".png");

  TH2F* tcmet_reliso_qcd30_e = h1->get2dHistogram((1ll << H_QCD30), "lep_tcMet_relIso", "", "e", 1, "_qcd30");
  tcmet_reliso_qcd30_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_qcd30_e->GetName()+".png");

  TH2F* tcmet_reliso_qcd80_e = h1->get2dHistogram((1ll << H_QCD80), "lep_tcMet_relIso", "", "e", 1, "_qcd80");
  tcmet_reliso_qcd80_e->Draw(opt);
  c1->SaveAs((TString)tcmet_reliso_qcd80_e->GetName()+".png");
  */
  //new Z
  TH2F* zll_pt_eta_all_e = getTH2F(h1, theSources, "dilep_ll_pt_eta", "", "ee", 1, "_all", opt); 
  TH2F* zll_pt_eta_z_e   = getTH2F(h1, sigdSources, "dilep_ll_pt_eta", "", "ee", 1, "_z", opt); 

  TH2F* zlt_pt_eta_all_e = getTH2F(h1, theSources, "dilep_lt_pt_eta", "", "ee", 1, "_all", opt); 
  //TH2F* zlt_pt_eta_z_e   = getTH2F(h1, sigdSources, "dilep_lt_pt_eta", "", "ee", 1, "_z", opt); 

  TH2F* zlt_ll20_pt_eta_all_e = getTH2F(h1, theSources, "dilep_lt_ll20_pt_eta", "", "ee", 1, "_all", opt); 
  //TH2F* zlt_ll20_pt_eta_z_e  = getTH2F(h1, sigdSources, "dilep_lt_ll20_pt_eta", "", "ee", 1, "_z", opt); 

  //TH2F* zlpt_reliso_all_e = getTH2F(h1, theSources, "dilep_lepPt_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zlpt_reliso_z_e = getTH2F(h1, sigdSources, "dilep_lepPt_relIso", "", "ee", 1, "_z", opt); 

  //TH2F* zlpt_reliso_scl_all_e = getTH2F(h1, theSources, "dilep_lepPt_Scl_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zlpt_reliso_scl_z_e = getTH2F(h1, sigdSources, "dilep_lepPt_Scl_relIso", "", "ee", 1, "_z", opt); 

  //TH2F* zmet_reliso_all_e = getTH2F(h1, theSources, "dilep_lepMet_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_z_e = getTH2F(h1, sigdSources, "dilep_lepMet_relIso", "", "ee", 1, "_z", opt); 

  //the ones that matter are below here
  TH2F* zmet_reliso_scl_all_e 			= getTH2F(h1, theSources, "dilep_lepMet_Scl_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_rscl_all_e 			= getTH2F(h1, theSources, "dilep_lepMet_rScl_relIso", "", "ee", 1, "_all", opt); 
  TH2F* zmet_reliso_scl_z_e 			= getTH2F(h1, sigdSources,"dilep_lepMet_Scl_relIso", "", "ee", 1, "_z", opt); 

  TH2F* zmet_reliso_sclt_all_e 			= getTH2F(h1, theSources, "dilep_lepMet_Scl_trth_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_rsclt_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_rScl_trth_relIso", "", "ee", 1, "_all", opt); 
  TH2F* zmet_reliso_sclt_z_e 			= getTH2F(h1, sigdSources, "dilep_lepMet_Scl_trth_relIso", "", "ee", 1, "_z", opt); 

  TH2F* zmet_reliso_scltm_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_Scl_tmas_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_rscltm_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_rScl_tmas_relIso", "", "ee", 1, "_all", opt); 
  TH2F* zmet_reliso_scltm_z_e 			= getTH2F(h1, sigdSources, "dilep_lepMet_Scl_tmas_relIso", "", "ee", 1, "_z", opt); 

  TH2F* zmet_reliso_scltmt_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_Scl_tmast_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_rscltmt_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_rScl_tmast_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_scltmt_z_e 			= getTH2F(h1, sigdSources, "dilep_lepMet_Scl_tmast_relIso", "", "ee", 1, "_z", opt); 

  TH2F* zmet_reliso_scltmtm_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_Scl_tmastmes_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_rscltmtm_all_e 		= getTH2F(h1, theSources, "dilep_lepMet_rScl_tmastmes_relIso", "", "ee", 1, "_all", opt); 
  //TH2F* zmet_reliso_scltmtm_z_e 		= getTH2F(h1, sigdSources, "dilep_lepMet_Scl_tmastmes_relIso", "", "ee", 1, "_z", opt); 


  TCanvas* c1 = new TCanvas();
  //save gen plots
  //dilep_0_genpt->Scale( lep_genpt->Integral()/dilep_0_genpt->Integral() );
  //dilep_0_genpt->SetLineColor(2); //red
  //dilep_0_genpt->Draw();
  lep_genpt->Draw();
  dilep_genpt->Scale( lep_genpt->Integral()/dilep_genpt->Integral() );
  dilep_genpt->SetLineColor(2); //red
  dilep_genpt->Draw("same");
  c1->SaveAs("Compare_genpts.png");

  //plot z lep pt from th2, met in w, zleppt+zmet.
  TH1D* proj_wmet 		= new TH1D( *(tcmet_reliso_sig_e		->ProjectionX( (TString)tcmet_reliso_sig_e->GetName()+"_projx", 0, 101)) );
  TH1D* proj_wmetnu 	= new TH1D( *(tcmetacc_reliso_sig_e		->ProjectionX( (TString)tcmetacc_reliso_sig_e->GetName()+"_projx", 0, 101)) );
  //TH1D* proj_zlpt 		= new TH1D( *(zlpt_reliso_all_e			->ProjectionX( (TString)zlpt_reliso_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmet 		= new TH1D( *(zmet_reliso_all_e			->ProjectionX( (TString)zmet_reliso_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zfakemet	= new TH1D( *(tcmet_reliso_zfake_e		->ProjectionX( (TString)tcmet_reliso_zfake_e->GetName()+"_projx", 0, 101)) );
  //the ones that matter are below here
  TH1D* proj_zmetsc 	= new TH1D( *(zmet_reliso_scl_all_e		->ProjectionX( (TString)zmet_reliso_scl_all_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsct 	= new TH1D( *(zmet_reliso_sclt_all_e	->ProjectionX( (TString)zmet_reliso_sclt_all_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsctm 	= new TH1D( *(zmet_reliso_scltm_all_e	->ProjectionX( (TString)zmet_reliso_scltm_all_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsctmt 	= new TH1D( *(zmet_reliso_scltmt_all_e	->ProjectionX( (TString)zmet_reliso_scltmt_all_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsctmtm	= new TH1D( *(zmet_reliso_scltmtm_all_e	->ProjectionX( (TString)zmet_reliso_scltmtm_all_e->GetName()+"_projx", 0, 101)) );  
  //same with 'refined' ('readjusted' 'retarded') scale
  //TH1D* proj_zmetrsc 	= new TH1D( *(zmet_reliso_rscl_all_e	->ProjectionX( (TString)zmet_reliso_rscl_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmetrsct 	= new TH1D( *(zmet_reliso_rsclt_all_e	->ProjectionX( (TString)zmet_reliso_rsclt_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmetrsctm 	= new TH1D( *(zmet_reliso_rscltm_all_e	->ProjectionX( (TString)zmet_reliso_rscltm_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmetrsctmt = new TH1D( *(zmet_reliso_rscltmt_all_e	->ProjectionX( (TString)zmet_reliso_rscltmt_all_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmetrsctmtm= new TH1D( *(zmet_reliso_rscltmtm_all_e->ProjectionX( (TString)zmet_reliso_rscltmtm_all_e->GetName()+"_projx", 0, 101)) );  

  //only in iso 0.1--iso signal region
  /*
  TH1D* proj_zlpt_is1	= new TH1D( *(zlpt_reliso_all_e			->ProjectionX( (TString)zlpt_reliso_all_e->GetName()+"_projx", 0, 10)) );  
  TH1D* proj_zmet_is1	= new TH1D( *(zmet_reliso_all_e			->ProjectionX( (TString)zmet_reliso_all_e->GetName()+"_projx", 0, 10)) );  
  TH1D* proj_zmetsc_is1	= new TH1D( *(zmet_reliso_scl_all_e		->ProjectionX( (TString)zmet_reliso_scl_all_e->GetName()+"_projx", 0, 10)) );  
  TH1D* proj_zmetsct_is1= new TH1D( *(zmet_reliso_sclt_all_e	->ProjectionX( (TString)zmet_reliso_sclt_all_e->GetName()+"_projx", 0, 10)) );  
  TH1D* proj_wmet_is1	= new TH1D( *(tcmet_reliso_sig_e		->ProjectionX( (TString)tcmet_reliso_sig_e->GetName()+"_projx", 0, 10)) );
  TH1D* proj_wmetnu_is1	= new TH1D( *(tcmetacc_reliso_sig_e		->ProjectionX( (TString)tcmetacc_reliso_sig_e->GetName()+"_projx", 0, 10)) );
  */
  //only signal for Zs
  //TH1D* proj_zlpt_sig	= new TH1D( *(zlpt_reliso_z_e			->ProjectionX( (TString)zlpt_reliso_z_e->GetName()+"_projx", 0, 101)) );  
  //TH1D* proj_zmet_sig	= new TH1D( *(zmet_reliso_z_e			->ProjectionX( (TString)zmet_reliso_z_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsc_sig	= new TH1D( *(zmet_reliso_scl_z_e		->ProjectionX( (TString)zmet_reliso_scl_z_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsct_sig= new TH1D( *(zmet_reliso_sclt_z_e  	->ProjectionX( (TString)zmet_reliso_sclt_z_e->GetName()+"_projx", 0, 101)) );  
  TH1D* proj_zmetsctm_sig= new TH1D( *(zmet_reliso_scltm_z_e  	->ProjectionX( (TString)zmet_reliso_scltm_z_e->GetName()+"_projx", 0, 101)) );  

  // ll pt vs eta, all samples -- ll eta when ll pt is 0-20
  TH1D* proj_ll_pteta_l20_alle = new TH1D( *(zll_pt_eta_all_e   ->ProjectionY( (TString)zll_pt_eta_all_e->GetName()+"_projy", 0, 20)) );
  proj_ll_pteta_l20_alle->Draw();
  c1->SaveAs( (TString)proj_ll_pteta_l20_alle->GetName()+".png" );

  // ll pt vs eta, z only -- ll eta when ll pt is 0-20
  TH1D* proj_ll_pteta_l20_sige = new TH1D( *(zll_pt_eta_z_e     ->ProjectionY( (TString)zll_pt_eta_z_e->GetName()+"_projy", 0, 20)) );
  proj_ll_pteta_l20_sige->Draw();
  c1->SaveAs( (TString)proj_ll_pteta_l20_sige->GetName()+".png" );

  // lt pt vs eta, all samples -- lt eta when lt pt is 0-20
  TH1D* proj_lt_pteta_l20_alle = new TH1D( *(zlt_pt_eta_all_e    ->ProjectionY( (TString)zlt_pt_eta_all_e->GetName()+"_projy", 0, 20)) );
  proj_lt_pteta_l20_alle->Draw();
  c1->SaveAs( (TString)proj_lt_pteta_l20_alle->GetName()+".png" );

  // lt pt vs eta, all samples -- lt pt all eta
  TH1D* proj_lt_pteta_alle = new TH1D( *(zlt_pt_eta_all_e    ->ProjectionX( (TString)zlt_pt_eta_all_e->GetName()+"_projx", 0, 101)) );
  proj_lt_pteta_alle->Draw();
  c1->SaveAs( (TString)proj_lt_pteta_alle->GetName()+".png" );

  // lt pt vs eta, all samples -- lt eta when ll and lt pt are 0-20
  TH1D* proj_lt_ll20_pteta_l20_alle = new TH1D( *(zlt_ll20_pt_eta_all_e->ProjectionY( (TString)zlt_ll20_pt_eta_all_e->GetName()+"_projy", 0, 20)) );
  proj_lt_ll20_pteta_l20_alle->Draw();
  c1->SaveAs( (TString)proj_lt_ll20_pteta_l20_alle->GetName()+".png" );
  
  // lt pt vs eta, all samples -- lt pt when ll pt is 0-20
  TH1D* proj_lt_ll20_pteta_alle = new TH1D( *(zlt_ll20_pt_eta_all_e->ProjectionX( (TString)zlt_ll20_pt_eta_all_e->GetName()+"_projx", 0, 101)) );
  proj_lt_ll20_pteta_alle->Draw();
  c1->SaveAs( (TString)proj_lt_ll20_pteta_alle->GetName()+".png" );

  //some simple hist stats
  //way this fn works is the last number printed is ratio of integral in bin range of last two args to total integral
  cout << "\nname\t\t\t\tmean\tintegral\t\% below 20\n";
  //printHistStats(proj_zlpt		,0,20);
  //printHistStats(proj_zmet		,0,20);
  //printHistStats(proj_zfakemet	,0,20);
  printHistStats(proj_zmetsc	,0,20);
  printHistStats(proj_zmetsct	,0,20);
  printHistStats(proj_zmetsctm	,0,20);
  printHistStats(proj_wmet		,0,20);
  printHistStats(proj_wmetnu	,0,20);
  printHistStats(dilep_genpt	,0,20);
  printHistStats(lep_genpt		,0,20);
  cout << endl << endl;

  //mc nu acceptance correction for z
  TH1F* zcorr = getAccCorr( proj_wmet, proj_wmetnu, proj_zmetsc );
  zcorr->Scale( proj_wmet->Integral()/zcorr->Integral() );
  zcorr->SetLineColor(2); //red
  zcorr->Draw();
  proj_wmet->Draw("same");
  c1->SaveAs("Compare_zcorr_wmet.png");

  //again with transverse mass cut on z
  TH1F* zcorrtm = getAccCorr( proj_wmet, proj_wmetnu, proj_zmetsctm );
  zcorrtm->SetLineColor(2); //red
  //use it before scaling/saving
  printHistStats(zcorrtm		,0,20);

  //plot correction hists
  TH1D* wcorrhist = new TH1D( *proj_wmet );
  wcorrhist->Divide( proj_wmetnu );
  wcorrhist->SetLineColor(2); //red
  wcorrhist->Draw();
  TH1F* wgencorrhist = new TH1F( *lep_genmet );
  wgencorrhist->Divide( lep_accgenmet );
  wgencorrhist->Draw("same");
  c1->SaveAs("Compare_corrfactors.png");

  //get th2s with correction--no transmass cut
  TH2F* zmetcorr_reliso_scl_all_e = corrTH2F( zmet_reliso_scl_all_e, proj_wmet, proj_wmetnu ); //the 2 last args are the hists for ratio
  zmetcorr_reliso_scl_all_e->Draw(opt);
  c1->SaveAs("dilep_lepmetcorr_Scl_reliso_ee_all.png"); //th2 alone
  //with tm cut
  TH2F* zmetcorr_reliso_scltm_all_e = corrTH2F( zmet_reliso_scltm_all_e, proj_wmet, proj_wmetnu );
  zmetcorr_reliso_scltm_all_e->Draw(opt);
  c1->SaveAs("dilep_lepmetcorr_Scl_tmas_reliso_ee_all.png"); //th2 alone
  //with tm from gen
  TH2F* zgenmetcorr_reliso_scltm_all_e = corrTH2F( zmet_reliso_scltm_all_e, lep_genmet, lep_accgenmet );
  zgenmetcorr_reliso_scltm_all_e->Draw(opt);
  c1->SaveAs("dilep_lepgenmetcorr_Scl_tmas_reliso_ee_all.png"); //th2 alone
  

  //compare projection after scaling to projection scaled--should be exact agreement (same scale for both)
  zmetcorr_reliso_scltm_all_e->ProjectionX( (TString)zmetcorr_reliso_scltm_all_e->GetName()+"_projx", 0, 101)->Draw();
  zcorrtm->Draw("same");
  c1->SaveAs("Compare_zcorrtmproj_wmet.png");

  //compare gen with w, z scaled--have the th2s, so can scale th1s freely
  TH1D* zgencorrtm = zgenmetcorr_reliso_scltm_all_e->ProjectionX( (TString)zgenmetcorr_reliso_scltm_all_e->GetName()+"_projx", 0, 101);
  cout << "zgencorrtm = " << zgencorrtm->GetName() << endl;
  printHistStats(zgencorrtm		,0,20);
  zgencorrtm->Scale( proj_wmet->Integral()/zgencorrtm->Integral() );
  zgencorrtm->Draw();
  //lep_genmet->Scale( proj_wmet->Integral()/lep_genmet->Integral() );
  //lep_genmet->SetLineColor( 6 ); //magenta
  //lep_genmet->Draw("same");
  //lep_accgenmet->Scale( proj_wmet->Integral()/lep_accgenmet->Integral() );
  //lep_accgenmet->SetLineColor( 7 ); //cyan (light blue)
  //lep_accgenmet->Draw("same");
  zcorrtm->Scale( proj_wmet->Integral()/zcorrtm->Integral() );
  zcorrtm->Draw("same");
  proj_wmet->SetLineColor(3); //green
  proj_wmet->Draw("same");
  c1->SaveAs("Compare_zcorrtm_wmet.png");




  //proj_zlpt->Scale( proj_wmet->Integral(0,101)/proj_zlpt->Integral(0,101) );
  //proj_zlpt->Scale( proj_wmetnu->Integral(0,101)/proj_zlpt->Integral(0,101) ); //scale to w dist with nu cut
  //proj_zlpt->SetStats(0);
  //proj_zlpt->Draw();
  //proj_zmet->Scale( proj_wmet->Integral()/proj_zmet->Integral() );
  //proj_zmet->Scale( proj_wmetnu->Integral()/proj_zmet->Integral() );
  //proj_zmet->SetLineColor(2); //red
  //proj_zmet->Draw("same");
  TH1D* proj_zmetsc2 = new TH1D( *proj_zmetsc );

  //1
  proj_zmetsc->Scale( proj_wmet->Integral()/proj_zmetsc->Integral() );
  proj_zmetsc->SetLineColor(2); //red
  proj_zmetsc->Draw();
  proj_wmet->SetLineColor(3); //green
  proj_wmet->Draw("same");
  c1->SaveAs("Compare_zmet_wmet.png");

  //2
  proj_zmetsc2->Scale( proj_wmetnu->Integral()/proj_zmetsc2->Integral() );
  proj_zmetsc2->SetLineColor(2); //red
  proj_zmetsc2->Draw();
  proj_wmetnu->SetLineColor(3); //green
  proj_wmetnu->Draw("same");
  c1->SaveAs("Compare_zmet_wmetnu.png");

  //3
  proj_zmetsctm->Scale( proj_wmet->Integral()/proj_zmetsctm->Integral() ); 
  proj_zmetsctm->SetLineColor(2); //red
  proj_zmetsctm->Draw();
  proj_wmet->Draw("same");
  c1->SaveAs("Compare_zmettm_wmet.png");

  //4
  proj_zmetsctmtm->Scale( proj_wmetnu->Integral()/proj_zmetsctmtm->Integral() ); 
  proj_zmetsctmtm->SetLineColor(2); //red
  proj_zmetsctmtm->Draw();
  proj_zmetsctmt->Scale( proj_wmetnu->Integral()/proj_zmetsctmt->Integral() ); 
  proj_zmetsctmt->SetLineColor(4); //blue
  proj_zmetsctmt->Draw("same");
  proj_wmetnu->Draw("same");
  c1->SaveAs("Compare_zmettmtm_wmetnu.png");

  /*
  proj_zmetsct->Scale( proj_wmetnu->Integral()/proj_zmetsct->Integral() ); 
  proj_zmetsct->Draw();
  proj_zmetsctm->Scale( proj_wmetnu->Integral()/proj_zmetsctm->Integral() );
  proj_zmetsctm->SetLineColor(4); //blue
  proj_zmetsctm->Draw("same");
  proj_zmetsc->Scale( proj_wmetnu->Integral()/proj_zmetsc->Integral() );
  proj_zmetsc->SetLineColor(2); //red
  proj_zmetsc->Draw("same");
  //proj_wmet->SetLineColor(4); //blue
  //proj_wmet->Draw("same");
  //proj_wmetnu->Scale( proj_wmet->Integral()/proj_wmetnu->Integral() ); //this is new baseline
  proj_wmetnu->SetLineColor(3); //green
  proj_wmetnu->Draw("same");
  //c1->SaveAs("Compare_zlpt_zmet_wmet.png");
  c1->SaveAs("Compare_zmett_zmet_wmet.png");
  */  

  //again, for signal samples only
  proj_zmetsct_sig->Scale( proj_wmetnu->Integral()/proj_zmetsct_sig->Integral() ); //new baseline
  proj_zmetsct_sig->Draw();
  proj_zmetsctm_sig->Scale( proj_wmetnu->Integral()/proj_zmetsctm_sig->Integral() ); //new baseline
  proj_zmetsctm_sig->SetLineColor(4); //blue
  proj_zmetsctm_sig->Draw();
  proj_zmetsc_sig ->Scale( proj_wmetnu->Integral()/proj_zmetsc_sig->Integral() );
  proj_zmetsc_sig ->SetLineColor(2); //red
  proj_zmetsc_sig ->Draw("same");
  proj_wmetnu->Draw("same"); //compare to same--this is just w signal
  c1->SaveAs("Compare_zsig_zmett_zmet_wmet.png");

  /*  //for iso 0.1
  proj_zmetsct_is1->Scale( proj_wmetnu->Integral()/proj_zmetsct_is1->Integral() ); //new baseline
  proj_zmetsct_is1->Draw();
  proj_zmetsc_is1 ->Scale( proj_wmetnu->Integral()/proj_zmetsc_is1->Integral() );
  proj_zmetsc_is1 ->SetLineColor(2); //red
  proj_zmetsc_is1 ->Draw("same");
  proj_wmetnu->Draw("same"); //compare to same--this is just w signal
  c1->SaveAs("Compare_iso1_zmett_zmet_wmet.png");
  */

  //ABCD results
  double metmax = 80.;

  TH2F* hists[] = {
	tcmet_reliso_all_e,
	tcmet_reliso_sig_e,
	tcmet_reliso_bkg_e,
	tcmet_reliso_zfake_e,
	//tcmet_reliso_mu15_e,
	//tcmet_reliso_qcdem_e,
	//tcmet_reliso_qcdbc_e,
	//tcmet_reliso_photn_e,
	//tcmet_reliso_ttbar_e,
	//tcmet_reliso_qcd30_e,
	//tcmet_reliso_qcd80_e,
  };
  const int N = sizeof(hists)/sizeof(TH2F*);

  /*  double x1 = 20.;
  double x2 = metmax;
  double x3 = 0.;
  double x4 = 20.;
  double y1 = 0.1;
  double y2 = 1.;
  double y3 = 0.;
  double y4 = 0.1; */
  //cout << abcd(tcmet_reliso_bkg, x1, x2, x3, x4, y1, y2, y3, y4) << endl;
  //cout << abcd(tcmet_reliso_bkg_e, x1, x2, x3, x4, y1, y2, y3, y4) << endl;
  //cout << abcd(tcmet_reliso_bkg_m, x1, x2, x3, x4, y1, y2, y3, y4) << endl;
  //cout << "Full range\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);
  
  //change x region only  //x4 = 15.;
  //cout << "5 border in x\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);

  //change y region only
  //x4 = 20.; //reset  //y1 = 0.2;  //cout << "0.1 border in y\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);

  //change x and y (y already changed)  //x4 = 15.;
  //cout << "5 border in x, 0.1 border in y\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);	

  //change met cut to 30  //x1 = 30.;  //x4 = 20.;  //y1 = 0.1;
  //cout << "met starts at 30, 10 border\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);

  //y2 = 0.2;
  //cout << "met starts at 30, 10 border, best y range\n";
  //Nabcd( hists, N, x1, x2, x3, x4, y1, y2, y3, y4);

  //lower iso upper limit, same met  
  cout << "Best iso range: 0-0.1, 0.1-0.2, full met no border\n";
  Nabcd( hists, N, 20., metmax, 0., 20., 0.1, 0.2, 0., 0.1);

  cout << "Iso range: 0-0.1, 0.3-0.4, increased met, 5 border\n";
  Nabcd( hists, N, 25., metmax, 0., 20., 0.3, 0.4, 0., 0.1);

  //cout << "exclude lowest iso bin: 0.01-0.1, 0.1-0.2, full met no border\n";
  //Nabcd( hists, N, 20., metmax, 0., 20., 0.1, 0.2, 0.01, 0.1);
  /*
  cout << "same width in iso as best, but with 0.1 border\n";
  Nabcd( hists, N, 20., metmax, 0., 20., 0.2, 0.3, 0., 0.1);

  cout << "Best iso range: 0-0.1, 0.1-0.2, change met range\n";
  Nabcd( hists, N, 20., metmax, 5., 15., 0.1, 0.2, 0., 0.1);

  cout << "Best iso range: 0-0.1, 0.1-0.2, change met range 2\n";
  Nabcd( hists, N, 20., metmax, 0., 15., 0.1, 0.2, 0., 0.1);

  cout << "Best iso range: 0-0.1, 0.1-0.2, change met range 3\n";
  Nabcd( hists, N, 20., metmax, 10., 20., 0.1, 0.2, 0., 0.1);

  cout << "iso border: 0-0.1, 0.2-0.3, and met range border\n";
  Nabcd( hists, N, 20., metmax, 5., 15., 0.2, 0.3, 0., 0.1);

  cout << "iso border: 0-0.1, 0.3-0.4, and met range border\n";
  Nabcd( hists, N, 20., metmax, 0., 10., 0.3, 0.4, 0., 0.1);

  cout << "lower iso, old met\n";
  Nabcd( hists, N, 20., metmax, 0., 20., 0.05, 0.1, 0., 0.05);

  cout << "lower iso w border, old met\n";
  Nabcd( hists, N, 20., metmax, 0., 20., 0.1, 0.15, 0., 0.05);

  cout << "lower iso w border, higher met\n";
  Nabcd( hists, N, 25., metmax, 15., 25., 0.1, 0.15, 0., 0.05);

  cout << "lower iso w border, higher met w border\n";
  Nabcd( hists, N, 25., metmax, 0., 15., 0.1, 0.15, 0., 0.05);
*/
  //projectX(tcmet_reliso_qcdem_e, 4);

  //cout << "\nUsing just lep pt as met\n";
  //aviCDtable( tcmet_reliso_all_e, zlpt_reliso_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  //make a third copy for each z category for cut on neutrino acceptance?

  cout << "\nUsing lep pt + met as met: baseline, orig regions\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_scl_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met: baseline, sig met 25, iso 0.3-0.4\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_scl_all_e, tcmet_reliso_sig_e, 25., metmax, -0.1, 20., 0.3, 0.4, 0., 0.1);
  /*
  cout << "\n******************************************\n";

  cout << "\nUsing lep pt + met corrected for nu: orig regions\n";
  aviCDtable( tcmet_reliso_all_e, zmetcorr_reliso_scl_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met corrected for nu: sig met 25, iso 0.3-0.4\n";
  aviCDtable( tcmet_reliso_all_e, zmetcorr_reliso_scl_all_e, tcmet_reliso_sig_e, 25., metmax, -0.1, 20., 0.3, 0.4, 0., 0.1);
  */
  cout << "\n******************************************\n";

  cout << "\nUsing lep pt + met corrected for nu w/ tmass: orig regions\n";
  aviCDtable( tcmet_reliso_all_e, zmetcorr_reliso_scltm_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met corrected for nu w/ tmass: sig met 25, iso 0.3-0.4\n";
  aviCDtable( tcmet_reliso_all_e, zmetcorr_reliso_scltm_all_e, tcmet_reliso_sig_e, 25., metmax, -0.1, 20., 0.3, 0.4, 0., 0.1);

  cout << "\n******************************************\n";

  cout << "\nUsing lep pt + genmet corrected for nu w/ tmass: orig regions\n";
  aviCDtable( tcmet_reliso_all_e, zgenmetcorr_reliso_scltm_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + genmet corrected for nu w/ tmass: sig met 25, iso 0.3-0.4\n";
  aviCDtable( tcmet_reliso_all_e, zgenmetcorr_reliso_scltm_all_e, tcmet_reliso_sig_e, 25., metmax, -0.1, 20., 0.3, 0.4, 0., 0.1);

  /*
  cout << "\nUsing lep pt + met as met: baseline rescaled\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_rscl_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z truth on ll\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_sclt_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z truth on ll rescale\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_rsclt_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_scltm_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut rescale\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_rscltm_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut + truth\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_scltmt_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut + truth rescale\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_rscltmt_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut + truth + measure pt to 10\%\n";
  aviCDtable( tcmet_reliso_all_e, zmet_reliso_scltmtm_all_e, tcmet_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);

  cout << "\nUsing lep pt + met as met, z tmass cut + truth + measure pt to 10\% rescale\n";
  aviCDtable( tcmetacc_reliso_all_e, zmet_reliso_rscltmtm_all_e, tcmetacc_reliso_sig_e, 20., metmax, -0.1, 20., 0.1, 0.2, 0., 0.1);
*/
}
//end plotResults

void aviCDtable( TH2F* data, TH2F* di, TH2F* ssig, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  double a = integrateTH2F(data, x1, x2, y3, y4);
  double b = integrateTH2F(data, x3, x4, y3, y4);
  double c = integrateTH2F(data, x1, x2, y1, y2);
  double d = integrateTH2F(data, x3, x4, y1, y2);
  double A = integrateTH2F(di, x1, x2, y3, y4);
  double B = integrateTH2F(di, x3, x4, y3, y4);
  double C = integrateTH2F(di, x1, x2, y1, y2);
  double D = integrateTH2F(di, x3, x4, y1, y2);
  double e = integrateTH2F(ssig, x1, x2, y3, y4);
  double f = integrateTH2F(ssig, x3, x4, y3, y4);
  double g = integrateTH2F(ssig, x1, x2, y1, y2);
  double h = integrateTH2F(ssig, x3, x4, y1, y2);
  //note: 80. is metmax above--same as upper edge of th2 in looper.cc
  double ditot = integrateTH2F(di, 0, 80., 0, 1);
  double sitot = integrateTH2F(ssig, 0, 80., 0, 1);

  double bpr = b - a*B/A;
  double cpr = c - a*C/A;
  double dpr = d - a*D/A;
  /*
  double bpr = b - a*B/ditot;
  double cpr = c - a*C/ditot;
  double dpr = d - a*D/ditot;
  */
  //next iteration--take the estimate of the bkg in a (bpr*cpr/dpr), subtract from a, and get bpr,cpr,dpr again, new ratio
  double a2 = a - bpr*cpr/dpr;
  double bpr2 = b - a2*B/A;
  double cpr2 = c - a2*C/A;
  double dpr2 = d - a2*D/A;

  cout << endl;
  cout <<        "& Data & Signal Only & Z Estimate & Data with subtraction & Subtraction iteration \\\\ \\hline \\hline" << endl;
  cout << " a & " << a << " & " << e << " & " << A << " & & "               << a2   << "\\\\\n";
  cout << " b & " << b << " & " << f << " & " << B << " & " << bpr << " & " << bpr2 << "\\\\\n";
  cout << " c & " << c << " & " << g << " & " << C << " & " << cpr << " & " << cpr2 << "\\\\\n";
  cout << " d & " << d << " & " << h << " & " << D << " & " << dpr << " & " << dpr2 << "\\\\\n";
  cout << "prediction for a & " << b*c/d << " & " << f*g/h << " & " << B*C/D << " & " << bpr*cpr/dpr << " & " << bpr2*cpr2/dpr2 << "\\\\" << endl << endl;

  //yanjun table
  cout << "sample    A/tot    B/tot   C/tot   D/tot\n";
  cout << "wjet     " << e/sitot << "  " << f/sitot << "  " << g/sitot << "  " << h/sitot << endl;
  cout << "zjet     " << A/ditot << "  " << B/ditot << "  " << C/ditot << "  " << D/ditot << endl << endl;
}
//end aviCDtable

template <class TH> TH2F* corrTH2F( TH2F* h, TH* met, TH* metacc) {
  TH2F* newh = new TH2F( *h );
  TH1F* ratio = new TH1F( *(TH1F*)met ); //make the correction hist
  ratio->Divide( metacc );
  for( int i=0;i<h->GetNbinsX();i++ ) { //x values: pt/met
	for( int j=0;j<h->GetNbinsY();j++ ) { //do for all y==iso
	  if( ratio->GetBinContent(i) < 5. ) //protect against large scale factors
		newh->SetBinContent( i, j, h->GetBinContent(i,j) * ratio->GetBinContent(i) );
	}
  }
  return newh;
}

template <class TH> TH1F* getAccCorr( TH* met, TH* metacc, TH* z ) { //returns th1f even though args are templated
//TH1F* accCorr( TH1F* met, TH1F* metacc, TH1F* z ) {
  TH1F* ratio = new TH1F( *(TH1F*)met );
  ratio->Divide( metacc );
  ratio->Multiply( z );
  return ratio;
}

void projectX(TH2F* h, const unsigned int n) {
  TCanvas* c = new TCanvas();
  TH1D* projs[n];
  int bins[5] = {0,1,10,20,100};
  //TH1D projs[n];
  //int nbins = h->GetNbinsY();
  cout << "i   mean      rms\n";
  for( unsigned int i=0; i<n; i++ ) {
	//projs[i] = new TH1D( *(h->ProjectionX( Form("%s_%i", h->GetName(), i), i*nbins/n, (i+1)*nbins/n)) );
	projs[i] = new TH1D( *(h->ProjectionX( Form("%s_%i", h->GetName(), i), bins[i], bins[i+1])) );
	cout << i << "  " << projs[i]->GetMean() << "  " << projs[i]->GetRMS() << endl;
	projs[i]->SetLineColor( i+1 );
	if( i == 0 )
	  projs[i]->Draw();
	else {
	  //scale to match 0th
	  projs[i]->Scale( projs[0]->Integral()/projs[i]->Integral() );
	  projs[i]->Draw("same");
	}
  }
  //projs[0]->Draw();
  c->SaveAs((TString)h->GetName()+"_proj.png");
}

void oneabcd(TH2F* h, TString title, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  cout << title << ":\n" ;
  double pred = abcd(h, x1, x2, x3, x4, y1, y2, y3, y4);
  double trth = integrateTH2F(h, x1, x2, y3, y4);
  cout << "b*c/d = " << pred << "   a = " << trth
	   << "    %diff = " << fabs(pred-trth)/trth << endl;
}

void Nabcd(TH2F** h, int n, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  printCoords(x1, x2, x3, x4, y1, y2, y3, y4);
  for( int i=0; i<n; i++ )
	oneabcd(h[i], h[i]->GetName(), x1, x2, x3, x4, y1, y2, y3, y4);
  //cout << "all  samples e:\nb*c/d = " << abcd(tcmet_reliso_all_e, x1, x2, x3, x4, y1, y2, y3, y4)
  //   << "   a = " << integrateTH2F(tcmet_reliso_all_e, x1, x2, y3, y4) << endl;
  cout << endl;
}


void printCoords(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  cout << "(x1, x2, x3, x4, y1, y2, y3, y4) = ("
	   << x1 << ", " << x2 << ", " << x3 << ", " << x4 << ", " << y1 << ", " << y2 << ", " << y3 << ", " << y4 << ")\n";
}

//currently, returns a = b*c/d, so call it such that this is what you want
double abcd(TH2F* h, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  //double a = integrateTH2F(h, x1, x2, y3, y4); //this is what you're predicting, dumb-ass
  double b = integrateTH2F(h, x3, x4, y3, y4);
  double c = integrateTH2F(h, x1, x2, y1, y2);
  double d = integrateTH2F(h, x3, x4, y1, y2);
  cout << b << "  " << c << "  " << d << "  ";
  return b*c/d;
}


double integrateTH2F(TH2F* h, double xlow, double xhgh, double ylow, double yhgh) {
  int binlowx = h->GetXaxis()->FindFixBin( xlow );
  int binhghx = h->GetXaxis()->FindFixBin( xhgh );
  int binlowy = h->GetYaxis()->FindFixBin( ylow );
  int binhghy = h->GetYaxis()->FindFixBin( yhgh );
  if( xhgh == 0. ) binhghx = h->GetNbinsX(); //default is max bin
  if( yhgh == 0. ) binhghy = h->GetNbinsY(); //default is max bin
  double integral = 0.;

  for( int j=binlowy; j < binhghy; j++ ) {
	for( int i=binlowx; i < binhghx; i++ ) { 
	  integral += h->GetBinContent(i,j);
	}
  }
  return integral;
}


//prints only when bin differences are found
//make sure h1, h2 have same NbinsX,Y
void compareTH2F(TH2F* h1, TH2F* h2) {
  //for( int i=0; i < h1->GetNbinsY(); i++ ) {
  for( int i=h1->GetNbinsY(); i >= 0; i-- ) { //start at top and go down
  //for( int i=h1->GetNbinsY(); i >= h1->GetNbinsY()-10; i-- ) { //start at top and go down, first 10 only
	for( int j=0; j < h1->GetNbinsX(); j++ ) {
	  if( h1->GetBinContent(j,i) - h2->GetBinContent(j,i) > 1.0e-5 )
		cout << "bin " << j << "," << i << "   first " << h1->GetBinContent(j,i) << "   second " << h2->GetBinContent(j,i) << endl;
	}
  }
}

void printTH2F(TH2F* h) {

  //for( int i=0; i < h->GetNbinsY(); i++ ) {
  for( int i=0; i < 10; i++ ) { //10 staring at bottom
  //for( int i=h->GetNbinsY(); i >= 0; i-- ) { //start at top and go down
  //for( int i=h->GetNbinsY(); i >= h->GetNbinsY()-10; i-- ) { //just 10 rows
	//for( int j=0; j < h->GetNbinsX(); j++ ) {
	for( int j=0; j < 30; j++ ) { //first 30 in x
	  printf("%4.1f ", h->GetBinContent(j,i) );
	}
	cout << endl;
  }

}

template <class TH> void printHistStats(TH* h, double xlow, double xhgh) {
  int hghbin = 0;
  if( xhgh == -1. ) hghbin = h->GetXaxis()->GetNbins() + 1; //include overflow
  else hghbin = h->GetXaxis()->FindFixBin( xhgh );
  cout << h->GetName() << " \t\t " << h->GetMean() << "  " << h->Integral()
	   << "  " << h->Integral( h->GetXaxis()->FindFixBin( xlow ), hghbin )/h->Integral() << endl;
}

//when shift is positive, the new one is old shifted up by shift
void HShift( TH1D*& h, int shift ) {
  TH1D* hnew = new TH1D( *h );
  int nbins = h->GetNbinsX();
  for( int i=0; i<=nbins; i++ ) {
	if( i-shift >= 0 ) //don't fill unless at underflow bin or higher
	  hnew->SetBinContent( i, h->GetBinContent(i-shift) );
	else
	  hnew->SetBinContent( i, 0 ); //zero these bins if they were filled at copy constructor call
  }
  //put content above bin (nbins - shift) into overflow
  hnew->SetBinContent( nbins+1, h->Integral(nbins-shift+1, nbins+1) );
  h = hnew;
}

TH2F* getTH2F( HistogramUtilities* h, sources_t theSources, TString var, TString var2, TString hyptyp, Int_t rebin, TString suffix, TString opt) {
  TCanvas* c = new TCanvas();
  TH2F* hnew = new TH2F(*(h->get2dHistogram(theSources, var, var2, hyptyp, rebin, suffix)));
  hnew->Draw(opt);
  c->SaveAs((TString)hnew->GetName()+".png");
  return hnew;
}

void saveStack(THStack* st, TLegend* leg, bool cpylog, double ymin, double ymax, TString name="", THStack* st2=0) {

  //double oldymax = st->GetMaximum();
  //double oldymin = st->GetMinimum();
  TCanvas *c = new TCanvas();
  TString usename;
  if( name != "" ) {
	usename = name;//+".png";
	//st->SetName(name);
	st->SetTitle(name);
  }
  else
	usename = (TString)st->GetName();//+".png";

  //for the y axis won't give me a max...just gives 1--don't use GetYaxis for stacks, use GetMaximum
  //cout << "ymin " << ymin << "  ymax " << ymax << "  st y max " << st->GetYaxis()->GetXmax() << "  max " << st->GetMaximum() << endl;
  if( ymax != -999 ) {
	//st->SetMinimum( ymin );
	st->SetMaximum( ymax );
  }
  //else
  //st->SetMinimum( ymin );
  //st->GetYaxis()->SetRangeUser( ymin, st->GetMaximum()+0.01*st->GetMaximum() );

  st->Draw("hist"); //must draw before the axis range is defined...(uh, don't ask me)
  if( st2 != 0 ) {
	st->SetMinimum( st2->GetMinimum() );
	st2->Draw("same,hist");
  }

  //if( verbose )
  //cout << "stack " << usename << "  max " << st->GetMaximum() << "  " << oldymax << "  min " << st->GetMinimum() << endl;
  //c->Update();
  //getc(stdin);
  leg->Draw();
  c->SaveAs(usename+".png");
  c->SaveAs(rootfilename);

  if(cpylog) {
	gPad->SetLogy();

	//st->SetMaximum( oldymax );
	//if( st2 == 0 ) { //don't reset min if st2 exists--this is wrong
	st->SetMinimum( ymin ); //need to reset regardless
	//cout << "set min to ymin " << ymin << endl;
	//}
  
	st->Draw("hist");
 	leg->Draw();
	c->Update();
	//c->SaveAs((TString)st->GetName()+"_log.png");
	c->SaveAs(usename+"_log.png");
	c->SaveAs(rootfilename);
  }
  
}

void makeStack(HistogramUtilities* h, TLegend* leg, sources_t theSources, TString title, TString subtitle, TString suffix, bool cpylog, double ymin, double ymax, TString name) {
  THStack *st = h->getStack(theSources, title, subtitle, suffix);
  saveStack( st, leg, cpylog, ymin, ymax, name );
}


  //Tests for ABCD
  //get this manually (not histogramUtilities way)
  /*
  TFile* file = new TFile("Results.root" , "READ");
  //TH2Fs
  //TH2F* h = (TH2F*)(file->Get(sources_[i].getName() + histNameSuffix)->Clone());
  TH2F* h = (TH2F*)(file->Get("wejetsAlpgen_lep_tcMet_relIso_all")->Clone());
  //printTH2F( h );
  cout << endl << endl;
  compareTH2F( h, h );

  TH2F* hdytt_a = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_all")->Clone());  
  TH2F* hdytt_m = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_m")->Clone());  
  //cout << endl << endl;
  //compareTH2F( hdytt_a, hdytt_m );

  TH2F* hdytt_e = (TH2F*)(file->Get("dyttAlpgen_lep_tcMet_relIso_e")->Clone());  
  //cout << endl << endl;
  hdytt_m->Add(hdytt_e);
  compareTH2F( hdytt_a, hdytt_m );
  */
  //cout << integrateTH2F( hdytt_a, 0, 50, 0.89, 1.0 ) << endl;
  //I expect this to be sum of:
  /*
bin 55,100  first 0.0754978   second 0
bin 60,98   first 0.00046883   second 0
bin 50,96   first 0.0008938   second 0
bin 21,95   first 0.0013502   second 0
bin 64,95   first 0.00787893   second 0
bin 66,94   first 0.00046883   second 0
bin 14,92   first 0.0000376945   3.76945e-05   second 0
bin 77,92   first 0.000284761   second 0
bin 81,92   first 0.00787893   second 0
bin 60,90   first 0.000927268   second 0

0.095687043
  */

  //same, but shifted--no, a shift is wrong, a scale (multiply x quantities by mw/mz) is needed in looper
  //cout << "int " << proj_wmet->Integral(0,100) << " undflw " << proj_wmet->GetBinContent(0) << " ovflw " << proj_wmet->GetBinContent(81) << endl;
  //HShift(proj_wmet, 5);
  //HShift(proj_wmetnu, 5);
  //cout << "int " << proj_wmet->Integral(0,100) << " undflw " << proj_wmet->GetBinContent(0) << " ovflw " << proj_wmet->GetBinContent(81) << endl;
  //cout << endl;
  //proj_zlpt->Draw();
  //proj_zmet->Draw("same");
  //proj_wmet->Draw("same");
  //proj_wmetnu->Draw("same");
  //c1->SaveAs("Compare_shift_zlpt_zmet_wmet.png");


/*
void aviCD( TH2F* single, TH2F* di, double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
  double a = integrateTH2F(single, x1, x2, y3, y4);
  double b = integrateTH2F(single, x3, x4, y3, y4);
  double c = integrateTH2F(single, x1, x2, y1, y2);
  double d = integrateTH2F(single, x3, x4, y1, y2);
  double A = integrateTH2F(di, x1, x2, y3, y4);
  double B = integrateTH2F(di, x3, x4, y3, y4);
  double C = integrateTH2F(di, x1, x2, y1, y2);
  double D = integrateTH2F(di, x3, x4, y1, y2);

  double bpr = b - a*B/A;
  double cpr = c - a*C/A;
  double dpr = d - a*D/A;
  cout << "prediction for a:  " << bpr*cpr/dpr << endl;
}
*/



