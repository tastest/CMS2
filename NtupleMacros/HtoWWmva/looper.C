#include "looper.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include <sstream>

#include "CORE/CMS2.h"
#include "CORE/trackSelections.h"
#include "CORE/metSelections.h"
#include "CORE/eventSelections.h"
#include "CORE/electronSelectionsParameters.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/mcSelections.h"
#include "Tools/goodrun.cc"
#include "CORE/utilities.cc"
#include "histtools.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"

char* iter          = "testFinal";
bool makebaby       = true;
bool makehist       = false;
bool maketext       = false;
bool debug          = false;
bool looseIdWJets   = false;

//--------------------------------------------------------------------
struct DorkyEventIdentifier {
  // this is a workaround for not having unique event id's in MC
  unsigned long int run, event,lumi;
  bool operator < (const DorkyEventIdentifier &) const;
  bool operator == (const DorkyEventIdentifier &) const;
};

//--------------------------------------------------------------------
bool DorkyEventIdentifier::operator < (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return run < other.run;
  if (event != other.event)
    return event < other.event;
  if(lumi != other.lumi)
    return lumi < other.lumi;
  return false;
}

//--------------------------------------------------------------------
bool DorkyEventIdentifier::operator == (const DorkyEventIdentifier &other) const
{
  if (run != other.run)
    return false;
  if (event != other.event)
    return false;
  return true;
}

//--------------------------------------------------------------------
std::set<DorkyEventIdentifier> already_seen;
bool is_duplicate (const DorkyEventIdentifier &id) {
  std::pair<std::set<DorkyEventIdentifier>::const_iterator, bool> ret =
    already_seen.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------


using namespace tas;
void looper::ScanChain (TChain* chain, const char* prefix, bool isData, int nEventsToProcess){

  TString prefixstr(prefix);

  if (makehist) bookHistos();

  std::vector<std::string> jetcorr_filenames_pf;
  jetcorr_filenames_pf.clear();
  jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L2Relative.txt");
  jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L3Absolute.txt");
  if (isData) 
    jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L2L3Residual.txt");
  jet_corrector_pf= makeJetCorrector(jetcorr_filenames_pf);    

  if (isData) {
    set_goodrun_file("goodruns.txt");
  }
  if (maketext) ofile.open( Form( "output/%s_%s_events.txt" , prefix , iter) );

  TH2F* frHisto = 0;
  TFile* rfFile = 0;
  if (prefixstr.Contains("WJetsToLNu") && looseIdWJets) {
    rfFile = TFile::Open("ww_el_fr_EGandEGMon.root");
    frHisto = (TH2F*) rfFile->Get("el_fr_v4_wwV1");
  }

  TObjArray *listOfFiles = chain->GetListOfFiles();

  unsigned int nEventsChain = 0;
  if(nEventsToProcess == -1) nEventsToProcess = chain->GetEntries();
  nEventsChain = nEventsToProcess;
  unsigned int nEventsTotal = 0;

  MakeBabyNtuple( Form( "output/%s_%s_baby.root" , prefix , iter) );

  if( debug ) cout << "Begin looping over files" << endl;

  // file loop
  TIter fileIter(listOfFiles);
  TFile* currentFile = 0;
  while ((currentFile = (TFile*)fileIter.Next())) {

    TFile f(currentFile->GetTitle());
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);

    // event loop
    unsigned int nEvents = tree->GetEntries();
    for (unsigned int event = 0; event < nEvents; ++event) {
      if( debug ) cout << "Event " << event << endl;
	  
      cms2.GetEntry(event);
      ++nEventsTotal;

      // progress feedback to user
      if (nEventsTotal % 1000 == 0) {
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal/(nEventsChain*0.01));
	  fflush(stdout);
	}
      }

      // skip duplicates 
      if( isData ) {
	DorkyEventIdentifier id = { evt_run(),evt_event(), evt_lumiBlock() };
	if (is_duplicate(id) ){
	  continue;
	}
      }
         
      if (prefixstr.Contains("VVJetsTo4L")) {
	if (!isWW()) continue;
      }

      //APPLY BASIC EVENT SELECTIONS 
      //TRIGGER, TRACKING AND VERTEX      
      if ( isData && !goodrun(cms2.evt_run(), cms2.evt_lumiBlock()))  continue;
      if( !cleaning_standardAugust2010( isData) )                     continue;

      if (nGoodVertex()<1) continue; 

      //trigger!!!!!!!!!!!!!     

      // N.B. BABY NTUPLE IS FILLED
      // FOR EACH EVENT
      InitBabyNtuple();

      if (prefixstr.Contains("VVJetsTo4L"))          sample_id_= 0;
      if (prefixstr.Contains("WWTo2L2Nu"))           sample_id_= 0;
      if (prefixstr.Contains("GluGluToWWTo4L"))      sample_id_= 1;
      if (prefixstr.Contains("TTJets"))              sample_id_= 2;
      if (prefixstr.Contains("WJetsToLNu"))          sample_id_= 3;
      if (prefixstr.Contains("DYToEEM10To20"))       sample_id_= 4;
      if (prefixstr.Contains("DYToEEM20"))           sample_id_= 4;
      if (prefixstr.Contains("DYToMuMuM10To20"))     sample_id_= 4;
      if (prefixstr.Contains("DYToMuMuM20"))         sample_id_= 4;
      if (prefixstr.Contains("DYToTauTauM10To20"))   sample_id_= 4;
      if (prefixstr.Contains("DYToTauTauM20"))       sample_id_= 4;
      if (prefixstr.Contains("WZ"))                  sample_id_= 5;
      if (prefixstr.Contains("ZZ"))                  sample_id_= 6;
      if (prefixstr.Contains("tW"))                  sample_id_= 7;
      if (prefixstr.Contains("HToWWTo2L2NuM130"))    sample_id_= 13;
      if (prefixstr.Contains("HToWWTo2Tau2NuM130"))  sample_id_= 13;
      if (prefixstr.Contains("HToWWToLNuTauNuM130")) sample_id_= 13;
      if (prefixstr.Contains("HToWWTo2L2NuM160"))    sample_id_= 16;
      if (prefixstr.Contains("HToWWTo2Tau2NuM160"))  sample_id_= 16;
      if (prefixstr.Contains("HToWWToLNuTauNuM160")) sample_id_= 16;
      if (prefixstr.Contains("HToWWTo2L2NuM200"))    sample_id_= 20;
      if (prefixstr.Contains("HToWWTo2Tau2NuM200"))  sample_id_= 20;
      if (prefixstr.Contains("HToWWToLNuTauNuM200")) sample_id_= 20;

      for(unsigned int i = 0; i < hyp_p4().size(); ++i) {

	//cout << "new event" << endl;
	int countWproducts = 0;
	if (!isData && i==0) {
	  for (unsigned int ig=0;ig<genps_id().size();++ig){
	    //cout << genps_id_mother().at(ig) << " " << genps_id().at(ig) << endl;
	    if (abs(genps_id_mother().at(ig))!=24) continue;
	    countWproducts++;
	    int id=genps_id().at(ig);
	    if (id==11||id==13||id==15) {
	      gen_lm_pt_= genps_p4().at(ig).pt();
	      gen_lm_eta_= genps_p4().at(ig).eta();
	      gen_lm_phi_= genps_p4().at(ig).phi();
	      gen_lm_id_= id;
	    }
	    else if (id==-11||id==-13||id==-15) {
	      gen_lp_pt_= genps_p4().at(ig).pt();
	      gen_lp_eta_= genps_p4().at(ig).eta();
	      gen_lp_phi_= genps_p4().at(ig).phi();
	      gen_lp_id_= id;
	    }
	    else if (id==12||id==14||id==16) {
	      gen_nm_pt_= genps_p4().at(ig).pt();
	      gen_nm_eta_= genps_p4().at(ig).eta();
	      gen_nm_phi_= genps_p4().at(ig).phi();
	      gen_nm_id_= id;
	    }
	    else if (id==-12||id==-14||id==-16) {
	      gen_np_pt_= genps_p4().at(ig).pt();
	      gen_np_eta_= genps_p4().at(ig).eta();
	      gen_np_phi_= genps_p4().at(ig).phi();
	      gen_np_id_= id;
	    }
	  }
	}        
	//if (i==0&&countWproducts!=4) {cout << "evt=" << evt_event() << " run=" << evt_run() << " nprod=" << countWproducts << endl;}

	//make selection for event baby
	event_type_ = hyp_type()[i];//0=mu+mu;1=el+mu;2=mu+el;3=el+el (where 2nd is lowerPt)

        //check that hyp leptons come from same vertex  
	//this does nothing on mc... check on data (present in ww official?)
        if( !hypsFromSameVtx( i ) ) continue;

        //OS, pt > (20,10) GeV 
        if( hyp_lt_id()[i] * hyp_ll_id()[i] > 0 )  continue;
	if (max(hyp_ll_p4()[i].pt(),hyp_lt_p4()[i].pt())<20) continue;
	if (min(hyp_ll_p4()[i].pt(),hyp_lt_p4()[i].pt())<10) continue;

	//store the info of who is the trailing lepton
	bool minLL = hyp_ll_p4()[i].pt() < hyp_lt_p4()[i].pt();

	//leptons have to pass WWv1 ID cuts (special treatment for trainling electrons in WJetsToLNu sample)
        if (abs(hyp_ll_id()[i])==13 && (! muonId(hyp_ll_index()[i] , NominalWWV1 ) ) )   continue;
        if (abs(hyp_lt_id()[i])==13 && (! muonId(hyp_lt_index()[i] , NominalWWV1 ) ) )   continue;          
	if (!(prefixstr.Contains("WJetsToLNu") && looseIdWJets)) {//standard
	  if (abs(hyp_ll_id()[i])==11 && (!pass_electronSelection(hyp_ll_index()[i], electronSelection_wwV1, false, false))) continue;
	  if (abs(hyp_lt_id()[i])==11 && (!pass_electronSelection(hyp_lt_index()[i], electronSelection_wwV1, false, false))) continue;
	} else { //WJetsToLNu, looseIdWJets==1 => apply FO selection on trailing electron
	  if (minLL) {
	    if (abs(hyp_ll_id()[i])==11 && (!pass_electronSelection(hyp_ll_index()[i], electronSelectionFO_el_wwV1_v4, false, false))) continue;
	    //if (abs(hyp_ll_id()[i])==11 && (!pass_electronSelection(hyp_ll_index()[i], electronSelection_wwV1_WP95, false, false))) continue;
	    if (abs(hyp_lt_id()[i])==11 && (!pass_electronSelection(hyp_lt_index()[i], electronSelection_wwV1, false, false))) continue;
	  } else {
	    if (abs(hyp_lt_id()[i])==11 && (!pass_electronSelection(hyp_lt_index()[i], electronSelectionFO_el_wwV1_v4, false, false))) continue;
	    //if (abs(hyp_lt_id()[i])==11 && (!pass_electronSelection(hyp_lt_index()[i], electronSelection_wwV1_WP95, false, false))) continue;
	    if (abs(hyp_ll_id()[i])==11 && (!pass_electronSelection(hyp_ll_index()[i], electronSelection_wwV1, false, false))) continue;
	  }
	}

	//reject events from low mass resonances (ee and mm only) 
	if( cms2.hyp_p4().at(i).mass2()<0 )  continue;
	if( (event_type_<0.5 || event_type_>2.5) && hyp_p4()[i].mass() < 12 )  continue;

	//reject events from Z peak (ee and mm only)  
	if( (event_type_<0.5 || event_type_>2.5) && fabs(hyp_p4()[i].mass()-91.1876)<15. )  continue;

	//event variables
        event_evt_= evt_event();
        event_run_= evt_run();
        event_lumi_= evt_lumiBlock();
	event_xsec_incl_= evt_xsec_incl();
	event_xsec_excl_= evt_xsec_excl();
	event_kfactor_= evt_kfactor();
	event_scale1fb_= evt_scale1fb();
	if (prefixstr.Contains("HToWWTo2L2Nu")) event_scale1fb_=event_scale1fb_*4./9.;//the cross section accounts for taus which are not in the sample
	else if (prefixstr.Contains("HToWWToLNuTauNu")) event_scale1fb_=event_scale1fb_*4./9.;//the cross section accounts for taus which are not in the sample
	else if (prefixstr.Contains("HToWWTo2Tau2Nu")) event_scale1fb_=event_scale1fb_*1./9.;//the cross section accounts for taus which are not in the sample
	else if (prefixstr.Contains("VVJetsTo4L")) event_scale1fb_ = 0.0060669;//magic number stolen from ww official code (processData.C), to be checked

	//dilepton variables
	dil_mass_= hyp_p4().at(i).mass();
	dil_pt_= hyp_p4().at(i).pt();
	dil_eta_= hyp_p4().at(i).eta();
	dil_phi_= hyp_p4().at(i).phi();
	dil_dphi_= looper::deltaPhi(hyp_ll_p4()[i].phi(),hyp_lt_p4()[i].phi());
	dil_metdphi_= looper::deltaPhi(hyp_p4()[i].phi(),evt_tcmetPhi());
	dil_deta_ = fabs(hyp_ll_p4()[i].eta()-hyp_lt_p4()[i].eta());
	dil_dr_ = sqrt(dil_dphi_*dil_dphi_+dil_deta_*dil_deta_);

	//hard lepton variables
        lephard_q_= minLL ? hyp_lt_charge()[i] : hyp_ll_charge()[i];
        lephard_id_= minLL ? hyp_lt_id()[i] : hyp_ll_id()[i];
        lephard_pt_= minLL ? hyp_lt_p4()[i].pt() : hyp_ll_p4()[i].pt();
        lephard_eta_= minLL ? hyp_lt_p4()[i].eta() : hyp_ll_p4()[i].eta();
        lephard_phi_= minLL ? hyp_lt_p4()[i].phi() : hyp_ll_p4()[i].phi();
	unsigned int lephard_index = minLL ? hyp_lt_index()[i] : hyp_ll_index()[i];
	if (abs(lephard_id_)==11) {
	  lephard_hOverE_= els_hOverE().at(lephard_index);
	  lephard_dEtaIn_= els_dEtaIn().at(lephard_index);
	  lephard_dPhiIn_= els_dPhiIn().at(lephard_index);
	  lephard_sigmaIEtaIEta_= els_sigmaIEtaIEta().at(lephard_index);
	  lephard_e2x5Max_= els_e2x5Max().at(lephard_index);
	  lephard_e1x5_= els_e1x5().at(lephard_index);
	  lephard_e5x5_= els_e5x5().at(lephard_index);
	  lephard_eSC_= els_eSC().at(lephard_index);
	  lephard_etaSC_= els_etaSC().at(lephard_index);
	  lephard_eOverPIn_= els_eOverPIn().at(lephard_index);
	  lephard_eOverPOut_= els_eOverPOut().at(lephard_index);
	  lephard_fbrem_= els_fbrem().at(lephard_index);
	  lephard_genId_= els_mc_id().at(lephard_index);
	  lephard_genMotherId_= els_mc_motherid().at(lephard_index);
	  if (event_type_>1.5 && event_type_<2.5) event_type_=1;//0=mu+mu;1=el+mu;2=mu+el;3=el+el (where 2nd is lowerPt)
	}

	//soft lepton variables
	lepsoft_q_= minLL ? hyp_ll_charge()[i] : hyp_lt_charge()[i];
        lepsoft_id_= minLL ? hyp_ll_id()[i] : hyp_lt_id()[i];
        lepsoft_pt_= minLL ? hyp_ll_p4()[i].pt() : hyp_lt_p4()[i].pt();
        lepsoft_eta_= minLL ? hyp_ll_p4()[i].eta() : hyp_lt_p4()[i].eta();
        lepsoft_phi_= minLL ? hyp_ll_p4()[i].phi() : hyp_lt_p4()[i].phi();
	unsigned int lepsoft_index = minLL ? hyp_ll_index()[i] : hyp_lt_index()[i];
	if (abs(lepsoft_id_)==11) {
	  lepsoft_hOverE_= els_hOverE().at(lepsoft_index);
	  lepsoft_dEtaIn_= els_dEtaIn().at(lepsoft_index);
	  lepsoft_dPhiIn_= els_dPhiIn().at(lepsoft_index);
	  lepsoft_sigmaIEtaIEta_= els_sigmaIEtaIEta().at(lepsoft_index);
	  lepsoft_e2x5Max_= els_e2x5Max().at(lepsoft_index);
	  lepsoft_e1x5_= els_e1x5().at(lepsoft_index);
	  lepsoft_e5x5_= els_e5x5().at(lepsoft_index);
	  lepsoft_eSC_= els_eSC().at(lepsoft_index);
	  lepsoft_etaSC_= els_etaSC().at(lepsoft_index);
	  lepsoft_eOverPIn_= els_eOverPIn().at(lepsoft_index);
	  lepsoft_eOverPOut_= els_eOverPOut().at(lepsoft_index);
	  lepsoft_fbrem_= els_fbrem().at(lepsoft_index);
	  lepsoft_genId_= els_mc_id().at(lepsoft_index);
	  lepsoft_genMotherId_= els_mc_motherid().at(lepsoft_index);
	  lepsoft_passTighterId_= 0;
	  if ( lepsoft_pt_>20 || (lepsoft_pt_<20 && (lepsoft_fbrem_>0.2 || (lepsoft_fbrem_<0.2&&fabs(lepsoft_etaSC_)<1.479 && lepsoft_eOverPIn_>0.95) ) ) 
	       ) lepsoft_passTighterId_= 1;
	  //set the fake rate in case of WJetsToLNu, looseIdWJets==1
	  if (prefixstr.Contains("WJetsToLNu") && looseIdWJets) lepsoft_fr_ = frHisto->GetBinContent(frHisto->GetXaxis()->FindBin(fabs(lepsoft_eta_)),
												     frHisto->GetYaxis()->FindBin(min(lepsoft_pt_,(static_cast<float>(34.9)))));
	  else lepsoft_fr_ = 1.;
	} else if (abs(lepsoft_id_)==13) {
	  lepsoft_passTighterId_= 1;
	  lepsoft_fr_ = 1.;
	}

	//met variables
	met_pt_= evt_tcmet();
	met_phi_= evt_tcmetPhi();
	met_sumet_=  evt_tcsumet();
	//projected met
	met_projpt_= evt_tcmet();
	double tightDPhi = fabs(cms2.hyp_lt_p4()[i].Phi() - met_phi_);
	tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
	double looseDPhi = fabs(cms2.hyp_ll_p4()[i].Phi() - met_phi_);
	looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
	double DeltaPhi =  TMath::Min(tightDPhi, looseDPhi);
	if (DeltaPhi < TMath::Pi()/2) met_projpt_= evt_tcmet()*TMath::Sin(DeltaPhi);
	//MT variables
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lepHardP4 = minLL ? hyp_lt_p4()[i] : hyp_ll_p4()[i];
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lepSoftP4 = minLL ? hyp_ll_p4()[i] : hyp_lt_p4()[i];
	ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metP4(met_pt_*cos(met_phi_),met_pt_*sin(met_phi_),0,met_pt_);
	mt_lephardmet_= (lepHardP4+metP4).mt();
	mt_lepsoftmet_= (lepSoftP4+metP4).mt();
	mt_dilmet_= (lepHardP4+lepSoftP4+metP4).mt();

	//jet and b-tag stuff
	vector<unsigned int> jetsForVeto = getJetIdVector(i, 25., 100000., 5.0);
	jets_num_=jetsForVeto.size();
	vector<unsigned int> jetsForBtags = getJetIdVector(i, 25., 100000., 2.5, true, 2.1);
	btags_num_=jetsForBtags.size();
	vector<unsigned int> jetsForLowPtBtags = getJetIdVector(i, 0., 25., 2.5, true, 2.1);
	lowptbtags_num_=jetsForLowPtBtags.size();
	vector<unsigned int> allJets = getJetIdVector(i, 0., 100000., 5.0);
	if (allJets.size()>0) {
	  unsigned int hardj_idx = getHardestJetId(allJets);
	  jethard_pt_ = getCorrectedJetPt(hardj_idx);
	  jethard_eta_ = pfjets_p4()[hardj_idx].eta();
	  jethard_phi_ = pfjets_p4()[hardj_idx].phi();
	  jethard_disc_ = pfjets_trackCountingHighEffBJetTag()[hardj_idx];
	}
	if (allJets.size()>1) {
	  unsigned int hardj2_idx = get2ndHardJetId(allJets);
	  jethard2_pt_ = getCorrectedJetPt(hardj2_idx);
	  jethard2_eta_ = pfjets_p4()[hardj2_idx].eta();
	  jethard2_phi_ = pfjets_p4()[hardj2_idx].phi();
	  jethard2_disc_ = pfjets_trackCountingHighEffBJetTag()[hardj2_idx];
	}
	vector<unsigned int> centralJets = getJetIdVector(i, 0., 100000., 2.5);
	if (centralJets.size()>0) {
	  unsigned int hardctrj_idx = getHardestJetId(centralJets);
	  ctrjethard_pt_ = getCorrectedJetPt(hardctrj_idx);
	  ctrjethard_eta_ = pfjets_p4()[hardctrj_idx].eta();
	  ctrjethard_phi_ = pfjets_p4()[hardctrj_idx].phi();
	  ctrjethard_disc_ = pfjets_trackCountingHighEffBJetTag()[hardctrj_idx];
	}
	float tmp_jets_discmax = -999;
	for (unsigned int j=0;j<centralJets.size();++j) {
	  unsigned int ijet = centralJets[j];
	  float disc = pfjets_trackCountingHighEffBJetTag()[ijet];
	  if (disc > tmp_jets_discmax) {
	    tmp_jets_discmax = disc;
	  }
	}
	jets_discmax_=tmp_jets_discmax;

	//veto on extra leptons and soft muons
	extralep_num_=getNExtraLeptons(i);
	softmu_num_=getNSoftMu(i,jetsForVeto);

	//sumpt's of jets, met and leptons
	llm_sumpt_ = lephard_pt_+lepsoft_pt_+met_pt_;
	llmj_sumpt_ = llm_sumpt_;
	if (jets_num_>0) llmj_sumpt_ += jethard_pt_;
	if (jets_num_>1) llmj_sumpt_ += jethard2_pt_;


	eventTree_->Fill();

	//make selection for histos
	if (makehist) {
	  // histos
	  fillUnderOverFlow( h_dPhi_ll, looper::deltaPhi(hyp_ll_p4()[i].phi(),hyp_lt_p4()[i].phi()), event_kfactor_*event_scale1fb_);
	  fillUnderOverFlow( h_dEta_ll, fabs(hyp_ll_p4()[i].eta()-hyp_lt_p4()[i].eta()), event_kfactor_*event_scale1fb_);
	  fillUnderOverFlow( h_mass_ll, hyp_p4()[i].mass(), event_kfactor_*event_scale1fb_);
	  fillUnderOverFlow( h_minPt_ll, min(hyp_ll_p4()[i].pt(),hyp_lt_p4()[i].pt()), event_kfactor_*event_scale1fb_);
	  fillUnderOverFlow( h_maxPt_ll, max(hyp_ll_p4()[i].pt(),hyp_lt_p4()[i].pt()), event_kfactor_*event_scale1fb_);
	}

      }
      
      //print info to output file
      if (maketext) ofile << "event=" << cms2.evt_event() << endl;

    } // end loop over events
  } // end loop over files
  
  if (nEventsChain != nEventsTotal)
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  
  CloseBabyNtuple();

  if (prefixstr.Contains("WJetsToLNu") && looseIdWJets) {
    rfFile->Close();
  }

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  if (makehist) {
    saveHist( Form( "output/%s_%s_histos.root" , prefix , iter ) );
    deleteHistos();
  }
  
} // end ScanChain

//--------------------------------------------------------------------

void looper::printEvent(  ostream& ostr ){
  ostr << cms2.evt_run() << " " << cms2.evt_lumiBlock() << " " << cms2.evt_event() << endl; 
}

//--------------------------------------------------------------------

void looper::InitBabyNtuple ()
{
  // event stuff
  sample_id_= -999999;
  event_evt_= -999999;
  event_run_= -999999;
  event_lumi_= -999999;
  event_type_= -999999;
  event_xsec_incl_= -999999;
  event_xsec_excl_= -999999;
  event_kfactor_= -999999;
  event_scale1fb_= -999999;
  dil_mass_= -999999;
  dil_pt_= -999999;
  dil_eta_= -999999;
  dil_phi_= -999999;
  dil_dphi_= -999999;
  lephard_q_= -999999;
  lephard_id_= -999999;
  lephard_pt_= -999999;
  lephard_eta_= -999999;
  lephard_phi_= -999999;
  lephard_hOverE_= -999999;
  lephard_dEtaIn_= -999999;
  lephard_dPhiIn_= -999999;
  lephard_sigmaIEtaIEta_= -999999;
  lephard_e2x5Max_= -999999;
  lephard_e1x5_= -999999;
  lephard_e5x5_= -999999;
  lephard_eSC_= -999999;
  lephard_etaSC_= -999999;
  lephard_eOverPIn_= -999999;
  lephard_eOverPOut_= -999999;
  lephard_fbrem_= -999999;
  lephard_genId_= -999999;
  lephard_genMotherId_= -999999;
  lepsoft_passTighterId_= -999999;
  lepsoft_q_= -999999;
  lepsoft_id_= -999999;
  lepsoft_pt_= -999999;
  lepsoft_eta_= -999999;
  lepsoft_phi_= -999999;
  lepsoft_hOverE_= -999999;
  lepsoft_dEtaIn_= -999999;
  lepsoft_dPhiIn_= -999999;
  lepsoft_sigmaIEtaIEta_= -999999;
  lepsoft_e2x5Max_= -999999;
  lepsoft_e1x5_= -999999;
  lepsoft_e5x5_= -999999;
  lepsoft_eSC_= -999999;
  lepsoft_etaSC_= -999999;
  lepsoft_eOverPIn_= -999999;
  lepsoft_eOverPOut_= -999999;
  lepsoft_fbrem_= -999999;
  lepsoft_genId_= -999999;
  lepsoft_genMotherId_= -999999;
  met_pt_= -999999;
  met_phi_= -999999;
  met_projpt_= -999999;
  jets_num_= -999999;
  lowptbtags_num_= -999999;
  btags_num_= -999999;
  extralep_num_= -999999;
  softmu_num_= -999999;
  jethard_pt_= -999999;
  jethard_eta_= -999999;
  jethard_phi_= -999999;
  jethard_disc_= -999999;
  jethard2_pt_= -999999;
  jethard2_eta_= -999999;
  jethard2_phi_= -999999;
  jethard2_disc_= -999999;
  ctrjethard_pt_= -999999;
  ctrjethard_eta_= -999999;
  ctrjethard_phi_= -999999;
  ctrjethard_disc_= -999999;
  
  jets_discmax_= -999999;

  met_sumet_= -999999;
  dil_metdphi_= -999999;
  dil_deta_= -999999;
  dil_dr_= -999999;

  mt_lephardmet_= -999999;
  mt_lepsoftmet_= -999999;
  mt_dilmet_= -999999;

  llm_sumpt_= -999999;
  llmj_sumpt_= -999999;
  
  gen_lp_pt_= -999999;
  gen_lp_eta_= -999999;
  gen_lp_phi_= -999999;
  gen_lp_id_= -999999;
  gen_lm_pt_= -999999;
  gen_lm_eta_= -999999;
  gen_lm_phi_= -999999;
  gen_lm_id_= -999999;
  gen_np_pt_= -999999;
  gen_np_eta_= -999999;
  gen_np_phi_= -999999;
  gen_np_id_= -999999;
  gen_nm_pt_= -999999;
  gen_nm_eta_= -999999;
  gen_nm_phi_= -999999;
  gen_nm_id_= -999999;

}

//--------------------------------------------------------------------

float looper::deltaPhi( float phi1 , float phi2 ){
  float dphi = fabs( phi1 - phi2 );
  if( dphi > TMath::Pi() ) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

//--------------------------------------------------------------------

void looper::bookHistos(){

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
  //float xmin = -100;

  h_dPhi_ll        = new TH1F("dPhi_ll","dPhi_ll",64,0,3.2);
  h_dEta_ll        = new TH1F("dEta_ll","dEta_ll",50,0,5);
  h_mass_ll        = new TH1F("mass_ll","mass_ll",150,0,300);
  h_minPt_ll        = new TH1F("minPt_ll","minPt_ll",100,0,200);
  h_maxPt_ll        = new TH1F("maxPt_ll","maxPt_ll",100,0,200);
}

//--------------------------------------------------------------------

void looper::MakeBabyNtuple (const char* babyFileName)
{
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();

  babyFile_ = new TFile(Form("%s", babyFileName), "RECREATE");

  babyFile_->cd();
  eventTree_ = new TTree("Events", "Events Tree");
  eventTree_->Branch("sample_id_"        , &sample_id_            , "sample_id/I");
  eventTree_->Branch("event_evt_"        , &event_evt_            , "event_evt/I");
  eventTree_->Branch("event_run_"        , &event_run_            , "event_run/I");
  eventTree_->Branch("event_lumi_"        , &event_lumi_            , "event_lumi/I");
  eventTree_->Branch("event_type_"        , &event_type_            , "event_type/F");
  eventTree_->Branch("event_xsec_incl_"        , &event_xsec_incl_            , "event_xsec_incl/F");
  eventTree_->Branch("event_xsec_excl_"        , &event_xsec_excl_            , "event_xsec_excl/F");
  eventTree_->Branch("event_kfactor_"        , &event_kfactor_            , "event_kfactor/F");
  eventTree_->Branch("event_scale1fb_"        , &event_scale1fb_            , "event_scale1fb/F");

  eventTree_->Branch("dil_mass_"        , &dil_mass_            , "dil_mass/F");
  eventTree_->Branch("dil_pt_"        , &dil_pt_            , "dil_pt/F");
  eventTree_->Branch("dil_eta_"        , &dil_eta_            , "dil_eta/F");
  eventTree_->Branch("dil_phi_"        , &dil_phi_            , "dil_phi/F");
  eventTree_->Branch("dil_dphi_"        , &dil_dphi_            , "dil_dphi/F");

  eventTree_->Branch("lephard_q_"        , &lephard_q_            , "lephard_q/I");
  eventTree_->Branch("lephard_id_"        , &lephard_id_            , "lephard_id/I");
  eventTree_->Branch("lephard_pt_"        , &lephard_pt_            , "lephard_pt/F");
  eventTree_->Branch("lephard_eta_"        , &lephard_eta_            , "lephard_eta/F");
  eventTree_->Branch("lephard_phi_"        , &lephard_phi_            , "lephard_phi/F");
  eventTree_->Branch("lephard_hOverE_"        , &lephard_hOverE_            , "lephard_hOverE/F");
  eventTree_->Branch("lephard_dEtaIn_"        , &lephard_dEtaIn_            , "lephard_dEtaIn/F");
  eventTree_->Branch("lephard_dPhiIn_"        , &lephard_dPhiIn_            , "lephard_dPhiIn/F");
  eventTree_->Branch("lephard_sigmaIEtaIEta_"        , &lephard_sigmaIEtaIEta_            , "lephard_sigmaIEtaIEta/F");
  eventTree_->Branch("lephard_e2x5Max_"        , &lephard_e2x5Max_            , "lephard_e2x5Max/F");
  eventTree_->Branch("lephard_e1x5_"        , &lephard_e1x5_            , "lephard_e1x5/F");
  eventTree_->Branch("lephard_e5x5_"        , &lephard_e5x5_            , "lephard_e5x5/F");
  eventTree_->Branch("lephard_eSC_"        , &lephard_eSC_            , "lephard_eSC/F");
  eventTree_->Branch("lephard_etaSC_"        , &lephard_etaSC_            , "lephard_etaSC/F");
  eventTree_->Branch("lephard_eOverPIn_"        , &lephard_eOverPIn_            , "lephard_eOverPIn/F");
  eventTree_->Branch("lephard_eOverPOut_"        , &lephard_eOverPOut_            , "lephard_eOverPOut/F");
  eventTree_->Branch("lephard_fbrem_"        , &lephard_fbrem_            , "lephard_fbrem/F");
  eventTree_->Branch("lephard_genId_"        , &lephard_genId_            , "lephard_genId/I");
  eventTree_->Branch("lephard_genMotherId_"        , &lephard_genMotherId_            , "lephard_genMotherId/I");
 
  eventTree_->Branch("lepsoft_q_"        , &lepsoft_q_            , "lepsoft_q/I");
  eventTree_->Branch("lepsoft_id_"        , &lepsoft_id_            , "lepsoft_id/I");
  eventTree_->Branch("lepsoft_pt_"        , &lepsoft_pt_            , "lepsoft_pt/F");
  eventTree_->Branch("lepsoft_eta_"        , &lepsoft_eta_            , "lepsoft_eta/F");
  eventTree_->Branch("lepsoft_phi_"        , &lepsoft_phi_            , "lepsoft_phi/F");
  eventTree_->Branch("lepsoft_hOverE_"        , &lepsoft_hOverE_            , "lepsoft_hOverE/F");
  eventTree_->Branch("lepsoft_dEtaIn_"        , &lepsoft_dEtaIn_            , "lepsoft_dEtaIn/F");
  eventTree_->Branch("lepsoft_dPhiIn_"        , &lepsoft_dPhiIn_            , "lepsoft_dPhiIn/F");
  eventTree_->Branch("lepsoft_sigmaIEtaIEta_"        , &lepsoft_sigmaIEtaIEta_            , "lepsoft_sigmaIEtaIEta/F");
  eventTree_->Branch("lepsoft_e2x5Max_"        , &lepsoft_e2x5Max_            , "lepsoft_e2x5Max/F");
  eventTree_->Branch("lepsoft_e1x5_"        , &lepsoft_e1x5_            , "lepsoft_e1x5/F");
  eventTree_->Branch("lepsoft_e5x5_"        , &lepsoft_e5x5_            , "lepsoft_e5x5/F");
  eventTree_->Branch("lepsoft_eSC_"        , &lepsoft_eSC_            , "lepsoft_eSC/F");
  eventTree_->Branch("lepsoft_etaSC_"        , &lepsoft_etaSC_            , "lepsoft_etaSC/F");
  eventTree_->Branch("lepsoft_eOverPIn_"        , &lepsoft_eOverPIn_            , "lepsoft_eOverPIn/F");
  eventTree_->Branch("lepsoft_eOverPOut_"        , &lepsoft_eOverPOut_            , "lepsoft_eOverPOut/F");
  eventTree_->Branch("lepsoft_fbrem_"        , &lepsoft_fbrem_            , "lepsoft_fbrem/F");
  eventTree_->Branch("lepsoft_genId_"        , &lepsoft_genId_            , "lepsoft_genId/I");
  eventTree_->Branch("lepsoft_genMotherId_"        , &lepsoft_genMotherId_            , "lepsoft_genMotherId/I");
  eventTree_->Branch("lepsoft_passTighterId_"        , &lepsoft_passTighterId_            , "lepsoft_passTighterId/I");
  eventTree_->Branch("lepsoft_fr_"        , &lepsoft_fr_            , "lepsoft_fr/F");

  eventTree_->Branch("met_pt_"        , &met_pt_            , "met_pt/F");
  eventTree_->Branch("met_phi_"        , &met_phi_            , "met_phi/F");
  eventTree_->Branch("met_projpt_"        , &met_projpt_            , "met_projpt/F");

  eventTree_->Branch("jets_num_"        , &jets_num_            , "jets_num/I");
  eventTree_->Branch("lowptbtags_num_"        , &lowptbtags_num_            , "lowptbtags_num/I");
  eventTree_->Branch("btags_num_"        , &btags_num_            , "btags_num/I");
  eventTree_->Branch("extralep_num_"        , &extralep_num_            , "extralep_num/I");
  eventTree_->Branch("softmu_num_"        , &softmu_num_            , "softmu_num/I");

  eventTree_->Branch("jethard_pt_"        , &jethard_pt_            , "jethard_pt/F");
  eventTree_->Branch("jethard_eta_"        , &jethard_eta_            , "jethard_eta/F");
  eventTree_->Branch("jethard_phi_"        , &jethard_phi_            , "jethard_phi/F");
  eventTree_->Branch("jethard_disc_"        , &jethard_disc_            , "jethard_disc/F");
  eventTree_->Branch("jethard2_pt_"        , &jethard2_pt_            , "jethard2_pt/F");
  eventTree_->Branch("jethard2_eta_"        , &jethard2_eta_            , "jethard2_eta/F");
  eventTree_->Branch("jethard2_phi_"        , &jethard2_phi_            , "jethard2_phi/F");
  eventTree_->Branch("jethard2_disc_"        , &jethard2_disc_            , "jethard2_disc/F");
  eventTree_->Branch("ctrjethard_pt_"        , &ctrjethard_pt_            , "ctrjethard_pt/F");
  eventTree_->Branch("ctrjethard_eta_"        , &ctrjethard_eta_            , "ctrjethard_eta/F");
  eventTree_->Branch("ctrjethard_phi_"        , &ctrjethard_phi_            , "ctrjethard_phi/F");
  eventTree_->Branch("ctrjethard_disc_"        , &ctrjethard_disc_            , "ctrjethard_disc/F");

  eventTree_->Branch("jets_discmax_"        , &jets_discmax_            , "jets_discmax/F");

  eventTree_->Branch("met_sumet_"        , &met_sumet_            , "met_sumet/F");
  eventTree_->Branch("dil_metdphi_"        , &dil_metdphi_            , "dil_metdphi/F");
  eventTree_->Branch("dil_deta_"        , &dil_deta_            , "dil_deta/F");
  eventTree_->Branch("dil_dr_"        , &dil_dr_            , "dil_dr/F");

  eventTree_->Branch("mt_lephardmet_"        , &mt_lephardmet_            , "mt_lephardmet/F");
  eventTree_->Branch("mt_lepsoftmet_"        , &mt_lepsoftmet_            , "mt_lepsoftmet/F");
  eventTree_->Branch("mt_dilmet_"        , &mt_dilmet_            , "mt_dilmet/F");

  eventTree_->Branch("llm_sumpt_"        , &llm_sumpt_            , "llm_sumpt/F");
  eventTree_->Branch("llmj_sumpt_"        , &llmj_sumpt_            , "llmj_sumpt/F");

  eventTree_->Branch("gen_lp_pt_"        , &gen_lp_pt_            , "gen_lp_pt/F");
  eventTree_->Branch("gen_lp_eta_"        , &gen_lp_eta_            , "gen_lp_eta/F");
  eventTree_->Branch("gen_lp_phi_"        , &gen_lp_phi_            , "gen_lp_phi/F");
  eventTree_->Branch("gen_lp_id_"        , &gen_lp_id_            , "gen_lp_id/I");
  eventTree_->Branch("gen_lm_pt_"        , &gen_lm_pt_            , "gen_lm_pt/F");
  eventTree_->Branch("gen_lm_eta_"        , &gen_lm_eta_            , "gen_lm_eta/F");
  eventTree_->Branch("gen_lm_phi_"        , &gen_lm_phi_            , "gen_lm_phi/F");
  eventTree_->Branch("gen_lm_id_"        , &gen_lm_id_            , "gen_lm_id/I");
  eventTree_->Branch("gen_np_pt_"        , &gen_np_pt_            , "gen_np_pt/F");
  eventTree_->Branch("gen_np_eta_"        , &gen_np_eta_            , "gen_np_eta/F");
  eventTree_->Branch("gen_np_phi_"        , &gen_np_phi_            , "gen_np_phi/F");
  eventTree_->Branch("gen_np_id_"        , &gen_np_id_            , "gen_np_id/I");
  eventTree_->Branch("gen_nm_pt_"        , &gen_nm_pt_            , "gen_nm_pt/F");
  eventTree_->Branch("gen_nm_eta_"        , &gen_nm_eta_            , "gen_nm_eta/F");
  eventTree_->Branch("gen_nm_phi_"        , &gen_nm_phi_            , "gen_nm_phi/F");
  eventTree_->Branch("gen_nm_id_"        , &gen_nm_id_            , "gen_nm_id/I");
}

//--------------------------------------------------------------------

void looper::fillUnderOverFlow(TH1F *h1, float value, float weight){

  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void looper::CloseBabyNtuple (){

  babyFile_->cd();
  eventTree_->Write();
  babyFile_->Close();
  
}

//--------------------------------------------------------------------

bool looper::isGoodTrack( int index ) {
  
     float corrected_d0 = trks_d0corr().at(index);

     if( trks_algo().at(index) < 8 ) {

       float d0cut = sqrt( pow(0.015,2) + pow(0.5/trks_trk_p4().at(index).pt(),2) );
       if( d0cut > 0.3 ) d0cut = 0.3;

       if( fabs( corrected_d0 ) > d0cut )                return false;
       //if( trks_nlayers().at(index) < nlayerscut_4567_ ) return false;
     }
     else {
       //if( trks_nlayers().at(index) < nlayerscut_89_ )                     return false;
       if( trks_validHits().at(index) < 9 )                                return false;
       if( trks_chi2().at(index) / trks_ndof().at(index) > 5. )            return false;
       if( trks_ptErr().at(index) / trks_trk_p4().at(index).pt() > 0.20 )  return false;
     }

     if( trks_validHits().at(index) < 6 )                               return false;
     if( trks_chi2().at(index) / trks_ndof().at(index) > 5 )            return false;
     if( fabs( trks_trk_p4().at(index).eta() ) > 2.65 )                 return false;
     if( trks_trk_p4().at(index).pt() > 100 )                           return false;
     if( trks_ptErr().at(index) / trks_trk_p4().at(index).pt() > 0.20 ) return false;
     if( !isTrackQuality( index, (1 << highPurity) ) )                  return false;

     if( trks_trk_p4().at(index).pt() > 0 && fabs(trks_trk_p4().at(index).eta()) > 2.5 ) return false;
 
     return true;
}

//--------------------------------------------------------------------

bool looper::isMuon( int index ) {

     for( unsigned int i = 0; i < mus_p4().size(); i++ ) {

	  if( mus_trkidx().at(i) == index ) return true;
     }

     return false;
}

//--------------------------------------------------------------------

bool looper::isElectron( int index ) {

     for( unsigned int i = 0; i < els_p4().size(); i++ ) {

	  if( els_trkidx().at(i) == index && els_hOverE().at(i) < 0.1 ) return true;
     }

     return false;
}

bool looper::isGoodVertex(size_t ivtx) {
    if (cms2.vtxs_isFake()[ivtx]) return false;
    if (cms2.vtxs_ndof()[ivtx] < 4.) return false;
    if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
    if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
    return true;
}
unsigned int looper::nGoodVertex() {
  unsigned int nVtx = 0;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    // if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    nVtx++;
  }
  return nVtx;
}
//--------------------------------------------------------------------

vector<unsigned int> looper::getJetIdVector(unsigned int i_hyp, float etMin, float etMax, float etaMax, bool doBtag, float discCut, float vetoCone) {

  vector<unsigned int> result;
  for ( unsigned int ij=0; ij < cms2.pfjets_p4().size(); ++ij) {        
    if ( TMath::Abs(cms2.pfjets_p4()[ij].eta()) > etaMax ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[ij])) < vetoCone ||
	 TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[ij])) < vetoCone ) continue;
        float pt = getCorrectedJetPt(ij);
    if ( pt < etMin ) continue;
    if ( pt > etMax ) continue;
    if (!doBtag) result.push_back(ij);
    else if (pfjets_trackCountingHighEffBJetTag()[ij]>discCut) result.push_back(ij);
  }
  return result;

}
//--------------------------------------------------------------------

unsigned int looper::getHardestJetId(vector<unsigned int> jetIdVec) {

  assert(jetIdVec.size()>0);

  unsigned int result = 0;
  float maxPt = 0;
  for (unsigned int j=0;j<jetIdVec.size();++j) {
    unsigned int ijet = jetIdVec[j];
    float pt = getCorrectedJetPt(ijet);
    if (pt > maxPt) {
      result=ijet;
      maxPt = pt;
    }
  }
  return result;
  
}
//--------------------------------------------------------------------

unsigned int looper::get2ndHardJetId(vector<unsigned int> jetIdVec) {

  assert(jetIdVec.size()>1);

  unsigned int result = 0;
  unsigned int maxId = 0;
  float maxPt = -1.;
  float secPt = -2.;
  for (unsigned int j=0;j<jetIdVec.size();++j) {
    unsigned int ijet = jetIdVec[j];
    float pt = getCorrectedJetPt(ijet);
    if (pt > maxPt) {
      result = maxId;
      maxId =ijet;
      secPt = maxPt;
      maxPt = pt;
    } else if (pt > secPt) {
      result=ijet;
      secPt = pt;      
    }
  }
  return result;

}
//--------------------------------------------------------------------
 
float looper::getCorrectedJetPt(unsigned int jetId){
  float jec = jetCorrection(cms2.pfjets_p4()[jetId], jet_corrector_pf);
    return cms2.pfjets_p4()[jetId].pt() * jec;
}
//--------------------------------------------------------------------

unsigned int looper::getNExtraLeptons(unsigned int i_hyp, double minLepPt) {

  unsigned int nMuons = 0;
  for (int ilep=0; ilep < int(cms2.mus_charge().size()); ++ilep) {
    // printf("Muon: %u, pt: %0.2f\n",i,cms2.mus_p4().at(ilep).pt());
    if ( cms2.mus_p4()[ilep].pt() < minLepPt ) continue;
    // printf("\tpassed minLepPt\n");
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == ilep ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == ilep ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( !muonId(ilep, NominalWWV1 ) ) continue;
    // printf("\tpassed all\n");
    ++nMuons;
  }
  unsigned int nElectrons = 0;
  for (int ilep=0; ilep < int(cms2.els_charge().size()); ++ilep) {
    // printf("Electron: %u, pt: %0.2f\n",i,cms2.els_p4().at(ilep).pt());
    if ( cms2.els_p4()[ilep].pt() < minLepPt ) continue;
    // printf("\tpassed minLepPt\n");
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.els_p4().at(ilep)) <0.1) ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.els_p4().at(ilep)) <0.1) ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( !pass_electronSelection(ilep, electronSelection_wwV1, false, false) ) continue;
    // printf("\tpassed all\n");
    ++nElectrons;
  }
  return nMuons+nElectrons;

}
//--------------------------------------------------------------------

unsigned int looper::getNSoftMu(unsigned int i_hyp, std::vector<unsigned int> jetIdVec) {

  unsigned int result = 0;
  bool nonisolated = true;
  for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
    // quality cuts
    if (  ((cms2.mus_goodmask()[imu]) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
    if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
    if ( TMath::Abs(cms2.mus_d0corr()[imu]) > 0.2) continue;
    //if ( TMath::Abs( fabs(cms2.mus_d0corr()[index])ww_mud0PV()) > 0.2) continue;//this is the official one... slightly different for pu
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    double sum = cms2.mus_iso03_sumPt().at(imu)+cms2.mus_iso03_emEt().at(imu)+cms2.mus_iso03_hadEt().at(imu);
    double pt  = cms2.mus_p4().at(imu).pt();
    float isoVal = sum/pt;
    if ( nonisolated && isoVal<0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
    bool skip = false;
    for (unsigned int j=0;j<jetIdVec.size();++j) {
      unsigned int ijet = jetIdVec[j];
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.pfjets_p4()[ijet],cms2.mus_p4()[imu])) < 0.3 ) skip=true;
    }
    if ( skip ) continue;
    ++result;
  }

  return result;

}
//--------------------------------------------------------------------
