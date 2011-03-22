/* Usage:
   root [0] .L ScanChain.C++
   root [1] TChain *chain = new TChain("Events")
   root [2] chain->Add("merged_ntuple.root")
   root [3] ScanChain(chain)
*/

// C++
#include <iostream>
#include <vector>

// ROOT
#include "TChain.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TAxis.h"

// CMS2
#include "branches.h"
#include "TBitSet.cc"
#include "CORE/CMS2.cc"
#include "CORE/electronSelections.cc"
#include "CORE/electronSelectionsParameters.cc"
#include "CORE/muonSelections.cc"
#include "CORE/jetSelections.cc"
#include "CORE/mcSelections.cc"
#include "MCUtil.h"

// ME related


using namespace tas;

bool applyJEC = true;

std::vector<std::string> jetcorr_filenames_jpt;
FactorizedJetCorrector *jet_corrector_jpt;

std::vector<std::string> jetcorr_filenames_pf;
FactorizedJetCorrector *jet_corrector_pf;

std::vector<std::string> jetcorr_filenames_calo;
FactorizedJetCorrector *jet_corrector_calo;

std::vector<std::string> jetcorr_filenames_trk;
FactorizedJetCorrector *jet_corrector_trk;


int getHypothesisType( int NTtype ){
  switch (NTtype){
  case ElEl:
    return ee;
    break;
  case MuMu:
    return mumu;
    break;
  case ElMu:
    return emu;
    break;
  case MuEl:
    return emu;
    break;
  }
  return emu; 
}


bool is_duplicate (const EventIdentifier &id){
     std::pair<std::set<EventIdentifier>::const_iterator, bool> ret = already_seen.insert(id);
     return !ret.second;
}



void progress( int nEventsTotal, int nEventsChain ){
  int period = 1000;
  if(nEventsTotal%1000 == 0) {
    // xterm magic from L. Vacavant and A. Cerri
    if (isatty(1)) {
      if( ( nEventsChain - nEventsTotal ) > period ){
        float frac = (float)nEventsTotal/(nEventsChain*0.01);
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
             "\033[0m\033[32m <---\033[0m\015", frac);
        fflush(stdout);
      }
      else {
        printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", 100.);
        cout << endl;
      }
    }
  }
}

void ScanChain(const char* process, TChain *chain, TFile *utilFile_,  int nEvents = -1, double IntLumi=100, double Xsect=1.0, int nProcessedEvents=-1, std::string skimFilePrefix="", bool realData=false, bool identifyEvents=false){

  _cutWord = new TBitSet(kNCuts);

  _cutMask = new TBitSet(kNCuts);
  _cutMask->SetAll();
  //_cutMask->SetFalse(kcut_zsel);

  already_seen.clear();

  for ( int count=0; count < kNdilepFlav; count++) {
  eventCount[count]=0;
  eventYield[count]=0;
  }

  // Apply the JEC
  if (applyJEC) {
    try { 
      jetcorr_filenames_jpt.clear();
      jetcorr_filenames_jpt.push_back("CORE/jetcorr/START38_V13_AK5JPT_L2Relative.txt");
      jetcorr_filenames_jpt.push_back("CORE/jetcorr/START38_V13_AK5JPT_L3Absolute.txt");
      if(realData)
	jetcorr_filenames_jpt.push_back("CORE/jetcorr/START38_V13_AK5JPT_L2L3Residual.txt");
      jet_corrector_jpt= makeJetCorrector(jetcorr_filenames_jpt);
      
      jetcorr_filenames_pf.clear();
      jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L2Relative.txt");
      jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L3Absolute.txt");
      if (realData) 
	jetcorr_filenames_pf.push_back("CORE/jetcorr/START38_V13_AK5PF_L2L3Residual.txt");
      jet_corrector_pf= makeJetCorrector(jetcorr_filenames_pf);
      
      jetcorr_filenames_calo.clear();
      jetcorr_filenames_calo.push_back("CORE/jetcorr/START38_V13_AK5Calo_L2Relative.txt");
      jetcorr_filenames_calo.push_back("CORE/jetcorr/START38_V13_AK5Calo_L3Absolute.txt");
      if(realData)
	jetcorr_filenames_calo.push_back("CORE/jetcorr/START38_V13_AK5Calo_L2L3Residual.txt");
      jet_corrector_calo= makeJetCorrector(jetcorr_filenames_calo);
      
      jetcorr_filenames_trk.clear();
      jetcorr_filenames_trk.push_back("CORE/jetcorr/START38_V13_AK5TRK_L2Relative.txt");
      jetcorr_filenames_trk.push_back("CORE/jetcorr/START38_V13_AK5TRK_L3Absolute.txt");
      jet_corrector_trk= makeJetCorrector(jetcorr_filenames_trk);
    } catch (...){
      cout << "\nFailed to setup correctors needed to get Jet Enetry Scale. Abort\n" << endl;
      assert(0);
    }
  }

  // File Loop
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  InitSkimmedTree(process);
  if(!realData) InitMCUtilHist(process, utilFile_);

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    TFile f( currentFile->GetTitle() );
    TTree *tree = (TTree*)f.Get("Events");
    cms2.Init(tree);
    
    // Event Loop
    unsigned int nEvents = tree->GetEntries();

    cout<<"Opened file "<<currentFile->GetTitle()<<" with events = "<<nEvents<<"\n";
    for( unsigned int event = 0; event < nEvents; ++event) {
    
      // Get Event Content
      cms2.GetEntry(event);
      
      // identifyEvents by MC truth
      if ( !realData && identifyEvents ){
	if(!isIdentified (process)) continue;
       }
      if (cms2.trks_d0().size() == 0) continue;  // needed to get rid of back Monte Carlo events in CMSSW_2_X analysis
      if (cms2.hyp_type().size() == 0) continue; // skip events without hypothesis
      EventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) continue;
      
      double mcweight;
      if ((process!="data")){
	mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	weight = IntLumi * mcweight * (Xsect>0?Xsect:cms2.evt_xsec_excl()*cms2.evt_kfactor()) /
	  (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
      } else weight =1;
      
      // cout << "weight = " << weight <<endl;
      
      // ==== select one particular event ====
      // if(cms2.evt_event()!=20791) continue;
      
      //      CalculateFakeRateProb();

      if(!realData) FillEffHist(TString(process), weight);

      unsigned int nHyps = cms2.hyp_type().size();

      for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) break;
	_cutWord->SetAllBitsFalse();
	ApplyEventSelection(i_hyp);
	
	int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);


	bool accept = true;
	for (int j = 0; j < kNCuts; j++) {
	  if (_cutMask->IsTrue(j)) accept &= _cutWord->IsTrue (j);
        }

	/*
	// Additional event cut for the HWW cut-based cross-check to Table~4 of EWK-10-009
	// ------------------------------------------
	// delta phi btw leptons
	// ------------------------------------------
	double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
	if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
	if(dphi > 60./180.*TMath::Pi()) accept = false;
	
	// ------------------------------------------
	// pt (30, 25) cuts
	// ------------------------------------------
	const double min1 = 30;
	const double min2 = 25;
	if ( !( cms2.hyp_lt_p4().at(i_hyp).pt() > min1 && cms2.hyp_ll_p4().at(i_hyp).pt() > min2 ) &&
	     !( cms2.hyp_lt_p4().at(i_hyp).pt() > min2 && cms2.hyp_ll_p4().at(i_hyp).pt() > min1 ) ) accept = false;
	*/

        if (accept){
	  if(!realData) FillKtHist(TString(process), weight);
	  FillNeededVariables(i_hyp);
	  if (getHypothesisType(cms2.hyp_type()[i_hyp]) == ee ) 
	    dilmass_ee_->Fill(cms2.hyp_p4().at(i_hyp).mass(), weight);
	  if (getHypothesisType(cms2.hyp_type()[i_hyp]) == emu ) 
	    dilmass_em_->Fill(cms2.hyp_p4().at(i_hyp).mass(), weight);
	  if (getHypothesisType(cms2.hyp_type()[i_hyp]) == mumu ) 
	    dilmass_mm_->Fill(cms2.hyp_p4().at(i_hyp).mass(), weight);
	  
	  eventCount[type]++;
	  eventYield[type]+=weight;
	  eventCount[all]++;
	  eventYield[all]+=weight;
	  outTree_->Fill();
	}
      }

      nEventsTotal++;
    
      // Progress
      progress( nEventsTotal, nEventsChain );

    }
  
    delete tree;
    f.Close();
  }
  if ( nEventsChain != nEventsTotal ) {
    cout<<" nEventsChain = "<<nEventsChain<<"    nEventsTotal="<<nEventsTotal<<"\n";
    std::cout << "ERROR: number of events from files is not equal to total number of events" << std::endl;
  }


  for ( int count=0; count < kNdilepFlav; count++) {
    cout<<" Total Count in "<<count<<" is equal "<<eventCount[count]<<"\n";
    cout<<" Total Yield in "<<count<<" is equal "<<eventYield[count]<<"\n";
  }

  Integral=eventYield[all];
  
  cout << "Total yield " << eventYield[all] << ": MM "<<  eventYield[mumu] 
       << "; EM "<< eventYield[emu]
       << "; EE "<< eventYield[ee] <<endl;
    

  outGlobalTree_->Fill();

  // std::cout << "tree is filled..." << std::endl;
  // std::cout << "listing utilFile_" << std::endl;
  // utilFile_->ls();

  outFile_->cd();
  outTree_->Write();
  outGlobalTree_->Write();
  outFile_->Close();
  
  // std::cout << "outFile_ is now closed" << std::endl;
  // std::cout << "now listing utilFile_ again" << std::endl;
  // utilFile_->ls();

  if(!realData) saveMCUtilOutput(utilFile_);
  
  //  samplehisto->Draw();
  
  cout<<"Total Events Before Selection "<<chain->GetEntries()<<"\n";
}




int ApplyEventSelection( unsigned int i_hyp){

  int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);

  // Trigger
  //if ( passedTriggerRequirements())
  _cutWord->SetTrue(kcut_Trigger);

  // OS
  if ( fast_sign (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] ) < 0)  _cutWord->SetTrue(kcut_OS);
  
  // Leading Lepton requirements:
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    if ( goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) _cutWord->SetTrue(kcut_LT);}
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    if ( goodElectronIsolated(cms2.hyp_lt_index()[i_hyp]) ) _cutWord->SetTrue(kcut_LT);}
  
  // Sub-Leading Lepton requirements:
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    if ( goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) _cutWord->SetTrue(kcut_LL);}
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    if ( goodElectronIsolated(cms2.hyp_ll_index()[i_hyp]) ) _cutWord->SetTrue(kcut_LL);}
  
  // Mll cut at mH = 160
  if ( cms2.hyp_p4()[i_hyp].mass() > 12 && cms2.hyp_p4()[i_hyp].mass() < 100 ) _cutWord->SetTrue(kcut_Mll); 
  
  // Z window veto
  if ( type == ee || type == mumu ) {
    if (!inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))   _cutWord->SetTrue(kcut_zsel);
  }
  if (type == emu) _cutWord->SetTrue(kcut_zsel);
  
  // MET cut
  double pMet = projectedMet(i_hyp);
  if  (( type == ee || type == mumu ) && ( pMet > 35 )) _cutWord->SetTrue(kcut_met);
  if  ((type == emu ) && ( pMet > 20 )) _cutWord->SetTrue(kcut_met);
  
  
  // Jet Veto
  int njets = getJets(pfJet, i_hyp, 25, 5.0, true, false).size(); 
  if ( njets == 0 ) _cutWord->SetTrue(kcut_jetveto);

  // ==Soft Muon Veto                                                                                                                                        
  if ( numberOfSoftMuons(i_hyp,true) == 0) _cutWord->SetTrue(kcut_softmuonveto);
  // == Extra Lepton Veto  
  
  if ( numberOfExtraLeptons(i_hyp,10) == 0) _cutWord->SetTrue(kcut_extraleptonveto);
  
  // Top tagging
  
  if ( !toptag(CaloJet,i_hyp,0) ) _cutWord->SetTrue(kcut_toptag);
  
  return 0;
}


int FillNeededVariables(unsigned int i_hyp){


  _evInfo.eventNumber = cms2.evt_event();
  _evInfo.runNumber = cms2.evt_run();
  _evInfo.lumiblockNumber = cms2.evt_lumiBlock();
  _evInfo.weight = weight;
  _evInfo.dilflavor = getHypothesisType(cms2.hyp_type()[i_hyp]);
  _evInfo.tcMet = cms2.evt_tcmet();
  _evInfo.tcMetSpec = projectedMet(i_hyp);
  _evInfo.tcMetX = cms2.evt_tcmet()* cos(cms2.evt_tcmetPhi());
  _evInfo.tcMetY = cms2.evt_tcmet()* sin(cms2.evt_tcmetPhi());
  _evInfo.Njets = 0;
  _evInfo.Qprod = fast_sign (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] );
  _evInfo.dimass = cms2.hyp_p4()[i_hyp].mass();
  _evInfo.dPhiLeptons =  fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
  if ( _evInfo.dPhiLeptons > TMath::Pi()) _evInfo.dPhiLeptons= 2*TMath::Pi() - _evInfo.dPhiLeptons;
  _evInfo.dEtaLeptons =  fabs(cms2.hyp_lt_p4()[i_hyp].eta() - cms2.hyp_ll_p4()[i_hyp].eta());
  _evInfo.dRLeptons =  sqrt(_evInfo.dPhiLeptons*_evInfo.dPhiLeptons+_evInfo.dEtaLeptons*_evInfo.dEtaLeptons);

  _evInfo.lep1_Type = TMath::Abs(cms2.hyp_lt_id()[i_hyp]);
  _evInfo.lep1_FakeType = 0;
  _evInfo.lep1_Px = cms2.hyp_lt_p4().at(i_hyp).px();
  _evInfo.lep1_Py = cms2.hyp_lt_p4().at(i_hyp).py();
  _evInfo.lep1_Pz = cms2.hyp_lt_p4().at(i_hyp).pz();
  _evInfo.lep1_E  = cms2.hyp_lt_p4().at(i_hyp).E();
  _evInfo.lep1_Charge = fast_sign (cms2.hyp_lt_id()[i_hyp] );

  _evInfo.lep2_Type = TMath::Abs(cms2.hyp_ll_id()[i_hyp]);
  _evInfo.lep2_FakeType = 0;
  _evInfo.lep2_Px = cms2.hyp_ll_p4().at(i_hyp).px();
  _evInfo.lep2_Py = cms2.hyp_ll_p4().at(i_hyp).py();
  _evInfo.lep2_Pz = cms2.hyp_ll_p4().at(i_hyp).pz();
  _evInfo.lep2_E  = cms2.hyp_ll_p4().at(i_hyp).E();
  _evInfo.lep2_Charge = fast_sign (cms2.hyp_ll_id()[i_hyp] );


  return 0;
}



//Utility fucntions


bool inZmassWindow(float mass){
  // return ( mass > 76. && mass < 106. );  
  return fabs(mass - 91.1876) < 15;
}


double projectedMet(int i_hyp){
  double DeltaPhi = nearestDeltaPhi(cms2.evt_tcmetPhi(),i_hyp);
  if (DeltaPhi < TMath::Pi()/2) return cms2.evt_tcmet()*TMath::Sin(DeltaPhi);
  return cms2.evt_tcmet();
}


double nearestDeltaPhi(double Phi, int i_hyp)
{
  double tightDPhi = fabs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi);
  tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
  double looseDPhi = fabs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi);
  looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
  return TMath::Min(tightDPhi, looseDPhi);
}


std::vector<LorentzVector> getJets(int type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag){
  std::vector<LorentzVector> jets;
  const double vetoCone = 0.3;
  double jec = 1.0;
  
  switch ( type ){
  case jptJet:
    for ( unsigned int i=0; i < cms2.jpts_p4().size(); ++i) {
      if(applyJEC)
	jec = jetCorrection(cms2.jpts_p4()[i], jet_corrector_jpt);
      if ( cms2.jpts_p4()[i].pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.jpts_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4()[i])) < vetoCone ) continue;
      jets.push_back(cms2.jpts_p4()[i] * jec);
    }
    break;
  case pfJet:
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
      if(applyJEC)
	jec = jetCorrection(cms2.pfjets_p4()[i], jet_corrector_pf);
      if ( cms2.pfjets_p4()[i].pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.pfjets_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
      jets.push_back(cms2.pfjets_p4()[i] * jec);
    }
    break;
  case GenJet:
    for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
      if ( cms2.genjets_p4()[i].pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.genjets_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4()[i])) < vetoCone ) continue;
      jets.push_back(cms2.genjets_p4()[i]);
    }
    break;
  case CaloJet:
    for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
      if ( cms2.jets_pat_jet_p4()[i].pt() < etThreshold ) continue; // note that this is already corrected
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.jets_pat_jet_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4()[i])) < vetoCone ) continue;
      jets.push_back(cms2.jets_pat_jet_p4()[i]);
    }
    break;
  case TrkJet:
    for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
      if(applyJEC)
	jec = jetCorrection(cms2.trkjets_p4()[i], jet_corrector_trk);
      if ( cms2.trkjets_p4()[i].pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.trkjets_p4()[i].eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[i])) < vetoCone ) continue;
      jets.push_back(cms2.trkjets_p4()[i] * jec);
    }
    break;
  default:
    std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
  }
  if ( sortJets ) std::sort(jets.begin(), jets.end(), comparePt);
  return jets;
}

bool comparePt(LorentzVector lv1, LorentzVector lv2) {
   return lv1.pt() > lv2.pt();
}

bool isGoodVertex(size_t ivtx) {
    if (cms2.vtxs_isFake()[ivtx]) return false;
    if (cms2.vtxs_ndof()[ivtx] < 4.) return false;
    if (cms2.vtxs_position()[ivtx].Rho() > 2.0) return false;
    if (fabs(cms2.vtxs_position()[ivtx].Z()) > 24.0) return false;
    return true;
}


double dzPV(const LorentzVector& vtx, const LorentzVector& p4, const LorentzVector& pv){
  return (vtx.z()-pv.z()) - ((vtx.x()-pv.x())*p4.x()+(vtx.y()-pv.y())*p4.y())/p4.pt() * p4.z()/p4.pt();
}





//
// Muon selections
//
bool goodMuonIsolated(unsigned int i){  
  bool ptcut = cms2.mus_p4().at(i).pt() >= 20.0;
  return ptcut && muonId(i, NominalWWV1);
}

bool ww_mud0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  double sumPtMax = -1;
  int iMax = -1;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  }
  if (iMax<0) return false;
  double dxyPV = cms2.mus_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.mus_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.mus_trk_p4()[index].phi());
  double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[iMax]);
  return fabs(dxyPV) < 0.02 && fabs(dzpv)<1.0;
}

double ww_muIsoVal(unsigned int index){
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/pt;
}

bool ww_muIso(unsigned int index){
  return ww_muIsoVal(index)<0.15;
}


bool ww_muId(unsigned int index){
  if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; //glb fit chisq  
  if (((cms2.mus_type().at(index)) & (1<<1)) == 0)    return false; // global muon
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
  if (cms2.mus_gfit_validSTAHits().at(index)==0 ) return false;
  if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
  if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
  if (cms2.mus_nmatches().at(index)<2) return false;
  return true;
}

//
// Electron selections
//

bool goodElectronIsolated(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 20.0;
  return ptcut && pass_electronSelection( i, electronSelection_wwV1);
}

bool ww_elId(unsigned int index){
  if( fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
      fabs(cms2.els_conv_dcot().at(index)) < 0.02) return false;
  if (! (electronId_VBTF(index, VBTF_35X_80) & (1<<ELEID_ID)) ) return false;
  if ( cms2.els_exp_innerlayers39X().at(index) > 0 ) return false;
  return true;
}

bool ww_eld0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  double sumPtMax = -1;
  int iMax = -1;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    if (!isGoodVertex(i)) continue;
    if ( cms2.vtxs_sumpt().at(i) > sumPtMax ){
      iMax = i;
      sumPtMax = cms2.vtxs_sumpt().at(i);
    }
  }
  if (iMax<0) return false;
  double dxyPV = cms2.els_d0()[index]-
    cms2.vtxs_position()[iMax].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.vtxs_position()[iMax].y()*cos(cms2.els_trk_p4()[index].phi());
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[iMax]);
  return fabs(dxyPV) < 0.02 && fabs(dzpv)<1.0;
}

double ww_elIsoVal(unsigned int index){
  float sum = cms2.els_tkIso().at(index);
  if ( fabs(cms2.els_etaSC()[index]) < 1.479 )
    // if ( fabs(cms2.els_p4().at(index).eta()) < 1.479) 
    sum += max(0., (cms2.els_ecalIso().at(index) -1.));
  else
    sum += cms2.els_ecalIso().at(index);
  sum += cms2.els_hcalIso().at(index);
  double pt = cms2.els_p4().at(index).pt();
  return sum/pt;
}

bool ww_elIso(unsigned int index){
  return ww_elIsoVal(index)<0.1;
}


unsigned int numberOfSoftMuons(int i_hyp, bool nonisolated,
			       const std::vector<LorentzVector>& vetojets)
{
  unsigned int nMuons = 0;
  for (int imu=0; imu < int(cms2.mus_charge().size()); ++imu) {
    // quality cuts
    // if (  ((cms2.mus_goodmask()[imu]) & (1<<14)) == 0 ) continue; // TMLastStationOptimizedLowPtTight
    if (  ((cms2.mus_goodmask()[imu]) & (1<<19)) == 0 ) continue; // TMLastStationAngTight
    if ( cms2.mus_p4()[imu].pt() < 3 ) continue;
    if ( TMath::Abs(ww_mud0PV(imu)) > 0.2) continue;
    if ( cms2.mus_validHits()[imu] < 11) continue;
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == imu ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == imu ) continue;
    if ( nonisolated && ww_muIsoVal(imu)<0.1 && cms2.mus_p4()[imu].pt()>20 ) continue;
    bool skip = false;
    for ( std::vector<LorentzVector>::const_iterator ijet = vetojets.begin();
	  ijet != vetojets.end(); ++ijet )
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(*ijet,cms2.mus_p4()[imu])) < 0.3 ) skip=true;
    if ( skip ) continue;
    ++nMuons;
  }
  return nMuons;
}

unsigned int numberOfExtraLeptons(int i_hyp, double minPt){
  unsigned int nMuons = 0;
  for (int i=0; i < int(cms2.mus_charge().size()); ++i) {
    // printf("Muon: %u, pt: %0.2f\n",i,cms2.mus_p4().at(i).pt());
    if ( cms2.mus_p4()[i].pt() < minPt ) continue;
    // printf("\tpassed minPt\n");
    if ( TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13 && cms2.hyp_lt_index()[i_hyp] == i ) continue;
    if ( TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13 && cms2.hyp_ll_index()[i_hyp] == i ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( ! (ww_mud0PV(i) && ww_muId(i) && ww_muIso(i)&&
	    fabs(cms2.mus_p4().at(i).eta()) <2.4) ) continue;
    // printf("\tpassed all\n");
    ++nMuons;
  }
  unsigned int nElectrons = 0;
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    // printf("Electron: %u, pt: %0.2f\n",i,cms2.els_p4().at(i).pt());
    if ( cms2.els_p4()[i].pt() < minPt ) continue;
    // printf("\tpassed minPt\n");
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.els_p4().at(i)) <0.1) ) continue;
    // printf("\tpassed hyp letpons\n");
    if ( !(ww_elId(i) && ww_eld0PV(i) && ww_elIso(i) && 
	   fabs(cms2.els_p4().at(i).eta()) < 2.5) ) continue;
    // printf("\tpassed all\n");
    ++nElectrons;
  }
  return nMuons+nElectrons;
}

bool toptag(int type, int i_hyp, double minPt){
  std::vector<LorentzVector> jets;
  const double vetoCone    = 0.3;

  switch ( type ){
  case pfJet:
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
      if ( cms2.pfjets_p4()[i].pt() < minPt ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4()[i])) < vetoCone ) continue;
      if ( cms2.pfjets_trackCountingHighEffBJetTag()[i]>2.1 ) return true;
    }
    break;
  case CaloJet:
    for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
      if ( cms2.jets_p4()[i].pt() < minPt ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[i])) < vetoCone ) continue;
      if ( cms2.jets_trackCountingHighEffBJetTag()[i]>2.1 ) return true;
    }
    break;
  default:
    std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
  }
  return false;
}


bool passedTriggerRequirements() {
  return true;
  bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele15_LW_L1R");
  //bool hlt_ele15_lw_l1r = cms2.passHLTTrigger("HLT_Ele17_SW_L1R");
  bool hltMu9           = cms2.passHLTTrigger("HLT_Mu9");
  // return hlt_ele15_lw_l1r || hltMu9;
  return true;
}

bool EventIdentifier::operator < (const EventIdentifier &other) const
{
     if (run != other.run)
	  return run < other.run;
     if (event != other.event)
	  return event < other.event;
     // the floating point numbers are not easy, because we're
     // comapring ones that are truncated (because they were written
     // to file and read back in) with ones that are not truncated.
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
	  return trks_d0 < other.trks_d0;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
	  return hyp_lt_pt < other.hyp_lt_pt;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
	  return hyp_lt_eta < other.hyp_lt_eta;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
	  return hyp_lt_phi < other.hyp_lt_phi;
     // if the records are exactly the same, then r1 is not less than
     // r2.  Duh!
     return false;
}

bool EventIdentifier::operator == (const EventIdentifier &other) const
{
     if (run != other.run)
	  return false;
     if (event != other.event)
	  return false;
     // the floating point numbers are not easy, because we're
     // comapring ones that are truncated (because they were written
     // to file and read back in) with ones that are not truncated.
     if (fabs(trks_d0 - other.trks_d0) > 1e-6 * trks_d0)
	  return false;
     if (fabs(hyp_lt_pt - other.hyp_lt_pt) > 1e-6 * hyp_lt_pt)
       return false;
     if (fabs(hyp_lt_eta - other.hyp_lt_eta) > 1e-6 * hyp_lt_eta)
       return false;
     if (fabs(hyp_lt_phi - other.hyp_lt_phi) > 1e-6 * hyp_lt_phi)
       return false;
     return true;
}


void CalculateFakeRateProb(){

  int size = cms2.genps_id().size();

  TDatabasePDG *pdg = new TDatabasePDG();
  // static TDatabasePDG *pdg = new TDatabasePDG();                                                                                          

  cout << "                " << "   pt    " << "  phi  " << "      eta   " << "    mass  "
       << "status " << "Mother  " << endl;
  std::cout << "---------------------------------------------------------------------" << std::endl;
  for (int j=0; j<size; j++) {
  float m2 = cms2.genps_p4().at(j).M2();
  float m = m2 >= 0 ? sqrt(m2) : 0.0;
  cout << setw(4) << left << j << " "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id().at(j))->GetName() << " "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).pt() << "  "
         << setw(7) << right << setprecision(4) << cms2.genps_p4().at(j).phi() << "  "
         << setw(10) << right << setprecision(4) << cms2.genps_p4().at(j).eta() << "  "
         << setw(7) << right << setprecision(4) << m << "  "
         << setw(4) << right << cms2.genps_status().at(j) << " "
         << setw(10) << left << pdg->GetParticle(cms2.genps_id_mother().at(j))->GetName()
         << " " << endl;
  }

}
double BTag(int type, unsigned int iJet){
  // find a jet in the jets list
  // that matches the current type jet
  LorentzVector jetP4;
  switch ( type ) {
  case jptJet:
    jetP4 = cms2.jpts_p4().at(iJet);
    break;
  case CaloJet:
    return cms2.jets_trackCountingHighEffBJetTag()[iJet];
    break;
  case pfJet:
    return cms2.pfjets_trackCountingHighEffBJetTag()[iJet];
    break;
  default:
    std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
    assert(0);
  }

  int refJet = -1;
  for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(jetP4,cms2.jets_p4()[i])) > 0.3 ) continue;
    refJet = i;
  }
  if (refJet == -1){
    // std::cout << "Warning: failed to find a matching jet for b-tagging." << std::endl; 
    return 0.0;
  }
  return cms2.jets_trackCountingHighEffBJetTag().at(refJet);
}

bool defaultBTag(int type, unsigned int iJet){
  return BTag(type,iJet)>2.1;
}

void FillEffHist(TString process, double weight) {
  
  if(!isIdentified(process)) return;
  // if(process != "WW") return;
  
  // Fill the lepton efficiency histograms
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    if (cms2.genps_status().at(i) != 3 || cms2.genps_id_mother().at(i) == 21212 ) continue;
    // consider only the leptons from the W and Z
    if (TMath::Abs(cms2.genps_id_mother().at(i)) != 24 && cms2.genps_id_mother().at(i) != 23) continue;
    // ==== Electron Efficiency
    if(TMath::Abs(cms2.genps_id().at(i)) == 11) {
      els_denom_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
      els_denom_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
      els_denom_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
      
      for(int i_els=0;i_els<cms2.els_charge().size();++i_els) {
	if( goodElectronIsolated(i_els) && cms2.els_mc3idx().at(i_els) == i) {
	  els_numer_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
	  els_numer_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
	  els_numer_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
	}
      }
    }    // ==== End of Electron Efficiency
    // ==== Muon Efficiency
    if(TMath::Abs(cms2.genps_id().at(i)) == 13) {
      mus_denom_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
      mus_denom_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
      mus_denom_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
      
      for(int i_mus=0;i_mus<cms2.mus_charge().size();++i_mus) {
	if( goodMuonIsolated(i_mus) && cms2.mus_mc3idx().at(i_mus) == i) {
	  mus_numer_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
	  mus_numer_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
	  mus_numer_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
	}
      }
    }   // ==== End of Muon Efficiency
  }
  // cout << "ScanChain:: " << __LINE__ << " els_denom_mc_pt->GetBinContent(10) =  " <<  els_denom_mc_pt->GetBinContent(10)  << endl;
  
}
void FillKtHist(TString process, double weight) {
  // if(!isIdentified(process)) return;
  // Currently only implemented for WW/HWW/WZ/ZZ
  LorentzVector systP4(0.,0.,0.,0.);
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    if (process=="WW" || process.Contains("HWW",TString::kExact)) {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24) 	systP4 += cms2.genps_p4().at(i);
    }
    else if (process=="WZ" || process == "ZZ") {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24 || cms2.genps_id().at(i) == 23 ) systP4 += cms2.genps_p4().at(i);
    }
  }
  kx_->Fill(systP4.Px(), weight);
  ky_->Fill(systP4.Py(), weight);
}

bool isIdentified( TString process) {
  // std::cout << "ScanChain::isIdentified()"<<endl;
  if(process == "WW" || process.Contains("HWW", TString::kExact)) 
    return getVVType()==0;
  else  if (process == "WZ")       return getVVType()==1;
  else  if (process == "ZZ")       return getVVType()==2;
  else  return true;
}



void ProcessSample(const char *process, std::vector<std::string> file_patterns, TFile *utilFile_,  int nEvents = -1, double IntLumi=100, double Xsect=1.0, int nProcessedEvents=-1, std::string skimFilePrefix="", bool realData=false, bool identifyEvents=false){

  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());
  
  ScanChain(process, tchain, utilFile_, nEvents, IntLumi, Xsect, nProcessedEvents, skimFilePrefix, realData, identifyEvents);

}

void ProcessSample(const char* process, std::string file_pattern, TFile *utilFile_, int nEvents = -1, double IntLumi=100, double Xsect=1.0, int nProcessedEvents=-1, std::string skimFilePrefix="", bool realData=false, bool identifyEvents=false){
  std::vector<std::string> vec;
  vec.push_back(file_pattern);
  ProcessSample(process, vec, utilFile_, nEvents, IntLumi, Xsect, nProcessedEvents, skimFilePrefix, realData, identifyEvents);
}

