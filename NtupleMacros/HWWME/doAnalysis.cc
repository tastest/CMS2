/* 
=== Description

  This is the looper for the H->WW analysis using the ME method, written by Sergo Yanyan and Kevin.
  The main goals of this looper are
  1. producing "Util.root" containing a bunch of histograms needed for ME calculation
  2. Producing Smurf ntuples (not synchronized with the main looper)

=== Notes
  1. This is the *slimmed* version of the WW looper
     http://www.t2.ucsd.edu/tastwiki/bin/view/CMS/HwwCode
  2. The jet-energy correction is simply the one from the CMS2 ntuple
*/

// C++
#include <iostream>
#include <vector>
#include "doAnalysis.h"

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
#include "Math/LorentzVector.h"
#include "TLorentzVector.h"
#include <vector>
#include "Math/VectorUtil.h"
#include "math.h"

// HWWME looper related
#include "doAnalysis.h"

// CMS2 related
#ifndef __CINT__
#include "TBitSet.hh"
#include "CORE/CMS2.h"
#include "CORE/electronSelections.h"
#include "CORE/muonSelections.h"
#include "CORE/jetSelections.h"
#include "CORE/metSelections.h"
#include "CORE/mcSelections.h"
#include "CORE/eventSelections.h"
#endif


using namespace tas;

bool useLHeleId = false;

TBitSet* _cutWord;
TBitSet* _cutMask;

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


static std::set<EventIdentifier> already_seen;

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



void ScanChain(const char* process, TChain *chain, TFile *utilFile_,  int nEvents, double IntLumi, double Xsect, int nProcessedEvents, std::string analysis, std::string skimFilePrefix, bool realData, bool identifyEvents){

  int eventCount[kNdilepFlav];
  double eventYield[kNdilepFlav];
  double weight;
  double Integral;
  
  _cutWord = new TBitSet(kNCuts);
  _cutMask = new TBitSet(kNCuts);
  _cutMask->SetAll();
  //_cutMask->SetFalse(kcut_zsel);
  
  already_seen.clear();
  
  for ( int count=0; count < kNdilepFlav; count++) {
    eventCount[count]=0;
    eventYield[count]=0;
  }
  
  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  if(!realData) InitMCUtilHist(process, utilFile_);
  
  // File Loop
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
      cms2.GetEntry(event);
      
      // if ( cms2.evt_event() != 100147) continue;
      // std::cout << "Run = " << cms2.evt_run() << ", lumi = " << cms2.evt_lumiBlock() << ", event = " << cms2.evt_event() << "\n";
      //dumpDocLines();
      
      // identifyEvents by MC truth if specified
      if ( (!realData && identifyEvents) && !isIdentified (process)) continue;

      // Calculate event weight for MC
      double mcweight;
      if (TString(process) != "data"){
	mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	weight = IntLumi * mcweight * (Xsect>0?Xsect:cms2.evt_xsec_excl()*cms2.evt_kfactor()) /
	  (nProcessedEvents>0 ? nProcessedEvents : cms2.evt_nEvts());
      } else weight = 1;
      
      // Get fake-rate related histograms
      if(!realData && TString(process) == "wjets") fillFOHist();
      
      // skip events without hypothesis
      unsigned int nHyps = cms2.hyp_type().size();
      if (nHyps == 0 ) continue;
      
      // hypothesis based duplicate event removal
      EventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) continue;
   
      
      // fill lepton efficiency
      if(!realData) 	fillEffHist(process, weight);
      
      // Start Event selection based on the hypothesis

      for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {

	if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) continue;
	if ( !isfinite(cms2.hyp_p4().at(i_hyp).mass2())) continue;

	_cutWord->SetAllBitsFalse();
	
	if (TString(analysis) == "HWW")
	  ApplyHWWEventSelection(i_hyp, realData);
	
	else if (TString(analysis) == "HZZ")
	  ApplyHZZEventSelection(i_hyp, realData);
	
	else
	  {
	    std::cout << "ERROR..analysis not being recognized...\n";
	    assert(0);
	  }

	int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
	bool accept = true;

	for (int j = 0; j < kNCuts; j++) {
	  if (_cutMask->IsTrue(j)) accept &= _cutWord->IsTrue (j);
        }

	// fill the system boost and smurf trees after the event selections
        if (accept){
	  const std::vector<LorentzVector>& jets =  getJets(pfJet, i_hyp, 30, 5.0, true, false);
	  if(!realData) fillKtHist(process, weight, jets.size());
	  eventCount[type]++;
	  eventYield[type]+=weight;
	  eventCount[all]++;
	  eventYield[all]+=weight;
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
    std::cout << "WARNING: number of events from files is not equal to total number of events" << std::endl;
  }


  for ( int count=0; count < kNdilepFlav; count++) {
    cout<<" Total Count in "<<count<<" is equal "<<eventCount[count]<<"\n";
    cout<<" Total Yield in "<<count<<" is equal "<<eventYield[count]<<"\n";
  }

  Integral=eventYield[all];
  
  cout << "Total yield " << eventYield[all] << ": MM "<<  eventYield[mumu] 
       << "; EM "<< eventYield[emu]
       << "; EE "<< eventYield[ee] <<endl;
    
  if(!realData) saveMCUtilOutput(process, utilFile_);
  cout<<"Total Events Before Selection "<<chain->GetEntries()<<"\n";
}


void ApplyHWWEventSelection( unsigned int i_hyp, bool realData){
  
  //std::cout << "doAnalysis::ApplyHWWEventSelection...\n";
  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return;
  if ( std::min(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<10 ) return;
  int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
  
  // no trigger requirements
  _cutWord->SetTrue(kcut_Trigger);

  // OS
  if ( fast_sign (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] ) < 0)  _cutWord->SetTrue(kcut_OS);

  // Require at least one good reconstructed primary vertex
  if (nGoodVertex() >= 1) _cutWord->SetTrue(kcut_GoodVertex);

  // di-lepton mass cut
  if ( cms2.hyp_p4()[i_hyp].mass() > 12 ) _cutWord->SetTrue(kcut_Mll); 

  // Z window veto
  if ( type == ee || type == mumu ) {
    if (!inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))   _cutWord->SetTrue(kcut_zsel);
  }
  if (type == emu) _cutWord->SetTrue(kcut_zsel);
  
  // == letpon ID and Isolation
  bool passedLTFinalRequirements = true;
  bool passedLLFinalRequirements = true;
  bool passedLTElFakableRequirements = true;
  bool passedLLElFakableRequirements = true;
  bool passedLTMuFakableRequirements = true;
  bool passedLLMuFakableRequirements = true;
  
  // muon selections

  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( !goodMuonIsolated(index))  passedLTFinalRequirements = false;
    if ( !fakableMuon(index) ) passedLTMuFakableRequirements = false;
    passedLTElFakableRequirements = false;
  }
  
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( !goodMuonIsolated(index) ) passedLLFinalRequirements = false;
    if ( !fakableMuon(index) ) passedLLMuFakableRequirements = false;
    passedLLElFakableRequirements = false;
  } 

  
  // electron selections
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( ! goodElectronIsolated(index) ) passedLTFinalRequirements = false;
    if ( !fakableElectron(index)) passedLTElFakableRequirements = false;
    passedLTMuFakableRequirements = false;
  }

  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( ! goodElectronIsolated(index) ) passedLLFinalRequirements = false;
    if ( !fakableElectron(index)) passedLLElFakableRequirements = false;
    passedLLMuFakableRequirements = false;
  }

  if ( passedLTFinalRequirements )     _cutWord->SetTrue(kcut_LT);
  if ( passedLLFinalRequirements )     _cutWord->SetTrue(kcut_LL);
  
  // MET cut
  if (passedMetRequirements(i_hyp)) _cutWord->SetTrue(kcut_met);
  
  // Jet Veto  
  const std::vector<LorentzVector>& jets =  getJets(pfJet, i_hyp, 30, 5.0, true, false);
  
  double dPhiDiLepJet1 = 0.0;
  if (jets.size() > 0)  
    dPhiDiLepJet1 = TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_p4().at(i_hyp),jets.at(0)));
  
  if ( jets.size() < 2 ) {
    if (type == emu || dPhiDiLepJet1 < 165.*TMath::Pi() / 180.)
      _cutWord->SetTrue(kcut_jetveto);
  }

  // ==Soft Muon Veto
  if ( numberOfSoftMuons(i_hyp,true) == 0) _cutWord->SetTrue(kcut_softmuonveto);
  
  // == Extra Lepton Veto  
  if ( numberOfExtraLeptons(i_hyp,10) == 0) _cutWord->SetTrue(kcut_extraleptonveto);
  
  // Top tagging
  if ( !toptag(CaloJet,i_hyp,0) ) _cutWord->SetTrue(kcut_toptag);
  
}

// HZZ selections

void ApplyHZZEventSelection( unsigned int i_hyp, bool realData){

  //std::cout << "doAnalysis::ApplyHZZEventSelection...\n";
  int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
  
  if ( type == emu ) return;
  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return;
  if ( std::min(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return;
  if ( cms2.hyp_p4().at(i_hyp).pt() < 25 ) return;

  // no trigger requirements
  _cutWord->SetTrue(kcut_Trigger);
  
  // OS
  if ( fast_sign (cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] ) < 0)  _cutWord->SetTrue(kcut_OS);
 
  // Require at least one good reconstructed primary vertex
  if (nGoodVertex() >= 1) _cutWord->SetTrue(kcut_GoodVertex);
  
  // di-lepton mass cut
  if ( cms2.hyp_p4()[i_hyp].mass() > 12 ) _cutWord->SetTrue(kcut_Mll); 
  
  // Z window veto
  if ( type == ee || type == mumu ) {
    if (inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))   _cutWord->SetTrue(kcut_zsel);
  }
  
  // == letpon ID and Isolation
  bool passedLTFinalRequirements = true;
  bool passedLLFinalRequirements = true;
  bool passedLTElFakableRequirements = true;
  bool passedLLElFakableRequirements = true;
  bool passedLTMuFakableRequirements = true;
  bool passedLLMuFakableRequirements = true;
  
  // muon selections

  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( !goodMuonIsolated(index))  passedLTFinalRequirements = false;
    if ( !fakableMuon(index) ) passedLTMuFakableRequirements = false;
    passedLTElFakableRequirements = false;
  }
  
  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 13){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( !goodMuonIsolated(index) ) passedLLFinalRequirements = false;
    if ( !fakableMuon(index) ) passedLLMuFakableRequirements = false;
    passedLLElFakableRequirements = false;
  } 

  
  // electron selections
  
  if (TMath::Abs(cms2.hyp_lt_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_lt_index()[i_hyp];
    if ( ! goodElectronIsolated(index) ) passedLTFinalRequirements = false;
    if ( !fakableElectron(index)) passedLTElFakableRequirements = false;
    passedLTMuFakableRequirements = false;
  }

  if (TMath::Abs(cms2.hyp_ll_id()[i_hyp]) == 11){
    unsigned int index = cms2.hyp_ll_index()[i_hyp];
    if ( ! goodElectronIsolated(index) ) passedLLFinalRequirements = false;
    if ( !fakableElectron(index)) passedLLElFakableRequirements = false;
    passedLLMuFakableRequirements = false;
  }

  if ( passedLTFinalRequirements )     _cutWord->SetTrue(kcut_LT);
  if ( passedLLFinalRequirements )     _cutWord->SetTrue(kcut_LL);

  // MET cut
  if ( metValue() > 60) _cutWord->SetTrue(kcut_met);
  
  // relax jet-veto  
  const std::vector<LorentzVector>& jets =  getJets(pfJet, i_hyp, 30, 5.0, true, false);
  if ( jets.size() < 2 ) 
    _cutWord->SetTrue(kcut_jetveto);
  
  // ==Soft Muon Veto
  if ( numberOfSoftMuons(i_hyp,true) == 0) _cutWord->SetTrue(kcut_softmuonveto);
  
  // == Extra Lepton Veto  
  if ( numberOfExtraLeptons(i_hyp,10) == 0) _cutWord->SetTrue(kcut_extraleptonveto);
  
  // Top tagging
  if ( !toptag(CaloJet,i_hyp,0) ) _cutWord->SetTrue(kcut_toptag);
  
}


// == Utility fucntions

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

double projectedMet(unsigned int i_hyp, double met, double phi)
{
  double DeltaPhi = nearestDeltaPhi(phi,i_hyp);
  if (DeltaPhi < TMath::Pi()/2) return met*TMath::Sin(DeltaPhi);
  return met;
}

double metValue(){    return cms2.evt_pfmet(); }
double metPhiValue(){ return cms2.evt_pfmetPhi(); }

bool passedMetRequirements(unsigned int i_hyp){

  int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
  metStruct trkMET = trackerMET(i_hyp,0.1); //,&jets);
  double pMet = std::min(projectedMet(i_hyp, metValue(), metPhiValue()),
			 projectedMet(i_hyp, trkMET.met, trkMET.metphi));
  if ( pMet < 20 ) return false;
  if (type == ee || type == mumu) {
    if ( pMet < 40 ) return false;
  }

  return true;
}

bool passedHZZMetRequirements(unsigned int i_hyp){
  metStruct trkMET = trackerMET(i_hyp,0.1); //,&jets);
  double met = std::min( metValue(), double(trkMET.met));
  if ( met < 50 ) return false;
  return true;
}




double nearestDeltaPhi(double Phi, int i_hyp)
{
  double tightDPhi = fabs(cms2.hyp_lt_p4()[i_hyp].Phi() - Phi);
  tightDPhi = std::min(2*TMath::Pi() - tightDPhi, tightDPhi);
  double looseDPhi = fabs(cms2.hyp_ll_p4()[i_hyp].Phi() - Phi);
  looseDPhi = std::min(2*TMath::Pi() - looseDPhi, looseDPhi);
  return TMath::Min(tightDPhi, looseDPhi);
}

// only developped for the pfjet and genjet
std::vector<LorentzVector> 
getJets(int type, int i_hyp, double etThreshold, double etaMax, bool sortJets, bool btag)
{
     std::vector<LorentzVector> jets;
     const double vetoCone = 0.3;
     double jec = 1.0;
     switch ( type ){
     case pfJet:
       for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {  
	 // take the default jec from the ntuple
	 jec = cms2.pfjets_cor()[i];
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
     default:
       std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
     }
     if ( sortJets ) std::sort(jets.begin(), jets.end(), comparePt);
     return jets;
}


bool comparePt(LorentzVector lv1, LorentzVector lv2) {
   return lv1.pt() > lv2.pt();
}

unsigned int nGoodVertex() {
  unsigned int nVtx = 0;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    // if (cms2.vtxs_isFake()[i]) continue;
    if (!isGoodVertex(i)) continue;
    nVtx++;
  }
  return nVtx;
}

int primaryVertex(){
  return 0;
}

//
// Muon selections
//

bool goodMuonIsolated(unsigned int i){  
  bool ptcut = cms2.mus_p4().at(i).pt() >= 10.0;
  //bool core = ptcut && muonId(i, NominalSmurfV5);
  bool internal = ww_muBase(i) && ww_mud0PV(i) && ww_mudZPV(i) && ww_muId(i) && ww_muIso(i); 
  //return ptcut && core;
  return ptcut && internal;
}

bool fakableMuon(unsigned int i){
  bool ptcut = cms2.mus_p4().at(i).pt() >= 10.0;
  // return ptcut && muonId(i, muonSelectionFO_mu_wwV1_iso10);
  return ptcut && ww_mudZPV(i) && muonId(i, muonSelectionFO_mu_smurf_04);
}

bool ww_muBase(unsigned int index){
  if (cms2.mus_p4().at(index).pt() < 10.0) return false;
  if (fabs(cms2.mus_p4().at(index).eta()) > 2.4) return false;
  if (cms2.mus_type().at(index) == 8) return false; // not STA
  return true;
}

double ww_mud0ValuePV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return 9999;
  double dxyPV = cms2.mus_d0()[index]-
    cms2.vtxs_position()[vtxIndex].x()*sin(cms2.mus_trk_p4()[index].phi())+
    cms2.vtxs_position()[vtxIndex].y()*cos(cms2.mus_trk_p4()[index].phi());
  return fabs(dxyPV);
}

bool ww_mud0PV(unsigned int index){
  if ( cms2.mus_p4().at(index).pt() < 20. ) return ww_mud0ValuePV(index) < 0.01;
  return ww_mud0ValuePV(index) < 0.02;
}

bool ww_mudZPV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  //double dzpv = cms2.mus_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.mus_vertex_p4()[index], cms2.mus_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
  return fabs(dzpv)<0.1;
}

bool ww_muId(unsigned int index){ 
  if (((cms2.mus_type().at(index)) & (1<<2)) == 0)    return false; // tracker muon
  if (cms2.mus_validHits().at(index) < 11)            return false; // # of tracker hits
  if (cms2.mus_ptErr().at(index)/cms2.mus_p4().at(index).pt()>0.1) return false;
  if (cms2.trks_valid_pixelhits().at(cms2.mus_trkidx().at(index))==0) return false;
  // global muon
  bool goodMuonGlobalMuon = false;
  if (((cms2.mus_type().at(index)) & (1<<1)) == (1<<1)){
    goodMuonGlobalMuon = true;
    if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) goodMuonGlobalMuon = false; //glb fit chisq
    if (cms2.mus_gfit_validSTAHits().at(index)==0 ) goodMuonGlobalMuon = false;
    if (cms2.mus_nmatches().at(index)<2) goodMuonGlobalMuon = false;
  }
  return goodMuonGlobalMuon || 
    cms2.mus_pid_TMLastStationTight().at(index) == 1; // TM id
}

double ww_muIsoVal(unsigned int index){
  double sum =  cms2.mus_iso03_sumPt().at(index) +
    cms2.mus_iso03_emEt().at(index)  +
    cms2.mus_iso03_hadEt().at(index);
  double pt  = cms2.mus_p4().at(index).pt();
  return sum/pt;
}

bool ww_muIso(unsigned int index){
  if (cms2.mus_p4().at(index).pt()>20) {
    if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
      return muonIsoValuePF(index,0,0.3) < 0.13;
    else 
      return muonIsoValuePF(index,0,0.3) < 0.09;
  } else {
    if (TMath::Abs(cms2.mus_p4()[index].eta())<1.479) 
      return muonIsoValuePF(index,0,0.3) < 0.06;
    else 
      return muonIsoValuePF(index,0,0.3) < 0.05;
  }
}


//
// Electron selections
//
bool goodElectronIsolated(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 10.0;
  bool internal = ww_elBase(i) && ww_elId(i) && ww_eld0PV(i) && ww_eldZPV(i) && ww_elIso(i);
  //return ptcut && pass_electronSelection( i, electronSelection_smurfV5);
  return ptcut && internal;
}

bool fakableElectron(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 10.0;
  // extrapolate in partial id, iso and d0
  //return ptcut && pass_electronSelection( i, electronSelectionFO_el_wwV1_v2);
  return ptcut && ww_eldZPV(i) && pass_electronSelection( i, electronSelectionFO_el_smurf_v4);
}
bool ww_elBase(unsigned int index){
  if (cms2.els_p4().at(index).pt() < 10.0) return false;
  if (fabs(cms2.els_p4().at(index).eta()) > 2.5) return false;
  return true;
}
bool ww_elId(unsigned int index){
  // if( fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
  //     fabs(cms2.els_conv_dcot().at(index)) < 0.02) return false;
  // if (! (electronId_VBTF(index, VBTF_35X_80) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_VBTF(index, VBTF_35X_70) & (1<<ELEID_ID)) ) return false;
  // if (! (electronId_CIC(index, 4, CIC_SUPERTIGHT) & (1<<ELEID_ID)) ) return false;

  if (useLHeleId) {
    if (cms2.els_p4().at(index).pt()>20 && (passLikelihoodId(index,cms2.els_lh().at(index),90) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false; 
    if (cms2.els_p4().at(index).pt()<20 && (passLikelihoodId(index,cms2.els_lh().at(index),80) & (1<<ELEID_ID))!=(1<<ELEID_ID) ) return false;
  } else {
    if (! pass_electronSelection(index, electronSelection_smurfV3_id, false, false) ) return false;
  }

  // MIT conversion
  if ( isFromConversionMIT(index) ) return false;

  // conversion rejection - hit based
  if ( cms2.els_exp_innerlayers().at(index) > 0 ) return false;
  // MIT conversion
  // if (! pass_electronSelection(index, (1ll<<ELENOTCONV_MIT), false, false) ) return false;
  // if ( cms2.els_exp_innerlayers39X().at(index) > 0 ) return false;
  //  int ctfIndex = cms2.els_trkidx().at(index);
  // if ( ctfIndex >=0 && 
  //     cms2.els_charge().at(index)!=cms2.trks_charge().at(ctfIndex) ) return false;
  // if ( !electronId_smurf_v2(index) ) return false;
  
  return true;
}

bool ww_eld0PV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  double dxyPV = cms2.els_d0()[index]-
    cms2.vtxs_position()[vtxIndex].x()*sin(cms2.els_trk_p4()[index].phi())+
    cms2.vtxs_position()[vtxIndex].y()*cos(cms2.els_trk_p4()[index].phi());
  return fabs(dxyPV) < 0.02;
}

bool ww_eldZPV(unsigned int index){
  int vtxIndex = primaryVertex();
  if (vtxIndex<0) return false;
  // double dzPV = cms2.els_z0corr()[index]-cms2.vtxs_position()[iMax].z();
  double dzpv = dzPV(cms2.els_vertex_p4()[index], cms2.els_trk_p4()[index], cms2.vtxs_position()[vtxIndex]);
  return fabs(dzpv)<0.1;
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
  return pass_electronSelection(index, electronSelection_smurfV5_iso);
  //return ww_elIsoVal(index)<0.1;
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
    if ( cms2.mus_p4().at(i).pt() < minPt ) continue;
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
    if ( cms2.els_p4().at(i).pt() < minPt ) continue;
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
      if ( cms2.pfjets_p4().at(i).pt() < minPt ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4().at(i))) < vetoCone ) continue;
      if ( cms2.pfjets_trackCountingHighEffBJetTag().at(i)>2.1 ) return true;
    }
    break;
  case CaloJet:
    for ( unsigned int i=0; i < cms2.jets_p4().size(); ++i) {
      if ( cms2.jets_p4().at(i).pt() < minPt ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4().at(i))) < vetoCone ) continue;
      if ( cms2.jets_trackCountingHighEffBJetTag().at(i)>2.1 ) return true;
    }
    break;
  default:
    std::cout << "ERROR: not supported jet type is requested: " << type << " FixIt!" << std::endl;
  }
  return false;
}



//
// Triger
//
bool passedTrigger(TString trigName) {
  if ( find(cms2.hlt_trigNames().begin(), cms2.hlt_trigNames().end(), trigName)
       == cms2.hlt_trigNames().end() ) return false;
  return cms2.passHLTTrigger(trigName);
}

bool passedTriggerRequirements() {
  if ( passedTrigger("HLT_Mu9") ) return true;
  if ( passedTrigger("HLT_Mu15_v1") ) return true;
  if ( passedTrigger("HLT_Ele10_LW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_LW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_SW_L1R") ) return true;
  if ( passedTrigger("HLT_Ele15_SW_CaloEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_CaloEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TightEleId_L1R") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v2") ) return true;
  if ( passedTrigger("HLT_Ele17_SW_TighterEleIdIsol_L1R_v3") ) return true;
  if ( passedTrigger("HLT_DoubleEle17_SW_L1R_v1") ) return true;
  if ( passedTrigger("HLT_DoubleEle15_SW_L1R_v1") ) return true;
  if ( passedTrigger("HLT_DoubleEle10_SW_L1R") ) return true;
  return false;
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
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(jetP4,cms2.jets_p4().at(i))) > 0.3 ) continue;
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



void findClosestEleFO(LorentzVector v_parton, double& minDR, int& idx_minDR) {
  minDR = 999.0;
  idx_minDR = -999;
  
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    if ( leptonIsFromW(i, 11, false) > 0) continue;
    if (cms2.els_p4().at(i).pt() < 10) continue; 
    if (!fakableElectron(i)) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton, cms2.els_p4().at(i))) < minDR ) {
      minDR = TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton,cms2.els_p4().at(i)));
      idx_minDR = i;
    }
  }
}

void findClosestMuFO(LorentzVector v_parton, double& minDR, int& idx_minDR) {
  minDR = 999.0;
  idx_minDR = -999;
  
  for (int i=0; i < int(cms2.mus_charge().size()); ++i) {  
    if ( leptonIsFromW(i, 13, false) > 0) continue;
    if (cms2.mus_p4().at(i).pt() < 10 ) continue; 
    if (!fakableMuon(i)) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton, cms2.mus_p4().at(i))) < minDR ) {
      minDR = TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton,cms2.mus_p4().at(i)));
      idx_minDR = i;
    }
  }
}

bool isIdentified(const char* process) {
  if(TString(process) == "ww")
    return isWW();
  else  if (TString(process) == "wz")       return isWZ();
  else  if (TString(process) == "zz")       return isZZ();
  else  return true;
}


void ProcessSample(const char* process, std::vector<std::string> file_patterns, TFile *utilFile_,  int nEvents, double IntLumi, double Xsect, int nProcessedEvents, std::string analysis, std::string skimFilePrefix, bool realData, bool identifyEvents){

  TChain *tchain = new TChain("Events");
  for ( std::vector<std::string>::const_iterator pattern = file_patterns.begin();
	pattern != file_patterns.end(); ++pattern )
    tchain->Add(pattern->c_str());
  ScanChain(process, tchain, utilFile_, nEvents, IntLumi, Xsect, nProcessedEvents, analysis, skimFilePrefix, realData, identifyEvents);

}

void ProcessSample(const char* process, std::string file_pattern, TFile *utilFile_, int nEvents, double IntLumi, double Xsect, int nProcessedEvents, std::string analysis, std::string skimFilePrefix, bool realData, bool identifyEvents){
  std::vector<std::string> vec;
  vec.push_back(file_pattern);
  ProcessSample(process, vec, utilFile_, nEvents, IntLumi, Xsect, nProcessedEvents, analysis, skimFilePrefix, realData, identifyEvents);
}

void findClosestGenPs(LorentzVector v_parton, double& minDR, int& idx_minDR) {
  minDR = 999.0;
  idx_minDR = -999;
  
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    if (cms2.genps_p4().at(i).pt() < 10) continue; 
    if ( TMath::Abs(cms2.genps_id().at(i)) == 11 || TMath::Abs(cms2.genps_id().at(i)) == 13 ||  TMath::Abs(cms2.genps_id().at(i)) == 15 || cms2.genps_id().at(i) == 22 ) {
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton, cms2.genps_p4().at(i))) < minDR ) {
	minDR = TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton,cms2.genps_p4().at(i)));
	idx_minDR = i;
      }
    }
  }
}


//
// Histograms
// 
//  Book histograms...


// ==  Histograms to calculate real lepton efficiencies
TH2F  *els_numer_mc_;
TH2F  *els_denom_mc_;
TH2F  *els_eff_mc_;

TH1F  *els_numer_mc_eta_;
TH1F  *els_denom_mc_eta_;
TH1F  *els_eff_mc_eta_;

TH1F  *els_numer_mc_pt_;
TH1F  *els_denom_mc_pt_;
TH1F  *els_eff_mc_pt_;

TH2F  *mus_numer_mc_;
TH2F  *mus_denom_mc_;
TH2F  *mus_eff_mc_;

TH1F  *mus_numer_mc_eta_;
TH1F  *mus_denom_mc_eta_;
TH1F  *mus_eff_mc_eta_;

TH1F  *mus_numer_mc_pt_;
TH1F  *mus_denom_mc_pt_;
TH1F  *mus_eff_mc_pt_;

// Histograms for System boost
// by default these correpsond to 0-jet
TH1F *kx_;
TH1F *ky_;
TH1F *kt_;
TH2F *boost_;

TH1F *kx_1jet_;
TH1F *ky_1jet_;
TH1F *kt_1jet_;
TH2F *boost_1jet_;

TH1F *kx_2jet_;
TH1F *ky_2jet_;
TH1F *kt_2jet_;
TH2F *boost_2jet_;


// Histograms to calculate the probablitiy of a parton -> lepton
TH2F *parton_;
TH2F *els_fake_;
TH2F *els_genfr_;
TH2F *mus_fake_;
TH2F *mus_genfr_;

TH1F *parton_eta_;
TH1F *els_fake_eta_;
TH1F *els_genfr_eta_;
TH1F *mus_fake_eta_;
TH1F *mus_genfr_eta_;

TH1F *parton_pt_;
TH1F *els_fake_pt_;
TH1F *els_genfr_pt_;
TH1F *mus_fake_pt_;
TH1F *mus_genfr_pt_;

// electron - letpon FO response 

TH2F *els_fo_parton_;
TH2F *mus_fo_parton_;

// Histograms to calculate the probablitiy of a lepton FO -> lepton

TH2F *els_fo_;
TH2F *els_good_;
TH2F *els_fr_;

TH2F *mus_fo_;
TH2F *mus_good_;
TH2F *mus_fr_;


// Not clean area, diagonistic histograms

TH1F *els_fo_parton_dEta_;
TH1F *els_fo_parton_dPhi_;
TH1F *els_fo_MCType_;
TH1F *els_fo_partonID_;
TH2F *els_fo_MCType_pt_;
TH2F *els_fo_partonID_pt_;

TH1F *mus_fo_parton_dEta_;
TH1F *mus_fo_parton_dPhi_;
TH1F *mus_fo_MCType_;
TH1F *mus_fo_partonID_;
TH2F *mus_fo_MCType_pt_;
TH2F *mus_fo_partonID_pt_;

TH1F *dR_parton_lepton_;

void InitMCUtilHist(const char* process, TFile *utilFile_) {
  cout << "InitMCUtilHist(): " << " process = " << process <<endl; 
  utilFile_->cd();  
  gROOT->cd();

  if( TString(process) == "ww" || TString(process) == "zz") {
  
    els_numer_mc_ = new TH2F(Form("%s_heleNumer", process), "electron numerators", 20, -2.5, 2.5, 20, 0, 100);
    els_denom_mc_ = new TH2F(Form("%s_heleDenom",process),  "electron denominators", 20, -2.5, 2.5, 20, 0, 100);
    els_eff_mc_ = new TH2F(Form("%s_heleEff",process),  "Electron Efficiency", 20, -2.5, 2.5, 20, 0, 100);
    
    els_numer_mc_eta_ = new TH1F(Form("%s_heleNumerEta",process), "electron numerators", 20, -2.5, 2.5);
    els_denom_mc_eta_ = new TH1F(Form("%s_heleDenomEta",process), "electron denominators", 20, -2.5, 2.5);
    els_eff_mc_eta_ = new TH1F(Form("%s_heleEffEta",process), "Electron Efficiency", 20, -2.5, 2.5);
    
    els_numer_mc_pt_ = new TH1F(Form("%s_heleNumerPt",process), "electron numerators", 20, 0, 100);
    els_denom_mc_pt_ = new TH1F(Form("%s_heleDenomPt",process), "electron denominators", 20, 0, 100);
    els_eff_mc_pt_ = new TH1F(Form("%s_heleEffPt",process), "Electron Efficiency", 20, 0, 100);
    
    mus_numer_mc_ = new TH2F(Form("%s_hmuNumer",process), "muon numerators", 20, -2.5, 2.5, 20, 0, 100);
    mus_denom_mc_ = new TH2F(Form("%s_hmuDenom",process), "muon denominators", 20, -2.5, 2.5, 20, 0, 100);
    mus_eff_mc_ = new TH2F(Form("%s_hmuEff",process), "Muon Efficiency", 20, -2.5, 2.5, 20, 0, 100);
    
    mus_numer_mc_eta_ = new TH1F(Form("%s_hmuNumerEta",process), "muon numerators", 20, -2.5, 2.5);
    mus_denom_mc_eta_ = new TH1F(Form("%s_hmuDenomEta",process), "muon denominators", 20, -2.5, 2.5);
    mus_eff_mc_eta_ = new TH1F(Form("%s_hmuEffEta",process),  "Muon Efficiency", 20, -2.5, 2.5);
    
    mus_numer_mc_pt_ = new TH1F(Form("%s_hmuNumerPt",process), "muon numerators", 20, 0, 100);
    mus_denom_mc_pt_ = new TH1F(Form("%s_hmuDenomPt",process), "muon denominators", 20, 0, 100);
    mus_eff_mc_pt_ = new TH1F(Form("%s_hmuEffPt",process), "Muon Efficiency", 20, 0, 100);
    
  }
  
  // boost in the 0-jet bin 
  kx_ = new TH1F(Form("%s_kx",process), "System Boost in X", 50, -50, 50);
  ky_ = new TH1F(Form("%s_ky",process), "System Boost in Y", 50, -50, 50);
  kt_ = new TH1F(Form("%s_kt",process), "System Boost", 50, 0, 50);
  boost_ = new TH2F(Form("%s_boost",process), "System Boost kX vs kY", 50, -50, 50, 50, -50, 50);
  
  // boost in the 1-jet bin
  kx_1jet_ = new TH1F(Form("%s_kx_1jet",process), "System Boost in X ( 1-Jet)", 50, -100, 100);
  ky_1jet_ = new TH1F(Form("%s_ky_1jet",process), "System Boost in Y (1-Jet)", 50, -100, 100);
  kt_1jet_ = new TH1F(Form("%s_kt_1jet",process), "System Boost (1-Jet)", 50, 0, 100);
  boost_1jet_ = new TH2F(Form("%s_boost_1jet",process), "System Boost (1-Jet)", 50, -100, 100, 50, -100, 100);

  // boost in the 1-jet bin
  kx_2jet_ = new TH1F(Form("%s_kx_2jet",process), "System Boost in X ( 2-Jet)", 50, -100, 100);
  ky_2jet_ = new TH1F(Form("%s_ky_2jet",process), "System Boost in Y (2-Jet)", 50, -100, 100);
  kt_2jet_ = new TH1F(Form("%s_kt_2jet",process), "System Boost (2-Jet)", 50, 0, 100);
  boost_2jet_ = new TH2F(Form("%s_boost_2jet",process), "System Boost (2-Jet)", 50, -100, 100, 50, -100, 100);


  if (TString(process) == "wjets") {
    const Double_t ptbins_genfr[9] = {10.,20.,25.,30.,35., 40., 50., 75., 100.};
    const int nptbins_genfr = 8;
    const Double_t etabins_genfr[5] = {0.0, 1.0, 1.479, 2.0, 2.5};
    const int netabins_genfr = 4;
     
    // parton to lepton FO

    parton_ = new TH2F(Form("%s_hparton",process), "Generator Parton", nptbins_genfr,ptbins_genfr,netabins_genfr,etabins_genfr);
    els_fake_ = new TH2F(Form("%s_heleFake",process), "FO Electrons",  nptbins_genfr,ptbins_genfr,netabins_genfr,etabins_genfr);
    mus_fake_ = new TH2F(Form("%s_hmuFake",process), "FO Muons", nptbins_genfr,ptbins_genfr,netabins_genfr,etabins_genfr);

    els_genfr_ = new TH2F(Form("%s_heleGenFR",process), "Parton to Electron Generator FR",  nptbins_genfr,ptbins_genfr,netabins_genfr,etabins_genfr);
    els_fo_parton_ = new TH2F(Form("%s_heleFOResponse",process), "matched parton pt / Electron FO pT", nptbins_genfr,ptbins_genfr, 10, 0.5, 4.5);
    
    mus_genfr_ = new TH2F(Form("%s_hmuGenFR",process), "Parton to Muon Generator FR", nptbins_genfr,ptbins_genfr,netabins_genfr,etabins_genfr);
    mus_fo_parton_ = new TH2F(Form("%s_hmuFOResponse",process), "matched parton pt / Muon FO pT", nptbins_genfr,ptbins_genfr, 10, 0.5, 4.5);
    

    // lepton FO to lepton 
    // due to the limited stat., truncate at 35
    const Double_t ptbins_fr[6] = {10.,15.,20.,25.,30.,35.};
    const int nptbins_fr = 5;

    els_fo_ = new TH2F(Form("%s_heleFO",process), "Electron FO", nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);
    els_good_ = new TH2F(Form("%s_heleGood",process), "Electron Passing Analysis Cuts", nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);
    els_fr_ = new TH2F(Form("%s_heleFR",process), "Electron FO to Good Electron FR",  nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);

    mus_fo_ = new TH2F(Form("%s_hmuFO",process), "Muon FO", nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);
    mus_good_ = new TH2F(Form("%s_hmuGood",process), "Muon Passing Analysis Cuts", nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);
    mus_fr_ = new TH2F(Form("%s_hmuFR",process), "Muon FO to Good Muon FR",  nptbins_fr,ptbins_fr,netabins_genfr,etabins_genfr);
	
    // 1-d parton - lepton FO histograms
    parton_eta_ = new TH1F(Form("%s_hpartonEta",process), "Generator Parton Eta", netabins_genfr, etabins_genfr);
    els_fake_eta_ = new TH1F(Form("%s_heleFakeEta",process), "FO electrons", netabins_genfr, etabins_genfr);
    els_genfr_eta_ = new TH1F(Form("%s_heleGenFREta",process), "Parton to Electron Generator FR", netabins_genfr, etabins_genfr);
    mus_fake_eta_ = new TH1F(Form("%s_hmuFakeEta",process), "FO muons", netabins_genfr, etabins_genfr);
    mus_genfr_eta_ = new TH1F(Form("%s_hmuGenFREta",process), "Parton to Muon Generator FR", netabins_genfr, etabins_genfr);
    
    parton_pt_ = new TH1F(Form("%s_hpartonPt",process), "Generator Parton Pt", netabins_genfr,etabins_genfr);
    els_fake_pt_ = new TH1F(Form("%s_heleFakePt",process), "FO electrons",netabins_genfr,etabins_genfr);
    els_genfr_pt_ = new TH1F(Form("%s_heleGenFRPt",process), "Parton to Electron Generator FR", netabins_genfr,etabins_genfr);
    mus_fake_pt_ = new TH1F(Form("%s_hmuFakePt",process), "FO muons", netabins_genfr,etabins_genfr);
    mus_genfr_pt_ = new TH1F(Form("%s_hmuGenFRPt",process), "Parton to Muon Generator FR", netabins_genfr,etabins_genfr);
     
    // Diagonistic histograms
    els_fo_parton_dEta_ = new TH1F(Form("%s_hdEtaEleFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    els_fo_parton_dPhi_ = new TH1F(Form("%s_hdPhiEleFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    els_fo_MCType_ = new TH1F(Form("%s_heleFakeMCType",process), "Electron FO MC Type", 7,-1,6);
    els_fo_partonID_ = new TH1F(Form("%s_helePartonID",process), "Electron FO Matched Parton ID",35, -10, 25);
    els_fo_MCType_pt_ = new TH2F(Form("%s_heleFakeMCTypeVsPt",process), "Electron FO MC Type",nptbins_genfr,ptbins_genfr, 7,-1,6);
    els_fo_partonID_pt_ = new TH2F(Form("%s_helePartonIDVsPt",process), "Electron FO Matched Parton ID",nptbins_genfr,ptbins_genfr, 35, -10, 25);

    mus_fo_parton_dEta_ = new TH1F(Form("%s_hdEtaMuFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    mus_fo_parton_dPhi_ = new TH1F(Form("%s_hdPhiMuFOParton",process), "#Delta#eta(Lepton FO, Parton)", 20,0,0.2);
    mus_fo_MCType_ = new TH1F(Form("%s_hmuFakeMCType",process), "Muctron FO MC Type", 7,-1,6);
    mus_fo_partonID_ = new TH1F(Form("%s_hmuPartonID",process), "Muctron FO Matched Parton ID",35, -10, 25);
    mus_fo_MCType_pt_ = new TH2F(Form("%s_hmuFakeMCTypeVsPt",process), "Muon FO MC Type",nptbins_genfr,ptbins_genfr, 7,-1,6);
    mus_fo_partonID_pt_ = new TH2F(Form("%s_hmuPartonIDVsPt",process), "Muon FO Matched Parton ID",nptbins_genfr,ptbins_genfr, 35, -10, 25);

    // mininum distrance between a status parton and status 3 leptons/photons
    dR_parton_lepton_ = new TH1F(Form("%s_hdRPartonLepton",process), "Min dR(Parton,lepton or photon)", 50, 0, 5);

  }
}

void fillKtHist(const char* process, double weight, const int njet) {
  //std::cout << "doAnalysis::fillKtHist()..\n";
  // define the system boost
  LorentzVector systP4(0.,0.,0.,0.);
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    //std::cout << "doAnalysis::fillKtHist() i = " << i << " at line "<< __LINE__ << ".\n";
    if (TString(process).Contains("ww",TString::kExact)) {
      //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
      if(TMath::Abs(cms2.genps_id().at(i)) == 24) 	systP4 += cms2.genps_p4().at(i);
      //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
    }
    else if (TString(process) == "wz") {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24 || cms2.genps_id().at(i) == 23 ) systP4 += cms2.genps_p4().at(i);
    }
    else if ( TString(process).Contains("zz", TString::kExact)) {
      if( cms2.genps_id().at(i) == 23 ) systP4 += cms2.genps_p4().at(i);
    }
    else if(TString(process).Contains("wjets", TString::kExact)) {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24) systP4 += cms2.genps_p4().at(i);
    }
  }
  //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
  // add the leading status(3) parton 4-momemtum as well for wjets
  if(TString(process) == "wjets") {
    unsigned int idx_maxpt = -1;
    double maxpt = -1.0;
    for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
      // selecting status = 3 colored objects
      if (cms2.genps_status().at(i) != 3 || (cms2.genps_id().at(i) != 21 && TMath::Abs(cms2.genps_id().at(i)) > 8 )) continue;
      if( maxpt < cms2.genps_p4().at(i).pt() ) {
	maxpt = cms2.genps_p4().at(i).pt();
	idx_maxpt = i;
      }
    }
    if( idx_maxpt >= 0) systP4 += cms2.genps_p4().at(idx_maxpt);
  }

  if ( njet == 0) {
    kx_->Fill(systP4.Px(), weight);
    ky_->Fill(systP4.Py(), weight);
    kt_->Fill(sqrt(systP4.Px()*systP4.Px()+systP4.Py()*systP4.Py()), weight);
    boost_->Fill(systP4.Px(), systP4.Py());
    //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
  }
  else if ( njet == 1) {
    kx_1jet_->Fill(systP4.Px(), weight);
    ky_1jet_->Fill(systP4.Py(), weight);
    kt_1jet_->Fill(sqrt(systP4.Px()*systP4.Px()+systP4.Py()*systP4.Py()), weight);
    boost_1jet_->Fill(systP4.Px(), systP4.Py());
    //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
  }
  else if ( njet >= 2) {
    kx_2jet_->Fill(systP4.Px(), weight);
    ky_2jet_->Fill(systP4.Py(), weight);
    kt_2jet_->Fill(sqrt(systP4.Px()*systP4.Px()+systP4.Py()*systP4.Py()), weight);
    boost_2jet_->Fill(systP4.Px(), systP4.Py());
    //std::cout << "doAnalysis::fillKtHist()" << __LINE__ << ".\n";
  }
}


void getEff(double & numer, double & denom, double & eff, double & efferr ) 
{
  if (denom == 0.0) return;
  eff = numer/denom;
  efferr = sqrt(eff*(1-eff)/denom);
}


void fill1DEffHist(TH1F* hist_numer, TH1F* hist_denom, TH1F* hist_eff) {
  // cout << "fill1DEffHist()" << endl;
  if(!hist_numer || !hist_denom ) return;
  // check the hist_numer and hist_denom have the same bin structures
  if(    hist_numer->GetNbinsX() != hist_denom->GetNbinsX() 
	 || hist_numer->GetXaxis()->GetXmin() != hist_denom->GetXaxis()->GetXmin()
	 || hist_numer->GetXaxis()->GetXmax() != hist_denom->GetXaxis()->GetXmax() ) return;
 
  for(int iX=0;iX<hist_numer->GetNbinsX()+2;iX++) {
    double numer = hist_numer->GetBinContent(iX);
    double denom = hist_denom->GetBinContent(iX);
    double eff(0.), efferr (0.); 
    if(denom!=0) 
      getEff(numer, denom, eff, efferr);

    hist_eff->SetBinContent(iX, eff);
    hist_eff->SetBinError(iX, efferr);
  }
}


void fill2DEffHist(TH2F* hist_numer, TH2F* hist_denom, TH2F* hist_eff) {
  cout << "fill2DEffHist()" << endl;
  if(!hist_numer || !hist_denom ) return;
  
  // check the hist_numer and hist_denom have the same bin structures
  if(   hist_numer->GetNbinsX() != hist_denom->GetNbinsX() 
	|| hist_numer->GetNbinsY() != hist_denom->GetNbinsY() ) return;
  if(   hist_numer->GetXaxis()->GetXmin() != hist_denom->GetXaxis()->GetXmin() 
	|| hist_numer->GetYaxis()->GetXmin() != hist_denom->GetYaxis()->GetXmin()) return;
  if(   hist_numer->GetXaxis()->GetXmax() != hist_denom->GetXaxis()->GetXmax() 
	|| hist_numer->GetYaxis()->GetXmax() != hist_denom->GetYaxis()->GetXmax()) return;
  
  // include the under/over-flow bins
  for(int iX=0;iX<hist_numer->GetNbinsX()+2;iX++) {
    for(int iY=0;iY<hist_numer->GetNbinsY()+2;iY++) {
      double numer = hist_numer->GetBinContent(iX, iY);
      double denom = hist_denom->GetBinContent(iX, iY);
      double eff(0.), efferr (0.); 
      if(denom!=0) 
	getEff(numer, denom, eff, efferr);
      hist_eff->SetBinContent(iX, iY, eff);
      hist_eff->SetBinError(iX, iY, efferr);
    }
  }
}

 
void fillEffHist(const char* process, double weight) {
  // std::cout << "doAnalysis::fillEffHist...\n";
  // only fill efficiency for the ww/zz processes..
  if(  TString(process) != "ww" && TString(process) != "zz" ) return;
  if (!isIdentified(process)) return;
  
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
      
      for(unsigned int i_els=0;i_els<cms2.els_charge().size();++i_els) {
	if( goodElectronIsolated(i_els) && cms2.els_mc3idx().at(i_els) == int(i)) {
	  els_numer_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
	  els_numer_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
	  els_numer_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
	}
      }
    } // ==== End of Electron Efficiency
    // ==== Muon Efficiency
    if(TMath::Abs(cms2.genps_id().at(i)) == 13) {
      mus_denom_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
      mus_denom_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
      mus_denom_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
      
      for(unsigned int i_mus=0;i_mus<cms2.mus_charge().size();++i_mus) {
	if( goodMuonIsolated(i_mus) && cms2.mus_mc3idx().at(i_mus) == int(i)) {
	  mus_numer_mc_->Fill(cms2.genps_p4().at(i).eta(),cms2.genps_p4().at(i).pt());
	  mus_numer_mc_eta_->Fill(cms2.genps_p4().at(i).eta());
	  mus_numer_mc_pt_->Fill(cms2.genps_p4().at(i).pt());
	}
      }
    } // ==== End of Muon Efficiency
  }
}
void fillFOHist() {
  
  // == Fill the parton-lepton FO related histograms

  // skipping the first 6 genparticles which are the incoming partons
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    // selecting status = 3 colored objects
    if (cms2.genps_status().at(i) != 3 || (cms2.genps_id().at(i) != 21 && TMath::Abs(cms2.genps_id().at(i)) > 8 )) continue;
    // apply acceptance cuts, which must be looser than the FO definition
    if (cms2.genps_p4().at(i).pt() < 10 || TMath::Abs(cms2.genps_p4().at(i).eta()) > 2.5 ) continue;

    parton_->Fill(genps_p4().at(i).pt(), TMath::Abs(genps_p4().at(i).eta()));
    parton_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
    parton_pt_->Fill(genps_p4().at(i).pt());

    double minDRGenPs = 999.0; 
    int idx_minDRGenPs = -1;
    findClosestGenPs( cms2.genps_p4().at(i), minDRGenPs, idx_minDRGenPs); 
    dR_parton_lepton_->Fill(minDRGenPs);

    // for electron FO
    double minDREleFO = 999.0;
    int idx_minDREleFO = -1;
    findClosestEleFO( cms2.genps_p4().at(i), minDREleFO, idx_minDREleFO); 
    
    if( minDREleFO < 0.2 ) {
      els_fake_->Fill(genps_p4().at(i).pt(), TMath::Abs(genps_p4().at(i).eta()));
      els_fake_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
      els_fake_pt_->Fill(genps_p4().at(i).pt());
      els_fo_parton_->Fill(els_p4().at(idx_minDREleFO).pt(), genps_p4().at(i).pt()/els_p4().at(idx_minDREleFO).pt());
      els_fo_parton_dEta_->Fill(TMath::Abs( cms2.els_p4().at(idx_minDREleFO).eta() - cms2.genps_p4().at(i).eta() ));
      els_fo_parton_dPhi_->Fill(TMath::Abs( cms2.els_p4().at(idx_minDREleFO).phi() - cms2.genps_p4().at(i).phi() ));
      els_fo_MCType_->Fill(elFakeMCCategory(idx_minDREleFO));
      els_fo_partonID_->Fill(genps_id().at(i));
      els_fo_MCType_pt_->Fill(els_p4().at(idx_minDREleFO).pt(), elFakeMCCategory(idx_minDREleFO));
      els_fo_partonID_pt_->Fill(els_p4().at(idx_minDREleFO).pt(), genps_id().at(i));
    }

    // for muon FO
    double minDRMuFO = 999.0;
    int idx_minDRMuFO = -1;
    findClosestMuFO( cms2.genps_p4().at(i), minDRMuFO, idx_minDRMuFO); 
    if( minDRMuFO < 0.2 ) {
      mus_fake_->Fill(genps_p4().at(i).pt(), TMath::Abs(genps_p4().at(i).eta()));
      mus_fake_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
      mus_fake_pt_->Fill(genps_p4().at(i).pt());
      mus_fo_parton_->Fill(mus_p4().at(idx_minDRMuFO).pt(), genps_p4().at(i).pt()/mus_p4().at(idx_minDRMuFO).pt());
      mus_fo_parton_dEta_->Fill(TMath::Abs( cms2.mus_p4().at(idx_minDRMuFO).eta() - cms2.genps_p4().at(i).eta() ));
      mus_fo_parton_dPhi_->Fill(TMath::Abs( cms2.mus_p4().at(idx_minDRMuFO).phi() - cms2.genps_p4().at(i).phi() ));
      mus_fo_MCType_->Fill(muFakeMCCategory(idx_minDRMuFO));
      mus_fo_partonID_->Fill(genps_id().at(i));
      mus_fo_MCType_pt_->Fill(mus_p4().at(idx_minDRMuFO).pt(), muFakeMCCategory(idx_minDRMuFO));
      mus_fo_partonID_pt_->Fill(mus_p4().at(idx_minDRMuFO).pt(), genps_id().at(i));
    }
  }
    
  // == Fill the lepton FO - lepton histograms
  // electrons
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    // veto the good leptons from the W decays
    // if ( TMath::Abs(cms2.els_mc3_motherid().at(i)) == 24) continue;
    // From UserCode/JRibnik/CMS2/NtupleMacros/CORE/mcSelections.cc
    if ( leptonIsFromW(i, 11, false) > 0) continue;
    if ( fakableElectron (i) ) {
      els_fo_->Fill(cms2.els_p4().at(i).pt(), TMath::Abs(cms2.els_p4().at(i).eta()));
      if( goodElectronIsolated(i))
	els_good_->Fill(cms2.els_p4().at(i).pt(), TMath::Abs(cms2.els_p4().at(i).eta()));
    }
  }

  // muons
  for (int i=0; i < int(cms2.mus_charge().size()); ++i) {
    // veto the good leptons from the W decays
    // if ( TMath::Abs(cms2.mus_mc3_motherid().at(i)) == 24) continue;
    if ( leptonIsFromW(i, 13, false) > 0) continue;
    if ( fakableMuon(i)) {
      mus_fo_->Fill(cms2.mus_p4().at(i).pt(), TMath::Abs(cms2.mus_p4().at(i).eta()));
      if( goodMuonIsolated(i))
	mus_good_->Fill(cms2.mus_p4().at(i).pt(), TMath::Abs(cms2.mus_p4().at(i).eta()));
    }
  }
}


void saveMCUtilOutput(const char* process, TFile *utilFile_)
{
  cout << "saveMCUtilHist()" << endl;
  utilFile_->cd();
  
  if(TString(process) == "ww" || TString(process) == "zz" ) {
    fill2DEffHist(els_numer_mc_, els_denom_mc_, els_eff_mc_);
    fill1DEffHist(els_numer_mc_eta_, els_denom_mc_eta_, els_eff_mc_eta_);
    fill1DEffHist(els_numer_mc_pt_, els_denom_mc_pt_, els_eff_mc_pt_);
    
    fill2DEffHist(mus_numer_mc_, mus_denom_mc_, mus_eff_mc_);
    fill1DEffHist(mus_numer_mc_eta_, mus_denom_mc_eta_, mus_eff_mc_eta_);
    fill1DEffHist(mus_numer_mc_pt_, mus_denom_mc_pt_, mus_eff_mc_pt_);
    
    els_denom_mc_->Write();
    els_numer_mc_->Write();
    els_eff_mc_->Write();
    
    els_denom_mc_eta_->Write();
    els_numer_mc_eta_->Write();
    els_eff_mc_eta_->Write();
    
    els_denom_mc_pt_->Write();
    els_numer_mc_pt_->Write();
    els_eff_mc_pt_->Write();
    
    mus_denom_mc_->Write();
    mus_numer_mc_->Write();
    mus_eff_mc_->Write();
    
    mus_denom_mc_eta_->Write();
    mus_numer_mc_eta_->Write();
    mus_eff_mc_eta_->Write();
    
    mus_denom_mc_pt_->Write();
    mus_numer_mc_pt_->Write();
    mus_eff_mc_pt_->Write();
  }
  
  kx_->Write();
  ky_->Write();
  kt_->Write();
  boost_->Write();

  kx_1jet_->Write();
  ky_1jet_->Write();
  kt_1jet_->Write();
  boost_1jet_->Write();
  
  kx_2jet_->Write();
  ky_2jet_->Write();
  kt_2jet_->Write();
  boost_2jet_->Write();



  // Fake rate related histograms filled only for wjets MC
  if(TString(process) == "wjets") {
    fill2DEffHist(els_fake_, parton_, els_genfr_);
    fill2DEffHist(mus_fake_, parton_, mus_genfr_);
    fill1DEffHist(els_fake_eta_, parton_eta_, els_genfr_eta_);
    fill1DEffHist(mus_fake_eta_, parton_eta_, mus_genfr_eta_);
    fill1DEffHist(els_fake_pt_, parton_pt_, els_genfr_pt_);
    fill1DEffHist(mus_fake_pt_, parton_pt_, mus_genfr_pt_);
   
    fill2DEffHist(els_good_, els_fo_, els_fr_);
    fill2DEffHist(mus_good_, mus_fo_, mus_fr_);
        
    parton_->Write();
    els_fake_->Write();
    els_genfr_->Write();
    mus_fake_->Write();
    mus_genfr_->Write();

    parton_eta_->Write();
    els_fake_eta_->Write();
    els_genfr_eta_->Write();
    mus_fake_eta_->Write();
    mus_genfr_eta_->Write();
    
    parton_pt_->Write();
    els_fake_pt_->Write();
    els_genfr_pt_->Write();
    mus_fake_pt_->Write();
    mus_genfr_pt_->Write();
    
    els_fo_->Write();
    els_good_->Write();
    els_fr_->Write();

    mus_fo_->Write();
    mus_good_->Write();
    mus_fr_->Write();
    
    els_fo_parton_->Write();
    mus_fo_parton_->Write();
    
    els_fo_parton_dEta_->Write();
    els_fo_parton_dPhi_->Write();
    els_fo_MCType_->Write();
    els_fo_partonID_->Write();
    els_fo_MCType_pt_->Write();
    els_fo_partonID_pt_->Write();

    mus_fo_parton_dEta_->Write();
    mus_fo_parton_dPhi_->Write();
    mus_fo_MCType_->Write();
    mus_fo_partonID_->Write();
    mus_fo_MCType_pt_->Write();
    mus_fo_partonID_pt_->Write(); 
    
    dR_parton_lepton_->Write();

  }
}

