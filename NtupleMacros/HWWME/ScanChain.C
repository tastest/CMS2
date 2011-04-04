/* 
=== Description

  This is the looper for the H->WW analysis using the ME method, written by Sergo Yanyan and Kevin.
  The main goals of this looper are
  1. producing "Util.root" containing a bunch of histograms needed for ME calculation
  2. Producing Smurf ntuples 

=== Usage:
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
using namespace tas;

// Smurf
const char* config_info = "Smurf HWW V1 selection; Winter10 FlatPU samples; 35.5/pb Run2010A (Sep17ReReco) & Run2010B (PromptReco)";
#include "../../../Smurf/Core/SmurfTree.h"


// JEC 
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

  if( nEvents == -1 ) nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  unsigned int nEventsTotal = 0;
  
  if(!realData) InitMCUtilHist(process, utilFile_);
  // make smurf ntuples
  SmurfTree smurfTree;
  smurfTree.CreateTree();
  smurfTree.tree_->SetDirectory(0);


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
      
      // identifyEvents by MC truth if specified
      if ( (!realData && identifyEvents) && !isIdentified (process)) continue;
      
      // Get fake-rate related histograms
      if(!realData && TString(process) == "wjets") fillFOHist();
      
      // skip events without hypothesis
      if (cms2.hyp_type().size() == 0) continue; 

      // duplicate event removal
      EventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.evt_lumiBlock(), cms2.trks_d0()[0], cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
      if (is_duplicate(id)) continue;

      // Calculate event weight for MC
      double mcweight;
      if ((process!="data")){
	mcweight = cms2.genps_weight() > 0.0 ? 1.0 : -1.0;
	weight = IntLumi * mcweight * (Xsect>0?Xsect:cms2.evt_xsec_excl()*cms2.evt_kfactor()) /
	  (nProcessedEvents>0?nProcessedEvents:cms2.evt_nEvts());
      } else weight =1;
      
      // fill lepton efficiency and system boost histograms which are needed for ME calculations
      if(!realData) 	fillEffHist(process, weight);
      
      // Start Event selection based on the hypothesis
      unsigned int nHyps = cms2.hyp_type().size();

      for( unsigned int i_hyp = 0; i_hyp < nHyps; ++i_hyp ) {
	if(cms2.hyp_p4().at(i_hyp).mass2() < 0 ) break;
	_cutWord->SetAllBitsFalse();
	ApplyEventSelection(i_hyp, realData);
	
	int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
	bool accept = true;
	for (int j = 0; j < kNCuts; j++) {
	  if (_cutMask->IsTrue(j)) accept &= _cutWord->IsTrue (j);
        }

        if (accept){
	  if(!realData) fillKtHist(process, weight);
	  FillSmurfNtuple(smurfTree,i_hyp,weight,process);
	  smurfTree.tree_->Fill();

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
    
  TFile* fSmurf = TFile::Open(Form("%s.root",process),"RECREATE");
  assert(fSmurf);
  smurfTree.tree_->Write();
  smurfTree.info_.SetTitle(config_info);
  smurfTree.info_.Write();
  fSmurf->Close();
  
  if(!realData) saveMCUtilOutput(process, utilFile_);
  
  cout<<"Total Events Before Selection "<<chain->GetEntries()<<"\n";
}


int ApplyEventSelection( unsigned int i_hyp, bool realData){

  if ( std::max(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<20 ) return false;
  if ( std::min(cms2.hyp_lt_p4().at(i_hyp).pt(),cms2.hyp_ll_p4().at(i_hyp).pt())<10 ) return false;

  int type =  getHypothesisType(cms2.hyp_type()[i_hyp]);
  
  // Apply trigger requirements in data
  _cutWord->SetTrue(kcut_Trigger);
  if (realData && ! passedTriggerRequirements())   _cutWord->SetFalse(kcut_Trigger);

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

void FillSmurfNtuple(SmurfTree& tree, unsigned int i_hyp, double weight, const char* process) {
  tree.InitVariables();
  tree.run_   = cms2.evt_run();
  tree.event_ = cms2.evt_event();
  tree.lumi_  = cms2.evt_lumiBlock();
  tree.nvtx_  = nGoodVertex();
  tree.scale1fb_ = weight;
  tree.met_    = metValue();
  tree.metPhi_ = metPhiValue();
  // sumet_;
  bool ltIsFirst = true;
  if ( cms2.hyp_lt_p4().at(i_hyp).pt()<cms2.hyp_ll_p4().at(i_hyp).pt() ) ltIsFirst = false;
  tree.type_ = SmurfTree::Type(cms2.hyp_type().at(i_hyp));
  if ( tree.type_ == SmurfTree::em || tree.type_ == SmurfTree::me ){
    if ( ltIsFirst )
      tree.type_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? SmurfTree::em : SmurfTree::me;
    else
      tree.type_ = abs(cms2.hyp_lt_id().at(i_hyp))==11 ? SmurfTree::me : SmurfTree::em;
  }
  tree.lep1_ = ltIsFirst ? cms2.hyp_lt_p4().at(i_hyp) : cms2.hyp_ll_p4().at(i_hyp);
  tree.lep2_ = ltIsFirst ? cms2.hyp_ll_p4().at(i_hyp) : cms2.hyp_lt_p4().at(i_hyp);
  tree.lq1_   = ltIsFirst ? cms2.hyp_lt_charge().at(i_hyp) : cms2.hyp_ll_charge().at(i_hyp);
  tree.lq2_   = ltIsFirst ? cms2.hyp_ll_charge().at(i_hyp) : cms2.hyp_lt_charge().at(i_hyp);
  tree.lid1_  = ltIsFirst ? cms2.hyp_lt_id().at(i_hyp) : cms2.hyp_ll_id().at(i_hyp);
  tree.lid2_  = ltIsFirst ? cms2.hyp_ll_id().at(i_hyp) : cms2.hyp_lt_id().at(i_hyp);
  const std::vector<LorentzVector>& jets = getJets(pfJet, i_hyp, 25, 5.0, true, false);
  if (jets.size()>0) tree.jet1_ = jets.at(0);
  if (jets.size()>1) tree.jet2_ = jets.at(1);
  // jet1_btag_;
  // jet2_btag_;
  tree.njets_ = jets.size();
  tree.evtype_ = SmurfTree::ZeroJet;
  tree.dilep_ = cms2.hyp_p4().at(i_hyp);
  tree.pmet_ = projectedMet(i_hyp);
  tree.dPhi_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4().at(i_hyp),cms2.hyp_ll_p4().at(i_hyp)));
  tree.dR_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4().at(i_hyp),cms2.hyp_ll_p4().at(i_hyp));
  if (jets.size()>0) {
    tree.dPhiLep1Jet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_lt_p4().at(i_hyp),jets.at(0)));
    tree.dRLep1Jet1_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4().at(i_hyp),jets.at(0));
    tree.dPhiLep2Jet1_ = fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_ll_p4().at(i_hyp),jets.at(0)));
    tree.dRLep2Jet1_   = ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4().at(i_hyp),jets.at(0));
    tree.dPhiDiLepJet1_= fabs(ROOT::Math::VectorUtil::DeltaPhi(cms2.hyp_p4().at(i_hyp),jets.at(0)));
    if (!ltIsFirst){
      std::swap(tree.dPhiLep1Jet1_,tree.dPhiLep2Jet1_);
      std::swap(tree.dRLep1Jet1_,tree.dRLep2Jet1_);
    }
  }
  tree.dPhiDiLepMET_ = acos(cos(cms2.hyp_p4().at(i_hyp).phi()-metPhiValue()));
  tree.dPhiLep1MET_ = acos(cos(tree.lep1_.phi()-metPhiValue()));
  tree.dPhiLep2MET_ = acos(cos(tree.lep2_.phi()-metPhiValue()));

  tree.mt_ = mt(tree.dilep_.pt(),tree.met_,tree.dPhiDiLepMET_);
  tree.mt1_ = mt(tree.lep1_.pt(),tree.met_,tree.dPhiLep1MET_);
  tree.mt2_ = mt(tree.lep2_.pt(),tree.met_,tree.dPhiLep2MET_);

  if (process!="data"){
    tree.genmet_ = cms2.gen_met();
    tree.genmetPhi_ = cms2.gen_metPhi();
    tree.lep1McId_ = ltIsFirst ? cms2.hyp_lt_mc_id().at(i_hyp) : cms2.hyp_ll_mc_id().at(i_hyp);
    tree.lep2McId_ = ltIsFirst ? cms2.hyp_ll_mc_id().at(i_hyp) : cms2.hyp_lt_mc_id().at(i_hyp);
  }

  tree.dstype_    = SmurfTree::other;
  if(TString(process) == "ww")  tree.dstype_ = SmurfTree::qqww;
  if(TString(process) == "wz")  tree.dstype_ = SmurfTree::wz;
  if(TString(process) == "zz")  tree.dstype_ = SmurfTree::zz;
  if(TString(process) == "wjets")  tree.dstype_ = SmurfTree::wjets;
  if(TString(process) == "dyee")  tree.dstype_ = SmurfTree::dyee;
  if(TString(process) == "dymm")  tree.dstype_ = SmurfTree::dymm;
  if(TString(process) == "dytt")  tree.dstype_ = SmurfTree::dytt;
  if(TString(process) == "ttbar")  tree.dstype_ = SmurfTree::ttbar;
  if(TString(process) == "tw")  tree.dstype_ = SmurfTree::tw;
  if(TString(process) == "qcd")  tree.dstype_ = SmurfTree::qcd;
  if(TString(process) == "data")  tree.dstype_ = SmurfTree::data;
  if(TString(process) == "hww120")  tree.dstype_ = SmurfTree::hww120;
  if(TString(process) == "hww130")  tree.dstype_ = SmurfTree::hww130;
  if(TString(process) == "hww140")  tree.dstype_ = SmurfTree::hww140;
  if(TString(process) == "hww150")  tree.dstype_ = SmurfTree::hww150;
  if(TString(process) == "hww160")  tree.dstype_ = SmurfTree::hww160;
  if(TString(process) == "hww170")  tree.dstype_ = SmurfTree::hww170;
  if(TString(process) == "hww180")  tree.dstype_ = SmurfTree::hww180;
  if(TString(process) == "hww190")  tree.dstype_ = SmurfTree::hww190;
  if(TString(process) == "hww200")  tree.dstype_ = SmurfTree::hww200;
  if(TString(process) == "hww210")  tree.dstype_ = SmurfTree::hww210;
  if(TString(process) == "hww220")  tree.dstype_ = SmurfTree::hww220;
  if(TString(process) == "hww230")  tree.dstype_ = SmurfTree::hww230;
  if(TString(process) == "hww250")  tree.dstype_ = SmurfTree::hww250;
  
}




// == Utility fucntions

double mt(double pt1, double pt2, double dphi){
  return 2*sqrt(pt1*pt2)*fabs(sin(dphi/2));
}

bool inZmassWindow(float mass){
  // return ( mass > 76. && mass < 106. );  
  return fabs(mass - 91.1876) < 15;
}

double projectedMet(int i_hyp)
{
  double DeltaPhi = nearestDeltaPhi(metPhiValue(),i_hyp);
  if (DeltaPhi < TMath::Pi()/2) return metValue()*TMath::Sin(DeltaPhi);
  return metValue();
}


double metValue(){    return cms2.evt_tcmet(); }
double metPhiValue(){ return cms2.evt_tcmetPhi(); }



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
	jec = jetCorrection(cms2.jpts_p4().at(i), jet_corrector_jpt);
      if ( cms2.jpts_p4().at(i).pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.jpts_p4().at(i).eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jpts_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jpts_p4().at(i))) < vetoCone ) continue;
      jets.push_back(cms2.jpts_p4().at(i) * jec);
    }
    break;
  case pfJet:
    for ( unsigned int i=0; i < cms2.pfjets_p4().size(); ++i) {
      if(applyJEC)
	jec = jetCorrection(cms2.pfjets_p4().at(i), jet_corrector_pf);
      if ( cms2.pfjets_p4().at(i).pt() * jec < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.pfjets_p4().at(i).eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.pfjets_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.pfjets_p4().at(i))) < vetoCone ) continue;
      jets.push_back(cms2.pfjets_p4().at(i) * jec);
    }
    break;
  case GenJet:
    for ( unsigned int i=0; i < cms2.genjets_p4().size(); ++i) {
      if ( cms2.genjets_p4().at(i).pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.genjets_p4().at(i).eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.genjets_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.genjets_p4().at(i))) < vetoCone ) continue;
      jets.push_back(cms2.genjets_p4().at(i));
    }
    break;
  case CaloJet:
    for ( unsigned int i=0; i < cms2.jets_pat_jet_p4().size(); ++i) {
      if ( cms2.jets_pat_jet_p4().at(i).pt() < etThreshold ) continue; // note that this is already corrected
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.jets_pat_jet_p4().at(i).eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.jets_pat_jet_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.jets_pat_jet_p4().at(i))) < vetoCone ) continue;
      jets.push_back(cms2.jets_pat_jet_p4().at(i));
    }
    break;
  case TrkJet:
    for ( unsigned int i=0; i < cms2.trkjets_p4().size(); ++i) {
      if(applyJEC)
	jec = jetCorrection(cms2.trkjets_p4().at(i), jet_corrector_trk);
      if ( cms2.trkjets_p4().at(i).pt() < etThreshold ) continue;
      if ( btag && !defaultBTag(type,i) ) continue;
      if ( TMath::Abs(cms2.trkjets_p4().at(i).eta()) > etaMax ) continue;
      if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4().at(i))) < vetoCone ||
	   TMath::Abs(ROOT::Math::VectorUtil::DeltaR(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4().at(i))) < vetoCone ) continue;
      jets.push_back(cms2.trkjets_p4().at(i) * jec);
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
    if (!isGoodVertex(i)) continue;
    nVtx++;
  }
  return nVtx;
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

bool fakableMuon(unsigned int i){
  bool ptcut = cms2.mus_p4().at(i).pt() >= 10.0;
  return ptcut && muonId(i, muonSelectionFO_mu_wwV1_iso10);
}


bool ww_mud0PV(unsigned int index){
  if ( cms2.vtxs_sumpt().empty() ) return false;
  double sumPtMax = -1;
  int iMax = -1;
  for ( unsigned int i = 0; i < cms2.vtxs_sumpt().size(); ++i ){
    if (cms2.vtxs_isFake().at(i)) continue;
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
  bool ptcut = cms2.els_p4().at(i).pt() >= 15.0;
  // return  ptcut && pass_electronSelection( i, electronSelection_wwV1);
  return  ptcut && pass_electronSelection( i, electronSelection_smurfV2);
}

bool fakableElectron(unsigned int i){
  bool ptcut = cms2.els_p4().at(i).pt() >= 15.0;
  // extrapolate in partial id, iso and d0
  return ptcut && pass_electronSelection( i, electronSelectionFO_el_wwV1_v2);
}

bool ww_elId(unsigned int index){
  if( fabs(cms2.els_conv_dist().at(index)) < 0.02 &&
      fabs(cms2.els_conv_dcot().at(index)) < 0.02) return false;
  if (! (electronId_VBTF(index, VBTF_35X_80) & (1<<ELEID_ID)) ) return false;
  if ( cms2.els_exp_innerlayers39X().at(index) > 0 ) return false;
  if ( !electronId_smurf_v2(index) ) return false;
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

void fillEffHist(const char* process, double weight) {
  if(TString(process) != "ww" || !isIdentified(process)) return;
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
    } // ==== End of Electron Efficiency
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
    } // ==== End of Muon Efficiency
  }
}

void findClosestEleFO(LorentzVector v_parton, double& minDR, int& idx_minDR) {
  minDR = 999.0;
  idx_minDR = -999;
  
  for (int i=0; i < int(cms2.els_charge().size()); ++i) {
    // if ( TMath::Abs(cms2.els_p4().at(i).eta()) > 2.0) continue; 
    if (cms2.els_p4().at(i).pt() < 10) continue; 
    if (goodElectronIsolated(i)) continue;
    if (! pass_electronSelection(i, electronSelectionFO_el_wwV1_v2))  continue;
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
    // if ( TMath::Abs(cms2.mus_p4().at(i).eta()) > 2.0) continue; 
    if (cms2.mus_p4().at(i).pt() < 10 ) continue; // just to av<oid the warning messages in the muonID code for low pt muons
    if (goodMuonIsolated(i)) continue;
    if (!muonId(i, muonSelectionFO_mu_wwV1_iso10)) continue;
    if ( TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton, cms2.mus_p4().at(i))) < minDR ) {
      minDR = TMath::Abs(ROOT::Math::VectorUtil::DeltaR(v_parton,cms2.mus_p4().at(i)));
      idx_minDR = i;
    }
  }
}

  

void fillKtHist(const char* process, double weight) {
  LorentzVector systP4(0.,0.,0.,0.);
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    if (TString(process).Contains("ww",TString::kExact)) {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24) 	systP4 += cms2.genps_p4().at(i);
    }
    else if (TString(process) == "wz" || TString(process) == "zz") {
      if(TMath::Abs(cms2.genps_id().at(i)) == 24 || cms2.genps_id().at(i) == 23 ) systP4 += cms2.genps_p4().at(i);
    }
  }
  kx_->Fill(systP4.Px(), weight);
  ky_->Fill(systP4.Py(), weight);
}

void fillFOHist() {
  // skipping the first 6 genparticles which are the incoming partons
  for (unsigned int i = 6;  i < cms2.genps_id().size(); ++i) {
    // selecting status = 3 colored objects
    if (cms2.genps_status().at(i) != 3 || (cms2.genps_id().at(i) != 21 && TMath::Abs(cms2.genps_id().at(i)) > 8 )) continue;
    // apply acceptance cuts, which must be looser than the FO definition
    if (cms2.genps_p4().at(i).pt() < 10 ) continue;

    parton_->Fill(TMath::Abs(genps_p4().at(i).eta()), genps_p4().at(i).pt());
    parton_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
    parton_pt_->Fill(genps_p4().at(i).pt());

    double minDREleFO = 999.0;
    int idx_minDREleFO = -1;
    findClosestEleFO( cms2.genps_p4().at(i), minDREleFO, idx_minDREleFO); 
    
    if( minDREleFO < 0.2 ) {
      els_fake_->Fill(TMath::Abs(genps_p4().at(i).eta()), genps_p4().at(i).pt());
      els_fake_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
      els_fake_pt_->Fill(genps_p4().at(i).pt());
      els_fo_parton_->Fill(els_p4().at(idx_minDREleFO).pt(), genps_p4().at(i).pt()/els_p4().at(idx_minDREleFO).pt());
      els_fo_parton_eta_->Fill(genps_p4().at(i).eta(), els_p4().at(idx_minDREleFO).eta());
      els_fo_parton_pt_->Fill(genps_p4().at(i).pt(), els_p4().at(idx_minDREleFO).pt());
      els_fo_MCType_->Fill(elFakeMCCategory(idx_minDREleFO));
      els_fo_partonID_->Fill(genps_id().at(i));
    }
    double minDRMuFO = 999.0;
    int idx_minDRMuFO = -1;
    findClosestMuFO( cms2.genps_p4().at(i), minDRMuFO, idx_minDRMuFO); 
    if( minDRMuFO < 0.2 ) {
      mus_fake_->Fill(TMath::Abs(genps_p4().at(i).eta()), genps_p4().at(i).pt());
      mus_fake_eta_->Fill(TMath::Abs(genps_p4().at(i).eta()));
      mus_fake_pt_->Fill(genps_p4().at(i).pt());
      mus_fo_parton_->Fill(mus_p4().at(idx_minDRMuFO).pt(), genps_p4().at(i).pt()/mus_p4().at(idx_minDRMuFO).pt());
      mus_fo_parton_eta_->Fill(genps_p4().at(i).eta(), mus_p4().at(idx_minDRMuFO).eta());
      mus_fo_parton_pt_->Fill(genps_p4().at(i).pt(), mus_p4().at(idx_minDRMuFO).pt());
      mus_fo_MCType_->Fill(muFakeMCCategory(idx_minDRMuFO));
      mus_fo_partonID_->Fill(genps_id().at(i));
    }
    
  }

}

bool isIdentified(const char* process) {
  // std::cout << "ScanChain::isIdentified()"<<endl;
  if(TString(process).Contains("ww", TString::kExact)) 
    return getVVType()==0;
  else  if (TString(process) == "wz")       return getVVType()==1;
  else  if (TString(process) == "zz")       return getVVType()==2;
  else  return true;
}


void ProcessSample(const char* process, std::vector<std::string> file_patterns, TFile *utilFile_,  int nEvents = -1, double IntLumi=100, double Xsect=1.0, int nProcessedEvents=-1, std::string skimFilePrefix="", bool realData=false, bool identifyEvents=false){

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

