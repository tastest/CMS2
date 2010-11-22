// -*- C++ -*-
//
// Package:    SParmMaker
// Class:      SParmMaker
// 
/**\class SParmMaker SParmMaker.cc CMS2/NtupleMaker/src/SParmMaker.cc

   Description: copy SUSY mSUGRA parameters into the EDM event tree

   Implementation:
   - extract and fill variables
*/
//
// Original Ben Hooberman
// Created:  Wed Mar  24 12:23:38 CDT 2010
// 
//
//


// system include files
#include <memory>
#include <vector>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/JetReco/interface/CaloJet.h"

#include "CMS2/NtupleMaker/interface/SParmMaker.h"

//
// class declaration
//

//
// constructors and destructor
//
using namespace edm;
using namespace std;

SParmMaker::SParmMaker(const edm::ParameterSet& iConfig) {

  aliasprefix_ = iConfig.getUntrackedParameter<std::string>("aliasPrefix");
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  // product of this EDProducer
  produces<float>               (branchprefix+"m0"         ).setBranchAlias(aliasprefix_+"_m0"       );
  produces<float>               (branchprefix+"m12"        ).setBranchAlias(aliasprefix_+"_m12"      );
  produces<float>               (branchprefix+"A"          ).setBranchAlias(aliasprefix_+"_A"        );
  produces<float>               (branchprefix+"mu"         ).setBranchAlias(aliasprefix_+"_mu"       );
  produces<float>               (branchprefix+"tanBeta"    ).setBranchAlias(aliasprefix_+"_tanBeta"  );
  produces<float>               (branchprefix+"xsec"       ).setBranchAlias(aliasprefix_+"_xsec"     );
  
  // parameters from configuration
  sparm_m0InputTag      = iConfig.getParameter<edm::InputTag>("sparm_m0InputTag"        );
  sparm_m12InputTag     = iConfig.getParameter<edm::InputTag>("sparm_m12InputTag"       );
  sparm_AInputTag       = iConfig.getParameter<edm::InputTag>("sparm_AInputTag"         );
  sparm_muInputTag      = iConfig.getParameter<edm::InputTag>("sparm_muInputTag"        );
  sparm_tanBetaInputTag = iConfig.getParameter<edm::InputTag>("sparm_tanBetaInputTag"   );
  sparm_xsecInputTag    = iConfig.getParameter<edm::InputTag>("sparm_xsecInputTag"      );

}

SParmMaker::~SParmMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SParmMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<float>  sparm_m0       (new float);
  auto_ptr<float>  sparm_m12      (new float);
  auto_ptr<float>  sparm_A        (new float);
  auto_ptr<float>  sparm_mu       (new float);
  auto_ptr<float>  sparm_tanBeta  (new float);
  auto_ptr<float>  sparm_xsec     (new float);
  
  edm::Handle<double> sparm_m0Handle;
  edm::Handle<double> sparm_m12Handle;
  edm::Handle<double> sparm_AHandle;
  edm::Handle<double> sparm_muHandle;
  edm::Handle<double> sparm_tanBetaHandle;
  edm::Handle<double> sparm_xsecHandle;

  iEvent.getByLabel(sparm_m0InputTag,      sparm_m0Handle); 
  iEvent.getByLabel(sparm_m12InputTag,     sparm_m12Handle); 
  iEvent.getByLabel(sparm_AInputTag,       sparm_AHandle); 
  iEvent.getByLabel(sparm_muInputTag,      sparm_muHandle); 
  iEvent.getByLabel(sparm_tanBetaInputTag, sparm_tanBetaHandle); 
  iEvent.getByLabel(sparm_xsecInputTag,    sparm_xsecHandle); 

  if( sparm_m0Handle.isValid() ){
    *sparm_m0 = (float)*(sparm_m0Handle.product());
  }
  else{
    *sparm_m0 = -9999.;
  }

  if( sparm_m12Handle.isValid() ){
    *sparm_m12 = (float)*(sparm_m12Handle.product());
  }
  else{
    *sparm_m12 = -9999.;
  }

  if( sparm_AHandle.isValid() ){
    *sparm_A = (float)*(sparm_AHandle.product());
  }
  else{
    *sparm_A = -9999.;
  }

  if( sparm_muHandle.isValid() ){
    *sparm_mu = (float)*(sparm_muHandle.product());
  }
  else{
    *sparm_mu = -9999.;
  }
  
  if( sparm_tanBetaHandle.isValid() ){
    *sparm_tanBeta = (float)*(sparm_tanBetaHandle.product());
  }
  else{
    *sparm_tanBeta = -9999.;
  }

  if( sparm_xsecHandle.isValid() ){
    *sparm_xsec = (float)*(sparm_xsecHandle.product());
  }
  else{
    *sparm_xsec = -9999.;
  }

  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(sparm_m0        ,branchprefix+"m0"        );
  iEvent.put(sparm_m12       ,branchprefix+"m12"       );
  iEvent.put(sparm_A         ,branchprefix+"A"         );
  iEvent.put(sparm_mu        ,branchprefix+"mu"        );
  iEvent.put(sparm_tanBeta   ,branchprefix+"tanBeta"   );
  iEvent.put(sparm_xsec      ,branchprefix+"xsec"      );
  
}

// ------------ method called once each job just before starting event loop  ------------
void SParmMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void SParmMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SParmMaker);
