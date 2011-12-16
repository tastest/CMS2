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
  produces<float>               (branchprefix+"mphi"         ).setBranchAlias(aliasprefix_+"_mphi"       );
  produces<float>               (branchprefix+"mxhi"         ).setBranchAlias(aliasprefix_+"_mxhi"       );
  produces<float>               (branchprefix+"xsec"         ).setBranchAlias(aliasprefix_+"_xsec"     );
 
 	
  // parameters from configuration
  sparm_mphiInputTag      = iConfig.getParameter<edm::InputTag>("sparm_mphiInputTag"        );
  sparm_mxhiInputTag      = iConfig.getParameter<edm::InputTag>("sparm_mxhiInputTag"        );
  sparm_xsecInputTag   	  = iConfig.getParameter<edm::InputTag>("sparm_xsecInputTag"      );


}

SParmMaker::~SParmMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SParmMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<float>  sparm_mphi       (new float);
  auto_ptr<float>  sparm_mxhi       (new float);
  auto_ptr<float>  sparm_xsec     (new float);
  
  edm::Handle<double> sparm_mphiHandle;
  edm::Handle<double> sparm_mxhiHandle;
  edm::Handle<double> sparm_xsecHandle;

  iEvent.getByLabel(sparm_mphiInputTag,      sparm_mphiHandle); 
  iEvent.getByLabel(sparm_mxhiInputTag,      sparm_mxhiHandle); 
  iEvent.getByLabel(sparm_xsecInputTag,    sparm_xsecHandle); 



  if( sparm_mphiHandle.isValid() ){
    *sparm_mphi = (float)*(sparm_mphiHandle.product());
  }
  else{
    *sparm_mphi = -9999.;
  }

  if( sparm_mxhiHandle.isValid() ){
    *sparm_mxhi = (float)*(sparm_mxhiHandle.product());
  }
  else{
    *sparm_mxhi = -9999.;
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

  iEvent.put(sparm_mphi        ,branchprefix+"mphi"        );
  iEvent.put(sparm_mxhi        ,branchprefix+"mxhi"        );
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
