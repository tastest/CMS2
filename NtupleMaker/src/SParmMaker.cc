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
  produces<float>               (branchprefix+"xsec"       ).setBranchAlias(aliasprefix_+"_xsec"     );
  produces<float>               (branchprefix+"m1"       	).setBranchAlias(aliasprefix_+"_m1"     	);
  produces<float>               (branchprefix+"m2"       	).setBranchAlias(aliasprefix_+"_m2"     	);
  produces<float>               (branchprefix+"m3"       	).setBranchAlias(aliasprefix_+"_m3"     	);
  produces<float>               (branchprefix+"m4"       	).setBranchAlias(aliasprefix_+"_m4"     	);
 
 	
  // parameters from configuration
  sparm_xsecInputTag    = iConfig.getParameter<edm::InputTag>("sparm_xsecInputTag"      );
  sparm_m1InputTag      = iConfig.getParameter<edm::InputTag>("sparm_m1InputTag"        );
  sparm_m2InputTag      = iConfig.getParameter<edm::InputTag>("sparm_m2InputTag"        );
  sparm_m3InputTag      = iConfig.getParameter<edm::InputTag>("sparm_m3InputTag"        );
  sparm_m4InputTag      = iConfig.getParameter<edm::InputTag>("sparm_m4InputTag"        );


}

SParmMaker::~SParmMaker() {}

//
// member functions
//

// ------------ method called to produce the data  ------------
void SParmMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<float>  sparm_xsec     (new float);
  auto_ptr<float>  sparm_m1      (new float);
  auto_ptr<float>  sparm_m2      (new float);
  auto_ptr<float>  sparm_m3      (new float);
  auto_ptr<float>  sparm_m4      (new float);
  
  edm::Handle<double> sparm_xsecHandle;
  edm::Handle<double> sparm_m1Handle;
  edm::Handle<double> sparm_m2Handle;
  edm::Handle<double> sparm_m3Handle;
  edm::Handle<double> sparm_m4Handle;

  iEvent.getByLabel(sparm_xsecInputTag,    sparm_xsecHandle); 
  iEvent.getByLabel(sparm_m1InputTag,      sparm_m1Handle); 
  iEvent.getByLabel(sparm_m2InputTag,      sparm_m2Handle); 
  iEvent.getByLabel(sparm_m3InputTag,      sparm_m3Handle); 
  iEvent.getByLabel(sparm_m4InputTag,      sparm_m4Handle); 



  if( sparm_xsecHandle.isValid() ){
    *sparm_xsec = (float)*(sparm_xsecHandle.product());
  }
  else{
    *sparm_xsec = -9999.;
  }
 
  if( sparm_m1Handle.isValid() ){
    *sparm_m1 = (float)*(sparm_m1Handle.product());
  }
  else{
    *sparm_m1 = -9999.;
  }

  if( sparm_m2Handle.isValid() ){
    *sparm_m2 = (float)*(sparm_m2Handle.product());
  }
  else{
    *sparm_m2 = -9999.;
  }

  if( sparm_m3Handle.isValid() ){
    *sparm_m3 = (float)*(sparm_m3Handle.product());
  }
  else{
    *sparm_m3 = -9999.;
  }

  if( sparm_m4Handle.isValid() ){
    *sparm_m4 = (float)*(sparm_m4Handle.product());
  }
  else{
    *sparm_m4 = -9999.;
  }

  // put containers into event
  std::string branchprefix = aliasprefix_;
  if(branchprefix.find("_") != std::string::npos) branchprefix.replace(branchprefix.find("_"),1,"");

  iEvent.put(sparm_xsec      ,branchprefix+"xsec"      );
  iEvent.put(sparm_m1        ,branchprefix+"m1"        );
  iEvent.put(sparm_m2        ,branchprefix+"m2"        );
  iEvent.put(sparm_m3        ,branchprefix+"m3"        );
  iEvent.put(sparm_m4        ,branchprefix+"m4"        );
  
}

// ------------ method called once each job just before starting event loop  ------------
void SParmMaker::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void SParmMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(SParmMaker);
