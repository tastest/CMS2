//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventMaker
// 
/**\class EventMaker EventMaker.cc CMS2/NtupleMakerMaker/src/EventMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventMaker.cc,v 1.16.2.1 2009/11/12 19:31:19 slava77 Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/TriggerNames.h"
#include "CMS2/NtupleMaker/interface/EventMaker.h"


#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"


#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "TString.h"

typedef math::XYZTLorentzVector LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

EventMaker::EventMaker(const edm::ParameterSet& iConfig) {

  
  produces<unsigned int>   ("evtrun"               ).setBranchAlias("evt_run"                  );
  produces<unsigned int>   ("evtevent"             ).setBranchAlias("evt_event"                );
  produces<unsigned int>   ("evtlumiBlock"         ).setBranchAlias("evt_lumiBlock"            );
  produces<TString>        ("evtdataset"           ).setBranchAlias("evt_dataset"              );
  produces<float>  ("evtbField"            ).setBranchAlias("evt_bField"               );
  produces<float>  ("evtweight"            ).setBranchAlias("evt_weight"               );
  produces<float>  ("evtxsecincl"          ).setBranchAlias("evt_xsec_incl"            );
  produces<float>  ("evtxsecexcl"          ).setBranchAlias("evt_xsec_excl"            );
  produces<float>  ("evtkfactor"           ).setBranchAlias("evt_kfactor"              );
  
  inclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("inclusiveCrossSection");
  exclusiveCrossSectionValue = iConfig.getUntrackedParameter<double>("exclusiveCrossSection");
  kfactorValue = iConfig.getUntrackedParameter<double>("kfactor");
  datasetName_ = iConfig.getParameter<std::string>("datasetName");
  
}


EventMaker::~EventMaker() {}

void  EventMaker::beginJob(const edm::EventSetup&) {
      
}

void EventMaker::endJob() {
  
}

// ------------ method called to produce the data  ------------
void EventMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  auto_ptr<unsigned int>      evt_run               (new unsigned int);
  auto_ptr<unsigned int>      evt_event             (new unsigned int);
  auto_ptr<unsigned int>      evt_lumiBlock         (new unsigned int);
  auto_ptr<TString>           evt_dataset           (new TString(datasetName_.c_str()));
  auto_ptr<float>    evt_bField            (new float);
  auto_ptr<float>    evt_weight            (new float);
  auto_ptr<float>    evt_xsec_incl         (new float);
  auto_ptr<float>    evt_xsec_excl         (new float);
  auto_ptr<float>    evt_kfactor           (new float);
  //auto_ptr<string>   evt_temp              (new string);
  
  *evt_run   = iEvent.id().run();
  *evt_event = iEvent.id().event();
  *evt_lumiBlock = iEvent.luminosityBlock();
  
 //need the magnetic field
  ESHandle<MagneticField> magneticField;
  iSetup.get<IdealMagneticFieldRecord>().get(magneticField);
  *evt_bField = magneticField->inTesla(GlobalPoint(0.,0.,0.)).z();
 
  
   //get the MC event weights
   //if weights do not exist (Pythia), default is weight of 1
  vector< Handle<HepMCProduct> > hepmc_vect;
  iEvent.getManyByType(hepmc_vect);
  HepMC::WeightContainer wc;
  if(hepmc_vect.size() != 0) { //found HepMC branch
    const HepMC::GenEvent *genEvt = hepmc_vect.at(0)->GetEvent();
     wc = genEvt->weights();
     float weight = -999.;
     if(wc.size() > 0 ) {
	weight = (float)wc[0];
	} 
     if(wc.size() == 0) weight = -999.;
     *evt_weight = weight;
  } else {
    try {
      Handle<double> evtwt;
      iEvent.getByLabel("genEventWeight", evtwt);
      *evt_weight = (float)*evtwt;
    } catch (edm::Exception const& x) {
      *evt_weight = 1.;
    }
  }   

  *evt_xsec_incl = inclusiveCrossSectionValue;
  *evt_xsec_excl = exclusiveCrossSectionValue;
  *evt_kfactor   = kfactorValue;
      
  iEvent.put(evt_run              ,"evtrun"             );
  iEvent.put(evt_event            ,"evtevent"           );
  iEvent.put(evt_lumiBlock        ,"evtlumiBlock"       );
  iEvent.put(evt_dataset          ,"evtdataset"         );
  iEvent.put(evt_bField           ,"evtbField"          );
  iEvent.put(evt_weight           ,"evtweight"          );
  iEvent.put(evt_xsec_incl        ,"evtxsecincl"        );
  iEvent.put(evt_xsec_excl        ,"evtxsecexcl"        );
  iEvent.put(evt_kfactor          ,"evtkfactor"         );
}


//define this as a plug-in
DEFINE_FWK_MODULE(EventMaker);
