// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      EventTrigMaker
// 
/**\class NtupleMaker NtupleMaker.cc CMS2/NtupleMaker/src/EventTrigMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  pts/4
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: EventTrigMaker.h,v 1.1.2.1 2009/11/12 19:31:19 slava77 Exp $
//
//
#ifndef NTUPLEMAKER_EVENTTRIGMAKER_H
#define NTUPLEMAKER_EVENTTRIGMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "TString.h"
//
// class decleration
//

class EventTrigMaker : public edm::EDProducer {
public:
     explicit EventTrigMaker (const edm::ParameterSet&);
      ~EventTrigMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

      // ----------member data ---------------------------
  
  
  void fillHLTInfo(const edm::Event&, 
		   int*, int*,int*, int*,
		   int*, int*,int*, int*, 
		   std::vector<TString>&);
  void fillL1Info(const edm::Event&, int*, 
		  int*, int*, int*, std::vector<TString>&,
		  const L1GtTriggerMenu* menu);
  
  bool haveTriggerInfo_;
};


#endif
