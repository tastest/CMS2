// -*- C++ -*-
//
// Package:    NtupleMaker
// Class:      L1DigiMaker
// 
/**\class L1DigiMaker.cc CMS2/NtupleMaker/src/L1DigiMaker.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: L1DigiMaker.h,v 1.1.6.1 2009/11/12 19:29:53 slava77 Exp $
//
//
#ifndef NTUPLEMAKER_L1DIGIMAKER_H
#define NTUPLEMAKER_L1DIGIMAKER_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//
// class decleration
//

class L1DigiMaker : public edm::EDProducer {
public:
     explicit L1DigiMaker (const edm::ParameterSet&);
      ~L1DigiMaker();

private:
     virtual void beginJob(const edm::EventSetup&) ;
     virtual void produce(edm::Event&, const edm::EventSetup&);
     virtual void endJob() ;

  
      // ----------member data ---------------------------
  std::string l1extraModName_;
};


#endif
