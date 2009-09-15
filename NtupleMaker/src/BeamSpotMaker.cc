//-*- C++ -*-
//
// Package:    NtupleMaker
// Class:      BeamSpotMaker
// 
/**\class BeamSpotMaker BeamSpotMaker.cc CMS2/NtupleMakerMaker/src/BeamSpotMaker.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Puneeth Kalavase
//         Created:  Fri Jun  6 11:07:38 CDT 2008
// $Id: BeamSpotMaker.cc,v 1.7 2009/09/15 14:16:22 fgolf Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CMS2/NtupleMaker/interface/BeamSpotMaker.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

typedef math::XYZTLorentzVectorF LorentzVector;
typedef math::XYZPoint Point;
using namespace reco;
using namespace edm;
using namespace std;

//
// constructors and destructor
//

BeamSpotMaker::BeamSpotMaker(const edm::ParameterSet& iConfig) {

  //p4 because we're not able to (yet) read XYZPointDs in bare root for some reason 
  //the 4th co-ordinate is 0
  produces<LorentzVector>   ("evtbsp4"         ).setBranchAlias("evt_bsp4"            );
  produces<int>             ("evtbsType"       ).setBranchAlias("evt_bsType"          );
  produces<float>           ("evtbsxErr"       ).setBranchAlias("evt_bs_xErr"         );
  produces<float>           ("evtbsyErr"       ).setBranchAlias("evt_bs_yErr"         );
  produces<float>           ("evtbszErr"       ).setBranchAlias("evt_bs_zErr"         );
  produces<float>           ("evtbssigmaZ"     ).setBranchAlias("evt_bs_sigmaZ"       );
  produces<float>           ("evtbssigmaZErr"  ).setBranchAlias("evt_bs_sigmaZErr"    );
  produces<float>           ("evtbsdxdz"       ).setBranchAlias("evt_bs_dxdz"         );
  produces<float>           ("evtbsdxdzErr"    ).setBranchAlias("evt_bs_dxdzErr"      );
  produces<float>           ("evtbsdydz"       ).setBranchAlias("evt_bs_dydz"         );
  produces<float>           ("evtbsdydzErr"    ).setBranchAlias("evt_bs_dydzErr"      );
  produces<float>           ("evtbsXwidth"     ).setBranchAlias("evt_bs_Xwidth"       );
  produces<float>           ("evtbsYwidth"     ).setBranchAlias("evt_bs_Ywidth"       );
  produces<float>           ("evtbsXwidthErr"  ).setBranchAlias("evt_bs_XwidthErr"    );
  produces<float>           ("evtbsYwidthErr"  ).setBranchAlias("evt_bs_YwidthErr"    );
  produces<vector <float> > ("evtcovMatrix"    ).setBranchAlias("evt_covMatrix"       ); 
  
  beamSpotInputTag = iConfig.getParameter<InputTag>("beamSpotInputTag");
  
}


BeamSpotMaker::~BeamSpotMaker() {}

void  BeamSpotMaker::beginJob(const edm::EventSetup&) {
}

void BeamSpotMaker::endJob() {
}


// ------------ method called to produce the data  ------------
void BeamSpotMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  auto_ptr<LorentzVector>  evt_bs_p4          (new LorentzVector);
  auto_ptr<int>            evt_bsType         (new int          );
  auto_ptr<float>          evt_bs_xErr        (new float        );
  auto_ptr<float>          evt_bs_yErr        (new float        );
  auto_ptr<float>          evt_bs_zErr        (new float        );
  auto_ptr<float>          evt_bs_sigmaZ      (new float        );
  auto_ptr<float>          evt_bs_sigmaZErr   (new float        );
  auto_ptr<float>          evt_bs_dxdz        (new float        );
  auto_ptr<float>          evt_bs_dxdzErr     (new float        );
  auto_ptr<float>          evt_bs_dydz        (new float        );
  auto_ptr<float>          evt_bs_dydzErr     (new float        );
  auto_ptr<float>          evt_bs_Xwidth      (new float        );
  auto_ptr<float>          evt_bs_Ywidth      (new float        );
  auto_ptr<float>          evt_bs_XwidthErr   (new float        );
  auto_ptr<float>          evt_bs_YwidthErr   (new float        );
  auto_ptr<vector<float> > evt_bs_covMatrix   (new vector<float>);
  
  Handle<BeamSpot> beamSpotH;
  iEvent.getByLabel(beamSpotInputTag, beamSpotH);

  bool haveBeamSpot = true;
  if(!beamSpotH.isValid() )
    haveBeamSpot = false;
  
  *evt_bs_p4         = haveBeamSpot ? LorentzVector(beamSpotH->position().x(),
						    beamSpotH->position().y(),
						    beamSpotH->position().z(),
						    0) : LorentzVector(0,0,0,0);

  *evt_bsType        = haveBeamSpot ? beamSpotH->type()           : -999;
  *evt_bs_xErr       = haveBeamSpot ? beamSpotH->x0Error()        : 0.0;
  *evt_bs_yErr       = haveBeamSpot ? beamSpotH->y0Error()        : 0.0;
  *evt_bs_zErr       = haveBeamSpot ? beamSpotH->z0Error()        : 0.0;
  *evt_bs_sigmaZ     = haveBeamSpot ? beamSpotH->sigmaZ()         : 0.0;
  *evt_bs_sigmaZErr  = haveBeamSpot ? beamSpotH->sigmaZ0Error()   : 0.0;
  *evt_bs_dxdz       = haveBeamSpot ? beamSpotH->dxdz()           : 0.0;
  *evt_bs_dxdzErr    = haveBeamSpot ? beamSpotH->dxdzError()      : 0.0;
  *evt_bs_dydz       = haveBeamSpot ? beamSpotH->dydz()           : 0.0;	
  *evt_bs_dydzErr    = haveBeamSpot ? beamSpotH->dydzError()      : 0.0;
  *evt_bs_Xwidth     = haveBeamSpot ? beamSpotH->BeamWidthX()     : 0.0;
  *evt_bs_Ywidth     = haveBeamSpot ? beamSpotH->BeamWidthY()     : 0.0;
  *evt_bs_XwidthErr  = haveBeamSpot ? beamSpotH->BeamWidthXError(): 0.0;
  *evt_bs_YwidthErr  = haveBeamSpot ? beamSpotH->BeamWidthYError(): 0.0;

  const unsigned int covMatrix_dim = 7;

  for( unsigned int i = 0; i < covMatrix_dim; i++ ) {
    for( unsigned int j = 0; j < covMatrix_dim; j++ ) {
      evt_bs_covMatrix->push_back( beamSpotH->covariance(i, j) );
    }
  }  


  iEvent.put(evt_bs_p4           , "evtbsp4"        );
  iEvent.put(evt_bsType          , "evtbsType"      );

  if(haveBeamSpot) {
    iEvent.put(evt_bs_xErr       , "evtbsxErr"      );
    iEvent.put(evt_bs_yErr       , "evtbsyErr"      );
    iEvent.put(evt_bs_zErr       , "evtbszErr"      );
    iEvent.put(evt_bs_sigmaZ     , "evtbssigmaZ"    );
    iEvent.put(evt_bs_sigmaZErr  , "evtbssigmaZErr" );
    iEvent.put(evt_bs_dxdz       , "evtbsdxdz"      );
    iEvent.put(evt_bs_dxdzErr    , "evtbsdxdzErr"   );
    iEvent.put(evt_bs_dydz       , "evtbsdydz"      );
    iEvent.put(evt_bs_dydzErr    , "evtbsdydzErr"   );
    iEvent.put(evt_bs_Xwidth     , "evtbsXwidth"    );
    iEvent.put(evt_bs_Ywidth     , "evtbsYwidth"    );
    iEvent.put(evt_bs_XwidthErr  , "evtbsXwidthErr" );
    iEvent.put(evt_bs_YwidthErr  , "evtbsYwidthErr" );
    iEvent.put(evt_bs_covMatrix  , "evtcovMatrix"   );
  }
  
}


//define this as a plug-in
DEFINE_FWK_MODULE(BeamSpotMaker);
