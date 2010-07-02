// -*- C++ -*-
//
// Package:    PATJetMaker
// Class:      PATJetMaker
// 
/**\class PATJetMaker PATJetMaker.cc CMS2/NtupleMaker/src/PATJetMaker.cc

Description: copy additional PAT jet variables in simple data structures into the EDM event tree

 Implementation:
     - take PAT jets
     - extract and fill variables
*/
//
// Original Author:  pts/4
// Thu Jun 12 22:55:46 UTC 2008
// $Id: PATJetMaker.cc,v 1.13 2010/05/09 19:34:48 kalavase Exp $
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CMS2/NtupleMaker/interface/PATJetMaker.h"

#include "DataFormats/Math/interface/LorentzVector.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Flags.h"

#include "CMS2/NtupleMaker/interface/MatchUtilities.h"
#include "CMS2/NtupleMaker/interface/MCUtilities.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "CMS2/NtupleMaker/interface/CommonUtils.h"
//#include "PhysicsTools/PatUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
typedef math::XYZTLorentzVectorF LorentzVector;

//
// class decleration
//

//
// constructors and destructor
//
PATJetMaker::PATJetMaker(const edm::ParameterSet& iConfig)
{
  using namespace std;
  // product of this EDProducer 
  produces<vector<int> >     ("jetspatgenPartonid"      ).setBranchAlias("jets_pat_genParton_id"); // PAT gen parton ID
				   
  produces<vector<int> >     ("jetspatgenPartonMotherid").setBranchAlias("jets_pat_genPartonMother_id"); // PAT gen parton mother ID
  produces<vector<int> >     ("jetspatpartonFlavour"    ).setBranchAlias("jets_pat_partonFlavour"); // PAT parton flavour
  produces<vector<uint32_t> >("jetspatflag"             ).setBranchAlias("jets_pat_flag"); //PAT flag
  
  produces<vector<float> >   ("jetspatnoCorrF"          ).setBranchAlias("jets_pat_noCorrF"); // PAT corr for corrjet->nocorr
  produces<vector<float> >   ("jetspatjetCharge"        ).setBranchAlias("jets_pat_jetCharge"); // PAT jet charge

  //btagging info
  produces<vector<float> >   ("jetspatcombinedSecondaryVertexBJetTag"     ).setBranchAlias("jets_pat_combinedSecondaryVertexBJetTag");
  produces<vector<float> >   ("jetspatcombinedSecondaryVertexMVABJetTag"  ).setBranchAlias("jets_pat_combinedSecondaryVertexMVABJetTag");
  produces<vector<float> >   ("jetspatconeIsolationTauJetTag"             ).setBranchAlias("jets_pat_coneIsolationTauJetTag");
  produces<vector<float> >   ("jetspatimpactParameterMVABJetTag"          ).setBranchAlias("jets_pat_impactParameterMVABJetTag");
  produces<vector<float> >   ("jetspatjetBProbabilityBJetTag"             ).setBranchAlias("jets_pat_jetBProbabilityBJetTag");
  produces<vector<float> >   ("jetspatjetProbabilityBJetTag"              ).setBranchAlias("jets_pat_jetProbabilityBJetTag");
  produces<vector<float> >   ("jetspatsimpleSecondaryVertexHighEffBJetTag").setBranchAlias("jets_pat_simpleSecondaryVertexHighEffBJetTag");
  produces<vector<float> >   ("jetspatsimpleSecondaryVertexHighPurBJetTag").setBranchAlias("jets_pat_simpleSecondaryVertexHighPurBJetTag");
  produces<vector<float> >   ("jetspatsoftElectronByIP3dBJetTag"          ).setBranchAlias("jets_pat_softElectronByIP3dBJetTag");
  produces<vector<float> >   ("jetspatsoftElectronByPtBJetTag"            ).setBranchAlias("jets_pat_softElectronByPtBJetTag");
  produces<vector<float> >   ("jetspatsoftMuonBJetTag"                    ).setBranchAlias("jets_pat_softMuonBJetTag");
  produces<vector<float> >   ("jetspatsoftMuonByIP3dBJetTag"              ).setBranchAlias("jets_pat_softMuonByIP3dBJetTag");
  produces<vector<float> >   ("jetspatsoftMuonByPtBJetTag"                ).setBranchAlias("jets_pat_softMuonByPtBJetTag");
  produces<vector<float> >   ("jetspattrackCountingHighEffBJetTag"        ).setBranchAlias("jets_pat_trackCountingHighEffBJetTag");
  produces<vector<float> >   ("jetspattrackCountingHighPurBJetTag"        ).setBranchAlias("jets_pat_trackCountingHighPurBJetTag");
  

  produces<vector<LorentzVector> > ("jetspatgenPartonp4"      ).setBranchAlias("jets_pat_genParton_p4"); // PAT gen parton p4
  produces<vector<LorentzVector> > ("jetspatgenPartonMotherp4").setBranchAlias("jets_pat_genPartonMother_p4"); // PAT gen parton mother p4
  produces<vector<LorentzVector> > ("jetspatgenJetp4"         ).setBranchAlias("jets_pat_genJet_p4"); // PAT gen jet p4
  produces<vector<LorentzVector> > ("jetspatjetp4"            ).setBranchAlias("jets_pat_jet_p4"); // PAT jet p4
  produces<vector<LorentzVector> > ("jetspatjetuncorp4"       ).setBranchAlias("jets_pat_jet_uncorp4"); // PAT jet p4

  //jetID flags
  produces<vector<int> > ("jetspatjetIDMinimal"     ).setBranchAlias("jets_pat_jetIDMinimal");
  produces<vector<int> > ("jetspatjetIDLooseAOD"    ).setBranchAlias("jets_pat_jetIDLooseAOD");
  produces<vector<int> > ("jetspatjetIDLoose"       ).setBranchAlias("jets_pat_jetIDLoose");
  produces<vector<int> > ("jetspatjetIDTight"       ).setBranchAlias("jets_pat_jetIDTight");
  // parameters from configuration
  patJetsInputTag_   = iConfig.getParameter<edm::InputTag>("patJetsInputTag"   );
  uncorRecoJetsTag_  = iConfig.getParameter<edm::InputTag>("uncorRecoJetsTag"  );

}


PATJetMaker::~PATJetMaker()
{
 
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void PATJetMaker::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace std;

  // get jet collection
  edm::Handle<vector<pat::Jet> > patJetsHandle;
  iEvent.getByLabel(patJetsInputTag_, patJetsHandle);
  vector<pat::Jet> v_patJets = *(patJetsHandle.product());

  edm::Handle<vector<reco::CaloJet> > uncorRecoJetsHandle;
  iEvent.getByLabel(uncorRecoJetsTag_, uncorRecoJetsHandle);
  vector<reco::CaloJet> v_uncorRecoJets = *(uncorRecoJetsHandle.product());

  MatchUtilities::alignRecoPatJetCollections(v_uncorRecoJets, v_patJets);

  // create containers
  auto_ptr<vector<int> >       jets_patgenParton_id(new vector<int>);
  auto_ptr<vector<int> >       jets_patgenPartonMother_id(new vector<int>);
  auto_ptr<vector<int> >       jets_patpartonFlavour(new vector<int>);
  auto_ptr<vector<uint32_t> >  jets_patflag(new vector<uint32_t>);

  auto_ptr<vector<float> >     jets_patnoCorrF(new vector<float>);
  auto_ptr<vector<float> >     jets_patjetCharge(new vector<float>);
  
  //btagging
  auto_ptr<vector<float> >   jets_pat_combinedSecondaryVertexBJetTag      (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_combinedSecondaryVertexMVABJetTag   (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_coneIsolationTauJetTag              (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_impactParameterMVABJetTag           (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_jetBProbabilityBJetTag              (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_jetProbabilityBJetTag               (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_simpleSecondaryVertexHighEffBJetTag (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_simpleSecondaryVertexHighPurBJetTag (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_softElectronByIP3dBJetTag           (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_softElectronByPtBJetTag             (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_softMuonBJetTag                     (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_softMuonByIP3dBJetTag               (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_softMuonByPtBJetTag                 (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_trackCountingHighEffBJetTag         (new vector<float>    );
  auto_ptr<vector<float> >   jets_pat_trackCountingHighPurBJetTag         (new vector<float>    );

  auto_ptr<vector<LorentzVector> > jets_patgenParton_p4(new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > jets_patgenPartonMother_p4(new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > jets_patgenJet_p4(new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > jets_patjet_p4(new vector<LorentzVector>);
  auto_ptr<vector<LorentzVector> > jets_patjet_uncorp4(new vector<LorentzVector>);

  auto_ptr<vector<int> > jets_pat_jetIDMinimal(new vector<int>);
  auto_ptr<vector<int> > jets_pat_jetIDLooseAOD(new vector<int>);
  auto_ptr<vector<int> > jets_pat_jetIDLoose(new vector<int>);
  auto_ptr<vector<int> > jets_pat_jetIDTight(new vector<int>);
    

  JetIDSelectionFunctor jetIDMinimal( JetIDSelectionFunctor::CRAFT08,
				      JetIDSelectionFunctor::MINIMAL );
  
  JetIDSelectionFunctor jetIDLooseAOD( JetIDSelectionFunctor::CRAFT08,
				       JetIDSelectionFunctor::LOOSE_AOD );
  
  JetIDSelectionFunctor jetIDLoose( JetIDSelectionFunctor::CRAFT08,
                                    JetIDSelectionFunctor::LOOSE );
  
  JetIDSelectionFunctor jetIDTight( JetIDSelectionFunctor::CRAFT08,
                                    JetIDSelectionFunctor::TIGHT );
  
  pat::strbitset ret = jetIDLoose.getBitTemplate();
  

  // loop over jets and fill containers
  vector<pat::Jet>::const_iterator patJetsEnd = v_patJets.end(); 
  for ( vector<pat::Jet>::const_iterator patJet = v_patJets.begin();
	patJet != patJetsEnd; 
	++patJet) {

    reco::GenParticle genParton = patJet->genParton() ? *patJet->genParton() : 
      reco::GenParticle(0, reco::Particle::LorentzVector(0, 0, 0, 0), reco::Particle::Point(0,0,0), 0, 0, true);

    const reco::GenParticle *mother = MCUtilities::motherID(genParton);

    
    jets_patgenParton_id->push_back(genParton.pdgId());
    jets_patgenPartonMother_id->push_back(mother ? mother->pdgId() : 0);
    jets_patpartonFlavour->push_back(patJet->partonFlavour());
    jets_patflag->push_back(patJet->status());

    
    float jetCorF  = patJet->p4().pt()/patJet->originalObject()->p4().pt(); 
           
    jets_patnoCorrF->push_back(1/jetCorF); 
    jets_patjetCharge->push_back(patJet->jetCharge());



    jets_pat_combinedSecondaryVertexBJetTag        ->push_back(CommonUtils::isinf(patJet->bDiscriminator("combinedSecondaryVertexBJetTags"))
							       ? -9999. : patJet->bDiscriminator("combinedSecondaryVertexBJetTags")        );
    jets_pat_combinedSecondaryVertexMVABJetTag     ->push_back(CommonUtils::isinf(patJet->bDiscriminator("combinedSecondaryVertexMVABJetTags"))
							       ? -9999. : patJet->bDiscriminator("combinedSecondaryVertexMVABJetTags")        );
    jets_pat_coneIsolationTauJetTag                ->push_back(CommonUtils::isinf(patJet->bDiscriminator("coneIsolationTauJetTags"))
							       ? -9999. : patJet->bDiscriminator("coneIsolationTauJetTags")        );
    jets_pat_impactParameterMVABJetTag             ->push_back(CommonUtils::isinf(patJet->bDiscriminator("impactParameterMVABJetTags"))
							       ? -9999. : patJet->bDiscriminator("impactParameterMVABJetTags")        );
    jets_pat_jetBProbabilityBJetTag                ->push_back(CommonUtils::isinf(patJet->bDiscriminator("jetBProbabilityBJetTags"))
							       ? -9999. : patJet->bDiscriminator("jetBProbabilityBJetTags")        );
    jets_pat_jetProbabilityBJetTag                 ->push_back(CommonUtils::isinf(patJet->bDiscriminator("jetProbabilityBJetTags"))
							       ? -9999. : patJet->bDiscriminator("jetProbabilityBJetTags")        );
    jets_pat_simpleSecondaryVertexHighEffBJetTag   ->push_back(CommonUtils::isinf(patJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"))
							       ? -9999. : patJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")        );
    jets_pat_simpleSecondaryVertexHighPurBJetTag   ->push_back(CommonUtils::isinf(patJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"))
							       ? -9999. : patJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags")        );
    jets_pat_softElectronByIP3dBJetTag             ->push_back(CommonUtils::isinf(patJet->bDiscriminator("softElectronByIP3dBJetTags"))
							       ? -9999. : patJet->bDiscriminator("softElectronByIP3dBJetTags")        );
    jets_pat_softElectronByPtBJetTag               ->push_back(CommonUtils::isinf(patJet->bDiscriminator("softElectronByPtBJetTags"))
							       ? -9999. : patJet->bDiscriminator("softElectronByPtBJetTags")        );
    jets_pat_softMuonBJetTag                       ->push_back(CommonUtils::isinf(patJet->bDiscriminator("softMuonBJetTags"))
							       ? -9999. : patJet->bDiscriminator("softMuonBJetTags")        );
    jets_pat_softMuonByIP3dBJetTag                 ->push_back(CommonUtils::isinf(patJet->bDiscriminator("softMuonByIP3dBJetTags"))
							       ? -9999. : patJet->bDiscriminator("softMuonByIP3dBJetTags")        );
    jets_pat_softMuonByPtBJetTag                   ->push_back(CommonUtils::isinf(patJet->bDiscriminator("softMuonByPtBJetTags"))
							       ? -9999. : patJet->bDiscriminator("softMuonByPtBJetTags")        );
    jets_pat_trackCountingHighEffBJetTag           ->push_back(CommonUtils::isinf(patJet->bDiscriminator("trackCountingHighEffBJetTags"))
							       ? -9999. : patJet->bDiscriminator("trackCountingHighEffBJetTags")        );
    jets_pat_trackCountingHighPurBJetTag           ->push_back(CommonUtils::isinf(patJet->bDiscriminator("trackCountingHighPurBJetTags"))
							       ? -9999. : patJet->bDiscriminator("trackCountingHighPurBJetTags")        );
    
    
    
    jets_patgenParton_p4->push_back( LorentzVector( genParton.p4() ) );
    jets_patgenPartonMother_p4->push_back(mother ? LorentzVector( mother->p4() ) : LorentzVector(0,0,0,0) );
    LorentzVector genJetP4 = patJet->genJet() ? LorentzVector( patJet->genJet()->p4() ) : LorentzVector(0, 0, 0, 0);
    jets_patgenJet_p4->push_back(genJetP4);
    jets_patjet_p4->push_back( LorentzVector( patJet->p4() ) );
    jets_patjet_uncorp4->push_back( LorentzVector( patJet->originalObject()->p4() ) );

  
    ret.set(false);
    jets_pat_jetIDMinimal->push_back(jetIDMinimal(*patJet, ret));
    ret.set(false);
    jets_pat_jetIDLooseAOD->push_back(jetIDLooseAOD(*patJet, ret));
    ret.set(false);
    jets_pat_jetIDLoose->push_back(jetIDLoose(*patJet, ret));
    ret.set(false);
    jets_pat_jetIDTight->push_back(jetIDTight(*patJet, ret));
    
  }
  iEvent.put(jets_patgenParton_id, "jetspatgenPartonid");
  iEvent.put(jets_patgenPartonMother_id, "jetspatgenPartonMotherid");
  iEvent.put(jets_patpartonFlavour, "jetspatpartonFlavour");
  iEvent.put(jets_patflag, "jetspatflag");     
  iEvent.put(jets_patnoCorrF, "jetspatnoCorrF");
  iEvent.put(jets_patjetCharge, "jetspatjetCharge");



  iEvent.put(jets_pat_combinedSecondaryVertexBJetTag,     "jetspatcombinedSecondaryVertexBJetTag");
  iEvent.put(jets_pat_combinedSecondaryVertexMVABJetTag,     "jetspatcombinedSecondaryVertexMVABJetTag");
  iEvent.put(jets_pat_coneIsolationTauJetTag,     "jetspatconeIsolationTauJetTag");
  iEvent.put(jets_pat_impactParameterMVABJetTag,     "jetspatimpactParameterMVABJetTag");
  iEvent.put(jets_pat_jetBProbabilityBJetTag,     "jetspatjetBProbabilityBJetTag");
  iEvent.put(jets_pat_jetProbabilityBJetTag,     "jetspatjetProbabilityBJetTag");
  iEvent.put(jets_pat_simpleSecondaryVertexHighEffBJetTag,     "jetspatsimpleSecondaryVertexHighEffBJetTag");
  iEvent.put(jets_pat_simpleSecondaryVertexHighPurBJetTag,     "jetspatsimpleSecondaryVertexHighPurBJetTag");
  iEvent.put(jets_pat_softElectronByIP3dBJetTag,     "jetspatsoftElectronByIP3dBJetTag");
  iEvent.put(jets_pat_softElectronByPtBJetTag,     "jetspatsoftElectronByPtBJetTag");
  iEvent.put(jets_pat_softMuonBJetTag,     "jetspatsoftMuonBJetTag");
  iEvent.put(jets_pat_softMuonByIP3dBJetTag,     "jetspatsoftMuonByIP3dBJetTag");
  iEvent.put(jets_pat_softMuonByPtBJetTag,     "jetspatsoftMuonByPtBJetTag");
  iEvent.put(jets_pat_trackCountingHighEffBJetTag,     "jetspattrackCountingHighEffBJetTag");
  iEvent.put(jets_pat_trackCountingHighPurBJetTag,     "jetspattrackCountingHighPurBJetTag");
  

  iEvent.put(jets_patgenParton_p4, "jetspatgenPartonp4");
  iEvent.put(jets_patgenPartonMother_p4, "jetspatgenPartonMotherp4");
  iEvent.put(jets_patgenJet_p4, "jetspatgenJetp4");
  iEvent.put(jets_patjet_p4, "jetspatjetp4");
  iEvent.put(jets_patjet_uncorp4, "jetspatjetuncorp4");

  iEvent.put(jets_pat_jetIDMinimal, "jetspatjetIDMinimal");
  iEvent.put(jets_pat_jetIDLooseAOD, "jetspatjetIDLooseAOD");
  iEvent.put(jets_pat_jetIDLoose, "jetspatjetIDLoose");
  iEvent.put(jets_pat_jetIDTight, "jetspatjetIDTight");
  

}

// ------------ method called once each job just before starting event loop  ------------
void 
PATJetMaker::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
PATJetMaker::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(PATJetMaker);