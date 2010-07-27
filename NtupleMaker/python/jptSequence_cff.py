import FWCore.ParameterSet.Config as cms

from RecoJets.JetAssociationProducers.trackExtrapolator_cfi import *
from RecoJets.JetAssociationProducers.ak5JTA_cff import *
ak5JetTracksAssociatorAtVertex.jets = cms.InputTag("ak5CaloJets")
ak5JetTracksAssociatorAtCaloFace.jets = cms.InputTag("ak5CaloJets")
ak5JetExtender.jets = cms.InputTag("ak5CaloJets")

from RecoJets.JetPlusTracks.JetPlusTrackCorrections_cff import *
JetPlusTrackZSPCorJetAntiKt5.src = cms.InputTag("ak5CaloJets")

JPTCorrections = cms.Sequence(trackExtrapolator * ak5JTA * JetPlusTrackCorrectionsAntiKt5)
