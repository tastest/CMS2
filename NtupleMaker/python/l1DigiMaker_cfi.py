import FWCore.ParameterSet.Config as cms

l1DigiMaker = cms.EDFilter("L1DigiMaker",
                           l1extraModName = cms.string("hltL1extraParticles"))
