import FWCore.ParameterSet.Config as cms

patMETMaker = cms.EDFilter("PATMETMaker",
    # met collection
    patMETsInputTag = cms.InputTag("selectedLayer1METs")
)


