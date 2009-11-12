import FWCore.ParameterSet.Config as cms

eventTrigMaker = cms.EDFilter("EventTrigMaker",
    haveTriggerInfo = cms.untracked.bool(True)
)


