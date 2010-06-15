import FWCore.ParameterSet.Config as cms

hltMaker = cms.EDProducer("HLTMaker",
    processName = cms.untracked.string("HLT"),
    aliasPrefix = cms.untracked.string("hlt"),                       
    fillTriggerObjects = cms.untracked.bool(True),
    prunedTriggerNames = cms.untracked.vstring(
        "HLT*Mu*",
        "HLT*Ele*",
        "HLT*EG*",
        "HLT*Photon*",
        "HLT*Jet*",
        "HLT*MET*",
    )
)
