import FWCore.ParameterSet.Config as cms

sParmMaker = cms.EDProducer("SParmMaker",
							sparm_m1InputTag        = cms.InputTag("susyScanMassP1"),
							sparm_m2InputTag        = cms.InputTag("susyScanMassP2"),
							sparm_m3InputTag        = cms.InputTag("susyScanMassP3"),
							sparm_m4InputTag        = cms.InputTag("susyScanMassP4"),
							sparm_xsecInputTag      = cms.InputTag("susyScanCrossSection"),
							aliasPrefix             = cms.untracked.string("sparm")
)
