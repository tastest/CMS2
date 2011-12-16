import FWCore.ParameterSet.Config as cms

sParmMaker = cms.EDProducer("SParmMaker",
							sparm_mphiInputTag        = cms.InputTag("susyScanMassPhi"),
							sparm_mxhiInputTag        = cms.InputTag("susyScanMassXhi"),
							sparm_xsecInputTag      = cms.InputTag("susyScanCrossSection"),
							aliasPrefix             = cms.untracked.string("sparm")
)
