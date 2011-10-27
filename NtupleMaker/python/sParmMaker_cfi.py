import FWCore.ParameterSet.Config as cms

sParmMaker = cms.EDProducer("SParmMaker",
							sparm_m0InputTag        = cms.InputTag("susyScanM0"),
							sparm_m12InputTag       = cms.InputTag("susyScanM12"),
							sparm_AInputTag         = cms.InputTag("susyScanA0"),
							sparm_muInputTag        = cms.InputTag("susyScanMu"),
							sparm_tanBetaInputTag   = cms.InputTag("susyScantanbeta"),
							sparm_xsecInputTag      = cms.InputTag("susyScanCrossSection"),
							sparm_susyScanBRInputTag       		= cms.InputTag("susyScanBR"),
							sparm_susyScanMassInputTag        	= cms.InputTag("susyScanMass"),
							sparm_susyScanRunInputTag        	= cms.InputTag("susyScanRun"),
							aliasPrefix             = cms.untracked.string("sparm")
)
