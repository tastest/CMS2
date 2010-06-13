import FWCore.ParameterSet.Config as cms

pfMuToMuAssMaker = cms.EDFilter("PFMuToMuAssMaker",
                                aliasPrefix = cms.untracked.string("pfmus"),                                
                                minDR = cms.double(0.1),
                                musInputTag = cms.InputTag("muonMaker"),
                                pfmusInputTag = cms.InputTag("pfMuonMaker")
)


