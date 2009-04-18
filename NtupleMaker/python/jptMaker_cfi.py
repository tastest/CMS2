import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
    # jpt collection
    #jptInputTag = cms.InputTag("JetPlusTrackZSPCorJetIcone5"),
    jptInputTag = cms.InputTag("JetPlusTrackZSPCorJetSisCone"),

    # calojet collection
#    ic5jetInputTag = cms.InputTag("iterativeCone5CaloJets")
    calojetInputTag = cms.InputTag("sisCone5CaloJets")
)


