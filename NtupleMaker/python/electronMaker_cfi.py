import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter(
	"ElectronMaker",
    # Electron collection
    electronsInputTag = cms.InputTag("gsfElectrons"),
    # Beamspot
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    # Isolation
    ecalIsoTag = cms.InputTag("eleIsoFromDepsEcalFromHitsCMS2"),
    hcalIsoTag = cms.InputTag("eleIsoFromDepsHcalFromTowersCMS2"),
    tkIsoTag = cms.InputTag("eleIsoFromDepsTkCMS2"),
    # egamma ID
    eidRobustLooseTag = cms.InputTag("eidRobustLooseCMS2"),
    eidRobustTightTag = cms.InputTag("eidRobustTightCMS2"),
    eidRobustHighEnergyTag = cms.InputTag("eidRobustHighEnergyCMS2"),
    eidLooseTag = cms.InputTag("eidLooseCMS2"),
    eidTightTag = cms.InputTag("eidTightCMS2")

)

