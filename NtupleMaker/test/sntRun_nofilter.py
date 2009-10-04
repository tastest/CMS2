import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1.2.3 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

# initialize MessageLogger and output report ----------------------------------

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'WARNING'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# conditions ------------------------------------------------------------------

process.load("Configuration.StandardSequences.MixingNoPileUp_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.GlobalTag.globaltag = cms.string('MC_31X_V3::All')

# source ----------------------------------------------------------------------

#-----------------------------------------------------------
# configure input data files and number of event to process
#-----------------------------------------------------------

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
	'file:test_input.root'
	)
)

#----------------------------------------------------
#CMS2 includes
#----------------------------------------------------
process.load("CMS2.NtupleMaker.aSkimFilter_cfi")
process.load("CMS2.NtupleMaker.beamSpotMaker_cfi")
process.load("CMS2.NtupleMaker.bTaggingSequence_cfi")
process.load("CMS2.NtupleMaker.bTaggingTrkSequence_cfi")
process.load("CMS2.NtupleMaker.bTagMaker_cfi")
process.load("CMS2.NtupleMaker.bTagTrkMaker_cfi")
process.load("CMS2.NtupleMaker.calotauMaker_cfi")
process.load("CMS2.NtupleMaker.candToGenAssMaker_cfi")
process.load("CMS2.NtupleMaker.conversionMaker_cfi")
process.load("CMS2.NtupleMaker.dilepGenFilter_cfi")
process.load("CMS2.NtupleMaker.elCaloIsoSequence_cff")
process.load("CMS2.NtupleMaker.electronMaker_cfi")
process.load("CMS2.NtupleMaker.electronSequence_cfi")
process.load("CMS2.NtupleMaker.elToJetAssMaker_cfi")
process.load("CMS2.NtupleMaker.elToMuAssMaker_cfi")
process.load("CMS2.NtupleMaker.eventMaker_cfi")
process.load("CMS2.NtupleMaker.flavorHistorySequence_cfi")
process.load("CMS2.NtupleMaker.genJetSequence_cff")

process.load("CMS2.NtupleMaker.genMaker_cfi")
process.load("CMS2.NtupleMaker.hltMaker_cff")
process.load("CMS2.NtupleMaker.hypDilepMaker_cfi")
process.load("CMS2.NtupleMaker.hypGenMaker_cfi")
process.load("CMS2.NtupleMaker.hypTrilepMaker_cfi")
process.load("CMS2.NtupleMaker.hypQuadlepMaker_cfi")
process.load("CMS2.NtupleMaker.hypIsoMaker_cfi")
process.load("CMS2.NtupleMaker.jetSequence_cff")
process.load("CMS2.NtupleMaker.jetMaker_cfi")
process.load("CMS2.NtupleMaker.jetToElAssMaker_cfi")
process.load("CMS2.NtupleMaker.jetToMuAssMaker_cfi")
process.load("CMS2.NtupleMaker.jptSequence_cff")
process.load("CMS2.NtupleMaker.l1Maker_cfi")
process.l1Maker.fillL1Particles = cms.untracked.bool(False)
process.load("CMS2.NtupleMaker.metSequence_cff")
process.load("CMS2.NtupleMaker.metMaker_cfi")
process.load("CMS2.NtupleMaker.muonMaker_cfi")
process.load("CMS2.NtupleMaker.muToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.muToJetAssMaker_cfi")
process.load("CMS2.NtupleMaker.patElectronMaker_cfi")
process.load("CMS2.NtupleMaker.patJetMaker_cfi")
process.load("CMS2.NtupleMaker.patMETMaker_cfi")
process.load("CMS2.NtupleMaker.patMuonMaker_cfi")
process.load("CMS2.NtupleMaker.patElectronMaker_cfi")
process.load("CMS2.NtupleMaker.patJetMaker_cfi")
process.load("CMS2.NtupleMaker.patMETMaker_cfi")
process.load("CMS2.NtupleMaker.patMuonMaker_cfi")
process.load("CMS2.NtupleMaker.pdfinfoMaker_cfi")
process.load("CMS2.NtupleMaker.pfJetMaker_cfi")
process.load("CMS2.NtupleMaker.pfmetMaker_cfi")
process.load("CMS2.NtupleMaker.pftauMaker_cfi")
process.load("CMS2.NtupleMaker.photonMaker_cfi")
process.load("CMS2.NtupleMaker.scMaker_cfi")
process.load("CMS2.NtupleMaker.tcmetMaker_cfi")
process.load("CMS2.NtupleMaker.trackMaker_cfi")
process.load("CMS2.NtupleMaker.trackToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.trackToMuonAssMaker_cfi")
process.load("CMS2.NtupleMaker.trkJetMaker_cfi")
process.load("CMS2.NtupleMaker.trkJetSequence_cfi")
process.load("CMS2.NtupleMaker.vertexMaker_cfi")

###Dilepton Filter
process.load("CMS2.NtupleMaker.theFilter_cfi")

#-------------------------------------------------
# JEC
#-------------------------------------------------

#############   Define the L2 correction service #####
process.L2RelativeJetCorrector = cms.ESSource("L2RelativeCorrectionService", 
     tagName = cms.string('Summer08Redigi_L2Relative_SC5Calo'),
     label = cms.string('L2RelativeJetCorrector')
)

#############   Define the L3 correction service #####
process.L3AbsoluteJetCorrector = cms.ESSource("L3AbsoluteCorrectionService", 
     tagName = cms.string('Summer08Redigi_L3Absolute_SC5Calo'),
     label = cms.string('L3AbsoluteJetCorrector')
)

#############   Define the chain corrector service - L2L3 ###
process.L2L3JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector','L3AbsoluteJetCorrector'),
    label = cms.string('L2L3JetCorrector')
)

# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3JetCorrector") 

#-------------------------------------------------
# PAT configuration
#-------------------------------------------------
# TQAF/PAT Layer 1 ------------------------------------------------------------

process.load("PhysicsTools.PatAlgos.patSequences_cff")

# add more jet collections to the tqaf layer1

from PhysicsTools.PatAlgos.tools.jetTools import *

# switch standard PAT jet collection to AntiKt5
process.prunedUncorrectedCMS2Jets.inputUncorrectedJetCollection = cms.InputTag('antikt5CaloJets')
process.candToGenAssMaker.genJetsInputTag = cms.InputTag("antikt5GenJets")
switchJetCollection(process,
                    cms.InputTag('prunedUncorrectedCMS2Jets'),
                    doJTA=True,
                    doBTagging=True,
                    jetCorrLabel=('AK5','Calo'),
                    doType1MET=True,
                    genJetCollection=cms.InputTag("antikt5GenJets")
                    )


process.load("RecoJets.JetProducers.antikt5GenJets_cff")
process.antikt5StGenJets = process.antikt5GenJets.clone()
process.antikt5StGenJets.src  = cms.InputTag("genParticlesAllStables")
process.antikt5StGenJets.jetPtMin       = cms.double(0.)
process.antikt5StGenJets.alias = "ANTIKT5StGenJet"
process.genJetMaker.genJetsInputTag = cms.InputTag("antikt5StGenJets")
process.bFlavorHistoryProducer.matchedSrc = cms.InputTag("antikt5StGenJets")
process.cFlavorHistoryProducer.matchedSrc = cms.InputTag("antikt5StGenJets")
process.genJetSequence.replace(process.sisCone5StGenJets, process.antikt5StGenJets)
process.jetMaker.L2L3JetCorrectorName =cms.string("L2L3JetCorrectorAK5Calo")
process.metMuonJESCorAK5CMS2 = process.metMuonJESCorSC5CMS2.clone()
process.metMuonJESCorAK5CMS2.corrector = "L2L3JetCorrectorAK5Calo"
process.metCorSequence += process.metMuonJESCorAK5CMS2
process.metMaker.MuonJEScorMET_tag_ = cms.InputTag("metMuonJESCorAK5CMS2")

#switchJetCollection(process, cms.InputTag('prunedUncorrectedCMS2Jets'), doJTA = True, doBTagging = True, jetCorrLabel = ('SC5', 'Calo'), doType1MET = True, genJetCollection = cms.InputTag("sisCone5GenJets") )

#-------------------------------------------------
# process output; first the event selection is
# defined: only those events that have passed the
# full production path are selected and written
# to file; the event content has been defined
# above
#-------------------------------------------------

## define event selection
process.EventSelection = cms.PSet(
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('p1', 'p2')
#    )
)

process.out_CMS2 = cms.OutputModule(
	"PoolOutputModule",
    process.EventSelection,
    verbose = cms.untracked.bool(True),
    dropMetaData = cms.untracked.string("NONE"),
    fileName = cms.untracked.string('ntuple.root')
)

process.out_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.out_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))

#-------------------------------------------------
# process paths;
#-------------------------------------------------

process.CMS2Reco      = cms.Sequence(process.egammaElectronIDCMS2 * process.cms2CaloJetSequence * process.cms2TrkJetSequence * process.metCorSequence * process.CMS2Btagging * process.CMS2TrkBtagging * process.metCorSequence)

process.eventmakers   = cms.Sequence(process.beamSpotMaker * process.vertexMaker * process.eventMaker )

process.trigmakers   = cms.Sequence(process.l1Maker * process.hltMakerSequence)

process.genmakers     = cms.Sequence(process.genMaker * process.pdfinfoMaker * process.genJetSequence * process.CMS2FlavorHistorySequence * process.candToGenAssMaker)

process.makers        = cms.Sequence(process.trackMaker * process.muonMaker * process.electronMaker * process.scMaker * process.jetMaker * process.JPTCorrections * process.trkJetMaker * process.metMaker * process.tcmetMaker * process.calotauMaker * process.photonMaker)

process.assmakers     = cms.Sequence(process.jetToMuAssMaker * process.jetToElAssMaker * process.muToElsAssMaker * process.muToJetAssMaker * process.elToMuAssMaker * process.elToJetAssMaker * process.trackToMuonAssMaker * process.trackToElsAssMaker)

process.hypmakers     = cms.Sequence(process.hypDilepMaker * process.hypTrilepMaker * process.hypQuadlepMaker * process.hypIsoMaker  * process.hypGenMaker)

process.othermakers   = cms.Sequence(process.elCaloIsoSequence * process.conversionMaker * process.bTagMaker * process.bTagTrkMaker )

process.pflowmakers   = cms.Sequence(process.pfmetMaker * process.pfJetMaker * process.pftauMaker)

process.patmakers     = cms.Sequence(process.patMuonMaker * process.patElectronMaker * process.patJetMaker * process.patMETMaker)

process.cms2          = cms.Sequence(process.eventmakers * process.trigmakers * process.makers * process.genmakers * process.assmakers * process.othermakers * process.hypmakers)

process.all           = cms.Sequence( process.CMS2Reco * process.cms2 * process.patDefaultSequence * process.patmakers * process.pflowmakers )

process.pAll          = cms.Path( process.all )
#process.p1            = cms.Path( process.all * process.theFilter )

#process.p2            = cms.Path( process.all * process.dilepGenFilter)

process.outpath       = cms.EndPath(process.out_CMS2)

