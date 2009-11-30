# Auto generated configuration file
# using: 
# Revision: 1.99.2.8 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: QCDForPF_cfi.py -s GEN,FASTSIM --pileup=NoPileUp --conditions=FrontierConditions_GlobalTag,IDEAL_V12::All --eventcontent=FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RECO -n 10 --relval 27000,1000 --no_exec --python_filename=fastsim_cms2_nofilt.py
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('FastSimulation/Configuration/RandomServiceInitialization_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator10TeV_cfi')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/Generator_cff')
process.load('FastSimulation/Configuration/FamosSequences_cff')
process.load('FastSimulation/Configuration/HLT_cff')
process.load('Configuration.StandardSequences.L1TriggerDefaultMenu_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('FastSimulation/Configuration/CommonInputs_cff')
process.load('FastSimulation/Configuration/EventContent_cff')

process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1.2.2 $'),
    annotation = cms.untracked.string('LHE nevts:10'),
    name = cms.untracked.string('PyCMS2')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)
process.options = cms.untracked.PSet(
    Rethrow = cms.untracked.vstring('ProductNotFound')
)
# Input source
process.source = cms.Source("PythiaSource",
    pythiaPylistVerbosity = cms.untracked.int32(0),
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(False),
    comEnergy = cms.untracked.double(10000.0),
    maxEventsToPrint = cms.untracked.int32(0),
    PythiaParameters = cms.PSet(
    pythiaUESettings = cms.vstring(
    'MSTJ(11)=3     ! Choice of the fragmentation function', 
    'MSTJ(22)=2     ! Decay those unstable particles', 
    'PARJ(71)=10 .  ! for which ctau  10 mm', 
    'MSTP(2)=1      ! which order running alphaS', 
    'MSTP(33)=0     ! no K factors in hard cross sections', 
    'MSTP(51)=10042 ! structure function chosen (external PDF CTEQ6L1)', 
    'MSTP(52)=2     ! work with LHAPDF', 
    'MSTP(81)=1     ! multiple parton interactions 1 is Pythia default', 
    'MSTP(82)=4     ! Defines the multi-parton model', 
    'MSTU(21)=1     ! Check on possible errors during program execution', 
    'PARP(82)=1.8387   ! pt cutoff for multiparton interactions', 
    'PARP(89)=1960. ! sqrts for which PARP82 is set', 
    'PARP(83)=0.5   ! Multiple interactions: matter distrbn parameter', 
    'PARP(84)=0.4   ! Multiple interactions: matter distribution parameter', 
    'PARP(90)=0.16  ! Multiple interactions: rescaling power', 
    'PARP(67)=2.5    ! amount of initial-state radiation', 
    'PARP(85)=1.0  ! gluon prod. mechanism in MI', 
    'PARP(86)=1.0  ! gluon prod. mechanism in MI', 
    'PARP(62)=1.25   ! ', 
    'PARP(64)=0.2    ! ', 
    'MSTP(91)=1      !', 
    'PARP(91)=2.1   ! kt distribution', 
    'PARP(93)=15.0  ! '),
    processParameters = cms.vstring(
    'MSEL=39                  ! All SUSY processes',
    'IMSS(1) = 11             ! Spectrum from external SLHA file',
    'IMSS(11) = 1             ! keeps gravitino mass from being overwritten',
    'IMSS(21) = 33            ! LUN number for SLHA File (must be 33)',
    'IMSS(22) = 33            ! Read-in SLHA decay table'), 
    parameterSets = cms.vstring('pythiaUESettings', 
                                'processParameters', 
                                'SLHAParameters'),
    SLHAParameters = cms.vstring('SLHAFILE = lesHouchesInput.lm1.out')    
    )
)

process.lepPt9Filter = cms.EDFilter("MCSingleParticleFilter",
                                    MinPt = cms.untracked.vdouble(     9,    9,    9,    9,    9,    9),
                                    MinEta = cms.untracked.vdouble( -2.7, -2.7, -2.7, -2.7, -2.7, -2.7),
                                    MaxEta = cms.untracked.vdouble(  2.7,  2.7,  2.7,  2.7,  2.7,  2.7),
                                    ParticleID = cms.untracked.vint32(11,  -11,   13,  -13,   15,  -15),
                                    Status =  cms.untracked.vint32(    3,    3,    3,    3,    3,    3)
                                    )


## Output definition
#process.output = cms.OutputModule("PoolOutputModule",
#    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#    fileName = cms.untracked.string('QCDForPF_cfi_py_GEN_FASTSIM.root'),
#    dataset = cms.untracked.PSet(
#        dataTier = cms.untracked.string('GEN-SIM-DIGI-RECO'),
#        filterName = cms.untracked.string('')
#    ),
#    SelectEvents = cms.untracked.PSet(
#        SelectEvents = cms.vstring('generation_step')
#    )
#)

# Additional output definition

# Other statements
process.famosPileUp.PileUpSimulator = process.PileUpSimulatorBlock.PileUpSimulator
process.famosPileUp.PileUpSimulator.averageNumber = 0
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.famosSimHits.ActivateDecays.comEnergy = 10000
process.simulation = cms.Sequence(process.simulationWithFamos)
process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)

# set correct vertex smearing
process.Early10TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Early10TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Early10TeVCollisionVtxSmearingParameters
# Apply Tracker misalignment
process.famosSimHits.ApplyAlignment = True
process.misalignedTrackerGeometry.applyAlignment = True

process.GlobalTag.globaltag = 'IDEAL_V12::All'

# Path and EndPath definitions
process.generation_step = cms.Path(cms.SequencePlaceholder("randomEngineStateProducer")+process.GeneInfo+process.genJetMET)
process.reconstruction = cms.Path(process.reconstructionWithFamos)
#process.out_step = cms.EndPath(process.output)

#----------------------------------------------------
#CMS2 includes
#----------------------------------------------------
process.load("CMS2.NtupleMaker.beamSpotMaker_cfi")
#CMS2Btagging
process.load("CMS2.NtupleMaker.bTaggingSequence_cfi")
#CMS2TrkBtagging
process.load("CMS2.NtupleMaker.bTaggingTrkSequence_cfi")
process.load("CMS2.NtupleMaker.bTagMaker_cfi")
process.load("CMS2.NtupleMaker.calotauMaker_cfi")
process.load("CMS2.NtupleMaker.candToGenAssMaker_cfi")
process.load("CMS2.NtupleMaker.conversionMaker_cfi")
#elCaloIsoSequence = cms.Sequence(egammaBasicClusterMerger*elCaloIsoMaker)
process.load("CMS2.NtupleMaker.elCaloIsoSequence_cff")
process.load("CMS2.NtupleMaker.electronMaker_cfi")
#electronSequence = cms.Sequence(uniqueElectrons*egammaIsolationSequenceCMS2*egammaElectronIDCMS2)
process.load("CMS2.NtupleMaker.electronSequence_cfi")
process.load("CMS2.NtupleMaker.elToJetAssMaker_cfi")
process.load("CMS2.NtupleMaker.elToMuAssMaker_cfi")
process.load("CMS2.NtupleMaker.elToTrackAssMaker_cfi")
process.load("CMS2.NtupleMaker.eventMaker_cfi")
process.load("CMS2.NtupleMaker.eventTrigMaker_cfi")
#gammaSequence = cms.Sequence(gamIsoDepositsCMS2 + gamIsoFromDepositsCMS2)
process.load("CMS2.NtupleMaker.gammaSequence_cfi") 
process.load("CMS2.NtupleMaker.genJetMaker_cfi")
process.load("CMS2.NtupleMaker.genMaker_cfi")
process.load("CMS2.NtupleMaker.hypDilepMaker_cfi")
process.load("CMS2.NtupleMaker.hypTrilepMaker_cfi")
process.load("CMS2.NtupleMaker.hypQuadlepMaker_cfi")
process.load("CMS2.NtupleMaker.jetMaker_cfi")
#cms2CaloJetSequence = cms.Sequence(prunedUncorrectedCMS2Jets)
process.load("CMS2.NtupleMaker.jetSequence_cff")
process.load("CMS2.NtupleMaker.jetToElAssMaker_cfi")
process.load("CMS2.NtupleMaker.jetToMuAssMaker_cfi")
#JPTCorrections = cms.Sequence(JetPlusTrackCorrections * L2L3CorJetSC5JPT * jptMaker)
process.load("CMS2.NtupleMaker.jptSequence_cff")
process.load("CMS2.NtupleMaker.l1DigiMaker_cfi")
process.load("CMS2.NtupleMaker.metMaker_cfi")
#metCorSequence
process.load("CMS2.NtupleMaker.metSequence_cff")
process.load("CMS2.NtupleMaker.muonMaker_cfi")
process.load("CMS2.NtupleMaker.muToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.muToJetAssMaker_cfi")
process.load("CMS2.NtupleMaker.muToTrackAssMaker_cfi")
process.load("CMS2.NtupleMaker.patElectronMaker_cfi")
process.load("CMS2.NtupleMaker.patJetMaker_cfi")
process.load("CMS2.NtupleMaker.patMETMaker_cfi")
process.load("CMS2.NtupleMaker.patMuonMaker_cfi")
process.load("CMS2.NtupleMaker.pdfinfoMaker_cfi")
process.load("CMS2.NtupleMaker.pfmetMaker_cfi")
process.load("CMS2.NtupleMaker.pftauMaker_cfi")
process.load("CMS2.NtupleMaker.photonMaker_cfi")
process.load("CMS2.NtupleMaker.scMaker_cfi")
process.load("CMS2.NtupleMaker.tcmetMaker_cfi")
process.load("CMS2.NtupleMaker.trackMaker_cfi")
process.load("CMS2.NtupleMaker.trackToElsAssMaker_cfi")
process.load("CMS2.NtupleMaker.trackToMuonAssMaker_cfi")
process.load("CMS2.NtupleMaker.triggerEventMaker_cfi")
process.load("CMS2.NtupleMaker.trkJetMaker_cfi")
process.load("CMS2.NtupleMaker.trkMuonFilter_cfi")
process.load("CMS2.NtupleMaker.vertexMaker_cfi")


###Dilepton Filter
process.load("CMS2.NtupleMaker.theFilter_cfi")

#-----------------------------------------------------------
# configure input data files and number of event to process
#-----------------------------------------------------------

#process.source = cms.Source("PoolSource",
#    skipEvents = cms.untracked.uint32(0),
#     fileNames = cms.untracked.vstring('file:/data/tmp/fgolf/1466166A-32CB-DD11-9779-001A92810AE2.root')
#)


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

#############   Define another L2 correction service #####
process.L2RelativeJetCorrector2 = cms.ESSource("L2RelativeCorrectionService", 
     tagName = cms.string('Summer08Redigi_L2Relative_SC5Calo'),
     label = cms.string('L2RelativeJetCorrector2')
)

#############   Define another L3 correction service #####
process.L3AbsoluteJetCorrector2 = cms.ESSource("L3AbsoluteCorrectionService", 
     tagName = cms.string('Summer08Redigi_L3Absolute_SC5Calo'),
     label = cms.string('L3AbsoluteJetCorrector2')
)

#############   Define the EMF correction service ####
process.L4EMFJetCorrector = cms.ESSource("L4EMFCorrectionService", 
    tagName = cms.string('CMSSW_152_L4EMF'),
    label = cms.string('L4EMFJetCorrector')
)

#############   Define the chain corrector service - L2L3 ###
process.L2L3JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector','L3AbsoluteJetCorrector'),
    label = cms.string('L2L3JetCorrector')
)

#############   Define the chain corrector service - L2L3L4 ###
process.L2L3L4JetCorrector = cms.ESSource("JetCorrectionServiceChain",  
    correctors = cms.vstring('L2RelativeJetCorrector2','L3AbsoluteJetCorrector2','L4EMFJetCorrector'),
    label = cms.string('L2L3L4JetCorrector')
)

# set the record's IOV. Must be defined once. Choose ANY correction service. #
process.prefer("L2L3L4JetCorrector") 

#-------------------------------------------------
# PAT configuration
#-------------------------------------------------

## change jet collection# PAT Layer 0+1
process.load("PhysicsTools.PatAlgos.patSequences_cff")
process.patDefaultSequence = cms.Sequence(
    process.beforeLayer1Objects *    # using '*', as the order is fixed.
    process.allLayer1Objects *
    process.selectedLayer1Objects
)

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process, 'prunedUncorrectedCMS2Jets', doJTA = True, doBTagging = True, jetCorrLabel = ('SC5', 'Calo'), doType1MET = True, genJetCollection = cms.InputTag("sisCone5GenJets") )

#new Flavor shit
process.load("PhysicsTools.HepMCCandAlgos.flavorHistoryPaths_cfi")

#-------------------------------------------------
# process output; first the event selection is
# defined: only those events that have passed the
# full production path are selected and written
# to file; the event content has been defined
# above
#-------------------------------------------------

## define event selection
process.EventSelection = cms.PSet(
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('path_cms2')
    )
)


process.out_CMS2 = cms.OutputModule("PoolOutputModule",
    process.EventSelection,
    verbose = cms.untracked.bool(True),
    dropMetaDataForDroppedData = cms.untracked.bool(True),
    fileName = cms.untracked.string('ntuple_lm1_1K.root')
)

process.out_CMS2.outputCommands = cms.untracked.vstring( 'drop *' )
process.out_CMS2.outputCommands.extend(cms.untracked.vstring('keep *_*Maker_*_*'))

#-------------------------------------------------
# process paths;
#-------------------------------------------------

process.CMS2Reco      = cms.Sequence(process.electronSequence * process.gammaSequence * process.caloJetMet*process.cms2CaloJetSequence * process.CMS2Btagging * process.metCorSequence)

process.eventmakers   = cms.Sequence(process.beamSpotMaker * process.vertexMaker * process.eventMaker * process.eventTrigMaker*process.pdfinfoMaker)
process.eventLevelMakers   = cms.Sequence(process.beamSpotMaker * process.vertexMaker  * process.eventMaker * process.pdfinfoMaker)
process.trigmmakers   = cms.Sequence(process.l1DigiMaker * process.triggerEventMaker)
process.genmakers     = cms.Sequence(process.genMaker * process.genjetmaker)
process.makers        = cms.Sequence(process.electronMaker * process.muonMaker * process.trackMaker * process.scMaker * process.jetMaker * process.JPTCorrections * process.trkmuonfilter * process.trkjetmaker * process.metMaker * process.tcmetMaker* process.recoJetAssociations*process.tautagging*process.calotauMaker*process.PFJetMet*process.PFTau*process.pftauMaker* process.CMS2TrkBtagging * process.photonMaker)
process.assmakers     = cms.Sequence(process.jetToMuAssMaker * process.jetToElAssMaker * process.muToElsAssMaker * process.sisCone5GenJetsPt10Seq*process.candToGenAssMaker * process.muToJetAssMaker * process.muToTrackAssMaker * process.elToTrackAssMaker * process.elToMuAssMaker *
                                     process.elToJetAssMaker * process.trackToMuonAssMaker * process.trackToElsAssMaker)
process.hypmakers     = cms.Sequence(process.hypDilepMaker * process.hypTrilepMaker * process.hypQuadlepMaker)
#process.flavorHistory = cms.Sequence(process.bFlavorHistoryProducer * process.cFlavorHistoryProducer)
process.othermakers   = cms.Sequence(process.elCaloIsoSequence * process.conversionMaker * process.bTagMaker * process.bTagTrkMaker)
process.pflowmakers   = cms.Sequence(process.pfmetMaker)
process.patmakers     = cms.Sequence(process.patMuonMaker * process.patElectronMaker * process.patJetMaker * process.patMETMaker)

process.cms2          = cms.Sequence(process.eventLevelMakers * process.trigmmakers * process.genmakers * process.makers * process.assmakers * process.othermakers * process.hypmakers)
#process.cms2          = cms.Sequence(process.trigmmakers * process.genmakers * process.makers * process.assmakers * process.othermakers * process.hypmakers)

process.load("CMS2.NtupleMaker.aSkimFilter_cfi")
process.path_cms2             = cms.Path(process.CMS2Reco * process.cms2 * process.pflowmakers
#                                         * process.caloJetMetGen
                                         * (process.genMETParticles+process.recoGenMET)
                                         * process.patDefaultSequence * process.patmakers * process.aSkimFilter)

process.outpath       = cms.EndPath(process.lepPt9Filter*process.eventTrigMaker*process.out_CMS2)


# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.HLTriggerFirstPath,process.HLT_L1Jet15,process.HLT_Jet30,process.HLT_Jet50,process.HLT_Jet80,process.HLT_Jet110,process.HLT_Jet180,process.HLT_Jet250,process.HLT_FwdJet20,process.HLT_DoubleJet150,process.HLT_DoubleJet125_Aco,process.HLT_DoubleFwdJet50,process.HLT_DiJetAve15,process.HLT_DiJetAve30,process.HLT_DiJetAve50,process.HLT_DiJetAve70,process.HLT_DiJetAve130,process.HLT_DiJetAve220,process.HLT_TripleJet85,process.HLT_QuadJet30,process.HLT_QuadJet60,process.HLT_SumET120,process.HLT_L1MET20,process.HLT_MET25,process.HLT_MET35,process.HLT_MET50,process.HLT_MET65,process.HLT_MET75,process.HLT_MET65_HT350,process.HLT_Jet180_MET60,process.HLT_Jet60_MET70_Aco,process.HLT_Jet100_MET60_Aco,process.HLT_DoubleJet125_MET60,process.HLT_DoubleFwdJet40_MET60,process.HLT_DoubleJet60_MET60_Aco,process.HLT_DoubleJet50_MET70_Aco,process.HLT_DoubleJet40_MET70_Aco,process.HLT_TripleJet60_MET60,process.HLT_QuadJet35_MET60,process.HLT_IsoEle15_L1I,process.HLT_IsoEle18_L1R,process.HLT_IsoEle15_LW_L1I,process.HLT_LooseIsoEle15_LW_L1R,process.HLT_Ele10_SW_L1R,process.HLT_Ele15_SW_L1R,process.HLT_Ele15_LW_L1R,process.HLT_EM80,process.HLT_EM200,process.HLT_DoubleIsoEle10_L1I,process.HLT_DoubleIsoEle12_L1R,process.HLT_DoubleIsoEle10_LW_L1I,process.HLT_DoubleIsoEle12_LW_L1R,process.HLT_DoubleEle5_SW_L1R,process.HLT_DoubleEle10_LW_OnlyPixelM_L1R,process.HLT_DoubleEle10_Z,process.HLT_DoubleEle6_Exclusive,process.HLT_IsoPhoton30_L1I,process.HLT_IsoPhoton10_L1R,process.HLT_IsoPhoton15_L1R,process.HLT_IsoPhoton20_L1R,process.HLT_IsoPhoton25_L1R,process.HLT_IsoPhoton40_L1R,process.HLT_Photon15_L1R,process.HLT_Photon25_L1R,process.HLT_DoubleIsoPhoton20_L1I,process.HLT_DoubleIsoPhoton20_L1R,process.HLT_DoublePhoton10_Exclusive,process.HLT_L1Mu,process.HLT_L1MuOpen,process.HLT_L2Mu9,process.HLT_IsoMu9,process.HLT_IsoMu11,process.HLT_IsoMu13,process.HLT_IsoMu15,process.HLT_Mu3,process.HLT_Mu5,process.HLT_Mu7,process.HLT_Mu9,process.HLT_Mu11,process.HLT_Mu13,process.HLT_Mu15,process.HLT_Mu15_Vtx2mm,process.HLT_DoubleIsoMu3,process.HLT_DoubleMu3,process.HLT_DoubleMu3_Vtx2mm,process.HLT_DoubleMu3_JPsi,process.HLT_DoubleMu3_Upsilon,process.HLT_DoubleMu7_Z,process.HLT_DoubleMu3_SameSign,process.HLT_DoubleMu3_Psi2S,process.HLT_BTagIP_Jet180,process.HLT_BTagIP_Jet120_Relaxed,process.HLT_BTagIP_DoubleJet120,process.HLT_BTagIP_DoubleJet60_Relaxed,process.HLT_BTagIP_TripleJet70,process.HLT_BTagIP_TripleJet40_Relaxed,process.HLT_BTagIP_QuadJet40,process.HLT_BTagIP_QuadJet30_Relaxed,process.HLT_BTagIP_HT470,process.HLT_BTagIP_HT320_Relaxed,process.HLT_BTagMu_DoubleJet120,process.HLT_BTagMu_DoubleJet60_Relaxed,process.HLT_BTagMu_TripleJet70,process.HLT_BTagMu_TripleJet40_Relaxed,process.HLT_BTagMu_QuadJet40,process.HLT_BTagMu_QuadJet30_Relaxed,process.HLT_BTagMu_HT370,process.HLT_BTagMu_HT250_Relaxed,process.HLT_DoubleMu3_BJPsi,process.HLT_DoubleMu4_BJPsi,process.HLT_TripleMu3_TauTo3Mu,process.HLT_IsoTau_MET65_Trk20,process.HLT_IsoTau_MET35_Trk15_L1MET,process.HLT_LooseIsoTau_MET30,process.HLT_LooseIsoTau_MET30_L1MET,process.HLT_DoubleIsoTau_Trk3,process.HLT_DoubleLooseIsoTau,process.HLT_IsoEle8_IsoMu7,process.HLT_IsoEle10_Mu10_L1R,process.HLT_IsoEle12_IsoTau_Trk3,process.HLT_IsoEle10_BTagIP_Jet35,process.HLT_IsoEle12_Jet40,process.HLT_IsoEle12_DoubleJet80,process.HLT_IsoEle5_TripleJet30,process.HLT_IsoEle12_TripleJet60,process.HLT_IsoEle12_QuadJet35,process.HLT_IsoMu14_IsoTau_Trk3,process.HLT_IsoMu7_BTagIP_Jet35,process.HLT_IsoMu7_BTagMu_Jet20,process.HLT_IsoMu7_Jet40,process.HLT_NoL2IsoMu8_Jet40,process.HLT_Mu14_Jet50,process.HLT_Mu5_TripleJet30,process.HLT_BTagMu_Jet20_Calib,process.HLT_ZeroBias,process.HLT_MinBias,process.HLT_MinBiasHcal,process.HLT_MinBiasEcal,process.HLT_MinBiasPixel,process.HLT_MinBiasPixel_Trk5,process.HLT_BackwardBSC,process.HLT_ForwardBSC,process.HLT_CSCBeamHalo,process.HLT_CSCBeamHaloOverlapRing1,process.HLT_CSCBeamHaloOverlapRing2,process.HLT_CSCBeamHaloRing2or3,process.HLT_TrackerCosmics,process.AlCa_IsoTrack,process.AlCa_EcalPhiSym,process.AlCa_EcalPi0,process.HLTriggerFinalPath,process.reconstruction,process.path_cms2,process.outpath)


process.eleIsoDepositTkCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.eleIsoDepositTkCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.eleIsoDepositEcalFromHitsCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.eleIsoDepositEcalFromHitsCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.eleIsoDepositHcalFromTowersCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.eleIsoDepositHcalFromTowersCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.gamIsoDepositTkCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.gamIsoDepositTkCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.gamIsoDepositEcalFromHitsCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.gamIsoDepositEcalFromHitsCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.gamIsoDepositHcalFromTowersCMS2.ExtractorPSet.barrelEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEB')
process.gamIsoDepositHcalFromTowersCMS2.ExtractorPSet.endcapEcalHits = cms.InputTag('caloRecHits:EcalRecHitsEE')
process.scMaker.ecalRecHitsInputTag_EB = cms.InputTag("caloRecHits","EcalRecHitsEB")
process.scMaker.ecalRecHitsInputTag_EE = cms.InputTag("caloRecHits","EcalRecHitsEE")
process.l1DigiMaker.l1extraModName = "l1extraParticles"        
#process.l1extraParticles.muonSource = cms.InputTag("l1ParamMuons","")
#process.eventMaker.haveTriggerInfo = cms.untracked.bool(False)

#process.Timing = cms.Service("Timing")
process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
process.options.wantSummary = cms.untracked.bool(True)

process.hypDilepMaker.TightLepton_PtCut = 10
process.hypDilepMaker.LooseLepton_PtCut =  5

for path in process.paths:
    getattr(process,path)._seq = process.lepPt9Filter*getattr(process,path)._seq
