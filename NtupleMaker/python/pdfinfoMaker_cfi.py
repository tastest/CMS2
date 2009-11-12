import FWCore.ParameterSet.Config as cms

pdfinfoMaker = cms.EDFilter(
	"PDFInfoMaker",
        sourceHepMCTag = cms.InputTag("source")
)


