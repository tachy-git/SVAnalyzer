import FWCore.ParameterSet.Config as cms

process = cms.Process("SVAnalysis")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Maximum number of events (-1 means all events in the file)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input files (will be overridden by CRAB)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

# TFileService for output
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('output.root')
)

# Analyzer
process.analyzer = cms.EDAnalyzer('SVAnalyzer',
    jets = cms.InputTag("slimmedJets"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    muons = cms.InputTag("slimmedMuons")
)

process.p = cms.Path(process.analyzer)
