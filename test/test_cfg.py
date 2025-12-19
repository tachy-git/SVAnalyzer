import FWCore.ParameterSet.Config as cms

process = cms.Process("SVAnalysis")

# Message logger
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

print(">>> Loading input files...")

# Max events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:1b799c1b-8ed6-4d56-a950-f6756e8e6273.root'
    )
)

print(">>> Setting TFileService")

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_test.root")
)

# Analyzer
process.analyzer = cms.EDAnalyzer("SVAnalyzer",
    jets = cms.InputTag("slimmedJets"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    secondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    muons = cms.InputTag("slimmedMuons")
)

print(">>> Adding analyzer to path")

process.p = cms.Path(process.analyzer)

print(">>> CMSSW config loaded successfully.")
