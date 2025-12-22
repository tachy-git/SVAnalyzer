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
        'file:/T2_KR_KISTI/store/mc/RunIII2024Summer24MiniAODv6/QCD_Bin-PT-470to600_TuneCP5_13p6TeV_pythia8/MINIAODSIM/150X_mcRun3_2024_realistic_v2-v2/120000/0000c969-5555-4938-8bf9-02c56eb1694f.root'
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
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    genParticles = cms.InputTag("prunedGenParticles"),
)

print(">>> Adding analyzer to path")

process.p = cms.Path(process.analyzer)

print(">>> CMSSW config loaded successfully.")
