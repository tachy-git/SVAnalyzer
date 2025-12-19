from CRABClient.UserUtilities import config
config = config()

# General settings
config.General.requestName = 'SVAnalysis_test_v1'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

# JobType settings
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'crab_cfg.py'  # Will use the CRAB-specific config
config.JobType.maxMemoryMB = 2500
config.JobType.numCores = 1

# Data settings
config.Data.inputDataset = '/QCD_Bin-PT-470to600_TuneCP5_13p6TeV_pythia8/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1  # 1 file per job for testing
config.Data.totalUnits = 5  # Process only 5 files for testing
config.Data.publication = False
config.Data.outputDatasetTag = 'SVAnalysis_test_v1'

# Site settings
config.Site.storageSite = 'T2_KR_KISTI'  # Change to your storage site
