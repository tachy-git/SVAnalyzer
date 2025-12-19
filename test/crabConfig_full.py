from CRABClient.UserUtilities import config
config = config()

# General settings
config.General.requestName = 'REQUESTNAME'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

# JobType settings
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'crab_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.numCores = 1
# config.JobType.allowUndistributedCMSSW = True  # Uncomment if using custom CMSSW

# Data settings
config.Data.inputDataset = 'INPUTDATASET'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1  # 2 files per job (adjust based on file size and processing time)
# config.Data.totalUnits = -1  # Process all files (default)
config.Data.publication = False  # Set to True if you want to publish output to DBS
config.Data.outputDatasetTag = 'OUTPUTDATASETTAG'

# Site settings
config.Site.storageSite = 'T2_KR_KISTI'  # Your storage site
# config.Site.whitelist = ['T2_KR_*', 'T2_US_*']  # Uncomment to specify allowed sites
# config.Site.blacklist = ['T2_XX_XXX']  # Uncomment to exclude specific sites

# Optional: Automatic splitting adjustment
# config.Data.splitting = 'Automatic'
