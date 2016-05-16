from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

runRange = '273017'

config.General.requestName = 'Charmonium_' + runRange 
config.General.workArea = 'NTuples_Prompt16'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../muonPogNtuples_cfg.py'
config.JobType.pyCfgParams = ['minMuPt=1.','minNMu=2','globalTag=80X_dataRun2_Prompt_v8','runOnMC=False','ntupleName=./muonPOGNtuple.root']

config.Data.inputDataset = '/Charmonium/Run2016B-PromptReco-v1/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/DCSOnly/json_DCSONLY.txt'
config.Data.runRange = runRange
config.Data.outLFNDirBase = '/store/user/battilan/MuonPOG/Ntuples/'
config.Data.publication = False
config.Data.outputDatasetTag = 'Charmonium_' + runRange

config.Site.storageSite = 'T2_CH_CERN' 

