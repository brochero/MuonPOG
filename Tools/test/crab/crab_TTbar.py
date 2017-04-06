from CRABClient.UserUtilities import config
config = config()

config.section_("General")
config.General.requestName = 'RelValTTbar_13TeV_HIPSIMU1p2'
config.General.workArea = 'NTuples_RelValTTbar_13TeV-NoMITI'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../muonPogNtuples_cfg.py'   
#config.JobType.pyCfgParams = ['hltPathFilter=IsoMu20']
config.JobType.priority = 1

config.section_("Data")
# HIP
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2-v1/GEN-SIM-RECO'
# No HIP
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_HIPoff-v1/GEN-SIM-RECO'

# HIPSIMU/HIPMITI
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip0p6-v1/GEN-SIM-RECO'
# HIPSIMU/NoHIPMITI
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip0p6_mtoff-v1/GEN-SIM-RECO'

# HIPSIMU/HIPMITI
#config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip0p8-v1/GEN-SIM-RECO'
# HIPSIMU/NoHIPMITI
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip0p8_mtoff-v1/GEN-SIM-RECO'

# HIPSIMU/HIPMITI
# config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip1p2-v1/GEN-SIM-RECO'
# HIPSIMU/NoHIPMITI
config.Data.inputDataset =  '/RelValTTbar_13/CMSSW_8_1_0_pre9-PU25ns_81X_mcRun2_asymptotic_v2_hip1p2_mtoff-v1/GEN-SIM-RECO'

config.Data.allowNonValidInputDataset = True
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.totalUnits = 10
config.Data.ignoreLocality  = True


config.section_("Site")
# config.Site.storageSite = 'T2_ES_IFCA'

#config.Data.outLFNDirBase    = '/store/group/phys_muon/Commissioning/Ntuples/Commissioning2016/mcSpring16/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
config.Site.storageSite = 'T2_CH_CERN'
