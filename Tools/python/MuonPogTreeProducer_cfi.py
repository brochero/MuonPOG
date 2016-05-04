import FWCore.ParameterSet.Config as cms

MuonPogTree = cms.EDAnalyzer("MuonPogTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                             TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),

                             TrigFilterCut = cms.untracked.string("all"),
                             TrigPathCut   = cms.untracked.string("all"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             PFMetTag         = cms.untracked.InputTag("pfMet"), 
                             PFChMetTag       = cms.untracked.InputTag("pfChMet"), 
                             CaloMetTag       = cms.untracked.InputTag("caloMet"),

                             ScalersTag = cms.untracked.InputTag("scalersRawToDigi"),
                             # not filled for now!
                             # l1MuonsTag = cms.untracked.InputTag("gmtStage2Digis:Muon:"),
                             l1MuonsTag = cms.untracked.InputTag("none"),

                             GenTag = cms.untracked.InputTag("prunedGenParticles"), # pruned
                             PileUpInfoTag = cms.untracked.InputTag("addPileupInfo"),
                             GenInfoTag = cms.untracked.InputTag("generator"),
                             MinMuPtCut = cms.untracked.double(0.),
                             MinNMuCut  = cms.untracked.int32(0)
                             )



