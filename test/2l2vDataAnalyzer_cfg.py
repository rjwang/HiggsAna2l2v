import FWCore.ParameterSet.Config as cms

process = cms.Process("PAT")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = gtag

process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),#True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )


## Input Module Configuration
process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring())
if(len(inputList)==0) : 
    if(isMC): 
	inputList = cms.untracked.vstring('/store/mc/Summer12_DR53X/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v3/0000/00277EBC-A2E0-E111-931D-00215E2223E2.root')
    else:
	inputList = cms.untracked.vstring('/store/data/Run2012D/DoubleMuParked/AOD/22Jan2013-v1/20000/002E1374-5F84-E211-83C4-20CF305616D0.root')
print inputList
process.source.fileNames=inputList

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) ) ## for testing


from CMGTools.HiggsAna2l2v.StandardSelections_cfi import getSelVersion
print 'Using the following global tag %s'%gtag


## Output Module Configuration
from CMGTools.HiggsAna2l2v.OutputConfiguration_cff import configureOutput
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands = cms.untracked.vstring('keep *'))
process.out.fileName = cms.untracked.string(outFile)




##-------------------- Import the JEC services -----------------------
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet = True
process.ak5PFJets.doAreaFastjet = True


################
# PRESELECTION #
################
#apply a good vertex selector and filter out scraping
from PhysicsTools.SelectorUtils.pvSelector_cfi import pvSelector
process.goodOfflinePrimaryVertices = cms.EDFilter( "PrimaryVertexObjectFilter",
                                                   filterParams = pvSelector.clone( minNdof = cms.double(4.0),
                                                                                    maxZ = cms.double(24.0),
                                                                                    maxd0 = cms.double(2.0)
                                                                                    ),
                                                   src=cms.InputTag('offlinePrimaryVertices')
                                                   )
process.goodVertexFilter = cms.EDFilter("GoodVertexFilter",
                                        vertexCollection = cms.InputTag('goodOfflinePrimaryVertices'),
                                        minimumNDOF = cms.uint32(4),
                                        maxAbsZ = cms.double(24),
                                        maxd0 = cms.double(2)
                                        )

process.noscraping = cms.EDFilter("FilterOutScraping",
                                  applyfilter = cms.untracked.bool(True),
                                  debugOn = cms.untracked.bool(False),
                                  numtrack = cms.untracked.uint32(10),
                                  thresh = cms.untracked.double(0.25)
                               )

# optional MET filters
# cf.https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters
#process.load('RecoMET.METFilters.metFilters_cff') ##missing package called RecoMET/METAnalyzer
#process.hcalLaserEventFilter.taggingMode   = cms.bool(True)
#process.EcalDeadCellTriggerPrimitiveFilter.taggingMode=cms.bool(True)
#process.eeBadScFilter.taggingMode           = cms.bool(True)
#process.ecalLaserCorrFilter.taggingMode     = cms.bool(True)
#process.trackingFailureFilter.VertexSource  = cms.InputTag('goodOfflinePrimaryVertices')
#process.trackingFailureFilter.taggingMode   = cms.bool(True)
#process.manystripclus53X.taggedMode         = cms.untracked.bool(True)
#process.manystripclus53X.forcedValue        = cms.untracked.bool(False)
#process.toomanystripclus53X.taggedMode      = cms.untracked.bool(True)
#process.toomanystripclus53X.forcedValue     = cms.untracked.bool(False)
#process.logErrorTooManyClusters.taggedMode  = cms.untracked.bool(True)
#process.logErrorTooManyClusters.forcedValue = cms.untracked.bool(False)

#process.metFilteringTaggers = cms.Sequence(process.HBHENoiseFilter*
#                                           process.hcalLaserEventFilter *
#                                           process.EcalDeadCellTriggerPrimitiveFilter *
#                                           process.eeBadScFilter *
#                                           process.ecalLaserCorrFilter *
#                                           process.trackingFailureFilter *
#                                           process.trkPOGFilters)

# optional MET filters : should add more? should run in tagging mode?
# cf.https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFilters                                                                                              
process.load('CommonTools.RecoAlgos.HBHENoiseFilter_cfi')                                                                                                           
process.metFilteringTaggers = cms.Sequence(process.HBHENoiseFilter)  


#PF2PAT
process.load("PhysicsTools.PatAlgos.patSequences_cff")
from PhysicsTools.PatAlgos.tools.pfTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *

postfix = "PFlow"
jetAlgo="AK5"
jecLevels=['L1FastJet', 'L2Relative', 'L3Absolute']
if(not isMC): jecLevels.append('L2L3Residual')

usePF2PAT(process,
          runPF2PAT=True,
          jetAlgo=jetAlgo,
          runOnMC=isMC,
          postfix=postfix,
          jetCorrections=('AK5PFchs', jecLevels),
          pvCollection=cms.InputTag('goodOfflinePrimaryVertices'),
          typeIMetCorrections=False)




#setup trigger matching
from CMGTools.HiggsAna2l2v.triggerMatching_cfg import *
addTriggerMatchingTo(process)




#custom electrons
useGsfElectrons(process,postfix=postfix,dR="03")
##process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
process.load('EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi')
process.eidMVASequence = cms.Sequence(  process.mvaTrigV0 + process.mvaNonTrigV0 )
process.patElectronsPFlow.electronIDSources.mvaTrigV0    = cms.InputTag("mvaTrigV0")
process.patElectronsPFlow.electronIDSources.mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0")
#from SHarper.HEEPAnalyzer.HEEPSelectionCuts_cfi import *
#process.selectedPatElectronsPFlowHeep = cms.EDProducer("HEEPAttStatusToPAT",
#                                                       eleLabel = cms.InputTag("selectedPatElectronsWithTrigger"),
#                                                       barrelCuts = cms.PSet(heepBarrelCuts),
#                                                       endcapCuts = cms.PSet(heepEndcapCuts),
#                                                       applyRhoCorrToEleIsol = cms.bool(True),
#                                                       eleIsolEffectiveAreas = cms.PSet (heepEffectiveAreas),
#                                                       eleRhoCorrLabel = cms.InputTag("kt6PFJets:rho"),
#                                                       verticesLabel = cms.InputTag("goodOfflinePrimaryVertices"),
#                                                       )

#custom muons
process.patMuonsPFlow.pfMuonSource = cms.InputTag("pfSelectedMuonsPFlow")
process.muonMatchPFlow.src = cms.InputTag("pfSelectedMuonsPFlow")

#custom jets for CHS
process.pfPileUpPFlow.checkClosestZVertex = cms.bool(False)
process.pfPileUpIsoPFlow.checkClosestZVertex = cms.bool(False)
getattr(process,"pfNoPileUp"+postfix).enable = True
getattr(process,"pfNoMuon"+postfix).enable = False     # to use muon-cleaned electron collection set to True (check iso)
getattr(process,"pfNoElectron"+postfix).enable = False # to use electron-cleaned tau collection set to True (check iso)
getattr(process,"pfNoTau"+postfix).enable = False      # to use tau-cleaned jet collection set to True (check what is a tau)
getattr(process,"pfNoJet"+postfix).enable = True       # this i guess it's for photons...      

#add q/g discriminator
#process.load('QuarkGluonTagger.EightTeV.QGTagger_RecoJets_cff')
#process.QGTagger.srcJets    = cms.InputTag("selectedPatJets"+postfix)
#process.QGTagger.isPatJet  = cms.untracked.bool(True)
#process.QGTagger.useCHS    = cms.untracked.bool(True)
#process.QGTagger.srcRho    = cms.InputTag('kt6PFJets','rho')
#process.QGTagger.srcRhoIso = cms.InputTag('kt6PFJetsCentral','rho')
#process.qgSequence=cms.Sequence(process.goodOfflinePrimaryVerticesQG+process.QGTagger)

#compute rho from central pf candidates only
from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
process.kt6PFJetsCentral = kt4PFJets.clone( rParam = cms.double(0.6),
                                            doAreaFastjet = cms.bool(True),
                                            doRhoFastjet = cms.bool(True),
                                            Rho_EtaMax = cms.double(2.5),
                                            Ghost_EtaMax = cms.double(2.5) )

#from CMGTools.HiggsAna2l2v.btvDefaultSequence_cff import *
#btvDefaultSequence(process,isMC,"selectedPatJets"+postfix,"goodOfflinePrimaryVertices")

# cf. https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag( cms.InputTag('pfMETcorrType0'),
                                                                 cms.InputTag('pfJetMETcorr', 'type1')
                                                                 )



######################################
# ANALYSIS                           #
######################################
#from CMGTools.HiggsAna2l2v.Analysis_cff import defineAnalysis
#defineAnalysis(process)
###### need to find out more details
    


#####################
#  PATH DEFINITION  #
#####################


#counters for specific filters
process.startCounter = cms.EDProducer("EventCountProducer")
process.scrapCounter = process.startCounter.clone()
process.vtxCounter   = process.startCounter.clone()
process.metCounter   = process.startCounter.clone()
process.endCounter = cms.EDProducer("EventCountProducer")

process.p = cms.Path( process.startCounter
                      *process.noscraping
                      *process.scrapCounter
                      *process.goodOfflinePrimaryVertices
                      *process.goodVertexFilter
                      *process.vtxCounter
                      *process.metFilteringTaggers
                      *process.metCounter
                      *process.eidMVASequence
                      *getattr(process,"patPF2PATSequence"+postfix)
                      #*process.btvSequence
                      *process.kt6PFJetsCentral
                      #*process.qgSequence
                      *process.type0PFMEtCorrection*process.producePFMETCorrections
                      *process.selectedPatElectronsWithTrigger#*process.selectedPatElectronsPFlowHeep
                      *process.selectedPatMuonsTriggerMatch
		      #*process.analysis
                      #*process.dataAnalyzer
		      #*process.endCounter*process.out
                      )



print '******************'
print process.out.fileName
print '*******************'
