import FWCore.ParameterSet.Config as cms

from CMGTools.External.pujetidproducer_cfi import pileupJetIdProducer,pileupJetIdProducerChs
from CMGTools.HiggsAna2l2v.TriggerSequences_cff import getTriggerPaths
import os

def getSelVersion():
    try:
        cmssw_version = os.environ["CMSSW_VERSION"].replace("CMSSW_","")
    except:
        cmssw_version = "5_X"
    selVersion=2012
    if cmssw_version.startswith("4"):   selVersion=2011
    return selVersion

selVersion=getSelVersion()
print 'CMSSW version %s - selection adapted for %d'%(os.environ['CMSSW_BASE'],selVersion)



#  _______   _                                
# |__   __| (_)                       
#    | |_ __ _  __ _  __ _  ___ _ __   
#    | | '__| |/ _` |/ _` |/ _ \ '__|  
#    | | |  | | (_| | (_| |  __/ |     
#    |_|_|  |_|\__, |\__, |\___|_|  
#               __/ | __/ |                                    
#              |___/ |___/                                  
# 
DoubleElectronTrigs, DoubleMuTrigs, MuEGTrigs, PhotonTrigs, SingleMuTrigs, SingleEleTrigs, mcTrigs = getTriggerPaths(version=selVersion)

# base values for trigger event
BaseTriggerSelection = cms.PSet( source = cms.InputTag("TriggerResults::HLT"),
                                 triggerPaths = cms.PSet( gamma=cms.vstring(PhotonTrigs),
                                                          ee=cms.vstring(DoubleElectronTrigs),
                                                          mumu=cms.vstring(DoubleMuTrigs),
                                                          emu=cms.vstring(MuEGTrigs),
                                                          singleMu=cms.vstring(SingleMuTrigs),
							  singleEle=cms.vstring(SingleEleTrigs)
                                                          )
                                 )
#met filters
BaseMetFilters = cms.PSet( metfilters=cms.vstring('HBHENoiseFilter',
						  'hcalLaserEventFilter',
						  'EcalDeadCellTriggerPrimitiveFilter',
						  'eeBadScFilter',
						  'ecalLaserCorrFilter',
						  'trackingFailureFilter')
			 )

# base values for the vertex selection ------------------------------------------
BaseGeneratorSelection = cms.PSet( source = cms.InputTag("genParticles"),
                                   filterId = cms.int32(25),
                                   genJets=cms.InputTag("selectedPatJets")
                                   )
#if(selVersion==2011): BaseGeneratorSelection.genJets=cms.InputTag("ak5GenJets")


# base values for the vertex selection ------------------------------------------
BaseVertexSelection = cms.PSet( source = cms.InputTag("goodOfflinePrimaryVertices"),
                                beamSpot = cms.InputTag("offlineBeamSpot"),
                                maxZ = cms.double(24),
                                maxRho = cms.double(2.0),
                                minNDOF = cms.int32(4)
                                )


#  __  __                    
# |  \/  |                  
# | \  / |_   _  ___  _ __  
# | |\/| | | | |/ _ \| '_ \
# | |  | | |_| | (_) | | | |
# |_|  |_|\__,_|\___/|_| |_|
#
# base values for loose muon selection ----------------------------------------------
BaseMuonsSelection = cms.PSet( source = cms.InputTag("selectedPatMuonsTriggerMatch"), 
                               sourceIsPF = cms.bool(False),
                               rho25Neut = cms.InputTag("kt6PFJetsCentralNeutral:rho"), #but using BaseJetSelection:rho
                               minPt = cms.double(10),
                               maxEta = cms.double(2.5),
                               id = cms.string("loose"),
                               vbtf2011 = cms.PSet( id = cms.string(""),
                                                    minValidMuonHits=cms.int32(1),
                                                    minMatchingMuonStations = cms.int32(2),
                                                    minValidTrackerHits = cms.int32(11),
                                                    minPixelHits = cms.int32(1),
                                                    maxTrackChi2 = cms.double(10),
                                                    maxRelPtUncertainty = cms.double(0.1),
                                                    maxD0=cms.double(0.02),
                                                    maxDz=cms.double(0.1),
                                                    maxRelIso = cms.double(0.15),
                                                    applySoftMuonIsolationVeto=cms.bool(False)
                                                    ),
                               usePFIso = cms.bool(True),
                               reComputePFIso = cms.bool(True),
			       #will do BetaCorrection at selection level
                               maxRelIso = cms.double(999999.),
                               doDeltaBetaCorrection = cms.bool(False)
                               )

# base values for soft muon selection ----------------------------------------------
BaseSoftMuonsSelection = BaseMuonsSelection.clone( minPt = cms.double(3),
                                                    id=cms.string("soft"),
                                                    vbtf2011 = cms.PSet( id = cms.string("TMLastStationAngTight"),
                                                                         minValidMuonHits=cms.int32(0),
                                                                         minMatchingMuonStations = cms.int32(0),
                                                                         minValidTrackerHits = cms.int32(11),
                                                                         minPixelHits = cms.int32(0),
                                                                         maxTrackChi2 = cms.double(9999.),
                                                                         maxRelPtUncertainty = cms.double(9999.),
                                                                         maxD0=cms.double(0.2),
                                                                         maxDz=cms.double(0.2),
                                                                         maxRelIso = cms.double(999999.),
                                                                         applySoftMuonIsolationVeto=cms.bool(True) )
						 ) 





#  ______ _           _                      
# |  ____| |         | |                    
# | |__  | | ___  ___| |_ _ __ ___  _ __   
# |  __| | |/ _ \/ __| __| '__/ _ \| '_ \   
# | |____| |  __/ (__| |_| | | (_) | | | |  
# |______|_|\___|\___|\__|_|  \___/|_| |_|
#  
# base values for electron selection ----------------------------------------------
BaseElectronsSelection = cms.PSet( source = cms.InputTag("selectedPatElectronsPFlowHeep"),
                                   id=cms.string("veto"),
                                   #cf. https://twiki.cern.ch/twiki/bin/view/CMS/RegressionSCCorrections
                                   scCorrector = cms.string("${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/data/EleEnRegress.root"),
                                   minPt = cms.double(20),
                                   maxEta = cms.double(2.5),
                                   vetoTransitionElectrons = cms.bool(True),
                                   minDeltaRtoMuons = cms.double(0.1),
                                   usePFIso = cms.bool(True),
                                   reComputePFIso = cms.bool(True),
				   #will do BetaCorrection at selection level
                                   doDeltaBetaCorrection = cms.bool(False),
				   maxRelIso    = cms.double(999999.)
                                   )

# base values for loose electron selection ----------------------------------------------
BaseLooseElectronsSelection = BaseElectronsSelection.clone(minPt = cms.double(8))




# base values for the dilepton System selection ------------------------------------------
BaseDileptonSelection = cms.PSet( minDileptonMass = cms.double(0),
                                  maxDileptonMass = cms.double(7000),
				  minLegPt = cms.double(20),
                                  maxDz = cms.double(1.0)
                                  )


#  ____   _                
# |  _  \| |          _ 
# | |_)  | |         | | 
# | ____/| |__   ___ | |_ ___  _ __    
# | |    |  _ \ / _ \| __/ _ \| '_ \  
# | |    | | | | (_) | |_ (_) | | | |
# |_|    |_| |_|\___/ \__\___/|_| |_|
#  
# base values for photon selection ----------------------------------------------
BasePhotonsSelection = cms.PSet( source = cms.InputTag("photons"),
                                 conversions=cms.InputTag("allConversions"),
                                 gsfElectrons = cms.InputTag("gsfElectrons"),
                                 #cf. https://twiki.cern.ch/twiki/bin/view/CMS/RegressionSCCorrections
                                 scCorrector = cms.string("${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/data/PhoEnRegress.root"),
                                 rho25 = cms.InputTag("kt6PFJetsForIso:rho"),
                                 ebrechits = cms.InputTag("reducedEcalRecHitsEB"),
                                 eerechits = cms.InputTag("reducedEcalRecHitsEE"),
                                 minEt = cms.double(5), 
                                 maxEta = cms.double(2.5),
                                 maxHoE = cms.double(0.05),
                                 minSipipEB = cms.double(0.001),
                                 minSihihEB = cms.double(0.001),
                                 maxSihihEB = cms.double(0.011),
                                 maxSihihEE = cms.double(0.03),
                                 trkIsoCoeffsEB = cms.vdouble(2.0,  0.001,  0.0167),
                                 trkIsoCoeffsEE = cms.vdouble(2.0,  0.001,  0.0032),
                                 ecalIsoCoeffsEB = cms.vdouble(4.2, 0.006,  0.183),
                                 ecalIsoCoeffsEE = cms.vdouble(4.2, 0.006,  0.090),
                                 hcalIsoCoeffsEB = cms.vdouble(2.2, 0.0025, 0.062),
                                 hcalIsoCoeffsEE = cms.vdouble(2.2, 0.0025,  0.180),
                                 trackSource = cms.InputTag("generalTracks"),
                                 gsfTrackSource = cms.InputTag("gsfElectronTracks")
                                 )


#      _      _                                
#     | |    | |   
#     | | ___| |_   
# _   | |/ _ \ __|   
#| |__| |  __/ |_   
# \____/ \___|\__|  
#                                        
#my base values for jet selection -------------------------------------------------
BaseJetSelection = cms.PSet( source = cms.InputTag("selectedPatJets"),#selectedPatJetsPFlow"), #should be normal PFlow 
                             rho = cms.InputTag("kt6PFJets:rho"),
                             minPt = cms.double(10),
                             maxEta = cms.double(5.0),
                             minDeltaRtoLepton = cms.double(0.4),
                             puJetIds = pileupJetIdProducer.algos,
                             #jetTags = cms.VInputTag("simpleInclusiveSecondaryVertexHighEffBJetTags",# mySimpleInclusiveSecondaryVertexHighEffBJetTags",
                             #                        "mySimpleInclusiveSecondaryVertexHighPurBJetTags",
                             #                        "combinedInclusiveSecondaryVertexPositiveBJetTags")

			     jetTags = cms.VInputTag( "trackCountingHighPurBJetTags",                  #1: tchp
  						      "jetProbabilityBJetTags",                        #2: jp
  						      "simpleSecondaryVertexHighEffBJetTags",          #3: ssvhe
  						      "simpleInclusiveSecondaryVertexHighEffBJetTags", #4: ivf
  						      "combinedSecondaryVertexBJetTags",               #5: origcsv
						      "combinedSecondaryVertexRetrainedBJetTags",      #6: csv
						      "combinedCSVJPBJetTags",                         #7: jpcsv
						      "combinedCSVSLBJetTags",                         #8: slcsv
						      "combinedCSVJPSLBJetTags",                       #9: supercsv
						      "trackCountingHighEffBJetTags",		      #10: tche
					       	      "simpleSecondaryVertexHighPurBJetTags"          #11: ssvhp
						    )
                             )
#CHS jets
AssocJetSelection = BaseJetSelection.clone(source = cms.InputTag("selectedPatJetsPFlow"),
                                           puJetIds = pileupJetIdProducerChs.algos
                                           )



#  __  __ ______ _______ 
# |  \/  |  ____|__   __|
# | \  / | |__     | |   
# | |\/| |  __|    | |   
# | |  | | |____   | |   
# |_|  |_|______|  |_|   
#                                   
# base values for met selection -----------------------------------------------------
BaseMetSelection = cms.PSet( source = cms.InputTag("patMETsPFlow"),
                             #trksource = cms.InputTag("trackMetProducer"),
                             mainSources = cms.VInputTag("pfMETPFlow","pfType1CorrectedMet","pfType1p2CorrectedMet"
							   #"ClusteredPFMetProducer:assoc",                #1
                                                           #"ClusteredPFMetProducer:standard",             #2  
                                                           #"ClusteredPFMetProducer:central",              #3
                                                           #"ClusteredPFMetProducer:cleaned",              #4 
                                                           #"ClusteredPFMetProducer:assocCharged",         #5 
                                                           #"ClusteredPFMetProducer:assocWithFwd",         #6
                                                           #"ClusteredPFMetProducer:mvaMET",               #7
                                                           #"ClusteredPFMetProducer:assocWithFwd",      #8  //to be replaced by something else
                                                           #"ClusteredPFMetProducer:assoc",             #9  //to be replaced by something else
                                                           #"ClusteredPFMetProducer:assocWithFwd",      #10 //to be replaced by something else
                                                           #"ClusteredPFMetProducer:assocBeta",            #11
                                                           #"ClusteredPFMetProducer:assocWithFwdBeta",    #12
							  ),
                             #pfCands = cms.InputTag("particleFlow"),
                             #pvAssocCandidatesSource = cms.InputTag("ClusteredPFMetProducer:pvAssocCandidates"),
                             #sumEtSources = cms.InputTag("ClusteredPFMetProducer:globalPfMetSums")
                             )

