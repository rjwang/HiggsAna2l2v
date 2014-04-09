import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.trigTools import *
from PhysicsTools.PatAlgos.producersLayer1.muonProducer_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.triggerLayer1.triggerMatchEmbedder_cfi import *
from FWCore.GuiBrowsers.ConfigToolBase import *
from PhysicsTools.PatAlgos.tools.helpers import *
from PhysicsTools.PatAlgos.patEventContent_cff import patTriggerL1RefsEventContent

"""
add all trigger matches to the leptons
"""
def addTriggerMatchingTo(process) :

    #triggers of interest
    pathTrigger_Mu_1  = 'path("HLT_Mu8_v*")'
    pathTrigger_Mu_2  = 'path("HLT_Mu12_v*")'
    pathTrigger_Mu_3  = 'path("HLT_Mu17_v*")'
    pathTrigger_Mu_4  = 'path("HLT_Mu24_*")'
    pathTrigger_Mu_5  = 'path("HLT_IsoMu24_*")'

    pathTrigger_Ele_1 = 'path("HLT_Ele8_*")'
    pathTrigger_Ele_2 = 'path("HLT_Ele8_CaloIdL_CaloIsoVL_v*")'
    pathTrigger_Ele_3 = 'path("HLT_Ele17_CaloIdL_CaloIsoVL_v*")'



    #muon trigger matching
    process.muonTriggerMatchHLTMu1 = cms.EDProducer("PATTriggerMatcherDRLessByR",
                                                            src     = cms.InputTag( 'selectedPatMuonsPFlow' ) ,
                                                            matched = cms.InputTag( 'patTrigger' ),    # selections of trigger objects ,
                                                            matchedCuts = cms.string( pathTrigger_Mu_1 ),    # selection of matches ,
                                                            maxDPtRel   = cms.double( 0.5 ), 
                                                            maxDeltaR   = cms.double( 0.3 ) ,
                                                            resolveAmbiguities    = cms.bool( True ) ,
                                                            resolveByMatchQuality = cms.bool( True ) 
                                                            )
    process.muonTriggerMatchHLTMu2  = process.muonTriggerMatchHLTMu1.clone( matchedCuts = cms.string( pathTrigger_Mu_2 ) )
    process.muonTriggerMatchHLTMu3  = process.muonTriggerMatchHLTMu1.clone( matchedCuts = cms.string( pathTrigger_Mu_3 ) )
    process.muonTriggerMatchHLTMu4  = process.muonTriggerMatchHLTMu1.clone( matchedCuts = cms.string( pathTrigger_Mu_4 ) )
    process.muonTriggerMatchHLTMu5  = process.muonTriggerMatchHLTMu1.clone( matchedCuts = cms.string( pathTrigger_Mu_5 ) )

    process.selectedPatMuonsTriggerMatch = cms.EDProducer( "PATTriggerMatchMuonEmbedder",
                                                           src     = cms.InputTag( "selectedPatMuonsPFlow" ),
                                                           matches = cms.VInputTag('muonTriggerMatchHLTMu1', 'muonTriggerMatchHLTMu2', 
                                                                                   'muonTriggerMatchHLTMu3', 'muonTriggerMatchHLTMu4',
                                                                                   'muonTriggerMatchHLTMu5'
                                                                                   )
                                                           )
    switchOnTriggerMatchEmbedding(process,
                                  triggerMatchers = [ 'muonTriggerMatchHLTMu1', 'muonTriggerMatchHLTMu2',
                                                      'muonTriggerMatchHLTMu3', 'muonTriggerMatchHLTMu4',
                                                      'muonTriggerMatchHLTMu5' ],
                                  sequence='patDefaultSequencePFlow')








    #electron trigger matching
    process.eleTriggerMatchHLTEle1 = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                                           src     = cms.InputTag( "selectedPatElectronsPFlow" ),
                                                           matched = cms.InputTag( "patTrigger"),
                                                           matchedCuts = cms.string(pathTrigger_Ele_1),
                                                           maxDPtRel = cms.double( 0.5 ),
                                                           maxDeltaR = cms.double( 0.3 ),
                                                           resolveAmbiguities    = cms.bool( True ),
                                                           resolveByMatchQuality = cms.bool( True )
                                                           )

    process.eleTriggerMatchHLTEle2 = process.eleTriggerMatchHLTEle1.clone( matchedCuts = cms.string(pathTrigger_Ele_2) )
    process.eleTriggerMatchHLTEle3 = process.eleTriggerMatchHLTEle1.clone( matchedCuts = cms.string(pathTrigger_Ele_3) )
    
    switchOnTriggerMatching( process,
                             ['eleTriggerMatchHLTEle1',
                              'eleTriggerMatchHLTEle2','eleTriggerMatchHLTEle3'],
                             sequence ='patDefaultSequencePFlow',
                             hltProcess = '*' )
    
    process.selectedPatElectronsWithTrigger = cms.EDProducer("PATTriggerMatchElectronEmbedder",
                                                             src     = cms.InputTag("selectedPatElectronsPFlow"),
                                                             matches = cms.VInputTag( cms.InputTag('eleTriggerMatchHLTEle1'),
                                                                                      cms.InputTag('eleTriggerMatchHLTEle2'),
                                                                                      cms.InputTag('eleTriggerMatchHLTEle3')
                                                                                      )
                                                             )
    
    removeCleaningFromTriggerMatching( process, sequence='patDefaultSequencePFlow' )
    process.patTrigger.processName    = cms.string( "*" )
    process.patTrigger.onlyStandAlone = cms.bool( False )
