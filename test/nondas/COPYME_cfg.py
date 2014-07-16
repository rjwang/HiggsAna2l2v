import FWCore.ParameterSet.Config as cms
import os,sys
isMC=True
gtag="START53_V23::All"

cfgFile=os.path.expandvars('${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/test/2l2vDataAnalyzer_nonDBS_cfg.py')
from CMGTools.HiggsAna2l2v.localPatTuples_cff import configureSourceFromCommandLine
castorDir, outFile, inputList = configureSourceFromCommandLine()
outFile='analysis.root'
inputList = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_0.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_1.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_2.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_3.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_4.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_5.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_6.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_7.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_8.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_9.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_10.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_11.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_12.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_13.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_14.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_15.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_16.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_17.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_18.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_19.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_20.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_21.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_22.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_23.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/EFT_Wimps_8TeV/COPYME_24.root',
) 
execfile(cfgFile)
