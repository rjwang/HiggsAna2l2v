import FWCore.ParameterSet.Config as cms
import os,sys
isMC=True
gtag="START53_V23::All"

cfgFile=os.path.expandvars('${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/test/2l2vDataAnalyzer_nonDBS_cfg.py')
from CMGTools.HiggsAna2l2v.localPatTuples_cff import configureSourceFromCommandLine
castorDir, outFile, inputList = configureSourceFromCommandLine()
outFile='analysis.root'
inputList = cms.untracked.vstring('root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_0.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_1.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_2.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_3.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_4.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_5.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_6.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_7.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_8.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_9.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_10.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_11.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_12.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_13.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_14.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_15.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_16.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_17.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_18.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_19.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_20.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_21.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_22.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_23.root',
'root://eoscms.cern.ch//eos/cms/store/user/rewang/AODSIM/MC8TeV_FermionWIMP_V_M50_FullSim/MC8TeV_FermionWIMP_V_M50_FullSim_24.root',
) 
execfile(cfgFile)
