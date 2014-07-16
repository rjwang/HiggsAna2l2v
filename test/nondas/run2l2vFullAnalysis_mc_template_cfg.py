import FWCore.ParameterSet.Config as cms
import os,sys
isMC=True
gtag="START53_V23::All"

cfgFile=os.path.expandvars('${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/test/2l2vDataAnalyzer_cfg.py')
from CMGTools.HiggsAna2l2v.localPatTuples_cff import configureSourceFromCommandLine
castorDir, outFile, inputList = configureSourceFromCommandLine()
outFile='analysis.root'
#inputList=HERE
execfile(cfgFile)
