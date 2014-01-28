import os,sys
runOnMC=False
runFull=True
runStd=True
cfgFile=os.path.expandvars('${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/test/2l2vPattuple_cfg.py')
from CMGTools.HiggsAna2l2v.localPatTuples_cff import configureSourceFromCommandLine
castorDir, outFile, inputList = configureSourceFromCommandLine()
outFile='analysis.root'
execfile(cfgFile)
