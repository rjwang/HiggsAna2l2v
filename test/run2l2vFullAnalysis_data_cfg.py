import os,sys
isMC=False
runFull=True
runStd=True
gtag="FT_53_V21_AN4::All"
#gtag="FT_53_V21_AN6::All"

cfgFile=os.path.expandvars('${CMSSW_BASE}/src/CMGTools/HiggsAna2l2v/test/2l2vDataAnalyzer_cfg.py')
from CMGTools.HiggsAna2l2v.localPatTuples_cff import configureSourceFromCommandLine
castorDir, outFile, inputList = configureSourceFromCommandLine()
#outFile='analysis.root'
execfile(cfgFile)
