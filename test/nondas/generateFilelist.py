#!/usr/bin/env python
import os,sys
import getopt
import commands

fromfile = 'run2l2vFullAnalysis_mc_template_cfg.py'
tofile = open('mc_cfg.py',"w")
dir = '/store/user/rewang/AODSIM/HZZd/HZZd_MZd-60GeV_mH-126GeV_8TeV_D3_evts10k_LSF/'


with open(fromfile) as fp:
	for line in fp:
		#print line
		if 'inputList=HERE' in line :
			status, allFiles = commands.getstatusoutput('cmsLs ' + dir + ' | grep root | awk \'{print $5}\'')
			tofile.writelines('inputList = cms.untracked.vstring(')
			for line in allFiles.split():
				#print line
				tofile.writelines('\'root://eoscms.cern.ch//eos/cms'+line+'\',\n')
			tofile.writelines(') \n')
		else :
			tofile.writelines(line)


tofile.close()
