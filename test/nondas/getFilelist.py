#!/usr/bin/env python
import os,sys
import getopt
import commands

tofile = open('nonDBS_FileList.txt',"w")
#dir = '/store/user/rewang/AODSIM/HZZd/HZZd_MZd-60GeV_mH-126GeV_8TeV_D3_evts10k_LSF/'
#dir = '/store/user/rewang/AODSIM/HZZd/HZZd_MZd-40GeV_mH-126GeV_8TeV_D5_evts10k'
dir = '/store/user/rewang/AODSIM/HZZd/MC8TeV_HZZd40D5_FullSim_evts50k/'

count=0
status, allFiles = commands.getstatusoutput('cmsLs ' + dir + ' | grep root | awk \'{print $5}\'')
for line in allFiles.split():
	#print line
	count = count+1
	tofile.writelines('root://eoscms.cern.ch//eos/cms'+line+', #FILE'+str(count)+'_'+'\n')

print 'Total Files: '+str(count)

tofile.close()
