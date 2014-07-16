#!/usr/bin/env python
import os,sys
import getopt
import commands

Njob='-1'

def help() :
   print 'generateSingleFileCfg.py -N give the nth input file'

#parse the options
try:
   # retrive command line options
   shortopts  = "N:?"
   opts, args = getopt.getopt( sys.argv[1:], shortopts )
except getopt.GetoptError:
   # print help information and exit:
   print "ERROR: unknown options in argument %s" % sys.argv[1:]
   help()
   sys.exit(1)

for o,a in opts:
   if o in("-?", "-h"):
      help()
      sys.exit(1)
   elif o in('-N'): Njob = a


if Njob == '-1':
   help()
   sys.exit(1)


fromfile = 'run2l2vFullAnalysis_mc_template_cfg.py'
AllList = 'nonDBS_FileList.txt'
tofile = open('mc_cfg.py',"w")


with open(fromfile) as fp:
	for line in fp:
		#print line
		if 'inputList=HERE' in line :

			with open(AllList) as ffp:
				for lline in ffp:
					if 'FILE'+Njob+'_' in lline :
						jobFile=lline.split(',')
						tofile.writelines('inputList = cms.untracked.vstring(')
						tofile.writelines('\''+jobFile[0]+'\',\n')
						tofile.writelines(') \n')
		else :
			tofile.writelines(line)


tofile.close()
