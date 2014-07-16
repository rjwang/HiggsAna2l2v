#!/usr/bin/env python
import os,sys
import getopt
import commands

SCRIPT = open('script_servertatus.sh',"w")
SCRIPT.writelines('#!bin/sh\n')
SERVER=['submit-4.t2.ucsd.edu','hcc-crabserver.unl.edu','submit-6.t2.ucsd.edu','vocms83.cern.ch','vocms20.cern.ch']

for ser in SERVER:
	#print ser
	status, output = commands.getstatusoutput('gsissh '+ser+' "du -sm *|sort -n"')
	if 'No such file or directory' in output: continue
	print '\n\033[91m'+output+'\033[0m'
	alljobs = output.split('\n')
	totalsize=0
	SCRIPT.writelines('\n')

	for ajob in alljobs:
		job = ajob.split('\t')
		size = int(job[0])
		totalsize = size + totalsize
		name = job[1]
		#SCRIPT.writelines('# '+name+': \t\t' + str(size/1000.) + " GB \n")
		SCRIPT.writelines('# gsissh '+ser+' "rm -rf '+name+'" ;\n')


	ratio=totalsize/1000.
	if ratio>50 :
		print '\n\033[35m>>> '+ser+' Total used size: '+str(totalsize/1000.)+'/(out of 100) GB <<<\033[0m \n'
	else :
		print '\n\033[92m>>> '+ser+' Total used size: '+str(totalsize/1000.)+'/(out of 100) GB <<<\033[0m \n'

	SCRIPT.writelines('\n# '+ser+' Total used size: '+str(totalsize/1000.)+'/(out of 100) GB \n')

SCRIPT.close()

#os.system('sh script_servertatus.sh')

