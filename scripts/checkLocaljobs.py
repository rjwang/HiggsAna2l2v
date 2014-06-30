#!/usr/bin/env python
import os,sys
import getopt
import commands


SCRIPT = open('script_resubmit2local.sh',"w")

status, output = commands.getstatusoutput('ps -fww')
alljobs = output.split('\n')

SCRIPT.writelines('#!bin/sh \n\n')
SCRIPT.writelines('set CWD = `pwd`; \n')
SCRIPT.writelines('cd $CMSSW_BASE/src/CMGTools/HiggsAna2l2v/; \n\n')

for line in alljobs:
	if 'run2014_dijetAnalysis' in line:
		print line
		ll=line.split()
		SCRIPT.writelines(ll[7]+' '+ll[8]+' & ; \n')



SCRIPT.close()



