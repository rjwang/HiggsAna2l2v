#!/usr/bin/env python
import os,sys
import getopt
import commands


allsubmits = 'SCRIPT_Local.sh'
count=0
with open(allsubmits) as fp:

    SCRIPT = open('localscript_0.sh',"w")
    SCRIPT.writelines('#!bin/sh \n\n')
    SCRIPT.writelines('set CWD = `pwd`; \n')
    SCRIPT.writelines('cd $CMSSW_BASE/src/CMGTools/HiggsAna2l2v/; \n\n')

    for line in fp:
	if 'run2014_' in line:			
		SCRIPT.writelines(line+'\n')
	if 'sleep 25' in line:
		count = count + 1
		SCRIPT.writelines('cd $CWD;')
		SCRIPT.close()

		SCRIPT = open('localscript_'+str(count)+'.sh',"w")
		SCRIPT.writelines('#!bin/sh \n\n')
		SCRIPT.writelines('set CWD = `pwd`; \n')
		SCRIPT.writelines('cd $CMSSW_BASE/src/CMGTools/HiggsAna2l2v/; \n\n')

SCRIPT.writelines('cd $CWD;')
SCRIPT.close()
