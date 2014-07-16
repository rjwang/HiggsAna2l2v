#!/usr/bin/env python
import os,sys
import getopt
import commands

SCRIPT = open('script_submitjobs.sh',"w")

status, output = commands.getstatusoutput('ls -d multicrab_*')
alljobs = output.split()
for ajob in alljobs:
    joblog = 'output_'+ajob+'.log'
    print 'running: '+joblog
    SCRIPT.writelines('\n###################' + '\n')
    SCRIPT.writelines('#' + joblog + '\n')
    SCRIPT.writelines('###################' + '\n\n')
    with open(joblog) as fp:

    	for line in fp:

	  if "Jobs with Wrapper Exit Code :" in line:
	      columns = line.split()
	      if columns[8] != "0" :
	    	SCRIPT.writelines('#' + line)
            	#print line
	    	nextline = fp.next()
	    	SCRIPT.writelines('#' + nextline + '\n')
	    	NLcolumns = nextline.split()
	    	#print NLcolumns[3]

	    	SCRIPT.writelines('multicrab -get '      + NLcolumns[3] + ' -c ' + ajob + ';\n')
	    	SCRIPT.writelines('multicrab -resubmit ' + NLcolumns[3] + ' -c ' + ajob + ';\n')

	  if "List of jobs Cancelled:" in line:
	      #print line
	      ccolumns = line.split()
	      #print ccolumns[4] 
	      SCRIPT.writelines('\n #' + line + '\n')
	      SCRIPT.writelines('multicrab -kill '     + ccolumns[4] + ' -c ' + ajob + ';\n')
	      SCRIPT.writelines('multicrab -resubmit ' + ccolumns[4] + ' -c ' + ajob + ';\n')

	  if "Jobs Aborted" in line:
	      #print line
	      nline = fp.next()
              nnline = fp.next()
	      #print nline
	      #print nnline
	      nncolumns = nnline.split()
	      #print nncolumns[3]
	      SCRIPT.writelines('\n #' + nline)
	      SCRIPT.writelines('#' + nnline + '\n')
	      SCRIPT.writelines('multicrab -resubmit ' + nncolumns[3] + ' -c ' + ajob + ';\n')

	  if "Jobs Created" in line:
	      #print line
	      SCRIPT.writelines('#' + line + '\n')
	      SCRIPT.writelines('multicrab -submit  -c ' + ajob + ';\n')


SCRIPT.close()

#os.system('more script_submitjobs.sh')

