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

	myfinishjobs='0'
	mysubmitjobs='0'
	myrunjobs='0'
	mycreatedjobs='0'
    	for line in fp:

	  #if "Jobs with Wrapper Exit Code :" in line:
	  #    columns = line.split()
	  #    if columns[8] != "0" :
	  #  	SCRIPT.writelines('#' + line)
          #  	#print line
	  #  	nextline = fp.next()
	  #  	SCRIPT.writelines('#' + nextline + '\n')
	  #  	NLcolumns = nextline.split()
	  #  	#print NLcolumns[3]
	  #  	SCRIPT.writelines('multicrab -get '      + NLcolumns[3] + ' -c ' + ajob + ';\n')
	  #  	SCRIPT.writelines('multicrab -resubmit ' + NLcolumns[3] + ' -c ' + ajob + ';\n')

	  #if "Log file" in line:
	  #    logfile = line.split()
	  #    llfile = logfile[3].split('/')
	  #    print llfile[16]

	  if "Jobs with Wrapper Exit Code : 0" in line:
	      #print line
	      finishjobs = line.split()
	      myfinishjobs = finishjobs[1]
	      #print finishjobs[1]

	  if "Total Jobs" in line:
	       #print line
	       totaljobs = line.split()
	       #print totaljobs[1]
	  #    SCRIPT.writelines('#' + line + '\n')
	  #    SCRIPT.writelines('multicrab -submit  -c ' + ajob + ';\n')

	  if "Jobs Submitted" in line:
	       submitjobs = line.split()
	       mysubmitjobs = submitjobs[1]

	  if "Jobs Running" in line:
	       runjobs = line.split()
	       myrunjobs = runjobs[1]

	  if "Jobs Created" in line:
	       createdjobs = line.split()
	       mycreatedjobs = createdjobs[1]

	print '# Total: '+totaljobs[1]+' >>> Finished: '+myfinishjobs+' >>> Submitted: '+mysubmitjobs+' >>> Running: '+myrunjobs+' >>> Created: '+mycreatedjobs+'\n'
	SCRIPT.writelines('# Total: '+totaljobs[1]+' >>> Finished: '+myfinishjobs+' >>> Submitted: '+mysubmitjobs+' >>> Running: '+myrunjobs+' >>> Created: '+mycreatedjobs+'\n')

	if int(totaljobs[1]) == int(mysubmitjobs) : continue
	if mycreatedjobs == '0' : continue

	if int(totaljobs[1]) < 500 :
	  SCRIPT.writelines('multicrab -submit 1-'+ totaljobs[1] +' -c ' + ajob + ';\n')

        elif int(totaljobs[1]) < 1000 :
	  SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 500-'+ totaljobs[1] +' -c ' + ajob + ';\n')

	elif int(totaljobs[1]) < 1500 :
	  SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 500-999 -c ' + ajob + ';\n')
	  if(int(totaljobs[1])!= 1000 ):
	  	SCRIPT.writelines('multicrab -submit 1000-'+ totaljobs[1] +' -c ' + ajob + ';\n')
	  else :
		SCRIPT.writelines('multicrab -submit 1000 -c ' + ajob + ';\n')
	
	elif int(totaljobs[1]) < 2000 :
          SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 500-999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1000-1499 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 1500-'+ totaljobs[1] +' -c ' + ajob + ';\n')

	elif int(totaljobs[1]) < 2500 :
          SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 500-999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1000-1499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1500-1999 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 2000-'+ totaljobs[1] +' -c ' + ajob + ';\n')

        elif int(totaljobs[1]) < 3000 :
          SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 500-999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1000-1499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1500-1999 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 2000-2499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 2500-'+ totaljobs[1] +' -c ' + ajob + ';\n')

        elif int(totaljobs[1]) < 3500 :
          SCRIPT.writelines('multicrab -submit 1-499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 500-999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1000-1499 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 1500-1999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 2000-2499 -c ' + ajob + ';\n')
	  SCRIPT.writelines('multicrab -submit 2500-2999 -c ' + ajob + ';\n')
          SCRIPT.writelines('multicrab -submit 3000-'+ totaljobs[1] +' -c ' + ajob + ';\n')


	else :
	  print 'need to fix the code'


SCRIPT.close()

#os.system('more script_completedjobs.sh')

