#!/usr/bin/env python
import os,sys
import getopt
import commands

dir='-1'

def help() :
   print 'changeXML.py -c [dir]'

#parse the options
try:
   # retrive command line options
   shortopts  = "c:?"
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
   elif o in('-c'): dir = a


if dir == '-1':
   help()
   sys.exit(1)


isFileExist = os.path.isfile(dir+'/share/arguments.xml.old')
if isFileExist : 
	print 'File: '+dir+'/share/arguments.xml.old'+' Exist!'
else : 
	os.system('cp '+dir+'/share/arguments.xml '+dir+'/share/arguments.xml.old')



fromfile = dir+'/share/arguments.xml.old'
#AllList = 'nonDBS_FileList.txt'
tofile = open(dir+'/share/arguments.xml',"w")

## change me
####################
NEvents=50000
MaxEvents=50
####################
print '++++++++++++++++++++++++++++++++++++++++++++'
print '  Please make sure to change the value:     '
print '	 NEvents: '+str(NEvents)
print '  MaxEvents: '+str(MaxEvents)
print '++++++++++++++++++++++++++++++++++++++++++++'


with open(fromfile) as fp:
	skipEvents=0
	njob=0
	for line in fp:
		#print line
		if '<Job' in line :
			jobids = line.split()
			lumiline = jobids[2]
			newline = jobids[3]
			jobid = jobids[3].split('"')[1]

			line = line.replace(lumiline,'SkipEvents="'+str(skipEvents)+'"')
			tofile.writelines(line)
			njob = njob+1
			skipEvents = MaxEvents*njob

		else :
			tofile.writelines(line)


tofile.close()
