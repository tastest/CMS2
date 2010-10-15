#! /usr/bin/env python

import string
import commands, os
import sys 


if len(sys.argv) != 3:
    print 'Usage: '
    print './retrieveCrabLogFiles.py -c <crab directory of task>'
    sys.exit

crabDir = ''
for i in range(0, len(sys.argv)):
    if sys.argv[i] == '-c':
        crabDir = sys.argv[i+1]


if commands.getoutput('echo $CRABDIR')=='':
    print 'crab is not available. Please source the right env vars. Exiting'
    sys.exit()

if commands.getoutput('echo $CMSSW_BASE')=='':
    print 'Please do eval `scramv1 runtime -(c)sh` in a CMSSW area. Exiting'
    sys.exit()

if os.path.isdir(crabDir) == False:
    print  crabDir + ' is not a directory. Exiting'
    sys.exit()

#get the task id
taskId = commands.getoutput('crab -c ' + crabDir + ' -printId | grep \'Task Id\'').replace('Task Id','').replace('=','').replace(' ', '')
numberOfJobs = commands.getoutput('crab -c ' + crabDir + ' -status | grep \'Total Jobs\'').replace('Total Jobs','').replace('crab:','').replace(' ', '')


print 'Getting log files for ' + numberOfJobs + ' jobs in Task: ' + taskId

taskURL = 'http\://glidein-2.t2.ucsd.edu/CSstoragePath/' + taskId + '/'
for i in range(1, int(numberOfJobs)+1):
    fname = 'CMSSW_' + str(i) + '.stdout'
    cmd = 'wget ' + taskURL + '/' + fname + ' -O ' + crabDir + '/res/' + fname
    status = commands.getstatusoutput(cmd)
    if status[0] != 0:
        print 'Had a problem getting file ' + fname + ' from server. Check URL:'
        print taskURL + '/' + fname

    fname = 'CMSSW_' + str(i) + '.stderr'
    cmd = 'wget ' + taskURL + '/' + fname + ' -O ' + crabDir + '/res/' + fname
    status = commands.getstatusoutput(cmd)
    if status[0] != 0:
        print 'Had a problem getting file ' + fname + ' from server. Check URL:'
        print taskURL + '/' + fname

    fname = 'crab_fjr_' + str(i) + '.xml'
    cmd = 'wget ' + taskURL + '/' + fname + ' -O ' + crabDir + '/res/' + fname
    status = commands.getstatusoutput(cmd)
    if status[0] != 0:
        print 'Had a problem getting file ' + fname + ' from server. Check URL:'
        print taskURL + '/' + fname



        
    









