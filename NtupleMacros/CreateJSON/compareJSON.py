#!/usr/bin/env python

import sys,os,json

def convertRange(runs):
    result = {}
    for run in runs.keys():
        tmp = []
        for lumirange in runs[run]:
            tmp.extend(range(lumirange[0],lumirange[1]+1))
        result[run] = tmp
    return result

one = json.load(open(sys.argv[1]))
two = json.load(open(sys.argv[2]))

sorted_runs_one = one.keys()
sorted_runs_one.sort()

sorted_runs_two = two.keys()
sorted_runs_two.sort()

runInBoth = []
runInNone = []

for run in sorted_runs_one:
    if run in sorted_runs_two:
        if run not in runInBoth:
            runInBoth.append(run)
    else :
        if run not in runInNone:
            runInNone.append(run)

for run in sorted_runs_two:
    if run in sorted_runs_one:
        if run not in runInBoth:
            runInBoth.append(run)
    else :
        if run not in runInNone:
            runInNone.append(run)
            
if len(runInNone) > 0 :
    print ''
    print 'Following runs are not in both files:'
    print ','.join(runInNone)
    print ''
else :
    print ''
    print 'Both files contain the same runs!'
    print ''
    
oneArray = convertRange(one)
twoArray = convertRange(two)

fileOneExcessLumi = {}
fileTwoExcessLumi = {}

for run in runInBoth:
    # print 'one',run,oneArray[run]
    # print 'two',run,twoArray[run]
    for lumi in oneArray[run] :
        if lumi not in twoArray[run] :
            if run in fileOneExcessLumi.keys():
                if lumi not in fileOneExcessLumi[run]:
                    fileOneExcessLumi[run].append(lumi)
            else :
                fileOneExcessLumi[run] = [lumi]
    for lumi in twoArray[run] :
        if lumi not in oneArray[run] :
            if run in fileTwoExcessLumi.keys():
                if lumi not in fileTwoExcessLumi[run]:
                    fileTwoExcessLumi[run].append(lumi)
            else :
                fileTwoExcessLumi[run] = [lumi]


fileOneExcessLumiSortedRuns = fileOneExcessLumi.keys()
fileOneExcessLumiSortedRuns.sort()                
if len(fileOneExcessLumiSortedRuns) > 0:
    print ''
    print 'File one contains these excess lumi sections:'
    for run in fileOneExcessLumiSortedRuns:
        print 'run:',run,'lumi sections:',fileOneExcessLumi[run]
    print ''
else:
    print ''
    print 'File one does not contain any excess lumi sections!'
    print ''
fileTwoExcessLumiSortedRuns = fileTwoExcessLumi.keys()
fileTwoExcessLumiSortedRuns.sort()                
if len(fileTwoExcessLumiSortedRuns) > 0:
    print ''
    print 'File two contains these excess lumi sections:'
    for run in fileTwoExcessLumiSortedRuns:
        print 'run:',run,'lumi sections:',fileTwoExcessLumi[run]
    print ''
else:
    print ''
    print 'File one does not contain any excess lumi sections!'
    print ''
    
                
