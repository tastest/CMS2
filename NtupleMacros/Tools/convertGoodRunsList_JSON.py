#! /usr/bin/env python

# usage to print out good run selection for goodrun.cc:
# convertGoodRunsList_JSON.py <json file>

import sys,json

runs = json.load(open(sys.argv[1],"r"))

for run in runs.keys():
    for lumiBlock in runs[run] :
        if len(lumiBlock) != 2 :
            print 'ERROR reading lumi block: run:',run,'lumiBlock:',lumiBlock
        else:
            # print 'run ' + str(run) + ', min lumi ' + str(lumiBlock[0]) + ', max lumi ' + str(lumiBlock[1])
            print run,lumiBlock[0],lumiBlock[1]

print ''

