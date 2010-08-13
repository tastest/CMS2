#!/usr/bin/env python

import sys,os,json

f = open(sys.argv[1])
runs = json.load(f)

sorted_runs = runs.keys()
sorted_runs.sort()

for run in sorted_runs:
    lumis = 0
    for range in runs[run]:
        lumis += len(range)
    print 'run:',run,'lumis:',lumis