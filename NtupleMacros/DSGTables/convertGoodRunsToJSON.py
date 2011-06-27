#! /usr/bin/env python

import string, random
import commands, re, os
import sys 

if len(sys.argv) != 2:
    print 'Usage: ./convertGoodRunsToJSON.py'
    sys.exit()

infilelines = open(sys.argv[1], 'r').read().split('\n')

outlines = [""]

for j in infilelines:
    #line = j.replace("  ", "")
    line = j.split(" ")
    if len(line) < 3:
        continue
    run = line[0]
    foundrun = False
    outlindex = 0
    for index,item in enumerate(outlines):
        if item.find(run) != -1:
            foundrun = True
            outlindex = index
    if foundrun == False:
        outlines.append(("\"" + run+"\": [["+line[1] + ", " + line[2] +"]]").strip('\n'))
        continue
    if foundrun == True:
        #print j
        #print line[1]
        #print line[2]
        outlines[outlindex] = outlines[outlindex].replace("]]", "], ") + "[" + line[1] +"," + line[2] + "]]"

print "{",
for i in range(len(outlines)):
    if outlines[i]=="":
        continue
    if i < len(outlines)-1:
        print (outlines[i]+",").strip("\n"),
    else:
        print(outlines[i]).strip("\n"),

print "}"



        

