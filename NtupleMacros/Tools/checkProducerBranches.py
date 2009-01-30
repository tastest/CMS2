#!/usr/bin/env python
####################################################################
#Simple script to check that for every branch you have declared,
#you also put it into the event
#Usage:
####################################################################
#system imports
import sys, commands, string, os


if len(sys.argv) != 2:
    print 'Usage: checkCMS2Producer.py [name of CMS Producer to check or src directory]'
    sys.exit()

#filelist = commands.getoutput('pwd ' + sys.argv[1] + '| ls *.cc').split('\n')

cmd = ''
if os.path.isdir(sys.argv[1]):
    cmd = "ls " + sys.argv[1] + "/*.cc"
if os.path.isfile(sys.argv[1]):
    cmd = "ls " + sys.argv[1]

filelist = commands.getoutput(cmd).split('\n')


if len(filelist) == 0:
    print 'No .cc files in directory'
    sys.exit()


for f in range(0, len(filelist) - 1):
    fname = filelist[f]
    print 'Reading File ' + fname + ':'
    if commands.getstatusoutput('ls ' +fname ) == 256:
        print 'The file does not exist'
        sys.exist()
    

    infile = open(fname, "r")
    lines = infile.readlines()
    infile.close()
    
    b_declared=[]
    b_put = []
    fileIsGood = True
    #get the branch names that are declared to go into
    #the event
    for line in range(0,len(lines)):
        if lines[line].find("produces") == -1:
            continue;
        temp = lines[line].replace(' ','')
        if temp.startswith('//'):
            continue
        blah = lines[line].split('\"')
        b_declared.append(blah[1])
        #print blah[1]
    

    ##get the names of the branches put into the event
    for line in range(0,len(lines)):
        if lines[line].find(".put(") == -1:
            continue;
        temp = lines[line].replace(' ','')
        if temp.startswith('//'): #commented out lines
            continue
        blah = lines[line].split('\"')
        if len(blah) > 1:
            b_put.append(blah[1])
    

    for branch in range(0, len(b_declared)):
        if b_put.count(b_declared[branch]) == 0:
            print b_declared[branch] + " not put in the event"
            fileIsGood = False
    
    for branch in range(0, len(b_put)):
        if b_declared.count(b_put[branch]) == 0:
            print b_put[branch] + " put into the event, but not declared"
            fileIsGood = False
    print fname + "checks out!"
    print '\n\n'




