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

for fname in filelist:    
    if commands.getstatusoutput('ls ' +fname ) == 256:
        print 'The file does not exist'
        sys.exist()
    
    print 'Reading File ' + fname + ':'
    infile = open(fname, "r")
    lines = infile.readlines()
    infile.close()

    ##check if this is an EDProducer
    isEDProducer = False
    for line in lines:
        if line.replace(' ','').find("::produce") != -1:
            isEDProducer = True
            break
        
    if isEDProducer == False:
        continue
    
    b_declared=[]
    b_put = []
    fileIsGood = True
    #get the branch names that are declared to go into
    #the event
    for line in lines:
        if line.find("produces") == -1:
            continue;
        temp = line.replace(' ','')
        if temp.startswith('//'):
            continue
        ##blah is the name of the branch
        blah = line.split('\"')
        if len(blah) < 4:
            continue
        b_declared.append(blah[1])
        ##check that the alias checks out
        if line.find("setBranchAlias") == -1:
            continue;
        if blah[1] != blah[3].replace('_', ''):
            print "Alias name and Branch Name don't match!"
            print "Branch Name is: " + blah[1] + ", and Alias is: " + blah[3]
            fileIsGood = False
    
            

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
            print "Branch " + b_declared[branch] + "  has been declared, but it has not been put into the Event!"
            fileIsGood = False
    
    for branch in range(0, len(b_put)):
        if b_declared.count(b_put[branch]) == 0:
            print "Branch " + b_put[branch] + " has been put into the event, but not declared"
            fileIsGood = False
    if(fileIsGood):
        print fname + ' checks out\n'
    else:
        print fname + ' Has a problem!\n'
        ##sys.exit()

    




