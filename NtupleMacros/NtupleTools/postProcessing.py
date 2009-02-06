#! /usr/bin/env python

import string
import commands, re
import sys 
import xml.dom.minidom
from xml.dom.minidom import Node

                      
################### Gets List of good XML files (files corresponding to jobs that have ###################
################### completed successfully                                             ###################

def getGoodXMLFiles(crabpath):
    global goodCrabXMLFiles
    # get list of Submission directories,
    cmd = 'ls -dt ' + crabpath + '/res/Submission_*'
    submissionDirs = []
    temp = []
    if commands.getstatusoutput(cmd)[0] != 256:
        submissionDirs = commands.getoutput(cmd).split('\n')
        cmd = 'cd ' + submissionDirs[0]+ '; ls *.xml' #start with the most recent directory
        temp = commands.getoutput(cmd).split('\n')
    
    
    for i in temp:
        print 'Parsing ' + i
        doc = xml.dom.minidom.parse(submissionDirs[0]+'/'+i) #read xml file to see if the job failed
        jobFailed = True
        for node in doc.getElementsByTagName("FrameworkJobReport"):
            key =  node.attributes.keys()[0].encode('ascii')
            if node.attributes[key].value == "Success":
                jobFailed = False
            if node.attributes[key].value == "Failed":
                print "Job " + i.split('/')[len(i.split('/'))-1].split('.')[0].split('_')[2] + " Failed!!!"
                jobFailed = True
        if jobFailed == False:
            goodCrabXMLFiles.append(submissionDirs[0] + '/' + i)
            
    
        
    for i in range(1,len(submissionDirs)):
        files = commands.getoutput('ls ' + submissionDirs[i] + '/crab_fjr_*.xml').split('\n')
        for j in files:
            duplicateFile = False
            for k in goodCrabXMLFiles:
                if k.split('/')[len(k.split('/'))-1] == j.split('/')[len(j.split('/'))-1]:
                    duplicateFile = True
                    break
            if duplicateFile == False:
                print 'Parsing + ' j
                doc = xml.dom.minidom.parse(j) #read xml file to see if the job failed
                jobFailed = True
                for node in doc.getElementsByTagName("FrameworkJobReport"):
                    key =  node.attributes.keys()[0].encode('ascii')
                    if node.attributes[key].value == "Success":
                        jobFailed = False
                    if node.attributes[key].value == "Failed":
                        print "Job " + j.split('/')[len(j.split('/'))-1].split('.')[0].split('_')[2] + " Failed!!!"
                        jobFailed = True
                if jobFailed == False:
                    goodCrabXMLFiles.append(j)
    
            
    ##now add the files in the top res directory that do not
    cmd = 'find ' + crabpath + '/res/ -iname *.xml'
    temp = commands.getoutput(cmd).split('\n')

    
    
    for j in temp:
        print 'Parsing ' + j
        duplicateFile = False
        for k in goodCrabXMLFiles:
            if k.split('/')[len(k.split('/'))-1] == j.split('/')[len(j.split('/'))-1]:
                duplicateFile = True
                break
        if duplicateFile == False:
            doc = xml.dom.minidom.parse(j) #read xml file to see if the job failed
            jobFailed = False
            for node in doc.getElementsByTagName("FrameworkJobReport"):
                key =  node.attributes.keys()[0].encode('ascii')
                if node.attributes[key].value == "Success":
                    jobFailed = False
                if node.attributes[key].value == "Failed":
                    print "Job " + j.split('/')[len(j.split('/'))-1].split('.')[0].split('_')[2] + " Failed!!!"
                    jobFailed = True                
            if jobFailed == False:
                goodCrabXMLFiles.append(j)


###############################################################################################################
################### Get Number of Events Run ###################

##eventsread = 0
##for i in range(0,len(xmlfiles)):
##    crabxmlfile = xmlfiles[i]
##    print 'reading '+crabxmlfile
##    #get the number of events processed
##    #parse the xml file
##    doc = xml.dom.minidom.parse(crabxmlfile)
##    for node in doc.getElementsByTagName("EventsRead"):
##        s = node.firstChild.data
##        #s is in $!#%^& unicode
##        #don't need to encode in ascii, but runs faster
##        s.encode('ascii') 
##        eventsread = eventsread+int(s)
        
##print str(eventsread)+' were processed \n'


def getNumEventsRun(crabpath):
    global totalNumEventsRun
    global goodCrabXMLFiles
    #get the job number from the crab file
    for i in goodCrabXMLFiles:
        #parse the crab file:
        print "Getting Number of Events run for Job: " + i 
        doc = xml.dom.minidom.parse(i)
        for node in doc.getElementsByTagName("EventsRead"):
            s = node.firstChild.data
            #s is in unicode
            #don't need to encode in ascii, but runs faster
            s.encode('ascii')
            totalNumEventsRun = totalNumEventsRun + int(s)

###############################################################################################################
################### Get List Of Root Files which are not corrupted and have the Events tree ###################
                
def getGoodRootFiles(datapath):
    global goodCrabXMLFiles
    global goodRootFiles
    #j = 0 # temp
    tempXMLFileList = []
    for i in goodCrabXMLFiles:
        #if j > 2:
            #break
        #j = j+1
        path = datapath + i.split('/')[len(i.split('/'))-1].replace('crab_fjr_', 'ntuple_').replace('.xml', '.root')
        print 'Checking File ' + path + ' for integrity'
        cmd = ""
        if path.find("pnfs") != -1:
            cmd = "/home/users/kalavase/crabTools/sweepRoot -o Events dcap://dcap-2.t2.ucsd.edu:22136/" + path + ' 2> /dev/null'
        else:
            cmd = "/home/users/kalavase/crabTools/sweepRoot -o Events " + path + ' 2> /dev/null'
        output = commands.getoutput(cmd).split('\n')
        for k in output:
            if k.find('SUMMARY') != -1:
                print k
            if k.find('SUMMARY: 1 bad, 0 good') != -1:
                print 'File: ' + path + 'does not seem to be good!!!!'
            elif k.find('SUMMARY: 0 bad, 1 good') != -1:
                print 'File: ' + path + ' looks just fine!'
                goodRootFiles.append(path)
                tempXMLFileList.append(i)
    
    goodCrabXMLFiles = tempXMLFileList


def makeRootMacros(outpath):
    global goodRootFiles
    global totalNumEventsRun
    #get the basic skeleton root script from my directory
    commands.getoutput("cp ~kalavase/crabTools/skelPostProcessingMacro.C postProcessingMacro.C")
    outFile = open("postProcessingMacro.C", "a")
    text = "\n\n\nvoid postProcess(Float_t xsec, Float_t kFactor, Float_t filterEfficiency=1) {\n"
    outFile.write(text)
    #write the 
    for i in goodRootFiles:
        fName = i.split('/')[len(i.split('/')) - 1]
        finPath  = "\"" + outpath + "preprocessing/" + fName + "\""
        foutPath = "\"" + outpath + "postprocessing/" + fName  + "\""
        temp = ", " + str(totalNumEventsRun)+ ", xsec, kFactor, filterEfficiency);\n"
        outFile.write("  doPostProcessing("+finPath + "," + foutPath+temp)
    outFile.write("\n}\n\n\n")
    outFile.write("void mergeAll() {\n")
    outFile.write("  TChain *chain = new TChain(\"Events\");\n\n")

    for i in goodRootFiles:
        fname = i.split("/")[len(i.split("/"))-1]
        outFile.write("  chain->Add(\"" + outpath + "postprocessing/" + fname+"\");\n")
    outFile.write("  chain->Merge(\"" + outpath + "postprocessing/merged_ntuple.root\",\"fast\");\n")
    outFile.write("\n}\n")
    outFile.close()
    cmd = "mv postProcessingMacro.C " + outpath
    commands.getoutput(cmd)



              
    

if( len(sys.argv)!=7 ):
    print 'Usage: postProcessing.py -c [name of crab directory]'
    print '                              -d [directory where root files are]'
    print '                              -o [where you would like to put your final merged, dieted file]'
    sys.exit()


crabpath = ''
datapath = ''
outpath = ''
for i in range (0, len(sys.argv)):
    if(sys.argv[i] == "-c"):
        crabpath   = sys.argv[i+1] + "/"
    if(sys.argv[i] == "-d"):
        datapath    = sys.argv[i+1] + "/"
    if(sys.argv[i] == "-o"):
        outpath     = sys.argv[i+1] + "/"

print 'Make sure that you have setup your CMSSW environment by doing eval `scramv1 runtime -sh` or eval `scramv1 runtime -csh in a CMSSW area!!!! \n\n'

if( commands.getstatusoutput('ls ' + crabpath)[0] == 256):
    print 'The crab path does not exist. Please supply a valid path'
    sys.exit()

if( commands.getstatusoutput('ls ' + datapath)[0] == 256):
    print 'The directory containing the root files does not exist. Please supply a valid path'
    sys.exit()

if( commands.getstatusoutput('ls ' + outpath)[0] == 256):
    print '*******************************************'
    print 'The directory where you want your final root files to end up does not exist'
    print 'Program will run, but keep in mind that ' + outpath + 'does not exist'


if datapath.find("pnfs") != -1:
    print "Files are on dcache. After some processing, will transfer them to " + outpath + "/preprocessing to speed up dieting/merging/weighting step"
    cmd = "mkdir " + outpath + "/preprocessing"
    if commands.getstatusoutput(cmd)[0] == 256 and commands.getstatusoutput("ls " + outpath + "/preprocessing")[1]!="":
        print "The directory " + outpath + "/preprocessing already exists and is not empty!. Exiting!"
        sys.exit()
    cmd = "mkdir " + outpath + "/postprocessing"
    if commands.getstatusoutput(cmd)[0] == 256 and commands.getstatusoutput("ls " + outpath + "/postprocessing")[1] != "":
        print "The directory " + outpath + "/postprocessing already exists and is not empty!. Exiting!"
        sys.exit()
    cmd = "voms-proxy-info"
    if commands.getstatusoutput(cmd)[0] == 256:
        print "Your environment is not set up correctly"
        print "Please get your voms-proxy and source the right CRAB scripts. Exiting!"
        sys.exit()

        

goodCrabXMLFiles = []
goodRootFiles = []
totalNumEventsRun = 0
rootFilesToMerge = []

getGoodXMLFiles(crabpath)
getGoodRootFiles(datapath)        
getNumEventsRun(crabpath)

print '+++++++++++++++++++++++++++++'
print 'Total number of events that were run over to produce ntuples: ' + str(totalNumEventsRun)

if datapath.find("pnfs") != -1:
    print "Moving files from dcache to " + outpath + "/preprocessing"
    for i in goodRootFiles:
        print i
        cmd = "dccp dcap://dcap-2.t2.ucsd.edu:22136/" + i + " " + outpath + "/preprocessing"
        print cmd
        print commands.getoutput(cmd)

    
print 'Now making root macros for dieting'
makeRootMacros(outpath)


    


    






    
    









    

    
    

