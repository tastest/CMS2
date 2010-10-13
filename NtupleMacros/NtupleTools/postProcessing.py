#! /usr/bin/env python

import string, random
import commands, re, os
import sys 
import xml.dom.minidom
from xml.dom.minidom import Node

                      
################### Gets List of good XML files (files corresponding to jobs that have ###################
################### completed successfully                                             ###################

def getGoodXMLFiles(crabpath):
    global goodCrabXMLFiles
    # get list of Submission directories,
    #cmd = 'ls -dt ' + crabpath + '/res/Submission_*'
    #submissionDirs = []
    #temp = []
    #if commands.getstatusoutput(cmd)[0] != 256:
    #    submissionDirs = commands.getoutput(cmd).split('\n')
    #    cmd = 'cd ' + submissionDirs[0]+ '; ls *.xml' #start with the most recent directory
    #    temp = commands.getoutput(cmd).split('\n')
    #
    #
    #for i in temp:
    #    print 'Parsing ' + i
    #    try:
    #        doc = xml.dom.minidom.parse(submissionDirs[0]+'/'+i) #read xml file to see if the job failed
    #    except:
    #        print 'FrameworkJobReport:',i,'could not be parsed and is skipped'
	#    continue
    #    jobFailed = True
    #    for node in doc.getElementsByTagName("FrameworkJobReport"):
    #        key =  node.attributes.keys()[0].encode('ascii')
    #        if node.attributes[key].value == "Success":
    #            jobFailed = False
    #        if node.attributes[key].value == "Failed":
    #            print "Job " + i.split('/')[len(i.split('/'))-1].split('.')[0].split('_')[2] + " Failed!!!"
    #            jobFailed = True
    #    if jobFailed == False:
    #        goodCrabXMLFiles.append(submissionDirs[0] + '/' + i)
    #        
    #
    #    
    #for i in range(1,len(submissionDirs)):
    #    files = commands.getoutput('ls ' + submissionDirs[i] + '/crab_fjr_*.xml').split('\n')
    #    for j in files:
    #        duplicateFile = False
    #        for k in goodCrabXMLFiles:
    #            if k.split('/')[len(k.split('/'))-1] == j.split('/')[len(j.split('/'))-1]:
    #                duplicateFile = True
    #                break
    #        if duplicateFile == False:
    #            print 'Parsing ' + j
    #            try:
    #                doc = xml.dom.minidom.parse(j) #read xml file to see if the job failed
    #            except:
    #                print 'FrameworkJobReport:',i,'could not be parsed and is skipped'
    #                continue
    #            jobFailed = True
    #            for node in doc.getElementsByTagName("FrameworkJobReport"):
    #                key =  node.attributes.keys()[0].encode('ascii')
    #                if node.attributes[key].value == "Success":
    #                    jobFailed = False
    #                if node.attributes[key].value == "Failed":
    #                    print "Job " + j.split('/')[len(j.split('/'))-1].split('.')[0].split('_')[2] + " Failed!!!"
    #                    jobFailed = True
    #            if jobFailed == False:
    #                goodCrabXMLFiles.append(j)
    ##now add the files in the top res directory that do not
    #cmd = 'find ' + crabpath + '/res/ -iname *.xml'
	
    #instead of all above, just use /res/*.xml
    cmd = 'ls ' + crabpath + '/res/*.xml'
    temp = commands.getoutput(cmd).split('\n')
    
    for j in temp:
        print 'Parsing ' + j
        duplicateFile = False
        for k in goodCrabXMLFiles:
            if k.split('/')[len(k.split('/'))-1] == j.split('/')[len(j.split('/'))-1]:
                duplicateFile = True
                break
        if duplicateFile == False:
            try:
                doc = xml.dom.minidom.parse(j) #read xml file to see if the job failed
            except:
                print 'FrameworkJobReport:',j,'could not be parsed and is skipped'
                continue
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

def getNumEventsRun(crabpath):
    global totalNumEventsRun
    global goodCrabXMLFiles
    #get the job number from the crab file
    for i in goodCrabXMLFiles:
        print "New: " + i
        jobNum = i.split("_")[len(i.split("_"))-1]
        jobNum = jobNum.split(".")[0]
        cmd = "ls " +  outpath + "/preprocessing/ntuple_" + jobNum + "_1.root"
        print cmd
        if commands.getstatusoutput(cmd)[0] == 256:
            continue
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
                
def getGoodRootFiles(datapath,outpath):
    global goodCrabXMLFiles
    global goodRootFiles
    global CMSSWpath
    global dcachePrefix
    tempXMLFileList = []
    for i in goodCrabXMLFiles:
        #no more resubmission dir. Use <LFN> tag, never assume ntuple file name.
        #many lfn's, but first should always be right, and start with '/store/user'
        doc = xml.dom.minidom.parse(i)
        path = '/cms'+doc.getElementsByTagName("LFN")[0].firstChild.data.strip().rstrip()
        fname = path.split('/')[len(path.split('/'))-1]
        
        #old
        #resubmissionNum = '1'
        #if i.find('Submission')!=-1:
        #    for j in i.split('/'):
        #        if j.find('Submission') !=-1:
        #            resubmissionNum = str(int(j.split('_')[len(j.split('_'))-1]) + 1)
        #path = ''
        #fname = i.split('/')[len(i.split('/'))-1].replace('crab_fjr_', 'ntuple_').replace('.xml', '_'+resubmissionNum+'.root')

        #if commands.getstatusoutput('echo $HOSTNAME')[1].find('ucsd') !=-1:
            #path = datapath + fname #no longer need user path, use LFN above
        #el
        if commands.getstatusoutput('echo $HOSTNAME')[1].find('fnal') !=-1:
            path = datapath.replace('pnfs', 'usr', 1) + fname
        
        if datapath.find('pnfs') != -1:
            print 'Copying ' + fname + ' from dcache to ' + outpath + '/temp'
            cmd = 'dccp ' + dcachePrefix + path + ' ' + outpath + '/temp'
            print commands.getoutput(cmd)
        if datapath.find('hadoop') != -1:
            print 'Copying ' + fname + ' from hadoop to ' + outpath + '/temp'
            #cmd = 'cp ' + dcachePrefix + path + ' ' + outpath + '/temp'
            cmd = 'hadoop fs -copyToLocal ' + dcachePrefix + path + ' ' + outpath + '/temp'
            commands.getoutput(cmd)
        
        print 'Checking File ' + outpath + '/temp/' + fname + ' for integrity'
        cmd = ""
        if datapath.find("pnfs") != -1:
            cmd = "./sweepRoot -o Events " + outpath + '/temp/' + fname +  ' 2> /dev/null'
            print cmd
        else:
            cmd = "./sweepRoot -o Events " + outpath + '/temp/' + fname + ' 2> /dev/null'

        output = commands.getoutput(cmd).split('\n')
        #print commands.getstatusoutput(cmd)
        for k in output:
            if k.find('SUMMARY') != -1:
                print k
            if k.find('SUMMARY: 1 bad, 0 good') != -1:
                print 'File: ' + outpath + '/temp/' + fname + ' does not seem to be good!!!!\n'
                commands.getoutput('rm ' + outpath+'/temp/' + fname)
            elif k.find('SUMMARY: 0 bad, 1 good') != -1:
                print 'File: ' + outpath + '/temp/' + fname + ' looks just fine!\n'
                commands.getoutput('mv ' + outpath + '/temp/' + fname + ' ' + outpath + '/preprocessing/')
                goodRootFiles.append(outpath + '/preprocessing/' + fname)
                tempXMLFileList.append(i)
    
    goodCrabXMLFiles = tempXMLFileList


###########################################################################
## Get the skeleton macro and use it to make the postProcessing file
############################################################################
def makeRootMacros(outpath):
    global goodRootFiles
    global totalNumEventsRun
    #get the basic skeleton root script from my directory
    commands.getoutput("cp ~kalavase/crabTools/skelPostProcessingMacro.C postProcessingMacro.C")
    outFile = open("postProcessingMacro.C", "a")
    text = "\n\n\nvoid postProcess(Float_t xsec, Float_t kFactor, Float_t filterEfficiency) {\n"
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
    outFile.write("  chain->Merge(\"" + outpath + "merged_ntuple.root\",\"fast\");\n")
    outFile.write("\n}\n")
    outFile.close()
    cmd = "mv postProcessingMacro.C " + outpath
    print outpath + '/postProcessingMacro.C written'
    commands.getoutput(cmd)



              
    

if( len(sys.argv)!=7 ):
    print 'Usage: postProcessing.py -c [name of crab directory]'
    print '                              -d [directory where root files are]'
    print '                              -o [where you would like to put your final merged, dieted file]'
    sys.exit()

##global variables here
crabpath = ''
datapath = ''
outpath = ''
dcachePrefix = ''     

    
CMSSWpath = commands.getstatusoutput('echo $CMSSW_BASE')[1]
for i in range (0, len(sys.argv)):
    if(sys.argv[i] == "-c"):
        crabpath   = sys.argv[i+1] + "/"
    if(sys.argv[i] == "-d"):
        datapath    = sys.argv[i+1] + "/"
    if(sys.argv[i] == "-o"):
        outpath     = sys.argv[i+1] + "/"
                         

if( commands.getstatusoutput('ls ' + crabpath)[0] == 256):
    print 'The crab path does not exist. Please supply a valid path'
    sys.exit()

if( commands.getstatusoutput('ls ' + datapath)[0] == 256):
    print 'The directory containing the root files does not exist. Please supply a valid path'
    sys.exit()

if( commands.getstatusoutput('ls ' + outpath)[0] == 256):
    print '*******************************************'
    print 'The directory where you want your final root files to end up does not exist'
    sys.exit()

if( CMSSWpath ==''):
    print '$CMSSW_BASE not set. Please do eval `scramv1 runtime -sh` (or -csh for cshell) in the release you created the ntuples in'
    sys.exit()
##check if files in NtupleMacros exists
makefilepath = CMSSWpath + '/src/CMS2/NtupleMacros/NtupleTools/Makefile'
macropath    = CMSSWpath + '/src/CMS2/NtupleTools/NtupleTools/sweepRoot.C'
if( commands.getstatusoutput('ls ' + makefilepath)[0] == 256 or os.path.isfile('ls ' + macropath) == 256):
    print CMSSWpath + '/src/CMS2/NtupleMacros/NtupleTools/Makefile or ' + CMSSWpath + '/src/CMS2/NtupleMacros/NtupleTools/sweepRoot.C do not exist. Please check those files out from CVS and re-try'
    sys.exit()


#do some gymanstics depending on whether or not we're at ucsd
#or FNAL
if datapath.find("pnfs") != -1:
    uname = commands.getstatusoutput('echo $USER')[1]
    if( commands.getstatusoutput('echo $HOSTNAME')[1].find('ucsd') !=-1):
        publicDoors = [22136, 22137, 22138, 22139, 22140]
        if(uname == 'kalavase' or uname == 'fgolf'):
            publicDoors += [22162]
        #if at UCSD, get the personal dcache door
        dcachePrefix = 'dcap://dcap-2.t2.ucsd.edu:' + str(publicDoors[random.randrange(0,len(publicDoors), 1)]) + '/'
    if( commands.getstatusoutput('echo $HOSTNAME')[1].find('fnal') != -1):
        dcachePrefix = 'dcap://cmsdca.fnal.gov:24125/pnfs/fnal.gov/'
        

if datapath.find("hadooop") != -1:
        dcachePrefix = ''
    

if datapath.find("pnfs") != -1 or datapath.find("hadoop") != -1:
    print "Files are in the SE. Will transfer them to " + outpath + "/preprocessing to speed up dieting/merging/weighting step"
    cmd = "mkdir " + outpath + "/preprocessing"
    if commands.getstatusoutput(cmd)[0] == 256 and commands.getstatusoutput("ls " + outpath + "/preprocessing")[1]!="":
        print "The directory " + outpath + "/preprocessing already exists and is not empty!. Exiting!"
        sys.exit()
    cmd = "mkdir " + outpath + "/postprocessing"
    if commands.getstatusoutput(cmd)[0] == 256 and commands.getstatusoutput("ls " + outpath + "/postprocessing")[1] != "":
        print "The directory " + outpath + "/postprocessing already exists and is not empty!. Exiting!"
        sys.exit()
    cmd = "mkdir " + outpath + "/temp"
    if commands.getstatusoutput(cmd)[0] == 256 and commands.getstatusoutput("ls " + outpath + "/temp")[1] != "":
        print "The directory " + outpath + "/temp which will be a temp area for ntuples already exists and is not empty!. Exiting!"
        sys.exit()

        

##now make the sweeproot macro
cmd = 'cp ' + CMSSWpath + '/src/CMS2/NtupleMacros/NtupleTools/Makefile .'
commands.getstatusoutput(cmd)
cmd = 'cp ' + CMSSWpath + '/src/CMS2/NtupleMacros/NtupleTools/sweepRoot.C .'
commands.getstatusoutput(cmd)
cmd = 'make '
commands.getstatusoutput(cmd)

goodCrabXMLFiles = []
goodRootFiles = []
totalNumEventsRun = 0
rootFilesToMerge = []

getGoodXMLFiles(crabpath)
getGoodRootFiles(datapath,outpath)        
getNumEventsRun(crabpath)
makeRootMacros(outpath)
##copy crab directory to output directory
commands.getstatusoutput('mkdir ' + outpath + '/crab_logs')
commands.getstatusoutput('cp -r ' + crabpath + '* ' + outpath + 'crab_logs')


print '+++++++++++++++++++++++++++++'
print 'Total number of events that were run over to produce ntuples: ' + str(totalNumEventsRun)



    
