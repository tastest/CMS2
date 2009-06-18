#! /usr/bin/env python

import string
import commands, os, re
import sys 
                      

cmsswSkelFile = ''
dataSet = ''
numEvtsTotal = -1
numEvtsPerJob = 5000
outNtupleName = 'ntuple.root'
storageElement = 'T2_US_UCSD'
tag = 'V01-02-06'
mode = 'glite'

def makeCrabConfig():
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2]
    outFile = open(outFileName + '.cfg', 'w')
    print 'Writing CRAB config file: ' + outFileName + '.cfg'
    outFile.write('[CRAB]\n')
    outFile.write('jobtype   = cmssw\n')
    outFile.write('scheduler = ' + mode + '\n\n')
    outFile.write('[CMSSW]\n')
    outFile.write('datasetpath            = ' + dataSet + '\n')
    outFile.write('pset                   = ' + outFileName + '_cfg.py \n')
    outFile.write('total_number_of_events = ' + str(numEvtsTotal) + '\n')
    outFile.write('events_per_job         = ' + str(numEvtsPerJob) + '\n')
    outFile.write('output_file            = ' + outNtupleName + '\n\n\n')
    outFile.write('[USER]\n')
    outFile.write('return_data            = 0\n')
    outFile.write('copy_data              = 1\n')
    outFile.write('storage_element        = ' + storageElement + '\n')
    outFile.write('ui_working_dir         = ' + outFileName + '\n')
    outFile.write('publish_data          = 1\n')
    outFile.write('publish_data_name     = CMS_' + tag + '\n')
    outFile.write('dbs_url_for_publication =  http://ming.ucsd.edu:8080/DBS1/servlet/DBSServlet\n\n')
    if ( mode == 'glite' ) :
        outFile.write('[GRID]\n')
        outFile.write('##Changing max tarball size to 20MB\n')
	outFile.write('maxtarballsize    =  20\n')
        outFile.write('##here are some default sites that we \n')
        outFile.write('##run at. Comment/Uncomment at will\n')
        outFile.write('##UCSD \n')
        outFile.write('CE_white_list = osg-gw-4.t2.ucsd.edu\n')
        outFile.write('SE_white_list = srm-3.t2.ucsd.edu\n')
        outFile.write('##WISC\n')
        outFile.write('#CE_white_list = cmsgrid02.hep.wisc.edu\n')
        outFile.write('#SE_white_list = cmssrm.hep.wisc.edu\n')
        outFile.write('##DESY\n')
        outFile.write('#CE_white_list = grid-ce3.desy.de\n')
        outFile.write('#SE_white_list = dcache-se-cms.desy.de\n')
        outFile.write('##Purdue\n')
        outFile.write('#CE_white_list = osg.rcac.purdue.edu \n')
        outFile.write('#SE_white_list = dcache.rcac.purdue.edu\n')
        outFile.write('##MIT\n')
        outFile.write('#CE_white_list = ce01.cmsaf.mit.edu\n')
        outFile.write('#SE_white_list = se01.cmsaf.mit.edu\n')
        outFile.write('##Nebraska\n')
        outFile.write('#CE_white_list = red.unl.edu,gpn-husker.unl.edu\n')
        outFile.write('#SE_white_list = srm.unl.edu\n')
        outFile.write('##IFCA\n')
        outFile.write('#CE_white_list = egeece01.ifca.es\n')
        outFile.write('#SE_white_list = storm.ifca.es\n')
        outFile.write('##Lyon\n')
        outFile.write('#CE_white_list = cclcgceli06.in2p3.fr or cclcgceli05.in2p3.fr \n')
        outFile.write('#SE_white_list = ccsrm.in2p3.fr\n')
        outFile.write('##CIEMAT\n')
        outFile.write('#CE_white_list = lcg02.ciemat.es\n')
        outFile.write('#SE_white_list = srm.ciemat.es\n')
        outFile.write('##IIHE\n')
        outFile.write('#CE_white_list = gridce.iihe.ac.be\n')
        outFile.write('#SE_white_list = maite.iihe.ac.be\n')
        outFile.write('##Aachen\n')
        outFile.write('#CE_white_list = grid-ce.physik.rwth-aachen.de\n')
        outFile.write('#SE_white_list = grid-srm.physik.rwth-aachen.de\n')

def makeCMSSWConfig(cmsswSkelFile):
    foundOutNtupleFile = False
    inFile = open(cmsswSkelFile, 'r').read().split('\n')
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2] + '_cfg.py'
    outFile = open(outFileName, 'w')
    print 'Writing CMS2 python config file : ' + outFileName
    for i in inFile:
        if i.find(outNtupleName) != -1:
            foundOutNtupleFile = True
        
        outFile.write(i+'\n')
        if i.find('process.load(\"CMS2.NtupleMaker.eventMaker_cfi\")') != -1:
            outFile.write('process.eventMaker.datasetName = cms.string(\"' +
                          dataSet+'\")\n')
            outFile.write('process.eventMaker.CMS2tag     = cms.string(\"' +
                          tag+'\")\n')
    if foundOutNtupleFile == False:
        print 'The root file you are outputting has not been found in your skeleton CMSSW config file'
        print 'Please check the name of the output root file in your PoolOutputModule, and try again'
        print 'Exiting!'
        sys.exit()
    outFile.close()

            
    
        


if len(sys.argv) < 5 :
    print 'Usage: makeCrabFiles.py [OPTIONS]'
    print 'Where the required options are: '
    print '\t-CMS2cfg name of the skeleton CMS2 config file '
    print '\t-d name of dataset'
    print '\t-t CMS2 tag, will be added to publish_data_name'
    print 'Optional arguments:'
    print '\t-strElem preferred storage element. Default is T2_US_UCSD if left unspecified'
    print '\t-nEvts Number of events you want to run on. Default is -1'
    print '\t-evtsPerJob Number of events per job. Default is 5000'
    print '\t-n Name of output Ntuple file. Default is ntuple.root'
    print '\t-m submission mode (possible: condor_g, condor, glite). Default is glite'
    sys.exit()


for i in range(0, len(sys.argv)):
    if sys.argv[i] == '-CMS2cfg':
        cmsswSkelFile = sys.argv[i+1]
    if sys.argv[i] == '-d':
        dataSet = sys.argv[i+1]
    if sys.argv[i] == '-nEvts':
        numEvtsTotal = sys.argv[i+1]
    if sys.argv[i] == '-evtsPerJob':
        numEvtsPerJob = sys.argv[i+1]
    if sys.argv[i] == '-strElem':
        storageElement = sys.argv[i+1]
    if sys.argv[i] == '-n':
        outNtupleName  = sys.argv[i+1]
    if sys.argv[i] == '-t':
        tag  = str(sys.argv[i+1])
    if sys.argv[i] == '-m':
        mode  = str(sys.argv[i+1])

##if os.path.exists(crabSkelFile) == False:
##    print 'Crab skeleton file does not exist. Exiting'
##    sys.exit()

if os.path.exists(cmsswSkelFile) == False:
    print 'CMSSW skeleton file does not exist. Exiting'
    sys.exit()

makeCMSSWConfig(cmsswSkelFile)
makeCrabConfig()
