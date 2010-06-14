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
mode = 'glidein'
dbs_url = 'http://ming.ucsd.edu:8080/DBS2/servlet/DBSServlet';
report_every = 1000;
global_tag_flag = '';


def makeCrabConfig():
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2]
    outFile = open(outFileName + '.cfg', 'w')
    print 'Writing CRAB config file: ' + outFileName + '.cfg'
    outFile.write('[CRAB]\n')
    outFile.write('jobtype   = cmssw\n')
    outFile.write('scheduler = ' + mode + '\n')
    outFile.write('use_server = ' + '1' + '\n')
    outFile.write('\n[CMSSW]\n')
    outFile.write('datasetpath             = ' + dataSet + '\n')
    outFile.write('pset                    = ' + outFileName + '_cfg.py \n')
    outFile.write('total_number_of_events  = ' + str(numEvtsTotal) + '\n')
    outFile.write('events_per_job          = ' + str(numEvtsPerJob) + '\n')
    outFile.write('output_file             = ' + outNtupleName + '\n\n\n')
    outFile.write('[USER]\n')
    outFile.write('return_data             = 0\n')
    outFile.write('copy_data               = 1\n')
    outFile.write('storage_element         = ' + storageElement + '\n')
    outFile.write('ui_working_dir          = ' + outFileName + '\n')
    outFile.write('user_remote_dir         = CMS2_' + tag + '/' + outFileName + '\n')
    outFile.write('publish_data            = 0\n')
    outFile.write('publish_data_name       = CMS2_' + tag + '\n')
    outFile.write('dbs_url_for_publication = ' + dbs_url + '\n\n')
    outFile.write('[GRID]\n')
    outFile.write('maxtarballsize = 20\n')
    
	 
    outFile.write('##here are some default sites that we \n')
    outFile.write('##run at. Comment/Uncomment at will\n')
    outFile.write('##UCSD \n')
    outFile.write('#SE_white_list = T2_US_UCSD\n')
    outFile.write('##WISC\n')
    outFile.write('#SE_white_list = T2_US_Wisconsin\n')
    outFile.write('##DESY\n')
    outFile.write('#SE_white_list = T2_DE_DESY\n')
    outFile.write('##Purdue\n')
    outFile.write('#SE_white_list = T2_US_Purdue\n')
    outFile.write('##MIT\n')
    outFile.write('#SE_white_list = T2_US_MIT\n')
    outFile.write('##Nebraska\n')
    outFile.write('#SE_white_list = T2_US_Nebraska\n')
    outFile.write('##IFCA\n')
    outFile.write('#SE_white_list = T2_ES_IFCA\n')
    outFile.write('##Lyon\n')
    outFile.write('#SE_white_list = T2_FR_CCIN2P3\n')
    outFile.write('##CIEMAT\n')
    outFile.write('#SE_white_list = T2_ES_CIEMAT\n')
    outFile.write('##IIHE\n')
    outFile.write('#SE_white_list = T2_BE_IIHE\n')
    outFile.write('##Aachen\n')
    outFile.write('#SE_white_list = T2_DE_RWTH\n')
#
def makeCMSSWConfig(cmsswSkelFile):
    foundOutNtupleFile = False
    foundreportEvery = False
    inFile = open(cmsswSkelFile, 'r').read().split('\n')
    for i in inFile:
        if i.find(outNtupleName) != -1:
            foundOutNtupleFile = True
        if i.find('reportEvery') != -1:
            foundOutNtupleFile = True
    if foundOutNtupleFile == False:
        print 'The root file you are outputting is not named ntuple.root as it should be for a CMS2 job.'
        print 'Please check the name of the output root file in your PoolOutputModule, and try again'
        print 'Exiting!'
        sys.exit()
    outFileName = dataSet.split('/')[1]+'_'+dataSet.split('/')[2] + '_cfg.py'
    print 'Writing CMS2 CMSSW python config file : ' + outFileName
    outFile = open(outFileName, 'w')
    for i in inFile:

        if i.find('reportEvery') != -1:
            outFile.write('process.MessageLogger.cerr.FwkReport.reportEvery = ' + str(report_every) + '\n'); continue

        if i.find('globaltag') != -1:
            outFile.write('process.GlobalTag.globaltag = "' + global_tag + '"\n'); continue

        outFile.write(i+'\n')
        
        if i.find('cms.Path') != -1:
            outFile.write('process.eventMaker.datasetName = cms.string(\"' +
                          dataSet+'\")\n')
            outFile.write('process.eventMaker.CMS2tag     = cms.string(\"' +
                          tag+'\")\n')

    outFile.close()



       


if len(sys.argv) < 5 :
    print 'Usage: makeCrabFiles.py [OPTIONS]'
    print '\nWhere the required options are: '
    print '\t-CMS2cfg\tname of the skeleton CMS2 config file '
    print '\t-d\t\tname of dataset'
    print '\t-t\t\tCMS2 tag, will be added to publish_data_name'
    print '\nOptional arguments:'
    print '\t-strElem\tpreferred storage element. Default is T2_US_UCSD if left unspecified'
    print '\t-nEvts\t\tNumber of events you want to run on. Default is -1'
    print '\t-evtsPerJob\tNumber of events per job. Default is 5000'
    #print '\t-n\t\tName of output Ntuple file. Default is ntuple.root'
    print '\t-m\t\tsubmission mode (possible: condor_g, condor, glite). Default is glidein'
    print '\t-dbs\t\tdbs url for publication. Default is http://ming.ucsd.edu:8080/DBS2/servlet/DBSServlet'
    print '\t-re\t\tMessage Logger modulus for error reporting. Default is 1000'
    print '\t-gtag\t\tglobal tag. Default is MC_31X_V3::All'
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
    #if sys.argv[i] == '-n':
    #    outNtupleName  = sys.argv[i+1]
    if sys.argv[i] == '-t':
        tag  = str(sys.argv[i+1])
    if sys.argv[i] == '-m':
        mode  = str(sys.argv[i+1])
    if sys.argv[i] == '-dbs':
        dbs_url = str(sys.argv[i+1])
    if sys.argv[i] == '-re':
        report_every = str(sys.argv[i+1])
    if sys.argv[i] == '-gtag':
        global_tag_flag = str(sys.argv[i+1])

if os.path.exists(cmsswSkelFile) == False:
    print 'CMSSW skeleton file does not exist. Exiting'
    sys.exit()


#print '\nGetting global tag from DBS...'
if( global_tag_flag != '' ):
	print '\nUsing \'' + global_tag_flag + '\' specified by -gtag flag.\n'
	global_tag = global_tag_flag
	makeCMSSWConfig(cmsswSkelFile)
	makeCrabConfig()
else :
    global_tag = '';
    dbs_result = '';
    command = 'dbsql find config.name,config.content where dataset=' + dataSet + '>config.content; while read line; do globaltag=`echo $line | sed -n \'s/^.*process.GlobalTag.globaltag = \([^p]*\).*$/\\1/p\'`; if [ "$globaltag" != "" ]; then echo $globaltag; break; fi; done <config.content; rm config.content';
    len = len( os.popen(command).readlines() )
    if( len > 0 ):
      lines = os.popen(command);
      for i in lines.readlines():
        dbs_result = re.sub('\n', '', i)
        global_tag = re.sub('#.+$', '', dbs_result)
        if( global_tag != '' and global_tag_flag == '' ):
            print '\nDBS Query results:\t\'' + dbs_result + '\' ?\n'
            print 'Use global tag from DBS:\t\'' + global_tag + '\' ?\n'
            answer = raw_input('[y/n]?')
            while(answer != 'y' and answer != 'n'): 
                print 'Please pick either \'y\' or \'n\''
                answer = raw_input('[y/n]?')
            if answer == 'n':
                print 'Enter alternative Global Tag:'
                global_tag = raw_input('new global tag:')
            makeCMSSWConfig(cmsswSkelFile)
            makeCrabConfig()
    else: 
      print '\nGlobal tag not found in DBS. Use -gtag to set global tag. Exiting...\n'
      sys.exit()
