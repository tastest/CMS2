

#! /usr/bin/python


# import xml.dom.minidom
import string
import subprocess
import sys, os
import argparse
import urllib2
import itertools

# print "bunch of stuff \n\
# here too."

# def checkLocalFiles(localFilePath):
    # p = subprocess.Popen(["ls", localfilepath ""], stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    # err, out = p.output()
    # print out
    # print err
    # print p.exitstatus()
    # cmd = 'ls ' + localFilePath + 'merged_ntuple*.xml'


#global vars 
successlist = {'fjr':[[],[],[],[]], 
               'stdout':[[],[],[],[]],
               'stderr':[[],[],[],[]],
               'stdlog':[[],[],[],[]],
               'outfiles':[[],[],[],[]],
               }
filepath = ''

# I borrowed this function from http://thoughtsbyclayg.blogspot.com/2008/10/parsing-list-of-numbers-in-python.html
def parseIntSet(nputstr=""):
  selection = set()
  invalid = set()
  # tokens are comma seperated values
  tokens = [x.strip() for x in nputstr.split(',')]
  for i in tokens:
     try:
        # typically tokens are plain old integers
        selection.add(int(i))
     except:
        # if not, then it might be a range
        try:
           token = [int(k.strip()) for k in i.split('-')]
           if len(token) > 1:
              token.sort()
              # we have items seperated by a dash
              # try to build a valid range
              first = token[0]
              last = token[len(token)-1]
              for x in range(first, last+1):
                 selection.add(x)
        except:
           # not an int and not a range...
           invalid.add(i)
  # Report invalid tokens before returning valid selection
  if len(invalid) > 0:         
      print "Invalid set: " + str(invalid)
      exit(1)
  return selection
# end parseIntSet


#Takes a url as input and checks to see that it is accessible
def isGoodURL(server, task = ''):
    try:
        f = urllib2.urlopen(urllib2.Request(server + task))
        goodLinkFound = True
    except:
        goodLinkFound = False
    return goodLinkFound

def getCrabLog(serveraddress,taskid,job):
    username = taskid.split("_")[0] # gets username appended to front of task
    # n = len(username)
    global filepath
    filepath = '%s/temp/' % (taskid[+len(username)+1:-7])
    ensure_dir(filepath) # makes sure path exists and if it does not, creates it.
    
    # print 'Attempting to retrieve files for job %s to store in %s' % (str(job), filepath)     
    # urlforfile = '%s%s/crab_fjr_%s.xml' % (serveraddress, taskid, str(job))

    try:
      crabfjr = urllib2.urlopen('%s/%s/crab_fjr_%s.xml' % (serveraddress, taskid, str(job)))
      output = open('%scrab_fjr_%s.xml' % (filepath, str(job)),'wb')
      output.write(crabfjr.read())
      output.close()
      successlist['fjr'][0].append(job)
    except:
      successlist['fjr'][1].append(job)

    try:
      CMSSWstderr = urllib2.urlopen('%s/%s/CMSSW_%s.stderr' % (serveraddress, taskid, str(job)))
      output = open('%sCMSSW_%s.stderr' % (filepath, str(job)),'wb')
      output.write(CMSSWstderr.read())
      output.close()
      successlist['stderr'][0].append(job)
    except:
      successlist['stderr'][1].append(job)

    try:
      CMSSWlog = urllib2.urlopen('%s/%s/CMSSW_%s.log' % (serveraddress, taskid, str(job)))
      output = open('%sCMSSW_%s.log' % (filepath, str(job)),'wb')
      output.write(CMSSWlog.read())
      output.close()
      successlist['stdlog'][0].append(job)
    except:
      successlist['stdlog'][1].append(job)

    try:
      CMSSWstdout = urllib2.urlopen('%s/%s/CMSSW_%s.stdout' % (serveraddress, taskid, str(job)))
      output = open('%sCMSSW_%s.stdout' % (filepath, str(job)),'wb')
      output.write(CMSSWstdout.read())
      output.close()
      successlist['stdout'][0].append(job)
    except:
      successlist['stdout'][1].append(job)

    try:
      outfile = urllib2.urlopen('%s/%s/out_files_%s.tgz' % (serveraddress, taskid, str(job)))
      output = open('%sout_files_%s.tgz' % (filepath, str(job)),'wb')
      output.write(outfile.read())
      output.close()
      successlist['outfiles'][0].append(job)
    except:
      successlist['outfiles'][1].append(job)

def getFiles(jobs, serveradd, taskID):
  # print jobs
    for jobid in jobs:
      # print 'Getting files for job %i' %(jobid)
      getCrabLog(serveradd, taskID, jobid)
    

def ensure_dir(f):
    d = os.path.dirname(f)
    # Print 'Crab files will be stored in: %s' %(f) 
    if not os.path.exists(d):
      print 'Directory %s does not exist. Making directory.' %(f)
      os.makedirs(d)

def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
      b = list(b)
      yield b[0][1], b[-1][1]

####ryan's code to parse input options
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# parser = OptionParser()
# parser.add_option("--nev"   , dest="num_events"  , default=-1          , help="The number of events to run (-1 for all)")
# parser.add_option("--lumi"  , dest="lumi"        , default=1.0         , help="luminosity in 1 fb^-1"                   )
# parser.add_option("--print" , dest="print_output", default="0"         , help="Print the histograms"                    )
# parser.add_option("--suffix", dest="suffix"      , default="png"       , help="The suffix for the histograms"           )
# parser.add_option("--data"  , dest="data"        , action="store_true" , default=False, help="is this for data?"        )
# parser.add_option("--submit", dest="submit"      , action="store_true" , default=False, help="submit the commands"      )
#
# (options, args) = parser.parse_args()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def main():
    #this parses input arguments
    parser = argparse.ArgumentParser(
        prog='getCrabLogs.py',
        description='This script manually gets files from crab servers using python. This tool is not the best way to get Crab logs; \
It is for the user who is stuck on the crab command line when crab is not responding. When using crab -status, every task \
sent to crab has a task name of the format username_Dataset_Name_AndOther_Things_randomStringOfDigits. This must be input by the user in \
its full form for this script to work. Also, It is best to check on the crab task page to see which jobs are complete and \
attempting to get the fjrs from those jobs with crab -get before using this script. This script will not tell you if a the files from a \
specific job do not exist, but it will give you an error if no file is found.')
    #add arguments
    parser.add_argument('-t', required=True, 
                        help='Required. Input the crab task name. this must include the string of 6 random digits tagged on at the end of the task name.')
    parser.add_argument('-s', #nargs='+', 
                        help='Optional. Input the server or servers where your crab task is running. \
                              Usage -s server1,server2 ...')
    parser.add_argument('-j', required=True, 
                        help='Required. Input the list of job numbers you need to retrieve. \
                              Usage -j 1,3-5,7,10-15 ...')
    #parse arguments
    args = parser.parse_args()

    print '\nAttempting to get crab logs for jobs: %s \n    Task: %s\n' %(args.j, args.t) 

    #makes a set of ints for the jobs to be looped over from the input of the same style of output that crab spits out.
    selection = parseIntSet(args.j)

    #makes a default list of servers
    servers = ['http://glidein-2.t2.ucsd.edu/CSstoragePath/','http://submit-2.t2.ucsd.edu/CSstoragePath/','http://submit-3.t2.ucsd.edu/CSstoragePath/','http://submit-4.t2.ucsd.edu/CSstoragePath/']
    #Adds user input servers to default list
    #make list unique
    if args.s != None:
        print '\nAdding servers to list of known servers... \n'
        servers = args.s.split(',')+servers

    #this finds the correct server where your job is located
    for server in servers:

        site = isGoodURL(server)
        if site == True:
            taskpage = isGoodURL(server, args.t)
            if taskpage == True:
                print 'Task found on server: %s Attempting to get crab logs... \n' % (server)
                #loop here
                getFiles(selection, server, args.t)
                break
            else:
                print 'Task not found on server: %s Checking next server... \n' % (server)
        else:
            print 'Server: %s not responding. Checking next server... \n' % (server)

    #Weird workaround I had to implement to make the output look the same as the userinput
    for key in successlist:
      for i in range(len(list(ranges(successlist[key][0])))):
        int1, int2 = list(ranges(successlist[key][0]))[i]
        if (int2-int1 != 0):
          successlist[key][2].append(str(int1)+'-'+str(int2))
        else:
          successlist[key][2].append(str(int1))
          
    for key in successlist:
      for i in range(len(list(ranges(successlist[key][1])))):
        int1, int2 = list(ranges(successlist[key][1]))[i]
        if (int2-int1 != 0):
          successlist[key][3].append(str(int1)+'-'+str(int2))
        else:
          successlist[key][3].append(str(int1))
    
    print '\nCompleted retrieval of crab files.\n'

    #here is where all the summary stuff is listed
    print '\nsummary:\n'
    print 'Task name: %s' %(args.t)
    print 'Output Directory: %s' %(filepath)
    print 'Server: %s' %(server)
    
    print '\nsuccessful jobs:'

    for key in successlist:
      print '%s successfully downloaded for jobs: %s' %(key, (','.join(map(str,successlist[key][2]))))

    print '\nunsuccessful jobs:'
    for key in successlist:
      print '%s failed to download for jobs: %s' %(key,(','.join(map(str,(successlist[key][3])))))
    print

if __name__ == "__main__":
    sys.exit(main())
        









