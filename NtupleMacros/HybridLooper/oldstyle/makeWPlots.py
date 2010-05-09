#!/usr/bin/env python

import sys
import optparse
import commands
import os
import glob
    
    
    
    
#######################
# Get options
#######################

parser = optparse.OptionParser("usage: %prog [options]\
<input directory> \n")
parser.add_option ('--o', dest='outdir', type='string',
                   default = './',
                   help="directory for output png")
parser.add_option ('--t', dest='title', type='string',
                   default = 'WFinder Selected Leptons',
                   help="Title of html page")
options, args = parser.parse_args()

if len(args) != 1:
    print "Please specify output dir.  Exiting."
    print "./makeWPlots.py outputdir"
    sys.exit()

indir  = args[0]+"/"
outdir = options.outdir+"/"
title = options.title


##############################################
# Check dir
##############################################
if not os.path.isdir(indir) :
    print "Cannot find %s.  Exiting." % infile
    sys.exit()

##############################################
# Run doData.C and makePlots.C script
##############################################

doData = "root -l -q -b doData.C"
#os.system(doData)

rmplotfile = "rm -f results/*.png; rm -f %s/*.png " % (indir)
os.system(rmplotfile)

makePlots = "root -l -q -b makePlots.C"
os.system(makePlots)


##############################################
# Copy files to the WEBPAGE
##############################################

# copy the files
cpplotfile = "cp results/histos_data_lin*_selected*.png results/histos_data*lin*_nm1*.png results/histos_data_lin*_antiselected*.png results/*FO*png %s" % (indir)
os.system(cpplotfile)

cpindexfile = "cp wplots.html %sindex.html" % (indir)
os.system(cpindexfile)



files = glob.glob(indir+'*png')

for file in files:
    base = os.path.splitext(file)[0]
    print file, base
    rmcmd = "rm -f %s_small.png" % (base)
    os.system(rmcmd)
    command = "convert %s.png -resize 250x250 %s_small.png" % (base, base)
    os.system(command)

