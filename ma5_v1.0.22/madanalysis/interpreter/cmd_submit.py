################################################################################
#  
#  Copyright (C) 2012 Eric Conte, Benjamin Fuks, Guillaume Serret
#  The MadAnalysis development team, email: <ma5team@iphc.cnrs.fr>
#  
#  This file is part of MadAnalysis 5.
#  Official website: <http://madanalysis.irmp.ucl.ac.be>
#  
#  MadAnalysis 5 is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  MadAnalysis 5 is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with MadAnalysis 5. If not, see <http://www.gnu.org/licenses/>
#  
################################################################################


import madanalysis.interpreter.cmd_base as CmdBase
from   madanalysis.IOinterface.job_writer import JobWriter
from   madanalysis.IOinterface.job_reader import JobReader
import logging
import glob
import os
import time

class CmdSubmit(CmdBase.CmdBase):
    """Command SUBMIT"""

    def __init__(self,main,resubmit=False):
        self.resubmit=resubmit
        if not resubmit:
            CmdBase.CmdBase.__init__(self,main,"submit")
        else:
            CmdBase.CmdBase.__init__(self,main,"resubmit")
        self.forbiddenpaths=[]
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/lib'))
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/bin'))
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/madanalysis'))

    @staticmethod
    def chronometer_display(diff):
        fill=False
        theLast=time.localtime(diff)
        theLastStr=""
        if theLast.tm_mday>1:
            theLastStr=time.strftime('%d days ',theLast)
            fill=True
        elif theLast.tm_mday==0:
            theLastStr=time.strftime('%d day ',theLast)
        if theLast.tm_hour>1:
            theLastStr=time.strftime('%H hours ',theLast)
            fill=True
        elif theLast.tm_hour==0 or fill:
            theLastStr=time.strftime('%H hour ',theLast)
        if theLast.tm_min>1:
            theLastStr=time.strftime('%M minutes ',theLast)
            fill=True
        elif theLast.tm_min==0 or fill:
            theLastStr=time.strftime('%M minute ',theLast)
        if theLast.tm_sec>1:
            theLastStr=time.strftime('%s seconds ',theLast)
        else:
            theLastStr=time.strftime('%s second ',theLast)
        return theLastStr    
       

    def do(self,args,history):
        if not self.resubmit:
            return self.do_submit(args,history)
        else:
            return self.do_resubmit(args,history)


    def do_resubmit(self,args,history):

        # Start time
        start_time=time.time()

        # Checking argument number
        if len(args)!=0:
            logging.warning("Command 'resubmit' takes no argument. Any argument will be skipped.")

        # Checking presence of a valid job
        if self.main.lastjob_name is "":
            logging.error("an analysis must be defined and ran before using the resubmit command.") 
            return False

        self.main.lastjob_status = False

        if not self.submit(self.main.lastjob_name,history):
            return
        if not self.extract(self.main.lastjob_name):
            return

        # Status = GOOD
        self.main.lastjob_status = True

        # End time 
        end_time=time.time()
           
        logging.info("   Well done ! Elapsed time = " + CmdSubmit.chronometer_display(end_time-start_time) )
        
    def do_submit(self,args,history):

        # Start time
        start_time=time.time()

        # Checking argument number
        if len(args)!=1:
            logging.error("wrong number of arguments for the command 'submit'.")
            self.help()
            return

        # Checking if a dataset has been defined
        if len(self.main.datasets)==0:
            logging.error("no dataset found; please define a dataset (via the command import).")
            logging.error("job submission aborted.")
            return

        # Checking if a selection item has been defined
        if len(self.main.selection)==0:
            logging.error("no analysis found. Please define an analysis (via the command plot).")
            logging.error("job submission aborted.")
            return

        # Treat the filename
        filename = os.path.expanduser(args[0])
        if not filename.startswith('/'):
            filename = self.main.currentdir + "/" + filename
        filename = os.path.normpath(filename)

        # Checking folder
        if filename in self.forbiddenpaths:
            logging.error("the folder '"+filename+"' is a MadAnalysis folder. " + \
                         "You cannot overwrite it. Please choose another folder.")
            return

        # Saving job name as global variable
        self.main.lastjob_name = filename
        self.main.lastjob_status = False


        if not self.submit(filename,history):
            return
        if not self.extract(filename):
            return

        # Status = GOOD
        self.main.lastjob_status = True

        # End time 
        end_time=time.time()
           
        logging.info("   Well done ! Elapsed time = " + CmdSubmit.chronometer_display(end_time-start_time) )




    def submit(self,dirname,history):

        # Initializing the JobWriter
        jobber = JobWriter(self.main.ma5dir,\
                           dirname,self.resubmit, self.main.libZIP)
        
        # Writing process
        if not self.resubmit:
            logging.info("   Creating folder '"+dirname+"'...")
        else:
            logging.info("   Checking the structure of the folder '"+dirname+"'...")
        if not jobber.Open():
            logging.error("job submission aborted.")
            return False
        
        if not self.resubmit:
            logging.info("   Copying 'SampleAnalyzer' source files...")
            if not jobber.CopyLHEAnalysis():
                logging.error("   job submission aborted.")
                return False

        logging.info("   Inserting your selection into 'SampleAnalyzer'...")
        if not jobber.WriteSelectionHeader(self.main):
            logging.error("job submission aborted.")
            return False
        if not jobber.WriteSelectionSource(self.main):
            logging.error("job submission aborted.")
            return False

        logging.info("   Writing the list of datasets...")
        for item in self.main.datasets:
            jobber.WriteDatasetList(item)

        logging.info("   Writing the command line history...")
        jobber.WriteHistory(history,self.main.firstdir)

        if not self.resubmit:
            logging.info("   Creating a 'Makefile'...")
            if not jobber.WriteMakefile():
                logging.error("job submission aborted.")
                return False

        logging.info("   Compiling 'SampleAnalyzer'...")
        if not jobber.CompileJob():
            logging.error("job submission aborted.")
            return False

        logging.info("   Linking 'SampleAnalyzer'...")
        if not jobber.LinkJob():
            logging.error("job submission aborted.")
            return False

        for item in self.main.datasets:
            logging.info("   Running 'SampleAnalyzer' over dataset '"
                         +item.name+"'...")
            logging.info("    *******************************************************")
            if not jobber.RunJob(item):
                logging.error("run over '"+item.name+"' aborted.")
            logging.info("    *******************************************************")
    
        return True   


    def extract(self,dirname):
        logging.info("   Checking SampleAnalyzer output...")
        jobber = JobReader(dirname)
        if not jobber.Open():
            logging.error("errors have occured during the analysis.")
            return False
        
        for item in self.main.datasets:
            if not jobber.CheckRootFile(item):
                logging.error("errors have occured during the analysis.")
                return False

        logging.info("   Extracting data from the output files...")
        for item in self.main.datasets:
            jobber.ExtractDatasetInfo(item)

        return True    
           

    def help(self):
        if not self.resubmit:
            logging.info("   Syntax: submit <dirname>")
            logging.info("   Performs an analysis over a list of datasets. Output is stored into the directory <dirname>.")
        else:
            logging.info("   Syntax: resubmit")
            logging.info("   Allows to update an analysis which has been already performed at least one time at.")


    def complete(self,text,line,begidx,endidx):

        #Resubmission case 
        if self.resubmit:
            return 

        #Getting back arguments
        args = line.split()
        nargs = len(args)
        if not text:
            nargs += 1
        
        #Checking number of arguments
        if nargs==2:
            output=[]
            for file in glob.glob(text+"*"):
                if os.path.isdir(file):
                    output.append(file)
            return self.finalize_complete(text,output)
        else:
            return []
