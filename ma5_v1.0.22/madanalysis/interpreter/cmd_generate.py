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
from madanalysis.layout.layout import Layout
from madanalysis.enumeration.report_format_type import ReportFormatType

import logging
import os
import glob

class CmdGenerate(CmdBase.CmdBase):
    """Command GENERATE"""

    def __init__(self,main,format):
        self.format  = format
        self.cmdname = ReportFormatType.convert2cmd(self.format)
        CmdBase.CmdBase.__init__(self,main,self.cmdname)
        self.forbiddenpaths=[]
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/lib'))
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/bin'))
        self.forbiddenpaths.append(os.path.normpath(self.main.ma5dir+'/madanalysis'))

    def do_possible(self):

        if self.main.lastjob_name is "":
            logging.error("an analysis must be defined and ran before using the '"+self.cmdname+"' command.") 
            return False
        if not self.main.lastjob_status:
            logging.error("errors occured during the analysis '"+self.main.lastjob_name+"'.")
            logging.error("Please resolve the problem before calling the '"+self.cmdname+"' command.")
            return False
        return True
        

    def do(self,args,history):

        # Checking argument number
        if len(args)!=1:
            logging.error("wrong number of arguments for the command '"+self.cmdname+"'.")
            self.help()
            return

        # Checking presence of a valid job
        if not self.do_possible():
            return

        # Getting filename
        filename = os.path.expanduser(args[0])
        if not filename.startswith('/'):
            filename = self.main.currentdir + "/" + filename
        filename = os.path.normpath(filename)

        # Checking folder
        if filename in self.forbiddenpaths:
            logging.error("the folder '"+filename+"' is MadAnalysis folder. " + \
                         "You cannot overwrite it. Please choose another folder.")
            return

        layout = Layout(self.main,filename,self.format)

        logging.info("   Preparing data for the report...")
        if not layout.Open():
            return

        if not layout.CopyLogo():
            return

        logging.info("   Saving all plots as image files...")
        if not layout.DoPlots():
            return
        
        logging.info("   Computing cut efficiencies...")
        if not layout.DoEfficiencies():
            return

        logging.info("   Generating the report...")
        if not layout.GenerateReport(history):
            return

        if self.format!=ReportFormatType.HTML:
            if not self.main.forced:
                logging.warning("Would you like to launch latex for creating a pdf file ? (Y/N)")
                
                allowed_answers=['n','no','y','yes']
                answer=""
                while answer not in  allowed_answers:
                   answer=raw_input("Answer : ")
                   answer=answer.lower()
                if answer=="no" or answer=="n":
                    return
                
            layout.CompileReport()
            

    def help(self):
        logging.info("   Syntax: "+self.cmdname+" <output dir>")
        logging.info("   Creates a "+ReportFormatType.convert2string(self.format)+ \
                     " report containing the results of the performed analysis.")

    def complete(self,text,line,begidx,endidx):

        #Getting back arguments
        args = line.split()
        nargs = len(args)
        if not text:
            nargs += 1
        
        #Checking number of arguments
        if nargs==2:
            output=[]
            for file in glob.glob(text+"*"):

                # directory presence 
                if not os.path.isdir(file):
                    continue
                
                output.append(file)

            return self.finalize_complete(text,output)
        else:
            return []

