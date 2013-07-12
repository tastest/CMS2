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


from madanalysis.interpreter.cmd_define           import CmdDefine
from madanalysis.enumeration.ma5_running_type     import MA5RunningType
from madanalysis.IOinterface.ufo_reader           import UFOReader
from madanalysis.IOinterface.job_writer           import JobWriter
from madanalysis.IOinterface.particle_reader      import ParticleReader
from madanalysis.IOinterface.multiparticle_reader import MultiparticleReader
from madanalysis.interpreter.cmd_define           import CmdDefine
from madanalysis.IOinterface.job_writer           import JobWriter
from madanalysis.IOinterface.job_reader           import JobReader
import madanalysis.interpreter.cmd_base as CmdBase
import logging
import glob
import os

class CmdImport(CmdBase.CmdBase):
    """Command IMPORT"""

    def __init__(self,main):
        CmdBase.CmdBase.__init__(self,main,"import")


    def do(self,args,myinterpreter):

        # Checking argument number
        if len(args)!=3 and len(args)!=1 :
            logging.error("wrong number of arguments for the command 'import'.")
            self.help()
            return

        # Getting dataset name
        if len(args)==3:
            if not args[1] == 'as':
                logging.error("syntax error with the command 'import'.")
                self.help()
                return

        # Getting filename
        filename = os.path.expanduser(args[0])

        # Normalize path
        filename = os.path.normpath(filename)

        # Checking if a directory
        if len(args)==3:
            if os.path.isdir(filename):
                filename += '/*'
                filename = os.path.normpath(filename)
            self.ImportDataset(filename,args[2])
              
            
        elif len(args)==1 and os.path.isdir(filename):
            if JobWriter.CheckJobStructureMute(filename):
                self.ImportJob(filename,myinterpreter)
                return
            elif UFOReader.CheckStructure(filename):
                self.ImportUFO(filename)
                return
            else:
                filename += '/*'
                filename = os.path.normpath(filename)
                self.ImportDataset(filename,"defaultset")
                return
        else:
            self.ImportDataset(filename,"defaultset")
            return
            

    def ImportUFO(self,filename):
        logging.info("UFO model folder is detected")

        # UFO mode is forbidden in RECO level
        if self.main.mode==MA5RunningType.RECO:
            logging.error("cannot import particles in MA5 reconstructed mode")
            return False
        
        # Other modes : parton and hadron
        logging.info("Import all particles defined in the model ...")
        cmd_define = CmdDefine(self.main)
        
        ufo = UFOReader(filename,cmd_define)
        if not ufo.OpenParticle():
            return False
        if not ufo.ReadParticle():
            return False
        if not ufo.CloseParticle():
            return False
        if not ufo.OpenParameter():
            return False
        if not ufo.ReadParameter():
            return False
        if not ufo.CloseParameter():
            return False

        if len(ufo.parts.parts)==0:
            logging.warning("UFO model contained no particles")
            return

        # Reseting only particles, kept multiparticles
        self.main.multiparticles.ResetParticles()

        # Use new particles
        ufo.CreateParticle()
        
        


    def ImportJob(self,filename,myinterpreter):
        logging.info("SampleAnalyzer job folder is detected")
        logging.info("Restore MadAnalysis configuration used for this job ...")

        # Ask question
        if not self.main.forced:
            logging.warning("You are going to reinitialize MadAnalysis 5. The current configuration will be lost.")
            logging.warning("Are you sure to do that ? (Y/N)")
            allowed_answers=['n','no','y','yes']
            answer=""
            while answer not in  allowed_answers:
               answer=raw_input("Answer : ")
               answer=answer.lower()
            if answer=="no" or answer=="n":
                return False

        self.main.datasets.Reset()

        # Reset selection
        self.main.selection.Reset()

        # Reset main
        self.main.ResetParameters()

        # Reset multiparticles
        self.main.multiparticles.Reset()

        # Opening a CmdDefine
        cmd_define = CmdDefine(self.main)

        # Loading particles
        input = ParticleReader(self.main.ma5dir,cmd_define,self.main.mode)
        input.Load()
        input = MultiparticleReader(self.main.ma5dir,cmd_define,self.main.mode)
        input.Load()

        # Reset history
        myinterpreter.history=[] 

        # Load script
        myinterpreter.load(filename+'/history.ma5')

        # Saving job name as global variable
        self.main.lastjob_name = filename
        self.main.lastjob_status = False

        # Extract info from ROOT file
        if not self.extract(filename):
            return

        # Initialize selection for generating report
        self.main.selection.RefreshStat() 
 
        # Status = GOOD
        self.main.lastjob_status = True


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
       
    def ImportDataset(self,filename,name):

        # Creating dataset if not exist
        newdataset=False
        if not self.main.datasets.Find(name):
            if not self.create(name):
                return
            else:
                newdataset=True

        # Calling fill
        if not self.fill(name,filename) and newdataset:
            
            # if the dataset has just been created but empty : remove it 
            self.main.datasets.Remove(name)


    def create(self,name):
        
        # Checking if the name is authorized
        if name in self.reserved_words:
            logging.error("name '" +name+ "' is a reserved keyword. Please choose another name.")
            return False

        # Checking if the name is authorized
        if not self.IsAuthorizedLabel(name):
            logging.error("syntax error with the name '" + name + "'.")
            logging.error("A correct name contains only characters being letters, digits or the '+', '-', '~' and '_' symbols.")
            logging.error("Moreover, a correct name  starts with a letter or the '_' symbol.")
            return False

        # Checking if no multiparticle with the same name has been defined
        if self.main.multiparticles.Find(name):
            logging.error("A (multi)particle '"+name+"' already exists. Please choose a different name.")
            return False

        # Creating dataset
        self.main.datasets.Add(name)
        return True

        
    def fill(self,name,filename):

        # Getting the dataset
        set = self.main.datasets.Get(name)
        
        # Getting all files corresponding to filename
        files=[]
        for file in glob.glob(filename):
            if os.path.isfile(file) and self.main.IsGoodFormat(file):
                   #( file.endswith(".lhe") or file.endswith(".lhe.gz")):
                files.append(file)
        

        # If no file
        if len(files)==0:
            logging.error("The dataset '"+filename+"' has not been found or has a unsupported format.")
            return False
            
        # Getting current dir
        theDir = os.getcwd()

        # Adding file
        for item in files:
            if item.startswith('/'):
                theFile = item
            else:    
                theFile = os.path.normpath(theDir+"/"+item)
            logging.info("   -> Storing the file '"+theFile+"' in the dataset '"+name+"'.")
            set.Add(theFile)

        return True    



    def help(self):
        logging.info("   Syntax: import <Sample file> as <dataset name>")
        logging.info("   Stores one or several data file(s) in a given dataset.")
        logging.info("   The supported event file formats are: ")
        logging.info("     - LHE (lhe or lhe.gz),")
        logging.info("     - HEPMC (hepmc or hepmc.gz),")
        logging.info("     - STDHEP (hep or hep.gz),") 
        logging.info("     - LHCO (lhco or lhco.gz).") 
        logging.info("   If the dataset does not exist, it is created.")



    def complete(self,text,line,begidx,endidx):

        #Getting back arguments
        if len(line)==0:
            return []
        
        args = line.split()
        nargs = len(args)
        if line[-1]==' ' or line[-1]=='\t':
            nargs += 1
        #if not text:
        #    nargs += 1
        #print "\n"
        #print args
        #print nargs
        #print text
        
        #Checking number of arguments
        if nargs==2:
            output=[]
            for file in glob.glob(text+"*"):
                if os.path.isfile(file):
                    if self.main.IsGoodFormat(file):
                        output.append(file)
                else:
                    output.append(file) # directory
            return self.finalize_complete(text,output)
        elif nargs==3:
            output=["as"]
            return self.finalize_complete(text,output)
        elif nargs==4:
            output=self.main.datasets.GetNames()
            return self.finalize_complete(text,output)
        else:
            return

