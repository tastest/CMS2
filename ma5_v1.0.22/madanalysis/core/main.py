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


from madanalysis.multiparticle.multiparticle_collection import MultiParticleCollection
from madanalysis.dataset.dataset_collection             import DatasetCollection
from madanalysis.selection.selection                    import Selection
from madanalysis.enumeration.uncertainty_type           import UncertaintyType
from madanalysis.enumeration.normalize_type             import NormalizeType
from madanalysis.enumeration.sb_ratio_type              import SBratioType
from madanalysis.interpreter.cmd_base                   import CmdBase
from madanalysis.core.linux_architecture                import LinuxArchitecture
from madanalysis.core.config_checker                    import ConfigChecker
from madanalysis.core.library_builder                   import LibraryBuilder
from madanalysis.IOinterface.library_writer             import LibraryWriter
from madanalysis.enumeration.ma5_running_type           import MA5RunningType
from madanalysis.enumeration.stacking_method_type       import StackingMethodType
from madanalysis.observable.observable_manager          import ObservableManager
import logging
import os
import sys

class IsolationSet():


    userVariables = { "algo"   : ["DELTAR","SUMPT"], \
                      "deltaR" : ["0.5"], \
                      "ET_PT"  : ["1."], \
                      "sumPT"  : ["1."] }

    def __init__(self):
        self.algo   = 'DELTAR'
        self.deltaR = 0.5
        self.ET_PT  = 1.
        self.sumPT  = 1.

    def user_GetParameters(self):
        return IsolationSet.userVariables.keys()

    def user_GetValues(self,variable):
        try:
            return IsolationSet.userVariables[variable]
        except:
            return []

        
class Main():

    userVariables = { "currentdir"  : [], \
                      "normalize"   : ["none","lumi","lumi_weight"], \
                      "lumi"        : [], \
                      "SBratio"     : ['"S/B"','"B/S"','"S/(S+B)"','"B/(B+S)"','"S/sqrt(S+B)"','"B/sqrt(B+S)"'],
                      "SBerror"     : [], \
                      "stacking_method" : ["stack","superimpose","normalize2one"]  }

    SBformula = { 'S/B'         : '1./(B**2)*sqrt(B**2*ES**2+S**2*EB**2)', \
                  'S/(S+B)'     : '1./(S+B)**2*sqrt(B**2*ES**2+S**2*EB**2)', \
                  'S/sqrt(S+B)' : '1./pow(S+B,3./2.)*sqrt((S+2*B)**2*ES**2+S**2*EB**2)' }

    forced = False

    def __init__(self):
        self.currentdir     = os.getcwd()
        self.firstdir       = os.getcwd()
        self.version        = ""
        self.date           = ""
        self.configLinux    = LinuxArchitecture()
        self.mode           = MA5RunningType.PARTON
        self.forced         = False
        self.multiparticles = MultiParticleCollection()
        self.datasets       = DatasetCollection()
        self.selection      = Selection()
        self.script         = False
        self.mg5            = False
        self.libZIP         = False
        self.FAC            = False
        self.observables    = ObservableManager(self.mode)
        self.ResetParameters()


    def ResetParameters(self):
        self.SBratio        = 'S/B'
        self.SBerror        = Main.SBformula['S/B']
        self.lumi           = 10
        self.lastjob_name   = ''
        self.lastjob_status = False
        self.stack          = StackingMethodType.STACK
        self.isolation      = IsolationSet()
        if self.mode==MA5RunningType.RECO:
            self.normalize = NormalizeType.NONE
        else:
            self.normalize = NormalizeType.LUMI_WEIGHT


    def InitObservables(self,mode):
        self.observables    = ObservableManager(mode)
        
        
    def IsGoodFormat(self,file):
        for item in self.GetSampleFormat():
            if file.endswith(item):
                return True
        return False 
        

    def GetSampleFormat(self):
        samples = []
        # Adding format according to MA5 level mode
        if self.mode in [MA5RunningType.PARTON,MA5RunningType.HADRON]:
            samples.append('.lhe')
            samples.append('.hep')
            samples.append('.hepmc')
        else:
            samples.append('.lhco')
            if self.FAC:
                samples.append('.root')

        # Adding gzip file
        if self.libZIP:
            zipsamples=[]
            for item in samples:
                zipsamples.append(item+'.gz')
            samples.extend(zipsamples)

        return samples     


    def Display(self):
        logging.info(" *********************************" )
        logging.info("            main program          " )
        logging.info(" *********************************" )
        self.user_DisplayParameter("currentdir")
        self.user_DisplayParameter("normalize")
        self.user_DisplayParameter("lumi")
        self.user_DisplayParameter("SBratio")
        self.user_DisplayParameter("SBerror")
        self.user_DisplayParameter("algo")
        if self.isolation.algo=='DELTAR':
            self.user_DisplayParameter("deltaR")
        else:
            self.user_DisplayParameter("ET_PT")
            self.user_DisplayParameter("sumPT")
        logging.info(" *********************************" )


    def user_DisplayParameter(self,parameter):
        if  parameter=="currentdir":
            logging.info(" currentdir = "+self.get_currentdir())
        elif parameter=="stacking_method":
            sentence = " stacking methode for histograms = "
            if self.stack==StackingMethodType.STACK:
                sentence+="stack"
            elif self.stack==StackingMethodType.SUPERIMPOSE:
                sentence+="superimpose"
            else:
                sentence+="normalize2one"
            logging.info(sentence)
        elif parameter=="normalize":
            word=""
            if self.normalize==NormalizeType.NONE:
                word="none"
            elif self.normalize==NormalizeType.LUMI:
                word="lumi"
            elif self.normalize==NormalizeType.LUMI_WEIGHT:
                word="lumi_weight"
            logging.info(" histogram normalization mode = " + word)
        elif parameter=="lumi":
            logging.info(" integrated luminosity = "+str(self.lumi)+" fb^{-1}" )
        elif parameter=="SBratio":
            logging.info(' S/B ratio formula = "' + self.SBratio + '"')
        elif parameter=="SBerror":
            logging.info(' S/B error formula = "' + self.SBerror + '"')
        elif parameter=='algo':
            logging.info(' Muon isolation | algo = "' + self.isolation.algo + "'")
        elif parameter=='deltaR':
            logging.info(' Muon isolation | deltaR = ' + str(self.isolation.deltaR) )
        elif parameter=='sumPT':
            logging.info(' Muon isolation | sumPT = ' + str(self.isolation.sumPT) )
        elif parameter=='ET_PT':
            logging.info(' Muon isolation | ET_PT = ' + str(self.isolation.ET_PT) )
        else:
            logging.error("'main' has no parameter called '"+parameter+"'")


    def user_GetValues(self,variable):
        if variable=="currentdir":
            return CmdBase.directory_complete()
        else:
            try:
                return Main.userVariables[variable]
            except:
                return []

    def user_GetParameters(self):
        return Main.userVariables.keys()

    def user_SetParameter(self,parameter,value):
        # currentdir
        if parameter=="currentdir":
            self.set_currentdir(value)

        # stacked
        elif parameter=="stacking_method":
            if value == "stack":
                self.stack=StackingMethodType.STACK
            elif value == "superimpose":
                self.stack=StackingMethodType.SUPERIMPOSE
            elif value == "normalize2one":
                self.stack=StackingMethodType.NORMALIZE2ONE
            else:
                logging.error("'stack' possible values are : 'stack', 'superimpose', 'normalize2one'")
                return False

        # normalize
        elif parameter=="normalize":
            if value == "none":
                self.normalize = NormalizeType.NONE
            elif value == "lumi":
                self.normalize = NormalizeType.LUMI
            elif value == "lumi_weight":
                self.normalize = NormalizeType.LUMI_WEIGHT
            else:
                logging.error("'normalize' possible values are : 'none', 'lumi', 'lumi_weight'")
                return False

        # lumi    
        elif (parameter=="lumi"):
            try:
                tmp = float(value)
            except:
                logging.error("'lumi' is a positive float value")
                return
            if (tmp>0):
                self.lumi=tmp
            else:
                logging.error("'lumi' is a positive float value")
                return

        # sbratio    
        elif (parameter=="SBratio"):
            quoteTag=False
            if value.startswith("'") and value.endswith("'"):
                quoteTag=True
            if value.startswith('"') and value.endswith('"'):
                quoteTag=True
            if quoteTag:
                value=value[1:-1]
            if Main.checkSBratio(value):
                self.SBratio=value
                self.suggestSBerror()
            else:
                logging.error("Specified formula is not correct.")
                return False

        # sberror
        elif (parameter=="SBerror"):
            quoteTag=False
            if value.startswith("'") and value.endswith("'"):
                quoteTag=True
            if value.startswith('"') and value.endswith('"'):
                quoteTag=True
            if quoteTag:
                value=value[1:-1]
            if Main.checkSBratio(value):
                self.SBerror=value
            else:
                logging.error("Specified formula is not correct.")
                return False
            
        # algo
        elif parameter=="algo":
            if value in ['DELTAR','SUMPT']:
                self.isolation.algo=value
            else:
                logging.error("Specified algo for muon isolation is not correct. "+\
                              "Only possible ones are : 'DELTAR' and 'SUMPT'.")
                return False
            
        # deltaR
        elif parameter=="deltaR":
            try:
                tmp=float(value)
                self.isolation.deltaR=tmp
            except:
                logging.error("'"+value+"' is not a float value.")
                return False

        # ET_PT
        elif parameter=="ET_PT":
            try:
                tmp=float(value)
                self.isolation.ET_PT=tmp
            except:
                logging.error("'"+value+"' is not a float value.")
                return False

        # sumPT
        elif parameter=="sumPT":
            try:
                tmp=float(value)
                self.isolation.sumPT=tmp
            except:
                logging.error("'"+value+"' is not a float value.")
                return False

        # other    
        else:
            logging.error("'main' has no parameter called '"+parameter+"'")

    @staticmethod
    def checkSBratio(text):
        logging.info("Checking the formula ...")
        text = text.replace("ES","z")
        text = text.replace("EB","t")
        text = text.replace("S","x")
        text = text.replace("B","y")
        from ROOT import TFormula
        formula = TFormula()
        test = formula.Compile(text)
        return (test==0)

    def suggestSBerror(self):
        # create a TFormula with the SBratio formula
        text = self.SBratio.replace("S","x")
        text = text.replace("B","y")
        from ROOT import TFormula
        ref = TFormula('SBratio',text)
        ref.Optimize()

        # Loop over SBerror formula and comparing
        for k, v in Main.SBformula.iteritems():
            text = k.replace("S","x")
            text = text.replace("B","y")
            error = TFormula('SBerror',text)
            error.Optimize()
            if ref.GetExpFormula()==error.GetExpFormula():
                logging.info("Formula corresponding to the uncertainty calculation has been found and set to the variable main.SBerror :")
                logging.info('  '+v)
                self.SBerror=v
                return True

        # Loop over SBerror formula and comparing
        # reverse S and B
        for k, v in Main.SBformula.iteritems():
            text = k.replace("S","y")
            text = text.replace("B","x")
            error = TFormula('SBerror',text)
            error.Optimize()
            if ref.GetExpFormula()==error.GetExpFormula():
                logging.info("Formula corresponding to the uncertainty calculation has been found and set to the variable main.SBerror :")
                v=v.replace('ES','ZZ')
                v=v.replace('EB','TT')
                v=v.replace('S','SS')
                v=v.replace('B','BB')
                v=v.replace('SS','B')
                v=v.replace('BB','S')
                v=v.replace('ZZ','EB')
                v=v.replace('TT','ES')
                logging.info('  '+v)
                self.SBerror=v
                return True


    def get_currentdir(self):
        return os.getcwd()

    def set_currentdir(self,dir):
        theDir = os.path.expanduser(dir)
        try:
            os.chdir(theDir)
        except:
            logging.error("Impossible to access the directory : "+theDir) 
        self.user_DisplayParameter("currentdir")
        
    currentdir = property(get_currentdir, set_currentdir)

    def CheckLinuxConfig(self):
        self.configLinux.ma5_version    = self.version
        self.configLinux.ma5_date       = self.date
        self.configLinux.python_version = sys.version.replace('\n','')
        self.configLinux.platform       = sys.platform

        checker = ConfigChecker(self.configLinux,self.ma5dir,self.script)
        
        if not checker.checkROOT():
            return False
        if not checker.checkGPP():
            return False
        self.libZIP = checker.checkZLIB()

        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + \
                                        ":" + self.ma5dir+'/lib/'
        os.environ['LIBRARY_PATH'] = os.environ['LIBRARY_PATH'] + \
                                        ":" + self.ma5dir+'/lib/'
        return True 
   
    def BuildLibrary(self):
        builder = LibraryBuilder(self.configLinux,self.ma5dir)
        UpdateNeed=False
        
        FirstUse = builder.checkMA5()
        if not FirstUse:
            UpdateNeed = not builder.compare()

        if not FirstUse and not UpdateNeed:
            return True

        # Compile library
        if FirstUse:
            logging.info("First use of MadAnalysis detected (or the library is missing)")
        elif UpdateNeed:
            logging.info("System configuration has changed since the last use. Need to rebuild the library")

        # Initializing the JobWriter
        compiler = LibraryWriter(self.ma5dir,'lib',self.libZIP,self.FAC)
    
        # Writing process
        logging.info("   Creating (or overwriting) folder 'lib'...")
        if not compiler.Open():
            logging.error("library building aborted.")
            sys.exit()
        if not os.path.isdir(self.ma5dir+'/lib'):
            logging.error('impossible to create the directory : ' + self.ma5dir+'/lib')
            logging.error("library building aborted.")
            sys.exit()

        # Dumping architecture
        if not self.configLinux.Export(self.ma5dir+'/lib/architecture.ma5'):
            sys.exit()

        logging.info("   Copying 'SampleAnalyzer' source files...")
        if not compiler.CopySampleAnalyzer():
            logging.error("   library building aborted.")
            sys.exit()

        logging.info("   Creating a 'Makefile'...")
        if not compiler.WriteMakefile('shared'):
            logging.error("library building aborted.")
            sys.exit()

        logging.info("   Compiling the MadAnalysis library...")
        if not compiler.Compile():
            logging.error("library building aborted.")
            sys.exit()

        logging.info("   Linking the MadAnalysis library...")
        if not compiler.Link():
            logging.error("library building aborted.")
            sys.exit()

        logging.info("   Checking the MadAnalysis library presence...")

        if not os.path.isfile(self.ma5dir+'/lib/SampleAnalyzer/libSampleAnalyzer.a'):
            logging.error("the library 'libSampleAnalyzer.a' is not produced.")
            sys.exit()

        # Moving the library into a different folder
        try:
            import shutil
            shutil.move(self.ma5dir+"/lib/SampleAnalyzer/libSampleAnalyzer.a",
            self.ma5dir+"/lib/libSampleAnalyzer.a")
        except:
            logging.error("impossible to move the SampleAnalyzer library into "\
	                 + " the folder '" + self.ma5dir+"/lib'")
            sys.exit()

        # Removing the SampleAnalyzer directory
        try: 
            shutil.rmtree(self.ma5dir+'/lib/SampleAnalyzer')
        except:
            logging.warning("impossible to remove the temporary folder "\
		    	    + "'"+ self.ma5dir + "/lib/SampleAnalyzer'")
    
        return True    
        
