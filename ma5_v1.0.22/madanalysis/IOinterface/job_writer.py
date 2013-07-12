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


from madanalysis.selection.instance_name      import InstanceName
from madanalysis.IOinterface.folder_writer    import FolderWriter
from madanalysis.enumeration.ma5_running_type import MA5RunningType
import logging
import shutil
import os
import commands

class JobWriter():

    def __init__(self,ma5dir,jobdir,resubmit=False,libZIP=True):
        self.ma5dir     = ma5dir
        self.path       = jobdir
        self.resubmit   = resubmit
        self.libZIP     = libZIP 

    @staticmethod     
    def CheckJobStructureMute(path):
        if not os.path.isdir(path):
            return False
        elif not os.path.isdir(path+"/SampleAnalyzer"):
            return False
        elif not os.path.isdir(path+"/SampleAnalyzer/Analysis"):
            return False
        elif not os.path.isdir(path+"/root"):
            return False
        elif not os.path.isdir(path+"/lists"):
            return False
        elif not os.path.isfile(path+"/history.ma5"):
            return False
        else:
            return True


    def CheckJobStructure(self):
        if not os.path.isdir(self.path):
            logging.error("folder '"+self.path+"' is not found")
            return False
        elif not os.path.isdir(self.path+"/SampleAnalyzer"):
            logging.error("folder '"+self.path+"/SampleAnalyzer' is not found")
            return False
        elif not os.path.isdir(self.path+"/SampleAnalyzer/Analysis"):
            logging.error("folder '"+self.path+"/SampleAnalyzer/Analysis' is not found")
            return False
        elif not os.path.isdir(self.path+"/root"):
            logging.error("folder '"+self.path+"/root' is not found")
            return False
        elif not os.path.isdir(self.path+"/lists"):
            logging.error("folder '"+self.path+"/lists' is not found")
            return False
        elif not os.path.isfile(self.path+"/history.ma5"):
            logging.error("file '"+self.path+"/history.ma5' is not found")
            return False
        else:
            return True

    def Open(self):
        if not self.resubmit:
            InstanceName.Clear()
            return FolderWriter.CreateDirectory(self.path,question=True)
        else:
            return self.CheckJobStructure()

    def CopyLHEAnalysis(self):
        try:
            os.mkdir(self.path+"/SampleAnalyzer")
            os.mkdir(self.path+"/SampleAnalyzer/Analysis")
            shutil.copyfile\
                      (\
                      self.ma5dir+"/madanalysis/SampleAnalyzer/newAnalysis.py",\
                      self.path+"/SampleAnalyzer/newAnalysis.py"\
                      )
        except:
            logging.error("An error occured during copying 'SampleAnalyzer'" +\
            "source files.")
            return False

        return True

    def WriteSelectionHeader(self,main):
        main.selection.RefreshStat();
        file = open(self.path+"/SampleAnalyzer/Analysis/user.h","w")
        import madanalysis.job.job_main as JobMain
        job = JobMain.JobMain(file,main)
        job.WriteHeader()
        file.close()
        return True

    def WriteSelectionSource(self,main):
        main.selection.RefreshStat();
        file = open(self.path+"/SampleAnalyzer/Analysis/user.cpp","w")
        import madanalysis.job.job_main as JobMain
        job = JobMain.JobMain(file,main)
        job.WriteSource()
        file.close()

        file = open(self.path+"/SampleAnalyzer/Analysis/analysisList.cpp","w")
        file.write('#include "Core/AnalysisManager.h"\n')
        file.write('#include "Analysis/user.h"\n')
        file.write('#include "Services/logger.h"\n')
        file.write('#include <stdlib.h>\n\n')
        file.write('// ------------------------------------------' +\
                   '-----------------------------------\n')
        file.write('// BuildTable\n')
        file.write('// ------------------------------------------' +\
                   '-----------------------------------\n')
        file.write('void AnalysisManager::BuildTable()\n')
        file.write('{\n')
        file.write('  Add(new user);\n')
        file.write('}\n')
        file.close()

        return True

    def WriteMakefile(self,option=""):
        file = open(self.path+"/SampleAnalyzer/Makefile","w")
        file.write('GCC = g++\n')
#        file.write('CXXFLAGS = -fPIC `root-config --cflags` -I./ -I' +\
#                   self.ma5dir + '/madanalysis/SampleAnalyzer/\n')
#        file.write('LIBFLAGS = `root-config --libs` -lboost_iostreams -lz' +\
#                   ' -L' + self.ma5dir + '/lib/ -lSampleAnalyzer\n\n')
        file.write('CXXFLAGS = -I./ -I' +\
                   self.ma5dir + '/madanalysis/SampleAnalyzer/ -lSampleAnalyzer\n')
#        file.write('LIBFLAGS = -lCore -lHist -lPhysics -lMathCore -lMatrix -lRIO -lboost_iostreams -lz' +\
#                   ' -lSampleAnalyzer\n\n')
        file.write('LIBFLAGS = -lGpad -lHist -lGraf -lGraf3d ' +\
                   '-lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore ' +\
                   '-lRIO -lNet -lThread -lCore -lCint -pthread -lm -ldl -rdynamic ' +\
                   '-lSampleAnalyzer')
        if self.libZIP:
            file.write(' -lz')
        file.write('\n\n')
        file.write('SRCS = $(wildcard */*.cpp)\n')
        file.write('OBJS = $(SRCS:.cpp=.o)\n')
        file.write('PROGRAM = SampleAnalyzer\n\n')
        file.write('all:\t compile link\n\n')
        file.write('compile:\t$(OBJS)\n\n')
        file.write('link:\t$(OBJS)\n')
        file.write('\t\t$(GCC) $(CXXFLAGS) $(OBJS) ')
        file.write('$(LIBFLAGS) -o $(PROGRAM)\n \n')
        file.write('clean:;\t@rm -f $(OBJS) $(PROGRAM) compilation.log linking.log *~ */*~ \n')
        file.close()

        file = open(self.path+"/SampleAnalyzer/setup.csh","w")
        file.write('#!/bin/csh -f\n')
        file.write('setenv LD_LIBRARY_PATH '    + os.environ['LD_LIBRARY_PATH']+'\n')
        file.write('setenv LIBRARY_PATH '       + os.environ['LD_LIBRARY_PATH']+'\n')
        file.write('setenv DYLD_LIBRARY_PATH '  + os.environ['DYLD_LIBRARY_PATH']+'\n')
        file.write('setenv CPLUS_INCLUDE_PATH ' + os.environ['CPLUS_INCLUDE_PATH']+'\n')
        file.close()
        
        file = open(self.path+"/SampleAnalyzer/setup.sh","w")
        file.write('#!/bin/sh\n')
        file.write('export LD_LIBRARY_PATH='    + os.environ['LD_LIBRARY_PATH']+'\n')
        file.write('export LIBRARY_PATH='       + os.environ['LD_LIBRARY_PATH']+'\n')
        file.write('export DYLD_LIBRARY_PATH='  + os.environ['DYLD_LIBRARY_PATH']+'\n')
        file.write('export CPLUS_INCLUDE_PATH=' + os.environ['CPLUS_INCLUDE_PATH']+'\n')
        file.close()

        return True

    def CompileJob(self):
        res=commands.getstatusoutput("cd "\
                                     +self.path+"/SampleAnalyzer/;"\
                                     +" make compile > compilation.log 2>&1")
        if res[0]==0:
            return True
        else:
            logging.error("errors occured during compilation. " +\
                          "For more details, see the file :")
            logging.error(" "+self.path+"/SampleAnalyzer/compilation.log")
            return False

    def LinkJob(self):
        res=commands.getstatusoutput("cd "\
                                     +self.path+"/SampleAnalyzer/;"\
                                     +" make link > linking.log 2>&1")
        if res[0]==0:
            return True
        else:
            logging.error("errors occured during compilation. For more details, see the file :")
            logging.error(" "+self.path+"/SampleAnalyzer/linking.log")
            return False

    def WriteHistory(self,history,firstdir):
        file = open(self.path+"/history.ma5","w")
        file.write('set main.currentdir = '+firstdir+'\n') 
        for line in history:
            items = line.split(';')
            for item in items :
                if item.startswith('help') or \
                   item.startswith('display') or \
                   item.startswith('generate') or \
                   item.startswith('history') or \
                   item.startswith('open') or \
                   item.startswith('preview') or \
                   item.startswith('resubmit') or \
                   item.startswith('shell') or \
                   item.startswith('!') or \
                   item.startswith('submit'):
                    pass
                else:
                    file.write(item)
                    file.write("\n")
        file.close()    

    def WriteDatasetList(self,dataset):
        if not os.path.isdir(self.path+"/lists"):
            os.mkdir(self.path+"/lists")
        name=InstanceName.Get(dataset.name)
        file = open(self.path+"/lists/"+name+".list","w")
        for item in dataset:
            file.write(item)
            file.write("\n")
        file.close()    

    def RunJob(self,dataset):
        if not os.path.isdir(self.path+"/root"):
            os.mkdir(self.path+"/root")
        name=InstanceName.Get(dataset.name)
        res=os.system('cd '\
                      +self.path+'/root;'\
                      +' ../SampleAnalyzer/'\
                      +'SampleAnalyzer --analysis=MadAnalysis5job ../lists/'+
                      name+'.list')
        return True

        
        
        
        
        
