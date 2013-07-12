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


from madanalysis.selection.instance_name import InstanceName
from madanalysis.IOinterface.folder_writer import FolderWriter
import logging
import shutil
import os
import commands

class LibraryWriter():

    def __init__(self,ma5dir,jobdir,libZIP,FAC):
        self.ma5dir     = ma5dir
        self.jobdir     = jobdir
        self.path       = os.path.normpath(ma5dir+"/"+jobdir)
        self.libZIP     = libZIP
        self.FAC        = FAC

    def Open(self):
        return FolderWriter.CreateDirectory(self.path,overwrite=True)

    def CopySampleAnalyzer(self):
        try:
            shutil.copytree\
                      (\
                      self.ma5dir+"/madanalysis/SampleAnalyzer",\
                      self.path+"/SampleAnalyzer/",\
                      ignore=shutil.ignore_patterns('*.pyc','*.o','*~','.*')\
                      )
        except:
            logging.error("An error occured during copying 'SampleAnalyzer'\
            source files.")
            return False

        return True

    def WriteMakefile(self,option=""):
        file = open(self.path + "/SampleAnalyzer/Makefile","w")
        file.write('GCC = g++\n')
        file.write('CXXFLAGS = `root-config --cflags` -I./')
        if self.libZIP:
            file.write(' -DZIP_USE')
        if self.FAC:
            file.write(' -DFAC_USE')
        file.write('\n')   
        #file.write('LIBFLAGS = `root-config --libs` -lboost_iostreams -lz\n\n')
        file.write('SRCS = $(wildcard */*.cpp)\n')
        file.write('OBJS = $(SRCS:.cpp=.o)\n')
        file.write('PROGRAM = SampleAnalyzer\n\n')
        file.write('all:\t precompile compile link\n\n')
        file.write('precompile:')
        if self.FAC:
            file.write('\n\t\trootcint -f "Reader/FACdict.cpp" -c -p Reader/FACdataformat.h Reader/FACLinkDef.h\n\n')
        else:
            file.write('\n\n')
        file.write('compile:\tprecompile $(OBJS)\n\n')
        file.write('link:\t$(OBJS)\n')
        file.write('\t\tar -ruc lib$(PROGRAM).a $(OBJS)\n')
        file.write('\t\tranlib lib$(PROGRAM).a\n\n')
        file.write('clean:;\t@rm -f $(OBJS) lib$(PROGRAM).a compilation.log linking.log *~ */*~ \n')
        file.close()
        return True


    def Compile(self):
        res=commands.getstatusoutput("cd "\
                                     +self.path+"/SampleAnalyzer/;"\
                                     +" make compile > compilation.log 2>&1")
        if res[0]==0:
            return True
        else:
            logging.error("errors occured during compilation. For more details, see the file :")
            logging.error(" "+self.path+"/SampleAnalyzer/compilation.log")
            return False

    def Link(self):
        res=commands.getstatusoutput("cd "\
                                     +self.path+"/SampleAnalyzer/;"\
                                     +" make link > linking.log 2>&1")
        if res[0]==0:
            return True
        else:
            logging.error("errors occured during compilation. For more details, see the file :")
            logging.error(" "+self.path+"/SampleAnalyzer/linking.log")
            return False


        
