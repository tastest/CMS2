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


import logging
import glob
import os
import commands
import sys

class ConfigChecker:

    @staticmethod
    def AddIfValid(path,container):
        dirs=glob.glob(path)
        for item in dirs:
            if not (item in container):
                container.append(item)
        

    def __init__(self,configLinux,ma5dir,script=False):

        # Getting paarmeter from the main program
        self.configLinux=configLinux
        self.ma5dir=ma5dir
        self.script=script

        self.libs = []
        # Filling container with paths included in LD_LIBRARY_PATH
        try:
            ld_library_path = os.environ['LD_LIBRARY_PATH'].split(':')
            for item in ld_library_path:
                ConfigChecker.AddIfValid(item,self.libs)
        except:
            os.environ['LD_LIBRARY_PATH']=''

        # Filling container with paths included in DYLD_LIBRARY_PATH
        try:
            ld_library_path = os.environ['DYLD_LIBRARY_PATH'].split(':')
            for item in ld_library_path:
                ConfigChecker.AddIfValid(item,self.libs)
        except:
            os.environ['DYLD_LIBRARY_PATH']=''
                
        # Filling container with paths included in LIBRARY_PATH
        try:
            library_path = os.environ['LIBRARY_PATH'].split(':')
            for item in library_path:
                ConfigChecker.AddIfValid(item,self.libs)
        except:
            os.environ['LIBRARY_PATH']=''

        # Filling container with standard library paths
        ConfigChecker.AddIfValid('/usr/lib*',self.libs)
        ConfigChecker.AddIfValid('/usr/local/lib*',self.libs)
        ConfigChecker.AddIfValid('/local/lib*',self.libs)
        ConfigChecker.AddIfValid('/opt/local/lib*',self.libs)

        self.includes = []

        # Filling container with paths included in CPLUS_INCLUDE_PATH
        try:
            cplus_include_path = os.environ['CPLUS_INCLUDE_PATH'].split(':')
            for item in cplus_include_path:
                ConfigChecker.AddIfValid(item,self.includes)
        except:
            os.environ['CPLUS_INCLUDE_PATH']=''

        # Filling container with standard include paths
        ConfigChecker.AddIfValid('/usr/include',self.includes)
        ConfigChecker.AddIfValid('/usr/local/include',self.includes)
        ConfigChecker.AddIfValid('/local/include',self.includes)
        ConfigChecker.AddIfValid('/opt/local/include',self.includes)


    def checkROOT(self):
        # Checking if ROOT is present
        logging.info("Checking ROOT libraries ...")

        # Trying to call root-config
        rootdirs = commands.getstatusoutput('root-config --libdir --incdir')
        if rootdirs[0]>0:
            logging.error('ROOT module called "root-config" is not detected.\n'\
		          +'Two explanations :n'\
		          +' - ROOT is not installed. You can download it '\
		          +'from http://root.cern.ch\n'\
		          +' - ROOT binary folder must be placed in the '\
                          +'global environment variable $PATH')
            return False

        # Extracting ROOT library and header path
        root_tmp = rootdirs[1].split() 
        self.includes.append(root_tmp[1])
        self.libs.append(root_tmp[0])
        os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + \
                    		        ":" + root_tmp[0]
        os.environ['DYLD_LIBRARY_PATH'] = os.environ['DYLD_LIBRARY_PATH'] + \
                    		        ":" + root_tmp[0]
        os.environ['LIBRARY_PATH'] = os.environ['LIBRARY_PATH'] + \
                    		        ":" + root_tmp[0]
        os.environ['CPLUS_INCLUDE_PATH'] = os.environ['CPLUS_INCLUDE_PATH'] + \
                    		        ":" + root_tmp[1]

        # Adding ROOT library path to Python path
        sys.path.append(root_tmp[0])

        # Looking for libPyROOT.so
        find=False
        for item in self.libs:
            files=glob.glob(item+"/libPyROOT.so")
            if len(files)!=0:
	        self.configLinux.libraries['PyROOT']=files[0]+":"+str(os.stat(files[0]).st_mtime)
    	        find=True
   	        break
        if not find:
	    logging.error("ROOT library called 'libPyROOT.so' is not found. Please check that ROOT is properly installed.")
            return False

        # Looking for ROOT.py
        find=False
        for item in self.libs:
            files=glob.glob(item+"/ROOT.py")
            if len(files)!=0:
	        self.configLinux.libraries['ROOT']=files[0]+":"+str(os.stat(files[0]).st_mtime)
	        find=True
	        break
        if not find:
	    logging.error("ROOT file called 'ROOT.py' is not found. Please check that ROOT is properly installed.")
            return False

        # Looking for TH1F.h
        find=False
        for item in self.includes:
            files=glob.glob(item+"/TH1F.h")
            if len(files)!=0:
	        self.configLinux.headers['ROOT']=files[0]+":"+str(os.stat(files[0]).st_mtime)
	        find=True
	        break
        if not find:
	    logging.error("ROOT headers are not found. " +\
		 "Please check that ROOT is properly installed.")
            return False

        # Loading ROOT library
        logging.info("Loading ROOT libraries ...")
        try :
	    from ROOT import gROOT
        except:
            logging.error("'root-config --libdir' indicates a wrong path for ROOT"\
	                  +" libraries. Please specify the ROOT library path"\
		          +" into the environnement variable $PYTHONPATH")
            return False

        # Setting ROOT batch mode
        if not self.script:
            from ROOT import TApplication
            from ROOT import gApplication
            TApplication.NeedGraphicsLibs()
            gApplication.InitializeGraphics()
        gROOT.SetBatch(True)
        

        # Checking ROOT release
        RootVersion = str(gROOT.GetVersionInt())
        if len(RootVersion)<3:
	    logging.error('Bad release of ROOT : '+gROOT.GetVersion()+\
                          '. MadAnalysis5 needs ROOT 5.27 or higher.\n Please upgrade your version of ROOT.')
            return False

        RootVersionA = int(RootVersion[0])
        RootVersionB = int(RootVersion[1]+RootVersion[2])
        if RootVersionA!=5 or RootVersionB<27:
	    logging.error('Bad release of ROOT : '+gROOT.GetVersion()+\
                          '. MadAnalysis5 needs ROOT 5.27 or higher.\n Please upgrade your version of ROOT.')
            return False

        self.configLinux.root_version   = RootVersion
        return True


    def checkGPP(self):
        # Checking g++ release
        logging.info("Checking g++ libraries ...")
        gcc_version = commands.getstatusoutput('g++ -dumpversion')
        if gcc_version[0]>0:
            logging.error('g++ compiler is not found. Please install it before ' + \
	             'using MadAnalysis 5')
            return False
        else:
            self.configLinux.gcc_version = gcc_version[1]
            return True


    def checkZLIB(self):

        logging.info("Checking zlib libraries ...")
        
        # Checking library libz.so
        find=False
        for item in self.libs:
            files=glob.glob(item+"/libz.so")
            files.extend(glob.glob(item+"/libz.a"))
            files.extend(glob.glob(item+"/libz.dylib"))
            if len(files)!=0:
	        self.configLinux.libraries['ZLib']=files[0]+":"+str(os.stat(files[0]).st_mtime)
   	        find=True
                os.environ['LD_LIBRARY_PATH'] = os.environ['LD_LIBRARY_PATH'] + \
                    		        ":" + item
                os.environ['DYLD_LIBRARY_PATH'] = os.environ['DYLD_LIBRARY_PATH'] + \
                    		        ":" + item
                os.environ['LIBRARY_PATH'] = os.environ['LIBRARY_PATH'] + \
                    		        ":" + item
	        break
        if not find:
	    logging.warning("Library called 'libz' is not found. Gzip format will be disabled.")
            logging.warning("To enable this format, please install 'zlib-devel' package.")
            return False

	# Checking header file filtreing_streambuf.hpp
        find=False
        for item in self.includes:
            files=glob.glob(item+"/zlib.h")
            if len(files)!=0:
  	        self.configLinux.headers['ZLib']=files[0]+":"+str(os.stat(files[0]).st_mtime)
	        find=True
                os.environ['CPLUS_INCLUDE_PATH'] = os.environ['CPLUS_INCLUDE_PATH'] + \
                    		        ":" + item
   	        break
        if not find:
	    logging.warning("Header file called 'zlib.h' is not found. Gzip format will be disabled.")
            logging.warning("To enable this format, please install 'zlib-devel' package.")

            return False

        return True


