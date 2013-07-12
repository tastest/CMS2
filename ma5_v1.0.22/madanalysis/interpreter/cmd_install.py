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


from madanalysis.interpreter.cmd_base import CmdBase
import logging
import os
import sys

class CmdInstall(CmdBase):
    """Command INSTALL"""

    def __init__(self,main):
        CmdBase.__init__(self,main,"install")

    @staticmethod
    def convert_bytes(bytes):
        bytes = float(bytes)
        if bytes >= 1099511627776:
            terabytes = bytes / 1099511627776
            size = '%.2fT' % terabytes
        elif bytes >= 1073741824:
            gigabytes = bytes / 1073741824
            size = '%.2fG' % gigabytes
        elif bytes >= 1048576:
            megabytes = bytes / 1048576
            size = '%.2fM' % megabytes
        elif bytes >= 1024:
            kilobytes = bytes / 1024
            size = '%.2fK' % kilobytes
        else:
            size = '%.2fb' % bytes
        return size

    @staticmethod
    def reporthook(numblocks,blocksize,filesize):
        try:
            step = int(filesize/(blocksize*10))
        except:
            step = 1

        if (numblocks+1)%step!=0:
            return
        try:
            percent = min(((numblocks+1)*blocksize*100)/filesize, 100)
        except:
            percent = 100
        theString="% 3.1f%%" % percent
        logging.info( theString + " of " + CmdInstall.convert_bytes(filesize) )
                            
    def do(self,args):

        # Checking argument number
        if len(args) != 1:
            logging.error("wrong number of arguments for the command 'install'.")
            self.help()
            return

        # Checking that arguments is 'samples'
        if args[0]!='samples':
            logging.error("the syntax is not correct.")
            self.help()
            return
        
        # Calling selection method
        return self.install_samples()

    def help(self):
        logging.info("   Syntax: install <component>")
        logging.info("   Download and install a MadAnalysis component from the official site.")
        logging.info("   List of available componentes : samples")

    def install_samples(self):

        # Vertical space
        logging.info("")

        # List of files
        files = { "ttbar_fh.lhe.gz" :   "http://madanalysis.irmp.ucl.ac.be/raw-attachment/wiki/samples/ttbar_fh.lhe.gz",\
                  "ttbar_sl_1.lhe.gz" : "http://madanalysis.irmp.ucl.ac.be/raw-attachment/wiki/samples/ttbar_sl_1.lhe.gz",\
                  "ttbar_sl_2.lhe.gz" : "http://madanalysis.irmp.ucl.ac.be/raw-attachment/wiki/samples/ttbar_sl_2.lhe.gz",\
                  "zz.lhe.gz" :         "http://madanalysis.irmp.ucl.ac.be/raw-attachment/wiki/samples/zz.lhe.gz" }

        # Checking connection with MA5 web site
        logging.info("Testing the access to MadAnalysis 5 website ...")
        import urllib
        try:
            urllib.urlopen('http://madanalysis.irmp.ucl.ac.be')
        except:
            logging.error("impossible to access MadAnalysis 5 website.")
            return False
    
        # Checking if the directory exits
        if os.path.isdir(self.main.ma5dir + '/samples'):
            logging.info("'samples' folder is already created")
        else:
            logging.info("Creating the 'samples' folder ...")
            try:
                os.mkdir(self.main.ma5dir + '/samples')
            except:
                logging.error("impossible to create the folder 'samples'")
                return False

        # Removing files
        for file in files.keys():

            if os.path.isfile(self.main.ma5dir + '/samples/' + file):
                logging.info("Removing file '" + file + "' ...")
                try:
                    os.remove(self.main.ma5dir + '/samples/' + file)
                except:
                    logging.error("impossible to remove the file 'samples/'"+file)
                    return False

        # Launching wget
        ind=0
        error=False
        
        for file,url in files.items():
            ind+=1
            logging.info(str(ind)+"/"+str(len(files.keys()))+" Downloading the file '"+file+"' ...")
            
#            try:
            urllib.urlretrieve(url,self.main.ma5dir+'/samples/'+file,CmdInstall.reporthook)
#            except:
#                logging.error("impossible to download '"+file+"'")
#                error=True

        # Result
        if error:
            logging.warning("Error(s) occured during the installation.")
        else:
            logging.info("Installation complete.")

        # Vertical space
        logging.info("")

        return True

    def complete(self,text,args,begidx,endidx):

        nargs = len(args)
        if not text:
            nargs +=1

        if nargs>2:
            return []
        else:
            output = ["samples"]
            return self.finalize_complete(text,output)
    


