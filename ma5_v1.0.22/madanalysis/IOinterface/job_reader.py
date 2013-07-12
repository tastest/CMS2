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
import logging
import shutil
import os
import commands

class JobReader():

    def __init__(self,jobdir):
        self.path       = jobdir
        self.rootdir    = os.path.normpath(self.path+"/root")

    def Open(self):
        if not os.path.isdir(self.path):
            logging.error("Directory called '"+self.path+"' is not found.")
            return False
        elif not os.path.isdir(self.rootdir):
            logging.error("Directory called '"+self.rootdir+"' is not found.")
            return False
        else:
            return True
                            
    def CheckRootFile(self,dataset):
        name=InstanceName.Get(dataset.name)
        if os.path.isfile(self.rootdir+"/"+name+".root"):
            return True
        else:
            logging.error("File called '"+self.rootdir+"/"+name+".root' is not found.")
            return False

    def ErrorMsg_BranchNotFound(self,branchname,filename):
        logging.error("branch '"+branchname+"' is not found in the file '"\
                      +filename+"'")

    def ErrorMsg_BranchEmpty(self,branchname,filename):
        logging.error("branch '"+branchname+"' is not empty in the file '"\
                      +filename+"'")

    def ExtractDatasetInfo(self,dataset):
        from ROOT import TFile
        name=InstanceName.Get(dataset.name)
        filename = self.rootdir+"/"+name+".root"
        rootfile = TFile(filename)
        if rootfile.IsZombie():
            logging.error("file called '"+self.rootdir+"/"+name+".root is not found")
            return


        # Getting data from ROOT file
        xsections = rootfile.Get("general/xsections")
        if not bool(xsections):
            ErrorMsg_BranchNotFound('general/xsections',filename)
            return

        xerrors = rootfile.Get("general/xerrors")
        if not bool(xerrors):
            ErrorMsg_BranchNotFound('general/xerrors',filename)
            return

        nevents = rootfile.Get("general/nevents")
        if not bool(nevents):
            ErrorMsg_BranchNotFound('general/nevents',filename)
            return


        # Checking indices

        if xsections.GetNoElements() is 0:
            ErrorMsg_BranchEmpty("branch 'general/xsections' is empty")
            return

        if xerrors.GetNoElements() is 0:
            ErrorMsg_BranchEmpty("branch 'general/xerrors' is empty")
            return

        if nevents.GetNoElements() is 0:
            ErrorMsg_BranchEmpty("branch 'general/nevents' is empty")
            return

        if xsections.GetNoElements()!=xerrors.GetNoElements() or \
           xerrors.GetNoElements()!=nevents.GetNoElements() or \
           nevents.GetNoElements()!=xsections.GetNoElements():
            logging.error("the 'general' branches have different size "\
                          "in the file '"+filename+"'")
            return

        if xsections.GetNoElements() is not (len(dataset)+1):
            logging.error("number of data files do not correspond in the file '"+filename+"'")
            return
        
        # Extracting data
        dataset.measured_xsection=xsections[xsections.GetNoElements()-1]
        dataset.measured_xerror=xerrors[xerrors.GetNoElements()-1]
        dataset.measured_n=int(nevents[nevents.GetNoElements()-1])
