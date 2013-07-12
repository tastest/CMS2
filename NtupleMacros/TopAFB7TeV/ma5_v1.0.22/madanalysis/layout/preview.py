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


from madanalysis.enumeration.report_format_type import ReportFormatType
from madanalysis.IOinterface.root_file_reader import RootFileReader
from madanalysis.IOinterface.folder_writer import FolderWriter
from madanalysis.selection.instance_name import InstanceName
from madanalysis.layout.plotflow import PlotFlow

from math import log10, floor
import os
import shutil

class Preview:

    def __init__(self,main):
        self.main = main
        self.input_path=self.main.lastjob_name
        self.files=[]
        self.plotflow = PlotFlow(self.main,self.files,"",ReportFormatType.HTML)


    def Open(self):

        # Checking input dir
        if not os.path.isdir(self.input_path):
            logging.error("no directory denoted by '"+self.input_path+"' found.")
            return False

        # Creating list of ROOT files
        for ind in range(0,len(self.main.datasets)):
            name=InstanceName.Get(self.main.datasets[ind].name)
            self.files.append(RootFileReader(os.path.normpath(self.input_path+"/root/"+name+".root")))

        # Trying to open each ROOT files
        for ind in range(0,len(self.files)):
            if not self.files[ind].Open():
                for i in range(0,ind):
                    self.files[i].Close()
                return False

        # Good end 
        return True     


    def DoThePlot(self,index):

        test = False

        # Checking negative index
        if index <= 0:
            logging.error("the index could not be negative or null")
            return False

        # Checking too big index
        if index > len(self.main.selection):
            logging.error("the index is "+str(index)+" whereas the number of plot/cut is : "+str(self.main.selection))
            return False
        
        # Looking for the plot 
        irelhisto=0
        for iabshisto in range(0,len(self.main.selection)):
            if self.main.selection[iabshisto].__class__.__name__!="Histogram":
                continue
            if index==(irelhisto+1):
                test=True
                break
            irelhisto+=1
        if not test:
            logging.error("the selection item is not a plot but a cut")
            return False

        # Initialize the plot flow
        if not self.plotflow.initialize():
            return False

        # Draw the corresponding histo
        self.plotflow.Preview(iabshisto,irelhisto)
        
        return True

                
