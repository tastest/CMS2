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
import madanalysis.dataset.dataset as Dataset

class DatasetCollection:

    def __init__(self):
        self.table = {}

    def __len__(self):
        return len(self.table)

    def __getitem__(self,i):
        return self.table.values()[i]

    def Display(self):
        logging.info(" ********* List of defined datasets *********" )
        sorted_keys = sorted(self.table.keys())
        for key in sorted_keys:
            logging.info(" "+key+" ("+self.table[key].GetStringTag()+")")
        logging.info(" ********************************************" )

    def Find(self,name):
        name.lower()
        if name in self.table.keys():
            return True
        return False

    def Add(self,name):
        name.lower()
        self.table[name]=Dataset.Dataset(name)

    def Get(self,name):
        name.lower()
        return self.table[name]

    def Remove(self,name):
        name.lower()
        if self.Find(name):
            del self.table[name]

    def Reset(self):
        for key in self.table.keys():
            del self.table[key]

    def GetNames(self):
        return sorted(self.table.keys())
