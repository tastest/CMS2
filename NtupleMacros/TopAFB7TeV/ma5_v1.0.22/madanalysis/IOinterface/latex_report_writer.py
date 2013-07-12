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


import madanalysis.IOinterface.text_file_writer as TextFileWriter
from madanalysis.enumeration.color_type import ColorType
from madanalysis.enumeration.font_type import FontType
from madanalysis.enumeration.script_type import ScriptType
from madanalysis.IOinterface.text_report import TextReport

import logging
import time
import os

class LATEXReportWriter(TextFileWriter.TextFileWriter):
    """Generate LaTeX report"""

    def __init__(self,filename,pdflatex=False):
        TextFileWriter.TextFileWriter.__init__(self,filename)
        self.pdflatex=pdflatex
        self.bullet=0
        self.table=0
        self.current_col=0
        self.number_col=0
        self.first_cell=True
        self.ext=''
        if pdflatex==False:
            self.ext='.eps'
        else:
            self.ext='.png'

    @staticmethod    
    def CheckStructure(dirname):
        if not os.path.isdir(dirname):
            return False
        if not os.path.isfile(dirname+'/main.tex'):
            return False
        if not os.path.isfile(dirname+'/logo.png') and \
           not os.path.isfile(dirname+'/logo.eps'):
            return False
        return True

    def WriteHeader(self):
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("%                     DOCUMENT'S CLASS & OPTIONS \n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("\\documentclass[10pt, a4paper, twoside]{article}\n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("%                     PACKAGES \n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("\\usepackage[english]{babel}\n")
        self.file.write("\\usepackage")
        if self.pdflatex==True:
            self.file.write("[pdftex]")
        self.file.write("{graphicx}\n")
        self.file.write("\\usepackage{geometry}\n")
        self.file.write("\\usepackage{verbatim}\n")
        self.file.write("\\usepackage{fancyhdr}\n")
        self.file.write("\\usepackage{multirow}\n")
        self.file.write("\\usepackage{subfigure}\n")
        self.file.write("\\usepackage{colortbl}\n")
        self.file.write("\\usepackage{hyperref}\n")
        self.file.write("\\usepackage{multicol}\n")
        self.file.write("\\usepackage{caption}\n")
        self.file.write("\\usepackage{alltt}\n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("%                      DOCUMENT CONFIGURATION\n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("\\geometry{verbose, a4paper, tmargin=2.5cm, bmargin=2.5cm, lmargin=1.5cm, rmargin=1.5cm}\n")
        self.file.write("\\definecolor{purple}{rgb}{0.62,0.12,0.94}")
        self.file.write("\\definecolor{grey}{rgb}{0.3,0.3,0.3}")
        self.file.write("\\definecolor{orange}{rgb}{1,0.5,0}")
        self.file.write("\\captionsetup{labelformat=empty,aboveskip=0pt}")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("%                      BEGIN DOCUMENT\n")
        self.file.write("%-------------------------------------------------------------------\n")
        self.file.write("\\begin{document}\n")

    def WriteTitle(self,title):
        self.file.write("\\begin{tabular}{cc}\n\\multirow{3}{*}{\\includegraphics[scale=")
        if self.pdflatex==True:
            self.file.write("0.5")
        else:
            self.file.write("0.4")
        self.file.write("]{logo"+self.ext+"}}&\\\\\n&\\textit{\\href{http://madanalysis.irmp.ucl.ac.be}{http://madanalysis.irmp.ucl.ac.be}}\\\\\n&\\\\\n\\end{tabular}\\\\\n")
        self.file.write("\\begin{center}\n\\begin{tabular}{m{18cm}}\n\\centering\n{\\huge{"+title+"}}\\\\\n\\textit{Created by \\textcolor{blue}{"+str(os.getlogin())+"}}\\\\\n\\textit{"+str(time.strftime("%A, %d %B %Y at  %H:%M:%S "))+"}\\\\\n\\end{tabular}\n\\end{center}\n")
        self.file.write("\\tableofcontents\n")

    def WriteSpacor(self):
        self.file.write("\n\\hrulefill\n")

    def WriteVspace(self):
        self.file.write("\\vspace{1cm}\n")

    def WriteSubTitle(self,subtitle):
#        if subtitle=="Histograms":
#           self.file.write("\\section{"+subtitle+"}\n")
#        else:    
        self.file.write("\\newpage \\section{"+subtitle+"}\n")

    def WriteSubSubTitle(self,subsubtitle):
        self.file.write("\\subsubsection{"+subsubtitle+"}")

    def WriteText(self,text):
         if self.bullet!=0:
            self.file.write("\\item")
         text.WriteLATEX(self.file)
                  
#    def NewLine(self):
#        self.file.write("\n\\newline\n")

    def OpenBullet(self):
        self.bullet=self.bullet+1
        self.file.write("\\begin{itemize}\n")

    def CloseBullet(self):
        self.bullet=self.bullet-1
        self.file.write("\\end{itemize}\n")

    def CreateTable(self,col):
        self.table=self.table+1
        self.number_col=len(col)
        self.file.write("\\begin{table}[!h]\n\\center\n\\begin{tabular}{|")
        for item in col:
            self.file.write("m{"+str(item)+"cm}|")
        self.file.write("}\n\\hline\n")

    def NewCell(self,color=ColorType.WHITE):
        self.current_col=self.current_col+1

        if  self.current_col>self.number_col:
            logging.warning("The number of the current column is larger than the total number of declared columns.")
        if self.first_cell==True:
            self.file.write("\\cellcolor{"+ColorType.convert2string(color)+"}")
            self.first_cell=False
        else:
            self.file.write("& \\cellcolor{"+ColorType.convert2string(color)+"}")
            
    def NewLine(self):
        self.current_col=0
        self.first_cell=True
        self.file.write("\\\\\n\\hline\n")
        
    def EndTable(self,caption):
        self.table=self.table-1
        if caption == "":
            self.file.write("\n\\hline\n\\end{tabular}\n")
        else :
            self.file.write("\n\\hline\n\\end{tabular}\n\\caption{")
            self.WriteText(caption)
            self.file.write("}\n")
        self.file.write("\\end{table}\n")


    def WriteFigure(self,caption,filename,scale):
        thefile = os.path.normpath(filename)
        if os.path.isfile(thefile+self.ext):
            if self.pdflatex:
                scale=0.60
            self.file.write("\\begin{figure}[!h]\n\center\n\\includegraphics[scale="+str(scale)+"]{"+\
                            os.path.basename(filename)+self.ext+"}\\\\\n\\caption{")
            self.WriteText(caption)
            self.file.write("}\n\\end{figure}\n")
            self.WriteSpacor()
            self.file.write("\\newpage")
        else:
            logging.warning(thefile+self.ext+" does not exist.")
        
    def WriteFoot(self):
        if self.bullet!=0:
            logging.warning("the number of 'OpenBullet()' and 'CloseBullet()' are different.")
        if self.table!=0:
            logging.warning("there is an open table. Please check for a missing 'EndTable()' tag.")
        self.file.write("\\end{document}\n")
