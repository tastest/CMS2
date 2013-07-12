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
import os
import logging
import time

class HTMLReportWriter(TextFileWriter.TextFileWriter):
    """Generate HTML report"""
    
    def __init__(self,filename):
        TextFileWriter.TextFileWriter.__init__(self,filename)
        self.page=[]
        self.section=[]
        self.sectionLevel=[]
        self.bullet=0
        self.table=0
        self.current_col=0
        self.number_col=0
        self.first_cell=True

    @staticmethod    
    def CheckStructure(dirname):
        if not os.path.isdir(dirname):
            return False
        if not os.path.isfile(dirname+'/index.html'):
            return False
        if not os.path.isfile(dirname+'/logo.png'):
            return False
        return True

    def WriteHeader(self):
        self.page.append("<HTML>\n")
        self.page.append("<HEAD><TITLE>MadAnalysis 5 HTML report</TITLE><HEAD>\n")
        self.page.append("<BODY>\n")
    
    def WriteTitle(self,title):
        self.page.append("<CENTER>")
        self.page.append("<table style=\"text-align: left; width: 839px; height: 71px;\" border=\"1\" cellpadding=\"2\" cellspacing=\"2\">\n")
        self.page.append("<tbody>\n<tr>\n<td style=\"vertical-align: top; text-align: center;\">")
        self.page.append("<img style=\"width: 182px; height: 53px;\" alt=\"\"src=\"logo.png\"><BR><SMALL><I>Please <A HREF=\"http://madanalysis.irmp.ucl.ac.be\">visit us.</I></SMALL></A></td>")
        self.page.append("<td style=\"vertical-align: center; text-align: center;\"><big><big><big>"+title+"</big></big></big><br><BR><I>Created by <FONT COLOR=\"#0000CC\">"+str(os.getlogin())+"</FONT></I>.</td>\n")
        self.page.append("</tr>\n</tbody>\n</table>\n")
        self.page.append("</CENTER>")
        self.page.append("<P ALIGN=CENTER><I>"+str(time.strftime("%A, %d %B %Y at  %H:%M:%S"))+"</I>.</P><BR>")
 
    def TableOfContents(self):
        contents = ""
        contents += "<H1>Table of contents.</H1>\n<P>\n"
        for i in range(len(self.section)):
            for j in range(self.sectionLevel[i]):
                contents += '&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;'
            contents += "<A href=\"#"+self.section[i]+"\">"
            contents += self.section[i]+"</A>.<BR>\n"
        return contents
    
    def WriteSpacor(self):
        self.page.append("<HR SIZE=\"4\" WIDTH=350 noshade>\n")

    def WriteVspace(self):
        self.page.append("<SPACER TYPE=vertical>")

    def WriteSubTitle(self,subtitle):
        self.section.append(subtitle)
        self.sectionLevel.append(1)
        self.page.append("<H2 ALIGN=LEFT><A name=\""+subtitle+"\"><U>"+subtitle+"</U></A>.</H2>\n")

    def WriteSubSubTitle(self,subtitle):
        self.section.append(subtitle)
        self.sectionLevel.append(2)
        self.page.append("<H3 ALIGN=LEFT><A name=\""+subtitle+"\"><U>"+subtitle+"</U></A>.</H3>\n")
        
    def WriteText(self,text):
         if self.bullet!=0:
            self.page.append("<LI>")
         text.WriteHTML(self.page)

    def NewLine(self):
        self.page.append("<BR>")

    def OpenBullet(self):
        self.bullet=self.bullet+1
        self.page.append("<UL>\n")

    def CloseBullet(self):
        self.bullet=self.bullet-1
        self.page.append("</UL>\n")
       
    def CreateTable(self,col):
        self.table=self.table+1
        self.number_col=len(col)
        self.page.append("\n<TABLE align=\"center\" frame=\"above\" border=\"5\" cellpadding=\"2\">\n<TR>\n")

    def NewCell(self,color=ColorType.WHITE):
        self.current_col=self.current_col+1
        
        if  self.current_col>self.number_col:
            logging.warning(" the number of the current column is bigger than the total number of declared columns.")
        if self.first_cell==True:
            self.page.append("<TD bgcolor=\""+ColorType.convert2hexa(color)+"\">")
        else:
            self.page.append("</TD> \n<TD bgcolor=\""+ColorType.convert2hexa(color)+"\">")
        self.first_cell=False
       
    def NewLine(self):
        self.current_col=0
        self.first_cell=True
        self.page.append("</TD> \n</TR>\n<TR>\n")
        
    def EndTable(self,caption):
        self.table=self.table-1
        if caption == "":
            self.page.append("</TD>\n</TR>\n</TBALE>\n")
        else :
            self.page.append("</TD>\n</TR>\n<CAPTION ALIGN=\"BOTTOM\">")
            self.WriteText(caption)
            self.page.append("</CAPTION>\n</TABLE>\n")
    

    def WriteFigure(self,caption,filename,scale):
        thefile = os.path.normpath(filename)
        from ROOT import TImage
        im = TImage.Open(thefile+".png",TImage.kPng)
        if not im.IsValid():
            logging.warning(" the picture "+ thefile+".png does not exist.")
        self.page.append("<CENTER>")
        self.page.append("<FIGURE>")
        self.page.append("<IMG ALIGN=\"center\" SRC=\""+os.path.basename(filename)+\
                         ".png\" WIDTH=\""+str(scale*im.GetWidth())+"\" HEIGHT=\""+str(scale*im.GetHeight())+"\">\n<FIGCAPTION>")
        self.WriteText(caption)
        self.page.append("</FIGCAPTION>\n</FIGURE>\n")
        self.page.append("</CENTER>\n")
        
    
    def WriteFoot(self):
        if self.bullet!=0:
            logging.warning(" the number of 'OpenBullet()' and 'CloseBullet()' are different.")
        if self.table!=0:
            logging.warning("open table found. Please check for a missing 'EndTable()'.")
        self.page.append("<HR SIZE=\"4\" WIDTH=500 noshade>\n")
        self.page.append("<CENTER><SMALL><I> Please visit us: </I><A HREF=\"http://madanalysis.irmp.ucl.ac.be\">http://madanalysis.irmp.ucl.ac.be</A>.</SMALL></CENTER>\n")
        self.page.append("</BODY>\n")
        self.page.append("</HTML>\n")
        self.page.insert(11,self.TableOfContents())
        for item in self.page:
            self.file.write(item)
            
