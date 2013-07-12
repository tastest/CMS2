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


from madanalysis.enumeration.sb_ratio_type       import SBratioType
from madanalysis.enumeration.color_type          import ColorType
from madanalysis.IOinterface.root_file_reader    import RootFileReader
from madanalysis.IOinterface.folder_writer       import FolderWriter
from madanalysis.selection.instance_name         import InstanceName
from madanalysis.enumeration.font_type           import FontType
from madanalysis.enumeration.script_type         import ScriptType
from madanalysis.IOinterface.text_report         import TextReport
from madanalysis.IOinterface.html_report_writer  import HTMLReportWriter
from madanalysis.IOinterface.latex_report_writer import LATEXReportWriter
from madanalysis.enumeration.report_format_type  import ReportFormatType
from madanalysis.enumeration.normalize_type      import NormalizeType
from madanalysis.enumeration.observable_type     import ObservableType
from madanalysis.layout.cutflow                  import CutFlow
from madanalysis.layout.plotflow                 import PlotFlow
from math                                        import log10, floor, ceil
import os
import shutil
import logging

class Layout:

    def __init__(self,main,path,mode):
        self.main=main
        self.input_path=self.main.lastjob_name
        self.output_path=path
        self.files=[]
        self.plots=[]
        self.mode=mode
        self.color=1
        self.cutflow  = [] 
        self.plotflow = PlotFlow(self.main,self.files,self.output_path,self.mode)
        self.isSignal = False
        self.isBackground = False

    @staticmethod
    def DisplayInteger(value):
        if type(value) is not int:
            return ""
        if value<0:
            return "-" + Layout.DisplayInteger(-value)
        elif value < 1000:
            return str(value)
        else:
            return Layout.DisplayInteger(value / 1000) +\
                   "," + '%03d' % (value % 1000)

    @staticmethod
    def Round_to_Ndigits(x,N):
        if N<1:
            return ""
        if x<(10**(N-1)):
            convert = '%.'+str(N)+'G'
            return '%s' % float(convert % x)
        else:
            tmp = '%s' % float('%.12G' % int(x))
            if len(tmp)>=3 and tmp.endswith('.0'):
                tmp = tmp[:-2]
            return tmp    

    @staticmethod
    def DisplayXsection(xsection,xerror):
        # xsection and xerror are null
        if xsection==0. and xerror==0.:
            return "0.0 +/- 0.0"
            
        # xsection is not null but xerror is
        # keep the 3 significative digit
        elif xerror==0:
            return Layout.Round_to_Ndigits(xsection,3)
        
        # error greater than xsection ? 
        elif xsection > xerror:
            string1 = Layout.Round_to_Ndigits(xerror,3)
            if 'e' in string1 or 'E' in string1:
                string2='%e' % xsection
            elif '.' in string1:
                convert='%.'+str(len(string1.split('.')[1]))+'f'
                string2=convert % xsection
            else:
                string2=str(int(xsection))
            return string2 + " +/- " + string1    

        else:
            string1 = Layout.Round_to_Ndigits(xsection,3)
            if 'e' in string1 or 'E' in string1:
                string2='%e' % xerror
            elif '.' in string1:
                convert='%.'+str(len(string1.split('.')[1]))+'f'
                string2=convert % xerror
            else:
                string2=str(int(xerror))
            return string1 + " +/- " + string2    

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

        # Creating production directory
        if not FolderWriter.CreateDirectory(self.output_path,True):
            return False

        # Creating cut flow for each file
        for ind in range(0,len(self.files)):
            self.cutflow.append(CutFlow(self.main.datasets[ind],\
                                        self.main.selection,\
                                        self.main.lumi,
                                        self.main))        
        # Good end 
        return True     


    def DoEfficiencies(self):

        if self.main.selection.Ncuts==0:
            return True

        for i in range(0,len(self.cutflow)):
            if not self.cutflow[i].initializeFromFile(self.files[i]):
                return False

            self.cutflow[i].calculate()

        self.signal     = CutFlow(self.main.datasets[0],\
                                  self.main.selection,\
                                  self.main.lumi,\
                                  self.main)
        self.background = CutFlow(self.main.datasets[0],\
                                  self.main.selection,\
                                  self.main.lumi,\
                                  self.main)

        signalvect     = []
        backgroundvect = []  

        for i in range(0,len(self.main.datasets)):
            if not self.main.datasets[i].background:
                signalvect.append(self.cutflow[i])
            else:
                backgroundvect.append(self.cutflow[i])

        if len(signalvect)!=0:
            self.signal.initializeFromCutflow(signalvect)
            self.isSignal=True
        if len(backgroundvect)!=0:
            self.background.initializeFromCutflow(backgroundvect)
            self.isBackground=True

        return True

    def DoPlots(self):

        if self.main.selection.Nhistos==0:
            return True

        if not self.plotflow.initialize():
            return False

        self.plotflow.DrawAll()
        
        return True

                

    def CopyLogo(self):
        
        # Filename
        filename = self.main.ma5dir+"/madanalysis/input/" + \
                   "logo." + \
                   ReportFormatType.convert2filetype(self.mode)

        # Checking file presence
        if not os.path.isfile(filename):
            logging.error("the image '" + \
                          filename + \
                          "' is not found.")
            return False

        # Copy file
        try :
            shutil.copy(filename,self.output_path)
            return True
        except:
            logging.error("Errors have occured during the copy of the file ")
            logging.error(" "+filename)
            return False

    def WriteDatasetTable(self,rootfile,report,xsection,weight):

        filenames = rootfile.Get("general/filenames","TClonesArray")
        if filenames is None:
            return False
                    
        xsections = rootfile.Get("general/xsections","TVectorT<float>")
        if xsections is None:
            return False
        
        xerrors   = rootfile.Get("general/xerrors","TVectorT<float>")
        if xerrors is None:
            return False
        
        nevents   = rootfile.Get("general/nevents","TVectorT<float>")
        if nevents is None:
            return False

        if filenames.GetEntries()!=xsections.GetNoElements() or \
           xsections.GetNoElements()!=xerrors.GetNoElements() or \
           xerrors.GetNoElements()!=nevents.GetNoElements() or \
           nevents.GetNoElements()!=filenames.GetEntries() :
            logging.error("the 'general' branches have different size "\
                          "in the file '"+filename+"'")
            return
        
       
        text=TextReport()

        text.Add('* Generation: ')
        text.SetColor(ColorType.BLUE)
        ngen = int(nevents[filenames.GetEntries()-1]);
        text.Add(str(ngen) + " ")
        text.SetColor(ColorType.BLACK)
        text.Add(' events.\n')
        if xsection != 0.0:
            text.Add('* Cross section imposed by the user: ')
            text.SetColor(ColorType.BLUE)
            text.Add(str(xsection))
            text.SetColor(ColorType.BLACK)
            text.Add(' pb.\n')
        if weight != 1.0:
            text.Add('* Event weight imposed by the user: ')
            text.SetColor(ColorType.BLUE)
            text.Add(str(weight))
            text.SetColor(ColorType.BLACK)
            text.Add('.\n')
        text.Add('* Normalization to ')
        text.Add(str(self.main.lumi))
        text.Add(' fb')
        text.SetScript(ScriptType.SUP)
        text.Add('-1')
        text.SetScript(ScriptType.none)
        text.Add(': ')
        text.SetColor(ColorType.BLUE)
        if xsection != 0.0:
            nlumi = int(xsection*1000*self.main.lumi)
        else:
            nlumi = int(xsections[filenames.GetEntries()-1]*1000*self.main.lumi*weight)
        text.Add(str(nlumi))
        text.Add(' +/- ')

        # round to the smallest integer greater than error
        if xsection != 0.0:
            elumi = 0.0
        else:
            elumi = ceil(xerrors[filenames.GetEntries()-1]*1000*self.main.lumi*weight)
        text.Add(str(int(elumi))+ " ")
        text.SetColor(ColorType.BLACK)
        text.Add(' events.\n')
        evw = float(nlumi)/float(ngen)*weight
        if evw > 1:
            text.SetColor(ColorType.RED)
        text.Add('* Ratio (event weight) =  ')
        if evw < 1:
            text.SetColor(ColorType.BLUE)
        text.Add(str(Layout.Round_to_Ndigits(evw,2)) + " ")
        if evw < 1:
            text.SetColor(ColorType.BLACK)
            text.Add('.')
        if evw > 1:
            text.Add(' - warning: please generate more events (weight larger than 1)!\n')
            text.SetColor(ColorType.BLACK)
        else:
            text.Add(' \n')
        if self.mode is ReportFormatType.HTML:
            text.Add(' \n')
        report.WriteText(text)


        report.CreateTable([11.5,2,3])
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        if filenames.GetEntries()>1:
            text.Add("Event files")
        else:
            text.Add("Event file")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Number of events")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Cross section (pb)")
        report.WriteText(text)
        report.NewLine()

        for ind in range(0,filenames.GetEntries()):
            color = ColorType.BLACK
            if ind is filenames.GetEntries()-1:
                color = ColorType.BLUE
            
            report.NewCell()
            text.Reset()
            text.SetColor(color)
            text.Add(str(filenames[ind]))
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            text.SetColor(color)
            text.Add(str(int(nevents[ind])))
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            text.SetColor(color)
            if xsection != 0.0:
                text.Add(str(xsection))
            else:
                text.Add(Layout.DisplayXsection(xsections[ind]*weight,xerrors[ind]*weight))
            report.WriteText(text)
            report.NewLine()
        text.Reset()
        report.EndTable(text)    
        text.Reset()

    # Writing Final Table
    def WriteFinalTable(self,report):

        # Caption
        text=TextReport()
        text.Reset()
        text.Add("Formula for signal(S)-background(B) comparison : ")
        text.SetColor(ColorType.BLUE)
        text.Add(self.main.SBratio+'\n')
        text.SetColor(ColorType.BLACK)
        report.WriteText(text)
        text.Reset()
        text.Add("Formula for uncertainty on signal(S)-background(B) comparison : ")
        text.SetColor(ColorType.BLUE)
        text.Add(self.main.SBerror+'\n')
        text.SetColor(ColorType.BLACK)
        report.WriteText(text)
        text.Reset()
        text.Add(' \n')
        report.WriteText(text)

        # Caption
        report.CreateTable([2.6,2.5,3.6,3.6,2.1])
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Cuts")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Signal (S)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Background (B)")# (+/- err)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add(" S vs B")
        report.WriteText(text)
        report.NewLine()
        text.Reset()

        # Initial
        report.NewCell()
        text.Reset()
        text.Add("Initial (no cut)")
        report.WriteText(text)
        report.NewCell()
        text.Reset()
        if self.isSignal:
            text.Add(Layout.DisplayXsection(self.signal.Ntotal.mean,\
                                            self.signal.Ntotal.error))
        else:
            text.Add("")
        report.WriteText(text)
        report.NewCell()
        text.Reset()
        if self.isBackground:
            text.Add(Layout.DisplayXsection(self.background.Ntotal.mean,\
                                            self.background.Ntotal.error))
        else:
            text.Add("")
        report.WriteText(text)
        report.NewCell()
        text.Reset()
        if self.isSignal and self.isBackground:
            value = self.signal.calculateBSratio(\
                self.background.Ntotal.mean,\
                self.background.Ntotal.error,\
                self.signal.Ntotal.mean,\
                self.signal.Ntotal.error)
            text.Add(Layout.DisplayXsection(value.mean,value.error))
        else:
            text.Add("")
        report.WriteText(text)
        report.NewLine()
        text.Reset()

        # Loop
        for ind in range(0,len(self.cutflow[0].Nselected)):
            report.NewCell()
            text.Reset()
            text.Add("cut " + str(ind+1))
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            if self.isSignal:
                text.Add(Layout.DisplayXsection(self.signal.Nselected[ind].mean,\
                                                self.signal.Nselected[ind].error))
            else:
                text.Add("")
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            if self.isBackground:
                text.Add(Layout.DisplayXsection(self.background.Nselected[ind].mean,\
                                                self.background.Nselected[ind].error))
            else:
                text.Add("")
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            if self.isSignal and self.isBackground:
                value = self.cutflow[0].calculateBSratio(\
                    self.background.Nselected[ind].mean,\
                    self.background.Nselected[ind].error,\
                    self.signal.Nselected[ind].mean,\
                    self.signal.Nselected[ind].error)
                text.Add(Layout.DisplayXsection(value.mean,value.error))
            else:
                text.Add("")
            report.WriteText(text)
            report.NewLine()
            text.Reset()


        text.Add("Signal and Background comparison")
        report.EndTable(text)    


    # Writing Efficiency Table
    def WriteEfficiencyTable(self,index,report):
        
        text=TextReport()
        report.CreateTable([2.6,2.5,3.6,3.6,2.1])
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Dataset")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Selected events (S)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Rejected events (R)")# (+/- err)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("S / (S + R)")# (+/- err)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("S / Initial")
        report.WriteText(text)
        report.NewLine()
        text.Reset()
        report.NewLine()
        for i in range(0,len(self.main.datasets)):
            # DatasetName
            report.NewCell()
            text.Reset()
            text.Add(self.main.datasets[i].name)
            report.WriteText(text)

            # SelectedEvents
            report.NewCell()
            text.Reset()
            text.Add(Layout.DisplayXsection(self.cutflow[i].Nselected[index].mean,self.cutflow[i].Nselected[index].error))
            report.WriteText(text)
            
            # RejectedEvents
            report.NewCell()
            text.Reset()
            text.Add(Layout.DisplayXsection(self.cutflow[i].Nrejected[index].mean,self.cutflow[i].Nrejected[index].error)) 
            report.WriteText(text)

            # Efficiency Events
            report.NewCell()
            text.Reset()
            text.Add(Layout.DisplayXsection(self.cutflow[i].eff[index].mean,self.cutflow[i].eff[index].error))
            report.WriteText(text)

            # Cumulative efficiency events
            report.NewCell()
            text.Reset()
            text.Add(Layout.DisplayXsection(self.cutflow[i].effcumu[index].mean,self.cutflow[i].effcumu[index].error))
            report.WriteText(text)
            
            report.NewLine()
            
        text.Reset()
        report.EndTable(text)    


    # Writing Statistics Table
    def WriteStatisticsTable(self,index,report):
        
        text=TextReport()
        report.CreateTable([2.6,2.5,3.6,3.6,3.6,2.1,2.1])
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Dataset")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Integral")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Entries / events")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Mean")# (+/- err)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("RMS")# (+/- err)")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Underflow")
        report.WriteText(text)
        report.NewCell(ColorType.YELLOW)
        text.Reset()
        text.Add("Overflow")
        report.WriteText(text)
        report.NewLine()
        
        histos=[]
        nevents=[]
        nentries=[]
        for iset in range(0,len(self.plotflow.plots)):
            histos.append(self.plotflow.plots[iset][index])
            nevents.append(self.plotflow.nevents[iset][index])
            nentries.append(self.plotflow.nentries[iset][index])
            
        # Looping on each histogram
            
        for item in range(0,len(histos)):
            
            # Getting the number of entries
            integral=0.
            nbinx = histos[item].GetNbinsX()+2
            for bin in range(0,nbinx):
                integral += histos[item].GetBinContent(bin)

            # Getting the number of events
            nevent = nevents[item]
                
            # Getting underflow and overflow - percentage    
            uflow = histos[item].GetBinContent(0)
            oflow  = histos[item].GetBinContent(nbinx)
           
            uflow_percent=0
            oflow_percent=0
            if integral!=0:
                uflow_percent = uflow*100/integral
                oflow_percent = oflow*100/integral
                
            # mean value + error
            mean=histos[item].GetMean(1)
            mean_error=histos[item].GetMeanError(1)
                
            # root mean square + error 
            rms=histos[item].GetRMS()
            rms_error=histos[item].GetRMSError()
                
            # writing the table
            report.NewCell()
            text.Reset()
            text.Add(self.main.datasets[item].name)
            report.WriteText(text)
            report.NewCell()

            # Nentries
            text.Reset()
            text.Add(Layout.DisplayXsection(integral,0))
            report.WriteText(text)
            report.NewCell()

            # Nentries / Nevents
            text.Reset()
            if nevents[item] !=0.:
                text.Add(str(Layout.Round_to_Ndigits(float(nentries[item])/float(nevents[item]),3)))
            else:
                text.Add("0.")
            report.WriteText(text)
            report.NewCell()

            # Mean value
            text.Reset()
            text.Add(str(Layout.Round_to_Ndigits(mean,6)))#+"(+/- "+str(Layout.Round_to_Ndigits(mean_error,3))+")")
            report.WriteText(text)
            report.NewCell()
            text.Reset()
            text.Add(str(Layout.Round_to_Ndigits(rms,4)))#+"(+/- "+str(Layout.Round_to_Ndigits(rms_error,3))+")")
            report.WriteText(text)
            if uflow_percent+oflow_percent<=5:
                report.NewCell(ColorType.GREEN)
            if uflow_percent+oflow_percent>5 and uflow_percent+oflow_percent<15:
                report.NewCell(ColorType.ORANGE)
            if uflow_percent+oflow_percent>15:
                report.NewCell(ColorType.RED)
            text.Reset()
            text.Add(str(Layout.Round_to_Ndigits(uflow_percent,4)))
            report.WriteText(text)
            if uflow_percent+oflow_percent<=5:
                report.NewCell(ColorType.GREEN)
            if uflow_percent+oflow_percent>5 and uflow_percent+oflow_percent<15:
                report.NewCell(ColorType.ORANGE)
            if uflow_percent+oflow_percent>15:
                report.NewCell(ColorType.RED)
            text.Reset()
            text.Add(str(Layout.Round_to_Ndigits(oflow_percent,4)))
            report.WriteText(text)
            report.NewLine()
        text.Reset()
        text.Add("Histogram number "+str(index+1)+" - Statistics")
        report.EndTable(text)

            
    def GenerateReport(self,history):

        # Defining report writing
        if self.mode is ReportFormatType.HTML:
            report = HTMLReportWriter(self.output_path+"/index.html")
        elif self.mode is ReportFormatType.LATEX:
            report = LATEXReportWriter(self.output_path+"/main.tex",False)
        else :
            report = LATEXReportWriter(self.output_path+"/main.tex",True)
                
        # Opening
        if not report.Open():
            return False

        # Create text
        text=TextReport()

        # Header
        report.WriteHeader()
        report.WriteTitle('MadAnalysis 5 report')
        report.WriteSpacor()

        # History of commands
        report.WriteSubTitle('Command history')
        text.Reset()
        text.SetFont(FontType.TT)
        for item in history:
            text.Add('ma5>'+ item+'\n')
        report.WriteText(text)

        # Configuration
        report.WriteSubTitle('Configuration')

        # Integrated luminosity 
        text.Reset()
        text.Add('MadAnalysis version ' + self.main.version + \
                 ' (' + self.main.date + ').\n')
        report.WriteText(text)

        # Integrated luminosity 
        text.Reset()

        # Normalization
        if self.main.normalize == NormalizeType.LUMI or \
           self.main.normalize == NormalizeType.LUMI_WEIGHT:
            text.Add('Histograms correspond to an integrated luminosity of ')
            text.SetColor(ColorType.BLUE)
            text.Add(str(self.main.lumi))
            text.SetColor(ColorType.BLUE)
            text.Add(' fb')
            text.SetScript(ScriptType.SUP)
            text.Add('-1')
            text.SetScript(ScriptType.none)
            text.Add('.\n')
        elif self.main.normalize == NormalizeType.NONE:
            text.Add('Histograms are not scaled.\n')
        
        report.WriteText(text)

        # Datasets
        report.WriteSubTitle('Datasets used')
        for ind in range(0,len(self.main.datasets)):
            datatype="signal"
            if self.main.datasets[ind].background:
                datatype="background"
            report.WriteSubSubTitle(
                self.main.datasets[ind].name + \
                ' (' + datatype + ')' )
            self.WriteDatasetTable(self.files[ind],report,self.main.datasets[ind].xsection,\
                                   self.main.datasets[ind].weight)
           

        # Plots display
        report.WriteSubTitle('Histograms / Cuts')
        
        # Plots
        ihisto=0
        icut=0
        for ind in range(0,len(self.main.selection)):
            if self.main.selection[ind].__class__.__name__=="Histogram":
                report.WriteSubSubTitle("Histogram number "+str(ihisto+1))
                if self.main.selection[ind].observable.name not in ['NPID','NAPID']:
                    self.WriteStatisticsTable(ihisto,report)
                text.Reset()
                text.Add("Histogram number "+str(ihisto+1))
                report.WriteFigure(text,self.output_path +"/"+
                           'selection_'+str(ihisto),\
                           1.0)
                ihisto+=1
            if self.main.selection[ind].__class__.__name__=="Cut":
                report.WriteSubSubTitle("Cut number "+str(icut+1))
                text.Reset()
                text.Add(self.main.selection[ind].GetStringDisplay()+'\n')
                report.WriteText(text)
                self.WriteEfficiencyTable(icut,report)
                icut+=1

        # Final table
        if self.main.selection.Ncuts!=0:
            report.WriteSubTitle('Signal and Background comparison')
            self.WriteFinalTable(report)
            
        # Foot
        report.WriteFoot()

        # Closing
        report.Close()

        return True
        

    def Close(self):
        for item in self.files:
            item.Close()

    @staticmethod
    def CheckLatexLog(file):
        if not os.path.isfile(file):
            return False
        for line in file:
            if line.startswith('!'):
                return False
        return True

    def CompileReport(self):
        
        # ---- LATEX MODE ----
        if self.mode==ReportFormatType.LATEX:

            # Launching latex and producing DVI file
            os.system('cd '+self.output_path+'; latex main.tex -interaction=nonstopmode > latex.log 2>&1')
            name=os.path.normpath(self.output_path+'/main.dvi')
            if not os.path.isfile(name):
                logging.error('DVI file cannot be produced')
                logging.error('Please have a look to the log file '+self.output_path+'/latex.log')
                return False
            
            # Checking latex log : are there errors
            if not Layout.CheckLatexLog(self.output_path+'/latex.log'):
                logging.error('some errors occured during LATEX compilation')
                logging.error('for more details, have a look to the log file : '+self.output_path+'/latex.log')
                return False
                
            # Converting DVI file to PDF file
            os.system('cd '+self.output_path+'; dvipdf main.dvi > dvipdf.log 2>&1')
            name=os.path.normpath(self.output_path+'/main.pdf')

            # Checking PDF file presence
            if not os.path.isfile(name):
                logging.error('PDF file cannot be produced')
                logging.error('Please have a look to the log file '+self.output_path+'/dvipdf.log')
                return False
                
        # ---- PDFLATEX MODE ----
        elif self.mode==ReportFormatType.PDFLATEX:

            # Launching latex and producing PDF file
            os.system('cd '+self.output_path+'; pdflatex -interaction=nonstopmode main.tex > latex.log 2>&1');

            # Checking latex log : are there errors
            if not Layout.CheckLatexLog(self.output_path+'/latex.log'):
                logging.error('some errors occured during LATEX compilation')
                logging.error('for more details, have a look to the log file : '+self.output_path+'/latex.log')
                return False
            
            # Checking PDF file presence
            name=os.path.normpath(self.output_path+'/main.pdf')
            if not os.path.isfile(name):
                logging.error('PDF file cannot be produced')
                logging.error('Please have a look to the log file '+self.output_path+'/latex2.log')
                return False
            
        
            
            
        
    
