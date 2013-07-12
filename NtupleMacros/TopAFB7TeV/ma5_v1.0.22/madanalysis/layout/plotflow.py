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


from madanalysis.enumeration.uncertainty_type     import UncertaintyType
from madanalysis.enumeration.normalize_type       import NormalizeType
from madanalysis.layout.root_config               import RootConfig
from madanalysis.enumeration.report_format_type   import ReportFormatType
from madanalysis.enumeration.observable_type      import ObservableType
from madanalysis.enumeration.color_type           import ColorType
from madanalysis.enumeration.linestyle_type       import LineStyleType
from madanalysis.enumeration.backstyle_type       import BackStyleType
from madanalysis.enumeration.stacking_method_type import StackingMethodType
from math import *


class PlotFlow:

    def __init__(self,main,files,output_path,mode):
        self.main        = main
        self.mode        = mode
        self.files       = files
        self.output_path = output_path
        self.plots       = []
        self.nevents     = []
        self.nentries    = []
        self.scales      = []
 
    def initialize(self):

        self.plots    = []
        self.nevents  = []
        self.nentries = []
        self.scales   = []

        # Getting histos from ROOT file    
        for item in self.files:
            myobject = item.Get("plots/plots_array","TClonesArray")
            if myobject is None:
                return False
            self.plots.append(myobject)
            myobject2 = item.Get("plots/nevents","TVectorT<float>")
            if myobject2 is None:
                return False
            self.nevents.append(myobject2)

            nentries=[]
            for iplot in range (0,len(myobject)):
                nentries.append(myobject[iplot].GetEntries())
            self.nentries.append(nentries)

        # Computing scales
        for iset in range(0,len(self.plots)):

            # Getting xsection
            xsection=self.main.datasets[iset].measured_xsection
            if self.main.datasets[iset].xsection!=0.:
                xsection=self.main.datasets[iset].xsection

            # No scale if no event
            if self.main.datasets[iset].measured_n==0:
                continue

            # Loop over plot
            scales=[]
            iplot=0
            for iabshisto in range(0,len(self.main.selection)):
                if self.main.selection[iabshisto].__class__.__name__!="Histogram":
                    continue
            
                scale=0.
                
                # Scaling histos
                scale2one = False
                if self.main.selection[iabshisto].stack==StackingMethodType.NORMALIZE2ONE or \
                   (self.main.stack==StackingMethodType.NORMALIZE2ONE and \
                   self.main.selection[iabshisto].stack==StackingMethodType.AUTO):
                    scale2one = True
                
                if scale2one or self.main.normalize == NormalizeType.NONE:
                    # Getting the number of entries
                    integral=0.
                    nbinx = self.plots[iset][iplot].GetNbinsX()+2
                    for bin in range(0,nbinx):
                        integral += self.plots[iset][iplot].GetBinContent(bin)
                    if scale2one and integral!=0.:
                        scale = 1./integral
                    elif self.main.normalize == NormalizeType.NONE:
                        scale = 1. 
                        
                elif self.main.normalize == NormalizeType.LUMI:
                    scale = xsection * self.main.lumi * 1000 / \
                            float(self.main.datasets[iset].measured_n)
                    
                elif self.main.normalize == NormalizeType.LUMI_WEIGHT:
                    scale = xsection * self.main.lumi * 1000 * \
                            self.main.datasets[iset].weight / \
                            float(self.main.datasets[iset].measured_n)

                scales.append(scale)

                iplot+=1

            self.scales.append(scales)

        # Reset Configuration
        RootConfig.Init()
        
        return True    

    def DrawAll(self):

        # Loop on each histo type
        irelhisto=0
        for iabshisto in range(0,len(self.main.selection)):
            if self.main.selection[iabshisto].__class__.__name__!="Histogram":
                continue
            self.color=1
            histos=[]
            scales=[]
            for iset in range(0,len(self.plots)):
                histos.append(self.plots[iset][irelhisto])
                scales.append(self.scales[iset][irelhisto])

            # Draw
            self.Draw(histos,scales,self.main.selection[iabshisto],irelhisto,preview=False)
                
            irelhisto+=1

        
    def Preview(self,iabs,irel):
        from ROOT import gROOT
        gROOT.SetBatch(False)
#        from ROOT import TCanvas
#        from ROOT import TApplication
#        from ROOT import gApplication
#        TApplication.NeedGraphicsLibs()
#        gApplication.InitializeGraphics()
#        eric = TCanvas("eric","eric",600,600)
        
#        answer=raw_input("Answer : ")

        self.color=1
        histos=[]
        ns=[]
        for iset in range(0,len(self.plots)):
            histos.append(self.plots[iset][irel])
            scales.append(self.scales[iset][irel])

        self.Draw(histos,scales,self.main.selection[iabs],irel,preview=True)

        gROOT.SetBatch(True)

    diconicetitle = {' ^ {':'^{', ' _ {':'_{', '\\\\':'#'}

    def NiceTitle(self,text):
        newtext=text 
        for i,j in self.diconicetitle.iteritems():
           newtext = newtext.replace(i,j)
        return newtext
        
    def Draw(self,histos,scales,ref,irelhisto,preview=False):

        from ROOT import TH1
        from ROOT import TH1F
        from ROOT import THStack
        from ROOT import TLegend
        from ROOT import TCanvas

        # Creating a canvas
        canvas = TCanvas("tempo","")

        # Special case : NPID, NAPID
        if ref.observable.name in ['NPID','NAPID']:

            # New collection with labels
            labels=[]
            
            # Loop over histos
            for ind in range(len(histos)):
                nBins = histos[ind].GetNbinsX()
                for b in range(1,nBins+1):
                    # Get the bin label
                    theLabel = histos[ind].GetXaxis().GetBinLabel(b)
                    # Add in the collection 
                    if theLabel not in labels:
                        labels.append(theLabel)
                        
            # Sorting labels (alphabetical order)
            labels = sorted(labels)

            # New histos
            newhistos = []
            for ind in range(len(histos)):

                nBins = histos[ind].GetNbinsX()

                # New histo
                newhisto = TH1F(histos[ind].GetName()+'_'+str(ind),\
                                histos[ind].GetTitle(),\
                                len(labels),
                                0,
                                len(labels))
                for ilabel in range(len(labels)):

                    # Set new labels
                    newhisto.GetXaxis().SetBinLabel(ilabel+1,labels[ilabel])
                    
                    fill=False 
                    for bin in range(1,nBins+1):
                        # Get the bin label
                        if histos[ind].GetXaxis().GetBinLabel(bin)==labels[ilabel]:
                            tmp = histos[ind].GetBinContent(bin)
                            newhisto.SetBinContent(ilabel+1,tmp)
                            fill=True
                            break
                    if not fill:
                        newhisto.SetBinContent(ilabel+1,0)
                        
                newhistos.append(newhisto)

            # overwrite histos
            histos=newhistos

        # Loop over datasets and histos
        for ind in range(0,len(histos)):

            if ref.observable.name in ['NPID','NAPID']:
                nBins = histos[ind].GetNbinsX()
                for b in range(1,nBins+1):
                    pid = int(histos[ind].GetXaxis().GetBinLabel(b))
                    if ref.observable.name=='NPID':
                        spid = self.main.multiparticles.GetName(pid)
                    else:
                        spid = self.main.multiparticles.GetAName(-pid,pid)
                    if spid!="":
                        histos[ind].GetXaxis().SetBinLabel(b,spid)

            # Scaling 
            histos[ind].Scale(scales[ind])
            
        # Stacking or superimposing histos ?
        stackmode = False
        if ref.stack==StackingMethodType.STACK or \
           ( ref.stack==StackingMethodType.AUTO and \
             self.main.stack==StackingMethodType.STACK ):
            stackmode=True

        # Setting AUTO settings
        if len(histos)==1:
            histos[0].SetLineColor(9)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
        elif len(histos)==2:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
        elif len(histos)==3:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
        elif len(histos)==4:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
        elif len(histos)==5:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
        elif len(histos)==6:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            histos[5].SetLineColor(2)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
                histos[5].SetFillColor(2)
                histos[5].SetFillStyle(3017)
        elif len(histos)==7:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            histos[5].SetLineColor(2)
            histos[6].SetLineColor(7)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
                histos[5].SetFillColor(2)
                histos[5].SetFillStyle(3017)
                histos[6].SetFillColor(7)
                histos[6].SetFillStyle(3022)
        elif len(histos)==8:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            histos[5].SetLineColor(2)
            histos[6].SetLineColor(7)
            histos[7].SetLineColor(3)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
                histos[5].SetFillColor(2)
                histos[5].SetFillStyle(3017)
                histos[6].SetFillColor(7)
                histos[6].SetFillStyle(3022)
                histos[7].SetFillColor(3)
                histos[7].SetFillStyle(3315)
        elif len(histos)==9:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            histos[5].SetLineColor(2)
            histos[6].SetLineColor(7)
            histos[7].SetLineColor(3)
            histos[8].SetLineColor(42)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
                histos[5].SetFillColor(2)
                histos[5].SetFillStyle(3017)
                histos[6].SetFillColor(7)
                histos[6].SetFillStyle(3022)
                histos[7].SetFillColor(3)
                histos[7].SetFillStyle(3315)
                histos[8].SetFillColor(42)
                histos[8].SetFillStyle(3351)
        elif len(histos)==10:
            histos[0].SetLineColor(9)
            histos[1].SetLineColor(46)
            histos[2].SetLineColor(8)
            histos[3].SetLineColor(4)
            histos[4].SetLineColor(6)
            histos[5].SetLineColor(2)
            histos[6].SetLineColor(7)
            histos[7].SetLineColor(3)
            histos[8].SetLineColor(42)
            histos[9].SetLineColor(48)
            if stackmode:
                histos[0].SetFillColor(9)
                histos[0].SetFillStyle(3004)
                histos[1].SetFillColor(46)
                histos[1].SetFillStyle(3005)
                histos[2].SetFillColor(8)
                histos[2].SetFillStyle(3006)
                histos[3].SetFillColor(4)
                histos[3].SetFillStyle(3007)
                histos[4].SetFillColor(6)
                histos[4].SetFillStyle(3013)
                histos[5].SetFillColor(2)
                histos[5].SetFillStyle(3017)
                histos[6].SetFillColor(7)
                histos[6].SetFillStyle(3022)
                histos[7].SetFillColor(3)
                histos[7].SetFillStyle(3315)
                histos[8].SetFillColor(42)
                histos[8].SetFillStyle(3351)
                histos[9].SetFillColor(48)
                histos[9].SetFillStyle(3481)
        else:
            histos[ind].SetLineColor(self.color)
            self.color += 1

        # Setting USER color
        for ind in range(0,len(histos)):

            # linecolor
            if self.main.datasets[ind].linecolor!=ColorType.AUTO:
                histos[ind].SetLineColor(ColorType.convert2root( \
                self.main.datasets[ind].linecolor,\
                self.main.datasets[ind].lineshade))

            # lineStyle
            histos[ind].SetLineStyle(LineStyleType.convert2code( \
                self.main.datasets[ind].linestyle))

            # linewidth
            histos[ind].SetLineWidth(self.main.datasets[ind].linewidth)

            # background color  
            if self.main.datasets[ind].backcolor!=ColorType.AUTO:
                histos[ind].SetFillColor(ColorType.convert2root( \
                self.main.datasets[ind].backcolor,\
                self.main.datasets[ind].backshade))

            # background color  
            if self.main.datasets[ind].backstyle!=BackStyleType.AUTO:
                histos[ind].SetFillStyle(BackStyleType.convert2code( \
                self.main.datasets[ind].backstyle))

        # Creating and filling the stack; computing the total number of events
        stack = THStack("mystack","")
        ntot = 0
        for item in histos:
            ntot+=item.GetEntries()
            stack.Add(item)

        # Drawing
        if stackmode:
            stack.Draw()
        else:
            stack.Draw("nostack")

        # Setting Y axis label
        axis_titleY = ref.GetYaxis()

        # Scale to one ?
        scale2one = False
        if ref.stack==StackingMethodType.NORMALIZE2ONE or \
           (self.main.stack==StackingMethodType.NORMALIZE2ONE and \
           ref.stack==StackingMethodType.AUTO):
            scale2one = True

        if scale2one:
            axis_titleY += " ( scaled to one )"
        elif self.main.normalize == NormalizeType.LUMI or \
           self.main.normalize == NormalizeType.LUMI_WEIGHT:
            axis_titleY += " ( L_{int} = " + str(self.main.lumi)+ " fb^{-1} )"
        elif self.main.normalize == NormalizeType.NONE:
            axis_titleY += " (not normalized)"

        if ref.titleY!="": 
            axis_titleY = self.NiceTitle(ref.titleY)

        stack.GetYaxis().SetTitle(axis_titleY)
        if(len(axis_titleY) > 35): 
           stack.GetYaxis().SetTitleSize(0.04)
        else:
           stack.GetYaxis().SetTitleSize(0.06)
        stack.GetYaxis().SetTitleFont(22)
        stack.GetYaxis().SetLabelSize(0.04)

        # Setting X axis label
        if ref.titleX=="": 
            axis_titleX = ref.GetXaxis()
        else:
            axis_titleX = self.NiceTitle(ref.titleX)
        
        # Setting X axis label
        stack.GetXaxis().SetTitle(axis_titleX)
        stack.GetXaxis().SetTitleSize(0.06)
        stack.GetXaxis().SetTitleFont(22)
        stack.GetXaxis().SetLabelSize(0.04)

        # Setting Log scale
        if ref.logX and ntot != 0:
            canvas.SetLogx()
        if ref.logY and ntot != 0:
            canvas.SetLogy()
        
        # Displaying a legend
        if len(self.main.datasets)>1:
            ymin_legend = .9-.055*len(histos)
            if ymin_legend<0.1:
                ymin_legend = 0.1
            legend = TLegend(.65,ymin_legend,.9,.9)
            legend.SetTextSize(0.05); 
            legend.SetTextFont(22); 
            for ind in range(0,len(histos)):
                legend.AddEntry(histos[ind],self.NiceTitle(self.main.datasets[ind].title))
            legend.SetFillColor(0)    
            legend.Draw()

        if not preview:
            
            # Save the canvas in the report format
            canvas.SaveAs(self.output_path+"/selection_"+str(irelhisto)+"."+\
                               ReportFormatType.convert2filetype(self.mode))

            # Save the canvas in the C format
            canvas.SaveAs(self.output_path+"/selection_"+str(irelhisto)+".C")
            
        else:
            # break 
            answer=raw_input("Press enter to continue : ")
  
            
                   
