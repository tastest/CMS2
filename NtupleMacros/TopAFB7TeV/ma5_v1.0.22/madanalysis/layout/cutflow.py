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


from madanalysis.enumeration.uncertainty_type import UncertaintyType
from madanalysis.enumeration.normalize_type import NormalizeType
from math import *

class Measure:

    def __init__(self):
        self.mean=0.
        self.error=0.


class CutFlow:

    def __init__(self,dataset,selection,lumi,main):
        self.dataset   = dataset
        self.selection = selection
        self.lumi      = lumi
        self.main      = main 

        self.Ntotal=0
        self.Nselected=[]
        self.Nrejected=[]
        self.eff=[]
        self.effcumu=[]
        self.formulaSBratio()


    def formulaSBratio(self):
        from ROOT import TFormula
        text = self.main.SBratio
        text = text.replace("S","x")
        text = text.replace("B","y")
        self.Mformula = TFormula("SBratio",text)
        text = self.main.SBerror
        text = text.replace("ES","z")
        text = text.replace("EB","t")
        text = text.replace("S","x")
        text = text.replace("B","y")
        self.Eformula = TFormula("SBerror",text)
    

    def calculateBSratio(self,B,eB,S,eS):
        value = Measure()
        value.mean  = self.Mformula.Eval(S,B)
        value.error = self.Eformula.Eval(S,B,eS,eB) 
        return value
                           
    @staticmethod
    def binomialNEventError(k,N):
        if N==0:
            return 0.
        else:
            return sqrt( float(k*abs(N-k)) / float(N) )

    @staticmethod
    def binomialError(k,N):
        if N==0:
            return 0.
        else:
            return sqrt( float(k*abs(N-k)) / float(N*N*N) )

    def initializeFromCutflow(self,cutflows):
        # Ntotal
        self.Ntotal=Measure()
        for item in cutflows:
            self.Ntotal.mean  += item.Ntotal.mean
            self.Ntotal.error += item.Ntotal.error**2
        self.Ntotal.error = sqrt(self.Ntotal.error)

        # Prepare vectors
        for i in range(0,len(cutflows[0].Nselected)):
            self.Nselected.append(Measure())
            self.Nrejected.append(Measure())
            self.eff.append(Measure())
            self.effcumu.append(Measure())

        # Fill selected and rejected
        for iset in range (0,len(cutflows)):
            for icut in range (0,len(cutflows[iset].Nselected)):
                self.Nselected[icut].mean  += cutflows[iset].Nselected[icut].mean
                self.Nrejected[icut].mean  += cutflows[iset].Nrejected[icut].mean
                self.Nselected[icut].error += cutflows[iset].Nselected[icut].error**2
                self.Nrejected[icut].error += cutflows[iset].Nrejected[icut].error**2
        self.Nselected[icut].error = sqrt(self.Nselected[icut].error)
        self.Nrejected[icut].error = sqrt(self.Nrejected[icut].error)

        # Compute efficiencies
        for i in range(0,len(self.eff)):
            if self.Ntotal.mean!=0:
                self.effcumu[i].mean=float(self.Nselected[i].mean)/float(self.Ntotal.mean)
                self.effcumu[i].error=CutFlow.binomialError(self.Nselected[i].mean,self.Ntotal.mean)
            if i==0:
                if self.Ntotal.mean!=0:
                    self.eff[i].mean=float(self.Nselected[i].mean)/float(self.Ntotal.mean)
                    self.eff[i].error=CutFlow.binomialError(self.Nselected[i].mean,self.Ntotal.mean)
            else:
                if self.Nselected[i-1].mean!=0:
                    self.eff[i].mean=float(self.Nselected[i].mean)/float(self.Nselected[i-1].mean)
                    self.eff[i].error=CutFlow.binomialError(self.Nselected[i].mean,self.Nselected[i-1].mean)
        
    
    def initializeFromFile(self,file):

        # Getting cuts from ROOT file    
        cut = file.Get("cuts/cuts","TVectorT<float>")
        if cut is None:
            return False

        # Preparing architecture for vectors
        self.Ntotal=Measure()
        for iabscut in range(0,len(self.selection)):
            if self.selection[iabscut].__class__.__name__!="Cut":
                continue
            self.Nselected.append(Measure())
            self.Nrejected.append(Measure())
            self.eff.append(Measure())
            self.effcumu.append(Measure())

        # Extracting Nselected information
        for icut in range(0,len(self.Nselected)):
            self.Nselected[icut].mean=cut[icut]

        # Extracting Ntotal
        self.Ntotal.mean=self.dataset.measured_n

        return True


    def calculate(self):
        
        # Calculating Nrejected
        for icut in range(0,len(self.Nselected)):
            if icut==0:
                self.Nrejected[icut].mean = self.Ntotal.mean - self.Nselected[icut].mean
            else:
                self.Nrejected[icut].mean = self.Nselected[icut-1].mean - self.Nselected[icut].mean

        # Calculating errors on Naccepted and Nrejected
        for icut in range(0,len(self.Nselected)):
            if icut==0:
                self.Nselected[icut].error = CutFlow.binomialNEventError(self.Nselected[icut].mean,self.Ntotal.mean)
                self.Nrejected[icut].error = CutFlow.binomialNEventError(self.Nrejected[icut].mean,self.Ntotal.mean)
            else:
                self.Nselected[icut].error = CutFlow.binomialNEventError(self.Nselected[icut].mean,self.Ntotal.mean)
                self.Nrejected[icut].error = CutFlow.binomialNEventError(self.Nrejected[icut].mean,self.Ntotal.mean)

        # efficiency calculation and its error
        for icut in range(0,len(self.Nselected)):
            
            if icut==0:
                if self.Ntotal.mean==0:
                    self.eff[icut].mean = 0
                else:                    
                    self.eff[icut].mean = float(self.Nselected[icut].mean) / \
                                                float(self.Ntotal.mean)
                self.eff[icut].error = CutFlow.binomialError(self.Nselected[icut].mean,self.Ntotal.mean)
            else:
                if self.Nselected[icut-1].mean==0:
                    self.eff[icut].mean = 0
                else:                    
                    self.eff[icut].mean = float(self.Nselected[icut].mean) / \
                                                float(self.Nselected[icut-1].mean)
                self.eff[icut].error = CutFlow.binomialError(self.Nselected[icut].mean,self.Nselected[icut-1].mean)
 
            if self.Ntotal.mean==0:
                self.effcumu[icut].mean=0
            else:
                self.effcumu[icut].mean = float(self.Nselected[icut].mean) / \
                                                float(self.Ntotal.mean)
            self.effcumu[icut].error = CutFlow.binomialError(self.Nselected[icut].mean,self.Ntotal.mean)

        # Getting xsection
        xsection = self.dataset.measured_xsection
        xerror  = self.dataset.measured_xerror
        if self.dataset.xsection!=0.:
            xsection=self.dataset.xsection
            xerror=0

        # Saving ntotal
        ntot = 0.+self.Ntotal.mean

        # Scaling Ntotal
        if self.main.normalize == NormalizeType.LUMI:
            self.Ntotal.error = xerror * self.lumi * 1000
            self.Ntotal.mean = xsection * self.lumi * 1000
        elif self.main.normalize == NormalizeType.LUMI_WEIGHT:
            self.Ntotal.error = xerror * self.lumi * 1000 * \
                                      self.dataset.weight
            self.Ntotal.mean = xsection * self.lumi * 1000 * \
                                      self.dataset.weight

        # Scaling Nselected
        for icut in range(0,len(self.Nselected)):
            if self.main.normalize == NormalizeType.LUMI:

                # error due to xsec 
                errXsec = xerror * self.lumi * 1000 * self.Nselected[icut].mean / ntot
                
                # scale factor
                factor = xsection * self.lumi * 1000 / ntot
                self.Nselected[icut].mean  *= factor
                self.Nselected[icut].error = CutFlow.binomialNEventError(self.Nselected[icut].mean,self.Ntotal.mean)

                # compute final error
                self.Nselected[icut].error = sqrt(\
                    self.Nselected[icut].error**2 + errXsec**2)
                
            elif self.main.normalize == NormalizeType.LUMI_WEIGHT:

                # error due to xsec 
                errXsec = xerror * self.lumi * 1000 * self.dataset.weight * self.Nselected[icut].mean / ntot
                
                # scale factor
                factor = xsection * self.lumi * self.dataset.weight * 1000 / ntot
                self.Nselected[icut].mean  *= factor
                self.Nselected[icut].error = CutFlow.binomialNEventError(self.Nselected[icut].mean,self.Ntotal.mean)
                
                # compute final error
                self.Nselected[icut].error = sqrt(\
                    self.Nselected[icut].error**2 + errXsec**2)

        # Scaling Nrejected
        for icut in range(0,len(self.Nrejected)):

            if self.main.normalize == NormalizeType.LUMI:

                # error due to xsec 
                errXsec = xerror * self.lumi * 1000 * self.Nrejected[icut].mean / ntot
                
                # scale factor
                factor = xsection * self.lumi * 1000 / ntot
                self.Nrejected[icut].mean  *= factor
                self.Nrejected[icut].error = CutFlow.binomialNEventError(self.Nrejected[icut].mean,self.Ntotal.mean)

                # compute final error
                self.Nrejected[icut].error = sqrt(\
                    self.Nrejected[icut].error**2 + errXsec**2)
                
            elif self.main.normalize == NormalizeType.LUMI_WEIGHT:

                # error due to xsec 
                errXsec = xerror * self.lumi * 1000 * self.dataset.weight * self.Nrejected[icut].mean / ntot
                
                # scale factor
                factor = xsection * self.lumi * self.dataset.weight * 1000 / ntot
                self.Nrejected[icut].mean  *= factor
                self.Nrejected[icut].error = CutFlow.binomialNEventError(self.Nrejected[icut].mean,self.Ntotal.mean)

                # compute final error
                self.Nrejected[icut].error = sqrt(\
                    self.Nrejected[icut].error**2 + errXsec**2)

                                            
                    
            
                                            

        
        
            
                
        
        

        
