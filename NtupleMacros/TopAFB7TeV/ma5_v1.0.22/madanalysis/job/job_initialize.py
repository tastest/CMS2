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


from madanalysis.enumeration.ma5_running_type import MA5RunningType
import logging

def WriteHadronicList(file,main):
    file.write('  // definition of the multiparticle "hadronic"\n')
    for item in main.multiparticles.Get("hadronic"):
        file.write('  PHYSICS->mcConfig().AddHadronicId('+str(item)+');\n')


def WriteInvisibleList(file,main):
    file.write('  // definition of the multiparticle "invisible"\n')
    for item in main.multiparticles.Get("invisible"):
        file.write('  PHYSICS->mcConfig().AddInvisibleId('+str(item)+');\n')


def WriteJobInitialize(file,main):

    # Function header
    file.write('void user::Initialize()\n{\n')

    # mcConfig initialization
    if main.mode!=MA5RunningType.RECO:
        file.write('  // Initializing PhysicsService for MC\n') 
        file.write('  PHYSICS->mcConfig().Reset();\n\n')
        WriteHadronicList(file,main)
        file.write('\n')
        WriteInvisibleList(file,main)
        file.write('\n')

    # recConfig initialization
    if main.mode==MA5RunningType.RECO:
        file.write('  // Initializing PhysicsService for RECO\n') 
        file.write('  PHYSICS->recConfig().Reset();\n\n')
        if main.isolation.algo=='DELTAR':
            file.write('  PHYSICS->recConfig().UseDeltaRIsolation('+\
                       str(main.isolation.deltaR) + ');\n')
        else:
            file.write('  PHYSICS->recConfig().UseSumPTIsolation('+\
                       str(main.isolation.sumPT) + ',' +\
                       str(main.isolation.ET_PT) + ');\n')
        file.write('\n')

    # Counting number of plots and cuts
    Nhistos = 0
    Ncuts   = 0
    for item in main.selection.table:
        if item.__class__.__name__=="Histogram":
            Nhistos+=1
        elif item.__class__.__name__=="Cut":
            Ncuts+=1

    # Declaring array of plots
    if Nhistos!=0:
        file.write('  // Initializing histogram array\n')
        file.write('  plots_array_=new TClonesArray("TH1F",' + \
                   str(Nhistos)+');\n\n')
        file.write('  plots_nevents_.resize('+str(Nhistos)+',0);\n\n')

    # Initializing array of cuts     
    if Ncuts!=0:
        file.write('  // Initializing cut array\n')
        file.write('  cuts_array_.resize('+str(Ncuts)+',0);\n\n')

    # Initializing each item
    ihisto  = 0
    icut    = 0
    iglobal = 0
    file.write('  // Initializing each selection item\n')
    for item in main.selection.table:

        # Histogram case
        if item.__class__.__name__=="Histogram":
            file.write('  new ((*plots_array_)['+str(ihisto)+']) ')
            file.write('TH1F("selection_'+str(ihisto)+'","",')       
            if item.observable in ['NPID','NAPID']:
                file.write('1,0,1);\n')
            else:
                file.write(str(item.nbins)+','+\
                           str(item.xmin)+','+\
                           str(item.xmax)+');\n')
            file.write('  P'+str(iglobal)+'_ = dynamic_cast<TH1F*>((*plots_array_)[' +\
                     str(ihisto)+']);\n')
            file.write('  NP'+str(iglobal)+'_=&(plots_nevents_['+str(ihisto)+']);\n')
            ihisto+=1

        # Cut case
        elif item.__class__.__name__=="Cut":
            file.write('  N'+str(iglobal)+'_=&(cuts_array_['+str(icut)+']);\n')
            icut+=1

        iglobal+=1

    # End    
    file.write('}\n\n')
