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
def WriteJobFinalize(file,main):
 
    # Function header
    file.write('void user::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)\n{\n')

    # Opening ROOT file
    file.write('  // Creating output root file\n')
    file.write('  TFile* output = new TFile((outputName_+".root").c_str(),"RECREATE");\n')

    # Creating subdirectories
    file.write('  output->mkdir("general");\n')
    file.write('  output->mkdir("plots");\n')
    file.write('  output->mkdir("cuts");\n\n')

    # Writing global information
    file.write('  // Saving global information\n')
    WriteGlobalInfo(file)
    file.write('\n')

    # Writing specific information
    WriteSpecificInfo(file,main)
    file.write('\n')

    # End
    file.write('  // Closing the output file\n')
    file.write('  delete output;\n')
    file.write('}\n')


def WriteGlobalInfo(file):
    file.write('  {\n')
    file.write('    TClonesArray filename_array("TObjString",files.size()+1);\n')
    file.write('    TVector xsection_array(files.size()+1);\n')
    file.write('    TVector xerror_array(files.size()+1);\n')
    file.write('    TVector nevent_array(files.size()+1);\n')
    file.write('    for (unsigned int i=0;i<files.size();i++)\n')
    file.write('    {\n')
    file.write('      new(filename_array[i]) TObjString(files[i].name().c_str());\n')
    file.write('      if (files[i].mc()!=0) {\n')
    file.write('        xsection_array[i] = files[i].mc()->xsection();\n')
    file.write('        xerror_array[i] = files[i].mc()->xsection_error();\n')
    file.write('        }\n')
    file.write('      else {\n')
    file.write('        xsection_array[i] = 0.;\n')
    file.write('        xerror_array[i] = 0.;\n')
    file.write('        }\n')
    file.write('      nevent_array[i] = files[i].nevents();\n')
    file.write('    }\n')
    file.write('    {\n')
    file.write('      new(filename_array[files.size()]) TObjString(summary.name().c_str());\n')
    file.write('      if (summary.mc()!=0) {\n')
    file.write('        xsection_array[files.size()] = summary.mc()->xsection();\n')
    file.write('        xerror_array[files.size()] = summary.mc()->xsection_error();\n')
    file.write('        }\n')
    file.write('      else {\n')
    file.write('        xsection_array[files.size()] = 0.;\n')
    file.write('        xerror_array[files.size()] = 0.;\n')
    file.write('        }\n')
    file.write('      nevent_array[files.size()] = summary.nevents();\n')
    file.write('    }\n')
    file.write('    output->cd("general");\n')
    file.write('    filename_array.Write("filenames",TObject::kSingleKey);\n')
    file.write('    xsection_array.Write("xsections");\n')
    file.write('    xerror_array.Write("xerrors");\n')
    file.write('    nevent_array.Write("nevents");\n')
    file.write('    output->cd("");\n')
    file.write('  }\n')


def WriteSpecificInfo(file,main):
    
    # Counting number of plots and cuts
    Nhistos = 0
    Ncuts   = 0
    PIDplot = False
    
    for item in main.selection.table:
        if item.__class__.__name__=="Histogram":
            Nhistos+=1
            if item.observable.name in ['NPID','NAPID']:
                PIDplot = True
        elif item.__class__.__name__=="Cut":
            Ncuts+=1

    # Layouting NPID/NAPID plots before saving
    if PIDplot:

        file.write('  // Finalizing special histograms\n')
        for ind in range(0,len(main.selection.table)):
            if main.selection.table[ind].__class__.__name__=="Histogram" and \
               main.selection.table[ind].observable.name in ['NPID','NAPID']:
                file.write('  C'+str(ind)+'_.FillHisto(P' +\
                   str(ind)+'_);\n')
        file.write('\n')      

    # Saving plots in a ROOT file
    if Nhistos!=0:
        file.write('  // Saving histogram\n')
        file.write('  {\n')
        file.write('    output->cd("plots");\n')
        file.write('    plots_array_->Write("plots_array",' +\
                   'TObject::kSingleKey);\n')
        file.write('    TVector plots2save_nevents(plots_nevents_.size());\n')
        file.write('    for (unsigned int i=0;i<plots_nevents_.size();i++)\n')
        file.write('      plots2save_nevents[i]=plots_nevents_[i];\n')
        file.write('    plots2save_nevents.Write("nevents");\n')
        file.write('    output->cd("");\n')
        file.write('  }\n')

    # Saving cuts in a ROOT file
    if Ncuts!=0:
        file.write('\n  // Saving efficiencies\n')
        file.write('  {\n')
        file.write('    TVector cuts2save_array(cuts_array_.size());\n')
        file.write('    for (unsigned int i=0;i<cuts_array_.size();i++)\n')
        file.write('      cuts2save_array[i]=cuts_array_[i];\n')
        file.write('    output->cd("cuts");\n')
        file.write('    cuts2save_array.Write("cuts");\n')
        file.write('    output->cd("");\n')
        file.write('  }\n')
