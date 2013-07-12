////////////////////////////////////////////////////////////////////////////////
//  
//  Copyright (C) 2012 Eric Conte, Benjamin Fuks, Guillaume Serret
//  The MadAnalysis development team, email: <ma5team@iphc.cnrs.fr>
//  
//  This file is part of MadAnalysis 5.
//  Official website: <http://madanalysis.irmp.ucl.ac.be>
//  
//  MadAnalysis 5 is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  MadAnalysis 5 is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with MadAnalysis 5. If not, see <http://www.gnu.org/licenses/>
//  
////////////////////////////////////////////////////////////////////////////////


// STL headers
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

// SampleHeader headers
#include "Core/SampleAnalyzer.h"
#include "Core/Configuration.h"
#include "Services/logger.h"
#include "Services/logger.h"
#include "Core/Counter.h"

// -----------------------------------------------------------------------
// printSyntax
// -----------------------------------------------------------------------
void printSyntax()
{
  INFO << std::endl;
  INFO << "Syntax : LHEanalyzer [option] <filelist>" << std::endl;
  INFO << " with <filelist> = txt file containing all LHE file names"
       << std::endl;
  INFO << " with [option] = --checkevent, "
       << "--analysis=\"name_of_the_analysis\"" << std::endl;
  INFO << std::endl;
}



// -----------------------------------------------------------------------
// execute function
// -----------------------------------------------------------------------
int execute(int argc, char *argv[])
{
  // -----------------------------------------------------------------------
  //                                HEADER
  // -----------------------------------------------------------------------
  INFO << "    * SampleAnalyzer 1.5 for MadAnalysis 5 - Welcome.";
  INFO << std::endl;

  // -----------------------------------------------------------------------
  //                      CHECKING NUMBER OF ARGUMENTS
  // -----------------------------------------------------------------------
  if (argc<2)
  {
    ERROR << "number of arguments is incorrect !!!" << std::endl;
    printSyntax();
    return 1;
  }

  // -----------------------------------------------------------------------
  //                        DECODING ARGUMENTS
  // -----------------------------------------------------------------------
  Configuration cfg;
  std::string filename = argv[argc-1];
  if (argc>=3) for (unsigned int i=1;i<(argc-1);i++)
  {
    // converting const char to string
    std::string option = std::string(argv[i]);

    // safety : skip empty string
    if (option.size()==0) continue;

    // check event
    if (option=="--checkevent") cfg.check_event = true;

    // analysis selection
    else if (option.find("--analysis=")==0)
    {
      cfg.analysis_selected = true;
      cfg.analysis_name = option.substr(11);
    }

    // unknown option
    else
    {
      ERROR << "argument '" << option << "' is unknown !!!" << std::endl;
      printSyntax();
      return 1;
    }
  }

  // -----------------------------------------------------------------------
  //                          DUMP ARGUMENTS
  // -----------------------------------------------------------------------
  INFO << "    * Option choices: ";
  if (!cfg.check_event && !cfg.analysis_selected) INFO << "none";
  if (cfg.check_event) INFO << "checking event mode ";
  if (cfg.check_event && cfg.analysis_selected) INFO << "+ ";
  if (cfg.analysis_selected) INFO << "selecting analysis '" << cfg.analysis_name << "'";
  INFO << "." << std::endl; 


  // -----------------------------------------------------------------------
  //                    EXTRACTING LIST OF SAMPLES
  // -----------------------------------------------------------------------
  INFO << "    * Extracting the following sample files:" << std::endl;

  std::ifstream input(filename.c_str());
  if (!input)
    {
      ERROR << "LHE files list called '"<< filename 
            << "' is not found !!!" << std::endl;
      return 1;
    }

  std::vector<std::string> filenames;
  std::string tmp;
  while(!input.eof() && !input.fail())
  {
    getline(input,tmp);
    std::stringstream str; str<<tmp; tmp=""; str>>tmp;
    if (!tmp.empty()) filenames.push_back(tmp);
  }

  input.close();

  if (filenames.size()==0)
  {
      ERROR << "LHE files list called '"<< filename 
            << "' contains no LHE file !!!" << std::endl;
      return 1;
  }
  
  // -----------------------------------------------------------------------
  //                      EXTRACTING NAME FOR ANALYSIS
  // -----------------------------------------------------------------------
  std::string::size_type pos=filename.rfind('.');
  if (pos!=std::string::npos) filename.resize(pos);
  pos=filename.rfind('/');
  if (pos!=std::string::npos) filename.erase(0,pos+1);


  // -----------------------------------------------------------------------
  //                          LAUNCHING ANALYZER       
  // -----------------------------------------------------------------------
  SampleAnalyzer analyzer(filename, cfg, filenames);
  analyzer.Run();

  // -----------------------------------------------------------------------
  //                               END                
  // -----------------------------------------------------------------------
  INFO << "    * Goodbye."<<std::endl;
  return 0;
}



// -----------------------------------------------------------------------
// main program
// -----------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Execute main program
  int result=execute(argc,argv);

  // Kill singletons
  PhysicsService::kill();
  Logger::kill();

  // Return error code (=0 if no error)
  return result;
}
