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


#include "Core/AnalysisManager.h"
#include "Services/logger.h"
#include <stdlib.h>


// -----------------------------------------------------------------------------
// Get
// -----------------------------------------------------------------------------
AnalysisBase* AnalysisManager::Get(const std::string& name)
{
  for (unsigned int i=0;i<stack_.size();i++)
    if (name==stack_[i]->name())
      return stack_[i];
  return 0;
} 


// -----------------------------------------------------------------------------
// ChoiceAnalysis
// -----------------------------------------------------------------------------
AnalysisBase* AnalysisManager::ChoiceAnalysis()
{
  // Display the list of analyses
  INFO << std::endl;
  Print();
  INFO << std::endl;

  // Choose an analysis
  INFO << "Please, choose an analysis (0.."
       << stack_.size()-1 << ") : ";
  unsigned int n=0;
  std::cin >> n;

  // Check the choice
  if (n<0 || n>=stack_.size())
    {
      ERROR << "wrong analysis" << std::endl;
      exit(1);
    }

  // Return the analysis
  return stack_[n];
}


// -----------------------------------------------------------------------------
// Print
// -----------------------------------------------------------------------------
void AnalysisManager::Print()
{
  INFO << "List of the implemented analyses :" << std::endl;
  INFO << "----------------------------------" << std::endl;

  //Display all analyses
  for (unsigned int i=0;i<stack_.size();i++)
    INFO << i << "\t" << stack_[i]->name() << std::endl;
}
