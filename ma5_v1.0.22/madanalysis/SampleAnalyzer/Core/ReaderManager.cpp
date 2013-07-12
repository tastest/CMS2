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


#include "Core/ReaderManager.h"
// STL headers
#include <typeinfo>
#include <functional>
#include <algorithm>

// SampleAnalyze headers
#include "Reader/LHEReader.h"
#include "Reader/LHCOReader.h"
#include "Reader/STDHEPreader.h"
#include "Reader/HEPMCReader.h"
#ifdef FAC_USE
#include "Reader/ROOTReader.h"
#endif
#include "Services/logger.h"


// -----------------------------------------------------------------------------
// BuildTable
// -----------------------------------------------------------------------------
void ReaderManager::BuildTable()
{
  // Adding LHE reader
  LHEReader* lhe = new LHEReader();
  Add("lhe",lhe);
#ifdef ZIP_USE
  Add("lhe.gz",lhe);
#endif

  LHCOReader* lhco = new LHCOReader();
  Add("lhco",lhco);
#ifdef ZIP_USE
  Add("lhco.gz",lhco);
#endif

  STDHEPreader* stdhep = new STDHEPreader();
  Add("hep",stdhep);
#ifdef ZIP_USE
  Add("hep.gz",stdhep);
#endif

  HEPMCReader* hepmc = new HEPMCReader();
  Add("hepmc",hepmc);
#ifdef ZIP_USE
  Add("hepmc.gz",hepmc);
#endif

#ifdef FAC_USE
  ROOTReader* root = new ROOTReader();
  Add("root",root);
#ifdef ZIP_USE
  Add("root.gz",root);
#endif
#endif

}


// -----------------------------------------------------------------------------
// Get
// -----------------------------------------------------------------------------
ReaderBase* ReaderManager::Get(std::string filename)
{
  // Set the extension in lower case
  std::transform(filename.begin(), filename.end(),
                 filename.begin(), std::ptr_fun<int, int>(std::tolower));

 // Loop over names
  for (std::map<std::string, unsigned int>::const_iterator
         it = names_.begin(); it != names_.end(); it++)
  {
    // easy case to reject
    if (filename.size()<it->first.size()) continue;

    // find pattern in filename
    if (it->first.compare(0,
                          it->first.size(),
                          filename,
                          filename.size()-it->first.size(),
                          it->first.size())==0)
    {
      return readers_[it->second];
    }
  }

  // No reader found : return null pointer
  return 0;
}


// -----------------------------------------------------------------------------
// Add
// -----------------------------------------------------------------------------
bool ReaderManager::Add(std::string extension, ReaderBase* reader)
{
  // Set the extension in lower case
  std::transform(extension.begin(), extension.end(),
                extension.begin(), std::ptr_fun<int, int>(std::tolower));

  // Insert extension in the data base
  std::pair<std::map<std::string,unsigned int>::iterator,bool>
    found = names_.insert(std::make_pair(extension,0));
  
  // Check if name insertion is failed
  if (!found.second) return false;

  // Look for reader in the data base
  for (unsigned int i=0;i<readers_.size();i++)
  {
    if (readers_[i]==reader)
    {
      found.first->second=i;
      return true;
    } 
  }

  // Case where the reader is not found in the data base
  readers_.push_back(reader);
  found.first->second=readers_.size()-1;
  return true;  
}


// -----------------------------------------------------------------------------
// Print
// -----------------------------------------------------------------------------
void ReaderManager::Print() const
{
  // Header
  INFO << "------------------------------------------" << std::endl;

  // Loop over names
  for (std::map<std::string, unsigned int>::const_iterator
         it = names_.begin(); it != names_.end(); it++)
  {
    INFO.width(10); 
    INFO << it->first;
    INFO << " : ";
    INFO << typeid(readers_[it->second]).name();
    INFO << " @ " ;
    INFO << readers_[it->second];
    INFO << std::endl; 
  }

  // Foot
  INFO << "------------------------------------------" << std::endl;
}


