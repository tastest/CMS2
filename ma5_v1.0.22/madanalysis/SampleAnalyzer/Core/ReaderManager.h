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


#ifndef READER_MANAGER_h
#define READER_MANAGER_h

// STL headers
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

// SampleAnalyzer headers
#include "Core/ReaderBase.h"


class ReaderManager
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 protected:

  //! Mapping between names and readers
  std::map<std::string, unsigned int> names_;

  //! List of readers
  std::vector<ReaderBase*> readers_;


  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  ReaderManager()
  { }

	//! Destructor
  ~ReaderManager()
  { 
    for (unsigned int i=0;i<readers_.size();i++)
      if (readers_[i]!=0) delete readers_[i];
  }

  //! Find a reader for a given filename
  ReaderBase* Get(std::string filename);

  //! Add a reader to the list
  bool Add(std::string extension, ReaderBase* reader);

  //! Display the content of the ReaderManager
  void Print() const;

  //! Build the table
  void BuildTable();

};

#endif
