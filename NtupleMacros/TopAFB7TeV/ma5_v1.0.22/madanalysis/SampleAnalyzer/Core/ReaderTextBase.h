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


#ifndef READER_TEXT_BASE_h
#define READER_TEXT_BASE_h

// STL headers
#include <fstream>
#include <iostream>
#include <sstream>

// RFIO headers
#ifdef RFIO_USE
   #include "Core/castorstream.h"
#endif

// SampleAnalyzer headers
#include "Core/ReaderBase.h"


class ReaderTextBase : public ReaderBase
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 protected:

  /// Streaming for reading input
  std::istream*  input_;


  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  /// Constructor without argument
  ReaderTextBase()
  {
    input_=0;
  }

	/// Destructor
  virtual ~ReaderTextBase()
  {
    if (input_ !=0) delete input_;
  }

  /// Initialize
  virtual bool Initialize(const std::string& rawfilename,
                          const Configuration& cfg);

  /// Read the sample (virtual pure)
  virtual bool ReadHeader(SampleFormat& mySample) = 0;

  /// Finalize the header (virtual pure)
  virtual bool FinalizeHeader(SampleFormat& mySample) = 0;

  /// Read the event (virtual pure)
  virtual bool ReadEvent(EventFormat& myEvent, SampleFormat& mySample) = 0;

  /// Finalize the event (virtual pure)
  virtual bool FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent) = 0;

  /// Finalize
  virtual bool Finalize();

  /// Read line text
  bool ReadLine(std::string& line);

};

#endif
