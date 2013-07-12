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


#ifndef LHE_READER_h
#define LHE_READER_h

// SampleAnalyzer headers
#include "Core/ReaderTextBase.h"

class LHEReader : public ReaderTextBase
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 protected:

  bool firstevent_;

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  //! Constructor without argument
  LHEReader()
  { firstevent_=false; }

	//! Destructor
  virtual ~LHEReader()
  { }

  //! Read the header
  virtual bool ReadHeader(SampleFormat& mySample);

  //! Finalize the header
  virtual bool FinalizeHeader(SampleFormat& mySample);

  //! Read the event
  virtual bool ReadEvent(EventFormat& myEvent, SampleFormat& mySample);

  //! Finalize the event
  virtual bool FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent);


 private:

  //! Fill the header from text line 
  void FillHeaderProcessLine(const std::string& line, SampleFormat& mySample);
  void FillHeaderInitLine   (const std::string& line, SampleFormat& mySample);

  //! Fill the event from text line 
  void FillEventInitLine(const std::string& line, EventFormat& myFormat);
  void FillEventParticleLine(const std::string& line, EventFormat& myFormat);

};

#endif
