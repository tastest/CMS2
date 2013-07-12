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


#ifndef PROCESS_FORMAT_H
#define PROCESS_FORMAT_H

// STL headers
#include <map>

// SampleAnalyzer headers
#include "Services/logger.h"

class LHEReader;
class LHCOReader;
class HEPMCReader;
class ROOTReader;

class ProcessFormat
{
  friend class LHEReader;
  friend class LHCOReader;
  friend class HEPMCReader;
  friend class ROOTReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:

  Float_t xsectionMean_;    /// cross-section (pb)
  Float_t xsectionError_;   /// statistical error on the cross-section
  Float_t weightMax_;       /// maximum weight encountered in the events
  UInt_t  processId_;       /// number identifying the process


  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public:

  /// Constructor without arguments
  ProcessFormat() {Reset();}

  /// Destructor
  ~ProcessFormat() {}

  /// Accessor to the mean value of the process cross section
  const Float_t&  xsection()      const {return xsectionMean_;  }

  /// Accessor to the mean value of the process cross section
  const Float_t&  xsectionMean()  const {return xsectionMean_;  }

  /// Accessor to the error value of the process cross section
  const Float_t&  xsectionError() const {return xsectionError_; }

  /// Accessor to the highest weight
  const Float_t&  weightMax()     const {return weightMax_;     }

  /// Accessor to the process identity
  const UInt_t&   processId()     const {return processId_;     }

  /// Clearing all information
  void Reset()
  {
     xsectionMean_=0.; xsectionError_=0.;
     weightMax_=0.;    processId_=0;
  }

  /// Displaying data member values
  void Print() const
  {
    INFO << "processId="        << processId_
         << " - xsectionMean="  << xsectionMean_
         << " - xsectionError=" << xsectionError_
         << " - weightMax="     << weightMax_ << std::endl;
  }

};


#endif
