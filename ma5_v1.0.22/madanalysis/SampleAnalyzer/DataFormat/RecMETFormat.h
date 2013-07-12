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


#ifndef RecMETFormat_h
#define RecMETFormat_h

// STL headers
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

// RecParticleFormat
#include "DataFormat/RecParticleFormat.h"
#include "Services/logger.h"

class LHCOReader;
class ROOTReader;

class RecMETFormat
{

  friend class LHCOReader;
  friend class ROOTReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------

 protected:
  TVector2 met_;     /// transverse missing energy

  // -------------------------------------------------------------
  //                        method members
  // -------------------------------------------------------------

 public:

  /// Constructor without arguments
  RecMETFormat()
  { Reset(); }

  /// Destructor                                                      
  ~RecMETFormat()
  {}

  /// Display information
  void Print() const
  {
    INFO << "MET = (" << met_.X() << " , " << met_.Y() << " )" ;
  }

  /// Clear all information
  void Reset()
  {
    met_.Set(0.,0.); 
  }

  /// Accessor to the MET vector
  const TVector2& vector()  const {return met_;}

  /// Accessor to the MET magnitude
  const Float_t magnitude() const {return met_.Mod(); }

  /// Accessor to the x-component MET
  const Float_t x()      const {return met_.X();}

  /// Accessor to the y-component MET
  const Float_t y()      const {return met_.Y();}

};

#endif
