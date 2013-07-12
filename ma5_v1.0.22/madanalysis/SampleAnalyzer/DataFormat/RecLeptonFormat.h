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


#ifndef RecLeptonFormat_h
#define RecLeptonFormat_h

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

class RecLeptonFormat : public RecParticleFormat
{

  friend class LHCOReader;
  friend class ROOTReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------             
 protected:

  Bool_t charge_;       /// charge of the particle 0 = -1, 1 = +1
  Float_t sumET_isol_;  /// sumET in an isolation cone
  Float_t sumPT_isol_;  /// sumPT in an isolation cone

  // -------------------------------------------------------------
  //                        method members
  // -------------------------------------------------------------             
 public:

  /// Constructor without arguments
  RecLeptonFormat()
  { Reset(); }

  /// Destructor
  virtual ~RecLeptonFormat()
  {}

  /// Dump information
  virtual void Print() const
  {
    INFO << "charge ="   << std::setw(8) << std::left << charge_  << ", "  
         << "sumET_isol_ = " << std::setw(8) << std::left << sumET_isol_ << ", "
         << "sumPT_isol_ = " << std::setw(8) << std::left << sumPT_isol_;

    RecParticleFormat::Print();
  }

  /// Clear all information
  virtual void Reset()
  {
    charge_=false;
    sumET_isol_=0.;
    sumPT_isol_=0.;
  }

  /// Accessor to the electric charge 
  virtual const Float_t charge() const
  { if (charge_) return +1.; else return -1; }

  /// Accessor to sumET_isol
  virtual const Float_t sumET_isol() const
  { return sumET_isol_; }

  /// Accessor to sumPT_isol
  virtual const Float_t sumPT_isol() const
  { return sumPT_isol_; }

  /// Accessor to ET_PT
  virtual const Float_t ET_PT_isol() const
  { if (sumPT_isol_!=0) return sumET_isol_/sumPT_isol_;
    else return 0; }

};

#endif
