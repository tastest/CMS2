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


#ifndef RecParticleFormat_h
#define RecParticleFormat_h

// STL headers
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

// MCParticleDataFormat
#include "DataFormat/ParticleBaseFormat.h"
#include "Services/logger.h"

class MCParticleFormat;
class LHEReader;
class LHCOReader;

class RecParticleFormat : public ParticleBaseFormat
{
  friend class LHEReader;
  friend class LHCOReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 protected:
   
  Float_t 	        HEoverEE_; /// hadronic energy over electromagnetic energy
  MCParticleFormat* mc_ ;      /// mother generated particle

  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public:

  /// Constructor without arguments
  RecParticleFormat()
  { Reset(); }

  /// Destructor
  virtual ~RecParticleFormat()
  {}

  /// Clear all information
  virtual void Reset()
  {
    momentum_.SetPxPyPzE(0.,0.,0.,0.);
    HEoverEE_=0.; 
    mc_=0;
  }

  /// Print particle informations
  virtual void Print() const
  {
    INFO << "momentum=(" << std::setw(8) << std::left << momentum_.Px()
         << ", "<<std::setw(8) << std::left << momentum_.Py()  
         << ", "<<std::setw(8) << std::left << momentum_.Pz() 
         << ", "<<std::setw(8) << std::left << momentum_.E() << ") - "
         << "EHoverEE=" << std::setw(8) << std::left << HEoverEE_
	       << " - ";

    if (mc_==0) ERROR << "NoMCmum" << " - ";
    else ERROR << "Mum1  " << " - ";
  }

  /// Accessor to matched Monte Carlo particle 
  const MCParticleFormat* mc() {return mc_;}

  /// Accessor to hadronic energy / electromagnetic energy ratio
  const Float_t& HEoverEE() const {return HEoverEE_;}

  /// Accessor to electromagnetic energy / hadronic energy ratio
  const Float_t EEoverHE() const 
  {
    if (HEoverEE_!=0) return 1./HEoverEE_; 
    else return 0.;
  }

  /// Accessor to the number of tracks
  virtual const UShort_t ntracks() const
  { return 0; }

  /// Accessor to the isolation tag
  virtual const Bool_t isolated() const
  { return false; }

  /// Accessor to the electric charge
  virtual const Float_t charge() const
  { return 0.; }

};

#endif
