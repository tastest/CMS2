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


#ifndef MCEventFormat_h
#define MCEventFormat_h

// STL headers
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

// SampleAnalyzer headers
#include "DataFormat/MCParticleFormat.h"
#include "Services/logger.h"

class LHEReader;
class LHCOReader;
class STDHEPreader;
class HEPMCReader;
class ROOTReader;


class MCEventFormat
{
  friend class LHEReader;
  friend class LHCOReader;
  friend class STDHEPreader;
  friend class HEPMCReader;
  friend class ROOTReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private : 

  UInt_t nparts_;       /// number of particles in the event
  UInt_t processId_;    /// identity of the current process
  Float_t weight_;      /// event weight
  Float_t scale_;       /// scale Q of the event
  Float_t alphaQED_;    /// ALPHA_em value used
  Float_t alphaQCD_;    /// ALPHA_s value used

  /// List of generated particles
  std::vector<MCParticleFormat> particles_;

  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public :

  /// Constructor withtout arguments
  MCEventFormat()
  { Reset(); }

  /// Destructor
  ~MCEventFormat()
  { }

  /// Accessor to the process identity
  const UInt_t&  processId() const {return processId_;}

  /// Accessor to the event weight
  const Float_t& weight()    const {return weight_;   }

  /// Accessor to the scale
  const Float_t& scale()     const {return scale_;    }

  /// Accessor to alpha_QED
  const Float_t& alphaQED()  const {return alphaQED_; }

  /// Accessor to alpha_QCD
  const Float_t& alphaQCD()  const {return alphaQCD_; }

  /// Accessor to the generated particle collection (read-only)
  const std::vector<MCParticleFormat>& particles() const {return particles_;}

  /// Accessor to the generated particle collection
  std::vector<MCParticleFormat>& particles() {return particles_;}

  /// Setting the process identity
  void setProcessId(UInt_t v)  {processId_=v;}

  /// Setting the event weight
  void setWeight   (Float_t v) {weight_=v;   }

  /// Setting the scale
  void setScale    (Float_t v) {scale_=v;    }

  /// Setting AlphaQED
  void setAlphaQED (Float_t v) {alphaQED_=v; }

  /// Setting AlphaQCD
  void setAlphaQCD (Float_t v) {alphaQCD_=v; }

  /// Clearing all information
  void Reset()
  { nparts_=0; processId_=0; weight_=0.;
    scale_=0.; alphaQED_=0.; alphaQCD_=0.;
    particles_.clear(); }

  /// Displaying data member values
  void Print() const
  {
    INFO << "nparts="      << nparts_
         << " - processId=" << processId_
         << " - weight="    << weight_
         << " - scale="     << scale_
         << " - alphaQED="  << alphaQED_
         << " - alphaQCD="  << alphaQCD_ << std::endl;
  }

  /// Giving a new particle
  MCParticleFormat* GetNewParticle()
  {
    particles_.push_back(MCParticleFormat());
    return &particles_.back();
  }
};


#endif
