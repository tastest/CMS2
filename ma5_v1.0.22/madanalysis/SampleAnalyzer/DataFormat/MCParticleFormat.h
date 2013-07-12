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


#ifndef MCParticleFormat_h
#define MCParticleFormat_h

// STL headers
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

// SampleAnalyzer
#include "DataFormat/ParticleBaseFormat.h"
#include "Services/logger.h"

class LHEReader;
class LHCOReader;
class STDHEPreader;
class HEPMCReader;
class ROOTReader;

class MCParticleFormat : public ParticleBaseFormat
{
  friend class LHEReader;
  friend class LHCOReader;
  friend class STDHEPreader;
  friend class HEPMCReader;
  friend class ROOTReader;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:
   
  Float_t 		    ctau_;	    /// proper lifetime ctau (in mm)
  Float_t 		    spin_;	    /// cosine of the angle btwn the spin vector and
                              /// its 3-momentum, in the lab frame
  Int_t	          pdgid_;		  /// PDG numbering of the particle
  Short_t	        statuscode_;/// status code (-1 for initial state, 
                              /// 2 intermediate state, 1 final state)
  UInt_t 	        mothup1_;   /// first mother index
  UInt_t 	        mothup2_;   /// second mother index
  Int_t           extra_;

  MCParticleFormat *mother1_ ;  // mother particle
  MCParticleFormat *mother2_ ;  // mother particle


  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public :

  /// Constructor without arguments
  MCParticleFormat()
  { Reset(); }

  /// Destructor
  virtual ~MCParticleFormat()
  {}

  /// Clear all information
  virtual void Reset()
  {
    momentum_.SetPxPyPzE(0.,0.,0.,0.);
    ctau_=0.; spin_=0.; pdgid_=0; 
    statuscode_=0; mothup1_=0; mothup2_=0; mother1_=0; mother2_=0; extra_=0;
  }

  /// Print particle informations
  virtual void Print() const
  {
    INFO << "momentum=(" << std::setw(8) << std::left << momentum_.Px()
         << ", "<<std::setw(8) << std::left << momentum_.Py()  
         << ", "<<std::setw(8) << std::left << momentum_.Pz() 
         << ", "<<std::setw(8) << std::left << momentum_.E() << ") - ";
    INFO << "ctau=" << std::setw(8) << std::left << ctau_ << " - "
         << "spin=" << std::setw(8) << std::left << spin_ << " - "
         << "PDGID=" << std::setw(8) << std::left << pdgid_ << " - "
         << "StatusCode=" << std::setw(3) << std::left 
         << static_cast<signed int>(statuscode_) << " - ";

    if (mother1_==0) ERROR << "NoMum1" << " - ";
    else ERROR << "Mum1  " << " - ";

    if (mother2_==0) ERROR << "NoMum2" << std::endl;
    else ERROR << "Mum2  " << std::endl;
  }


  const Float_t& ctau() const {return ctau_;}
  const Float_t& spin() const {return spin_;}
  const Int_t& pdgid() const {return pdgid_;}
  const Short_t& statuscode() const {return statuscode_;}
  const MCParticleFormat* mother1() const {return mother1_;}
  const MCParticleFormat* mother2() const {return mother2_;}

  // mutators
  void setCtau(Float_t v)  {ctau_=v;}
  void setSpin(Float_t v)  {spin_=v;}
  void setPdgid(Int_t v)   {pdgid_=v;}
  void setStatuscode(Short_t v)  {statuscode_=v;}
  void setMomentum(const TLorentzVector& v)  {momentum_=v;}
  void setMothUp1(UInt_t v) {mothup1_=v;}
  void setMothUp2(UInt_t v) {mothup2_=v;}

};

#endif
