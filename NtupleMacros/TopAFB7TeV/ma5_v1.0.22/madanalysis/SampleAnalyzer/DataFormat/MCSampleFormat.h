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


#ifndef MCSAMPLE_DATAFORMAT_H
#define MCSAMPLE_DATAFORMAT_H

// STL headers
#include <map>
#include <iostream>
#include <vector>

// SampleAnalyzer headers
#include "DataFormat/MCProcessFormat.h"
#include "Services/logger.h" 

class LHEReader;
class LHCOReader;
class HEPMCReader;
class ROOTReader;
class SampleAnalyzer;

class MCSampleFormat
{
  friend class LHEReader;
  friend class LHCOReader;
  friend class HEPMCReader;
  friend class ROOTReader;
  friend class SampleAnalyzer;

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private:

  // ---------------------- physics info -------------------------
  std::pair<Int_t,Int_t>      beamPDGID_;    
  std::pair<Float_t,Float_t>  beamE_;        
  std::pair<UInt_t,UInt_t>    beamPDFauthor_;
  std::pair<UInt_t,UInt_t>    beamPDFID_;
  Int_t                       weightMode_;
  UInt_t                      nProcesses_;
  std::vector<ProcessFormat>  processes_;

  // ----------------------- file info ---------------------------
  Bool_t      MadgraphTag_;       // Is the file produced by Madgraph
  Bool_t      MadgraphPythiaTag_; // Is the file produced by MadgraphPythia
  Float_t     xsection_;          // cross-section (fb^{-1})
  Float_t     xsection_error_;    // cross-section error (fb^{-1})

  // -------------------------------------------------------------
  //                      method members
  // -------------------------------------------------------------
 public :

  /// Constructor withtout arguments
  MCSampleFormat()
  { Reset(); }

  /// Destructor
  ~MCSampleFormat()
  { }

  /// Clear all the content
  void Reset()
  {
    beamPDGID_         = std::make_pair(0,0); 
    beamE_             = std::make_pair(0,0); 
    beamPDFauthor_     = std::make_pair(0,0); 
    beamPDFID_         = std::make_pair(0,0);
    weightMode_        = 0.; 
    MadgraphPythiaTag_ = false; 
    MadgraphTag_       = false;
    xsection_          = 0.; 
    xsection_error_    = 0.;
    processes_.clear();
  }
  
  /// Accessor to PDG ID of the intial partons
  const std::pair<Int_t,Int_t>& beamPDGID() const
  { return beamPDGID_; }

  /// Accessor to the beam energy
  const std::pair<Float_t,Float_t>& beamE() const
  { return beamE_; } 

  /// Accessor to the PDF authors
  const std::pair<UInt_t,UInt_t>& beamPDFauthor() const
  { return beamPDFauthor_; }

  /// Accessor to the PDF identity
  const std::pair<UInt_t,UInt_t>& beamPDFID() const
  { return beamPDFID_; }

  /// Accessor to the weight mode
  const Int_t& weightMode() const
  { return weightMode_; }

  /// Accessor to the MadGraph tag
  const Bool_t&       MadgraphTag()      const {return MadgraphTag_; }

  /// Accessor to the MadGraphPythia tag
  const Bool_t&       MadgraphPythiaTag()const {return MadgraphPythiaTag_;}

  /// Accessor to the sample cross section mean value
  const Float_t&      xsection()         const {return xsection_; }

  /// Accessor to the sample cross section error
  const Float_t&      xsection_error()   const {return xsection_error_; }

  /// Accessor to the process collection (read-only)
  const std::vector<ProcessFormat>& processes() const
  { return processes_; }

  /// Accessor to the process collection
  std::vector<ProcessFormat>& processes()
  { return processes_; }

  /// Set the PDG ID of the intial partons
  void setBeamPDGID(Int_t a, Int_t b) 
  {beamPDGID_=std::make_pair(a,b); }

  /// Set the beam energy
  void setBeamE(Float_t a, Float_t b)
  {beamE_=std::make_pair(a,b); }

  /// Set the PDF authors
  void setBeamPDFauthor(UInt_t a, UInt_t b)
  {beamPDFauthor_=std::make_pair(a,b); }

  /// Set the the PDF identity
  void setBeamPDFid(UInt_t a, UInt_t b)
  {beamPDFID_=std::make_pair(a,b); }

  /// Set the weight mode
  void setWeightMode(Int_t v) 
  {weightMode_=v;}

  /// Enable MadgraphTag
  void enableMadgraphTag()        
  {MadgraphTag_=true;}

  /// Disable MadgraphTag
  void disableMadgraphTag()       
  {MadgraphTag_=false;}

  /// Enable MadgraphTag
  void enableMadgraphPythiaTag()
  {MadgraphPythiaTag_=true;}

  /// Disable MadgraphTag
  void disableMadgraphPythiaTag()
  {MadgraphPythiaTag_=false;}

  /// Set event cross section mean value
  void set_xsection(float value)
  { xsection_=value;
    if (MadgraphPythiaTag_) xsection_*=1e9;}

  /// Set event cross section error
  void set_xsection_error(float value)
  { xsection_error_=value;
    if (MadgraphPythiaTag_) xsection_error_*=1e9;}

  /// Displaying subtitle for file
  void printSubtitle() const
  {
    INFO << "    * \t => sample produced by ";
    if (MadgraphPythiaTag_) INFO << "MadGraph + Pythia.";
    else if (MadgraphTag_) INFO << "MadGraph.";
    else INFO << "an unknown generator (cross section assumed in pb).";
    INFO<<std::endl;
  }

  /// Giving a new process entry
  ProcessFormat* GetNewProcess()
  {
    processes_.push_back(ProcessFormat());
    return &processes_.back();
  }

};


#endif
