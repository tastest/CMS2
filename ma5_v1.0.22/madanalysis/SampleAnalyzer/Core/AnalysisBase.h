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


#ifndef ANALYSISBASE_h
#define ANALYSISBASE_h

// SampleAnalyzer headers
#include "DataFormat/EventFormat.h"
#include "DataFormat/SampleFormat.h"
#include "Core/Counter.h"
#include "Services/Physics.h"
#include "Services/logger.h"

// STL headers
#include <set>
#include <string>
#include <cmath>

// ROOT headers
#include <TTree.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TFile.h>
#include <TROOT.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TVector.h>
#include <TClonesArray.h>


// initializing MACRO 
#define INIT_ANALYSIS(CLASS,NAME) public: CLASS() {setName(NAME);} virtual ~CLASS() {} private:

class AnalysisBase
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 public :

  /// name of the analysis
  std::string name_;

  /// output ROOT file name
  std::string outputName_;

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public :

  /// Constructor without argument 
  AnalysisBase()
  { name_="unknown"; }

  /// Destructor
  virtual ~AnalysisBase()
  { }

  /// Initialize (common part to all analyses)
  void PreInitialize(const std::string& outputName)
  { outputName_ = outputName; }

  /// Initialize (specific to the analysis)
  virtual void Initialize()=0;

  /// Finalize
  void PreFinalize(const SampleFormat& summary, 
                   const std::vector<SampleFormat>& samples)
  { }
  virtual void Finalize(const SampleFormat& summary, 
                        const std::vector<SampleFormat>& samples)=0;

  /// Execute
  void PreExecute(const SampleFormat& mySample,
                  const EventFormat& myEvent)
  { PHYSICS->SetFinalState(myEvent.mc());
    PHYSICS->SetInitialState(myEvent.mc());
 }

  virtual void Execute(const SampleFormat& mySample,
                       const EventFormat& myEvent)=0;

  /// Accessor to analysis name
  const std::string name() const {return name_;}

  /// Mutator to analysis name
  void setName(const std::string& Name) {name_=Name;}

 protected :


};



#endif
