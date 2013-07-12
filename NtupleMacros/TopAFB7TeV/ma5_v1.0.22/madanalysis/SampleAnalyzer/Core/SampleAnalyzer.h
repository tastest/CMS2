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


#ifndef SAMPLE_ANALYZER_H
#define SAMPLE_ANALYZER_H

// STL headers
#include <iostream>
#include <string>
#include <vector>


// SampleAnalyzer headers
#include "Core/Configuration.h"
#include "DataFormat/EventFormat.h"
#include "DataFormat/SampleFormat.h"
#include "Core/ReaderManager.h"
#include "Core/AnalysisManager.h"


class SampleAnalyzer
{
 private :
 
  const std::string&              analysisName_; 
  const Configuration&            cfg_;
  const std::vector<std::string>& inputs_;

  
 public:

  /// Constructor with arguments
  SampleAnalyzer(const std::string& analysisName,
                 const Configuration& cfg,
                 const std::vector<std::string>& inputs): 
    analysisName_(analysisName), cfg_(cfg), inputs_(inputs)
  { }

  /// Run the analysis of the files
  bool Run() const;

  /// Run the analysis of an individual file
  ULong_t IndividualRun(ReaderBase* & lhe,
                        AnalysisBase* & theanalysis,
                        SampleFormat& mySample) const;

  /// Produce ROOT output file
  bool ProduceOutput(SampleFormat& summary, std::vector<SampleFormat>& samples) const;
};


#endif
