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


#include "Core/SampleAnalyzer.h"
#include "Services/logger.h"


// ----------------------------------------------------------------------------
// Run function
// ----------------------------------------------------------------------------
bool SampleAnalyzer::Run() const
{
  // -----------------------------------------------------------------------
  //                   GETTING THE PROPER ANALYSIS
  // -----------------------------------------------------------------------

  // Building the list of analyses
  AnalysisManager analyses;
  analyses.BuildTable();

  // Getting the analysis
  AnalysisBase* myAnalysis = 0;
  if (cfg_.analysis_selected) myAnalysis = analyses.Get(cfg_.analysis_name);
  else myAnalysis = analyses.ChoiceAnalysis();

  if (myAnalysis==0)
    {
      ERROR << "analysis called '" << cfg_.analysis_name << "' is not found" 
            << std::endl;
      return false;
    }

  // Initialize (common part to all analyses)
  myAnalysis->PreInitialize(analysisName_);
  ULong_t ncounter=0;

  // Initialize (specific to the analysis)
  myAnalysis->Initialize();

  // Declaring format
  std::vector<SampleFormat> samples;

  // -----------------------------------------------------------------------
  //                          LOOP ON LHE SAMPLES
  // -----------------------------------------------------------------------
  for (unsigned int i=0;i<inputs_.size();i++)
  {

    // Progression bar
    INFO << "    * " << i+1 <<"/" << inputs_.size() << " ";  
    INFO << " " << inputs_[i] << std::endl;

    // Data format
    SampleFormat mySample;
    mySample.setName(inputs_[i]);

    // Building the list of readers
    ReaderManager readers;
    readers.BuildTable();

    // Find an appropiate reader for the file
    ReaderBase* myReader = readers.Get(inputs_[i]);
    if (myReader==0)
    {
      ERROR << "the format of the input file is not supported. Sorry."
            << std::endl;
      continue;
    }

    // Initialize the reader
    myReader->Initialize(inputs_[i], cfg_);

    // Read the event
    ncounter+=IndividualRun(myReader,myAnalysis,mySample);

    // Finalize the reader
    myReader->Finalize();

    // Saving information related to the sample
    samples.push_back(mySample);
  }

  // -----------------------------------------------------------------------
  //                    ADDING MEAN RESULTS TO FILES
  // -----------------------------------------------------------------------
  SampleFormat summary;
  ProduceOutput(summary,samples);

  // -----------------------------------------------------------------------
  //                    PRODUCING OUTPUT ROOT FILE
  // -----------------------------------------------------------------------
  INFO << "    * Creating a root file..." << std::endl;
  myAnalysis->PreFinalize(summary,samples);
  myAnalysis->Finalize(summary,samples);

  // -----------------------------------------------------------------------
  //                      DUMP NUMBER OF EVENT
  // -----------------------------------------------------------------------
  INFO << "    * Total number of processed events: ";
  INFO << ncounter << "." << std::endl;

  for (unsigned int i=0;i<samples.size();i++)
  { samples[i].Delete(); }


}
 

// ----------------------------------------------------------------------------
// IndividualRun function
// ----------------------------------------------------------------------------
ULong_t SampleAnalyzer::IndividualRun(ReaderBase* &   myReader,
                                      AnalysisBase* & myAnalysis,
                                      SampleFormat&   mySample) const
{ 
  bool test=true;
  ULong_t ncounter=0;

  // Read the header block
  test=myReader->ReadHeader(mySample);
  if (!test)
  {
    ERROR << "No header has been found." << std::endl;
    exit(1);
  }

  // Finalize the header block
  myReader->FinalizeHeader(mySample);

  // Dump the header block
  if(mySample.mc()!=0)
  {
    mySample.mc()->printSubtitle();
  }

  // Process the event block 
  do
  {
    // Creating an event
    EventFormat  myEvent;

    // Read an event
    test=myReader->ReadEvent(myEvent, mySample);

    // Finalize the event
    myReader->FinalizeEvent(mySample,myEvent);

    if (test) 
    {
      // Analyze the event
      myAnalysis->PreExecute(mySample,myEvent);  
      myAnalysis->Execute(mySample,myEvent);
      ncounter++;
    }
    myEvent.Delete();
  }
  while(test);

  // Dump number of read events
  INFO << "    * \t => Number of processed events: ";
  INFO << ncounter << "." << std::endl;

  // Saving number of read events
  mySample.setNEvents(ncounter);
  
  // Return number of events
  return ncounter;
}


// ----------------------------------------------------------------------------
// IndividualRun function
// ----------------------------------------------------------------------------
bool SampleAnalyzer::ProduceOutput(SampleFormat& summary,
                                   std::vector<SampleFormat>& samples ) const
{ 
  // Create a SampleFormat container for summary info
  summary.setName("FINAL");
  summary.InitializeMC();
  summary.InitializeRec();

  // Loop over samples
  for (unsigned int i=0;i<samples.size();i++)
  {
    summary.nevents_              += samples[i].nevents_;
    if(samples[i].mc()==0) {continue ;}
    summary.mc()->xsection_       += samples[i].mc()->xsection_ * 
                                   samples[i].nevents_;
    summary.mc()->xsection_error_ += samples[i].mc()->xsection_error_ * 
                                   samples[i].mc()->xsection_error_ *
                                   samples[i].nevents_ *
                                   samples[i].nevents_;
  }
  if (samples.size()!=0)
  {
    summary.mc()->xsection_       /= summary.nevents_;
    summary.mc()->xsection_error_  = sqrt(summary.mc()->xsection_error_)
                                   / summary.nevents_;
  }

  // Normal end
  return true;
}




