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


//STL headers
#include <sstream>

//SampleHeader headers
#include "Reader/HEPMCReader.h"
#include "Services/logger.h"


// -----------------------------------------------------------------------------
// ReadHeader
// -----------------------------------------------------------------------------
bool HEPMCReader::ReadHeader(SampleFormat& mySample)
{
  // Reset the saved line 
  savedline_="";
  
  // Initialize MC
  mySample.InitializeMC();

  // Skipping header line until first event line
  std::string firstWord;
  std::string line;

  while(firstWord!="E")
  {
    // Getting the next non-empty line
    if (!ReadLine(line)) return false;

    // Splitting the line in words
    std::stringstream str;
    str << line;

    // Extracting the first word
    str >> firstWord;
  }
  
  savedline_  = line;

  // Normal end
  return true;
}


// -----------------------------------------------------------------------------
// FinalizeHeader
// -----------------------------------------------------------------------------
bool HEPMCReader::FinalizeHeader(SampleFormat& mySample)
{
  return true;
}

// -----------------------------------------------------------------------------
// ReadEvent
// -----------------------------------------------------------------------------
bool HEPMCReader::ReadEvent(EventFormat& myEvent, SampleFormat& mySample)
{
  // Initializing MC event
  myEvent.InitializeMC();

  Bool_t eventOnGoing=false;

  // Read the saved line
  if (savedline_!="") 
  {
    FillEvent(savedline_, myEvent, mySample);
    eventOnGoing=true;
    savedline_="";
  }

  bool endEvent=false;
  
  // Loop over particle
  while(!endEvent)
  {
    std::string line;

    // Getting a line from the file
    if (!ReadLine(line))
    {
      if (eventOnGoing) return true; else return false;
    }

    // Splitting the line in words
    std::stringstream str;
    str << line;

    // Extracting the first word
    std::string firstWord;
    str >> firstWord;

    // Is next event ?
    if (firstWord=="E")
    {
      savedline_  = line;
      return true;
    }
    else
    {
      // Decoding the line
      endEvent=!FillEvent(line, myEvent, mySample);
      eventOnGoing=true;
    }
  }

  // Normal end 
  return true; 
}


// -----------------------------------------------------------------------------
// FinalizeEvent
// -----------------------------------------------------------------------------
bool HEPMCReader::FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent)
{

  double value=0;
  double err_value=0;
  if(nevents_!=0) 
  {
    value=(event_xsection_/static_cast<double>(nevents_));
    err_value=(event_xsection_err_/static_cast<double>(nevents_));
  }

  mySample.mc()->set_xsection(value);
  mySample.mc()->set_xsection_error(err_value);
  // Normal end 

  return true; 
}

//------------------------------------------------------------------------------
// FillEventHeader
//------------------------------------------------------------------------------
Bool_t HEPMCReader::FillEvent(const std::string& line,
                              EventFormat& myEvent, 
                              SampleFormat& mySample)
{
  // Splitting line in words
  std::stringstream str;
  str << line ;

  // Getting the first word
  std::string firstWord;
  str >> firstWord;

  // Action to do according to the first word
  if(firstWord=="E")
  {
    FillEventInformations(line, myEvent);
  }
  else if (firstWord=="N")
  {
    // Weight names
    // WARNING << "HEPMC weight names not considered." << std::endl;
  }
  else if (firstWord=="U")
  {
    // Event units
    FillEventUnits(line);
  }
  else if (firstWord=="C")
  {
    // Cross section
    FillEventXS(line);
  }
  else if (firstWord=="H")
  {
    // HeavyIon line
    // WARNING << "HEPMC Heavy Ions line not considered." << std::endl;
  }
  else if (firstWord=="F")
  {
    // PDF Info
    FillEventPDFInfo(line,mySample);
  }
  else if (firstWord=="V")
  {
    // Vertex line 
    FillEventVertexLine(line,myEvent);
  }
  else if (firstWord=="P")
  {
    // Particle Line
    FillEventParticleLine(line,myEvent);
  }
  else if (firstWord=="HepMC::IO_GenEvent-END_EVENT_LISTING")
  {
    return false;
  }
  else
  {
    // ignore other cases
    WARNING << "HEPMC linecode unknown" << std::endl;
  }

  // Normal end
  return true;
}

// -----------------------------------------------------------------------------
// FillEventInformations
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventInformations(const std::string& line,
                                  EventFormat& myEvent)
{
  std::stringstream str;
  str << line;
  std::string firstc;
  int tmp=0;

  str >> firstc;
  str >> tmp;                    //event number
  str >> myEvent.mc()->scale_;
  str >> myEvent.mc()->alphaQCD_;
  str >> myEvent.mc()->alphaQED_;
  str >> myEvent.mc()->processId_;
  str >> tmp;                    // process vertex
  str >> tmp;                    // vertices number
  str >> myEvent.mc()->nparts_;
  str >> myEvent.mc()->weight_;
}

// -----------------------------------------------------------------------------
// FillEventUnits
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventUnits(const std::string& line)
{
  std::stringstream str;
  str << line;
  
  std::string tmp;
  
  str >> tmp;
  str >> tmp;
  
  if (tmp=="GEV") units_=1;
  else if (tmp=="MEV") units_=0.001;
  else if (tmp=="KEV") units_=0.000001;
}


// -----------------------------------------------------------------------------
// FillEventXS
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventXS(const std::string& line)
{
  std::stringstream str;
  str << line;
  double xsectmp=0;
  double xsectmp_err=0;
  std::string firstc;
  
  nevents_++;
  
  str >> firstc;
  str >> xsectmp;
  event_xsection_+=xsectmp;
  str >> event_xsection_err_;
  event_xsection_err_+=xsectmp_err;
}

// -----------------------------------------------------------------------------
// FillEventPDFInfo
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventPDFInfo(const std::string& line, SampleFormat& mySample)
{
  std::stringstream str;
  str << line;
  std::string firstc;
  int tmp=0;
  str >> firstc;
  str >> mySample.mc()->beamPDGID_.first;
  str >> mySample.mc()->beamPDGID_.second;
  str >> mySample.mc()->beamE_.first;
  str >> mySample.mc()->beamE_.first;
  str >> tmp;
  str >> tmp;
  str >> tmp;
  str >> mySample.mc()->beamPDFID_.first;
  str >> mySample.mc()->beamPDFID_.first;
}

// -----------------------------------------------------------------------------
// FillEventParticleLine
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventParticleLine(const std::string& line,
                                      EventFormat& myEvent)
{
  std::stringstream str;
  str << line;

  double   	tmp;    // temporary variable to fill in LorentzVector

  // Get a new particle
  MCParticleFormat * part = myEvent.mc()->GetNewParticle();
  char linecode;
  str >> linecode;
  UInt_t partnum;
  str >> partnum;
  str >> part->pdgid_;
  str >> tmp; part->momentum_.SetPx(tmp*units_);
  str >> tmp; part->momentum_.SetPy(tmp*units_);
  str >> tmp; part->momentum_.SetPz(tmp*units_);
  str >> tmp; part->momentum_.SetE(tmp*units_);
  str >> tmp; 
  str >> part->statuscode_;
  str >> tmp; 
  str >> tmp; 
  str >> part->extra_;

  SetMother(part, myEvent);
}

// -----------------------------------------------------------------------------
// FillEventVertexLine
// -----------------------------------------------------------------------------
void HEPMCReader::FillEventVertexLine(const std::string& line, EventFormat& myEvent)
{
  std::stringstream str;
  str << line;

  int tmp=0;
  
  str >> current_vertex_.barcode_;
  str >> tmp;
  str >> tmp;
  str >> tmp;
  str >> current_vertex_.ctau_;
}


//--------------------------------------------------------------------------
// SetMother
//--------------------------------------------------------------------------
void HEPMCReader::SetMother(MCParticleFormat* part, EventFormat& myEvent)
{
  if (myEvent.mc()->particles().size()>1) return;
  
  unsigned int nmother=0;
  
  for (unsigned int i =0; i<(myEvent.mc()->particles().size()-1);i++)
  {
    if(myEvent.mc()->particles()[i].extra_==current_vertex_.barcode_)
    {
      nmother++;
      if(nmother==1)      part->mother1_=&(myEvent.mc()->particles()[i]);
      else if(nmother==2) part->mother2_=&(myEvent.mc()->particles()[i]);
      else WARNING << "Number of mothers greather than 2" << std::endl;
    }
  }
}
