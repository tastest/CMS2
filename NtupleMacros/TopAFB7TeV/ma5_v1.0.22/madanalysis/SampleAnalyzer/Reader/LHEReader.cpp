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


// STL headers
#include <sstream>

// SampleHeader headers
#include "Reader/LHEReader.h"
#include "Services/logger.h"

// -----------------------------------------------------------------------------
// ReadHeader
// -----------------------------------------------------------------------------
bool LHEReader::ReadHeader(SampleFormat& mySample)
{
  // Initiliaze MC
  mySample.InitializeMC();

  // Declaring a new string for line
  std::string line;

  // Read line by line the file until tag <init>
  bool EndOfLoop=false;
  do 
  { 
    if (!ReadLine(line)) return false;
    if ((line.find("<MGGenerationInfo>")!=std::string::npos) ||
         (line.find("<mgversion>")!=std::string::npos) )
      mySample.mc()->enableMadgraphTag();
    if ( (line.find("<MGPythiaCard>")!=std::string::npos) ||
         (line.find("<mgpythiacard>")!=std::string::npos) ) 
      mySample.mc()->enableMadgraphPythiaTag();
    EndOfLoop = (line.find("<init>")!=std::string::npos);
  }
  while(!EndOfLoop);

  // Read line by line the file until tag </init>
  EndOfLoop=false;
  bool first=true;
  do 
  { 
    if (!ReadLine(line)) return false;
    EndOfLoop = (line.find("</init>")!=std::string::npos);
    if (!EndOfLoop)
    {
      if (first) FillHeaderInitLine(line,mySample);
      else FillHeaderProcessLine(line,mySample);
    }
    first=false;
  }
  while(!EndOfLoop);

  // Read line by line the file until tag <event>
  EndOfLoop=false;
  do 
  { 
    if (!ReadLine(line)) return false;
    if ((line.find("<MGGenerationInfo>")!=std::string::npos)||
         (line.find("<mgversion>")!=std::string::npos) )
      mySample.mc()->enableMadgraphTag();
    if ((line.find("<MGPythiaCard>")!=std::string::npos) ||
         (line.find("<mgpythiacard>")!=std::string::npos) )
      mySample.mc()->enableMadgraphPythiaTag();
    EndOfLoop = (line.find("<event>")!=std::string::npos);
  }
  while(!EndOfLoop);

  // Normal end
  firstevent_=true;
  return true;
}


// -----------------------------------------------------------------------------
// FinalizeHeader
// -----------------------------------------------------------------------------
bool LHEReader::FinalizeHeader(SampleFormat& mySample)
{
  // Computing xsection an its error for the sample
  float xsection = 0;
  float xerror = 0;
  for (unsigned int i=0;i<mySample.mc()->processes().size();i++)
  {
    xsection += mySample.mc()->processes()[i].xsectionMean();
    xerror   += mySample.mc()->processes()[i].xsectionError();
  }

  // Filling xsection and its error
  mySample.mc()->set_xsection( xsection );
  mySample.mc()->set_xsection_error( xerror );

  // Normal end 
  return true;
}


// -----------------------------------------------------------------------------
// ReadEvent
// -----------------------------------------------------------------------------
bool LHEReader::ReadEvent(EventFormat& myEvent, SampleFormat& mySample)
{
  // Initiliaze MC
  myEvent.InitializeMC();

  // Declaring a new string for line
  std::string line;
  bool EndOfLoop=false;

  // Read line by line the file until tag <event>
  if (!firstevent_)
  {
    EndOfLoop=false;
    do 
    { 
      if (!ReadLine(line)) return false;
      EndOfLoop = (line.find("<event>")!=std::string::npos);
    }
    while(!EndOfLoop);
  }

  // Read the particles
  EndOfLoop=false;
  firstevent_=false;
  bool first=true;
  do 
  { 
    if (!ReadLine(line)) return false;
    EndOfLoop = (line.find("</event>")!=std::string::npos);
    if (!EndOfLoop)
    {
      if (first) FillEventInitLine(line,myEvent);
      else FillEventParticleLine(line,myEvent);
    }
    first=false;
  }
  while(!EndOfLoop);

  // Normal end
  return true; 
}


// -----------------------------------------------------------------------------
// FinalizeEvent
// -----------------------------------------------------------------------------
bool LHEReader::FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent)
{
  // Mother pointer assignment
  for (unsigned int i=0; i<myEvent.mc()->particles_.size();i++)
  {
    unsigned int index1=myEvent.mc()->particles_[i].mothup1_;
    unsigned int index2=myEvent.mc()->particles_[i].mothup2_;
    if (index1!=0 && index2!=0)
    {
      if (index1>=myEvent.mc()->particles_.size() ||
          index2>=myEvent.mc()->particles_.size())
      {
        /*        for (unsigned int num=0;num<myEvent.mc()->particles_.size();num++)
        {
          std::cout << "num=" << num+1 << " | id=" << 
                       myEvent.mc()->particles_[num].pdgid_ 
                    << " ; mothup1_=" << myEvent.mc()->particles_[num].mothup1_ 
                    << " ; mothup2_=" << myEvent.mc()->particles_[num].mothup2_ << std::endl;
                    }*/
        
        ERROR << "mother index is greater to nb of particles" << std::endl;
        ERROR << " - index1 = " << index1 << std::endl;
        ERROR << " - index2 = " << index2 << std::endl;
        ERROR << " - particles.size() " << myEvent.mc()->particles_.size() << std::endl;
        exit(1);
      }

      myEvent.mc()->particles_[i].mother1_ = &myEvent.mc()->particles_[index1-1];
      myEvent.mc()->particles_[i].mother2_ = &myEvent.mc()->particles_[index2-1];
    }
  }

  // Normal end
  return true; 
}


// -----------------------------------------------------------------------------
// FillHeaderInitLine
// -----------------------------------------------------------------------------
void LHEReader::FillHeaderInitLine(const std::string& line, 
                                   SampleFormat& mySample)
{
  std::stringstream str;
  str << line;

  str >> mySample.mc()->beamPDGID_.first;
  str >> mySample.mc()->beamPDGID_.second;
  str >> mySample.mc()->beamE_.first;
  str >> mySample.mc()->beamE_.second;
  str >> mySample.mc()->beamPDFauthor_.first;
  str >> mySample.mc()->beamPDFauthor_.second;
  str >> mySample.mc()->beamPDFID_.first;
  str >> mySample.mc()->beamPDFID_.second;
  str >> mySample.mc()->weightMode_;
  str >> mySample.mc()->nProcesses_;
}


// -----------------------------------------------------------------------------
// FillHeaderProcessLine
// -----------------------------------------------------------------------------
void LHEReader::FillHeaderProcessLine(const std::string& line,
                                      SampleFormat& mySample)
{
  std::stringstream str;
  str << line;
  float tmp;

  // Get a new process
  ProcessFormat * proc = mySample.mc()->GetNewProcess();

  str >> proc->xsectionMean_;
  str >> proc->xsectionError_;
  str >> proc->weightMax_;
  str >> proc->processId_;
}


// -----------------------------------------------------------------------------
// FillEventInitLine
// -----------------------------------------------------------------------------
void LHEReader::FillEventInitLine(const std::string& line,
                                  EventFormat& myEvent)
{
  std::stringstream str;
  str << line;

  str >> myEvent.mc()->nparts_;
  str >> myEvent.mc()->processId_;
  str >> myEvent.mc()->weight_;
  str >> myEvent.mc()->scale_;
  str >> myEvent.mc()->alphaQED_;
  str >> myEvent.mc()->alphaQCD_;
}


// -----------------------------------------------------------------------------
// FillEventParticleLine
// -----------------------------------------------------------------------------
void LHEReader::FillEventParticleLine(const std::string& line,
                                      EventFormat& myEvent)
{
  std::stringstream str;
  str << line;

  signed int 	color1;	// color 1 not stored 
  signed int	color2;	// color 2 not stored
  double   		tmp;	  // temporary variable to fill in LorentzVector

  // Get a new particle
  MCParticleFormat * part = myEvent.mc()->GetNewParticle();

  str >> part->pdgid_;
  str >> part->statuscode_;
  str >> part->mothup1_;
  str >> part->mothup2_;
  str >> color1;
  str >> color2;
  str >> tmp; part->momentum_.SetPx(tmp);
  str >> tmp; part->momentum_.SetPy(tmp);
  str >> tmp; part->momentum_.SetPz(tmp);
  str >> tmp; part->momentum_.SetE(tmp);
  str >> tmp;
  str >> part->ctau_;
  str >> part->spin_;
}

