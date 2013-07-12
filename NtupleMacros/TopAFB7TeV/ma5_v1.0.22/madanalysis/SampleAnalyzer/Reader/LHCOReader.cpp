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


#include "Reader/LHCOReader.h"
#include <sstream>
#include <cmath>

// -----------------------------------------------------------------------------
// ReadHeader
// -----------------------------------------------------------------------------
bool LHCOReader::ReadHeader(SampleFormat& mySample)
{
  EndOfFile_=false;

  mySample.InitializeRec();

  firstevent_=true;
  saved_=false;

  // Normal end
  return true;
}

// -----------------------------------------------------------------------------
// FinalizeHeader
// -----------------------------------------------------------------------------
bool LHCOReader::FinalizeHeader(SampleFormat& mySample)
{
  // Normal end
  return true;
}


// -----------------------------------------------------------------------------
// ReadEvent
// -----------------------------------------------------------------------------
bool LHCOReader::ReadEvent(EventFormat& myEvent, SampleFormat& mySample)
{
  if(EndOfFile_) return false;

  myEvent.InitializeRec();

  bool EndOfLoop = false;
  
// Declarging a new string for line
  std::string line;
  std::string tmp;
     
  do
    {
      std::stringstream str;
      str.str("");
      
      
      if (!ReadLine(line))
	{
	  EndOfFile_=true;
	  return true;
	}
      str << line;
      str >> tmp;

      if (saved_)
	{
	  FillEventInitLine(savedline_,myEvent);
	  saved_=false;
	}

      if(tmp=="0")
	{
	  if(firstevent_ )
	    {
	      FillEventInitLine(line,myEvent);
	    }
	  else 
	    {
	      EndOfLoop = true;
	      savedline_=line;
	      saved_=true;
	      continue;
	    }
	}

      FillEventParticleLine(line,myEvent);
      
      firstevent_=false;
    }
  while(!EndOfLoop);
  
  // Normal end
  return true; 
}


// -----------------------------------------------------------------------------
// FinalizeEvent
// -----------------------------------------------------------------------------
bool LHCOReader::FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent)
{
  // Normal end
  return true; 
}


// -----------------------------------------------------------------------------                                     
// FillEventInitLine                                     
// -----------------------------------------------------------------------------     
void LHCOReader::FillEventInitLine(const std::string& line, EventFormat& myEvent)
{
  std::stringstream str;
  std::string muf;
  str << line;
  str >> muf;
}
// -----------------------------------------------------------------------------                                     
// FillEventParticleLine                                     
// -----------------------------------------------------------------------------     
void LHCOReader::FillEventParticleLine(const std::string& line, EventFormat& myEvent)
{
  std::stringstream str;
  std::string firstchar;
  std::string muf;
  str << line;
  str >> firstchar;

  Double_t                tmp;      // temporary variable to fill in LorentzVector
  Double_t                eta;
  Double_t                phi;      // to define the MET
  Double_t                pt;
  Double_t                mass;

  str >> muf;
  if (muf=="0")
  {

  }
  else if (muf=="1")
  {
    RecLeptonFormat * elec = myEvent.rec()->GetNewElectron();
    str >> eta; 
    str >> phi; 
    str >> pt; 
    str >> mass; elec->momentum_.SetPtEtaPhiM(pt,eta,phi,mass);
    str >> tmp; 
    if(tmp<0) elec->charge_ = false;
    else elec->charge_= true ;
    str << tmp;
    str >> elec->HEoverEE_;
      
  }
  else if (muf=="2")
  {
    RecLeptonFormat * muon = myEvent.rec()->GetNewMuon();
    str >> eta; 
    str >> phi; 
    str >> pt; 
    str >> mass; muon->momentum_.SetPtEtaPhiM(pt,eta,phi,mass);
    str >> tmp; 
    if(tmp<0) muon->charge_ = false;
    else muon->charge_= true ;
    str >> tmp;
    str >> tmp;
    muon->sumPT_isol_=std::floor(tmp);
    Float_t ET_PT=tmp-muon->sumPT_isol_;
    Bool_t test=false;
    for (unsigned int j=0;j<5;j++)
    {
      ET_PT*=10;
      if (std::floor(ET_PT)==ET_PT)
      {
        test=true;
        break;
      }
    }
    if (test) muon->sumET_isol_=std::floor(ET_PT)*muon->sumPT_isol_;
    else muon->sumET_isol_=0;
  }
  else if (muf=="3")
  {
    RecTauFormat * tau = myEvent.rec()->GetNewTau();
    str >> eta; 
    str >> phi; 
    str >> pt; 
    str >> mass; tau->momentum_.SetPtEtaPhiM(pt,eta,phi,mass);
    str >> tmp; 
    if(tmp<0) tau->charge_ = false;
    else tau->charge_= true ;
    str >> tmp;
    str >> tmp; tau->HEoverEE_=tmp;
  }
  else if (muf=="4")
  {
    RecJetFormat * jet = myEvent.rec()->GetNewJet(); 
    str >> eta; 
    str >> phi; 
    str >> pt; 
    str >> mass; jet->momentum_.SetPtEtaPhiM(pt,eta,phi,mass);
    str >> tmp; jet->ntracks_=tmp; 
    str >> tmp;
    if ( tmp == 1. || tmp ==2.) jet->btag_=true;
    else jet->btag_ =false;
    str >> jet->HEoverEE_;
  }
  else if (muf=="6")
  {
    str >> tmp; 
    str >> phi; 
    str >> tmp;
    myEvent.rec()->MET_.met_.Set(tmp*cos(phi),tmp*sin(phi));
  }
  else if (firstchar!="0")
  {
    WARNING << "Unknown type of object : " << muf << std::endl;
  }
}
