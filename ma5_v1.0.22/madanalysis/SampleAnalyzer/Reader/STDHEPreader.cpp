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


#include "Reader/STDHEPreader.h"
#include "Services/logger.h"

// -----------------------------------------------------------------------------
// Initialize
// -----------------------------------------------------------------------------
bool STDHEPreader::Initialize(const std::string& rawfilename,
                              const Configuration& cfg)
{
  if (!ReaderTextBase::Initialize(rawfilename, cfg)) return false;
  xdrinput_=new xdr_istream(*input_);

  return true;
}


// -----------------------------------------------------------------------------
// ReadHeader
// -----------------------------------------------------------------------------
void STDHEPreader::Reset()
{
  nevhept_=0;
  nhept_=0;
  isthept_.clear();
  idhept_.clear();
  jmohept_.clear();
  jdahept_.clear();
  phept_.clear();
  vhept_.clear();
}

// -----------------------------------------------------------------------------
// ReadHeader
// -----------------------------------------------------------------------------
bool STDHEPreader::ReadHeader(SampleFormat& mySample)
{
  // Initiliaze MC
  mySample.InitializeMC();

  if (!DecodeFileHeader(mySample)) return false;

  return true;
}


// -----------------------------------------------------------------------------
// ReadEvent
// -----------------------------------------------------------------------------
bool STDHEPreader::ReadEvent(EventFormat& myEvent, SampleFormat& mySample)
{
  // Initiliaze MC
  myEvent.InitializeMC();

  bool eventRead=false;
  while(!eventRead)
  {
    // Read blockid
    Int_t blockid=0;
    *xdrinput_ >> blockid;
    if (xdrinput_->eof()) return false;

    Int_t ntot=0;
    *xdrinput_ >> ntot;

    std::string version;
    *xdrinput_ >> version;

    if      (blockid==EVENTTABLE )       DecodeEventTable (version);
    else if (blockid==EVENTHEADER)       DecodeEventHeader(version);
    else if (blockid==MCFIO_STDHEPBEG ||
             blockid==MCFIO_STDHEPEND)   {DecodeSTDCM1 (version,mySample); }
    else if (blockid==MCFIO_STDHEP)      {DecodeEventData(version, myEvent);
                                          eventRead=true;}
    else
    {
      ERROR << "Block with the ID=" << blockid 
            << " is not managed by SampleAnalyzer" << std::endl;
      exit(1);
    }
  }

  return true;
}


// -----------------------------------------------------------------------------
// DecodeFileHeader
// -----------------------------------------------------------------------------
bool STDHEPreader::DecodeFileHeader(SampleFormat& mySample)
{
  // temporary variables used for reading the xdr format file
  std::string  tmps;
  Int_t  tmpi = 0;
  UInt_t tmpui = 0;

  // BlockID
  *xdrinput_ >> tmpi;
  if (tmpi != FILEHEADER)
  {
    ERROR << "header block not found" << std::endl;
    return false;
  }

  // Ntot
  *xdrinput_ >> tmpi;

  // STDHEP version
  *xdrinput_ >> tmps;
  SetVersion(tmps);
  if (version_==UNKNOWN) 
  {
    ERROR << "stdhep version unknown : '" << tmps << std::endl;
    return false;
  }

  // Title
  *xdrinput_ >> tmps;
  //std::cout << "title=" << tmps << std::endl;

  // Comment
  *xdrinput_ >> tmps;
  //std::cout << "comment=" << tmps << std::endl; 

  // Creation date
  *xdrinput_ >> tmps;
  //std::cout << "date=" << tmps << std::endl;

  // Closing date (only in version 2.01)
  if (version_=V21)
  {
    *xdrinput_ >> tmps;
    //std::cout << "cdate=" << tmps << std::endl;
  }

  // Expected number of events
  *xdrinput_ >> tmpui;
  //std::cout << "Nevents = " << tmpui << std::endl;

  // Number of events
  *xdrinput_ >> tmpui;
  //std::cout << "Nevents = " << tmpui << std::endl;

  // First table
  *xdrinput_ >> tmpui;
  //std::cout << "First table=" << tmpui << std::endl;

  // Dim Table
  UInt_t dimTable;
  *xdrinput_ >> dimTable;
  //std::cout << "Dim table=" << dimTable << std::endl;

  // Number of blocks
  UInt_t nBlocks;
  *xdrinput_ >> nBlocks;
  //std::cout << "N blocks = " << nBlocks << std::endl;

  // Number of NTuples
  UInt_t nNTuples = 0;
  if (version_!=V1)
  {
    *xdrinput_ >> nNTuples;
    //std::cout << "Nb NTuples = " << nNTuples << std::endl;
  }

  // Processing blocks extraction 
  if (nBlocks!=0)
  {
    // Extracting blocks
    std::vector<Int_t> blocks;
    *xdrinput_ >> blocks;  
 
    // Extracting block names
    for (UInt_t i=0;i<blocks.size();i++)
    {
      *xdrinput_ >> tmps;
      //std::cout << "Block " << i << " = " << tmps; 
    }
  }

  // Processing ntuple extraction (only in version 2)
  /*if (version_!=V1 && nNTuples!=0)
  {
    // Loop over ntuple
    for (unsigned int i=0;i<nNTuples;i++)
    {
      // Number of characters in title
      if ( !xdr_int(&xdr_, &tmpi) ) return false;
      unsigned int nctitle = tmpi;

      // Number of characters in category
      if ( !xdr_int(&xdr_, &tmpi) ) return false;
      unsigned int nccategory = tmpi;

      // IdRef
      if ( !xdr_int(&xdr_, &tmpi) ) return false;

      // uid
      if ( !xdr_int(&xdr_, &tmpi) ) return false;

      // Title of the Ntuple
      if ( !xdr_string(&xdr_, &tmps,nctitle) ) return false;

      // Category of the Ntuple
      if ( !xdr_string(&xdr_, &tmps,nccategory) ) return false;

      INFO << "ntu ... to finish" << std::endl;

    }
  }
  */

  return true;
}


// -----------------------------------------------------------------------------
// FinalizeHeader
// -----------------------------------------------------------------------------
bool STDHEPreader::FinalizeHeader(SampleFormat& mySample)
{


  return true;
}


// -----------------------------------------------------------------------------
// DecodeEventTable
// -----------------------------------------------------------------------------
bool STDHEPreader::DecodeEventTable(const std::string& evt_version)
{
  // Decoding the event
  Int_t idat=0;
  *xdrinput_ >> idat;

  UInt_t uidat=0;
  *xdrinput_ >> uidat; 

  // Extracting evtnums
  std::vector<Int_t> evtnums;
  *xdrinput_ >> evtnums;

  // Extracting storenums
  std::vector<Int_t> storenums;
  *xdrinput_ >> storenums;

  // Extracting runnums
  std::vector<Int_t> runnums;
  *xdrinput_ >> runnums;

  // Extracting trigMasks
  std::vector<UInt_t> NtrigMasks;
  *xdrinput_ >> NtrigMasks;

  // Extracting prtEvents
  std::vector<UInt_t> NptrEvents;
  *xdrinput_ >> NptrEvents;

  return true;
}


// -----------------------------------------------------------------------------
// DecodeEventHeader
// -----------------------------------------------------------------------------
bool STDHEPreader::DecodeEventHeader(const std::string& evt_version)
{
  Int_t evtnum=0;
  *xdrinput_ >> evtnum;

  Int_t storenums=0;
  *xdrinput_ >> storenums;

  Int_t runnum=0;
  *xdrinput_ >> runnum;

  Int_t trigMask=0;
  *xdrinput_ >> trigMask;

  UInt_t nBlocks=0;
  *xdrinput_ >> nBlocks;

  UInt_t dimBlocks=0;
  *xdrinput_ >> dimBlocks;

  UInt_t nNTuples=0;
  UInt_t dimNTuples=0;

  if (evt_version=="2.00")
  {
    *xdrinput_ >> nNTuples;
    *xdrinput_ >> dimNTuples;
  }

  // Processing blocks extraction 
  if (dimBlocks>0)
  {
    // Extracting blocks
    std::vector<Int_t> blocks;
    *xdrinput_ >> blocks;

    // Extracting blocks
    std::vector<UInt_t> ptrBlocks;
    *xdrinput_ >> ptrBlocks;
  }

  // Processing blocks extraction 
  if (dimNTuples>0 && evt_version=="2.00")
  {
    // Extracting blocks
    std::vector<Int_t> nTupleIds;
    *xdrinput_ >> nTupleIds;

    // Extracting blocks
    std::vector<UInt_t> ptrNTuples;
    *xdrinput_ >> ptrNTuples;
  }

  return true;
}


// -----------------------------------------------------------------------------
// DecodeSTDCM1
// -----------------------------------------------------------------------------
bool STDHEPreader::DecodeSTDCM1(const std::string& version, SampleFormat& mySample)
{

  Int_t nevtreq; 
  *xdrinput_ >> nevtreq;

  Int_t nevtgen; 
  *xdrinput_ >> nevtgen;

  Int_t nevtwrt;
  *xdrinput_ >> nevtwrt;

  Float_t stdecom;
  *xdrinput_ >> stdecom; 

  Float_t stdxsec;
  *xdrinput_ >> stdxsec;
  if (stdxsec!=0) mySample.mc()->set_xsection(stdxsec);

  Double_t stdseed1; 
  *xdrinput_ >> stdseed1;

  Double_t stdseed2;
  *xdrinput_ >> stdseed2;

  if (version.find("1.")==0 || version.find("2.")==0 || 
      version.find("3.")==0 || version.find("4.")==0 || 
      version.find("5.00")==0) return true;

  std::string tmps;
  *xdrinput_ >> tmps;
  *xdrinput_ >> tmps;

  if (version.find("5.00")==0 || version.find("5.01")==0 )
   return true;

  Int_t nevtlh=0;
  *xdrinput_ >> nevtlh;

  return true;

}


// -----------------------------------------------------------------------------
// DecodeEventFormat
// -----------------------------------------------------------------------------
bool STDHEPreader::DecodeEventData(const std::string& version, 
                                   EventFormat& myEvent)
{
  Reset();

  // Extracting the event number
  *xdrinput_ >> nevhept_;  

  // Extracting the number of particles
  *xdrinput_ >> nhept_;

  // Extracting isthept
  *xdrinput_ >> isthept_;

  // Extracting idhept
  *xdrinput_ >> idhept_;

  // Extracting jmohept
  *xdrinput_ >>  jmohept_;

  // Extracting jdahept
  *xdrinput_ >> jdahept_;

  // Extracting 
  *xdrinput_ >> phept_;

  // Extracting 
  *xdrinput_ >> vhept_;

  // Check the size of the collections
  if ( nhept_<0 ||
       nhept_!=isthept_.size()       || nhept_!=idhept_.size()        || 
       (2*nhept_)!=jmohept_.size()   || (2*nhept_)!=jdahept_.size()   ||
       (5*nhept_)!=phept_.size()     || (4*nhept_)!=vhept_.size() )
  {
    ERROR << "Inconsistent size of collections. "
          << "The file is probably corrupted" << std::endl;
    return false;
  }

  // Loop over particles
  for (unsigned int i=0;i<static_cast<unsigned int>(nhept_);i++)
  {
    // Get a new particle
    MCParticleFormat * part = myEvent.mc()->GetNewParticle();

    // Fill the data format
    part->pdgid_=idhept_[i];
    part->statuscode_=isthept_[i];
    part->mothup1_=jmohept_[2*i];
    part->mothup2_=jmohept_[2*i+1];
    part->momentum_.SetPx(phept_[5*i]);
    part->momentum_.SetPy(phept_[5*i+1]);
    part->momentum_.SetPz(phept_[5*i+2]);
    part->momentum_.SetE (phept_[5*i+3]);
  }

  return true;
}


// -----------------------------------------------------------------------------
// FinalizeEvent
// -----------------------------------------------------------------------------
bool STDHEPreader::FinalizeEvent(SampleFormat& mySample, EventFormat& myEvent)
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
        ERROR << "mother index is greater to nb of particles"
              << std::endl
              << " - index1 = " << index1 << std::endl
              << " - index2 = " << index2 << std::endl
              << " - particles.size() " << myEvent.mc()->particles_.size()
              << std::endl;
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
// Finalize
// -----------------------------------------------------------------------------
bool STDHEPreader::Finalize()
{
  if (!ReaderTextBase::Finalize()) return false;
  if (xdrinput_!=0) delete xdrinput_;
  return true;  
}


// -----------------------------------------------------------------------------
// SetVersion
// -----------------------------------------------------------------------------
void STDHEPreader::SetVersion(const std::string& version)
{
  if (version.size()<2)       version_=UNKNOWN;
  else if (version[0]==1)     version_=V1;
  else if (version=="2.01")   version_=V21;
  else if (version[0]==2)     version_=V2;
  else version_=UNKNOWN; 
}
