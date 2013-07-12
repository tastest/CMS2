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


#ifndef FAC_FORMAT_H
#define FAC_FORMAT_H

// ROOT header
#include <TChain.h>
#include <TLorentzVector.h>
#include <TObject.h>
#include <TVector2.h>

// STL header
#include <iostream>
#include <vector>


namespace FAC
{
  class MCLeptonFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  private:
    TLorentzVector p4_;
    std::vector<Int_t> mothers_;


    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------
  public:

    MCLeptonFormat()
    {
      reset();
    }

    MCLeptonFormat(double px, double py, double pz, double e)
    {
      p4_.SetPxPyPzE(px,py,pz,e);
    }

    virtual ~MCLeptonFormat()
    { }

    void reset()
    {
      p4_.SetPxPyPzE(0.,0.,0.,0.);
      mothers_.clear();
    }

    const TLorentzVector& momentum() const
    { return p4_; }

    void setMomentum(double px, double py, double pz, double e)
    { p4_.SetPxPyPzE(px,py,pz,e); }

    const std::vector<Int_t>& mothers() const
    { return mothers_; }

    std::vector<Int_t>& mothers()
    { return mothers_; }

    void addMother(Int_t mum)
    { mothers_.push_back(mum); }

    virtual void Print(Option_t *option="") const
    {
      std::cout << "           qgen(" << p4_.Px() << ", " << p4_.Py() << ", "
                << p4_.Pz() << ", " << p4_.E() << ")" << std::endl;
      std::cout << "          ";
      for (unsigned int i=0;i<mothers_.size();i++)
      {
        std::cout << " -> " << mothers_[i];
      }
      std::cout << std::endl;
 
    }

    ClassDef(MCLeptonFormat,1)

  };

  class ElectronFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  private:
    TLorentzVector p4_;
    Bool_t isolated_;
    Bool_t charge_;
    Float_t EHoverEE_;

  public:
    MCLeptonFormat mc;


    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------
  public:

    ElectronFormat()
    {
      reset();
    }

    ElectronFormat(double px, double py, double pz, double e)
    {
      reset();
      p4_.SetPxPyPzE(px,py,pz,e);
    }

    virtual ~ElectronFormat()
    { }

    void reset()
    {
      p4_.SetPxPyPzE(0.,0.,0.,0.);
      isolated_ = false;
      charge_ = false;
      EHoverEE_ = 0.;
      mc.reset();
    }

    const TLorentzVector& momentum() const
    { return p4_; }

    void setMomentum(double px, double py, double pz, double e)
    { p4_.SetPxPyPzE(px,py,pz,e); }

    const Bool_t& isolated() const
    { return isolated_; }

    void setIsolated(Bool_t isolated)
    { isolated_=isolated; }

    Float_t charge() const
    { if (charge_) return +1.;
      else return -1.; }

    void setCharge(Float_t value)
    { if (value>0) charge_=true;
      else charge_=false; }

    const Float_t& ehoveree() const
    { return EHoverEE_; }

    void setEhoveree(Float_t EHoverEE)
    { EHoverEE_=EHoverEE; }

    virtual void Print(Option_t *option="") const
    {
      std::cout << "Electron : q(" << p4_.Px() << ", " << p4_.Py() << ", " 
                << p4_.Pz() << ", " << p4_.E() << ")" << std::endl;
      std::cout << "           isolated="<<isolated_<<", charge=" << charge_
                << ", EHoverEE="<<EHoverEE_ << std::endl;
      mc.Print();
    }

    ClassDef(ElectronFormat,1)

  };


  class GenJetFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  protected:
    TLorentzVector p4_;
    Char_t partonFlavour_;

    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------

  public :

    GenJetFormat()
    {
      reset();
    }

    GenJetFormat(double px, double py, double pz, double e)
    {
      reset();
      p4_.SetPxPyPzE(px,py,pz,e);
    }

    virtual ~GenJetFormat()
    { }

    virtual void reset()
    {
      p4_.SetPxPyPzE(0.,0.,0.,0.);
      partonFlavour_=0;
    }

    const TLorentzVector& momentum() const
    { return p4_; }

    void setMomentum(double px, double py, double pz, double e)
    { p4_.SetPxPyPzE(px,py,pz,e); }


    const Char_t& partonFlavour() const
    { return partonFlavour_; }

    void setPartonFlavour(Char_t d)
    { partonFlavour_=d; }

    ClassDef(GenJetFormat,1)

  };

  class GenParticleFormat : public TObject
  {

  private:
    TLorentzVector p4_;
    Int_t pid_;
    std::pair<Int_t,Int_t> mothers_;

  public:

    GenParticleFormat()
    {
      p4_.SetPxPyPzE(0,0,0,0);
      pid_=0;
      mothers_=std::make_pair(0,0);
    }

    virtual ~GenParticleFormat()
    { }

    GenParticleFormat(const double& px, const double& py, const double& pz,
                      const double& e,  const signed int& pid)
    {
      p4_.SetPxPyPzE(px,py,pz,e);
      pid_=pid;
    }

    const Int_t& pid() const
    { return pid_; }

    void setPid(Int_t pid)
    { pid_=pid; }

    const TLorentzVector& momentum() const
    { return p4_; }

    void setMomentum(double px, double py, double pz, double e)
    { p4_.SetPxPyPzE(px,py,pz,e); }
  
    const std::pair<Int_t, Int_t>& mothers() const
    { return mothers_; }

    void setMothers(Int_t a, Int_t b)
    { mothers_=std::make_pair(a,b); }

    virtual void Print(Option_t *option="") const
    {
    }

    ClassDef(GenParticleFormat,1)

  };


  class JetFormat : public GenJetFormat
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  private:
    Bool_t  bTag_;
    UInt_t  NTracks_;
    UInt_t  NCalo_;
    Float_t EHoverEE_;
    Char_t  MCID_;

    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------

  public :

    JetFormat()
    { reset(); }

    JetFormat(double px, double py, double pz, double e):GenJetFormat(px,py,pz,e)
    { }

    virtual ~JetFormat()
    { }

    virtual void reset()
    {
      bTag_ = false;
      NTracks_ = 0;
      NCalo_ = 0;
      EHoverEE_ = 0.;
      MCID_ = 0;
    }

    const Bool_t& btag() const
    { return bTag_; }

    void setBtag(Bool_t bTag)
    { bTag_=bTag; }

    const UInt_t& ntracks() const
    { return NTracks_; }

    void setNtracks(UInt_t Ntracks)
    { NTracks_=Ntracks; }

    const UInt_t& ncalo() const
    { return NCalo_; }

    void setNcalo(UInt_t NCalo)
    { NCalo_=NCalo; }

    const Float_t& ehoveree() const
    { return EHoverEE_; }

    void setEhoveree(Float_t EHoverEE)
    { EHoverEE_=EHoverEE; }

    virtual void Print(Option_t *option="") const
    {}

    void setMCID(Char_t value)
    { MCID_=value; }

    const Char_t& mcID() const
    { return MCID_; }

    ClassDef(JetFormat,1)

  };


  class METFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  protected:
    TVector2 p_;


    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------

  public:

    METFormat()
    { reset(); } 

    virtual ~METFormat()
    { } 

    double met() const 
    { return p_.Mod(); }

    double px() const 
    { return p_.Px(); }

    double py() const 
    { return p_.Py(); }

    void set(double px, double py)
    { p_.Set(px,py); }

    double phi() const 
    { return TVector2::Phi_mpi_pi(p_.Phi()); }

    virtual void reset()
    { p_.Set(0.,0.); }

    virtual void Print(Option_t *option="") const
    {
    }

    ClassDef(METFormat,1)

  };
  

  class MuonFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  private:
    TLorentzVector p4_;
    Bool_t isolated_;
    Bool_t charge_;

  public:
    MCLeptonFormat mc;

    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------

  public :

    MuonFormat()
    {
      reset();
    }

    MuonFormat(double px, double py, double pz, double e)
    {
      reset();
      p4_.SetPxPyPzE(px,py,pz,e);
    }

    virtual ~MuonFormat()
    { }

    void reset()
    {
      p4_.SetPxPyPzE(0.,0.,0.,0.);
      isolated_ = false;
      charge_ = false;
      mc.reset();
    }

    const TLorentzVector& momentum() const
    { return p4_; }

    void setMomentum(double px, double py, double pz, double e)
    { p4_.SetPxPyPzE(px,py,pz,e); }

    const Bool_t& isolated() const
    { return isolated_; }

    void setIsolated(Bool_t isolated)
    { isolated_=isolated; }

    float charge() const
    { if (charge_) return +1.;
      else return -1.; }

    void setCharge(float value)
    { if (value>0) charge_=true;
      else charge_=false; }

    virtual void Print(Option_t *option="") const
    {
    }

    ClassDef(MuonFormat,1)

  };



  class EventFormat : public TObject
  {

    // -------------------------------------------------------------
    //                        data members
    // -------------------------------------------------------------

  public:
    std::vector<GenParticleFormat> mcparticles;
    std::vector<JetFormat>         jets;
    std::vector<GenJetFormat>      mcjets;
    std::vector<MuonFormat>        muons;
    std::vector<ElectronFormat>    electrons;
    METFormat                      met;
    METFormat                      mcmet;


    // -------------------------------------------------------------
    //                        method members
    // -------------------------------------------------------------

  public :
    EventFormat()
    { }

    virtual ~EventFormat()
    {}

    void reset()
    {
      mcparticles.clear(); jets.clear(); muons.clear();
      electrons.clear(); met.reset(); mcmet.reset(); mcjets.clear();
    }

    virtual void Print(Option_t *option="") const
    {
    }

    ClassDef(EventFormat,1)

  };

}

#endif

