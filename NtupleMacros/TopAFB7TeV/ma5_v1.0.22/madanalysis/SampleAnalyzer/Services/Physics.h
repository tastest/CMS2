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


#ifndef PHYSICS_SERVICE_h
#define PHYSICS_SERVICE_h

// STL headers
#include <set>
#include <string>
#include <algorithm>

// SampleAnalyzer headers
#include "Services/MCconfig.h"
#include "Services/RECconfig.h"
#include "DataFormat/MCEventFormat.h"
#include "DataFormat/RecEventFormat.h"


#define PHYSICS PhysicsService::getInstance()

enum OrderingObservable{Eordering, Pordering, PTordering, 
                        ETordering, PXordering, PYordering,
                        PZordering, ETAordering};



struct PointerComparison
{
  template<typename T>
  static bool ESortPredicate(T* part1, 
                             T* part2)
  { return part1->e() > part2->e(); }

  template<typename T>
  static bool ETSortPredicate(T* part1, 
                              T* part2)
  { return part1->et() > part2->et(); }

  template<typename T>
  static bool PSortPredicate(T* part1, 
                             T* part2)
  { return part1->p() > part2->p(); }

  template<typename T>
  static bool PTSortPredicate(T* part1, 
                              T* part2)
  { return part1->pt() > part2->pt(); }

  template<typename T>
  static bool ETASortPredicate(T* part1, 
                               T* part2)
  { return part1->eta() > part2->eta(); }

  template<typename T>
  static bool PXSortPredicate(T* part1, 
                              T* part2)
  { return part1->px() > part2->px(); }

  template<typename T>
  static bool PYSortPredicate(T* part1, 
                       T* part2)
  { return part1->py() > part2->py(); }

  template<typename T>
  static bool PZSortPredicate(T* part1, 
                              T* part2)
  { return part1->pz() > part2->pz(); }

};


class PhysicsService
{

  // -------------------------------------------------------------
  //                       data members
  // -------------------------------------------------------------
 protected:

  MCconfig mcConfig_;
  RECconfig recConfig_;
  Int_t finalstate_;
  Int_t initialstate_;
  static PhysicsService* service_;

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public:

  /// GetInstance
  static PhysicsService* getInstance()
  {
    if (service_==0) service_ = new PhysicsService;
    return service_;
  }

  /// Kill
  static void kill()
  {
    if (service_!=0) delete service_;
    service_=0;
  }


  /// Is Initial State
  Bool_t IsInitialState(const MCParticleFormat& part) const
  {
    return (part.statuscode()==initialstate_);
  }

  /// Is Final State
  Bool_t IsFinalState(const MCParticleFormat& part) const
  {
    return (part.statuscode()==finalstate_);
  }

  /// Is Final State
  Bool_t IsInterState(const MCParticleFormat& part) const
  {
    return (part.statuscode()!=finalstate_ && part.statuscode()!=initialstate_);
  }

  /// Is Initial State
  Bool_t IsInitialState(const MCParticleFormat* part) const
  {
    return (part->statuscode()==initialstate_);
  }

  /// Is Final State
  Bool_t IsFinalState(const MCParticleFormat* part) const
  {
    return (part->statuscode()==finalstate_);
  }

  /// Is Final State
  Bool_t IsInterState(const MCParticleFormat* part) const
  {
    return (part->statuscode()!=finalstate_ && part->statuscode()!=initialstate_);
  }

  /// Set Initial State
  void SetInitialState(const MCEventFormat* myEvent)
  {
    if (myEvent==0) return; 
    if (myEvent->particles().empty()) return;
    initialstate_=myEvent->particles()[myEvent->particles().size()-1].statuscode(); 
  }

  /// Set Final State
  void SetFinalState(const MCEventFormat* myEvent)
  {
   
    if (myEvent==0) return; 
    if (myEvent->particles().empty()) return;
    finalstate_=myEvent->particles()[myEvent->particles().size()-1].statuscode(); 
  }

  /// Get MCconfig
  const MCconfig& mcConfig() const  
  { return mcConfig_; }
  MCconfig& mcConfig()
  { return mcConfig_; }

  /// Get RECconfig
  const RECconfig& recConfig() const  
  { return recConfig_; }
  RECconfig& recConfig()
  { return recConfig_; }

  /// Is hadronic ?
  inline bool IsHadronic(const RecParticleFormat* part) const
  {
    if (dynamic_cast<const RecJetFormat*>(part)==0) return false;
    else return true;
  }

  /// Is invisible ?
  inline bool IsInvisible(const RecParticleFormat* part) const
  {
    return false;
  }

  inline bool IsHadronic(const RecParticleFormat& part) const
  {
    return IsHadronic(&part);
  }

  /// Is invisible ?
  inline bool IsInvisible(const RecParticleFormat& part) const
  {
    return IsInvisible(&part);
  }

  /// Is hadronic ?
  inline bool IsHadronic(const MCParticleFormat& part) const
  {
    std::set<Int_t>::iterator found = mcConfig_.hadronic_ids_.find(part.pdgid());
    if (found==mcConfig_.hadronic_ids_.end()) return false; else return true;
  }

  /// Is invisible ?
  inline bool IsInvisible(const MCParticleFormat& part) const
  {
    std::set<Int_t>::iterator found = mcConfig_.invisible_ids_.find(part.pdgid());
    if (found==mcConfig_.invisible_ids_.end()) return false; else return true;
  }

  /// Is hadronic ?
  inline bool IsHadronic(const MCParticleFormat* part) const
  {
    if (part==0) return false;
    return IsHadronic(*part);
  }

  /// Is invisible ?
  inline bool IsInvisible(const MCParticleFormat* part) const
  {
    if (part==0) return false;
    return IsInvisible(*part);
  }

  /// Compute the total transverse energy
  inline double EventTET(const MCEventFormat* event) const
  { 
    double energy=0;
    for (unsigned int i=0;i<event->particles().size();i++)
    {
      if (event->particles()[i].statuscode()==finalstate_)
          energy+=event->particles()[i].pt();
    }
    return energy;
  }

  /// Compute the missing transverse energy
  inline double EventMET(const MCEventFormat* event) const
  {
    TLorentzVector q(0.,0.,0.,0.);
    for (unsigned int i=0;i<event->particles().size();i++)
    {
      if (event->particles()[i].statuscode()==finalstate_)
        if (!IsInvisible(event->particles()[i]))
          q+=event->particles()[i].momentum();
    }
    return q.Perp();
  }

  /// Compute the total hadronic transverse energy
  inline double EventTHT(const MCEventFormat* event) const
  {
    double energy=0;
    for (unsigned int i=0;i<event->particles().size();i++)
    {
      if (event->particles()[i].statuscode()==finalstate_)
        if (IsHadronic(event->particles()[i]))
          energy+=event->particles()[i].pt();
    }
    return energy;
  }

  /// Compute the missing hadronic transverse energy
  inline double EventMHT(const MCEventFormat* event) const
  {
    TLorentzVector q(0.,0.,0.,0.);
    for (unsigned int i=0;i<event->particles().size();i++)
    {
      if (event->particles()[i].statuscode()==finalstate_)
        if (!IsInvisible(event->particles()[i]) &&
            IsHadronic(event->particles()[i]))
          q+=event->particles()[i].momentum();
    }
    return q.Perp();
  }

  /// Compute the total transverse energy
  inline double EventTET(const RecEventFormat* event) const
  { 
    double energy=0;

    for (unsigned int i=0;i<event->jets().size();i++)
      energy+=event->jets()[i].et();
    for (unsigned int i=0;i<event->electrons().size();i++)
      energy+=event->electrons()[i].et();
    for (unsigned int i=0;i<event->muons().size();i++)
      energy+=event->muons()[i].et();
    for (unsigned int i=0;i<event->taus().size();i++)
      energy+=event->taus()[i].et();
  
    return energy;
  }

  /// Compute the missing transverse energy
  inline double EventMET(const RecEventFormat* event) const
  {
    return event->MET().magnitude();
  }

  /// Compute the total hadronic transverse energy
  inline double EventTHT(const RecEventFormat* event) const
  {
    double energy=0;
    for (unsigned int i=0;i<event->jets().size();i++)
    {
      energy+=event->jets()[i].et();
    }
    return energy;
  }

  /// Compute the missing hadronic transverse energy
  inline double EventMHT(const RecEventFormat* event) const
  {
    TLorentzVector q(0.,0.,0.,0.);
    for (unsigned int i=0;i<event->jets().size();i++)
    {
      q+=event->jets()[i].momentum();
    }
    return q.Et(); 
  }

  /// rank filter
  void rankFilter(std::vector<const MCParticleFormat*>& parts,Short_t rank,
                    OrderingObservable obs=PTordering) const
  {
    // rank equal to zero
    if (rank==0)
    {
      WARNING << "Rank equal to 0 is not possible. Allowed values are 1,2,3,... and -1,-2,-3,..." << std::endl;
      parts.clear();
      return;
    }

    // Number of particle is not correct
    if ((static_cast<Int_t>(parts.size()) - std::abs(static_cast<Int_t>(rank)))<0 ) 
    { 
      parts.clear();
      return;
    }

    // Sorting particle collection
    sort(parts,obs);

    // Keeping the only particle
    std::vector<const MCParticleFormat*> tmp(1);  
    if (rank>0) tmp[0]=parts[rank-1];
    else tmp[0]=parts[parts.size()+rank];

    // Saving tmp
    parts = tmp;
  }

  /// rank filter
  void rankFilter(std::vector<const RecParticleFormat*>& parts,Short_t rank,
                    OrderingObservable obs=PTordering) const
  {
    // rank equal to zero
    if (rank==0)
    {
      WARNING << "Rank equal to 0 is not possible. Allowed values are 1,2,3,... and -1,-2,-3,..." << std::endl;
      parts.clear();
      return;
    }

    // Number of particle is not correct
    if ((static_cast<Int_t>(parts.size()) - std::abs(static_cast<Int_t>(rank)))<0 ) 
    { 
      parts.clear();
      return;
    }

    // Sorting particle collection
    sort(parts,obs);

    // Keeping the only particle
    std::vector<const RecParticleFormat*> tmp(1);  
    if (rank>0) tmp[0]=parts[rank-1];
    else tmp[0]=parts[parts.size()+rank];

    // Saving tmp
    parts = tmp;
  }


  /// sort particle
  void sort(std::vector<const MCParticleFormat*>& parts,
            OrderingObservable obs=PTordering) const
  {
    if (obs==PTordering) 
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PTSortPredicate<const MCParticleFormat>);
    else if (obs==ETordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ETSortPredicate<const MCParticleFormat>);
    else if (obs==Eordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ESortPredicate<const MCParticleFormat>);
    else if (obs==ETAordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ETASortPredicate<const MCParticleFormat>);
    else if (obs==PXordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PXSortPredicate<const MCParticleFormat>);
    else if (obs==PYordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PYSortPredicate<const MCParticleFormat>);
    else if (obs==PZordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PZSortPredicate<const MCParticleFormat>);
  }

  /// sort particle
  void sort(std::vector<const RecParticleFormat*>& parts,
            OrderingObservable obs=PTordering) const
  {
    if (obs==PTordering) 
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PTSortPredicate<const RecParticleFormat>);
    else if (obs==ETordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ETSortPredicate<const RecParticleFormat>);
    else if (obs==Eordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ESortPredicate<const RecParticleFormat>);
    else if (obs==ETAordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::ETASortPredicate<const RecParticleFormat>);
    else if (obs==PXordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PXSortPredicate<const RecParticleFormat>);
    else if (obs==PYordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PYSortPredicate<const RecParticleFormat>);
    else if (obs==PZordering)
        std::sort(parts.begin(),parts.end(),
                  PointerComparison::PZSortPredicate<const RecParticleFormat>);
  }


  void ToRestFrame(MCParticleFormat& part, const MCParticleFormat* boost) const
  {
    if (boost==0) return;
    ToRestFrame(part,*boost);
  }

  void ToRestFrame(MCParticleFormat& part, const MCParticleFormat& boost) const
  {
    TVector3 b = -1. * boost.momentum().BoostVector();
    part.momentum().Boost(b);
  }

  /// Compute srqt(S)
  inline double SqrtS(const MCEventFormat* event) const
  {
    TLorentzVector q(0.,0.,0.,0.);
    for (UInt_t i=0;i<event->particles().size();i++)
    {
      if ( event->particles()[i].statuscode() == initialstate_ )
        q += event->particles()[i].momentum();
    }
    return sqrt(q.Mag2());
  }

  /// Muon isolation
  Bool_t IsIsolatedMuon(const RecLeptonFormat* muon,
                        const RecEventFormat* event) const
  {
    // Safety
    if (muon==0 || event==0) return false;

    // Method : DeltaR
    if (recConfig_.deltaRalgo_)
    {
      // Loop over jets
      for (unsigned int i=0;i<event->jets().size();i++)
      {
        if ( muon->dr(event->jets()[i]) < recConfig_.deltaR_ ) return false;
      }
      return true;
    }

    // Method : SumPT
    else
    {
      return ( muon->sumPT_isol() < recConfig_.sumPT_ && 
               muon->ET_PT_isol() < recConfig_.ET_PT_  );
    }

    return true;
  } 

  /// Muon isolation
  Bool_t IsIsolatedMuon(const RecLeptonFormat& part,
                        const RecEventFormat* event) const
  {
    return IsIsolatedMuon(&part,event);
  }

 private:

  /// Constructor
  PhysicsService()  
  { initialstate_=-1; finalstate_=1; }

  /// Destructor
  ~PhysicsService()
  {}
};

#endif
