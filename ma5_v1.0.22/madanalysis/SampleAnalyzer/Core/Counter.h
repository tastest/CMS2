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


#ifndef COUNTER_h
#define COUNTER_h

// STL headers
#include <map>
#include <string>
#include <sstream>

// ROOT headers
#include <TH1F.h>

template <typename T> 
class Counter
{

  // -------------------------------------------------------------
  //                        data members
  // -------------------------------------------------------------
 private :

  //! Collection of observables
  std::map<T,ULong_t> stack_;

  // -------------------------------------------------------------
  //                       method members
  // -------------------------------------------------------------
 public :

  //! Constructor without argument
  Counter()
  { }

  //! Destructor
  ~Counter()
  { }

  typedef typename std::map<T,ULong_t>::iterator       iterator; 
  typedef typename std::map<T,ULong_t>::const_iterator const_iterator; 
  typedef typename std::map<T,ULong_t>::size_type      size_type;

  //! Adding an entry for a given observable
  void Add(const T& obs)
  {
    iterator it = stack_.find(obs);
    
    if (it==stack_.end()) stack_[obs]=1; // Obs not found
    else stack_[obs]+=1;                 // Obs found  
  }

  //! Returning the size of the stack
  size_type size() const
  { return stack_.size(); }

  //! Return a pointer to the first element of the stack
  const_iterator begin() const
  { return stack_.begin(); }

  //! Return a pointer to the first element of the stack
  const_iterator end() const
  { return stack_.end(); }

  //! Creating a ROOT histo
  void FillHisto(TH1F* histo)
  {

    histo -> SetBins(size(),0, size());
    
    unsigned int i=0;
    for (const_iterator it=begin();it!=end();it++)
    {
      histo->SetBinContent(i+1,it->second);
      std::string tmp;
      std::stringstream str;
      str << it->first;
      str >> tmp;
      histo->GetXaxis()->SetBinLabel(i+1,tmp.c_str());
      i++;
    }
  }

};

#endif
