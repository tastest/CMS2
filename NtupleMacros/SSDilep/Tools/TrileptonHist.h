// -*- C++ -*-
#ifndef TRILEPTONHIST_H
#define TRILEPTONHIST_H

#include <string>
#include "TH1.h"
#include "TrileptonHypType.h"
#include "HistParm.h"

// wrapper around TH1F to handle the trilepton hypothesis types
//
// - use like a TH1F, but supply the hypo type when you fill
//
// - supply the sample name in the constructor, since the histo name
//   is of the form samplename_variablename_hyponame
//
// - these histograms are added to the gDirectory, so the
//   saving-everything-to-file feature works

class Sample;
class HistParm;
class TrileptonHist {
public:
     typedef TH1F H_t; // just in case we need to switch to TH1D eventually
     
public:
     TrileptonHist (const Sample &, const std::string &name, int bins, 
		   double min, double max);
     TrileptonHist (const Sample &, const HistParm &);
     virtual ~TrileptonHist ();
     H_t 	&Hist (enum TrileptonHypType i) { return *histos[i]; }
     H_t	&operator [] (enum TrileptonHypType i) { return *histos[i]; }
     int	Fill (enum TrileptonHypType, double, double w = 1);
     unsigned int Entries (enum TrileptonHypType i) const { return entries[i]; }
     double	Integral (enum TrileptonHypType i) const { return integral[i]; }

protected:
     void 	build (const Sample &, const std::string &name, int bins, 
		       double min, double max);
     H_t	*histos[21]; // need to be pointers, since that's the root way
     unsigned int	entries[21];
     double		integral[21];
};

#endif
