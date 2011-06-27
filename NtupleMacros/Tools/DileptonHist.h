// -*- C++ -*-
#ifndef DILEPTONHIST_H
#define DILEPTONHIST_H

#include <string>
#include "TH1.h"
#include "DileptonHypType.h"
#include "HistParm.h"

// wrapper around TH1F to handle the dilepton hypothesis types
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
class DileptonHist {
public:
     typedef TH1F H_t; // just in case we need to switch to TH1D eventually
     
public:
     DileptonHist (const Sample &, const std::string &name, int bins, 
		   double min, double max);
     DileptonHist (const Sample &, const HistParm &);
     virtual ~DileptonHist ();
     H_t 	&Hist (enum DileptonHypType i) { return *histos[i]; }
     H_t	&operator [] (enum DileptonHypType i) { return *histos[i]; }
     int	Fill (enum DileptonHypType, double, double w = 1);
     unsigned int Entries (enum DileptonHypType i) const { return entries[i]; }
     double	Integral (enum DileptonHypType i) const { return integral[i]; }

protected:
     void 	build (const Sample &, const std::string &name, int bins, 
		       double min, double max);
     H_t	*histos[4]; // need to be pointers, since that's the root way
     unsigned int	entries[4];
     double		integral[4];
};

#endif
