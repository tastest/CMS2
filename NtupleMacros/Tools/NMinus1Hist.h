// -*- C++ -*-
#ifndef NMINUS1HIST_H
#define NMINUS1HIST_H

#include <string>
#include <vector>
#include "TH1.h"
#include "DileptonHist.h"
#include "LooperBase.h"

// extension of the dilepton hypo hist:
//
// - takes a set of cuts (the "N cuts", default event selection)
//
// - plus a cut mask (the "1" in N - 1)
//
// - for fancy people, there is also an option to specify a whole
//   vector of cut masks (useful for plotting MET without MET or MET
//   special cut, for example)
// 
// - accumulates a histo of events that pass the N - 1 cuts (as well
//   as the events that pass all N)
//

class Sample;
class HistParm;
class NMinus1Hist {
public:
     NMinus1Hist (const Sample &, const std::string &name, 
		  int bins, double min, double max,
		  cuts_t cuts, cuts_t cut_mask);
     NMinus1Hist (const Sample &, const std::string &name, 
		  int bins, double min, double max,
		  cuts_t cuts, std::vector<cuts_t> cut_mask, 
		  std::vector<std::string> cut_names = std::vector<std::string>());
     virtual ~NMinus1Hist () { }
     void		Fill (cuts_t cuts_passed, enum DileptonHypType, double, 
			      double w = 1);
     DileptonHist	&N () 		{ return h_N; }
     DileptonHist	&NMinus1 ()	{ return *h_NMinus1[0]; }
     operator const DileptonHist & () const { return h_N; }

protected:
     cuts_t		cuts;
     std::vector<cuts_t> cut_mask;
     std::vector<cuts_t> mask;
     DileptonHist	h_N;
     std::vector<DileptonHist *> h_NMinus1;
};

#endif
