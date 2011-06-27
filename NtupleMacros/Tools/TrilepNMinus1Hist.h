// -*- C++ -*-
#ifndef TRILEPNMINUS1HIST_H
#define TRILEPNMINUS1HIST_H

#include <string>
#include <vector>
#include "TH1.h"
#include "TrileptonHist.h"
#include "LooperBase.h"
#include "LooperTypes.h"

// extension of the trilepton hypo hist:
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
class TrilepNMinus1Hist {
public:
     TrilepNMinus1Hist (const Sample &, const std::string &name, 
		  int bins, double min, double max,
		  cuts_t cuts, cuts_t cut_mask);
     TrilepNMinus1Hist (const Sample &, const std::string &name, 
		  int bins, double min, double max,
		  cuts_t cuts, std::vector<cuts_t> cut_mask, 
		  std::vector<std::string> cut_names = std::vector<std::string>());
     virtual ~TrilepNMinus1Hist () { }
     void		Fill (cuts_t cuts_passed, enum TrileptonHypType, double, 
			      double w = 1);
     TrileptonHist	&N () 		{ return h_N; }
     TrileptonHist	&NMinus1 ()	{ return *h_NMinus1[0]; }
     operator const TrileptonHist & () const { return h_N; }

protected:
     cuts_t		cuts;
     std::vector<cuts_t> cut_mask;
     std::vector<cuts_t> mask;
     TrileptonHist	h_N;
     std::vector<TrileptonHist *> h_NMinus1;
};

#endif
