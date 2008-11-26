// -*- C++ -*-

#ifndef LOOPERBASE_H
#define LOOPERBASE_H

#include <stdio.h>
#include <math.h>
#include "Rtypes.h"
#include "LooperTypes.h"
#include "Sample.h"
#include "DileptonHist.h"
#include "NMinus1Hist.h"
#include "TH1.h"
#include "TH2.h"

// Base class for loopers.  Takes care of the looping, calls the
// virtual functions for each event / dilep candidate / trilep
// candidate.  

class LooperBase {
public:
     // constructor; tell the looper what sample to loop on (see
     // Tools/Sample.h), what cuts candidates need to pass, and a file
     // name for dumping log messages
     LooperBase (Sample, cuts_t, const char *logfilename = 0);

protected:
     // book your histograms here
     virtual void	BookHistos	()		{ }

     // filter out this event.  If FilterEvent returns true, no
     // further processing is done on this event
     virtual bool	FilterEvent	() 		{ return false; }
     
     // The *Select() functions return a bitmask that says which cuts
     // a candidate passes.  "Candidate" means event, dilepton
     // candidate, trilepton candidate or quadlepton candidate, as
     // appropriate for your analysis
     virtual cuts_t	EventSelect 	()		{ return 0; 	}
     virtual cuts_t	DilepSelect 	(int idx)	{ return 0; 	}
     virtual cuts_t	TrilepSelect 	(int idx) 	{ return 0; 	}
     virtual cuts_t	QuadlepSelect 	(int idx) 	{ return 0; 	}

     // These functions are called for every event, dilepton/
     // trilepton/ quadlepton candidate.  Fill your histograms in
     // these functions.  It is up to you to check that the candidate
     // passes your cuts.
     virtual void	FillEventHistos 	()		{ return; 	}
     virtual void	FillDilepHistos 	(int idx)	{ return; 	}
     virtual void	FillTrilepHistos 	(int idx)	{ return; 	}
     virtual void	FillQuadlepHistos 	(int idx)	{ return; 	}

     // Begin() is called once before looping stars, End() is called
     // once after looping ends.  If you need extra initializations
     // etc, put them here.
     virtual void	Begin ()		{ }
     virtual void	End () 			{ }

     // get the event weight
     // - by default, this returns evt_scale1fb * kFactor
     //   (the kFactor is set when the samples are defined, see Sample.h)
     // - sometimes, specific applications might require a fancier weight
     virtual double	Weight		(int idx);

     //------------------------------------------------------------
     // everything below here is implementation details ...
     //------------------------------------------------------------

public:
     virtual ~LooperBase ();
     virtual uint64 	Loop ();
     const std::string 	&SampleName () const { return sample.name; }
     const bool 	&SampleSM () const { return sample.sm; }
     
protected:
     Sample		sample;
     cuts_t		cuts;
     FILE		*logfile;
     double		hypos_total_weight[4];
     uint64		hypos_total_n[4];
     int                duplicates_total_n_; // for duplicate filter
     double             duplicates_total_weight_; // for duplicate filter
};

inline cuts_t CUT_BIT (int i)
{
     return 1ll << i;
}

#endif
