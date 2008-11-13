// -*- C++ -*-

#ifndef LOOPERBASE_H
#define LOOPERBASE_H

#include "Rtypes.h"
#include "Sample.h"

typedef UInt_t 		uint32;
typedef ULong64_t 	uint64;

typedef uint64	cuts_t;

// Base class for loopers.  Takes care of the looping, calls the
// virtual functions for each event / dilep candidate / trilep
// candidate.  

class LooperBase {
public:
     LooperBase (Sample);
     virtual ~LooperBase ();
     virtual uint64 	Loop ();
     const std::string 	&SampleName () const { return sample.name; }
     const bool 	&SampleSM () const { return sample.sm; }
     
protected:
     virtual cuts_t	EventSelect 	()		{ return 1; 	}
     virtual void	Event 		()		{ return; 	}
     virtual cuts_t	DilepSelect 	(int idx)	{ return 1; 	}
     virtual void	Dilep 		(int idx)	{ return; 	}
     virtual cuts_t	TrilepSelect 	(int idx) 	{ return 1; 	}
     virtual void	Trilep 		(int idx)	{ return; 	}
     virtual cuts_t	QuadlepSelect 	(int idx) 	{ return 1; 	}
     virtual void	Quadlep 	(int idx)	{ return; 	}

     virtual void	Begin ();
     virtual void	End () 			{ }

     Sample		sample;
     double		hypos_total_weight[4];
     uint64		hypos_total_n[4];
};

#endif
