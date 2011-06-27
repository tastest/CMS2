// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
     CUT_SAME_SIGN,
     CUT_OPP_SIGN,
     CUT_LT_GOOD,
     CUT_LL_GOOD,
};

//----------------------------------------------------------------------
// Cut combinations for selections.  These are examples that are used
// for various tables in the WW analysis.
// ----------------------------------------------------------------------

// Note: some quick reminders for bitwise operations.  
//
// - To turn the enum into a cut bit, use the CUT_BIT(x) function.
//   For example, to set the bit for CUT_LT_PT, use CUT_BIT(CUT_LT_PT).
//
// - Combine bits with the bitwise OR operator (|).  For example, to
//   turn on the bits for CUT_LT_PT and CUT_LL_PT, use
//   CUT_BIT(CUT_LT_PT) | CUT_BIT(CUT_LL_PT).  For another example,
//   making a new set of cuts that is the same as an existing set with
//   one cut added (the hypothetical CUT_EXTRA_CUT), use the
//   following: 
//   cuts_t baseline_cuts_with_extra_cut = baseline_cuts | CUT_BIT(CUT_EXTRA_CUT);
//
// - To turn off a bit (useful when defining a set of cuts that is the
//   same as another set of cuts with one cut removed), use the
//   following: 
//   cuts_t baseline_cuts_without_lt_pt = baseline_cuts & ~CUT_BIT(CUT_LT_PT);
 
// baseline cuts
const static cuts_t baseline_cuts = 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) ;   

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

// Looper for an analysis.  
//
// - switching between files, removing duplicates and other technical stuff is handled by LooperBase
// - analysis-specific stuff is defined here: 
//    * declaring histograms and booking them
//    * naming cuts (using the enum from above) and checking which cuts are passed
//    * filling histograms
//    * counting passing candidates
class Looper : public LooperBase {

public:
     // constructor; tell the looper what sample to loop on (see
     // Tools/Sample.h), what cuts candidates need to pass, and a file
     // name for dumping log messages
     Looper (Sample, cuts_t cuts, const char *logfilename = 0);
     virtual ~Looper () { }

protected:
     // this is where we book our histograms
     virtual void	BookHistos ();
     virtual cuts_t	DilepSelect 	(int idx);
     virtual void	FillDilepHistos (int idx);
     virtual void	End		();

protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------

     TH1F	*hdilMass_;
};
#endif
