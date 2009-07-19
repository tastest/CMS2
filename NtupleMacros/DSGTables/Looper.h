// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "DSGTable.h"

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
     CUT_LT_PT,
     CUT_LL_PT,
     CUT_LT_GOOD,
     CUT_LL_GOOD,
     CUT_LT_CALOISO,
     CUT_LL_CALOISO,
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
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_CALOISO)	) |  
     (CUT_BIT(CUT_LL_CALOISO)	);   

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

     // filter out this event.  If FilterEvent returns true, no
     // further processing is done on this event
     virtual bool	FilterEvent();
     virtual cuts_t	DilepSelect (int idx);
     virtual void	FillDilepHistos (int idx);
 
protected:
     // functions to classify a dilepton candidate
     int		Zcat (int) const;
     int		METcat (int) const;
     int                SumJetcat (int) const;
     int		Jetcat (int) const;
     int		Bucket (int) const;
     
public:
     // these functions are called by the table-printing code
//      virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
//      virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
//      virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

public:
     DSGTable		dsgTable;
};
#endif
