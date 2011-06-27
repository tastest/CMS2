// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "TDatabasePDG.h"

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
  CUT_MU_NUMERATOR,
  CUT_MU_DENOMINATOR,
  CUT_QCD_BIN_UPPER_PTHAT,
  CUT_EVEN,
  CUT_ODD,
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
 
// define useful cut combinations here

const static cuts_t observation_cuts = 
  (CUT_BIT(CUT_MU_NUMERATOR));   

const static cuts_t observation_cuts_qcd_bins = 
  observation_cuts | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));   

const static cuts_t prediction_cuts = (CUT_BIT(CUT_MU_DENOMINATOR));

const static cuts_t prediction_cuts_qcd_bins = prediction_cuts | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));

const static cuts_t observation_cuts_even = 
  (CUT_BIT(CUT_MU_NUMERATOR)) | (CUT_BIT(CUT_EVEN));   

const static cuts_t observation_cuts_qcd_bins_even = 
  observation_cuts_even | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));   

const static cuts_t prediction_cuts_even = 
  (CUT_BIT(CUT_MU_DENOMINATOR)) | (CUT_BIT(CUT_EVEN));

const static cuts_t prediction_cuts_qcd_bins_even = 
  prediction_cuts_even | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));
						 
const static cuts_t observation_cuts_odd = 
  (CUT_BIT(CUT_MU_NUMERATOR)) | (CUT_BIT(CUT_ODD));   

const static cuts_t observation_cuts_qcd_bins_odd = 
  observation_cuts_odd | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));   

const static cuts_t prediction_cuts_odd = 
  (CUT_BIT(CUT_MU_DENOMINATOR)) | (CUT_BIT(CUT_ODD));

const static cuts_t prediction_cuts_qcd_bins_odd = 
  prediction_cuts_odd | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));
						 
//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

// WW looper: make all the standard histograms and yields
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

  // we define an analysis-specific EventSelect(), DilepSelect(),
  // TrilepSelect() and QuadlepSelect() that check which cuts the
  // event, dilepton/trilepton/quadlepton candidate passes
  virtual cuts_t	EventSelect	();
  virtual cuts_t	DilepSelect 	(int idx);
  virtual cuts_t	TrilepSelect 	(int idx);
  virtual cuts_t	QuadlepSelect 	(int idx);
  // we define an analysis-specific set of FillEventHistos(),
  // FillDilepHistos(), FillTrilepHistos() and FillQuadlepHistos()
  // that fill our histograms.  
  // 
  // the framework calls our FillEventHistos() function for every event
  virtual void	FillEventHistos ();
  // the framework calls our FillDilepHistos() function for every
  // dilepton candidate; the argument is the index of the candidate
  // in the dilepton block
  virtual void	FillDilepHistos (int idx);
  // the framework calls our FillTrilepHistos() function for every
  // trilepton candidate; the argument is the index of the candidate
  // in the trilepton block
  virtual void	FillTrilepHistos (int idx);
  // the framework calls our FillQuadlepHistos() function for every
  // quadlepton candidate; the argument is the index of the candidate
  // in the quadlepton block
  virtual void	FillQuadlepHistos (int idx);
  // at the end of the loop, we get a callback to do things like
  // printing a status message
  virtual void	End		();

  // weight for dilepton candidate
  virtual double 	Weight (int);

public:
  // these functions are called by the table-printing code
  virtual double	CandsPassing () const { return cands_passing_; }
  virtual int       	CandsCount () const { return cands_count_; }
  virtual double	RMS () const { return sqrt(cands_passing_w2_); }

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  TH1F	*hmuPt_;
  TH1F	*hmuEta_;
  
protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double	cands_passing_;
  double	cands_passing_w2_;
  unsigned int	cands_count_;
  unsigned int  events_;

};

#endif
