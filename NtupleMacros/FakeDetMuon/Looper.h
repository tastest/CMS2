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
  CUT_PT_LEADING_JET,
  CUT_QCD_BIN_UPPER_PTHAT,
  CUT_NO_CUT,
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
const static cuts_t muo_fakes_cuts = 
  (CUT_BIT(CUT_NO_CUT));   

const static cuts_t muo_fakes_wo_trigger_jet_cuts = 
  (CUT_BIT(CUT_PT_LEADING_JET));   

const static cuts_t muo_fakes_cuts_using_bins = 
  muo_fakes_cuts | (CUT_BIT(CUT_QCD_BIN_UPPER_PTHAT));

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

public:
  // these functions are called by the table-printing code
  virtual double	EventsPassing () const { return events_passing_; }
  virtual int       	EventsCount () const { return events_count_; }
  virtual double	RMS () const { return sqrt(events_passing_w2_); }

  virtual double calculateFakeRateError(unsigned int nNum, unsigned int nDen);
  virtual double calculateWeightedFakeRateError(double nNum, double nNumErr2, double nDen, double nDenErr2);

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  TH1F *pt_num_mll_;
  TH1F *pt_num_mlt_;
  TH1F *pt_den_muo_;

  TH1F *eta_num_mll_;
  TH1F *eta_num_mlt_;
  TH1F *eta_den_muo_;

  TH2F *num_mll_;
  TH2F *num_mlt_;
  TH2F *den_muo_;

  TH1F *pt_num_wo_leading_mll_;
  TH1F *pt_num_wo_leading_mlt_;
  TH1F *pt_den_wo_leading_muo_;

  TH1F *eta_num_wo_leading_mll_;
  TH1F *eta_num_wo_leading_mlt_;
  TH1F *eta_den_wo_leading_muo_;

  TH2F *num_wo_leading_mll_;
  TH2F *num_wo_leading_mlt_;
  TH2F *den_wo_leading_muo_;

  TH1F *pt_num_wo_second_leading_mll_;
  TH1F *pt_num_wo_second_leading_mlt_;
  TH1F *pt_den_wo_second_leading_muo_;

  TH1F *eta_num_wo_second_leading_mll_;
  TH1F *eta_num_wo_second_leading_mlt_;
  TH1F *eta_den_wo_second_leading_muo_;

  TH2F *num_wo_second_leading_mll_;
  TH2F *num_wo_second_leading_mlt_;
  TH2F *den_wo_second_leading_muo_;
	
  TH1F *muo_n_;

  TH1F *njets_;

  TH1F *jetpt_;

  TH1F *jeteta_;

  TH1F *deltaR_;

  TH1F *tkIso_;
  TH1F *tkIso_uncut_;

  TH1F *EOverp_;
  TH1F *EOverp_uncut_;

  TH1F *HOverE_;
  TH1F *HOverE_uncut_;

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  unsigned int          events_;
  double                events_weighted_;
  double		events_passing_;
  double		events_passing_w2_;
  unsigned int		events_count_;
  unsigned int          nDenominator_;
  unsigned int          nDenominator_wo_leading_;
  unsigned int          nDenominator_wo_second_leading_;
  unsigned int          nNumerator_ll_;
  unsigned int          nNumerator_ll_wo_leading_;
  unsigned int          nNumerator_ll_wo_second_leading_;
  unsigned int          nNumerator_lt_;
  unsigned int          nNumerator_lt_wo_leading_;
  unsigned int          nNumerator_lt_wo_second_leading_;

  double                nDenominator_weighted_;
  double                nDenominator_wo_leading_weighted_;
  double                nDenominator_wo_second_leading_weighted_;
  double                nNumerator_ll_weighted_;
  double                nNumerator_ll_wo_leading_weighted_;
  double                nNumerator_ll_wo_second_leading_weighted_;
  double                nNumerator_lt_weighted_;
  double                nNumerator_lt_wo_leading_weighted_;
  double                nNumerator_lt_wo_second_leading_weighted_;

  double                nDenominator_weighted_w2_;
  double                nDenominator_wo_leading_weighted_w2_;
  double                nDenominator_wo_second_leading_weighted_w2_;
  double                nNumerator_ll_weighted_w2_;
  double                nNumerator_ll_wo_leading_weighted_w2_;
  double                nNumerator_ll_wo_second_leading_weighted_w2_;
  double                nNumerator_lt_weighted_w2_;
  double                nNumerator_lt_wo_leading_weighted_w2_;
  double                nNumerator_lt_wo_second_leading_weighted_w2_;

  TDatabasePDG *pdg;

};
#endif
