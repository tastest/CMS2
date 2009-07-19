// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include <map>

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
  CUT_ALLLEP_PT20,
  CUT_FIRSTLEP_PT20,
  CUT_SECONDLEP_PT20,
  CUT_THIRDLEP_PT20,
  CUT_HIGHEST_PT_LEP_PT20,
  CUT_SECOND_HIGHEST_PT_LEP_PT20,
  CUT_THIRD_HIGHEST_PT_LEP_PT20,
  CUT_ALLLEP_GOOD,
  CUT_FIRSTLEP_GOOD,
  CUT_SECONDLEP_GOOD,
  CUT_THIRDLEP_GOOD,
  CUT_ALLLEP_ISO,
  CUT_FIRSTLEP_ISO,
  CUT_SECONDLEP_ISO,
  CUT_THIRDLEP_ISO,
  CUT_HIGHEST_PT_LEP_ISO,
  CUT_SECOND_HIGHEST_PT_LEP_ISO,
  CUT_THIRD_HIGHEST_PT_LEP_ISO,
  CUT_PRIM_Z,
  CUT_TCMET_15,
  CUT_PTCMET_15,
  CUT_NO_SEC_Z,
  CUT_ADD_ELECTRONS_VETO_CUT,
  CUT_ADD_MUONS_VETO_CUT,
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
  (CUT_BIT(CUT_HIGHEST_PT_LEP_PT20)) |
  (CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_PT20)) |
  (CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_PT20)) |
  (CUT_BIT(CUT_ALLLEP_GOOD)) |
  (CUT_BIT(CUT_HIGHEST_PT_LEP_ISO)) |
  (CUT_BIT(CUT_SECOND_HIGHEST_PT_LEP_ISO)) |
  (CUT_BIT(CUT_THIRD_HIGHEST_PT_LEP_ISO)) |
  (CUT_BIT(CUT_PRIM_Z)) |
  (CUT_BIT(CUT_TCMET_15)) |
  (CUT_BIT(CUT_NO_SEC_Z)) |
  (CUT_BIT(CUT_ADD_ELECTRONS_VETO_CUT)) |
  (CUT_BIT(CUT_ADD_MUONS_VETO_CUT)) ;

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
  virtual double	CandsPassing (enum TrileptonHypType i) const { return cands_passing_[i]; }
  virtual int       	CandsCount (enum TrileptonHypType i) const { return cands_count_[i]; }
  virtual double	RMS (enum TrileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  // NMinus1Hists take care of N - 1 plots and splitting by hypothesis automatically
  TrilepNMinus1Hist	*htcmet_;
  TrilepNMinus1Hist	*hptcmet_;
  TrilepNMinus1Hist	*h_highest_lep_pt_;
  TrilepNMinus1Hist	*h_second_highest_lep_pt_;
  TrilepNMinus1Hist	*h_third_highest_lep_pt_;
  TrilepNMinus1Hist	*h_highest_lep_iso_;
  TrilepNMinus1Hist	*h_second_highest_lep_iso_;
  TrilepNMinus1Hist	*h_third_highest_lep_iso_;
  TrilepNMinus1Hist	*h_highest_iso_lep_iso_;
  TrilepNMinus1Hist	*h_second_highest_iso_lep_iso_;
  TrilepNMinus1Hist	*h_third_highest_iso_lep_iso_;

//   TrilepNMinus1Hist	*h_additional_muon_pt_;
//   TrilepNMinus1Hist	*h_additional_muon_iso_;
//   TrilepNMinus1Hist	*h_additional_electron_pt_;
//   TrilepNMinus1Hist	*h_additional_electron_iso_;

  TrilepNMinus1Hist	*h_counter_electrons_;
  TrilepNMinus1Hist	*h_counter_muons_;

  TrilepNMinus1Hist	*h_njets_;

  TrilepNMinus1Hist	*h_DeltaPhiMETNearestLepton_;

  TrilepNMinus1Hist	*h_primZMass_;
  TrilepNMinus1Hist	*h_addZMass_;

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double		cands_passing_[21];
  double		cands_passing_w2_[21];
  unsigned int	cands_count_[21];
  
  std::multimap<float,float,std::greater<float> > trileptonPt_;
  std::multimap<float,float,std::greater<float> > trileptonIso_;

  unsigned int addElectronsCounter_;
  unsigned int addMuonsCounter_;

  double primZMass_;
  double addZMass_;
};
#endif
