// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "Math/LorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

// List of all cuts that can be applied.  The cuts are handled as a
// bitfield; these labels define which bit corresponds to which cut.
// The cut are tested and the corresponding bits are set for each
//  - event in EventSelect()
//  - dilepton candidate in DilepSelect()
//  - trilepton candidate in TrilepSelect()
//  - quadlepton candidate in QuadlepSelect().
enum {
  CUT_MORE_THAN_TWO_TRACKS,
  CUT_NL,
  CUT_E,
  CUT_M,
  CUT_EE,
  CUT_EM,
  CUT_MM,
  CUT_ML,
  CUT_OS,
  CUT_ZMASS,
  CUT_ANTI_ZMASS,
  CUT_MET,
  CUT_ANTI_MET,
  CUT_MT,
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
 
//single electron cuts
const static cuts_t baseline_single_electron_cuts = 
//   (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_E) ) |
  (CUT_BIT(CUT_ANTI_ZMASS) ) |
  (CUT_BIT(CUT_MET) ) |
  (CUT_BIT(CUT_MT) );

//single muon cuts
const static cuts_t baseline_single_muon_cuts = 
//   (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_M) ) |
  (CUT_BIT(CUT_ANTI_ZMASS) ) |
  (CUT_BIT(CUT_MET) ) |
  (CUT_BIT(CUT_MT) );

//dielectron cuts
const static cuts_t baseline_dielectron_cuts = 
//   (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_EE) ) |
  (CUT_BIT(CUT_ZMASS) ) |
  (CUT_BIT(CUT_ANTI_MET) ) |
  (CUT_BIT(CUT_OS) );

//dimuon cuts
const static cuts_t baseline_dimuon_cuts = 
//   (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_MM) ) |
  (CUT_BIT(CUT_ZMASS) ) |
  (CUT_BIT(CUT_ANTI_MET) ) |
  (CUT_BIT(CUT_OS) );

//EM cuts
const static cuts_t baseline_emu_cuts = 
//   (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_EM) ) |
  (CUT_BIT(CUT_ZMASS) ) |
  (CUT_BIT(CUT_ANTI_MET) ) |
  (CUT_BIT(CUT_OS) );

//single electron test cuts
const static cuts_t test_single_electron_cuts = 
  (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_E) );

//single muon test cuts
const static cuts_t test_single_muon_cuts = 
  (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_M) );

//dielectron test cuts
const static cuts_t test_dielectron_cuts = 
  (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_OS) ) |
  (CUT_BIT(CUT_EE) );

//dimuon test cuts
const static cuts_t test_dimuon_cuts = 
  (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_OS) ) |
  (CUT_BIT(CUT_MM) );

//emu test cuts
const static cuts_t test_emu_cuts = 
  (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS) ) |
  (CUT_BIT(CUT_EM) );

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

  // variables
  static const int mu_shift = 0x8000000;
  enum Type {NL, E, M, EE, EM, MM, ML};

  float wmt_;

  LorentzVector lep1;
  LorentzVector lep2;
  LorentzVector boson;
  LorentzVector wp4_;
  LorentzVector genp4_;

public:
  // these functions are called by the table-printing code
  virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
  virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
  virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  // a simple TH1/TH2, times four to split by hypothesis type:
  TH1F		*helPt_[4];
  TH1F		*hmuPt_[4];
  TH2F		*hCaloEtaPt_[4];

  // NMinus1Hists take care of N - 1 plots and splitting by hypothesis automatically
  NMinus1Hist	*hmt_;
  NMinus1Hist	*hmet_;
     NMinus1Hist	*hnjets_;
     NMinus1Hist	*hnjptjets_;
  NMinus1Hist	*hdilMass_;
  NMinus1Hist	*hLepMetMass_;
  NMinus1Hist   *hGenLepEta_;
  NMinus1Hist   *hGenLepPt_;

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double	cands_passing_[4];
  double	cands_passing_w2_[4];
  unsigned int	cands_count_[4];
};
#endif
