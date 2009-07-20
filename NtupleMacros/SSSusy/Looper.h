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
// #define PRETTY_PRINT(args ...) args
// PRETTY_PRINT (
     enum {
       CUT_TTBAR_TYPE_WW,
       CUT_TTBAR_TYPE_WO,
       CUT_TTBAR_TYPE_OO,
       CUT_TRUE_MU_FROM_W,
       CUT_TRUE_EL_FROM_W,
       CUT_NOT_TRUE_GAMMA_FROM_MUON,
       CUT_MAX_PT,
       CUT_MIN_PT,
       CUT_MU_PT,
       CUT_SAME_SIGN,
       CUT_OPP_SIGN,
       // 	  CUT_PASS2_MET,
       // 	  CUT_PASS4_MET,
       // 	  CUT_PASS2_TCMET,
       // 	  CUT_PASS4_TCMET,
       CUT_TCMET,
       CUT_LT_GOOD,
       CUT_LL_GOOD,
       CUT_LT_ISO,
       CUT_LL_ISO,
       CUT_MU_GOOD,
       CUT_MU_ISO,
       CUT_PASS_ZVETO,
       CUT_IN_Z_WINDOW,
       CUT_PASS_ADDZVETO,
       CUT_ELFAKE_FAKEABLE_OBJECT, // these are
       CUT_ELFAKE_NUMERATOR,	      // here for
       CUT_ELFAKE_NOT_NUMERATOR,   // applying it to emu only first
       //	  CUT_MORE_THAN_TWO_TRACKS,
       CUT_PASS_TRIGGER,
     };
//      );

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

// this is the current baseline set of cuts
const static cuts_t baseline_cuts = 
  (CUT_BIT(CUT_MIN_PT)		) | 
  (CUT_BIT(CUT_MAX_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)        ) | 
  //  (CUT_BIT(CUT_TCMET)		) |   //temp out 090714
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_ISO)	        ) |  
  (CUT_BIT(CUT_LL_ISO)        	) |
  (CUT_BIT(CUT_MU_PT)		) |  // currently always require the muon to have pt>20
  ( CUT_BIT(CUT_NOT_TRUE_GAMMA_FROM_MUON) ) |
  (CUT_BIT(CUT_TRUE_MU_FROM_W)	) //|
  //  (CUT_BIT(CUT_TTBAR_TYPE_WO)   )   //temp out 090714
// |  
//     (CUT_BIT(CUT_PASS_ZVETO)	) | 
//     (CUT_BIT(CUT_PASS_TRIGGER))
  ;   

// denominator object cuts for the fake rate prediction 
const static cuts_t fakerate_denominator_cuts = (baseline_cuts & 
						 ~(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
						   CUT_BIT(CUT_LT_ISO)  | CUT_BIT(CUT_LL_ISO)  | 
                                                   CUT_BIT(CUT_TTBAR_TYPE_WO)  )) |
  CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);

// numerator object cuts for the fake rate prediction 
const static cuts_t fakerate_numerator_cuts = 
  //  fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NUMERATOR) | CUT_BIT(CUT_TTBAR_TYPE_WO); // remove - being suspicious of the WO cut 090715_17:39
  fakerate_denominator_cuts | 
  (CUT_BIT(CUT_ELFAKE_NUMERATOR)) |
  (CUT_BIT(CUT_MU_PT)		) | 
  (CUT_BIT(CUT_MU_GOOD)		) | 
  (CUT_BIT(CUT_MU_ISO)	        ) ;

// denominator and not numerator (this is the yield that should be
// multiplied by FR / (1 - FR))
const static cuts_t fakerate_denominator_not_numerator_cuts = 
  fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);

static const cuts_t fakerate_ss_numerator_cuts = 
  (fakerate_numerator_cuts & ~CUT_BIT(CUT_OPP_SIGN)) | CUT_BIT(CUT_SAME_SIGN);

static const cuts_t fakerate_ss_denominator_not_numerator_cuts = 
  (fakerate_denominator_not_numerator_cuts & ~CUT_BIT(CUT_OPP_SIGN) )
  | CUT_BIT(CUT_SAME_SIGN);

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

     // count candidates passing various cuts
     virtual void	CountCuts (cuts_t, double);
public:
     // these functions are called by the table-printing code
     virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
     virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------

     NMinus1Hist 	*hnJet;
     NMinus1Hist	*hnCaloJet;
     NMinus1Hist	*hnTrackJet;
     NMinus1Hist	*hnJPTJet;
     NMinus1Hist	*hcaloJetPt;
     NMinus1Hist	*htrackJetPt;
     NMinus1Hist	*hJPTJetPt;
     NMinus1Hist	*hminLepPt;
     NMinus1Hist	*hmaxLepPt;
     NMinus1Hist	*hltPt;
     NMinus1Hist	*hllPt;
     NMinus1Hist	*helPt;
     NMinus1Hist	*hmuPt;
     NMinus1Hist	*helEta;
     NMinus1Hist	*hmuEta;
     NMinus1Hist	*hdphiLep;
     NMinus1Hist	*hdilMass;
     NMinus1Hist	*hdilPt;
     NMinus1Hist	*hmet;
     NMinus1Hist	*hmetSpec;
     NMinus1Hist	*hmetTrkCorr;
     NMinus1Hist	*hptJet1;
     NMinus1Hist	*hptJet2;
     NMinus1Hist	*hptJet3;
     NMinus1Hist	*hptJet4;
     NMinus1Hist	*hetaJet1;
     NMinus1Hist	*hetaJet2;
     NMinus1Hist	*hetaJet3;
     NMinus1Hist	*hetaJet4;
     NMinus1Hist	*hnumTightLep;
     NMinus1Hist	*heleRelIso;
     NMinus1Hist	*heleRelIsoTrk;
     NMinus1Hist	*hmuRelIso;
     NMinus1Hist	*hminRelIso;
     NMinus1Hist	*hminRelIso_withCalo;
     NMinus1Hist	*htagMuPt;
     NMinus1Hist	*htagMuRelIso;
     // mc matches for fake studies
     NMinus1Hist 	*hmuPdgId;
     NMinus1Hist	*hmuMoPdgId;
     NMinus1Hist	*helPdgId;
     NMinus1Hist	*helMoPdgId;
     // electron id variables
     NMinus1Hist	*helEop;
     NMinus1Hist	*held0;
     NMinus1Hist	*helfbrem;
     NMinus1Hist	*helHE;
     NMinus1Hist	*helsee;
     NMinus1Hist	*helsppEB;
     NMinus1Hist	*helsppEE;
     NMinus1Hist	*heldphiin;
     NMinus1Hist	*heldetain;
     NMinus1Hist	*helEseedopin;
     // for conversion killing
     NMinus1Hist	*helConvDeltaPhi_ss;
     NMinus1Hist	*helConvDeltaPhi_os;
     // for random other things
     TH2F		*held0vsRelIso; 
     TH2F		*heldphiinvsRelIso;
     TH2F		*held0vsRelIsoMCgamma; 
     TH2F		*heldphiinvsRelIsoMCgamma;
     TH2F		*htrkCalodRvsPtSum;
     TH2F		*hCaloEtaPt;

protected:
     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing_[4];
     double		cands_passing_w2_[4];
     unsigned int	cands_count_[4];
     double		count_cuts_[64];
     double		count_correlation_[64][64];
};

// background estimate for W+jets from fake rates (only electrons in
// emu for now)
class FakeRateLooper : public Looper {
public:
     FakeRateLooper (Sample s, cuts_t cuts, const char *fname = 0);
     virtual double	CandsPassingSystHi (enum DileptonHypType i) const { return cands_passing_syst_hi[i]; }
     virtual double	CandsPassingSystLo (enum DileptonHypType i) const { return cands_passing_syst_lo[i]; }
//      virtual double	FakeSyst (enum DileptonHypType i) const;
protected:
     virtual void	BookHistos 	();
     virtual cuts_t	DilepSelect 	(int idx);
     virtual void	FillDilepHistos (int idx);
//      virtual double	Weight		(int idx);
//      virtual double	Weight		(int idx, int n_sig_syst);

protected:
     double		cands_passing_syst_hi[4];
     double		cands_passing_syst_lo[4];
     TH2F		*fake_syst;
};

class EventCountingLooper : public Looper {
public:
     EventCountingLooper (Sample s, cuts_t cuts, const char *fname = 0);
protected:
     virtual void	End		();
     virtual void	BeforeDilepHistos	();
     virtual void	AfterDilepHistos	();
protected:
     double		events_passing_[4];
     double		events_passing_w2_[4];
     double		cands_passing_prev_[4];
};     

#endif
