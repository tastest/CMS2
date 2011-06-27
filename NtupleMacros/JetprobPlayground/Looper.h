// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "signedImpact.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

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
  CUT_LT_PT,
  CUT_LL_PT,
  CUT_SAME_SIGN,
  CUT_OPP_SIGN,
  CUT_PASS2_TCMET,
  CUT_PASS4_TCMET,
  CUT_LT_GOOD,
  CUT_LL_GOOD,
  CUT_LT_CALOISO,
  CUT_LL_CALOISO,
  CUT_PASS_ZVETO,
  CUT_PASS_JETVETO_JPT20,
  CUT_PASS_JETVETO_JPT25,
  CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT,	
  CUT_MORE_THAN_TWO_TRACKS,
  CUT_PASS_TRIGGER,
  CUT_PASS_JETVETO_SIP,
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

// this is the current baseline set of cuts
const static cuts_t baseline_cuts = 
  (CUT_BIT(CUT_LT_PT)		) | 
  (CUT_BIT(CUT_LL_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)		) | 
  (CUT_BIT(CUT_PASS4_TCMET)		) |  
  (CUT_BIT(CUT_PASS2_TCMET)		) |  
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_CALOISO)	) |  
  (CUT_BIT(CUT_LL_CALOISO)	) |  
  (CUT_BIT(CUT_PASS_ZVETO)	) | 
  (CUT_BIT(CUT_PASS_JETVETO_JPT20)	) |  
  (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT))	|
  (CUT_BIT(CUT_PASS_TRIGGER)) ;

const static cuts_t jetprob_study_cuts = 
  (CUT_BIT(CUT_LT_PT)		) | 
  (CUT_BIT(CUT_LL_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)		) | 
  (CUT_BIT(CUT_PASS4_TCMET)		) |  
  (CUT_BIT(CUT_PASS2_TCMET)		) |  
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_CALOISO)	) |  
  (CUT_BIT(CUT_LL_CALOISO)	) |  
  (CUT_BIT(CUT_PASS_ZVETO)	) | 
  (CUT_BIT(CUT_PASS_JETVETO_JPT20)	) |  
  //(CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT))	|
  (CUT_BIT(CUT_PASS_TRIGGER)) ;
  //  (CUT_BIT(CUT_PASS_JETVETO_SIP) );

const static cuts_t new_baseline_cuts = 
  (CUT_BIT(CUT_LT_PT)		) | 
  (CUT_BIT(CUT_LL_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)		) | 
  (CUT_BIT(CUT_PASS4_TCMET)		) |  
  (CUT_BIT(CUT_PASS2_TCMET)		) |  
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_CALOISO)	) |  
  (CUT_BIT(CUT_LL_CALOISO)	) |  
  (CUT_BIT(CUT_PASS_ZVETO)	) | 
  (CUT_BIT(CUT_PASS_JETVETO_JPT20)	) |  
  (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT))	|
  (CUT_BIT(CUT_PASS_TRIGGER)) |
  (CUT_BIT(CUT_PASS_JETVETO_SIP) );

const static cuts_t mu_trkprob_cuts = 
  (CUT_BIT(CUT_LT_PT)		) | 
  (CUT_BIT(CUT_LL_PT)		) | 
  (CUT_BIT(CUT_OPP_SIGN)		) | 
  (CUT_BIT(CUT_LT_GOOD)		) | 
  (CUT_BIT(CUT_LL_GOOD)		) | 
  (CUT_BIT(CUT_LT_CALOISO)	) |  
  (CUT_BIT(CUT_LL_CALOISO)	) ;

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

  // for signed impact parameter study
  NMinus1Hist *hnJPTJet;
  NMinus1Hist *hJPTJetPt;
  NMinus1Hist *hntracks;
  NMinus1Hist *htrd0; 
  NMinus1Hist *htrd0sig; 
  NMinus1Hist *htrd0ByPt[4];      
  NMinus1Hist *htrd0sigByPt[4];   
  NMinus1Hist *htrd0ByNtrks[4];   
  NMinus1Hist *htrd0sigByNtrks[4];
  NMinus1Hist *htrd0Strange; // mother is a Ks or Lambda
  NMinus1Hist *htrd0sigStrange; // mother is a Ks or Lambda
  NMinus1Hist *hjetpt; 
  NMinus1Hist *hjetntr;  
  NMinus1Hist *htrkprob;
  NMinus1Hist *hjetprob;
  NMinus1Hist *hlogjetprob;
  NMinus1Hist *hjetprobByPt[4];
  NMinus1Hist *hjetprobByNtrks[4];
  NMinus1Hist *hminjetprob;
  NMinus1Hist *hmaxptjetprob;
  NMinus1Hist *hmaxptjetprobByPt[4];
  NMinus1Hist *htrkd0[3];
  NMinus1Hist *htrknchi2[3];
  NMinus1Hist *htrkvalidhits[3];
  NMinus1Hist *hmaxptjetprobNtrks5[5];
  NMinus1Hist *hjetprobNtrks2Pt10;
  NMinus1Hist *hmaxntrksPt10;
  TH2D        *hntrksvsjptpt;

  // use muons to study track probability
  NMinus1Hist *hmutrkprob;
  NMinus1Hist *hmud0sig;

  // track jets: a track jet is a pair<LorentzVector, vector<track indices>>
  std::vector<std::pair<LorentzVector, std::vector<unsigned int> > > &TrackJets () { return trackjets_; }
  std::vector<double>	&TrackJetProbs () { return trackjet_jp_; }
  virtual void	MakeTrackJets (int i_hyp);
  virtual void	MakeTrackJetProbs ();

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double		cands_passing_[4];
  double		cands_passing_w2_[4];
  unsigned int	cands_count_[4];
  double		count_cuts_[64];
  double		count_correlation_[64][64];

  std::vector<std::pair<LorentzVector, std::vector<unsigned int> > > trackjets_;
  std::vector<double>	trackjet_jp_;
};

// background estimate for W+jets from fake rates (only electrons in emu for now)
class FakeRateLooper : public Looper {
public:
  FakeRateLooper (Sample s, cuts_t cuts, const char *fname = 0);
  virtual double	CandsPassingSystHi (enum DileptonHypType i) const { return cands_passing_syst_hi[i]; }
  virtual double	CandsPassingSystLo (enum DileptonHypType i) const { return cands_passing_syst_lo[i]; }

protected:
  virtual void	BookHistos 	();
  virtual cuts_t	DilepSelect 	(int idx);
  virtual void	FillDilepHistos (int idx);

protected:
  double		cands_passing_syst_hi[4];
  double		cands_passing_syst_lo[4];
  TH2F		*fake_syst;
};

#endif
