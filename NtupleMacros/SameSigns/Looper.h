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
  CUT_LT_PT,
  CUT_LL_PT,
  CUT_LT_GOOD,
  CUT_LL_GOOD,
  CUT_LT_CALOISO,
  CUT_LL_CALOISO,
  CUT_SAME_SIGN,
  CUT_OPP_SIGN,
  CUT_PASS2_MET,
  CUT_PASS4_MET,
  CUT_PASS2_METCORR,
  CUT_PASS4_METCORR,
  CUT_PASS_MET_10,
  CUT_PASS_MET_1,
  CUT_PASS_METCORR_10,
  CUT_PASS_METCORR_1,
  CUT_PASS_SUMET_10,
  CUT_PASS_SUMET_1,
  CUT_PASS_ZVETO,
  CUT_IN_Z_WINDOW,
  CUT_PASS_ADDZVETO,
  CUT_Z_TRACK_VETO_HYP,
  CUT_Z_TRACK_VETO_TRACKS,
  CUT_Z_TRACK_VETO_HYP_OR_TRACKS,
  CUT_PASS_MUON_B_VETO,	
  CUT_MUON_TAGGED,
  CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT,	
  CUT_MUON_TAGGED_WITHOUT_PTCUT,
  CUT_ELFAKE_FAKEABLE_OBJECT,
  CUT_ELFAKE_NUMERATOR,
  CUT_ELFAKE_NOT_NUMERATOR,
  CUT_MORE_THAN_TWO_TRACKS,
  CUT_LJETS,
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
// define useful cut combinations here
const static cuts_t samesigns_baseline_cuts = (
					 (CUT_BIT(CUT_LJETS)) |
					 (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS)) |
					 (CUT_BIT(CUT_LT_PT)		) | 
					 (CUT_BIT(CUT_LL_PT)		) |  
					 (CUT_BIT(CUT_LT_GOOD)	) | 
					 (CUT_BIT(CUT_LL_GOOD)	) | 
					 (CUT_BIT(CUT_LT_CALOISO)	) |  
					 (CUT_BIT(CUT_LL_CALOISO)	) | 
                                         //                                         (CUT_BIT(CUT_OPP_SIGN)	) | // TEMP out - need to write a dedicated looper for SJ cuts
					 (CUT_BIT(CUT_PASS4_METCORR)	) |  
					 (CUT_BIT(CUT_PASS2_METCORR)	) 
					 );   


// const static cuts_t samesigns_susybaseline_cuts = (
// 					 (CUT_BIT(CUT_LJETS)) |
// 					 (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS)) |
//                                          // 					 (CUT_BIT(CUT_LT_PT)		) | 
//                                          // 					 (CUT_BIT(CUT_LL_PT)		) |  
// 					 (CUT_BIT(CUT_LT_SUSYGOOD)	) | 
// 					 (CUT_BIT(CUT_LL_SUSYGOOD)	) | 
// 					 (CUT_BIT(CUT_LT_SUSYCALOISO)	) |  
// 					 (CUT_BIT(CUT_LL_SUSYCALOISO)	) | 
//                                          //					 (CUT_BIT(CUT_OPP_SIGN)	) | //WARTNING SS is in TOO!! TEMP 090629, IBL
// 					 (CUT_BIT(CUT_SUSYMETCORR)	)  
// 					 );   
		  
const static cuts_t samesigns_met_10_cuts = samesigns_baseline_cuts | CUT_BIT(CUT_PASS_METCORR_10);

const static cuts_t samesigns_met_1_cuts = samesigns_baseline_cuts | CUT_BIT(CUT_PASS_METCORR_1);

const static cuts_t samesigns_sumet_10_cuts = samesigns_baseline_cuts | CUT_BIT(CUT_PASS_SUMET_10);

const static cuts_t samesigns_sumet_1_cuts = samesigns_baseline_cuts | CUT_BIT(CUT_PASS_SUMET_1);

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
  virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
  virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
  virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
  //----------------------------------------------------------------------
  // declare your histograms here:
  //----------------------------------------------------------------------

  // NMinus1Hists take care of N - 1 plots and splitting by hypothesis automatically
  NMinus1Hist *hnJet_;
  NMinus1Hist *hnCaloJet_;
  NMinus1Hist *hnTrackJet_;
  NMinus1Hist *hltPt_;
  NMinus1Hist *hllPt_;
  NMinus1Hist *hltEta_;
  NMinus1Hist *hllEta_;
  NMinus1Hist *hltCaloIso_;
  NMinus1Hist *hllCaloIso_;
  NMinus1Hist *helPt_;
  NMinus1Hist *hmuPt_;
  NMinus1Hist *helEta_;
  NMinus1Hist *hmuEta_;
  NMinus1Hist *helCaloIso_;
  NMinus1Hist *hmuCaloIso_;
  NMinus1Hist *hnTrack_;
  NMinus1Hist *hdphiLep_;
  NMinus1Hist *hdilMass_;
  NMinus1Hist *hdilPt_;
  NMinus1Hist *hmet_;
  NMinus1Hist *hmetSpec_;
  NMinus1Hist *hmetTrkCorr_;
  NMinus1Hist *hsumet_;
  NMinus1Hist *hthrust_;
  NMinus1Hist *hjptthrust_;
  NMinus1Hist *hmeff_;
  NMinus1Hist *hmeffcorr_;

  TDatabasePDG *pdg;

  Double_t thrust;
  TVector3 thrust_axis;
  Double_t thrust_axis_phi, thrust_axis_theta;
  Double_t jptthrust;
  TVector3 jptthrust_axis;
  Double_t jptthrust_axis_phi, jptthrust_axis_theta;

  bool isWW_type1;
  bool isWOth_type2;
  bool isOthOth_type3;

protected:
  // count the (weighted and unweighted) number of candidates passing our cuts
  double		cands_passing_[4];
  double		cands_passing_w2_[4];
  unsigned int	        cands_count_[4];
  double                sumEt_;

  Double_t              calc_thrust(Double_t * thrust_axis_phi, Double_t * thrust_axis_theta, TString EFO="tracks");
  Double_t              find_thrust_axis(Double_t * thrust_axis_phi, Double_t * thrust_axis_theta, Double_t *thrust, TString EFO="tracks");
  bool                  idIsCharm(int id);
  bool                  idIsBeauty(int id);

};
#endif
