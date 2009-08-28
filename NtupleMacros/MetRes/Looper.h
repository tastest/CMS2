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
     CUT_LT_PT,
     CUT_LL_PT,
	 CUT_HYP_PT,
     //CUT_SAME_SIGN,
     CUT_OPP_SIGN,
     //CUT_PASS2_MET,
     //CUT_PASS4_MET,
     //CUT_PASS2_METCORR,
     //CUT_PASS4_METCORR,
	 CUT_TCMET,
     CUT_LT_GOOD,
     CUT_LL_GOOD,
     //CUT_EL_GOOD,
     //CUT_EL_GOOD_NO_D0,
     //CUT_LT_TIGHT_DPHIIN,
     //CUT_LL_TIGHT_DPHIIN,
     //CUT_MU_GOOD,
     //CUT_ONE_SUPERTIGHT,
     //CUT_TWO_SUPERTIGHT,
     CUT_LT_ISO,
     CUT_LL_ISO,
     //CUT_ONE_ISO,
     //CUT_TWO_ISO,
     //CUT_EL_ISO,
     //CUT_MU_ISO,
     //CUT_LT_CALOISO,
     //CUT_LL_CALOISO,
     //CUT_ONE_CALOISO,
     //CUT_TWO_CALOISO,
     //CUT_EL_CALOISO,
     //CUT_PASS_ZVETO,
     CUT_IN_Z_WINDOW,
     CUT_IN_OUTER_Z,
	 CUT_Z_SF,
	 CUT_SUSY_TRIGGER,
	 CUT_SUMJETPT200,
	 CUT_NJETS2
};


const static cuts_t baseline_susy = 
 CUT_BIT(CUT_HYP_PT) |
 CUT_BIT(CUT_OPP_SIGN) |
 CUT_BIT(CUT_TCMET) |
 CUT_BIT(CUT_LT_GOOD) |
 CUT_BIT(CUT_LL_GOOD) |
 CUT_BIT(CUT_LT_ISO) |
 CUT_BIT(CUT_LL_ISO) |
//CUT_BIT(CUT_IN_Z_WINDOW) |
 CUT_BIT(CUT_Z_SF) |
 CUT_BIT(CUT_SUSY_TRIGGER) |
 CUT_BIT(CUT_SUMJETPT200) |
 CUT_BIT(CUT_NJETS2);

const static cuts_t baseline_metres = 
//CUT_BIT(CUT_IN_Z_WINDOW) |
 CUT_BIT(CUT_HYP_PT) |
 CUT_BIT(CUT_OPP_SIGN) |
 CUT_BIT(CUT_SUSY_TRIGGER) |
 CUT_BIT(CUT_LT_GOOD) |
 CUT_BIT(CUT_LL_GOOD) |
 CUT_BIT(CUT_LT_ISO) |
 CUT_BIT(CUT_LL_ISO);

const static cuts_t inz_metres = baseline_metres | CUT_BIT(CUT_IN_Z_WINDOW);
const static cuts_t outerz_metres = baseline_metres | CUT_BIT(CUT_IN_OUTER_Z);

     //CUT_PASS_ADDZVETO,
     //CUT_PASS_JETVETO_CALO,
     //CUT_PASS_JETVETO_TRACKJETS,
     //CUT_PASS_JETVETO_CALOTRACKJETS_COMBO,
     //CUT_PASS_EXTRALEPTON_VETO,
     //CUT_EL_BARREL,
     //CUT_ELFAKE_FAKEABLE_OBJECT,
     //CUT_ELFAKE_NUMERATOR,
     //CUT_ELFAKE_NOT_NUMERATOR,

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

  virtual void NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max);
  virtual void NewProf(TProfile*& h, char* name, char* title, int bins, double min, double max);

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

  double sumjetpt;
  int sjpbin; //from 0 to nsjpbins-1
  double counttt;
  bool printevt;
  
  // TH1/TH2s, times four to split by hypothesis type:
  TH1F* htcmet[4];
  TH1F* htcmetinz[4];
  TH1F* htcmetouz[4];
  //components of tcmet
  TH1F* htcmetx[4];
  TH1F* htcmetinzx[4];
  TH1F* htcmetouzx[4];
  TH1F* htcmety[4];
  TH1F* htcmetinzy[4];
  TH1F* htcmetouzy[4];
  //both components together (not sum, two entries per hyp)
  TH1F* htcmetxy[4];
  TH1F* htcmetouzxy[4];
  TH1F* htcmetinzxy[4];

  //number sumjetpt bins: 0,30-100, 100-200, 200+ (note, nothing btwn 0 and 30 b'c min jet pt is 30)
#define nsjpbins 4 
  TH1F* htcmetxy_sjp[4][nsjpbins];
  TH1F* htcmetouzxy_sjp[4][nsjpbins];
  TH1F* htcmetinzxy_sjp[4][nsjpbins];

  TH1F* hdphi_sjp[4][nsjpbins];
  TH1F* hdphiouz_sjp[4][nsjpbins];
  TH1F* hdphiinz_sjp[4][nsjpbins];

  TH1F* hmetvmll_sjp[4][nsjpbins];
  TH1F* hmetvmllouz_sjp[4][nsjpbins];
  TH1F* hmetvmllinz_sjp[4][nsjpbins];
  TH1F* hmetvmllden_sjp[4][nsjpbins];
  TH1F* hmetvmllouzden_sjp[4][nsjpbins];
  TH1F* hmetvmllinzden_sjp[4][nsjpbins];

  TProfile* pmetvmll_sjp[4][nsjpbins];
  TProfile* pmetvmllouz_sjp[4][nsjpbins];
  TProfile* pmetvmllinz_sjp[4][nsjpbins];

  //2d
  TH2F* htcmetxvy[4];

  TH1F* hsumjetpt[4];
  TH1F* hsumjetptinZ[4];
  TH1F* hsumjetptouZ[4];

  TH1F* helPt_[4];
  TH1F* hmuPt_[4];

  // NMinus1Hists take care of N - 1 plots and splitting by hypothesis automatically
  //NMinus1Hist	*hltPt_;
  //NMinus1Hist	*hllPt_;
  //NMinus1Hist	*hdilMass_;

protected:
     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing_[4];
     double		cands_passing_w2_[4];
     unsigned int	cands_count_[4];
};
#endif
