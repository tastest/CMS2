// -*- C++ -*-

#ifndef DSGLOOPER_H
#define DSGLOOPER_H

#include <math.h>
#include <stdio.h>
#include "TH2.h"
#include "LooperBase.h"
#include "DileptonHist.h"
#include "NMinus1Hist.h"
#include "../CORE/CMS2.h"

// enums for cuts
enum {
  DSG_LT_PT,
  DSG_LL_PT,
  DSG_SAME_SIGN,
  DSG_OPP_SIGN,
  DSG_PASS2_MET,
  DSG_PASS4_MET,
  DSG_PASS_MET_10,
  DSG_PASS_MET_1,
  DSG_PASS_SUMET_10,
  DSG_PASS_SUMET_1,
  DSG_PASS2_METCORR,
  DSG_PASS4_METCORR,
  DSG_LT_GOOD,
  DSG_LL_GOOD,
  DSG_EL_GOOD,
  DSG_EL_GOOD_NO_D0,
  DSG_LT_TIGHT_DPHIIN,
  DSG_LL_TIGHT_DPHIIN,
  DSG_MU_GOOD,
  DSG_ONE_SUPERTIGHT,
  DSG_TWO_SUPERTIGHT,
  DSG_LT_ISO,
  DSG_LL_ISO,
  DSG_ONE_ISO,
  DSG_TWO_ISO,
  DSG_EL_ISO,
  DSG_MU_ISO,
  DSG_LT_CALOISO,
  DSG_LL_CALOISO,
  DSG_ONE_CALOISO,
  DSG_TWO_CALOISO,
  DSG_EL_CALOISO,
  DSG_PASS_ZVETO,
  DSG_IN_Z_WINDOW,
  DSG_PASS_ADDZVETO,
  DSG_PASS_JETVETO_CALO,
  DSG_PASS_JETVETO_TRACKJETS,
  DSG_PASS_JETVETO_CALOTRACKJETS_COMBO,
  DSG_PASS_MUON_B_VETO,	
  DSG_MUON_TAGGED,
  DSG_PASS_MUON_B_VETO_WITHOUT_PTCUT,	
  DSG_MUON_TAGGED_WITHOUT_PTCUT,
  DSG_PASS_EXTRALEPTON_VETO,
  DSG_EL_BARREL,
  DSG_ELFAKE_FAKEABLE_OBJECT,
  DSG_ELFAKE_NUMERATOR,
  DSG_ELFAKE_NOT_NUMERATOR,
  DSG_MORE_THAN_TWO_TRACKS,
  DSG_Z_TRACK_VETO_HYP,
  DSG_Z_TRACK_VETO_TRACKS,
  DSG_Z_TRACK_VETO_HYP_OR_TRACKS
};

inline cuts_t CUT_BIT (int i)
{
  return 1ll << i;
}

// define useful cut combinations here
const static cuts_t dsg_baseline_cuts = 
		(CUT_BIT(DSG_MORE_THAN_TWO_TRACKS) |
		 // 		(CUT_BIT(DSG_PASS_ZVETO)	) | 
		 // 		(CUT_BIT(DSG_PASS_ADDZVETO)	) | 
		 // 		(CUT_BIT(DSG_PASS_MUON_B_VETO_WITHOUT_PTCUT) | 
		 //             (CUT_BIT(DSG_PASS_MUON_B_VETO)) |
// 		 (CUT_BIT(DSG_Z_TRACK_VETO_HYP_OR_TRACKS)) |
// 		 (CUT_BIT(DSG_Z_TRACK_VETO_HYP)) |
// 		 (CUT_BIT(DSG_Z_TRACK_VETO_TRACKS)) |
		 (CUT_BIT(DSG_LT_PT)		) | 
		 (CUT_BIT(DSG_LL_PT)		) | 
		 (CUT_BIT(DSG_OPP_SIGN)		) | 
		 (CUT_BIT(DSG_PASS4_MET)		) |  
		 (CUT_BIT(DSG_PASS2_MET)		) |  
		 (CUT_BIT(DSG_LT_GOOD)		) | 
		 (CUT_BIT(DSG_LL_GOOD)		) | 
		 (CUT_BIT(DSG_LT_CALOISO)	) |  
		 (CUT_BIT(DSG_LL_CALOISO)	) );   
		  
const static cuts_t dsg_met_10_cuts = dsg_baseline_cuts | CUT_BIT(DSG_PASS_MET_10);

const static cuts_t dsg_met_1_cuts = dsg_baseline_cuts | CUT_BIT(DSG_PASS_MET_1);

const static cuts_t dsg_sumet_10_cuts = dsg_baseline_cuts | CUT_BIT(DSG_PASS_SUMET_10);

const static cuts_t dsg_sumet_1_cuts = dsg_baseline_cuts | CUT_BIT(DSG_PASS_SUMET_1);

const static cuts_t dsg_fakerate_cuts = dsg_baseline_cuts & 
		~(CUT_BIT(DSG_LT_GOOD) | CUT_BIT(DSG_LL_GOOD) |
		  CUT_BIT(DSG_LT_ISO) | CUT_BIT(DSG_LL_ISO) |
		  CUT_BIT(DSG_LT_CALOISO) | CUT_BIT(DSG_LL_CALOISO));

class DSGLooperBase : public LooperBase {
  // DSG looper base class.  
  //
  // - provides unified cut bit pattern (but actually cutting is up to
  //   the derived classes)
  //
  // - provides function to fill all the histograms
public:
  DSGLooperBase (Sample, cuts_t cuts, const char *logfilename = 0);
  virtual ~DSGLooperBase () { }
  virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing[i]; }
  virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count[i]; }
  virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2[i]); }

protected:
  virtual cuts_t	DilepSelect 	(int idx);
  virtual void	FillHists 	(int idx);
  virtual double	Weight		(int idx);
  virtual void	End		();

protected:
  cuts_t		cuts;
  NMinus1Hist 	hnJet,
    hnCaloJet,
    hnTrackJet,
    hcaloJetPt,
    htrackJetPt,
    hminLepPt,
    hmaxLepPt,
    hltPt,
    hllPt,
    helPt,
    hmuPt,
    helEta,
    hmuEta,
    hdphiLep,
    hdilMass,
    hdilPt,
    hmet,
    hmetSpec,
    hmetTrkCorr,
    hptJet1,
    hptJet2,
    hptJet3,
    hptJet4,
    hetaJet1,
    hetaJet2,
    hetaJet3,
    hetaJet4,
    hnumTightLep,
    heleRelIso,
    heleRelIsoTrk,
    hmuRelIso,
    hminRelIso,
    hminRelIso_withCalo,
    htagMuPt,
    htagMuRelIso;
  // mc matches for fake studies
  NMinus1Hist 	hmuPdgId,
    hmuMoPdgId,
    helPdgId,
    helMoPdgId;
  // electron id variables
  NMinus1Hist	helEop,
    held0,
    helfbrem,
    helHE,
    helsee,
    helsppEB,
    helsppEE,
    heldphiin,
    heldetain,
    helEseedopin;
  // for conversion killing
  NMinus1Hist	helConvDeltaPhi_ss,
    helConvDeltaPhi_os;
  TH2F	held0vsRelIso, heldphiinvsRelIso,
    held0vsRelIsoMCgamma, heldphiinvsRelIsoMCgamma;
  TH2F	htrkCalodRvsPtSum;
  TH2F	hCaloEtaPt;

protected:
  cuts_t		cuts_passed; 
  FILE		*logfile;
  double		cands_passing[4];
  double		cands_passing_w2[4];
  unsigned int	cands_count[4];
};

// - stack plots (including N-1)
// - tables
class DSGResultsLooper : public DSGLooperBase {
public:
  DSGResultsLooper (Sample, 
		    cuts_t cuts = dsg_baseline_cuts, 
		    const char *fname = 0);
protected:
  virtual void	Dilep 		(int idx);
};

// - background estimate for W+jets from fake rates (only electrons for now)
class DSGFakeRateLooper : public DSGLooperBase {
public:
  DSGFakeRateLooper (Sample s, cuts_t cuts = dsg_fakerate_cuts, 
		     const char *fname = 0);
  virtual double	CandsPassingSystHi (enum DileptonHypType i) const { return cands_passing_syst_hi[i]; }
  virtual double	CandsPassingSystLo (enum DileptonHypType i) const { return cands_passing_syst_lo[i]; }
  virtual double	FakeSyst (enum DileptonHypType i) const;
protected:
  virtual void	FillHists 	(int idx);
  virtual cuts_t	DilepSelect	(int idx);
  virtual void	Dilep 		(int idx);
  virtual double	Weight		(int idx);
  virtual double	Weight		(int idx, int n_sig_syst);

protected:
  double		cands_passing_syst_hi[4];
  double		cands_passing_syst_lo[4];
  TH2F		fake_syst;
};

#endif
