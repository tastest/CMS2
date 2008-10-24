// -*- C++ -*-

#ifndef WWLOOPER_H
#define WWLOOPER_H

#include <math.h>
#include <stdio.h>
#include "TH2.h"
#include "LooperBase.h"
#include "DileptonHist.h"
#include "NMinus1Hist.h"
#include "CMS2.h"

// enums for cuts
enum {
     WW_LT_PT,
     WW_LL_PT,
     WW_SAME_SIGN,
     WW_OPP_SIGN,
     WW_PASS2_MET,
     WW_PASS4_MET,
     WW_LT_GOOD,
     WW_LL_GOOD,
     WW_EL_GOOD,
     WW_MU_GOOD,
     WW_ONE_SUPERTIGHT,
     WW_TWO_SUPERTIGHT,
     WW_LT_ISO,
     WW_LL_ISO,
     WW_ONE_ISO,
     WW_TWO_ISO,
     WW_EL_ISO,
     WW_MU_ISO,
     WW_LT_CALOISO,
     WW_LL_CALOISO,
     WW_ONE_CALOISO,
     WW_TWO_CALOISO,
     WW_EL_CALOISO,
     WW_PASS_ZVETO,
     WW_PASS_ADDZVETO,
     WW_PASS_JETVETO_CALO,
     WW_PASS_JETVETO_TRACKJETS,
     WW_PASS_MUON_B_VETO,	
     WW_MUON_TAGGED,
     WW_PASS_MUON_B_VETO_WITHOUT_PTCUT,	
     WW_MUON_TAGGED_WITHOUT_PTCUT,
     WW_PASS_EXTRALEPTON_VETO,
     WW_EL_BARREL,
     WW_ELFAKE_FAKEABLE_OBJECT,
     WW_ELFAKE_NUMERATOR,
     WW_ELFAKE_NOT_NUMERATOR,
};

inline cuts_t CUT_BIT (int i)
{
     return 1ll << i;
}

// define useful cut combinations here
const static cuts_t ww_baseline_cuts = 
     (CUT_BIT(WW_LT_PT)		) | 
     (CUT_BIT(WW_LL_PT)		) | 
     (CUT_BIT(WW_OPP_SIGN)		) | 
     (CUT_BIT(WW_PASS4_MET)		) |  
     (CUT_BIT(WW_PASS2_MET)		) |  
     (CUT_BIT(WW_LT_GOOD)		) | 
     (CUT_BIT(WW_LL_GOOD)		) | 
     (CUT_BIT(WW_LT_CALOISO)	) |  
     (CUT_BIT(WW_LL_CALOISO)	) |  
     (CUT_BIT(WW_PASS_ADDZVETO)	) | 
     (CUT_BIT(WW_PASS_JETVETO_CALO)	) |
     (CUT_BIT(WW_PASS_JETVETO_TRACKJETS)	) |  
     (CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT)	);   
/*      (CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT)	) | */
/*      (CUT_BIT(WW_PASS_EXTRALEPTON_VETO)	);    */

// these cuts are used to measure the mu tagging efficiency for top
const static cuts_t ww_baseline_mu_tageff_cuts = ww_baseline_cuts & ~((CUT_BIT(WW_PASS_JETVETO_CALO)	) | 
								      (CUT_BIT(WW_PASS_JETVETO_TRACKJETS)	) | 
								      (CUT_BIT(WW_PASS_EXTRALEPTON_VETO)	) | 
								      (CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT)	) | 
								      (CUT_BIT(WW_PASS_MUON_B_VETO)));

const static cuts_t ww_baseline_no_trackjets_cuts = ww_baseline_cuts & 
     ~(CUT_BIT(WW_PASS_JETVETO_TRACKJETS)); 

const static cuts_t ww_baseline_no_btags_cuts = ww_baseline_cuts & 
     ~(CUT_BIT(WW_PASS_MUON_B_VETO) | CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT)); 

const static cuts_t ww_baseline_no_caloiso_cuts = (ww_baseline_cuts & ~((CUT_BIT(WW_LT_CALOISO)	) | (CUT_BIT(WW_LL_CALOISO)	))) | 
     (CUT_BIT(WW_LT_ISO)		) | (CUT_BIT(WW_LL_ISO)		);

const static cuts_t ww_old_baseline_cuts = 
     (CUT_BIT(WW_LT_PT)		) | 
     (CUT_BIT(WW_LL_PT)		) | 
     (CUT_BIT(WW_OPP_SIGN)		) | 
     (CUT_BIT(WW_PASS4_MET)		) |  
     (CUT_BIT(WW_PASS2_MET)		) | 
     (CUT_BIT(WW_LT_GOOD)		) | 
     (CUT_BIT(WW_LL_GOOD)		) | 
     (CUT_BIT(WW_LT_ISO)		) | 
     (CUT_BIT(WW_LL_ISO)		) | 
     (CUT_BIT(WW_PASS_ZVETO)	) | 
     (CUT_BIT(WW_PASS_ADDZVETO)	) | 
     (CUT_BIT(WW_PASS_JETVETO_CALO)	); 

const static cuts_t ww_ss_baseline_cuts = (ww_baseline_cuts & ~(CUT_BIT(WW_OPP_SIGN))) | (CUT_BIT(WW_SAME_SIGN));

const static cuts_t ww_dumbo_cuts = (ww_baseline_cuts & ~(CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO) |
							  CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO))) | // relax iso cuts
     (CUT_BIT(WW_MU_ISO)); // so only muon has to be isolated

const static cuts_t ww_ss_dumbo_cuts = (ww_dumbo_cuts & ~(CUT_BIT(WW_OPP_SIGN))) | (CUT_BIT(WW_SAME_SIGN));

const static cuts_t ww_fakerate_cuts = ww_baseline_cuts & 
     ~(CUT_BIT(WW_LT_GOOD) | CUT_BIT(WW_LL_GOOD) |
       CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO) |
       CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO));

/* const static cuts_t ww_fakerate_cuts =  */
/* /\*      (CUT_BIT(WW_LT_PT)) |  *\/ */
/* /\*      (CUT_BIT(WW_LL_PT)) |  *\/ */
/*      (CUT_BIT(WW_OPP_SIGN)); */

const static cuts_t ww_ss_fakerate_cuts = (ww_fakerate_cuts & ~(CUT_BIT(WW_OPP_SIGN))) | (CUT_BIT(WW_SAME_SIGN));

const static cuts_t ww_baseline_cuts_nomet_nozveto = ww_baseline_cuts & 
     ~(CUT_BIT(WW_PASS2_MET) | CUT_BIT(WW_PASS4_MET) |
       CUT_BIT(WW_PASS_JETVETO_TRACKJETS) |
       CUT_BIT(WW_PASS_ZVETO) | 
       CUT_BIT(WW_PASS_ADDZVETO));

const static cuts_t ww_baseline_cuts_nomet = ww_baseline_cuts & 
     ~(CUT_BIT(WW_PASS2_MET) | CUT_BIT(WW_PASS4_MET));

class WWLooperBase : public LooperBase {
// WW looper base class.  
//
// - provides unified cut bit pattern (but actually cutting is up to
//   the derived classes)
//
// - provides function to fill all the histograms
public:
     WWLooperBase (Sample, cuts_t cuts, const char *logfilename = 0);
     virtual ~WWLooperBase () { }
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
     TH2F	held0vsRelIso, heldphiinvsRelIso;

protected:
     cuts_t		cuts_passed; 
     FILE		*logfile;
     double		cands_passing[4];
     double		cands_passing_w2[4];
     unsigned int	cands_count[4];
};

// - stack plots (including N-1)
// - tables
class WWResultsLooper : public WWLooperBase {
public:
     WWResultsLooper (Sample, 
		      cuts_t cuts = ww_baseline_cuts, 
		      const char *fname = 0);
protected:
     virtual void	Dilep 		(int idx);
};

// - background estimate for W+jets from Dumbo method
class WWDumboLooper : public WWLooperBase {
public:
     WWDumboLooper (Sample s, cuts_t cuts = ww_dumbo_cuts, const char *fname = 0);
protected:
     virtual void	Dilep 		(int idx);
};

// - background estimate for W+jets from fake rates (only electrons for now)
class WWFakeRateLooper : public WWLooperBase {
public:
     WWFakeRateLooper (Sample s, cuts_t cuts = ww_fakerate_cuts, 
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

// - stack plots and tables for same-sign control sample
class WWSSResultsLooper : public WWLooperBase {
public:
     WWSSResultsLooper (Sample s, cuts_t cuts = ww_ss_baseline_cuts, 
			const char *fname = 0);
protected:
     virtual void	Dilep 		(int idx);
};

// - mu tag efficiency estimate
class WWMuTagEffLooper : public WWLooperBase {
public:
     WWMuTagEffLooper (Sample s, cuts_t cuts,
		       cuts_t mutag_cuts,
		       const char *fname = 0);
     virtual double	CandTwoJet (enum DileptonHypType i) const { return cand_twojet[i]; }
     virtual double	CandTwoJetMuTagged (enum DileptonHypType i) const { return cand_twojet_mutagged[i]; }
     virtual double	CandTwoJetRMS (enum DileptonHypType i) const { return sqrt(cand_twojet_w2[i]); }
     virtual double	CandTwoJetMuTaggedRMS (enum DileptonHypType i) const { return sqrt(cand_twojet_mutagged_w2[i]); }
     WWMuTagEffLooper	operator + (const WWMuTagEffLooper &) const;
     const double	*MuTagEff () const { return mutag_eff; }
protected:
     virtual void	Dilep 		(int idx);
     virtual void	End ();
     cuts_t		mutag_cuts;
     double		cand_twojet[4];
     double		cand_twojet_mutagged[4];
     double		cand_twojet_w2[4];
     double		cand_twojet_mutagged_w2[4];
     double		mutag_eff[4];
};

// - top background estimate
class WWTopEstimateLooper : public WWLooperBase {
public:
     WWTopEstimateLooper (Sample s, const double mutag_eff[4], 
			  cuts_t cuts,
			  cuts_t mutag_cuts,
			  const char *fname = 0);
     virtual double	CandsPassing (enum DileptonHypType i) const { return top_prediction[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(top_prediction_w2[i]); }
protected:
     virtual void	Dilep 		(int idx);
     cuts_t		mutag_cuts;
     double		mutag_eff[4];
     double		top_prediction[4];
     double		top_prediction_w2[4];
};

class WWDYInEMuLooper : public WWLooperBase {
public:
     WWDYInEMuLooper (Sample s, cuts_t cuts = ww_baseline_cuts, 
		      const char *fname = 0)
	  : WWLooperBase(s, cuts, fname),
	    hDphi		(s,  "Dphi", 50, 0, M_PI, cuts, CUT_BIT(WW_PASS_MUON_B_VETO))
	  { }
protected:
     void		Dilep (int i_hyp)
	  {
	       cuts_t cuts_passed = DilepSelect(i_hyp);
	       const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
	       const double weight = Weight(i_hyp);
	       double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
	       if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
	       hDphi.Fill(cuts_passed, myType, dphi, weight);
	  }

protected:
     NMinus1Hist	hDphi;
};
#endif
