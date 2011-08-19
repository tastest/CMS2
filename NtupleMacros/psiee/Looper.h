// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "CORE/selections.h"
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
     CUT_SAME_SIGN,
     CUT_OPP_SIGN,
     CUT_PASS2_MET,
     CUT_PASS4_MET,
     CUT_PASS2_METCORR,
     CUT_PASS4_METCORR,
     CUT_LT_GOOD,
     CUT_LL_GOOD,
     CUT_EL_GOOD,
     CUT_EL_GOOD_NO_D0,
     CUT_LT_TIGHT_DPHIIN,
     CUT_LL_TIGHT_DPHIIN,
     CUT_MU_GOOD,
     CUT_ONE_SUPERTIGHT,
     CUT_TWO_SUPERTIGHT,
     CUT_LT_ISO,
     CUT_LL_ISO,
     CUT_ONE_ISO,
     CUT_TWO_ISO,
     CUT_EL_ISO,
     CUT_MU_ISO,
     CUT_LT_CALOISO,
     CUT_LL_CALOISO,
     CUT_ONE_CALOISO,
     CUT_TWO_CALOISO,
     CUT_EL_CALOISO,
     CUT_PASS_ZVETO,
     CUT_IN_Z_WINDOW,
     CUT_PASS_ADDZVETO,
     CUT_PASS_JETVETO_CALO,
     CUT_PASS_JETVETO_TRACKJETS,
     CUT_PASS_JETVETO_CALOTRACKJETS_COMBO,
     CUT_PASS_EXTRALEPTON_VETO,
     CUT_EL_BARREL,
     CUT_ELFAKE_FAKEABLE_OBJECT,
     CUT_ELFAKE_NUMERATOR,
     CUT_ELFAKE_NOT_NUMERATOR,
     CUT_MORE_THAN_TWO_TRACKS,
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
     (CUT_BIT(CUT_MORE_THAN_TWO_TRACKS)) |
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_PASS4_METCORR)	) |  
     (CUT_BIT(CUT_PASS2_METCORR)	) |  
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_CALOISO)	) |  
     (CUT_BIT(CUT_LL_CALOISO)	) |  
     (CUT_BIT(CUT_PASS_ZVETO)	) | 
     (CUT_BIT(CUT_PASS_ADDZVETO)	) | 
     (CUT_BIT(CUT_PASS_JETVETO_CALO)	) |
  (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)	) ;

// baseline + tight delta phi in
const static cuts_t baseline_tight_dphiin_cuts = baseline_cuts | CUT_BIT(CUT_LT_TIGHT_DPHIIN) | CUT_BIT(CUT_LL_TIGHT_DPHIIN);

const static cuts_t baseline_metcorr_cuts = (baseline_cuts & ~(CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET)))
		   | CUT_BIT(CUT_PASS4_METCORR) | CUT_BIT(CUT_PASS2_METCORR);

// these cuts are used to measure the mu tagging efficiency for top
const static cuts_t baseline_mu_tageff_cuts = baseline_cuts & ~((CUT_BIT(CUT_PASS_JETVETO_CALO)	) | 
								(CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)	) | 
								(CUT_BIT(CUT_PASS_EXTRALEPTON_VETO)	) );

const static cuts_t baseline_no_trackjets_cuts = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)); 

const static cuts_t baseline_no_caloiso_cuts = (baseline_cuts & ~((CUT_BIT(CUT_LT_CALOISO)	) | (CUT_BIT(CUT_LL_CALOISO)	))) | 
     (CUT_BIT(CUT_LT_ISO)		) | (CUT_BIT(CUT_LL_ISO)		);

const static cuts_t old_baseline_cuts = 
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_PASS4_MET)		) |  
     (CUT_BIT(CUT_PASS2_MET)		) | 
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_ISO)		) | 
     (CUT_BIT(CUT_LL_ISO)		) | 
     (CUT_BIT(CUT_PASS_ZVETO)	) | 
     (CUT_BIT(CUT_PASS_ADDZVETO)	) | 
     (CUT_BIT(CUT_PASS_JETVETO_CALO)	); 

const static cuts_t ss_baseline_cuts = (baseline_cuts & ~(CUT_BIT(CUT_OPP_SIGN))) | (CUT_BIT(CUT_SAME_SIGN));

const static cuts_t dumbo_cuts = 
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_PASS4_MET)		) |  
     (CUT_BIT(CUT_PASS2_MET)		) |  
     (CUT_BIT(CUT_MU_GOOD)		) | 
     (CUT_BIT(CUT_EL_GOOD)		) | 
     (CUT_BIT(CUT_MU_ISO)	) |  
     (CUT_BIT(CUT_PASS_ADDZVETO)	) | 
     (CUT_BIT(CUT_PASS_JETVETO_CALO)	) |
								(CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)	) ;


const static cuts_t dumbo_cuts_nod0 = (dumbo_cuts & ~CUT_BIT(CUT_EL_GOOD)) | CUT_BIT(CUT_EL_GOOD_NO_D0);

const static cuts_t ss_dumbo_cuts = (dumbo_cuts & ~(CUT_BIT(CUT_OPP_SIGN))) | (CUT_BIT(CUT_SAME_SIGN));

const static cuts_t fakerate_cuts = baseline_cuts & 
     ~(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
       CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) |
       CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO));

const static cuts_t ss_fakerate_cuts = (fakerate_cuts & ~(CUT_BIT(CUT_OPP_SIGN))) | (CUT_BIT(CUT_SAME_SIGN));

const static cuts_t baseline_cuts_nomet_nozveto = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_MET) |
       CUT_BIT(CUT_PASS_JETVETO_TRACKJETS) |
       CUT_BIT(CUT_PASS_ZVETO) | 
       CUT_BIT(CUT_PASS_ADDZVETO));

const static cuts_t baseline_cuts_nomet = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_MET));

const static cuts_t baseline_cuts_zwindow = (baseline_cuts & 
						~(CUT_BIT(CUT_PASS_ZVETO) | CUT_BIT(CUT_PASS_ADDZVETO))) |
		   CUT_BIT(CUT_IN_Z_WINDOW);

const static cuts_t baseline_cuts_zwindow_trkcorr = (baseline_cuts_zwindow & 
							~(CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_MET))) |
		   CUT_BIT(CUT_PASS2_METCORR) | CUT_BIT(CUT_PASS4_METCORR);

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

     // the framework calls our FillEventHistos() function for every event
     virtual void	FillEventHistos ();
     virtual void 	End ();

     static double	DeltaR (int itrk, const LorentzVector &v);


protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------
     TH1F		*h_m_els;
     TH1F		*h_m_eltrks;
     TH1F		*h_m_sctrks;

     TH1F		*h_m_elsctrk;

     TH1F 		*h_m_els_pt_ee;
     TH1F		*h_m_els_pt_eb;
     TH1F               *h_m_els_et_ee;
     TH1F               *h_m_els_et_eb;

     TH1F               *h_m_elspf_pt_ee;
     TH1F               *h_m_elspf_pt_eb;
     TH1F               *h_m_elspf_et_ee;
     TH1F               *h_m_elspf_et_eb;

	TH1F	*h1_mass_eltrkeltrk_;
        TH1F    *h1_mass_sctrkeltrk_;
        TH1F    *h1_mass_sctrksctrk_;

};
#endif
