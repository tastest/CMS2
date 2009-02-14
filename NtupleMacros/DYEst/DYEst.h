// -*- C++ -*-

#ifndef DYEST_H
#define DYEST_H

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
	  CUT_PASS2_TCMET,
	  CUT_PASS4_TCMET,
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
	  CUT_PASS_JETVETO_JPT20,
	  CUT_PASS_JETVETO_JPT25,
	  CUT_PASS_MUON_B_VETO,	
	  CUT_MUON_TAGGED,
	  CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT,	
	  CUT_MUON_TAGGED_WITHOUT_PTCUT,
	  CUT_PASS_EXTRALEPTON_VETO,
	  CUT_EL_BARREL,
	  CUT_ELFAKE_FAKEABLE_OBJECT,
	  CUT_ELFAKE_NUMERATOR,
	  CUT_ELFAKE_NOT_NUMERATOR,
	  CUT_MORE_THAN_TWO_TRACKS,

     CUT_MUON_RECO_CLEANING,
     CUT_MUON_RECO_CLEANING20,
     CUT_MET_SIMPLE20,
     CUT_MET_SIMPLE35,
     CUT_MET_SIMPLE45,
     CUT_MET_BALLANCE,
     CUT_MET_PROJECTED,
     CUT_PASS5_MET,


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
//      (CUT_BIT(CUT_PASS_ADDZVETO)	) | 
//      (CUT_BIT(CUT_PASS_JETVETO_CALO)	) |
//      (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)	) |  
     (CUT_BIT(CUT_PASS_JETVETO_JPT20)	) |  
     (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)	);   


// current baseline cuts with pass5 met

// this is the current baseline set of cuts
const static cuts_t baseline_cuts_pass5 = 
     (CUT_BIT(CUT_LT_PT)                ) |
     (CUT_BIT(CUT_LL_PT)                ) |
     (CUT_BIT(CUT_OPP_SIGN)             ) |
//     (CUT_BIT(CUT_PASS4_TCMET)          ) |
//     (CUT_BIT(CUT_PASS2_TCMET)          ) |
     (CUT_BIT(CUT_PASS5_MET)          ) |
//          (CUT_BIT(CUT_MET_SIMPLE45)   ) |
//          (CUT_BIT(CUT_MET_BALLANCE)   ) |
//          (CUT_BIT(CUT_MET_PROJECTED)   ) |
     (CUT_BIT(CUT_LT_GOOD)              ) |
     (CUT_BIT(CUT_LL_GOOD)              ) |
     (CUT_BIT(CUT_LT_CALOISO)   ) |
     (CUT_BIT(CUT_LL_CALOISO)   ) |
     (CUT_BIT(CUT_PASS_ZVETO)   ) |
//      (CUT_BIT(CUT_PASS_ADDZVETO)     ) | 
//      (CUT_BIT(CUT_PASS_JETVETO_CALO) ) |
//      (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)    ) |  
     (CUT_BIT(CUT_PASS_JETVETO_JPT20)   ) |  
     (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)       );


const static cuts_t jet_z_veto_cuts =
        (CUT_BIT(CUT_PASS_ZVETO)) |
        (CUT_BIT(CUT_PASS_JETVETO_JPT20));

const static cuts_t simple_met_cuts =
        (CUT_BIT(CUT_MET_SIMPLE20)) |
        (CUT_BIT(CUT_MET_SIMPLE45));

// cuts used in the Feb 08 presentation by fkw
const static cuts_t feb_baseline_cuts = 
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

// + fix for broken CSA07 alpgen events 
const static cuts_t feb_baseline_with_ntrks_cuts = 
     feb_baseline_cuts | CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);

// + tcmet instead of hyp_met
const static cuts_t feb_baseline_with_tcmet_cuts = (feb_baseline_with_ntrks_cuts & 
					     ~(CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET)))
     | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET);

// + trkjet veto
const static cuts_t feb_baseline_with_trackjets_cuts = 
     feb_baseline_with_ntrks_cuts | CUT_BIT(CUT_PASS_JETVETO_TRACKJETS); 

// + extra muon veto 
const static cuts_t feb_baseline_with_btags_cuts =
     feb_baseline_with_ntrks_cuts | CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT); 

// use trk+calo instead of track-based rel iso for electrons
const static cuts_t feb_baseline_with_caloiso_cuts = (feb_baseline_with_ntrks_cuts & 
						      ~(CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO))) | 
     CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO);

// skip pass4 met cut (in emu, this is the same as not making a pMET cut)
const static cuts_t feb_baseline_no_pass4met_cuts = feb_baseline_with_ntrks_cuts & 
     ~(CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS4_TCMET));

// cuts used in the Oct 08 presentations by J.Mü Dima (except
// that we used a rel iso > 0.9 cut for the electron CUT_L*_CALOISO)
const static cuts_t oct_baseline_cuts = 
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_PASS4_MET)		) |  
     (CUT_BIT(CUT_PASS2_MET)		) |  
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_CALOISO)	) |  
     (CUT_BIT(CUT_LL_CALOISO)	) |  
     (CUT_BIT(CUT_PASS_ZVETO)	) | 
     (CUT_BIT(CUT_PASS_ADDZVETO)	) | 
     (CUT_BIT(CUT_PASS_JETVETO_CALO)	) |
     (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)	) |  
     (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)	);   

// replace tcmet with hyp_met for met cuts
const static cuts_t baseline_no_tcmet_cuts = (baseline_cuts & 
					     ~(CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)))
     | CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET);

// skip trkjet veto
const static cuts_t baseline_no_trackjets_cuts = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)); 

// skip extra muon veto 
const static cuts_t baseline_no_btags_cuts = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS_MUON_B_VETO) | CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)); 

// use track-based rel iso for electrons instead of trk+calo
const static cuts_t baseline_no_caloiso_cuts = (baseline_cuts & 
						~(CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO))) | 
     CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO);

// skip pass4 met cut (in emu, this is the same as not making a pMET cut)
const static cuts_t baseline_no_pass4met_cuts = baseline_cuts & 
     ~(CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS4_TCMET));

// skip ntrks > 2 cut
const static cuts_t baseline_no_ntrks_cuts = baseline_cuts & ~CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);

// denominator object cuts for the fake rate prediction 
const static cuts_t fakerate_denominator_cuts = (baseline_cuts & 
						 ~(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
						   CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) |
						   CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO))) |
     CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
						 
// numerator object cuts for the fake rate prediction 
const static cuts_t fakerate_numerator_cuts = 
     fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NUMERATOR);

// denominator and not numerator (this is the yield that should be
// multiplied by FR / (1 - FR))
const static cuts_t fakerate_denominator_not_numerator_cuts = 
     fakerate_denominator_cuts | CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);

const static cuts_t oingo_cuts = 
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_PASS4_MET)		) |  
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_CALOISO)	) |  
     (CUT_BIT(CUT_LL_CALOISO)	);   

const static cuts_t fakerate_histat_cuts = 
     (CUT_BIT(CUT_LT_PT)		) | 
     (CUT_BIT(CUT_LL_PT)		) | 
     (CUT_BIT(CUT_OPP_SIGN)		) | 
     (CUT_BIT(CUT_LT_GOOD)		) | 
     (CUT_BIT(CUT_LL_GOOD)		) | 
     (CUT_BIT(CUT_LT_CALOISO)	) |  
     (CUT_BIT(CUT_LL_CALOISO)	);

const static cuts_t fakerate_histat_denominator_cuts = (fakerate_histat_cuts & 
							~(CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
							  CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) |
							  CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO))) |
     CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
						 
// denominator and not numerator (this is the yield that should be
// multiplied by FR / (1 - FR))
const static cuts_t fakerate_histat_denominator_not_numerator_cuts = 
     fakerate_histat_denominator_cuts | CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);

static const cuts_t fakerate_ss_numerator_cuts = 
     (fakerate_numerator_cuts & ~CUT_BIT(CUT_OPP_SIGN)) | CUT_BIT(CUT_SAME_SIGN);

static const cuts_t fakerate_ss_denominator_not_numerator_cuts = 
     (fakerate_denominator_not_numerator_cuts & ~CUT_BIT(CUT_OPP_SIGN)) 
     | CUT_BIT(CUT_SAME_SIGN);


//----------------------------------------------------------------------
// DYEsts 
//----------------------------------------------------------------------

// DYEst for a dilepton analysis.  
//
// - switching between files, removing duplicates and other technical stuff is handled by DYEstBase
// - analysis-specific stuff is defined here: 
//    * declaring histograms and booking them
//    * naming cuts (using the enum from above) and checking which cuts are passed
//    * filling histograms
//    * counting passing candidates
class DYEst : public LooperBase {

public:
     // constructor; tell the looper what sample to loop on (see
     // Tools/Sample.h), what cuts candidates need to pass, and a file
     // name for dumping log messages
     DYEst (Sample, cuts_t cuts, const char *logfilename = 0);
     virtual ~DYEst () { }

protected:
     // this is where we book our histograms
     virtual void	BookHistos ();
     // we define an analysis-specific DilepSelect() that checks which
     // cuts a dilepton candidate passed
     virtual cuts_t	DilepSelect 	(int idx);
     // we define an analysis-specific FillDilepHistos() that fills our histograms
     virtual void	FillDilepHistos	(int idx);
     // at the end of the loop, we want to print a quick status message
     virtual void	End		();

public:
     // these functions are called by the table-printing code
     virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing[i]; }
     virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2[i]); }

protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------

     NMinus1Hist	*hnm1_mll_0j_;
     NMinus1Hist        *hnm1_mll_1j_;
     NMinus1Hist        *hnm1_mll_2j_;

     NMinus1Hist        *hnm1_met_0j_in_;
     NMinus1Hist        *hnm1_met_1j_in_;
     NMinus1Hist        *hnm1_met_2j_in_;

     NMinus1Hist        *hnm1_met_0j_out_;
     NMinus1Hist        *hnm1_met_1j_out_;
     NMinus1Hist        *hnm1_met_2j_out_;

	Long_t n_WZ_ee_;
	Long_t n_ZZ_ee_;
	Long_t n_Total_ee_;

        Long_t n_WZ_mm_;
        Long_t n_ZZ_mm_;
        Long_t n_Total_mm_;

	float test_;

protected:
     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing[4];
     double		cands_passing_w2[4];
     unsigned int	cands_count[4];
};
#endif
