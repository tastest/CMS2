// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "EffH1F.h"
#include "EffH2F.h"

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
	 CUT_PT20,
	 CUT_ETA24,
	 CUT_MU_GLOBAL,
	 CUT_EL_DR,
	 CUT_MU_DR,
};
 //my dR matching--in place of mc_el
	 //CUT_MC_EL,
	 //CUT_MC3_EL,
	 //CUT_MC_MU,
	 //CUT_MC3_MU,


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

//const static cuts_t el_base =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC_EL) ) |
//( CUT_BIT(CUT_MC3_EL) ) |
//( CUT_BIT(CUT_EL_GOOD) ) |
//( CUT_BIT(CUT_EL_ISO) ) ;

//const static cuts_t mu_base =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC_MU) ) |
//( CUT_BIT(CUT_MC3_MU) ) |
//( CUT_BIT(CUT_MU_GOOD) ) |
//( CUT_BIT(CUT_MU_ISO) ) ;

const static cuts_t stat1cuts =
( CUT_BIT(CUT_ETA24) ) |
( CUT_BIT(CUT_PT20) ) |
( CUT_BIT(CUT_IN_Z_WINDOW) ) ;

//ISO CUTS
//const static cuts_t els_iso_denom =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC_EL) ) |
//( CUT_BIT(CUT_EL_GOOD) );

//const static cuts_t els_iso3_denom =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC3_EL) ) |
//( CUT_BIT(CUT_EL_GOOD) );

//const static cuts_t els_iso_numer = ( els_iso_denom | CUT_BIT(CUT_EL_ISO) );
//const static cuts_t els_iso3_numer = ( els_iso3_denom | CUT_BIT(CUT_EL_ISO) );

//const static cuts_t mus_iso_denom =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC_MU) ) |
//( CUT_BIT(CUT_MU_GOOD) );

//const static cuts_t mus_iso3_denom =
//( CUT_BIT(CUT_OPP_SIGN) ) |
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) |
//( CUT_BIT(CUT_MC3_MU) ) |
//( CUT_BIT(CUT_MU_GOOD) );

//const static cuts_t mus_iso_numer = mus_iso_denom | CUT_BIT(CUT_MU_ISO);
//const static cuts_t mus_iso3_numer = mus_iso3_denom | CUT_BIT(CUT_MU3_ISO);

//RECO CUTS
//no need b'c it's just the stat1 plus dr match
//const static cuts_t els_reco_denom =
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) );
//
//const static cuts_t els_reco_numer = els_reco_denom | CUT_BIT(CUT_MC_EL) ;
//
//const static cuts_t els_reco3_denom = els_reco_denom;
//
//const static cuts_t els_reco3_numer = els_reco3_denom | CUT_BIT(CUT_MC3_EL);
//
//const static cuts_t mus_reco_denom =
//( CUT_BIT(CUT_ETA24) ) |
//( CUT_BIT(CUT_PT20) ) |
//( CUT_BIT(CUT_IN_Z_WINDOW) ) ;
//
//const static cuts_t mus_reco_numer = mus_reco_denom | CUT_BIT(CUT_MC_MU);
//
//const static cuts_t mus_reco3_denom = mus_reco_denom;
//
//const static cuts_t mus_reco3_numer = mus_reco3_denom | CUT_BIT(CUT_MC3_MU);


struct counts {
  double total; //all events in sample
  double denom; //all denominator events
  double numer; //all numerator events
  double onelep; //events w/ only 1 lepton passing numerator
  double failnumer; //events which pass denom but fail numer (0 pass)
  double geneta; //num evts which don't have 2 leptons in eta 2.5
  double gencuts; //num evts which fail the gen cut
  //double opp_sign; //num which fail numerator because of this (but i'm not cutting on it, so this is useless)
  double same_flv;
  double pt20; //ie, this is num events which aren't in numerator b'c of pt cut
  double el_good;
  double bad_mom; //mc_motherid != 23
  double no_match_z; // (num leptons with mc_motherid == 23) != 2
  //double multihyp; //evt has more that two pairs of leptons passing denom--no longer possible
  double dupematch; //both stat 1 match the same reco
  //void fill //not sure about this: don't want a separate fill for each data member, but passing in a data member is dumb?
};


class Looper : public LooperBase {

public:
  // constructor; tell the looper what sample to loop on (see
  // Tools/Sample.h), what cuts candidates need to pass, and a file
  // name for dumping log messages
  Looper (Sample, cuts_t cuts, const char *logfilename = 0, bool usew = true);
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
  virtual cuts_t LepSelect(int i, int flv); //flv=0 for el, =1 for mu
  //virtual void PairSelect( int flv );
  virtual cuts_t Stat1Select(vector<int>);
  virtual cuts_t EventSelect();
  //virtual cuts_t CutCount(cuts_t, cuts_t, cuts_t, int);

  virtual void FillEventHistos();
  //virtual void FillStat1Histos( int , vector<int>, vector<int> ); 
  //virtual void FillStat1Histos();
  virtual void End		();

public:
  // these functions are called by the table-printing code
  virtual double CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
  virtual int CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
  virtual double RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:

  //TH1F* hgen_z_mass;
  //TH1F* hgen_z_p;
  //TH1F* hgen_z_pt;
  //TH1F* hgen_z_eta;
  //TH1F* hgen3lep_z_mass; //status 3 just z
  //const unsigned int nlepplots = 2;
  //#define nlepplots 2
  //TH1F* hgen1_z_eta;
  //TH1F* hgen1_lep_mass;
  //TH1F* hgen1_lep_pt[nlepplots];
  //TH1F* hgen3_lep_pt[nlepplots];
  //TH1F* hgen1_lep_eta[nlepplots];
  //
  //TH1F* hels_size;
  //TH1F* hmus_size;
  //TH1F* hmus_type;
  //
  //TH1F* hels_iso;
  //TH1F* hels_chg;
  //TH1F* hmus_iso;

  //kinematic (not efficiency) plots
  TH1F* hels_dr_recostat1[2]; //one for each lepton from Z
  TH1F* e_eta[2];
  TH1F* e_dr; //dr of stat1
  TH2F* e_dr_zpt; //dr of stat1 vs zpt

  //efficiency plots
  //reco
  EffH1F* eff_ep_reco[2]; //eff as fn of leading, second e p
  EffH1F* eff_ept_reco[2]; //pt 1,2
  EffH1F* eff_eeta_reco[2]; //eff as fn of leading, second e eta
  EffH1F* eff_edr_reco; //eff versus dr btwn two stat1 leps
  EffH1F* eff_zp_ind_reco; // Z plots
  EffH1F* eff_zpt_ind_reco;
  EffH1F* eff_zp_pr_reco; // Z plots
  EffH1F* eff_zpt_pr_reco;
  //unrepeated reco
  EffH1F* ineff_ep_reco[2]; //ineffienciey is ones which fail
  EffH1F* ineff_ept_reco[2];
  EffH1F* eff_ep_sum_reco; //sum of e1, e2
  EffH1F* eff_ept_sum_reco;

  //id
  EffH1F* eff_ep_id[2]; //eff as fn of leading, second e p
  EffH1F* eff_ept_id[2]; //pt 1,2
  EffH1F* eff_eeta_id[2]; //eff as fn of leading, second e eta
  EffH1F* eff_edr_id; //eff versus dr btwn two stat1 leps
  EffH1F* eff_zp_ind_id;
  EffH1F* eff_zpt_ind_id;
  EffH1F* eff_zp_pr_id; 
  EffH1F* eff_zpt_pr_id;
  
  //iso
  EffH1F* eff_ep_iso[2]; //eff as fn of leading, second e p
  EffH1F* eff_ept_iso[2]; //pt 1,2
  EffH1F* eff_eeta_iso[2]; //eff as fn of leading, second e eta
  EffH1F* eff_edr_iso; //eff versus dr btwn two stat1 leps
  EffH1F* eff_zp_ind_iso;
  EffH1F* eff_zpt_ind_iso;
  EffH1F* eff_zp_pr_iso;
  EffH1F* eff_zpt_pr_iso;

  //id+iso
  EffH1F* eff_ep_id_iso[2]; //eff as fn of leading, second e p
  EffH1F* eff_ept_id_iso[2]; //pt 1,2
  EffH1F* eff_eeta_id_iso[2]; //eff as fn of leading, second e eta
  EffH1F* eff_edr_id_iso; //eff versus dr btwn two stat1 leps
  EffH1F* eff_zp_ind_id_iso;
  EffH1F* eff_zpt_ind_id_iso;
  EffH1F* eff_zp_pr_id_iso; 
  EffH1F* eff_zpt_pr_id_iso;

  //EffH2F* eff_p_eta_reco;
  //EffH2F* eff_pt_eta_reco;
  //
  //#define netabins 3
  //EffH1F* eff_p_bineta_reco[netabins];
  //EffH1F* eff_pt_bineta_reco[netabins];
  
  //jet EffH's
  //EffH1F* eff_njets_reco;
  //EffH1F* eff_jetEt_reco;
  //EffH2F* eff_pt_njets_reco;
  //EffH2F* eff_pt_jetEt_reco;

protected:
  bool useweight;

  //event_cuts evt_cuts[2]; //struct defined above, one for e, mu
  //#define lepeffs 2
  //int elsidx[lepeffs]; //idx of passing numerator, denominator--REUSE for different cuts!!! 

#define ncounts 1  //should be n cuts (not counting e,mu diffs)
  counts count[ncounts]; //my struct, declared above Looper
  double weight;
  int denomitr; 

  // count the (weighted and unweighted) number of candidates passing our cuts
  double		cands_passing_[4];
  double		cands_passing_w2_[4];
  unsigned int	cands_count_[4];
};
#endif

