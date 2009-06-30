// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "EffH1F.h"

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


const static cuts_t stat1cuts =
( CUT_BIT(CUT_ETA24) ) |
( CUT_BIT(CUT_PT20) ) |
( CUT_BIT(CUT_IN_Z_WINDOW) ) ;


class Looper : public LooperBase {

public:
  // constructor; tell the looper what sample to loop on (see
  // Tools/Sample.h), what cuts candidates need to pass, and a file
  // name for dumping log messages
  //Looper (Sample, cuts_t cuts, const char *logfilename = 0);
  Looper (Sample, cuts_t cuts, const char *logfilename = 0, bool usew = true);
  virtual ~Looper () { }

protected:

  virtual void	BookHistos ();

  virtual bool	FilterEvent();
  virtual cuts_t LepSelect(int i, int flv); //flv=0 for el, =1 for mu
  virtual cuts_t Stat1Select(vector<int>);
  virtual cuts_t EventSelect	();
  virtual void FillEventHistos ();

  virtual void	End		();

public:
  // these functions are called by the table-printing code
  virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
  virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
  virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:

  TH1F* e_hcal_iso;
  TH1F* e_ecal_iso;
  TH1F* e_trck_iso;
  TH1F* e_hcal_iso_dr05_1;
  TH1F* e_ecal_iso_dr05_1;
  TH1F* e_trck_iso_dr05_1;

  //eff versus dr btwn two stat1 leps
  EffH1F* eff_edr_hcal_iso_dr05_1;
  EffH1F* eff_edr_ecal_iso_dr05_1;
  EffH1F* eff_edr_trck_iso_dr05_1;
  

protected:

  bool useweight;
  double weight;
  
  // count the (weighted and unweighted) number of candidates passing our cuts
  double		cands_passing_[4];
  double		cands_passing_w2_[4];
  unsigned int	cands_count_[4];
};
#endif
