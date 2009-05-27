// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#define MAXELE 50

#include "Tools/LooperBase.h"

// root includes
#include "TTree.h"
#include "TFile.h"

// List of all cuts that can be applied.  The cuts are handled as a
enum {
     CUT_LT_PT,
};

const static cuts_t baseline_cuts = (CUT_BIT(CUT_LT_PT));

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

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
     // the framework calls our FillEventHistos() function for every event
     virtual void	FillEventHistos ();
     virtual void	End		();

	// matching
	float match(unsigned int leptonIndex);


public:
     // these functions are called by the table-printing code
     virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
     virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:

     // file and tree (necessary)
     TFile *outFile_;
     TTree *outTree_;

     // variables (you choose)
     Int_t sampleId_;

     // event properties
	Float_t z_pt_;
	Float_t z_p_;

     // electron variables
     Int_t ele_count_;
     Float_t ele_pt_[MAXELE];
     Float_t ele_p_[MAXELE];
     Float_t ele_eta_[MAXELE];
     Float_t ele_etaDet_[MAXELE];

     Float_t ele_hOverE_[MAXELE];
     Float_t ele_dPhiIn_[MAXELE];
     Float_t ele_dEtaIn_[MAXELE];
     Float_t ele_sigmaIEtaIEta_[MAXELE];
     Float_t ele_tkIso_[MAXELE];

	Int_t ele_matchMC_[MAXELE];


     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing_[4];
     double		cands_passing_w2_[4];
     unsigned int	cands_count_[4];
};
#endif
