// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"
#include "Cuts.h"

//----------------------------------------------------------------------
// Loopers 
//----------------------------------------------------------------------

class Looper : public LooperBase {

public:
     Looper (Sample, cuts_t cuts, const char *logfilename = 0);
     virtual ~Looper () { }

protected:
     // this is where we book our histograms
     virtual void	BookHistos ();

     // filter out this event.  If FilterEvent returns true, no
     // further processing is done on this event
     virtual bool	FilterEvent();
     virtual cuts_t	EventSelect	();

	void WEvent();
	void ZEvent();

     virtual void	FillEventHistos ();
     virtual void	End		();

	// do stuff with histogram
	void FormatHist(TH1* hist);



public:
     // these functions are called by the table-printing code
     virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
     virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------

     	TH1F	*h1_lep_pt_[3];
	TH1F	*h1_lep_met_[3];
        TH1F    *h1_lep_tkIso_[3];

	TH1F	*h1_dilep_0_pt_[4];
	TH1F	*h1_dilep_1_pt_[4];
	TH1F	*h1_dilep_nhyp_;

protected:
     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing_[4];
     double		cands_passing_w2_[4];
     unsigned int	cands_count_[4];
};
#endif

