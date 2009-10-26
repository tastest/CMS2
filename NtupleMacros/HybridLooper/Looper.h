// -*- C++ -*-

#ifndef LOOPER_H
#define LOOPER_H

#include "Tools/LooperBase.h"

#include "Cuts.h"
#include "EffMulti.h"


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
     // we define an analysis-specific set of FillEventHistos(),
     // FillDilepHistos(), FillTrilepHistos() and FillQuadlepHistos()
     // that fill our histograms.  
     // 
     // the framework calls our FillEventHistos() function for every event
     virtual void	FillEventHistos ();
     virtual void	End		();

	// get subdetector for histogram filling
	int getSubdet(int eleIndex);

	// study of track isolation
	//	- why is the tail so long for signal?
	void trackIsolationStudy(int eleIndex, int det);

	// function to recompute track iso
	float recomputeTrackIsolation(int eleIndex, float strip, float dRIn, float dROut, float &shCutSum);

	// w efficiency studies
	void wEfficiency();

	// results of AN2009_98 comparison
	// well... related
	void AN2009_98();

	// do stuff with histogram
	void Format2DHist(TH2F** hist, std::string name, Int_t nx, Float_t minx, Float_t maxx,
                                Int_t ny, Float_t miny, Float_t maxy);
	void FormatHist(TH1F** hist, std::string name, Int_t n, Float_t min, Float_t max);
	void FormatEffHist(EffMulti** hist, bool lessThan, 
				float thresholdEB, float ThresholdEE, std::string name);

public:
     // these functions are called by the table-printing code
     virtual double	CandsPassing (enum DileptonHypType i) const { return cands_passing_[i]; }
     virtual int       	CandsCount (enum DileptonHypType i) const { return cands_count_[i]; }
     virtual double	RMS (enum DileptonHypType i) const { return sqrt(cands_passing_w2_[i]); }

     virtual double     CandsPassingW (unsigned int i) const { return wEvents_passing_[i]; }
     virtual int        CandsCountW (unsigned int i) const { return wEvents_count_[i]; }
     virtual double     RMSW (unsigned int i) const { return sqrt(wEvents_passing_w2_[i]); }

     virtual double     CandsPassingAN2009_98 (unsigned int i) const { return AN2009_98Events_passing_[i]; }
     virtual int        CandsCountAN2009_98 (unsigned int i) const { return AN2009_98Events_count_[i]; }
     virtual double     RMSAN2009_98 (unsigned int i) const { return sqrt(AN2009_98Events_passing_w2_[i]); }


protected:
     //----------------------------------------------------------------------
     // declare your histograms here:
     //----------------------------------------------------------------------

	// general
	//
        TH1F    *h1_pt_[2];
        TH1F    *h1_eta_[2];
	TH1F	*h1_phi_[2];
        TH1F    *h1_wwIsoAll_[2];

	// isolation
	//
	TH1F 	*h1_wwIso_[2];
	TH1F	*h1_tkIso03All_[2];
	TH1F	*h1_ecalIso03All_[2];
	TH1F	*h1_hcalIso03All_[2];
	TH1F	*h1_caloIso03All_[2];

	TH2F	*h2_tkIso03All_[2];
        TH2F    *h2_ecalIso03All_[2];
        TH2F    *h2_hcalIso03All_[2];
	TH2F	*h2_caloIso03All_[2];

	// N-1
        TH1F    *h1_tkIso03AllNM1_[2];
        TH1F    *h1_tkIso03AllIDNM1_[2];
        TH1F    *h1_tkIso03AllConvNM1_[2];
        TH1F    *h1_tkIso03AllConvIDNM1_[2];

        TH1F    *h1_ecalIso03AllNM1_[2];
        TH1F    *h1_hcalIso03AllNM1_[2];
	TH1F	*h1_tkIso03AllReJura01In015NM1_[2];
        TH1F    *h1_tkIso03AllReJura01In015IDNM1_[2];
        TH1F    *h1_tkIso03AllReJura01In015ConvNM1_[2];
        TH1F    *h1_tkIso03AllReJura01In015ConvIDNM1_[2];
	TH1F	*h1_tkIso03AllReShCutNM1_[2];

	// track isolation investigation
        TH1F    *h1_tkIso03AllRe_[2];
	TH1F	*h1_tkIso03AllReShVeto_[2];
        TH1F    *h1_tkIso03AllReRel_[2];
	TH1F	*h1_tkIso03Alld0corr_[2];

        TH1F    *h1_tkIso03AllReJura01In015_[2];

	TH1F 	*h1_tkIso03AllReJura01_[2];
        TH1F    *h1_tkIso03AllReJura02_[2];
        TH1F    *h1_tkIso03AllReJura03_[2];


	TH1F	*h1_tkIso03AllRedEta_[2];
	TH1F    *h1_tkIso03AllRedPhi_[2];
	TH2F	*h2_tkIso03AllRedR2D_[2];

	// electron ID related
	//
	TH1F	*h1_dEtaIn_[2];
	TH1F	*h1_dPhiIn_[2];
	TH1F	*h1_dPhiInSigned_[2];
	TH1F	*h1_hoe_[2];
	TH1F	*h1_sigmaIEtaIEta_[2];
        TH1F    *h1_sigmaIPhiIPhi_[2];
	TH1F	*h1_E2x5Norm5x5_[2];
	TH1F    *h1_E1x5Norm5x5_[2];
	TH1F	*h1_eopIn_[2];
	TH1F 	*h1_d0corr_[2];
	TH1F	*h1_closestMuon_[2];

	//EffMulti	*em_tasElectronV1_[2];

	//
	// W effciency studies related
	//
        TH1F *h1_weff_pt_[2];
	TH1F *h1_weff_iso_[2];
	TH1F *h1_weff_tcmet_[2];
	TH1F *h1_weff_jptpt_[2];
	TH1F *h1_weff_tcmet_after_iso_[2];
	TH1F *h1_weff_jptpt_after_iso_[2];
	TH1F *h1_weff_leadjptphi_after_iso_[2];
        TH1F *h1_weff_jptphimax_after_iso_[2];
        TH1F *h1_weff_d0corr_after_iso_[2];

        TH1F *h1_weff_tcmet_after_iso_jpt_[2];
        TH1F *h1_weff_tcmet_after_iso_jpt_conv_[2];


        TH1F *h1_weff_leadjptphi_after_iso_jpt_[2];
        TH1F *h1_weff_d0corr_after_iso_jpt_[2];
        TH1F *h1_weff_leadjptphi_after_iso_jpt_tcmet_[2];
	TH1F *h1_weff_leastemjpt_after_iso_[2];

	// the actual distributions
	TH1F *h1_weffs_sigmaIEtaIEta_[2];
        TH1F *h1_weffbg_sigmaIEtaIEta_[2];

        EffMulti *em_dEtaIn_[2];
        EffMulti *em_dPhiIn_[2];
        EffMulti *em_hoe_[2];           
        EffMulti *em_sieie_[2];

        EffMulti *em_classBasedTight_[2];
        EffMulti *em_robustTight_[2];

        EffMulti *em_eopInLT30_[2];
        EffMulti *em_eopInGT05_[2];

	EffMulti *em_tasElectronV1_[2];

	//
	// AN2009-98 studies related
	//

	TH1F *h1_AN2009_098_pt2_[2];
        TH1F *h1_AN2009_098_eta1_[2];

        TH1F *h1_AN2009_098_tkIso_[2];
        TH1F *h1_AN2009_098_ecalIso_[2];
        TH1F *h1_AN2009_098_hcalIso_[2];
        TH1F *h1_AN2009_098_tcmet_[2];

        TH1F *h1_AN2009_098_tcmet_after_selection_[2];	

protected:

     unsigned int       AN2009_98Events_count_[3];
     double             AN2009_98Events_passing_[3];
     double             AN2009_98Events_passing_w2_[3];


     unsigned int	wEvents_count_[3];
     double 		wEvents_passing_[3];
     double		wEvents_passing_w2_[3];

     // count the (weighted and unweighted) number of candidates passing our cuts
     double		cands_passing_[4];
     double		cands_passing_w2_[4];
     unsigned int	cands_count_[4];
};
#endif

