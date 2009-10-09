#include <math.h>
#include "TVector3.h"
//#include "Math/VectorUtil.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Looper.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(dcands_passing_	    , 0, sizeof(dcands_passing_       ));
     memset(dcands_passing_w2_	, 0, sizeof(dcands_passing_w2_    ));
     memset(dcands_count_		, 0, sizeof(dcands_count_         ));
     memset(scands_passing_	    , 0, sizeof(scands_passing_       ));
     memset(scands_passing_w2_	, 0, sizeof(scands_passing_w2_    ));
     memset(scands_count_		, 0, sizeof(scands_count_         ));

	 //initialize data members
	 transmass = 0;
	 elidxs[0] = -1;
	 elidxs[1] = -1;
	 muidxs[0] = -1;
	 muidxs[1] = -1;

}

void Looper::FormatHist(TH1* hist)
{
	hist->SetFillColor(sample_.histo_color);
}

void Looper::NewHist(TH1F*& h, char* name, char* title, int bins, double min, double max) {
  h = new TH1F(name, title, bins, min, max);
  h->Sumw2();
  h->SetFillColor(sample_.histo_color);
  h->SetLineColor(sample_.histo_color);
}


//to add: mass, transverse mass, njets
void Looper::BookHistos ()
{

  // single lepton histograms (two + 1 types)
  for (unsigned int i = 0; i < 3; ++i) {
	std::string hyp = "e";
	if (i == 1) hyp = "m";
	else if (i == 2) hyp = "all";

	h1_lep_Highpt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_pt", hyp.c_str()), 
								  "Highlep_pt", 100, 0.0, 100.0);
	FormatHist(h1_lep_Highpt_[i]);
        
	h1_lep_HighptMet_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_Met", hyp.c_str()), 
									 "Highlep_Met", 100, 0.0, 100.0);
	FormatHist(h1_lep_HighptMet_[i]);
        
	h1_lep_HighptRelIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_RelIso", hyp.c_str()), 
										"Highlep_RelIso", 120, -0.1, 1.1);
	FormatHist(h1_lep_HighptRelIso_[i]);
        
	h1_lep_HighptRelIsoPtLg20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "Highlep_RelIsoPtLg20", hyp.c_str()), 
											  "Highlep_RelIsoPtLg20", 120, -0.1, 1.1);
	FormatHist(h1_lep_HighptRelIsoPtLg20_[i]);

	////
	h1_lep_Lowpt_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_pt", hyp.c_str()), 
								 "LowLep_pt", 100, 0.0, 100.0);
	FormatHist(h1_lep_Lowpt_[i]);
        
	h1_lep_LowptMet_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_Met", hyp.c_str()), 
									"LowLep_Met", 100, 0.0, 100.0);
	FormatHist(h1_lep_LowptMet_[i]);
        
	h1_lep_LowptRelIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_RelIso", hyp.c_str()), 
									   "LowLep_RelIso", 120, -0.1, 1.1);
	FormatHist(h1_lep_LowptRelIso_[i]);
        
	h1_lep_LowptRelIsoPtLg20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_RelIsoPtLg20", hyp.c_str()), 
											 "LowLep_RelIsoPtLg20", 120, -0.1, 1.1);
	FormatHist(h1_lep_LowptRelIsoPtLg20_[i]);

	h1_lep_LowptNLepGt10Lt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt10Lt20", hyp.c_str()), 
											 "LowLep_NLepGt10Lt20", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt10Lt20_[i]);

	h1_lep_LowptNLepGt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20", hyp.c_str()), 
										 "LowLep_NLepGt20", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20_[i]);

	h1_lep_LowptNLepGt20tightIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20tightIDIso0_1", hyp.c_str()), 
													  "LowLep_NLepGt20tightIDIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20tightIDIso0_1_[i]);

	h1_lep_LowptNLepGt20looseIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20looseIDIso0_1", hyp.c_str()), 
													  "LowLep_NLepGt20looseIDIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20looseIDIso0_1_[i]);

	h1_lep_LowptNLepGt20vlooseIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20vlooseIDIso0_1", hyp.c_str()), 
													   "LowLep_NLepGt20vlooseIDIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20vlooseIDIso0_1_[i]);

	h1_lep_LowptNLepGt20NOIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20NOIDIso0_1", hyp.c_str()), 
												   "LowLep_NLepGt20NOIDIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20NOIDIso0_1_[i]);

        h1_lep_LowptNLepGt20tightIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20tightIDNoIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20tightIDNoIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20tightIDNoIso_[i]);

        h1_lep_LowptNLepGt20looseIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20looseIDNoIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20looseIDNoIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20looseIDNoIso_[i]);

        h1_lep_LowptNLepGt20vlooseIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20vlooseIDNoIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20vlooseIDNoIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20vlooseIDNoIso_[i]);

        h1_lep_LowptNLepGt20tightIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20tightIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20tightIDConvRIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20tightIDConvRIso0_1_[i]);

        h1_lep_LowptNLepGt20looseIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20looseIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20looseIDConvRIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20looseIDConvRIso0_1_[i]);

        h1_lep_LowptNLepGt20vlooseIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20vlooseIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20vlooseIDConvRIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20vlooseIDConvRIso0_1_[i]);

        h1_lep_LowptNLepGt20NOIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20NOIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20NOIDConvRIso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20NOIDConvRIso0_1_[i]);


//         h1_lep_LowptNLepGt20NOIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20NOIDIso0_1", hyp.c_str()), 
//                                                           "LowLep_NLepGt20NOIDIso0_1", 21, -0.5, 20.5);
// 	FormatHist(h1_lep_LowptNLepGt20NOIDIso0_1_[i]);

        h1_lep_LowptNLepGt20ID0Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID0Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID0Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID0Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID1Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID1Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID1Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID1Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID2Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID2Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID2Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID2Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID3Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID3Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID3Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID3Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID4Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID4Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID4Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID4Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID5Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID5Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID5Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID5Iso0_1_[i]);
        h1_lep_LowptNLepGt20ID6Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID6Iso0_1", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID6Iso0_1", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID6Iso0_1_[i]);

        h1_lep_LowptNLepGt20ID0NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID0NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID0NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID0NOIso_[i]);
        h1_lep_LowptNLepGt20ID1NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID1NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID1NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID1NOIso_[i]);
        h1_lep_LowptNLepGt20ID2NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID2NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID2NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID2NOIso_[i]);
        h1_lep_LowptNLepGt20ID3NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID3NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID3NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID3NOIso_[i]);
        h1_lep_LowptNLepGt20ID4NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID4NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID4NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID4NOIso_[i]);
        h1_lep_LowptNLepGt20ID5NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID5NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID5NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID5NOIso_[i]);
        h1_lep_LowptNLepGt20ID6NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_NLepGt20ID6NOIso", hyp.c_str()), 
                                                          "LowLep_NLepGt20ID6NOIso", 21, -0.5, 20.5);
	FormatHist(h1_lep_LowptNLepGt20ID6NOIso_[i]);


        // ibl paste begin;
        //wo ID/iso
	h1_lep_LowptmLepGt10Lt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt10Lt20", hyp.c_str()), 
                                                 "LowLep_mLepGt10Lt20", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt10Lt20_[i]);
	h1_lep_LowptmLepGt20_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20", hyp.c_str()), 
                                             "LowLep_mLepGt20", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20_[i]);
        //w ID/iso
        h1_lep_LowptmLepGt20NOIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20NOIDIso0_1", hyp.c_str()), 
                                                       "LowLep_mLepGt20NOIDIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20NOIDIso0_1_[i]);
        h1_lep_LowptmLepGt20ID0Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID0Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID0Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID0Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID1Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID1Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID1Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID1Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID2Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID2Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID2Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID2Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID3Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID3Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID3Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID3Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID4Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID4Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID4Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID4Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID5Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID5Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID5Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID5Iso0_1_[i]);
        h1_lep_LowptmLepGt20ID6Iso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID6Iso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID6Iso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID6Iso0_1_[i]);

        h1_lep_LowptmLepGt20ID0NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID0NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID0NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID0NOIso_[i]);
        h1_lep_LowptmLepGt20ID1NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID1NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID1NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID1NOIso_[i]);
        h1_lep_LowptmLepGt20ID2NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID2NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID2NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID2NOIso_[i]);
        h1_lep_LowptmLepGt20ID3NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID3NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID3NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID3NOIso_[i]);
        h1_lep_LowptmLepGt20ID4NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID4NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID4NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID4NOIso_[i]);
        h1_lep_LowptmLepGt20ID5NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID5NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID5NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID5NOIso_[i]);
        h1_lep_LowptmLepGt20ID6NOIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20ID6NOIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20ID6NOIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20ID6NOIso_[i]);

        // ibl paste end;
        //ibl paste resume
        h1_lep_LowptmLepGt20tightIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20tightIDIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20tightIDIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20tightIDIso0_1_[i]);

        h1_lep_LowptmLepGt20looseIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20looseIDIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20looseIDIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20looseIDIso0_1_[i]);

        h1_lep_LowptmLepGt20vlooseIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20vlooseIDIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20vlooseIDIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20vlooseIDIso0_1_[i]);

//         h1_lep_LowptmLepGt20NOIDIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20NOIDIso0_1", hyp.c_str()), 
//                                                           "LowLep_mLepGt20NOIDIso0_1", 200, 0., 200.);
// 	FormatHist(h1_lep_LowptmLepGt20NOIDIso0_1_[i]);

        h1_lep_LowptmLepGt20tightIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20tightIDNoIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20tightIDNoIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20tightIDNoIso_[i]);

        h1_lep_LowptmLepGt20looseIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20looseIDNoIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20looseIDNoIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20looseIDNoIso_[i]);

        h1_lep_LowptmLepGt20vlooseIDNoIso_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20vlooseIDNoIso", hyp.c_str()), 
                                                          "LowLep_mLepGt20vlooseIDNoIso", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20vlooseIDNoIso_[i]);

        h1_lep_LowptmLepGt20tightIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20tightIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20tightIDConvRIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20tightIDConvRIso0_1_[i]);

        h1_lep_LowptmLepGt20looseIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20looseIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20looseIDConvRIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20looseIDConvRIso0_1_[i]);

        h1_lep_LowptmLepGt20vlooseIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20vlooseIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20vlooseIDConvRIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20vlooseIDConvRIso0_1_[i]);

        h1_lep_LowptmLepGt20NOIDConvRIso0_1_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_mLepGt20NOIDConvRIso0_1", hyp.c_str()), 
                                                          "LowLep_mLepGt20NOIDConvRIso0_1", 200, 0., 200.);
	FormatHist(h1_lep_LowptmLepGt20NOIDConvRIso0_1_[i]);

        //ibl paste really end



	h1_lep_LowpthOverE_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_hOverE", hyp.c_str()), "LowLep_hOverE", 150, -0.1, .2);
	FormatHist( h1_lep_LowpthOverE_[i]);
	h1_lep_lowpteOverPIn_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_eOverPIn", hyp.c_str()), "LowLep_eOverPIn", 160, -0.1, 1.5);
	FormatHist( h1_lep_lowpteOverPIn_[i]);
	//	      h1_lep_lowpteSeedOverPOut_[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
	//	      h1_lep_lowpteSeedOverPIn_[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
	h1_lep_lowptfBrem_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_fBrem", hyp.c_str()), "LowLep_fBrem", 200, -2.0, 2.0);
	FormatHist(      h1_lep_lowptfBrem_[i]);
	h1_lep_lowptdEtaIn_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_dEtaIn", hyp.c_str()), "LowLep_dEtaIn", 200, -0.1, 0.1);
	FormatHist( h1_lep_lowptdEtaIn_[i]);
	h1_lep_lowptdEtaOut_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_dEtaOut", hyp.c_str()), "LowLep_dEtaOut", 200, -0.1, 0.1);
	FormatHist( h1_lep_lowptdEtaOut_[i]);
	h1_lep_lowptdPhiIn_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_dPhiIn", hyp.c_str()), "LowLep_dPhiIn", 100, -0.5, 0.5);
	FormatHist( h1_lep_lowptdPhiIn_[i]);
	h1_lep_lowptdPhiInPhiOut_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_dPhiInPhiOut", hyp.c_str()), "LowLep_dPhiInPhiOut", 100, -0.5, 0.5);
	FormatHist(  h1_lep_lowptdPhiInPhiOut_[i]);
	h1_lep_lowptdPhiOut_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_dPhiOut", hyp.c_str()), "LowLep_dPhiOut", 100, -0.5, 0.5);
	FormatHist( h1_lep_lowptdPhiOut_[i]);
        
	h1_lep_lowptsigmaPhiPhi_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_sigmaPhiPhi", hyp.c_str()), "LowLep_sigmaPhiPhi", 200, -0.2, 0.2);
	FormatHist( h1_lep_lowptsigmaPhiPhi_[i]);
	h1_lep_lowptsigmaIPhiIPhi_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_sigmaIPhiIPhi", hyp.c_str()), "LowLep_sigmaIPhiIPhi",  200, -0.2, 0.2);
	FormatHist( h1_lep_lowptsigmaIPhiIPhi_[i]);
	h1_lep_lowptsigmaEtaEta_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_sigmaEtaEta", hyp.c_str()), "LowLep_sigmaEtaEta", 200, -0.2, 0.2);
	FormatHist(  h1_lep_lowptsigmaEtaEta_[i]);
	h1_lep_lowptsigmaIEtaIEta_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_sigmaIEtaIEta", hyp.c_str()), "LowLep_sigmaIEtaIEta", 200, -0.2, 0.2);
	FormatHist(  h1_lep_lowptsigmaIEtaIEta_[i]);
        
	h1_lep_lowptegamma_robustLooseId_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_egammarobustLooseId", hyp.c_str()), "LowLep_egammarobustLooseId", 2, -0.5, 1.5);
	FormatHist(   h1_lep_lowptegamma_robustLooseId_[i]);
	h1_lep_lowptegamma_robustTightId_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_egammarobustTightId", hyp.c_str()), "LowLep_egammarobustTightId", 2, -0.5, 1.5);
	FormatHist(  h1_lep_lowptegamma_robustTightId_[i]);
	h1_lep_lowptegamma_looseId_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_egammalooseId", hyp.c_str()), "LowLep_egammalooseId", 2, -0.5, 1.5);
	FormatHist( h1_lep_lowptegamma_looseId_[i]);
	h1_lep_lowptegamma_tightId_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_egammatightId", hyp.c_str()), "LowLep_egammatightId", 2, -0.5, 1.5);
	FormatHist( h1_lep_lowptegamma_tightId_[i]);
	h1_lep_lowptegamma_robustHighEnergy_[i] = new TH1F( Form("%s_%s_%s", SampleName().c_str(), "LowLep_egammarobustHighEnergy", hyp.c_str()), "LowLep_egammarobustHighEnergy",  2, -0.5, 1.5);
	FormatHist(   h1_lep_lowptegamma_robustHighEnergy_[i]);

	double d0max = 0.1;
        
	NewHist( hlep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "lep_pt", hyp.c_str()), "lep_pt", 100, 0.0, 100.0);
	NewHist( hlep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "lep_transmass", hyp.c_str()), "lep_transmass", 200, 0.0, 200.0);
	NewHist( hlep_tcmet[i], Form("%s_%s_%s", SampleName().c_str(), "lep_tcmet", hyp.c_str()), "lep_met", 100, 0.0, 100.0);
	NewHist( hlep_clmumet[i], Form("%s_%s_%s", SampleName().c_str(), "lep_calomet_muon", hyp.c_str()), "lep_calomet_muon", 100, 0.0, 100.0);
	NewHist( hlep_met_dphi[i], Form("%s_%s_%s", SampleName().c_str(), "lep_met_dphi", hyp.c_str()), "lep_met_dphi", 100, 0, 2 * 3.14159);
	NewHist( hlep_trckIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_trckIso", hyp.c_str()), "lep_trckIso", 100, 0.0, 10.0);
	NewHist( hlep_ecalIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_ecalIso", hyp.c_str()), "lep_ecalIso", 100, 0.0, 10.0);
	NewHist( hlep_relIso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_relIso", hyp.c_str()), "lep_relIso", 100, 0.0, 1.0);
	NewHist( hlep_d0[i], Form("%s_%s_%s", SampleName().c_str(), "lep_d0", hyp.c_str()), "lep_d0", 100, 0.0, d0max);

	//for nlep, fill before cutting on it--same for W,Z
	NewHist( hlep_nlep_nod0iso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep_nod0iso", hyp.c_str()), "lep_nlep_nod0iso", 10, -0.5, 9.5);
	NewHist( hlep_nlep_nometiso[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep_nometiso", hyp.c_str()), "lep_nlep_nometiso", 10, -0.5, 9.5);
	NewHist( hlep_nlep[i], Form("%s_%s_%s", SampleName().c_str(), "lep_nlep", hyp.c_str()), "lep_nlep", 10, -0.5, 9.5);
	NewHist( hlep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet20", hyp.c_str()), "lep_njet20", 10, -0.5, 9.5);
	NewHist( hlep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "lep_njet30", hyp.c_str()), "lep_njet30", 10, -0.5, 9.5);
	NewHist( hlep_conv[i], Form("%s_%s_%s", SampleName().c_str(), "lep_conversions", hyp.c_str()), "lep_conversions", 2, -0.5, 1.5);

	//TH2F's for ABCD
	hlep_d0_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_trckIso", hyp.c_str()), "lep_d0_trckIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_ecalIso", hyp.c_str()), "lep_d0_ecalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_hcalIso", hyp.c_str()), "lep_d0_hcalIso", 100, 0.0, d0max, 100, 0.0, 10.0);
	hlep_d0_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_d0_relIso",  hyp.c_str()), "lep_d0_relIso", 100, 0.0, d0max, 100, 0.0, 1.0);
	
	hlep_met_trckIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_trckIso", hyp.c_str()), "lep_tcMet_trckIso", 100, 0.0, 50.0, 100, 0.0, 10.0);
	hlep_met_ecalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_ecalIso", hyp.c_str()), "lep_tcMet_ecalIso", 100, 0.0, 50.0, 100, 0.0, 10.0);
	hlep_met_hcalIso[i] = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_hcalIso", hyp.c_str()), "lep_tcMet_hcalIso", 100, 0.0, 50.0, 100, 0.0, 10.0);
	hlep_met_relIso[i]  = new TH2F( Form("%s_%s_%s", SampleName().c_str(), "lep_tcMet_relIso",  hyp.c_str()), "lep_tcMet_relIso",  100, 0.0, 50.0, 100, 0.0, 1.0);
  }

  // di-lepton histograms (three + 1 types)
  for (unsigned int i = 0; i < 4; ++i) {
	NewHist( hdilep_0_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_0_pt", dilepton_hypo_names[i]), "dilep_0_pt", 100, 0.0, 100.0);
	NewHist( hdilep_1_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_1_pt", dilepton_hypo_names[i]), "dilep_1_pt", 100, 0.0, 100.0);
	NewHist( hdilep_pt[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_pt", dilepton_hypo_names[i]), "dilep_pt", 100, 0.0, 100.0);
	NewHist( hdilep_mass[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_mass", dilepton_hypo_names[i]), "dilep_mass", 200, 0.0, 200.0);
	NewHist( hdilep_tcmet[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_tcmet", dilepton_hypo_names[i]), "dilep_tcmet", 100, 0.0, 100.0);
	NewHist( hdilep_clmumet[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_calomet_muon", dilepton_hypo_names[i]), "dilep_calomet_muon", 100, 0.0, 100.0);
	NewHist( hdilep_njet20[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet20", dilepton_hypo_names[i]), "dilep_njet20", 10, -0.5, 9.5);
	NewHist( hdilep_njet30[i], Form("%s_%s_%s", SampleName().c_str(), "dilep_njet30", dilepton_hypo_names[i]), "dilep_njet30", 10, -0.5, 9.5);

  }
	
  // event level histograms
  //NewHist( hdilep_nhyp, Form("%s_%s_%s", SampleName().c_str(), "dilep_nhyp", "all"), "dilep_nhyp", 10, -0.5, 9.5);

}

cuts_t Looper::DilepSelect() //(int i_hyp), no hyp, just idxs
{
  cuts_t ret = 0;
  float ptcut = 20.0;
  //int idx1 = (elidxs[0] != -1 ? elidxs[0] : muidxs[0]);
  //int idx2 = (elidxs[1] != -1 ? elidxs[1] : muidxs[1]);
  //int idx1, idx2;

  if( elidxs[0] != -1 && elidxs[1] != -1 ) {
	 
	//if (cms2.hyp_lt_p4()[i_hyp].pt() > ptcut && cms2.hyp_ll_p4()[i_hyp].pt() > ptcut)
	if( cms2.els_p4()[elidxs[0]].pt() >= ptcut && cms2.els_p4()[elidxs[1]].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	//if( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 )
	if( cms2.els_charge()[elidxs[0]] * cms2.els_charge()[elidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  else if( muidxs[0] != -1 && muidxs[1] != -1 ) {

	if( cms2.mus_p4()[muidxs[0]].pt() >= ptcut && cms2.mus_p4()[muidxs[1]].pt() >= ptcut)
	  ret |= CUT_BIT(DILEP_PT);

	if( cms2.mus_charge()[muidxs[0]] * cms2.mus_charge()[muidxs[1]] < 0 )
	  ret |= CUT_BIT(DILEP_OS);

	double mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
	if( mass > 76. && mass < 106. )
	  ret |= CUT_BIT(DILEP_MASS);

  }
  else
	cout << "BAD DILEPSELECT CALL" << endl;

  return ret;
}

cuts_t Looper::LepSelect(int lep_type, int i)
{
  cuts_t ret = 0;
  float ptcut = 20.0;
  LorentzVector lep_p4;
  
  // e
  if( lep_type == 0 ) {

	lep_p4 = cms2.els_p4()[i];

	if (cms2.els_p4()[i].pt() >= ptcut)
	  ret |= CUT_BIT(LEP_PT);

	if( GoodSusyElectronWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyElectronIsolation(i, true) ) //bool is for use calo
	  ret |= CUT_BIT(LEP_ISO);

	//put in all cuts from GoodSusyElectronWithoutIsolation
	//if ( cms2.els_egamma_tightId().at(index)     !=  1) return false; 
	//if ( fabs(cms2.els_d0corr().at(index)) >= 0.02)   return false;
	//if ( cms2.els_closestMuon().at(index) != -1) return false; 
	//if ( TMath::Abs(cms2.els_p4()[index].eta()) > 2.4) return false;

	if ( cms2.els_egamma_tightId().at(i) ==  1
		 && cms2.els_closestMuon().at(i) == -1
		 && TMath::Abs(cms2.els_p4()[i].eta()) < 2.4 )
	  ret |= CUT_BIT(LEP_GOOD_NOD0);

	if ( fabs(cms2.els_d0corr().at(i)) <= 0.02)
	  ret |= CUT_BIT(LEP_D0);

  }

  // m
  else if( lep_type == 1 ) {

	lep_p4 = cms2.mus_p4()[i];

	if (cms2.mus_p4()[i].pt() >= ptcut)
	  ret |= CUT_BIT(LEP_PT);

	if( GoodSusyMuonWithoutIsolation(i) )
	  ret |= CUT_BIT(LEP_GOOD);

	if( PassSusyMuonIsolation(i) )
	  ret |= CUT_BIT(LEP_ISO);

	//put in all cuts from GoodSusyMuonWithoutIsolation
	//if (((cms2.mus_type().at(index)) & (1<<1)) == 0) return false; // global muon
	//if (((cms2.mus_type().at(index)) & (1<<2)) == 0) return false; // tracker muon
	//if (cms2.mus_validHits().at(index) < 11)    return false; 
	//if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) return false; 
	//if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) return false;
	//if (cms2.mus_pat_ecalvetoDep().at(index) >= 4) return false; // ECalE < 4 
	//if (cms2.mus_pat_hcalvetoDep().at(index) >= 6) return false; // HCalE < 6 
	//if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) return false;
	
	if( (cms2.mus_type().at(i) & (1<<1)) == (1<<1) // global muon
		&& (cms2.mus_type().at(i) & (1<<2)) == (1<<2) // tracker muon
		&& cms2.mus_validHits().at(i) >= 11
		&& cms2.mus_gfit_chi2().at(i)/cms2.mus_gfit_ndof().at(i) < 10
		&& cms2.mus_pat_ecalvetoDep().at(i) < 4 // ECalE < 4 
		&& cms2.mus_pat_hcalvetoDep().at(i) < 6 // HCalE < 6 
		&& TMath::Abs(cms2.mus_p4()[i].eta()) < 2.4)
	  ret |= CUT_BIT(LEP_GOOD_NOD0);
	
	if (fabs(cms2.mus_d0corr().at(i)) <= 0.02) 
	  ret |= CUT_BIT(LEP_D0);
	
  }

  //transverse mass
  double dphi = ROOT::Math::VectorUtil::DeltaPhi( LorentzVector( cms2.evt_tcmet()*cos(cms2.evt_tcmetPhi()),
																 cms2.evt_tcmet()*sin(cms2.evt_tcmetPhi()),
																 0, cms2.evt_tcmet())
												  , lep_p4 );
  double masst = sqrt( 2 * cms2.evt_tcmet() * lep_p4.Et() * ( 1 - cos(dphi) ) );
  
  //check yanjun's code:
  TVector3 tcMET;
  tcMET.SetPtEtaPhi(cms2.evt_tcmet(), 0, cms2.evt_tcmetPhi());
  double massyj = sqrt( ( tcMET.Pt() + lep_p4.Et())*( tcMET.Pt() + lep_p4.Et())
					   - ( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )*( tcMET.Pt()*cos(tcMET.Phi()) + lep_p4.Px() )
					   - ( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() )*( tcMET.Pt()*sin(tcMET.Phi()) + lep_p4.Py() ) );

  if( fabs( tcMET.Phi() - cms2.evt_tcmetPhi() ) > 0.01 )
	cout << "Phi error " << tcMET.Phi() << "   " << cms2.evt_tcmetPhi() << endl;

  if( fabs( massyj - masst ) > 0.1 && masst > 2 && massyj > 2 )
	cout << "Mass disagreement " << masst << "   " << massyj << endl;

  //set looper member
  transmass = masst;

  if( masst > 40 && masst < 100 )
	ret |= CUT_BIT(TMASS);


  return ret;
}

void Looper::FillEventHistos ()
{

  // a first loop over the most energetic 
  // muon and electron to determine a good cut 
  // in pT and iso
  int   hiPtIdx    = -1;
  float hiPtmax    = -1.;

  // get the event weight
  if( sample_.kFactor != 1 ) cout << "kFactor non-unity " << sample_.kFactor << endl;
  double weight = Looper::Weight();
  
  // have a look at the highest pt electron
  for(int ele = 0; ele < (int)cms2.els_p4().size(); ++ele) {
    if(cms2.els_p4()[ele].pt() > hiPtmax && GoodSusyElectronWithoutIsolation(ele) ) hiPtIdx = ele;
  }

  // histogram indices are e, m, all (0, 1, 2)
  if(hiPtIdx != -1)  {
    h1_lep_Highpt_[0]                                                    ->Fill(cms2.els_p4()[hiPtIdx].pt(), weight);
    h1_lep_HighptMet_[0]                                                 ->Fill(cms2.evt_tcmet(), weight);
    h1_lep_HighptRelIso_[0]                                              ->Fill(inv_el_relsusy_iso(hiPtIdx, true), weight);
    if(cms2.els_p4()[hiPtIdx].pt() > 20. ) h1_lep_HighptRelIsoPtLg20_[0] ->Fill(inv_el_relsusy_iso(hiPtIdx, true), weight);
  }

 // have a look at all but the highest pt electron
  uint nEleGt10Lt20 = 0;
  uint nEleGt20     = 0;
  uint nEleGt20tightIDIso0_1     = 0;
  uint nEleGt20looseIDIso0_1     = 0;
  uint nEleGt20vlooseIDIso0_1     = 0;
  uint nEleGt20NOIDIso0_1     = 0;
  uint nEleGt20tightIDNoIso     = 0;
  uint nEleGt20looseIDNoIso     = 0;
  uint nEleGt20vlooseIDNoIso     = 0;

  uint nEleGt20tightIDConvRIso0_1     = 0;
  uint nEleGt20looseIDConvRIso0_1     = 0;
  uint nEleGt20vlooseIDConvRIso0_1     = 0;
  uint nEleGt20NOIDConvRIso0_1     = 0;

  double iblmass = -1.0;

  for(int ele = 0; ele <  (int)cms2.els_p4().size(); ++ele) {
    if(hiPtIdx != -1 && hiPtIdx != ele) {

      iblmass = (cms2.els_p4()[hiPtIdx] + cms2.els_p4()[ele]).mass();

      h1_lep_Lowpt_[0]                                                ->Fill(cms2.els_p4()[ele].pt(), weight);
      h1_lep_LowptMet_[0]                                             ->Fill(cms2.evt_tcmet(), weight);
      h1_lep_LowptRelIso_[0]                                          ->Fill(inv_el_relsusy_iso(ele, true), weight);
      if(cms2.els_p4()[ele].pt() > 20. ) h1_lep_LowptRelIsoPtLg20_[0] ->Fill(inv_el_relsusy_iso(ele, true), weight); // was buggy until 090928

      if(cms2.els_p4()[ele].pt() > 10. && 
         (cms2.els_p4()[ele].pt() < 20. )) {       
        ++nEleGt10Lt20;
        h1_lep_LowptmLepGt10Lt20_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. ) {
        ++nEleGt20;
        h1_lep_LowptmLepGt20_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_tightId()[ele] && 
         inv_el_relsusy_iso(ele, true) < 0.1 ) {
        ++nEleGt20tightIDIso0_1;
        h1_lep_LowptmLepGt20tightIDIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_looseId()[ele] && 
         inv_el_relsusy_iso(ele, true) < 0.1 ) {
        ++nEleGt20looseIDIso0_1;
        h1_lep_LowptmLepGt20looseIDIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt()        > 20.  && 
         fabs(cms2.els_dEtaOut()[ele])  < 0.02 && 
         cms2.els_eOverPIn()[ele]       > 0.4  &&  
         fabs(cms2.els_hOverE()[ele])   < 0.03 && 
         inv_el_relsusy_iso(ele, true)  < 0.1 ) {
        ++nEleGt20vlooseIDIso0_1;
        h1_lep_LowptmLepGt20vlooseIDIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt()        > 20.  && 
         inv_el_relsusy_iso(ele, true)  < 0.1 ) {
        ++nEleGt20NOIDIso0_1;
        h1_lep_LowptmLepGt20NOIDIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_tightId()[ele] ) {
        ++nEleGt20tightIDNoIso;
        h1_lep_LowptmLepGt20tightIDNoIso_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_looseId()[ele] ) {
        ++nEleGt20looseIDNoIso;
        h1_lep_LowptmLepGt20looseIDNoIso_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt()        > 20.  && 
         fabs(cms2.els_dEtaOut()[ele])  < 0.02 && 
         cms2.els_eOverPIn()[ele]       > 0.4  &&  
         fabs(cms2.els_hOverE()[ele])   < 0.03 ) {
        ++nEleGt20vlooseIDNoIso;
        h1_lep_LowptmLepGt20vlooseIDNoIso_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_tightId()[ele] &&
         !conversionElectron(ele) && 
         inv_el_relsusy_iso(ele, true) < 0.1 ) {
        ++nEleGt20tightIDConvRIso0_1;
        h1_lep_LowptmLepGt20tightIDConvRIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt() > 20. && 
         cms2.els_egamma_looseId()[ele] && 
         !conversionElectron(ele) && 
         inv_el_relsusy_iso(ele, true) < 0.1 ) {
        ++nEleGt20looseIDConvRIso0_1;
        h1_lep_LowptmLepGt20looseIDConvRIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt()        > 20.  && 
         fabs(cms2.els_dEtaOut()[ele])  < 0.02 && 
         cms2.els_eOverPIn()[ele]       > 0.4  &&  
         fabs(cms2.els_hOverE()[ele])   < 0.03 && 
         !conversionElectron(ele)              &&
         inv_el_relsusy_iso(ele, true)  < 0.1 ) {
        ++nEleGt20vlooseIDConvRIso0_1;
        h1_lep_LowptmLepGt20vlooseIDConvRIso0_1_[0]->Fill(iblmass, weight); 
      }
      if(cms2.els_p4()[ele].pt()        > 20.  && 
         !conversionElectron(ele)              &&
         inv_el_relsusy_iso(ele, true)  < 0.1 ) {
        ++nEleGt20NOIDConvRIso0_1;
        h1_lep_LowptmLepGt20NOIDConvRIso0_1_[0]->Fill(iblmass, weight); 
      }

      h1_lep_LowpthOverE_[0]->Fill(cms2.els_hOverE()[ele], weight);
      h1_lep_lowpteOverPIn_[0]->Fill(cms2.els_eOverPIn()[ele], weight);
//       h1_lep_lowpteSeedOverPOut_[0]->Fill(cms2.els_eSeedOverPOut()[ele], weight);
//       h1_lep_lowpteSeedOverPIn_[0]->Fill(cms2.els_eSeedOverPIn()[ele], weight);
      h1_lep_lowptfBrem_[0]->Fill(cms2.els_fBrem()[ele], weight);
      h1_lep_lowptdEtaIn_[0]->Fill(cms2.els_dEtaIn()[ele], weight);
      h1_lep_lowptdEtaOut_[0]->Fill(cms2.els_dEtaOut()[ele], weight);
      h1_lep_lowptdPhiIn_[0]->Fill(cms2.els_dPhiIn()[ele], weight);
      h1_lep_lowptdPhiInPhiOut_[0]->Fill(cms2.els_dPhiInPhiOut()[ele], weight);
      h1_lep_lowptdPhiOut_[0]->Fill(cms2.els_dPhiOut()[ele], weight);

      h1_lep_lowptsigmaPhiPhi_[0]->Fill(cms2.els_sigmaPhiPhi()[ele], weight);  
      h1_lep_lowptsigmaIPhiIPhi_[0]->Fill(cms2.els_sigmaIPhiIPhi()[ele], weight);  
      h1_lep_lowptsigmaEtaEta_[0]->Fill(cms2.els_sigmaEtaEta()[ele], weight);  
      h1_lep_lowptsigmaIEtaIEta_[0]->Fill(cms2.els_sigmaIEtaIEta()[ele], weight);  

      h1_lep_lowptegamma_robustLooseId_[0]->Fill(cms2.els_egamma_robustLooseId()[ele], weight);  
      h1_lep_lowptegamma_robustTightId_[0]->Fill(cms2.els_egamma_robustTightId()[ele], weight);  
      h1_lep_lowptegamma_looseId_[0]->Fill(cms2.els_egamma_looseId()[ele], weight);  
      h1_lep_lowptegamma_tightId_[0]->Fill(cms2.els_egamma_tightId()[ele], weight);  
      h1_lep_lowptegamma_robustHighEnergy_[0]->Fill(cms2.els_egamma_robustHighEnergy()[ele], weight);  
 
    }
  }
  h1_lep_LowptNLepGt10Lt20_[0]->Fill(nEleGt10Lt20, weight);
  h1_lep_LowptNLepGt20_[0]    ->Fill(nEleGt20, weight);
  h1_lep_LowptNLepGt20tightIDIso0_1_[0]    ->Fill(nEleGt20tightIDIso0_1, weight);
  h1_lep_LowptNLepGt20looseIDIso0_1_[0]    ->Fill(nEleGt20looseIDIso0_1, weight);
  h1_lep_LowptNLepGt20vlooseIDIso0_1_[0]   ->Fill(nEleGt20vlooseIDIso0_1, weight);
  h1_lep_LowptNLepGt20NOIDIso0_1_[0]       ->Fill(nEleGt20NOIDIso0_1, weight);
  h1_lep_LowptNLepGt20tightIDNoIso_[0]    ->Fill(nEleGt20tightIDNoIso, weight);
  h1_lep_LowptNLepGt20looseIDNoIso_[0]    ->Fill(nEleGt20looseIDNoIso, weight);
  h1_lep_LowptNLepGt20vlooseIDNoIso_[0]   ->Fill(nEleGt20vlooseIDNoIso, weight);

  h1_lep_LowptNLepGt20tightIDConvRIso0_1_[0]    ->Fill(nEleGt20tightIDConvRIso0_1, weight);
  h1_lep_LowptNLepGt20looseIDConvRIso0_1_[0]    ->Fill(nEleGt20looseIDConvRIso0_1, weight);
  h1_lep_LowptNLepGt20vlooseIDConvRIso0_1_[0]   ->Fill(nEleGt20vlooseIDConvRIso0_1, weight);
  h1_lep_LowptNLepGt20NOIDConvRIso0_1_[0]       ->Fill(nEleGt20NOIDConvRIso0_1, weight);

  // have a look at the highest pt muon
  // reset highpt index
  hiPtIdx    = -1;
  hiPtmax    = -1.;
  for(int muo = 0; muo <  (int)cms2.mus_p4().size(); ++muo) {
    if(cms2.mus_p4()[muo].pt() > hiPtmax && GoodSusyMuonWithoutIsolation(muo) ) hiPtIdx = muo;
  }
  
  // histogram indices are e, m, all (0, 1, 2)
  if(hiPtIdx != -1)  {
    h1_lep_Highpt_[1]                                                    ->Fill(cms2.mus_p4()[hiPtIdx].pt(), weight);
    h1_lep_HighptMet_[1]                                                 ->Fill(cms2.evt_tcmet(), weight);
    h1_lep_HighptRelIso_[1]                                              ->Fill(inv_mu_relsusy_iso(hiPtIdx), weight);
    if(cms2.mus_p4()[hiPtIdx].pt() > 20. ) h1_lep_HighptRelIsoPtLg20_[1] ->Fill(inv_mu_relsusy_iso(hiPtIdx), weight);
  }

  // have a look at all but the highest pt muoctron ;)
  uint nMuoGt10Lt20 = 0;
  uint nMuoGt20     = 0;

  uint nMuoGt20NOIDIso0_1     = 0;

  uint nMuoGt20ID0Iso0_1     = 0;
  uint nMuoGt20ID1Iso0_1     = 0;
  uint nMuoGt20ID2Iso0_1     = 0;
  uint nMuoGt20ID3Iso0_1     = 0;
  uint nMuoGt20ID4Iso0_1     = 0;
  uint nMuoGt20ID5Iso0_1     = 0;
  uint nMuoGt20ID6Iso0_1     = 0;

  uint nMuoGt20ID0NOIso     = 0;
  uint nMuoGt20ID1NOIso     = 0;
  uint nMuoGt20ID2NOIso     = 0;
  uint nMuoGt20ID3NOIso     = 0;
  uint nMuoGt20ID4NOIso     = 0;
  uint nMuoGt20ID5NOIso     = 0;
  uint nMuoGt20ID6NOIso     = 0;

  iblmass = -1.0;

  for(int muo = 0; muo <  (int)cms2.mus_p4().size(); ++muo) {
    if(hiPtIdx != -1 && hiPtIdx != muo) {
      
      iblmass = (cms2.mus_p4()[hiPtIdx] + cms2.mus_p4()[muo]).mass();

      h1_lep_Lowpt_[1]                                                ->Fill(cms2.mus_p4()[muo].pt(), weight);
      h1_lep_LowptMet_[1]                                             ->Fill(cms2.evt_tcmet(), weight);
      h1_lep_LowptRelIso_[1]                                          ->Fill(inv_mu_relsusy_iso(muo), weight);
      if(cms2.mus_p4()[muo].pt() > 20. ) h1_lep_LowptRelIsoPtLg20_[1] ->Fill(inv_mu_relsusy_iso(muo), weight);
      if(cms2.mus_p4()[muo].pt() > 10. && 
         (cms2.mus_p4()[muo].pt() < 20. ))       {
        ++nMuoGt10Lt20;
        h1_lep_LowptmLepGt10Lt20_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. ) {
        ++nMuoGt20;
        h1_lep_LowptmLepGt20_[1]->Fill(iblmass, weight); 
      }

      // GoodTMTestMuonWithoutIsolation(int index, int mode)
      //   //modes:
      // mode = 1: TMLS tight, d0, ecal, hcal, eta
      // mode = 2: TMLS tight, d0, eta
      // mode = 3: TMLS loose, d0, ecal, hcal, eta
      // mode = 4: TMLS loose, d0, eta
      // mode = 5: TMC2D tight, d0, eta
      // mode = 6: TMC2D loose, d0, eta

      // require any Global or Tracker muon
      bool isGlobalMu  = true;
      bool isTrackerMu = true;
      if (((cms2.mus_type().at(muo)) & (1<<1)) == 0) isGlobalMu  = false; // global muon
      if (((cms2.mus_type().at(muo)) & (1<<2)) == 0) isTrackerMu = false; // tracker muon

      if(cms2.mus_p4()[muo].pt() > 20. && 
         ( isGlobalMu || isTrackerMu ) &&
         inv_mu_relsusy_iso(muo) < 0.1 ) {        
        ++nMuoGt20NOIDIso0_1;
        h1_lep_LowptmLepGt20NOIDIso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodSusyMuonWithoutIsolation(muo) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID0Iso0_1;
        h1_lep_LowptmLepGt20ID0Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 1) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID1Iso0_1;
      h1_lep_LowptmLepGt20ID1Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 2) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID2Iso0_1;
      h1_lep_LowptmLepGt20ID2Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 3) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID3Iso0_1;
      h1_lep_LowptmLepGt20ID3Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 4) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID4Iso0_1;
      h1_lep_LowptmLepGt20ID4Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 5) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID5Iso0_1;
      h1_lep_LowptmLepGt20ID5Iso0_1_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 6) && 
         inv_mu_relsusy_iso(muo) < 0.1 ) {
        ++nMuoGt20ID6Iso0_1;
      h1_lep_LowptmLepGt20ID6Iso0_1_[1]->Fill(iblmass, weight); 
      }

      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodSusyMuonWithoutIsolation(muo) ) {
        ++nMuoGt20ID0NOIso;
        h1_lep_LowptmLepGt20ID0NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 1) ) {
        ++nMuoGt20ID1NOIso;
        h1_lep_LowptmLepGt20ID1NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 2) ) {
        ++nMuoGt20ID2NOIso;
        h1_lep_LowptmLepGt20ID2NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 3) ) {
        ++nMuoGt20ID3NOIso;
        h1_lep_LowptmLepGt20ID3NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 4) ) {
        ++nMuoGt20ID4NOIso;
        h1_lep_LowptmLepGt20ID4NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 5) ) {
        ++nMuoGt20ID5NOIso;
        h1_lep_LowptmLepGt20ID5NOIso_[1]->Fill(iblmass, weight); 
      }
      if(cms2.mus_p4()[muo].pt() > 20. && 
         GoodTMTestMuonWithoutIsolation(muo, 6) ) {
        ++nMuoGt20ID6NOIso;
        h1_lep_LowptmLepGt20ID6NOIso_[1]->Fill(iblmass, weight); 
      }
     }
  }
  h1_lep_LowptNLepGt10Lt20_[1]->Fill(nMuoGt10Lt20, weight);
  h1_lep_LowptNLepGt20_[1]    ->Fill(nMuoGt20, weight);

  h1_lep_LowptNLepGt20NOIDIso0_1_[1]    ->Fill(nMuoGt20NOIDIso0_1, weight);

  h1_lep_LowptNLepGt20ID0Iso0_1_[1]    ->Fill(nMuoGt20ID0Iso0_1, weight);
  h1_lep_LowptNLepGt20ID1Iso0_1_[1]    ->Fill(nMuoGt20ID1Iso0_1, weight);
  h1_lep_LowptNLepGt20ID2Iso0_1_[1]    ->Fill(nMuoGt20ID2Iso0_1, weight);
  h1_lep_LowptNLepGt20ID3Iso0_1_[1]    ->Fill(nMuoGt20ID3Iso0_1, weight);
  h1_lep_LowptNLepGt20ID4Iso0_1_[1]    ->Fill(nMuoGt20ID4Iso0_1, weight);
  h1_lep_LowptNLepGt20ID5Iso0_1_[1]    ->Fill(nMuoGt20ID5Iso0_1, weight);
  h1_lep_LowptNLepGt20ID6Iso0_1_[1]    ->Fill(nMuoGt20ID6Iso0_1, weight);

  h1_lep_LowptNLepGt20ID0NOIso_[1]    ->Fill(nMuoGt20ID0NOIso, weight);
  h1_lep_LowptNLepGt20ID1NOIso_[1]    ->Fill(nMuoGt20ID1NOIso, weight);
  h1_lep_LowptNLepGt20ID2NOIso_[1]    ->Fill(nMuoGt20ID2NOIso, weight);
  h1_lep_LowptNLepGt20ID3NOIso_[1]    ->Fill(nMuoGt20ID3NOIso, weight);
  h1_lep_LowptNLepGt20ID4NOIso_[1]    ->Fill(nMuoGt20ID4NOIso, weight);
  h1_lep_LowptNLepGt20ID5NOIso_[1]    ->Fill(nMuoGt20ID5NOIso, weight);
  h1_lep_LowptNLepGt20ID6NOIso_[1]    ->Fill(nMuoGt20ID6NOIso, weight);
  
  // need to determine if this is a di-lepton
  // or a single lepton event
  int nels = 0, nmus = 0;
  int nels_nod0iso = 0, nels_nometiso = 0;
  int nmus_nod0iso = 0, nmus_nometiso = 0;
  int nels_nopt = 0, nmus_nopt = 0;
  cuts_t pass_lep_cut = 0; //cuts of lepton which passes lepcuts--only good for single lepton
  elidxs[0] = elidxs[1] = -1;
  muidxs[0] = muidxs[1] = -1;
  int elidxs_nopt[] = {-1, -1}; //this one can be local--still just 2
  int muidxs_nopt[] = {-1, -1};
  //DEFINE CUTS
  //cut for counting good, isolated leptons: no mass, met reqmnts
  cuts_t lepcuts = (CUT_BIT(LEP_PT)
					| CUT_BIT(LEP_GOOD)
					| CUT_BIT(LEP_ISO)
					);

  cuts_t w_evt_sel = ( CUT_BIT(TCMET) | CUT_BIT(TMASS) );
  //cut for pt plot, also no mass, met cuts
  cuts_t lepcuts_nopt  = ( CUT_BIT(LEP_GOOD) | CUT_BIT(LEP_ISO) );
  //for (d0,iso) plots, relax (d0,iso), apply everything else including tmass, tcmet
  cuts_t lepcuts_nod0  = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD_NOD0) | w_evt_sel );
  cuts_t lepcuts_noiso = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | w_evt_sel );
  cuts_t lepcuts_notm  = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_ISO) | CUT_BIT(LEP_GOOD) | CUT_BIT(TCMET) );

  //cuts for ABCD: CUT ON TMASS, MET
  cuts_t lepcuts_nod0iso = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD_NOD0) | w_evt_sel );
  cuts_t lepcuts_nometiso = ( CUT_BIT(LEP_PT) | CUT_BIT(LEP_GOOD) | CUT_BIT(TMASS) );

  //check tcMET here b'c no point doing for every lepton--don't use for Z
  cuts_t tcmetcut = 0;
  if( cms2.evt_tcmet() > 20 )
	tcmetcut |= CUT_BIT(TCMET);

  //select els
  for( unsigned int i=0; i<cms2.els_p4().size(); i++ ) {
	transmass = 0;
	cuts_t elcut = LepSelect(0, i); //0 for els
	elcut |= tcmetcut;
	if( (elcut & lepcuts) == lepcuts ) { //all cuts
	  nels++;
	  pass_lep_cut = elcut;
	  if( elidxs[0] == -1 )
		elidxs[0] = i;
	  else if( elidxs[1] == -1 )
		elidxs[1] = i;
	  //if > 2 els, ignore the rest--should be sorted by pt already
	}
	if( (elcut & lepcuts_nopt) == lepcuts_nopt ) { //no pt cut
	  nels_nopt++;
	  if( elidxs_nopt[0] == -1 )
		elidxs_nopt[0] = i;
	  else if( elidxs_nopt[1] == -1 )
		elidxs_nopt[1] = i;
	}
	if( (elcut & lepcuts_noiso) == lepcuts_noiso ) {
	  hlep_trckIso[0]->Fill( cms2.els_tkIso()[i], weight );
	  hlep_trckIso[2]->Fill( cms2.els_tkIso()[i], weight );
	  hlep_ecalIso[0]->Fill( cms2.els_pat_ecalIso()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.els_pat_ecalIso()[i], weight );
	  hlep_relIso[0]->Fill( inv_el_relsusy_iso(i, true), weight );
	  hlep_relIso[2]->Fill( inv_el_relsusy_iso(i, true), weight );
	}
	if( (elcut & lepcuts_nod0) == lepcuts_nod0 ) {
	  hlep_d0[0]->Fill( cms2.els_d0corr()[i], weight );
	  hlep_d0[2]->Fill( cms2.els_d0corr()[i], weight );
	}
	//for mass plots, check all but mass
	if( (elcut & lepcuts_notm) == lepcuts_notm ) {
	  hlep_mass[0]->Fill( transmass, weight ); //check all cuts but mass for mass plot (n-1)
	  hlep_mass[2]->Fill( transmass, weight );
	}
	//ABCD in d0, Iso
	if( (elcut & lepcuts_nod0iso) == lepcuts_nod0iso ) {
	  nels_nod0iso++;
	  hlep_d0_trckIso[0]->Fill( cms2.els_d0corr()[i], cms2.els_tkIso()[i], weight );
	  hlep_d0_trckIso[2]->Fill( cms2.els_d0corr()[i], cms2.els_tkIso()[i], weight );
	  hlep_d0_ecalIso[0]->Fill( cms2.els_d0corr()[i], cms2.els_pat_ecalIso()[i], weight );
	  hlep_d0_ecalIso[2]->Fill( cms2.els_d0corr()[i], cms2.els_pat_ecalIso()[i], weight );
	  hlep_d0_hcalIso[0]->Fill( cms2.els_d0corr()[i], cms2.els_pat_hcalIso()[i], weight );
	  hlep_d0_hcalIso[2]->Fill( cms2.els_d0corr()[i], cms2.els_pat_hcalIso()[i], weight );
	  hlep_d0_relIso[0] ->Fill( cms2.els_d0corr()[i], inv_el_relsusy_iso(i, true), weight );
	  hlep_d0_relIso[2] ->Fill( cms2.els_d0corr()[i], inv_el_relsusy_iso(i, true), weight );
	}
	//ABCD in tcmet, Iso
	if( (elcut & lepcuts_nometiso) == lepcuts_nometiso ) {
	  nels_nometiso++;
	  hlep_met_trckIso[0]->Fill( cms2.evt_tcmet(), cms2.els_tkIso()[i], weight );		
	  hlep_met_trckIso[2]->Fill( cms2.evt_tcmet(), cms2.els_tkIso()[i], weight );		
	  hlep_met_ecalIso[0]->Fill( cms2.evt_tcmet(), cms2.els_pat_ecalIso()[i], weight );	
	  hlep_met_ecalIso[2]->Fill( cms2.evt_tcmet(), cms2.els_pat_ecalIso()[i], weight );	
	  hlep_met_hcalIso[0]->Fill( cms2.evt_tcmet(), cms2.els_pat_hcalIso()[i], weight );	
	  hlep_met_hcalIso[2]->Fill( cms2.evt_tcmet(), cms2.els_pat_hcalIso()[i], weight );	
	  hlep_met_relIso[0] ->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(i, true), weight );
	  hlep_met_relIso[2] ->Fill( cms2.evt_tcmet(), inv_el_relsusy_iso(i, true), weight );
	}
  }

  //select mus
  for( unsigned int i=0; i<cms2.mus_p4().size(); i++ ) {
	transmass = 0;
	cuts_t mucut = LepSelect(1, i); //1 for mus
	mucut |= tcmetcut;
	if( (mucut & lepcuts) == lepcuts ) { //all cuts
	  nmus++;
	  pass_lep_cut = mucut;
	  if( muidxs[0] == -1 )
		muidxs[0] = i;
	  else if( muidxs[1] == -1 )
		muidxs[1] = i;
	}
	if( (mucut & lepcuts_nopt) == lepcuts_nopt ) { //no pt cut
	  nmus_nopt++;
	  if( muidxs_nopt[0] == -1 )
		muidxs_nopt[0] = i;
	  else if( muidxs_nopt[1] == -1 )
		muidxs_nopt[1] = i;
	}
	if( (mucut & lepcuts_noiso) == lepcuts_noiso ) {
	  hlep_trckIso[1]->Fill( cms2.mus_iso03_sumPt()[i], weight);
	  hlep_trckIso[2]->Fill( cms2.mus_iso03_sumPt()[i], weight);
	  hlep_ecalIso[1]->Fill( cms2.mus_iso03_emEt()[i], weight );
	  hlep_ecalIso[2]->Fill( cms2.mus_iso03_emEt()[i], weight );
	  hlep_relIso[1]->Fill( inv_mu_relsusy_iso(i), weight );
	  hlep_relIso[2]->Fill( inv_mu_relsusy_iso(i), weight );
	}
	if( (mucut & lepcuts_nod0) == lepcuts_nod0 ) {
	  hlep_d0[1]->Fill( cms2.mus_d0corr()[i], weight );
	  hlep_d0[2]->Fill( cms2.mus_d0corr()[i], weight );
	}
	//for mass plots, check all but mass
	if( (mucut & lepcuts_notm) == lepcuts_notm ) {
	  hlep_mass[1]->Fill( transmass, weight ); //check all cuts but mass for mass plot (n-1)
	  hlep_mass[2]->Fill( transmass, weight );
	}
	//ABCD in d0, Iso
	if( (mucut & lepcuts_nod0iso) == lepcuts_nod0iso ) {
	  nmus_nod0iso++;
	  hlep_d0_trckIso[1]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_sumPt()[i], weight );
	  hlep_d0_trckIso[2]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_sumPt()[i], weight );
	  hlep_d0_ecalIso[1]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_emEt()[i], weight );
	  hlep_d0_ecalIso[2]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_emEt()[i], weight );
	  hlep_d0_hcalIso[1]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_hadEt()[i], weight );
	  hlep_d0_hcalIso[2]->Fill( cms2.mus_d0corr()[i], cms2.mus_iso03_hadEt()[i], weight );
	  hlep_d0_relIso[1] ->Fill( cms2.mus_d0corr()[i], inv_mu_relsusy_iso(i), weight );
	  hlep_d0_relIso[2] ->Fill( cms2.mus_d0corr()[i], inv_mu_relsusy_iso(i), weight );
	}
	//ABCD in tcmet, Iso
	if( (mucut & lepcuts_nometiso) == lepcuts_nometiso ) {
	  nmus_nometiso++;
	  hlep_met_trckIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_sumPt()[i], weight );		
	  hlep_met_trckIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_sumPt()[i], weight );		
	  hlep_met_ecalIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_emEt()[i], weight );	
	  hlep_met_ecalIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_emEt()[i], weight );	
	  hlep_met_hcalIso[1]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_hadEt()[i], weight );	
	  hlep_met_hcalIso[2]->Fill( cms2.evt_tcmet(), cms2.mus_iso03_hadEt()[i], weight );	
	  hlep_met_relIso[1] ->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(i), weight );
	  hlep_met_relIso[2] ->Fill( cms2.evt_tcmet(), inv_mu_relsusy_iso(i), weight );
	}
  }

  //fill nlep cut before checking nels, nmus
  //0=el, 1=mu, 2=all
  hlep_nlep[0]->Fill( nels, weight );
  hlep_nlep[1]->Fill( nmus, weight );
  hlep_nlep[2]->Fill( nels+nmus, weight );

  hlep_nlep_nod0iso[0]->Fill( nels_nod0iso, weight );
  hlep_nlep_nod0iso[1]->Fill( nmus_nod0iso, weight );
  hlep_nlep_nod0iso[2]->Fill( nels_nod0iso + nmus_nod0iso, weight );
  hlep_nlep_nometiso[0]->Fill( nels_nometiso, weight );
  hlep_nlep_nometiso[1]->Fill( nmus_nometiso, weight );
  hlep_nlep_nometiso[2]->Fill( nels_nometiso + nmus_nometiso, weight );

  //fill lep pt hists based on n nopt--single leps (each flavor independent of other)
  if( nels_nopt == 1 ) { //allow mus
	hlep_pt[0]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hlep_pt[2]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
  }
  if( nmus_nopt == 1 ) {
	hlep_pt[1]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hlep_pt[2]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
  }
  //now for N leps  //DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE
  if( nels_nopt > 1 ) { //allow mus
	hdilep_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_0_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_0_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[0]].pt(), weight );
	hdilep_1_pt[DILEPTON_EE]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
	hdilep_1_pt[DILEPTON_ALL]->Fill( cms2.els_p4()[elidxs_nopt[1]].pt(), weight );
  }
  if( nmus_nopt > 1 ) {
	hdilep_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_0_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_0_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[0]].pt(), weight );
	hdilep_1_pt[DILEPTON_MUMU]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
	hdilep_1_pt[DILEPTON_ALL]->Fill( cms2.mus_p4()[muidxs_nopt[1]].pt(), weight );
  }
  
  //enforce exactly two leptons and SF requirement
  //remember, nels, nmus require NO MET, NOR MASS CUTS
  if( (nels == 2 && nmus == 0)
	  || (nmus == 2 && nels == 0) )
	ZEvent();

  else if( (nels == 1 && nmus == 0) || (nmus == 1 && nels == 0) ) {
	unsigned int lep_type = (nels == 1 ? 0 : 1);
	//for met plots, only need to pass tmass (all lep cuts are passed via nels, nmus if)
	if( (pass_lep_cut & CUT_BIT(TMASS)) == CUT_BIT(TMASS) ) {
	  hlep_tcmet[lep_type]->Fill(cms2.evt_tcmet(), weight);
	  hlep_tcmet[2]->Fill(cms2.evt_tcmet(), weight);
	  hlep_clmumet[lep_type]->Fill(cms2.evt_metMuonCorr(), weight);
	  hlep_clmumet[2]->Fill(cms2.evt_metMuonCorr(), weight);
	}
	//for WEvent, need to pass tmass, met
	if( (pass_lep_cut & w_evt_sel) == w_evt_sel ) 
	  WEvent();
  }

}
//end FillEventHistos

void Looper::WEvent() {
  
  double weight = Looper::Weight();

  // histogram indices are e, m, all (0, 1, 2)
  // get lep_type, lep_p4
  unsigned int lep_type = 0; //default el
  LorentzVector lep_p4;
  if( elidxs[0] == -1 ) { //no els
	lep_type = 1;
	lep_p4 = cms2.mus_p4()[muidxs[0]];
  }
  else 
	lep_p4 = cms2.els_p4()[elidxs[0]];

  //jet vars--single lepton
  int njets_20 = 0;
  int njets_30 = 0;
  //this code ripped from selections.cc->getCaloJets
  //vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets;
  //calo_jets.clear();
  
  for( unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
    //if( dRbetweenVectors(lep_p4, cms2.jets_p4()[jj]) < 0.4 )
	if(  ROOT::Math::VectorUtil::DeltaR(lep_p4, cms2.jets_p4()[jj]) < 0.4 )
		//(dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.jets_p4()[jj]) < 0.4)
	  continue;
    if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4)
	  continue;
    //if (cms2.jets_emFrac()[jj] < 0.1) continue;
	//count
    if (cms2.jets_p4()[jj].pt() > 20)
	  njets_20++;
    if (cms2.jets_p4()[jj].pt() > 30)
	  njets_30++;
    //calo_jets.push_back(cms2.jets_p4()[jj]);
  }
  
  //if (calo_jets.size() > 1) {
  //sort(calo_jets.begin(), calo_jets.end(),  comparePt);
  //}
  //return calo_jets;

  hlep_njet20[lep_type]->Fill(njets_20, weight);
  hlep_njet20[2]->Fill(njets_20, weight);
  hlep_njet30[lep_type]->Fill(njets_30, weight);
  hlep_njet30[2]->Fill(njets_30, weight);


  if( lep_type == 0 ) {
	//conversion plot has ALL cuts applied
	hlep_conv[lep_type]->Fill( conversionElectron(elidxs[0]), weight );
	hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );

	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.els_p4()[elidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);
  }
  else if (lep_type == 1) {
	float dphi = acos(cos(cms2.evt_tcmetPhi() - cms2.mus_p4()[muidxs[0]].Phi() ));
	hlep_met_dphi[lep_type]->Fill(dphi, weight);
	hlep_met_dphi[2]->Fill(dphi, weight);
  }

  scands_passing_[lep_type] += weight;
  scands_passing_w2_[lep_type] += weight * weight;
  scands_count_[lep_type]++;
  scands_passing_[2] += weight;
  scands_passing_w2_[2] += weight * weight;
  scands_count_[2]++;

  //} //end if passed cuts

}

void Looper::ZEvent ()
{
  double weight = Looper::Weight();

  //hdilep_nhyp_->Fill(cms2.hyp_p4().size(), weight);

  // define the cuts to be used
  cuts_t cuts = (CUT_BIT(DILEP_PT)
				 | CUT_BIT(DILEP_OS)
				 | CUT_BIT(DILEP_MASS)
				 //| CUT_BIT(LEP_GOOD) //these already checked
				 //| CUT_BIT(LEP_ISO)
				 );
  cuts_t cuts_nomass = (CUT_BIT(DILEP_PT) | CUT_BIT(DILEP_OS));

  //get leptons from indicies, not hyps
  //for( unsigned int i = 0; i < cms2.hyp_p4().size(); ++i) {
  // get what type of di-lepton hypothesis this is
  //const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i]);
  //DILEPTON_ALL,     DILEPTON_MUMU,     DILEPTON_EMU,     DILEPTON_EE

  const enum DileptonHypType myType = (elidxs[0] != -1 ? DILEPTON_EE : DILEPTON_MUMU);
  //int myType = (elidxs[0] != -1 ? 3 : 1);

  // does this hypothesis pass the required cuts?
  cuts_t cuts_passed = DilepSelect();

  //these already checked
  //require both hyps to pass lep select (1st arg:0=e,1=m)
  //cuts_passed |= LepSelect( abs(cms2.hyp_lt_id()[i]) == 11 ? 0 : 1, cms2.hyp_lt_index()[i] );
  //cuts_passed |= LepSelect( abs(cms2.hyp_ll_id()[i]) == 11 ? 0 : 1, cms2.hyp_ll_index()[i] );

  //lep vars
  LorentzVector lep1_p4, lep2_p4;	
  double pt1=0, pt2=0, mass=0;
  if( elidxs[0] != -1 && elidxs[1] != -1 ) {
	lep1_p4 = cms2.els_p4()[elidxs[0]];
	lep2_p4 = cms2.els_p4()[elidxs[1]];
	pt1 = cms2.els_p4()[elidxs[0]].pt();
	pt2 = cms2.els_p4()[elidxs[1]].pt();
	mass = (cms2.els_p4()[elidxs[0]] + cms2.els_p4()[elidxs[1]]).mass();
  }
  else if( muidxs[0] != -1 && muidxs[1] != -1 ) {
	lep1_p4 = cms2.mus_p4()[muidxs[0]];
	lep2_p4 = cms2.mus_p4()[muidxs[1]];
	pt1 = cms2.mus_p4()[muidxs[0]].pt();
	pt2 = cms2.mus_p4()[muidxs[1]].pt();
	mass = (cms2.mus_p4()[muidxs[0]] + cms2.mus_p4()[muidxs[1]]).mass();
  }
  else {
	cout << "BAD LEP INDXS IN ZEVENT\n" << endl;
	exit(1);
  }

  //fill pt hists before/during lep selection
  // fill histograms for the 0th lepton
  //hdilep_0_pt[myType]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
  //hdilep_0_pt[DILEPTON_ALL]->Fill(cms2.hyp_lt_p4()[i].pt(), weight);
  //hdilep_0_pt[myType]->Fill(pt1, weight);
  //hdilep_0_pt[DILEPTON_ALL]->Fill(pt1, weight);
	
  // fill histograms for the 1th lepton
  //hdilep_1_pt[myType]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
  //hdilep_1_pt[DILEPTON_ALL]->Fill(cms2.hyp_ll_p4()[i].pt(), weight);
  //hdilep_1_pt[myType]->Fill(pt2, weight);
  //hdilep_1_pt[DILEPTON_ALL]->Fill(pt2, weight);
	
  // fill hypothesis level histograms
  //hdilep_mass[myType]->Fill(cms2.hyp_p4()[i].mass(), weight);
  //hdilep_mass[DILEPTON_ALL]->Fill(cms2.hyp_p4()[i].mass(), weight);
  if( (cuts_passed & cuts_nomass) == cuts_nomass ) {
	hdilep_mass[myType]->Fill(mass, weight);
	hdilep_mass[DILEPTON_ALL]->Fill(mass, weight);
  }

  if( (cuts_passed & cuts) == cuts ) {

	//jet vars--dilep case
	int njets_20 = 0;
	int njets_30 = 0;
	//this code ripped from selections.cc->getCaloJets
	//vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > calo_jets;
	//calo_jets.clear();
  
	for( unsigned int jj=0; jj < cms2.jets_p4().size(); ++jj) {
	  //if( dRbetweenVectors(lep1_p4, cms2.jets_p4()[jj]) < 0.4 ||
	  //  dRbetweenVectors(lep2_p4, cms2.jets_p4()[jj]) < 0.4 )
	  if(  ROOT::Math::VectorUtil::DeltaR(lep1_p4, cms2.jets_p4()[jj]) < 0.4 ||
		   ROOT::Math::VectorUtil::DeltaR(lep2_p4, cms2.jets_p4()[jj]) < 0.4 )
		  continue;
	  if (fabs(cms2.jets_p4()[jj].Eta()) > 2.4)
		continue;
	  //if (cms2.jets_emFrac()[jj] < 0.1) continue;
	  //count
	  if (cms2.jets_p4()[jj].pt() > 20)
		njets_20++;
	  if (cms2.jets_p4()[jj].pt() > 30)
		njets_30++;
	  //calo_jets.push_back(cms2.jets_p4()[jj]);
	}
  
	//if (calo_jets.size() > 1) {
	//sort(calo_jets.begin(), calo_jets.end(),  comparePt);
	//}
	//return calo_jets;

	if( myType == DILEPTON_EE ) {
	  hlep_conv[0]->Fill( conversionElectron(elidxs[0]), weight );
	  hlep_conv[2]->Fill( conversionElectron(elidxs[0]), weight );
	}
	
	hdilep_njet20[myType]->Fill(njets_20, weight);
	hdilep_njet20[DILEPTON_ALL]->Fill(njets_20, weight);
	hdilep_njet30[myType]->Fill(njets_30, weight);
	hdilep_njet30[DILEPTON_ALL]->Fill(njets_30, weight);

	//no met cut, so require all cuts for met plot
	hdilep_tcmet[myType]->Fill(cms2.evt_tcmet(), weight);
	hdilep_tcmet[DILEPTON_ALL]->Fill(cms2.evt_tcmet(), weight);
	hdilep_clmumet[myType]->Fill(cms2.evt_metMuonCorr(), weight);
	hdilep_clmumet[DILEPTON_ALL]->Fill(cms2.evt_metMuonCorr(), weight);

	//all other hists need to fill before checking cuts

	dcands_passing_[myType] += weight;
	dcands_passing_w2_[myType] += weight * weight;
	dcands_count_[myType]++;
	dcands_passing_[DILEPTON_ALL] += weight;
	dcands_passing_w2_[DILEPTON_ALL] += weight * weight;
	dcands_count_[DILEPTON_ALL]++;

  }
  //} //end loop on hyp

}


bool Looper::GoodTMTestMuonWithoutIsolation(int index, int mode) {
  //modes:
  // mode = 1: TMLS tight, d0, ecal, hcal, eta
  // mode = 2: TMLS tight, d0, eta
  // mode = 3: TMLS loose, d0, ecal, hcal, eta
  // mode = 4: TMLS loose, d0, eta
  // mode = 5: TMC2D tight, d0, eta
  // mode = 6: TMC2D loose, d0, eta

  if (((cms2.mus_type().at(index)) & (1<<2)) == 0) return false; // not a tracker muon - skip further consideration

  bool isGoodMu = true;
  if(mode == 1) {
    if (!cms2.mus_pid_TMLastStationTight().at(index))     isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
    if (cms2.mus_pat_ecalvetoDep().at(index) >= 4)   isGoodMu = false; // ECalE < 4
    if (cms2.mus_pat_hcalvetoDep().at(index) >= 6)   isGoodMu = false; // HCalE < 6
   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }
  else if(mode == 2) {
    if (!cms2.mus_pid_TMLastStationTight().at(index))     isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }
  else if(mode == 3) {
    if (!cms2.mus_pid_TMLastStationLoose().at(index))     isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
    if (cms2.mus_pat_ecalvetoDep().at(index) >= 4)   isGoodMu = false; // ECalE < 4
    if (cms2.mus_pat_hcalvetoDep().at(index) >= 6)   isGoodMu = false; // HCalE < 6
   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }
  else if(mode == 4) {
    if (!cms2.mus_pid_TMLastStationLoose().at(index))     isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }
  else if(mode == 5) {
    if (!cms2.mus_pid_TM2DCompatibilityTight().at(index)) isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }
  else if(mode == 6) {
    if (!cms2.mus_pid_TM2DCompatibilityLoose().at(index)) isGoodMu = false;
    if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
    if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;
  }

  else isGoodMu = false;

  return isGoodMu;


// mus_type
// mus_goodmask
// mus_p4
// mus_trk_p4
// mus_d0
// mus_z0
// mus_d0corr
// mus_z0corr
// mus_vertexphi
// mus_chi2
// mus_ndof
// mus_validHits
// mus_lostHits
// mus_d0Err
// mus_z0Err
// mus_ptErr
// nchAlias(
// mus_phiErr
// mus_charge
// mus_trk_charge
  
// mus_qoverp
  
// mus_qoverpError
  
// mus_outerPhi
// mus_outerEta
// mus_trkrefkey

// mus_nmatches
// mus_e_em
// mus_e_had
// mus_e_ho
// mus_e_emS9
// mus_e_hadS9
// mus_e_hoS9
// mus_iso
// mus_iso03_sumPt
// mus_iso03_emEt
// mus_iso03_hadEt
// mus_iso03_hoEt
// mus_iso03_ntrk
// mus_iso05_sumPt
// mus_iso05_emEt
// mus_iso05_hadEt
// mus_iso05_hoEt
// mus_iso05_ntrk
        
// mus_gfit_chi2
// mus_gfit_ndof
// mus_gfit_validHits
// mus_pid_TMLastStationLoose
// mus_pid_TMLastStationTight
// mus_pid_TM2DCompatibilityLoose
// mus_pid_TM2DCompatibilityTight
// mus_caloCompatibility
// mus_vertex_p4
// mus_gfit_outerPos_p4

// mus_pid_TMLastStationLoose
// mus_pid_TMLastStationTight
// mus_pid_TM2DCompatibilityLoose
// mus_pid_TM2DCompatibilityTight
// mus_caloCompatibility
//   if (((cms2.mus_type().at(index)) & (1<<1)) == 0) isGoodMu = false; // global muon
//   if (((cms2.mus_type().at(index)) & (1<<2)) == 0) isGoodMu = false; // tracker muon
//   if (cms2.mus_validHits().at(index) < 11)    isGoodMu = false;
//   if (cms2.mus_gfit_chi2().at(index)/cms2.mus_gfit_ndof().at(index) >= 10) isGoodMu = false;
//   if (fabs(cms2.mus_d0corr().at(index))   >= 0.02) isGoodMu = false;
//   if (cms2.mus_pat_ecalvetoDep().at(index) >= 4) isGoodMu = false; // ECalE < 4
//   if (cms2.mus_pat_hcalvetoDep().at(index) >= 6) isGoodMu = false; // HCalE < 6
//   if ( TMath::Abs(cms2.mus_p4()[index].eta()) > 2.4) isGoodMu = false;

}


void Looper::End ()
{
  /*
  int ret = fprintf(logfile_, 
					"Sample %10s: Total candidate count (ee em mm all): %8u %8u %8u %8u."
					" Total weight %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f %10.1f +- %10.1f\n",   
					sample_.name.c_str(),
					CandsCount(DILEPTON_EE), CandsCount(DILEPTON_EMU), CandsCount(DILEPTON_MUMU), CandsCount(DILEPTON_ALL), 
					CandsPassing(DILEPTON_EE)  , RMS(DILEPTON_EE),  
					CandsPassing(DILEPTON_EMU) , RMS(DILEPTON_EMU),  
					CandsPassing(DILEPTON_MUMU), RMS(DILEPTON_MUMU), 
					CandsPassing(DILEPTON_ALL) , RMS(DILEPTON_ALL));
  if (ret < 0)
	perror("writing to log file");
  */
}

