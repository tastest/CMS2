#include <math.h>
#include "TVector3.h"
#include "CORE/selections.h"
#include "CORE/utilities.h"
#include "CORE/CMS2.h"
#include "Tools/tools.h"
#include "Tools/fakerates.h"
#include "Looper.h"

static const double d0_bins[] = { 0, 0.025, 0.03, 0.055 };
static const double dphiin_bins[] = { -0.04, 0, 0.04, 0.045, 0.085 };
static const double iso_bins[] = { 0, 0.82, 0.9, 0.92, 1.0001};

bool passTrkjetCuts (int i_trk);
void calculateJetProb (const std::vector<std::pair<LorentzVector, std::vector<unsigned int> > > &trackJets,
		       vector<double> *jetprobs);

Looper::Looper (Sample s, cuts_t c, const char *fname) 
     : LooperBase(s, c, fname)
{
     // zero out the candidate counters (don't comment this out)
     memset(cands_passing_	, 0, sizeof(cands_passing_       ));
     memset(cands_passing_w2_	, 0, sizeof(cands_passing_w2_    ));
     memset(cands_count_		, 0, sizeof(cands_count_         ));
}

void Looper::BookHistos ()
{
       hnJet		= new NMinus1Hist(sample_, "nJet"            ,	 6	, -0.5, 5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
       hnCaloJet	= new NMinus1Hist(sample_, "nCaloJet"        ,	 6	, -0.5, 5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS))	);
       hnTrackJet	= new NMinus1Hist(sample_, "nTrackJet"       ,	 6	, -0.5, 5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) | (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
       hcaloJetPt	= new NMinus1Hist(sample_, "caloJetPt"       ,	 6	, -0.5, 5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_CALO)) 	);
//        hjpJetPt		= new NMinus1Hist(sample_, "jpJetPt"         ,	 100	, 0, 100	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_SIP)) 	);
       hjpJetNtr	= new NMinus1Hist(sample_, "jpJetNtr"        ,	 11	, -0.5, 10.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_SIP)) 	);
       hnjpJets		= new NMinus1Hist(sample_, "njpJets"         ,	 6	, -0.5, 5.5	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_SIP)) 	);
       hntracks		= new NMinus1Hist(sample_, "ntracks"         ,	 101	, -0.5, 100.5	, cuts_, 0);
       hminLepPt	= new NMinus1Hist(sample_, "minLepPt"        ,	 150	, 0, 150	, cuts_, CUT_BIT(CUT_LL_PT)	);
       hmaxLepPt	= new NMinus1Hist(sample_, "maxLepPt"        ,	 150	, 0, 150	, cuts_, CUT_BIT(CUT_LL_PT)  	);
       hltPt		= new NMinus1Hist(sample_, "ltPt"            ,	 150	, 0, 150	, cuts_, (CUT_BIT(CUT_LT_PT))	);
       hllPt		= new NMinus1Hist(sample_, "llPt"            ,	 150	, 0, 150	, cuts_, (CUT_BIT(CUT_LL_PT))	);
       helPt		= new NMinus1Hist(sample_, "elPt"            ,	 16	, 0, 160	, cuts_, CUT_BIT(CUT_LL_PT)	);
       hmuPt		= new NMinus1Hist(sample_, "muPt"            ,	 16	, 0, 160	, cuts_, CUT_BIT(CUT_LL_PT)	);
       helEta		= new NMinus1Hist(sample_, "elEta"           ,	 12	, -3, 3		, cuts_, 0);
       hmuEta		= new NMinus1Hist(sample_, "muEta"           ,	 12	, -3, 3		, cuts_, 0);
       hdphiLep		= new NMinus1Hist(sample_, "dphiLep"         ,	 50	, 0, M_PI	, cuts_, 0	);
       hdilMass		= new NMinus1Hist(sample_, "dilMass"         ,	 100	, 0, 300	, cuts_, CUT_BIT(CUT_PASS_ZVETO) | CUT_BIT(CUT_PASS_ADDZVETO) | CUT_BIT(CUT_IN_Z_WINDOW));
       hdilPt		= new NMinus1Hist(sample_, "dilPt"           ,	 100	, 0, 300	, cuts_, 0	);
       hmet		= new NMinus1Hist(sample_, "met"             ,	 100	, 0, 200	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)	);
       hmetSpec		= new NMinus1Hist(sample_, "metSpec"         ,	 100	, 0, 200	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)  );
       hmetTrkCorr	= new NMinus1Hist(sample_, "metTrkCorr"      ,	 100	, 0, 200	, cuts_, CUT_BIT(CUT_PASS4_MET) | CUT_BIT(CUT_PASS2_MET) | CUT_BIT(CUT_PASS4_TCMET) | CUT_BIT(CUT_PASS2_TCMET)	);
       hptJet1		= new NMinus1Hist(sample_, "ptJet1"          ,	 100	, 0, 300	, cuts_, 0		);
       hptJet2		= new NMinus1Hist(sample_, "ptJet2"          ,	 100	, 0, 300	, cuts_, 0       	);
       hptJet3		= new NMinus1Hist(sample_, "ptJet3"          ,	 100	, 0, 300	, cuts_, 0       	);
       hptJet4		= new NMinus1Hist(sample_, "ptJet4"          ,	 100	, 0, 300	, cuts_, 0       	);
       hetaJet1		= new NMinus1Hist(sample_, "etaJet1"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet2		= new NMinus1Hist(sample_, "etaJet2"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet3		= new NMinus1Hist(sample_, "etaJet3"         ,	 50	, -4, 4		, cuts_, 0       	);
       hetaJet4		= new NMinus1Hist(sample_, "etaJet4"         ,	 50	, -4, 4		, cuts_, 0       	);
       hnumTightLep	= new NMinus1Hist(sample_, "numTightLep"     ,	 6	, -0.5, 5.5	, cuts_, 0             	);
       heleRelIso	= new NMinus1Hist(sample_, "eleRelIso"       ,	 101	, 0, 1.01	, cuts_, 0		);
       heleRelIsoTrk	= new NMinus1Hist(sample_, "eleRelIsoTrk"    ,	 101	, 0, 1.01	, cuts_, 0		);
       hmuRelIso	= new NMinus1Hist(sample_, "muRelIso"        ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO))	);
       hminRelIso	= new NMinus1Hist(sample_, "minRelIso"       ,	 101	, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO))	);
       hminRelIso_withCalo = new NMinus1Hist	(sample_, "minRelIso_withCalo", 101, 0, 1.01	, cuts_, (CUT_BIT(CUT_LT_ISO)) | (CUT_BIT(CUT_LL_ISO))	);
       htagMuPt		= new NMinus1Hist(sample_, "tagMuPt"	      ,	 100	, 0, 100	, cuts_ & ~((CUT_BIT(CUT_PASS_MUON_B_VETO)) | (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO))), (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)));
       htagMuRelIso	= new NMinus1Hist(sample_, "tagMuRelIso"     ,	 101	, 0, 1.01	, cuts_ & ~((CUT_BIT(CUT_PASS_MUON_B_VETO)) | (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO))), (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)));
       hmuPdgId		= new NMinus1Hist(sample_, "muPdgId", 		2301, -0.5, 2300.5, cuts_, 0);
       hmuMoPdgId	= new NMinus1Hist(sample_, "muMoPdgId", 	2301, -0.5, 2300.5, cuts_, 0);
       helPdgId		= new NMinus1Hist(sample_, "elPdgId", 		2301, -0.5, 2300.5, cuts_, 0);
       helMoPdgId	= new NMinus1Hist(sample_, "elMoPdgId", 	2301, -0.5, 2300.5, cuts_, 0);
       helEop      	= new NMinus1Hist(sample_, "elEop"	      ,  10	, 0, 10	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       held0    	= new NMinus1Hist(sample_, "eld0"	      ,  50	, 0, 0.1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helfbrem    	= new NMinus1Hist(sample_, "elfbrem"	      ,  11	, -0.1, 1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helHE       	= new NMinus1Hist(sample_, "elHE"	      ,  11	, -0.03, 0.3	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsee      	= new NMinus1Hist(sample_, "elsee"	      ,  50	, 0, 0.05	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsppEB  	= new NMinus1Hist(sample_, "elsppEB"	      ,  50	, 0, 0.0	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helsppEE	  	= new NMinus1Hist(sample_, "elsppEC"	      ,  50	, 0, 0.0	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       heldphiin   	= new NMinus1Hist(sample_, "eldphiin"	      ,  10	, -0.1, 0.1	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       heldetain   	= new NMinus1Hist(sample_, "eldetain"	      ,  10	, -0.02, 0.02	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helEseedopin	= new NMinus1Hist(sample_, "elEseedopin"     ,  10	, 0, 20	, cuts_, (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_LL_CALOISO)));
       helConvDeltaPhi_ss = new NMinus1Hist   	(sample_, "elConvDeltaPhi_ss"	    ,  10	, 0, 0.1, cuts_, 0);
       helConvDeltaPhi_os = new NMinus1Hist	(sample_, "elConvDeltaPhi_os"     ,  10	, 0, 0.1, cuts_, 0);
       held0vsRelIso 		= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "d0vsRelIso"), ";d0;rel iso", sizeof(d0_bins) / sizeof(double) - 1, d0_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       heldphiinvsRelIso 	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "dphiinvsRelIso"), ";dphiin;rel iso", sizeof(dphiin_bins) / sizeof(double) - 1, dphiin_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       held0vsRelIsoMCgamma 	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "d0vsRelIsoMCgamma"), ";d0;rel iso", sizeof(d0_bins) / sizeof(double) - 1, d0_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       heldphiinvsRelIsoMCgamma = new TH2F(Form("%s_%s_em", sample_.name.c_str(), "dphiinvsRelIsoMCgamma"), ";dphiin;rel iso", sizeof(dphiin_bins) / sizeof(double) - 1, dphiin_bins, sizeof(iso_bins) / sizeof(double) - 1, iso_bins);
       htrkCalodRvsPtSum	= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "trkCalodRvsPtSum"), ";pt sum;#DeltaR", 10, 0, 30, 60, 0, M_PI);
       hCaloEtaPt		= new TH2F(Form("%s_%s_em", sample_.name.c_str(), "CaloEtaPt"), ";pt;#eta", 10, 0, 30, 10, -5, 5);
       hntrklj = new NMinus1Hist(sample_, "ntrklj" , 100, 0, 100, cuts_, 0 );
       hntrkbj = new NMinus1Hist(sample_, "ntrkbj" , 100, 0, 100, cuts_, 0 );
       hntrknm = new NMinus1Hist(sample_, "ntrknm" , 100, 0, 100, cuts_, 0 );
       hjplj = new NMinus1Hist(sample_, "jplj" , 40, 0, 20, cuts_, 0 );
       hjpbj = new NMinus1Hist(sample_, "jpbj" , 40, 0, 20, cuts_, 0 );
       hjpnm = new NMinus1Hist(sample_, "jpnm" , 40, 0, 20, cuts_, 0 );
       hsiplj = new NMinus1Hist(sample_, "siplj" , 100, -0.5, 0.5, cuts_, 0 );
       hsipbj = new NMinus1Hist(sample_, "sipbj" , 100, -0.5, 0.5, cuts_, 0 );
       hsipnm = new NMinus1Hist(sample_, "sipnm" , 100, -0.5, 0.5, cuts_, 0 );
       hsipsiglj = new NMinus1Hist(sample_, "sipsiglj" , 100, -10, 10, cuts_, 0 );
       hsipsigbj = new NMinus1Hist(sample_, "sipsigbj" , 100, -10, 10, cuts_, 0 );
       hsipsignm = new NMinus1Hist(sample_, "sipsignm" , 100, -10, 10, cuts_, 0 );
       
       htrd0 			= new NMinus1Hist(sample_, "trd0"                    ,	100, -0.1, 0.1	, cuts_, 0);
       htrd0sig 		= new NMinus1Hist(sample_, "trd0sig"                 ,	100, -10, 10  	, cuts_, 0);
       htrd0Strange 		= new NMinus1Hist(sample_, "trd0Strange"             ,	100, -0.1, 0.1	, cuts_, 0); // mother is a Ks or Lambda
       htrd0sigStrange		= new NMinus1Hist(sample_, "trd0sigStrange" 	     ,	100, -10, 10  	, cuts_, 0); // mother is a Ks or Lambda
       hjetpt 			= new NMinus1Hist(sample_, "jetpt"                   ,	100, 0, 100	, cuts_, 0);
       hjetntr	 		= new NMinus1Hist(sample_, "jetntr"                ,	21, -0.5, 20.5 	, cuts_, 0);
       hjetd0 			= new NMinus1Hist(sample_, "jetd0"                   ,	100, -0.1, 0.1	, cuts_, 0);
       hjetd0sig 		= new NMinus1Hist(sample_, "jetd0sig"                ,	100, -10, 10  	, cuts_, 0);
       hjetmaxd0sig 		= new NMinus1Hist(sample_, "jetmaxd0sig"             ,	100, 0, 10  	, cuts_, 0);
       htrkprob			= new NMinus1Hist(sample_, "trkprob"                 ,	100, 0, 1	, cuts_, 0);
       hjetprob			= new NMinus1Hist(sample_, "jetprob"                 ,	100, 0, 1	, cuts_, 0);
       hlogjetprob		= new NMinus1Hist(sample_, "logjetprob"              ,	21, -20, 1	, cuts_, 0);
       hminjetprob		= new NMinus1Hist(sample_, "minjetprob"              ,	120, 0, 1.2	, cuts_, 0);
       hmaxptjetprob		= new NMinus1Hist(sample_, "maxptjetprob"            ,	120, 0, 1.2	, cuts_, 0);
       for (int i = 0; i < 4; ++i) {
	    htrackJetPt[i]	= new NMinus1Hist(sample_, Form("%s%d", "trackJetPt"	 , i)	      ,	100, 0, 100	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS)) 	);
	    hjpJetPt[i]		= new NMinus1Hist(sample_, Form("%s%d", "jpJetPt"	 , i)	      ,	100, 0, 100	, cuts_, (CUT_BIT(CUT_PASS_JETVETO_SIP)) 	);
	    htrd0ByPt 	    [i]	= new NMinus1Hist(sample_, Form("%s%d", "trd0ByPt"       , i)         ,	100, -0.1, 0.1	, cuts_, 0);
	    htrd0sigByPt    [i]	= new NMinus1Hist(sample_, Form("%s%d", "trd0sigByPt"    , i)         ,	100, -10, 10  	, cuts_, 0);
	    htrd0ByNtrks    [i]	= new NMinus1Hist(sample_, Form("%s%d", "trd0ByNtrks"    , i)         ,	100, -0.1, 0.1	, cuts_, 0);
	    htrd0sigByNtrks [i]	= new NMinus1Hist(sample_, Form("%s%d", "trd0sigByNtrks" , i)         ,	100, -10, 10  	, cuts_, 0);
	    hjetd0ByPt 	    [i]	= new NMinus1Hist(sample_, Form("%s%d", "jetd0ByPt"      , i)         ,	100, -0.1, 0.1	, cuts_, 0);
	    hjetd0sigByPt   [i]	= new NMinus1Hist(sample_, Form("%s%d", "jetd0sigByPt"   , i)         ,	100, -10, 10  	, cuts_, 0);
	    hjetmaxd0sigByPt[i]	= new NMinus1Hist(sample_, Form("%s%d", "jetmaxd0sigByPt", i)         ,	100, 0, 10  	, cuts_, 0);
	    hjetd0ByNtrks   [i]	= new NMinus1Hist(sample_, Form("%s%d", "jetd0ByNtrks"   , i)         ,	100, -0.1, 0.1	, cuts_, 0);
	    hjetd0sigByNtrks[i]	= new NMinus1Hist(sample_, Form("%s%d", "jetd0sigByNtrks", i)         ,	100, -10, 10  	, cuts_, 0);
	    hjetmaxd0sigByNtrks[i]	= new NMinus1Hist(sample_, Form("%s%d", "jetmaxd0sigByNtrks", i)         ,	100, 0, 10  	, cuts_, 0);
	    hjetprobByPt[i]	= new NMinus1Hist(sample_, Form("%s%d", "jetprobByPt", i)      ,	100, 0, 1	, cuts_, 0);
	    hjetprobByNtrks[i]	= new NMinus1Hist(sample_, Form("%s%d", "jetprobByNtrks", i)      ,	100, 0, 1	, cuts_, 0);
	    hmaxptjetprobByPt[i]= new NMinus1Hist(sample_, Form("%s%d", "maxptjetprobByPt", i)      ,	100, 0, 1	, cuts_, 0);
       }
       hhypDeltaz0sig 	= new NMinus1Hist(sample_, "hypDeltaz0sig" , 100, -10, 10, cuts_, 0 );
       for (int i = 0; i < 3; ++i) {
	    htrkd0	   [i] = new NMinus1Hist(sample_, Form("%s%d", "trkd0"        	, i), 100, -1, 1,	cuts_, 0);
	    htrkDeltaz0	   [i] = new NMinus1Hist(sample_, Form("%s%d", "trkDeltaz0"   	, i), 100, -1, 1,	cuts_, 0);
	    htrkDeltaz0sig [i] = new NMinus1Hist(sample_, Form("%s%d", "trkDeltaz0sig"	, i), 100, -10, 10,	cuts_, 0);
	    htrknchi2	   [i] = new NMinus1Hist(sample_, Form("%s%d", "trknchi2"     	, i), 100, 0, 10,	cuts_, 0);
	    htrkvalidhits  [i] = new NMinus1Hist(sample_, Form("%s%d", "trkvalidhits" 	, i), 21, -0.5, 20.5,	cuts_, 0);
       }
       hjpJetPt0Pt1	= new TH2D(Form("%s_jpJetPt0Pt1_all", sample_.name.c_str()), Form("%s_jpJetPt0Pt1_all", sample_.name.c_str()), 
				   102, -2, 100, 102, -2, 100);
       hjpJetPt0Pt1->SetLineColor(sample_.histo_color);
       hjpJetPt0Ntr0	= new TH2D(Form("%s_jpJetPt0Ntr0_all", sample_.name.c_str()), Form("%s_jpJetPt0Ntr0_all", sample_.name.c_str()), 
				   102, -2, 100, 21, -0.5, 20.5);
       hjpJetPt0Ntr0->SetLineColor(sample_.histo_color);
       hjpJetPt0trJetPt0	= new TH2D(Form("%s_jpJetPt0trJetPt0_all", sample_.name.c_str()), Form("%s_jpJetPt0trJetPt0_all", sample_.name.c_str()), 
				   102, -2, 100, 102, -2, 100);
       hjpJetPt0trJetPt0->SetLineColor(sample_.histo_color);
}

bool Looper::FilterEvent()
{ 
     //
     // duplicate filter, based on trk information and dilepton hyp
     //
     if (cms2.trks_d0corr().size() == 0)
	  return true;
     DorkyEventIdentifier id = { cms2.evt_run(), cms2.evt_event(), cms2.trks_d0corr()[0], 
				 cms2.hyp_lt_p4()[0].pt(), cms2.hyp_lt_p4()[0].eta(), cms2.hyp_lt_p4()[0].phi() };
     return is_duplicate(id); 
}

cuts_t Looper::EventSelect ()
{
     //------------------------------------------------------------
     // In an event-based analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::DilepSelect (int i_hyp)
{
     cuts_t ret = 0;
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);

     // enough tracks?
     if (cms2.trks_trk_p4().size() > 2)
	  ret |= CUT_BIT(CUT_MORE_THAN_TWO_TRACKS);
     // pt cuts
     if (cms2.hyp_lt_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LT_PT));
     if (cms2.hyp_ll_p4()[i_hyp].pt() > 20.0) 
	  ret |= (CUT_BIT(CUT_LL_PT));
     // sign cuts
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] < 0 ) 
	  ret |= (CUT_BIT(CUT_OPP_SIGN));
     if ( cms2.hyp_lt_id()[i_hyp] * cms2.hyp_ll_id()[i_hyp] > 0 ) 
	  ret |= (CUT_BIT(CUT_SAME_SIGN));
     // track corrected MET
     const TVector3 trkCorr = correctMETforTracks();
     if (pass4Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS4_TCMET));
     if (pass2Met(i_hyp, trkCorr))
	  ret |= (CUT_BIT(CUT_PASS2_TCMET));
     // MET
     if (pass4Met(i_hyp, TVector3()))
	  ret |= (CUT_BIT(CUT_PASS4_MET));
     if (pass2Met(i_hyp, TVector3()))
	  ret |= (CUT_BIT(CUT_PASS2_MET));
     // muon quality
     int n_iso_mu = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) ) 
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_MU_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_TIGHT_DPHIIN)  | CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
	  n_iso_mu++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_TIGHT_DPHIIN)  | CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_MU_GOOD) | CUT_BIT(CUT_MU_ISO);
	  n_iso_mu++;
     }
     // electron quality
     int n_iso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_lt_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_EL_GOOD);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolation(cms2.hyp_ll_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_EL_GOOD);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_lt_index()[i_hyp], false)) {
	  ret |= CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_EL_ISO);
	  n_iso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && passElectronIsolation(cms2.hyp_ll_index()[i_hyp], false)) {
	  ret |= CUT_BIT(CUT_LL_ISO) | CUT_BIT(CUT_EL_ISO);
	  n_iso_el++;
     }     
     if (n_iso_mu + n_iso_el >= 1)
	  ret |= (CUT_BIT(CUT_ONE_ISO));
     if (n_iso_mu + n_iso_el >= 2)
	  ret |= (CUT_BIT(CUT_TWO_ISO));
     // electrons without d0
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_lt_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronWithoutIsolationWithoutd0(cms2.hyp_ll_index()[i_hyp]) )
	  ret |= CUT_BIT(CUT_EL_GOOD_NO_D0);
     // supertight cuts (only for electrons)
     int n_supertight_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (supertightElectron(cms2.hyp_lt_index()[i_hyp])) {
// 	       ret |= (CUT_BIT(CUT_LT_SUPERTIGHT));
	       n_supertight_el++;
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (supertightElectron(cms2.hyp_ll_index()[i_hyp])) {
// 	       ret |= (CUT_BIT(CUT_LL_SUPERTIGHT));
	       n_supertight_el++;
	  }
     }
     if (n_supertight_el >= 1)
	  ret |= CUT_BIT(CUT_ONE_SUPERTIGHT);
     if (n_supertight_el >= 2)
   	  ret |= CUT_BIT(CUT_TWO_SUPERTIGHT);
     // supertight dphiin cut
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (deltaPhiInElectron(cms2.hyp_lt_index()[i_hyp])) {
	       ret |= CUT_BIT(CUT_LT_TIGHT_DPHIIN);
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (deltaPhiInElectron(cms2.hyp_ll_index()[i_hyp])) {
	       ret |= CUT_BIT(CUT_LL_TIGHT_DPHIIN);
	  }
     }
     // barrel
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (fabs(cms2.els_p4()[cms2.hyp_lt_index()[i_hyp]].eta()) < 1.479)
 	       ret |= (CUT_BIT(CUT_EL_BARREL));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  if (fabs(cms2.els_p4()[cms2.hyp_ll_index()[i_hyp]].eta()) < 1.479)
 	       ret |= (CUT_BIT(CUT_EL_BARREL));
     }
     // calo iso
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_lt_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LT_GOOD)) | (CUT_BIT(CUT_LT_CALOISO));
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13 && goodMuonIsolated(cms2.hyp_ll_index()[i_hyp]) ) {
	  ret |= (CUT_BIT(CUT_LL_GOOD)) | (CUT_BIT(CUT_LL_CALOISO));
     }
     int n_caloiso_el = 0;
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_lt_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
	  n_caloiso_el++;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && goodElectronIsolated(cms2.hyp_ll_index()[i_hyp], true)) {
	  ret |= CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_CALOISO) | CUT_BIT(CUT_EL_CALOISO);
	  n_caloiso_el++;
     }     
     if (n_iso_mu + n_caloiso_el >= 1)
 	  ret |= (CUT_BIT(CUT_ONE_CALOISO));
     if (n_iso_mu + n_caloiso_el >= 2)
 	  ret |= (CUT_BIT(CUT_TWO_CALOISO));
     // jet veto
     if (cms2.hyp_njets()[i_hyp] == 0)
	  ret |= (CUT_BIT(CUT_PASS_JETVETO_CALO));
     // track jets
     if (passTrkJetVeto(i_hyp))
	  ret |= (CUT_BIT(CUT_PASS_JETVETO_TRACKJETS));
     // muon b tag, with 20 GeV upper cut on the muon pt
//      if (passMuonBVeto(i_hyp, true))
// 	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO));
//      else ret |= (CUT_BIT(CUT_MUON_TAGGED));
     // muon b tag, with no upper cut on the muon pt
//      if (passMuonBVeto(i_hyp, false))
// 	  ret |= (CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT));
//      else ret |= (CUT_BIT(CUT_MUON_TAGGED_WITHOUT_PTCUT));
     if (cms2.hyp_type()[i_hyp] == 1 || cms2.hyp_type()[i_hyp] == 2)
	  ret |= (CUT_BIT(CUT_PASS_ZVETO));
     else if (not inZmassWindow(cms2.hyp_p4()[i_hyp].mass()))
	  ret |= (CUT_BIT(CUT_PASS_ZVETO));
     else ret |= (CUT_BIT(CUT_IN_Z_WINDOW));
     // Z veto using additional leptons in the event
     if (not additionalZveto())
	  ret |= (CUT_BIT(CUT_PASS_ADDZVETO));
     // any additional high-pt, isolated leptons?
     if (passTriLepVeto(i_hyp))
	  ret |= (CUT_BIT(CUT_PASS_EXTRALEPTON_VETO));

     //*****************************************************************
     // track jet veto with signed impact parameter 
     //*****************************************************************
     bool pass_sip = true;
     for (unsigned int itrkjet = 0; itrkjet < TrackJets().size(); ++itrkjet) {
	  if (TrackJets()[itrkjet].first.pt() < 10)
	       break;
	  if (TrackJetProbs()[itrkjet] < 1e-2) {
	       pass_sip = false;
	       break;
	  }
     }
     if (pass_sip)
	  ret |= CUT_BIT(CUT_PASS_JETVETO_SIP);
	  
     //*****************************************************************
     // special handling for the fake rate cuts for now, because they
     // only work for emu
     //*****************************************************************
     if (myType != DILEPTON_EMU)
	  return ret;
     // in addition, for the muons, check that they pass tight+iso
     // (since the fake rate is electron only right now)
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  if ((ret & (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO))) != 
	      (CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LT_ISO)))
	       return ret;
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  if ((ret & (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO))) != 
	      (CUT_BIT(CUT_LL_GOOD) | CUT_BIT(CUT_LL_ISO)))
	       return ret;
     }
     // now set the fake flags for the electron
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  if (isFakeable(cms2.hyp_lt_index()[i_hyp]))
	       ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
	  if (isNumeratorElectron(cms2.hyp_lt_index()[i_hyp]))
	       ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
	  else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
     } else {
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	       if (isFakeable(cms2.hyp_ll_index()[i_hyp]))
		    ret |= CUT_BIT(CUT_ELFAKE_FAKEABLE_OBJECT);
	       if (isNumeratorElectron(cms2.hyp_ll_index()[i_hyp]))
		    ret |= CUT_BIT(CUT_ELFAKE_NUMERATOR);
	       else ret |= CUT_BIT(CUT_ELFAKE_NOT_NUMERATOR);
	  }
     }

     return ret;
}

double Looper::Weight (int)
{
     return cms2.evt_scale1fb() * sample_.kFactor;
}

cuts_t Looper::TrilepSelect (int i_hyp)
{
     //------------------------------------------------------------
     // In a trilepton analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

cuts_t Looper::QuadlepSelect (int i_hyp)
{
     //------------------------------------------------------------
     // In a quadlepton analysis, you would make your cuts here
     //------------------------------------------------------------

     cuts_t ret = 0;
     return ret;
}

void Looper::FillEventHistos ()
{

}

bool compareJetpt (const std::pair<LorentzVector, std::vector<unsigned int> > &j1, 
		   const std::pair<LorentzVector, std::vector<unsigned int> > &j2) 
{
     return j1.first.pt() > j2.first.pt();
}

void Looper::MakeTrackJets (int i_hyp)
{
     trackjets_.clear();
     for ( unsigned int itrkjet = 0; itrkjet < cms2.trkjets_p4().size(); ++itrkjet) {
	  // get rid of trkjets that are hypothesis electrons
	  if ((abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4) ||
	      (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4))
	       continue;
	  // why not muons?
	  TLorentzVector trkJet(cms2.trkjets_p4()[itrkjet].Px(), cms2.trkjets_p4()[itrkjet].Py(),
				cms2.trkjets_p4()[itrkjet].Pz(), cms2.trkjets_p4()[itrkjet].E());

	  vector<unsigned int> trkidx;
	  // Loop over all tracks and match to the track jet
	  for (unsigned int trkIter = 0; trkIter < cms2.trks_trk_p4().size(); ++trkIter)
	  {
	       TLorentzVector trk(cms2.trks_trk_p4()[trkIter].Px(), cms2.trks_trk_p4()[trkIter].Py(),
				  cms2.trks_trk_p4()[trkIter].Pz(), cms2.trks_trk_p4()[trkIter].E() );
	       
	       if (trkJet.DeltaR(trk) < 0.5) {
		    // This track matches the track jet, make some quality cuts
		    if (passTrkjetCuts(trkIter)) {
			 trkidx.push_back(trkIter);
		    }
	       }
	  }
	  std::pair<LorentzVector, std::vector<unsigned int> > trkjet(cms2.trkjets_p4()[itrkjet],
								      trkidx);
	  trackjets_.push_back(trkjet);
     }
     // sort track jets by pt
     sort(trackjets_.begin(), trackjets_.end(), compareJetpt);
}

void Looper::MakeTrackJetProbs ()
{
     calculateJetProb(TrackJets(), &trackjet_jp_);
}
    
void Looper::FillDilepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // Example dilepton histo filling; edit for your application
     //------------------------------------------------------------
     
     // make trackjets (from the ntuple or (by virtual function call) by kt clustering)
     MakeTrackJets(i_hyp);
     // and cache jet probabilities
     MakeTrackJetProbs();

     // every histogram needs to know what hypothesis he is 
     const enum DileptonHypType myType = hyp_typeToHypType(cms2.hyp_type()[i_hyp]);
     // and what the event weight is 
     const double weight = Weight(i_hyp);
     
     // these are the cuts that the candidate passes:
     cuts_t cuts_passed = DilepSelect(i_hyp);
     
     if ((cuts_passed & cuts_) == cuts_) {
	  cands_passing_[myType] += weight;
	  cands_passing_w2_[myType] += weight * weight;
	  cands_count_[myType]++;
	  cands_passing_[DILEPTON_ALL] += weight;
	  cands_passing_w2_[DILEPTON_ALL] += weight * weight;
	  cands_count_[DILEPTON_ALL]++;
     }

     hntracks->Fill(cuts_passed, myType, cms2.trks_trk_p4().size(), weight);

     // sip jet properties
     int nsipjets = 0;
     double sippts[4] = { -1, -1, -1, -1 };
     int sipntrs[4] = { -1, -1, -1, -1 };
     double trjpts[4] = { -1, -1, -1, -1 };
     for (unsigned int itrkjet = 0; itrkjet < TrackJets().size(); ++itrkjet) {
	  if (itrkjet < 4)
	       trjpts[itrkjet] = TrackJets()[itrkjet].first.pt();
	  if (TrackJetProbs()[itrkjet] > 1e-2) 
	       continue;
	  if (nsipjets < 4) {
	       hjpJetPt[nsipjets]->Fill(cuts_passed, myType, TrackJets()[itrkjet].first.pt(), weight);
	       sippts[nsipjets] = TrackJets()[itrkjet].first.pt();
	       sipntrs[nsipjets] = TrackJets()[itrkjet].second.size();
	  }
	  hjpJetNtr->Fill(cuts_passed, myType, TrackJets()[itrkjet].second.size(), weight);
	  nsipjets++;
     }
     hnjpJets->Fill(cuts_passed, myType, nsipjets, weight);
     if ((cuts_passed & cuts_) == cuts_) {
	  hjpJetPt0Pt1->Fill(sippts[0], sippts[1], weight);
	  hjpJetPt0Ntr0->Fill(sippts[0], sipntrs[0], weight);
	  hjpJetPt0trJetPt0->Fill(sippts[0], trjpts[0], weight);
     }

     // signed impact parameter study (AY && fkw)
#if 0
     if ((cuts_passed & cuts_) == cuts_) {
	  for ( unsigned int itrkjet = 0; itrkjet < cms2.trkjets_p4().size(); ++itrkjet) {

	       // get rid of trkjets that are hypothesis electrons
	       if ( (abs(cms2.hyp_lt_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_lt_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4) ||
		    (abs(cms2.hyp_ll_id()[i_hyp]) == 11 && dRbetweenVectors(cms2.hyp_ll_p4()[i_hyp],cms2.trkjets_p4()[itrkjet]) < 0.4)
		    ) continue;

	       int iGen = match4vector(cms2.trkjets_p4()[itrkjet],cms2.genps_p4(),0.3); // match track jet to generator particle
	       vector<unsigned int> idxIntrkJet = idxInCone( cms2.trkjets_p4()[itrkjet],cms2.trks_trk_p4()); // create vector of indices of all tracks in the track jet

	       std::vector<double> sImpactVec;
	       std::vector<double> sImpactErrVec;
	       std::vector<unsigned int>    trackIndexVec;
	       getSignedImpact(itrkjet, sImpactVec, sImpactErrVec, trackIndexVec); // calculate the signed impact parameter for all tracks in the track jet
	       double jetProb = getJetProb(sImpactVec, sImpactErrVec, 0.1, 2.0); // calculate the jet probability for each track jet

	       if ( abs(cms2.genps_id()[iGen]) == 5) {

		    hntrkbj->Fill(cuts_passed, myType, idxIntrkJet.size(), weight);
		    hjpbj->Fill(cuts_passed, myType, jetProb, weight);

		    for ( unsigned int it = 0; it < sImpactVec.size(); ++it) {
			 hsipbj->Fill(cuts_passed, myType, sImpactVec[it], weight);
			 hsipsigbj->Fill(cuts_passed, myType, (sImpactVec[it] / sImpactErrVec[it]), weight);
		    }

	       }

	       else if ( abs(cms2.genps_id()[iGen]) < 5) {

		    hntrklj->Fill(cuts_passed, myType, idxIntrkJet.size(), weight);
		    hjplj->Fill(cuts_passed, myType, jetProb, weight);

		    for ( unsigned int it = 0; it < sImpactVec.size(); ++it) {
			 hsiplj->Fill(cuts_passed, myType, sImpactVec[it], weight);
			 hsipsiglj->Fill(cuts_passed, myType, (sImpactVec[it] / sImpactErrVec[it]), weight);
		    }

	       }

	       else {
		    hntrknm->Fill(cuts_passed, myType, idxIntrkJet.size(), weight);
		    hjpnm->Fill(cuts_passed, myType, jetProb, weight);

		    for ( unsigned int it = 0; it < sImpactVec.size(); ++it) {
			 hsipnm->Fill(cuts_passed, myType, sImpactVec[it], weight);
			 hsipsignm->Fill(cuts_passed, myType, (sImpactVec[it] / sImpactErrVec[it]), weight);
		    }
	       }
	  }
     }
#endif
#if 1
     // signed impact parameter study (FG3 && J.Mü.)
     double min_jetprob = 1.1;
     double max_ptjet = 0;
     double jetprob_max_ptjet = 1.1;
     for (unsigned int itrkjet = 0; itrkjet < TrackJets().size(); ++itrkjet) {
	  int i_pt;
	  if (TrackJets()[itrkjet].first.pt() > 15)
	       i_pt = 0;
	  else if (TrackJets()[itrkjet].first.pt() > 10)
	       i_pt = 1;
	  else if (TrackJets()[itrkjet].first.pt() > 5)
	       i_pt = 2;
	  else i_pt = 3;
	  hjetpt->Fill(cuts_passed, myType, TrackJets()[itrkjet].first.pt(), weight);
	  
	  int i_ntr;
	  if (TrackJets()[itrkjet].second.size() > 5)
	       i_ntr = 0;
	  else if (TrackJets()[itrkjet].second.size() > 3)
	       i_ntr = 1;
	  else if (TrackJets()[itrkjet].second.size() > 1)
	       i_ntr = 2;
	  else i_ntr = 3;
	  hjetntr->Fill(cuts_passed, myType, TrackJets()[itrkjet].second.size(), weight);
	  
	  double jetprob = 1;
	  int n_jetprob_trks = 0;
	  double jet_sip = 0;
	  double jet_sipsig = 0;
	  double jet_maxsipsig = -999;
	  for (unsigned int i = 0; i < TrackJets()[itrkjet].second.size(); ++i) {
	       std::pair<double, double> sip = 
		    signedImpact(TrackJets()[itrkjet].first, TrackJets()[itrkjet].second[i]);
	       double d0 = sip.first;
	       double d0err = sip.second;
	       const double sipsig = d0 / d0err;
	       htrd0	->Fill(cuts_passed, myType, d0, weight);
	       htrd0sig	->Fill(cuts_passed, myType, sipsig, weight);
	       htrd0ByPt	[i_pt]->Fill(cuts_passed, myType, d0, weight);
	       htrd0sigByPt	[i_pt]->Fill(cuts_passed, myType, sipsig, weight);
	       htrd0ByNtrks	[i_ntr]->Fill(cuts_passed, myType, d0, weight);
	       htrd0sigByNtrks	[i_ntr]->Fill(cuts_passed, myType, sipsig, weight);
	       unsigned int trkidx = TrackJets()[itrkjet].second[i]; 
	       int momcid = cms2.trk_mc_motherid()[trkidx];
	       // see if we need to put a cutoff for strange hadrons
	       switch (abs(momcid)) {
	       case 130: case 310: case 311: // Kl/Ks/K0
	       case 3122: // Lambda
	       case 3222: case 3112: // Sigma+, Sigma-
	       case 3322: case 3312: // Cascade0, Cascade-
	       case 3334: // Omega-
		    htrd0Strange->Fill(cuts_passed, myType, d0, weight);
		    htrd0sigStrange->Fill(cuts_passed, myType, sipsig, weight);
		    break;
	       default:
		    break;
	       }
	       if (sipsig > 2) {
		    jet_sip += d0;
		    jet_sipsig += sipsig;
	       }
	       if (sipsig > jet_maxsipsig)
		    jet_maxsipsig = sipsig;
	       if (sipsig > 0) {
		    jetprob *= exp(-0.5 * sipsig * sipsig);
		    n_jetprob_trks++;
		    if (TrackJets()[itrkjet].second.size() > 3 && 
			TrackJets()[itrkjet].first.pt() > 10) {
			 htrkprob->Fill(cuts_passed, myType, exp(-0.5 * sipsig * sipsig), weight);
		    }
	       }
	  }
	  hjetd0		       ->Fill(cuts_passed, myType, jet_sip, 	weight);
	  hjetd0sig	       ->Fill(cuts_passed, myType, jet_sipsig,	weight);
	  hjetmaxd0sig	       ->Fill(cuts_passed, myType, jet_maxsipsig,	weight);
	  hjetd0ByPt	[i_pt] ->Fill(cuts_passed, myType, jet_sip,     weight);
	  hjetd0sigByPt	[i_pt] ->Fill(cuts_passed, myType, jet_sipsig,  weight);
	  hjetmaxd0sigByPt	[i_pt] ->Fill(cuts_passed, myType, jet_maxsipsig,  weight);
	  hjetd0ByNtrks	[i_ntr]->Fill(cuts_passed, myType, jet_sip,     weight);
	  hjetd0sigByNtrks	[i_ntr]->Fill(cuts_passed, myType, jet_sipsig,  weight);
	  hjetmaxd0sigByNtrks	[i_ntr]->Fill(cuts_passed, myType, jet_maxsipsig,  weight);
	  hjetprobByPt[i_pt]	->Fill(cuts_passed, myType, jetprob,  weight);
	  hjetprobByNtrks[i_ntr]->Fill(cuts_passed, myType, jetprob,  weight);
	  if (TrackJets()[itrkjet].second.size() > 3 && 
	      TrackJets()[itrkjet].first.pt() > 10) {
// 	       double jp = pow(jetprob, 1. / n_jetprob_trks);
	       double jp = jetprob;
	       hjetprob->Fill(cuts_passed, myType, jp, weight);
	       hlogjetprob->Fill(cuts_passed, myType, log(jp), weight);
	       if (jp < min_jetprob)
		    min_jetprob = jp;
	  }
	  if (TrackJets()[itrkjet].second.size() > 3 && 
	      TrackJets()[itrkjet].first.pt() > max_ptjet) {
	       max_ptjet = TrackJets()[itrkjet].first.pt();
	       jetprob_max_ptjet = jetprob;
	  }
     }
     hminjetprob->Fill(cuts_passed, myType, min_jetprob, weight);
     hmaxptjetprob->Fill(cuts_passed, myType, jetprob_max_ptjet, weight);
     int i_pt;
     if (max_ptjet > 15)
	  i_pt = 0;
     else if (max_ptjet > 10)
	  i_pt = 1;
     else if (max_ptjet > 5)
	  i_pt = 2;
     else i_pt = 3;
     hmaxptjetprobByPt[i_pt]->Fill(cuts_passed, myType, jetprob_max_ptjet, weight);
#endif
     // jet count
     hnJet->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     hnCaloJet	->Fill(cuts_passed, myType, cms2.hyp_njets()[i_hyp], weight);
     hnTrackJet	->Fill(cuts_passed, myType, nTrkJets(i_hyp), weight);
     // jet momentum
     for (unsigned int i = 0; i < TrackJets().size() && i < 4; ++i) {
	  htrackJetPt[i]->Fill(cuts_passed, myType, TrackJets()[i].first.pt(), weight);
     }
     // track quality
     const double hyp_delta_z0 = cms2.hyp_lt_z0()[i_hyp] - cms2.hyp_ll_z0()[i_hyp];
     const double hyp_z0 = 0.5 * (cms2.hyp_lt_z0()[i_hyp] + cms2.hyp_ll_z0()[i_hyp]);
     const double hyp_delta_z0Err = sqrt(cms2.hyp_lt_z0Err()[i_hyp] * cms2.hyp_lt_z0Err()[i_hyp] + 
					 cms2.hyp_ll_z0Err()[i_hyp] * cms2.hyp_ll_z0Err()[i_hyp]);
     hhypDeltaz0sig->Fill(cuts_passed, myType, hyp_delta_z0 / hyp_delta_z0Err, weight);
     for (unsigned int i = 0; i < cms2.trks_trk_p4().size(); ++i) {
	  // quality index: 
	  // 0 for tracks with |ip sig| < 2; 
	  // 1 for 2 < // |ip sig| < 5; 
	  // 2 for |ip sig| > 5
	  int i_qual = 0;
	  if (fabs(cms2.trks_d0corr()[i] / cms2.trks_d0Err()[i]) > 2)
	       i_qual = 1;
	  if (fabs(cms2.trks_d0corr()[i] / cms2.trks_d0Err()[i]) > 5)
	       i_qual = 2;
	  enum { D0, Z0, NCHI2, HITS };
	  unsigned int trkcuts = 0;
	  if (fabs(cms2.trks_d0corr()[i]) < 0.1)
	       trkcuts |= 1 << D0;
	  if (fabs(cms2.trks_z0()[i] - hyp_z0) < 0.3)
	       trkcuts |= 1 << Z0;
	  if (cms2.trks_chi2()[i] / cms2.trks_ndof()[i] < 4)
	       trkcuts |= 1 << NCHI2;
	  if (cms2.trks_validHits()[i] > 6)
	       trkcuts |= 1 << HITS;
	  const unsigned int all = 1 << D0 | 1 << Z0 | 1 << NCHI2 | 1 << HITS;
	  if (((trkcuts | 1 << D0) & all) == all)
	       htrkd0	[i_qual]->Fill(cuts_passed, myType, cms2.trks_d0corr()[i]	, weight);
	  if (((trkcuts | 1 << Z0) & all) == all)
	       htrkDeltaz0	[i_qual]->Fill(cuts_passed, myType, cms2.trks_z0()[i] - hyp_z0	, weight);
	  if (((trkcuts | 1 << Z0) & all) == all)
	       htrkDeltaz0sig[i_qual]->Fill(cuts_passed, myType, 
					    (cms2.trks_z0()[i] - hyp_z0) / 
					    sqrt(0.25 * hyp_delta_z0Err * hyp_delta_z0Err + 
						 cms2.trks_z0Err()[i] * cms2.trks_z0Err()[i]), weight);
	  if (((trkcuts | 1 << NCHI2) & all) == all)
	       htrknchi2	[i_qual]->Fill(cuts_passed, myType, cms2.trks_chi2()[i] / cms2.trks_ndof()[i], weight);
	  if (((trkcuts | 1 << HITS) & all) == all)
	       htrkvalidhits	[i_qual]->Fill(cuts_passed, myType, cms2.trks_validHits()[i], weight);
     }

     // lepton pt's
     hminLepPt->Fill(cuts_passed, myType, std::min(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hmaxLepPt->Fill(cuts_passed, myType, std::max(cms2.hyp_ll_p4()[i_hyp].pt(), cms2.hyp_lt_p4()[i_hyp].pt()), weight);
     hltPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
     hllPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  helPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  helEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].pt(), weight);
	  hmuEta->Fill(cuts_passed, myType, cms2.hyp_lt_p4()[i_hyp].eta(), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  helPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  helEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     } else {
	  hmuPt->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].pt(), weight);
	  hmuEta->Fill(cuts_passed, myType, cms2.hyp_ll_p4()[i_hyp].eta(), weight);
     }
    
     // dilepton mass
     hdilMass->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].mass(), weight);
    
     // delta phi btw leptons
     double dphi = fabs(cms2.hyp_lt_p4()[i_hyp].phi() - cms2.hyp_ll_p4()[i_hyp].phi());
     if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
     hdphiLep->Fill(cuts_passed, myType, dphi, weight);
    
     // Relative isolation... muons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	  hmuRelIso->Fill(cuts_passed, myType, reliso_lt(i_hyp), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	  hmuRelIso->Fill(cuts_passed, myType, reliso_ll(i_hyp), weight);
     }

     // Relative isolation... electrons
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  heleRelIso->Fill(cuts_passed, myType, reliso_lt(i_hyp, true), weight);
	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_lt(i_hyp, false), weight);
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  heleRelIso->Fill(cuts_passed, myType, reliso_ll(i_hyp, true), weight);
	  heleRelIsoTrk->Fill(cuts_passed, myType, reliso_ll(i_hyp, false), weight);
     }

     // lower of the two isolations, regardless of species (used for Dumbo method)
     hminRelIso->Fill(cuts_passed, myType, std::min(reliso_lt(i_hyp), reliso_ll(i_hyp)), weight);
     // lower of the two isolations, regardless of species (used for
     // Dumbo method) --- using calo iso for electrons
     hminRelIso_withCalo->Fill(cuts_passed, myType, 
			      std::min(reliso_lt(i_hyp, true), reliso_ll(i_hyp, true)), 
			      weight);

     // dilepton pt
     hdilPt->Fill(cuts_passed, myType, cms2.hyp_p4()[i_hyp].pt(), weight);
    
     // Met and Met special
     //     hmet->Fill(cuts_passed, myType, cms2.hyp_met()[i_hyp], weight);      
     hmet->Fill(cuts_passed, myType, cms2.evt_metMuonJESCorr(), weight);      
//     hmetSpec->Fill(cuts_passed, myType, 
//		   MetSpecial(cms2.hyp_met()[i_hyp], cms2.hyp_metPhi()[i_hyp], i_hyp),
//		   weight);
     hmetSpec->Fill(cuts_passed, myType, 
		   MetSpecial(cms2.evt_metMuonJESCorr(), cms2.evt_metMuonJESCorrPhi(), i_hyp),
		   weight);
     // track correction to the met
//     const TVector3 trkCorr = correctMETforTracks();
//     TVector3 hyp_met;
//     hyp_met.SetPtEtaPhi(cms2.hyp_met()[i_hyp], 0, cms2.hyp_metPhi()[i_hyp]);
//     hyp_met += trkCorr;
     hmetTrkCorr->Fill(cuts_passed, myType, cms2.evt_tcmet(), weight);

     // tag muon pt and iso
     htagMuPt->Fill(cuts_passed, myType, tagMuonPt(i_hyp), weight);
     htagMuRelIso->Fill(cuts_passed, myType, tagMuonRelIso(i_hyp), weight);

     if (myType == DILEPTON_EMU) {
	  int mu_idx = -1;
	  int el_idx = -1;
	  if (abs(cms2.hyp_lt_id()[i_hyp]) == 13) {
	       mu_idx = cms2.hyp_lt_index()[i_hyp];
	       el_idx = cms2.hyp_ll_index()[i_hyp];
	  }
	  if (abs(cms2.hyp_ll_id()[i_hyp]) == 13) {
	       mu_idx = cms2.hyp_ll_index()[i_hyp];
	       el_idx = cms2.hyp_lt_index()[i_hyp];
	  }
	  assert(el_idx != -1 && mu_idx != -1);
	  hmuPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_id()[mu_idx]), weight);
	  hmuMoPdgId->Fill(cuts_passed, myType, abs(cms2.mus_mc_motherid()[mu_idx]), weight);
	  helPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_id()[el_idx]), weight);
	  helMoPdgId->Fill(cuts_passed, myType, abs(cms2.els_mc_motherid()[el_idx]), weight);
	  helEop      ->Fill(cuts_passed, myType, cms2.els_eOverPIn	()[el_idx], weight);
	  held0	      ->Fill(cuts_passed, myType, fabs(cms2.els_d0		()[el_idx]), weight);
	  if ((cuts_passed & cuts_) == cuts_) {
	       held0vsRelIso->Fill(fabs(cms2.els_d0		()[el_idx]), el_rel_iso(el_idx, true), weight);
	       heldphiinvsRelIso->Fill(cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], el_rel_iso(el_idx, true), weight);
	       if (abs(cms2.els_mc_id()[el_idx]) == 22) {
		    held0vsRelIsoMCgamma->Fill(fabs(cms2.els_d0		()[el_idx]), el_rel_iso(el_idx, true), weight);
		    heldphiinvsRelIsoMCgamma->Fill(cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], el_rel_iso(el_idx, true), weight);
	       }
	  }
	  helfbrem    ->Fill(cuts_passed, myType, cms2.els_fBrem		()[el_idx], weight);
	  helHE       ->Fill(cuts_passed, myType, cms2.els_hOverE	()[el_idx], weight);
	  helsee      ->Fill(cuts_passed, myType, cms2.els_sigmaEtaEta	()[el_idx], weight);
	  if (cuts_passed & CUT_BIT(CUT_EL_BARREL))
	       helsppEB      ->Fill(cuts_passed, myType, cms2.els_sigmaPhiPhi	()[el_idx], weight);
	  else helsppEE      ->Fill(cuts_passed, myType, cms2.els_sigmaPhiPhi	()[el_idx], weight);
	  heldphiin   ->Fill(cuts_passed, myType, cms2.els_charge	()[el_idx] * cms2.els_dPhiIn	()[el_idx], weight);
	  heldetain   ->Fill(cuts_passed, myType, cms2.els_dEtaIn	()[el_idx], weight);
	  helEseedopin->Fill(cuts_passed, myType, cms2.els_eSeedOverPOut	()[el_idx], weight);
// 	  if (cms2.els_mc_id()[el_idx] == 22 && (cuts_passed & cuts) == cuts)
// 	       printf("run %10u, event %10u: weight %f\n", cms2.evt_run(), cms2.evt_event(),
// 		      cms2.evt_scale1fb());
     }

     // coversions
     if (abs(cms2.hyp_lt_id()[i_hyp]) == 11) {
	  const int conv = conversionPartner(cms2.hyp_lt_index()[i_hyp]);
	  if (conv != -1) {
	       const double dphi = conversionDeltaPhi(conv, cms2.hyp_lt_index()[i_hyp]);
	       if (cms2.hyp_lt_id()[i_hyp] * cms2.trks_charge()[conv] > 0)
		    helConvDeltaPhi_os->Fill(cuts_passed, myType, dphi, weight);
	       else helConvDeltaPhi_ss->Fill(cuts_passed, myType, dphi, weight);
	  }
     }
     if (abs(cms2.hyp_ll_id()[i_hyp]) == 11) {
	  const int conv = conversionPartner(cms2.hyp_ll_index()[i_hyp]);
	  if (conv != -1) {
	       const double dphi = conversionDeltaPhi(conv, cms2.hyp_ll_index()[i_hyp]);
	       if (cms2.hyp_ll_id()[i_hyp] * cms2.trks_charge()[conv] > 0)
		    helConvDeltaPhi_os->Fill(cuts_passed, myType, dphi, weight);
	       else helConvDeltaPhi_ss->Fill(cuts_passed, myType, dphi, weight);
	  }
     }
}

void Looper::FillTrilepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // In a trilepton analysis, you would fill your histos here
     //------------------------------------------------------------
}

void Looper::FillQuadlepHistos (int i_hyp)
{
     //------------------------------------------------------------
     // In a quadlepton analysis, you would fill your histos here
     //------------------------------------------------------------
}

void Looper::End ()
{
     //------------------------------------------------------------
     //Example status message at the end of a looper; edit for your
     //application
     //------------------------------------------------------------

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
}
