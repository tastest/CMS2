#include <unistd.h>
#include <string>
#include "TDirectory.h"
#include "WWLooper.h"
#include "Sample.h"
#include "utilities.h"
#include "printWWTable.h"

enum {
     LOOP_WW	,
     LOOP_WZ	,
     LOOP_ZZ	,
     LOOP_WJETS	,
     LOOP_DYEE	,
     LOOP_DYMM	,
     LOOP_DYTT	,
     LOOP_TTBAR	,
     LOOP_TW	,
};

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     using std::string;
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     saveHist(hist.c_str());
     const WWLooperBase *loopers[] = { 
	  &looper_ww          ,
	  &looper_wz          ,
	  &looper_zz          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printWWTable(loopers, sizeof(loopers) / sizeof(WWLooperBase *), tbl.c_str(), which_ones);
     return 0;
}

int WW_OldResults_trkjets ()
{
     return run<WWResultsLooper>(ww_old_baseline_cuts | CUT_BIT(WW_PASS_JETVETO_TRACKJETS), "WW_OldResults_trkjets");
}

int WW_OldResults_muveto ()
{
     return run<WWResultsLooper>(ww_old_baseline_cuts | CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT), "WW_OldResults_muveto");
}

int WW_OldResults_caloiso ()
{
     return run<WWResultsLooper>((ww_old_baseline_cuts & ~(CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO))) | 
				 CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO), 	
				 "WW_OldResults_caloiso");
}

int WW_OldResults_caloiso_trkjets ()
{
     return run<WWResultsLooper>((ww_old_baseline_cuts & ~(CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO))) | 
				 CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO) | CUT_BIT(WW_PASS_JETVETO_TRACKJETS), 	
				 "WW_OldResults_caloiso_trkjets");
}

int WW_OldResults_caloiso_muveto ()
{
     return run<WWResultsLooper>((ww_old_baseline_cuts & ~(CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO))) | 
				 CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO) | CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT), 	
				 "WW_OldResults_caloiso_muveto");
}

int WW_Results_softmu ()
{
     return run<WWResultsLooper>((ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO_WITHOUT_PTCUT))) | 
				 CUT_BIT(WW_PASS_MUON_B_VETO), 	
				 "WW_Results_softmu");
}

int WW_NoLLPtCut ()
{
     return run<WWResultsLooper>(ww_baseline_cuts & ~(CUT_BIT(WW_LL_PT)), 	
				 "WW_NoLLPtCut");
}

int WW_Results ()
{
     return run<WWResultsLooper>(ww_baseline_cuts, "WW_Results");
}

int WW_NoMuTag ()
{
     return run<WWResultsLooper>(ww_baseline_no_btags_cuts, "WW_NoMuTag");
}

int WW_NoTrackJets ()
{
     return run<WWResultsLooper>(ww_baseline_no_trackjets_cuts, "WW_NoTrackJets");
}

int WW_NoCaloIso ()
{
     return run<WWResultsLooper>(ww_baseline_no_caloiso_cuts, "WW_NoCaloIso");
}

int WW_OldResults ()
{
     return run<WWResultsLooper>(ww_old_baseline_cuts, "WW_OldResults");
}

int WW_TopEstimate ()
{
     WWMuTagEffLooper looper_eff_ww	(fWW()		, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_ww          .Loop();	
     WWMuTagEffLooper looper_eff_wjets	(fWjets()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_wjets       .Loop();	
     WWMuTagEffLooper looper_eff_dyee	(fDYee()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_dyee        .Loop();	
     WWMuTagEffLooper looper_eff_dymm	(fDYmm()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_dymm        .Loop();	
     WWMuTagEffLooper looper_eff_dytt	(fDYtt()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_dytt        .Loop();	
     WWMuTagEffLooper looper_eff_ttbar	(fttbar()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_ttbar       .Loop();	
     WWMuTagEffLooper looper_eff_tw	(ftW()		, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED));	looper_eff_tw          .Loop();	
     double mu_tag_eff[4];
     memcpy(mu_tag_eff, (looper_eff_ww    + 
			 looper_eff_wjets + 
			 looper_eff_dyee  + 
			 looper_eff_dymm  + 
			 looper_eff_dytt  + 
			 looper_eff_ttbar + 
			 looper_eff_tw    ).MuTagEff(), sizeof(mu_tag_eff)); 
     printf("mu tag efficiency: %f %f %f %f\n", mu_tag_eff[0], mu_tag_eff[1], 
	    mu_tag_eff[2], mu_tag_eff[3]);
     const WWLooperBase *loopers_eff[] = { 
	  &looper_eff_ww          ,
	  &looper_eff_wjets       ,
	  &looper_eff_dyee        ,
	  &looper_eff_dymm        ,
	  &looper_eff_dytt        ,
	  &looper_eff_ttbar       ,
	  &looper_eff_tw          ,
     };
     printWWTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagCount.tbl");
     printWWMuTagTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagEff.tbl");
     histo_directory->Clear();
     WWTopEstimateLooper looper_ww	(fWW()		, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_ww          .Loop();
     WWTopEstimateLooper looper_wjets	(fWjets()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_wjets       .Loop();
     WWTopEstimateLooper looper_dyee	(fDYee()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_dyee        .Loop();
     WWTopEstimateLooper looper_dymm	(fDYmm()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_dymm        .Loop();
     WWTopEstimateLooper looper_dytt	(fDYtt()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_dytt        .Loop();
     WWTopEstimateLooper looper_ttbar	(fttbar()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_ttbar       .Loop();
     WWTopEstimateLooper looper_tw	(ftW()		, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED));	looper_tw          .Loop();
     const WWLooperBase *loopers[] = { 
	  &looper_ww          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printWWTable(loopers, sizeof(loopers) / sizeof(WWResultsLooper *), "WW_TopEstimate.tbl");
     return 0;
}

int WW_TopEstimate_noptcut ()
{
     WWMuTagEffLooper looper_eff_ww	(fWW()		, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_ww          .Loop();	
     WWMuTagEffLooper looper_eff_wjets	(fWjets()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_wjets       .Loop();	
     WWMuTagEffLooper looper_eff_dyee	(fDYee()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dyee        .Loop();	
     WWMuTagEffLooper looper_eff_dymm	(fDYmm()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dymm        .Loop();	
     WWMuTagEffLooper looper_eff_dytt	(fDYtt()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dytt        .Loop();	
     WWMuTagEffLooper looper_eff_ttbar	(fttbar()	, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_ttbar       .Loop();	
     WWMuTagEffLooper looper_eff_tw	(ftW()		, ww_baseline_mu_tageff_cuts, CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_tw          .Loop();	
     double mu_tag_eff[4];
     memcpy(mu_tag_eff, (looper_eff_ww    + 
			 looper_eff_wjets + 
			 looper_eff_dyee  + 
			 looper_eff_dymm  + 
			 looper_eff_dytt  + 
			 looper_eff_ttbar + 
			 looper_eff_tw    ).MuTagEff(), sizeof(mu_tag_eff)); 
     printf("mu tag efficiency: %f %f %f %f\n", mu_tag_eff[0], mu_tag_eff[1], 
	    mu_tag_eff[2], mu_tag_eff[3]);
     const WWLooperBase *loopers_eff[] = { 
	  &looper_eff_ww          ,
	  &looper_eff_wjets       ,
	  &looper_eff_dyee        ,
	  &looper_eff_dymm        ,
	  &looper_eff_dytt        ,
	  &looper_eff_ttbar       ,
	  &looper_eff_tw          ,
     };
     printWWTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagCount_noptcut.tbl");
     printWWMuTagTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagEff_noptcut.tbl");
     histo_directory->Clear();
     WWTopEstimateLooper looper_ww	(fWW()		, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_ww          .Loop();
     WWTopEstimateLooper looper_wjets	(fWjets()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_wjets       .Loop();
     WWTopEstimateLooper looper_dyee	(fDYee()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dyee        .Loop();
     WWTopEstimateLooper looper_dymm	(fDYmm()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dymm        .Loop();
     WWTopEstimateLooper looper_dytt	(fDYtt()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dytt        .Loop();
     WWTopEstimateLooper looper_ttbar	(fttbar()	, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_ttbar       .Loop();
     WWTopEstimateLooper looper_tw	(ftW()		, mu_tag_eff, ww_baseline_cuts & ~(CUT_BIT(WW_PASS_MUON_B_VETO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_tw          .Loop();
     const WWLooperBase *loopers[] = { 
	  &looper_ww          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printWWTable(loopers, sizeof(loopers) / sizeof(WWResultsLooper *), "WW_TopEstimate_noptcut.tbl");
     return 0;
}

int WW_TopEstimate_oldcuts ()
{
     WWMuTagEffLooper looper_eff_ww	(fWW()		, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_ww          .Loop();	
     WWMuTagEffLooper looper_eff_wjets	(fWjets()	, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_wjets       .Loop();	
     WWMuTagEffLooper looper_eff_dyee	(fDYee()	, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dyee        .Loop();	
     WWMuTagEffLooper looper_eff_dymm	(fDYmm()	, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dymm        .Loop();	
     WWMuTagEffLooper looper_eff_dytt	(fDYtt()	, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_dytt        .Loop();	
     WWMuTagEffLooper looper_eff_ttbar	(fttbar()	, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_ttbar       .Loop();	
     WWMuTagEffLooper looper_eff_tw	(ftW()		, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_eff_tw          .Loop();	
     double mu_tag_eff[4];
     memcpy(mu_tag_eff, (looper_eff_ww    + 
			 looper_eff_wjets + 
			 looper_eff_dyee  + 
			 looper_eff_dymm  + 
			 looper_eff_dytt  + 
			 looper_eff_ttbar + 
			 looper_eff_tw    ).MuTagEff(), sizeof(mu_tag_eff)); 
     printf("mu tag efficiency: %f %f %f %f\n", mu_tag_eff[0], mu_tag_eff[1], 
	    mu_tag_eff[2], mu_tag_eff[3]);
     const WWLooperBase *loopers_eff[] = { 
	  &looper_eff_ww          ,
	  &looper_eff_wjets       ,
	  &looper_eff_dyee        ,
	  &looper_eff_dymm        ,
	  &looper_eff_dytt        ,
	  &looper_eff_ttbar       ,
	  &looper_eff_tw          ,
     };
     printWWTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagCount_oldcuts.tbl");
     printWWMuTagTable(loopers_eff, sizeof(loopers_eff) / sizeof(WWResultsLooper *), "WW_TopMuTagEff_oldcuts.tbl");
     histo_directory->Clear();
     WWTopEstimateLooper looper_ww	(fWW()		, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_ww          .Loop();
     WWTopEstimateLooper looper_wjets	(fWjets()	, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_wjets       .Loop();
     WWTopEstimateLooper looper_dyee	(fDYee()	, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dyee        .Loop();
     WWTopEstimateLooper looper_dymm	(fDYmm()	, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dymm        .Loop();
     WWTopEstimateLooper looper_dytt	(fDYtt()	, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_dytt        .Loop();
     WWTopEstimateLooper looper_ttbar	(fttbar()	, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_ttbar       .Loop();
     WWTopEstimateLooper looper_tw	(ftW()		, mu_tag_eff, ww_old_baseline_cuts & ~(CUT_BIT(WW_PASS_JETVETO_CALO)), CUT_BIT(WW_MUON_TAGGED_WITHOUT_PTCUT));	looper_tw          .Loop();
     const WWLooperBase *loopers[] = { 
	  &looper_ww          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printWWTable(loopers, sizeof(loopers) / sizeof(WWResultsLooper *), "WW_TopEstimate_oldcuts.tbl");
     return 0;
}

int WW_SS_Results ()
{
     return run<WWResultsLooper>(ww_ss_baseline_cuts, "WW_SS_Results");
}

int WW_SS_NoCaloIso ()
{
     return run<WWResultsLooper>((ww_ss_baseline_cuts & 
				  ~(CUT_BIT(WW_LT_CALOISO) | CUT_BIT(WW_LL_CALOISO))) | 
				 CUT_BIT(WW_LT_ISO) | CUT_BIT(WW_LL_ISO), 
				 "WW_SS_NoCaloIso");
}

int WW_Wjets_Dumbo ()
{
     return run<WWResultsLooper>(ww_dumbo_cuts, "WW_Wjets_Dumbo", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_Dumbo_nod0 ()
{
     return run<WWResultsLooper>(ww_dumbo_cuts_nod0, "WW_Wjets_Dumbo_nod0", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_SS_Dumbo ()
{
     WWDumboLooper looper_ww		(fWW()		, ww_ss_dumbo_cuts, "WW_SS_Dumbo.log");	looper_ww          .Loop();
     WWDumboLooper looper_wjets		(fWjets()	, ww_ss_dumbo_cuts, "WW_SS_Dumbo.log");	looper_wjets       .Loop();
     saveHist("WWSSDumboResults.root");
     return 0;
}

int WW_Wjets_Numerator ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR), "WW_Wjets_Numerator",
				 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff);
}

int WW_Wjets_Numerator_barrel ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL), 
				 "WW_Wjets_Numerator_barrel", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_Numerator_barrel_supertight ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL)
				 | CUT_BIT(WW_ONE_SUPERTIGHT), 
				 "WW_Wjets_Numerator_barrel_supertight", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_Numerator_barrel_supertight_caloiso ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL)
				 | CUT_BIT(WW_ONE_SUPERTIGHT) | CUT_BIT(WW_TWO_CALOISO), 
				 "WW_Wjets_Numerator_barrel_supertight_caloiso", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_Fakerates ()
{
     return run<WWFakeRateLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR), 
				  "WW_Wjets_Fakerates", 1 | LOOP_WW | 1 << LOOP_WJETS | 0xffffffff );
}

int WW_Wjets_Fakerates_barrel ()
{
     return run<WWFakeRateLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR) | CUT_BIT(WW_EL_BARREL), 
				  "WW_Wjets_Fakerates_barrel", 1 | LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_FOs ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT), 
				 "WW_Wjets_FOs", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_FOs_supertight ()
{
     printf("CUT_BIT(WW_ELFAKE_NOT_NUMERATOR) == %llx\n", cuts_t(1) << 35);
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT)  | CUT_BIT(WW_ONE_SUPERTIGHT), 
				 "WW_Wjets_FOs_supertight", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_FOs_not_numerator ()
{
     return run<WWResultsLooper>(ww_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR), 
				  "WW_Wjets_FOs_not_numerator");
}

int WW_Wjets_SS_Numerator ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR), "WW_Wjets_SS_Numerator", 1 << LOOP_WW | 1 << LOOP_WJETS & 0xffffffff);
}

int WW_Wjets_SS_Numerator_barrel ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL), 
				 "WW_Wjets_SS_Numerator_barrel", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_SS_Numerator_barrel_supertight ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL) | CUT_BIT(WW_ONE_SUPERTIGHT), 
				 "WW_Wjets_SS_Numerator_barrel_supertight", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_SS_Numerator_barrel_supertight_caloiso ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_NUMERATOR) | CUT_BIT(WW_EL_BARREL) | CUT_BIT(WW_ONE_SUPERTIGHT) | CUT_BIT(WW_TWO_CALOISO), 
				 "WW_Wjets_SS_Numerator_barrel_supertight_caloiso", 1 << LOOP_WW | 1 << LOOP_WJETS);
}

int WW_Wjets_SS_Fakerates ()
{
     return run<WWFakeRateLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR), 
				  "WW_Wjets_SS_Fakerates", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff);
}

int WW_Wjets_SS_Fakerates_barrel ()
{
     return run<WWFakeRateLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR) | CUT_BIT(WW_EL_BARREL), 
				  "WW_Wjets_SS_Fakerates_barrel", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_SS_FOs ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT), 
				 "WW_Wjets_SS_FOs", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_SS_FOs_supertight ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ONE_SUPERTIGHT), 
				 "WW_Wjets_SS_FOs_supertight", 1 << LOOP_WW | 1 << LOOP_WJETS | 0xffffffff & 0);
}

int WW_Wjets_SS_FOs_not_numerator ()
{
     return run<WWResultsLooper>(ww_ss_fakerate_cuts | CUT_BIT(WW_ELFAKE_FAKEABLE_OBJECT) | CUT_BIT(WW_ELFAKE_NOT_NUMERATOR), 
				  "WW_Wjets_SS_FOs_not_numerator");
}

int WW_DY_in_emu ()
{
     WWDYInEMuLooper l(fDYmm());
     l.Loop();
     saveHist("WWDYInEMu.root");
     return 0;
}

int WW_TrkCorrMET_NoMetNoZVeto ()
{
     return run<WWResultsLooper>(ww_baseline_cuts_nomet_nozveto, "WW_TrkCorrMET_NoMetNoZVeto", 1 << LOOP_WW | 1 << LOOP_DYEE); //  | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

int WW_TrkCorrMET_NoMet ()
{
     return run<WWResultsLooper>(ww_baseline_cuts_nomet, "WW_TrkCorrMET_NoMet", 1 << LOOP_WW | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

int WW_TrkCorrMET ()
{
     return run<WWResultsLooper>(ww_baseline_cuts 
        	                        & ~(CUT_BIT(WW_PASS2_MET))
	                                & ~(CUT_BIT(WW_PASS4_MET))
					| (CUT_BIT(WW_PASS2_METCORR))
      					| (CUT_BIT(WW_PASS4_METCORR)),
		"WW_TrkCorrMET", 1 << LOOP_WW | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

int WW_DYEstimate_Results ()
{

        // cuts as in WW_Results but remove the Z VETO
     return run<WWDYEstimateLooper>(ww_baseline_cuts
                                & ~(CUT_BIT(WW_PASS_ADDZVETO))
                                & ~(CUT_BIT(WW_PASS2_MET))
                                & ~(CUT_BIT(WW_PASS4_MET))
				& ~(CUT_BIT(WW_PASS_JETVETO_CALO))
     				& ~(CUT_BIT(WW_PASS_JETVETO_TRACKJETS)),
                                 "WW_DYEstimate_Results", 
				1 << LOOP_WW | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}


int WW_In_Zwindow ()
{
     return run<WWResultsLooper>(ww_baseline_cuts_zwindow, "WW_In_Zwindow");
}

int WW_In_Zwindow_TrkCorr ()
{
     return run<WWResultsLooper>(ww_baseline_cuts_zwindow_trkcorr, "WW_In_Zwindow_TrkCorr");
}

