#include <unistd.h>
#include <string>
#include "DYEst.h"
#include "Tools/Sample.h"
#include "Tools/tools.h"

using std::string;

// this enum says which samples should actually be used (to shorten
// looping time if you only care about the yields for one or two
// samples)
enum {
     LOOP_WW	,
     LOOP_WZ	,
     LOOP_ZZ	,
     LOOP_WWINCL    , 
     LOOP_WZINCL    ,
     LOOP_ZZINCL    ,
     LOOP_WJETS	,
     LOOP_DYEE	,
     LOOP_DYMM	,
     LOOP_DY20MM,
     LOOP_DYTT	,
     LOOP_TTBAR	,
     LOOP_TW	,
};

// helper function used to print yield tables
void printTable (const DYEst **hists, int n, const char *fname, 
		 uint32 which_ones)
{
     FILE *f = 0;
     if (fname == 0 || strlen(fname) == 0)
	  f = stdin;
     else f = fopen(fname, "w");
     if (f == 0) {
	  perror("printing table");
	  return;
     }
     fprintf(f, "| %10s", "");
     for (int j = 0; j < n; ++j) {
	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
     }
     fprintf(f, "|%30s  |\n", "total");
     for (int i = 0; i < 4; ++i) {
	  fprintf(f, "|%10s  ", dilepton_hypo_names[i]);
	  double cands = 0;
	  double w2 = 0;
	  for (int j = 0; j < n; ++j) {
	       fprintf(f, "|  %10.1f &plusmn; %10.1f", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
	       cands += hists[j]->CandsPassing(DileptonHypType(i));
	       w2 += hists[j]->RMS(DileptonHypType(i)) * 
		    hists[j]->RMS(DileptonHypType(i));
	  }
	  fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
     }
     if (f != stdin) 
	  fclose(f);
}

// run a looper on each sample and produce a yield table; arguments:
//
// class Looper: which type of looper to run (usually: Looper)
// cuts: cut definition from Looper.h (usually: baseline_cuts)
// name: name for the output files (usually: "Results", which produces Results.tbl, Results.root, Results.log)
// which_ones: which samples to run (usually: all); to run only WW and ttbar, use: (1 << LOOP_WW) | (1 << LOOP_TTBAR)
//
// examples:
// run<Looper>(baseline_cuts, "Results", 1 << LOOP_WW)				// produce table with default cuts, WW only
// run<Looper>(baseline_cuts, "Results", 1 << LOOP_WW | 1 << LOOP_WJETS)	// produce table with default cuts, WW and Wjets only
// run<Looper>(baseline_cuts, "Results")					// produce table with default cuts, all samples
template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();

     Looper looper_wwincl       (fWW_incl()     , cuts, log.c_str());   if (which_ones & (1 << LOOP_WWINCL)) looper_wwincl      .Loop();
     Looper looper_wzincl       (fWZ_incl()     , cuts, log.c_str());   if (which_ones & (1 << LOOP_WZINCL)) looper_wzincl      .Loop();
     Looper looper_zzincl       (fZZ_incl()     , cuts, log.c_str());   if (which_ones & (1 << LOOP_ZZINCL)) looper_zzincl      .Loop();

     Looper looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dy20mm       (fDY20mm()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_DY20MM)) looper_dy20mm      .Loop();

     Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();

     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     // then we collect them all and print a table
     const Looper *loopers[] = { 
	  &looper_ww          ,
	  &looper_wz          ,
	  &looper_zz          ,
          &looper_wwincl      ,
          &looper_wzincl      ,
          &looper_zzincl      ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
     };
     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
     return 0;
}

int DYEstResults_ForWW_MET20 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE20))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET20");
}

int DYEstResults_ForWW_MET25 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE25))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET25");
}

int DYEstResults_ForWW_MET30 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE30))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET30");
}


int DYEstResults_ForWW_MET35 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE35))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET35");
}


int DYEstResults_ForWW_MET40 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE40))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET35");
}

int DYEstResults_ForWW_MET45 ()
{
        // make sure that the 'pass5' is removed and
        // the same set of cuts is applied to each hyp_type
        // as far as the met is concerned (this is needed
        // to get the non peaking bg estimate from emu
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE45))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET45",
                        1 << LOOP_WW | 1 << LOOP_WZ | 1 << LOOP_ZZ |
                        1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
                        1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

int DYEstResults_ForWW_MET45_DY20MM ()
{
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE45))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET45_DY20MM",
//                        1 << LOOP_WWINCL | 1 << LOOP_WZINCL | 1 << LOOP_ZZINCL |
//                        1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
//                        1 << LOOP_DYEE | 1 << LOOP_DY20MM | 1 << LOOP_DYTT);
//			  1 << LOOP_WWINCL | 
			  1 << LOOP_DY20MM);
}

// use the inclusive WW, ZZ and WZ samples
int DYEstResults_ForWW_MET45_INCL ()
{
     return run<DYEst>(baseline_cuts_nomet
                                     | (CUT_BIT(CUT_MET_SIMPLE45))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET45_INCL",
                        1 << LOOP_WWINCL | 1 << LOOP_WZINCL | 1 << LOOP_ZZINCL |
                        1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
                        1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

// use the inclusive WW, ZZ and WZ samples 
// AND the trigger selection
int DYEstResults_ForWW_MET45_INCL_TRIG ()
{
     return run<DYEst>(baseline_cuts_nomet
//					| (CUT_BIT(CUT_PASS_TRIGGER))
                                     | (CUT_BIT(CUT_MET_SIMPLE45))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)),
                        "DYEstResults_ForWW_MET45_INCL_TRIG",
                        1 << LOOP_WWINCL | 1 << LOOP_WZINCL | 1 << LOOP_ZZINCL |
                        1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
                        1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}


int DYEstResults_ForWW22X_ValStandard ()
{
	// standard new baseline cuts for 22X
     return run<DYEst>(baseline_cuts,
                        "DYEstResults_ForWW22X_ValStandard");
                        //1 << LOOP_WW);// | 1 << LOOP_TTBAR);

}

int DYEstResults_ForWW22X_ValPass5 ()
{
        // standard new baseline cuts for 22X
     return run<DYEst>(baseline_cuts_pass5,
                        "DYEstResults_ForWW22X_ValPass5");
                      //  1 << LOOP_WW);// | 1 << LOOP_TTBAR);

}

// needed to get the difference in efficiency
// between e and mu without any met cuts applied
// for use in estimating the background with the emu
// final state.
int DYEstResults_GetEMuEff ()
{
     return run<DYEst>(baseline_cuts_pass5
                                     & ~(CUT_BIT(CUT_PASS5_MET)),
                        "DYEstResults_GetEMuEff");

}

// inclusive samples 
// AND the trigger
int DYEstResults_GetEMuEff_TRIG ()
{
     return run<DYEst>(		(baseline_cuts_pass5
                                        | (CUT_BIT(CUT_PASS_TRIGGER)))
                                     & ~(CUT_BIT(CUT_PASS5_MET)),
                        "DYEstResults_GetEMuEff_TRIG",
                        1 << LOOP_WWINCL | 1 << LOOP_WZINCL | 1 << LOOP_ZZINCL |
                        1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
                        1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

// use the inclusive WW, ZZ and WZ samples 
// AND the trigger selection
int DYEstResults_ForWW_MET45_INCL_TRIG_WZStudy ()
{
     return run<DYEst>((baseline_cuts_nomet
//                                      | (CUT_BIT(CUT_PASS_TRIGGER))
                                     | (CUT_BIT(CUT_MET_SIMPLE45))
                                     | (CUT_BIT(CUT_MET_BALLANCE))
                                     | (CUT_BIT(CUT_MET_PROJECTED)))

					& ~(CUT_BIT(CUT_PASS_ZVETO)),


                        "DYEstResults_ForWW_MET45_INCL_TRIG_WZStudy",
                        //1 << LOOP_WWINCL | 
			1 << LOOP_WZINCL);
			//| 1 << LOOP_ZZINCL |
                        //1 << LOOP_TTBAR | 1 << LOOP_WJETS | 1 << LOOP_TW |
                        //1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT);
}

// test
int test ()
{
     return run<DYEst>(baseline_cuts_pass5, "test");

}


