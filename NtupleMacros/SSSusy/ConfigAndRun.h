// -*- C++ -*-

#include <unistd.h>
#include <string>
#include "Looper.h"
#include "Tools/Sample.h"
#include "Tools/tools.h"

using std::string;

// this enum says which samples should actually be used (to shorten
// looping time if you only care about the yields for one or two
// samples)
enum {
     LOOP_WW	,
     LOOP_WW_EXCL	,
     LOOP_WZ	,
     LOOP_ZZ	,
     LOOP_WJETS	,
     LOOP_WJETS_AND_FRIENDS	,
     LOOP_DYEE	,
     LOOP_DYMM	,
     LOOP_DYTT	,
     LOOP_DY_AND_FRIENDS	,
     LOOP_WGAMMA,
     LOOP_ZGAMMA,
     LOOP_TTBAR	,
     LOOP_TTBAR_TAUOLA	,
     LOOP_TW	,
     LOOP_TW_AND_FRIENDS	,
     LOOP_LM0	,
     LOOP_LM1	,
     LOOP_LM2	,
     LOOP_LM3	,
     LOOP_LM4	,
     LOOP_LM5	,
     LOOP_LM6	,
     LOOP_LM7	,
     LOOP_LM8	,
     LOOP_LM9	,
     LOOP_LM10	,
     LOOP_LM11	,
};

uint32 default_samples = (1 <<      LOOP_WW)	|
     (1 << LOOP_WZ	)	|
     (1 << LOOP_ZZ	)	|
     (1 << LOOP_WJETS	)	|
     (1 << LOOP_DYEE	)	|
     (1 << LOOP_DYMM	)	|
     (1 << LOOP_DYTT	)	|
     (1 << LOOP_TTBAR	)	|
     (1 << LOOP_TW	)       |
     (1 << LOOP_LM0     )       |
     (1 << LOOP_LM1     )       |
     (1 << LOOP_LM2     )       |
     (1 << LOOP_LM3     )       |
     (1 << LOOP_LM4     )       |
     (1 << LOOP_LM5     )       |
     (1 << LOOP_LM6     )       |
     (1 << LOOP_LM7     )       |
     (1 << LOOP_LM8     )       |
     (1 << LOOP_LM9     )       |
     (1 << LOOP_LM10    )       |
     (1 << LOOP_LM11    )      
  ;

uint32 eff_samples = default_samples | (1 << LOOP_WW_EXCL);
// uint32 eff_samples = (default_samples & ~(1 << LOOP_WW)) | (1 << LOOP_WW_EXCL);

// #define TWIKI_OUTPUT
#define LATEX_OUTPUT
//#define SUMMARY_OUTPUT

// helper function used to print yield tables
void printTable (const Looper **hists, int n, const char *fname, 
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
#if defined(TWIKI_OUTPUT)
     fprintf(f, "| %10s", "");
     for (int j = 0; j < n; ++j) {
	  if (not hists[j]->HasRun())
	       continue;
	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
     }
     fprintf(f, "|*%30s*  |\n", "total bg");
#else 
#if defined(LATEX_OUTPUT)
     fprintf(f, "\\hline\\hline\n%10s", "");
     for (int j = 0; j < n; ++j) {
	  if (not hists[j]->HasRun())
	       continue;
	  fprintf(f, "&  \\%-30s  ", hists[j]->SampleName().c_str());
     }
     fprintf(f, "\\\\\\hline\n");
#endif
#endif
     for (int i = 0; i < 4; ++i) {
#if defined(TWIKI_OUTPUT)
	  fprintf(f, "|%10s  ", dilepton_hypo_names[i]);
#else 
#if defined(LATEX_OUTPUT)
	  fprintf(f, "\\%-10s  ", dilepton_hypo_names[i]);
#endif
#endif
	  double cands = 0;
	  double w2 = 0;
	  int is_background = 0;
	  for (int j = 0; j < n; ++j) {
	       if (not hists[j]->HasRun())
		    continue;
#if defined(TWIKI_OUTPUT)
	       fprintf(f, "|  %10.1f &plusmn; %10.1f", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
#else 
#if defined(LATEX_OUTPUT)
	       fprintf(f, "&  %10.1f $\\pm$ %10.1f", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
#endif
#endif
	       if (is_background) {
		    cands += hists[j]->CandsPassing(DileptonHypType(i));
		    w2 += hists[j]->RMS(DileptonHypType(i)) * 
			 hists[j]->RMS(DileptonHypType(i));
	       }
	       is_background++;
#if 0
	       const FakeRateLooper *looper = 
		    dynamic_cast<const FakeRateLooper *>(hists[j]);
	       if (looper != 0) {
		    fprintf(f, "(stat) &plusmn; %5.1f (fake)", 
			    looper->FakeSyst(DileptonHypType(i)));
	       }
#endif
	  }
#if defined(TWIKI_OUTPUT)
	  fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
#else
#if defined(LATEX_OUTPUT)
	  fprintf(f, "\\\\\n");
#endif
#endif
     }
#if defined(LATEX_OUTPUT)
	  fprintf(f, "\\hline\\hline\n");
#endif
     if (f != stdin) 
	  fclose(f);
}

void printTableVertically (const Looper **hists, int n, const char *fname, 
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
     double cands[4] = {0, 0, 0, 0};
     double w2[4] = {0, 0, 0, 0};
     int is_background[4] = {0, 0, 0, 0};
#if defined(LATEX_OUTPUT) && !defined(SUMMARY_OUTPUT)
     fprintf(f, "\\hline\\hline\n");
#endif
     for (int j = 0; j < n; ++j) {
	  if (not hists[j]->HasRun())
	       continue;
#if defined(TWIKI_OUTPUT)
	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
#else 
#if defined(LATEX_OUTPUT)
	  fprintf(f, "\\%-30s  ", hists[j]->SampleName().c_str());
#endif
#endif
	  for (int i = 0; i < 4; ++i) {
#if defined(TWIKI_OUTPUT)
	       fprintf(f, "|  %10.1f &plusmn; %10.1f  |\n", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
#else 
#if defined(LATEX_OUTPUT)
	       const double n = hists[j]->CandsPassing(DileptonHypType(i));
#if defined(SUMMARY_OUTPUT)
	       if (n < 1000000)
		    fprintf(f, "& %18.1f", n);
	       else {
		    const double expo = log(n) / log(10);
		    const int log10 = (int)floor(expo);
		    const int exp10 = pow(10.0, log10);
		    fprintf(f, "& %3.1f$\\cdot$10$^{%d}$\\\\\n", n / exp10, log10);
	       }
#else
	       fprintf(f, "&  %10.2f & %5.2f", 
		       hists[j]->CandsPassing(DileptonHypType(i)),
		       hists[j]->RMS(DileptonHypType(i)));
#endif
#endif
#endif
	       if (is_background[i]) {
		    cands[i] += hists[j]->CandsPassing(DileptonHypType(i));
		    w2[i] += hists[j]->RMS(DileptonHypType(i)) * 
			 hists[j]->RMS(DileptonHypType(i));
	       }
	       is_background[i]++;
#if 0
	       const FakeRateLooper *looper = 
		    dynamic_cast<const FakeRateLooper *>(hists[j]);
	       if (looper != 0) {
		    fprintf(f, "(stat) &plusmn; %5.1f (fake)", 
			    looper->FakeSyst(DileptonHypType(i)));
	       }
#endif	       
	  }
#if defined(LATEX_OUTPUT) 
	  fprintf(f, "\\\\\n", n);
#endif
     }
#if defined(TWIKI_OUTPUT)
     fprintf(f, "|  *%30s*  ", "total");
     for (int i = 0; i < 4; ++i) {
	  fprintf(f, "|  %10.1f &plusmn; %10.1f  |\n", cands[i], sqrt(w2[i]));
     }
#else
#if defined(LATEX_OUTPUT)
     fprintf(f, "\\hline\n %-30s  ", "total");
     for (int i = 0; i < 4; ++i) {
	  const double n = cands[i];
#if defined(SUMMARY_OUTPUT)
	  if (n < 1000000)
	       fprintf(f, "& %18.1f", n);
	  else {
	       const double expo = log(n) / log(10);
	       const int log10 = (int)floor(expo);
	       const int exp10 = pow(10.0, log10);
	       fprintf(f, "& %3.1f$\\cdot$10$^{%d}$\\\\\n", n / exp10, log10);
	  }
#else
	  fprintf(f, "&  %10.2f & %5.2f", cands[i], sqrt(w2[i]));
#endif
     }
     fprintf(f, "\\\\\n", n);
#endif
#endif
#if defined(LATEX_OUTPUT) && !defined(SUMMARY_OUTPUT)
     fprintf(f, "\\hline\\hline\n");
#endif
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
template <class L> int run (cuts_t cuts, const string &name, uint32 which_ones = default_samples,
			    void (*print)(const Looper **, int, const char *, uint32) = printTableVertically)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     L looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
//      L looper_ww_excl		(fWW_excl()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW_EXCL    )) looper_ww_excl          .Loop();
     L looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
//      L looper_wz_incl		(fWZ_incl()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz_incl     .Loop();
     L looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     L looper_wjetsAlpgen	(fWjetsAlpgenSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjetsAlpgen       .Loop();
     L looper_wjets		(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     //     L looper_wjets		(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     L looper_wc		(fWc()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS_AND_FRIENDS )) looper_wc       .Loop();
     L looper_vlqq		(fVlqq()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS_AND_FRIENDS )) looper_vlqq       .Loop();
     L looper_dyee		(fDY20ee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     L looper_dymm		(fDY20mm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     L looper_dytt		(fDY20tt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     L looper_astar		(fAstar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY_AND_FRIENDS  )) looper_astar        .Loop();
     L looper_dy20tt		(fDY20tt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY_AND_FRIENDS  )) looper_dy20tt        .Loop();
     L looper_dy20mm		(fDY20mm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY_AND_FRIENDS  )) looper_dy20mm        .Loop();
//      L looper_wgamma		(fWgamma()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WGAMMA  )) looper_wgamma        .Loop();
//      L looper_zgamma		(fZgamma()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZGAMMA  )) looper_zgamma        .Loop();
     L looper_ttbar		(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     L looper_ttbar_tauola	(fttbar_taula()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR_TAUOLA )) looper_ttbar_tauola.Loop();
     L looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     L looper_singletop_tchan	(fSingleTop_tChannel()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW_AND_FRIENDS    )) looper_singletop_tchan          .Loop();
     L looper_singletop_schan	(fSingleTop_sChannel()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW_AND_FRIENDS    )) looper_singletop_schan          .Loop();
     L looper_lm0		(fLM0()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM0   )) looper_lm0         .Loop();
     L looper_lm1		(fLM1()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM1   )) looper_lm1         .Loop();
     L looper_lm2		(fLM2()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM2   )) looper_lm2         .Loop();
     L looper_lm3		(fLM3()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM3   )) looper_lm3         .Loop();
     L looper_lm4		(fLM4()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM4   )) looper_lm4         .Loop();
     L looper_lm5		(fLM5()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM5   )) looper_lm5         .Loop();
     L looper_lm6		(fLM6()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM6   )) looper_lm6         .Loop();
     L looper_lm7		(fLM7()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM7   )) looper_lm7         .Loop();
     L looper_lm8		(fLM8()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM8   )) looper_lm8         .Loop();
     L looper_lm9		(fLM9()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM9   )) looper_lm9         .Loop();
     L looper_lm10		(fLM10()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM10  )) looper_lm10        .Loop();
     L looper_lm11		(fLM11()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_LM11  )) looper_lm11        .Loop();
     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     // then we collect them all and print a table
     const Looper *loopers[] = { 
	  &looper_ww          ,
// 	  &looper_ww_excl     ,
	  &looper_wz          ,
// 	  &looper_wz_incl     ,
	  &looper_zz          ,
 	  &looper_wjetsAlpgen ,
 	  &looper_wjets       ,
 	  &looper_wc       ,
 	  &looper_vlqq       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_astar        ,
	  &looper_dy20tt        ,
	  &looper_dy20mm        ,
// 	  &looper_wgamma        ,
// 	  &looper_zgamma        ,
	  &looper_ttbar       ,
	  &looper_ttbar_tauola,
	  &looper_tw          ,
	  &looper_singletop_tchan          ,
	  &looper_singletop_schan          ,
          &looper_lm0         ,
          &looper_lm1         ,
          &looper_lm2         ,
          &looper_lm3         ,
          &looper_lm4         ,
          &looper_lm5         ,
          &looper_lm6         ,
          &looper_lm7         ,
          &looper_lm8         ,
          &looper_lm9         ,
          &looper_lm10        ,
          &looper_lm11        ,
     };
     print(loopers, sizeof(loopers) / sizeof(L *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(baseline_cuts, "Results");
}

int  OSSUSY()
{
     return run<Looper>(baselineOSSUSY_cuts, "OSSUSY");
}

int  SSSUSY()
{
     return run<Looper>(baselineSSSUSY_cuts, "SSSUSY");
}

int Wjets_Numerator ()
{
     return run<Looper>(fakerate_numerator_cuts, "Wjets_Numerator");
}

int Wjets_FOs_Not_Numerator ()
{
     return run<Looper>(fakerate_denominator_not_numerator_cuts, "Wjets_FOs_Not_Numerator");
}

int Wjets_Fakerate ()
{
  // had this crap til 090715 18:24
//      return run<FakeRateLooper>(baseline_cuts & 
// 				~(CUT_BIT(CUT_PASS_TRIGGER) |
// 				  CUT_BIT(CUT_LT_GOOD) | CUT_BIT(CUT_LL_GOOD) |
// 				  CUT_BIT(CUT_LT_ISO) | CUT_BIT(CUT_LL_ISO)), 
// 				"Wjets_Fakerate");
  return run<FakeRateLooper>(fakerate_denominator_not_numerator_cuts,
				"Wjets_Fakerate");
}

int Wjets_SS_Numerator ()
{
  //  return run<Looper>(fakerate_ss_numerator_cuts, "Wjets_SS_Numerator", 1 << LOOP_TTBAR);
  return run<Looper>(fakerate_ss_numerator_cuts, "Wjets_SS_Numerator");
}

int Wjets_SS_FOs_Not_Numerator ()
{
     return run<Looper>(fakerate_ss_denominator_not_numerator_cuts, "Wjets_SS_FOs_Not_Numerator");
}

int Wjets_SS_Fakerate ()
{
  //     return run<FakeRateLooper>(fakerate_ss_denominator_not_numerator_cuts, "Wjets_SS_Fakerate", 1 << LOOP_TTBAR);
     return run<FakeRateLooper>(fakerate_ss_denominator_not_numerator_cuts, "Wjets_SS_Fakerate");
}
