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
     LOOP_WZ	,
     LOOP_ZZ	,
     LOOP_WJETS	,
     LOOP_DYEE	,
     LOOP_DYMM	,
     LOOP_DYTT	,
     LOOP_TTBAR	,
     LOOP_TW	,
};

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
template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     Looper looper_ttbar	(fttbarSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
	 //Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
	 //Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
	 //Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_dyee		(fZeejetsAlpgenSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm		(fZmmjetsAlpgenSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dytt		(fZttjetsAlpgenSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_wjets	(fWjetsAlpgenSingle()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
	 //Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     // then we collect them all and print a table
     const Looper *loopers[] = { 
// 	  &looper_ww          ,
// 	  &looper_wz          ,
// 	  &looper_zz          ,
	  &looper_ttbar       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
 	  &looper_dytt        ,
	  &looper_wjets       ,
// 	  &looper_tw          ,
     };
     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results_se ()
{
  return run<Looper>(baseline_single_electron_cuts, "Results_se", 1 << LOOP_WJETS | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT | 1 << LOOP_TTBAR);
}

int Results_sm ()
{
  return run<Looper>(baseline_single_muon_cuts, "Results_sm", 1 << LOOP_WJETS | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT | 1 << LOOP_TTBAR);
}

int Results_de ()
{
  return run<Looper>(baseline_dielectron_cuts, "Results_de", 1 << LOOP_WJETS | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT | 1 << LOOP_TTBAR);
}

int Results_dm ()
{
  return run<Looper>(baseline_dimuon_cuts, "Results_dm", 1 << LOOP_WJETS | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT | 1 << LOOP_TTBAR);
}

int Results_em ()
{
  return run<Looper>(baseline_emu_cuts, "Results_em", 1 << LOOP_WJETS | 1 << LOOP_DYEE | 1 << LOOP_DYMM | 1 << LOOP_DYTT | 1 << LOOP_TTBAR);
}

