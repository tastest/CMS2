#include <unistd.h>
#include <string>
#include "TDirectory.h"
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

// helper function used to print yield tables
void printTable (const Looper **hists, int n, const char *fname, 
		 uint32 which_ones)
{
/*      FILE *f = 0; */
/*      if (fname == 0 || strlen(fname) == 0) */
/* 	  f = stdin; */
/*      else f = fopen(fname, "w"); */
/*      if (f == 0) { */
/* 	  perror("printing table"); */
/* 	  return; */
/*      } */
/*      fprintf(f, "| %10s", ""); */
/*      for (int j = 0; j < n; ++j) { */
/* 	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str()); */
/*      } */
/*      fprintf(f, "|%30s  |\n", "total"); */
/*      for (int i = 0; i < 4; ++i) { */
/* 	  fprintf(f, "|%10s  ", dilepton_hypo_names[i]); */
/* 	  double cands = 0; */
/* 	  double w2 = 0; */
/* 	  for (int j = 0; j < n; ++j) { */
/* 	       fprintf(f, "|  %10.1f &plusmn; %10.1f",  */
/* 		       hists[j]->CandsPassing(DileptonHypType(i)), */
/* 		       hists[j]->RMS(DileptonHypType(i))); */
/* 	       cands += hists[j]->CandsPassing(DileptonHypType(i)); */
/* 	       w2 += hists[j]->RMS(DileptonHypType(i)) *  */
/* 		    hists[j]->RMS(DileptonHypType(i)); */
/* 	  } */
/* 	  fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2)); */
/*      } */
/*      if (f != stdin)  */
/* 	  fclose(f); */
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

//     Looper looper_psi		(fFreeForm("/data/tmp/cms2/JPsiEE_Summer09-DESIGN_3X_V8A_2360GeV-v1/V03-00-23/merged_ntuple.root"), 0, log.c_str());
     gDirectory = new TDirectory;
     Looper looper_psi                (fFreeForm("/Users/dlevans/Documents/UCSD/psiee/JPsiEE_Summer09-DESIGN_3X_V8A_2360GeV-v1/V03-00-23/merged_ntuple.root"), 0, log.c_str());
     gDirectory = new TDirectory;


     looper_psi.Loop();
     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(baseline_cuts, "Results");
}

int SS_Results ()
{
     return run<Looper>(ss_baseline_cuts, "SS_Results");
}

int In_Zwindow ()
{
     return run<Looper>(baseline_cuts_zwindow, "In_Zwindow");
}
