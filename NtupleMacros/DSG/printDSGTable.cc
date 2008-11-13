#include <stdio.h>
#include "DSGLooper.h"
#include "DileptonHypType.h"

void printDSGTable (const DSGLooperBase **hists, int n, const char *fname, 
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
  fprintf(f, "|%30s  |\n", "SM total");
  for (int i = 0; i < 4; ++i) {
    fprintf(f, "|%10s  ", dilepton_hypo_names[i]);
    double cands = 0;
    double w2 = 0;
    for (int j = 0; j < n; ++j) {
      fprintf(f, "|  %10.1f &plusmn; %10.1f", 
	      hists[j]->CandsPassing(DileptonHypType(i)),
	      hists[j]->RMS(DileptonHypType(i)));
      // only sum up SM samples
      if ( hists[j]->SampleSM() ) {
	cands += hists[j]->CandsPassing(DileptonHypType(i));
	w2 += hists[j]->RMS(DileptonHypType(i)) * 
	  hists[j]->RMS(DileptonHypType(i));
      }
      if (not (which_ones & 1 << j))
	continue;
      const DSGFakeRateLooper *looper = 
	dynamic_cast<const DSGFakeRateLooper *>(hists[j]);
      if (looper != 0) {
#if 0
	fprintf(f, " + %5.1f &minus; %5.1f", 
		looper->CandsPassingSystHi(DileptonHypType(i)) 
		- looper->CandsPassing(DileptonHypType(i)),
		looper->CandsPassing(DileptonHypType(i)) 
		- looper->CandsPassingSystLo(DileptonHypType(i)));
#else
	fprintf(f, "(stat) &plusmn; %5.1f (fake)", 
		looper->FakeSyst(DileptonHypType(i)));
#endif
      }

    }
    fprintf(f, "|  %10.1f &plusmn; %10.1f|\n", cands, sqrt(w2));
  }
  if (f != stdin) 
    fclose(f);
}
