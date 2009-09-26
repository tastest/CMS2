#include <unistd.h>
#include <string>
#include "Looper.h"
#include "Tools/Sample.h"
#include "LocalSample.h"

#include "Tools/tools.h"

#include "Cuts.h"

using std::string;

enum {

	LOOP_VALIDATION,

	LOOP_WENU,
	LOOP_EM30_80,
	LOOP_BC30_80,

	LOOP_ZZ_7TeV,
	
	// 2_1_X
	LOOP_QCD30,
	LOOP_QCD80,
	LOOP_WJET_ALP
};

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

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";

     Looper looper_validation(fValidation(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_VALIDATION)) looper_validation.Loop();

     Looper looper_wenu(fWenu(), cuts, log.c_str());
	if (which_ones & (1 << LOOP_WENU)) looper_wenu.Loop();

     Looper looper_bc30_80(fBC30_80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_BC30_80)) looper_bc30_80.Loop();

     Looper looper_em30_80(fEM30_80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_EM30_80)) looper_em30_80.Loop();

     Looper looper_ZZ_7TeV(fZZ_7TeV(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_ZZ_7TeV)) looper_ZZ_7TeV.Loop();


	// 2_2_1
     Looper looper_qcd30(fQCDpt30(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_QCD30)) looper_qcd30.Loop();
     Looper looper_qcd80(fQCDpt80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_QCD80)) looper_qcd80.Loop();
     Looper looper_wjet_alp(fWjetsAlpgenSingle(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_WJET_ALP)) looper_wjet_alp.Loop();


     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());

     // then we collect them all and print a table
     const Looper *loopers[] = { 
	&looper_validation,
	  &looper_wenu,
	&looper_em30_80,
	&looper_bc30_80,
	&looper_qcd30,
	&looper_qcd80,
	&looper_wjet_alp,

	&looper_ZZ_7TeV
     };

     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(event_cuts, "Results", 1 << LOOP_ZZ_7TeV);
//LOOP_QCD30 | 1 << LOOP_QCD80 | 1 << LOOP_WJET_ALP);
}

