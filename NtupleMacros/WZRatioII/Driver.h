#include <unistd.h>
#include <string>
#include "Looper.h"
#include "Tools/Sample.h"
#include "LocalSample.h"

#include "Tools/tools.h"

#include "Cuts.h"

using std::string;

enum {
	LOOP_WENU,
	LOOP_EM30_80,
	LOOP_BC30_80,
	
	// 2_1_X
	LOOP_QCD30,
	LOOP_QCD80,
	LOOP_WJET_ALP,
        LOOP_ZEEJET_ALP,
        LOOP_ZMMJET_ALP,
        LOOP_ZTTJET_ALP,

	LOOP_MU15_SINGLE,

	LOOP_Z_0JET
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

     Looper looper_wenu(fWenu(), cuts, log.c_str());
	if (which_ones & (1 << LOOP_WENU)) looper_wenu.Loop();

     Looper looper_bc30_80(fBC30_80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_BC30_80)) looper_bc30_80.Loop();

     Looper looper_em30_80(fEM30_80(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_EM30_80)) looper_em30_80.Loop();

	// 2_2_1
     Looper looper_qcd30(fQCDpt30(), cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD30)) looper_qcd30.Loop();
     Looper looper_qcd80(fQCDpt80(), cuts, log.c_str()); if (which_ones & (1 << LOOP_QCD80)) looper_qcd80.Loop();
     Looper looper_wjet_alp(fWjetsAlpgenSingle(), cuts, log.c_str()); if (which_ones & (1 << LOOP_WJET_ALP)) looper_wjet_alp.Loop();
     Looper looper_zeejet_alp(fZjetsAlpgenSingle(), cuts, log.c_str()); if (which_ones & (1 << LOOP_ZEEJET_ALP)) looper_zeejet_alp.Loop();
     Looper looper_zmmjet_alp(fZmmjetsAlpgenSingle(), cuts, log.c_str());if (which_ones & (1 << LOOP_ZMMJET_ALP)) looper_zmmjet_alp.Loop();
     Looper looper_zttjet_alp(fZttjetsAlpgenSingle(), cuts, log.c_str()); if (which_ones & (1 << LOOP_ZTTJET_ALP)) looper_zttjet_alp.Loop();
     Looper looper_mu15_alp(fInclusiveMuPt15Single(), cuts, log.c_str()); if (which_ones & (1 << LOOP_MU15_SINGLE)) looper_mu15_alp.Loop();


     Looper looper_z_0jet(fZ_0Jet(), cuts, log.c_str());
        if (which_ones & (1 << LOOP_Z_0JET)) looper_z_0jet.Loop();


     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());

     // then we collect them all and print a table
     const Looper *loopers[] = { 
	  &looper_wenu,
	&looper_em30_80,
	&looper_bc30_80,
	&looper_qcd30,
	&looper_qcd80,
	&looper_wjet_alp,
        &looper_zeejet_alp,
        &looper_zmmjet_alp,
        &looper_zttjet_alp,
        &looper_mu15_alp,
	&looper_z_0jet
     };

     printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(event_cuts, "Results", 1 << LOOP_WJET_ALP | 1 << LOOP_QCD30 |
				1 << LOOP_ZEEJET_ALP | 1 << LOOP_ZMMJET_ALP | 1 << LOOP_ZTTJET_ALP);
		//1 << LOOP_WENU);
		// 1 << LOOP_QCD30 | 1 << LOOP_WJET_ALP);
}

