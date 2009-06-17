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
	  if (not (which_ones & 1 << j))
	       continue;
	  fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str());
     }
     fprintf(f, "|%30s  |\n", "total");
#else 
#if defined(LATEX_OUTPUT)
     fprintf(f, "\\hline\\hline\n%10s", "");
     for (int j = 0; j < n; ++j) {
	  if (not (which_ones & 1 << j))
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
	  for (int j = 0; j < n; ++j) {
	       if (not (which_ones & 1 << j))
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
	       cands += hists[j]->CandsPassing(DileptonHypType(i));
	       w2 += hists[j]->RMS(DileptonHypType(i)) * 
		    hists[j]->RMS(DileptonHypType(i));
	       const FakeRateLooper *looper = 
		    dynamic_cast<const FakeRateLooper *>(hists[j]);
	       if (looper != 0) {
		    fprintf(f, "(stat) &plusmn; %5.1f (fake)", 
			    looper->FakeSyst(DileptonHypType(i)));
	       }
	       
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
template <class L> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     L looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     L looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     L looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     L looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     L looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     L looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     L looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     L looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     L looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     // when all the loopers are done, we save the histograms to file
     saveHist(hist.c_str());
     // then we collect them all and print a table
     const Looper *loopers[] = { 
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
     printTable(loopers, sizeof(loopers) / sizeof(L *), tbl.c_str(), which_ones);
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(baseline_cuts, "Results");
}

int Feb_Results ()
{
     return run<Looper>(feb_baseline_cuts, "Feb_Results");
}

int Feb_Results_With_Ntrks ()
{
     return run<Looper>(feb_baseline_with_ntrks_cuts, "Feb_Results_With_Ntrks");
}

int Feb_Results_With_TcMET ()
{
     return run<Looper>(feb_baseline_with_tcmet_cuts, "Feb_Results_With_TcMET");
}

int Feb_Results_With_TrackJets ()
{
     return run<Looper>(feb_baseline_with_trackjets_cuts, "Feb_Results_With_TrackJets");
}

int Feb_Results_With_Btags ()
{
     return run<Looper>(feb_baseline_with_btags_cuts, "Feb_Results_With_Btags");
}

int Feb_Results_With_CaloIso ()
{
     return run<Looper>(feb_baseline_with_caloiso_cuts, "Feb_Results_With_CaloIso");
}

int Feb_Results_NoPass4MET ()
{
     return run<Looper>(feb_baseline_no_pass4met_cuts, "Feb_Results_NoPass4MET");
}

int Oct_Results ()
{
     return run<Looper>(oct_baseline_cuts, "Oct_Results");
}

int Results_NoTcMET ()
{
     return run<Looper>(baseline_no_tcmet_cuts, "Results_NoTcMET");
}

int Results_NoTrackJets ()
{
     return run<Looper>(baseline_no_trackjets_cuts, "Results_NoTrackJets", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int Results_NoSipJets ()
{
     return run<Looper>(baseline_no_sipjets_cuts, "Results_NoSipJets", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int Results_NoBtags ()
{
     return run<Looper>(baseline_no_btags_cuts, "Results_NoBtags");
}

int Results_NoCaloIso ()
{
     return run<Looper>(baseline_no_caloiso_cuts, "Results_NoCaloIso");
}

int Results_NoPass4MET ()
{
     return run<Looper>(baseline_no_pass4met_cuts, "Results_NoPass4MET");
}

int Results_NoNtrks ()
{
     return run<Looper>(baseline_no_ntrks_cuts, "Results_NoNtrks");
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
     return run<FakeRateLooper>(fakerate_denominator_not_numerator_cuts, "Wjets_Fakerate");
}

int Sip_Results ()
{
  return run<Looper>(sip_cuts, "Sip_Results", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int Sip_PlusCaloJetVeto_Results ()
{
  return run<Looper>(sip_plus_calojetveto_cuts, "Sip_PlusCaloJetVeto_Results", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int JetVetoSip ()
{
  return run<Looper>(baseline_plus_sip_cuts, "JetVetoSip", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int KtResults ()
{
     return run<KtSipLooper>(baseline_no_trackjets_cuts & ~CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT), "KtResults", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int ConeResults ()
{
     return run<Looper>(baseline_no_trackjets_cuts & ~CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT), "ConeResults", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int ConeWithSipResults ()
{
     return run<Looper>((baseline_no_trackjets_cuts & ~CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | CUT_BIT(CUT_PASS_JETVETO_SIP), 
			"ConeWithSipResults", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

int KtWithSipResults ()
{
     return run<KtSipLooper>((baseline_no_trackjets_cuts & ~CUT_BIT(CUT_PASS_MUON_B_VETO_WITHOUT_PTCUT)) | CUT_BIT(CUT_PASS_JETVETO_SIP), 
			     "KtWithSipResults", 1 << LOOP_WW | 1 << LOOP_TTBAR | 1 << LOOP_TW);
}

