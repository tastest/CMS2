// -*- C++ -*-

#include <unistd.h>
#include <string>
#include "Looper.h"
#include "FakeRateLooper.h"
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
  LOOP_InclusiveMu5Pt50,
  LOOP_InclusiveMuPt15,
  LOOP_QCDBCtoEPt20to30,
  LOOP_QCDBCtoEPt30to80,
  LOOP_QCDBCtoEPt80to170,
  LOOP_QCDEMenrichedPt20to30,
  LOOP_QCDEMenrichedPt30to80,
  LOOP_QCDEMenrichedPt80to170,
  LOOP_QCDpt30,
  LOOP_QCDpt30to80,
  LOOP_QCDpt80to170,
  LOOP_QCDpt170to300,
  LOOP_QCDpt300to470,
  LOOP_QCDpt470to800,
  LOOP_QCDpt800toInf,
};

// helper function used to print yield tables
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
  for (int j = 0; j < n; ++j) { 
    fprintf(f, "|  *%30s*  ", hists[j]->SampleName().c_str()); 
  } 
  fprintf(f, "|%30s  |\n", "total"); 
  for (int j = 0; j < n; ++j) { 
    fprintf(f, "|  %10.1f &plusmn; %10.1f",  
            hists[j]->CandsPassing(), 
            hists[j]->RMS()); 
  } 
  fprintf(f, "\n"); 
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
  L looper_InclusiveMu5Pt50		(fInclusiveMu5Pt50()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_InclusiveMu5Pt50          )) looper_InclusiveMu5Pt50       .Loop();
  L looper_InclusiveMuPt15		        (fInclusiveMuPt15()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_InclusiveMuPt15           )) looper_InclusiveMuPt15        .Loop();
  L looper_QCDBCtoEPt20to30		(fQCDBCtoEPt20to30()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDBCtoEPt20to30          )) looper_QCDBCtoEPt20to30       .Loop();
  L looper_QCDBCtoEPt30to80		(fQCDBCtoEPt30to80()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDBCtoEPt30to80          )) looper_QCDBCtoEPt30to80       .Loop();
  L looper_QCDBCtoEPt80to170		(fQCDBCtoEPt80to170()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDBCtoEPt80to170         )) looper_QCDBCtoEPt80to170      .Loop();
  L looper_QCDEMenrichedPt20to30		(fQCDEMenrichedPt20to30()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDEMenrichedPt20to30     )) looper_QCDEMenrichedPt20to30  .Loop();
  L looper_QCDEMenrichedPt30to80		(fQCDEMenrichedPt30to80()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDEMenrichedPt30to80     )) looper_QCDEMenrichedPt30to80  .Loop();
  L looper_QCDEMenrichedPt80to170	        (fQCDEMenrichedPt80to170()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_QCDEMenrichedPt80to170    )) looper_QCDEMenrichedPt80to170 .Loop();
  L looper_QCDpt30          (fQCDpt30()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt30    )) looper_QCDpt30 .Loop();
  L looper_QCDpt30to80          (fQCDpt30to80()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt30to80    )) looper_QCDpt30to80 .Loop();
  L looper_QCDpt80to170          (fQCDpt80to170()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt80to170    )) looper_QCDpt80to170 .Loop();
  L looper_QCDpt170to300          (fQCDpt170to300()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt170to300    )) looper_QCDpt170to300 .Loop();
  L looper_QCDpt300to470          (fQCDpt300to470()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt300to470    )) looper_QCDpt300to470 .Loop();
  L looper_QCDpt470to800          (fQCDpt470to800()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt470to800    )) looper_QCDpt470to800 .Loop();
  L looper_QCDpt800toInf          (fQCDpt800toInf()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt800toInf    )) looper_QCDpt800toInf .Loop();
  // when all the loopers are done, we save the histograms to file
  saveHist(hist.c_str());
  // then we collect them all and print a table
  const Looper *loopers[] = { 
    &looper_InclusiveMu5Pt50,
    &looper_InclusiveMuPt15,
    &looper_QCDBCtoEPt20to30,
    &looper_QCDBCtoEPt30to80,
    &looper_QCDBCtoEPt80to170,
    &looper_QCDEMenrichedPt20to30,
    &looper_QCDEMenrichedPt30to80,
    &looper_QCDEMenrichedPt80to170,
    &looper_QCDEMenrichedPt80to170,
    &looper_QCDpt30,
    &looper_QCDpt30to80,
    &looper_QCDpt80to170,
    &looper_QCDpt170to300,
    &looper_QCDpt300to470,
    &looper_QCDpt470to800,
    &looper_QCDpt800toInf,
  };
  printTable(loopers, sizeof(loopers) / sizeof(L *), tbl.c_str(), which_ones);
  return 0;
}

int Observation ()
{
  return run<Looper>(observation_cuts, "Observation", 1 << LOOP_QCDpt30);
}

int ObservationQCDBins ()
{
  return run<Looper>(observation_cuts_qcd_bins, "ObservationQCDBins", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

int Prediction ()
{
  return run<FakeRateLooper>(prediction_cuts, "Prediction", 1 << LOOP_QCDpt30);
}

int PredictionQCDBins ()
{
  return run<FakeRateLooper>(prediction_cuts_qcd_bins, "PredictionQCDBins", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

int ObservationQCDBinsEven ()
{
  return run<Looper>(observation_cuts_qcd_bins_even, "ObservationQCDBinsEven", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

int PredictionQCDBinsEven ()
{
  return run<FakeRateLooper>(prediction_cuts_qcd_bins_even, "PredictionQCDBinsEven", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

int ObservationQCDBinsOdd ()
{
  return run<Looper>(observation_cuts_qcd_bins_odd, "ObservationQCDBinsOdd", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

int PredictionQCDBinsOdd ()
{
  return run<FakeRateLooper>(prediction_cuts_qcd_bins_odd, "PredictionQCDBinsOdd", 1 << LOOP_QCDpt30to80 | 1 << LOOP_QCDpt80to170 | 1 << LOOP_QCDpt170to300 | 1 << LOOP_QCDpt300to470 | 1 << LOOP_QCDpt470to800 | 1 << LOOP_QCDpt800toInf);
}

