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
  LOOP_InclusiveMu5Pt50, 
  LOOP_InclusiveMuPt15, 
  LOOP_QCDBCtoEPt20to30, 
  LOOP_QCDBCtoEPt30to80, 
  LOOP_QCDBCtoEPt80to170, 
  LOOP_QCDEMenrichedPt20to30, 
  LOOP_QCDEMenrichedPt30to80, 
  LOOP_QCDEMenrichedPt80to170, 
  LOOP_QCDpt30, 
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
template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
  const string hist = name + ".root";
  const string tbl = name + ".tbl";
  const string log = name + ".log";
  // by default, we run this list of samples; if we're told by the
  // which_ones bit field to skip a sample, we skip it
//   Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
//   Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
//   Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
//   Looper looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
//   Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
//   Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
//   Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
//   Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
//   Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
  Looper looper_InclusiveMu5Pt50                (fInclusiveMu5Pt50()            , cuts, log.c_str());   if (which_ones & (1 << LOOP_InclusiveMu5Pt50 )) looper_InclusiveMu5Pt50       .Loop(); 
  Looper looper_InclusiveMuPt15                 (fInclusiveMuPt15()             , cuts, log.c_str());   if (which_ones & (1 << LOOP_InclusiveMuPt15  )) looper_InclusiveMuPt15        .Loop(); 
  Looper looper_QCDBCtoEPt20to30                (fQCDBCtoEPt20to30()            , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDBCtoEPt20to30 )) looper_QCDBCtoEPt20to30       .Loop(); 
  Looper looper_QCDBCtoEPt30to80                (fQCDBCtoEPt30to80()            , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDBCtoEPt30to80 )) looper_QCDBCtoEPt30to80       .Loop(); 
  Looper looper_QCDBCtoEPt80to170               (fQCDBCtoEPt80to170()           , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDBCtoEPt80to170)) looper_QCDBCtoEPt80to170      .Loop(); 
  Looper looper_QCDEMenrichedPt20to30           (fQCDEMenrichedPt20to30()       , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDEMenrichedPt20to30     )) looper_QCDEMenrichedPt20to30  .Loop(); 
  Looper looper_QCDEMenrichedPt30to80           (fQCDEMenrichedPt30to80()       , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDEMenrichedPt30to80     )) looper_QCDEMenrichedPt30to80  .Loop(); 
  Looper looper_QCDEMenrichedPt80to170          (fQCDEMenrichedPt80to170()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDEMenrichedPt80to170    )) looper_QCDEMenrichedPt80to170 .Loop(); 
  Looper looper_QCDpt30          (fQCDpt30()      , cuts, log.c_str());   if (which_ones & (1 << LOOP_QCDpt30    )) looper_QCDpt30 .Loop(); 
  // when all the loopers are done, we save the histograms to file
  saveHist(hist.c_str());
  // then we collect them all and print a table
  const Looper *loopers[] = { 
//     &looper_ww          ,
//     &looper_wz          ,
//     &looper_zz          ,
//     &looper_wjets       ,
//     &looper_dyee        ,
//     &looper_dymm        ,
//     &looper_dytt        ,
//     &looper_ttbar       ,
//     &looper_tw          ,
    &looper_InclusiveMu5Pt50, 
    &looper_InclusiveMuPt15, 
    &looper_QCDBCtoEPt20to30, 
    &looper_QCDBCtoEPt30to80, 
    &looper_QCDBCtoEPt80to170, 
    &looper_QCDEMenrichedPt20to30, 
    &looper_QCDEMenrichedPt30to80, 
    &looper_QCDEMenrichedPt80to170, 
    &looper_QCDpt30, 
  };
  printTable(loopers, sizeof(loopers) / sizeof(Looper *), tbl.c_str(), which_ones);
  return 0;
}

// default yield table
int Results ()
{
  return run<Looper>(baseline_cuts, "Results");
}
