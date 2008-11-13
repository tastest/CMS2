#include <unistd.h>
#include <string>
#include "DSGLooper.h"
#include "Sample.h"
#include "utilities.h"
#include "printDSGTable.h"

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
     LOOP_LM1   ,
     LOOP_LM2   ,
     LOOP_LM4   ,
};

template <class Looper> int runDSG (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     using std::string;
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     Looper looper_lm1               (fLM1()         , cuts, log.c_str());   if (which_ones & (1 << LOOP_LM1   )) looper_lm1         .Loop();
     Looper looper_lm2               (fLM2()         , cuts, log.c_str());   if (which_ones & (1 << LOOP_LM2   )) looper_lm2         .Loop();
     Looper looper_lm4               (fLM4()         , cuts, log.c_str());   if (which_ones & (1 << LOOP_LM4   )) looper_lm4         .Loop();
     saveHist(hist.c_str());
     const DSGLooperBase *loopers[] = { 
	  &looper_ww          ,
	  &looper_wz          ,
	  &looper_zz          ,
	  &looper_wjets       ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw          ,
          &looper_lm1         ,
          &looper_lm2         ,
          &looper_lm4         ,
     };
     printDSGTable(loopers, sizeof(loopers) / sizeof(DSGLooperBase *), tbl.c_str(), which_ones);
     return 0;
}

int DSG_BaseLine ()
{
     return runDSG<DSGResultsLooper>(dsg_baseline_cuts, "DSG_BaseLine");
}

int DSG_MET_10 ()
{
  return runDSG<DSGResultsLooper>(dsg_met_10_cuts, "DSG_MET_10");
}

int DSG_MET_1 ()
{
     return runDSG<DSGResultsLooper>(dsg_met_1_cuts, "DSG_MET_1");
}
