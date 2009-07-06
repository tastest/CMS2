// -*- C++ -*-

#include <unistd.h>
#include <string>
#include "TFile.h"
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

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     Sample data = 
	  fWW()          +
	  fWZ()          +
	  fZZ()          +
	  fWjets()       +
	  fDYee()        +
	  fttbar()       +
	  fLM8();         
     data.name = "data";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
     Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
     Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
     Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_wjets	(fWjets()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wjets       .Loop();
     Looper looper_dyee		(fDYee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYEE  )) looper_dyee        .Loop();
     Looper looper_dymm		(fDYmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYMM  )) looper_dymm        .Loop();
     Looper looper_dytt		(fDYtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DYTT  )) looper_dytt        .Loop();
     Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW    )) looper_tw          .Loop();
     Looper looper_data		(data		, cuts, log.c_str());	looper_data          .Loop();
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
	  &looper_data        ,
     };
     // when all the loopers are done, we save the histograms to file
//      saveHist(hist.c_str());
     TFile *f = TFile::Open(hist.c_str(), "recreate");
//      DSGTables *tables = new DSGTables("tables");
     for (unsigned int i = 0; i < sizeof(loopers) / sizeof(Looper *); ++i) {
// 	  tables->tables_.push_back(loopers[i]->dsgTable);
	  loopers[i]->dsgTable.Write();
     }
//      tables->Write();
     f->Close();
     return 0;
}

// default yield table
int Results ()
{
     return run<Looper>(baseline_cuts, "Results");
}

void DSGDisplay ();
