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
     LOOP_VV	,
     LOOP_WJETS	,
     LOOP_DY	,
     LOOP_TTBAR	,
     LOOP_TW	,
};

template <class Looper> int run (cuts_t cuts, const string &name, uint32 which_ones = 0xffffffff)
{
     const string hist = name + ".root";
     const string tbl = name + ".tbl";
     const string log = name + ".log";
     Sample data = fdata();
// 	  fWW()          +
//  	  fWZ()          +
//  	  fZZ()          +
//  	  fWe()       +
//  	  fWmu()       +
//  	  fWtau()       +
//  	  fZee()        +
//  	  fZmm()        +
//  	  fZtt()        +
//  	  fttbar()       +
// 	  fLM8();         
//      data.name = "data";
     // by default, we run this list of samples; if we're told by the
     // which_ones bit field to skip a sample, we skip it
//      Looper looper_ww		(fWW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WW    )) looper_ww          .Loop();
//      Looper looper_wz		(fWZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WZ    )) looper_wz          .Loop();
//      Looper looper_zz		(fZZ()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_ZZ    )) looper_zz          .Loop();
     Looper looper_vv		(fVV()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_VV    )) looper_vv          .Loop();
     Looper looper_we		(fWe()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_we	        .Loop();
     Looper looper_wmu		(fWmu()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wmu	        .Loop();
     Looper looper_wtau		(fWtau()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_WJETS )) looper_wtau        .Loop();
     Looper looper_dyee		(fZee()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY  )) looper_dyee        .Loop();
     Looper looper_dymm		(fZmm()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY  )) looper_dymm        .Loop();
     Looper looper_dytt		(fZtt()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_DY  )) looper_dytt        .Loop();
     Looper looper_ttbar	(fttbar()	, cuts, log.c_str());	if (which_ones & (1 << LOOP_TTBAR )) looper_ttbar       .Loop();
     Looper looper_tw		(ftW()		, cuts, log.c_str());	if (which_ones & (1 << LOOP_TW )) looper_tw       .Loop();
     Looper looper_data		(data		, cuts, log.c_str());	looper_data          .Loop();
     // then we collect them all and print a table
     const Looper *loopers[] = { 
// 	  &looper_ww          ,
// 	  &looper_wz          ,
// 	  &looper_zz          ,
 	  &looper_vv          ,
	  &looper_we       ,
	  &looper_wmu       ,
	  &looper_wtau      ,
	  &looper_dyee        ,
	  &looper_dymm        ,
	  &looper_dytt        ,
	  &looper_ttbar       ,
	  &looper_tw       ,
	  &looper_data        ,
     };
     // when all the loopers are done, we save the histograms to file
//      saveHist(hist.c_str());
     TFile *f = TFile::Open(hist.c_str(), "recreate");
//      DSGTables *tables = new DSGTables("tables");
     for (unsigned int i = 0; i < sizeof(loopers) / sizeof(Looper *); ++i) {
// 	  tables->tables_.push_back(loopers[i]->dsgTable);
	  loopers[i]->dsgTable.Write();
	  loopers[i]->dsgTable_emu.Write();
     }
//      tables->Write();
     f->Close();
     return 0;
}

// default yield table, using everything but ttbar
int Results ()
{
     return run<Looper>(baseline_cuts, "Results", ~(1 << LOOP_TTBAR));
}

void DSGDisplay ();
