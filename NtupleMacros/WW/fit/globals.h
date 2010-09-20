#ifndef GLOBALS_H
#define GLOBALS_H

#include "RooRealVar.h"
class RooAbsData;
class RooAbsPdf;

//
// variables
//

RooRealVar* var_pt = new RooRealVar("pt","pt",20.,500.);

// ROOT 5.26 workaround for bug when
// using parameterless likelihoods
RooRealVar* var_dummy = new RooRealVar("var_dummy","var_dummy",0);

//
// datasets
//

RooAbsData *ds_ww_pt;
RooAbsData *ds_tt_pt;
RooAbsData *ds_lz_pt[20];
// njv := no jet veto
//RooAbsData *ds_ww_pt_njv;
//RooAbsData *ds_tt_pt_njv;
//RooAbsData *ds_lz_pt_njv[20];

void setDataSets();
void clearDataSets();

//
// pdfs
//

RooAbsPdf *pdf_ww_pt;
RooAbsPdf *pdf_tt_pt;
RooAbsPdf *pdf_lz_pt[20];
// njv := no jet veto
//RooAbsPdf *pdf_ww_pt_njv;
//RooAbsPdf *pdf_tt_pt_njv;
//RooAbsPdf *pdf_lz_pt_njv[20];

#endif
