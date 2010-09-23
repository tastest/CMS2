#ifndef GLOBALS_H
#define GLOBALS_H

#include "RooRealVar.h"
class RooAbsData;
class RooAbsPdf;

//
// variables
//

RooRealVar* var_deta  = new RooRealVar("deta","deta",0.,5.);
RooRealVar* var_dilpt = new RooRealVar("dilpt","dilpt",0.,500.);
RooRealVar* var_dphi  = new RooRealVar("dphi","dphi",0.,3.15);
RooRealVar* var_mass  = new RooRealVar("mass","mass",12.,800.);
RooRealVar* var_met   = new RooRealVar("met","met",20.,500.);
RooRealVar* var_pt1    = new RooRealVar("pt1","pt1",20.,500.);
RooRealVar* var_pt2    = new RooRealVar("pt2","pt2",20.,500.);

// ROOT 5.26 workaround for bug when
// using parameterless likelihoods
RooRealVar* var_dummy = new RooRealVar("var_dummy","var_dummy",0);

//
// datasets
//

RooAbsData *ds_ww_deta;
RooAbsData *ds_tt_deta;
RooAbsData *ds_lz_deta[20];

RooAbsData *ds_ww_dilpt;
RooAbsData *ds_tt_dilpt;
RooAbsData *ds_lz_dilpt[20];

RooAbsData *ds_ww_dphi;
RooAbsData *ds_tt_dphi;
RooAbsData *ds_lz_dphi[20];

RooAbsData *ds_ww_mass;
RooAbsData *ds_tt_mass;
RooAbsData *ds_lz_mass[20];

RooAbsData *ds_ww_met;
RooAbsData *ds_tt_met;
RooAbsData *ds_lz_met[20];

RooAbsData *ds_ww_pt1;
RooAbsData *ds_tt_pt1;
RooAbsData *ds_lz_pt1[20];
// njv := no jet veto
//RooAbsData *ds_ww_pt1_njv;
//RooAbsData *ds_tt_pt1_njv;
//RooAbsData *ds_lz_pt1_njv[20];

RooAbsData *ds_ww_pt2;
RooAbsData *ds_tt_pt2;
RooAbsData *ds_lz_pt2[20];
// njv := no jet veto
//RooAbsData *ds_ww_pt2_njv;
//RooAbsData *ds_tt_pt2_njv;
//RooAbsData *ds_lz_pt2_njv[20];

void setDataSets();
void clearDataSets();

//
// pdfs
//

RooAbsPdf *pdf_ww_deta;
RooAbsPdf *pdf_tt_deta;
RooAbsPdf *pdf_lz_deta[20];

RooAbsPdf *pdf_ww_dilpt;
RooAbsPdf *pdf_tt_dilpt;
RooAbsPdf *pdf_lz_dilpt[20];

RooAbsPdf *pdf_ww_dphi;
RooAbsPdf *pdf_tt_dphi;
RooAbsPdf *pdf_lz_dphi[20];

RooAbsPdf *pdf_ww_mass;
RooAbsPdf *pdf_tt_mass;
RooAbsPdf *pdf_lz_mass[20];

RooAbsPdf *pdf_ww_met;
RooAbsPdf *pdf_tt_met;
RooAbsPdf *pdf_lz_met[20];

RooAbsPdf *pdf_ww_pt1;
RooAbsPdf *pdf_tt_pt1;
RooAbsPdf *pdf_lz_pt1[20];
// njv := no jet veto
//RooAbsPdf *pdf_ww_pt1_njv;
//RooAbsPdf *pdf_tt_pt1_njv;
//RooAbsPdf *pdf_lz_pt1_njv[20];

RooAbsPdf *pdf_ww_pt2;
RooAbsPdf *pdf_tt_pt2;
RooAbsPdf *pdf_lz_pt2[20];
// njv := no jet veto
//RooAbsPdf *pdf_ww_pt2_njv;
//RooAbsPdf *pdf_tt_pt2_njv;
//RooAbsPdf *pdf_lz_pt2_njv[20];

#endif
