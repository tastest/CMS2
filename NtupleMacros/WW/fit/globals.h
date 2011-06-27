#ifndef GLOBALS_H
#define GLOBALS_H

class RooAbsData;
class RooAbsPdf;
class RooRealVar;

//
// variables
//

RooRealVar* var_deta;
RooRealVar* var_dilpt;
RooRealVar* var_dphi;
RooRealVar* var_mass;
RooRealVar* var_met;
RooRealVar* var_pt1;
RooRealVar* var_pt2;
// ROOT 5.26 workaround for bug when
// using parameterless likelihoods
RooRealVar* var_dummy;

void setVars();

//
// datasets
//

RooAbsData *ds_ww;
RooAbsData *ds_ww_deta;
RooAbsData *ds_ww_dilpt;
RooAbsData *ds_ww_dphi;
RooAbsData *ds_ww_mass;
RooAbsData *ds_ww_met;
RooAbsData *ds_ww_pt1;
RooAbsData *ds_ww_pt2;

RooAbsData *ds_ww_pt1vdeta;
RooAbsData *ds_ww_pt1vdilpt;
RooAbsData *ds_ww_pt1vdphi;
RooAbsData *ds_ww_pt1vmass;
RooAbsData *ds_ww_pt1vmet;
RooAbsData *ds_ww_pt1vpt2;

RooAbsData *ds_tt;
RooAbsData *ds_tt_deta;
RooAbsData *ds_tt_dilpt;
RooAbsData *ds_tt_dphi;
RooAbsData *ds_tt_mass;
RooAbsData *ds_tt_met;
RooAbsData *ds_tt_pt1;
RooAbsData *ds_tt_pt2;

RooAbsData *ds_tt_pt1vdeta;
RooAbsData *ds_tt_pt1vdilpt;
RooAbsData *ds_tt_pt1vdphi;
RooAbsData *ds_tt_pt1vmass;
RooAbsData *ds_tt_pt1vmet;
RooAbsData *ds_tt_pt1vpt2;

RooAbsData *ds_lz[20];
RooAbsData *ds_lz_deta[20];
RooAbsData *ds_lz_dilpt[20];
RooAbsData *ds_lz_dphi[20];
RooAbsData *ds_lz_mass[20];
RooAbsData *ds_lz_met[20];
RooAbsData *ds_lz_pt1[20];
RooAbsData *ds_lz_pt2[20];

RooAbsData *ds_lz_pt1vdeta[20];
RooAbsData *ds_lz_pt1vdilpt[20];
RooAbsData *ds_lz_pt1vdphi[20];
RooAbsData *ds_lz_pt1vmass[20];
RooAbsData *ds_lz_pt1vmet[20];
RooAbsData *ds_lz_pt1vpt2[20];

void setDataSets();
void clearDataSets();

//
// pdfs
//

RooAbsPdf *pdf_ww_deta;
RooAbsPdf *pdf_ww_dilpt;
RooAbsPdf *pdf_ww_dphi;
RooAbsPdf *pdf_ww_mass;
RooAbsPdf *pdf_ww_met;
RooAbsPdf *pdf_ww_pt1;
RooAbsPdf *pdf_ww_pt2;

RooAbsPdf *pdf_ww_pt1vdeta;
RooAbsPdf *pdf_ww_pt1vdilpt;
RooAbsPdf *pdf_ww_pt1vdphi;
RooAbsPdf *pdf_ww_pt1vmass;
RooAbsPdf *pdf_ww_pt1vmet;
RooAbsPdf *pdf_ww_pt1vpt2;

RooAbsPdf *pdf_tt_deta;
RooAbsPdf *pdf_tt_dilpt;
RooAbsPdf *pdf_tt_dphi;
RooAbsPdf *pdf_tt_mass;
RooAbsPdf *pdf_tt_met;
RooAbsPdf *pdf_tt_pt1;
RooAbsPdf *pdf_tt_pt2;

RooAbsPdf *pdf_tt_pt1vdeta;
RooAbsPdf *pdf_tt_pt1vdilpt;
RooAbsPdf *pdf_tt_pt1vdphi;
RooAbsPdf *pdf_tt_pt1vmass;
RooAbsPdf *pdf_tt_pt1vmet;
RooAbsPdf *pdf_tt_pt1vpt2;

RooAbsPdf *pdf_lz_deta[20];
RooAbsPdf *pdf_lz_dilpt[20];
RooAbsPdf *pdf_lz_dphi[20];
RooAbsPdf *pdf_lz_mass[20];
RooAbsPdf *pdf_lz_met[20];
RooAbsPdf *pdf_lz_pt1[20];
RooAbsPdf *pdf_lz_pt2[20];

RooAbsPdf *pdf_lz_pt1vdeta[20];
RooAbsPdf *pdf_lz_pt1vdilpt[20];
RooAbsPdf *pdf_lz_pt1vdphi[20];
RooAbsPdf *pdf_lz_pt1vmass[20];
RooAbsPdf *pdf_lz_pt1vmet[20];
RooAbsPdf *pdf_lz_pt1vpt2[20];

#endif
