#ifndef CRUIJFFPDFS_H
#define CRUIJFFPDFS_H

class RooRealVar;

RooRealVar* var_cru_m0_ww_pt1;
RooRealVar* var_cru_sl_ww_pt1;
RooRealVar* var_cru_sr_ww_pt1;
RooRealVar* var_cru_al_ww_pt1;
RooRealVar* var_cru_ar_ww_pt1;
RooRealVar* var_cru_m0_tt_pt1;
RooRealVar* var_cru_sl_tt_pt1;
RooRealVar* var_cru_sr_tt_pt1;
RooRealVar* var_cru_al_tt_pt1;
RooRealVar* var_cru_ar_tt_pt1;
RooRealVar* var_cru_m0_lz_pt1[20];
RooRealVar* var_cru_sl_lz_pt1[20];
RooRealVar* var_cru_sr_lz_pt1[20];
RooRealVar* var_cru_al_lz_pt1[20];
RooRealVar* var_cru_ar_lz_pt1[20];

RooRealVar* var_cru_m0_ww_pt2;
RooRealVar* var_cru_sl_ww_pt2;
RooRealVar* var_cru_sr_ww_pt2;
RooRealVar* var_cru_al_ww_pt2;
RooRealVar* var_cru_ar_ww_pt2;
RooRealVar* var_cru_m0_tt_pt2;
RooRealVar* var_cru_sl_tt_pt2;
RooRealVar* var_cru_sr_tt_pt2;
RooRealVar* var_cru_al_tt_pt2;
RooRealVar* var_cru_ar_tt_pt2;
RooRealVar* var_cru_m0_lz_pt2[20];
RooRealVar* var_cru_sl_lz_pt2[20];
RooRealVar* var_cru_sr_lz_pt2[20];
RooRealVar* var_cru_al_lz_pt2[20];
RooRealVar* var_cru_ar_lz_pt2[20];

void setPDFs();
void clearPDFs();

#endif
