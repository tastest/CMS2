#ifndef CRUIJFFPDFS_H
#define CRUIJFFPDFS_H

class RooRealVar;

RooRealVar* var_cru_m0_ww_pt;
RooRealVar* var_cru_sl_ww_pt;
RooRealVar* var_cru_sr_ww_pt;
RooRealVar* var_cru_al_ww_pt;
RooRealVar* var_cru_ar_ww_pt;
RooRealVar* var_cru_m0_tt_pt;
RooRealVar* var_cru_sl_tt_pt;
RooRealVar* var_cru_sr_tt_pt;
RooRealVar* var_cru_al_tt_pt;
RooRealVar* var_cru_ar_tt_pt;
RooRealVar* var_cru_m0_lz_pt[20];
RooRealVar* var_cru_sl_lz_pt[20];
RooRealVar* var_cru_sr_lz_pt[20];
RooRealVar* var_cru_al_lz_pt[20];
RooRealVar* var_cru_ar_lz_pt[20];

void setPDFs();
void clearPDFs();

#endif
