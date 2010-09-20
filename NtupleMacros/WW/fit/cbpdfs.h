#ifndef CBPDFS_H
#define CBPDFS_H

class RooRealVar;

RooRealVar* var_cb_alpha_ww_pt;
RooRealVar* var_cb_m0_ww_pt;
RooRealVar* var_cb_n_ww_pt;
RooRealVar* var_cb_sigma_ww_pt;
RooRealVar* var_cb_alpha_tt_pt;
RooRealVar* var_cb_m0_tt_pt;
RooRealVar* var_cb_n_tt_pt;
RooRealVar* var_cb_sigma_tt_pt;
RooRealVar* var_cb_alpha_lz_pt[20];
RooRealVar* var_cb_m0_lz_pt[20];
RooRealVar* var_cb_n_lz_pt[20];
RooRealVar* var_cb_sigma_lz_pt[20];

void setPDFs();
void clearPDFs();

#endif
