#ifndef HISTPDFS_H
#define HISTPDFS_H

class RooDataHist;

RooDataHist* dh_ww_pt;
RooDataHist* dh_tt_pt;
RooDataHist* dh_lz_pt[20];
RooDataHist* dh_ww_pt_njv;
RooDataHist* dh_tt_pt_njv;
RooDataHist* dh_lz_pt_njv[20];

void setPDFs();
void clearPDFs();

#endif
