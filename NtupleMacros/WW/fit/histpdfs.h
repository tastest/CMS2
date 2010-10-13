#ifndef HISTPDFS_H
#define HISTPDFS_H

class RooDataHist;

RooDataHist* dh_ww_deta;
RooDataHist* dh_tt_deta;
RooDataHist* dh_lz_deta[20];

RooDataHist* dh_ww_dilpt;
RooDataHist* dh_tt_dilpt;
RooDataHist* dh_lz_dilpt[20];

RooDataHist* dh_ww_dphi;
RooDataHist* dh_tt_dphi;
RooDataHist* dh_lz_dphi[20];

RooDataHist* dh_ww_mass;
RooDataHist* dh_tt_mass;
RooDataHist* dh_lz_mass[20];

RooDataHist* dh_ww_met;
RooDataHist* dh_tt_met;
RooDataHist* dh_lz_met[20];

RooDataHist* dh_ww_pt1;
RooDataHist* dh_tt_pt1;
RooDataHist* dh_lz_pt1[20];
//RooDataHist* dh_ww_pt1_njv;
//RooDataHist* dh_tt_pt1_njv;
//RooDataHist* dh_lz_pt1_njv[20];

RooDataHist* dh_ww_pt2;
RooDataHist* dh_tt_pt2;
RooDataHist* dh_lz_pt2[20];
//RooDataHist* dh_ww_pt2_njv;
//RooDataHist* dh_tt_pt2_njv;
//RooDataHist* dh_lz_pt2_njv[20];

void setPDFs();
void clearPDFs();

#endif
