#ifndef HISTPDFS_H
#define HISTPDFS_H

class RooDataHist;

RooDataHist* dh_ww_deta;
RooDataHist* dh_ww_dilpt;
RooDataHist* dh_ww_dphi;
RooDataHist* dh_ww_mass;
RooDataHist* dh_ww_met;
RooDataHist* dh_ww_pt1;
RooDataHist* dh_ww_pt2;

RooDataHist* dh_ww_pt1vdeta;
RooDataHist* dh_ww_pt1vdilpt;
RooDataHist* dh_ww_pt1vdphi;
RooDataHist* dh_ww_pt1vmass;
RooDataHist* dh_ww_pt1vmet;
RooDataHist* dh_ww_pt1vpt2;

RooDataHist* dh_tt_deta;
RooDataHist* dh_tt_dilpt;
RooDataHist* dh_tt_dphi;
RooDataHist* dh_tt_mass;
RooDataHist* dh_tt_met;
RooDataHist* dh_tt_pt1;
RooDataHist* dh_tt_pt2;

RooDataHist* dh_tt_pt1vdeta;
RooDataHist* dh_tt_pt1vdilpt;
RooDataHist* dh_tt_pt1vdphi;
RooDataHist* dh_tt_pt1vmass;
RooDataHist* dh_tt_pt1vmet;
RooDataHist* dh_tt_pt1vpt2;

RooDataHist* dh_lz_deta[20];
RooDataHist* dh_lz_dilpt[20];
RooDataHist* dh_lz_dphi[20];
RooDataHist* dh_lz_mass[20];
RooDataHist* dh_lz_met[20];
RooDataHist* dh_lz_pt1[20];
RooDataHist* dh_lz_pt2[20];

RooDataHist* dh_lz_pt1vdeta[20];
RooDataHist* dh_lz_pt1vdilpt[20];
RooDataHist* dh_lz_pt1vdphi[20];
RooDataHist* dh_lz_pt1vmass[20];
RooDataHist* dh_lz_pt1vmet[20];
RooDataHist* dh_lz_pt1vpt2[20];

void setPDFs();
void clearPDFs();

#endif
