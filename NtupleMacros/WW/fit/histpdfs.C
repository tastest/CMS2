#include "globals.h"
#include "histpdfs.h"

#include "RooDataHist.h"
#include "RooHistPdf.h"

void setPDFs() {
    dh_ww_pt = new RooDataHist("dh_ww_pt","dh_ww_pt",*var_pt,*ds_ww_pt);
    dh_tt_pt = new RooDataHist("dh_tt_pt","dh_tt_pt",*var_pt,*ds_tt_pt);
    //dh_ww_pt_njv = new RooDataHist("dh_ww_pt_njv","dh_ww_pt_njv",*var_pt,*ds_ww_pt_njv);
    //dh_tt_pt_njv = new RooDataHist("dh_tt_pt_njv","dh_tt_pt_njv",*var_pt,*ds_tt_pt_njv);

    pdf_ww_pt = new RooHistPdf("pdf_ww_pt","pdf_ww_pt",*var_pt,*dh_ww_pt);
    pdf_tt_pt = new RooHistPdf("pdf_tt_pt","pdf_tt_pt",*var_pt,*dh_tt_pt);
    //pdf_ww_pt_njv = new RooHistPdf("pdf_ww_pt_njv","pdf_ww_pt_njv",*var_pt,*dh_ww_pt_njv);
    //pdf_tt_pt_njv = new RooHistPdf("pdf_tt_pt_njv","pdf_tt_pt_njv",*var_pt,*dh_tt_pt_njv);

    for(unsigned int i = 0; i < 20; ++i) {
        dh_lz_pt[i] = new RooDataHist(Form("dh_lz%i_pt",i),Form("dh_lz%i_pt",i),*var_pt,*ds_lz_pt[i]);
        //dh_lz_pt_njv[i] = new RooDataHist(Form("dh_lz%i_pt_njv",i),Form("dh_lz%i_pt_njv",i),*var_pt,*ds_lz_pt_njv[i]);
        pdf_lz_pt[i] = new RooHistPdf(Form("pdf_lz%i_pt",i),Form("pdf_lz%i_pt",i),*var_pt,*dh_lz_pt[i]);
        //pdf_lz_pt_njv[i] = new RooHistPdf(Form("pdf_lz%i_pt_njv",i),Form("pdf_lz%i_pt_njv",i),*var_pt,*dh_lz_pt_njv[i]);
    }
}

void clearPDFs() {
    delete dh_ww_pt;
    delete dh_tt_pt;
    //delete dh_ww_pt_njv;
    //delete dh_tt_pt_njv;
    delete (RooHistPdf*)pdf_ww_pt;
    delete (RooHistPdf*)pdf_tt_pt;
    //delete (RooHistPdf*)pdf_ww_pt_njv;
    //delete (RooHistPdf*)pdf_tt_pt_njv;

    for(unsigned int i = 0; i < 20; ++i) {
        delete dh_lz_pt[i];
        //delete dh_lz_pt_njv[i];
        delete (RooHistPdf*)pdf_lz_pt[i];
        //delete (RooHistPdf*)pdf_lz_pt_njv[i];
    }
}
