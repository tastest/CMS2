#include "globals.h"
#include "cruijffpdfs.h"

#include "../../Tools/CruijffPdf.h"
#include "RooRealVar.h"

void setPDFs() {
    var_cru_m0_ww_pt1 = new RooRealVar("var_cru_m0_ww_pt1","var_cru_m0_ww_pt1",0.,100.);
    var_cru_sl_ww_pt1 = new RooRealVar("var_cru_sl_ww_pt1","var_cru_sl_ww_pt1",0.,100.);
    var_cru_sr_ww_pt1 = new RooRealVar("var_cru_sr_ww_pt1","var_cru_sr_ww_pt1",0.,100.);
    var_cru_al_ww_pt1 = new RooRealVar("var_cru_al_ww_pt1","var_cru_al_ww_pt1",0.,100.);
    var_cru_ar_ww_pt1 = new RooRealVar("var_cru_ar_ww_pt1","var_cru_ar_ww_pt1",0.,100.);
    pdf_ww_pt1 = new CruijffPdf("pdf_ww_pt1","pdf_ww_pt1",*var_pt1,*var_cru_m0_ww_pt1,*var_cru_sl_ww_pt1,*var_cru_sr_ww_pt1,*var_cru_al_ww_pt1,*var_cru_ar_ww_pt1);
    ((CruijffPdf*)pdf_ww_pt1)->fitTo(*ds_ww_pt1);

    var_cru_m0_tt_pt1 = new RooRealVar("var_cru_m0_tt_pt1","var_cru_m0_tt_pt1",0.,100.);
    var_cru_sl_tt_pt1 = new RooRealVar("var_cru_sl_tt_pt1","var_cru_sl_tt_pt1",0.,100.);
    var_cru_sr_tt_pt1 = new RooRealVar("var_cru_sr_tt_pt1","var_cru_sr_tt_pt1",0.,100.);
    var_cru_al_tt_pt1 = new RooRealVar("var_cru_al_tt_pt1","var_cru_al_tt_pt1",0.,100.);
    var_cru_ar_tt_pt1 = new RooRealVar("var_cru_ar_tt_pt1","var_cru_ar_tt_pt1",0.,100.);
    pdf_tt_pt1 = new CruijffPdf("pdf_tt_pt1","pdf_tt_pt1",*var_pt1,*var_cru_m0_tt_pt1,*var_cru_sl_tt_pt1,*var_cru_sr_tt_pt1,*var_cru_al_tt_pt1,*var_cru_ar_tt_pt1);
    ((CruijffPdf*)pdf_tt_pt1)->fitTo(*ds_tt_pt1);

    var_cru_m0_ww_pt2 = new RooRealVar("var_cru_m0_ww_pt2","var_cru_m0_ww_pt2",0.,100.);
    var_cru_sl_ww_pt2 = new RooRealVar("var_cru_sl_ww_pt2","var_cru_sl_ww_pt2",0.,100.);
    var_cru_sr_ww_pt2 = new RooRealVar("var_cru_sr_ww_pt2","var_cru_sr_ww_pt2",0.,100.);
    var_cru_al_ww_pt2 = new RooRealVar("var_cru_al_ww_pt2","var_cru_al_ww_pt2",0.,100.);
    var_cru_ar_ww_pt2 = new RooRealVar("var_cru_ar_ww_pt2","var_cru_ar_ww_pt2",0.,100.);
    pdf_ww_pt2 = new CruijffPdf("pdf_ww_pt2","pdf_ww_pt2",*var_pt2,*var_cru_m0_ww_pt2,*var_cru_sl_ww_pt2,*var_cru_sr_ww_pt2,*var_cru_al_ww_pt2,*var_cru_ar_ww_pt2);
    ((CruijffPdf*)pdf_ww_pt2)->fitTo(*ds_ww_pt2);

    var_cru_m0_tt_pt2 = new RooRealVar("var_cru_m0_tt_pt2","var_cru_m0_tt_pt2",0.,100.);
    var_cru_sl_tt_pt2 = new RooRealVar("var_cru_sl_tt_pt2","var_cru_sl_tt_pt2",0.,100.);
    var_cru_sr_tt_pt2 = new RooRealVar("var_cru_sr_tt_pt2","var_cru_sr_tt_pt2",0.,100.);
    var_cru_al_tt_pt2 = new RooRealVar("var_cru_al_tt_pt2","var_cru_al_tt_pt2",0.,100.);
    var_cru_ar_tt_pt2 = new RooRealVar("var_cru_ar_tt_pt2","var_cru_ar_tt_pt2",0.,100.);
    pdf_tt_pt2 = new CruijffPdf("pdf_tt_pt2","pdf_tt_pt2",*var_pt2,*var_cru_m0_tt_pt2,*var_cru_sl_tt_pt2,*var_cru_sr_tt_pt2,*var_cru_al_tt_pt2,*var_cru_ar_tt_pt2);
    ((CruijffPdf*)pdf_tt_pt2)->fitTo(*ds_tt_pt2);

    for(unsigned int i = 0; i < 20; ++i) {
        var_cru_m0_lz_pt1[i] = new RooRealVar(Form("var_cru_m0_lz%i_pt1",i),Form("var_cru_m0_lz%i_pt1",i),0.,100.);
        var_cru_sl_lz_pt1[i] = new RooRealVar(Form("var_cru_sl_lz%i_pt1",i),Form("var_cru_sl_lz%i_pt1",i),0.,100.);
        var_cru_sr_lz_pt1[i] = new RooRealVar(Form("var_cru_sr_lz%i_pt1",i),Form("var_cru_sr_lz%i_pt1",i),0.,100.);
        var_cru_al_lz_pt1[i] = new RooRealVar(Form("var_cru_al_lz%i_pt1",i),Form("var_cru_al_lz%i_pt1",i),0.,100.);
        var_cru_ar_lz_pt1[i] = new RooRealVar(Form("var_cru_ar_lz%i_pt1",i),Form("var_cru_ar_lz%i_pt1",i),0.,100.);
        pdf_lz_pt1[i] = new CruijffPdf(Form("pdf_lz%i_pt1",i),Form("pdf_lz%i_pt1",i),*var_pt1,*var_cru_m0_lz_pt1[i],*var_cru_sl_lz_pt1[i],*var_cru_sr_lz_pt1[i],*var_cru_al_lz_pt1[i],*var_cru_ar_lz_pt1[i]);
        ((CruijffPdf*)pdf_lz_pt1[i])->fitTo(*ds_lz_pt1[i]);

        var_cru_m0_lz_pt2[i] = new RooRealVar(Form("var_cru_m0_lz%i_pt2",i),Form("var_cru_m0_lz%i_pt2",i),0.,100.);
        var_cru_sl_lz_pt2[i] = new RooRealVar(Form("var_cru_sl_lz%i_pt2",i),Form("var_cru_sl_lz%i_pt2",i),0.,100.);
        var_cru_sr_lz_pt2[i] = new RooRealVar(Form("var_cru_sr_lz%i_pt2",i),Form("var_cru_sr_lz%i_pt2",i),0.,100.);
        var_cru_al_lz_pt2[i] = new RooRealVar(Form("var_cru_al_lz%i_pt2",i),Form("var_cru_al_lz%i_pt2",i),0.,100.);
        var_cru_ar_lz_pt2[i] = new RooRealVar(Form("var_cru_ar_lz%i_pt2",i),Form("var_cru_ar_lz%i_pt2",i),0.,100.);
        pdf_lz_pt2[i] = new CruijffPdf(Form("pdf_lz%i_pt2",i),Form("pdf_lz%i_pt2",i),*var_pt2,*var_cru_m0_lz_pt2[i],*var_cru_sl_lz_pt2[i],*var_cru_sr_lz_pt2[i],*var_cru_al_lz_pt2[i],*var_cru_ar_lz_pt2[i]);
        ((CruijffPdf*)pdf_lz_pt2[i])->fitTo(*ds_lz_pt2[i]);
    }
}

void clearPDFs() {
    delete var_cru_m0_ww_pt1;
    delete var_cru_sl_ww_pt1;
    delete var_cru_sr_ww_pt1;
    delete var_cru_al_ww_pt1;
    delete var_cru_ar_ww_pt1;
    delete var_cru_m0_tt_pt1;
    delete var_cru_sl_tt_pt1;
    delete var_cru_sr_tt_pt1;
    delete var_cru_al_tt_pt1;
    delete var_cru_ar_tt_pt1;
    delete (CruijffPdf*)pdf_ww_pt1;
    delete (CruijffPdf*)pdf_tt_pt1;
    delete var_cru_m0_ww_pt2;
    delete var_cru_sl_ww_pt2;
    delete var_cru_sr_ww_pt2;
    delete var_cru_al_ww_pt2;
    delete var_cru_ar_ww_pt2;
    delete var_cru_m0_tt_pt2;
    delete var_cru_sl_tt_pt2;
    delete var_cru_sr_tt_pt2;
    delete var_cru_al_tt_pt2;
    delete var_cru_ar_tt_pt2;
    delete (CruijffPdf*)pdf_ww_pt2;
    delete (CruijffPdf*)pdf_tt_pt2;

    for(unsigned int i = 0; i < 20; ++i) {
        delete var_cru_m0_lz_pt1[i];
        delete var_cru_sl_lz_pt1[i];
        delete var_cru_sr_lz_pt1[i];
        delete var_cru_al_lz_pt1[i];
        delete var_cru_ar_lz_pt1[i];
        delete (CruijffPdf*)pdf_lz_pt1[i];
        delete var_cru_m0_lz_pt2[i];
        delete var_cru_sl_lz_pt2[i];
        delete var_cru_sr_lz_pt2[i];
        delete var_cru_al_lz_pt2[i];
        delete var_cru_ar_lz_pt2[i];
        delete (CruijffPdf*)pdf_lz_pt2[i];
    }
}
