#include "globals.h"
#include "cruijffpdfs.h"

#include "../../Tools/CruijffPdf.h"
#include "RooRealVar.h"

void setPDFs() {
    var_cru_m0_ww_pt = new RooRealVar("var_cru_m0_ww_pt","var_cru_m0_ww_pt",0.,100.);
    var_cru_sl_ww_pt = new RooRealVar("var_cru_sl_ww_pt","var_cru_sl_ww_pt",0.,100.);
    var_cru_sr_ww_pt = new RooRealVar("var_cru_sr_ww_pt","var_cru_sr_ww_pt",0.,100.);
    var_cru_al_ww_pt = new RooRealVar("var_cru_al_ww_pt","var_cru_al_ww_pt",0.,100.);
    var_cru_ar_ww_pt = new RooRealVar("var_cru_ar_ww_pt","var_cru_ar_ww_pt",0.,100.);
    pdf_ww_pt = new CruijffPdf("pdf_ww_pt","pdf_ww_pt",*var_pt,*var_cru_m0_ww_pt,*var_cru_sl_ww_pt,*var_cru_sr_ww_pt,*var_cru_al_ww_pt,*var_cru_ar_ww_pt);
    ((CruijffPdf*)pdf_ww_pt)->fitTo(*ds_ww_pt);

    var_cru_m0_tt_pt = new RooRealVar("var_cru_m0_tt_pt","var_cru_m0_tt_pt",0.,100.);
    var_cru_sl_tt_pt = new RooRealVar("var_cru_sl_tt_pt","var_cru_sl_tt_pt",0.,100.);
    var_cru_sr_tt_pt = new RooRealVar("var_cru_sr_tt_pt","var_cru_sr_tt_pt",0.,100.);
    var_cru_al_tt_pt = new RooRealVar("var_cru_al_tt_pt","var_cru_al_tt_pt",0.,100.);
    var_cru_ar_tt_pt = new RooRealVar("var_cru_ar_tt_pt","var_cru_ar_tt_pt",0.,100.);
    pdf_tt_pt = new CruijffPdf("pdf_tt_pt","pdf_tt_pt",*var_pt,*var_cru_m0_tt_pt,*var_cru_sl_tt_pt,*var_cru_sr_tt_pt,*var_cru_al_tt_pt,*var_cru_ar_tt_pt);
    ((CruijffPdf*)pdf_tt_pt)->fitTo(*ds_tt_pt);

    for(unsigned int i = 0; i < 20; ++i) {
        var_cru_m0_lz_pt[i] = new RooRealVar(Form("var_cru_m0_lz%i_pt",i),Form("var_cru_m0_lz%i_pt",i),0.,100.);
        var_cru_sl_lz_pt[i] = new RooRealVar(Form("var_cru_sl_lz%i_pt",i),Form("var_cru_sl_lz%i_pt",i),0.,100.);
        var_cru_sr_lz_pt[i] = new RooRealVar(Form("var_cru_sr_lz%i_pt",i),Form("var_cru_sr_lz%i_pt",i),0.,100.);
        var_cru_al_lz_pt[i] = new RooRealVar(Form("var_cru_al_lz%i_pt",i),Form("var_cru_al_lz%i_pt",i),0.,100.);
        var_cru_ar_lz_pt[i] = new RooRealVar(Form("var_cru_ar_lz%i_pt",i),Form("var_cru_ar_lz%i_pt",i),0.,100.);
        pdf_lz_pt[i] = new CruijffPdf(Form("pdf_lz%i_pt",i),Form("pdf_lz%i_pt",i),*var_pt,*var_cru_m0_lz_pt[i],*var_cru_sl_lz_pt[i],*var_cru_sr_lz_pt[i],*var_cru_al_lz_pt[i],*var_cru_ar_lz_pt[i]);
        ((CruijffPdf*)pdf_lz_pt[i])->fitTo(*ds_lz_pt[i]);
    }
}

void clearPDFs() {
    delete var_cru_m0_ww_pt;
    delete var_cru_sl_ww_pt;
    delete var_cru_sr_ww_pt;
    delete var_cru_al_ww_pt;
    delete var_cru_ar_ww_pt;
    delete var_cru_m0_tt_pt;
    delete var_cru_sl_tt_pt;
    delete var_cru_sr_tt_pt;
    delete var_cru_al_tt_pt;
    delete var_cru_ar_tt_pt;
    delete (CruijffPdf*)pdf_ww_pt;
    delete (CruijffPdf*)pdf_tt_pt;

    for(unsigned int i = 0; i < 20; ++i) {
        delete var_cru_m0_lz_pt[i];
        delete var_cru_sl_lz_pt[i];
        delete var_cru_sr_lz_pt[i];
        delete var_cru_al_lz_pt[i];
        delete var_cru_ar_lz_pt[i];
        delete (CruijffPdf*)pdf_lz_pt[i];
    }
}
