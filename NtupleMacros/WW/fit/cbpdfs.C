#include "globals.h"
#include "cbpdfs.h"

#include "RooCBShape.h"
#include "RooRealVar.h"

void setPDFs() {
    var_cb_alpha_ww_pt = new RooRealVar("var_cb_alpha_ww_pt","var_cb_alpha_ww_pt",-10.,0.);
    var_cb_m0_ww_pt    = new RooRealVar("var_cb_m0_ww_pt","var_cb_m0_ww_pt",0.,100.);
    var_cb_n_ww_pt     = new RooRealVar("var_cb_n_ww_pt","var_cb_n_ww_pt",0.,10.);
    var_cb_sigma_ww_pt = new RooRealVar("var_cb_sigma_ww_pt","var_cb_sigma_ww_pt",0.,100.);
    pdf_ww_pt = new RooCBShape("pdf_ww_pt","pdf_ww_pt",*var_pt,*var_cb_m0_ww_pt,*var_cb_sigma_ww_pt,*var_cb_alpha_ww_pt,*var_cb_n_ww_pt);
    ((RooCBShape*)pdf_ww_pt)->fitTo(*ds_ww_pt);

    var_cb_alpha_tt_pt = new RooRealVar("var_cb_alpha_tt_pt","var_cb_alpha_tt_pt",-10.,0.);
    var_cb_m0_tt_pt    = new RooRealVar("var_cb_m0_tt_pt","var_cb_m0_tt_pt",0.,100.);
    var_cb_n_tt_pt     = new RooRealVar("var_cb_n_tt_pt","var_cb_n_tt_pt",0.,10.);
    var_cb_sigma_tt_pt = new RooRealVar("var_cb_sigma_tt_pt","var_cb_sigma_tt_pt",0.,100.);
    pdf_tt_pt = new RooCBShape("pdf_tt_pt","pdf_tt_pt",*var_pt,*var_cb_m0_tt_pt,*var_cb_sigma_tt_pt,*var_cb_alpha_tt_pt,*var_cb_n_tt_pt);
    ((RooCBShape*)pdf_tt_pt)->fitTo(*ds_tt_pt);

    for(unsigned int i = 0; i < 20; ++i) {
        var_cb_alpha_lz_pt[i] = new RooRealVar(Form("var_cb_alpha_lz%i_pt",i),Form("var_cb_alpha_lz%i_pt",i),-10.,0.);
        var_cb_m0_lz_pt[i]    = new RooRealVar(Form("var_cb_m0_lz%i_pt",i),Form("var_cb_m0_lz%i_pt",i),0.,100.);
        var_cb_n_lz_pt[i]     = new RooRealVar(Form("var_cb_n_lz%i_pt",i),Form("var_cb_n_lz%i_pt",i),0.,10.);
        var_cb_sigma_lz_pt[i] = new RooRealVar(Form("var_cb_sigma_lz%i_pt",i),Form("var_cb_sigma_lz%i_pt",i),0.,100.);
        pdf_lz_pt[i] = new RooCBShape(Form("pdf_lz%i_pt",i),Form("pdf_lz%i_pt",i),*var_pt,*var_cb_m0_lz_pt[i],*var_cb_sigma_lz_pt[i],*var_cb_alpha_lz_pt[i],*var_cb_n_lz_pt[i]);
        ((RooCBShape*)pdf_lz_pt[i])->fitTo(*ds_lz_pt[i]);
    }
}

void clearPDFs() {
    delete var_cb_alpha_ww_pt;
    delete var_cb_m0_ww_pt;
    delete var_cb_n_ww_pt;
    delete var_cb_sigma_ww_pt;
    delete var_cb_alpha_tt_pt;
    delete var_cb_m0_tt_pt;
    delete var_cb_n_tt_pt;
    delete var_cb_sigma_tt_pt;
    delete (RooCBShape*)pdf_ww_pt;
    delete (RooCBShape*)pdf_tt_pt;

    for(unsigned int i = 0; i < 20; ++i) {
        delete var_cb_alpha_lz_pt[i];
        delete var_cb_m0_lz_pt[i];
        delete var_cb_n_lz_pt[i];
        delete var_cb_sigma_lz_pt[i];
        delete (RooCBShape*)pdf_lz_pt[i];
    }
}
