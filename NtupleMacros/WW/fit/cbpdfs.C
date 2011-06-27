#include "globals.h"
#include "cbpdfs.h"

#include "RooCBShape.h"
#include "RooRealVar.h"

void setPDFs() {
    var_cb_alpha_ww_pt1 = new RooRealVar("var_cb_alpha_ww_pt1","var_cb_alpha_ww_pt1",-10.,0.);
    var_cb_m0_ww_pt1    = new RooRealVar("var_cb_m0_ww_pt1","var_cb_m0_ww_pt1",0.,100.);
    var_cb_n_ww_pt1     = new RooRealVar("var_cb_n_ww_pt1","var_cb_n_ww_pt1",0.,10.);
    var_cb_sigma_ww_pt1 = new RooRealVar("var_cb_sigma_ww_pt1","var_cb_sigma_ww_pt1",0.,100.);
    pdf_ww_pt1 = new RooCBShape("pdf_ww_pt1","pdf_ww_pt1",*var_pt1,*var_cb_m0_ww_pt1,*var_cb_sigma_ww_pt1,*var_cb_alpha_ww_pt1,*var_cb_n_ww_pt1);
    ((RooCBShape*)pdf_ww_pt1)->fitTo(*ds_ww_pt1);

    var_cb_alpha_tt_pt1 = new RooRealVar("var_cb_alpha_tt_pt1","var_cb_alpha_tt_pt1",-10.,0.);
    var_cb_m0_tt_pt1    = new RooRealVar("var_cb_m0_tt_pt1","var_cb_m0_tt_pt1",0.,100.);
    var_cb_n_tt_pt1     = new RooRealVar("var_cb_n_tt_pt1","var_cb_n_tt_pt1",0.,10.);
    var_cb_sigma_tt_pt1 = new RooRealVar("var_cb_sigma_tt_pt1","var_cb_sigma_tt_pt1",0.,100.);
    pdf_tt_pt1 = new RooCBShape("pdf_tt_pt1","pdf_tt_pt1",*var_pt1,*var_cb_m0_tt_pt1,*var_cb_sigma_tt_pt1,*var_cb_alpha_tt_pt1,*var_cb_n_tt_pt1);
    ((RooCBShape*)pdf_tt_pt1)->fitTo(*ds_tt_pt1);

    var_cb_alpha_ww_pt2 = new RooRealVar("var_cb_alpha_ww_pt2","var_cb_alpha_ww_pt2",-10.,0.);
    var_cb_m0_ww_pt2    = new RooRealVar("var_cb_m0_ww_pt2","var_cb_m0_ww_pt2",0.,100.);
    var_cb_n_ww_pt2     = new RooRealVar("var_cb_n_ww_pt2","var_cb_n_ww_pt2",0.,10.);
    var_cb_sigma_ww_pt2 = new RooRealVar("var_cb_sigma_ww_pt2","var_cb_sigma_ww_pt2",0.,100.);
    pdf_ww_pt2 = new RooCBShape("pdf_ww_pt2","pdf_ww_pt2",*var_pt2,*var_cb_m0_ww_pt2,*var_cb_sigma_ww_pt2,*var_cb_alpha_ww_pt2,*var_cb_n_ww_pt2);
    ((RooCBShape*)pdf_ww_pt2)->fitTo(*ds_ww_pt2);

    var_cb_alpha_tt_pt2 = new RooRealVar("var_cb_alpha_tt_pt2","var_cb_alpha_tt_pt2",-10.,0.);
    var_cb_m0_tt_pt2    = new RooRealVar("var_cb_m0_tt_pt2","var_cb_m0_tt_pt2",0.,100.);
    var_cb_n_tt_pt2     = new RooRealVar("var_cb_n_tt_pt2","var_cb_n_tt_pt2",0.,10.);
    var_cb_sigma_tt_pt2 = new RooRealVar("var_cb_sigma_tt_pt2","var_cb_sigma_tt_pt2",0.,100.);
    pdf_tt_pt2 = new RooCBShape("pdf_tt_pt2","pdf_tt_pt2",*var_pt2,*var_cb_m0_tt_pt2,*var_cb_sigma_tt_pt2,*var_cb_alpha_tt_pt2,*var_cb_n_tt_pt2);
    ((RooCBShape*)pdf_tt_pt2)->fitTo(*ds_tt_pt2);

    for(unsigned int i = 0; i < 20; ++i) {
        var_cb_alpha_lz_pt1[i] = new RooRealVar(Form("var_cb_alpha_lz%i_pt1",i),Form("var_cb_alpha_lz%i_pt1",i),-10.,0.);
        var_cb_m0_lz_pt1[i]    = new RooRealVar(Form("var_cb_m0_lz%i_pt1",i),Form("var_cb_m0_lz%i_pt1",i),0.,100.);
        var_cb_n_lz_pt1[i]     = new RooRealVar(Form("var_cb_n_lz%i_pt1",i),Form("var_cb_n_lz%i_pt1",i),0.,10.);
        var_cb_sigma_lz_pt1[i] = new RooRealVar(Form("var_cb_sigma_lz%i_pt1",i),Form("var_cb_sigma_lz%i_pt1",i),0.,100.);
        pdf_lz_pt1[i] = new RooCBShape(Form("pdf_lz%i_pt1",i),Form("pdf_lz%i_pt1",i),*var_pt1,*var_cb_m0_lz_pt1[i],*var_cb_sigma_lz_pt1[i],*var_cb_alpha_lz_pt1[i],*var_cb_n_lz_pt1[i]);
        ((RooCBShape*)pdf_lz_pt1[i])->fitTo(*ds_lz_pt1[i]);

        var_cb_alpha_lz_pt2[i] = new RooRealVar(Form("var_cb_alpha_lz%i_pt2",i),Form("var_cb_alpha_lz%i_pt2",i),-10.,0.);
        var_cb_m0_lz_pt2[i]    = new RooRealVar(Form("var_cb_m0_lz%i_pt2",i),Form("var_cb_m0_lz%i_pt2",i),0.,100.);
        var_cb_n_lz_pt2[i]     = new RooRealVar(Form("var_cb_n_lz%i_pt2",i),Form("var_cb_n_lz%i_pt2",i),0.,10.);
        var_cb_sigma_lz_pt2[i] = new RooRealVar(Form("var_cb_sigma_lz%i_pt2",i),Form("var_cb_sigma_lz%i_pt2",i),0.,100.);
        pdf_lz_pt2[i] = new RooCBShape(Form("pdf_lz%i_pt2",i),Form("pdf_lz%i_pt2",i),*var_pt2,*var_cb_m0_lz_pt2[i],*var_cb_sigma_lz_pt2[i],*var_cb_alpha_lz_pt2[i],*var_cb_n_lz_pt2[i]);
        ((RooCBShape*)pdf_lz_pt2[i])->fitTo(*ds_lz_pt2[i]);
    }
}

void clearPDFs() {
    delete var_cb_alpha_ww_pt1;
    delete var_cb_m0_ww_pt1;
    delete var_cb_n_ww_pt1;
    delete var_cb_sigma_ww_pt1;
    delete var_cb_alpha_tt_pt1;
    delete var_cb_m0_tt_pt1;
    delete var_cb_n_tt_pt1;
    delete var_cb_sigma_tt_pt1;
    delete (RooCBShape*)pdf_ww_pt1;
    delete (RooCBShape*)pdf_tt_pt1;
    delete var_cb_alpha_ww_pt2;
    delete var_cb_m0_ww_pt2;
    delete var_cb_n_ww_pt2;
    delete var_cb_sigma_ww_pt2;
    delete var_cb_alpha_tt_pt2;
    delete var_cb_m0_tt_pt2;
    delete var_cb_n_tt_pt2;
    delete var_cb_sigma_tt_pt2;
    delete (RooCBShape*)pdf_ww_pt2;
    delete (RooCBShape*)pdf_tt_pt2;

    for(unsigned int i = 0; i < 20; ++i) {
        delete var_cb_alpha_lz_pt1[i];
        delete var_cb_m0_lz_pt1[i];
        delete var_cb_n_lz_pt1[i];
        delete var_cb_sigma_lz_pt1[i];
        delete (RooCBShape*)pdf_lz_pt1[i];
        delete var_cb_alpha_lz_pt2[i];
        delete var_cb_m0_lz_pt2[i];
        delete var_cb_n_lz_pt2[i];
        delete var_cb_sigma_lz_pt2[i];
        delete (RooCBShape*)pdf_lz_pt2[i];
    }
}
