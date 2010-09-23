#include "globals.h"
#include "histpdfs.h"

#include "RooDataHist.h"
#include "RooHistPdf.h"

void setPDFs() {
    dh_ww_deta = new RooDataHist("dh_ww_deta","dh_ww_deta",*var_deta,*ds_ww_deta);
    dh_tt_deta = new RooDataHist("dh_tt_deta","dh_tt_deta",*var_deta,*ds_tt_deta);

    dh_ww_dilpt = new RooDataHist("dh_ww_dilpt","dh_ww_dilpt",*var_dilpt,*ds_ww_dilpt);
    dh_tt_dilpt = new RooDataHist("dh_tt_dilpt","dh_tt_dilpt",*var_dilpt,*ds_tt_dilpt);

    dh_ww_dphi = new RooDataHist("dh_ww_dphi","dh_ww_dphi",*var_dphi,*ds_ww_dphi);
    dh_tt_dphi = new RooDataHist("dh_tt_dphi","dh_tt_dphi",*var_dphi,*ds_tt_dphi);

    dh_ww_mass = new RooDataHist("dh_ww_mass","dh_ww_mass",*var_mass,*ds_ww_mass);
    dh_tt_mass = new RooDataHist("dh_tt_mass","dh_tt_mass",*var_mass,*ds_tt_mass);

    dh_ww_met = new RooDataHist("dh_ww_met","dh_ww_met",*var_met,*ds_ww_met);
    dh_tt_met = new RooDataHist("dh_tt_met","dh_tt_met",*var_met,*ds_tt_met);

    dh_ww_pt1 = new RooDataHist("dh_ww_pt1","dh_ww_pt1",*var_pt1,*ds_ww_pt1);
    dh_tt_pt1 = new RooDataHist("dh_tt_pt1","dh_tt_pt1",*var_pt1,*ds_tt_pt1);
    //dh_ww_pt1_njv = new RooDataHist("dh_ww_pt1_njv","dh_ww_pt1_njv",*var_pt1,*ds_ww_pt1_njv);
    //dh_tt_pt1_njv = new RooDataHist("dh_tt_pt1_njv","dh_tt_pt1_njv",*var_pt1,*ds_tt_pt1_njv);

    dh_ww_pt2 = new RooDataHist("dh_ww_pt2","dh_ww_pt2",*var_pt2,*ds_ww_pt2);
    dh_tt_pt2 = new RooDataHist("dh_tt_pt2","dh_tt_pt2",*var_pt2,*ds_tt_pt2);
    //dh_ww_pt2_njv = new RooDataHist("dh_ww_pt2_njv","dh_ww_pt2_njv",*var_pt2,*ds_ww_pt2_njv);
    //dh_tt_pt2_njv = new RooDataHist("dh_tt_pt2_njv","dh_tt_pt2_njv",*var_pt2,*ds_tt_pt2_njv);

    pdf_ww_deta = new RooHistPdf("pdf_ww_deta","pdf_ww_deta",*var_deta,*dh_ww_deta);
    pdf_tt_deta = new RooHistPdf("pdf_tt_deta","pdf_tt_deta",*var_deta,*dh_tt_deta);

    pdf_ww_dilpt = new RooHistPdf("pdf_ww_dilpt","pdf_ww_dilpt",*var_dilpt,*dh_ww_dilpt);
    pdf_tt_dilpt = new RooHistPdf("pdf_tt_dilpt","pdf_tt_dilpt",*var_dilpt,*dh_tt_dilpt);

    pdf_ww_dphi = new RooHistPdf("pdf_ww_dphi","pdf_ww_dphi",*var_dphi,*dh_ww_dphi);
    pdf_tt_dphi = new RooHistPdf("pdf_tt_dphi","pdf_tt_dphi",*var_dphi,*dh_tt_dphi);

    pdf_ww_mass = new RooHistPdf("pdf_ww_mass","pdf_ww_mass",*var_mass,*dh_ww_mass);
    pdf_tt_mass = new RooHistPdf("pdf_tt_mass","pdf_tt_mass",*var_mass,*dh_tt_mass);

    pdf_ww_met = new RooHistPdf("pdf_ww_met","pdf_ww_met",*var_met,*dh_ww_met);
    pdf_tt_met = new RooHistPdf("pdf_tt_met","pdf_tt_met",*var_met,*dh_tt_met);

    pdf_ww_pt1 = new RooHistPdf("pdf_ww_pt1","pdf_ww_pt1",*var_pt1,*dh_ww_pt1);
    pdf_tt_pt1 = new RooHistPdf("pdf_tt_pt1","pdf_tt_pt1",*var_pt1,*dh_tt_pt1);
    //pdf_ww_pt1_njv = new RooHistPdf("pdf_ww_pt1_njv","pdf_ww_pt1_njv",*var_pt1,*dh_ww_pt1_njv);
    //pdf_tt_pt1_njv = new RooHistPdf("pdf_tt_pt1_njv","pdf_tt_pt1_njv",*var_pt1,*dh_tt_pt1_njv);

    pdf_ww_pt2 = new RooHistPdf("pdf_ww_pt2","pdf_ww_pt2",*var_pt2,*dh_ww_pt2);
    pdf_tt_pt2 = new RooHistPdf("pdf_tt_pt2","pdf_tt_pt2",*var_pt2,*dh_tt_pt2);
    //pdf_ww_pt2_njv = new RooHistPdf("pdf_ww_pt2_njv","pdf_ww_pt2_njv",*var_pt2,*dh_ww_pt2_njv);
    //pdf_tt_pt2_njv = new RooHistPdf("pdf_tt_pt2_njv","pdf_tt_pt2_njv",*var_pt2,*dh_tt_pt2_njv);

    for(unsigned int i = 0; i < 20; ++i) {
        dh_lz_deta[i] = new RooDataHist(Form("dh_lz%i_deta",i),Form("dh_lz%i_deta",i),*var_deta,*ds_lz_deta[i]);
        dh_lz_dilpt[i] = new RooDataHist(Form("dh_lz%i_dilpt",i),Form("dh_lz%i_dilpt",i),*var_dilpt,*ds_lz_dilpt[i]);
        dh_lz_dphi[i] = new RooDataHist(Form("dh_lz%i_dphi",i),Form("dh_lz%i_dphi",i),*var_dphi,*ds_lz_dphi[i]);
        dh_lz_mass[i] = new RooDataHist(Form("dh_lz%i_mass",i),Form("dh_lz%i_mass",i),*var_mass,*ds_lz_mass[i]);
        dh_lz_met[i] = new RooDataHist(Form("dh_lz%i_met",i),Form("dh_lz%i_met",i),*var_met,*ds_lz_met[i]);
        dh_lz_pt1[i] = new RooDataHist(Form("dh_lz%i_pt1",i),Form("dh_lz%i_pt1",i),*var_pt1,*ds_lz_pt1[i]);
        //dh_lz_pt1_njv[i] = new RooDataHist(Form("dh_lz%i_pt1_njv",i),Form("dh_lz%i_pt1_njv",i),*var_pt1,*ds_lz_pt1_njv[i]);
        dh_lz_pt2[i] = new RooDataHist(Form("dh_lz%i_pt2",i),Form("dh_lz%i_pt2",i),*var_pt2,*ds_lz_pt2[i]);
        //dh_lz_pt2_njv[i] = new RooDataHist(Form("dh_lz%i_pt2_njv",i),Form("dh_lz%i_pt2_njv",i),*var_pt2,*ds_lz_pt2_njv[i]);

        pdf_lz_deta[i] = new RooHistPdf(Form("pdf_lz%i_deta",i),Form("pdf_lz%i_deta",i),*var_deta,*dh_lz_deta[i]);
        pdf_lz_dilpt[i] = new RooHistPdf(Form("pdf_lz%i_dilpt",i),Form("pdf_lz%i_dilpt",i),*var_dilpt,*dh_lz_dilpt[i]);
        pdf_lz_dphi[i] = new RooHistPdf(Form("pdf_lz%i_dphi",i),Form("pdf_lz%i_dphi",i),*var_dphi,*dh_lz_dphi[i]);
        pdf_lz_mass[i] = new RooHistPdf(Form("pdf_lz%i_mass",i),Form("pdf_lz%i_mass",i),*var_mass,*dh_lz_mass[i]);
        pdf_lz_met[i] = new RooHistPdf(Form("pdf_lz%i_met",i),Form("pdf_lz%i_met",i),*var_met,*dh_lz_met[i]);
        pdf_lz_pt1[i] = new RooHistPdf(Form("pdf_lz%i_pt1",i),Form("pdf_lz%i_pt1",i),*var_pt1,*dh_lz_pt1[i]);
        //pdf_lz_pt1_njv[i] = new RooHistPdf(Form("pdf_lz%i_pt1_njv",i),Form("pdf_lz%i_pt1_njv",i),*var_pt1,*dh_lz_pt1_njv[i]);
        pdf_lz_pt2[i] = new RooHistPdf(Form("pdf_lz%i_pt2",i),Form("pdf_lz%i_pt2",i),*var_pt2,*dh_lz_pt2[i]);
        //pdf_lz_pt2_njv[i] = new RooHistPdf(Form("pdf_lz%i_pt2_njv",i),Form("pdf_lz%i_pt2_njv",i),*var_pt2,*dh_lz_pt2_njv[i]);
    }
}

void clearPDFs() {
    delete dh_ww_deta;
    delete dh_tt_deta;
    delete dh_ww_dilpt;
    delete dh_tt_dilpt;
    delete dh_ww_dphi;
    delete dh_tt_dphi;
    delete dh_ww_mass;
    delete dh_tt_mass;
    delete dh_ww_met;
    delete dh_tt_met;
    delete dh_ww_pt1;
    delete dh_tt_pt1;
    //delete dh_ww_pt1_njv;
    //delete dh_tt_pt1_njv;
    delete dh_ww_pt2;
    delete dh_tt_pt2;
    //delete dh_ww_pt2_njv;
    //delete dh_tt_pt2_njv;
    delete (RooHistPdf*)pdf_ww_deta;
    delete (RooHistPdf*)pdf_tt_deta;
    delete (RooHistPdf*)pdf_ww_dilpt;
    delete (RooHistPdf*)pdf_tt_dilpt;
    delete (RooHistPdf*)pdf_ww_dphi;
    delete (RooHistPdf*)pdf_tt_dphi;
    delete (RooHistPdf*)pdf_ww_mass;
    delete (RooHistPdf*)pdf_tt_mass;
    delete (RooHistPdf*)pdf_ww_met;
    delete (RooHistPdf*)pdf_tt_met;
    delete (RooHistPdf*)pdf_ww_pt1;
    delete (RooHistPdf*)pdf_tt_pt1;
    //delete (RooHistPdf*)pdf_ww_pt1_njv;
    //delete (RooHistPdf*)pdf_tt_pt1_njv;
    delete (RooHistPdf*)pdf_ww_pt2;
    delete (RooHistPdf*)pdf_tt_pt2;
    //delete (RooHistPdf*)pdf_ww_pt2_njv;
    //delete (RooHistPdf*)pdf_tt_pt2_njv;

    for(unsigned int i = 0; i < 20; ++i) {
        delete dh_lz_deta[i];
        delete dh_lz_dilpt[i];
        delete dh_lz_dphi[i];
        delete dh_lz_mass[i];
        delete dh_lz_met[i];
        delete dh_lz_pt1[i];
        //delete dh_lz_pt1_njv[i];
        delete dh_lz_pt2[i];
        //delete dh_lz_pt2_njv[i];
        delete (RooHistPdf*)pdf_lz_deta[i];
        delete (RooHistPdf*)pdf_lz_dilpt[i];
        delete (RooHistPdf*)pdf_lz_dphi[i];
        delete (RooHistPdf*)pdf_lz_mass[i];
        delete (RooHistPdf*)pdf_lz_met[i];
        delete (RooHistPdf*)pdf_lz_pt1[i];
        //delete (RooHistPdf*)pdf_lz_pt1_njv[i];
        delete (RooHistPdf*)pdf_lz_pt2[i];
        //delete (RooHistPdf*)pdf_lz_pt2_njv[i];
    }
}
