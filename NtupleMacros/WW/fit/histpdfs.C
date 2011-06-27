#include "globals.h"
#include "histpdfs.h"

#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooNDKeysPdf.h"
#include "RooRealVar.h"
#include "TH1F.h"
#include "TH2F.h"

void setPDFs() {
    // build hist pdfs from keys pdfs
    // probably there is a better way to do this...
    RooNDKeysPdf keys_ww_deta("keys_ww_deta","keys_ww_deta",*var_deta,*((RooDataSet*)ds_ww_deta));
    RooNDKeysPdf keys_ww_dilpt("keys_ww_dilpt","keys_ww_dilpt",*var_dilpt,*((RooDataSet*)ds_ww_dilpt));
    RooNDKeysPdf keys_ww_dphi("keys_ww_dphi","keys_ww_dphi",*var_dphi,*((RooDataSet*)ds_ww_dphi));
    RooNDKeysPdf keys_ww_mass("keys_ww_mass","keys_ww_mass",*var_mass,*((RooDataSet*)ds_ww_mass));
    RooNDKeysPdf keys_ww_met("keys_ww_met","keys_ww_met",*var_met,*((RooDataSet*)ds_ww_met));
    RooNDKeysPdf keys_ww_pt1("keys_ww_pt1","keys_ww_pt1",*var_pt1,*((RooDataSet*)ds_ww_pt1));
    RooNDKeysPdf keys_ww_pt2("keys_ww_pt2","keys_ww_pt2",*var_pt2,*((RooDataSet*)ds_ww_pt2));

    TH1F* h_keys_ww_deta = (TH1F*)keys_ww_deta.createHistogram("h_keys_ww_deta",*var_deta);
    TH1F* h_keys_ww_dilpt = (TH1F*)keys_ww_dilpt.createHistogram("h_keys_ww_dilpt",*var_dilpt);
    TH1F* h_keys_ww_dphi = (TH1F*)keys_ww_dphi.createHistogram("h_keys_ww_dphi",*var_dphi);
    TH1F* h_keys_ww_mass = (TH1F*)keys_ww_mass.createHistogram("h_keys_ww_mass",*var_mass);
    TH1F* h_keys_ww_met = (TH1F*)keys_ww_met.createHistogram("h_keys_ww_met",*var_met);
    TH1F* h_keys_ww_pt1 = (TH1F*)keys_ww_pt1.createHistogram("h_keys_ww_pt1",*var_pt1);
    TH1F* h_keys_ww_pt2 = (TH1F*)keys_ww_pt2.createHistogram("h_keys_ww_pt2",*var_pt2);

    dh_ww_deta = new RooDataHist("dh_ww_deta","dh_ww_deta",*var_deta,h_keys_ww_deta);
    dh_ww_dilpt = new RooDataHist("dh_ww_dilpt","dh_ww_dilpt",*var_dilpt,h_keys_ww_dilpt);
    dh_ww_dphi = new RooDataHist("dh_ww_dphi","dh_ww_dphi",*var_dphi,h_keys_ww_dphi);
    dh_ww_mass = new RooDataHist("dh_ww_mass","dh_ww_mass",*var_mass,h_keys_ww_mass);
    dh_ww_met = new RooDataHist("dh_ww_met","dh_ww_met",*var_met,h_keys_ww_met);
    dh_ww_pt1 = new RooDataHist("dh_ww_pt1","dh_ww_pt1",*var_pt1,h_keys_ww_pt1);
    dh_ww_pt2 = new RooDataHist("dh_ww_pt2","dh_ww_pt2",*var_pt2,h_keys_ww_pt2);

    pdf_ww_deta = new RooHistPdf("pdf_ww_deta","pdf_ww_deta",*var_deta,*dh_ww_deta);
    pdf_ww_dilpt = new RooHistPdf("pdf_ww_dilpt","pdf_ww_dilpt",*var_dilpt,*dh_ww_dilpt);
    pdf_ww_dphi = new RooHistPdf("pdf_ww_dphi","pdf_ww_dphi",*var_dphi,*dh_ww_dphi);
    pdf_ww_mass = new RooHistPdf("pdf_ww_mass","pdf_ww_mass",*var_mass,*dh_ww_mass);
    pdf_ww_met = new RooHistPdf("pdf_ww_met","pdf_ww_met",*var_met,*dh_ww_met);
    pdf_ww_pt1 = new RooHistPdf("pdf_ww_pt1","pdf_ww_pt1",*var_pt1,*dh_ww_pt1);
    pdf_ww_pt2 = new RooHistPdf("pdf_ww_pt2","pdf_ww_pt2",*var_pt2,*dh_ww_pt2);

    delete h_keys_ww_deta;
    delete h_keys_ww_dilpt;
    delete h_keys_ww_dphi;
    delete h_keys_ww_mass;
    delete h_keys_ww_met;
    delete h_keys_ww_pt1;
    delete h_keys_ww_pt2;

    RooNDKeysPdf keys_ww_pt1vdeta("keys_ww_pt1vdeta","keys_ww_pt1vdeta",RooArgList(*var_pt1,*var_deta),*((RooDataSet*)ds_ww_pt1vdeta));
    RooNDKeysPdf keys_ww_pt1vdilpt("keys_ww_pt1vdilpt","keys_ww_pt1vdilpt",RooArgList(*var_pt1,*var_dilpt),*((RooDataSet*)ds_ww_pt1vdilpt));
    RooNDKeysPdf keys_ww_pt1vdphi("keys_ww_pt1vdphi","keys_ww_pt1vdphi",RooArgList(*var_pt1,*var_dphi),*((RooDataSet*)ds_ww_pt1vdphi));
    RooNDKeysPdf keys_ww_pt1vmass("keys_ww_pt1vmass","keys_ww_pt1vmass",RooArgList(*var_pt1,*var_mass),*((RooDataSet*)ds_ww_pt1vmass));
    RooNDKeysPdf keys_ww_pt1vmet("keys_ww_pt1vmet","keys_ww_pt1vmet",RooArgList(*var_pt1,*var_met),*((RooDataSet*)ds_ww_pt1vmet));
    RooNDKeysPdf keys_ww_pt1vpt2("keys_ww_pt1vpt2","keys_ww_pt1vpt2",RooArgList(*var_pt1,*var_pt2),*((RooDataSet*)ds_ww_pt1vpt2));

    TH2F* h_keys_ww_pt1vdeta = (TH2F*)keys_ww_pt1vdeta.createHistogram("h_keys_ww_pt1vdeta",*var_pt1,RooFit::YVar(*var_deta));
    TH2F* h_keys_ww_pt1vdilpt = (TH2F*)keys_ww_pt1vdilpt.createHistogram("h_keys_ww_pt1vdilpt",*var_pt1,RooFit::YVar(*var_dilpt));
    TH2F* h_keys_ww_pt1vdphi = (TH2F*)keys_ww_pt1vdphi.createHistogram("h_keys_ww_pt1vdphi",*var_pt1,RooFit::YVar(*var_dphi));
    TH2F* h_keys_ww_pt1vmass = (TH2F*)keys_ww_pt1vmass.createHistogram("h_keys_ww_pt1vmass",*var_pt1,RooFit::YVar(*var_mass));
    TH2F* h_keys_ww_pt1vmet = (TH2F*)keys_ww_pt1vmet.createHistogram("h_keys_ww_pt1vmet",*var_pt1,RooFit::YVar(*var_met));
    TH2F* h_keys_ww_pt1vpt2 = (TH2F*)keys_ww_pt1vpt2.createHistogram("h_keys_ww_pt1vpt2",*var_pt1,RooFit::YVar(*var_pt2));

    dh_ww_pt1vdeta = new RooDataHist("dh_ww_pt1vdeta","dh_ww_pt1vdeta",RooArgSet(*var_pt1,*var_deta),h_keys_ww_pt1vdeta);
    dh_ww_pt1vdilpt = new RooDataHist("dh_ww_pt1vdilpt","dh_ww_pt1vdilpt",RooArgSet(*var_pt1,*var_dilpt),h_keys_ww_pt1vdilpt);
    dh_ww_pt1vdphi = new RooDataHist("dh_ww_pt1vdphi","dh_ww_pt1vdphi",RooArgSet(*var_pt1,*var_dphi),h_keys_ww_pt1vdphi);
    dh_ww_pt1vmass = new RooDataHist("dh_ww_pt1vmass","dh_ww_pt1vmass",RooArgSet(*var_pt1,*var_mass),h_keys_ww_pt1vmass);
    dh_ww_pt1vmet = new RooDataHist("dh_ww_pt1vmet","dh_ww_pt1vmet",RooArgSet(*var_pt1,*var_met),h_keys_ww_pt1vmet);
    dh_ww_pt1vpt2 = new RooDataHist("dh_ww_pt1vpt2","dh_ww_pt1vpt2",RooArgSet(*var_pt1,*var_pt2),h_keys_ww_pt1vpt2);

    pdf_ww_pt1vdeta = new RooHistPdf("pdf_ww_pt1vdeta","pdf_ww_pt1vdeta",RooArgSet(*var_pt1,*var_deta),*dh_ww_pt1vdeta);
    pdf_ww_pt1vdilpt = new RooHistPdf("pdf_ww_pt1vdilpt","pdf_ww_pt1vdilpt",RooArgSet(*var_pt1,*var_dilpt),*dh_ww_pt1vdilpt);
    pdf_ww_pt1vdphi = new RooHistPdf("pdf_ww_pt1vdphi","pdf_ww_pt1vdphi",RooArgSet(*var_pt1,*var_dphi),*dh_ww_pt1vdphi);
    pdf_ww_pt1vmass = new RooHistPdf("pdf_ww_pt1vmass","pdf_ww_pt1vmass",RooArgSet(*var_pt1,*var_mass),*dh_ww_pt1vmass);
    pdf_ww_pt1vmet = new RooHistPdf("pdf_ww_pt1vmet","pdf_ww_pt1vmet",RooArgSet(*var_pt1,*var_met),*dh_ww_pt1vmet);
    pdf_ww_pt1vpt2 = new RooHistPdf("pdf_ww_pt1vpt2","pdf_ww_pt1vpt2",RooArgSet(*var_pt1,*var_pt2),*dh_ww_pt1vpt2);

    delete h_keys_ww_pt1vdeta;
    delete h_keys_ww_pt1vdilpt;
    delete h_keys_ww_pt1vdphi;
    delete h_keys_ww_pt1vmass;
    delete h_keys_ww_pt1vmet;
    delete h_keys_ww_pt1vpt2;

    RooNDKeysPdf keys_tt_deta("keys_tt_deta","keys_tt_deta",*var_deta,*((RooDataSet*)ds_tt_deta));
    RooNDKeysPdf keys_tt_dilpt("keys_tt_dilpt","keys_tt_dilpt",*var_dilpt,*((RooDataSet*)ds_tt_dilpt));
    RooNDKeysPdf keys_tt_dphi("keys_tt_dphi","keys_tt_dphi",*var_dphi,*((RooDataSet*)ds_tt_dphi));
    RooNDKeysPdf keys_tt_mass("keys_tt_mass","keys_tt_mass",*var_mass,*((RooDataSet*)ds_tt_mass));
    RooNDKeysPdf keys_tt_met("keys_tt_met","keys_tt_met",*var_met,*((RooDataSet*)ds_tt_met));
    RooNDKeysPdf keys_tt_pt1("keys_tt_pt1","keys_tt_pt1",*var_pt1,*((RooDataSet*)ds_tt_pt1));
    RooNDKeysPdf keys_tt_pt2("keys_tt_pt2","keys_tt_pt2",*var_pt2,*((RooDataSet*)ds_tt_pt2));

    TH1F* h_keys_tt_deta = (TH1F*)keys_tt_deta.createHistogram("h_keys_tt_deta",*var_deta);
    TH1F* h_keys_tt_dilpt = (TH1F*)keys_tt_dilpt.createHistogram("h_keys_tt_dilpt",*var_dilpt);
    TH1F* h_keys_tt_dphi = (TH1F*)keys_tt_dphi.createHistogram("h_keys_tt_dphi",*var_dphi);
    TH1F* h_keys_tt_mass = (TH1F*)keys_tt_mass.createHistogram("h_keys_tt_mass",*var_mass);
    TH1F* h_keys_tt_met = (TH1F*)keys_tt_met.createHistogram("h_keys_tt_met",*var_met);
    TH1F* h_keys_tt_pt1 = (TH1F*)keys_tt_pt1.createHistogram("h_keys_tt_pt1",*var_pt1);
    TH1F* h_keys_tt_pt2 = (TH1F*)keys_tt_pt2.createHistogram("h_keys_tt_pt2",*var_pt2);

    dh_tt_deta = new RooDataHist("dh_tt_deta","dh_tt_deta",*var_deta,h_keys_tt_deta);
    dh_tt_dilpt = new RooDataHist("dh_tt_dilpt","dh_tt_dilpt",*var_dilpt,h_keys_tt_dilpt);
    dh_tt_dphi = new RooDataHist("dh_tt_dphi","dh_tt_dphi",*var_dphi,h_keys_tt_dphi);
    dh_tt_mass = new RooDataHist("dh_tt_mass","dh_tt_mass",*var_mass,h_keys_tt_mass);
    dh_tt_met = new RooDataHist("dh_tt_met","dh_tt_met",*var_met,h_keys_tt_met);
    dh_tt_pt1 = new RooDataHist("dh_tt_pt1","dh_tt_pt1",*var_pt1,h_keys_tt_pt1);
    dh_tt_pt2 = new RooDataHist("dh_tt_pt2","dh_tt_pt2",*var_pt2,h_keys_tt_pt2);

    pdf_tt_deta = new RooHistPdf("pdf_tt_deta","pdf_tt_deta",*var_deta,*dh_tt_deta);
    pdf_tt_dilpt = new RooHistPdf("pdf_tt_dilpt","pdf_tt_dilpt",*var_dilpt,*dh_tt_dilpt);
    pdf_tt_dphi = new RooHistPdf("pdf_tt_dphi","pdf_tt_dphi",*var_dphi,*dh_tt_dphi);
    pdf_tt_mass = new RooHistPdf("pdf_tt_mass","pdf_tt_mass",*var_mass,*dh_tt_mass);
    pdf_tt_met = new RooHistPdf("pdf_tt_met","pdf_tt_met",*var_met,*dh_tt_met);
    pdf_tt_pt1 = new RooHistPdf("pdf_tt_pt1","pdf_tt_pt1",*var_pt1,*dh_tt_pt1);
    pdf_tt_pt2 = new RooHistPdf("pdf_tt_pt2","pdf_tt_pt2",*var_pt2,*dh_tt_pt2);

    delete h_keys_tt_deta;
    delete h_keys_tt_dilpt;
    delete h_keys_tt_dphi;
    delete h_keys_tt_mass;
    delete h_keys_tt_met;
    delete h_keys_tt_pt1;
    delete h_keys_tt_pt2;

    RooNDKeysPdf keys_tt_pt1vdeta("keys_tt_pt1vdeta","keys_tt_pt1vdeta",RooArgList(*var_pt1,*var_deta),*((RooDataSet*)ds_tt_pt1vdeta));
    RooNDKeysPdf keys_tt_pt1vdilpt("keys_tt_pt1vdilpt","keys_tt_pt1vdilpt",RooArgList(*var_pt1,*var_dilpt),*((RooDataSet*)ds_tt_pt1vdilpt));
    RooNDKeysPdf keys_tt_pt1vdphi("keys_tt_pt1vdphi","keys_tt_pt1vdphi",RooArgList(*var_pt1,*var_dphi),*((RooDataSet*)ds_tt_pt1vdphi));
    RooNDKeysPdf keys_tt_pt1vmass("keys_tt_pt1vmass","keys_tt_pt1vmass",RooArgList(*var_pt1,*var_mass),*((RooDataSet*)ds_tt_pt1vmass));
    RooNDKeysPdf keys_tt_pt1vmet("keys_tt_pt1vmet","keys_tt_pt1vmet",RooArgList(*var_pt1,*var_met),*((RooDataSet*)ds_tt_pt1vmet));
    RooNDKeysPdf keys_tt_pt1vpt2("keys_tt_pt1vpt2","keys_tt_pt1vpt2",RooArgList(*var_pt1,*var_pt2),*((RooDataSet*)ds_tt_pt1vpt2));

    TH2F* h_keys_tt_pt1vdeta = (TH2F*)keys_tt_pt1vdeta.createHistogram("h_keys_tt_pt1vdeta",*var_pt1,RooFit::YVar(*var_deta));
    TH2F* h_keys_tt_pt1vdilpt = (TH2F*)keys_tt_pt1vdilpt.createHistogram("h_keys_tt_pt1vdilpt",*var_pt1,RooFit::YVar(*var_dilpt));
    TH2F* h_keys_tt_pt1vdphi = (TH2F*)keys_tt_pt1vdphi.createHistogram("h_keys_tt_pt1vdphi",*var_pt1,RooFit::YVar(*var_dphi));
    TH2F* h_keys_tt_pt1vmass = (TH2F*)keys_tt_pt1vmass.createHistogram("h_keys_tt_pt1vmass",*var_pt1,RooFit::YVar(*var_mass));
    TH2F* h_keys_tt_pt1vmet = (TH2F*)keys_tt_pt1vmet.createHistogram("h_keys_tt_pt1vmet",*var_pt1,RooFit::YVar(*var_met));
    TH2F* h_keys_tt_pt1vpt2 = (TH2F*)keys_tt_pt1vpt2.createHistogram("h_keys_tt_pt1vpt2",*var_pt1,RooFit::YVar(*var_pt2));

    dh_tt_pt1vdeta = new RooDataHist("dh_tt_pt1vdeta","dh_tt_pt1vdeta",RooArgSet(*var_pt1,*var_deta),h_keys_tt_pt1vdeta);
    dh_tt_pt1vdilpt = new RooDataHist("dh_tt_pt1vdilpt","dh_tt_pt1vdilpt",RooArgSet(*var_pt1,*var_dilpt),h_keys_tt_pt1vdilpt);
    dh_tt_pt1vdphi = new RooDataHist("dh_tt_pt1vdphi","dh_tt_pt1vdphi",RooArgSet(*var_pt1,*var_dphi),h_keys_tt_pt1vdphi);
    dh_tt_pt1vmass = new RooDataHist("dh_tt_pt1vmass","dh_tt_pt1vmass",RooArgSet(*var_pt1,*var_mass),h_keys_tt_pt1vmass);
    dh_tt_pt1vmet = new RooDataHist("dh_tt_pt1vmet","dh_tt_pt1vmet",RooArgSet(*var_pt1,*var_met),h_keys_tt_pt1vmet);
    dh_tt_pt1vpt2 = new RooDataHist("dh_tt_pt1vpt2","dh_tt_pt1vpt2",RooArgSet(*var_pt1,*var_pt2),h_keys_tt_pt1vpt2);

    pdf_tt_pt1vdeta = new RooHistPdf("pdf_tt_pt1vdeta","pdf_tt_pt1vdeta",RooArgSet(*var_pt1,*var_deta),*dh_tt_pt1vdeta);
    pdf_tt_pt1vdilpt = new RooHistPdf("pdf_tt_pt1vdilpt","pdf_tt_pt1vdilpt",RooArgSet(*var_pt1,*var_dilpt),*dh_tt_pt1vdilpt);
    pdf_tt_pt1vdphi = new RooHistPdf("pdf_tt_pt1vdphi","pdf_tt_pt1vdphi",RooArgSet(*var_pt1,*var_dphi),*dh_tt_pt1vdphi);
    pdf_tt_pt1vmass = new RooHistPdf("pdf_tt_pt1vmass","pdf_tt_pt1vmass",RooArgSet(*var_pt1,*var_mass),*dh_tt_pt1vmass);
    pdf_tt_pt1vmet = new RooHistPdf("pdf_tt_pt1vmet","pdf_tt_pt1vmet",RooArgSet(*var_pt1,*var_met),*dh_tt_pt1vmet);
    pdf_tt_pt1vpt2 = new RooHistPdf("pdf_tt_pt1vpt2","pdf_tt_pt1vpt2",RooArgSet(*var_pt1,*var_pt2),*dh_tt_pt1vpt2);

    delete h_keys_tt_pt1vdeta;
    delete h_keys_tt_pt1vdilpt;
    delete h_keys_tt_pt1vdphi;
    delete h_keys_tt_pt1vmass;
    delete h_keys_tt_pt1vmet;
    delete h_keys_tt_pt1vpt2;

    /*
    for(unsigned int i = 0; i < 20; ++i) {
        RooNDKeysPdf keys_lz_deta("keys_lz_deta","keys_lz_deta",*var_deta,*((RooDataSet*)ds_lz_deta[i]));
        RooNDKeysPdf keys_lz_dilpt("keys_lz_dilpt","keys_lz_dilpt",*var_dilpt,*((RooDataSet*)ds_lz_dilpt[i]));
        RooNDKeysPdf keys_lz_dphi("keys_lz_dphi","keys_lz_dphi",*var_dphi,*((RooDataSet*)ds_lz_dphi[i]));
        RooNDKeysPdf keys_lz_mass("keys_lz_mass","keys_lz_mass",*var_mass,*((RooDataSet*)ds_lz_mass[i]));
        RooNDKeysPdf keys_lz_met("keys_lz_met","keys_lz_met",*var_met,*((RooDataSet*)ds_lz_met[i]));
        RooNDKeysPdf keys_lz_pt1("keys_lz_pt1","keys_lz_pt1",*var_pt1,*((RooDataSet*)ds_lz_pt1[i]));
        RooNDKeysPdf keys_lz_pt2("keys_lz_pt2","keys_lz_pt2",*var_pt2,*((RooDataSet*)ds_lz_pt2[i]));

        TH1F* h_keys_lz_deta = (TH1F*)keys_lz_deta.createHistogram("h_keys_lz_deta",*var_deta);
        TH1F* h_keys_lz_dilpt = (TH1F*)keys_lz_dilpt.createHistogram("h_keys_lz_dilpt",*var_dilpt);
        TH1F* h_keys_lz_dphi = (TH1F*)keys_lz_dphi.createHistogram("h_keys_lz_dphi",*var_dphi);
        TH1F* h_keys_lz_mass = (TH1F*)keys_lz_mass.createHistogram("h_keys_lz_mass",*var_mass);
        TH1F* h_keys_lz_met = (TH1F*)keys_lz_met.createHistogram("h_keys_lz_met",*var_met);
        TH1F* h_keys_lz_pt1 = (TH1F*)keys_lz_pt1.createHistogram("h_keys_lz_pt1",*var_pt1);
        TH1F* h_keys_lz_pt2 = (TH1F*)keys_lz_pt2.createHistogram("h_keys_lz_pt2",*var_pt2);

        dh_lz_deta[i] = new RooDataHist(Form("dh_lz%i_deta",i),Form("dh_lz%i_deta",i),*var_deta,h_keys_lz_deta);
        dh_lz_dilpt[i] = new RooDataHist(Form("dh_lz%i_dilpt",i),Form("dh_lz%i_dilpt",i),*var_dilpt,h_keys_lz_dilpt);
        dh_lz_dphi[i] = new RooDataHist(Form("dh_lz%i_dphi",i),Form("dh_lz%i_dphi",i),*var_dphi,h_keys_lz_dphi);
        dh_lz_mass[i] = new RooDataHist(Form("dh_lz%i_mass",i),Form("dh_lz%i_mass",i),*var_mass,h_keys_lz_mass);
        dh_lz_met[i] = new RooDataHist(Form("dh_lz%i_met",i),Form("dh_lz%i_met",i),*var_met,h_keys_lz_met);
        dh_lz_pt1[i] = new RooDataHist(Form("dh_lz%i_pt1",i),Form("dh_lz%i_pt1",i),*var_pt1,h_keys_lz_pt1);
        dh_lz_pt2[i] = new RooDataHist(Form("dh_lz%i_pt2",i),Form("dh_lz%i_pt2",i),*var_pt2,h_keys_lz_pt2);

        pdf_lz_deta[i] = new RooHistPdf(Form("pdf_lz%i_deta",i),Form("pdf_lz%i_deta",i),*var_deta,*dh_lz_deta[i]);
        pdf_lz_dilpt[i] = new RooHistPdf(Form("pdf_lz%i_dilpt",i),Form("pdf_lz%i_dilpt",i),*var_dilpt,*dh_lz_dilpt[i]);
        pdf_lz_dphi[i] = new RooHistPdf(Form("pdf_lz%i_dphi",i),Form("pdf_lz%i_dphi",i),*var_dphi,*dh_lz_dphi[i]);
        pdf_lz_mass[i] = new RooHistPdf(Form("pdf_lz%i_mass",i),Form("pdf_lz%i_mass",i),*var_mass,*dh_lz_mass[i]);
        pdf_lz_met[i] = new RooHistPdf(Form("pdf_lz%i_met",i),Form("pdf_lz%i_met",i),*var_met,*dh_lz_met[i]);
        pdf_lz_pt1[i] = new RooHistPdf(Form("pdf_lz%i_pt1",i),Form("pdf_lz%i_pt1",i),*var_pt1,*dh_lz_pt1[i]);
        pdf_lz_pt2[i] = new RooHistPdf(Form("pdf_lz%i_pt2",i),Form("pdf_lz%i_pt2",i),*var_pt2,*dh_lz_pt2[i]);

        delete h_keys_lz_deta;
        delete h_keys_lz_dilpt;
        delete h_keys_lz_dphi;
        delete h_keys_lz_mass;
        delete h_keys_lz_met;
        delete h_keys_lz_pt1;
        delete h_keys_lz_pt2;

        RooNDKeysPdf keys_lz_pt1vdeta("keys_lz_pt1vdeta","keys_lz_pt1vdeta",RooArgList(*var_pt1,*var_deta),*((RooDataSet*)ds_lz_pt1vdeta[i]));
        RooNDKeysPdf keys_lz_pt1vdilpt("keys_lz_pt1vdilpt","keys_lz_pt1vdilpt",RooArgList(*var_pt1,*var_dilpt),*((RooDataSet*)ds_lz_pt1vdilpt[i]));
        RooNDKeysPdf keys_lz_pt1vdphi("keys_lz_pt1vdphi","keys_lz_pt1vdphi",RooArgList(*var_pt1,*var_dphi),*((RooDataSet*)ds_lz_pt1vdphi[i]));
        RooNDKeysPdf keys_lz_pt1vmass("keys_lz_pt1vmass","keys_lz_pt1vmass",RooArgList(*var_pt1,*var_mass),*((RooDataSet*)ds_lz_pt1vmass[i]));
        RooNDKeysPdf keys_lz_pt1vmet("keys_lz_pt1vmet","keys_lz_pt1vmet",RooArgList(*var_pt1,*var_met),*((RooDataSet*)ds_lz_pt1vmet[i]));
        RooNDKeysPdf keys_lz_pt1vpt2("keys_lz_pt1vpt2","keys_lz_pt1vpt2",RooArgList(*var_pt1,*var_pt2),*((RooDataSet*)ds_lz_pt1vpt2[i]));

        TH2F* h_keys_lz_pt1vdeta = (TH2F*)keys_lz_pt1vdeta.createHistogram("h_keys_lz_pt1vdeta",*var_pt1,RooFit::YVar(*var_deta));
        TH2F* h_keys_lz_pt1vdilpt = (TH2F*)keys_lz_pt1vdilpt.createHistogram("h_keys_lz_pt1vdilpt",*var_pt1,RooFit::YVar(*var_dilpt));
        TH2F* h_keys_lz_pt1vdphi = (TH2F*)keys_lz_pt1vdphi.createHistogram("h_keys_lz_pt1vdphi",*var_pt1,RooFit::YVar(*var_dphi));
        TH2F* h_keys_lz_pt1vmass = (TH2F*)keys_lz_pt1vmass.createHistogram("h_keys_lz_pt1vmass",*var_pt1,RooFit::YVar(*var_mass));
        TH2F* h_keys_lz_pt1vmet = (TH2F*)keys_lz_pt1vmet.createHistogram("h_keys_lz_pt1vmet",*var_pt1,RooFit::YVar(*var_met));
        TH2F* h_keys_lz_pt1vpt2 = (TH2F*)keys_lz_pt1vpt2.createHistogram("h_keys_lz_pt1vpt2",*var_pt1,RooFit::YVar(*var_pt2));

        dh_lz_pt1vdeta[i] = new RooDataHist(Form("dh_lz%i_pt1vdeta",i),Form("dh_lz%i_pt1vdeta",i),RooArgSet(*var_pt1,*var_deta),h_keys_lz_pt1vdeta);
        dh_lz_pt1vdilpt[i] = new RooDataHist(Form("dh_lz%i_pt1vdilpt",i),Form("dh_lz%i_pt1vdilpt",i),RooArgSet(*var_pt1,*var_dilpt),h_keys_lz_pt1vdilpt);
        dh_lz_pt1vdphi[i] = new RooDataHist(Form("dh_lz%i_pt1vdphi",i),Form("dh_lz%i_pt1vdphi",i),RooArgSet(*var_pt1,*var_dphi),h_keys_lz_pt1vdphi);
        dh_lz_pt1vmass[i] = new RooDataHist(Form("dh_lz%i_pt1vmass",i),Form("dh_lz%i_pt1vmass",i),RooArgSet(*var_pt1,*var_mass),h_keys_lz_pt1vmass);
        dh_lz_pt1vmet[i] = new RooDataHist(Form("dh_lz%i_pt1vmet",i),Form("dh_lz%i_pt1vmet",i),RooArgSet(*var_pt1,*var_met),h_keys_lz_pt1vmet);
        dh_lz_pt1vpt2[i] = new RooDataHist(Form("dh_lz%i_pt1vpt2",i),Form("dh_lz%i_pt1vpt2",i),RooArgSet(*var_pt1,*var_pt2),h_keys_lz_pt1vpt2);

        pdf_lz_pt1vdeta[i] = new RooHistPdf(Form("pdf_lz%i_pt1vdeta",i),Form("pdf_lz%i_pt1vdeta",i),RooArgSet(*var_pt1,*var_deta),*dh_lz_pt1vdeta[i]);
        pdf_lz_pt1vdilpt[i] = new RooHistPdf(Form("pdf_lz%i_pt1vdilpt",i),Form("pdf_lz%i_pt1vdilpt",i),RooArgSet(*var_pt1,*var_dilpt),*dh_lz_pt1vdilpt[i]);
        pdf_lz_pt1vdphi[i] = new RooHistPdf(Form("pdf_lz%i_pt1vdphi",i),Form("pdf_lz%i_pt1vdphi",i),RooArgSet(*var_pt1,*var_dphi),*dh_lz_pt1vdphi[i]);
        pdf_lz_pt1vmass[i] = new RooHistPdf(Form("pdf_lz%i_pt1vmass",i),Form("pdf_lz%i_pt1vmass",i),RooArgSet(*var_pt1,*var_mass),*dh_lz_pt1vmass[i]);
        pdf_lz_pt1vmet[i] = new RooHistPdf(Form("pdf_lz%i_pt1vmet",i),Form("pdf_lz%i_pt1vmet",i),RooArgSet(*var_pt1,*var_met),*dh_lz_pt1vmet[i]);
        pdf_lz_pt1vpt2[i] = new RooHistPdf(Form("pdf_lz%i_pt1vpt2",i),Form("pdf_lz%i_pt1vpt2",i),RooArgSet(*var_pt1,*var_pt2),*dh_lz_pt1vpt2[i]);

        delete h_keys_lz_pt1vdeta;
        delete h_keys_lz_pt1vdilpt;
        delete h_keys_lz_pt1vdphi;
        delete h_keys_lz_pt1vmass;
        delete h_keys_lz_pt1vmet;
        delete h_keys_lz_pt1vpt2;
    }
*/
}

void clearPDFs() {
    delete dh_ww_deta;
    delete dh_ww_dilpt;
    delete dh_ww_dphi;
    delete dh_ww_mass;
    delete dh_ww_met;
    delete dh_ww_pt1;
    delete dh_ww_pt2;
    delete dh_ww_pt1vdeta;
    delete dh_ww_pt1vdilpt;
    delete dh_ww_pt1vdphi;
    delete dh_ww_pt1vmass;
    delete dh_ww_pt1vmet;
    delete dh_ww_pt1vpt2;
    delete dh_tt_deta;
    delete dh_tt_dilpt;
    delete dh_tt_dphi;
    delete dh_tt_mass;
    delete dh_tt_met;
    delete dh_tt_pt1;
    delete dh_tt_pt2;
    delete dh_tt_pt1vdeta;
    delete dh_tt_pt1vdilpt;
    delete dh_tt_pt1vdphi;
    delete dh_tt_pt1vmass;
    delete dh_tt_pt1vmet;
    delete dh_tt_pt1vpt2;
    delete (RooHistPdf*)pdf_ww_deta;
    delete (RooHistPdf*)pdf_ww_dilpt;
    delete (RooHistPdf*)pdf_ww_dphi;
    delete (RooHistPdf*)pdf_ww_mass;
    delete (RooHistPdf*)pdf_ww_met;
    delete (RooHistPdf*)pdf_ww_pt1;
    delete (RooHistPdf*)pdf_ww_pt2;
    delete (RooHistPdf*)pdf_ww_pt1vdeta;
    delete (RooHistPdf*)pdf_ww_pt1vdilpt;
    delete (RooHistPdf*)pdf_ww_pt1vdphi;
    delete (RooHistPdf*)pdf_ww_pt1vmass;
    delete (RooHistPdf*)pdf_ww_pt1vmet;
    delete (RooHistPdf*)pdf_ww_pt1vpt2;
    delete (RooHistPdf*)pdf_tt_deta;
    delete (RooHistPdf*)pdf_tt_dilpt;
    delete (RooHistPdf*)pdf_tt_dphi;
    delete (RooHistPdf*)pdf_tt_mass;
    delete (RooHistPdf*)pdf_tt_met;
    delete (RooHistPdf*)pdf_tt_pt1;
    delete (RooHistPdf*)pdf_tt_pt2;
    delete (RooHistPdf*)pdf_tt_pt1vdeta;
    delete (RooHistPdf*)pdf_tt_pt1vdilpt;
    delete (RooHistPdf*)pdf_tt_pt1vdphi;
    delete (RooHistPdf*)pdf_tt_pt1vmass;
    delete (RooHistPdf*)pdf_tt_pt1vmet;
    delete (RooHistPdf*)pdf_tt_pt1vpt2;

    for(unsigned int i = 0; i < 20; ++i) {
        delete dh_lz_deta[i];
        delete dh_lz_dilpt[i];
        delete dh_lz_dphi[i];
        delete dh_lz_mass[i];
        delete dh_lz_met[i];
        delete dh_lz_pt1[i];
        delete dh_lz_pt2[i];
        delete dh_lz_pt1vdeta[i];
        delete dh_lz_pt1vdilpt[i];
        delete dh_lz_pt1vdphi[i];
        delete dh_lz_pt1vmass[i];
        delete dh_lz_pt1vmet[i];
        delete dh_lz_pt1vpt2[i];
        delete (RooHistPdf*)pdf_lz_deta[i];
        delete (RooHistPdf*)pdf_lz_dilpt[i];
        delete (RooHistPdf*)pdf_lz_dphi[i];
        delete (RooHistPdf*)pdf_lz_mass[i];
        delete (RooHistPdf*)pdf_lz_met[i];
        delete (RooHistPdf*)pdf_lz_pt1[i];
        delete (RooHistPdf*)pdf_lz_pt2[i];
        delete (RooHistPdf*)pdf_lz_pt1vdeta[i];
        delete (RooHistPdf*)pdf_lz_pt1vdilpt[i];
        delete (RooHistPdf*)pdf_lz_pt1vdphi[i];
        delete (RooHistPdf*)pdf_lz_pt1vmass[i];
        delete (RooHistPdf*)pdf_lz_pt1vmet[i];
        delete (RooHistPdf*)pdf_lz_pt1vpt2[i];
    }
}
