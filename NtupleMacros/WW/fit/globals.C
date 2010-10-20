#include "globals.h"

#include "RooDataSet.h"
#include "RooRealVar.h"
#include "TFile.h"

void setVars() {
    var_deta  = new RooRealVar("deta","deta",0.,4.5);
    var_dilpt = new RooRealVar("dilpt","dilpt",0.,500.);
    var_dphi  = new RooRealVar("dphi","dphi",0.,3.15);
    var_mass  = new RooRealVar("mass","mass",12.,800.);
    var_met   = new RooRealVar("met","met",20.,500.);
    var_pt1   = new RooRealVar("pt1","pt1",20.,500.);
    var_pt2   = new RooRealVar("pt2","pt2",20.,400.);
    var_dummy = new RooRealVar("var_dummy","var_dummy",0);

    var_deta->setBins(40);
    var_dilpt->setBins(40);
    var_dphi->setBins(40);
    var_mass->setBins(40);
    var_met->setBins(40);
    var_pt1->setBins(40);
    var_pt2->setBins(40);
}

void setDataSets() {
    TFile fww("/tas/jribnik/processed_wwsmp1j.root");
    TFile ftt("/tas/jribnik/dataset_ttbarll.root");
    TFile flz("/tas/jribnik/processed_lambdaz.root");

    ds_ww = (RooAbsData*)fww.Get("ww");
    ds_ww->SetName("ds_ww");

    ds_ww_deta = ds_ww->reduce(*var_deta,"selected==1&&unique==1");
    ds_ww_deta->SetName("ds_ww_deta");
    ds_ww_dilpt = ds_ww->reduce(*var_dilpt,"selected==1&&unique==1");
    ds_ww_dilpt->SetName("ds_ww_dilpt");
    ds_ww_dphi = ds_ww->reduce(*var_dphi,"selected==1&&unique==1");
    ds_ww_dphi->SetName("ds_ww_dphi");
    ds_ww_mass = ds_ww->reduce(*var_mass,"selected==1&&unique==1");
    ds_ww_mass->SetName("ds_ww_mass");
    ds_ww_met = ds_ww->reduce(*var_met,"selected==1&&unique==1");
    ds_ww_met->SetName("ds_ww_met");
    ds_ww_pt1 = ds_ww->reduce(*var_pt1,"selected==1&&unique==1");
    ds_ww_pt1->SetName("ds_ww_pt1");
    ds_ww_pt2 = ds_ww->reduce(*var_pt2,"selected==1&&unique==1");
    ds_ww_pt2->SetName("ds_ww_pt2");

    ds_ww_pt1vdeta = ds_ww->reduce(RooArgSet(*var_pt1,*var_deta),"selected==1&&unique==1");
    ds_ww_pt1vdeta->SetName("ds_ww_pt1vdeta");
    ds_ww_pt1vdilpt = ds_ww->reduce(RooArgSet(*var_pt1,*var_dilpt),"selected==1&&unique==1");
    ds_ww_pt1vdilpt->SetName("ds_ww_pt1vdilpt");
    ds_ww_pt1vdphi = ds_ww->reduce(RooArgSet(*var_pt1,*var_dphi),"selected==1&&unique==1");
    ds_ww_pt1vdphi->SetName("ds_ww_pt1vdphi");
    ds_ww_pt1vmass = ds_ww->reduce(RooArgSet(*var_pt1,*var_mass),"selected==1&&unique==1");
    ds_ww_pt1vmass->SetName("ds_ww_pt1vmass");
    ds_ww_pt1vmet = ds_ww->reduce(RooArgSet(*var_pt1,*var_met),"selected==1&&unique==1");
    ds_ww_pt1vmet->SetName("ds_ww_pt1vmet");
    ds_ww_pt1vpt2 = ds_ww->reduce(RooArgSet(*var_pt1,*var_pt2),"selected==1&&unique==1");
    ds_ww_pt1vpt2->SetName("ds_ww_pt1vpt2");

    ds_tt = (RooAbsData*)ftt.Get("ttbar");
    ds_tt->SetName("ds_tt");

    ds_tt_deta = ds_tt->reduce(*var_deta,"refv0==1");
    ds_tt_deta->SetName("ds_tt_deta");
    ds_tt_dilpt = ds_tt->reduce(*var_dilpt,"refv0==1");
    ds_tt_dilpt->SetName("ds_tt_dilpt");
    ds_tt_dphi = ds_tt->reduce(*var_dphi,"refv0==1");
    ds_tt_dphi->SetName("ds_tt_dphi");
    ds_tt_mass = ds_tt->reduce(*var_mass,"refv0==1");
    ds_tt_mass->SetName("ds_tt_mass");
    ds_tt_met = ds_tt->reduce(*var_met,"refv0==1");
    ds_tt_met->SetName("ds_tt_met");
    ds_tt_pt1 = ds_tt->reduce(*var_pt1,"refv0==1");
    ds_tt_pt1->SetName("ds_tt_pt1");
    ds_tt_pt2 = ds_tt->reduce(*var_pt2,"refv0==1");
    ds_tt_pt2->SetName("ds_tt_pt2");

    ds_tt_pt1vdeta = ds_tt->reduce(RooArgSet(*var_pt1,*var_deta),"refv0==1");
    ds_tt_pt1vdeta->SetName("ds_tt_pt1vdeta");
    ds_tt_pt1vdilpt = ds_tt->reduce(RooArgSet(*var_pt1,*var_dilpt),"refv0==1");
    ds_tt_pt1vdilpt->SetName("ds_tt_pt1vdilpt");
    ds_tt_pt1vdphi = ds_tt->reduce(RooArgSet(*var_pt1,*var_dphi),"refv0==1");
    ds_tt_pt1vdphi->SetName("ds_tt_pt1vdphi");
    ds_tt_pt1vmass = ds_tt->reduce(RooArgSet(*var_pt1,*var_mass),"refv0==1");
    ds_tt_pt1vmass->SetName("ds_tt_pt1vmass");
    ds_tt_pt1vmet = ds_tt->reduce(RooArgSet(*var_pt1,*var_met),"refv0==1");
    ds_tt_pt1vmet->SetName("ds_tt_pt1vmet");
    ds_tt_pt1vpt2 = ds_tt->reduce(RooArgSet(*var_pt1,*var_pt2),"refv0==1");
    ds_tt_pt1vpt2->SetName("ds_tt_pt1vpt2");

    /*
    // lambdaz points
    RooDataSet* ds_lzs = (RooDataSet*)flz.Get("ww");
    ds_lzs->SetName("ds_lzs");

    float lambdazs[20] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5};

    for(unsigned int i = 0; i < 20; ++i) {
        float lo = lambdazs[i]-0.0001;
        float hi = lambdazs[i]+0.0001;
        ds_lz[i] = ds_lzs->reduce(Form("lambdaz>%f&&lambdaz<%f",lo,hi));
        ds_lz[i]->SetName(Form("ds_lz%i",i));

        ds_lz_deta[i] = ds_lzs->reduce(*var_deta,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_deta[i]->SetName(Form("ds_lz%i_deta",i)); 
        ds_lz_dilpt[i] = ds_lzs->reduce(*var_dilpt,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_dilpt[i]->SetName(Form("ds_lz%i_dilpt",i)); 
        ds_lz_dphi[i] = ds_lzs->reduce(*var_dphi,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_dphi[i]->SetName(Form("ds_lz%i_dphi",i)); 
        ds_lz_mass[i] = ds_lzs->reduce(*var_mass,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_mass[i]->SetName(Form("ds_lz%i_mass",i)); 
        ds_lz_met[i] = ds_lzs->reduce(*var_met,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_met[i]->SetName(Form("ds_lz%i_met",i)); 
        ds_lz_pt1[i] = ds_lzs->reduce(*var_pt1,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1[i]->SetName(Form("ds_lz%i_pt1",i)); 
        ds_lz_pt2[i] = ds_lzs->reduce(*var_pt2,Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt2[i]->SetName(Form("ds_lz%i_pt2",i)); 

        ds_lz_pt1vdeta[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_deta),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vdeta[i]->SetName(Form("ds_lz%i_pt1vdeta",i)); 
        ds_lz_pt1vdilpt[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_dilpt),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vdilpt[i]->SetName(Form("ds_lz%i_pt1vdilpt",i)); 
        ds_lz_pt1vdphi[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_dphi),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vdphi[i]->SetName(Form("ds_lz%i_pt1vdphi",i)); 
        ds_lz_pt1vmass[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_mass),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vmass[i]->SetName(Form("ds_lz%i_pt1vmass",i)); 
        ds_lz_pt1vmet[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_met),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vmet[i]->SetName(Form("ds_lz%i_pt1vmet",i)); 
        ds_lz_pt1vpt2[i] = ds_lzs->reduce(RooArgSet(*var_pt1,*var_pt2),Form("lambdaz>%f&&lambdaz<%f&&selected==1&&unique==1",lo,hi));
        ds_lz_pt1vpt2[i]->SetName(Form("ds_lz%i_pt1vpt2",i)); 
    }
    */
}

void clearDataSets() {
    delete ds_ww;
    delete ds_ww_deta;
    delete ds_ww_dilpt;
    delete ds_ww_dphi;
    delete ds_ww_mass;
    delete ds_ww_met;
    delete ds_ww_pt1;
    delete ds_ww_pt2;
    delete ds_ww_pt1vdeta;
    delete ds_ww_pt1vdilpt;
    delete ds_ww_pt1vdphi;
    delete ds_ww_pt1vmass;
    delete ds_ww_pt1vmet;
    delete ds_ww_pt1vpt2;
    delete ds_tt;
    delete ds_tt_deta;
    delete ds_tt_dilpt;
    delete ds_tt_dphi;
    delete ds_tt_mass;
    delete ds_tt_met;
    delete ds_tt_pt1;
    delete ds_tt_pt2;
    delete ds_tt_pt1vdeta;
    delete ds_tt_pt1vdilpt;
    delete ds_tt_pt1vdphi;
    delete ds_tt_pt1vmass;
    delete ds_tt_pt1vmet;
    delete ds_tt_pt1vpt2;

    for(unsigned int i = 0; i < 20; ++i) {
        delete ds_lz[i];
        delete ds_lz_deta[i];
        delete ds_lz_dilpt[i];
        delete ds_lz_dphi[i];
        delete ds_lz_mass[i];
        delete ds_lz_met[i];
        delete ds_lz_pt1[i];
        delete ds_lz_pt2[i];
        delete ds_lz_pt1vdeta[i];
        delete ds_lz_pt1vdilpt[i];
        delete ds_lz_pt1vdphi[i];
        delete ds_lz_pt1vmass[i];
        delete ds_lz_pt1vmet[i];
        delete ds_lz_pt1vpt2[i];
    }
}
