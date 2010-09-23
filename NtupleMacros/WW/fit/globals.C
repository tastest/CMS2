#include "globals.h"

#include "RooDataSet.h"
#include "TFile.h"

void setDataSets() {
    TFile fww("/tas/jribnik/processed_wwsmp1j.root");
    TFile ftt("/tas/jribnik/dataset_ttbarll.root");
    TFile flz("/tas/jribnik/processed_lambdaz.root");

    ds_ww_deta = ((RooAbsData*)fww.Get("ww"))->reduce(*var_deta,"selected==1");
    ds_ww_deta->SetName("ds_ww_deta");
    ds_tt_deta = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_deta,"refv0==1");
    ds_tt_deta->SetName("ds_tt_deta");

    ds_ww_dilpt = ((RooAbsData*)fww.Get("ww"))->reduce(*var_dilpt,"selected==1");
    ds_ww_dilpt->SetName("ds_ww_dilpt");
    ds_tt_dilpt = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_dilpt,"refv0==1");
    ds_tt_dilpt->SetName("ds_tt_dilpt");

    ds_ww_dphi = ((RooAbsData*)fww.Get("ww"))->reduce(*var_dphi,"selected==1");
    ds_ww_dphi->SetName("ds_ww_dphi");
    ds_tt_dphi = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_dphi,"refv0==1");
    ds_tt_dphi->SetName("ds_tt_dphi");

    ds_ww_mass = ((RooAbsData*)fww.Get("ww"))->reduce(*var_mass,"selected==1");
    ds_ww_mass->SetName("ds_ww_mass");
    ds_tt_mass = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_mass,"refv0==1");
    ds_tt_mass->SetName("ds_tt_mass");

    ds_ww_met = ((RooAbsData*)fww.Get("ww"))->reduce(*var_met,"selected==1");
    ds_ww_met->SetName("ds_ww_met");
    ds_tt_met = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_met,"refv0==1");
    ds_tt_met->SetName("ds_tt_met");

    ds_ww_pt1 = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt1,"selected==1");
    ds_ww_pt1->SetName("ds_ww_pt1");
    ds_tt_pt1 = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt1,"refv0==1");
    ds_tt_pt1->SetName("ds_tt_pt1");
    //ds_ww_pt1_njv = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt1);
    //ds_ww_pt1_njv->SetName("ds_ww_pt1_njv");
    //ds_tt_pt1_njv = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt1);
    //ds_tt_pt1_njv->SetName("ds_tt_pt1_njv");

    ds_ww_pt2 = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt2,"selected==1");
    ds_ww_pt2->SetName("ds_ww_pt2");
    ds_tt_pt2 = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt2,"refv0==1");
    ds_tt_pt2->SetName("ds_tt_pt2");
    //ds_ww_pt2_njv = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt2);
    //ds_ww_pt2_njv->SetName("ds_ww_pt2_njv");
    //ds_tt_pt2_njv = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt2);
    //ds_tt_pt2_njv->SetName("ds_tt_pt2_njv");

    // lambdaz points
    RooDataSet* ds_lz = (RooDataSet*)flz.Get("ww");

    float lambdazs[20] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5};

    for(unsigned int i = 0; i < 20; ++i) {
        float lo = lambdazs[i]-0.0001;
        float hi = lambdazs[i]+0.0001;
        ds_lz_deta[i] = ds_lz->reduce(*var_deta,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_deta[i]->SetName(Form("ds_lz%i_deta",i)); 
        ds_lz_dilpt[i] = ds_lz->reduce(*var_dilpt,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_dilpt[i]->SetName(Form("ds_lz%i_dilpt",i)); 
        ds_lz_dphi[i] = ds_lz->reduce(*var_dphi,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_dphi[i]->SetName(Form("ds_lz%i_dphi",i)); 
        ds_lz_mass[i] = ds_lz->reduce(*var_mass,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_mass[i]->SetName(Form("ds_lz%i_mass",i)); 
        ds_lz_met[i] = ds_lz->reduce(*var_met,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_met[i]->SetName(Form("ds_lz%i_met",i)); 
        ds_lz_pt1[i] = ds_lz->reduce(*var_pt1,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_pt1[i]->SetName(Form("ds_lz%i_pt1",i)); 
        //ds_lz_pt1_njv[i] = ds_lz->reduce(*var_pt1,Form("lambdaz>%f&&lambdaz<%f",lo,hi));
        //ds_lz_pt1_njv[i]->SetName(Form("ds_lz%i_pt1_njv",i)); 
        ds_lz_pt2[i] = ds_lz->reduce(*var_pt2,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_pt2[i]->SetName(Form("ds_lz%i_pt2",i)); 
        //ds_lz_pt2_njv[i] = ds_lz->reduce(*var_pt2,Form("lambdaz>%f&&lambdaz<%f",lo,hi));
        //ds_lz_pt2_njv[i]->SetName(Form("ds_lz%i_pt2_njv",i)); 
    }
}

void clearDataSets() {
    delete ds_ww_deta;
    delete ds_tt_deta;
    delete ds_ww_dilpt;
    delete ds_tt_dilpt;
    delete ds_ww_dphi;
    delete ds_tt_dphi;
    delete ds_ww_mass;
    delete ds_tt_mass;
    delete ds_ww_met;
    delete ds_tt_met;
    delete ds_ww_pt1;
    delete ds_tt_pt1;
    //delete ds_ww_pt1_njv;
    //delete ds_tt_pt1_njv;
    delete ds_ww_pt2;
    delete ds_tt_pt2;
    //delete ds_ww_pt2_njv;
    //delete ds_tt_pt2_njv;

    for(unsigned int i = 0; i < 20; ++i) {
        delete ds_lz_deta[i];
        delete ds_lz_dilpt[i];
        delete ds_lz_dphi[i];
        delete ds_lz_mass[i];
        delete ds_lz_met[i];
        delete ds_lz_pt1[i];
        //delete ds_lz_pt1_njv[i];
        delete ds_lz_pt2[i];
        //delete ds_lz_pt2_njv[i];
    }
}
