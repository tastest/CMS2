#include "globals.h"

#include "RooDataSet.h"
#include "TFile.h"

void setDataSets() {
    TFile fww("/tas/jribnik/processed_wwsmp1j.root");
    TFile ftt("/tas/jribnik/dataset_ttbarll.root");
    TFile flz("/tas/jribnik/processed_lambdaz.root");

    ds_ww_pt = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt,"selected==1");
    ds_ww_pt->SetName("ds_ww_pt");
    ds_tt_pt = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt,"refv0==1");
    ds_tt_pt->SetName("ds_tt_pt");
    //ds_ww_pt_njv = ((RooAbsData*)fww.Get("ww"))->reduce(*var_pt);
    //ds_ww_pt_njv->SetName("ds_ww_pt_njv");
    //ds_tt_pt_njv = ((RooAbsData*)ftt.Get("ttbar"))->reduce(*var_pt);
    //ds_tt_pt_njv->SetName("ds_tt_pt_njv");

    // lambdaz points
    RooDataSet* ds_lz = (RooDataSet*)flz.Get("ww");

    float lambdazs[20] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,
        -0.05,-0.1,-0.15,-0.2,-0.25,-0.3,-0.35,-0.4,-0.45,-0.5};

    for(unsigned int i = 0; i < 20; ++i) {
        float lo = lambdazs[i]-0.0001;
        float hi = lambdazs[i]+0.0001;
        ds_lz_pt[i] = ds_lz->reduce(*var_pt,Form("lambdaz>%f&&lambdaz<%f&&selected==1",lo,hi));
        ds_lz_pt[i]->SetName(Form("ds_lz%i_pt",i)); 
        //ds_lz_pt_njv[i] = ds_lz->reduce(*var_pt,Form("lambdaz>%f&&lambdaz<%f",lo,hi));
        //ds_lz_pt_njv[i]->SetName(Form("ds_lz%i_pt_njv",i)); 
    }
}

void clearDataSets() {
    delete ds_ww_pt;
    delete ds_tt_pt;
    //delete ds_ww_pt_njv;
    //delete ds_tt_pt_njv;

    for(unsigned int i = 0; i < 20; ++i) {
        delete ds_lz_pt[i];
        //delete ds_lz_pt_njv[i];
    }
}
