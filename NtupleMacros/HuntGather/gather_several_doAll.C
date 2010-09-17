#include "BabySample.h"
#include "cuts.h"
#include "gather.h"
#include "goodrun.h"

#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"

#include <iostream>

void gather_several_doAll( TCut thisBaseSel = inclusivez_dilep )
{
    gROOT->SetStyle("Plain");
    gStyle->SetHistMinimumZero();
    //gStyle->SetOptStat("nemoui");
    gStyle->SetOptStat(0);

    // this sets the json file obviously; it's
    // easiest to just  soft link to  the real
    // file as json.txt
    std::cout << "Using json.txt for goodruns\n";
    set_goodrun_file_json("json.txt");

    // calculate  integrated luminosity  in /fb
    // note that I input /nb and so get out /nb
    // then scale to /fb which  is what  I need
    // 1e-3*GetIntLumi(/pb)  is the  same thing
    // doesn't matter so long as it corresponds
    // to the lumi of the json file and that it
    // is correctly scaled to /fb afterward
    float f_intlumifb = 1e-6*GetIntLumi(2790);
    std::cout << "Integrated luminosity: " << f_intlumifb << "/fb\n";

    //
    // data babies
    //

    BabySample *bs_data_emu    = new BabySample("data","/tas05/disk00/jribnik/hunt/emu_baby/*.root","",1.,true);
    BabySample *bs_data_dilep  = new BabySample("data","/tas05/disk00/jribnik/hunt/dilep_baby/*.root","",1.,true);
    BabySample *bs_data_trilep = new BabySample("data","/tas05/disk00/jribnik/hunt/trilep_baby/*.root","",1.,true);

    //
    // mc babies
    //

    float kttbarjets = 157.5/165.;
    float ksingletop = 1.;
    float kvvjets    = 1.;
    float kwjets     = 31314./28049.;
    float kzjets     = 3048./2400.;
    float kzll       = 1666./1300.;
    float kdyll      = 3457./2659.;

    // emu
    BabySample *bs_zjets_emu       = new BabySample("zjets","/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/emu_baby/*.root","",kzjets,false,kAzure-2,1001);

    // dilep
    BabySample *bs_ttbarjets_dilep = new BabySample("ttbar","/tas05/disk00/jribnik/huntmc/TTbarJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kttbarjets,false,kRed+1,1001);
    BabySample *bs_singletop_dilep = new BabySample("tW","/tas05/disk00/jribnik/huntmc/SingleTop_tWChannel-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",ksingletop,false,kMagenta,1001);

    BabySample *bs_vvjets_dilep    = new BabySample("vvjets","/tas05/disk00/jribnik/huntmc/VVJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kvvjets,false,10,1001);
    BabySample *bs_wjets_dilep     = new BabySample("wjets","/tas05/disk00/jribnik/huntmc/WJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kwjets,false,kGreen-3,1001);
    BabySample *bs_ztautau_dilep   = new BabySample("ztautau","/tas05/disk00/jribnik/huntmc/Ztautau_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","mass<50",kzll,false,kAzure+8,1001);

    // Note that a common prefix means
    // a common histogram when used in
    // the same DrawAll, i.e. the five
    // samples below are combined
    BabySample *bs_zjets_dilep     = new BabySample("zll","/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kzjets,false,kAzure-2,1001);
    BabySample *bs_zee_dilep       = new BabySample("zll","/tas05/disk00/jribnik/huntmc/Zee_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","mass<50",kzll,false,kAzure-2,1001);
    BabySample *bs_zmumu_dilep     = new BabySample("zll","/tas05/disk00/jribnik/huntmc/Zmumu_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","mass<50",kzll,false,kAzure-2,1001);
    BabySample *bs_dyee_dilep      = new BabySample("zll","/tas05/disk00/jribnik/huntmc/DYee_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kdyll,false,kAzure-2,1001);
    BabySample *bs_dymumu_dilep    = new BabySample("zll","/tas05/disk00/jribnik/huntmc/DYmumu_M10to20_Spring10-START3X_V26_S09-v1/dilep_baby/*.root","",kdyll,false,kAzure-2,1001);

    // trilep
    BabySample *bs_zjets_trilep    = new BabySample("zjets","/tas05/disk00/jribnik/huntmc/ZJets-madgraph_Spring10-START3X_V26_S09-v1/trilep_baby/*.root","",kzjets,false,1001);

    //
    // time to plot!
    //
    TString filename = "";

    /*
       TCanvas *c_dileptonictopv4_dilep_mnjets = DrawAll( "njets", dileptonictopv4_dilep_mnjets, f_intlumifb, 11, -0.5, 10.5, false,
       bs_data_dilep      ,
       bs_ttbarjets_dilep ,
       bs_singletop_dilep ,
       bs_dyee_dilep      ,
       bs_dymumu_dilep    ,
       bs_vvjets_dilep    ,
       bs_wjets_dilep     ,
       bs_zjets_dilep     ,
       bs_zee_dilep       ,
       bs_zmumu_dilep     ,
       bs_ztautau_dilep  
       );
       filename = "c_dileptonictopv4_dilep_njets.png"
       c_dileptonictopv4_dilep_mnjets->SaveAs(filename);
     */

    
    TCanvas *c1 = DrawAll( "mass", thisBaseSel, f_intlumifb, 100, 70, 170, false,
            bs_data_dilep      ,
            bs_ttbarjets_dilep ,
            bs_singletop_dilep ,
            bs_dyee_dilep      ,
            bs_dymumu_dilep    ,
            bs_vvjets_dilep    ,
            bs_wjets_dilep     ,
            bs_zjets_dilep     ,
            bs_zee_dilep       ,
            bs_zmumu_dilep     ,
            bs_ztautau_dilep  
            );
    filename = thisBaseSel.GetName();
    filename.Append("_mass.png");
    c1->SaveAs(filename);
    
     TCanvas *c1a = DrawAll( "mass", thisBaseSel, f_intlumifb, 100, 70, 170, true,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
    filename = thisBaseSel.GetName();
    filename.Append("_mass_int.png");
     c1a->SaveAs(filename);

     TCanvas *c2 = DrawAll( "mt2", thisBaseSel, f_intlumifb, 100, 0, 120, true,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
     c2->SetLogy(1);
    filename = thisBaseSel.GetName();
    filename.Append("_mt2_int.png");
     c2->SaveAs(filename);

     TCanvas *c3 = DrawAll( "tcmet", thisBaseSel, f_intlumifb, 100, 0, 400, true,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
     c3->SetLogy(1);
    filename = thisBaseSel.GetName();
    filename.Append("_tcmet_int.png");
    c3->SaveAs(filename);

     // Histogram names do not currently distinguish between
     // standard and integrated histograms
     // Names must be unique before we can uncomment this as
     // we support filling existing histograms
     /*
     TCanvas *c4 = DrawAll( "tcmet", thisBaseSel, f_intlumifb, 100, 0, 400, false,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
             );
     c4->SetLogy(1);
    filename = thisBaseSel.GetName();
    filename.Append("_tcmet.png");
     c4->SaveAs(filename);
     */

     TCanvas *c5 = DrawAll( "njetsClean", thisBaseSel, f_intlumifb, 15, 0, 15, false,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
     c5->SetLogy(1);
    filename = thisBaseSel.GetName();
    filename.Append("_njetsClean.png");
     c5->SaveAs(filename);

     TCanvas *c6 = DrawAll( "dilpt", thisBaseSel, f_intlumifb, 200, 0, 200, false,
             bs_data_dilep      ,
             bs_ttbarjets_dilep ,
             bs_singletop_dilep ,
             bs_dyee_dilep      ,
             bs_dymumu_dilep    ,
             bs_vvjets_dilep    ,
             bs_wjets_dilep     ,
             bs_zjets_dilep     ,
             bs_zee_dilep       ,
             bs_zmumu_dilep     ,
             bs_ztautau_dilep  
             );
     c6->SetLogy(1);
    filename = thisBaseSel.GetName();
    filename.Append("_dilpt.png");
    c6->SaveAs(filename);
}
