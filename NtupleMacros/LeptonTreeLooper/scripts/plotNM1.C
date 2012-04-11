
{

    //
    // load lins
    //

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gROOT->ProcessLine(".L ../libLeptonTreeLooper.so");

    gROOT->ProcessLine(".L ~/tdrStyle.C");
    gROOT->ProcessLine("setTDRStyle()");
    gStyle->SetOptTitle(1);
    gROOT->ForceStyle();

    //
    // open file
    //

    TFile *f = new TFile("../histos.root", "READ");
    gROOT->cd();

    //
    // make plots
    //

    bool drawNormalized = true;
    TString prefix = "test";

    // general plots
    TCanvas *c_nvtx_tp = ComparePlots(f, "data_2011B_h1_tp_nvtx", "data_2012A_h1_tp_nvtx", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_tagAndProbeMass_tp = ComparePlots(f, "data_2011B_h1_tp_tagAndProbeMass", "data_2012A_h1_tp_tagAndProbeMass", "2011B", "2012A", 1, drawNormalized);
    c_nvtx_tp->SaveAs("../plots/"+prefix+"_nvtx.png");
    c_tagAndProbeMass_tp->SaveAs("../plots/"+prefix+"_tagAndProbeMass.png");

    // variables
    TCanvas *c_detain_tp = ComparePlots(f, "data_2011B_h1_tp_detain", "data_2012A_h1_tp_detain", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_dphiin_tp = ComparePlots(f, "data_2011B_h1_tp_dphiin", "data_2012A_h1_tp_dphiin", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_sieie_tp = ComparePlots(f, "data_2011B_h1_tp_sieie", "data_2012A_h1_tp_sieie", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_hoe_tp = ComparePlots(f, "data_2011B_h1_tp_hoe", "data_2012A_h1_tp_hoe", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_ooemoop_tp = ComparePlots(f, "data_2011B_h1_tp_ooemoop", "data_2012A_h1_tp_ooemoop", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_d0vtx_tp = ComparePlots(f, "data_2011B_h1_tp_d0vtx", "data_2012A_h1_tp_d0vtx", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_dzvtx_tp = ComparePlots(f, "data_2011B_h1_tp_dzvtx", "data_2012A_h1_tp_dzvtx", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_vfitprob_tp = ComparePlots(f, "data_2011B_h1_tp_vfitprob", "data_2012A_h1_tp_vfitprob", "2011B", "2012A", 1, drawNormalized);
    TCanvas *c_mhit_tp = ComparePlots(f, "data_2011B_h1_tp_mhit", "data_2012A_h1_tp_mhit", "2011B", "2012A", 1, drawNormalized);

    c_detain_tp->SaveAs("../plots/"+prefix+"_detain.png");
    c_dphiin_tp->SaveAs("../plots/"+prefix+"_dphiin.png");
    c_sieie_tp->SaveAs("../plots/"+prefix+"_sieie.png");
    c_hoe_tp->SaveAs("../plots/"+prefix+"_hoe.png");
    c_ooemoop_tp->SaveAs("../plots/"+prefix+"_ooemoop.png");
    c_d0vtx_tp->SaveAs("../plots/"+prefix+"_d0vtx.png");
    c_dzvtx_tp->SaveAs("../plots/"+prefix+"_dzvtx.png");
    c_vfitprob_tp->SaveAs("../plots/"+prefix+"_vfitprob.png");
    c_mhit_tp->SaveAs("../plots/"+prefix+"_mhit.png");
    

}


