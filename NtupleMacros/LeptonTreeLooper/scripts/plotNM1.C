
void plotNM1 (TString inputFile, TString proc1, TString proc2, unsigned int opt)
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
    //gStyle->SetOptTitle(1);
    gROOT->ForceStyle();

    //
    // open file
    //

    TFile *f = new TFile("../"+inputFile+".root", "READ");
    TString prefix = inputFile + "_" + proc1 + "_" + proc2;
    gROOT->cd();

    //
    // make plots
    //

    bool drawNormalized = true;

    // general plots
    TCanvas *c_tagAndProbeMass_tp = ComparePlots(f, proc1 + "_h1_tp_tagAndProbeMass", proc2 + "_h1_tp_tagAndProbeMass", proc1, proc2, 1, drawNormalized, false, opt);
    c_tagAndProbeMass_tp->SaveAs("../plots/"+prefix+"_tp_tagAndProbeMass.png");

    // variables in tag and probe
    TCanvas *c_nvtxraw_tp = ComparePlots(f, proc1 + "_h1_tp_nvtxraw", proc2 + "_h1_tp_nvtxraw", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_nvtx_tp = ComparePlots(f, proc1 + "_h1_tp_nvtx", proc2 + "_h1_tp_nvtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_detain_tp = ComparePlots(f, proc1 + "_h1_tp_detain", proc2 + "_h1_tp_detain", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_dphiin_tp = ComparePlots(f, proc1 + "_h1_tp_dphiin", proc2 + "_h1_tp_dphiin", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_sieie_tp = ComparePlots(f, proc1 + "_h1_tp_sieie", proc2 + "_h1_tp_sieie", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_hoe_tp = ComparePlots(f, proc1 + "_h1_tp_hoe", proc2 + "_h1_tp_hoe", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_hoetow_tp = ComparePlots(f, proc1 + "_h1_tp_hoetow", proc2 + "_h1_tp_hoetow", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ooemoop_tp = ComparePlots(f, proc1 + "_h1_tp_ooemoop", proc2 + "_h1_tp_ooemoop", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_d0vtx_tp = ComparePlots(f, proc1 + "_h1_tp_d0vtx", proc2 + "_h1_tp_d0vtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_dzvtx_tp = ComparePlots(f, proc1 + "_h1_tp_dzvtx", proc2 + "_h1_tp_dzvtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_vfitprob_tp = ComparePlots(f, proc1 + "_h1_tp_vfitprob", proc2 + "_h1_tp_vfitprob", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_mhit_tp = ComparePlots(f, proc1 + "_h1_tp_mhit", proc2 + "_h1_tp_mhit", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfchiso_tp = ComparePlots(f, proc1 + "_h1_tp_pfchiso", proc2 + "_h1_tp_pfchiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfemiso_tp = ComparePlots(f, proc1 + "_h1_tp_pfemiso", proc2 + "_h1_tp_pfemiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfnhiso_tp = ComparePlots(f, proc1 + "_h1_tp_pfnhiso", proc2 + "_h1_tp_pfnhiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfchisomvacut_tp = ComparePlots(f, proc1 + "_h1_tp_pfchisomvacut", proc2 + "_h1_tp_pfchisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfemisomvacut_tp = ComparePlots(f, proc1 + "_h1_tp_pfemisomvacut", proc2 + "_h1_tp_pfemisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfnhisomvacut_tp = ComparePlots(f, proc1 + "_h1_tp_pfnhisomvacut", proc2 + "_h1_tp_pfnhisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ecaliso_tp = ComparePlots(f, proc1 + "_h1_tp_ecaliso", proc2 + "_h1_tp_ecaliso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ecalisoval_tp = ComparePlots(f, proc1 + "_h1_tp_ecalisoval", proc2 + "_h1_tp_ecalisoval", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_hcaliso_tp = ComparePlots(f, proc1 + "_h1_tp_hcaliso", proc2 + "_h1_tp_hcaliso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_trackiso_tp = ComparePlots(f, proc1 + "_h1_tp_trackiso", proc2 + "_h1_tp_trackiso", proc1, proc2, 1, drawNormalized, true, opt);
    c_nvtxraw_tp->SaveAs("../plots/"+prefix+"_tp_nvtxraw.png");
    c_nvtx_tp->SaveAs("../plots/"+prefix+"_tp_nvtx.png");
    c_detain_tp->SaveAs("../plots/"+prefix+"_tp_detain.png");
    c_dphiin_tp->SaveAs("../plots/"+prefix+"_tp_dphiin.png");
    c_sieie_tp->SaveAs("../plots/"+prefix+"_tp_sieie.png");
    c_hoe_tp->SaveAs("../plots/"+prefix+"_tp_hoe.png");
    c_hoetow_tp->SaveAs("../plots/"+prefix+"_tp_hoetow.png");
    c_ooemoop_tp->SaveAs("../plots/"+prefix+"_tp_ooemoop.png");
    c_d0vtx_tp->SaveAs("../plots/"+prefix+"_tp_d0vtx.png");
    c_dzvtx_tp->SaveAs("../plots/"+prefix+"_tp_dzvtx.png");
    c_vfitprob_tp->SaveAs("../plots/"+prefix+"_tp_vfitprob.png");
    c_mhit_tp->SaveAs("../plots/"+prefix+"_tp_mhit.png");
    c_pfchiso_tp->SaveAs("../plots/"+prefix+"_tp_pfchiso.png");
    c_pfemiso_tp->SaveAs("../plots/"+prefix+"_tp_pfemiso.png");
    c_pfnhiso_tp->SaveAs("../plots/"+prefix+"_tp_pfnhiso.png");
    c_pfchisomvacut_tp->SaveAs("../plots/"+prefix+"_tp_pfchisomvacut.png");
    c_pfemisomvacut_tp->SaveAs("../plots/"+prefix+"_tp_pfemisomvacut.png");
    c_pfnhisomvacut_tp->SaveAs("../plots/"+prefix+"_tp_pfnhisomvacut.png");
    c_ecaliso_tp->SaveAs("../plots/"+prefix+"_tp_ecaliso.png");
    c_ecalisoval_tp->SaveAs("../plots/"+prefix+"_tp_ecalisoval.png");
    c_hcaliso_tp->SaveAs("../plots/"+prefix+"_tp_hcaliso.png");
    c_trackiso_tp->SaveAs("../plots/"+prefix+"_tp_trackiso.png");
 
    // variables in FR
    TCanvas *c_nvtxraw_fr = ComparePlots(f, proc1 + "_h1_fr_nvtxraw", proc2 + "_h1_fr_nvtxraw", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_nvtx_fr = ComparePlots(f, proc1 + "_h1_fr_nvtx", proc2 + "_h1_fr_nvtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_detain_fr = ComparePlots(f, proc1 + "_h1_fr_detain", proc2 + "_h1_fr_detain", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_dphiin_fr = ComparePlots(f, proc1 + "_h1_fr_dphiin", proc2 + "_h1_fr_dphiin", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_sieie_fr = ComparePlots(f, proc1 + "_h1_fr_sieie", proc2 + "_h1_fr_sieie", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_hoe_fr = ComparePlots(f, proc1 + "_h1_fr_hoe", proc2 + "_h1_fr_hoe", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_hoetow_fr = ComparePlots(f, proc1 + "_h1_fr_hoetow", proc2 + "_h1_fr_hoetow", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ooemoop_fr = ComparePlots(f, proc1 + "_h1_fr_ooemoop", proc2 + "_h1_fr_ooemoop", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_d0vtx_fr = ComparePlots(f, proc1 + "_h1_fr_d0vtx", proc2 + "_h1_fr_d0vtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_dzvtx_fr = ComparePlots(f, proc1 + "_h1_fr_dzvtx", proc2 + "_h1_fr_dzvtx", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_vfitprob_fr = ComparePlots(f, proc1 + "_h1_fr_vfitprob", proc2 + "_h1_fr_vfitprob", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_mhit_fr = ComparePlots(f, proc1 + "_h1_fr_mhit", proc2 + "_h1_fr_mhit", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfchiso_fr = ComparePlots(f, proc1 + "_h1_fr_pfchiso", proc2 + "_h1_fr_pfchiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfemiso_fr = ComparePlots(f, proc1 + "_h1_fr_pfemiso", proc2 + "_h1_fr_pfemiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfnhiso_fr = ComparePlots(f, proc1 + "_h1_fr_pfnhiso", proc2 + "_h1_fr_pfnhiso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfchisomvacut_fr = ComparePlots(f, proc1 + "_h1_fr_pfchisomvacut", proc2 + "_h1_fr_pfchisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfemisomvacut_fr = ComparePlots(f, proc1 + "_h1_fr_pfemisomvacut", proc2 + "_h1_fr_pfemisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_pfnhisomvacut_fr = ComparePlots(f, proc1 + "_h1_fr_pfnhisomvacut", proc2 + "_h1_fr_pfnhisomvacut", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ecaliso_fr = ComparePlots(f, proc1 + "_h1_fr_ecaliso", proc2 + "_h1_fr_ecaliso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_ecalisoval_fr = ComparePlots(f, proc1 + "_h1_fr_ecalisoval", proc2 + "_h1_fr_ecalisoval", proc1, proc2, 1, drawNormalized, false, opt);
    TCanvas *c_hcaliso_fr = ComparePlots(f, proc1 + "_h1_fr_hcaliso", proc2 + "_h1_fr_hcaliso", proc1, proc2, 1, drawNormalized, true, opt);
    TCanvas *c_trackiso_fr = ComparePlots(f, proc1 + "_h1_fr_trackiso", proc2 + "_h1_fr_trackiso", proc1, proc2, 1, drawNormalized, true, opt);
    c_nvtxraw_fr->SaveAs("../plots/"+prefix+"_fr_nvtxraw.png");
    c_nvtx_fr->SaveAs("../plots/"+prefix+"_fr_nvtx.png");
    c_detain_fr->SaveAs("../plots/"+prefix+"_fr_detain.png");
    c_dphiin_fr->SaveAs("../plots/"+prefix+"_fr_dphiin.png");
    c_sieie_fr->SaveAs("../plots/"+prefix+"_fr_sieie.png");
    c_hoe_fr->SaveAs("../plots/"+prefix+"_fr_hoe.png");
    c_hoetow_fr->SaveAs("../plots/"+prefix+"_fr_hoetow.png");
    c_ooemoop_fr->SaveAs("../plots/"+prefix+"_fr_ooemoop.png");
    c_d0vtx_fr->SaveAs("../plots/"+prefix+"_fr_d0vtx.png");
    c_dzvtx_fr->SaveAs("../plots/"+prefix+"_fr_dzvtx.png");
    c_vfitprob_fr->SaveAs("../plots/"+prefix+"_fr_vfitprob.png");
    c_mhit_fr->SaveAs("../plots/"+prefix+"_fr_mhit.png");
    c_pfchiso_fr->SaveAs("../plots/"+prefix+"_fr_pfchiso.png");
    c_pfemiso_fr->SaveAs("../plots/"+prefix+"_fr_pfemiso.png");
    c_pfnhiso_fr->SaveAs("../plots/"+prefix+"_fr_pfnhiso.png");
    c_pfchisomvacut_fr->SaveAs("../plots/"+prefix+"_fr_pfchisomvacut.png");
    c_pfemisomvacut_fr->SaveAs("../plots/"+prefix+"_fr_pfemisomvacut.png");
    c_pfnhisomvacut_fr->SaveAs("../plots/"+prefix+"_fr_pfnhisomvacut.png");
    c_ecaliso_fr->SaveAs("../plots/"+prefix+"_fr_ecaliso.png");
    c_ecalisoval_fr->SaveAs("../plots/"+prefix+"_fr_ecalisoval.png");
    c_hcaliso_fr->SaveAs("../plots/"+prefix+"_fr_hcaliso.png");
    c_trackiso_fr->SaveAs("../plots/"+prefix+"_fr_trackiso.png");

    //
    // tidy up (a bit)
    //

    f->Close();

}


