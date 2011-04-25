//#include "cuts.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TIterator.h"

#include "SampleType.h"
#include "cuts.h"

#include <vector>

enum Mode {
    EXPRESS,
    PROMPT,
    DEBUG
};

void makeGatherPlots(TString base, Mode mode = PROMPT) {

    gROOT->ProcessLine(".L tdrstyle.C");
    gROOT->ProcessLine("setTDRStyle()");

    gSystem->Load("libTree.so");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libEG.so");
    gSystem->Load("libMathCore.so");
    gROOT->ProcessLine(".L ../libs/libHuntGather2011Plotter.so");
    gROOT->ProcessLine(".L ../libs/libCMS2NtupleMacrosCORE.so");
    gROOT->ProcessLine(".L ../libs/libCMS2NtupleMacrosTools.so");
    gROOT->ProcessLine(".L ../libs/libHuntGather2011Babymaker.so");

    gROOT->ProcessLine(".L makeGatherPlotsElectrons.C");
    gROOT->ProcessLine(".L makeGatherPlotsMuons.C");
    gROOT->ProcessLine(".L makeGatherTriggerMonitor.C");
    gROOT->ProcessLine(".L makeGatherPlotsValidation.C");
    gROOT->ProcessLine(".L makeGatherPlotsHiggs.C");
    gROOT->ProcessLine(".L makeGatherPlotsOS.C");
    gROOT->ProcessLine(".L makeGatherPlotsZMet.C");
    gROOT->ProcessLine(".L makeGatherPlotsSS.C");
    gROOT->ProcessLine(".L makeGatherPlotsST.C");
    gROOT->ProcessLine(".L makeGatherPlotsExotica.C");
    gROOT->ProcessLine(".L makeGatherPlotsExpress.C");

    //
    // define samples
    //

    TCut cut_notau("cut_notau", "ngentaus==0");
    TCut cut_tau("cut_tau", "ngentaus==2");
    TCut cut_dyee("cut_dyee", "ngenels==2");
    TCut cut_dymm("cut_dymm", "ngenmus==2");

    float k_ww = 1.0;
    float k_wz = 1.0;
    float k_zz = 1.0;    
    float k_dy = 1.0;
    float k_vgammajets = 1.0;
    float k_ttbar = 1.0;
    float k_tw = 1.0;
    float k_wjets = 1.0;
/*
    BabySample *bs_dilep_wz  = new BabySample("wz", "mc",                   base+"/mc/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_wz, BACKGROUND, kGray, 1001);
    BabySample *bs_dilep_zz  = new BabySample("zz", "mc",                   base+"/mc/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_zz, BACKGROUND, 10, 1001);
    BabySample *bs_dilep_dyeemm  = new BabySample("dyeemm", "mc",           base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", cut_notau, k_dy, BACKGROUND, kAzure-2, 1001);
    bs_dilep_dyeemm->add(   base+"/mc/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilep2010_ZMassLessThan50Skim/baby_gather.root");
    bs_dilep_dyeemm->add(   base+"/mc/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");
    BabySample *bs_dilep_dytt  = new BabySample("dytt", "mc",               base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", cut_tau, k_dy, BACKGROUND, kCyan, 1001);
    bs_dilep_dytt->add(base+"/mc/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");
    BabySample *bs_dilep_vgammajets  = new BabySample("vgammajets", "mc",   base+"/mc/PhotonVJets_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_vgammajets, BACKGROUND, kOrange-3, 1001);
    BabySample *bs_dilep_ttbar  = new BabySample("ttbar", "mc",             base+"/mc/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_ttbar, BACKGROUND, kRed+1, 1001);
    BabySample *bs_dilep_tw  = new BabySample("tw", "mc",                   base+"/mc/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_tw, BACKGROUND, kMagenta, 1001);
    BabySample *bs_dilep_wjets  = new BabySample("wjets", "mc",             base+"/mc/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", "", k_dy, BACKGROUND, kGreen-3, 1001);
    BabySample *bs_dilep_ww  = new BabySample("ww", "mc",                   base+"/mc/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/baby_gather.root", "", k_ww, BACKGROUND, kGray+1, 1001);

    // special for Tag and Probe
    BabySample *bs_dilep_dyee  = new BabySample("dyee", "mc",               base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", cut_dymm, k_dy, BACKGROUND, kAzure-2, 1001);
    BabySample *bs_dilep_dymm  = new BabySample("dymm", "mc",               base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", cut_dymm, k_dy, BACKGROUND, kAzure-2, 1001);
    // end special for Tag and Probe

    //
    // BSM
    //

    float k_hww160 = 1.0;
    float k_hww130 = 1.0;
    float k_hww200 = 1.0;
    float k_lm0 = 1.0;

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc",           base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", "", k_hww160, SIGNAL, kBlack, kSolid);
    BabySample *bs_dilep_hww130  = new BabySample("hww130", "mc",           base+"/mc/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", "", k_hww130, SIGNAL, kGray, kSolid);
    BabySample *bs_dilep_hww200  = new BabySample("hww200", "mc",           base+"/mc/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root",
            "", k_hww200, SIGNAL, kCyan, kSolid);

    BabySample *bs_dilep_lm0  = new BabySample("LM0", "mc",                 base+"/mc/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/baby_gather.root", "", k_lm0, SIGNAL, kGray, kSolid);
*/

    BabySample *bs_dilep_dyeemm  = new BabySample("dyeemm", "mc",           base+"/mctmp/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", cut_notau, k_dy, BACKGROUND, kAzure-2, 1001);
    //bs_dilep_dyeemm->add(   base+"/mc/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilep2010_ZMassLessThan50Skim/baby_gather.root");
    //bs_dilep_dyeemm->add(   base+"/mc/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");
    BabySample *bs_dilep_dytt  = new BabySample("dytt", "mc",               base+"/mctmp/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", cut_tau, k_dy, BACKGROUND, kCyan, 1001);
    //bs_dilep_dytt->add(base+"/mc/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");

    BabySample *bs_dilep_wz  = new BabySample("wz", "mc",                   base+"/mctmp/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_wz, BACKGROUND, kGray, 1001);
    BabySample *bs_dilep_zz  = new BabySample("zz", "mc",                   base+"/mctmp/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_zz, BACKGROUND, 10, 1001);
    BabySample *bs_dilep_vgammajets  = new BabySample("vgammajets", "mc",   base+"/mctmp/PhotonVJets_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_vgammajets, BACKGROUND, kOrange-3, 1001);
    BabySample *bs_dilep_ttbar  = new BabySample("ttbar", "mc",             base+"/mctmp/TTJets_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_ttbar, BACKGROUND, kRed+1, 1001);
    BabySample *bs_dilep_tw  = new BabySample("tw", "mc",                   base+"/mctmp/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_tw, BACKGROUND, kMagenta, 1001);
    BabySample *bs_dilep_wjets  = new BabySample("wjets", "mc",             base+"/mctmp/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_dy, BACKGROUND, kGreen-3, 1001);
    BabySample *bs_dilep_ww  = new BabySample("ww", "mc",                   base+"/mctmp/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_ww, BACKGROUND, kGray+1, 1001);
    
    // special for Tag and Probe
    BabySample *bs_dilep_dyee  = new BabySample("dyee", "mc",               base+"/mctmp/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", cut_dymm, k_dy, BACKGROUND, kAzure-2, 1001);
    BabySample *bs_dilep_dymm  = new BabySample("dymm", "mc",               base+"/mctmp/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", cut_dymm, k_dy, BACKGROUND, kAzure-2, 1001);
    // end special for Tag and Probe

    //
    // BSM
    //

    float k_hww160 = 1.0;
    float k_hww130 = 1.0;
    float k_hww200 = 1.0;
    float k_lm0 = 1.0;

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc",           base+"/mctmp/GluGluToHToWWToLNuTauNu_M-130_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_hww160, SIGNAL, kBlack, kSolid);
    BabySample *bs_dilep_hww130  = new BabySample("hww130", "mc",           base+"/mctmp/GluGluToHToWWToLNuTauNu_M-160_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_hww130, SIGNAL, kGray, kSolid);
    BabySample *bs_dilep_hww200  = new BabySample("hww200", "mc",           base+"/mctmp/GluGluToHToWWToLNuTauNu_M-200_7TeV-powheg-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root",
            "", k_hww200, SIGNAL, kCyan, kSolid);

    BabySample *bs_dilep_lm0  = new BabySample("LM0", "mc",                 base+"/mctmp/LM0_SUSY_sftsht_7TeV-pythia6_Spring11-PU_S1_START311_V1G1-v1/V04-01-01/baby_gather.root", "", k_lm0, SIGNAL, kGray, kSolid);


    //
    // Data
    //

    TCut c_remove_end2010bad("remove_end2010bad", "run <= 149294 || run >= 160325");
    TCut c_remove_2T2011runs("remove_2T2011runs", "run != 162713");
    // set up the duplicate removal cut and
    // set up the preselection cut for data
    TCut c_notduplicate("! is_duplicate(run,evt,ls,pt1,pt2)");
    TCut c_datapresel = c_notduplicate + c_remove_end2010bad + c_remove_2T2011runs;

    // ALL
    BabySample *bs_data = new BabySample("data", "data", base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root", c_datapresel, 1.0, DATA, kBlack);
    bs_data->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleMu_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    bs_data->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/MuEG_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    // add in 2010B
    bs_data->add(base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");
    bs_data->add(base+"/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root");
    // add in 2010A
    bs_data->add(base+"/data/EG_Run2010A-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");
    bs_data->add(base+"/data/Mu_Run2010A-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");

    // 2011 only
    BabySample *bs_data2011 = new BabySample("data", "data", base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root", c_datapresel, 1.0, DATA, kBlack);
    bs_data2011->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleMu_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    bs_data2011->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/MuEG_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    bs_data2011->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleElectron_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");
    bs_data2011->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleMu_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");
    bs_data2011->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/MuEG_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");

    // 2010 only
    BabySample *bs_data2010 = new BabySample("data", "data", base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root", c_datapresel, 1.0, DATA, kBlack);
    // add in 2010B
    bs_data2010->add(base+"/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root");
    // add in 2010A
    bs_data2010->add(base+"/data/EG_Run2010A-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");
    bs_data2010->add(base+"/data/Mu_Run2010A-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root");

    // 2011 Express Pre Tech Stop
    BabySample *bs_data_express = new BabySample("express", "data", base+"mc/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/baby_gather.root", c_datapresel, 1.0, DATA, kRed);

    // 2011 Express Post Tech Stop
    BabySample *bs_data_express_post_tech_stop = new BabySample("express", "data", base+"data/ExpressPhysics_Run2011A-Express-v2_FEVT/V04-01-02/baby_gather_merged_ntuple*.root", c_datapresel, 1.0, DATA, kRed);

    // 2011 Max (prompt reco + post tech stop express)
    //BabySample *bs_data2011max = new BabySample("data", "data", base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleElectron_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root", c_datapresel, 1.0, DATA, kBlack);
    //bs_data2011max->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/DoubleMu_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    //bs_data2011max->add(base+"/mctmp/CMSSW_4_1_2_patch1_V04-00-13/MuEG_Run2011A-PromptReco-v1_AOD/CMSSW_4_1_2_patch1_V04-00-13_merged/V04-00-13/baby_gather*.root");
    //bs_data2011max->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleElectron_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");
    //bs_data2011max->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/DoubleMu_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");
    //bs_data2011max->add(base+"/data/CMSSW_4_1_2_patch1_V04-01-03/MuEG_Run2011A-PromptReco-v2_AOD/CMSSW_4_1_2_patch1_V04-01-03_merged/V04-01-03/baby_gather*.root");
    //bs_data2011max->add(base+"/data/ExpressPhysics_Run2011A-Express-v2_FEVT/V04-01-02/baby_gather_merged_ntuple*.root");

    BabySample *bs_data2011max = new BabySample("data", "data", base+"mctmp/ExpressPhysicsRun2011A-Express-v1FEVT/V04-00-08/baby_gather.root", c_datapresel, 1.0, DATA, kBlack);
    bs_data2011max->add(base+"/data/ExpressPhysics_Run2011A-Express-v2_FEVT/V04-01-02/baby_gather_merged_ntuple*.root");

    //
    // Define the mixtures of signals, background 
    // and data that can be plotted
    //

    std::vector<BabySample*> babyVectorTPee;
    babyVectorTPee.push_back(bs_data_express);
    babyVectorTPee.push_back(bs_dilep_dyee);
    std::vector<BabySample*> babyVectorTPmm;
    babyVectorTPmm.push_back(bs_data_express);
    babyVectorTPmm.push_back(bs_dilep_dymm);

    //
    // Standard Model
    //

    std::vector<BabySample*> babyVectorSM;
    babyVectorSM.push_back(bs_data);
    babyVectorSM.push_back(bs_dilep_ww);
    babyVectorSM.push_back(bs_dilep_wz);
    babyVectorSM.push_back(bs_dilep_zz);
    babyVectorSM.push_back(bs_dilep_dyeemm);
    babyVectorSM.push_back(bs_dilep_dytt);
    babyVectorSM.push_back(bs_dilep_vgammajets);
    babyVectorSM.push_back(bs_dilep_ttbar);
    babyVectorSM.push_back(bs_dilep_tw);
    babyVectorSM.push_back(bs_dilep_wjets);

    std::vector<BabySample*> babyVectorSM_express;
    babyVectorSM_express.push_back(bs_data_express);
    babyVectorSM_express.push_back(bs_dilep_ww);
    babyVectorSM_express.push_back(bs_dilep_wz);
    babyVectorSM_express.push_back(bs_dilep_zz);
    babyVectorSM_express.push_back(bs_dilep_dyeemm);
    babyVectorSM_express.push_back(bs_dilep_dytt);
    babyVectorSM_express.push_back(bs_dilep_vgammajets);
    babyVectorSM_express.push_back(bs_dilep_ttbar);
    babyVectorSM_express.push_back(bs_dilep_tw);
    babyVectorSM_express.push_back(bs_dilep_wjets);

    std::vector<BabySample*> babyVectorSM2011PostTechStopExpress;
    babyVectorSM2011PostTechStopExpress.push_back(bs_data_express_post_tech_stop);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_ww);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_wz);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_zz);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_dyeemm);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_dytt);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_vgammajets);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_ttbar);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_tw);
    babyVectorSM2011PostTechStopExpress.push_back(bs_dilep_wjets);

    std::vector<BabySample*> babyVectorSM2011;
    babyVectorSM2011.push_back(bs_data2011);
    babyVectorSM2011.push_back(bs_dilep_ww);
    babyVectorSM2011.push_back(bs_dilep_wz);
    babyVectorSM2011.push_back(bs_dilep_zz);
    babyVectorSM2011.push_back(bs_dilep_dyeemm);
    babyVectorSM2011.push_back(bs_dilep_dytt);
    babyVectorSM2011.push_back(bs_dilep_vgammajets);
    babyVectorSM2011.push_back(bs_dilep_ttbar);
    babyVectorSM2011.push_back(bs_dilep_tw);
    babyVectorSM2011.push_back(bs_dilep_wjets);

    std::vector<BabySample*> babyVectorSM2011max;
    babyVectorSM2011max.push_back(bs_data2011max);
    babyVectorSM2011max.push_back(bs_dilep_ww);
    babyVectorSM2011max.push_back(bs_dilep_wz);
    babyVectorSM2011max.push_back(bs_dilep_zz);
    babyVectorSM2011max.push_back(bs_dilep_dyeemm);
    babyVectorSM2011max.push_back(bs_dilep_dytt);
    babyVectorSM2011max.push_back(bs_dilep_vgammajets);
    babyVectorSM2011max.push_back(bs_dilep_ttbar);
    babyVectorSM2011max.push_back(bs_dilep_tw);
    babyVectorSM2011max.push_back(bs_dilep_wjets);

    std::vector<BabySample*> babyVectorSM2010;
    babyVectorSM2010.push_back(bs_data2010);
    babyVectorSM2010.push_back(bs_dilep_ww);
    babyVectorSM2010.push_back(bs_dilep_wz);
    babyVectorSM2010.push_back(bs_dilep_zz);
    babyVectorSM2010.push_back(bs_dilep_dyeemm);
    babyVectorSM2010.push_back(bs_dilep_dytt);
    babyVectorSM2010.push_back(bs_dilep_vgammajets);
    babyVectorSM2010.push_back(bs_dilep_ttbar);
    babyVectorSM2010.push_back(bs_dilep_tw);
    babyVectorSM2010.push_back(bs_dilep_wjets);


    //
    // SUSY
    //

    std::vector<BabySample*> babyVectorSusy;
    babyVectorSusy.push_back(bs_dilep_ww);
    babyVectorSusy.push_back(bs_dilep_wz);
    babyVectorSusy.push_back(bs_dilep_zz);
    babyVectorSusy.push_back(bs_dilep_dyeemm);
    babyVectorSusy.push_back(bs_dilep_dytt);
    babyVectorSusy.push_back(bs_dilep_vgammajets);
    babyVectorSusy.push_back(bs_dilep_ttbar);
    babyVectorSusy.push_back(bs_dilep_tw);
    babyVectorSusy.push_back(bs_dilep_wjets);
    babyVectorSusy.push_back(bs_dilep_lm0);
    babyVectorSusy.push_back(bs_data);

    std::vector<BabySample*> babyVectorSusy2011;
    babyVectorSusy2011.push_back(bs_dilep_ww);
    babyVectorSusy2011.push_back(bs_dilep_wz);
    babyVectorSusy2011.push_back(bs_dilep_zz);
    babyVectorSusy2011.push_back(bs_dilep_dyeemm);
    babyVectorSusy2011.push_back(bs_dilep_dytt);
    babyVectorSusy2011.push_back(bs_dilep_vgammajets);
    babyVectorSusy2011.push_back(bs_dilep_ttbar);
    babyVectorSusy2011.push_back(bs_dilep_tw);
    babyVectorSusy2011.push_back(bs_dilep_wjets);
    babyVectorSusy2011.push_back(bs_dilep_lm0);
    babyVectorSusy2011.push_back(bs_data2011);

    std::vector<BabySample*> babyVectorSusy2011max;
    babyVectorSusy2011max.push_back(bs_dilep_ww);
    babyVectorSusy2011max.push_back(bs_dilep_wz);
    babyVectorSusy2011max.push_back(bs_dilep_zz);
    babyVectorSusy2011max.push_back(bs_dilep_dyeemm);
    babyVectorSusy2011max.push_back(bs_dilep_dytt);
    babyVectorSusy2011max.push_back(bs_dilep_vgammajets);
    babyVectorSusy2011max.push_back(bs_dilep_ttbar);
    babyVectorSusy2011max.push_back(bs_dilep_tw);
    babyVectorSusy2011max.push_back(bs_dilep_wjets);
    babyVectorSusy2011max.push_back(bs_dilep_lm0);
    babyVectorSusy2011max.push_back(bs_data2011max);


    //
    // Luminosity determination
    //

    const char *goodrunlist = "../runlists/top_nov22.txt";
    float goodruns_lumi = 36.1; // (I missed 2010A for the time being)
    std::cout << "[The Gathering] Determining luminosity" << std::endl;
    std::cout << "[The Gathering] Z Norm determined using: " << goodrunlist << std::endl;
    set_goodrun_file(goodrunlist);
    float zPerPb = GetZPerPb(bs_data->chain(), goodruns_lumi);

    // now set the dcs good run list for 2011A
    const char *goodrunlist = "../runlists/dcs_jmu.txt";
    set_goodrun_file(goodrunlist);
    std::cout << "[The Gathering] Good run list set: " << goodrunlist << std::endl;

    //
    // Make the plots
    //

    if (mode == EXPRESS) {
        float est_extra_lumi_express = GetAllLumi(bs_data2011max->chain(), zPerPb);
        std::cout << "[The Gathering] Estimated L in EXPRESS = " << est_extra_lumi_express << std::endl;
        makeGatherPlotsExpress("post_techstop", babyVectorSM2011max, goodruns_lumi, est_extra_lumi_express);
    }

    if (mode == DEBUG) {
        float est_extra_lumi = GetNewLumi(bs_data2011->chain(), zPerPb);
        std::cout << "[The Gathering] Estimated L in 2011 Prompt = " << est_extra_lumi << std::endl;
        //makeGatherPlotsElectrons(babyVectorTPee, est_lumi);
        //makeGatherPlotsMuons(babyVectorTP, est_lumi);
        //makeGatherTriggerMonitor(babyVectorSM, est_lumi);
        //makeGatherPlotsValidation("express", babyVectorSM_express, goodruns_lumi, est_extra_lumi);
        makeGatherPlotsExotica("run2011", babyVectorSM2011, est_extra_lumi);
        makeGatherPlotsExotica("run2010", babyVectorSM2010, goodruns_lumi);
        //makeGatherPlotsExotica("all", babyVectorSM, est_lumi);
    }

    if (mode == PROMPT) {

        float est_extra_lumi = GetAllLumi(bs_data2011->chain(), zPerPb);
        std::cout << "[The Gathering] Estimated L in 2011 Prompt = " << est_extra_lumi << std::endl;
        float est_extra_lumi_2011max = GetAllLumi(bs_data2011max->chain(), zPerPb);
        std::cout << "[The Gathering] Estimated L in 2011max = " << est_extra_lumi_2011max << std::endl;

        // validation for all
        //makeGatherPlotsValidation("all", babyVectorSM, goodruns_lumi, est_extra_lumi);
        makeGatherPlotsValidation("run2011", babyVectorSM2011, goodruns_lumi, est_extra_lumi);
        makeGatherPlotsValidation("max2011", babyVectorSM2011max, goodruns_lumi, est_extra_lumi_2011max);

        // higgs for all data and for 2010 + 2011
        //makeGatherPlotsHiggs("all", babyVectorSM, est_lumi);
        makeGatherPlotsHiggs("run2011", babyVectorSM2011, est_extra_lumi);
        makeGatherPlotsHiggs("max2011", babyVectorSM2011max, est_extra_lumi_2011max);

        // OS for all data and for 2010 + 2011
        //makeGatherPlotsOS("all", babyVectorSusy, est_lumi);
        makeGatherPlotsOS("run2011", babyVectorSusy2011, est_extra_lumi);
        makeGatherPlotsOS("max2011", babyVectorSusy2011max, est_extra_lumi_2011max);

        // ZMet for all data and for 2010 + 2011
        //makeGatherPlotsZMet("all", babyVectorSusy, est_lumi);
        makeGatherPlotsZMet("run2011", babyVectorSusy2011, est_extra_lumi);
        makeGatherPlotsZMet("max2011", babyVectorSusy2011max, est_extra_lumi_2011max);

        // SS plots for all data and for 2010 + 2011
        //makeGatherPlotsSS("all", babyVectorSusy, est_lumi);
        makeGatherPlotsSS("run2011", babyVectorSusy2011, est_extra_lumi);
        makeGatherPlotsSS("max2011", babyVectorSusy2011max, est_extra_lumi_2011max);

        // ST plots for all data and for 2010 + 2011
        //makeGatherPlotsST("all", babyVectorSusy, est_lumi);
        makeGatherPlotsST("run2011", babyVectorSusy2011, est_extra_lumi);
        makeGatherPlotsST("max2011", babyVectorSusy2011max, est_extra_lumi_2011max);
        
        // Exotica (in this context... Misc) for all data and for 2010 + 2011
        //makeGatherPlotsExotica("all", babyVectorSM, est_lumi);
        makeGatherPlotsExotica("run2011", babyVectorSM2011, est_extra_lumi);
        makeGatherPlotsExotica("max2011", babyVectorSM2011max, est_extra_lumi_2011max);

    } 

    //
    // Save the plots
    //

    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {

        // lin
        c1->Print(Form("../output/%s.png", c1->GetName()));
        c1->Print(Form("../output/%s.root", c1->GetName()));

        // log
        TPad *p_main = 0;
        if ((p_main = (TPad*)c1->FindObject("p_main"))) p_main->SetLogy(1);
        else c1->SetLogy(1);
        c1->Print(Form("../output/%s_log.png", c1->GetName()));

    }

    //
    // Tidy up  
    //

    delete bs_dilep_ww;
    delete bs_dilep_wz;
    delete bs_dilep_zz;
    delete bs_dilep_dyeemm;
    delete bs_dilep_dytt;
    delete bs_dilep_vgammajets;
    delete bs_dilep_ttbar;
    delete bs_dilep_tw;
    delete bs_dilep_wjets;
    delete bs_dilep_hww160;
    delete bs_dilep_lm0;
    delete bs_data;
    delete bs_data_express;
    delete bs_data_express_post_tech_stop;
    delete bs_data2011;
    delete bs_dilep_dyee;
    delete bs_dilep_dymm;
    delete bs_data2010;

}
