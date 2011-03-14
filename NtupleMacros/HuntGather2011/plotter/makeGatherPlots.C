//#include "cuts.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TIterator.h"

#include "SampleType.h"
#include "cuts.h"

#include <vector>

void makeGatherPlots(TString base, bool debug = false) {

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

    //
    // define samples
    //

    TCut cut_notau("cut_notau", "ngentaus==0");
    TCut cut_tau("cut_tau", "ngentaus==2");

    float k_ww = 1.0;
    float k_wz = 1.0;
    float k_zz = 1.0;    
    float k_dy = 1.0;
    float k_vgammajets = 1.0;
    float k_ttbar = 1.0;
    float k_tw = 1.0;
    float k_wjets = 1.0;

    BabySample *bs_dilep_wz  = new BabySample("wz", "mc", 
            base+"/mc/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
            "", k_wz, BACKGROUND, kGray, 1001);

    BabySample *bs_dilep_zz  = new BabySample("zz", "mc", 
            base+"/mc/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
            "", k_zz, BACKGROUND, 10, 1001);

    BabySample *bs_dilep_dyeemm  = new BabySample("dyeemm", "mc", 
            base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
            cut_notau, k_dy, BACKGROUND, kAzure-2, 1001);
    bs_dilep_dyeemm->add(base+"/mc/DYToEE_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/dilep2010_ZMassLessThan50Skim/baby_gather.root");
    bs_dilep_dyeemm->add(base+"/mc/DYToMuMu_M-20_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");

    BabySample *bs_dilep_dytt  = new BabySample("dytt", "mc",
            base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root",
            cut_tau, k_dy, BACKGROUND, kCyan, 1001);
    bs_dilep_dytt->add(base+"/mc/DYToTauTau_M-20_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/diLep2010_ZMassLessThan50Skim/baby_gather.root");

    BabySample *bs_dilep_vgammajets  = new BabySample("vgammajets", "mc", 
            base+"/mc/PhotonVJets_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
            "", k_vgammajets, BACKGROUND, kOrange-3, 1001);

    BabySample *bs_dilep_ttbar  = new BabySample("ttbar", "mc", 
            base+"/mc/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
            "", k_ttbar, BACKGROUND, kRed+1, 1001);

    BabySample *bs_dilep_tw  = new BabySample("tw", "mc",
            base+"/mc/TToBLNu_TuneZ2_tW-channel_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root",
            "", k_tw, BACKGROUND, kMagenta, 1001);

    BabySample *bs_dilep_wjets  = new BabySample("wjets", "mc", 
            base+"/mc/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", 
            "", k_dy, BACKGROUND, kGreen-3, 1001);

    BabySample *bs_dilep_ww  = new BabySample("ww", "mc",
            base+"/mc/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/baby_gather.root",
            "", k_ww, BACKGROUND, kGray+1, 1001);

    //
    // BSM
    //

    float k_hww160 = 1.0;
    float k_hww130 = 1.0;
    float k_hww200 = 1.0;
    float k_lm0 = 1.0;

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc", 
            base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", 
            "", k_hww160, SIGNAL, kBlack, kSolid);
    BabySample *bs_dilep_hww130  = new BabySample("hww130", "mc",
            base+"/mc/GluGluToHToWWTo2L2Nu_M-130_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root",
            "", k_hww130, SIGNAL, kGray, kSolid);
    BabySample *bs_dilep_hww200  = new BabySample("hww200", "mc",
            base+"/mc/GluGluToHToWWTo2L2Nu_M-200_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root",
            "", k_hww200, SIGNAL, kCyan, kSolid);

    BabySample *bs_dilep_lm0  = new BabySample("LM0", "mc",
            base+"/mc/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/baby_gather.root",
            "", k_lm0, SIGNAL, kGray, kSolid);

    //
    // Data
    //

    // set up the good run list cut
    int brun = min_run();
    int bls  = min_run_min_lumi();
    int erun = max_run();
    int els  = max_run_max_lumi();
    TCut c_goodrunplus(Form("(((run>%i&&run<%i)||(run==%i&&ls>=%i)||(run==%i&&ls<=%i))&&goodrun(run,ls))||(run>%i||(run==%i&&ls>%i))",
                brun, erun, brun, bls, erun, els, erun, erun, els));

    // set up the duplicate removal cut and
    // set up the preselection cut for data
    TCut c_notduplicate = ("! is_duplicate(run,evt,ls,pt1,pt2)");
    TCut c_datapresel = c_goodrunplus + c_notduplicate;

    BabySample *bs_data = new BabySample("data", "data", 
            base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather.root",
            c_datapresel, 1.0, DATA);
    bs_data->add(base+"/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root");

    //
    // Define the mixtures of signals, background 
    // and data that can be plotted
    //

    std::vector<BabySample*> babyVectorTP;
    babyVectorTP.push_back(bs_data);
    babyVectorTP.push_back(bs_dilep_dyeemm);
    babyVectorTP.push_back(bs_dilep_dytt);
    babyVectorTP.push_back(bs_dilep_ttbar);
    babyVectorTP.push_back(bs_dilep_wjets);

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

    std::vector<BabySample*> babyVectorHiggs;
    babyVectorHiggs.push_back(bs_dilep_ww);
    babyVectorHiggs.push_back(bs_dilep_wz);
    babyVectorHiggs.push_back(bs_dilep_zz);
    babyVectorHiggs.push_back(bs_dilep_dyeemm);
    babyVectorHiggs.push_back(bs_dilep_dytt);
    babyVectorHiggs.push_back(bs_dilep_vgammajets);
    babyVectorHiggs.push_back(bs_dilep_ttbar);
    babyVectorHiggs.push_back(bs_dilep_tw);
    babyVectorHiggs.push_back(bs_dilep_wjets);
    babyVectorHiggs.push_back(bs_data);
    babyVectorHiggs.push_back(bs_dilep_hww160);
    babyVectorHiggs.push_back(bs_dilep_hww130);
    babyVectorHiggs.push_back(bs_dilep_hww200);

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
    babyVectorSusy.push_back(bs_data);
    babyVectorSusy.push_back(bs_dilep_lm0);

    //
    // Luminosity determination
    //
    const char *goodrunlist = "../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt";
    float goodruns_lumi = 35.0;

    std::cout << "[The Gathering] Determining luminosity" << std::endl;
    std::cout << "[The Gathering] " << goodrunlist << std::endl;
    set_goodrun_file(goodrunlist);
    float est_lumi = GetIntLumi(bs_data, goodruns_lumi);
    float est_newruns_lumi = est_lumi - goodruns_lumi;
    std::cout << "[The Gathering] Estimated L = " << est_lumi << std::endl;
    std::cout << std::endl;

    //
    // Make the plots
    //

    if (debug) {
        //makeGatherPlotsElectrons(babyVectorTP, est_lumi);
        //makeGatherPlotsMuons(babyVectorTP, est_lumi);
        makeGatherTriggerMonitor(babyVectorSM, est_lumi);
    }
    else {
        makeGatherPlotsValidation(babyVectorSM, goodruns_lumi, est_newruns_lumi);
        makeGatherTriggerMonitor(babyVectorSM, est_lumi);
        makeGatherPlotsHiggs(babyVectorHiggs, est_lumi);
        makeGatherPlotsOS(babyVectorSusy, est_lumi);
        makeGatherPlotsZMet(babyVectorSusy, est_lumi);
        makeGatherPlotsSS(babyVectorSusy, est_lumi);
        makeGatherPlotsST(babyVectorSusy, est_lumi);
        makeGatherPlotsExotica(babyVectorSM, est_lumi);
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
        c1->SetLogy(1);
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

}
