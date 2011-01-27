//#include "cuts.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TIterator.h"

#include "SampleType.h"
#include "cuts.h"

#include <vector>

void makeGatherPlots() {

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

    //
    // define samples
    //

    float k_wz = 1.0;
    float k_zz = 1.0;    
    float k_dy = 1.0;
    float k_gammajets = 1.0;
    float k_ttbar = 1.0;
    float k_wjets = 1.0;

    TString base = "/nfs-3/userdata/cms2/gather/";

    BabySample *bs_dilep_wz  = new BabySample("wz", "mc", 
        base+"/mc/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
        "", k_wz, BACKGROUND, kGray, 1001);

    BabySample *bs_dilep_zz  = new BabySample("zz", "mc", 
        base+"/mc/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
        "", k_zz, BACKGROUND, 10, 1001);

    BabySample *bs_dilep_dy  = new BabySample("dy", "mc", 
        base+"/mc/DYJetsToLL_TuneD6T_M-50_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
        "", k_dy, BACKGROUND, kAzure-2, 1001);

    BabySample *bs_dilep_gammajets  = new BabySample("gammajets", "mc", 
        base+"/mc/PhotonVJets_7TeV-madgraph_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
        "", k_gammajets, BACKGROUND, kOrange-3, 1001);

    BabySample *bs_dilep_ttbar  = new BabySample("ttbar", "mc", 
        base+"/mc/TTJets_TuneD6T_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", 
        "", k_ttbar, BACKGROUND, kRed+1, 1001);

    BabySample *bs_dilep_wjets  = new BabySample("wjets", "mc", 
        base+"/mc/WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", 
        "", k_dy, BACKGROUND, kGreen-3, 1001);

    //
    // BSM
    //

    float k_hww160 = 1.0;

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc", 
        base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", 
        "", k_hww160, SIGNAL, kBlack, kDashed);
    BabySample *bs_dilep_hww160x50  = new BabySample("hww160x50", "mc",
        base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root",
        "", k_hww160*50.0, SIGNAL, kBlack, kSolid);

    //
    // Data
    //

    BabySample *bs_data = new BabySample("electrons", "data", 
        base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple*.root",
        "", 1.0, DATA);

    //
    // Luminosity determination
    //

    TChain *chain_all_data = new TChain("tree");
    chain_all_data->Add(base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple_*.root");

    const char *goodrunlist = "../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt";
    float goodrunlumi = 35.0;
    std::cout << "Using " << goodrunlist << " for goodruns\n";
    set_goodrun_file(goodrunlist);
    //unsigned int lastgoodrun = max_run();
    //unsigned int lastgoodlumi = max_run_max_lumi();
    float est_lumi = 35.0; //GetIntLumi(chain_all_data, goodrunlumi);
    std::cout << "Integrated luminosity: " << est_lumi << std::endl;

    //
    // Define the mixtures of signals, background 
    // and data that can be plotted
    //

    std::vector<BabySample*> babyVector;
    babyVector.push_back(bs_dilep_wz);
    babyVector.push_back(bs_dilep_zz);
    babyVector.push_back(bs_dilep_dy);
    babyVector.push_back(bs_dilep_gammajets);
    babyVector.push_back(bs_dilep_ttbar);
    babyVector.push_back(bs_dilep_wjets);
    babyVector.push_back(bs_dilep_hww160);
    babyVector.push_back(bs_dilep_hww160x50);
    babyVector.push_back(bs_data);    

    //
    // Make the plots
    //

    TCut validation_ee ("validation_ee", base_dilep+ee_dilep);

    TCanvas *c1 = DrawAll("mass","test", validation_ee, "", est_lumi, 40, 0., 500., false, babyVector);
    c1->Draw();

    //
    // Save the plots
    //

    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {
        c1->Print(Form("%s.png", c1->GetName()));
        c1->SetLogy(1);
        c1->Print(Form("%s_log.png", c1->GetName()));
    }

    //
    // Tidy up  
    //

    delete bs_dilep_wz;
    delete bs_dilep_zz;
    delete bs_dilep_dy;
    delete bs_dilep_gammajets;
    delete bs_dilep_ttbar;
    delete bs_dilep_wjets;
    delete bs_dilep_hww160;
    delete bs_dilep_hww160x50;
    delete bs_data;

}
