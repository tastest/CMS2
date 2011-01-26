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
    float k_hww160 = 1.0;
    float k_data = 1.0;

    BabySample *bs_dilep_wz  = new BabySample("wz", "mc", "/tas/cms2/gather/mc/WZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_wz, BACKGROUND, kAzure-2, 1001);

    BabySample *bs_dilep_zz  = new BabySample("zz", "mc", "/tas/cms2/gather/mc/ZZtoAnything_TuneZ2_7TeV-pythia6-tauola_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-17/baby_gather.root", "", k_zz, BACKGROUND, kAzure, 1001);

    //
    // BSM
    //

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc", "/tas/cms2/gather/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", "", k_hww160, SIGNAL, kGreen, 10);
    BabySample *bs_dilep_hww160x10  = new BabySample("hww160x10", "mc", "/tas/cms2/gather/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", "", k_hww160*10, SIGNAL, kGreen+2, 10);

    //
    // Data
    //

    BabySample *bs_data = new BabySample("electrons", "data", "/tas/cms2/gather/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_nt
uple_*.root", "", k_data, DATA);

    TChain *chain_all_data = new TChain("tree");
    chain_all_data->Add("/tas/cms2/gather/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple_*.root");

    const char *goodrunlist = "../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt";
    float goodrunlumi = 35.0;
    std::cout << "Using " << goodrunlist << " for goodruns\n";
    set_goodrun_file(goodrunlist);
    //unsigned int lastgoodrun = max_run();
    //unsigned int lastgoodlumi = max_run_max_lumi();

    float est_lumi = GetIntLumi(chain_all_data, goodrunlumi);
    std::cout << "Integrated luminosity: " << est_lumi << std::endl;

    BabySample *bs_data = new BabySample("electrons", "data", "/tas/cms2/gather/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple_*.root", "", k_data, DATA);

    //
    // Standard Model MC
    //

    std::vector<BabySample*> babyVector;
    babyVector.push_back(bs_dilep_wz);
    babyVector.push_back(bs_dilep_zz);
    babyVector.push_back(bs_dilep_hww160);
    babyVector.push_back(bs_dilep_hww160x10);
    babyVector.push_back(bs_data);    

    //
    // BSM MC
    //

    TCanvas *c1 = DrawAll("mass","os_of_mass", "pt1>20.0&&pt2>20.0&&iso1<0.2&&iso2<0.2", "", est_lumi, 40, 0., 500., false, babyVector);
    c1->Draw();


/*
    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {
        c1->Print(Form("%s.png", c1->GetName()));
        c1->SetLogy(1);
        c1->Print(Form("%s_log.png", c1->GetName()));
    }
*/

    delete bs_dilep_wz;
    delete bs_dilep_zz;
    delete bs_data;

}
