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

    float k_ww = 1.0;
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

    BabySample *bs_dilep_ww  = new BabySample("ww", "mc",
        base+"/mc/WWTo2L2Nu_TuneZ2_7TeV-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-14/baby_gather.root",
        "", k_ww, BACKGROUND, kGray+1, 1001);

    //
    // BSM
    //

    float k_hww160 = 1.0;
    float k_lm0 = 1.0;

    BabySample *bs_dilep_hww160  = new BabySample("hww160", "mc", 
        base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root", 
        "", k_hww160, SIGNAL, kBlack, kDashed);
    BabySample *bs_dilep_hww160x50  = new BabySample("hww160x50", "mc",
        base+"/mc/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6_Fall10-E7TeV_ProbDist_2010Data_BX156_START38_V12-v1/V03-06-18/baby_gather.root",
        "", k_hww160*50.0, SIGNAL, kBlack, kSolid);

    BabySample *bs_dilep_lm0  = new BabySample("LM0", "mc",
        base+"/mc/LM0_SUSY_sftsht_7TeV-pythia6_Fall10-START38_V12-v1/V03-06-18/baby_gather.root",
        "", k_lm0, SIGNAL, kGray, kSolid);

    //
    // Data
    //

    BabySample *bs_data_el2010b = new BabySample("data", "data", 
        base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple*.root",
        "", 1.0, DATA);

    BabySample *bs_data_mu2010b = new BabySample("data", "data",
        base+"/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root",
        "", 1.0, DATA);


    //
    // Luminosity determination
    //

    TChain *chain_all_data = new TChain("tree");
    chain_all_data->Add(base+"/data/Electron_Run2010B-Nov4ReReco_v1_RECO/V03-06-16/diLepPt1020Skim/baby_gather_skimmed_ntuple_*.root");
    chain_all_data->Add(base+"/data/Mu_Run2010B-Nov4ReReco_v1_RECO/V03-06-17/diLepPt1020Skim/baby_gather.root");
    const char *goodrunlist = "../runlists/Cert_TopNov5_Merged_135821-149442_allPVT.txt";
    float goodruns_lumi = 35.0;
    std::cout << "Using " << goodrunlist << " for goodruns\n";
    set_goodrun_file(goodrunlist);
    unsigned int lastgoodrun = max_run();
    unsigned int lastgoodlumi = max_run_max_lumi();
    float est_lumi = GetIntLumi(chain_all_data, goodruns_lumi);
    float est_newruns_lumi = est_lumi - goodruns_lumi;
    std::cout << "Integrated luminosity total estimate: " << est_lumi << std::endl;

    //
    // Define the mixtures of signals, background 
    // and data that can be plotted
    //

    std::vector<BabySample*> babyVectorSM;
    babyVectorSM.push_back(bs_dilep_ww);
    babyVectorSM.push_back(bs_dilep_wz);
    babyVectorSM.push_back(bs_dilep_zz);
    babyVectorSM.push_back(bs_dilep_dy);
    babyVectorSM.push_back(bs_dilep_gammajets);
    babyVectorSM.push_back(bs_dilep_ttbar);
    babyVectorSM.push_back(bs_dilep_wjets);
    babyVectorSM.push_back(bs_data_el2010b);
    babyVectorSM.push_back(bs_data_mu2010b);

    std::vector<BabySample*>babyVectorHiggs = babyVectorSM;
    babyVectorHiggs.push_back(bs_dilep_hww160);    
    babyVectorHiggs.push_back(bs_dilep_hww160x50);

    std::vector<BabySample*>babyVectorSusy = babyVectorSM;
    babyVectorHiggs.push_back(bs_dilep_lm0);

    //
    // Make the plots
    //
/*
    TCut validation_ee ("validation_ee", base_dilep+ee_dilep);
    TCut validation_mm ("validation_mm", base_dilep+mm_dilep);

    DrawAll("mass", "validation_mass_goodruns_ee", validation_ee, 
        Form("!isdata||(run < %i || (run == %i && ls <= %i))", 
        lastgoodrun, lastgoodrun, lastgoodlumi), goodruns_lumi, 50,0., 200., 0, babyVectorSM);
    DrawAll("mass", "validation_mass_newruns_ee", validation_ee, 
        Form("!isdata||(run > %i || (run == %i && ls > %i))", 
        lastgoodrun, lastgoodrun, lastgoodlumi), est_newruns_lumi, 50,0., 200., 0, babyVectorSM);
    DrawAll("mass", "validation_mass_goodruns_mm", validation_mm, 
        Form("!isdata||(run < %i || (run == %i && ls <= %i))", 
        lastgoodrun, lastgoodrun, lastgoodlumi), goodruns_lumi, 50,0., 200., 0, babyVectorSM);
    DrawAll("mass", "validation_mass_newruns_mm", validation_mm, 
        Form("!isdata||(run > %i || (run == %i && ls > %i))", 
        lastgoodrun, lastgoodrun, lastgoodlumi), est_newruns_lumi, 50,0., 200., 0, babyVectorSM);
*/
    //
    // OS PLOTS
    //
    std::cout << "Making OS plots...\n";

    TCut osanal_dilep   ("osanal_dilep"   ,base_dilep+os_dilep+"tcmet>50.&&sumjetpt>150.");
    TCut osanal_of_dilep("osanal_of_dilep",osanal_dilep+of_dilep);
    TCut osanal_sf_dilep("osanal_sf_dilep",osanal_dilep+sf_dilep);

    DrawAll("mass","os_of_mass",osanal_of_dilep,"",est_lumi,40,0.,500.,0, babyVectorSusy);
    DrawAll("mass","os_sf_mass",osanal_sf_dilep,"",est_lumi,40,0.,500.,0, babyVectorSusy);
    DrawAll("sumjetpt","os_sumjetpt_int",osanal_dilep,"",est_lumi,40,0.,800.,1, babyVectorSusy);
    DrawAll("tcmet","os_tcmet_int",osanal_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
/*
    //
    // SS PLOTS
    //

    std::cout << "Making SS plots...\n";

    TCut ssanal_lepid1 = "(abs(eormu1)==11&&e1_ctfCharge==e1_scCharge&&e1_ctfCharge==e1_gsfCharge&&e1_vbtf70) || abs(eormu1)==13";
    TCut ssanal_lepid2 = "(abs(eormu2)==11&&e2_ctfCharge==e2_scCharge&&e2_ctfCharge==e2_gsfCharge&&e2_vbtf70) || abs(eormu2)==13";
    TCut ssanal_lepid  = ssanal_lepid1+ssanal_lepid2;
    TCut ssanal_dilep   ("ssanal_dilep"   ,base_dilep+ss_dilep+ssanal_lepid);
    TCut ssanal_ee_dilep("ssanal_ee_dilep",ssanal_dilep+ee_dilep);
    TCut ssanal_mm_dilep("ssanal_mm_dilep",ssanal_dilep+mm_dilep);

    DrawAll("mass","ss_ee_mass",ssanal_ee_dilep,"",est_lumi,40,0.,500.,0, babyVectorSusy);
    DrawAll("mass","ss_mm_mass",ssanal_mm_dilep,"",est_lumi,40,0.,500.,0, babyVectorSusy);
    DrawAll("sumjetpt","ss_sumjetpt_int",ssanal_dilep,"",est_lumi,40,0.,800.,1, babyVectorSusy);
    DrawAll("tcmet","ss_tcmet_int",ssanal_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);

    //
    // Z+MET
    //

    std::cout << "Making Z+MET plots...\n";
    TCut zmet_os_0j_dilep("zmet_os_0j_dilep",inclusivez_dilep+os_dilep+"njets==0");
    TCut zmet_os_1j_dilep("zmet_os_1j_dilep",inclusivez_dilep+os_dilep+"njets==1");
    TCut zmet_os_2j_dilep("zmet_os_2j_dilep",inclusivez_dilep+os_dilep+"njets>=2");
    TCut zmet_sumjetptgt100_os_ee_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+ee_dilep+"sumjetpt>100.");
    TCut zmet_sumjetptgt100_os_mm_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+mm_dilep+"sumjetpt>100.");
    DrawAll("tcmet","zmet_os_0j_tcmet_int",zmet_os_0j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("tcmet","zmet_os_1j_tcmet_int",zmet_os_1j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("tcmet","zmet_os_2j_tcmet_int",zmet_os_2j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("pfmet","zmet_os_0j_pfmet_int",zmet_os_0j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("pfmet","zmet_os_1j_pfmet_int",zmet_os_1j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("pfmet","zmet_os_2j_pfmet_int",zmet_os_2j_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("tcmet","zmet_sumjetptgt100_os_ee_tcmet_int",zmet_sumjetptgt100_os_ee_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("pfmet","zmet_sumjetptgt100_os_ee_pfmet_int",zmet_sumjetptgt100_os_ee_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("tcmet","zmet_sumjetptgt100_os_mm_tcmet_int",zmet_sumjetptgt100_os_mm_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);
    DrawAll("pfmet","zmet_sumjetptgt100_os_mm_pfmet_int",zmet_sumjetptgt100_os_mm_dilep,"",est_lumi,40,0.,300.,1, babyVectorSusy);

    //
    // Effective Mass
    //

    std::cout << "Making other plots...\n";

    TCut base_tcmeffgt400_dilep  ("base_tcmeffgt400_dilep",base_dilep+"tcmeff>400.");
    TCut base_tcmetgt50_dilep    ("base_tcmetgt50_dilep",base_dilep+"tcmet>50.");
    TCut base_sumjetptgt200_dilep("base_sumjetptgt200_dilep",base_dilep+"sumjetpt>200.");
    TCut base_dilptgt100_dilep   ("base_dilptgt100_dilep",base_dilep+"dilpt>100.");
    DrawAll("njets","meff_njetsclean",base_dilep,"",est_lumi,10,-0.5,9.5,0, babyVectorSusy);
    DrawAll("njets","meff_tcmeffgt400_njetsclean",base_tcmeffgt400_dilep,"",est_lumi,10,-0.5,9.5,0, babyVectorSusy);
    DrawAll("tcmeff","meff_tcmeff_int",base_dilep,"",est_lumi,40,0.,1000.,1, babyVectorSusy);
    DrawAll("tcmeff","meff_tcmetgt50_tcmeff_int",base_tcmetgt50_dilep,"",est_lumi,40,0.,1000.,1, babyVectorSusy);
    DrawAll("tcmeff","meff_sumjetptgt200_tcmeff_int",base_sumjetptgt200_dilep,"",est_lumi,40,0.,1000.,1, babyVectorSusy);
    DrawAll("tcmeff","meff_dilptgt100_tcmeff_int",base_dilptgt100_dilep,"",est_lumi,40,0.,1000.,1, babyVectorSusy);

    //
    // Exotica 
    //

    std::cout << "Making additional plots...\n";

    TCut cut_z_dijets         ("z_dijets", inclusivez_dilep+"njets>=2");
    TCut cut_z_highptdijets   ("z_highptdijets", inclusivez_dilep+"njets>=2&&jet1pt>150.&&jet2pt>150.");
    DrawAll("jetmass","exotica_z_highptdijets",cut_z_highptdijets,"",est_lumi,40,0.,2000.,0, babyVectorSM);
    DrawAll("jetmass","exotica_z_dijets",cut_z_dijets,"",est_lumi,40,0.,2000.,0, babyVectorSM);
*/

    //
    // Save the plots
    //

    TSeqCollection *list = gROOT->GetListOfCanvases();
    TIterator *iter = list->MakeIterator();
    TCanvas *c1 = 0;

    while ((c1 = (TCanvas*)iter->Next())) {
        c1->Print(Form("../output/%s.png", c1->GetName()));
        c1->SetLogy(1);
        c1->Print(Form("../output/%s_log.png", c1->GetName()));
    }

    //
    // Tidy up  
    //

    delete bs_dilep_ww;
    delete bs_dilep_wz;
    delete bs_dilep_zz;
    delete bs_dilep_dy;
    delete bs_dilep_gammajets;
    delete bs_dilep_ttbar;
    delete bs_dilep_wjets;
    delete bs_dilep_hww160;
    delete bs_dilep_hww160x50;
    delete bs_dilep_lm0;
    delete bs_data_el2010b;
    delete bs_data_mu2010b;

}
