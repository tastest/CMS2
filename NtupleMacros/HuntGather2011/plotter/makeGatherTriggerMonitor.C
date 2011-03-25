#include "cuts.h"

void makeGatherTriggerMonitor(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Trigger Monitor...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep);
    const BabySample *data = babyVector[0]; // fix me

    //
    // Define stuff
    //

    Int_t run_start = 160000;
    Int_t run_end = 170000;

    // ok, so it's actually 4 runs per bin
    // but configure as needed
    Int_t n_run_bins = (run_end - run_start)/4;

    // event selections
    TCut cut_trg_zee("trg_zee", inclusivez_dilep+ee_dilep);
    TCut cut_trg_zmm("trg_zmm", inclusivez_dilep+mm_dilep);
    TCut cut_trg_ttof("trg_ttof", inclusivenonz_dilep1+base_dilep+os_dilep+of_dilep);

    // the single object triggers
    // note that trg_single_x is a mask
    // bit zero means LT passed and bit one means LL passed
    TCut cut_trg_single_e("trig_single_e", "(trg_single_e & (1<<0)) || (trg_single_e & (1<<1))");
    TCut cut_trg_double_single_e("trig_double_single_e", "(trg_single_e & (1<<0)) && (trg_single_e & (1<<1))");
    TCut cut_trg_single_mu("trig_single_mu", "(trg_single_mu & (1<<0)) || (trg_single_mu & (1<<1))");
    TCut cut_trg_double_single_mu("trig_double_single_mu", "(trg_single_mu & (1<<0)) && (trg_single_mu & (1<<1))");

    // the double object triggers
    // electrons
    TCut cut_trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1("trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", "(trg_double_e & (1<<0)) == (1<<0)");
    TCut cut_trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1("trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1", "(trg_double_e & (1<<1)) == (1<<1)");
    TCut cut_trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1("trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1", "(trg_double_e & (1<<2)) == (1<<2)");
    TCut cut_trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1("trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", "(trg_double_e & (1<<3)) == (1<<3)");

    // muons
    TCut cut_trg_DoubleMu7("trg_DoubleMu7", "(trg_double_mu & (1<<0)) == (1<<0)");

    // cross
    TCut cut_trg_Mu17_Ele8_CaloIdL_v1("trg_Mu17_Ele8_CaloIdL_v1", "(trg_cross_emu & (1<<0)) == (1<<0)");
    TCut cut_trg_Mu8_Ele17_CaloIdL_v1("trg_Mu8_Ele17_CaloIdL_v1", "(trg_cross_emu & (1<<1)) == (1<<1)");
    TCut cut_trg_Mu8_Photon20_CaloIdVT_IsoT_v2("trg_Mu8_Photon20_CaloIdVT_IsoT_v2", "(trg_cross_emu & (1<<2)) == (1<<2)");

    //
    // Make the plots
    //

    // monitor denominator distributions
/*
    DrawAll("njets", "trigger_cross_emu_njets", inclusivenonz_dilep1+base_dilep+os_dilep+of_dilep, luminosity, 5, -0.5, 4.5, 0, babyVector);
    DrawAll("nvtx", "trigger_single_nvtx", inclusivez_dilep+os_dilep+sf_dilep, luminosity, 20, -0.5, 19.5, 0, babyVector);

    // single object triggers.  Efficiency is ~ unbiassed tag and probe
    TriggerMonitor("trigger_single_e", "run","e1_scet>30&&e2_scet>30"+cut_trg_zee+cut_trg_single_e, cut_trg_double_single_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_single_mu", "run", "pt1>20&&pt2>20"+cut_trg_zmm+cut_trg_single_mu, cut_trg_double_single_mu, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_single_e_vtx", "nvtx", "e1_scet>30&&e2_scet>30"+cut_trg_zee+cut_trg_single_e, cut_trg_double_single_e, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_single_mu_vtx", "nvtx", "pt1>20&&pt2>20"+cut_trg_zmm+cut_trg_single_mu, cut_trg_double_single_mu, luminosity, 20, -0.5, 19.5, false, data);

    //
    // cross
    //

    // e denom
    // run
    TriggerMonitor("trigger_singleedenom_cross_emu_trg_Mu17_Ele8_CaloIdL_v1", "run", 
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu17_Ele8_CaloIdL_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singleedenom_cross_emu_trg_cut_trg_Mu8_Ele17_CaloIdL_v1", "run",
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu8_Ele17_CaloIdL_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singleedenom_cross_emu_trg_Mu8_Photon20_CaloIdVT_IsoT_v2", "run",
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30)) && ((abs(eormu1)==13&&pt1>10)||(abs(eormu2)==13&&pt2>10))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu8_Photon20_CaloIdVT_IsoT_v2, luminosity, n_run_bins, run_start, run_end, false, data);

    // nvtx
    TriggerMonitor("trigger_nvtx_singleedenom_cross_emu_trg_Mu17_Ele8_CaloIdL_v1", "nvtx",
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu17_Ele8_CaloIdL_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singleedenom_cross_emu_trg_cut_trg_Mu8_Ele17_CaloIdL_v1", "nvtx",
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu8_Ele17_CaloIdL_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singleedenom_cross_emu_trg_Mu8_Photon20_CaloIdVT_IsoT_v2", "nvtx",
            "((abs(eormu1)==11&&e1_scet>30)||(abs(eormu2)==11&&e2_scet>30)) && ((abs(eormu1)==13&&pt1>10)||(abs(eormu2)==13&&pt2>10))"+cut_trg_ttof+cut_trg_single_e, cut_trg_Mu8_Photon20_CaloIdVT_IsoT_v2, luminosity, 20, -0.5, 19.5, false, data);

    // mu denom
    // run
    TriggerMonitor("trigger_singlemudenom_cross_emu_trg_Mu17_Ele8_CaloIdL_v1", "run",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu17_Ele8_CaloIdL_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singlemudenom_cross_emu_trg_cut_trg_Mu8_Ele17_CaloIdL_v1", "run",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20)) && ((abs(eormu1)==11&&e1_scet>20)||(abs(eormu2)==11&&e2_scet>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu8_Ele17_CaloIdL_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singlemudenom_cross_emu_trg_Mu8_Photon20_CaloIdVT_IsoT_v2", "run",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20)) && ((abs(eormu1)==11&&e1_scet>20)||(abs(eormu2)==11&&e2_scet>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu8_Photon20_CaloIdVT_IsoT_v2, luminosity, n_run_bins, run_start, run_end, false, data);

    // nvtx
    TriggerMonitor("trigger_nvtx_singlemudenom_cross_emu_trg_Mu17_Ele8_CaloIdL_v1", "nvtx",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu17_Ele8_CaloIdL_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singlemudenom_cross_emu_trg_cut_trg_Mu8_Ele17_CaloIdL_v1", "nvtx",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20))&& ((abs(eormu1)==11&&e1_scet>20)||(abs(eormu2)==11&&e2_scet>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu8_Ele17_CaloIdL_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singlemudenom_cross_emu_trg_Mu8_Photon20_CaloIdVT_IsoT_v2", "nvtx",
            "((abs(eormu1)==13&&pt1>20)||(abs(eormu2)==13&&pt2>20))&& ((abs(eormu1)==11&&e1_scet>20)||(abs(eormu2)==11&&e2_scet>20))"+cut_trg_ttof+cut_trg_single_mu, cut_trg_Mu8_Photon20_CaloIdVT_IsoT_v2, luminosity, 20, -0.5, 19.5, false, data);

    //
    // double muons
    // 
    // run
    TriggerMonitor("trigger_singledenom_trg_DoubleMu7", "run", 
            "(pt1>20||pt2>20)" + cut_trg_zmm+cut_trg_single_mu, cut_trg_DoubleMu7, luminosity, n_run_bins, run_start, run_end, false, data);

    // nvtx
    TriggerMonitor("trigger_nvtx_singledenom_trg_DoubleMu7", "nvtx",
            "(pt1>20||pt2>20)" + cut_trg_zmm+cut_trg_single_mu, cut_trg_DoubleMu7, luminosity, 20, -0.5, 19.5, false, data);

    //
    // double electrons
    //
    // run
    TriggerMonitor("trigger_singledenom_trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", "run", 
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singledenom_trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1", "run",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singledenom_trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1", "run",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30) && TMath::Min(e1_scet,e2_scet)>20"+cut_trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singledenom_trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", "run",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1, luminosity, n_run_bins, run_start, run_end, false, data);

    // nvtx
    TriggerMonitor("trigger_nvtx_singledenom_trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", "nvtx",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singledenom_trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1", "nvtx",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_nvtx_singledenom_trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1", "nvtx",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)&& TMath::Min(e1_scet,e2_scet)>20"+cut_trg_Ele32_CaloIdL_CaloIsoVL_SC17_v1, luminosity, 20, -0.5, 19.5, false, data);
   TriggerMonitor("trigger_nvtx_singledenom_trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", "nvtx",
            cut_trg_zee+cut_trg_single_e, "(e1_scet>30||e2_scet>30)"+cut_trg_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1, luminosity, 20, -0.5, 19.5, false, data);
*/
}

