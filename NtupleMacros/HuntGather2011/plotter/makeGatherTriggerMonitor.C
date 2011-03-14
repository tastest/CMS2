#include "cuts.h"

//void makeGatherTriggerMonitor(const std::vector<BabySample*> &babyVector, const BabySample *data, const float &luminosity)
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

    Int_t run_start = 146000;
    Int_t run_end = 150000;

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
    TCut cut_trg_double_e("trg_double_e", "trg_double_e == 1");
    TCut cut_trg_double_mu("trg_double_mu", "trg_double_mu == 1");
    TCut cut_trg_cross_emu("trg_cross_emu", "trg_cross_emu == 1");

    //
    // Make the plots
    //

    // monitor denominator distributions

    DrawAll("njets", "trigger_cross_emu_njets", inclusivenonz_dilep1+base_dilep+os_dilep+of_dilep, luminosity, 5, -0.5, 4.5, 0, babyVector);

    // single object triggers.  Efficiency is ~ unbiassed tag and probe
    TriggerMonitor("trigger_single_e", "run", cut_trg_zee+cut_trg_single_e, cut_trg_double_single_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_single_mu", "run", cut_trg_zmm+cut_trg_single_mu, cut_trg_double_single_mu, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_single_e_vtx", "nvtx", cut_trg_zee+cut_trg_single_e, cut_trg_double_single_e, luminosity, 20, -0.5, 19.5, false, data);
    TriggerMonitor("trigger_single_mu_vtx", "nvtx", cut_trg_zmm+cut_trg_single_mu, cut_trg_double_single_mu, luminosity, 20, -0.5, 19.5, false, data);

    // double object triggers wrt single
    TriggerMonitor("trigger_singledenom_double_e", "run", cut_trg_zee+cut_trg_single_e, cut_trg_double_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singledenom_double_mu", "run", cut_trg_zmm+cut_trg_single_mu, cut_trg_double_mu, luminosity, n_run_bins, run_start, run_end, false, data);

    // double object triggers wrt reco... this is just fraction of events that passed a double trigger
    // however since we monitor for bin by bin changes, rather than absolute "efficiency" it's ok
    TriggerMonitor("trigger_recodenom_double_e", "run", cut_trg_zee, cut_trg_double_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_recodenom_double_mu", "run", cut_trg_zmm, cut_trg_double_mu, luminosity, n_run_bins, run_start, run_end, false, data);

    // now the e-mu wrt reco denominator and single e and single mu
    TriggerMonitor("trigger_recodenom_cross_emu", "run", cut_trg_ttof, cut_trg_cross_emu, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singleedenom_cross_emu", "run", cut_trg_ttof+cut_trg_single_e, cut_trg_cross_emu, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_singlemudenom_cross_emu", "run", cut_trg_ttof+cut_trg_single_mu, cut_trg_cross_emu, luminosity, n_run_bins, run_start, run_end, false, data);

}

