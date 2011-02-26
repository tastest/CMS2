#include "cuts.h"

void makeGatherTriggerMonitor(const std::vector<BabySample*> &babyVector, const BabySample *data, const float &luminosity)
{

    std::cout << "Making Trigger Monitor...\n";

    //
    // Define specific cuts for trigger monitoring
    //


    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    //PreselectBabies(babyVector, os_dilep_all);

    //
    // Make the plots
    //

    Int_t run_start = 145000;
    Int_t run_end = 150000;
    Int_t n_run_bins = (run_end - run_start)/4;

    TCut cut_trg_zee = ("trg_zee", inclusivez_dilep+ee_dilep);
    cut_trg_zee.SetName("trg_zee");
    TCut cut_trg_zmm = ("trg_zmm", inclusivez_dilep+mm_dilep);
    cut_trg_zmm.SetName("trg_zmm");

    TCut cut_trg_single_e = ("trg_single_e", "trg_single_e == 1");
    cut_trg_single_e.SetName("trg_single_e");

    TCut cut_trg_double_e = ("trg_double_e", "trg_double_e == 1");
    cut_trg_double_e.SetName("trg_double_e");

    TCut cut_trg_double_mu = ("trg_double_mu", "trg_double_mu == 1");
    cut_trg_double_mu.SetName("trg_double_mu");

    TriggerMonitor("trigger_double_e", cut_trg_zee, cut_trg_double_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_double_mu", cut_trg_zmm, cut_trg_double_mu, luminosity, n_run_bins, run_start, run_end, false, data);

}

