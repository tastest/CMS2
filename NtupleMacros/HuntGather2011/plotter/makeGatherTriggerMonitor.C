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

    Int_t run_start = 146000;
    Int_t run_end = 150000;
    Int_t n_run_bins = (run_end - run_start)/4;

    // event selections
    TCut cut_trg_zee("trg_zee", inclusivez_dilep+ee_dilep);
    TCut cut_trg_zmm("trg_zmm", inclusivez_dilep+mm_dilep);
    TCut cut_trg_ttof("trg_ttof", base_dilep+"(abs(eormu1)!=abs(eormu2))");

    // trigger selections
    TCut cut_trg_single_e("trg_single_e", "trg_single_e == 1");
    TCut cut_trg_single_mu("trg_single_mu", "trg_single_mu == 1");
    TCut cut_trg_double_e("trg_double_e", "trg_double_e == 1");
    TCut cut_trg_double_mu("trg_double_mu", "trg_double_mu == 1");
    TCut cut_trg_cross_emu("trg_cross_emu", "trg_cross_emu == 1");

    TriggerMonitor("trigger_double_e", cut_trg_zee+cut_trg_single_e, cut_trg_double_e, luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor("trigger_double_mu", cut_trg_zmm+cut_trg_single_mu, cut_trg_double_mu, luminosity, n_run_bins, run_start, run_end, false, data);

    DrawAll("ntchelbtags", "trigger_cross_emu_ntchelbtags", cut_trg_ttof, luminosity, 5, -0.5, 4.5, 0, babyVector);
    TriggerMonitor("trigger_cross_emu", cut_trg_ttof+"ntchelbtags>0", cut_trg_cross_emu, luminosity, n_run_bins, run_start, run_end, false, data);

}

