#include "cuts.h"

void makeGatherMETMonitor(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making MET Monitor...\n";

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
    Int_t n_run_bins = (run_end - run_start);

    //
    // Make the plots
    //

    TriggerMonitor(prefix+"_tcmetmon", "run", inclusivez_dilep, TCut("metmon_tcmet", "tcmet > 40"), luminosity, n_run_bins, run_start, run_end, false, data);
    TriggerMonitor(prefix+"_pfmetmon", "run", inclusivez_dilep, TCut("metmon_pfmet", "pfmet > 40"), luminosity, n_run_bins, run_start, run_end, false, data);

}

