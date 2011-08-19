#include "cuts.h"

void makeGatherPlotsValidation(TString prefix, const std::vector<BabySample*> &babyVector, const float &goodruns_lumi, const float &est_newruns_lumi)
{

    std::cout << "Making validation plots...\n";

    //
    // Define specific cuts for validation plots
    //

    unsigned int lastgoodrun = max_run();
    unsigned int lastgoodlumi = max_run_max_lumi();

    TCut validation_ee ("validation_ee", base_dilep+ee_dilep);
    TCut validation_mm ("validation_mm", base_dilep+mm_dilep);
    TCut validation_newrun("validation_newrun", Form("!isdata||(run > %i || (run == %i && ls > %i))", lastgoodrun, lastgoodrun, lastgoodlumi));
    TCut validation_goodrun("validation_goodrun", Form("!isdata||(run < %i || (run == %i && ls <= %i))", lastgoodrun, lastgoodrun, lastgoodlumi));

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep);

    //
    // Make the plots
    //

    // mass in new and goodruns
    //DrawAll("mass", prefix+"_validation_mass_goodruns_ee", validation_ee+validation_goodrun, goodruns_lumi, 50,0., 200., 0, babyVector);
    //DrawAll("mass", prefix+"_validation_mass_newruns_ee", validation_ee+validation_newrun, est_newruns_lumi, 50,0., 200., 0, babyVector);
    //DrawAll("mass", prefix+"_validation_mass_goodruns_mm", validation_mm+validation_goodrun, goodruns_lumi, 50,0., 200., 0, babyVector);
    //DrawAll("mass", prefix+"_validation_mass_newruns_mm", validation_mm+validation_newrun, est_newruns_lumi, 50,0., 200., 0, babyVector);

    DrawAll("mass", prefix+"_validation_mass_newruns_ee", validation_ee, est_newruns_lumi, 50,0., 200., 0, babyVector);
    DrawAll("mass", prefix+"_validation_mass_newruns_mm", validation_mm, est_newruns_lumi, 50,0., 200., 0, babyVector);

    // hyps with ttbar and dy like selection
    DrawAll("hyp_type", prefix+"_validation_hyp_type_top", dileptonictopv4_dilep, goodruns_lumi+est_newruns_lumi, 7, -0.5, 7.5, 0, babyVector);
    DrawAll("hyp_type", prefix+"_validation_hyp_type_base", base_dilep, goodruns_lumi+est_newruns_lumi, 7, -0.5, 6.5, 0, babyVector);

}

