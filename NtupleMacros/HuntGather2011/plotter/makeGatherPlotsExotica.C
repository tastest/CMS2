#include "cuts.h"

void makeGatherPlotsExotica(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Exotica plots...\n";

    //
    // Define specific cuts for Exotica plots
    //

    TCut cut_z_dijets         ("z_dijets", inclusivez_dilep+"njets>=2");
    TCut cut_z_highptdijets   ("z_highptdijets", inclusivez_dilep+"njets>=2&&jet1pt>150.&&jet2pt>150.");

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, inclusivez_dilep);

    //
    // Make the plots
    //

    DrawAll("jetmass",prefix+"_exotica_z_highptdijets",cut_z_highptdijets,luminosity,40,0.,2000.,0, babyVector);
    DrawAll("jetmass",prefix+"_exotica_z_dijets",cut_z_dijets,luminosity,40,0.,2000.,0, babyVector);

}

