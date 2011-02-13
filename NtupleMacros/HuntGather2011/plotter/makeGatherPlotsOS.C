#include "cuts.h"

void makeGatherPlotsOS(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making OS plots...\n";

    //
    // Define specific cuts for OS plots
    //

    TCut osanal_dilep   ("osanal_dilep"   ,inclusivenonz_dilep1+base_dilep+os_dilep+"tcmet>50.&&sumjetpt>100.");
    TCut osanal_of_dilep("osanal_of_dilep",osanal_dilep+of_dilep);
    TCut osanal_sf_dilep("osanal_sf_dilep",osanal_dilep+sf_dilep);

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, osanal_dilep);

    //
    // Make the plots
    //

    DrawAll("mass", "os_of_mass", osanal_of_dilep, luminosity, 40, 0., 500., 0, babyVector);
    DrawAll("mass", "os_sf_mass", osanal_sf_dilep, luminosity,40, 0., 500.,0, babyVector);
    DrawAll("sumjetpt", "os_sumjetpt_int", osanal_dilep, luminosity,40, 0.,800., 1, babyVector);
    DrawAll("tcmet", "os_tcmet_int", osanal_dilep, luminosity, 40, 0.,300., 1, babyVector);

}

