#include "cuts.h"

void makeGatherPlotsZMet(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Z+Met plots...\n";

    //
    // Define specific cuts for ZMet plots
    //

    TCut zmet_os_0j_dilep("zmet_os_0j_dilep",inclusivez_dilep+os_dilep+"njets==0");
    TCut zmet_os_1j_dilep("zmet_os_1j_dilep",inclusivez_dilep+os_dilep+"njets==1");
    TCut zmet_os_2j_dilep("zmet_os_2j_dilep",inclusivez_dilep+os_dilep+"njets>=2");
    TCut zmet_sumjetptgt100_os_ee_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+ee_dilep+"sumjetpt>100.");
    TCut zmet_sumjetptgt100_os_mm_dilep("zmet_sumjetptgt100_os_ee_dilep",inclusivez_dilep+os_dilep+mm_dilep+"sumjetpt>100.");

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, inclusivez_dilep+os_dilep);

    //
    // Make the plots
    //

    DrawAll("tcmet","zmet_os_0j_tcmet_int",zmet_os_0j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("tcmet","zmet_os_1j_tcmet_int",zmet_os_1j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("tcmet","zmet_os_2j_tcmet_int",zmet_os_2j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("pfmet","zmet_os_0j_pfmet_int",zmet_os_0j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("pfmet","zmet_os_1j_pfmet_int",zmet_os_1j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("pfmet","zmet_os_2j_pfmet_int",zmet_os_2j_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("tcmet","zmet_sumjetptgt100_os_ee_tcmet_int",zmet_sumjetptgt100_os_ee_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("pfmet","zmet_sumjetptgt100_os_ee_pfmet_int",zmet_sumjetptgt100_os_ee_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("tcmet","zmet_sumjetptgt100_os_mm_tcmet_int",zmet_sumjetptgt100_os_mm_dilep,luminosity,40,0.,300.,1, babyVector);
    DrawAll("pfmet","zmet_sumjetptgt100_os_mm_pfmet_int",zmet_sumjetptgt100_os_mm_dilep,luminosity,40,0.,300.,1, babyVector);

}

