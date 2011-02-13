#include "cuts.h"

void makeGatherPlotsSS(const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making SS plots...\n";

    //
    // Define specific cuts for OS plots
    //

    TCut ssanal_lepid1 = "(abs(eormu1)==11&&e1_ctfCharge==e1_scCharge&&e1_ctfCharge==e1_gsfCharge&&e1_vbtf70) || abs(eormu1)==13";
    TCut ssanal_lepid2 = "(abs(eormu2)==11&&e2_ctfCharge==e2_scCharge&&e2_ctfCharge==e2_gsfCharge&&e2_vbtf70) || abs(eormu2)==13";
    TCut ssanal_lepid  = ssanal_lepid1+ssanal_lepid2;
    TCut ssanal_dilep   ("ssanal_dilep"   ,base_dilep+ss_dilep+ssanal_lepid);
    TCut ssanal_ee_dilep("ssanal_ee_dilep",ssanal_dilep+ee_dilep);
    TCut ssanal_mm_dilep("ssanal_mm_dilep",ssanal_dilep+mm_dilep);

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, ssanal_dilep);

    //
    // Make the plots
    //

    DrawAll("mass","ss_ee_mass",ssanal_ee_dilep,luminosity,40,0.,500.,0, babyVector);
    DrawAll("mass","ss_mm_mass",ssanal_mm_dilep,luminosity,40,0.,500.,0, babyVector);
    DrawAll("sumjetpt","ss_sumjetpt_int",ssanal_dilep,luminosity,40,0.,800.,1, babyVector);
    DrawAll("tcmet","ss_tcmet_int",ssanal_dilep,luminosity,40,0.,300.,1, babyVector);

}

