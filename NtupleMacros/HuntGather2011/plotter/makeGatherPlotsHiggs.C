#include "cuts.h"

void makeGatherPlotsHiggs(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making Higgs plots...\n";

    //
    // Define specific cuts for OS plots
    //

    TCut base_hwwpt("(pt1 > 20. && pt2 > 20.) || (pt2 > 20.0 && pt1 > 20.)");
    TCut base_hwwid1("(abs(eormu1) == 11 && e1_vbtf90full && e1_vbtf80) || (abs(eormu1) == 13 && mu1_muonidfull && ! mu1_cosmic)");
    TCut base_hwwid2("(abs(eormu2) == 11 && e2_vbtf90full && e2_vbtf80) || (abs(eormu2) == 13 && mu2_muonidfull && ! mu2_cosmic)");
    TCut base_hwwclean("evt_clean082010 == 1");
    TCut base_hwwmet("((abs(eormu1)==abs(eormu2) && proj_tcmet > 35) || (abs(eormu1)!=abs(eormu2) && proj_tcmet > 20))");
    TCut base_hwwnjets("njets == 0");
    TCut base_hwwnoz("(abs(eormu1)!=abs(eormu2)) || (mass < 76.0 || mass > 106.0)");
    TCut hww_incldilep("hww_incldilep", base_hwwpt+base_hwwid1+base_hwwid2+base_hwwclean+base_hwwmet+base_hwwnoz);
    TCut hww_excldilep("hww_excldilep", base_hwwpt+base_hwwid1+base_hwwid2+base_hwwclean+base_hwwmet+base_hwwnoz+base_hwwnjets);

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, hww_incldilep);

    //
    // Make the plots
    //

    DrawAll("njets",prefix+"_hww_njets", hww_incldilep, luminosity, 10, -0.5, 9.5, 0, babyVector);
    DrawAll("mass",prefix+"_hww_mass", hww_excldilep, luminosity, 20, 0.0, 400.0, 0, babyVector);

}

