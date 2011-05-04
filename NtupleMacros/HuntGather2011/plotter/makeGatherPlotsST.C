#include "cuts.h"

void makeGatherPlotsST(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making ST plots...\n";

    //
    // Define specific cuts for ST plots
    //

    TCut base_tcmeffgt400_dilep  ("base_tcmeffgt400_dilep",base_dilep+"tcmeff>400.");
    TCut base_tcmetgt50_dilep    ("base_tcmetgt50_dilep",base_dilep+"tcmet>50.");
    TCut base_sumjetptgt200_dilep("base_sumjetptgt200_dilep",base_dilep+"sumjetptSS>200.");
    TCut base_dilptgt100_dilep   ("base_dilptgt100_dilep",base_dilep+"dilpt>100.");

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep);

    //
    // Make the plots
    //

    DrawAll("njetsSS", prefix+"_meff_njetsclean",base_dilep,luminosity,10,-0.5,9.5,0, babyVector);
    DrawAll("njetsSS", prefix+"_meff_tcmeffgt400_njetsclean",base_tcmeffgt400_dilep,luminosity,10,-0.5,9.5,0, babyVector);
    DrawAll("tcmeff", prefix+"_meff_tcmeff_int",base_dilep,luminosity,40,0.,1000.,1, babyVector);
    DrawAll("tcmeff", prefix+"_meff_tcmetgt50_tcmeff_int",base_tcmetgt50_dilep,luminosity,40,0.,1000.,1, babyVector);
    DrawAll("tcmeff", prefix+"_meff_sumjetptgt200_tcmeff_int",base_sumjetptgt200_dilep,luminosity,40,0.,1000.,1, babyVector);
    DrawAll("tcmeff", prefix+"_meff_dilptgt100_tcmeff_int",base_dilptgt100_dilep,luminosity,40,0.,1000.,1, babyVector);

}

