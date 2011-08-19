#include "cuts.h"

void makeGatherPlotsST(TString prefix, const std::vector<BabySample*> &babyVector, const float &luminosity)
{

    std::cout << "Making ST plots...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep+"tcmeff>250");

    //
    // Make the plots
    //

    DrawAll("tcmeff", prefix+"_tcmeff_all",TCut("base_dilep", base_dilep),luminosity,40,250.,2250., false, babyVector);
    DrawAll("njetsSS", prefix+"_njets_tcmeff250_all",TCut("base_dilep", base_dilep+"tcmeff>250"),luminosity,10,-0.5,9.5,0, babyVector);

    DrawAll("tcmeff", prefix+"_tcmeff_met50_all",TCut("tcmet>50", base_dilep + "tcmet>50"),luminosity,40,250.,2250., false, babyVector);
    DrawAll("tcmeff", prefix+"_tcmeff_met200_all",TCut("tcmet>200", base_dilep + "tcmet>200"),luminosity,40,250.,2250., false, babyVector);

    DrawAll("tcmeff", prefix+"_tcmeff_lep50_all",TCut("dilepton pt > 50", base_dilep + "dilpt>50"),luminosity,40,250.,2250., false, babyVector);
    DrawAll("tcmeff", prefix+"_tcmeff_lep200_all",TCut("dilepton pt > 200", base_dilep + "dilpt>200"),luminosity,40,250.,2250., false, babyVector);

    DrawAll("tcmeff", prefix+"_tcmeff_jet100_all", TCut("sumjetpt > 100", base_dilep + "sumjetptSS>100"),luminosity,40,250.,2250., false, babyVector);
    DrawAll("tcmeff", prefix+"_tcmeff_jet400_all", TCut("sumjetpt > 400", base_dilep + "sumjetptSS>400"),luminosity,40,250.,2250., false, babyVector);

}

