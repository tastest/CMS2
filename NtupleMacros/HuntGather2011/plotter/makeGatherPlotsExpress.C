#include "cuts.h"

void makeGatherPlotsExpress(TString prefix, const std::vector<BabySample*> &babyVector, const float &goodruns_lumi, const float &est_newruns_lumi)
{

    std::cout << "Making express plots...\n";

    //
    // Apply preselection cuts to the samples in the baby vector
    // These preselection cuts will apply to all plots!
    //

    PreselectBabies(babyVector, base_dilep);

    //
    // Make the plots
    //


    //
    // Z region
    //

    TCut express_ee ("express_ee", base_dilep+ee_dilep);
    TCut express_mm ("express_mm", base_dilep+mm_dilep);
    TCut express_eemm ("express_eemm", base_dilep+sf_dilep);
//    TCut express_zjet("express_zjet", "acos(cos(dilphi - jet1phi)) < 0.25 && njets == 1");
    TCut express_zjet("express_zjet", "njets == 1");
    TCut express_jetbal("express_jetbal", base_dilep+express_zjet);

    DrawAll("mass", prefix+"_express_mass_ee", express_ee, est_newruns_lumi, 50,0., 200., 0, babyVector);
    DrawAll("mass", prefix+"_express_mass_mm", express_mm, est_newruns_lumi, 50,0., 200., 0, babyVector);
    DrawAll("pfmet", prefix+"_express_pfmet_ee", express_ee, est_newruns_lumi, 50,0., 200., 0, babyVector);
    DrawAll("pfmet", prefix+"_express_pfmet_mm", express_mm, est_newruns_lumi, 50,0., 200., 0, babyVector);
    DrawAll("nvtx", prefix+"_express_nvtx_eemm", express_eemm, est_newruns_lumi, 20, -0.5, 19.5, 0, babyVector);
    DrawAll("jet1pt", prefix+"_express_jet1pt_eemm", express_eemm, est_newruns_lumi, 50, 0.0, 200., 0, babyVector);
   
 
    TCut var_jetbal("Z jet pT balance", "(dilpt-jet1pt)/dilpt");
    DrawAll(var_jetbal, prefix+"_express_jetbal_eemm", express_jetbal, est_newruns_lumi, 50, -5.0, 5.0, 0, babyVector);

    //
    // Top region
    //

    TCut express_ttof("express_ttof", inclusivenonz_dilep1+base_dilep+os_dilep+of_dilep);
    TCut express_ttofjj("express_ttofjj", express_ttof+"njets>=2");

    DrawAll("njets", prefix+"_express_njets_ttof", express_ttof, est_newruns_lumi, 10, -0.5, 9.5, 0, babyVector);
    DrawAll("ntchelbtags", prefix+"_express_nbjets_ttof", express_ttof, est_newruns_lumi, 10, -0.5, 9.5, 0, babyVector);
    DrawAll("pfmet", prefix+"_express_pfmet_ttofjj", express_ttofjj, est_newruns_lumi, 20, 0.0, 200., 0, babyVector);

}

